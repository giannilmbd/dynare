function dime(objective_function, init_x, mh_bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, dr , steady_state, exo_steady_state, exo_det_steady_state)
%dime(objective_function, init_x, mh_bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, dr , steady_state, exo_steady_state, exo_det_steady_state)
% Differential-Independence Mixture Ensemble ("DIME") MCMC sampling
% as proposed in "Ensemble MCMC Sampling for Robust Bayesian Inference"
% (Gregor Boehl, 2022, SSRN No. 4250395):
%
%   https://gregorboehl.com/live/dime_mcmc_boehl.pdf
%
% INPUTS
% - objective_function      [char]     string specifying the name of the objective function (posterior kernel).
% - init_x                  [double]   p×1 vector of parameters to be estimated (initial values, not used).
% - mh_bounds               [double]   p×2 matrix defining lower and upper bounds for the parameters.
% - dataset_                [dseries]  sample
% - dataset_info            [struct]   information about the dataset
% - options_                [struct]   Dynare's options
% - M_                      [struct]   model description
% - estim_params_           [struct]   estimated parameters
% - bayestopt_              [struct]   describing the priors
% - steady_state            [double]   steady state of endogenous variables
% - exo_steady_state        [double]   steady state of exogenous variables
% - exo_det_steady_state    [double]   steady state of exogenous deterministic variables
%
% SPECIAL REQUIREMENTS
% None.

% Copyright © 2024-2026 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

    % initialize
    ndim = length(init_x);
    prop_cov = eye(ndim);
    prop_mean = zeros(1,ndim);
    naccepted = 1;
    cumlweight = -inf;
    fixPSD = 1e-16*prop_cov;

    % get options & assign problem specific default value
    opts = options_.posterior_sampler_options.current_options;
    if isfield(opts,'gamma')
        g0 = opts.gamma;
    else
        g0 = 2.38 / sqrt(2 * ndim);
    end
    nchain = opts.nchain;
    dft = opts.df_proposal_dist;
    tune = opts.tune;
    if tune > opts.niter
        error('dime: number of ensemble iterations to keep (tune=%d) exceeds number of iterations (niter=%d)!', tune, opts.niter)
    end

    % temporary workaround to deal with parallelization
    if opts.parallel
        % check if parallel is possible
        installed_toolboxes = ver;
        toolbox_installed = any(strcmp("Parallel Computing Toolbox", {installed_toolboxes.Name}));
        if ~toolbox_installed || feature('numcores') == 1
            error('dime: parallel processing option is chosen but Parallel Computing Toolbox is not installed or machine has only one core')
        end
    end

    % Set location for the simulated particles.
    SimulationFolder = CheckPath('dime', M_.dname);
    % Define prior distribution
    Prior = dprior(bayestopt_, options_.prior_trunc);
    bounds = check_bounds(Prior, mh_bounds);

    % Set function handle for the objective
    eval(sprintf('%s = @(x) %s(x, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, mh_bounds, dr , steady_state, exo_steady_state, exo_det_steady_state, []);', 'funobj', func2str(objective_function)));

    % Initialization of the sampler (draws from the prior distribution with finite logged likelihood)
    t0 = tic;
    [x, ~, lprob] = ...
        smc_samplers_initialization(funobj, 'dime', nchain, Prior, SimulationFolder, opts.niter, options_.DynareRandomStreams.seed, options_.parallel_info);
    x = ptransform(x', bounds, true);

    disp_verbose(sprintf('Estimation:dime: log-posterior standard deviation (post.std) of a multivariate normal would be %.2f.\n', sqrt(0.5*ndim)), options_.verbosity);

    % preallocate
    chains = zeros(opts.niter, nchain, ndim);
    lprobs = zeros(opts.niter, nchain);
    isplit = fix(nchain/2) + 1;

    for iter = 1:opts.niter

        for complementary_ensemble = [false true]

            % define current ensemble
            if complementary_ensemble
                idcur = 1:isplit;
                idref = (isplit+1):nchain;
            else
                idcur = (isplit+1):nchain;
                idref = 1:isplit;
            end
            cursize = length(idcur);
            refsize = nchain - cursize;

            % get differential evolution proposal
            % draw the indices of the complementary chains
            i1 = (0:cursize-1) + randi([1 cursize-1], 1, cursize);
            i2 = (0:cursize-1) + randi([1 cursize-2], 1, cursize);
            i2(i2 >= i1) = i2(i2 >= i1) + 1;

            % add small noise and calculate proposal
            f = opts.sigma*randn(cursize, 1);
            q = x(idcur,:) + g0 .* (x(idref(mod(i1, refsize) + 1),:) - x(idref(mod(i2, refsize) + 1),:)) + f;
            factors = zeros(cursize,1);

            % get AIMH proposal
            xchnge = rand(cursize,1) <= opts.aimh_prob;

            % draw alternative candidates and calculate their proposal density
            prop_cov_chol = chol(prop_cov*(dft - 2)/dft + fixPSD);
            xcand = rand_multivariate_student(prop_mean, prop_cov_chol, dft, sum(xchnge));
            lprop_old = log(multivariate_student_pdf(x(idcur(xchnge),:), prop_mean, prop_cov_chol, dft));
            lprop_new = log(multivariate_student_pdf(xcand, prop_mean, prop_cov_chol, dft));

            % update proposals and factors
            q(xchnge,:) = xcand;
            factors(xchnge) = lprop_old - lprop_new;

            % Metropolis-Hastings
            newlprob = log_prob_fun(funobj, bounds, opts.parallel, q);
            lnpdiff = factors + newlprob - lprob(idcur);
            accepted = lnpdiff > log(rand(cursize,1));
            naccepted = naccepted + sum(accepted);

            % update chains
            x(idcur(accepted),:) = q(accepted,:);
            lprob(idcur(accepted)) = newlprob(accepted);
        end

        % store
        chains(iter,:,:) = ptransform(x, bounds, false);
        lprobs(iter,:) = lprob;

        % log weight of current ensemble
        lweight = logsumexp(lprob) + log(naccepted) - log(nchain);

        % calculate stats for current ensemble
        ncov = cov(x);
        nmean = mean(x);

        % update AIMH proposal distribution
        newcumlweight = logsumexp([cumlweight lweight]);
        prop_cov = exp(cumlweight - newcumlweight) .* prop_cov + exp(lweight - newcumlweight) .* ncov;
        prop_mean = exp(cumlweight - newcumlweight) .* prop_mean + exp(lweight - newcumlweight) .* nmean;
        cumlweight = newcumlweight + log(opts.rho);

        % be informative
        tt = toc(t0);
        if iter == 1 || ~mod(iter,15)
            disp_verbose('     #iter.      post.mode     post.std    %accept     lapsed', options_.verbosity)
            disp_verbose('     ------      ---------     --------    -------     ------', options_.verbosity)
        end
        counter = sprintf('%u/%u', iter, opts.niter);
        disp_verbose(sprintf('%11s     %10.3f     %1.2e     %5.2f%%     %5.2fs', counter, max(lprob), std(lprob), 100*naccepted/nchain, tt), options_.verbosity)
        naccepted = 0;
    end

    save(sprintf('%s%schains.mat', SimulationFolder, filesep()), 'chains', 'lprobs', 'tune')
end

function lprobs = log_prob_fun(loglikefun, bounds, runs_in_parallel, x)

    dim = length(x);
    x = ptransform(x, bounds, false);
    lprobs = zeros(dim,1);
    if runs_in_parallel
        parfor i=1:dim
            lprobs(i) = -loglikefun(x(i,:)');
        end
    else
        for i=1:dim
            lprobs(i) = -loglikefun(x(i,:)');
        end
    end
end

function nbounds = check_bounds(Prior, bounds)

    % ensure that prior is continuous at bounds
    nbounds = bounds;
    prior_at_mean = Prior.density(Prior.mean);
    dim = length(Prior.mean);
    for j = 1:dim
        if ~isinf(bounds.ub(j))
            x = Prior.mean;
            x(j) = bounds.ub(j);
            if Prior.density(x - 1e-16) - prior_at_mean > log(1e-16)
                nbounds.ub(j) = nbounds.ub(j) + 1e-2;
            end
            x = Prior.mean;
            x(j) = bounds.lb(j);
            if Prior.density(x + 1e-16) - prior_at_mean > log(1e-16)
                nbounds.lb(j) = nbounds.lb(j) - 1e-2;
            end
        end
    end
end

function x = ptransform(x, bounds, con2unc)

    % parameter transformation
    lb = bounds.lb;
    ub = bounds.ub;
    dim = size(x,2);
    for j = 1:dim
        % one-sided
        if ~isinf(lb(j)) && isinf(ub(j))
            if con2unc
                x(:,j) = log(x(:,j) - lb(j));
            else
                x(:,j) = lb(j) + exp(x(:,j));
            end
        % two-sided
        elseif ~isinf(lb(j)) && ~isinf(ub(j))
            if con2unc
                y = (x(:,j) - lb(j))/(ub(j) - lb(j));
                x(:,j) = log(y ./ (1 - y));
            else
                x(:,j) = lb(j) + (ub(j) - lb(j)) ./ (1 + exp(-x(:,j)));
            end
        end
    end
end
