function mdd = dsmh(objective_function, mh_bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, oo_)
% function mdd = dsmh(objective_function, mh_bounds, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, oo_)
% Dynamic Striated Metropolis-Hastings algorithm.
% based on Waggoner/Wu/Zha (2016): "Striated Metropolis–Hastings sampler for high-dimensional models,"
% Journal of Econometrics, 192(2): 406-420, https://doi.org/10.1016/j.jeconom.2016.02.007
%
% INPUTS
%   o objective_function  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters.
%   o dataset_              data structure
%   o dataset_info          dataset info structure
%   o options_              options structure
%   o M_                    model structure
%   o estim_params_         estimated parameters structure
%   o bayestopt_            estimation options structure
%   o oo_                   outputs structure
%
% SPECIAL REQUIREMENTS
%   None.
% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code suitable to be executed in
% parallel on multi core or cluster machine (in general a 'for' cycle)
% has been removed from this function and been placed in the posterior_sampler_core.m function.
%
% The DYNARE parallel packages comprise a i) set of pairs of MATLAB functions that can be executed in
% parallel and called name_function.m and name_function_core.m and ii) a second set of functions used
% to manage the parallel computations.
%
% This function was the first function to be parallelized. Later, other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management functions.
% Copyright © 2022-2024 Dynare Team
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

opts = options_.posterior_sampler_options.current_options;

% Set location for the simulated particles.
SimulationFolder = CheckPath('dsmh', M_.dname);
%delete old stale files before creating new ones
delete_stale_file(sprintf('%s%sparticles-*.mat', SimulationFolder,filesep))

% Define prior distribution
Prior = dprior(bayestopt_, options_.prior_trunc);

% Set function handle for the objective
eval(sprintf('%s = @(x) %s(x, dataset_, dataset_info, options_, M_, estim_params_, bayestopt_, mh_bounds, oo_.dr , oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);', 'funobj', func2str(objective_function)));

lambda = exp(bsxfun(@minus,opts.H,1:1:opts.H)/(opts.H-1)*log(opts.lambda1));
c = 0.055 ;
MM = int64(opts.N*opts.G/10) ;

% Step 0: Initialization of the sampler
[param, tlogpost_iminus1, loglik] = ...
    smc_samplers_initialization(funobj, 'dsmh', opts.particles, Prior, SimulationFolder, opts.H) ;

ESS = zeros(opts.H,1) ;
zhat = 0 ;

% The DSMH starts here
dprintf('#Iter.       lambda         ESS                 c     Accept. rate       scale      resample    seconds')
for i=2:opts.H
    t0 = tic;
    % Step 1: sort the densities and compute IS weights
    [tlogpost_iminus1,loglik,param] = sort_matrices(tlogpost_iminus1, loglik,param) ;
    [tlogpost_i,weights,zhat,ESS,Omegachol] = compute_IS_weights_and_moments(param, tlogpost_iminus1, loglik, lambda, i, zhat, ESS) ;
    % Step 2: tune c_i
    [c,acpt] = tune_c(funobj, param, tlogpost_i, lambda, i, c, Omegachol, weights, mh_bounds, opts, Prior) ;
    % Step 3: Metropolis step
    [param,tlogpost_iminus1,loglik] = mutation_DSMH(funobj, param, tlogpost_i, tlogpost_iminus1, loglik, lambda, i, c, MM, Omegachol, weights, mh_bounds, opts, Prior) ;
    tt = toc(t0) ;
    dprintf('%3u          %5.4f     %9.5E         %5.4f        %5.4f        %+5.4f        %3s       %5.2f', i, lambda(i), ESS(i), c, acpt, zhat, 'no', tt)
    %save(sprintf('%s%sparticles-%u-%u.mat', SimulationFolder, filesep(), i, opts.H), 'param', 'tlogpost', 'loglik')
end

tlogpost = tlogpost_iminus1 + loglik*(lambda(end)-lambda(end-1));
weights = exp(loglik*(lambda(end)-lambda(end-1)));
weights = weights/sum(weights);
iresample = kitagawa(weights);
particles = param(:,iresample);
tlogpostkernel = tlogpost(iresample);
loglikelihood = loglik(iresample);
save(sprintf('%s%sparameters_particles_final.mat', SimulationFolder, filesep()), 'particles', 'tlogpostkernel', 'loglikelihood')

mdd = zhat;

function [tlogpost_iminus1,loglik,param] = sort_matrices(tlogpost_iminus1, loglik,param)
[~,indx_ord] = sortrows(tlogpost_iminus1);
tlogpost_iminus1 = tlogpost_iminus1(indx_ord);
param = param(:,indx_ord);
loglik = loglik(indx_ord);

function [tlogpost_i,weights,zhat,ESS,Omegachol] = compute_IS_weights_and_moments(param, tlogpost_iminus1, loglik, lambda, i, zhat, ESS)
if i==1
    tlogpost_i = tlogpost_iminus1 + loglik*lambda(i);
else
    tlogpost_i = tlogpost_iminus1 + loglik*(lambda(i)-lambda(i-1));
end
weights = exp(tlogpost_i-tlogpost_iminus1);
zhat = zhat + log(mean(weights)) ;
weights = weights/sum(weights);
ESS(i) = 1/sum(weights.^2);
% estimates of mean and variance
mu = param*weights;
z = bsxfun(@minus,param,mu);
Omega = z*diag(weights)*z';
Omegachol = chol(Omega)';

function [c,acpt] = tune_c(funobj, param, tlogpost_i, lambda, i, c, Omegachol, weights, mh_bounds, opts, Prior)
%        disp('tuning c_i...');
%        disp('Initial value =');
%        disp(c) ;
npar = size(param,1);
lower_prob = (.5*(opts.alpha0+opts.alpha1))^5;
upper_prob = (.5*(opts.alpha0+opts.alpha1))^(1/5);
stop=0 ;
outer_iter=1;
while stop==0 && outer_iter<200
    acpt = 0.0;
    indx_resmpl = kitagawa(weights,rand(1,1),opts.G);
    param0 = param(:,indx_resmpl);
    tlogpost0 = tlogpost_i(indx_resmpl);
    for j=1:opts.G
        for l=1:opts.K
            validate = 0;
            l_iter=1;
            while validate == 0 && l_iter<200
                candidate = param0(:,j) + sqrt(c)*Omegachol*randn(npar,1);
                l_iter=l_iter+1;
                if all(candidate >= mh_bounds.lb) && all(candidate <= mh_bounds.ub)
                    [tlogpostx,loglikx] = tempered_likelihood(funobj, candidate, lambda(i), Prior);
                    if isfinite(loglikx) % if returned log-density is not Inf or Nan (penalized value)
                        validate = 1;
                        if rand(1,1)<exp(tlogpostx-tlogpost0(j)) % accept
                            acpt = acpt + 1/(opts.G*opts.K);
                            param0(:,j)= candidate;
                            tlogpost0(j) = tlogpostx;
                        end
                    end
                end
            end
            if l_iter==200
                error('dsmh: Inner loop reached maximum of iterations.')
            end
        end
    end
    %           disp('Acceptation rate =') ;
    %           disp(acpt) ;
    if opts.alpha0<=acpt && acpt<=opts.alpha1
        %                disp('done!');
        stop=1;
    else
        if acpt<lower_prob
            c = c/5;
        elseif lower_prob<=acpt && acpt<=upper_prob
            c = c*log(.5*(opts.alpha0+opts.alpha1))/log(acpt);
        else
            c = 5*c;
        end
        %                disp('Trying with c= ') ;
        %                disp(c)
    end
    outer_iter=outer_iter+1;
end
if outer_iter==200
    error('dsmh: Outer loop reached maximum of iterations.')
end

function [out_param,out_tlogpost_iminus1,out_loglik] = mutation_DSMH(funobj, param, tlogpost_i, tlogpost_iminus1, loglik, lambda, i, c, MM, Omegachol, weights, mh_bounds, opts, Prior)
indx_levels = (1:1:MM-1)*opts.N*opts.G/MM;
npar = size(param,1) ;
p = 1/(10*opts.tau);
%        disp('Metropolis step...');
% build the dynamic grid of levels
levels = [0.0;tlogpost_iminus1(indx_levels)];
% initialize the outputs
out_param = param;
out_tlogpost_iminus1 = tlogpost_i;
out_loglik = loglik;
% resample and initialize the starting groups
indx_resmpl = kitagawa(weights,rand(1,1),opts.G);
param0 = param(:,indx_resmpl);
tlogpost_iminus10 = tlogpost_iminus1(indx_resmpl);
tlogpost_i0 = tlogpost_i(indx_resmpl);
loglik0 = loglik(indx_resmpl);
% Start the Metropolis
for l=1:opts.N*opts.tau
    for j=1:opts.G
        u1 = rand(1,1);
        u2 = rand(1,1);
        if u1<p
            k=1 ;
            for m=1:MM-1
                if levels(m)<=tlogpost_iminus10(j) && tlogpost_iminus10(j)<levels(m+1)
                    k = m+1;
                    break
                end
            end
            indx = floor( (k-1)*opts.N*opts.G/MM+1 + u2*(opts.N*opts.G/MM-1) );
            if i==1
                alp = (loglik(indx)-loglik0(j))*lambda(i);
            else
                alp = (loglik(indx)-loglik0(j))*(lambda(i)-lambda(i-1));
            end
            if u2<exp(alp)
                param0(:,j) = param(:,indx);
                tlogpost_i0(j) = tlogpost_i(indx);
                loglik0(j) = loglik(indx);
                tlogpost_iminus10(j) = tlogpost_iminus1(indx);
            end
        else
            validate= 0;
            while validate==0
                candidate = param0(:,j) + sqrt(c)*Omegachol*randn(npar,1);
                if all(candidate(:) >= mh_bounds.lb) && all(candidate(:) <= mh_bounds.ub)
                    [tlogpostx, loglikx] = tempered_likelihood(funobj, candidate, lambda(i), Prior);
                    if isfinite(loglikx) % if returned log-density is not Inf or Nan (penalized value)
                        validate = 1;
                        if u2<exp(tlogpostx-tlogpost_i0(j)) % accept
                            param0(:,j) = candidate;
                            tlogpost_i0(j) = tlogpostx;
                            loglik0(j) = loglikx;
                            if i==1
                                tlogpost_iminus10(j) = tlogpostx-loglikx*lambda(i);
                            else
                                tlogpost_iminus10(j) = tlogpostx-loglikx*(lambda(i)-lambda(i-1));
                            end
                        end
                    end
                end
            end
        end
    end
    if mod(l,opts.tau)==0
        rang = (l/opts.tau-1)*opts.G+1:l*opts.G/opts.tau;
        out_param(:,rang) = param0;
        out_tlogpost_iminus1(rang) = tlogpost_i0;
        out_loglik(rang) = loglik0;
    end
end
