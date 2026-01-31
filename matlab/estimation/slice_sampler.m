function [theta, fxsim, neval, sampler_options] = slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin)
% [theta, fxsim, neval, sampler_options] = slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin)
% ----------------------------------------------------------
% UNIVARIATE SLICE SAMPLER - stepping out (Neal, 2003)
% W: optimal value in the range (3,10)*std(x)
%    - see C.Planas and A.Rossi (2014)
% objective_function(theta,varargin): -log of any unnormalized pdf
% with varargin (optional) a vector of auxiliary parameters
% to be passed to f( ).
% ----------------------------------------------------------
%
% INPUTS
%   objective_function:       objective function (expressed as minus the log of a density)
%   theta:                    last value of theta
%   thetaprior:               bounds of the theta space
%   sampler_options:          posterior sampler options
%   varargin:                 optional input arguments to objective function
%
% OUTPUTS
%   theta:       new theta sample
%   fxsim:       value of the objective function for the new sample
%   neval:       number of function evaluations
%   sampler_options:          posterior sampler options
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2015-2025 Dynare Team
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

endo_init_state = false;
draw_endo_init_state_from_smoother=false;
draw_endo_init_state_with_rotated_slice = false;
if isfield(sampler_options,'draw_init_state_with_rotated_slice') && sampler_options.draw_init_state_with_rotated_slice
    endo_init_state = true;
    draw_endo_init_state_with_rotated_slice = true;
end
if isfield(sampler_options,'draw_init_state_from_smoother') && sampler_options.draw_init_state_from_smoother
    endo_init_state = true;
    draw_endo_init_state_from_smoother=true;
end

if endo_init_state
    [index_init_state, IS, index_deep_parameters] = get_init_state_estim_params(varargin{4}, varargin{6}, varargin{8});
    thetaprior(index_init_state,1) = theta(index_init_state);
    thetaprior(index_init_state,2) = theta(index_init_state);
end

if sampler_options.rotated %&& ~isempty(sampler_options.V1),
    sampler_options.endo_init_state.status = endo_init_state;
    if endo_init_state
        sampler_options.endo_init_state.IB = index_init_state;
        sampler_options.endo_init_state.IP = index_deep_parameters;
    end
    [theta, fxsim, neval] = rotated_slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin{:});
    if endo_init_state
        % draw initial states
        if draw_endo_init_state_from_smoother
            [theta, fxsim] = draw_init_state_from_smoother([false 5],sampler_options,theta,fxsim,thetaprior,varargin{:});
        elseif draw_endo_init_state_with_rotated_slice
            [V, D]=get_init_state_prior(theta,varargin{3:end});
            % take eigenvectors of state priors and set zero wieghts for other
            % params
            nslice = size(V,2);
            V=V(IS,:);
            V1 = zeros(length(theta),size(V,2));
            V1(index_init_state,:) = V;
            sampler_options.V1=V1;
            stderr = sqrt(diag(D));
            sampler_options.WR=stderr*3;
            for k=1:nslice
                bounds.lb(k) = norminv(1e-10, 0, stderr(k));
                bounds.ub(k) = norminv(1-1e-10, 0, stderr(k));
            end
            sampler_options.rthetaprior=[bounds.lb(:) bounds.ub(:)];
             % here params are fixed, so no need to account for change in
             % state prior!
            sampler_options.endo_init_state.status = false;
            %         sampler_options.WR=sampler_options.initial_step_size*(bounds.ub-bounds.lb);
            [theta, fxsim, neval1] = rotated_slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin{:});
            [~, icheck]=set_init_state(theta,varargin{3:end});
            neval(index_init_state(1:nslice)) = neval1(1:nslice);
        end
    end
    if isempty(sampler_options.mode) % jumping
        return
    else
        nevalR=sum(neval);
    end
end

thetaprior0=thetaprior;
if isfield(sampler_options,'fast_likelihood_evaluation_for_rejection') && sampler_options.fast_likelihood_evaluation_for_rejection
    fast_likelihood_evaluation_for_rejection = true;
    rejection_penalty=sampler_options.fast_likelihood_evaluation_for_rejection_penalty;
else
    fast_likelihood_evaluation_for_rejection = false;
end

use_prior_draws = false;
if isfield(sampler_options,'use_prior_draws') && sampler_options.use_prior_draws.status
    use_prior_draws = sampler_options.use_prior_draws.status;
    try_prior_draws = sampler_options.use_prior_draws.mh_blck(sampler_options.curr_block);
end

theta=theta(:);
npar = length(theta);
W1 = sampler_options.W1;
neval = zeros(npar,1);

fname = [ int2str(sampler_options.curr_block)];

Prior = dprior(varargin{6},varargin{3}.prior_trunc);

if use_prior_draws && try_prior_draws
    fxsim = sampler_options.last_posterior;
    Z1 = fxsim + log(rand(1,1));
    ilogpo2=-inf;
    nattempts=0;
    while ilogpo2<Z1 && nattempts<10
        nattempts=nattempts+1;
        validate=false;
        while not(validate)
            candidate = Prior.draw();
            if all(candidate >= thetaprior0(:,1)) && all(candidate <= thetaprior0(:,2))
                if fast_likelihood_evaluation_for_rejection
                    itest = -rejection_objective_function(objective_function,theta,Z1-rejection_penalty,varargin{:});
                else
                    itest = -feval(objective_function,candidate,varargin{:});
                end
                if isfinite(itest)
                    validate=true;
                end
            end
        end
        if itest>ilogpo2
            ilogpo2= itest;
            best_candidate = candidate;
        end
    end
    if ilogpo2>=Z1
        theta= best_candidate;
        fxsim = ilogpo2;
    else
        % prior draw is worse than current draw, so I stop drawing from
        % prior in this chain
        sampler_options.use_prior_draws.mh_blck(sampler_options.curr_block)=false;
    end
end


it=0;
islow=false(npar,1);
while it<npar
    it=it+1;
    neval(it) = 0;
    W = W1(it);
    xold  = theta(it);
    theta0=theta;
    XLB   = thetaprior(it,1);
    XUB   = thetaprior(it,2);
    if XLB==XUB
        continue
    end


    % -------------------------------------------------------
    % 1. DRAW Z = ln[f(X0)] - EXP(1) where EXP(1)=-ln(U(0,1))
    %    THIS DEFINES THE SLICE S={x: z < ln(f(x))}
    % -------------------------------------------------------
    fxold = -feval(objective_function,theta,varargin{:});
    if endo_init_state
        ys0 = get_steady_state(theta,varargin{3:end});
    end
    if ~isfinite(fxold)
        disp(['slice_sampler:: Iteration ' int2str(it) ' started with bad parameter set (fval is inf or nan)'])
        icount=0;
        while ~isfinite(fxold) && icount<1000
            icount=icount+1;
            theta = Prior.draw();
            if all(theta >= thetaprior(:,1)) && all(theta <= thetaprior(:,2))
                fxold = -feval(objective_function,theta,varargin{:});
            end
        end
        % restart from 1
        it = 1;
        neval(it) = 0;
        W = W1(it);
        xold  = theta(it);
        XLB   = thetaprior(it,1);
        XUB   = thetaprior(it,2);
    end
    neval(it) = neval(it) + 1;
    Z = fxold + log(rand(1,1));
    % -------------------------------------------------------------
    % 2. FIND I=(L,R) AROUND X0 THAT CONTAINS S AS MUCH AS POSSIBLE
    %    STEPPING-OUT PROCEDURE
    % -------------------------------------------------------------
    u = rand(1,1);
    L = max(XLB,xold-W*u);
    R = min(XUB,L+W);
    mytxt{it,1} = '';
    while(L > XLB)
        xsim = L;
        theta(it) = xsim;
        if endo_init_state
            theta0(it) = xsim;
            [theta, icheck]=set_init_state(theta0,ys0,varargin{3:end});
        end
        if fast_likelihood_evaluation_for_rejection
            fxl = -rejection_objective_function(objective_function,theta,Z-rejection_penalty,varargin{:});
         else
            fxl = -feval(objective_function,theta,varargin{:});
        end
        neval(it) = neval(it) + 1;
        if (fxl <= Z)
            break
        end
        L = max(XLB,L-W);
        if neval(it)>30
            L=XLB;
            xsim = L;
            theta(it) = xsim;
            if endo_init_state
                theta0(it) = xsim;
                [theta, icheck]=set_init_state(theta0,varargin{3:end});
            end
            fxl = -feval(objective_function,theta,varargin{:});
            icount = 0;
            while (isinf(fxl) || isnan(fxl)) && icount<300
                icount = icount+1;
                L=L+sqrt(eps);
                xsim = L;
                theta(it) = xsim;
                if endo_init_state
                    theta0(it) = xsim;
                    [theta, icheck]=set_init_state(theta0,varargin{3:end});
                end
                fxl = -feval(objective_function,theta,varargin{:});
            end
            mytxt{it,1} = sprintf('Getting L for [%s] is taking too long.', varargin{6}.name{it});
            if sampler_options.save_iter_info_file
                save([varargin{4}.dname filesep 'metropolis/slice_iter_info_' fname],'mytxt','neval','it','theta','fxl')
            end
        end
    end
    neval1 = neval(it);
    mytxt{it,2} = '';
    while(R < XUB)
        xsim = R;
        theta(it) = xsim;
        if endo_init_state
            theta0(it) = xsim;
            [theta, icheck]=set_init_state(theta0,ys0,varargin{3:end});
        end
        if fast_likelihood_evaluation_for_rejection
            fxr = -rejection_objective_function(objective_function,theta,Z-rejection_penalty,varargin{:});
        else
            fxr = -feval(objective_function,theta,varargin{:});
        end
        neval(it) = neval(it) + 1;
        if (fxr <= Z)
            break
        end
        R = min(XUB,R+W);
        if neval(it)>(neval1+30)
            R=XUB;
            xsim = R;
            theta(it) = xsim;
            if endo_init_state
                theta0(it) = xsim;
                [theta, icheck]=set_init_state(theta0,ys0,varargin{3:end});
            end
            fxr = -feval(objective_function,theta,varargin{:});
            icount = 0;
            while (isinf(fxr) || isnan(fxr)) && icount<300
                icount = icount+1;
                R=R-sqrt(eps);
                xsim = R;
                theta(it) = xsim;
                if endo_init_state
                    theta0(it) = xsim;
                    [theta, icheck]=set_init_state(theta0,ys0,varargin{3:end});
                end
                fxr = -feval(objective_function,theta,varargin{:});
            end
            mytxt{it,2} = sprintf('Getting R for [%s] is taking too long.', varargin{6}.name{it});
            if sampler_options.save_iter_info_file
                save([varargin{4}.dname filesep 'metropolis/slice_iter_info_' fname],'mytxt','neval','it','theta','fxr')
            end
        end
    end
    % ------------------------------------------------------
    % 3. SAMPLING FROM THE SET A = (I INTERSECT S) = (LA,RA)
    % ------------------------------------------------------
    fxsim = Z-1;
    mytxt{it,3} = '';
    while (fxsim < Z)
        u = rand(1,1);
        xsim = L + u*(R - L);
        theta(it) = xsim;
        if endo_init_state
            theta0(it) = xsim;
            [theta, icheck]=set_init_state(theta0,ys0,varargin{3:end});
        end
        if fast_likelihood_evaluation_for_rejection
            fxsim = -rejection_objective_function(objective_function,theta,Z-rejection_penalty,varargin{:});
        else
            fxsim = -feval(objective_function,theta,varargin{:});
        end
        neval(it) = neval(it) + 1;
        if (xsim > xold)
            R = xsim;
        else
            L = xsim;
        end
        if (R-L)<1.e-6 %neval(it)>(30+neval2)
            fprintf('The sampling for parameter [%s] is taking too long as the sampling set is too tight. Check the prior.\n', varargin{6}.name{it})
            mytxt{it,3} = sprintf('Sampling [%s] is taking too long.', varargin{6}.name{it});
            if sampler_options.save_iter_info_file
                save([varargin{4}.dname filesep 'metropolis/slice_iter_info_' fname],'mytxt','neval','it')
            end
            islow(it)=true;
            theta(it) = xold;
            fxsim = fxold;
            break
        end
    end
    if sampler_options.save_iter_info_file
        save([varargin{4}.dname filesep 'metropolis/slice_iter_info_' fname],'mytxt','neval','it','theta','fxsim')
    end
    if isinf(fxsim) || isnan(fxsim)
        theta(it) = xold;
        fxsim = fxold;
        disp('SLICE: posterior density is infinite. Reset values at initial ones.')
    end
    if endo_init_state && icheck %(icheck || islow)
        [theta, fxsim] = draw_init_state_from_smoother([false 1],sampler_options,theta,fxsim,thetaprior,varargin{:});
    end
    if isinf(fxsim) || isnan(fxsim)
        theta(it) = xold;
        fxsim = fxold;
        disp('SLICE: posterior density is infinite. Reset values at initial ones.')
    end
    if sampler_options.save_iter_info_file
        save([varargin{4}.dname filesep 'metropolis/slice_iter_info_' fname],'mytxt','neval','it','theta','fxsim')
    end
end

if any(islow)
    fxsim = -feval(objective_function,theta,varargin{:});
    Z1 = fxsim + log(rand(1,1));
    ilogpo2=-inf;
    nattempts=0;
    while ilogpo2<Z1 && nattempts<10
        nattempts=nattempts+1;


        validate=false;
        while not(validate)
            candidate = theta;
            my_candidate = Prior.draw();
            candidate(islow) = my_candidate(islow);

            if all(candidate(islow) >= thetaprior(islow,1)) && all(candidate(islow) <= thetaprior(islow,2))
                if endo_init_state
                    init = true;
                    [candidate] = draw_init_state_from_smoother(init,sampler_options,candidate,nan,thetaprior,varargin{:});
                end
                itest = -feval(objective_function,candidate,varargin{:});
                if isfinite(itest)
                    validate=true;
                end

            end
        end
        if itest>ilogpo2
            ilogpo2= itest;
            best_candidate = candidate;
        end
    end
    theta= best_candidate;
    fxsim = ilogpo2;
end
if endo_init_state
    % draw initial states
    if draw_endo_init_state_from_smoother
        [theta, fxsim] = draw_init_state_from_smoother([false 5],sampler_options,theta,fxsim,thetaprior,varargin{:});
    elseif draw_endo_init_state_with_rotated_slice
        [V, D]=get_init_state_prior(theta,varargin{3:end});
        % take eigenvectors of state priors and set zero wieghts for other
        % params
        nslice = size(V,2);
        V=V(IS,:);
        V1 = zeros(length(theta),size(V,2));
        V1(index_init_state,:) = V;
        sampler_options.V1=V1;
        stderr = sqrt(diag(D));
        sampler_options.WR=stderr*3;
        for k=1:nslice
            bounds.lb(k) = norminv(1e-10, 0, stderr(k));
            bounds.ub(k) = norminv(1-1e-10, 0, stderr(k));
        end
        sampler_options.rthetaprior=[bounds.lb(:) bounds.ub(:)];
        %         sampler_options.WR=sampler_options.initial_step_size*(bounds.ub-bounds.lb);
        [theta, fxsim, neval1] = rotated_slice_sampler(objective_function,theta,thetaprior,sampler_options,varargin{:});
        [~, icheck]=set_init_state(theta,ys0,varargin{3:end});
        neval(index_init_state(1:nslice)) = neval1(1:nslice);
    end
    save([varargin{4}.dname filesep 'metropolis/slice_iter_info_' fname],'mytxt','neval','it','theta','fxsim')
end

if sampler_options.rotated && ~isempty(sampler_options.mode) % jumping
    neval=sum(neval)+nevalR;
end


function ys = get_steady_state(xparam1, options_,M_,estim_params_,~,~,~, endo_steady_state, exo_steady_state, exo_det_steady_state)
% wrapper function to get steady state 

M_ = set_all_parameters(xparam1,estim_params_,M_);
ys = evaluate_steady_state(endo_steady_state,[exo_steady_state; exo_det_steady_state],M_,options_,true);
