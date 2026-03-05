function [xparam1, logpost0, mh_bounds, M_, neval] = draw_init_state_from_smoother(init,sampler_options,xparam1,logpost0,mh_bounds, ...
    dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info)
% [xparam1, logpost0, mh_bounds, M_, neval] = draw_init_state_from_smoother(init,sampler_options,xparam1,logpost0,mh_bounds, ...
%     dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info)
% Draws the Kalman filter initial state from the DSGE smoother.
%
% The function runs an unconditional state smoother (OccBin or linear)
% to obtain the smoothed state mean and covariance, then proposes an
% initial state draw. If init is false, it performs a short Metropolis–
% Hastings step to update the endogenous initial state parameters in
% xparam1 (independent MH centered at the smoothed state when
% options_.estimate_initial_states_endogenous_prior is true; otherwise a RW-MH around
% the current state), honoring parameter bounds. If init is true, it
% directly sets the initial state parameters in xparam1 from the smoothed
% draw without acceptance testing.
%
% INPUTS
% - init                [logical|double] when logical, selects mode:
%                        true  → initialize from smoothed draw without MH;
%                        false → sample via MH around smoothed/current state.
%                        When a 2-element vector [flag, target], uses
%                        flag as above and sets target_accepted=target.
% - sampler_options     [structure]     MCMC/sampler options; uses fields
%                                       .bounds.lb/.ub and optionally
%                                       .fast_likelihood_evaluation_for_rejection.
% - xparam1             [double]        current values for the estimated parameters
%                                       (includes endogenous initial state parameters).
% - logpost0            [double]        current log-posterior at xparam1.
% - mh_bounds           [structure]     optional; when struct, may be tightened
%                                       to pin initial state parameters after sampling.
% - dataset_            [structure]     dataset after transformations
% - dataset_info        [structure]     information about the sample (missing data, indices)
% - options_            [structure]     options; toggles lik_init and smoother controls
% - M_                  [structure]     model structure (updated endo_initial_state)
% - estim_params_       [structure]     parameters to be estimated
% - bayestopt_          [structure]     priors and state indexing (mf0)
% - BoundsInfo          [structure]     prior bounds info
% - dr                  [structure]     reduced-form decision rules
% - endo_steady_state   [vector]        steady state for endogenous variables
% - exo_steady_state    [vector]        steady state for exogenous variables
% - exo_det_steady_state [vector]       steady state for deterministic exogenous variables
% - derivatives_info    [structure]     derivative info for identification
%
% OUTPUTS
% - xparam1             [double]        updated parameter vector with initial state entries set.
% - logpost0            [double]        updated log-posterior at returned xparam1.
% - mh_bounds           [structure]     possibly tightened bounds for initial state parameters.
% - M_                  [structure]     updated model with endo_initial_state values.
% - neval               [double]        number of objective function evaluations.
%
% SEE ALSO
%   get_init_state_prior - computes Pstar-based constraints for initial-state consistency.
%   set_init_state - enforces null-space projection on init-state parameters.

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

% varargin
% 1        2            3        4  5             6          7          8   9                  10                11                   12
% dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info

options_.noprint = ~options_.debug;

neval = 0;

if options_.occbin.likelihood.status
    options_.occbin.smoother.status = true;
end

if nargin<17
    derivatives_info=[];
end

target_accepted = 5;
if not(islogical(init))
    target_accepted = init(2);
    init=logical(init(1));
end
if init
    target_accepted = 1;
end

% here I run unconditional smoother, so I need to undo the init state
% estimation setup and set lik_init = 1
M_ = set_all_parameters(xparam1,estim_params_,M_);
ys0 = evaluate_steady_state(endo_steady_state,[exo_steady_state; exo_det_steady_state],M_,options_,true);
store_endo_initial_state=M_.endo_initial_state;
M_.endo_initial_state.status = false;
error_flag=0;
options_.lik_init=1;
if init
    [Pstar, info]=get_pstar(xparam1,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
    if info(1)
        return
    end
    state_uncertainty0=zeros(M_.endo_nbr,M_.endo_nbr);
    state_uncertainty0(dr.restrict_var_list,dr.restrict_var_list)=Pstar;
    alphahat0=zeros(M_.endo_nbr,1);
end

options_.smoother_redux=true;
options_.smoothed_state_uncertainty = true;
options_.occbin.smoother.debug = false;
options_.occbin.smoother.plot = false;
options_.occbin.smoother.store_results = false;
options_.occbin.smoother.waitbar = false;
dr.ys = evaluate_steady_state(endo_steady_state,[exo_steady_state; exo_det_steady_state],M_,options_,true);
oo_.dr = dr;
oo_.steady_state= endo_steady_state;
oo_.exo_steady_state = exo_steady_state;
oo_.exo_det_steady_state = exo_det_steady_state;
options_.verbosity=false;
if options_.occbin.smoother.status
    if not(init)
        % check first that PKF with latent states provides sensible
        % likelihood
        options_.estimate_initial_states_endogenous_prior=false;
        logpost2  = -rejection_objective_function(@dsge_likelihood,xparam1,logpost0-10,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
        neval = neval + 1;
        options_.estimate_initial_states_endogenous_prior=true;
        if (logpost0-logpost2)<1.e3
            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,oo_,bayestopt_.mf,alphahat0,state_uncertainty0] = occbin.DSGE_smoother(xparam1,M_,oo_,options_,bayestopt_,estim_params_,dataset_,dataset_info,true);
        else
            oo_.occbin.smoother.error_flag = 313;
            alphahat0 = [];
        end
    end
    if init || (oo_.occbin.smoother.error_flag && isempty(alphahat0))
        % use linear smoother to initialize or if error
        options_.occbin.smoother.status=false;
        [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,bayestopt_.mf,alphahat0,state_uncertainty0] = DsgeSmoother(xparam1,dataset_,dataset_info,M_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,options_,bayestopt_,estim_params_);
        options_.occbin.smoother.status=true;
    end
else
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,bayestopt_.mf,alphahat0,state_uncertainty0] = DsgeSmoother(xparam1,dataset_,dataset_info,M_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,options_,bayestopt_,estim_params_);
end
% end unconditional smoother to get mean (alphahat0) and covariance (state_uncertainty0) of the proposal for init state
% now I reset init state estimation stuff
M_.endo_initial_state = store_endo_initial_state;
options_.lik_init=2;
if error_flag==0
    % smoother converged!
    % draw initial state from smoothed distribution
    [U,X] = svd(0.5*(state_uncertainty0(dr.restrict_var_list(bayestopt_.mf0),dr.restrict_var_list(bayestopt_.mf0))+state_uncertainty0(dr.restrict_var_list(bayestopt_.mf0),dr.restrict_var_list(bayestopt_.mf0))'));
    is = find(diag(X)>options_.kalman_tol);
    StateVectorVarianceSquareRoot = chol(X(is,is))';

    % Get the rank of StateVectorVarianceSquareRoot
    state_variance_rank = size(StateVectorVarianceSquareRoot,2);
    U = U(:,is);

    % store current init state values
    IB = startsWith(bayestopt_.name, 'init ');
    M_.endo_initial_state.values(dr.state_var) = xparam1(IB);

    fast_likelihood_evaluation_for_rejection = false;
    if isfield(sampler_options,'fast_likelihood_evaluation_for_rejection') && sampler_options.fast_likelihood_evaluation_for_rejection
        fast_likelihood_evaluation_for_rejection = true;
    end

    if not(init)
        if options_.estimate_initial_states_endogenous_prior
            % independent MH
            % check probability of smoothed in t=0 (which is the mode for
            % linear case, but how about occbin?)
            yhat = alphahat0(dr.restrict_var_list(bayestopt_.mf0));
            % build declaration order init states alphahat01
            alphahat01=alphahat0;
            alphahat01(dr.restrict_var_list(bayestopt_.mf0))=yhat;
            alphahat01 = alphahat01(dr.inv_order_var);

            xproposal=xparam1;
            if options_.loglinear
                xproposal(IB) = exp(alphahat01(dr.state_var)).*dr.ys(dr.state_var);
            else
                xproposal(IB) = alphahat01(dr.state_var)+dr.ys(dr.state_var);
            end
            if not(all(xproposal(:)>=sampler_options.bounds.lb) && all(xproposal(:)<=sampler_options.bounds.ub))
                xproposal(xproposal(:)<sampler_options.bounds.lb)=sampler_options.bounds.lb(xproposal(:)<sampler_options.bounds.lb)+sqrt(eps);
                xproposal(xproposal(:)>sampler_options.bounds.ub)=sampler_options.bounds.ub(xproposal(:)>sampler_options.bounds.ub)-sqrt(eps);
            end

            is_smoothed_state_optimal=true;
            [xcheck, icheck]=set_init_state(xproposal,ys0,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
            if icheck
                xproposal=xcheck;
            end
            logpost1 = -dsge_likelihood(xproposal,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
            neval = neval + 1;
            if logpost1<logpost0
                is_smoothed_state_optimal=false;
            end
            logpostSMO = logpost1;
        end
        if not(options_.estimate_initial_states_endogenous_prior) || not(is_smoothed_state_optimal)
            % use previous draw RW MH
            if options_.loglinear
                alphahat0=log(store_endo_initial_state.values)-log(dr.ys);
            else
                alphahat0=store_endo_initial_state.values-dr.ys;
            end
            alphahat0=alphahat0(dr.order_var); % decision rule order
        end
    end
    naccepted=0;
    nattempts=0;
    disp_verbose('draw_init_state_from_smoother: Starting MH sampling for initial states', options_.debug);
    while naccepted<target_accepted && nattempts<10
        niter = 0;
        nattempts = nattempts+1;
        disp_verbose(sprintf('draw_init_state_from_smoother: Outer loop (variance scaling attempt) - attempt %d/10 (%.1f%%), accepted %d/%d (%.1f%%)', nattempts, 100*nattempts/10, naccepted, target_accepted, 100*naccepted/target_accepted), options_.debug);
        while naccepted<target_accepted && niter<20
            niter = niter+1;
            if mod(niter, 10) == 0
                disp_verbose(sprintf('draw_init_state_from_smoother:   MH iteration - iteration %d/20 (%.1f%%)', niter, 100*niter/20), options_.debug);
            end
            new_draw_out_of_bounds= true;
            icount = 0;
            while new_draw_out_of_bounds && icount<10
                icount = icount+1;
                draw_base = randn(state_variance_rank,1);
                if not(isempty(is))
                    yhat = U(:,is)*StateVectorVarianceSquareRoot*draw_base+alphahat0(dr.restrict_var_list(bayestopt_.mf0));
                else
                    yhat = alphahat0(dr.restrict_var_list(bayestopt_.mf0));
                end
                alphahat01=alphahat0;
                alphahat01(dr.restrict_var_list(bayestopt_.mf0))=yhat;
                alphahat01 = alphahat01(dr.inv_order_var);

                M_=update_parameters_filter_initial_state(M_,alphahat01,dr,options_);
                xproposal=xparam1;
                xproposal(IB) = M_.endo_initial_state.values(dr.state_var);
                if all(xproposal(:)>=sampler_options.bounds.lb) && all(xproposal(:)<=sampler_options.bounds.ub)
                    new_draw_out_of_bounds = false;
                end
            end
            if new_draw_out_of_bounds
                xproposal(xproposal(:)<sampler_options.bounds.lb)=sampler_options.bounds.lb(xproposal(:)<sampler_options.bounds.lb)+sqrt(eps);
                xproposal(xproposal(:)>sampler_options.bounds.ub)=sampler_options.bounds.ub(xproposal(:)>sampler_options.bounds.ub)-sqrt(eps);
            end
            [xcheck, icheck]=set_init_state(xproposal,ys0,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
            if icheck
                % if init states have been modified to match the null space
                % of Pstar
                xproposal=xcheck;
            end
            if init
                naccepted = 1;
                xparam1=xproposal;
                disp_verbose(sprintf('draw_init_state_from_smoother initialization:     ACCEPTED - accepted %d/%d (%.1f%%)', naccepted, target_accepted, 100*naccepted/target_accepted), options_.debug);
                break
            end
            lnrand = log(rand);
            if fast_likelihood_evaluation_for_rejection
                fval=lnrand+logpost0-10;
                logpost1  = -rejection_objective_function(@dsge_likelihood,xproposal,fval,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
                neval = neval + 1;
                if (logpost1 >= fval)
                    logcheck = -dsge_likelihood(xproposal,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
                    neval = neval + 1;
                    logpost1 = logcheck;
                end
            else
                logpost1 = -dsge_likelihood(xproposal,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
                neval = neval + 1;
            end
            if logpost1<logpost0
                r = logpost1-logpost0;
                if (logpost1 > -inf) && (lnrand < r)
                    accepted = 1;
                    xparam1 = xproposal;
                else
                    accepted = 0;
                    M_.endo_initial_state = store_endo_initial_state;
                    logpost1 = logpost0;
                end
            else
                xparam1 = xproposal;
                accepted = 1;
            end
            if accepted
                logpost0 = logpost1;
                naccepted = naccepted+1;
                disp_verbose(sprintf('draw_init_state_from_smoother:     ACCEPTED - accepted %d/%d (%.1f%%)', naccepted, target_accepted, 100*naccepted/target_accepted), options_.debug);
                M_.endo_initial_state.values(dr.state_var) = xparam1(IB);
                store_endo_initial_state = M_.endo_initial_state;
                if options_.estimate_initial_states_endogenous_prior && logpostSMO<logpost0
                    % switch from independent to RW Metropolis
                    is_smoothed_state_optimal=false;
                end

                if not(options_.estimate_initial_states_endogenous_prior && is_smoothed_state_optimal) %nattempts==1)
                    alphahat0=store_endo_initial_state.values-dr.ys;
                    alphahat0=alphahat0(dr.order_var); % decision rule order
                end
            end
        end
        if init
            break
        end
        if naccepted==0 && options_.estimate_initial_states_endogenous_prior
            % try reducing variance of state uncertainty in the
            % proposal and continue with MH centered on alphahat0
            StateVectorVarianceSquareRoot = StateVectorVarianceSquareRoot*0.66;
        end
    end
    if naccepted < target_accepted
        disp_verbose(sprintf('draw_init_state_from_smoother: WARNING - MH sampling terminated without reaching acceptance criterion - accepted %d/%d (%.1f%%), attempts %d/10\n', naccepted, target_accepted, 100*naccepted/target_accepted, nattempts), options_.debug);
    end
end
if not(init)
    if isstruct(mh_bounds)
        % set ub=lb
        mh_bounds.lb(IB)= xparam1(IB);
        mh_bounds.ub(IB)= xparam1(IB);
    end
    disp_verbose(['draw_init_state_from_smoother: Final - accepted=' int2str(naccepted) ' | iterations=' int2str(niter) ' | attempts=' int2str(nattempts)],options_.debug)
end

%% Local helper function
function M_=update_parameters_filter_initial_state(M_,alphahat01,dr,options_)
% Updates M_.endo_initial_state.values from state deviations.
%
% Given a state vector `alphahat01` in declaration order (representing
% deviations from steady state), this routine updates
% `M_.endo_initial_state.values` by adding the deviations to the steady
% state `ys`. In log-linear mode it uses exponential transformation;
% otherwise a direct level addition.

if options_.loglinear && ~options_.logged_steady_state
    M_.endo_initial_state.values(dr.state_var) = exp(log(dr.ys(dr.state_var))+alphahat01(dr.state_var));
elseif ~options_.loglinear && ~options_.logged_steady_state
    M_.endo_initial_state.values(dr.state_var)= dr.ys(dr.state_var)+alphahat01(dr.state_var);
else
    error('The steady state is logged. This should not happen. Please contact the developers')
end
