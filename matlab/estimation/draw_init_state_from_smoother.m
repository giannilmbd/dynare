function [xparam1, logpost0, mh_bounds, M_] = draw_init_state_from_smoother(init,sampler_options,xparam1,logpost0,mh_bounds, ...
    dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info)
% [xparam1, logpost0, mh_bounds, M_] = draw_init_state_from_smoother(init,sampler_options,xparam1,logpost0,mh_bounds, ...
%     dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info)
% Computes the endogenous log prior addition to the initial prior
%
% INPUTS
% - xparam1             [double]        current values for the estimated parameters.
% - mh_bounds           [structure]     containing mh_bounds
% - dataset_            [structure]     dataset after transformations
% - dataset_info        [structure]     storing informations about the
%                                       sample; not used but required for interface
% - options_            [structure]     Matlab's structure describing the current options
% - M_                  [structure]     Matlab's structure describing the model
% - estim_params_       [structure]     characterizing parameters to be estimated
% - bayestopt_          [structure]     describing the priors
% - BoundsInfo          [structure]     containing prior bounds
% - dr                  [structure]     Reduced form model.
% - endo_steady_state   [vector]        steady state value for endogenous variables
% - exo_steady_state    [vector]        steady state value for exogenous variables
% - exo_det_steady_state [vector]       steady state value for exogenous deterministic variables
% - derivatives_info    [structure]     derivative info for identification
%
% OUTPUTS
% - xparam1             [double]        current values for the estimated parameters.
% - mh_bounds           [structure]     containing mh_bounds
% - M_                  [structure]     Matlab's structure describing the model

% Copyright © 2024-2025 Dynare Team
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

filter_initial_state=M_.filter_initial_state;
M_.filter_initial_state = [];
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
data = dataset_.data;
data_index = dataset_info.missing.aindex;
missing_value = dataset_info.missing.state;

% Set number of observations
gend = dataset_.nobs;
options_.smoother_redux=true;
options_.smoothed_state_uncertainty = true;
options_.occbin.smoother.debug = false;
options_.occbin.smoother.plot = false;
options_.occbin.smoother.store_results = false;
options_.occbin.smoother.waitbar = false;
M_ = set_all_parameters(xparam1,estim_params_,M_);
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
        options_.init_state_endogenous_prior=false;
        logpost2  = -rejection_objective_function(@dsge_likelihood,xparam1,logpost0-10,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);        options_.init_state_endogenous_prior=true;
        if (logpost0-logpost2)<1.e3
            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,oo_,bayestopt_.mf,alphahat0,state_uncertainty0] = occbin.DSGE_smoother(xparam1,gend,transpose(data),data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_,dataset_,dataset_info);
        else
            oo_.occbin.smoother.error_flag = 313;
            alphahat0 = [];
        end
    end
    if init || (oo_.occbin.smoother.error_flag && isempty(alphahat0))
        options_.occbin.smoother.status=false;
        [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,bayestopt_.mf,alphahat0,state_uncertainty0] = DsgeSmoother(xparam1,gend,transpose(data),data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
        options_.occbin.smoother.status=true;
    end
else
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,bayestopt_,alphahat0,state_uncertainty0] = DsgeSmoother(xparam1,gend,transpose(data),data_index,missing_value,M_,oo_,options_,bayestopt_,estim_params_);
end
% end
M_.filter_initial_state = filter_initial_state;
options_.lik_init=2;
if error_flag==0
    % smoother converged!
    params0 = M_.params;
    % draw initial state from smoothed distribution
    [U,X] = svd(0.5*(state_uncertainty0(dr.restrict_var_list(bayestopt_.mf0),dr.restrict_var_list(bayestopt_.mf0))+state_uncertainty0(dr.restrict_var_list(bayestopt_.mf0),dr.restrict_var_list(bayestopt_.mf0))'));
    is = find(diag(X)>options_.kalman_tol);
    StateVectorVarianceSquareRoot = chol(X(is,is))';

    % Get the rank of StateVectorVarianceSquareRoot
    state_variance_rank = size(StateVectorVarianceSquareRoot,2);
    U = U(:,is);

    pvec=[];
    for ii=1:size(M_.filter_initial_state,1)
        if ~isempty(M_.filter_initial_state{ii,1})
            tmp1 = strrep(M_.filter_initial_state{ii,2},');','');
            tmp1 = strrep(tmp1,'M_.params(','');
            pvec = [pvec eval(tmp1)];
        end
    end
    [~,~,IB] = intersect(M_.param_names(pvec),bayestopt_.name,'stable');
    M_.params(pvec) = xparam1(IB);

    fast_likelihood_evaluation_for_rejection = false;
    if isfield(sampler_options,'fast_likelihood_evaluation_for_rejection') && sampler_options.fast_likelihood_evaluation_for_rejection
        fast_likelihood_evaluation_for_rejection = true;
    end

    alphahat00 = alphahat0;
    if not(init)
        if options_.init_state_endogenous_prior
            % independent MH
            % check probability of smoothed in t=0 (which is the mode for
            % linear case, but how about occbin?)
            yhat = alphahat0(dr.restrict_var_list(bayestopt_.mf0));
            alphahat01=alphahat00;
            alphahat01(dr.restrict_var_list(bayestopt_.mf0))=yhat;
            alphahat01 = alphahat01(dr.inv_order_var);

            M_=update_parameters_filter_initial_state(M_,alphahat01,dr.ys,options_);
            xparam11=xparam1;
            xparam11(IB) = M_.params(pvec);
            if not(all(xparam11(:)>=sampler_options.bounds.lb) && all(xparam11(:)<=sampler_options.bounds.ub))
                xparam11(xparam11(:)<sampler_options.bounds.lb)=sampler_options.bounds.lb(xparam11(:)<sampler_options.bounds.lb)+sqrt(eps);
                xparam11(xparam11(:)>sampler_options.bounds.ub)=sampler_options.bounds.ub(xparam11(:)>sampler_options.bounds.ub)-sqrt(eps);
            end

            is_smoothed_state_optimal=true;
            [xcheck, icheck]=set_init_state(xparam11,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
            if icheck
                xparam11=xcheck;
            end
            logpost1 = -dsge_likelihood(xparam11,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
            if logpost1<logpost0
                is_smoothed_state_optimal=false;
            end
            logpostSMO = logpost1;
        end
        if not(options_.init_state_endogenous_prior) || not(is_smoothed_state_optimal)
            % use previous draw RW MH
            for ii=1:size(M_.filter_initial_state,1)
                if ~isempty(M_.filter_initial_state{ii,1})
                    if options_.loglinear && ~options_.logged_steady_state
                        eval(['alphahat0(ii) = log(' strrep(M_.filter_initial_state{ii,2},';','') ') - log(dr.ys(ii));']);
                    elseif ~options_.loglinear && ~options_.logged_steady_state
                        eval(['alphahat0(ii) = ' strrep(M_.filter_initial_state{ii,2},';','') '- dr.ys(ii);'])
                    else
                        error('The steady state is logged. This should not happen. Please contact the developers')
                    end
                end
            end
            alphahat0=alphahat0(dr.order_var);
        end
    end
    naccepted=0;
    nattempts=0;
    while naccepted<target_accepted && nattempts<10
        niter = 0;
        nattempts = nattempts+1;
        while naccepted<target_accepted && niter<20
            niter = niter+1;
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
                alphahat01=alphahat00;
                alphahat01(dr.restrict_var_list(bayestopt_.mf0))=yhat;
                alphahat01 = alphahat01(dr.inv_order_var);

                M_=update_parameters_filter_initial_state(M_,alphahat01,dr.ys,options_);
                xparam11=xparam1;
                xparam11(IB) = M_.params(pvec);
                if all(xparam11(:)>=sampler_options.bounds.lb) && all(xparam11(:)<=sampler_options.bounds.ub)
                    new_draw_out_of_bounds = false;
                end
            end
            if new_draw_out_of_bounds
                xparam11(xparam11(:)<sampler_options.bounds.lb)=sampler_options.bounds.lb(xparam11(:)<sampler_options.bounds.lb)+sqrt(eps);
                xparam11(xparam11(:)>sampler_options.bounds.ub)=sampler_options.bounds.ub(xparam11(:)>sampler_options.bounds.ub)-sqrt(eps);
            end
            [xcheck, icheck]=set_init_state(xparam11,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
            if icheck
                xparam11=xcheck;
            end
            if init
                xparam1=xparam11;
                break
            end
            lnrand = log(rand);
            if fast_likelihood_evaluation_for_rejection
                fval=lnrand+logpost0-10;
                logpost1  = -rejection_objective_function(@dsge_likelihood,xparam11,fval,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
                if (logpost1 >= fval)
                    logcheck = -dsge_likelihood(xparam11,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
                    logpost1 = logcheck;
                end
            else
                logpost1 = -dsge_likelihood(xparam11,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info);
            end
            if logpost1<logpost0
                r = logpost1-logpost0;
                if (logpost1 > -inf) && (lnrand < r)
                    accepted = 1;
                    xparam1 = xparam11;
                else
                    accepted = 0;
                    M_.params = params0;
                    logpost1 = logpost0;
                end
            else
                xparam1 = xparam11;
                accepted = 1;
            end
            if accepted
                logpost0 = logpost1;
                naccepted = naccepted+1;
                M_.params(pvec) = xparam1(IB);
                if options_.init_state_endogenous_prior && logpostSMO<logpost0
                    is_smoothed_state_optimal=false;
                end

                if not(options_.init_state_endogenous_prior && is_smoothed_state_optimal) %nattempts==1)
                    for ii=1:size(M_.filter_initial_state,1)
                        if ~isempty(M_.filter_initial_state{ii,1})
                            if options_.loglinear && ~options_.logged_steady_state
                                eval(['alphahat0(ii) = log(' strrep(M_.filter_initial_state{ii,2},';','') ') - log(dr.ys(ii));']);
                            elseif ~options_.loglinear && ~options_.logged_steady_state
                                eval(['alphahat0(ii) = ' strrep(M_.filter_initial_state{ii,2},';','') '- dr.ys(ii);'])
                            else
                                error('The steady state is logged. This should not happen. Please contact the developers')
                            end
                        end
                    end
                    alphahat0=alphahat0(dr.order_var);
                end
            end
        end
        if init
            break
        end
        if naccepted==0 && options_.init_state_endogenous_prior
            if not(is_smoothed_state_optimal)
                M_.params(pvec) = xparam1(IB);
                alphahat00 = alphahat0;

                for ii=1:size(M_.filter_initial_state,1)
                    if ~isempty(M_.filter_initial_state{ii,1})
                        if options_.loglinear && ~options_.logged_steady_state
                            eval(['alphahat0(ii) = log(' strrep(M_.filter_initial_state{ii,2},';','') ') - log(dr.ys(ii));']);
                        elseif ~options_.loglinear && ~options_.logged_steady_state
                            eval(['alphahat0(ii) = ' strrep(M_.filter_initial_state{ii,2},';','') '- dr.ys(ii);'])
                        else
                            error('The steady state is logged. This should not happen. Please contact the developers')
                        end
                    end
                end
                alphahat0=alphahat0(dr.order_var);
            else
                StateVectorVarianceSquareRoot = StateVectorVarianceSquareRoot*0.66;
            end
        end
    end
end
if not(init)
    if isstruct(mh_bounds)
        mh_bounds.lb(IB)= xparam1(IB);
        mh_bounds.ub(IB)= xparam1(IB);
    end
    disp_verbose(['naccepted=' int2str(naccepted) '| niter=' int2str(niter) '| nattempts=' int2str(nattempts)],options_.debug)
end

