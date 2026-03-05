function [alphahat,etahat,epsilonhat,ahat0,SteadyState,trend_coeff,aKK,T0,R0,P,PKK,decomp,Trend,state_uncertainty,oo_,mf,alphahat0,state_uncertainty0] = DSGE_smoother(xparam1,M_,oo_,options_,bayestopt_,estim_params_,dataset_, dataset_info, store_linear_smoother)
% [alphahat,etahat,epsilonhat,ahat0,SteadyState,trend_coeff,aKK,T0,R0,P,PKK,decomp,Trend,state_uncertainty,oo_,bayestopt_,alphahat0,state_uncertainty0] = DSGE_smoother(xparam1,M_,oo_,options_,bayestopt_,estim_params_,dataset_, dataset_info, store_linear_smoother)
% Runs a DSGE smoother with occasionally binding constraints
%
% INPUTS
% - xparam1       [double]        (p*1) vector of (estimated) parameters.
% - M_            [structure]     MATLAB's structure describing the model (M_).
% - oo_           [structure]     MATLAB's structure containing the results (oo_).
% - options_      [structure]     MATLAB's structure describing the current options (options_).
% - bayestopt_    [structure]     describing the priors
% - estim_params_ [structure]     characterizing parameters to be estimated
% - dataset_      [structure]     the dataset after required transformation
% - dataset_info  [structure]     Various information about the dataset (descriptive statistics and missing observations)
% - store_linear_smoother [Boolean]    indicator whether linear smoother results should be written to oo_
%
% OUTPUTS
% - alphahat      [double]  (m*T) matrix, smoothed endogenous variables (a_{t|T})  (decision-rule order)
% - etahat        [double]  (r*T) matrix, smoothed structural shocks (r>=n is the number of shocks).
% - epsilonhat    [double]  (n*T) matrix, smoothed measurement errors.
% - ahat0         [double]  (m*T) matrix, updated (endogenous) variables (a_{t|t}) (decision-rule order)
% - SteadyState   [double]  (m*1) vector specifying the steady state level of each endogenous variable (declaration order)
% - trend_coeff   [double]  (n*1) vector, parameters specifying the slope of the trend associated to each observed variable.
% - aKK           [double]  (K,n,T+K) array, k (k=1,...,K) steps ahead
%                                   filtered (endogenous) variables  (decision-rule order)
% - T0 and R0     [double]  Matrices defining the state equation (T is the (m*m) transition matrix).
% - P:            (m*m*(T+1)) 3D array of one-step ahead forecast error variance
%                       matrices (decision-rule order)
% - PKK           (K*m*m*(T+K)) 4D array of k-step ahead forecast error variance
%                       matrices (meaningless for periods 1:d) (decision-rule order)
% - decomp        (K*m*r*(T+K)) 4D array of shock decomposition of k-step ahead
%                       filtered variables (decision-rule order)
% - Trend         [double] (n*T) pure trend component; stored in options_.varobs order
% - state_uncertainty [double] (K,K,T) array, storing the uncertainty
%                                   about the smoothed state (decision-rule order)
% - M_            [structure] describing the model
% - oo_           [structure] storing the results
% - options_      [structure] describing the options
% - bayestopt_    [structure] describing the priors
% - alphahat0     [double]  (m*1) array, smoothed endogenous variables in period 0 (a_{0|T})  (decision-rule order)
% - state_uncertainty0 [double] (K,K,1) array, storing the uncertainty in period 0

% Copyright © 2021-2026 Dynare Team
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

smoother_field_list = {'SmoothedVariables', 'UpdatedVariables', 'SmoothedShocks'};

% get dataset information

if not(isempty(xparam1))
    M_ = set_all_parameters(xparam1,estim_params_,M_);
end
if options_.occbin.smoother.linear_smoother
    %% linear smoother
    options_.occbin.smoother.status=false;
    [alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,T0,R0,P,PK,decomp,Trend,state_uncertainty,oo_.dr,mf,alphahat0,state_uncertainty0] = ...
        DsgeSmoother(xparam1,dataset_,dataset_info,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,options_,bayestopt_,estim_params_);
    bayestopt_.mf=mf;
    if store_linear_smoother
        tmp_smoother=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,...
            aK,P,PK,decomp,Trend,state_uncertainty,alphahat0,state_uncertainty0);
        for jf=1:length(smoother_field_list)
            oo_.occbin.linear_smoother.(smoother_field_list{jf}) = tmp_smoother.(smoother_field_list{jf});
        end
    end
    oo_.occbin.linear_smoother.alphahat=alphahat;
    oo_.occbin.linear_smoother.etahat=etahat;
    oo_.occbin.linear_smoother.epsilonhat=epsilonhat;
    oo_.occbin.linear_smoother.ahat=ahat;
    oo_.occbin.linear_smoother.SteadyState=SteadyState;
    oo_.occbin.linear_smoother.trend_coeff=trend_coeff;
    oo_.occbin.linear_smoother.aK=aK;
    oo_.occbin.linear_smoother.T0=T0;
    oo_.occbin.linear_smoother.R0=R0;
    oo_.occbin.linear_smoother.decomp=decomp;
    oo_.occbin.linear_smoother.alphahat0=alphahat0;
    oo_.occbin.linear_smoother.state_uncertainty0=state_uncertainty0;

    disp_verbose('OccBin: linear smoother done.',options_.verbosity)
    options_.occbin.smoother.status=true;
end
% if init_mode
%% keep these commented lines for the moment, should some nasty problem in the future call back for this option
%     shocks1     = etahat(:,2:end)';
%     clear regime_history ;
%
%     for ii = 1:gend-1
%
%         fprintf('Period number %d out of %d \n',ii,gend-1)
%         opts_simul.SHOCKS = shocks1(ii,:);
%         opts_simul.endo_init = alphahat(oo_.dr.inv_order_var,ii);
%         out =  occbin.runsim_fn(opts_simul, M_, oo_, options_);
%         regime_history(ii) = out.regime_history;
%     end
%
%     [TT,RR] = dynare_resolve(M_,options_,oo_);
%     CC = zeros([size(RR,1),gend ]);
% else

opts_simul = options_.occbin.simul;
opts_simul.curb_retrench = options_.occbin.smoother.curb_retrench;
opts_simul.maxit = options_.occbin.smoother.maxit;
opts_simul.waitbar = options_.occbin.smoother.waitbar;
opts_simul.periods = options_.occbin.smoother.periods;
opts_simul.check_ahead_periods = options_.occbin.smoother.check_ahead_periods;
opts_simul.max_check_ahead_periods = options_.occbin.smoother.max_check_ahead_periods;
opts_simul.periodic_solution = options_.occbin.smoother.periodic_solution;
opts_simul.full_output = options_.occbin.smoother.full_output;
opts_simul.piecewise_only = options_.occbin.smoother.piecewise_only;
occbin_options = struct();

occbin_options.first_period_occbin_update = options_.occbin.smoother.first_period_occbin_update;
occbin_options.opts_simul = opts_simul; % this builds the opts_simul options field needed by occbin.solver
occbin_options.opts_regime.binding_indicator = options_.occbin.smoother.init_binding_indicator;
occbin_options.opts_regime.regime_history=options_.occbin.smoother.init_regime_history;

options_.noprint = true;

is_realtime_smoother_successful = true;

[alphahat,etahat,epsilonhat,ahat,SteadyState,trend_coeff,aK,T0,R0,P,PK,decomp,Trend,state_uncertainty,oo_.dr,mf,alphahat0,state_uncertainty0,~,error_indicator,oo_.occbin.smoother.regime_history] = ...
    DsgeSmoother(xparam1,dataset_,dataset_info,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,options_,bayestopt_,estim_params_,occbin_options);%     T1=TT;
bayestopt_.mf=mf;

if error_indicator(1) || isempty(alphahat0)
    is_realtime_smoother_successful = false;
    if ~options_.occbin.smoother.linear_smoother || ~store_linear_smoother %make sure linear smoother results are set before using them
        options_.occbin.smoother.status=false;
        [~,etahat,~,~,~,~,~,~,~,~,~,~,~,~,~,~,alphahat0] = ...
            DsgeSmoother(xparam1,dataset_,dataset_info,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,options_,bayestopt_,estim_params_);
         options_.occbin.smoother.status=true;
    else
        etahat= oo_.occbin.linear_smoother.etahat;
        alphahat0= oo_.occbin.linear_smoother.alphahat0;
        state_uncertainty0 = oo_.occbin.linear_smoother.state_uncertainty0;
    end
    base_regime = struct();
    if M_.occbin.constraint_nbr==1
        base_regime.regime = 0;
        base_regime.regimestart = 1;
    else
        base_regime.regime1 = 0;
        base_regime.regimestart1 = 1;
        base_regime.regime2 = 0;
        base_regime.regimestart2 = 1;
    end
    oo_.occbin.smoother.regime_history = [];
    for jper=1:size(alphahat,2)+1
        if jper == 1
            oo_.occbin.smoother.regime_history = base_regime;
        else
            oo_.occbin.smoother.regime_history(jper) = base_regime;
        end
    end
end

oo_.occbin.smoother.realtime_regime_history = oo_.occbin.smoother.regime_history;
regime_history = oo_.occbin.smoother.regime_history;
opts_regime.regime_history = oo_.occbin.smoother.regime_history;

ahat0 = ahat;
aKK=aK;
PKK=PK;
clear aK PK;

opts_regime.binding_indicator=[];
regime_history0 = regime_history;

disp_verbose('OccBin smoother iteration 1.',options_.verbosity)
opts_simul.SHOCKS = [etahat(:,1:end)'; zeros(1,M_.exo_nbr)];
opts_simul.exo_pos = 1:M_.exo_nbr;
opts_simul.endo_init = alphahat0(oo_.dr.inv_order_var,1);
opts_simul.init_regime=regime_history; % use realtime regime for guess, to avoid multiple solution issues!
opts_simul.periods = size(opts_simul.SHOCKS,1);
options_.occbin.simul=opts_simul;
occbin_smoother_debug=options_.occbin.smoother.debug;
[~, out, ss] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
oo_.occbin.smoother.error_flag=0;
is_last_simulation_converged = true;
if out.error_flag
    is_last_simulation_converged = false;
    disp_verbose('OccBin smoother:: simulation within smoother did not converge.',options_.verbosity)
    oo_.occbin.smoother.error_flag=321;
    % store regimes consistent with the last smoother run
    out.regime_history = regime_history ;
elseif not(isequal(out.regime_history(1:options_.occbin.likelihood.first_period_binding_regime_allowed-1),regime_history(1:options_.occbin.likelihood.first_period_binding_regime_allowed-1)))
    disp_verbose('Occbin smoother:: simulation violates first_period_binding_regime_allowed.',~options_.noprint)
    oo_.occbin.smoother.error_flag=322;
    % store regimes consistent with the last smoother run
    out.regime_history = regime_history ;
end
regime_history = out.regime_history;
if options_.smoother_redux
    occbin_options.opts_simul.restrict_state_space =1;
    [T0,R0] = dynare_resolve(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    oo_.occbin.linear_smoother.T0=T0;
    oo_.occbin.linear_smoother.R0=R0;
end
TT = ss.T(oo_.dr.order_var,oo_.dr.order_var,:);
RR = ss.R(oo_.dr.order_var,:,:);
CC = ss.C(oo_.dr.order_var,:);

opts_regime.regime_history = regime_history;
opts_regime.binding_indicator = [];
[TT, RR, CC, regime_history] = occbin.check_regimes(TT, RR, CC, opts_regime, M_, options_ , oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
is_changed = ~isequal(regime_history0,regime_history);
if isempty(regime_history0)
    regime_history0 = regime_history;
end
iter=1;

is_periodic = 0;
is_changed_start = 0;
maxiter = options_.occbin.smoother.max_number_of_iterations;
occbin_smoother_fast = options_.occbin.smoother.fast;

sto_alphahat=alphahat;
sto_etahat={etahat};
sto_CC = CC;
sto_RR = RR;
sto_TT = TT;
sto_eee=NaN(size(TT,1),size(TT,3));
for k=1:size(TT,3)
    sto_eee(:,k) = eig(TT(:,:,k));
end

is_conditional_smoother_converged = false;
is_realtime_smoother_converged = true;

if (is_changed || oo_.occbin.smoother.error_flag)
    is_realtime_smoother_converged = false;
end
if is_realtime_smoother_successful
    sto_realtime.alphahat=alphahat;
    sto_realtime.etahat=etahat;
    sto_realtime.epsilonhat=epsilonhat;
    sto_realtime.SteadyState=SteadyState;
    sto_realtime.trend_coeff=trend_coeff;
    sto_realtime.T0=T0;
    sto_realtime.R0=R0;
    sto_realtime.P=P;
    sto_realtime.decomp=decomp;
    sto_realtime.Trend=Trend;
    sto_realtime.state_uncertainty=state_uncertainty;
    sto_realtime.dr=oo_.dr;
    sto_realtime.alphahat0=alphahat0;
    sto_realtime.state_uncertainty0=state_uncertainty0;
    sto_realtime.opts_simul=opts_simul;
    sto_realtime.regime_history=regime_history;
end

if ~is_realtime_smoother_converged && not(options_.lik_init==2 && options_.Harvey_scale_factor==0)
    % try conditional smoother, where state uncertainty is smaller [unless too many shocks in excess]
    disp_verbose(sprintf('OccBin: try conditional smoother iteration.'),options_.verbosity)
    occbin_options.opts_regime.regime_history=regime_history;
    %%%% resort to conditional smoother
    alphahat1 = alphahat0(oo_.dr.inv_order_var)+SteadyState;
    alphahat1 = alphahat1(oo_.dr.state_var);
    M_local = M_;
    % set direct assignment of initial states
    M_local.endo_initial_state.status = true;
    M_local.endo_initial_state.values = zeros(M_.endo_nbr,1);
    opts_local1 = options_;
    opts_local1.lik_init = 2;
    opts_local1.Harvey_scale_factor = 0;
    M_local.endo_initial_state.values(oo_.dr.state_var) = alphahat1;
    occbin_options.first_period_occbin_update = 1;
    [c.alphahat,c.etahat,c.epsilonhat,~,~,~,~,c.T0,c.R0,c.P,~,~,~,~,~,~,c.alphahat0,~,~,c.error_indicator,c.regime_history] = ...
        DsgeSmoother(xparam1,dataset_,dataset_info,M_local,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,opts_local1,bayestopt_,estim_params_,occbin_options);%     T1=TT;
    if c.error_indicator(1)
        disp_verbose('OccBin smoother:: there was an error in running conditional smoother.',options_.verbosity)
        oo_.occbin.smoother.error_flag=322;
    else
        %%%% compute state uncertainty consistent with conditional smoother
        %%%% regime sequence
        regime_history = c.regime_history;
        opts_simul.SHOCKS = [c.etahat(:,1:end)'; zeros(1,M_.exo_nbr)];
        opts_simul.endo_init = alphahat0(oo_.dr.inv_order_var,1);
        options_.occbin.simul=opts_simul;
        out0 = out;
        ss0=ss;
        [~, out, ss] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
        oo_.occbin.smoother.error_flag=0;
        if out.error_flag
            disp_verbose('OccBin smoother:: simulation within conditional smoother did not converge.',options_.verbosity)
            oo_.occbin.smoother.error_flag=321;
            % use regimes consistent with the last smoother run
            out = out0;
            out.regime_history = regime_history;
            ss=ss0;
        elseif not(isequal(out.regime_history(1:options_.occbin.likelihood.first_period_binding_regime_allowed-1),regime_history(1:options_.occbin.likelihood.first_period_binding_regime_allowed-1)))
            disp_verbose('Occbin smoother:: simulation violates first_period_binding_regime_allowed.',~options_.noprint)
            oo_.occbin.smoother.error_flag=322;
            % use regimes consistent with the last smoother run
            out.regime_history = regime_history ;
        end
        c.TT = ss.T(oo_.dr.order_var,oo_.dr.order_var,:);
        c.RR = ss.R(oo_.dr.order_var,:,:);
        c.CC = ss.C(oo_.dr.order_var,:);

        opts_regime.regime_history = out.regime_history;
        [c.TT, c.RR, c.CC, regime_history] = occbin.check_regimes(c.TT, c.RR, c.CC, opts_regime, M_, options_ , oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
        is_changed = ~isequal(c.regime_history,regime_history);
        if not(is_changed || oo_.occbin.smoother.error_flag)
            is_conditional_smoother_converged = true;
            occbin_options.opts_regime.regime_history=c.regime_history;
            occbin_options.first_period_occbin_update = inf;
            % compute uncertainty consistent with sequence of regimes of conditional smoother
            CC = c.CC;
            TT = c.TT;
            RR = c.RR;
            [~,~,~,~,~,~,~,~,~,P,~,decomp,~,state_uncertainty,~,~,~,state_uncertainty0]...
                = DsgeSmoother(xparam1,dataset_,dataset_info,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,options_,bayestopt_,estim_params_,occbin_options,TT,RR,CC);
            alphahat = c.alphahat;
            etahat = c.etahat;
            epsilonhat = c.epsilonhat;
            T0 = c.T0;
            R0 = c.R0;
            alphahat0 = c.alphahat0;
            oo_.occbin.smoother.regime_history = c.regime_history;
        else
            disp_verbose(sprintf('OccBin: conditional smoother did not converge.'),options_.verbosity)
        end
    end
end


occbin_options.first_period_occbin_update = inf;

% try consistent estimation of shocks and regime sequence
while is_changed && maxiter>iter && ~is_periodic && is_last_simulation_converged
    if iter==1
        regime_history = sto_realtime.regime_history;
        oo_.occbin.smoother.error_flag=0;
    end
    iter=iter+1;
    disp_verbose(sprintf('OccBin smoother iteration %u.', iter),options_.verbosity)
    occbin_options.opts_regime.regime_history=regime_history;
    [alphahat,etahat,epsilonhat,~,SteadyState,trend_coeff,~,T0,R0,P,~,decomp,Trend,state_uncertainty,oo_.dr,~,alphahat0,state_uncertainty0]...
        = DsgeSmoother(xparam1,dataset_,dataset_info,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,options_,bayestopt_,estim_params_,occbin_options,TT,RR,CC);
    sto_etahat(iter)={etahat};
    regime_history0(iter,:) = regime_history;
    if occbin_smoother_debug
        save('Occbin_smoother_debug_regime_history','regime_history0');
    end

    sto_CC = CC;
    sto_RR = RR;
    sto_TT = TT;

    opts_simul.SHOCKS = [etahat(:,1:end)'; zeros(1,M_.exo_nbr)];
    opts_simul.endo_init = alphahat0(oo_.dr.inv_order_var,1);
    options_.occbin.simul=opts_simul;
    ss0=ss;
    [~, out, ss] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    if out.error_flag
        disp_verbose('OccBin smoother:: simulation within smoother did not converge.',options_.verbosity)
        oo_.occbin.smoother.error_flag=321;
        % use regimes consistent with the last smoother run
        out.regime_history = regime_history;
        break
    elseif not(isequal(out.regime_history(1:options_.occbin.likelihood.first_period_binding_regime_allowed-1),regime_history(1:options_.occbin.likelihood.first_period_binding_regime_allowed-1)))
        disp_verbose('Occbin smoother:: simulation violates first_period_binding_regime_allowed.',~options_.noprint)
        oo_.occbin.smoother.error_flag=322;
        % use regimes consistent with the last smoother run
        out.regime_history = regime_history ;
        break
    end
    regime_history = out.regime_history;
    TT = ss.T(oo_.dr.order_var,oo_.dr.order_var,:);
    RR = ss.R(oo_.dr.order_var,:,:);
    CC = ss.C(oo_.dr.order_var,:);

    opts_regime.regime_history = regime_history;
    [TT, RR, CC, regime_history] = occbin.check_regimes(TT, RR, CC, opts_regime, M_, options_ , oo_.dr, oo_.steady_state, oo_.exo_steady_state, oo_.exo_det_steady_state);
    is_changed = ~isequal(regime_history0(iter,:),regime_history);
    isdiff_regime=NaN(size(regime_history0,2),M_.occbin.constraint_nbr);
    isdiff_start=NaN(size(isdiff_regime));
    isdiff_=NaN(size(isdiff_regime));
    if M_.occbin.constraint_nbr==2
        for k=1:size(regime_history0,2)
            isdiff_regime(k,1) = ~isequal(regime_history0(end,k).regime1,regime_history(k).regime1);
            isdiff_start(k,1) = ~isequal(regime_history0(end,k).regimestart1,regime_history(k).regimestart1);
            isdiff_(k,1) =  isdiff_regime(k,1) ||  isdiff_start(k,1);
            isdiff_regime(k,2) = ~isequal(regime_history0(end,k).regime2,regime_history(k).regime2);
            isdiff_start(k,2) = ~isequal(regime_history0(end,k).regimestart2,regime_history(k).regimestart2);
            isdiff_(k,2) =  isdiff_regime(k,2) ||  isdiff_start(k,2);
        end
        is_changed_regime = any(isdiff_regime(:,1)) || any(isdiff_regime(:,2));
        is_changed_start = any(isdiff_start(:,1)) || any(isdiff_start(:,2));
    else
        for k=1:size(regime_history0,2)
            isdiff_regime(k,1) = ~isequal(regime_history0(end,k).regime,regime_history(k).regime);
            isdiff_start(k,1) = ~isequal(regime_history0(end,k).regimestart,regime_history(k).regimestart);
            isdiff_(k,1) =  isdiff_regime(k,1) ||  isdiff_start(k,1);
        end
        is_changed_regime = any(isdiff_regime(:,1));
        is_changed_start = any(isdiff_start(:,1));
    end
    if occbin_smoother_fast
        is_changed = is_changed_regime;
    end
    if iter>1
        for kiter=1:iter-1
            is_tmp = isequal(regime_history0(kiter,:),regime_history);
            if is_tmp
                break
            end
        end
        is_periodic = (is_changed && is_tmp);
    end

    if is_changed
        eee=NaN(size(TT,1),size(TT,3));
        for k=1:size(TT,3)
            eee(:,k) = eig(TT(:,:,k));
        end
    end

    if occbin_smoother_debug || (is_periodic && ~options_.noprint)
        regime_ = cell(0);
        regime_new = regime_;
        start_ = regime_;
        start_new = regime_;
        if M_.occbin.constraint_nbr==2
            indx_init_1 = find(isdiff_(:,1));
            if ~isempty(indx_init_1)
                qq={regime_history0(end,indx_init_1).regime1}';
                for j=1:length(qq), regime_(j,1) = {int2str(qq{j})}; end
                qq={regime_history0(end,indx_init_1).regimestart1}';
                for j=1:length(qq), start_(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_1).regime1}';
                for j=1:length(qq), regime_new(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_1).regimestart1}';
                for j=1:length(qq), start_new(j,1) = {int2str(qq{j})}; end
                disp('Time points where regime 1 differs')
                if ~isoctave
                    disp(table(indx_init_1, regime_, start_, regime_new, start_new))
                else % The table() function is not implemented in Octave, print something more or less equivalent (though much less readable)
                    disp(vertcat({'indx_init_1', 'regime_', 'start_', 'regime_new', 'start_new'}, ...
                                 horzcat(num2cell(indx_init_1), regime_, start_, regime_new, start_new)))
                end
            end

            indx_init_2 = find(isdiff_(:,2));
            if ~isempty(indx_init_2)
                regime_ = cell(0);
                regime_new = regime_;
                start_ = regime_;
                start_new = regime_;
                qq={regime_history0(end,indx_init_2).regime2}';
                for j=1:length(qq), regime_(j,1) = {int2str(qq{j})}; end
                qq={regime_history0(end,indx_init_2).regimestart2}';
                for j=1:length(qq), start_(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_2).regime2}';
                for j=1:length(qq), regime_new(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_2).regimestart2}';
                for j=1:length(qq), start_new(j,1) = {int2str(qq{j})}; end
                disp('Time points where regime 2 differs ')
                if ~isoctave
                    disp(table(indx_init_2, regime_, start_, regime_new, start_new))
                else % The table() function is not implemented in Octave, print something more or less equivalent (though much less readable)
                    disp(vertcat({'indx_init_2', 'regime_', 'start_', 'regime_new', 'start_new'}, ...
                                 horzcat(num2cell(indx_init_2), regime_, start_, regime_new, start_new)))
                end
            end
        else
            indx_init_1 = find(isdiff_(:,1));
            if ~isempty(indx_init_1)
                qq={regime_history0(end,indx_init_1).regime}';
                for j=1:length(qq), regime_(j,1) = {int2str(qq{j})}; end
                qq={regime_history0(end,indx_init_1).regimestart}';
                for j=1:length(qq), start_(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_1).regime}';
                for j=1:length(qq), regime_new(j,1) = {int2str(qq{j})}; end
                qq={regime_history(indx_init_1).regimestart}';
                for j=1:length(qq), start_new(j,1) = {int2str(qq{j})}; end
                disp('Time points where regime differs')
                if ~isoctave
                    disp(table(indx_init_1, regime_, start_, regime_new, start_new))
                else % The table() function is not implemented in Octave, print something more or less equivalent (though much less readable)
                    disp(vertcat({'indx_init_1', 'regime_', 'start_', 'regime_new', 'start_new'}, ...
                                 horzcat(num2cell(indx_init_1), regime_, start_, regime_new, start_new)))
                end
            end

        end
    end
    if is_changed
        sto_alphahat=alphahat;
        sto_eee = eee;
    end
end
regime_history0(max(iter+1,1),:) = regime_history;
oo_.occbin.smoother.regime_history=regime_history0(end,:);
oo_.occbin.smoother.regime_history_iter=regime_history0;
if occbin_smoother_debug
    save('Occbin_smoother_debug_regime_history','regime_history0')
end

is_smoother_converged=false;
oo_.occbin.smoother.warning_flag = 0;

if ( (maxiter==iter || ~is_last_simulation_converged) && is_changed) || is_periodic || oo_.occbin.smoother.error_flag
    disp_verbose('occbin.DSGE_smoother: smoother did not converge.',options_.verbosity)
    disp_verbose('occbin.DSGE_smoother: The algorithm did not reach a fixed point for the smoothed regimes.',options_.verbosity)
    if is_periodic
        is_smoother_converged=true;
        oo_.occbin.smoother.error_flag=0;
        oo_.occbin.smoother.warning_flag = 321;
        disp_verbose('occbin.DSGE_smoother: For the periods indicated above, regimes loops between the "regime_" and the "regime_new_" pattern displayed above.',options_.verbosity)
        disp_verbose('occbin.DSGE_smoother: We provide smoothed shocks consistent with "regime_" in oo_.',options_.verbosity)
    elseif ~is_realtime_smoother_successful
        disp_verbose('occbin.DSGE_smoother: The realtime smoother was not successful either.',options_.verbosity)
        disp_verbose('occbin.DSGE_smoother: The respective fields in oo_ will be left empty.',options_.verbosity)
        oo_.occbin.smoother=[];
        oo_.occbin.smoother.error_flag=322;
    end
else
    disp_verbose('occbin.DSGE_smoother: smoother converged.',options_.verbosity)
    is_smoother_converged=true;
    oo_.occbin.smoother.error_flag=0;
    if ~is_conditional_smoother_converged && ~is_realtime_smoother_converged
        oo_.occbin.smoother.warning_flag = 320;
        disp_verbose('occbin.DSGE_smoother: WARNING: algorithm converged to a different regime wrt realtime regime',options_.verbosity)
        disp_verbose('occbin.DSGE_smoother: WARNING: this usually indicates that for some period there is no unique regime dominance.',options_.verbosity)
    elseif occbin_smoother_fast && is_changed_start
        oo_.occbin.smoother.warning_flag = 323;
        disp_verbose('occbin.DSGE_smoother: WARNING: fast algo is used, regime duration was not forced to converge',options_.verbosity)
        disp_verbose('occbin.DSGE_smoother: WARNING: this usually indicates that for some period there is no unique regime dominance.',options_.verbosity)
    end
end

if is_realtime_smoother_successful && ~is_conditional_smoother_converged && ~is_realtime_smoother_converged && ~is_smoother_converged
    % fallback solution
    disp_verbose('occbin.DSGE_smoother: We provide results consistent with realtime regime',options_.verbosity)
    disp_verbose('occbin.DSGE_smoother: WARNING: this usually indicates that for some period there is no unique regime dominance.',options_.verbosity)
    oo_.occbin.smoother.warning_flag = 322;
    oo_.occbin.smoother.error_flag = 0;
    options_.occbin.simul=sto_realtime.opts_simul;
    options_.occbin.simul.maxit=1;
    [~, out] = occbin.solver(M_,options_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    oo_.occbin.smoother.regime_history = oo_.occbin.smoother.realtime_regime_history;
    is_changed=false;
    alphahat=sto_realtime.alphahat;
    etahat=sto_realtime.etahat;
    epsilonhat=sto_realtime.epsilonhat;
    SteadyState=sto_realtime.SteadyState;
    trend_coeff=sto_realtime.trend_coeff;
    T0=sto_realtime.T0;
    R0=sto_realtime.R0;
    P=sto_realtime.P;
    decomp=sto_realtime.decomp;
    Trend=sto_realtime.Trend;
    state_uncertainty=sto_realtime.state_uncertainty;
    oo_.dr=sto_realtime.dr;
    alphahat0=sto_realtime.alphahat0;
    state_uncertainty0=sto_realtime.state_uncertainty0;
end


if (~is_changed || (occbin_smoother_debug && iter>1)) && store_linear_smoother
    if is_changed
        % this can happen when realtime smoother did not work and
        % iterations done starting from linear smoother did not converge
        CC = sto_CC;
        RR = sto_RR;
        TT = sto_TT;
        oo_.occbin.smoother.regime_history=regime_history0(end-1,:);
    end
    if options_.occbin.smoother.store_results
        tmp_smoother=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,alphahat,etahat,epsilonhat,ahat0,SteadyState,trend_coeff,aKK,P,PKK,decomp,Trend,state_uncertainty,alphahat0,state_uncertainty0);
        for jf=1:length(smoother_field_list)
            oo_.occbin.smoother.(smoother_field_list{jf}) = tmp_smoother.(smoother_field_list{jf});
        end
        oo_.occbin.smoother.alphahat=alphahat;
        oo_.occbin.smoother.etahat=etahat;
        oo_.occbin.smoother.epsilonhat=epsilonhat;
        oo_.occbin.smoother.ahat=ahat0;
        oo_.occbin.smoother.SteadyState=SteadyState;
        oo_.occbin.smoother.trend_coeff=trend_coeff;
        oo_.occbin.smoother.aK=aKK;
        oo_.occbin.smoother.T0=TT;
        oo_.occbin.smoother.R0=RR;
        oo_.occbin.smoother.C0=CC;
        if isfield(out,'piecewise')
            oo_.occbin.smoother.simul.piecewise = out.piecewise(1:end-1,:);
        end
        if ~options_.occbin.simul.piecewise_only &&  isfield(out,'linear')
            oo_.occbin.smoother.simul.linear = out.linear(1:end-1,:);
        end
    end
    if options_.occbin.smoother.plot
        GraphDirectoryName = CheckPath('graphs',M_.fname);
        latexFolder = CheckPath('latex',M_.dname);

        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([latexFolder filesep M_.fname '_OccBin_smoother_plots.tex'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by occbin.DSGE_smoother.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
            fprintf(fidTeX,' \n');
        end
        j1=0;
        ifig=0;
        for j=1:M_.exo_nbr
            if max(abs(oo_.occbin.smoother.etahat(j,:)))>1.e-8
                j1=j1+1;
                if mod(j1,9)==1
                    hh_fig = dyn_figure(options_.nodisplay,'name','OccBin smoothed shocks');
                    ifig=ifig+1;
                    isub=0;
                end
                isub=isub+1;
                subplot(3,3,isub)
                if  options_.occbin.smoother.linear_smoother
                    plot(oo_.occbin.linear_smoother.etahat(j,:)','linewidth',2)
                    hold on,
                end
                plot(oo_.occbin.smoother.etahat(j,:)','r--','linewidth',2)
                hold on, plot([0 options_.nobs],[0 0],'k--')
                set(gca,'xlim',[0 options_.nobs])
                if options_.TeX
                    title(['$' M_.exo_names_tex{j,:} '$'],'interpreter','latex')
                else
                    title(M_.exo_names{j,:},'interpreter','none')
                end

                if mod(j1,9)==0
                    if  options_.occbin.smoother.linear_smoother
                        annotation('textbox', [0.1,0,0.35,0.05],'String', 'Linear','Color','Blue','horizontalalignment','center','interpreter','none');
                    end
                    annotation('textbox', [0.55,0,0.35,0.05],'String', 'Piecewise','Color','Red','horizontalalignment','center','interpreter','none');
                    dyn_saveas(gcf,[GraphDirectoryName filesep M_.fname,'_smoothedshocks_occbin',int2str(ifig)],options_.nodisplay,options_.graph_format);
                    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                        % TeX eps loader file
                        fprintf(fidTeX,'\\begin{figure}[H]\n');
                        fprintf(fidTeX,'\\centering \n');
                        fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_smoothedshocks_occbin%s}\n',options_.figures.textwidth*min(j1/3,1),[GraphDirectoryName '/' M_.fname],int2str(ifig)); % don't use filesep as it will create issues with LaTeX on Windows
                        fprintf(fidTeX,'\\caption{OccBin smoothed shocks.}');
                        fprintf(fidTeX,'\\label{Fig:smoothedshocks_occbin:%s}\n',int2str(ifig));
                        fprintf(fidTeX,'\\end{figure}\n');
                        fprintf(fidTeX,' \n');
                    end
                end
            end
        end

        if mod(j1,9)~=0 && j==M_.exo_nbr
                annotation('textbox', [0.1,0,0.35,0.05],'String', 'Linear','Color','Blue','horizontalalignment','center','interpreter','none');
                annotation('textbox', [0.55,0,0.35,0.05],'String', 'Piecewise','Color','Red','horizontalalignment','center','interpreter','none');
                dyn_saveas(hh_fig,[GraphDirectoryName filesep M_.fname,'_smoothedshocks_occbin',int2str(ifig)],options_.nodisplay,options_.graph_format);
            if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                % TeX eps loader file
                fprintf(fidTeX,'\\begin{figure}[H]\n');
                fprintf(fidTeX,'\\centering \n');
                fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s_smoothedshocks_occbin%s}\n',options_.figures.textwidth*min(j1/3,1),[GraphDirectoryName '/' M_.fname],int2str(ifig)); % don't use filesep as it will create issues with LaTeX on Windows
                fprintf(fidTeX,'\\caption{OccBin smoothed shocks.}');
                fprintf(fidTeX,'\\label{Fig:smoothedshocks_occbin:%s}\n',int2str(ifig));
                fprintf(fidTeX,'\\end{figure}\n');
                fprintf(fidTeX,' \n');
            end
        end
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fclose(fidTeX);
        end
    end
end
