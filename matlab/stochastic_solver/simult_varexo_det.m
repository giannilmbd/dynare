function [y_,yf_CI,yf_CI_ME]=simult_varexo_det(y0,ex,ex_det,iorder,var_list,M_,dr,exo_det_steady_state,options_)
%function [y_,yf_CI,yf_CI_ME]=simult_varexo_det(y0,ex,ex_det, iorder,var_list,M_,dr,exo_det_steady_state,options_)
%
% Simulates a stochastic model in the presence of deterministic exogenous
% shocks; allows bootstrapping standard errors for forecasting
%
% INPUTS:
%    y0:        initial values, of length M_.maximum_lag
%    ex:        matrix of stochastic exogenous shocks, starting at period 1
%    ex_det:    matrix of deterministic exogenous shocks, starting at period 1-M_.maximum_lag
%    iorder:    order of approximation
%    var_list:  list of endogenous variables to simulate
%    M_:        Dynare model structure
%    dr:        decision rules
%    exo_det_steady_state_  steady state of varexo_det
%    options_:  Dynare options structure
%
% OUTPUTS:
%   yf:          mean forecast
%   int_width:   distance between upper bound and
%                mean forecast
%   yf_CI:       upper and low bound for forecast
%   yf_CI_ME:    upper and low bound for forecast when
%                   considering measurement error
%
%
% The forecast horizon is equal to size(ex, 1).
% The condition size(ex,1)+M_.maximum_lag=size(ex_det,1) must be verified
%  for consistency.

% Copyright © 2008-2025 Dynare Team
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

if options_.order>2
    error('simult_varexo_det.m: Forecasting with varexo_det does not support order>2.')
elseif options_.order==2 && options_.pruning
    error('simult_varexo_det.m: Forecasting with varexo_det does not support pruning.')
end
horizon = size(ex,1);
if size(ex_det, 1) ~= horizon+M_.maximum_lag
    error('simult_varexo_det: Size mismatch: number of forecasting periods for stochastic exogenous and deterministic exogenous don''t match')
end
y_ = zeros(size(y0,1),horizon+M_.maximum_lag);
y_(:,1:M_.maximum_lag) = y0;
k1 = M_.maximum_lag:-1:1;

%% get indices of requested variables
nvar = length(var_list);
if nvar == 0
    nvar = M_.endo_nbr;
    ivar = 1:nvar;
else
    ivar=zeros(nvar,1);
    for i=1:nvar
        i_tmp = strmatch(var_list{i}, M_.endo_names, 'exact');
        if isempty(i_tmp)
            error ('simult_varexo_det: the requested variable %s does not exist',var_list{i}) ;
        else
            ivar(i) = i_tmp;
        end
    end
end
if nargout==3 % add measurement error variance uncertainty for observables
    [loc_H, loc_varlist] = ismember(options_.varobs, var_list);
    loc_varlist(loc_varlist==0) = []; %location of observables in variable list
end

fact = norminv((1-options_.forecasts.conf_sig)/2,0,1);
if iorder == 1
    for i = M_.maximum_lag+1: horizon+M_.maximum_lag
        tempx1 = y_(dr.order_var,k1);
        tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,M_.maximum_lag);
        tempx = tempx2(M_.nstatic+(1:M_.nspred));
        y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghx*tempx+dr.ghu*ex(i-M_.maximum_lag,:)';
        for j=1:min(M_.maximum_lag+dr.exo_det_length+1-i,dr.exo_det_length)
            y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*(ex_det(i+j-1,:)'-exo_det_steady_state);
        end

        k1 = k1+1;
    end
    y_=y_(ivar,:); %sort according to requested var_list_
    if nargout>1
        % Define union of requested and state variables
        [var_pos]=union(dr.inv_order_var(dr.state_var),dr.inv_order_var(ivar),'stable');
        % set state_pos to positions of state variables for inner part of
        % variance loop
        [~,state_pos] = ismember(dr.inv_order_var(dr.state_var),var_pos);
        % set var_list_pos to positions of var_list_ for reading out requested
        % variances
        [~,var_list_pos] = ismember(dr.inv_order_var(ivar),var_pos);
        ghx1 = dr.ghx(var_pos,:);
        ghu1 = dr.ghu(var_pos,:);

        sigma_y1 = zeros(size(ghx1,1),size(ghx1,1)); %no initial uncertainty about the states
        var_yf = zeros(nvar,1+horizon); %initialize
        for horizon_iter=1:horizon
        % variance recursion
            sigma_y1 = ghx1*sigma_y1(state_pos,state_pos)*ghx1'+ghu1*M_.Sigma_e*ghu1';
            var_yf(:,1+horizon_iter) = diag(sigma_y1(var_list_pos,var_list_pos))'; % first period is initial condition without uncertainty
        end
        int_width = -fact*sqrt(var_yf);
        yf_CI=cat(3,y_-int_width,y_+int_width);
        if nargout==3 % add measurement error variance uncertainty for observables
            var_yf_ME=var_yf;
            if ~isempty(loc_varlist)
                var_yf_ME(loc_varlist,:) = var_yf(loc_varlist,:)+repmat(diag(M_.H(loc_H,loc_H)), 1,horizon+1);
            end
            int_width_ME = -fact*sqrt(var_yf_ME);
            yf_CI_ME=cat(3,y_-int_width_ME,y_+int_width_ME); %yf is already in correct order
        end
    end
elseif iorder == 2
    y_=simulate_varexo_det_order_2(M_,dr,y_,k1,ex,ex_det,exo_det_steady_state,horizon);
    if nargout>1
        %% generate confidence bands via simulation
        % eliminate shocks with 0 variance
        i_exo_var = setdiff(1:M_.exo_nbr,find(diag(M_.Sigma_e) == 0));
        nxs = length(i_exo_var);
        exo_simul = zeros(horizon,M_.exo_nbr);
        chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));
        y_simul=NaN(M_.endo_nbr,M_.maximum_lag+horizon,options_.forecast_replic);
        for forecast_iter=1:options_.forecast_replic
            if ~isempty(M_.Sigma_e)
                % we fill the shocks row wise to have the same values
                % independently of the length of the simulation
                exo_simul(:,i_exo_var) = randn(nxs,horizon)'*chol_S;
            end
            y_simul(:,:,forecast_iter) = simulate_varexo_det_order_2(M_,dr,y_(:,1:M_.maximum_lag),k1,exo_simul,ex_det,exo_det_steady_state,horizon);
        end
        yf_CI=quantile(y_simul(ivar,:,:),[(1-options_.forecasts.conf_sig)/2 1-(1-options_.forecasts.conf_sig)/2],3);
    end
    if nargout==3 % add measurement error variance uncertainty for observables
        if ~isempty(loc_varlist)
            % eliminate shocks with 0 variance
            i_endo_var = setdiff(1:length(options_.varobs),find(diag(M_.H) == 0));
            endo_simul = zeros(horizon+1,length(options_.varobs));
            chol_S = chol(M_.H(i_endo_var,i_endo_var));
            for forecast_iter=1:options_.forecast_replic
                if ~isempty(M_.H)
                    % we fill the shocks row wise to have the same values
                    % independently of the length of the simulation
                    endo_simul(:,i_endo_var,forecast_iter) = randn(length(i_endo_var),horizon+1)'*chol_S;
                end
            end
            y_simul(options_.varobs_id,:,:) = y_simul(options_.varobs_id,:,:)+permute(endo_simul,[2,1,3]);%overwrite old variable, not needed anymore
        else
            %no reason to add ME as no variable was requested
        end
        yf_CI_ME=quantile(y_simul(ivar,:,:),[(1-options_.forecasts.conf_sig)/2 1-(1-options_.forecasts.conf_sig)/2],3);
    end
else
    error('simultxdet.m: order>2 not supported.')
end


function y_=simulate_varexo_det_order_2(M_,dr,y_,k1,exo,exo_det,exo_det_steady_state,horizon)
for i = M_.maximum_lag+1:horizon+M_.maximum_lag
    tempx1 = y_(dr.order_var,k1);
    tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,M_.maximum_lag);
    tempx = tempx2(M_.nstatic+(1:M_.nspred));
    tempu = exo(i-M_.maximum_lag,:)';
    tempuu = kron(tempu,tempu);
    tempxx = kron(tempx,tempx);
    tempxu = kron(tempx,tempu);
    y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghs2/2+dr.ghx*tempx+ ...
        dr.ghu*tempu+0.5*(dr.ghxx*tempxx+dr.ghuu*tempuu)+dr.ghxu* ...
        tempxu;
    for j=1:min(M_.maximum_lag+dr.exo_det_length+1-i,dr.exo_det_length)
        tempud = exo_det(i+j-1,:)'-exo_det_steady_state;
        tempudud = kron(tempud,tempud);
        tempxud = kron(tempx,tempud);
        tempuud = kron(tempu,tempud);
        y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*tempud + ...
            dr.ghxud{j}*tempxud + dr.ghuud{j}*tempuud + ...
            0.5*dr.ghudud{j,j}*tempudud;
        for k=1:j-1
            tempudk = exo_det(i+k-1,:)'-exo_det_steady_state;
            tempududk = kron(tempudk,tempud);
            y_(dr.order_var,i) = y_(dr.order_var,i) + ...
                dr.ghudud{k,j}*tempududk;
        end
    end
    k1 = k1+1;
end