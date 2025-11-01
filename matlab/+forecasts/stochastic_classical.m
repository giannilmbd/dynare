function [yf,yf_CI,yf_CI_ME]=stochastic_classical(dr,y0,horizon,var_list,M_,options_)
% function [yf,yf_CI,yf_CI_ME]=stochastic_classical(dr,y0,horizon,var_list,M_,options_)
%   computes mean forecast for a given value of the parameters
%   computes also confidence band for the forecast
%
% INPUTS:
%   dr:          structure containing decision rules
%   y0:          initial values
%   horizon:     nbr of periods to forecast
%   var_list:    list of variables (character matrix)
%   M_:          Dynare model structure
%   options_:    Dynare options structure

% OUTPUTS:
%   yf:          mean forecast, ivar by (horizon+1),includes starting value
%   yf_CI:       upper and low bound for forecast
%   yf_CI_ME:    upper and low bound for forecast when
%                   considering measurement error
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2003-2024 Dynare Team
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

%% get variable list and indices
if isempty(var_list)
    var_list = M_.endo_names(1:M_.orig_endo_nbr);
end
nvar = length(var_list);
ivar = zeros(nvar, 1);
for i=1:nvar
    i_tmp = strmatch(var_list{i}, M_.endo_names, 'exact');
    if isempty(i_tmp)
        disp(var_list{i});
        error ('One of the variable specified does not exist') ;
    else
        ivar(i) = i_tmp;
    end
end

fact = norminv((1-options_.forecasts.conf_sig)/2, 0, 1);
%compute point forecast
yf = simult_(M_,options_,y0,dr,zeros(horizon,M_.exo_nbr),options_.order);
yf = yf(ivar,:); %bring into requested order before adding measurement error
if nargout==3 % add measurement error variance uncertainty for observables
    [loc_H, loc_varlist] = ismember(options_.varobs, var_list);
    loc_varlist(loc_varlist==0) = []; %location of observables in variable list
end
if options_.order==1 % use analytical formulas, no need to simulate
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

    %initialize recursion
    sigma_y1 = zeros(size(ghx1,1),size(ghx1,1)); %no initial uncertainty about the states
    var_yf = zeros(nvar,1+horizon); %initialize
    for horizon_iter = 1:horizon
        % variance recursion
        sigma_y1 = ghx1*sigma_y1(state_pos,state_pos)*ghx1'+ghu1*M_.Sigma_e*ghu1';
        var_yf(:,1+horizon_iter) = diag(sigma_y1(var_list_pos,var_list_pos))'; % first period is initial condition without uncertainty
    end

    int_width = -fact*sqrt(var_yf);
    yf_CI=cat(3,yf-int_width,yf+int_width);
    if nargout==3 % add measurement error variance uncertainty for observables
        var_yf_ME=var_yf;
        if ~isempty(loc_varlist)
            var_yf_ME(loc_varlist,:) = var_yf(loc_varlist,:)+repmat(diag(M_.H(loc_H,loc_H)), 1,horizon+1);
        end
        int_width_ME = -fact*sqrt(var_yf_ME);
        yf_CI_ME=cat(3,yf-int_width_ME,yf+int_width_ME); %yf is already in correct order
    end
else %use simulations
    % eliminate shocks with 0 variance
    i_exo_var = setdiff(1:M_.exo_nbr,find(diag(M_.Sigma_e) == 0));
    nxs = length(i_exo_var);
    exo_simul = zeros(horizon,M_.exo_nbr);
    chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));
    yf_simul=NaN(M_.endo_nbr,1+horizon,options_.forecast_replic);
    for forecast_iter=1:options_.forecast_replic
        if ~isempty(M_.Sigma_e)
            % we fill the shocks row wise to have the same values
            % independently of the length of the simulation
            exo_simul(:,i_exo_var) = randn(nxs,horizon)'*chol_S;
        end
        yf_simul(:,:,forecast_iter) = simult_(M_,options_,y0,dr,exo_simul,options_.order);
    end
    yf_CI=quantile(yf_simul(ivar,:,:),[(1-options_.forecasts.conf_sig)/2 1-(1-options_.forecasts.conf_sig)/2],3);
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
            yf_simul(options_.varobs_id,:,:) = yf_simul(options_.varobs_id,:,:)+permute(endo_simul,[2,1,3]);%overwrite old variable, not needed anymore
        else
            %no reason to add ME as no variable was requested
        end
        yf_CI_ME=quantile(yf_simul(ivar,:,:),[(1-options_.forecasts.conf_sig)/2 1-(1-options_.forecasts.conf_sig)/2],3);
    end
end