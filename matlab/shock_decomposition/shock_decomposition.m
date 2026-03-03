function [oo_,M_] = shock_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% shock_decomposition(M_,oo_,options_,varlist,bayestopt_,estim_params_)
% Computes shocks contribution to a simulated trajectory. The field sets
% are oo_.shock_decomposition, oo_.forecast_shock_decomposition.
% Subfields are arrays n_var by nshock+2 by nperiods. The first nshock columns store
% the respective shock contributions, column n+1 stores the role of the initial
% conditions, while column n+2 stores the value of the smoothed variables.
% Both the variables and shocks are stored in the order of declaration,
% i.e. M_.endo_names and M_.exo_names, respectively.
%
% INPUTS
%    M_:          [structure]  Definition of the model
%    oo_:         [structure]  Storage of results
%    options_:    [structure]  Options
%    varlist:     [char]       List of variables
%    bayestopt_:  [structure]  describing the priors
%    estim_params_: [structure] characterizing parameters to be estimated
%
% OUTPUTS
%    oo_:         [structure]  Storage of results
%    M_:          [structure]  Definition of the model; makes sure that
%                               M_.params is correctly updated
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2009-2026 Dynare Team
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

% indices of endogenous variables

if isfield(oo_,'shock_decomposition_info') && isfield(oo_.shock_decomposition_info,'i_var')
    if isfield (oo_,'realtime_conditional_shock_decomposition') ...
            || isfield (oo_,'realtime_forecast_shock_decomposition') ...
            || isfield (oo_,'realtime_shock_decomposition') ...
            || isfield (oo_,'conditional_shock_decomposition') ...
            || isfield (oo_,'initval_decomposition')
        error('shock_decomposition::squeezed shock decompositions are already stored in oo_')
    end
end

if isempty(varlist)
    varlist = M_.endo_names(1:M_.orig_endo_nbr);
end

[i_var,~,index_uniques] = varlist_indices(varlist,M_.endo_names);
varlist=varlist(index_uniques);

% number of shocks
nshocks = M_.exo_nbr;

% parameter set
parameter_set = options_.parameter_set;
if isempty(parameter_set)
    if isfield(oo_,'posterior_mean')
        parameter_set = 'posterior_mean';
    elseif isfield(oo_,'mle_mode')
        parameter_set = 'mle_mode';
    elseif isfield(oo_,'posterior')
        parameter_set = 'posterior_mode';
    else
        error(['shock_decomposition: option parameter_set is not specified ' ...
               'and posterior mode is not available'])
    end
    options_.parameter_set=parameter_set; %store local copy to make sure subsequently called routines use same value
end

if ~isempty(options_.shock_decomp.forecast_type)
    switch options_.shock_decomp.forecast_type
        case 'unconditional'
            if strcmp(parameter_set,'calibration')
                if ~isfield(oo_, 'forecast')
                    error('Can''t find unconditional forecasts in oo_.forecasts. Did you run the forecast command?');
                else
                    uf = oo_.forecast;
                end
            elseif strcmp(parameter_set,'mle_mode')
                if ~isfield(oo_, 'forecast')
                    error('Can''t find unconditional forecasts in oo_.forecasts. Did you run ML estimation with the forecast option?');
                else
                    uf = oo_.forecast;
                end
            else
                if ~isfield(oo_, 'MeanForecast')
                    error('Can''t find unconditional forecasts in oo_.MeanForecast. Did you run Bayesian estimation with the forecast option?');
                else
                    uf = oo_.MeanForecast;
                end
            end
        case 'conditional'
            if ~isfield(oo_, 'conditional_forecast')
                error('Can''t find conditional forecasts');
            end
            cf = oo_.conditional_forecast;
    end
end

options_.selected_variables_only = 0; %make sure all variables are stored
options_.nograph=true;
options_.plot_priors=0;
[oo_temp, ~, ~, ~, Smoothed_Variables_deviation_from_mean, initial_date] = evaluate_smoother(parameter_set, varlist, M_, oo_, options_, bayestopt_, estim_params_);

options_.initial_date=initial_date;
% reduced form
dr = oo_temp.dr;

% initialization
gend = size(oo_temp.SmoothedShocks.(M_.exo_names{1}),1);

if ~isempty(options_.shock_decomp.forecast_type)
    smoothed_data_present = false;
    if ~isempty(options_.dataset) || ~isempty(options_.datafile)
        smoothed_data_present = true;
    else
        gend = 0;
    end

    if strcmp(options_.shock_decomp.forecast_type, 'unconditional')
        options_.plot_shock_decomp.forecast_length = length(uf.Mean.(M_.endo_names{i_var(1)}));
        gend = gend + options_.plot_shock_decomp.forecast_length;
    else
        options_.plot_shock_decomp.forecast_length = length(cf.cond.Mean.(M_.endo_names{i_var(1)}));
        gend = gend + options_.plot_shock_decomp.forecast_length;
    end

    if smoothed_data_present && strcmp(options_.shock_decomp.forecast_type, 'conditional')
        % initial condition is in the conditional forecast data
        gend = gend - 1;
    end

    % initialization
    epsilon=zeros(nshocks,gend);
    z = zeros(M_.endo_nbr,nshocks+2,gend);

    if smoothed_data_present
        smoothed_periods = size(Smoothed_Variables_deviation_from_mean,2);
        for i=1:nshocks
            epsilon(i,1:smoothed_periods) = oo_temp.SmoothedShocks.(M_.exo_names{i});
        end
        z(:,end,1:smoothed_periods) = Smoothed_Variables_deviation_from_mean;
    end

    if strcmp(options_.shock_decomp.forecast_type, 'unconditional')
        for i=1:size(i_var, 1)
            z(i_var(i),end,smoothed_periods+1:end) = ...
                uf.Mean.(M_.endo_names{i_var(i)})(:) - dr.ys(i_var(i));
        end
    else
        for i=1:size(i_var, 1)
            z(i_var(i),end,smoothed_periods+1:end) = ...
                cf.cond.Mean.(M_.endo_names{i_var(i)})(2:end) - dr.ys(i_var(i));
        end

        conditional_periods = length(cf.controlled_exo_variables.Mean.(M_.exo_names{1}));
        for i=1:nshocks
            epsilon(i,smoothed_periods+(1:conditional_periods)) = ...
                cf.controlled_exo_variables.Mean.(M_.exo_names{i});
        end
    end
else
    epsilon=NaN(nshocks,gend);
    z = zeros(M_.endo_nbr,nshocks+2,gend);

    for i=1:nshocks
        epsilon(i,:) = oo_temp.SmoothedShocks.(M_.exo_names{i});
    end
    z(:,end,:) = Smoothed_Variables_deviation_from_mean;
end


i_state = dr.order_var(M_.nstatic+(1:M_.nspred));
for i=1:gend
    if i > 1 && i <= M_.maximum_lag+1
        lags = min(i-1,M_.maximum_lag):-1:1;
    end

    if i > 1
        tempx = permute(z(:,1:nshocks,lags),[1 3 2]);
        m = min(i-1,M_.maximum_lag);
        tempx = [reshape(tempx,M_.endo_nbr*m,nshocks); zeros(M_.endo_nbr*(M_.maximum_lag-i+1),nshocks)];
        z(:,1:nshocks,i) = dr.ghx(dr.inv_order_var,:)*tempx(i_state,:);
        lags = lags+1;
    end

    if i > options_.shock_decomp.init_state
        z(:,1:nshocks,i) = z(:,1:nshocks,i) + dr.ghu(dr.inv_order_var,:).*repmat(epsilon(:,i)',M_.endo_nbr,1);
    end
    z(:,nshocks+1,i) = z(:,nshocks+2,i) - sum(z(:,1:nshocks,i),2);
end

if options_.shock_decomp.with_epilogue
    [z, oo_.shock_decomposition_info.epilogue_steady_state] = epilogue_shock_decomposition(z, M_, oo_temp);
end

if ~isempty(options_.shock_decomp.forecast_type)
    oo_.forecast_shock_decomposition.(options_.shock_decomp.forecast_type) = z;
    oo_.shock_decomposition_info.forecast.(options_.shock_decomp.forecast_type).parameter_set=options_.parameter_set;
else
    oo_.shock_decomposition = z;
    oo_.shock_decomposition_info.parameter_set=options_.parameter_set;
end

if ~options_.no_graph.shock_decomposition
    oo_ = plot_shock_decomposition(M_,oo_,options_,varlist);
end

oo_.gui.ran_shock_decomposition = true;
