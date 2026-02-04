function oo_ = simulate_stochastic_shocks(M_, options_, oo_, var_list, shock_list)
% Simulate stochastic paths for heterogeneous-agent models
%
% SYNTAX:
%   oo_ = heterogeneity.simulate_stochastic_shocks(M_, options_, oo_, var_list, shock_list)
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure (uses options_.periods)
%   oo_          [struct]  Dynare results structure (must contain oo_.heterogeneity.dr.G)
%   var_list     [cell]    List of variable names to plot
%   shock_list   [cell]    List of shock names
%
% OUTPUTS:
%   oo_          [struct]  Updated results with oo_.endo_simul [n_endo × periods]
%
% NOTE:
%   This is a lightweight function with minimal validation. For full input
%   validation and automatic mode detection, use heterogeneity.simulate() instead.

% Copyright © 2025 Dynare Team
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
%
% Original author: Normann Rion <normann@dynare.org>

% Check simulation horizon vs drop (following stoch_simul.m pattern)
% Commented out as it is not useful yet.
% if options_.periods <= options_.drop
%     error('heterogeneity:simulate_stochastic_shocks:PeriodsTooShort', ...
%           'The horizon of simulation (periods=%d) is shorter than the number of observations to be dropped (drop=%d).', ...
%           options_.periods, options_.drop);
% end

% Get shock covariance matrix from M_.Sigma_e
if ~isfield(M_, 'Sigma_e') || isempty(M_.Sigma_e)
    error('heterogeneity:simulate_stochastic_shocks:NoShockVariance', ...
          'M_.Sigma_e is empty. Declare shock variances in shocks block.');
end

% Get Jacobian horizon
truncation_horizon = options_.heterogeneity.solve.truncation_horizon;

% Eliminate shocks with 0 variance (following simult.m pattern)
i_exo_var = setdiff(1:M_.exo_nbr, find(diag(M_.Sigma_e) == 0));
nxs = length(i_exo_var);

% Draw random shocks: eps ~ N(0, Sigma_e)
% Following SSJ approach: draw (truncation_horizon - 1) pre-sample shocks
% plus options_.periods simulation shocks to avoid startup effects
set_dynare_seed_local_options([],false,'default');
total_shocks = options_.periods + truncation_horizon - 1;
eps_draws = zeros(M_.exo_nbr, total_shocks);
if nxs > 0
    chol_S = chol(M_.Sigma_e(i_exo_var, i_exo_var));
    eps_draws(i_exo_var, :) = (randn(nxs, total_shocks)' * chol_S)';
end

% Initialize oo_.endo_simul: [n_endo × options_.periods]
% Following stoch_simul.m: no initial steady state column
oo_.endo_simul = zeros(M_.endo_nbr, options_.periods);

% For each variable, compute path as linear combination of shock responses
oo_het = oo_.heterogeneity;
for i_var = 1:M_.endo_nbr
    var_name = M_.endo_names{i_var};
    path = zeros(options_.periods, 1);

    % Sum over all shocks
    for i_shock = 1:length(shock_list)
        shock_name = shock_list{i_shock};
        shock_idx = find(strcmp(shock_name, M_.exo_names));

        % Error if G matrix doesn't exist
        if ~isfield(oo_.heterogeneity.dr.G, var_name) || ...
           ~isfield(oo_.heterogeneity.dr.G.(var_name), shock_name)
            error('heterogeneity:simulate_stochastic_shocks:NoJacobian', ...
                  'No Jacobian found for variable "%s" and shock "%s". Run heterogeneity_solve first.', ...
                  var_name, shock_name);
        end

        % Get unit IRF (response to unit shock at period 1)
        G = oo_.heterogeneity.dr.G.(var_name).(shock_name);
        unit_irf = G(1:truncation_horizon, 1);

        % Convolve unit IRF with shock sequence
        % Formula: Y_t = sum_{s=1}^T unit_irf(s) * eps_draws(t + T - s)
        for t = 1:options_.periods
            for s = 1:truncation_horizon
                path(t) = path(t) + unit_irf(s) * eps_draws(shock_idx, t + truncation_horizon - s);
            end
        end
    end

    % oo_.endo_simul contains levels (not deviations), consistent with stoch_simul
    oo_.endo_simul(i_var, :) = oo_het.mat.y(i_var) + path';
end

% Plot if requested
if ~options_.nograph
    config = struct();
    config.figure_title = 'Heterogeneity: Stochastic Simulation';
    config.filename_suffix = '_simulation';
    config.tex_caption = 'Stochastic simulation.';
    config.tex_label = 'simulation';
    heterogeneity.plot_simulation(M_, options_, oo_.endo_simul, var_list, config);
end

end
