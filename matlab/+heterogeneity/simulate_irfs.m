function oo_irfs = simulate_irfs(M_, options_, oo_, var_list, shock_list)
% Compute impulse response functions for heterogeneous-agent models
%
% SYNTAX:
%   oo_irfs = heterogeneity.simulate_irfs(M_, options_, oo_, var_list, shock_list)
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure
%   oo_          [struct]  Dynare results structure (must contain oo_.heterogeneity.dr.G)
%   var_list     [cell]    List of variable names to compute IRFs for
%   shock_list   [cell]    List of shock names to compute IRFs for
%
% OUTPUTS:
%   oo_irfs          [struct]  Results with oo_irfs.(var)_(shock) for each pair
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

% Initialization
oo_irfs = struct();
irf_horizon = options_.irf;

% Adjust if requested T differs from Jacobian T
truncation_horizon = options_.heterogeneity.solve.truncation_horizon;
if irf_horizon > truncation_horizon
    error('heterogeneity:simulate_irfs:HorizonTooLong', ...
          'Requested irf=%d exceeds Jacobian truncation horizon %d. Increase truncation_horizon in heterogeneity_solve or reduce irf.', ...
          irf_horizon, truncation_horizon);
end

% Get Cholesky decomposition of shock covariance matrix
if isfield(M_, 'Sigma_e') && ~isempty(M_.Sigma_e)
    cs = get_lower_cholesky_covariance(M_.Sigma_e, options_.add_tiny_number_to_cholesky);
else
    cs = [];
end

% For each shock, compute IRFs for all variables
for i_shock = 1:length(shock_list)
    shock_name = shock_list{i_shock};

    % Get shock standard deviation from Cholesky decomposition
    shock_idx = find(strcmp(shock_name, M_.exo_names));
    if ~isempty(cs) && size(cs, 1) >= shock_idx
        shock_std = cs(shock_idx, shock_idx);
    else
        shock_std = 0;
    end

    % Shock path: one standard deviation shock at period 1 (or unit for relative_irf)
    shock_path = zeros(truncation_horizon, 1);
    if options_.relative_irf
        % For relative_irf: use unit shock (100 %).
        shock_path(1) = 100;
    else
        % Default: one standard deviation shock
        shock_path(1) = shock_std;
    end

    % Compute IRF for each variable
    for i_var = 1:length(var_list)
        var_name = var_list{i_var};

        % Check if variable exists
        if ~any(strcmp(var_name, M_.endo_names))
            warning('heterogeneity:simulate_irfs:UnknownVariable', ...
                    'Variable "%s" not found, skipping', var_name);
            continue;
        end

        % Check if G matrix exists for this variable-shock pair
        if ~isfield(oo_.heterogeneity.dr.G, var_name) || ...
           ~isfield(oo_.heterogeneity.dr.G.(var_name), shock_name)
            warning('heterogeneity:simulate_irfs:NoJacobian', ...
                    'No Jacobian found for variable "%s" and shock "%s", skipping', var_name, shock_name);
            continue;
        end

        % Get Jacobian G^{var, shock}
        G_block = oo_.heterogeneity.dr.G.(var_name).(shock_name);  % [truncation_horizon × truncation_horizon]

        % Compute IRF: y(t) = sum_{s=1}^t G(t,s) * epsilon(s)
        irf_var = G_block * shock_path;  % [truncation_horizon × 1]

        % Adjust to requested horizon irf_horizon
        irf_var = irf_var(1:irf_horizon);

        % Store in oo_.irfs with Dynare naming convention: varname_shockname
        field_name = [var_name '_' shock_name];
        oo_irfs.(field_name) = irf_var;
    end
end

% Plot if requested
if ~options_.nograph
    heterogeneity.plot_irfs(M_, options_, oo_irfs, shock_list, var_list);
end

end
