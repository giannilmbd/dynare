function [oo_, options_] = simulate(M_, options_, oo_, var_list)
% Compute IRFs, deterministic simulations, or stochastic simulations for models
% featuring rich heterogeneity
%
% SYNTAX:
%   [oo_, options_] = heterogeneity.simulate(M_, options_, oo_, var_list)
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure with optional fields:
%                            - irf: IRF horizon (default: 40)
%                            - periods: simulation horizon (default: 0 = IRF mode)
%                            - irf_shocks: cell array of shock names for IRFs (default: all shocks)
%                            - relative_irf: plot IRFs as % deviations from steady state (default: false)
%                            - nograph, nodisplay, graph_format: graph display options
%   oo_          [struct]  Dynare results structure (must contain oo_.heterogeneity.dr.G)
%   var_list     [cell array] List of variable names to display/compute (empty cell array {} for all variables)
%
% SIMULATION MODES (automatically detected):
%   1. Stochastic simulation (no M_.det_shocks):
%      Computes responses to unanticipated shocks drawn from distributions in shocks block.
%      - IRFs: Always computed (responses to one-time shocks)
%      - Simulation paths: Computed when periods > 0 (random shocks each period)
%      Both can be requested together in a single command.
%   2. News shock sequence (M_.det_shocks exists):
%      Triggered when mod file contains shocks block with 'periods' and 'values' keywords.
%      Simulates response to anticipated shocks (news shocks) known at date 0.
%      Following Auclert et al. (2021): households learn at t=0 about shocks at various future dates.
%      Mutually exclusive with stochastic simulation options.
%
% OUTPUTS:
%   oo_          [struct] Updated results structure with:
%   options_     [struct] Updated options structure with graph-related settings
%
%                  Stochastic simulation mode:
%                    oo_.irfs.(var)_(shock): [T × 1] IRF for each variable-shock pair (always computed)
%                    oo_.endo_simul: [n_vars × (periods+1)] simulated paths in levels (only when periods > 0)
%
%                  News shock sequence mode:
%                    oo_.endo_simul: [n_vars × (periods+1)] simulated paths (levels)
%                    oo_.exo_simul: [total_periods × n_exo] news shock paths (populated by make_ex_)
%
% EXAMPLES:
%   % Compute IRFs for all shocks and all variables
%   oo_ = heterogeneity.simulate(M_, options_, oo_, {});
%
%   % Compute IRFs only for TFP shock
%   options_.irf_shocks = {'Z'};
%   oo_ = heterogeneity.simulate(M_, options_, oo_, {});
%
%   % Stochastic simulation over 1000 periods
%   options_.periods = 1000;
%   oo_ = heterogeneity.simulate(M_, options_, oo_, {});
%
%   % Compute BOTH IRFs and stochastic simulation
%   options_.irf = 40;
%   options_.periods = 1000;
%   oo_ = heterogeneity.simulate(M_, options_, oo_, {});
%
%   % News shock sequence via shocks block in mod file:
%   % In mod file (anticipated shocks known at t=0):
%   %   shocks;
%   %   var Z; periods 1:10; values 0.01;
%   %   var Z; periods 20; values -0.005;
%   %   end;
%   % Then call:
%   options_.periods = 40;
%   oo_ = heterogeneity.simulate(M_, options_, oo_, {});
%
%   % Compute only consumption and output IRFs
%   oo_ = heterogeneity.simulate(M_, options_, oo_, {'C', 'Y'});
%
%   % Access stored IRFs
%   oo_ = heterogeneity.simulate(M_, options_, oo_, {});
%   plot(oo_.irfs.Y_Z);

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

%% Step 1: Input Validation and Parsing

% Default to all endogenous variables if var_list is empty
if isempty(var_list)
    var_list = M_.endo_names;
end

% Check that solution exists
if ~isfield(oo_, 'heterogeneity') || ~isfield(oo_.heterogeneity, 'dr') || ~isfield(oo_.heterogeneity.dr, 'G')
    error('No solution found. Run heterogeneity_solve first.');
end

% Determine shock list for IRFs
if isfield(options_, 'irf_shocks') && ~isempty(options_.irf_shocks)
    shock_list = options_.irf_shocks;
    % Validate shock names
    for i = 1:length(shock_list)
        if ~any(strcmp(shock_list{i}, M_.exo_names))
            error('heterogeneity:simulate:UnknownShock', ...
                  'Shock "%s" not found in M_.exo_names', shock_list{i});
        end
    end
else
    % Default: all aggregate shocks
    shock_list = M_.exo_names;
end

% Validate option compatibility with simulation mode
% News shock mode (det_shocks) is incompatible with stochastic IRF options
if isfield(M_, 'det_shocks') && ~isempty(M_.det_shocks)
    % Check if any stochastic IRF options have non-default values
    stochastic_options_set = {};

    if options_.irf ~= 40
        stochastic_options_set{end+1} = sprintf('irf=%d (default: 40)', options_.irf);
    end
    if options_.periods ~= 0
        stochastic_options_set{end+1} = sprintf('periods=%d (default: 0)', options_.periods);
    end
    if isfield(options_, 'irf_shocks') && ~isempty(options_.irf_shocks)
        stochastic_options_set{end+1} = 'irf_shocks';
    end
    if options_.relative_irf ~= false
        stochastic_options_set{end+1} = 'relative_irf=true (default: false)';
    end
    if options_.impulse_responses.plot_threshold ~= 1e-10
        stochastic_options_set{end+1} = sprintf('irf_plot_threshold=%g (default: 1e-10)', ...
            options_.impulse_responses.plot_threshold);
    end

    if ~isempty(stochastic_options_set)
        error('heterogeneity:simulate:IncompatibleOptions', ...
              ['News shock simulation (shocks block with periods/values) is incompatible with ' ...
               'stochastic IRF options.\nThe following non-default options were detected:\n  - %s\n' ...
               'Remove these options or use a shocks block without periods/values for stochastic simulation.'], ...
              strjoin(stochastic_options_set, '\n  - '));
    end
end

%% Step 3: Determine simulation mode

% Two mutually exclusive modes:
%
% 1. STOCHASTIC SIMULATION (no shocks block):
%    - Computes IRFs to unanticipated shocks
%    - Can also compute stochastic simulation paths (periods > 0)
%    - Uses linear G matrices to compute responses
%    - Follows stoch_simul.m patterns
%
% 2. NEWS SHOCK SEQUENCE (shocks block with periods/values):
%    - Simulates response to anticipated shocks known at date 0
%    - Uses make_ex_ to parse shocks block (convenient syntax)
%    - Does NOT use perfect foresight solver (histval/endval/boundary conditions)
%    - Treats as linear news shocks from steady state
%    - Uses G matrices (same as stochastic mode, different shock sequences)

if isfield(M_, 'det_shocks') && ~isempty(M_.det_shocks)
    % NEWS SHOCK SEQUENCE MODE
    oo_ = heterogeneity.simulate_news_shocks(M_, options_, oo_, var_list, shock_list);
else
    % STOCHASTIC SIMULATION MODE
    % Validate periods option
    if options_.periods < 0
        error('heterogeneity:simulate:NegativePeriods', ...
              'Option periods cannot be negative.');
    end

    % Always compute IRFs
    oo_.irfs = heterogeneity.simulate_irfs(M_, options_, oo_, var_list, shock_list);

    % Compute stochastic simulation paths if requested
    if options_.periods > 0
        oo_ = heterogeneity.simulate_stochastic_shocks(M_, options_, oo_, var_list, shock_list);
    end
end

end  % main function