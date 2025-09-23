function oo_ = simulate_news_shocks(M_, options_, oo_, var_list, shock_list)
% Simulate responses to anticipated (news) shocks for heterogeneous-agent models
%
% SYNTAX:
%   oo_ = heterogeneity.simulate_news_shocks(M_, options_, oo_, var_list, shock_list)
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure
%   oo_          [struct]  Dynare results structure (must contain oo_.heterogeneity.dr.G and oo_.exo_simul)
%   var_list     [cell]    List of variable names to plot
%   shock_list   [cell]    List of shock names
%
% OUTPUTS:
%   oo_          [struct]  Updated results with:
%                            oo_.endo_simul [n_endo × periods] simulated paths
%                            oo_.exo_simul  [total_periods × n_exo] shock paths (from make_ex_)
%
% DESCRIPTION:
%   This function simulates responses to a sequence of ANTICIPATED shocks (news shocks)
%   specified via the shocks block with periods and values syntax.
%
%   APPROACH:
%   - Uses make_ex_ to parse the shocks block (convenient periods/values syntax)
%   - Does NOT use perfect foresight machinery (nonlinear solver, boundary conditions)
%   - Treats shocks as news known at date 0, computes linear responses via G matrices
%   - Economy starts at steady state (no histval/endval/complex initial conditions)
%
%   UNSUPPORTED PERFECT FORESIGHT FEATURES:
%   - histval (initial conditions ≠ steady state)
%   - endval (terminal conditions ≠ steady state)
%   - initval_series (time-varying initial/terminal paths)
%   These features are part of perfect_foresight_setup but not relevant for news shocks.
%
% NOTE:
%   This is a lightweight function with minimal validation. For full input
%   validation and automatic mode detection, use heterogeneity.simulate() instead.
%   The caller must have called make_ex_(M_, options_, oo_) before calling this function
%   to populate oo_.exo_simul.

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

% Call make_ex_ to parse shocks block and populate oo_.exo_simul
oo_ = make_ex_(M_, options_, oo_);

% Validate that perfect foresight features are not used
if ~isempty(M_.endo_histval)
    error('heterogeneity:simulate_news_shocks:HistvalNotSupported', ...
          'histval not supported for news shock simulation. Simulation starts from steady state.');
end
if ~isempty(oo_.initial_steady_state) || ~isempty(oo_.initial_exo_steady_state)
    error('heterogeneity:simulate_news_shocks:EndvalNotSupported', ...
          'endval not supported for news shock simulation. Simulation starts from steady state.');
end
if ~isempty(oo_.initval_series)
    error('heterogeneity:simulate_news_shocks:InitvalSeriesNotSupported', ...
          'initval_series not supported for news shock simulation. Simulation starts from steady state.');
end

% Validate that shock sequence was populated by make_ex_
if ~isfield(oo_, 'exo_simul') || isempty(oo_.exo_simul)
    error('heterogeneity:simulate_news_shocks:NoExoSimul', ...
          'No shock sequence found. Specify shocks using a shocks block with periods and values.');
end

% Infer number of simulation periods from shock sequence created by make_ex_
% For news shocks, periods is not set in options_.periods but from shock sequence
periods = size(oo_.exo_simul, 1) - M_.maximum_lag - M_.maximum_lead;
if periods <= 0
    error('heterogeneity:simulate_news_shocks:InvalidPeriods', ...
          'Invalid simulation horizon: periods must be positive.');
end

% Validate against Jacobian horizon
truncation_horizon = options_.heterogeneity.solve.truncation_horizon;
if periods > truncation_horizon
    error('heterogeneity:simulate_news_shocks:SequenceTooLong', ...
          'Simulation horizon periods=%d exceeds Jacobian truncation horizon %d. Reduce periods or increase truncation_horizon in heterogeneity_solve.', ...
          periods, truncation_horizon);
end

% Extract shock values for simulation periods (ignore lags/leads from make_ex_)
% make_ex_ creates: [max_lag + periods + max_lead] × n_exo
% We only need the periods simulation periods: rows [max_lag+1 : max_lag+periods]
start_idx = M_.maximum_lag + 1;
end_idx = M_.maximum_lag + periods;
eps_seq = oo_.exo_simul(start_idx:end_idx, :)';  % [n_exo × periods]

% Initialize oo_.endo_simul: [n_endo × periods], starting from steady state
% Following stoch_simul pattern: periods 1 to periods, no initial column
oo_.endo_simul = zeros(M_.endo_nbr, periods);

% Compute the aggregate steady state
oo_het = oo_.heterogeneity;
steady_state = heterogeneity.internal.compute_agg_steady_state(M_, oo_het.ss.agg, oo_het.mat.pol.x_bar_dash, oo_het.mat.d.Phi, oo_het.mat.d.hist, oo_het.indices.Ix.in_x);

% Compute linear responses using impulse response Jacobians from heterogeneity_solve
for i_var = 1:M_.endo_nbr
    var_name = M_.endo_names{i_var};
    deviation = zeros(periods, 1);  % Deviation from steady state

    % Sum contributions from all shocks
    for i_shock = 1:length(shock_list)
        shock_name = shock_list{i_shock};
        shock_idx = find(strcmp(shock_name, M_.exo_names));

        % Error if no Jacobian available for this variable-shock pair
        if ~isfield(oo_.heterogeneity.dr.G, var_name) || ...
           ~isfield(oo_.heterogeneity.dr.G.(var_name), shock_name)
            error('heterogeneity:simulate_news_shocks:NoJacobian', ...
                  'No Jacobian found for variable "%s" and shock "%s". Run heterogeneity_solve first.', ...
                  var_name, shock_name);
        end

        % Get impulse response Jacobian: G(t,s) = response at t to unit shock at s
        G = oo_.heterogeneity.dr.G.(var_name).(shock_name);  % [truncation_horizon × truncation_horizon]

        % Convolve with shock sequence: response(t) = sum_s G(t,s) * shock(s)
        deviation = deviation + G(1:periods,1:periods)*eps_seq(shock_idx, :)';
    end

    % Store as levels (steady state + deviation), consistent with stoch_simul
    oo_.endo_simul(i_var, :) = steady_state(i_var) + deviation';
end

% Note: oo_.exo_simul already populated by make_ex_, no modification needed

% Plot if requested
if ~options_.nograph
    config = struct();
    config.figure_title = 'Heterogeneity: Deterministic Shock Sequence';
    config.filename_suffix = '_det_simulation';
    config.tex_caption = 'Deterministic shock sequence.';
    config.tex_label = 'det_simulation';
    heterogeneity.plot_simulation(M_, options_, oo_.endo_simul, var_list, config);
end

end
