function [oo_, ts] = perfect_foresight_with_expectation_errors_solver(M_, options_, oo_)
% INPUTS
%   M_                  [structure] describing the model
%   options_            [structure] describing the options
%   oo_                 [structure] storing the results
%
% OUTPUTS
%   oo_                 [structure] storing the results

% Copyright © 2021-2025 Dynare Team
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

if options_.pfwee.constant_simulation_length && ~isempty(options_.simul.last_simulation_period)
    error('Options constant_simulation_length and last_simulation_period cannot be used together')
end

if options_.pfwee.constant_simulation_length && ~all(cellfun(@(x) isempty(x), {oo_.pfwee.controlled_paths_by_period.endogenize_id}))
    error('Options constant_simulation_length is incompatible with perfect_foresight_controlled_paths block')
end

[periods, first_simulation_period] = get_simulation_periods(options_);

% Retrieve initial paths built by pfwee_setup
% (the versions in oo_ will be truncated before calling perfect_foresight_solver)
endo_simul = oo_.endo_simul;
exo_simul = oo_.exo_simul;

% Enforce the endval steady option in pf_solver
if M_.maximum_lead > 0
    options_.simul.endval_steady = true;
end

% When there is no regular endval block (or endval(learnt_in=1) block), ensure
% that the initial steady state is stored, so that the correct “base” scenario is
% constructed in the homotopy loop (see perfect_foresight_solver.m)
if isempty(oo_.initial_steady_state)
    oo_.initial_steady_state = oo_.steady_state;
    oo_.initial_exo_steady_state = oo_.exo_steady_state;
end

% Start main loop around informational periods
info_period = 1;
increment = 0;
while info_period <= periods
    if ~options_.noprint
        fprintf('perfect_foresight_with_expectations_errors_solver: computing solution for information available at period %d\n', info_period)
    end

    % Compute terminal steady state as anticipated
    oo_.exo_steady_state = oo_.pfwee.terminal_info(:, info_period);
    % oo_.steady_state will be updated by pf_solver (since endval_steady=true)

    if options_.pfwee.constant_simulation_length && increment > 0
        % Use previous terminal steady state as guess value for: simulation periods that don’t yet have an initial guess (i.e. are NaNs at this point); and also for the terminal steady state
        endo_simul = [ endo_simul repmat(oo_.steady_state, 1, increment)];
        exo_simul = [ exo_simul; NaN(increment, M_.exo_nbr)];
    end

    oo_.endo_simul = endo_simul(:, info_period:end); % Take initial conditions + guess values from previous simulation
    if options_.pfwee.constant_simulation_length
        sim_length = periods;
    else
        sim_length = periods - info_period + 1;
    end

    oo_.exo_simul = exo_simul(info_period:end, :);
    oo_.exo_simul(M_.maximum_lag+(1:periods-info_period+1), :) = oo_.pfwee.shocks_info(:, info_period:end, info_period)';
    oo_.exo_simul(M_.maximum_lag+periods-info_period+2:end, :) = repmat(oo_.exo_steady_state', sim_length+M_.maximum_lead-(periods-info_period+1), 1);

    oo_.deterministic_simulation.controlled_paths_by_period = oo_.pfwee.controlled_paths_by_period(info_period:end, info_period);
    if all(cellfun(@(x) isempty(x), {oo_.deterministic_simulation.controlled_paths_by_period.endogenize_id}))
        oo_.deterministic_simulation.controlled_paths_by_period = []; % Expected by pf_solver when there is no controlled var
    end

    options_.periods = sim_length;
    % The following two options are reset to empty, so as to avoid an inconsistency with periods
    options_.simul.first_simulation_period = dates();
    options_.simul.last_simulation_period = dates();

    if info_period > 1 && homotopy_completion_share < 1 && options_.simul.homotopy_marginal_linearization_fallback > 0
        marginal_linearization_previous_raw_sims.sim1.endo_simul = oo_.deterministic_simulation.sim1.endo_simul(:, increment+1:end);
        marginal_linearization_previous_raw_sims.sim1.exo_simul = oo_.deterministic_simulation.sim1.exo_simul(increment+1:end, :);
        marginal_linearization_previous_raw_sims.sim1.homotopy_completion_share = oo_.deterministic_simulation.sim1.homotopy_completion_share;
        marginal_linearization_previous_raw_sims.sim2.endo_simul = oo_.deterministic_simulation.sim2.endo_simul(:, increment+1:end);
        marginal_linearization_previous_raw_sims.sim2.exo_simul = oo_.deterministic_simulation.sim2.exo_simul(increment+1:end, :);
        marginal_linearization_previous_raw_sims.sim2.homotopy_completion_share = oo_.deterministic_simulation.sim2.homotopy_completion_share;
    else
        marginal_linearization_previous_raw_sims = [];
    end

    oo_= perfect_foresight_solver(M_, options_, oo_, marginal_linearization_previous_raw_sims);

    if ~oo_.deterministic_simulation.status
        error('perfect_foresight_with_expectation_errors_solver: failed to compute solution for information available at period %d\n', info_period)
    end

    if info_period == 1
        homotopy_completion_share = oo_.deterministic_simulation.homotopy_completion_share;
        options_.simul.homotopy_max_completion_share = homotopy_completion_share;
    elseif oo_.deterministic_simulation.homotopy_completion_share ~= homotopy_completion_share
        %% NB: in the following message, we don’t use the %.3f formatter because it may round the
        %% share to a greater number, which would lead to an incorrect suggestion
        error('perfect_foresight_solver_with_expectation_errors: could not find a solution for information available at period %d with the same homotopy completion share as period 1; a possible solution is to retry the simulation with homotopy_max_completion_share=%s\n', ...
              info_period, num2str(floor(oo_.deterministic_simulation.homotopy_completion_share*1000)/1000))
    end

    endo_simul(:, info_period:end) = oo_.endo_simul;
    exo_simul(info_period:end, :) = oo_.exo_simul;

    % Increment info_period (as much as possible, if information set does not change for some time)
    increment = 1;
    while info_period+increment <= periods && ...
          all(oo_.pfwee.terminal_info(:, info_period) == oo_.pfwee.terminal_info(:, info_period+increment)) && ...
          all(all(oo_.pfwee.shocks_info(:, info_period+increment:end, info_period) == oo_.pfwee.shocks_info(:, info_period+increment:end, info_period+increment))) && ...
          isequal(oo_.pfwee.controlled_paths_by_period(info_period+increment:end, info_period), oo_.pfwee.controlled_paths_by_period(info_period+increment:end, info_period+increment))
        increment = increment + 1;
    end
    info_period = info_period + increment;
end

% Set final paths
oo_.endo_simul = endo_simul;
oo_.exo_simul = exo_simul;

if nargout > 1
    ts = construct_simulation_dseries(oo_, M_, first_simulation_period);
end
