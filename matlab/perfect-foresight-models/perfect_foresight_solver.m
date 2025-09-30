function [oo_, ts]=perfect_foresight_solver(M_, options_, oo_, marginal_linearization_previous_raw_sims)
% Computes deterministic simulations
%
% INPUTS
%   M_                  [structure] describing the model
%   options_            [structure] describing the options
%   oo_                 [structure] storing the results
%   marginal_linearization_previous_raw_sims [struct, optional]
%       if not empty, contains the two simulations used to compute the extrapolation by marginal
%       linearization in a previous informational period, in the context of
%       perfect_foresight_with_expectation_errors in combination with homotopy and marginal
%       linearization
%
% OUTPUTS
%   oo_                 [structure] storing the results
%   ts                  [dseries]   final simulation paths
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 1996-2025 Dynare Team
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

check_input_arguments(options_, M_, oo_);

if nargin < 4
    marginal_linearization_previous_raw_sims = [];
end

[periods, first_simulation_period, last_simulation_period] = get_simulation_periods(options_);

if options_.debug
    model_static = str2func([M_.fname,'.sparse.static_resid']);
    for ii=1:size(oo_.exo_simul,1)
        [residual(:,ii)] = model_static(oo_.steady_state, oo_.exo_simul(ii,:),M_.params);
    end
    problematic_periods=find(any(isinf(residual)) | any(isnan(residual)))-M_.maximum_endo_lag;
    if ~isempty(problematic_periods)
        period_string=num2str(problematic_periods(1));
        for ii=2:length(problematic_periods)
            period_string=[period_string, ', ', num2str(problematic_periods(ii))];
        end
        fprintf('\n\nWARNING: Value for the exogenous variable(s) in period(s) %s inconsistent with the static model.\n',period_string);
        fprintf('WARNING: Check for division by 0.\n')
    end
end

% Various sanity checks

if options_.no_homotopy && (options_.simul.homotopy_initial_step_size ~= 1 ...
                            || options_.simul.homotopy_max_completion_share ~= 1 ...
                            || options_.simul.homotopy_linearization_fallback ...
                            || options_.simul.homotopy_marginal_linearization_fallback ~= 0)
    error('perfect_foresight_solver: the no_homotopy option is incompatible with homotopy_initial_step_size, homotopy_max_completion_share, homotopy_linearization_fallback and homotopy_marginal_linearization_fallback options')
end

if options_.simul.homotopy_initial_step_size > 1 || options_.simul.homotopy_initial_step_size < 0
    error('perfect_foresight_solver: The value given to homotopy_initial_step_size option must be in the [0,1] interval')
end

if options_.simul.homotopy_min_step_size > 1 || options_.simul.homotopy_min_step_size < 0
    error('perfect_foresight_solver: The value given to homotopy_min_step_size option must be in the [0,1] interval')
end

if options_.simul.homotopy_max_completion_share > 1 || options_.simul.homotopy_max_completion_share < 0
    error('perfect_foresight_solver: The value given to homotopy_max_completion_share option must be in the [0,1] interval')
end

if options_.simul.homotopy_marginal_linearization_fallback > 1 || options_.simul.homotopy_marginal_linearization_fallback < 0
    error('perfect_foresight_solver: The value given to homotopy_marginal_linearization_fallback option must be in the [0,1] interval')
end

if options_.simul.homotopy_initial_step_size < options_.simul.homotopy_min_step_size
    error('perfect_foresight_solver: The value given to homotopy_initial_step_size option must be greater or equal to that given to homotopy_min_step_size option')
end

if options_.simul.homotopy_linearization_fallback && options_.simul.homotopy_marginal_linearization_fallback > 0
    error('perfect_foresight_solver: Options homotopy_linearization_fallback and homotopy_marginal_linearization_fallback cannot be used together')
end

if options_.simul.homotopy_max_completion_share < 1 && ~options_.simul.homotopy_linearization_fallback && options_.simul.homotopy_marginal_linearization_fallback == 0
    error('perfect_foresight_solver: Option homotopy_max_completion_share has a value less than 1, so you must also specify either homotopy_linearization_fallback or homotopy_marginal_linearization_fallback')
end

if ~isempty(oo_.deterministic_simulation.controlled_paths_by_period)
    for p = 1:length(oo_.deterministic_simulation.controlled_paths_by_period)
        if isempty(oo_.deterministic_simulation.controlled_paths_by_period(p).exogenize_id)
            continue
        end
        exo_ids = oo_.deterministic_simulation.controlled_paths_by_period(p).endogenize_id;
        if ~isempty(options_.simul.homotopy_exclude_varexo)
            [is_excluded, excluded_exo_ids] = ismember(options_.simul.homotopy_exclude_varexo, M_.exo_names(exo_ids));
            if any(is_excluded)
                error('Exogenous %s cannot be in the homotopy_exclude_varexo option and a perfect_foresight_controlled_paths block at the same time', M_.exo_names{excluded_exo_ids(1)})
            end
        end
    end
end

initperiods = 1:M_.maximum_lag;
simperiods = M_.maximum_lag+(1:periods);
lastperiods = M_.maximum_lag+periods+(1:M_.maximum_lead);

% Create base scenario for homotopy, which corresponds to the initial steady
% state (i.e. a known solution to the perfect foresight problem, assuming that
% oo_.steady_state/oo_.initial_steady_state effectively contains a steady state)
if isempty(oo_.initial_steady_state)
    endobase = repmat(oo_.steady_state, 1,M_.maximum_lag+periods+M_.maximum_lead);
    exobase = repmat(oo_.exo_steady_state',M_.maximum_lag+periods+M_.maximum_lead,1);
else
    endobase = repmat(oo_.initial_steady_state, 1, M_.maximum_lag+periods+M_.maximum_lead);
    exobase = repmat(oo_.initial_exo_steady_state', M_.maximum_lag+periods+M_.maximum_lead, 1);
end

% Determine whether to recompute the final steady state (either because
% option “endval_steady” was passed, or because there is an “endval” block and the
% terminal condition is a steady state)
if options_.simul.endval_steady
    recompute_final_steady_state = true;
elseif ~isempty(oo_.initial_steady_state)
    recompute_final_steady_state = true;
    for j = lastperiods
        endval_resid = evaluate_static_model(oo_.endo_simul(:,j), oo_.exo_simul(j,:)', M_.params, M_, options_);
        if norm(endval_resid, 'Inf') > options_.simul.steady_tolf
            recompute_final_steady_state = false;
            break
        end
    end
else
    recompute_final_steady_state = false;
end

% Perform the homotopy loop
if isempty(marginal_linearization_previous_raw_sims)
    shareorig = 1;
    endoorig = oo_.endo_simul;
    exoorig = oo_.exo_simul;
else
    shareorig = marginal_linearization_previous_raw_sims.sim1.homotopy_completion_share;
    endoorig = marginal_linearization_previous_raw_sims.sim1.endo_simul;
    exoorig = shareorig*oo_.exo_simul + (1-shareorig)*exobase;
    if ~isempty(oo_.deterministic_simulation.controlled_paths_by_period)
        for p = 1:length(oo_.deterministic_simulation.controlled_paths_by_period)
            if isempty(oo_.deterministic_simulation.controlled_paths_by_period(p).exogenize_id)
                continue
            end
            exo_ids = oo_.deterministic_simulation.controlled_paths_by_period(p).endogenize_id;
            exoorig(length(initperiods)+p,exo_ids) = marginal_linearization_previous_raw_sims.sim1.exo_simul(length(initperiods)+p,exo_ids);
        end
    end
end
[completed_share, endo_simul, exo_simul, steady_state, exo_steady_state, iteration, maxerror, solver_iter, per_block_status] = homotopy_loop(M_,options_,oo_,options_.simul.homotopy_max_completion_share, shareorig, endoorig, exoorig, endobase, exobase, initperiods, simperiods, lastperiods, recompute_final_steady_state, oo_.steady_state, oo_.exo_steady_state);

% Do linearization if needed and requested, and put results and solver status information in oo_
if completed_share == 1
    oo_.endo_simul = endo_simul;
    if options_.simul.endval_steady
        oo_.steady_state = steady_state;
    end
    % NB: no need to modify oo_.exo_simul and oo_.exo_steady_state, since we simulated 100% of the shock, except if controlled paths
    if ~isempty(oo_.deterministic_simulation.controlled_paths_by_period)
        oo_.exo_simul = exo_simul;
    end

    if ~options_.noprint
        fprintf('Perfect foresight solution found.\n\n')
    end
    oo_.deterministic_simulation.status = true;
elseif options_.simul.homotopy_linearization_fallback && completed_share > 0
    oo_.deterministic_simulation.sim1.endo_simul = endo_simul;
    oo_.deterministic_simulation.sim1.exo_simul = exo_simul;
    oo_.deterministic_simulation.sim1.steady_state = steady_state;
    oo_.deterministic_simulation.sim1.exo_steady_state = exo_steady_state;
    oo_.deterministic_simulation.sim1.homotopy_completion_share = completed_share;

    oo_.endo_simul = endobase + (endo_simul - endobase)/completed_share;
    if options_.simul.endval_steady
        % The following is needed for the STEADY_STATE() operator to work properly,
        % and thus must come before computing the maximum error.
        % This is not a true steady state, but it is the closest we can get to
        oo_.steady_state = oo_.endo_simul(:, end);
    end
    % NB: no need to modify oo_.exo_simul and oo_.exo_steady_state, since we simulated 100% of the shock (although with an approximation), except if controlled paths
    if ~isempty(oo_.deterministic_simulation.controlled_paths_by_period)
        oo_.exo_simul = exobase + (exo_simul - exobase)/completed_share;
    end

    maxerror = compute_maxerror(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);

    if ~options_.noprint
        fprintf('Perfect foresight solution found for %.1f%% of the shock, then extrapolation was performed using linearization\n\n', completed_share*100)
    end
    oo_.deterministic_simulation.status = true;
    oo_.deterministic_simulation.homotopy_linearization = true;
elseif options_.simul.homotopy_marginal_linearization_fallback > 0 && completed_share > options_.simul.homotopy_marginal_linearization_fallback
    oo_.deterministic_simulation.sim1.endo_simul = endo_simul;
    oo_.deterministic_simulation.sim1.exo_simul = exo_simul;
    oo_.deterministic_simulation.sim1.steady_state = steady_state;
    oo_.deterministic_simulation.sim1.exo_steady_state = exo_steady_state;
    oo_.deterministic_simulation.sim1.homotopy_completion_share = completed_share;

    % Now compute extra simulation. First try using the first simulation as guess value.
    if isempty(marginal_linearization_previous_raw_sims)
        shareorig = 1;
        endoorig = oo_.endo_simul;
        exoorig = oo_.exo_simul;
    else
        shareorig = marginal_linearization_previous_raw_sims.sim2.homotopy_completion_share;
        endoorig = marginal_linearization_previous_raw_sims.sim2.endo_simul;
        exoorig = shareorig*oo_.exo_simul + (1-shareorig)*exobase;
        if ~isempty(oo_.deterministic_simulation.controlled_paths_by_period)
            for p = 1:length(oo_.deterministic_simulation.controlled_paths_by_period)
                if isempty(oo_.deterministic_simulation.controlled_paths_by_period(p).exogenize_id)
                    continue
                end
                exo_ids = oo_.deterministic_simulation.controlled_paths_by_period(p).endogenize_id;
                exoorig(length(initperiods)+p,exo_ids) = marginal_linearization_previous_raw_sims.sim2.exo_simul(length(initperiods)+p,exo_ids);
            end
        end
    end
    extra_share = completed_share - options_.simul.homotopy_marginal_linearization_fallback;
    if ~options_.noprint
        fprintf('Only %.1f%% of the shock could be simulated. Since marginal linearization was requested as a fallback, now running an extra simulation for %.1f%% of the shock\n\n', completed_share*100, extra_share*100)
        fprintf('%s\n\n', repmat('*', 1, 80))
    end
    extra_simul_time_counter = tic;
    [extra_success, extra_endo_simul, extra_exo_simul, extra_steady_state, extra_exo_steady_state, extra_controlled_paths_by_period] = create_scenario(M_,options_,oo_,extra_share, shareorig, endoorig, exoorig, endobase, exobase, initperiods, lastperiods, recompute_final_steady_state, endo_simul, exo_simul, steady_state, exo_steady_state);
    if extra_success
        try
            [extra_endo_simul, extra_success, ~, ~, ~, extra_exo_simul] = perfect_foresight_solver_core(extra_endo_simul, extra_exo_simul, extra_steady_state, extra_exo_steady_state, extra_controlled_paths_by_period, M_, options_);
        catch ME
            if is_numerical_exception(ME)
                extra_success = false;
            else
                rethrow(ME)
            end
        end
    end
    if ~extra_success
        if ~options_.noprint
            fprintf('The extra simulation for %.1f%% of the shock did not run when using the first simulation as a guess value. Now trying a full homotopy loop to get that extra simulation working\n\n', extra_share*100)
            fprintf('%s\n\n', repmat('*', 1, 80))
        end
        [extra_completed_share, extra_endo_simul, extra_exo_simul, extra_steady_state, extra_exo_steady_state] = homotopy_loop(M_,options_,oo_,extra_share, shareorig, endoorig, exoorig, endobase, exobase, initperiods, simperiods, lastperiods, recompute_final_steady_state, oo_.steady_state, oo_.exo_steady_state);
        extra_success = (extra_completed_share == extra_share);
    end
    extra_simul_time_elapsed = toc(extra_simul_time_counter);
    if extra_success
        oo_.deterministic_simulation.sim2.endo_simul = extra_endo_simul;
        oo_.deterministic_simulation.sim2.exo_simul = extra_exo_simul;
        oo_.deterministic_simulation.sim2.steady_state = extra_steady_state;
        oo_.deterministic_simulation.sim2.exo_steady_state = extra_exo_steady_state;
        oo_.deterministic_simulation.sim2.homotopy_completion_share = extra_share;

        oo_.endo_simul = endo_simul + (endo_simul - extra_endo_simul)*(1-completed_share)/options_.simul.homotopy_marginal_linearization_fallback;
        if options_.simul.endval_steady
            % The following is needed for the STEADY_STATE() operator to work properly,
            % and thus must come before computing the maximum error.
            % This is not a true steady state, but it is the closest we can get to
            oo_.steady_state = oo_.endo_simul(:, end);
        end
        % NB: no need to modify oo_.exo_simul and oo_.exo_steady_state, since we simulated 100% of the shock (although with an approximation), except if controlled paths
        if ~isempty(oo_.deterministic_simulation.controlled_paths_by_period)
            oo_.exo_simul = exo_simul + (exo_simul - extra_exo_simul)*(1-completed_share)/options_.simul.homotopy_marginal_linearization_fallback;
        end

        maxerror = compute_maxerror(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, M_, options_);

        if ~options_.noprint
            fprintf('Perfect foresight solution found for %.1f%% of the shock, then extrapolation was performed using marginal linearization (extra simulation took %.1f seconds)\n\n', completed_share*100, extra_simul_time_elapsed)
        end
        oo_.deterministic_simulation.homotopy_marginal_linearization = true;
    else
        % Set oo_ values to the partial simulation
        oo_.endo_simul = endo_simul;
        oo_.exo_simul = exo_simul;
        oo_.steady_state = steady_state;
        oo_.exo_steady_state = exo_steady_state;
        fprintf('perfect_foresight_solver: marginal linearization failed, unable to find solution for %.1f%% of the shock (extra simulation took %.1f seconds). Try to modify the value of homotopy_marginal_linearization_fallback option\n\n', extra_share*100, extra_simul_time_elapsed)
    end
    oo_.deterministic_simulation.status = extra_success;
else
    % Set oo_ values to the partial simulation
    oo_.endo_simul = endo_simul;
    oo_.exo_simul = exo_simul;
    oo_.steady_state = steady_state;
    oo_.exo_steady_state = exo_steady_state;
    fprintf('Failed to solve perfect foresight model\n\n')
    oo_.deterministic_simulation.status = false;
end

oo_.deterministic_simulation.error = maxerror;
oo_.deterministic_simulation.homotopy_completion_share = completed_share;
oo_.deterministic_simulation.homotopy_iterations = iteration;
if ~isempty(solver_iter)
    oo_.deterministic_simulation.iterations = solver_iter;
end
if ~isempty(per_block_status)
    oo_.deterministic_simulation.block = per_block_status;
end

if nargout > 1
    ts = construct_simulation_dseries(oo_, M_, first_simulation_period);
end

oo_.gui.ran_perfect_foresight = oo_.deterministic_simulation.status;


function [completed_share, endo_simul, exo_simul, steady_state, exo_steady_state, iteration, maxerror, solver_iter, per_block_status] = homotopy_loop(M_,options_,oo_,max_share, shareorig, endoorig, exoorig, endobase, exobase, initperiods, simperiods, lastperiods, recompute_final_steady_state, steady_state, exo_steady_state)
% INPUTS
%   M_               [structure] describing the model
%   options_         [structure] describing the options
%   share            [double]    the share of the shock that we want to simulate
%   simperiods       [vector]    period indices of simulation periods (between initial and terminal conditions)
%   endoorig         [matrix]    path of endogenous corresponding to shareorig of the shock (also possibly used as guess value for first iteration if relevant)
%   exoorig          [matrix]    path of exogenous corresponding to shareorig of the shock (also possibly used as guess value for first iteration if controlled paths)
%   …                            other inputs have the same meaning as in the create_scenario function
%
% OUTPUTS
%   completed_share  [double]    the share that has been successfully computed
%   endo_simul       [matrix]    path of endogenous corresponding to completed share
%   exo_simul        [matrix]    path of exogenous corresponding to completed share
%   steady_state     [vector]    steady state of endogenous corresponding to the completed share (equal to the input if terminal steady state not recomputed)
%   exo_steady_state [vector]    steady state of exogenous corresponding to the completed share (equal to the input if terminal steady state not recomputed)
%   iteration        [integer]   number of homotopy iterations performed
%   maxerror         [double]    as returned by perfect_foresight_solver_core
%   solver_iter      [integer]   corresponds to iter as returned by perfect_foresight_solver_core
%   per_block_status [struct]    as returned by perfect_foresight_solver_core


completed_share = 0;  % Share of shock successfully completed so far
step = min(options_.simul.homotopy_initial_step_size, max_share);
success_counter = 0;
iteration = 0;

endo_simul = endoorig;
exo_simul = exoorig;
periods = get_simulation_periods(options_);

while step > options_.simul.homotopy_min_step_size

    iteration = iteration+1;

    saved_endo_simul = endo_simul;
    saved_exo_simul = exo_simul;

    new_share = completed_share + step; % Try this share, and see if it succeeds

    if new_share > max_share
        new_share = max_share; % Don't go beyond target point
        step = new_share - completed_share;
    end

    iter_time_counter = tic;

    [steady_success, endo_simul, exo_simul, steady_state, exo_steady_state, controlled_paths_by_period] = create_scenario(M_,options_,oo_,new_share, shareorig, endoorig, exoorig, endobase, exobase, initperiods, lastperiods, recompute_final_steady_state, endo_simul, exo_simul, steady_state, exo_steady_state);

    if steady_success
        % At the first iteration, use the initial guess given by
        % perfect_foresight_setup or the user (but only if new_share==shareorig, otherwise it
        % does not make much sense). Afterwards, until a converging iteration has been obtained,
        % use the rescaled terminal condition (or, if there is no lead, the base
        % scenario / initial steady state).
        if completed_share == 0
            if iteration == 1 && new_share == shareorig
                % Nothing to do, at this point endo_simul(:, simperiods) == endoorig(:, simperiods)
            elseif M_.maximum_lead > 0
                endo_simul(:, simperiods) = repmat(endo_simul(:, lastperiods(1)), 1, periods);
            else
                endo_simul(:, simperiods) = endobase(:, simperiods);
            end

            % Also handle initial guess in exo_simul when there are controlled paths:
            % first try the initial guess given by perfect_foresight_setup or the user,
            % afterwards use the base scenario
            if ~isempty(oo_.deterministic_simulation.controlled_paths_by_period)
                if iteration == 1 && new_share == shareorig
                    % Nothing to do, at this point exo_simul(:, simperiods) == exoorig(:, simperiods)
                else
                    exo_simul(simperiods,:) = exobase(simperiods,:);
                end
            end
        end

        % Solve for the paths of the endogenous variables.
        try
            [endo_simul, success, maxerror, solver_iter, per_block_status, exo_simul] = perfect_foresight_solver_core(endo_simul, exo_simul, steady_state, exo_steady_state, controlled_paths_by_period, M_, options_);
        catch ME
            % Catch some numerical-related exceptions, and treat them as failure for the current homotopy step
            if is_numerical_exception(ME)
                success = false;
                maxerror = NaN;
                solver_iter = [];
                per_block_status = [];
            else
                rethrow(ME)
            end
        end
    else
        success = false;
        maxerror = NaN;
        solver_iter = [];
        per_block_status = [];
    end

    iter_time_elapsed = toc(iter_time_counter);

    if options_.no_homotopy || (iteration == 1 && success && new_share == 1)
        % Skip homotopy
        if success
            completed_share = new_share;
        end
        break
    end

    if iteration == 1 && ~options_.noprint
        fprintf('\nEntering the homotopy method iterations...\n')
        iter_summary_table = { sprintf('\nIter. \t | Share \t | Status \t | Max. residual\t | Duration (sec)\n'),
                               sprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n') };
    end

    if success
        % Successful step
        if ~options_.noprint
            iter_summary_table{end+1} = sprintf('%i \t | %1.5f \t | %s \t | %e \t | %.1f\n', iteration, new_share, 'succeeded', maxerror, iter_time_elapsed);
        end
        completed_share = new_share;
        if completed_share >= max_share
            % Print the iterations summary table for the last time, to show convergence
            fprintf('%s', iter_summary_table{:})
            break
        end
        success_counter = success_counter + 1;
        if options_.simul.homotopy_step_size_increase_success_count > 0 ...
           && success_counter >= options_.simul.homotopy_step_size_increase_success_count
            success_counter = 0;
            step = step * 2;
        end
    else
        endo_simul = saved_endo_simul;
        exo_simul = saved_exo_simul;
        success_counter = 0;
        step = step / 2;
        if ~options_.noprint
            if ~steady_success
                iter_summary_table{end+1} = sprintf('%i \t | %1.5f \t | failed (in endval steady) \t\t | %.1f\n', iteration, new_share, iter_time_elapsed);
            elseif isreal(maxerror)
                iter_summary_table{end+1} = sprintf('%i \t | %1.5f \t | failed \t | %e \t | %.1f\n', iteration, new_share, maxerror, iter_time_elapsed);
            else
                iter_summary_table{end+1} = sprintf('%i \t | %1.5f \t | failed \t | Complex \t | %.1f\n', iteration, new_share, iter_time_elapsed);
            end
        end
    end

    % Print the iterations summary table at every iteration
    fprintf('%s', iter_summary_table{:})
end

%If simulated paths are complex, take real part and recompute the residuals to check whether this is actually a solution
if ~isreal(endo_simul(:)) % cannot happen with bytecode or the perfect_foresight_problem DLL
    real_simul = real(endo_simul);
    real_maxerror = compute_maxerror(real_simul, exo_simul, steady_state, M_, options_);
    if real_maxerror <= options_.dynatol.f
        endo_simul = real_simul;
        maxerror = real_maxerror;
    else
        completed_share = 0;
        disp('Simulation terminated with imaginary parts in the residuals or endogenous variables.')
    end
end

fprintf('\n')


function [steady_success, endo_simul, exo_simul, steady_state, exo_steady_state, controlled_paths_by_period] = create_scenario(M_,options_,oo_,share, shareorig, endoorig, exoorig, endobase, exobase, initperiods, lastperiods, recompute_final_steady_state, endo_simul, exo_simul, steady_state, exo_steady_state)
% For a given share, computes the exogenous path and also the initial and
% terminal conditions for the endogenous path (but do not modify the initial
% guess for endogenous)
%
% INPUTS
%   M_               [structure] describing the model
%   options_         [structure] describing the options
%   oo_              [structure] storing the results
%   share            [double]    the share of the shock that we want to simulate
%   shareorig        [double]    the share to which endoorig and exoorig correspond (typically 100%, except for perfect_foresight_with_expectation_errors_solver with homotopy and marginal linearization)
%   endoorig         [matrix]    path of endogenous corresponding to shareorig of the shock (only initial and terminal conditions are used, and the latter only if recompute_final_steady_state=false)
%   exoorig          [matrix]    path of exogenous corresponding to shareorig of the shock
%   endobase         [matrix]    path of endogenous corresponding to 0% of the shock (only initial and terminal conditions are used, except if oo_.deterministic_simulation.controlled_paths_by_period is not empty)
%   exobase          [matrix]    path of exogenous corresponding to 0% of the shock
%   initperiods      [vector]    period indices of initial conditions
%   lastperiods      [vector]    period indices of terminal conditions
%   recompute_final_steady_state [boolean] self-explanatory
%   endo_simul       [matrix]    path of endogenous, used to construct the guess values (initial condition not used; terminal condition used as guess value iff recompute_final_steady_state=true)
%   exo_simul        [matrix]    path of exogenous, used to construct the guess values (only if oo_.deterministic_simulation.controlled_paths_by_period is not empty)
%   steady_state     [vector]    steady state of endogenous, only used if terminal steady state is *not* recomputed by the function
%   exo_steady_state [vector]    steady state of exogenous, only used if terminal steady state is *not* recomputed by the function
%
% OUTPUTS
%   steady_success   [boolean]   whether the recomputation of the steady state was successful (always true if no recomputation was tried)
%   endo_simul       [matrix]    path of endogenous corresponding to the scenario
%   exo_simul        [matrix]    path of exogenous corresponding to the scenario
%   steady_state     [vector]    steady state of endogenous corresponding to the scenario (equal to the input if terminal steady state not recomputed)
%   exo_steady_state [vector]    steady state of exogenous corresponding to the scenario (equal to the input if terminal steady state not recomputed)
%   controlled_paths_by_period [struct] flips between endos and exos corresponding to the scenario

saved_exo_simul = exo_simul; % For guess value of controlled paths

% Compute convex combination for the path of exogenous
exo_simul = exoorig*share/shareorig + exobase*(1-share/shareorig);
if ~isempty(options_.simul.homotopy_exclude_varexo)
    [is_exo, excluded_exo_ids] = ismember(options_.simul.homotopy_exclude_varexo, M_.exo_names);
    if ~all(is_exo)
        error('Option homotopy_exclude_varexo must be given exogenous variable names')
    end
    exo_simul(:, excluded_exo_ids) = exoorig(:, excluded_exo_ids);
end

% Compute convex combination for the initial condition
% In most cases, the initial condition is a steady state and this does nothing
% This is for cases when the initial condition is out of equilibrium
endo_simul(:, initperiods) = share/shareorig*endoorig(:, initperiods)+(1-share/shareorig)*endobase(:, initperiods);

% If there is a permanent shock, ensure that the rescaled terminal condition is
% a steady state (if the user asked for this recomputation, or if the original
% terminal condition is a steady state)
steady_success = true;
if recompute_final_steady_state
    % Set “local” options for steady state computation (after saving the global values)
    saved_steady_solve_algo = options_.solve_algo;
    options_.solve_algo = options_.simul.steady_solve_algo;
    saved_steady_maxit = options_.steady.maxit;
    options_.steady.maxit = options_.simul.steady_maxit;
    saved_steady_tolf = options_.solve_tolf;
    options_.solve_tolf = options_.simul.steady_tolf;
    saved_steady_tolx = options_.solve_tolx;
    options_.solve_tolx = options_.simul.steady_tolx;
    saved_steady_markowitz = options_.markowitz;
    options_.markowitz = options_.simul.steady_markowitz;

    saved_ss = endo_simul(:, lastperiods);
    % Effectively compute the terminal steady state
    for j = lastperiods
        % First use the terminal steady of the previous homotopy iteration as guess value (or the contents of the endval block if this is the first iteration)
        [endo_simul(:, j), ~, info] = evaluate_steady_state(endo_simul(:, j), exo_simul(j, :)', M_, options_, ~options_.simul.endval_steady_nocheck);
        if info(1)
            % If this fails, then try again using the initial steady state as guess value
            if isempty(oo_.initial_steady_state)
                guess_value = oo_.steady_state;
            else
                guess_value = oo_.initial_steady_state;
            end
            [endo_simul(:, j), ~, info] = evaluate_steady_state(guess_value, exo_simul(j, :)', M_, options_, ~options_.simul.endval_steady_nocheck);
            if info(1)
                % If this fails again, give up and restore last periods in endo_simul
                endo_simul(:, lastperiods) = saved_ss;
                steady_success = false;
                break;
            end
        end
    end

    % The following is needed for the STEADY_STATE() operator to work properly
    steady_state = endo_simul(:, end);

    exo_steady_state = exo_simul(end, :)';

    options_.solve_algo = saved_steady_solve_algo;
    options_.steady.maxit = saved_steady_maxit;
    options_.solve_tolf = saved_steady_tolf;
    options_.solve_tolx = saved_steady_tolx;
    options_.markowitz = saved_steady_markowitz;
else
    % The terminal condition is not a steady state, compute a convex combination
    endo_simul(:, lastperiods) = share/shareorig*endoorig(:, lastperiods)+(1-share/shareorig)*endobase(:, lastperiods);
end

controlled_paths_by_period = oo_.deterministic_simulation.controlled_paths_by_period;
if ~isempty(controlled_paths_by_period)
    for p = 1:length(controlled_paths_by_period)
        if isempty(controlled_paths_by_period(p).exogenize_id)
            continue
        end
        controlled_paths_by_period(p).values = controlled_paths_by_period(p).values*share + endobase(controlled_paths_by_period(p).exogenize_id,p)'*(1-share);

        % Handle guess values for endogenized exos
        exo_ids = controlled_paths_by_period(p).endogenize_id;
        exo_simul(length(initperiods)+p,exo_ids) = saved_exo_simul(length(initperiods)+p,exo_ids);
    end
end


function check_input_arguments(options_, M_, oo_)

periods = get_simulation_periods(options_);

if options_.stack_solve_algo < 0 || options_.stack_solve_algo > 7
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: stack_solve_algo must be between 0 and 7')
end

if ~options_.block && ~options_.bytecode && ~ismember(options_.stack_solve_algo, [0:3 6 7])
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: you must use stack_solve_algo={0,1,2,3,6,7} when not using block nor bytecode option')
end

if options_.block && ~options_.bytecode && options_.stack_solve_algo == 5
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: you can''t use stack_solve_algo = 5 without bytecode option')
end

if isempty(oo_.endo_simul) || any(size(oo_.endo_simul) ~= [ M_.endo_nbr, M_.maximum_lag+periods+M_.maximum_lead ])

    if options_.initval_file
        fprintf('PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Check whether your initval-file provides %d periods.',M_.maximum_endo_lag+periods+M_.maximum_endo_lead)
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    else
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    end
end

if (M_.exo_nbr > 0) && ...
        (isempty(oo_.exo_simul) || any(size(oo_.exo_simul) ~= [ M_.maximum_lag+periods+M_.maximum_lead, M_.exo_nbr ]))
    if options_.initval_file
        fprintf('PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size. Check whether your initval-file provides %d periods.',M_.maximum_endo_lag+periods+M_.maximum_endo_lead)
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size.')
    else
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    end
end


function r = is_numerical_exception(ME)
% Deals with:
% - erf and erfc (under MATLAB, because under Octave they accept complex inputs)
% - normpdf and normcdf (under MATLAB, from missing/stats/; and under Octave from the statistics package)

r = (~isoctave && (strcmp(ME.identifier, 'MATLAB:erf:notFullReal') ...
                   || strcmp(ME.identifier, 'MATLAB:erfc:notFullReal'))) ...
    || (isoctave && (strcmp(ME.message, 'normcdf: X, MU, and SIGMA must not be complex.') ...
                     || strcmp(ME.message, 'normpdf: X, MU, and SIGMA must not be complex.')));
