function [y, success, maxerror, iter, per_block_status, exo_simul] = perfect_foresight_solver_core(y, exo_simul, steady_state, exo_steady_state, controlled_paths_by_period, M_, options_)

% Core function calling solvers for perfect foresight model
%
% INPUTS
% - y                   [matrix] initial path of endogenous (typically oo_.endo_simul)
% - exo_simul           [matrix] path of exogenous
% - steady_state        [vector] steady state of endogenous variables
% - exo_steady_state    [vector] steady state of exogenous variables
% - controlled_paths_by_period [struct] describes flips between endos and exos, typically
%                          extracted from perfect_foresight_controlled_paths block;
%                          Can be an empty array if no such flip is requested
% - M_                  [struct] contains a description of the model.
% - options_            [struct] contains various options.
%
% OUTPUTS
% - y                   [double array] path for the endogenous variables (solution)
% - success             [logical] Whether a solution was found
% - maxerror            [double] contains the maximum absolute error
% - iter                [integer] Number of iterations of the underlying nonlinear solver (empty for non-iterative methods)
% - per_block_status    [struct] In the case of block decomposition, provides per-block solver status information (empty if no block decomposition)

% Copyright © 2015-2025 Dynare Team
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

if options_.lmmcp.status
    options_.stack_solve_algo=7;
    options_.solve_algo = 10;
elseif options_.stack_solve_algo==7 && options_.solve_algo == 11
    options_.lmmcp.status = 1; %Path solver
end

periods = get_simulation_periods(options_);

if options_.linear_approximation
    if ~ismember(options_.stack_solve_algo, [0 7])
        error('perfect_foresight_solver: option linear_approximation is only available with stack_solve_algo option equal to 0 or 7')
    end
    if options_.block
        error('perfect_foresight_solver: option linear_approximation is not available with the block option')
    end
    if options_.bytecode
        error('perfect_foresight_solver: option linear_approximation is not available with the bytecode option')
    end
    if M_.maximum_endo_lead == 0 || M_.maximum_endo_lag == 0
        error('perfect_foresight_solver: option linear_approximation is not available with purely backward, purely forward or static models')
    end
end

if options_.endogenous_terminal_period
    if options_.stack_solve_algo ~= 0
        error('perfect_foresight_solver: option endogenous_terminal_period is only available with stack_solve_algo option equal to 0')
    end
    if options_.block
        error('perfect_foresight_solver: option endogenous_terminal_period is not available with the block option')
    end
    if options_.bytecode
        error('perfect_foresight_solver: option endogenous_terminal_period is not available with the bytecode option')
    end
    if M_.maximum_endo_lead == 0 || M_.maximum_endo_lag == 0
        error('perfect_foresight_solver: option endogenous_terminal_period is not available with purely backward, purely forward or static models')
    end
end

if ~isempty(controlled_paths_by_period)
    assert(nargout >= 6); % Ensure modified exos are used
    if ~ismember(options_.stack_solve_algo, [0 1 2 3 6 7])
        error('perfect_foresight_controlled_paths is only available with stack_solve_algo option equal to 0, 1, 2, 3, 6 or 7')
    end
    if options_.stack_solve_algo == 7 && ismember(options_.solve_algo, [10, 11])
        error('perfect_foresight_controlled_paths is not available for mixed-complementarity problems (LMMCP or PATH solvers)')
    end
    if options_.bytecode
        error('perfect_foresight_controlled_paths is not available with the bytecode option')
    end
    if options_.block
        error('perfect_foresight_controlled_paths is not available with the block option')
    end
    if M_.maximum_endo_lead == 0 || M_.maximum_endo_lag == 0
        error('perfect_foresight_controlled_paths is not available with purely backward, purely forward or static models')
    end
end


if options_.linear && ismember(options_.stack_solve_algo, [0 7]) && ~options_.block ...
        && ~options_.bytecode && M_.maximum_endo_lead > 0 && M_.maximum_endo_lag > 0
    options_.linear_approximation = true;
end

maxerror = [];
iter = [];
per_block_status = [];

if options_.block
    if M_.block_structure.time_recursive
        error('Internal error: can''t perform stacked perfect foresight simulation with time-recursive block decomposition')
    end
    if options_.bytecode && ~ismember(options_.stack_solve_algo, [1 6])
        try
            y = bytecode('dynamic', 'block_decomposed', M_, options_, y, exo_simul, M_.params, steady_state, periods);
            success = true;
        catch ME
            if options_.verbosity >= 1
                disp(ME.message)
            end
            success = false;
        end
    else
        [y, success, maxerror, per_block_status] = solve_block_decomposed_problem(y, exo_simul, steady_state, options_, M_);
    end
else
    if options_.bytecode && ~ismember(options_.stack_solve_algo, [1 6])
        try
            y = bytecode('dynamic', M_, options_, y, exo_simul, M_.params, repmat(steady_state, 1, periods+2), periods);
            success = true;
        catch ME
            if options_.verbosity >= 1
                disp(ME.message)
            end
            success = false;
        end
    else
        if M_.maximum_endo_lead == 0 && M_.maximum_endo_lag>0 && ~options_.lmmcp.status % Purely backward model
            [y, success] = sim1_purely_backward(y, exo_simul, steady_state, M_, options_);
        elseif M_.maximum_endo_lag == 0 && M_.maximum_endo_lead>0 && ~options_.lmmcp.status % Purely forward model
            [y, success] = sim1_purely_forward(y, exo_simul, steady_state, M_, options_);
        elseif M_.maximum_endo_lag == 0 && M_.maximum_endo_lead == 0 && ~options_.lmmcp.status % Purely static model
            [y, success] = sim1_purely_static(y, exo_simul, steady_state, M_, options_);
        else % General case
            switch options_.stack_solve_algo
              case {0 2 3}
                if options_.linear_approximation
                    [y, success, maxerror, exo_simul] = sim1_linear(y, exo_simul, steady_state, exo_steady_state, controlled_paths_by_period, M_, options_);
                else
                    [y, success, maxerror, iter, exo_simul] = sim1(y, exo_simul, steady_state, controlled_paths_by_period, M_, options_);
                end
              case {1 6}
                [y, success, maxerror, iter, exo_simul] = sim1_lbj(y, exo_simul, steady_state, controlled_paths_by_period, M_, options_);
              case 7
                if options_.linear_approximation
                    if options_.solve_algo == 10 || options_.solve_algo == 11
                        warning('Since you are requesting an MCP solver, you should not specify your model as model(linear) or use the linear_approximation option!')
                    end
                    [y, success] = solve_stacked_linear_problem(y, exo_simul, steady_state, exo_steady_state, M_, options_);
                else
                    [y, success, maxerror, exo_simul] = solve_stacked_problem(y, exo_simul, steady_state, controlled_paths_by_period, M_, options_);
                end
              otherwise
                error('Invalid value of stack_solve_algo option!')
            end
        end
    end
end

% Some solvers do not compute the maximum error, so do it here if needed
if nargout > 2 && isempty(maxerror)
    if options_.bytecode
        residuals = bytecode('dynamic', 'evaluate', M_, options_, y, exo_simul, M_.params, steady_state, periods);
    else
        ny = size(y, 1);
        if M_.maximum_lag > 0
            y0 = y(:, M_.maximum_lag);
        else
            y0 = NaN(ny, 1);
        end
        if M_.maximum_lead > 0
            yT = y(:, M_.maximum_lag+periods+1);
        else
            yT = NaN(ny, 1);
        end
        yy = y(:,M_.maximum_lag+(1:periods));

        residuals = perfect_foresight_problem(yy(:), y0, yT, exo_simul, M_.params, steady_state, periods, M_, options_);
    end
    maxerror = max(max(abs(residuals)));
end
