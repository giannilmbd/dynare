function [endogenousvariables, success, maxerror, exogenousvariables] = solve_stacked_problem(endogenousvariables, exogenousvariables, steadystate, controlled_paths_by_period, M_, options_)
% [endogenousvariables, success, maxerror] = solve_stacked_problem(endogenousvariables, exogenousvariables, steadystate, M_, options_)
% Solves the perfect foresight model using dynare_solve
%
% INPUTS
% - endogenousvariables [double] N*(T+M_.maximum_lag+M_.maximum_lead) array, paths for the endogenous variables (initial guess).
% - exogenousvariables  [double] (T+M_.maximum_lag+M_.maximum_lead)*M array, paths for the exogenous variables.
% - steadystate         [double] N*1 array, steady state for the endogenous variables.
%   - controlled_paths_by_period [struct] data from perfect_foresight_controlled_paths block
% - M_                   [struct] contains a description of the model.
% - options_             [struct] contains various options.
%
% OUTPUTS
% - endogenousvariables [double] N*(T+M_.maximum_lag+M_.maximum_lead) array, paths for the endogenous variables (solution of the perfect foresight model).
% - success             [logical] Whether a solution was found
% - maxerror            [double] 1-norm of the residual
% - exogenousvariables  [double] (T+M_.maximum_lag+M_.maximum_lead)*M array, paths for the exogenous variables
%                                (may be modified if perfect_foresight_controlled_paths present)

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

periods = get_simulation_periods(options_);

if M_.maximum_lag > 0
    y0 = endogenousvariables(:, M_.maximum_lag);
else
    y0 = NaN(M_.endo_nbr, 1);
end
if M_.maximum_lead > 0
    yT = endogenousvariables(:, M_.maximum_lag+periods+1);
else
    yT = NaN(M_.endo_nbr, 1);
end
z = endogenousvariables(:,M_.maximum_lag+(1:periods));

if (options_.solve_algo == 10 || options_.solve_algo == 11)% mixed complementarity problem
    assert(isempty(controlled_paths_by_period))
    [lb, ub] = feval(sprintf('%s.dynamic_complementarity_conditions', M_.fname), M_.params);
    if options_.linear_approximation
        lb = lb - steadystate_y;
        ub = ub - steadystate_y;
    end
    if options_.solve_algo == 10
        options_.lmmcp.lb = repmat(lb,periods,1);
        options_.lmmcp.ub = repmat(ub,periods,1);
    elseif options_.solve_algo == 11
        options_.mcppath.lb = repmat(lb,periods,1);
        options_.mcppath.ub = repmat(ub,periods,1);
    end
    dynamic_resid_function = str2func([M_.fname,'.sparse.dynamic_resid']);
    dynamic_g1_function = str2func([M_.fname,'.sparse.dynamic_g1']);
    [y, check, res, ~, errorcode] = dynare_solve(@perfect_foresight_mcp_problem, z(:), ...
                                                 options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                                 options_, ...
                                                 dynamic_resid_function, dynamic_g1_function, y0, yT, ...
                                                 exogenousvariables, M_.params, steadystate, ...
                                                 M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr, ...
                                                 M_.maximum_lag, periods, M_.endo_nbr, ...
                                                 M_.dynamic_mcp_equations_reordering);
    eq_to_ignore=find(isfinite(lb) | isfinite(ub));

elseif isempty(controlled_paths_by_period)
    [y, check, res, ~, errorcode] = dynare_solve(@perfect_foresight_problem, z(:), ...
                                               options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                               options_, y0, yT, exogenousvariables, M_.params, steadystate, periods, M_, options_);
else % controlled paths
    for p = 1:periods
        if isempty(controlled_paths_by_period(p).exogenize_id)
            continue
        end
        z(controlled_paths_by_period(p).exogenize_id,p) = exogenousvariables(p+M_.maximum_lag,controlled_paths_by_period(p).endogenize_id);
    end

    [y, check, res, ~, errorcode] = dynare_solve(@controlled_paths_wrapper, z(:), ...
                                                 options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                                 options_, y0, yT, exogenousvariables, M_.params, steadystate, periods, controlled_paths_by_period, M_, options_);

    for p = 1:periods
        if isempty(controlled_paths_by_period(p).exogenize_id)
            continue
        end
        exogenousvariables(p+M_.maximum_lag,controlled_paths_by_period(p).endogenize_id) = y(controlled_paths_by_period(p).exogenize_id+(p-1)*M_.endo_nbr);
        y(controlled_paths_by_period(p).exogenize_id+(p-1)*M_.endo_nbr) = controlled_paths_by_period(p).values;
    end
end

if all(imag(y)<.1*options_.dynatol.x)
    if ~isreal(y)
        y = real(y);
    end
else
    check = 1;
end

endogenousvariables(:, M_.maximum_lag+(1:periods)) = reshape(y, M_.endo_nbr, periods);
residuals=zeros(size(endogenousvariables));
residuals(:, M_.maximum_lag+(1:periods)) = reshape(res, M_.endo_nbr, periods);
if (options_.solve_algo == 10 || options_.solve_algo == 11)% mixed complementarity problem
    residuals(eq_to_ignore,bsxfun(@le, endogenousvariables(eq_to_ignore,:), lb(eq_to_ignore)+eps) | bsxfun(@ge,endogenousvariables(eq_to_ignore,:),ub(eq_to_ignore)-eps))=0;
end
maxerror = max(max(abs(residuals)));
success = ~check;

if ~success && options_.debug
    dprintf('solve_stacked_problem: Nonlinear solver routine failed with errorcode=%i.', errorcode)
end

end


function [res, A] = controlled_paths_wrapper(z, y0, yT, exogenousvariables, params, steadystate, periods, controlled_paths_by_period, M_, options_)

y = z;

for p = 1:periods
    if isempty(controlled_paths_by_period(p).exogenize_id)
        continue
    end
    exogenousvariables(p+M_.maximum_lag,controlled_paths_by_period(p).endogenize_id) = z(controlled_paths_by_period(p).exogenize_id+(p-1)*M_.endo_nbr);
    y(controlled_paths_by_period(p).exogenize_id+(p-1)*M_.endo_nbr) = controlled_paths_by_period(p).values;
end

[res, A] = perfect_foresight_problem(y, y0, yT, exogenousvariables, params, steadystate, periods, M_, options_);

A = controlled_paths_substitute_stacked_jacobian(A, y, y0, yT, exogenousvariables, steadystate, controlled_paths_by_period, M_);

end
