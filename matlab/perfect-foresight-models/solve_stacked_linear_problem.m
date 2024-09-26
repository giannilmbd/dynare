function [endogenousvariables, success] = solve_stacked_linear_problem(endogenousvariables, exogenousvariables, steadystate_y, steadystate_x, M_, options_)

% Copyright © 2015-2024 Dynare Team
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

if M_.maximum_lag > 0
    y0 = endogenousvariables(:, M_.maximum_lag);
else
    y0 = NaN(M_.endo_nbr, 1);
end
if M_.maximum_lead > 0
    yT = endogenousvariables(:, M_.maximum_lag+options_.periods+1);
else
    yT = NaN(M_.endo_nbr, 1);
end
z = endogenousvariables(:,M_.maximum_lag+(1:options_.periods));

% Evaluate the residuals and Jacobian of the dynamic model at the deterministic steady state.
y3n = repmat(steadystate_y, 3, 1);
[d1, TT_order, TT] = feval([M_.fname,'.sparse.dynamic_resid'], y3n, steadystate_x', M_.params, ...
                           steadystate_y);
jacobian = feval([M_.fname,'.sparse.dynamic_g1'], y3n, steadystate_x', M_.params, steadystate_y, ...
                 M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, ...
                 M_.dynamic_g1_sparse_colptr, TT_order, TT);

% Check that the dynamic model was evaluated at the steady state.
if ~options_.steadystate.nocheck && max(abs(d1))>1e-12
    error('Jacobian is not evaluated at the steady state!')
end

z = bsxfun(@minus, z, steadystate_y);
x = bsxfun(@minus, exogenousvariables, steadystate_x');

[y, check, ~, ~, errorcode] = dynare_solve(@linear_perfect_foresight_problem, z(:), ...
                                           options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, ...
                                           options_, ...
                                           jacobian, y0-steadystate_y, yT-steadystate_y, ...
                                           x, M_.params, steadystate_y, ...
                                           M_.maximum_lag, options_.periods, M_.endo_nbr);

if all(imag(y)<.1*options_.dynatol.x)
    if ~isreal(y)
        y = real(y);
    end
else
    check = 1;
end

endogenousvariables = [y0 bsxfun(@plus,reshape(y,M_.endo_nbr,options_.periods), steadystate_y) yT];

success = ~check;

if ~success && options_.debug
    dprintf('solve_stacked_linear_problem: Nonlinear solver routine failed with errorcode=%i.', errorcode)
end
