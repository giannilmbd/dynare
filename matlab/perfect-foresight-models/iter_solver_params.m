function [iter_tol, iter_maxit, gmres_restart] = iter_solver_params(options_, A, b)
% Computes tolerance used for gmres or bicgstab functions in perfect foresight solver.
% b is the RHS of the linear equation to solve. Also ensures that maxit and restart
% are within acceptable range.
%
% The tolerance passed to gmres and bicgstab is a relative one (‖Ax−b‖/‖b‖).
% However, we test the convergence of algorithms via options_.dynatol.f, which is an
% absolute error (‖Ax−b‖). Hence the need to rescale by the norm of the RHS (‖b‖).

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

iter_tol = options_.simul.iter_tol;
if isempty(iter_tol)
    iter_tol = options_.dynatol.f / max(abs(b)) / 10;
end

iter_maxit = min(options_.simul.iter_maxit, size(A, 1));
gmres_restart = min(options_.simul.gmres_restart, size(A, 1));
