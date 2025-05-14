function iter_solver_error_flag(flag)
% Interprets error flag from gmres or bicgstab functions

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

switch flag
    case 1
        error('Maximum number of iterations exceeded in GMRES/BiCGStab. You may want to increase the iter_maxit option or the gmres_restart option (if using stack_solve_algo=2 for the latter), or decrease iter_tol.')
    case 2
        error('The preconditioner matrix is ill conditioned in GMRES/BiCGStab. You should try another preconditioner.')
    case 3
        error('No progress between two successive iterations in GMRES/BiCGStab. You may want to increase the iter_tol option.')
    case 4
        error('One of the scalar quantities calculated by GMRES/BiCGStab became too small or too large. Try another preconditioner or algorithm.')
end
