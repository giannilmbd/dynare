function result = compute_nullspace(matrix, jacobian_tolerance)
% result = compute_nullspace(matrix, jacobian_tolerance)
% Helper function to compute null space with appropriate method
% based on MATLAB/Octave version and tolerance settings

% Copyright © 1996-2026 Dynare Team
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

if isempty(jacobian_tolerance) && isoctave
    result = null(matrix);
elseif ~isempty(jacobian_tolerance) && (isoctave || ~matlab_ver_less_than('9.12'))
    result = null(matrix, jacobian_tolerance);
else %use rational basis in MATLAB if no tolerance specified or Matlab version is too old
    if matlab_ver_less_than('9.12')
        result = null(matrix, 'r'); %old syntax for rational basis
    else
        result = null(matrix, "rational");
    end
end