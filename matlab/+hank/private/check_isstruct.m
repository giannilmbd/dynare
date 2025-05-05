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
% Check that an input is a structure.
%
% INPUTS
% - f [any type]: variable to check
% - f_name [char]: name to display in the error message
%
% OUTPUTS
% - (none) This function throws an error if the input is not a structure.
%
% DESCRIPTION
% Checks whether `f` is a structure.
% If not, throws an error indicating that the field `f_name` should be a structure.
function check_isstruct(f, f_name)
   if ~isstruct(f)
      error('Misspecified steady-state input `ss`: the `%s` field should be a structure.', f_name);
   end
end