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
% Check consistency between two sets of variable names.
%
% INPUTS
% - f1 [cell array of strings]: list of variable names to be checked for consistency
% - f2 [cell array of strings]: reference list of valid variable names
% - f1_name [char]: name of the first set (for error messages)
% - f2_name [char]: name of the second set (for error messages)
%
% OUTPUTS
% - (none) This function throws an error if there are elements in `f1` that are not present in `f2`.
%
% DESCRIPTION
% Checks that all elements of `f1` are included in `f2`. If not, throws a descriptive error 
% mentioning the missing variables and the corresponding structure names for easier debugging.
function check_consistency(f1, f2, f1_name, f2_name)
   f1_in_f2 = ismember(f1,f2);
   if ~all(f1_in_f2)
      error('Misspecified steady-state input `ss`: the following variables are mentioned in `%s`, but are missing in `%s`: %s', f1_name, f2_name, strjoin(f1(~f1_in_f2)));
   end
end