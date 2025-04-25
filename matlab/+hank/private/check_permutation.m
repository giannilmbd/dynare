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
% Check and returns a permutation order for variables.
%
% INPUTS
% - s [structure]: structure containing the permutation information
% - f_name [char]: name of the field within `s` that should contain the permutation list
% - s_name [char]: full name of the structure `s` for error messages
% - symbs [cell array of char]: list of expected variable names
%
% OUTPUTS
% - out_order [cell array or string array]: permutation order of variables
%
% DESCRIPTION
% Checks that the field `f_name` exists in structure `s` and contains a valid
% permutation (i.e., a cell array of character vectors or a string array).
% If the field is missing, returns the default order given by `symbs`.
% Throws an error if the specified variables are not consistent with `symbs`.
function [out_order] = check_permutation(s, f_name, s_name, symbs)
   if ~isfield(s, f_name)
      out_order = symbs;
   else
      f = s.(f_name);
      if ~isstring(f) && ~iscellstr(f)
         error('Misspecified steady-state input `ss`: the `%s.%s` field should be a cell array of character vectors or a string array.', s_name, f_name);
      end
      err_var = setdiff(f, symbs);
      if ~isempty(err_var)
         error('Misspecified steady-state input `ss`: the set of variables of the `%s.%s` field is not consistent with the information in M_ and the other fields in the steady-state structure `ss`. Problematic variables: %s', s_name, f_name, strjoin(err_var));
      end
      out_order = f;
   end
end