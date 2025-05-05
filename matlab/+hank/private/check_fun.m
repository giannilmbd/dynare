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
% Apply a validation function to a set of fields within a structure.
%
% INPUTS
% - f [structure]: structure containing fields to be validated
% - f_name [char]: name of the structure (for error messages)
% - symbs [cell array of strings]: list of field names to validate
% - fun [function handle]: validation function applied to each field
% - text [char]: descriptive text for the error message in case of validation failure
% - eq [optional, scalar]: if provided, the output of `fun` is compared to `eq`
%
% OUTPUTS
% - (none) This function throws an error if any field fails the validation.
%
% DESCRIPTION
% For each symbol in `symbs`, applies the validation function `fun` to the
% corresponding field in `f`. If an optional `eq` argument is provided, the
% output of `fun` is compared to `eq`. Throws a descriptive error listing all
% fields that do not satisfy the validation criteria.
function check_fun(f, f_name, symbs, fun, text, eq)
   elt_check = cellfun(@(s) fun(f.(s)), symbs);
   if nargin == 6
      elt_check = elt_check == eq;
   end
   if ~all(elt_check)
      error('Misspecified steady-state input `ss`: the following fields in `%s` %s: %s.', f_name, text, strjoin(symbs(~elt_check)));
   end
end