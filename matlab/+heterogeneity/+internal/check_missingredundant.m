function check_missingredundant(f, f_name, symbs, no_warning_redundant)
% Check for missing and redundant fields in a structure.
%
% INPUTS
% - f [structure]: structure to check
% - f_name [char]: name to display in error/warning messages
% - symbs [cell array of char]: list of expected field names
% - no_warning_redundant [logical]: if true, suppresses warnings about redundant fields
%
% OUTPUTS
% - (none) This function throws an error if any expected field is missing,
%   and issues a warning if redundant fields are found (unless no_warning_redundant is true).
%
% DESCRIPTION
% Verifies that all fields listed in `symbs` are present in the structure `f`.
% Throws an error if any expected field is missing.
% If `no_warning_redundant` is false, checks if `f` contains fields not listed
% in `symbs`, and issues a warning if redundant fields are found.

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
% Original author: Normann Rion <normann@dynare.org>
   fields = fieldnames(f);
   symbs_in_fields = ismember(symbs, fields);
   if ~all(symbs_in_fields)
      error('Misspecified steady-state input `ss`. The following fields are missing in `%s`: %s.', f_name, strjoin(symbs(~symbs_in_fields)));
   end
   if ~no_warning_redundant
      fields_in_symbs = ismember(fields, symbs);
      if ~all(fields_in_symbs)
         warning('Steady-state input `ss`. The following fields are redundant in `%s`: %s.', f_name, strjoin(fields(~fields_in_symbs)));
      end
   end
end