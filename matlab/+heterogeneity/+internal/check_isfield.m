function check_isfield(f, s, f_name, details)
% Check the existence of a field in a structure.
%
% INPUTS
% - f [char]: name of the field to check
% - s [structure]: structure in which the field should be present
% - f_name [optional, char]: name to display in the error message (defaults to `f`)
% - details [optional, char]: additional details to append to the error message
%
% OUTPUTS
% - (none) This function throws an error if the specified field is missing.
%
% DESCRIPTION
% Checks whether the field `f` exists in the structure `s`.
% If not, throws an error indicating the missing field.
% If provided, `details` is appended to the error message for additional context.

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
   if nargin < 4
      details = '';
   end
   if nargin < 3
      f_name = f;
   end
   if ~isfield(s, f)
      error('Misspecified steady-state input `ss`: the `%s` field is missing.%s', f_name, details);
   end
end