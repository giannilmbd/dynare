function dr = populate_dr_timing_fields_transposed(dr, fieldname, output, data, indices, N_Y, j)
% Populates decision rule fields with transposed structure (input.output instead of output.input).
%
% INPUTS
% - dr [structure]: decision rule structure to populate
% - fieldname [string]: name of field to populate (e.g., 'curlyYs')
% - output [string]: output variable name
% - data [array]: data matrix (T × N_Y × N_Ix or similar)
% - indices [structure]: indexing structure containing Y.names
% - N_Y [structure]: structure with timing counts (lag, current, lead, exo)
% - j [scalar]: index for third dimension of data
%
% OUTPUTS
% - dr [structure]: updated decision rule structure

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

in_Y = 0;

% Process lag timing
for i=1:N_Y.lag
    in_Y = in_Y+1;
    input = sprintf("%s_lag", indices.Y.names.lag{i});
    dr.(fieldname).(input).(output) = data(:,in_Y,j);
end

% Process current timing
for i=1:N_Y.current
    in_Y = in_Y+1;
    input = indices.Y.names.current{i};
    dr.(fieldname).(input).(output) = data(:,in_Y,j);
end

% Process lead timing
for i=1:N_Y.lead
    in_Y = in_Y+1;
    input = sprintf("%s_lead", indices.Y.names.lead{i});
    dr.(fieldname).(input).(output) = data(:,in_Y,j);
end

% Process exogenous timing
for i=1:N_Y.exo
    in_Y = in_Y+1;
    input = indices.Y.names.exo{i};
    dr.(fieldname).(input).(output) = data(:,in_Y,j);
end

end