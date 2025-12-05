function dr = populate_dr_timing_fields(dr, fieldnames, output, data, indices, N_Y)
% Populates decision rule fields across all timing categories (lag, current, lead, exo).
%
% INPUTS
% - dr [structure]: decision rule structure to populate
% - fieldnames [cell array]: names of fields to populate (e.g., {'J_ha', 'F'})
% - output [string]: output variable name
% - data [cell array]: cell array of data matrices corresponding to fieldnames
% - indices [structure]: indexing structure containing Y.names
% - N_Y [structure]: structure with timing counts (lag, current, lead, exo)
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
    for f=1:numel(fieldnames)
        dr.(fieldnames{f}).(output).(input) = data{f}(:,:,in_Y);
    end
end

% Process current timing
for i=1:N_Y.current
    in_Y = in_Y+1;
    input = indices.Y.names.current{i};
    for f=1:numel(fieldnames)
        dr.(fieldnames{f}).(output).(input) = data{f}(:,:,in_Y);
    end
end

% Process lead timing
for i=1:N_Y.lead
    in_Y = in_Y+1;
    input = sprintf("%s_lead", indices.Y.names.lead{i});
    for f=1:numel(fieldnames)
        dr.(fieldnames{f}).(output).(input) = data{f}(:,:,in_Y);
    end
end

% Process exogenous timing
for i=1:N_Y.exo
    in_Y = in_Y+1;
    input = indices.Y.names.exo{i};
    for f=1:numel(fieldnames)
        dr.(fieldnames{f}).(output).(input) = data{f}(:,:,in_Y);
    end
end

end