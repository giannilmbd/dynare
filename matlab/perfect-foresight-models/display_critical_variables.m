function display_critical_variables(y, endo_names, caller, noprint, period_offset, iteration)
% display_critical_variables(y, endo_names, caller, noprint,period_offset,iteration)
% INPUTS
% - y                      [double]    nvars×periods array
% - endo_names             [string]    cell array of variable names
% - caller                 [string]    name of caller function
% - noprint                [Boolean]   indicator whether output should be printed
% - period_offset          [integer]   offset for period display (to start counting with 0)
% - iteration              [integer]   current iteration
%
% OUTPUTS
%  None
%
% ALGORITHM


% Copyright © 2025-2026 Dynare Team
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


if noprint
    return
end
if nargin<5
    period_offset=0;
end
if nargin<6
    iteration=[];
end

if any(any(isnan(y)))
    [variable_index,period_index] = find(isnan(y));
    if isempty(iteration)
        fprintf('%s: Simulation provided NaN for the following variables:\n',caller)
    else
        fprintf('%s: Iteration %u provided NaN for the following variables:\n',caller,iteration)
    end
    for var_iter=1:length(variable_index)
        fprintf('%s in period %u\n', endo_names{variable_index(var_iter)},period_index(var_iter)-period_offset)
    end
    fprintf('\n'),
end
if any(any(isinf(y)))
    [variable_index,period_index] = find(isinf(y));
    if isempty(iteration)
        fprintf('%s: Simulation diverged (Inf) for the following variables:\n',caller)
    else
        fprintf('%s: Iteration %u diverged (Inf) for the following variables:\n',caller,iteration)
    end
    for var_iter=1:length(variable_index)
        fprintf('%s in period %u\n', endo_names{variable_index(var_iter)},period_index(var_iter)-period_offset)
    end
    fprintf('\n'),
end

return % --*-- Unit tests --*--

%@test:1
y=[NaN,Inf,1;Inf, 1, NaN; 1 1 NaN];
endo_names={'y','c','k'};
caller='unit test'
noprint=false;
period_offset=0;
try
    display_critical_variables(y, endo_names, caller, noprint);
    t(1)=true;
catch
    t(1)=false
end

period_offset=1;
iteration=1;
try
    display_critical_variables(y, endo_names, caller, noprint,period_offset,iteration);
    t(2)=true;
catch
    t(2)=false
end

T = all(t);
%@eof:1