function [periods, first_simulation_period, last_simulation_period] = get_simulation_periods(options_)
% Given the options_ structure, returns the number of simulation periods, and
% also possibly the first and last simulation periods (as dates objects).
% Also verify the consistency of the three options if all are given by the user.

% Copyright © 2024 Dynare Team
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

first_simulation_period = options_.simul.first_simulation_period;
last_simulation_period = options_.simul.last_simulation_period;

if options_.periods > 0
    periods = options_.periods;
    if ~isempty(first_simulation_period) && ~isempty(last_simulation_period)
        if periods ~= last_simulation_period - first_simulation_period + 1
            error('Options first_simulation_period, last_simulation_period and periods are inconsistent')
        end
    elseif ~isempty(first_simulation_period) && isempty(last_simulation_period)
        last_simulation_period = first_simulation_period + periods - 1;
    elseif isempty(first_simulation_period) && ~isempty(last_simulation_period)
        first_simulation_period = last_simulation_period - periods + 1;
    end
else
    if isempty(options_.simul.first_simulation_period) || isempty(options_.simul.last_simulation_period)
        error('Dynare:periodsNotSet', 'Options first_simulation_period and last_simulation_period must be set when periods is not provided')
    end
    periods = last_simulation_period - first_simulation_period + 1;
end
