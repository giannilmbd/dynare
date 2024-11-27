function ts = construct_simulation_dseries(oo_, M_, first_simulation_period)
% Returns a dseries object for the perfect foresight simulation.
% Merges it with the contents of the datafile option on the initval_file command if provided.
% first_simulation_period may be empty, in which case some default value will be used.

% Copyright © 2021-2024 Dynare Team
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

if isempty(first_simulation_period)
    if isfield(oo_, 'initval_series') && ~isempty(oo_.initval_series)
        first_simulation_period = oo_.initval_series.dates(1)+(M_.orig_maximum_lag-1);
    else
        first_simulation_period = dates(1,1);
    end
end

ts = dseries([transpose(oo_.endo_simul(1:M_.orig_endo_nbr,:)), oo_.exo_simul], ...
             first_simulation_period - M_.maximum_lag, [M_.endo_names(1:M_.orig_endo_nbr); M_.exo_names]);

if isfield(oo_, 'initval_series') && ~isempty(oo_.initval_series)
    names = ts.name;
    ts = merge(oo_.initval_series{names{:}}, ts);
end
