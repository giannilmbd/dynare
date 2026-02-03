function v = extract_from_database(database_id, varname, period, first_simulation_period, M_)
% Helper for the “shock_paths” block, used when variables are referenced from a
% database namespace.

% Copyright © 2026 Dynare Team
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

dbname = M_.database{database_id};
db = evalin('base', dbname);

if isa(db, 'dseries')
    if isempty(first_simulation_period)
        error('%s is a dseries database is used but neither first_simulation_period nor last_simulation_period option was passed', dbname)
    end
    v = db{varname}(first_simulation_period+period-1).data;
elseif isa(db, 'table')
    if size(db, 1) == 1
        % If single-row table, ignore the period information
        v = db{1, varname};
    else
        v = db{period, varname};
    end
else
    error('Unsupported database type for %s', dbname)
end
