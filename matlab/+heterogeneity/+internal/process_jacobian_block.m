function dr = process_jacobian_block(dr, row_idx, col_idx, values, M_, timing, var_names)
% Processes sparse Jacobian block from mat.dG and populates dr.J structure.
%
% INPUTS
% - dr [structure]: decision rule structure to populate
% - row_idx [vector]: row indices from find()
% - col_idx [vector]: column indices from find()
% - values [vector]: values from find()
% - M_ [structure]: Dynare model structure
% - timing [scalar]: timing indicator (-1=lag, 0=current, 1=lead, NaN=exo)
% - var_names [cell array]: variable names to use (M_.endo_names or M_.exo_names)
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

for i=1:numel(row_idx)
    % Get equation name
    if row_idx(i) <= M_.orig_endo_nbr
        eq = regexprep(['res_' M_.equations_tags{row_idx(i),3}], '\W', '_');
    end

    % Get variable name
    var = var_names{col_idx(i)};

    % Build entry [timing, value]
    entry = [timing values(i)];

    % Assign to dr.J structure
    % For current timing (timing=-1), need to check if field exists and append
    if timing == -1
        if isfield(dr.J.(eq), var)
            dr.J.(eq).(var) = [dr.J.(eq).(var) ; entry];
        else
            dr.J.(eq).(var) = entry;
        end
    else
        % For other timings, direct assignment
        dr.J.(eq).(var) = entry;
    end
end

end