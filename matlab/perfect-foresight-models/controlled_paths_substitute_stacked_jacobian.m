function A = controlled_paths_substitute_stacked_jacobian(A, y, y0, yT, exogenousvariables, steadystate, controlled_paths_by_period, M_)
% Given a stacked Jacobian for the perfect foresight problem, modifies some columns
% according to the perfect_foresight_controlled_paths block.

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

dynamic_g1 = str2func([M_.fname '.dynamic_g1']);

periods = numel(controlled_paths_by_period);

for p = 1:periods
    nsubst = length(controlled_paths_by_period(p).exogenize_id);
    if nsubst > 0
        if p == 1
            y3n = [ y0; y(1:2*M_.endo_nbr) ];
        elseif p == periods
            y3n = [ y(end-2*M_.endo_nbr+1:end); yT ];
        else
            y3n = y((p-2)*M_.endo_nbr+(1:3*M_.endo_nbr));
        end
        g1 = dynamic_g1(y3n, exogenousvariables(p+M_.maximum_lag,:), M_.params, steadystate, M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr);

        subst = [ zeros((p-1)*M_.endo_nbr, nsubst);
                  g1(:,3*M_.endo_nbr+controlled_paths_by_period(p).endogenize_id);
                  zeros((periods-p)*M_.endo_nbr, nsubst) ];
        A(:,(p-1)*M_.endo_nbr+controlled_paths_by_period(p).exogenize_id) = subst;
    end
end
