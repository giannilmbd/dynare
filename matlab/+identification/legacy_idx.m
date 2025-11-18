function idx_old = legacy_idx(idx_new, M_)
% Converts a dynamic Jacobian index column from the new (sparse) to the old
% (legacy) representation (see issue #1859 for the historical background)

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

tmp = idx_new - 3 * M_.endo_nbr;
if tmp < 1
    lli = M_.lead_lag_incidence';
    if M_.maximum_lag == 0
        idx_new = idx_new - M_.endo_nbr;
    end
    idx_old = lli(idx_new);
else
    idx_old = tmp + M_.nspred + M_.endo_nbr + M_.nsfwrd;
end
