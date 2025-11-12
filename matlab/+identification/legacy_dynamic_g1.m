function g1_legacy = legacy_dynamic_g1(g1, M_)
% Given g1 as returned by dynamic_g1.m, construct the g1 matrix in the legacy
% representation (see issue #1859 for the historical background). Ideally this
% file should go away when the identification code is adapted to the new model
% representation.

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

idx = find(M_.lead_lag_incidence')';
if M_.maximum_lag == 0
    idx = M_.endo_nbr + idx;
end

g1_legacy = full(g1(:,[idx 3*M_.endo_nbr+(1:(M_.exo_nbr+M_.exo_det_nbr))]));
end
