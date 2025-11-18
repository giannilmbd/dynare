function g4_legacy = legacy_dynamic_g4(g4_v, M_)
% Given g4 as returned by dynamic_g4.m, construct the g4 matrix in the legacy
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

nvar = M_.nspred + M_.endo_nbr + M_.nsfwrd + M_.exo_nbr + M_.exo_det_nbr;
nnz = size(M_.dynamic_g4_sparse_indices, 1);

g4_i = int32(zeros(nnz, 1));
g4_j = int32(zeros(nnz, 1));

for k = 1:length(g4_v)
    eq = M_.dynamic_g4_sparse_indices(k,1);
    var1 = identification.legacy_idx(M_.dynamic_g4_sparse_indices(k,2), M_) - 1;
    var2 = identification.legacy_idx(M_.dynamic_g4_sparse_indices(k,3), M_) - 1;
    var3 = identification.legacy_idx(M_.dynamic_g4_sparse_indices(k,4), M_) - 1;
    var4 = identification.legacy_idx(M_.dynamic_g4_sparse_indices(k,5), M_) - 1;

    g4_i(k) = eq;
    g4_j(k) = ((var1 * nvar + var2) * nvar + var3) * nvar + var4 + 1;
end

g4_legacy = sparse(g4_i, g4_j, g4_v, M_.endo_nbr, nvar ^ 4);
