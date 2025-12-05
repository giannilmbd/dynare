function shock_jacobian = compute_shock_jacobian(endo_nbr, exo_nbr, dG, heterogeneous_jacobian, T, sizes, indices)
% Constructs the shock Jacobian matrix for the model featuring rich
% heterogeneity.
%
% INPUTS
% - endo_nbr [integer]: number of aggregate endogenous variables
% - exo_nbr [integer]: number of aggregate exogenous variables
% - dG [matrix]: aggregate Jacobian
% - M_ [structure]: Dynare's model structure
% - mat [structure]: matrices from steady-state computation
% - heterogeneous_jacobian [4-D array]: heterogeneous agent Jacobian (T × T × N_Y × N_Ix)
% - T [scalar]: number of time periods
% - sizes [structure]: size information for various dimensions
% - indices [structure]: permutation and indexing information
%
% OUTPUTS
% - shock_jacobian [sparse matrix]: shock Jacobian matrix (N*T × M*T)
%
% DESCRIPTION
% Constructs the sparse shock Jacobian matrix that maps exogenous shocks
% to endogenous variables in the linearized model.

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

d = dG(:, 3*endo_nbr+1:end);
[d_i, d_j, d_v] = find(d);
d_nnz = nnz(d);
N_Ix = sizes.N_Ix;
N_Y = sizes.N_Y;

exact_nnz = d_nnz * T + N_Y.exo * N_Ix * T * T;

rows = zeros(exact_nnz, 1);
cols = zeros(exact_nnz, 1);
vals = zeros(exact_nnz, 1);
idx = 0;

if d_nnz ~= 0
    one_to_T = (1:T)';
    n_identity = d_nnz * T;
    idx_range = 1:n_identity;
    rows(idx_range) = repelem((d_i-1)*T, T, 1) + repmat(one_to_T, d_nnz, 1);
    cols(idx_range) = repelem((d_j-1)*T, T, 1) + repmat(one_to_T, d_nnz, 1);
    vals(idx_range) = repelem(d_v, T, 1);
    idx = idx+n_identity;
end

one_to_T = (1:T).';
[jc, jr] = meshgrid(one_to_T, one_to_T);
for i=1:N_Ix
    row_base = (indices.Ix.in_endo(i)-1)*T;
    in_Y = N_Y.all-N_Y.exo;
    for j=1:N_Y.exo
        in_Y = in_Y+1;
        col_base = (indices.Y.ind_in_dG.exo(j)-1) * T;
        J_ij = heterogeneous_jacobian(:,:,in_Y,i);
        n_entries = T*T;
        idx_range = idx+(1:n_entries);
        rows(idx_range) = row_base + jr(:);
        cols(idx_range) = col_base + jc(:);
        vals(idx_range) = -J_ij(:);
        idx = idx + n_entries;
    end
end

shock_jacobian = sparse(rows, cols, vals, endo_nbr*T, exo_nbr*T);

end