function state_jacobian = compute_state_jacobian(endo_nbr, state_var, dG, heterogeneous_jacobian, T, sizes, indices)
% Constructs the state jacobian matrix
%
% INPUTS
% - endo_nbr [integer]: number of aggregate endogenous variables
% - state_var [vector]: vector of state indices
% - dG [matrix]: aggregate Jacobian
% - heterogeneous_jacobian [4-D array]: heterogeneous agent jacobian (T × T × N_Y × N_Ix)
% - T [scalar]: number of time periods
% - sizes [structure]: size information for various dimensions
% - indices [structure]: permutation and indexing information
%
% OUTPUTS
% - state_jacobian [sparse matrix]: state jacobian matrix (endo_nbr*T × endo_nbr*T)
%
% DESCRIPTION
% Constructs the sparse state jacobian matrix that includes both aggregate
% and heterogeneous agent dynamics for the linearized model.

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

N_Ix = sizes.N_Ix;
N_Y = sizes.N_Y;

a = dG(:, endo_nbr+(1:endo_nbr));
b = dG(:, (1:endo_nbr));
c = dG(:, 2*endo_nbr+(1:endo_nbr));

[a_i, a_j, a_v] = find(a);
[b_i, b_j, b_v] = find(b);
[c_i, c_j, c_v] = find(c);
a_nnz = nnz(a);
b_nnz = nnz(b);
c_nnz = nnz(c);

exact_nnz = a_nnz*T + b_nnz*(T-1) + c_nnz*(T-1) + N_Ix*((N_Y.lag+N_Y.lead)*T*(T-1)+N_Y.current*T*T);

rows = zeros(exact_nnz, 1);
cols = zeros(exact_nnz, 1);
vals = zeros(exact_nnz, 1);
idx = 0;
one_to_T = (1:T).';
two_to_T = one_to_T(2:end);
one_to_Tm1 = one_to_T(1:end-1);

if a_nnz ~= 0
    n_identity = a_nnz * T;
    idx_range = (idx+1):(idx+n_identity);
    rows(idx_range) = repelem((a_i-1)*T, T, 1) + repmat(one_to_T, a_nnz, 1);
    cols(idx_range) = repelem((a_j-1)*T, T, 1) + repmat(one_to_T, a_nnz, 1);
    vals(idx_range) = repelem(a_v, T, 1);
    idx = idx + n_identity;
end

if b_nnz ~= 0
    n_lag = b_nnz * (T-1);
    idx_range = (idx+1):(idx+n_lag);
    rows(idx_range) = repelem((b_i-1)*T, T-1, 1) + repmat(two_to_T, b_nnz, 1);
    cols(idx_range) = repelem((b_j-1)*T, T-1, 1) + repmat(one_to_Tm1, b_nnz, 1);
    vals(idx_range) = repelem(b_v, T-1, 1);
    idx = idx + n_lag;
end

if c_nnz ~= 0
    n_lead = c_nnz * (T-1);
    idx_range = (idx+1):(idx+n_lead);
    rows(idx_range) = repelem((c_i-1)*T, T-1, 1) + repmat(one_to_Tm1, c_nnz, 1);
    cols(idx_range) = repelem((c_j-1)*T, T-1, 1) + repmat(two_to_T, c_nnz, 1);
    vals(idx_range) = repelem(c_v, T-1, 1);
    idx = idx + n_lead;
end

if N_Y.lag > 0
    [jc_lag, jr_lag] = meshgrid(one_to_Tm1, one_to_T);
end
if N_Y.current > 0
    [jc_current, jr_current] = meshgrid(one_to_T, one_to_T);
end
if N_Y.lead > 0
    [jc_lead, jr_lead] = meshgrid(two_to_T, one_to_T);
end
for i = 1:N_Ix
    row_base = (indices.Ix.in_endo(i)-1)*T;
    in_Y = 0;
    for j=1:N_Y.lag
        in_Y = in_Y+1;
        if ismember(indices.Y.ind_in_dG.lag(j), state_var)
            col_base = (indices.Y.ind_in_dG.lag(j)-1) * T;
            J_ij = heterogeneous_jacobian(:,two_to_T,in_Y,i);
            n_entries = T*(T-1);
            idx_range = idx+(1:n_entries);
            rows(idx_range) = row_base + jr_lag(:);
            cols(idx_range) = col_base + jc_lag(:);
            vals(idx_range) = -J_ij(:);
            idx = idx + n_entries;
        end
    end
    for j=1:N_Y.current
        in_Y = in_Y+1;
        col_base = (indices.Y.ind_in_dG.current(j)-1) * T;
        J_ij = heterogeneous_jacobian(:,:,in_Y,i);
        n_entries = T*T;
        idx_range = idx+(1:n_entries);
        rows(idx_range) = row_base + jr_current(:);
        cols(idx_range) = col_base + jc_current(:);
        vals(idx_range) = -J_ij(:);
        idx = idx + n_entries;
    end
    for j=1:N_Y.lead
        in_Y = in_Y+1;
        col_base = (indices.Y.ind_in_dG.lead(j)-1) * T;
        J_ij = heterogeneous_jacobian(:,one_to_Tm1,in_Y,i);
        n_entries = T*(T-1);
        idx_range = idx+(1:n_entries);
        rows(idx_range) = row_base + jr_lead(:);
        cols(idx_range) = col_base + jc_lead(:);
        vals(idx_range) = -J_ij(:);
        idx = idx + n_entries;
    end
end

state_jacobian = sparse(rows, cols, vals, endo_nbr*T, endo_nbr*T);

end