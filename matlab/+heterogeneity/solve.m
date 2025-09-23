function dr = solve(M_, options_solve, oo_het, varargin)
% Computes a first-order linearization around the model's aggregate steady state
%
% INPUTS
% - M_ [structure]: Dynare's model structure
% - options_solve [structure]: Options for the solver
% - oo_het [structure]: Heterogeneity results structure (from heterogeneity_load_steady_state)
% - varargin [optional]: custom impulse responses for testing purposes
%
% OUTPUTS
% - dr [structure]: Decision rule structure containing:
%   * G: structured general equilibrium Jacobians
%   * J_ha: heterogeneous agent Jacobians
%   * F: transition matrices
%   * J: aggregate model Jacobians
%   * curlyYs: aggregate impact coefficients
%   * curlyDs: distribution derivative coefficients
%
% DESCRIPTION
% Computes the linearized solution of a model featuring rich heterogeneity by
% constructing impulse response functions, distribution derivatives, and
% Jacobian matrices.

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
T = options_solve.truncation_horizon;
if nargin > 3 && ~isempty(varargin{1})
    impulse_responses = varargin{1};
else
    impulse_responses = heterogeneity.internal.compute_impulse_responses(T, M_, oo_het);
end

mat = oo_het.mat;
sizes = oo_het.sizes;
indices = oo_het.indices;
N_Y = sizes.N_Y;
N_Ix = sizes.N_Ix;
N_x = sizes.N_x;
N_sp = sizes.N_sp;

% N_x x N_Y x N_sp x T 
impulse_responses_perm = permute(impulse_responses, [1 2 4 3]);
% N_x x N_Y x T x N_sp
x_hat = reshape(impulse_responses_perm, [], N_sp) * mat.d.Phi;
% (N_x ⋅ N_Y ⋅ T) x N_om
x_hat = reshape(x_hat, N_x, N_Y.all, T, []);
% N_x x N_Y x T x N_om
x_hat = permute(x_hat, [4 3 2 1]);
% N_om x T x N_Y x N_x
px_hat = x_hat(:,:,:,indices.x.ind.states);
% N_om x T x N_Y x n

curlyDs = compute_curlyDs(mat.d.ind, mat.d.w, mat.d.inv_h, px_hat, mat.d.hist, mat.Mu, indices.states, sizes.d.states);
% N_om x (T*N_Y)
curlyDs = reshape(curlyDs, size(px_hat,1), T, N_Y.all);
% N_om x T x N_Y

curlyYs = reshape(sum(x_hat(:,:,:,indices.Ix.in_x).*reshape(mat.d.hist,[],1,1,1), 1), T, N_Y.all, N_Ix);

expectations = heterogeneity.internal.compute_expected_pol(mat, sizes, indices, T, N_Ix);

transition_matrices = compute_transition_matrices(curlyYs, expectations, curlyDs);

heterogeneous_jacobian = heterogeneity.internal.compute_heterogeneous_jacobian(transition_matrices);

state_jacobian = heterogeneity.internal.compute_state_jacobian(M_.endo_nbr, M_.state_var, mat.dG, heterogeneous_jacobian, T, sizes, indices);

shock_jacobian = heterogeneity.internal.compute_shock_jacobian(M_.endo_nbr, M_.exo_nbr, mat.dG, heterogeneous_jacobian, T, sizes, indices);

G = -full(state_jacobian) \ full(shock_jacobian);

dr = struct;

for i=1:M_.endo_nbr
    var = M_.endo_names{i};
    for j=1:M_.exo_nbr
        inn = M_.exo_names{j};
        dr.G.(var).(inn) = G((i-1)*T+(1:T), (j-1)*T+(1:T));
    end
end

H_ = M_.heterogeneity(1);
for j=1:N_Ix
    output = [upper(M_.heterogeneity_aggregates{j,1}) '_' H_.endo_names{M_.heterogeneity_aggregates{j,3}}];
    dr = heterogeneity.internal.populate_dr_timing_fields(dr, {'J_ha', 'F'}, output, ...
        {heterogeneous_jacobian(:,:,:,j), transition_matrices(:,:,:,j)}, indices, N_Y);
end

[a_i, a_j, a_v] = find(mat.dG(:,M_.endo_nbr+(1:M_.endo_nbr)));
[b_i, b_j, b_v] = find(mat.dG(:,(1:M_.endo_nbr)));
[c_i, c_j, c_v] = find(mat.dG(:,2*M_.endo_nbr+(1:M_.endo_nbr)));
[d_i, d_j, d_v] = find(mat.dG(:,3*M_.endo_nbr+1:end));

% Process Jacobian blocks: lead (timing=0), current (timing=-1), lag (timing=1), exo (timing=0)
dr = heterogeneity.internal.process_jacobian_block(dr, a_i, a_j, a_v, M_, 0, M_.endo_names);
dr = heterogeneity.internal.process_jacobian_block(dr, b_i, b_j, b_v, M_, -1, M_.endo_names);
dr = heterogeneity.internal.process_jacobian_block(dr, c_i, c_j, c_v, M_, 1, M_.endo_names);
dr = heterogeneity.internal.process_jacobian_block(dr, d_i, d_j, d_v, M_, 0, M_.exo_names);

for j=1:N_Ix
    output = indices.x.names{indices.Ix.in_x(j)};
    dr = heterogeneity.internal.populate_dr_timing_fields_transposed(dr, 'curlyYs', output, curlyYs, indices, N_Y, j);
end

curlyDs = permute(curlyDs, [2 1 3]);
n_dims = numel(indices.order);
n_shocks_endo = numel(indices.shocks.endo);
n_shocks_exo = numel(indices.shocks.exo);
n_states = numel(indices.states);
dims = zeros(1, n_dims+1);
dims(1) = T;
for i=1:n_shocks_endo
    dims(1+i) = sizes.shocks.(indices.shocks.endo{i});
end
for i=1:n_shocks_exo
    dims(1+n_shocks_endo+i) = sizes.shocks.(indices.shocks.exo{i});
end
for i=1:n_states
    dims(1+n_shocks_endo+n_shocks_exo+i) = sizes.d.states.(indices.states{i});
end

dr = heterogeneity.internal.populate_dr_timing_fields_reshape_2d(dr, 'curlyDs', curlyDs, dims, indices, N_Y);

end