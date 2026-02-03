function dr = solve(M_, options_solve, oo_het, varargin)
% Computes a first-order linearization around the model's aggregate steady state
%
% INPUTS
% - M_ [structure]: Dynare's model structure
% - options_solve [structure]: Options for the solver containing:
%     .truncation_horizon [integer]: Time horizon T for sequence-space Jacobians
% - oo_het [structure]: Heterogeneity results structure (from load_steady_state)
%     containing .mat, .sizes, .indices fields
% - varargin [optional]: custom impulse responses for testing purposes
%
% OUTPUTS
% - dr [structure]: Decision rule structure containing:
%   * G.(var).(shock) [T×T matrix]: General equilibrium Jacobians mapping
%       shock sequences to variable responses
%   * J_ha.(output).(input) [T×T matrix]: Heterogeneous agent Jacobians
%   * F.(output).(input) [T×T matrix]: Transition matrices
%   * J.(eq).(var) [vector]: Aggregate model Jacobian entries [timing, value]
%   * curlyYs.(input).(output) [T×1 vector]: Aggregate impact coefficients
%   * curlyDs.(input) [T×N_om array]: Distribution derivative coefficients
%
% DESCRIPTION
% Computes the linearized solution of a model featuring rich heterogeneity by
% constructing impulse response functions, distribution derivatives, and
% Jacobian matrices following the sequence-space Jacobian (SSJ) method.
%
% HELPER FUNCTIONS (defined at end of file)
% - compute_expected_pol: Computes expectation matrices for forward-looking dynamics
% - compute_heterogeneous_jacobian: Transforms transition matrices to HA Jacobian
% - compute_impulse_responses: Computes IRFs for heterogeneous agents
% - compute_shock_jacobian: Constructs shock Jacobian matrix (sparse)
% - compute_state_jacobian: Constructs state Jacobian matrix (sparse)
% - process_jacobian_block: Processes sparse Jacobian blocks into dr.J structure
% - populate_dr_timing_fields: Populates dr fields across timing categories
% - populate_dr_timing_fields_transposed: Populates dr with transposed structure
% - populate_dr_timing_fields_reshape_2d: Populates dr with reshape for curlyDs

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
T = int32(options_solve.truncation_horizon);
H_ = M_.heterogeneity(1);
if nargin > 3 && ~isempty(varargin{1})
    impulse_responses = varargin{1};
else
    impulse_responses = compute_impulse_responses(T, H_, oo_het);
end

mat = oo_het.mat;
sizes = oo_het.sizes;
indices = oo_het.indices;
N_Y = sizes.N_Y;
N_Ix = sizes.N_Ix;
N_sp = sizes.N_sp;

% H_.endo_nbr x N_Y x N_sp x T 
impulse_responses_perm = permute(impulse_responses, [1 2 4 3]);
% H_.endo_nbr x N_Y x T x N_sp
x_hat = reshape(impulse_responses_perm, [], N_sp) * mat.d.Phi;
% (H_.endo_nbr ⋅ N_Y ⋅ T) x N_om
x_hat = reshape(x_hat, H_.endo_nbr, N_Y.all, T, []);
% H_.endo_nbr x N_Y x T x N_om
x_hat = permute(x_hat, [4 3 2 1]);
% N_om x T x N_Y x H_.endo_nbr
px_hat = x_hat(:,:,:,H_.state_var);
% N_om x T x N_Y x n

curlyDs = compute_curlyDs(mat.d.ind, mat.d.w, mat.d.inv_h, px_hat, mat.d.hist, mat.Mu, mat.d.dims);
% N_om x (T*N_Y)
curlyDs = reshape(curlyDs, size(px_hat,1), T, N_Y.all);
% N_om x T x N_Y

curlyYs = reshape(sum(x_hat(:,:,:,indices.Ix.in_het).*reshape(mat.d.hist,[],1,1,1), 1), T, N_Y.all, N_Ix);

expectations = compute_expected_pol(mat, sizes, indices, T, N_Ix);

transition_matrices = compute_transition_matrices(curlyYs, expectations, curlyDs);

heterogeneous_jacobian = compute_heterogeneous_jacobian(transition_matrices);

state_var = find(M_.lead_lag_incidence(1,:));
state_jacobian = compute_state_jacobian(M_.endo_nbr, state_var, mat.dG, heterogeneous_jacobian, T, sizes, indices);

shock_jacobian = compute_shock_jacobian(M_.endo_nbr, M_.exo_nbr, mat.dG, heterogeneous_jacobian, T, sizes, indices);

G = -full(state_jacobian) \ full(shock_jacobian);

dr = struct;

for i=1:M_.endo_nbr
    var = M_.endo_names{i};
    for j=1:M_.exo_nbr
        inn = M_.exo_names{j};
        dr.G.(var).(inn) = G((i-1)*T+(1:T), (j-1)*T+(1:T));
    end
end

for j=1:N_Ix
    output = [upper(M_.heterogeneity_aggregates{j,1}) '_' H_.endo_names{M_.heterogeneity_aggregates{j,3}}];
    dr = populate_dr_timing_fields(dr, {'J_ha', 'F'}, output, ...
        {heterogeneous_jacobian(:,:,:,j), transition_matrices(:,:,:,j)}, indices, N_Y);
end

[a_i, a_j, a_v] = find(mat.dG(:,M_.endo_nbr+(1:M_.endo_nbr)));
[b_i, b_j, b_v] = find(mat.dG(:,(1:M_.endo_nbr)));
[c_i, c_j, c_v] = find(mat.dG(:,2*M_.endo_nbr+(1:M_.endo_nbr)));
[d_i, d_j, d_v] = find(mat.dG(:,3*M_.endo_nbr+1:end));

% Process Jacobian blocks: lead (timing=0), current (timing=-1), lag (timing=1), exo (timing=0)
dr = process_jacobian_block(dr, a_i, a_j, a_v, M_, 0, M_.endo_names);
dr = process_jacobian_block(dr, b_i, b_j, b_v, M_, -1, M_.endo_names);
dr = process_jacobian_block(dr, c_i, c_j, c_v, M_, 1, M_.endo_names);
dr = process_jacobian_block(dr, d_i, d_j, d_v, M_, 0, M_.exo_names);

for j=1:N_Ix
    output = H_.endo_names{indices.Ix.in_het(j)};
    dr = populate_dr_timing_fields_transposed(dr, 'curlyYs', output, curlyYs, indices, N_Y, j);
end

curlyDs = permute(curlyDs, [2 1 3]);
n_dims = numel(indices.order);
n_shocks = numel(indices.shocks);
n_states = numel(indices.states);
dims = zeros(1, n_dims+1);
dims(1) = T;
for i=1:n_shocks
    dims(1+i) = sizes.shocks.(indices.shocks{i});
end
for i=1:n_states
    dims(1+n_shocks+i) = sizes.d.states.(indices.states{i});
end

dr = populate_dr_timing_fields_reshape_2d(dr, 'curlyDs', curlyDs, dims, indices, N_Y);

end

function expectations = compute_expected_pol(mat, sizes, indices, T, N_Ix)
% Computes expectation matrices for heterogeneous agent dynamics.
%
% INPUTS
% - mat [structure]: matrices from steady-state computation
% - sizes [structure]: size information for various dimensions
% - indices [structure]: permutation and indexing information
% - T [scalar]: number of time periods
% - N_Ix [scalar]: number of aggregated heterogeneous variables
%
% OUTPUTS
% - expectations [4-D array]: expectation matrices (N_e × N_a × N_Ix × T-1)
%
% DESCRIPTION
% Computes the expectation operators needed for the forward-looking
% components of the heterogeneous agent dynamics.
N_om = sizes.d.N_a * sizes.N_e;
expectations = zeros(N_om, N_Ix, T-1);

for i = 1:N_Ix
    expectations(:,i,1) = mat.d.x_bar(indices.Ix.in_het(i),:);
end
expectations(:,:,1) = expectations(:,:,1) - sum(expectations(:,:,1), 1) / N_om;

for s = 2:T-1
    expectations(:,:,s) = compute_expected_y(mat.d.ind, mat.d.w, expectations(:,:,s-1), mat.Mu, mat.d.dims);
    expectations(:,:,s) = expectations(:,:,s) - sum(expectations(:,:,s), 1) / N_om;
end

end

function heterogeneous_jacobian = compute_heterogeneous_jacobian(transition_matrices)
% Computes heterogeneous agent Jacobian from transition matrices.
%
% INPUTS
% - transition_matrices [4-D array]: F matrices (T × T × N_Y × N_Ix)
%
% OUTPUTS
% - heterogeneous_jacobian [4-D array]: heterogeneous agent Jacobian (T × T × N_Y × N_Ix)
%
% DESCRIPTION
% Transforms the transition matrices into the heterogeneous agent Jacobian
% by accumulating effects over time periods.
T = size(transition_matrices, 1);
heterogeneous_jacobian = transition_matrices;
for t = 2:T
    heterogeneous_jacobian(2:T,t,:,:) = heterogeneous_jacobian(2:T,t,:,:) + heterogeneous_jacobian(1:T-1,t-1,:,:);
end

end

function impulse_responses = compute_impulse_responses(T, H_, oo_het)
% Computes impulse response functions for heterogeneous agents.
%
% INPUTS
% - T      [integer]:   Truncation horizon
% - H_     [structure]: Heterogeneous model structure
% - oo_het [structure]: Heterogeneous results structure (must contain steady state)
%
% OUTPUTS
% - impulse_responses [4-D array]: impulse responses (H_.endo_nbr × N_Y × N_sp × T)
%
% DESCRIPTION
% Computes the initial impulse response (period 0) and iteratively builds
% the transmission mechanism for subsequent periods using policy matrices.
mat = oo_het.mat;
sizes = oo_het.sizes;
indices = oo_het.indices;

N_sp = sizes.N_sp;
N_Y = sizes.N_Y.all;

impulse_responses = zeros(H_.endo_nbr*N_Y, N_sp, T);

ind_x = H_.endo_nbr + (1:H_.endo_nbr);
ind_xe = H_.endo_nbr + ind_x;

Ex = compute_expected_dx(mat.pol.x_bar, mat.pol.ind, mat.pol.w, mat.pol.inv_h, mat.Mu, mat.pol.dims);
Ae3 = pagemtimes(mat.dF(:, ind_xe, :), Ex);
f3 = mat.dF(:, ind_x, :);
col = H_.state_var;
f3(:, col, :) = f3(:, col, :) + Ae3;

RHS3 = cat(2, mat.dF(:, ind_xe, :), mat.dF(:, indices.Y.ind_in_dF, :));
Z3 = -pagemldivide(f3, RHS3);
cFx_next = Z3(:, 1:H_.endo_nbr, :);
x = Z3(:, H_.endo_nbr+1:end, :);

impulse_responses(:,:,1) = ((reshape(x, [], N_sp) / mat.pol.U) / mat.pol.L) * mat.pol.P;

for s = 2:T
    Ex_dash = impulse_responses(:,:,s-1) * mat.pol.Phi_e;
    x = pagemtimes(cFx_next, reshape(Ex_dash, H_.endo_nbr, N_Y, N_sp));
    impulse_responses(:,:,s) = ((reshape(x, [], N_sp) / mat.pol.U) / mat.pol.L) * mat.pol.P;
end

impulse_responses = reshape(impulse_responses, H_.endo_nbr, N_Y, N_sp, T);

end

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
d = dG(:, 3*endo_nbr+1:end);
[d_i, d_j, d_v] = find(d);
d_i = int32(d_i); d_j = int32(d_j);
d_nnz = nnz(d);
N_Ix = sizes.N_Ix;
N_Y = sizes.N_Y;

exact_nnz = d_nnz * T + N_Y.exo * N_Ix * T * T;

rows = zeros(exact_nnz, 1);
cols = zeros(exact_nnz, 1);
vals = zeros(exact_nnz, 1);
idx = 0;
one_to_T = int32((1:T).');

if d_nnz ~= 0
    n_identity = d_nnz * T;
    idx_range = 1:n_identity;
    rows(idx_range) = repelem((d_i-1)*T, T, 1) + repmat(one_to_T, d_nnz, 1);
    cols(idx_range) = repelem((d_j-1)*T, T, 1) + repmat(one_to_T, d_nnz, 1);
    vals(idx_range) = repelem(d_v, T, 1);
    idx = idx+n_identity;
end

[jc, jr] = meshgrid(one_to_T, one_to_T);
for i=1:N_Ix
    row_base = (indices.Ix.in_agg(i)-1)*T;
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
N_Ix = sizes.N_Ix;
N_Y = sizes.N_Y;

a = dG(:, endo_nbr+(1:endo_nbr));
b = dG(:, (1:endo_nbr));
c = dG(:, 2*endo_nbr+(1:endo_nbr));

[a_i, a_j, a_v] = find(a);
[b_i, b_j, b_v] = find(b);
[c_i, c_j, c_v] = find(c);
a_i = int32(a_i); a_j = int32(a_j);
b_i = int32(b_i); b_j = int32(b_j);
c_i = int32(c_i); c_j = int32(c_j);
a_nnz = nnz(a);
b_nnz = nnz(b);
c_nnz = nnz(c);

exact_nnz = a_nnz*T + b_nnz*(T-1) + c_nnz*(T-1) + N_Ix*((N_Y.lag+N_Y.lead)*T*(T-1)+N_Y.current*T*T);

rows = zeros(exact_nnz, 1, "int32");
cols = zeros(exact_nnz, 1, "int32");
vals = zeros(exact_nnz, 1);
idx = int32(0);
one_to_T = int32((1:T).');
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
    row_base = (indices.Ix.in_agg(i)-1)*T;
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

function dr = populate_dr_timing_fields_reshape_2d(dr, fieldname, data, dims, indices, N_Y)
% Populates decision rule fields with reshape operation for 3D curlyDs data.
%
% INPUTS
% - dr [structure]: decision rule structure to populate
% - fieldname [string]: name of field to populate (e.g., 'curlyDs')
% - data [array]: data matrix (T × N_om × N_Y.all)
% - dims [vector]: dimensions for reshaping
% - indices [structure]: indexing structure containing Y.names
% - N_Y [structure]: structure with timing counts (lag, current, lead, exo)
%
% OUTPUTS
% - dr [structure]: updated decision rule structure
in_Y = 0;

% Process lag timing
for i=1:N_Y.lag
    in_Y = in_Y+1;
    input = sprintf("%s_lag", indices.Y.names.lag{i});
    dr.(fieldname).(input) = reshape(data(:,:,in_Y), dims);
end

% Process current timing
for i=1:N_Y.current
    in_Y = in_Y+1;
    input = indices.Y.names.current{i};
    dr.(fieldname).(input) = reshape(data(:,:,in_Y), dims);
end

% Process lead timing
for i=1:N_Y.lead
    in_Y = in_Y+1;
    input = sprintf("%s_lead", indices.Y.names.lead{i});
    dr.(fieldname).(input) = reshape(data(:,:,in_Y), dims);
end

% Process exogenous timing
for i=1:N_Y.exo
    in_Y = in_Y+1;
    input = indices.Y.names.exo{i};
    dr.(fieldname).(input) = reshape(data(:,:,in_Y), dims);
end

end

function dr = populate_dr_timing_fields_transposed(dr, fieldname, output, data, indices, N_Y, j)
% Populates decision rule fields with transposed structure (input.output instead of output.input).
%
% INPUTS
% - dr [structure]: decision rule structure to populate
% - fieldname [string]: name of field to populate (e.g., 'curlyYs')
% - output [string]: output variable name
% - data [array]: data matrix (T × N_Y × N_Ix or similar)
% - indices [structure]: indexing structure containing Y.names
% - N_Y [structure]: structure with timing counts (lag, current, lead, exo)
% - j [scalar]: index for third dimension of data
%
% OUTPUTS
% - dr [structure]: updated decision rule structure
in_Y = 0;

% Process lag timing
for i=1:N_Y.lag
    in_Y = in_Y+1;
    input = sprintf("%s_lag", indices.Y.names.lag{i});
    dr.(fieldname).(input).(output) = data(:,in_Y,j);
end

% Process current timing
for i=1:N_Y.current
    in_Y = in_Y+1;
    input = indices.Y.names.current{i};
    dr.(fieldname).(input).(output) = data(:,in_Y,j);
end

% Process lead timing
for i=1:N_Y.lead
    in_Y = in_Y+1;
    input = sprintf("%s_lead", indices.Y.names.lead{i});
    dr.(fieldname).(input).(output) = data(:,in_Y,j);
end

% Process exogenous timing
for i=1:N_Y.exo
    in_Y = in_Y+1;
    input = indices.Y.names.exo{i};
    dr.(fieldname).(input).(output) = data(:,in_Y,j);
end

end

function dr = populate_dr_timing_fields(dr, fieldnames, output, data, indices, N_Y)
% Populates decision rule fields across all timing categories (lag, current, lead, exo).
%
% INPUTS
% - dr [structure]: decision rule structure to populate
% - fieldnames [cell array]: names of fields to populate (e.g., {'J_ha', 'F'})
% - output [string]: output variable name
% - data [cell array]: cell array of data matrices corresponding to fieldnames
% - indices [structure]: indexing structure containing Y.names
% - N_Y [structure]: structure with timing counts (lag, current, lead, exo)
%
% OUTPUTS
% - dr [structure]: updated decision rule structure
in_Y = 0;

% Process lag timing
for i=1:N_Y.lag
    in_Y = in_Y+1;
    input = sprintf("%s_lag", indices.Y.names.lag{i});
    for f=1:numel(fieldnames)
        dr.(fieldnames{f}).(output).(input) = data{f}(:,:,in_Y);
    end
end

% Process current timing
for i=1:N_Y.current
    in_Y = in_Y+1;
    input = indices.Y.names.current{i};
    for f=1:numel(fieldnames)
        dr.(fieldnames{f}).(output).(input) = data{f}(:,:,in_Y);
    end
end

% Process lead timing
for i=1:N_Y.lead
    in_Y = in_Y+1;
    input = sprintf("%s_lead", indices.Y.names.lead{i});
    for f=1:numel(fieldnames)
        dr.(fieldnames{f}).(output).(input) = data{f}(:,:,in_Y);
    end
end

% Process exogenous timing
for i=1:N_Y.exo
    in_Y = in_Y+1;
    input = indices.Y.names.exo{i};
    for f=1:numel(fieldnames)
        dr.(fieldnames{f}).(output).(input) = data{f}(:,:,in_Y);
    end
end

end