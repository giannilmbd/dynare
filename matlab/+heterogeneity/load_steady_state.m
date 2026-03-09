function oo_het = load_steady_state(M_, options_het, oo_het, steady_state, flag_initial_guess)
% Wrapper to validate, complete, and store steady-state input for models
% featuring rich heterogeneity.
%
% This function performs validation of the user-provided steady-state structure `steady_state`, constructs
% interpolation objects and basis matrices, computes missing policy functions for complementarity
% multipliers (if not supplied), and populates the global `oo_het` structure with the results.
%
% INPUTS
%   M_       [struct] : Dynare model structure
%   options_het [struct] : Heterogeneity-specific options structure. When steady_state
%                        is not passed as 4th argument:
%                        - .steady_state_file_name: path to .mat file. If empty,
%                          the variable is loaded from the base workspace instead.
%                        - .steady_state_variable_name: variable name in the .mat
%                          file (when file_name is non-empty) or in the base
%                          workspace (when file_name is empty).
%   oo_het      [struct] : Heterogeneity-specific output structure to which results will be written
%   steady_state (optional) [struct] : User-provided steady-state structure. When provided,
%                        file_name and variable_name options are ignored.
%   flag_initial_guess (optional) [bool] : If true, computes only policy-related matrices
%                        and returns early. Skips distribution-dependent computations
%                        (mat.d.ind/w/inv_h, Ix, mat.G/dG, mat.F/dF). Default: false.
%
% OUTPUTS
%   oo_het      [struct] : Updated output structure containing:
%                        - oo_het.steady_state    : validated and completed steady-state structure
%                        - oo_het.sizes : grid size structure
%                        - oo_het.mat : interpolation and policy function matrices

% Copyright © 2025-2026 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify it under the terms of the
% GNU General Public License as published by the Free Software Foundation, either version 3 of
% the License, or (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with Dynare. If not,
% see <https://www.gnu.org/licenses/>.
%
% Original author: Normann Rion <normann@dynare.org>
   assert(isscalar(M_.heterogeneity), 'Heterogeneous-agent models with more than one heterogeneity dimension are not allowed yet!');

   if nargin < 4
      file_name = options_het.steady_state_file_name;
      variable_name = options_het.steady_state_variable_name;

      if isempty(file_name)
         % No file specified — load variable from base workspace
         if ~evalin('base', ['exist(''' variable_name ''', ''var'')'])
            error('heterogeneity_load_steady_state: variable ''%s'' not found in workspace. Provide either ''filename'' or ''variable'' option.', variable_name);
         end
         steady_state = evalin('base', variable_name);
      else
         % File specified — load from .mat file
         [filepath, basename, extension] = fileparts(file_name);

         % Auto-detect extension if not provided
         if isempty(extension)
            if ~isempty(filepath)
               base_path = fullfile(filepath, basename);
            else
               base_path = basename;
            end

            if isfile([base_path '.mat'])
               extension = '.mat';
            else
               error('heterogeneity_load_steady_state: Cannot find file: %s.mat', base_path);
            end
            file_name_to_load = [base_path extension];
         else
            file_name_to_load = file_name;
         end

         if ~isfile(file_name_to_load)
            error('heterogeneity_load_steady_state: File not found: %s', file_name_to_load);
         end

         loaded_data = load(file_name_to_load);
         if ~isfield(loaded_data, variable_name)
            error('heterogeneity_load_steady_state: Variable ''%s'' not found in file ''%s''.', ...
                  variable_name, file_name_to_load);
         end
         steady_state = loaded_data.(variable_name);
      end
   end

   if nargin < 5
      flag_initial_guess = false;
   end

   % Check steady-state input
   [steady_state, sizes, indices] = heterogeneity.check_steady_state_input(M_, options_het.check, steady_state, flag_initial_guess);
   % Compute useful indices and names lists
   [indices, sizes] = compute_model_indices(M_, indices, sizes);
   % Compute Φ and ̃Φ with useful basis matrices
   mat = struct;
   Phi = speye(sizes.N_e);
   Phi_tilde = speye(sizes.N_sp);
   % Allocate policy interpolation matrices
   mat.pol.ind = zeros(sizes.N_sp, sizes.n_a, 'int32');
   mat.pol.w = zeros(sizes.N_sp, sizes.n_a);
   mat.pol.inv_h = zeros(sizes.N_sp, sizes.n_a);
   mat.d.ind = zeros(sizes.d.N_a, sizes.n_a, 'int32');
   mat.d.w = zeros(sizes.d.N_a, sizes.n_a);
   n_repmat = sizes.d.N_a;
   n_repelem = 1;
   for i=1:sizes.n_a
      state = indices.states{i};
      grid = steady_state.pol.grids.(state);
      [ind, w] = find_bracket_linear_weight(grid, steady_state.d.grids.(state));
      n_repmat = n_repmat/sizes.d.states.(state);
      mat.d.ind(:,i) = repmat(repelem(ind, n_repelem, 1), n_repmat, 1);
      mat.d.w(:,i) = repmat(repelem(w, n_repelem, 1), n_repmat, 1);
      n_repelem = n_repelem*sizes.d.states.(state);
      Phi = kron(build_linear_splines_basis_matrix(ind, w, sizes.pol.states.(state)), Phi);
      [ind, w] = find_bracket_linear_weight(grid, steady_state.pol.values.(state)(:));
      mat.pol.ind(:, i) = ind;
      mat.pol.w(:, i) = w;
      mat.pol.inv_h(:, i) = 1 ./ (grid(ind+1) - grid(ind));
   end
   % Create dims vectors for MEX functions
   mat.pol.dims = int32(cellfun(@(s) sizes.pol.states.(s), indices.states));
   mat.d.dims = int32(cellfun(@(s) sizes.d.states.(s), indices.states));
   mat.d.Phi = Phi;
   mat.pol.Phi = Phi_tilde;
   % Compute shock transition matrix
   H_ = M_.heterogeneity(1);
   mat.Mu = eye(1);
   for i=sizes.n_e:-1:1
      mat.Mu = kron(mat.Mu, steady_state.shocks.Pi.(H_.exo_names{i}));
   end
   % Compute state and shock grids for the policy state grids
   mat.pol.sm = set_state_matrix(steady_state.pol.grids, steady_state.shocks.grids, sizes, indices.shocks, indices.states);
   % LU decomposition of Phi_tilde
   [mat.pol.L, mat.pol.U, mat.pol.P, mat.pol.Q] = lu(Phi_tilde, 'vector');
   % Computation of Phi_tilde_e
   [I_mex, J_mex, V_mex] = compute_Phi_tilde_e(mat.pol.ind, mat.pol.w, mat.pol.dims, mat.Mu);
   mat.pol.Phi_e = sparse(I_mex, J_mex, V_mex, sizes.N_sp, sizes.N_sp);
   % Compute policy matrices without the auxiliary variables
   [x_bar, x_bar_dash] = compute_pol_matrices(steady_state.pol.values, H_.orig_endo_nbr, sizes.N_sp, mat.pol.U, mat.pol.L, mat.pol.P, mat.pol.Q, H_.endo_names(1:H_.orig_endo_nbr));
   % Compute y, x, yh, xh once for reuse throughout the function
   % - Aggregate endogenous variable vector at the steady state - %
   y = NaN(M_.endo_nbr,1);
   for i=1:M_.orig_endo_nbr
      var = M_.endo_names{i};
      y(i) = steady_state.agg.(var);
   end
   % - Exogenous variables - %
   mat.x = zeros(M_.exo_nbr,1);
   % - Auxiliary variables - %
   yagg = NaN(sizes.N_Ix,1);
   if M_.set_auxiliary_variables
      fun_set_auxiliary_variables = str2func([M_.fname '.set_auxiliary_variables']);
      y = fun_set_auxiliary_variables(y, mat.x, yagg, M_.params);
   end
   mat.y = repmat(y, 3, 1);
   % Build yh and xh arrays
   yh = NaN(3*H_.endo_nbr, sizes.N_sp);
   xh = NaN(H_.exo_nbr, sizes.N_sp);
   pol_ind = 1:H_.orig_endo_nbr;
   % t-1: states from state matrix
   yh(H_.state_var, :) = mat.pol.sm(sizes.n_e + (1:sizes.n_a), :);
   % t: current policy values
   yh(H_.endo_nbr + pol_ind, :) = x_bar;
   % t+1: expected future values
   yh(2*H_.endo_nbr + pol_ind, :) = x_bar_dash * mat.pol.Phi_e;
   % Exogenous shocks
   if ~isempty(indices.shocks)
      xh = mat.pol.sm(1:sizes.n_e, :);
   end
   % Compute auxiliary policy values using the preprocessor-generated function
   % Level-by-level computation with Phi_e between levels ensures correct
   % expectations for chain auxiliaries (e.g., AUX2(t) = AUX1(t+1))
   if H_.set_auxiliary_variables
      n_aux = H_.endo_nbr - H_.orig_endo_nbr;
      aux_range = H_.orig_endo_nbr+1:H_.endo_nbr;
      set_aux_fn = str2func([M_.fname '.dynamic_het1_set_auxiliary_variables']);
      pol_shape = size(steady_state.pol.values.(H_.endo_names{1}));
      n_levels = numel(H_.het_aux_levels);

      % Process each topological level
      for lv = 0:n_levels-1
         level_vars = H_.het_aux_levels{lv+1};  % 1-based MATLAB indexing

         % 1. Compute this level's aux at t (pointwise over all grid points)
         for j = 1:sizes.N_sp
            yh(:, j) = set_aux_fn(mat.y, mat.x, M_.params, [], yh(:, j), xh(:, j), [], lv);
         end

         % 2. Build policy for this level's computed variables
         for i = 1:length(level_vars)
            idx = level_vars(i);
            aux_name = H_.endo_names{idx};
            steady_state.pol.values.(aux_name) = reshape(yh(H_.endo_nbr + idx, :), pol_shape);
         end

         % 3. Apply Phi_e to get E[aux(t+1)] for this level's computed variables
         [lv_x_bar, lv_x_bar_dash] = compute_pol_matrices(...
            steady_state.pol.values, length(level_vars), sizes.N_sp, ...
            mat.pol.U, mat.pol.L, mat.pol.P, mat.pol.Q, H_.endo_names(level_vars));
         yh(H_.endo_nbr + level_vars, :) = lv_x_bar;
         yh(2*H_.endo_nbr + level_vars, :) = lv_x_bar_dash * mat.pol.Phi_e;
      end

      % Build final aux_x_bar and aux_x_bar_dash for all aux variables
      [aux_x_bar, aux_x_bar_dash] = compute_pol_matrices(...
         steady_state.pol.values, n_aux, sizes.N_sp, ...
         mat.pol.U, mat.pol.L, mat.pol.P, mat.pol.Q, H_.endo_names(aux_range));

      % Concatenate x_bar and x_bar_dash
      x_bar = [x_bar ; aux_x_bar];
      x_bar_dash = [x_bar_dash ; aux_x_bar_dash];
   end
   % Store the final x_bar and x_bar_dash in mat.pol
   mat.pol.x_bar = x_bar;
   mat.pol.x_bar_dash = x_bar_dash;
   % Extract state grids
   mat.pol.grids_array = cell(1, sizes.n_a);
   mat.d.grids_array = cell(1, sizes.n_a);
   for k = 1:sizes.n_a
      mat.pol.grids_array{k} = steady_state.pol.grids.(indices.states{k});
      mat.d.grids_array{k} = steady_state.d.grids.(indices.states{k});
   end
   % Early return for initial guess mode - skip distribution-dependent computations
   % When flag_initial_guess is true:
   %   - mat.pol.* (policy interpolation matrices) are computed
   %   - mat.Mu (shock transition matrix) is computed
   %   - mat.d.Phi (distribution interpolation basis) is computed
   %   - SKIPS mat.d.hist-dependent computations (mat.d.ind/w/inv_h, Ix, mat.G/dG, mat.F/dF)
   if flag_initial_guess
      oo_het.steady_state = steady_state;
      oo_het.sizes = sizes;
      oo_het.mat = mat;
      oo_het.indices = indices;
      return
   end
   % Generate relevant distribution-related objects
   mat.d.x_bar = mat.pol.x_bar_dash * mat.d.Phi;
   N_om = sizes.N_e * sizes.d.N_a;
   % Allocate distribution interpolation matrices
   mat.d.ind = zeros(N_om, sizes.n_a, 'int32');
   mat.d.w = zeros(N_om, sizes.n_a);
   mat.d.inv_h = zeros(N_om, sizes.n_a);
   for i = 1:sizes.n_a
      state = indices.states{i};
      grid = steady_state.d.grids.(state);
      [ind, w] = find_bracket_linear_weight(grid, mat.d.x_bar(H_.state_var(i), :));
      mat.d.ind(:, i) = ind;
      mat.d.w(:, i) = w;
      inv_h = 1 ./ (grid(ind+1) - grid(ind));
      mat.d.inv_h(:, i) = inv_h(:);
   end
   % Get inputs for Dynare residual and Jacobian functions
   mat.d.hist = reshape(steady_state.d.hist, [], 1);
   % Aggregate inputs with actual Ix values
   Ix = mat.pol.x_bar_dash * mat.d.Phi * mat.d.hist;
   yagg = Ix(indices.Ix.in_het);
   % Update the auxiliary variables accordingly
   if M_.set_auxiliary_variables
      fun_set_auxiliary_variables = str2func([M_.fname '.set_auxiliary_variables']);
      y = fun_set_auxiliary_variables(mat.y(M_.endo_nbr+(1:M_.endo_nbr)), mat.x, yagg, M_.params);
   end
   y = repmat(y, 3, 1);
   mat.y = y;
   % Call to the aggregate residual and Jacobian functions
   if M_.orig_endo_nbr > 0
      mat.G = feval([M_.fname '.dynamic_resid'], mat.y, mat.x, M_.params, [], yagg);
      mat.dG = feval([M_.fname '.dynamic_g1'], mat.y, mat.x, M_.params, [], yagg, ...
                      M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, ...
                      M_.dynamic_g1_sparse_colptr);
   end
   % Simple loop-based residual computation (for comparison with MEX)
   mat.F = NaN(H_.endo_nbr, sizes.N_sp);
   het_resid_fn = str2func([M_.fname '.dynamic_het1_resid']);
   for j = 1:sizes.N_sp
      mat.F(:, j) = het_resid_fn(mat.y, mat.x, M_.params, [], yh(:, j), xh(:, j), []);
   end
   % Call the heterogeneous Jacobian function (still using the old loop for now)
   mat.dF = NaN(H_.endo_nbr, 3*H_.endo_nbr+H_.exo_nbr+3*M_.endo_nbr+M_.exo_nbr, sizes.N_sp);
   het_jac_fn = str2func([M_.fname '.dynamic_het1_g1']);
   for j = 1:sizes.N_sp
      mat.dF(:,:,j) = feval(het_jac_fn, mat.y, mat.x, M_.params, [], yh(:,j), xh(:,j), [], H_.dynamic_g1_sparse_rowval, H_.dynamic_g1_sparse_colval, H_.dynamic_g1_sparse_colptr);
   end

   % Store the results in oo_
   oo_het.steady_state = steady_state;
   oo_het.sizes = sizes;
   oo_het.mat = mat;
   oo_het.indices = indices;
end

function sm = set_state_matrix(pol_grids, shocks_grids, sizes, shocks, states)
% Builds the full tensor-product grid of state and shock variables.
%
% This function constructs a matrix containing all combinations of state and shock values
% used for evaluating policy functions or distributions. The result is a matrix of size
% (n_a + n_e) × (N_a * N_e), where each column corresponds to one grid point in the joint
% state space.
%
% INPUTS
%   pol_grids    [struct] : Steady-state structure containing states grids for policy functions
%   shocks_grids  [struct] : Steady-state structure containing shocks grids
%   sizes  [struct] : Structure containing sizes of endogenous states and exogenous shocks
%   field  [char]   : Field name of `ss` ('pol' or 'd'), indicating whether to construct
%                     the grid for policy functions or for the distribution
%
% OUTPUT
%   sm     [matrix] : Matrix of size (n_a + n_e) × N_sp, where each column contains
%                     one combination of state and shock values, ordered lexicographically.
%
% Notes:
% - The first `n_e` rows correspond to shocks; the next `n_a` to endogenous states.
% - Ordering is compatible with Kronecker-product-based interpolation logic.
% - Used in residual evaluation, interpolation, and expectation computation routines.
   % Setting the states matrix
   sm = zeros(sizes.N_sp, sizes.n_a+sizes.n_e);

   % Filling the state matrix
   n_repmat = sizes.N_sp;
   n_repelem = 1;
   for j=1:sizes.n_e
      shock = shocks{j};
      n_repmat = n_repmat/sizes.shocks.(shock);
      sm(:,j) = repmat(repelem(shocks_grids.(shock), n_repelem, 1), n_repmat, 1);
      n_repelem = n_repelem*sizes.shocks.(shock);
   end
   for j=1:sizes.n_a
      state = states{j};
      n_repmat = n_repmat/sizes.pol.states.(state);
      sm(:,sizes.n_e+j) = repmat(repelem(pol_grids.(state), n_repelem, 1), n_repmat, 1);
      n_repelem = n_repelem*sizes.pol.states.(state);
   end
   sm = sm';
end

function B = build_linear_splines_basis_matrix(ind, w, m)
% Constructs a sparse linear interpolation matrix for 1D bracketed data.
%
% This function builds a sparse matrix `B` of size m×n that maps `n` query points,
% located between known grid points, to the surrounding `m`-dimensional grid
% using linear interpolation weights.
%
% It is typically used in Dynare heterogeneity-specific routines after calling
% `find_bracket_linear_weight` MEX to form interpolation
% matrices over endogenous grids.
%
% INPUTS
%   ind [n×1 int32]   : Lower bracket indices for each query (i.e., index `ilow` s.t.
%                       x(ind(i)) ≤ xq(i) ≤ x(ind(i)+1))
%
%   w   [n×1 double]  : Linear interpolation weights (relative position in interval):
%                       w(i) = (x(ind(i)+1) - xq(i)) / (x(ind(i)+1) - x(ind(i)))
%
%   m   [int scalar]  : Size of the full grid (number of basis functions, i.e., length of `x`)
%
% OUTPUT
%   B   [m×n sparse]  : Sparse interpolation matrix such that:
%                         f_interp = B * f_grid
%                       where `f_grid` is a column vector of function values on the grid.
% NOTES
% - Each column of `B` contains at most two nonzero entries: `w(i)` and `1 - w(i)`
% - When w == 1, only the lower index entry is included (upper weight is zero)
% - Assumes linear interpolation over uniform or non-uniform 1D grids
    n = length(ind);
    % Identify interior points (w < 1) that need both entries
    interior = w < 1;
    % Build sparse triplets: all lower entries + upper entries for interior points
    i = [ind; ind(interior)+1];
    j = [int32(1:n)'; int32(find(interior))];
    v = [w; 1-w(interior)];
    B = sparse(i, j, v, m, n);
end

function [indices, sizes] = compute_model_indices(M_, indices, sizes)
% COMPUTE_MODEL_INDICES Compute all indices and sizes for heterogeneous-agent models
%
% This function computes comprehensive index structures for aggregate variables (Y),
% heterogeneous variables (x), and aggregation operators (Ix). It consolidates
% the index computation logic for use across multiple heterogeneity module functions.
%
% INPUTS
%   M_       [struct] : Dynare model structure
%   indices  [struct] : Partial indices structure containing:
%                       .states - State variable names (cell array)
%                       .shocks - Shock variable names (cell array)
%   sizes    [struct] : Partial sizes structure containing:
%                       .n_e - Number of shocks
%                       .n_a - Number of states
%
% OUTPUTS
%   indices  [struct] : Complete indices structure containing:
%       .Y.ind_in_dG.all      [vector] : All aggregate variable indices in dG
%       .Y.ind_in_dG.exo      [vector] : Exogenous aggregate variable indices
%       .Y.ind_in_dG.lead     [vector] : Lead aggregate variable indices
%       .Y.ind_in_dG.current  [vector] : Current aggregate variable indices
%       .Y.ind_in_dG.lag      [vector] : Lag aggregate variable indices
%       .Y.names.exo          [cell]   : Exogenous variable names
%       .Y.names.lead         [cell]   : Lead variable names
%       .Y.names.current      [cell]   : Current variable names
%       .Y.names.lag          [cell]   : Lag variable names
%       .Y.ind_in_dF          [vector] : Aggregate variable indices in dF
%       .het_states           [vector] : State variable indices
%       .Ix.in_agg            [vector] : Aggregated variable indices in M_.endo_names
%       .Ix.in_het            [vector] : Aggregated variable indices in H_.endo_names
%   sizes    [struct] : Size structure containing:
%       .N_Y.all              [scalar] : Total number of aggregate variables
%       .N_Y.exo              [scalar] : Number of exogenous aggregate variables
%       .N_Y.lead             [scalar] : Number of lead aggregate variables
%       .N_Y.current          [scalar] : Number of current aggregate variables
%       .N_Y.lag              [scalar] : Number of lag aggregate variables
%       .N_Ix                 [scalar] : Number of aggregation operators
    H_ = M_.heterogeneity(1);

    % ========================================================================
    % PART 1: Aggregate Variables (Y)
    % ========================================================================
    % Identify which aggregate variables appear in the heterogeneous block's
    % Jacobian (dF). These are the aggregate variables that affect household
    % decisions.

    n = 3*H_.endo_nbr + H_.exo_nbr;  % Offset for aggregate variables in dF

    % Find columns in dF that correspond to aggregate variables
    ind_Y_in_dF = ismember(H_.dynamic_g1_sparse_colval, n+(1:(3*M_.endo_nbr+M_.exo_nbr)));
    ind_Y_in_dF = H_.dynamic_g1_sparse_colval(ind_Y_in_dF);
    ind_Y_in_dF = unique(ind_Y_in_dF);
    ind_Y_in_dF = double(ind_Y_in_dF);

    % Convert to indices in dG (aggregate block's Jacobian)
    ind_Y_in_dG_all = ind_Y_in_dF - n;

    % Classify aggregate variables by timing (lag/current/lead/exo)
    % Layout in dG: [Y(-1), Y(0), Y(+1), exo]
    %               [1:M_.endo_nbr, M_.endo_nbr+1:2*M_.endo_nbr,
    %                2*M_.endo_nbr+1:3*M_.endo_nbr, 3*M_.endo_nbr+1:end]

    indices.Y.ind_in_dG.all = ind_Y_in_dG_all;
    indices.Y.ind_in_dG.exo = ind_Y_in_dG_all(ind_Y_in_dG_all > 3*M_.endo_nbr) - 3*M_.endo_nbr;
    indices.Y.ind_in_dG.lead = ind_Y_in_dG_all((ind_Y_in_dG_all > 2*M_.endo_nbr) & ...
                                                (ind_Y_in_dG_all <= 3*M_.endo_nbr)) - 2*M_.endo_nbr;
    indices.Y.ind_in_dG.current = ind_Y_in_dG_all((ind_Y_in_dG_all > M_.endo_nbr) & ...
                                                   (ind_Y_in_dG_all <= 2*M_.endo_nbr)) - M_.endo_nbr;
    indices.Y.ind_in_dG.lag = ind_Y_in_dG_all(ind_Y_in_dG_all <= M_.endo_nbr);

    % Get variable names for each timing
    indices.Y.names.exo = M_.exo_names(indices.Y.ind_in_dG.exo);
    if M_.orig_endo_nbr > 0
        indices.Y.names.lead = M_.endo_names(indices.Y.ind_in_dG.lead);
        indices.Y.names.current = M_.endo_names(indices.Y.ind_in_dG.current);
        indices.Y.names.lag = M_.endo_names(indices.Y.ind_in_dG.lag);
    end

    % Store original indices in dF
    indices.Y.ind_in_dF = ind_Y_in_dF;

    % Count variables by category
    sizes.N_Y.all = numel(ind_Y_in_dF);
    sizes.N_Y.exo = numel(indices.Y.names.exo);
    if M_.orig_endo_nbr
        sizes.N_Y.lead = numel(indices.Y.names.lead);
        sizes.N_Y.current = numel(indices.Y.names.current);
        sizes.N_Y.lag = numel(indices.Y.names.lag);
    end

    % ========================================================================
    % PART 2: Aggregation Operators (Ix)
    % ========================================================================
    % Identify which heterogeneous variables are aggregated using SUM(...)
    % operators and appear in the aggregate model block.
    % Extract aggregation operators from aux_vars
    N_Ix = size(M_.heterogeneity_aggregates, 1);
    iIx = zeros(N_Ix, 1);
    iIx_in_endo = zeros(N_Ix, 1);
    N_aux = numel(M_.aux_vars);
    j = 0;

    for i = 1:N_aux
        if M_.aux_vars(i).type == 14  % Heterogeneity aggregation auxiliary variable
            j = j + 1;
            if strcmp(M_.heterogeneity_aggregates{j, 1}, 'sum')
                iIx(j) = M_.aux_vars(i).endo_index;  % Index in M_.endo_names
                iIx_in_endo(j) = M_.heterogeneity_aggregates{j, 3};  % Index in H_.endo_names
            end
        end
    end

    % Store aggregation indices
    indices.Ix.in_agg = int32(iIx);  % Indices in M_.endo_names (aggregate auxiliary variables)
    indices.Ix.in_het = int32(iIx_in_endo);  % Indices in H_.endo_names (heterogeneous variables)
    sizes.N_Ix = N_Ix;

    % ========================================================================
    % PART 3: Heterogeneous MCP multipliers
    % ========================================================================
    % Extract multiplier indices for Fischer-Burmeister solver
    % Multiplier variables (MULT_L_*, MULT_U_*) are created by the preprocessor
    % when unfolding complementarity conditions. Type 15 = heterogeneousMultiplier
    mult_in_het = [];
    if isfield(H_, 'aux_vars') && ~isempty(H_.aux_vars)
        for i = 1:length(H_.aux_vars)
            if H_.aux_vars(i).type == 15  % Type 15 = heterogeneousMultiplier
                mult_in_het(end+1) = H_.aux_vars(i).endo_index;
            end
        end
    end
    indices.mult.in_het = int32(mult_in_het);

end

function [x_bar, x_bar_dash] = compute_pol_matrices(pol_values, N_x, N_sp, U, L, P, Q, x_names)
% Construct interpolated policy function matrices for state transitions.
%
% Given discretized policy functions, this function builds:
% - `x_bar`: matrix of policy function values (flattened over the joint state space)
% - `x_bar_dash`: transformed version of `x_bar` projected onto the basis used for expectations
%
% INPUTS
%   pol_values  [struct]   : Structure containing policy function values for each variable
%   N_x         [integer]  : Number of policy variables
%   N_sp        [integer]  : Number of grid points in the joint state space
%   U           [matrix]   : Upper triangular matrix from basis decomposition
%   L           [matrix]   : Lower triangular matrix from basis decomposition
%   P           [vector]   : Row permutation vector from basis decomposition
%   Q           [vector]   : Column permutation vector from basis decomposition
%   x_names     [cell]     : Cell array of policy variable names
%
% OUTPUTS
%   x_bar       [matrix]   : Matrix of policy function values (N_x × N_sp)
%   x_bar_dash  [matrix]   : Projected matrix in basis representation for expectations (N_x × N_sp)
   % Computing x_bar
   x_bar = NaN(N_x, N_sp);
   for i=1:N_x
      x = x_names{i};
      if isfield(pol_values, x)
         x_bar(i,:) = reshape(pol_values.(x),1,[]);
      end
   end
   % Computing x_bar^#
   x_bar_dash = zeros(N_x, N_sp);
   x_bar_dash(:,P) = (x_bar(:,Q) / U) / L;
end
