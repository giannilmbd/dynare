function [out_ss, sizes, indices] = check_steady_state_input(M_, options_het, ss)
% Validate and process steady-state input for heterogeneous-agent models
% Supports AR(1) discretization only. I.i.d shocks are handled as AR(1) with zero persistence.
%
% INPUTS
% - M_          [structure] MATLAB's structure describing the model
% - options_het [structure] MATLAB's structure describing the current options
% - ss          [structure] MATLAB's structure with the model steady-state
%                       information. It contains the following fields:
%    - pol [structure] MATLAB's structure containing the policy functions
%                      discretization
%       - pol.grids [structure]: MATLAB's structure containing the nodes of the
%                                state grids as column vectors.
%       - pol.values [structure]: MATLAB's structure containing the policy
%                                 function as matrices. Row indices and column
%                                 indices follow the lexicographic order
%                                 specified in pol.shocks and pol.states
%       - pol.order [array]: a list containing the order of
%                             variables used for the dimensions of pol.values. If not
%                             specified, it follows the declaration order in the
%                             var(heterogeneity=) statement (resp.
%                             var(heterogeneity=) and varexo(heterogeneity=)
%                             statements) for discretized AR(1) processes (resp.
%                             discretized gaussian i.i.d innovations).
%    - shocks [structure]: Matlab's structure describing the discretization of
%                          individual shocks or innovations:
%       - shocks.grids [structure]: MATLAB's structure containing the nodes of
%                                   the shock grids
%       - shocks.Pi [structure]: Matlab's structure containing the Markov
%                                matrices for AR(1) processes. For i.i.d shocks
%                                (in H_.exo_names), transition matrices are
%                                automatically generated using Rouwenhorst
%                                discretization with zero persistence.
%    - d [structure]: Matlab's structure describing the steady-state
%                     distribution
%       - d.grids [structure]: structure containing the states grids as column
%                              vectors. If one of the states grid is not
%                              specified, it is assumed to be identical to
%                              its policy-function value stored in pol.grids
%       - d.hist [array]: the histogram of the distribution as an array.
%                         Dimension indices follow the lexicographic order
%                         specified in d.order.
%       - d.order [array]: a list containing the order of variables used for the
%                          dimensions of d.hist. If it is not specified, it
%                          falls back to the value induced by pol.order.
%    - agg [structure]: Matlab's tructure containing the steady-state values of
%                       aggregate variables
% OUTPUTS
% - out_ss [structure]: validated and normalized steady-state input structure.
%                       It has a similar structure as the ss structure in inputs.
% - sizes [structure]: structure containing sizes information
%    - n_e [integer]: number of shocks
%    - N_e [integer]: total number of nodes in the shocks grids
%    - shocks [structure]: number of nodes for each shock grid
%    - n_a [integer]: number of states
%    - pol [structure]:
%       - pol.N_a [integer]: total number of nodes in the policy state grids
%       - pol.states [structure]: number of nodes for each state grid
%    - d [structure]:
%       - d.states [structure]: number of nodes for each distribution state grid
%       - d.N_a [integer]: total number of nodes in the distribution grids
%    - n_pol [integer]: number of policy variables
%    - agg [integer]: number of aggregate variables

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
   %% Checks
   % Retrieve variables information from M_
   H_ = M_.heterogeneity(1);
   state_symbs = H_.endo_names(H_.state_var);
   inn_symbs = H_.exo_names;
   agg_symbs = M_.endo_names(1:M_.orig_endo_nbr);
   required_pol_symbs = H_.endo_names(1:H_.orig_endo_nbr);
   mult_symbs = cell(0,1);
   for i = 1:numel(H_.aux_vars)
      if (H_.aux_vars(i).type == 15)
         mult_symbs{end+1,1} = H_.endo_names{H_.aux_vars(i).endo_index};
      end
   end
   % Initialize output variables
   sizes = struct;
   out_ss = struct;
   if ~isstruct(ss)
      error('Misspecified steady-state input `ss`: the steady-state input `ss` is not a structure.');
   end
   %% Determine variables used in policy functions and classify shocks
   % First, check that `ss.pol` exists to determine what variables are used
   heterogeneity.internal.check_isfield('pol', ss, 'ss.pol');
   pol = ss.pol;
   heterogeneity.internal.check_isstruct(ss.pol, 'ss.pol');

   % The ordering of `ss.pol.values` dimensions is the fundamental element
   heterogeneity.internal.check_isfield('order', ss.pol, 'ss.pol.order');
   pol_order = check_permutation(pol, 'order', 'ss.pol', ss.pol.order);

   % Filter variables in pol_order by their origin
   % pol_vars_from_endo = intersect(pol_order, required_pol_symbs); % From H_.endo_names - COMMENTED: only used by Method 1
   pol_vars_from_exo = intersect(pol_order, inn_symbs);           % From H_.exo_names

   %% Method 1: REMOVED - AR(1) state variables in H_.endo_names with proper equations
   %
   % REASON FOR REMOVAL:
   % Method 1 is not currently used in practice. Users typically specify all shocks
   % via varexo(heterogeneity=) (Method 2) rather than via AR(1) state variables in
   % var(heterogeneity=). This code is kept here commented for reference and potential
   % internal testing purposes only.
   %
   % Original Method 1 detection code (DISABLED):
   method1_ar1_shocks = {};
   ind_ar1 = [];
   % if ~isempty(pol_vars_from_endo)
   %    % Build Jacobian sparsity pattern to detect AR(1) equations
   %    J = sparse(H_.dynamic_g1_sparse_rowval, H_.dynamic_g1_sparse_colval, ones(numel(H_.dynamic_g1_sparse_colval),1), ...
   %               H_.endo_nbr, 3*H_.endo_nbr+H_.exo_nbr+3*M_.endo_nbr+M_.exo_nbr);
   %
   %    % Extract relevant parts of Jacobian
   %    J_lagged = J(:, 1:H_.endo_nbr);                              % Lagged hetero endo
   %    J_current = J(:, H_.endo_nbr+1:2*H_.endo_nbr);               % Current hetero endo
   %    J_lead = J(:, 2*H_.endo_nbr+1:3*H_.endo_nbr);               % Lead hetero endo
   %    J_exo = J(:, 3*H_.endo_nbr+1:3*H_.endo_nbr+H_.exo_nbr);     % Hetero innovations
   %
   %    % Find variables that appear in both lagged and current
   %    common_vars = J_lagged .* J_current;  % Element-wise multiply
   %
   %    % For each equation, check AR(1) pattern
   %    for eq_idx = 1:H_.endo_nbr
   %       if any(common_vars(eq_idx, :))  % Has variables in both lagged and current
   %          % Count variables in common between lagged and current
   %          common_count = sum(common_vars(eq_idx, :));
   %
   %          % Check AR(1) conditions:
   %          % 1. Exactly one variable appears in both lagged and current
   %          % 2. No lead variables
   %          % 3. Exactly one innovation
   %          % 4. All other lagged and current variables are zero
   %          if common_count == 1 && ...
   %             sum(J_lead(eq_idx, :)) == 0 && ...
   %             sum(J_exo(eq_idx, :)) == 1 && ...
   %             sum(J_lagged(eq_idx, :)) == 1 && ...
   %             sum(J_current(eq_idx, :)) == 1
   %
   %             % Find which variable this is
   %             var_idx = find(common_vars(eq_idx, :));
   %             var_name = H_.endo_names{var_idx};
   %
   %             % Check if this variable is in our policy function variables from H_.endo_names
   %             if ismember(var_name, pol_vars_from_endo)
   %                method1_ar1_shocks{end+1} = var_name;
   %                ind_ar1(end+1) = eq_idx;
   %             end
   %          end
   %       end
   %    end
   % end
   
   %% Method 2: Exogenous shocks in H_.exo_names
   method2a_ar1_shocks = {}; % Exogenous shocks with provided discretization
   % NOTE: method2b (automatic IID discretization) has been removed
   % Reason: Users provide policy function values on tensored grids, so all shock
   % discretizations must be user-provided to ensure consistency

   % Check shocks structure
   heterogeneity.internal.check_isfield('shocks', ss, 'ss.shocks')
   shocks = ss.shocks;
   heterogeneity.internal.check_isstruct(shocks, 'ss.shocks')
   heterogeneity.internal.check_isfield('Pi',shocks);
   heterogeneity.internal.check_isstruct(shocks.Pi, 'ss.shocks.Pi');
   heterogeneity.internal.check_isfield('grids', shocks);
   heterogeneity.internal.check_isstruct(shocks.grids, 'ss.shocks.grids')

   for i = 1:numel(pol_vars_from_exo)
      s = pol_vars_from_exo{i};
      % Check if discretization is provided
      has_grid = isfield(shocks.grids, s);
      has_pi = isfield(shocks.Pi, s);

      if has_grid && has_pi
         % Method 2a: Exogenous shocks with provided discretization
         method2a_ar1_shocks{end+1} = s;
      else
         % Both grids and Pi must be provided for all shocks
         error('Misspecified steady-state input `ss`: shock %s requires both `ss.shocks.grids.%s` and `ss.shocks.Pi.%s` to be provided. Automatic discretization is not supported as it would be incompatible with user-provided policy function values.', s, s, s);
      end
   end
   
   %% Generate indices structure in canonical order
   indices = struct();
   
   % Build canonical order: shocks first (H_.endo_names then H_.exo_names), then states
   % Method 1 AR(1) shocks in H_.endo_names declaration order
   method1_ordered = {};
   for i = 1:numel(H_.endo_names)
      if ismember(H_.endo_names{i}, method1_ar1_shocks)
         method1_ordered{end+1} = H_.endo_names{i};
      end
   end
   
   % Method 2a shocks in H_.exo_names declaration order
   method2_ordered = {};
   for i = 1:numel(H_.exo_names)
      if ismember(H_.exo_names{i}, method2a_ar1_shocks)
         method2_ordered{end+1} = H_.exo_names{i};
      end
   end
   
   % States in H_.endo_names declaration order (excluding Method 1 AR(1) shocks)
   states_ordered = {};
   for i = 1:numel(H_.endo_names)
      if ismember(H_.endo_names{i}, state_symbs) && ~ismember(H_.endo_names{i}, method1_ar1_shocks)
         states_ordered{end+1} = H_.endo_names{i};
      end
   end
   
   % Set indices structure
   indices.shocks.endo = method1_ordered';  % Method 1 AR(1) shocks from H_.endo_names
   indices.shocks.exo = method2_ordered';   % Method 2a/2b shocks from H_.exo_names  
   indices.shocks.all = [indices.shocks.endo; indices.shocks.exo];  % Combined list
   indices.states = states_ordered';
   indices.order = [indices.shocks.all; indices.states];
   
   % Generate indices.x.ind.eq (indices of non-AR(1) equations)
   ind_eq = 1:H_.endo_nbr;
   ind_eq = ind_eq(~ismember(ind_eq, ind_ar1));
   indices.x.ind.eq = ind_eq;
   
   %% Process shock discretizations
   % Initialize output shocks structure
   out_ss.shocks.grids = struct();
   out_ss.shocks.Pi = struct();

   % All shocks for final processing (method2b removed - all shocks must be user-provided)
   all_shocks = [method1_ar1_shocks, method2a_ar1_shocks];
   shock_symbs = all_shocks;
   %% Method 1: REMOVED - AR(1) state variables in H_.endo_names
   %
   % REASON FOR REMOVAL:
   % Method 1 is not currently used in practice. Users typically specify all shocks
   % via varexo(heterogeneity=) (Method 2) rather than via AR(1) state variables in
   % var(heterogeneity=). This code is kept here commented for reference and potential
   % internal testing purposes only.
   %
   % Original Method 1 code (DISABLED):
   % for i = 1:numel(method1_ar1_shocks)
   %    s = method1_ar1_shocks{i};
   %    % User must provide discretization for Method 1 AR(1) shocks
   %    has_grid = isfield(shocks, 'grids') && isstruct(shocks.grids) && isfield(shocks.grids, s);
   %    has_pi = isfield(shocks, 'Pi') && isstruct(shocks.Pi) && isfield(shocks.Pi, s);
   %
   %    if ~has_grid && ~has_pi
   %       error('Misspecified steady-state input `ss`: AR(1) process %s (Method 1) requires discretization in either `ss.shocks.grids.%s` or `ss.shocks.Pi.%s`.', s, s, s);
   %    end
   %
   %    if has_grid
   %       % Validate and copy grid
   %       if ~(isnumeric(shocks.grids.(s)) && isreal(shocks.grids.(s)) && isvector(shocks.grids.(s)) && ~issparse(shocks.grids.(s)))
   %          error('Misspecified steady-state input `ss`: `ss.shocks.grids.%s` is not a dense real vector.', s);
   %       end
   %       if isrow(shocks.grids.(s))
   %          out_ss.shocks.grids.(s) = shocks.grids.(s)';
   %       else
   %          out_ss.shocks.grids.(s) = shocks.grids.(s);
   %       end
   %    end
   %
   %    if has_pi
   %       % Validate and copy Pi matrix
   %       if ~(ismatrix(shocks.Pi.(s)) && isnumeric(shocks.Pi.(s)) && isreal(shocks.Pi.(s)) && ~issparse(shocks.Pi.(s)))
   %          error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` is not a dense real matrix.', s);
   %       end
   %       if has_grid && size(shocks.Pi.(s), 1) ~= numel(out_ss.shocks.grids.(s))
   %          error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` row count does not match `ss.shocks.grids.%s` size.', s, s);
   %       end
   %       if has_grid && size(shocks.Pi.(s), 2) ~= numel(out_ss.shocks.grids.(s))
   %          error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` column count does not match `ss.shocks.grids.%s` size.', s, s);
   %       end
   %       % Validate Markov matrix properties
   %       if any(any(shocks.Pi.(s) < 0))
   %          error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` contains negative elements.', s);
   %       end
   %       if any(abs(sum(shocks.Pi.(s), 2) - 1) > options_het.check.tol_check_sum)
   %          error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` has row sums different from 1.', s);
   %       end
   %       out_ss.shocks.Pi.(s) = shocks.Pi.(s);
   %    end
   % end
   
   %% Process Method 2a: Exogenous AR(1) shocks with provided discretization
   for i = 1:numel(method2a_ar1_shocks)
      s = method2a_ar1_shocks{i};
      % Process provided discretization
      has_grid = isfield(shocks, 'grids') && isstruct(shocks.grids) && isfield(shocks.grids, s);
      has_pi = isfield(shocks, 'Pi') && isstruct(shocks.Pi) && isfield(shocks.Pi, s);
      
      if has_grid
         % Validate and copy grid
         if ~(isnumeric(shocks.grids.(s)) && isreal(shocks.grids.(s)) && isvector(shocks.grids.(s)) && ~issparse(shocks.grids.(s)))
            error('Misspecified steady-state input `ss`: `ss.shocks.grids.%s` is not a dense real vector.', s);
         end
         if isrow(shocks.grids.(s))
            out_ss.shocks.grids.(s) = shocks.grids.(s)';
         else
            out_ss.shocks.grids.(s) = shocks.grids.(s);
         end
         % Check for NaN values
         temp_struct.(s) = out_ss.shocks.grids.(s);
         heterogeneity.internal.check_fun(temp_struct, 'ss.shocks.grids', {s}, @(x) ~any(isnan(x)), 'contain NaN elements');
         % Check for strict monotonicity
         heterogeneity.internal.check_fun(temp_struct, 'ss.shocks.grids', {s}, @(x) all(diff(x) > 0), 'are not strictly increasing');
      end
      
      if has_pi
         % Validate and copy Pi matrix
         if ~(ismatrix(shocks.Pi.(s)) && isnumeric(shocks.Pi.(s)) && isreal(shocks.Pi.(s)) && ~issparse(shocks.Pi.(s)))
            error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` is not a dense real matrix.', s);
         end
         if has_grid && size(shocks.Pi.(s), 1) ~= numel(out_ss.shocks.grids.(s))
            error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` row count does not match `ss.shocks.grids.%s` size.', s, s);
         end
         if has_grid && size(shocks.Pi.(s), 2) ~= numel(out_ss.shocks.grids.(s))
            error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` column count does not match `ss.shocks.grids.%s` size.', s, s);
         end
         % Validate Markov matrix properties
         if (any(any(isnan(shocks.Pi.(s)))))
            error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` contains NaN elements.', s);
         end
         if any(any(shocks.Pi.(s) < 0))
            error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` contains negative elements.', s);
         end
         if any(abs(sum(shocks.Pi.(s), 2) - 1) > options_het.check.tol_check_sum)
            error('Misspecified steady-state input `ss`: `ss.shocks.Pi.%s` has row sums different from 1.', s);
         end
         out_ss.shocks.Pi.(s) = shocks.Pi.(s);
      end
   end
   
   %% Method 2b: REMOVED - Automatic IID discretization no longer supported
   %
   % REASON FOR REMOVAL:
   % Users provide policy function values on a tensored grid of (shocks × states).
   % Automatic discretization of shocks would be incompatible with user-provided
   % policy function values, as the grids would not necessarily match.
   %
   % Therefore, all shock discretizations (both grids and Pi) must be user-provided
   % to ensure consistency with the policy function values.
   %
   % This code is kept here commented for reference and potential internal testing
   % purposes only. It should never be reached by end users.
   %
   % Original Method 2b code (DISABLED):
   % for i = 1:numel(method2b_iid_shocks)
   %    s = method2b_iid_shocks{i};
   %    N = options_het.rouwenhorst.grid_size;
   %    shock_idx = find(strcmp(inn_symbs, s));
   %    if isempty(shock_idx)
   %       error('Misspecified steady-state input `ss`: cannot find i.i.d shock %s in H_.exo_names.', s);
   %    end
   %    sigma = sqrt(H_.Sigma_e(shock_idx, shock_idx));
   %    tol = options_het.check.tol_check_sum;
   %    max_iter = options_het.rouwenhorst.max_iter;
   %    [y_grid, ~, Pi_mat] = rouwenhorst(0.0, sigma, N, tol, max_iter);
   %    out_ss.shocks.grids.(s) = y_grid;
   %    out_ss.shocks.Pi.(s) = Pi_mat;
   %    if ~options_het.check.no_warning_redundant
   %       warning('Steady-state input `ss`: i.i.d shock %s has been automatically discretized using Rouwenhorst method with zero persistence (N=%d, sigma=%.4f).', s, N, sigma);
   %    end
   % end
   
   % Store shock information
   sizes.n_e = numel(shock_symbs);
   grid_nb = cellfun(@(s) numel(out_ss.shocks.grids.(s)), shock_symbs); 
   sizes.N_e = prod(grid_nb);
   for i=1:numel(shock_symbs)
      s = shock_symbs{i};
      sizes.shocks.(s) = numel(out_ss.shocks.grids.(s));
   end
   
   % Update state and policy variables based on new classification
   final_state_symbs = indices.states;  % States excluding AR(1) shocks
   
   % Policy variables that should be present (excluding AR(1) shocks but including multipliers)
   required_pol_symbs_updated = setdiff(required_pol_symbs, method1_ar1_shocks);
   relevant_pol_symbs = [required_pol_symbs_updated; mult_symbs];
   
   %% Policy functions
   % pol structure was already checked above, now validate grids
   % Check that the `ss.pol.grids` field exists
   heterogeneity.internal.check_isfield('grids', ss.pol, 'ss.pol.grids');
   % Check that `ss.pol.grids` is a structure
   heterogeneity.internal.check_isstruct(ss.pol.grids, 'ss.pol.grids');
   % Check that the state grids specification and the var(heterogeneity=)
   % statement are compatible
   heterogeneity.internal.check_missingredundant(pol.grids, 'ss.pol.grids', final_state_symbs, options_het.check.no_warning_redundant);
   % Check that `ss.pol.grids` values are dense real vectors
   heterogeneity.internal.check_fun(pol.grids, 'ss.pol.grids', final_state_symbs, @(x) isnumeric(x) && isreal(x) && isvector(x) && ~issparse(x), 'are not dense real vectors');
   % Check for NaN values
   heterogeneity.internal.check_fun(pol.grids, 'ss.pol.grids', final_state_symbs, @(x) ~any(isnan(x)), 'contain NaN elements');
   % Check for strict monotonicity
   heterogeneity.internal.check_fun(pol.grids, 'ss.pol.grids', final_state_symbs, @(x) all(diff(x) > 0), 'are not strictly increasing');
   % Store:
   %    - the number of states
   sizes.n_a = numel(final_state_symbs);
   %    - the size of the tensor product of the states grids 
   grid_nb = cellfun(@(s) numel(pol.grids.(s)), final_state_symbs);
   sizes.pol.N_a = prod(grid_nb);
   %    - the size of each state grid
   for i=1:numel(final_state_symbs)
      s = final_state_symbs{i};
      sizes.pol.states.(s) = numel(pol.grids.(s));
   end
   % Make `ss.pol.grids` values column vector in the output steady-state
   % structure
   for i=1:numel(final_state_symbs)
      s = final_state_symbs{i};
      if isrow(ss.pol.grids.(s))
         out_ss.pol.grids.(s) = ss.pol.grids.(s)';
      else
         out_ss.pol.grids.(s) = ss.pol.grids.(s);
      end
   end
   % Check that the field `ss.pol.values` exists 
   heterogeneity.internal.check_isfield('values', pol, 'ss.pol.values');
   % Check that `ss.pol.values` is a struct
   heterogeneity.internal.check_isstruct(pol.values, 'ss.pol.values');
   % Check the missing and redundant variables in `ss.pol.values`
   heterogeneity.internal.check_missingredundant(pol.values, 'ss.pol.values', required_pol_symbs_updated, true);
   pol_values = fieldnames(pol.values);
   pol_values_in_relevant_pol_symbs = ismember(pol_values, relevant_pol_symbs);
   if ~options_het.check.no_warning_redundant
      if ~all(pol_values_in_relevant_pol_symbs)
         warning('Steady-state input `ss`. The following fields are redundant in `%s`: %s.', 'ss.pol.values', strjoin(pol_values(~pol_values_in_relevant_pol_symbs)));
      end
   end
   provided_relevant_pol_symbs = pol_values(pol_values_in_relevant_pol_symbs);
   % Check that `ss.pol.values` values are dense real matrices
   heterogeneity.internal.check_fun(pol.values, 'ss.pol.values', provided_relevant_pol_symbs, @(x) isnumeric(x) && isreal(x) && ~issparse(x), 'are not dense real arrays');
   % Check for NaN values
   heterogeneity.internal.check_fun(pol.values, 'ss.pol.values', provided_relevant_pol_symbs, @(x) ~any(isnan(x(:))), 'contain NaN elements');
   % Check that policy functions for state variables stay within grid bounds
   for i = 1:numel(final_state_symbs)
      s = final_state_symbs{i};
      pol_min = min(pol.values.(s)(:));
      pol_max = max(pol.values.(s)(:));
      grid_min = min(pol.grids.(s));
      grid_max = max(pol.grids.(s));

      if pol_min < grid_min
         error('Misspecified steady-state input `ss`: policy function for state variable `%s` has minimum value %.6f which is below the grid minimum %.6f.', s, pol_min, grid_min);
      end

      if pol_max > grid_max
         error('Misspecified steady-state input `ss`: policy function for state variable `%s` has maximum value %.6f which exceeds the grid maximum %.6f.', s, pol_max, grid_max);
      end
   end
   % Check the number of dimensions of `ss.pol.values` fields
   sizes.n_pol = H_.endo_nbr;
   heterogeneity.internal.check_fun(pol.values, 'ss.pol.values', provided_relevant_pol_symbs, @(x) ndims(x), 'have a number of dimensions that is not consistent with the number of states and shocks', sizes.n_a+sizes.n_e);
   % Validate that user's pol.order matches our canonical indices.order
   out_ss.pol.order = check_permutation(pol, 'order', 'ss.pol', indices.order);
   % Check the internal size of dimensions of `ss.pol.values` fields 
   isshock = ismember(out_ss.pol.order, indices.shocks.all);
   ind_shock = find(isshock);
   ind_state = find(~isshock);
   for i=1:numel(ind_shock)
      d = ind_shock(i);
      var = out_ss.pol.order{d};
      heterogeneity.internal.check_fun(pol.values, 'ss.pol.values', provided_relevant_pol_symbs, @(x) size(x, d), sprintf('have a dimension-%s size that is not compatible with the size ss.shocks.grids.%s.', num2str(d), var), sizes.shocks.(var));
   end
   for i=1:numel(ind_state)
      d = ind_state(i);
      var = out_ss.pol.order{d};
      heterogeneity.internal.check_fun(pol.values, 'ss.pol.values', provided_relevant_pol_symbs, @(x) size(x, d), sprintf('have a dimension-%s size that is not compatible with the size ss.pol.grids.%s.', num2str(d), var), sizes.pol.states.(var));
   end
   % Copy `ss.pol.values` in `out_ss`
   out_ss.pol.values = ss.pol.values;
   %% Distribution
   % Check that the field `ss.d` exists 
   heterogeneity.internal.check_isfield('d', ss, 'ss.d');
   d = ss.d;
   % Check that the field `ss.d.hist` exists
   heterogeneity.internal.check_isfield('hist', ss.d, 'ss.d.hist');
   % Check the type of `ss.d.hist`
   if ~(isnumeric(d.hist) && isreal(d.hist) && ~issparse(d.hist))
      error('Misspecified steady-state input `ss`: `ss.d.hist` is not a dense real array.');
   end
   % Check for NaN values
   if any(isnan(d.hist(:)))
      error('Misspecified steady-state input `ss`: `ss.d.hist` contains NaN elements.');
   end
   % Check the consistency of `ss.d.grids` 
   if ~isfield(d, 'grids')
      if ~options_het.check.no_warning_d_grids
         warning('In the steady-state input `ss.d.grids`, no distribution-specific grid is set for states %s. The policy grids in `ss.pol.grids` shall be used.' , strjoin(final_state_symbs));
      end
      % Copy the relevant sizes from the pol field
      sizes.d = sizes.pol;
      out_ss.d.grids = out_ss.pol.grids;
   else
      % Check that `ss.d.grids` is a struct
      heterogeneity.internal.check_isstruct(d.grids, 'ss.d.grids');
      % Check redundant variables in `ss.d.grids`
      d_grids_symbs = fieldnames(d.grids);
      if isempty(d_grids_symbs)
         if ~options_het.check.no_warning_redundant
            warning('In the steady-state input `ss.d.grids`, no distribution-specific grid is set for states %s. The policy grids in `ss.pol.grids` shall be used.' , strjoin(final_state_symbs));
         end
         % Copy the relevant sizes from the pol field
         sizes.d = sizes.pol;
         out_ss.d.grids = out_ss.pol.grids;
      else
         d_grids_symbs_in_states = ismember(d_grids_symbs, final_state_symbs);
         if ~all(d_grids_symbs_in_states)
            if ~options_het.check.no_warning_redundant
               warning('In the steady-state input `ss.d.states`, the following specification for the states grids in the distribution structure are not useful: %s.', strjoin(d_grids_symbs(~d_grids_symbs_in_states)));
            end
            d_grids_symbs = d_grids_symbs(d_grids_symbs_in_states);
         end
         % Check the types of `ss.d.grids` elements
         heterogeneity.internal.check_fun(pol.grids, 'ss.d.grids', d_grids_symbs, @(x) isnumeric(x) && isreal(x) && isvector(x) && ~issparse(x), 'are not dense real vectors');
         % Check for NaN values
         heterogeneity.internal.check_fun(d.grids, 'ss.d.grids', d_grids_symbs, @(x) ~any(isnan(x)), 'contain NaN elements');
         % Check for strict monotonicity
         heterogeneity.internal.check_fun(d.grids, 'ss.d.grids', d_grids_symbs, @(x) all(diff(x) > 0), 'are not strictly increasing');
         % Check that policy functions for state variables stay within d.grids bounds
         for i = 1:numel(d_grids_symbs)
            s = d_grids_symbs{i};
            pol_min = min(pol.values.(s)(:));
            pol_max = max(pol.values.(s)(:));
            grid_min = min(d.grids.(s));
            grid_max = max(d.grids.(s));

            if pol_min < grid_min
               error('Misspecified steady-state input `ss`: policy function for state variable `%s` has minimum value %.6f which is below the d.grids minimum %.6f.', s, pol_min, grid_min);
            end

            if pol_max > grid_max
               error('Misspecified steady-state input `ss`: policy function for state variable `%s` has maximum value %.6f which exceeds the d.grids maximum %.6f.', s, pol_max, grid_max);
            end
         end
         % Store the size of the states grids for the distribution
         for i=1:numel(d_grids_symbs)
            s = d_grids_symbs{i};
            sizes.d.states.(s) = numel(d.grids.(s));
            out_ss.d.grids.(s) = d.grids.(s);
         end
         states_out_of_d = setdiff(final_state_symbs, d_grids_symbs);
         if ~isempty(states_out_of_d)
            if ~options_het.check.no_warning_d_grids
               warning('heterogeneity.bbeg.het_stoch_simul: in the steady-state structure, no distribution-specific grid set for states %s. The policy grids specified in pol.grids shall be used.', strjoin(states_out_of_d));
            end
            % Store the sizes of the unspecified states grids from the pol field
            for i=1:numel(states_out_of_d)
               s = states_out_of_d{i};
               sizes.d.states.(s) = sizes.pol.states.(s);
               out_ss.d.grids.(s) = pol.grids.(s);
            end
         end
      end
   end
   % Check the permutation of variables in the distribution histogram
   out_ss.d.order = check_permutation(d, 'order', 'ss.d', out_ss.pol.order);
   % Check the internal size of dimensions of `ss.d.hist` 
   isshock = ismember(out_ss.d.order, indices.shocks.all);
   ind_shock = find(isshock);
   ind_state = find(~isshock);
   for i=1:numel(ind_shock)
      dim = ind_shock(i);
      var = out_ss.d.order{dim};
      if size(d.hist,dim) ~= sizes.shocks.(var)
         error('Misspecified steady-state input `ss`: dimension %s of `ss.d.hist` and `ss.shocks.grids.%s` have incompatible sizes.', num2str(dim), var);
      end
   end
   for i=1:numel(ind_state)
      dim = ind_state(i);
      var = out_ss.d.order{dim};
      if size(d.hist,dim) ~= sizes.d.states.(var)
         error('Misspecified steady-state input `ss`: dimension %s of `ss.d.hist` and `ss.pol.grids.%s` have incompatible sizes.', num2str(dim), var);
      end
   end
   % Copy `ss.d.hist` in `out_ss` 
   out_ss.d.hist = ss.d.hist;
   %% Aggregate variables
   % Check that the field `ss.agg` exists
   heterogeneity.internal.check_isfield('agg', ss, 'ss.agg');
   agg = ss.agg;
   % Check that `ss.agg` is a structure
   heterogeneity.internal.check_isstruct(agg, 'ss.agg');
   % Check that the aggregate variables specification and the var statement are
   % compatible
   heterogeneity.internal.check_missingredundant(agg, 'ss.agg', agg_symbs, options_het.check.no_warning_redundant);
   % Store the number of aggregate variables
   sizes.agg = numel(fieldnames(agg));
   % Check the types of `ss.agg` values
   heterogeneity.internal.check_fun(agg, 'ss.agg', agg_symbs, @(x) isreal(x) && isscalar(x), 'are not real scalars');
   % Check for NaN values
   heterogeneity.internal.check_fun(agg, 'ss.agg', agg_symbs, @(x) ~isnan(x), 'contain NaN values');
   % Copy `ss.agg` into `out_ss`
   out_ss.agg = agg;
end

function out_order = check_permutation(s, f_name, s_name, symbs)
% Check and returns a permutation order for variables.
%
% INPUTS
% - s [structure]: structure containing the permutation information
% - f_name [char]: name of the field within `s` that should contain the permutation list
% - s_name [char]: full name of the structure `s` for error messages
% - symbs [cell array of char]: list of expected variable names
%
% OUTPUTS
% - out_order [cell array or string array]: permutation order of variables
%
% DESCRIPTION
% Checks that the field `f_name` exists in structure `s` and contains a valid
% permutation (i.e., a cell array of character vectors or a string array).
% If the field is missing, returns the default order given by `symbs`.
% Throws an error if the specified variables are not consistent with `symbs`.
   if ~isfield(s, f_name)
      out_order = symbs;
   else
      f = s.(f_name);
      if ~isstring(f) && ~iscellstr(f)
         error('Misspecified steady-state input `ss`: the `%s.%s` field should be a cell array of character vectors or a string array.', s_name, f_name);
      end
      err_var = setdiff(f, symbs);
      if ~isempty(err_var)
         error('Misspecified steady-state input `ss`: the set of variables of the `%s.%s` field is not consistent with the information in M_ and the other fields in the steady-state structure `ss`. Problematic variables: %s.', s_name, f_name, strjoin(err_var));
      end
      out_order = f;
   end
end