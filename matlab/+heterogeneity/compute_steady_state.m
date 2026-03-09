function [oo_het, M_params] = compute_steady_state(M_, options_het, oo_het, steady_state)
% Compute steady state for heterogeneous-agent DSGE models
%
% This function computes the steady state of a heterogeneous-agent model by
% iterating on policy functions (backward iteration) and the distribution
% (forward iteration) until convergence. It also calibrates free parameters
% to satisfy calibration target equations.
%
% INPUTS
% - M_           [structure] Dynare model structure containing:
%                  .heterogeneity - heterogeneity dimension specification
%                  .fname - model filename for dynamic functions
%                  .param_names, .params - parameter names and values
%                  .endo_nbr, .endo_names - endogenous variable info
%                  .aux_vars - auxiliary variable definitions
% - options_het  [structure] Heterogeneity-specific options:
%                  .steady_state_file_name - path to .mat file with initial guess.
%                    If empty and steady_state is not passed as 4th argument,
%                    the variable named by .steady_state_variable_name is
%                    loaded from the base workspace via evalin.
%                  .steady_state_variable_name - variable name in file (when
%                    .steady_state_file_name is non-empty) or in the base
%                    workspace (when .steady_state_file_name is empty)
%                  .check - validation options (passed to check_steady_state_input)
%                  .forward.max_iter, .forward.tol, .forward.check_every - distribution iteration
%                  .calibration.tolf, .calibration.verbosity - calibration options
% - oo_het       [structure] Heterogeneity results structure to populate
% - steady_state [structure] (optional) User-provided initial guess structure:
%                  .pol.grids - state grids for policy functions
%                  .pol.values - initial policy function values
%                  .pol.order - dimension ordering
%                  .shocks.grids, .shocks.Pi - shock discretization
%                  .agg - aggregate variable steady-state values
%                  .free_parameters.(param).initial_guess - parameters to calibrate
%                  .free_parameters.(param).lower_bound, .upper_bound - optional bounds
%
% OUTPUTS
% - oo_het       [structure] Updated heterogeneity results containing:
%                  .steady_state - validated steady-state with converged policies
%                  .sizes - grid dimension information
%                  .mat - interpolation and transition matrices
%                  .indices - variable and equation indices

% Copyright © 2025-2026 Dynare Team
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
   assert(isscalar(M_.heterogeneity), 'Heterogeneous-agent models with more than one heterogeneity dimension are not allowed yet!');

   if nargin < 4
      file_name = options_het.steady_state_file_name;
      variable_name = options_het.steady_state_variable_name;

      if isempty(file_name)
         % No file specified — load variable from base workspace
         if ~evalin('base', ['exist(''' variable_name ''', ''var'')'])
            error('heterogeneity_compute_steady_state: variable ''%s'' not found in workspace. Provide either ''filename'' or ''variable'' option.', variable_name);
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
               error('heterogeneity_compute_steady_state: Cannot find file: %s.mat', base_path);
            end
            file_name_to_load = [base_path extension];
         else
            file_name_to_load = file_name;
         end

         if ~isfile(file_name_to_load)
            error('heterogeneity_compute_steady_state: File not found: %s', file_name_to_load);
         end

         loaded_data = load(file_name_to_load);
         if ~isfield(loaded_data, variable_name)
            error('heterogeneity_compute_steady_state: Variable ''%s'' not found in file ''%s''.', ...
                  variable_name, file_name_to_load);
         end
         steady_state = loaded_data.(variable_name);
      end
   end

   % Load initial guess and prepare matrices via load_steady_state ===
   % This computes: mat.pol.*, mat.Mu, mat.d.Phi, sizes.*, indices.*
   % Skips: mat.d.hist-dependent computations
   oo_het = heterogeneity.load_steady_state(M_, options_het, oo_het, steady_state, true);

   % Extract prepared structures
   mat = oo_het.mat;
   sizes = oo_het.sizes;
   indices = oo_het.indices;
   steady_state = oo_het.steady_state;
   H_ = M_.heterogeneity(1);

   % Extract complementarity conditions (bounds on endogenous variables)
   [mat.pol.bounds.lower_bounds, mat.pol.bounds.upper_bounds] = ...
       feval(sprintf('%s.dynamic_het1_complementarity_conditions', M_.fname), M_.params);

   % Deal with calibration parameters: initial values and bounds
   % mat is already populated by Phase 1 - add free_parameters fields
   % free_parameters is optional - if not provided, we compute residuals only (no calibration)
   if ~isfield(steady_state, 'free_parameters')
      steady_state.free_parameters = struct;
   end
   free_parameters_names = fieldnames(steady_state.free_parameters);
   M_params_in_free_parameters = ismember(M_.param_names, free_parameters_names);
   free_parameters_in_M_params = ismember(free_parameters_names, M_.param_names);
   indices.free_parameters.ind = int32(find(M_params_in_free_parameters));
   indices.free_parameters.names = M_.param_names(indices.free_parameters.ind);
   if ~all(free_parameters_in_M_params)
      error('Misspecified steady-state input `steady_state.free_parameters`. The following parameters are missing in `M_.param_names`: %s.', strjoin(indices.free_parameters.names(~free_parameters_in_M_params)));
   end
   n_free_parameters = numel(free_parameters_names);
   mat.free_parameters.bounds = [-Inf(1, n_free_parameters); +Inf(1, n_free_parameters)];
   mat.free_parameters.initial_values = NaN(n_free_parameters, 1);
   for i = 1:n_free_parameters
      param = indices.free_parameters.names{i};
      check_isfield('initial_guess', steady_state.free_parameters.(param), sprintf('steady_state.free_parameters.%s.initial_guess', param));
      mat.free_parameters.initial_values(i) = steady_state.free_parameters.(param).initial_guess;
      if isfield(steady_state.free_parameters.(param), 'lower_bound')
         mat.free_parameters.bounds(1,i) = steady_state.free_parameters.(param).lower_bound;
      end
      if isfield(steady_state.free_parameters.(param), 'upper_bound')
         mat.free_parameters.bounds(2,i) = steady_state.free_parameters.(param).upper_bound;
      end
      if isfinite(mat.free_parameters.bounds(1,i)) && isfinite(mat.free_parameters.bounds(2,i)) && mat.free_parameters.bounds(1,i) > mat.free_parameters.bounds(2,i)
         error('Misspecified steady-state input: incompatible values for `steady_state.free_parameters.%s.lower_bound` and `steady_state.free_parameters.%s.upper_bound`.', param, param);
      end
   end

   %% Determine calibration target equations
   % Either user-specified via options or auto-detected from SUM operators
   user_target_eqs = options_het.calibration.target_equations;
   if ~isempty(user_target_eqs)
      % User-specified equations
      if ~iscell(user_target_eqs)
         user_target_eqs = num2cell(user_target_eqs);  % Convert numeric array to cell
      end

      indices.target_equations = zeros(1, numel(user_target_eqs), 'int32');
      for i = 1:numel(user_target_eqs)
         eq = user_target_eqs{i};
         if isnumeric(eq)
            % Equation number (1-based index)
            if eq < 1 || eq > M_.orig_endo_nbr
               error('Invalid calibration target equation index %d. Must be between 1 and %d.', eq, M_.orig_endo_nbr);
            end
            indices.target_equations(i) = int32(eq);
         elseif ischar(eq) || isstring(eq)
            % Equation name - find matching index in M_.equations_tags
            eq_idx = get_equation_number_by_tag(eq, M_);
            if eq_idx == 0
               error('Equation name ''%s'' not found in M_.equations_tags.', eq);
            end
            if eq_idx > M_.orig_endo_nbr
               error('Equation ''%s'' (index %d) is not an aggregate equation (must be <= %d).', eq, eq_idx, M_.orig_endo_nbr);
            end
            indices.target_equations(i) = int32(eq_idx);
         else
            error('Calibration target equation must be a number or string, got %s.', class(eq));
         end
      end
      indices.target_equations = int32(unique(indices.target_equations));  % Remove duplicates, sort
      target_eq_source = 'user-specified';
   else
      % Auto-detect equations with SUM operators (backward compatible default)
      % Step 1: Find indices of SUM auxiliary variables (type 14)
      sum_var_indices = [];
      if isfield(M_, 'aux_vars') && ~isempty(M_.aux_vars)
         for i = 1:length(M_.aux_vars)
            if M_.aux_vars(i).type == 14  % Type 14 = SUM aggregation
               sum_var_indices(end+1) = M_.aux_vars(i).endo_index;
            end
         end
      end

      % Step 2: Find which aggregate equations reference these SUM variables at time t
      if ~isempty(sum_var_indices)
         % Use pre-computed sparsity pattern from preprocessor
         % Dynamic Jacobian columns: [y(-1), y(0), y(+1), x] where y(0) is at time t
         y_t_cols = M_.endo_nbr + sum_var_indices;  % Columns for SUM variables at time t

         % Find rows (equations) that have non-zero entries in these columns
         mask = ismember(M_.dynamic_g1_sparse_colval, y_t_cols);
         rows = M_.dynamic_g1_sparse_rowval(mask);
         rows = rows(rows <= M_.orig_endo_nbr);  % Only aggregate equations
         indices.target_equations = int32(unique(rows));
      else
         indices.target_equations = int32([]);
      end
      target_eq_source = 'auto-detected (SUM operators)';
   end

   % Validate: number of target equations must match number of free parameters
   if n_free_parameters > 0
      if numel(indices.target_equations) ~= n_free_parameters
         error(['Number of calibration target equations (%d) must match number of unknown ' ...
                'parameters (%d).\nTarget equations: [%s]\nUnknown parameters: [%s]'], ...
                numel(indices.target_equations), n_free_parameters, ...
                strjoin(arrayfun(@num2str, indices.target_equations, 'UniformOutput', false), ', '), ...
                strjoin(indices.free_parameters.names, ', '));
      end
   end

   % Warn if calibration tolerance is tighter than time iteration tolerance
   if n_free_parameters > 0 && options_het.calibration.tolf < options_het.time_iteration.tol
      warning('heterogeneity:tolerance_mismatch', ...
              ['Calibration tolerance (calibration_tolf=%.2e) is smaller than time iteration ' ...
               'tolerance (time_iteration_tol=%.2e).\nThe calibration step may attempt to match ' ...
               'noise in the inner solver. Consider increasing calibration_tolf or decreasing ' ...
               'time_iteration_tol.'], ...
              options_het.calibration.tolf, options_het.time_iteration.tol);
   end

   % Display selected equations when verbosity == 1
   if options_het.calibration.verbosity == 1 && n_free_parameters > 0
      fprintf('\n=== Calibration Target Equations ===\n');
      for i = 1:numel(indices.target_equations)
         eq_idx = indices.target_equations(i);
         eq_name = get_equation_name_by_number(eq_idx, M_);
         if ~isempty(eq_name)
            eq_name = sprintf(' (%s)', eq_name);
         end
         fprintf('  Equation %3d%s\n', eq_idx, eq_name);
      end
      fprintf('Source: %s\n\n', target_eq_source);
   end

   %% Build equation names array for MEX (name if available, otherwise index as string)
   equation_names = cell(M_.orig_endo_nbr, 1);
   for i = 1:M_.orig_endo_nbr
      name = get_equation_name_by_number(i, M_);
      if isempty(name)
         equation_names{i} = sprintf('%d', i);
      else
         equation_names{i} = name;
      end
   end

   %% Call MEX function
   if M_.orig_endo_nbr == 0
      M_.dynamic_tmp_nbr = [];
   end
   % Tensor product grid - call compute_steady_state_tensor
   output = compute_steady_state_tensor(M_.fname, equation_names, M_.params, H_.orig_endo_nbr, H_.set_auxiliary_variables, H_.het_aux_levels, int32(H_.dynamic_mcp_equations_reordering), int32(H_.state_var), H_.dynamic_g1_sparse_rowval, H_.dynamic_g1_sparse_colval, H_.dynamic_g1_sparse_colptr, options_het, mat, indices);

   %% Process MEX output and update steady_state structure
    % 1. Check convergence
    if ~output.converged
        error(['Steady-state computation failed to converge after %d iterations.\n' ...
               'Final residual norm: %.6e\n' ...
               'Time iteration: %s (iter=%d, resid=%.3e)\n' ...
               'Distribution: %s (iter=%d, resid=%.3e)'], ...
               output.iterations, output.residual_norm, ...
               mat2str(output.time_iteration.converged), output.time_iteration.iterations, ...
               output.time_iteration.residual_norm, ...
               mat2str(output.distribution.converged), output.distribution.iterations, ...
               output.distribution.residual_norm);
    end

    % 2. Update M_.params with calibrated parameter values
    if n_free_parameters > 0 && ~isempty(output.params)
        M_.params(indices.free_parameters.ind) = output.params;

        % Print calibrated parameters
        if options_het.calibration.verbosity >= 1
            fprintf('\n=== Calibrated Parameters ===\n');
            for i = 1:numel(output.param_names)
                fprintf('  %-20s = %12.6f\n', output.param_names{i}, output.params(i));
            end
            fprintf('\n');
        end
    end

    % 3. Reshape policies from [n_het_endo × N_sp] to multi-dimensional arrays
    % Grid sizes: [N_e1, N_e1, ... , N_a1, N_a2, ...]
    pol_grid_sizes = [];
    for i = 1:sizes.n_e
        shock = indices.shocks{i};
        pol_grid_sizes = [pol_grid_sizes, sizes.shocks.(shock)];
    end
    d_grid_sizes = pol_grid_sizes;
    for i = 1:sizes.n_a
        state = indices.states{i};
        pol_grid_sizes = [pol_grid_sizes, sizes.pol.states.(state)];
        d_grid_sizes = [d_grid_sizes, sizes.d.states.(state)];
    end

    for i = 1:H_.orig_endo_nbr
        var_name = H_.endo_names{i};
        % Extract row i from policies matrix and reshape
        steady_state.pol.values.(var_name) = reshape(output.time_iteration.policies(i,:), pol_grid_sizes);
    end

    % 4. Reshape distribution from [N_sp × 1] to multi-dimensional array
    steady_state.d.hist = reshape(output.distribution.hist, d_grid_sizes);

    % 5. Store aggregate values (Ix) - these are integrated heterogeneous variables
    % Note: mapping from Ix to aggregate variables depends on model structure
    % For now, store in output structure for diagnostic purposes
    output.aggregates_named = struct();
    for i = 1:numel(indices.Ix.in_het)
        het_var_idx = indices.Ix.in_het(i);
        var_name = H_.endo_names{het_var_idx};
        output.aggregates_named.(var_name) = output.aggregates.Ix(het_var_idx);
    end

    % 6. Print convergence diagnostics
    if options_het.calibration.verbosity == 1
        if n_free_parameters > 0
            % Full calibration was performed
            fprintf('\n=== Steady-State Convergence ===\n');
            fprintf('Overall: %s after %d calibration iterations\n', ...
                    mat2str(output.converged), output.iterations);
            fprintf('Target equations residual norm: %.6e\n\n', output.residual_norm);
        else
            % No calibration - just residual computation
            fprintf('\n=== Heterogeneous Block Computation ===\n');
            fprintf('Aggregate residual norm: %.6e\n\n', output.residual_norm);
        end

        fprintf('Time Iteration (Backward):\n');
        fprintf('  Converged: %s\n', mat2str(output.time_iteration.converged));
        fprintf('  Iterations: %d\n', output.time_iteration.iterations);
        fprintf('  Residual norm (sup): %.6e\n\n', output.time_iteration.residual_norm);

        fprintf('Distribution (Forward):\n');
        fprintf('  Converged: %s\n', mat2str(output.distribution.converged));
        fprintf('  Iterations: %d\n', output.distribution.iterations);
        fprintf('  Residual norm (L1): %.6e\n\n', output.distribution.residual_norm);
    end

    % 7. Check calibration target equation residuals against tolerance (only when calibrating)
    if n_free_parameters > 0
        tol = options_het.calibration.tolf;
        if output.residual_norm > tol
            warning(['Target equations residual (%.6e) exceeds tolerance (%.6e).\n' ...
                     'Steady state may not be accurate.'], ...
                     output.residual_norm, tol);
        end
    end

    % 8. Print detailed residuals if verbosity == 1
    if options_het.calibration.verbosity == 1
        fprintf('=== Aggregate Equation Residuals ===\n');
        if n_free_parameters > 0
            % Only show target equations when calibrating
            for i = 1:numel(indices.target_equations)
                eq_idx = indices.target_equations(i);
                eq_name = get_equation_name_by_number(eq_idx, M_);
                if ~isempty(eq_name)
                    eq_name = sprintf(' (%s)', eq_name);
                end
                fprintf('  Equation %3d%s: %12.6e\n', eq_idx, eq_name, ...
                        output.aggregates.residuals(eq_idx));
            end
        else
            % Show all residuals when not calibrating
            for i = 1:numel(output.aggregates.residuals)
                eq_name = get_equation_name_by_number(i, M_);
                if ~isempty(eq_name)
                    eq_name = sprintf(' (%s)', eq_name);
                end
                fprintf('  Equation %3d%s: %12.6e\n', i, eq_name, ...
                        output.aggregates.residuals(i));
            end
        end
        fprintf('\n');
    end

    % Finalize via load_steady_state with computed steady_state ===
    % This computes: mat.d.x_bar, mat.d.ind/w/inv_h, mat.d.hist, Ix, mat.y, mat.G/dG, mat.F/dF
    oo_het = heterogeneity.load_steady_state(M_, options_het, oo_het, steady_state, false);
    M_params = M_.params;

end

function check_isfield(f, s, f_name, details)
% Check the existence of a field in a structure.
%
% INPUTS
% - f [char]: name of the field to check
% - s [structure]: structure in which the field should be present
% - f_name [optional, char]: name to display in the error message (defaults to `f`)
% - details [optional, char]: additional details to append to the error message
%
% OUTPUTS
% - (none) This function throws an error if the specified field is missing.
%
% DESCRIPTION
% Checks whether the field `f` exists in the structure `s`.
% If not, throws an error indicating the missing field.
% If provided, `details` is appended to the error message for additional context.
   if nargin < 4
      details = '';
   end
   if nargin < 3
      f_name = f;
   end
   if ~isfield(s, f)
      error('Misspecified steady-state input `steady_state`: the `%s` field is missing.%s', f_name, details);
   end
end
