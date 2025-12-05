function oo_het = load_steady_state(M_, options_het, oo_het, ss)
% Wrapper to validate, complete, and store steady-state input for models
% featuring rich heterogeneity.
%
% This function performs validation of the user-provided steady-state structure `ss`, constructs
% interpolation objects and basis matrices, computes missing policy functions for complementarity
% multipliers (if not supplied), and populates the global `oo_het` structure with the results.
%
% INPUTS
%   M_       [struct] : Dynare model structure
%   options_het [struct] : Heterogeneity-specific options structure
%   oo_het      [struct] : Heterogeneity-specific output structure to which results will be written
%   ss (optional) [struct] : User-provided steady-state structure
%
% OUTPUTS
%   oo_het      [struct] : Updated output structure containing:
%                        - oo_het.ss    : validated and completed steady-state structure
%                        - oo_het.sizes : grid size structure
%                        - oo_het.mat   : interpolation and policy function matrices

% Copyright © 2025 Dynare Team
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

      % Check that filename was provided
      if isempty(file_name)
         error('heterogeneity_load_steady_state: filename option is required');
      end

      % Parse file path to separate directory, basename, and extension
      [filepath, basename, extension] = fileparts(file_name);

      % Auto-detect extension if not provided
      if isempty(extension)
         % Reconstruct potential file names
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
         % Extension provided - use as-is
         file_name_to_load = file_name;
      end

      % Verify file exists
      if ~isfile(file_name_to_load)
         error('heterogeneity_load_steady_state: File not found: %s', file_name_to_load);
      end

      % Load the file
      loaded_data = load(file_name_to_load);
      if ~isfield(loaded_data, variable_name)
         error('heterogeneity_load_steady_state: Variable ''%s'' not found in file ''%s''.', ...
               variable_name, file_name_to_load);
      end
      ss = loaded_data.(variable_name);
   end

   % Check steady-state input
   [out_ss, sizes, indices] = heterogeneity.check_steady_state_input(M_, options_het, ss);
   % Compute useful indices and names lists
   H_ = M_.heterogeneity(1);
   % --- Aggregate variables ---
   n = 3*H_.endo_nbr+H_.exo_nbr;
   ind_Y_in_dF = ismember(H_.dynamic_g1_sparse_colval, n+(1:(3*M_.endo_nbr+M_.exo_nbr)));
   ind_Y_in_dF = H_.dynamic_g1_sparse_colval(ind_Y_in_dF);
   ind_Y_in_dF = unique(ind_Y_in_dF);
   ind_Y_in_dF = double(ind_Y_in_dF);
   ind_Y_in_dG_all = ind_Y_in_dF-n;
   indices.Y.ind_in_dG.all = ind_Y_in_dG_all;
   indices.Y.ind_in_dG.exo = ind_Y_in_dG_all(ind_Y_in_dG_all > 3*M_.endo_nbr) - 3*M_.endo_nbr;
   indices.Y.ind_in_dG.lead = ind_Y_in_dG_all((ind_Y_in_dG_all > 2*M_.endo_nbr) & (ind_Y_in_dG_all <= 3*M_.endo_nbr)) - 2*M_.endo_nbr;
   indices.Y.ind_in_dG.current = ind_Y_in_dG_all((ind_Y_in_dG_all > M_.endo_nbr) & (ind_Y_in_dG_all <= 2*M_.endo_nbr)) - M_.endo_nbr;
   indices.Y.ind_in_dG.lag = ind_Y_in_dG_all(ind_Y_in_dG_all <= M_.endo_nbr);
   indices.Y.names.exo = M_.exo_names(indices.Y.ind_in_dG.exo);
   indices.Y.names.lead = M_.endo_names(indices.Y.ind_in_dG.lead);
   indices.Y.names.current = M_.endo_names(indices.Y.ind_in_dG.current);
   indices.Y.names.lag = M_.endo_names(indices.Y.ind_in_dG.lag);
   indices.Y.ind_in_dF = ind_Y_in_dF;
   sizes.N_Y.all = numel(ind_Y_in_dF);
   sizes.N_Y.exo = numel(indices.Y.names.exo);
   sizes.N_Y.lead = numel(indices.Y.names.lead);
   sizes.N_Y.current = numel(indices.Y.names.current);
   sizes.N_Y.lag = numel(indices.Y.names.lag);
   % --- Heterogeneous variables ---
   x_names = H_.endo_names;
   ind_x_in_endo = 1:H_.endo_nbr;
   mask = ismember(x_names, indices.shocks.endo);
   ind_x_in_endo = ind_x_in_endo(~mask);
   ind_a = cellfun(@(x) find(strcmp(x,x_names)), indices.states);
   x_names(mask) = [];
   indices.x.names = x_names;
   indices.x.ind.states = ind_a;
   indices.x.ind.in_endo = ind_x_in_endo;
   sizes.N_x = numel(x_names);
   % --- Aggregated heterogeneous variables ---
   ind_endo_in_x = zeros(H_.endo_nbr,1);
   for i=1:H_.endo_nbr
      ind = find(strcmp(H_.endo_names{i}, x_names), 1, "first");
      if ~isempty(ind)
         ind_endo_in_x(i) = ind;
      end
   end
   N_Ix = size(M_.heterogeneity_aggregates,1);
   sizes.N_Ix = N_Ix;
   iIx = zeros(N_Ix,1);
   iIx_in_endo = zeros(N_Ix,1);
   N_aux = numel(M_.aux_vars);
   j = 0;
   for i=1:N_aux
      if M_.aux_vars(i).type == 14
         j = j + 1;
         if (M_.heterogeneity_aggregates{j,1} == 'sum')
            iIx(j) = M_.aux_vars(i).endo_index;
            iIx_in_endo(j) = M_.heterogeneity_aggregates{j,3};
         end
      end
   end
   indices.Ix.in_endo = iIx;
   indices.Ix.in_x = ind_endo_in_x(iIx_in_endo);
   % Compute Φ and ̃Φ with useful basis matrices
   mat = struct;
   Phi = speye(1);
   N_sp = sizes.N_e*sizes.pol.N_a;
   sizes.N_sp = N_sp;
   Phi_tilde = speye(N_sp);
   for i=sizes.n_a:-1:1
      state = indices.states{i};
      grid = out_ss.pol.grids.(state);
      [ind, w] = find_bracket_linear_weight(grid, out_ss.d.grids.(state));
      Phi = kron(Phi, heterogeneity.internal.build_linear_splines_basis_matrix(ind, w, sizes.pol.states.(state)));
      [ind, mat.pol.w.(state)] = find_bracket_linear_weight(grid, out_ss.pol.values.(state)(:));
      mat.pol.ind.(state) = ind;
      mat.pol.inv_h.(state) = 1 ./ (grid(ind+1) - grid(ind));
   end
   Phi = kron(Phi, speye(sizes.N_e));
   mat.d.Phi = Phi;
   mat.pol.Phi = Phi_tilde;
   % Compute Mu
   Mu = eye(1);
   if isfield(out_ss.shocks, 'Pi')
      for i=sizes.n_e:-1:1
         e = indices.shocks.all{i};
         Mu = kron(Mu, out_ss.shocks.Pi.(e));
      end
   end
   mat.Mu = Mu;
   % Compute state and shock grids for the policy state grids
   mat.pol.sm = heterogeneity.internal.set_state_matrix(out_ss.pol.grids, out_ss.shocks.grids, sizes, indices);
   % LU decomposition of Phi_tilde
   [mat.pol.L, mat.pol.U, mat.pol.P] = lu(Phi_tilde);
   % Computation of Phi_tilde_e
   mat.pol.Phi_e = heterogeneity.internal.build_Phi_tilde_e(sizes, indices, mat);
   % Compute policy matrices without the complementarity conditions multipliers
   [mat.pol.x_bar, mat.pol.x_bar_dash] = heterogeneity.internal.compute_pol_matrices(out_ss, sizes, mat, indices);
   % Compute the policy values of complementarity conditions multipliers
   mult_values = heterogeneity.internal.compute_pol_mcp_mul(M_, out_ss, sizes, mat, indices);
   % Update ss.pol.values accordingly
   f = fieldnames(mult_values);
   for i=1:numel(f)
     out_ss.pol.values.(f{i}) = mult_values.(f{i});
   end
   % Update the policy matrices
   [mat.pol.x_bar, mat.pol.x_bar_dash] = heterogeneity.internal.compute_pol_matrices(out_ss, sizes, mat, indices);
   % Generate relevant distribution-related objects
   mat.d.x_bar = mat.pol.x_bar_dash*mat.d.Phi;
   for i=1:sizes.n_a
      state = indices.states{i};
      grid = out_ss.d.grids.(state);
      [ind, w] = find_bracket_linear_weight(grid, mat.d.x_bar(ind_a(i),:));
      mat.d.ind.(state) = reshape(ind, sizes.N_e, sizes.d.N_a);
      mat.d.w.(state) = reshape(w, sizes.N_e, sizes.d.N_a);
      inv_h = 1 ./ (grid(ind+1) - grid(ind));
      mat.d.inv_h.(state) = reshape(inv_h, sizes.N_e, sizes.d.N_a);
   end
   % Get inputs for Dynare residual and Jacobian functions 
   [y,x,yagg,yh,xh,mat.d.hist] = heterogeneity.internal.compute_input(M_, out_ss, sizes, mat, indices);
   % Call to the aggregate residual and Jacobian functions
   mat.G = feval([M_.fname '.dynamic_resid'], y, x, M_.params, [], yagg);
   mat.dG = feval([M_.fname '.dynamic_g1'], y, x, M_.params, [], yagg, ...
                   M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, ...
                   M_.dynamic_g1_sparse_colptr);
   % Call the heterogeneous residual and Jacobian functions
   N_sp = sizes.N_e*sizes.pol.N_a;
   F = NaN(H_.endo_nbr, N_sp);
   dF = NaN(H_.endo_nbr, 3*H_.endo_nbr+H_.exo_nbr+3*M_.endo_nbr+M_.exo_nbr, N_sp);
   for j=1:N_sp
      F(:,j) = feval([M_.fname '.dynamic_het1_resid'], y, x, M_.params, [], yh(:,j), xh(:,j), []);
      dF(:,:,j) = feval([M_.fname '.dynamic_het1_g1'], y, x, M_.params, [], yh(:,j), xh(:,j), [], H_.dynamic_g1_sparse_rowval, H_.dynamic_g1_sparse_colval, H_.dynamic_g1_sparse_colptr);
   end
   mat.F = F;
   mat.dF = dF;
   % Store the results in oo_
   oo_het.ss = out_ss;
   oo_het.sizes = sizes;
   oo_het.mat = mat;
   oo_het.indices = indices;
end