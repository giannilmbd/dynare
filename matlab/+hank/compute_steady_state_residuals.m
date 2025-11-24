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
% COMPUTE_STEADY_STATE_RESIDUALS Computes residuals for a heterogeneous-agent model at the steady state
%
% This function evaluates the residuals of both aggregate and individual (heterogeneous-agent) 
% equilibrium conditions given a steady-state object stored in `oo_.hank`. It uses previously 
% computed policy functions and basis matrices to evaluate the model’s residuals pointwise.
%
% INPUTS
%   M_   [struct] : Dynare model structure
%   oo_  [struct] : Dynare output structure. Must contain the fields:
%                    - oo_.hank.ss     : steady-state structure (see initialize_steady_state)
%                    - oo_.hank.sizes  : tensor-product grid sizes
%                    - oo_.hank.mat    : basis and interpolation matrices
%
% OUTPUTS
%   F [matrix] : Residuals of the heterogeneous-agent equations (size: H_.endo_nbr × N_sp),
%                where N_sp is the number of points in the full state space.
%   G [vector] : Residuals of the aggregate equilibrium equations (size: M_.endo_nbr × 1)
%
% NOTES
% - For the aggregate block:
%     - The individual state distribution `ss.d.hist` is reordered to match the policy grid.
%     - Aggregate variables defined as expectations over the distribution are recomputed
%       using `mat.pol.x_bar_dash` and the `Phi` basis matrix.
%     - The residuals are evaluated using the sparse dynamic function `dynamic_resid`.
%
% - For the heterogeneous block:
%     - Constructs the full stacked endogenous (`yh`) and exogenous (`xh`) state vectors
%       at times t-1, t, and t+1 using the stored policy functions and interpolation structure.
%     - Calls `dynamic_het1_resid` at each point in the tensor-product state space.
%
%           dynamic_resid, dynamic_het1_resid
function [F,G] = compute_steady_state_residuals(M_, oo_)
   ss = oo_.hank.ss;
   sizes = oo_.hank.sizes;
   mat = oo_.hank.mat;
   % Rearrange ss.d.hist so that the state and shock orders is similar to the one
   % used for policy functions. order_shocks(i) contains the position of
   % ss.pol.shocks(i) in ss.d.shocks. In other words, ss.d.shocks(order_shocks)
   % equals ss.pol.shocks. A similar logic applies to order_states.
   order_shocks = cellfun(@(x) find(strcmp(x,ss.d.shocks)), ss.pol.shocks);
   order_states = cellfun(@(x) find(strcmp(x,ss.d.states)), ss.pol.states);
   dim_states = cellfun(@(s) sizes.d.states.(s), ss.d.states);
   dim_shocks = cellfun(@(s) sizes.shocks.(s), ss.d.shocks);
   hist = reshape(ss.d.hist, [dim_shocks dim_states]);
   hist = permute(hist, [order_shocks sizes.n_e+order_states]);
   % Compute aggregated heterogeneous variables
   Ix = mat.pol.x_bar_dash*mat.Phi*reshape(hist,[],1);
   % Computing aggregate residuals
   % - Extended aggregate endogenous variable vector at the steady state - %
   y = NaN(3*M_.endo_nbr,1);
   for i=1:M_.orig_endo_nbr
      var = M_.endo_names{i};
      y(i) = ss.agg.(var); 
      y(M_.endo_nbr+i) = ss.agg.(var); 
      y(2*M_.endo_nbr+i) = ss.agg.(var); 
   end
   % - Variables resulting from an aggregation operator - %
   iIx = zeros(size(M_.heterogeneity_aggregates,1),1);
   ix = zeros(size(M_.heterogeneity_aggregates,1),1);
   for i=1:size(M_.heterogeneity_aggregates,1)
      if (M_.heterogeneity_aggregates{i,1} == 'sum')
         iIx(i) = M_.aux_vars(M_.heterogeneity_aggregates{i,2}).endo_index;
         ix(i) = M_.heterogeneity_aggregates{i,3};
      end
   end
   yagg = Ix(ix);
   y(iIx) = yagg;
   y(M_.endo_nbr+iIx) = yagg;
   y(2*M_.endo_nbr+iIx) = yagg;
   % - Exogenous variables - %
   x = zeros(M_.exo_nbr,1);
   % - Call to the residual function - %
   G = feval([M_.fname '.dynamic_resid'], y, x, M_.params, [], yagg);
   % Compute heterogeneous residuals
   N_sp = sizes.N_e*sizes.pol.N_a;
   H_ = M_.heterogeneity(1);
   % Heterogeneous endogenous variable vector
   yh = NaN(3*H_.endo_nbr, N_sp);
   % Heterogeneous exogenous variable vector
   xh = NaN(H_.exo_nbr, N_sp);
   % Get the H_.endo_names indices for the steady-state ordering ss.pol.shocks and ss.pol.states 
   if isfield(ss.shocks, 'Pi')
      order_shocks = cellfun(@(x) find(strcmp(x,H_.endo_names)), ss.pol.shocks);
   else
      order_shocks = cellfun(@(x) find(strcmp(x,H_.exo_names)), ss.pol.shocks);
   end
   % Get the H_.endo_names indices for the steady-state ordering ss.pol.shocks and ss.pol.states 
   order_states = cellfun(@(x) find(strcmp(x,H_.endo_names)), ss.pol.states);
   % Get the indices of originally declared pol values in H_.endo_names
   pol_names = H_.endo_names;
   pol_names(ismember(H_.endo_names, ss.pol.shocks)) = [];
   pol_ind = cellfun(@(x) find(strcmp(x,H_.endo_names)), pol_names);
   % t-1
   yh(order_states,:) = mat.pol.sm(1:sizes.n_a,:);
   % t
   yh(H_.endo_nbr+pol_ind,:) = mat.pol.x_bar;
   if isfield(ss.shocks, 'Pi')
      yh(H_.endo_nbr+order_shocks,:) = mat.pol.sm(sizes.n_a+1:end,:);
   else
      xh(order_shocks,:) = mat.pol.sm(sizes.n_a+1:end,:);
   end
   % t+1
   yh(2*H_.endo_nbr+pol_ind,:) = mat.pol.x_bar_dash*mat.Phi_tilde_e;
   % Call the residual function
   F = NaN(H_.endo_nbr, N_sp);
   for j=1:N_sp
      F(:,j) = feval([M_.fname '.dynamic_het1_resid'], y, x, M_.params, [], yh(:,j), xh(:,j), []);
   end
end
