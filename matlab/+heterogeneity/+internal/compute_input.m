function [y,x,yagg,yh,xh,hist] = compute_input(M_, ss, sizes, mat, indices)
% Computes model input arrays for heterogeneous agent dynamics.
%
% INPUTS
% - M_ [structure]: Dynare's model structure
% - ss [structure]: steady-state structure
% - sizes [structure]: size information for various dimensions
% - mat [structure]: matrices from steady-state computation
% - indices [structure]: permutation and indexing information
%
% OUTPUTS
% - y [array]: aggregate variable values
% - x [array]: heterogeneous variable values
% - yagg [array]: aggregate variable array
% - yh [array]: heterogeneous variable array
% - xh [array]: heterogeneous state array
% - hist [array]: steady-state distribution
%
% DESCRIPTION
% Constructs input arrays needed for model evaluation and residual computation
% in heterogeneous agent models.

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
   % Rearrange ss.d.hist so that the state and shock orders is similar to the one
   % used for policy functions. order_shocks(i) contains the position of
   % ss.pol.shocks(i) in ss.d.shocks. In other words, ss.d.shocks(order_shocks)
   % equals ss.pol.shocks. A similar logic applies to order_states.
   order = cellfun(@(x) find(strcmp(x,ss.d.order)), indices.order);
   hist = reshape(permute(ss.d.hist, order), [], 1);
   % Aggregate inputs:
   [y, yagg, x] = heterogeneity.internal.compute_agg_steady_state(M_, ss.agg, mat.pol.x_bar_dash, mat.d.Phi, hist, indices.Ix.in_x);
   y = repmat(y, 3, 1);
   % Heterogeneous inputs:
   N_sp = sizes.N_e*sizes.pol.N_a;
   H_ = M_.heterogeneity(1);
   % - Heterogeneous endogenous variable vector - %
   yh = NaN(3*H_.endo_nbr, N_sp);
   % - Heterogeneous exogenous variable vector - %
   xh = NaN(H_.exo_nbr, N_sp);
   % Get indices for shocks - Method 1 AR(1) are in H_.endo_names, Method 2a/2b are in H_.exo_names
   order_shocks_endo = cellfun(@(x) find(strcmp(x,H_.endo_names)), indices.shocks.endo);
   order_shocks_exo = cellfun(@(x) find(strcmp(x,H_.exo_names)), indices.shocks.exo);
   % Get the H_.endo_names indices for the steady-state ordering ss.pol.shocks and ss.pol.states 
   order_states = cellfun(@(x) find(strcmp(x,H_.endo_names)), indices.states);
   % Get the indices of originally declared pol values in H_.endo_names
   pol_names = H_.endo_names;
   pol_names(ismember(H_.endo_names, indices.shocks.all)) = [];
   pol_ind = cellfun(@(x) find(strcmp(x,H_.endo_names)), pol_names);
   % t-1
   yh(order_states,:) = mat.pol.sm(sizes.n_e+(1:sizes.n_a),:);
   % t
   yh(H_.endo_nbr+pol_ind,:) = mat.pol.x_bar;
   % Method 1 AR(1) shocks go to yh (endogenous), Method 2a/2b go to xh (exogenous)
   if ~isempty(indices.shocks.endo)
      yh(H_.endo_nbr+order_shocks_endo,:) = mat.pol.sm(1:numel(indices.shocks.endo),:);
   end
   if ~isempty(indices.shocks.exo)
      xh(order_shocks_exo,:) = mat.pol.sm(numel(indices.shocks.endo)+(1:numel(indices.shocks.exo)),:);
   end
   % t+1
   yh(2*H_.endo_nbr+pol_ind,:) = mat.pol.x_bar_dash*mat.pol.Phi_e;
end