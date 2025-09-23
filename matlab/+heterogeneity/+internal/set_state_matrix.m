function sm = set_state_matrix(pol_grids, shocks_grids, sizes, indices)
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
%   sm     [matrix] : Matrix of size (n_a + n_e) × (N_a * N_e), where each column contains
%                     one combination of state and shock values, ordered lexicographically.
%
% Notes:
% - The first `n_a` rows correspond to endogenous states; the next `n_e` to exogenous shocks.
% - Ordering is compatible with Kronecker-product-based interpolation logic.
% - Used in residual evaluation, interpolation, and expectation computation routines.

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
   % Useful size vectors and variables
   N = sizes.pol.N_a*sizes.N_e;
   
   % Setting the states matrix
   sm = zeros(N, sizes.n_a+sizes.n_e);
   
   % Filling the state matrix
   n_repmat = N;
   n_repelem = 1;
   for j=1:sizes.n_e
      shock = indices.shocks.all{j};
      n_repmat = n_repmat/sizes.shocks.(shock);
      sm(:,j) = repmat(repelem(shocks_grids.(shock), n_repelem, 1), n_repmat, 1);
      n_repelem = n_repelem*sizes.shocks.(shock);
   end
   for j=1:sizes.n_a
      state = indices.states{j};
      n_repmat = n_repmat/sizes.pol.states.(state);
      sm(:,sizes.n_e+j) = repmat(repelem(pol_grids.(state), n_repelem, 1), n_repmat, 1);
      n_repelem = n_repelem*sizes.pol.states.(state);
   end
   sm = sm';
end