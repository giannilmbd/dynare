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
% SET_STATE_MATRIX Builds the full tensor-product grid of state and shock variables.
%
% This function constructs a matrix containing all combinations of state and shock values 
% used for evaluating policy functions or distributions. The result is a matrix of size 
% (n_a + n_e) × (N_a * N_e), where each column corresponds to one grid point in the joint 
% state space.
%
% INPUTS
%   ss     [struct] : Steady-state structure containing grids and state/shock names
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
function sm = set_state_matrix(ss, sizes, field)
   % Useful size vectors and variables
   f = ss.(field);
   N = sizes.(field).N_a*sizes.N_e;
   
   % Setting the states matrix
   sm = zeros(N, sizes.n_a+sizes.n_e);
   
   % Filling the state matrix
   prod_a_1_to_jm1 = 1;
   prod_a_jp1_to_n_a = N;
   for j=1:sizes.n_a
      state = f.states{j};
      prod_a_jp1_to_n_a = prod_a_jp1_to_n_a/sizes.(field).states.(state);
      sm(:,j) = kron(ones(prod_a_1_to_jm1,1), kron(f.grids.(state), ones(prod_a_jp1_to_n_a,1)));
      prod_a_1_to_jm1 = prod_a_1_to_jm1*sizes.(field).states.(state);
   end
   prod_theta_1_to_jm1 = sizes.(field).N_a;
   prod_theta_jp1_to_n_theta = sizes.N_e;
   for j=1:sizes.n_e
      shock = f.shocks{j};
      prod_theta_jp1_to_n_theta = prod_theta_jp1_to_n_theta/sizes.shocks.(shock);
      sm(:,sizes.n_a+j) = kron(ones(prod_theta_1_to_jm1,1), kron(ss.shocks.grids.(shock), ones(prod_theta_jp1_to_n_theta,1)));
      prod_theta_1_to_jm1 = prod_theta_1_to_jm1*sizes.shocks.(shock);
   end
   sm = sm';
end