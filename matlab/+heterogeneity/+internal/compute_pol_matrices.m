function [x_bar, x_bar_dash] = compute_pol_matrices(ss, sizes, mat, indices)
% Construct interpolated policy function matrices for state transitions.
%
% Given a steady-state structure `ss` with discretized policy functions, this function builds:
% - `x_bar`: matrix of policy function values (flattened over the joint state space)
% - `x_bar_dash`: transformed version of `x_bar` projected onto the basis used for expectations
%
% INPUTS
%   ss      [struct] : Validated steady-state structure
%   sizes   [struct] : Structure with grid sizes
%   mat     [struct] : Structure containing transformation matrices (e.g., `U_Phi_tilde`)
%   indices [struct] : Index structure with variable names and mappings
%
% OUTPUTS
%   x_bar       [matrix] : Matrix of policy function values (N_x × N_sp)
%   x_bar_dash  [matrix] : Projected matrix in basis representation for expectations (N_x × N_sp)

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
   % Computing x_bar
   x_bar = NaN(sizes.N_x, sizes.N_sp);
   for i=1:sizes.N_x
      x = indices.x.names{i};
      if isfield(ss.pol.values, x)
         x_bar(i,:) = reshape(ss.pol.values.(x),1,[]);
      end
   end
   % Computing x_bar^#
   x_bar_dash = ((x_bar / mat.pol.U) / mat.pol.L) * mat.pol.P;
end 