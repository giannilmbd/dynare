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
% INITIALIZE_STEADY_STATE Wrapper to validate, complete, and store steady-state input for HANK models.
%
% This function performs validation of the user-provided steady-state structure `ss`, constructs
% interpolation objects and basis matrices, computes missing policy functions for complementarity 
% multipliers (if not supplied), and populates the global `oo_.hank` structure with the results.
%
% INPUTS
%   M_       [struct] : Dynare model structure
%   options_ [struct] : Dynare options structure
%   oo_      [struct] : Dynare output structure to which results will be written
%   ss       [struct] : User-supplied steady-state structure (partial or complete)
% OUTPUTS
%   oo_      [struct] : Updated Dynare output structure containing:
%                        - oo_.hank.ss    : validated and completed steady-state structure
%                        - oo_.hank.sizes : grid size structure
%                        - oo_.hank.mat   : interpolation and policy function matrices
function oo_ = initialize_steady_state(M_, options_, oo_, ss)
   assert(isscalar(M_.heterogeneity), 'Heterogeneous-agent models with more than one heterogeneity dimension are not allowed yet!');
   % Check steady-state input
   [out_ss, sizes] = hank.check_steady_state_input(M_, options_, ss);
   % Build useful basis and interpolation matrices
   H_ = M_.heterogeneity(1);
   mat = build_grids_and_interpolation_matrices(out_ss, sizes);
   % Compute policy matrices without the complementarity conditions multipliers
   [mat.pol.x_bar, mat.pol.x_bar_dash] = compute_pol_matrices(H_, out_ss, sizes, mat);
   % Compute the policy values of complementarity conditions multipliers
   mult_values = compute_pol_mcp_mul(M_, out_ss, sizes, mat);
   % Update ss.pol.values accordingly
   f = fieldnames(mult_values);
   for i=1:numel(f)
     out_ss.pol.values.(f{i}) = mult_values.(f{i});
   end
   % Update the policy matrices
   [mat.pol.x_bar, mat.pol.x_bar_dash] = compute_pol_matrices(H_, out_ss, sizes, mat);
   % Store the results in oo_
   oo_.hank.ss = out_ss;
   oo_.hank.sizes = sizes;
   oo_.hank.mat = mat;
end