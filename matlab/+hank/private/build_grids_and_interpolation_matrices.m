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
% BUILD_GRIDS_AND_INTERPOLATION_MATRICES Constructs interpolation and basis matrices
% for evaluating policy functions and distributions in heterogeneous-agent models.
%
% This function builds all numerical structures needed to evaluate policy residuals,
% compute expectations, and perform interpolation over the state space. It constructs:
% - tensor-product basis matrices (`Phi`)
% - expectation operators (`Phi_tilde_e`)
% - grid point matrices (`pol.sm`, `d.sm`)
% - LU decomposition of the interpolation basis for efficient linear solves.
%
% INPUTS
%   ss    [struct] : Steady-state structure validated by `check_steady_state_input`
%   sizes [struct] : Structure containing grid dimensions and counts
%
% OUTPUT
%   mat   [struct] : Structure containing matrices used in solution and simulation:
%       - Phi            : basis matrix for evaluating interpolated functions over distribution grid
%       - Phi_tilde      : expanded identity matrix for state variables (used in solving)
%       - Phi_tilde_e    : expectation operator incorporating both interpolation and transition shocks
%       - Mu             : Kronecker product of exogenous transition or quadrature weights
%       - pol.sm         : full tensor grid of state and shock combinations for policy evaluation
%       - d.sm           : full tensor grid of state and shock combinations for distribution evaluation
%       - L_Phi_tilde,
%         U_Phi_tilde,
%         P_Phi_tilde    : LU decomposition of `Phi_tilde` (used to project policy values)
%
% Notes:
% - The structure `ss` must contain valid grids and values in `ss.pol.values`, `ss.pol.states`,
%   `ss.pol.shocks`, `ss.d.grids`, and `ss.shocks`.
% - Assumes linear interpolation and separable basis construction across state variables.
% - `Phi_tilde_e` applies the expected operator `E_t[f(x')]` for the simulation step.
function mat = build_grids_and_interpolation_matrices(ss, sizes)
   mat = struct;
   % Compute Φ and ̃Φ with useful basis matrices
   Phi = speye(sizes.N_e);
   Phi_tilde = speye(sizes.N_e);
   B = struct;
   B_tilde_bar = struct;
   for i=sizes.n_a:-1:1
      s = ss.pol.states{i};
      [ind, w] = find_bracket_linear_weight(ss.pol.grids.(s), ss.d.grids.(s));
      B.(s) = lin_spli(ind, w, sizes.pol.states.(s));
      [ind, w] = find_bracket_linear_weight(ss.pol.grids.(s), reshape(ss.pol.values.(s), 1, []));
      B_tilde_bar.(s) = lin_spli(ind, w, sizes.pol.states.(s));
      Phi = kron(B.(s), Phi);
      Phi_tilde = kron(speye(sizes.pol.states.(s)), Phi_tilde);
   end
   mat.Phi = Phi;
   mat.Phi_tilde = Phi_tilde;
   % Compute ̃Φᵉ and Mu
   N_sp = sizes.pol.N_a*sizes.N_e;
   hadamard_prod = ones(sizes.pol.N_a,N_sp);
   prod_a_1_to_jm1 = 1;
   prod_a_jp1_to_n_a = sizes.pol.N_a;
   for j=1:sizes.n_a
      s = ss.pol.states{j};
      prod_a_jp1_to_n_a = prod_a_jp1_to_n_a/sizes.pol.states.(s);
      hadamard_prod = hadamard_prod .* kron(kron(ones(prod_a_1_to_jm1,1), B_tilde_bar.(s)), ones(prod_a_jp1_to_n_a,1));
      prod_a_1_to_jm1 = prod_a_1_to_jm1*sizes.pol.states.(s);
   end
   Mu = eye(1);
   if isfield(ss.shocks, 'Pi')
      for i=sizes.n_e:-1:1
         e = ss.pol.shocks{i};
         Mu = kron(ss.shocks.Pi.(e), Mu);
      end
      Phi_tilde_e = kron(hadamard_prod, ones(sizes.N_e,1)).*kron(ones(sizes.pol.N_a), Mu');
   else
      for i=sizes.n_e:-1:1
         e = ss.pol.shocks{i};
         Mu = kron(ss.shocks.w.(e), Mu);
      end
      Phi_tilde_e = kron(hadamard_prod, Mu);
   end
   mat.Mu = Mu;
   mat.Phi_tilde_e = Phi_tilde_e;
   % Compute state and shock grids for the policy and distribution state grids
   mat.pol.sm = set_state_matrix(ss, sizes, "pol");
   mat.d.sm = set_state_matrix(ss, sizes, "d");
   % LU decomposition of Phi_tilde 
   [L_Phi_tilde,U_Phi_tilde,P_Phi_tilde] = lu(Phi_tilde);
   mat.L_Phi_tilde = L_Phi_tilde;
   mat.U_Phi_tilde = U_Phi_tilde;
   mat.P_Phi_tilde = P_Phi_tilde;
end