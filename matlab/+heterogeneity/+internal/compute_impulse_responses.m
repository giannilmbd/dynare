function impulse_responses = compute_impulse_responses(T, M_, oo_het)
% Computes impulse response functions for heterogeneous agents.
%
% INPUTS
% - T [integer]: Truncation horizon
% - M_ [structure]: Dynare's model structure
% - oo_ [structure]: Dynare's results structure (must contain initialized steady state)
%
% OUTPUTS
% - impulse_responses [4-D array]: impulse responses (N_x × N_Y × N_sp × T)
%
% DESCRIPTION
% Computes the initial impulse response (period 0) and iteratively builds
% the transmission mechanism for subsequent periods using policy matrices.

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

mat = oo_het.mat;
sizes = oo_het.sizes;
indices = oo_het.indices;

N_x = sizes.N_x;
N_sp = sizes.N_sp;
N_Y = sizes.N_Y.all;
n = sizes.n_a;

impulse_responses = zeros(N_x*N_Y, N_sp, T);

H_ = M_.heterogeneity(1);
ind_x = H_.endo_nbr + indices.x.ind.in_endo;
ind_xe = H_.endo_nbr + ind_x;

Ex = compute_expected_dx(mat.pol.x_bar, mat.pol.ind, mat.pol.w, mat.pol.inv_h, mat.Mu, indices.states, sizes.pol.states);
Ae3 = pagemtimes(mat.dF(indices.x.ind.eq, ind_xe, :), Ex);
f3 = mat.dF(indices.x.ind.eq, ind_x, :);
col = indices.x.ind.states;
f3(:, col, :) = f3(:, col, :) + Ae3;

RHS3 = cat(2, mat.dF(indices.x.ind.eq, ind_xe, :), mat.dF(indices.x.ind.eq, indices.Y.ind_in_dF, :));
Z3 = -pagemldivide(f3, RHS3);
cFx_next = Z3(:, 1:N_x, :);
x = Z3(:, N_x+1:end, :);

impulse_responses(:,:,1) = ((reshape(x, [], N_sp) / mat.pol.U) / mat.pol.L) * mat.pol.P;

for s = 2:T
    Ex_dash = impulse_responses(:,:,s-1) * mat.pol.Phi_e;
    x = pagemtimes(cFx_next, reshape(Ex_dash, N_x, N_Y, N_sp));
    impulse_responses(:,:,s) = ((reshape(x, [], N_sp) / mat.pol.U) / mat.pol.L) * mat.pol.P;
end

impulse_responses = reshape(impulse_responses, N_x, N_Y, N_sp, T);

end