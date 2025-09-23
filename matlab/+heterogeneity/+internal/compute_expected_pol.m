function expectations = compute_expected_pol(mat, sizes, indices, T, N_Ix)
% Computes expectation matrices for heterogeneous agent dynamics.
%
% INPUTS
% - mat [structure]: matrices from steady-state computation
% - sizes [structure]: size information for various dimensions
% - indices [structure]: permutation and indexing information
% - T [scalar]: number of time periods
% - N_Ix [scalar]: number of aggregated heterogeneous variables
%
% OUTPUTS
% - expectations [4-D array]: expectation matrices (N_e × N_a × N_Ix × T-1)
%
% DESCRIPTION
% Computes the expectation operators needed for the forward-looking
% components of the heterogeneous agent dynamics.

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

N_om = sizes.d.N_a * sizes.N_e;
expectations = zeros(N_om, N_Ix, T-1);

for i = 1:N_Ix
    expectations(:,i,1) = mat.d.x_bar(indices.Ix.in_x(i),:);
end
expectations(:,:,1) = expectations(:,:,1) - sum(expectations(:,:,1), 1) / N_om;

for s = 2:T-1
    expectations(:,:,s) = compute_expected_y(mat.d.ind, mat.d.w, expectations(:,:,s-1), mat.Mu, indices.states, sizes.d.states);
    expectations(:,:,s) = expectations(:,:,s) - sum(expectations(:,:,s), 1) / N_om;
end

end