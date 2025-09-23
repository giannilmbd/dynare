function heterogeneous_jacobian = compute_heterogeneous_jacobian(transition_matrices)
% Computes heterogeneous agent Jacobian from transition matrices.
%
% INPUTS
% - transition_matrices [4-D array]: F matrices (T × T × N_Y × N_Ix)
%
% OUTPUTS
% - heterogeneous_jacobian [4-D array]: heterogeneous agent Jacobian (T × T × N_Y × N_Ix)
%
% DESCRIPTION
% Transforms the transition matrices into the heterogeneous agent Jacobian
% by accumulating effects over time periods.

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

T = size(transition_matrices, 1);
heterogeneous_jacobian = transition_matrices;
for t = 2:T
    heterogeneous_jacobian(2:T,t,:,:) = heterogeneous_jacobian(2:T,t,:,:) + heterogeneous_jacobian(1:T-1,t-1,:,:);
end

end