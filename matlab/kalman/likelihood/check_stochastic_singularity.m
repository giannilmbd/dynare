function null_space = check_stochastic_singularity(F, d_index, varobs, t, matrix_title, F_tolerance, null_space_tolerance)
% Checks and reports observables involved in stochastic singularity.
%
% This function inspects the prediction-error covariance matrix
% F. If it is rank-deficient, it computes a basis for its null
% space and prints the implied rank-deficient linear combinations of
% observables.
%
% INPUTS
% - F                     [double]      [n x n] prediction-error covariance matrix at time t
% - d_index               [vector]      [n x 1] indices of observables used at time t (indices in varobs)
% - varobs                [cell]        [p x 1] cell array with names of observed variables
% - t                     [integer]     time index (typically diffuse period)
% - matrix_title          [string]      matrix title for printout
% - F_tolerance           [double]      (optional) threshold used to zero-out spurious entries in F (default: 1e-10)
% - null_space_tolerance  [double]      (optional) reporting threshold for null-space coefficients (default: 1e-8)
%
% OUTPUTS
% - null_space            [double]      basis vectors (columns) spanning the null space of cleaned F
%
% This function is called by: none
% This function calls: null, isoctave, matlab_ver_less_than

% Copyright © 2004-2026 Dynare Team
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

if nargin < 5 || isempty(matrix_title)
    matrix_title = 'forecast-error variance matrix';
end

if nargin < 6 || isempty(F_tolerance)
    F_tolerance = 1e-10;
end

if nargin < 7 || isempty(null_space_tolerance)
    null_space_tolerance = 1e-8;
end

if size(F, 1) ~= size(F, 2)
    error('check_stochastic_singularity: F must be a square matrix.');
end

if size(F, 1) ~= length(d_index)
    error('check_stochastic_singularity: length(d_index) must match size(F,1).');
end

if ~iscell(varobs)
    error('check_stochastic_singularity: varobs must be a cell array of variable names.');
end

if any(d_index < 1) || any(d_index > length(varobs))
    error('check_stochastic_singularity: d_index contains indices outside varobs.');
end

F(abs(F) < F_tolerance) = 0;
null_space = compute_nullspace(F, []);

fprintf('\nSingularity detected in %s at t=%d.\n',matrix_title, t);
fprintf('Rank-deficient linear combination(s) of observables:\n');
if isempty(null_space)
    fprintf('  None found (null space is empty at current numerical tolerance).\n');
    return
end

for j = 1:size(null_space, 2)
    fprintf('  Combination %d: ', j);
    for i = 1:length(d_index)
        if abs(null_space(i,j)) > null_space_tolerance
            fprintf('%s (%.4f)  ', varobs{d_index(i)}, null_space(i,j));
        end
    end
    fprintf('\n');
end