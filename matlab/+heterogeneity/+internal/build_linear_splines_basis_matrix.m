function B = build_linear_splines_basis_matrix(ind, w, m)
% Constructs a sparse linear interpolation matrix for 1D bracketed data.
%
% This function builds a sparse matrix `B` of size m×n that maps `n` query points,
% located between known grid points, to the surrounding `m`-dimensional grid
% using linear interpolation weights.
%
% It is typically used in Dynare heterogeneity-specific routines after calling
% `find_bracket_linear_weight` MEX to form interpolation
% matrices over endogenous grids.
%
% INPUTS
%   ind [n×1 int32]   : Lower bracket indices for each query (i.e., index `ilow` s.t.
%                       x(ind(i)) ≤ xq(i) ≤ x(ind(i)+1))
%
%   w   [n×1 double]  : Linear interpolation weights (relative position in interval):
%                       w(i) = (x(ind(i)+1) - xq(i)) / (x(ind(i)+1) - x(ind(i)))
%
%   m   [int scalar]  : Size of the full grid (number of basis functions, i.e., length of `x`)
%
% OUTPUT
%   B   [m×n sparse]  : Sparse interpolation matrix such that:
%                         f_interp = B * f_grid
%                       where `f_grid` is a column vector of function values on the grid.
% NOTES
% - Each column of `B` contains exactly two nonzero entries: `w(i)` and `1 - w(i)`
% - Assumes linear interpolation over uniform or non-uniform 1D grids

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
    i = [ind;ind+1];
    n = length(ind);
    j = repmat(int32(1:n)',2,1);
    v = [w;1-w];
    B = sparse(i,j,v,m,n);
end