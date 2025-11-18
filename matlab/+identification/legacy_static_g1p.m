function g1p_legacy = legacy_static_g1p(g1p, M_)
% Given g1p as returned by static_params_derivs.m, construct the g1p matrix in
% the legacy representation (see issue #1859 for the historical background).
% Ideally this file should go away when the identification code is adapted to
% the new model representation.

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

g1p_legacy = zeros(M_.eq_nbr, M_.endo_nbr, M_.param_nbr);

for i = 1:size(g1p, 1)
    g1p_legacy(g1p(i,1), g1p(i,2), g1p(i,3)) = g1p(i,4);
end
