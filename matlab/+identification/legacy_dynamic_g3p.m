function g3p_legacy = legacy_dynamic_g3p(g3p, M_)
% Given g3p as returned by dynamic_params_derivs.m, construct the g3p matrix in
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

g3p_legacy = g3p;
for i = 1:size(g3p, 1)
    g3p_legacy(i,2) = identification.legacy_idx(g3p(i,2), M_);
    g3p_legacy(i,3) = identification.legacy_idx(g3p(i,3), M_);
    g3p_legacy(i,4) = identification.legacy_idx(g3p(i,4), M_);
end
