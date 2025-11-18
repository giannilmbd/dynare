function g1pp_legacy = legacy_dynamic_g1pp(g1pp, M_)
% Given g1pp as returned by dynamic_params_derivs.m, construct the g1pp matrix
% in the legacy representation (see issue #1859 for the historical background).
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

g1pp_legacy = g1pp;
next_sym_idx = size(g1pp, 1) + 1;
for i = 1:size(g1pp, 1)
    g1pp_legacy(i,2) = identification.legacy_idx(g1pp(i,2), M_);

    %% Add symmetric elements, which are not included in the sparse representation
    if g1pp_legacy(i,3) ~= g1pp_legacy(i,4)
        g1pp_legacy(next_sym_idx,:) = g1pp_legacy(i,:);
        g1pp_legacy(next_sym_idx,3) = g1pp_legacy(i,4);
        g1pp_legacy(next_sym_idx,4) = g1pp_legacy(i,3);
        next_sym_idx = next_sym_idx + 1;
    end
end
