function [L, U, P, Q] = iterstack_preconditioner(A, options_)
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

ny = size(A, 1) / options_.periods;

if options_.simul.iterstack_nperiods > 0
    if options_.simul.iterstack_nperiods > options_.periods
        error('iterstack preconditionner: iterstack_nperiods option is greater than the periods options')
    end
    lusize = ny * options_.simul.iterstack_nperiods;
elseif options_.simul.iterstack_nlu ~= 0
    if options_.simul.iterstack_nlu < 0 % To avoid misleading TROLL users
        error('iterstack preconditionner: negative value of iterstack_nlu option is not supported')
    end
    lusize = floor(floor(size(A, 1) / options_.simul.iterstack_nlu) / ny) * ny;
    if lusize == 0
        error('iterstack preconditionner: iterstack_nlu option is too large')
    end
else
    if size(A, 1) < options_.simul.iterstack_maxlu
        error('iterstack preconditionner: too small problem; either decrease iterstack_maxlu option, or use iterstack_nperiods or iterstack_nlu options')
    end
    if ny > options_.simul.iterstack_maxlu
        error('iterstack preconditionner: too large problem; either increase iterstack_maxlu option, or use iterstack_nperiods or iterstack_nlu options')
    end
    lusize = floor(options_.simul.iterstack_maxlu / ny) * ny;
end

if options_.simul.iterstack_relu < 0 || options_.simul.iterstack_relu > 1
    error('iterstack preconditionner: iterstack_relu option must be between 0 and 1')
end

luidx = floor((floor(size(A, 1) / lusize) - 1) * options_.simul.iterstack_relu) * lusize + (1:lusize);
[L1, U1, P1, Q1] = lu(A(luidx,luidx));
eyek = speye(floor(size(A, 1) / lusize));
L = kron(eyek, L1);
U = kron(eyek, U1);
P = kron(eyek, P1);
Q = kron(eyek, Q1);
r = rem(size(A, 1), lusize);
if r > 0 % Compute additional smaller LU for remainder, if any
    [L2, U2, P2, Q2] = lu(A(end-r+1:end,end-r+1:end));
    L = blkdiag(L, L2);
    U = blkdiag(U, U2);
    P = blkdiag(P, P2);
    Q = blkdiag(Q, Q2);
end

if options_.debug
    fprintf('iterstack preconditionner: main LU size = %d, remainder LU size = %d\n', lusize, r)
end
