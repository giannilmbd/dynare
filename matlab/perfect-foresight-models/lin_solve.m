function [x, first_iter_lu] = lin_solve(A, b, options_, first_iter_lu, force_lu)
% [x, first_iter_lu] = lin_solve(A, b, options_, first_iter_lu, force_lu)
% Solves the linear system A·x=b. Used at the heart of the perfect foresight solver when
% stack_solve_algo equals 0 (LU), 2 (GMRES) or 3 (BiCGStab).
%
% If force_lu is true, then the value of options_.stack_solve_algo is ignored and a LU is used.
%
% first_iter_lu corresponds to the preconditioner used when preconditioner=first_iter_lu.
% If empty on input, then the routine computes the preconditioner and returns on output.
% If not empty, use that preconditioner without recomputing it, and pass it unmodified on output.

% Copyright © 1996-2025 Dynare Team
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

if nargin < 4
    first_iter_lu = [];
end
if nargin < 5
    force_lu = false;
end

if norm(b) < sqrt(eps) % then x = 0 is a solution
    x = 0;
    return
end

if options_.stack_solve_algo == 0 || force_lu
    x = A\b;
else % Iterative algorithm
    if strcmp(options_.simul.preconditioner, 'first_iter_lu')
        if isempty(first_iter_lu)
            [L, U, P, Q] = lu(A);
        else
            L = first_iter_lu.L;
            U = first_iter_lu.U;
            P = first_iter_lu.P;
            Q = first_iter_lu.Q;
        end
    elseif strcmp(options_.simul.preconditioner, 'block_diagonal_lu')
        [L, U, P, Q] = block_diagonal_lu_preconditioner(A, options_);
    elseif strcmp(options_.simul.preconditioner, 'incomplete_lu')
        [L, U, P] = ilu(A, options_.simul.incomplete_lu);
        Q = speye(size(A));
    end

    if strcmp(options_.simul.preconditioner, 'first_iter_lu') && isempty(first_iter_lu)
        z = L\(P*b);
        y = U\z;
        first_iter_lu = struct('L', L, 'U', U, 'P', P, 'Q', Q);
    else
        iter_tol = options_.simul.iter_tol;
        if isempty(iter_tol)
            % The tolerance passed to gmres and bicgstab is a relative one (‖Ax−b‖/‖b‖).
            % However, we test the convergence of algorithms via options_.dynatol.f, which is an
            % absolute error (‖Ax−b‖). Hence the need to rescale by the norm of the RHS (‖b‖).
            iter_tol = options_.dynatol.f / max(abs(b)) / 10;
        end

        iter_maxit = min(options_.simul.iter_maxit, size(A, 1));

        if options_.stack_solve_algo == 2
            gmres_restart = min(options_.simul.gmres_restart, size(A, 1));
            [y, flag] = gmres(P*A*Q, P*b, gmres_restart, iter_tol, iter_maxit, L, U);
        elseif options_.stack_solve_algo == 3
            [y, flag] = bicgstab(P*A*Q, P*b, iter_tol, iter_maxit, L, U);
        else
            error('Unsupported value for options_.stack_solve_algo')
        end
        iter_solver_error_flag(flag)
    end

    x = Q*y;
end
if ~options_.simul.allow_nonfinite_values
    x(~isfinite(x)) = 0; %prevent non-finite values from propagating, see #1975
end


function [L, U, P, Q] = block_diagonal_lu_preconditioner(A, options_)

ny = size(A, 1) / options_.periods;

if options_.simul.block_diagonal_lu_nperiods > 0
    if options_.simul.block_diagonal_lu_nperiods > options_.periods
        error('block_diagonal_lu preconditioner: block_diagonal_lu_nperiods option is greater than the periods options')
    end
    lusize = ny * options_.simul.block_diagonal_lu_nperiods;
elseif options_.simul.block_diagonal_lu_nlu ~= 0
    if options_.simul.block_diagonal_lu_nlu < 0 % To avoid misleading TROLL users
        error('block_diagonal_lu preconditioner: negative value of block_diagonal_lu_nlu option is not supported')
    end
    lusize = floor(floor(size(A, 1) / options_.simul.block_diagonal_lu_nlu) / ny) * ny;
    if lusize == 0
        error('block_diagonal_lu preconditioner: block_diagonal_lu_nlu option is too large')
    end
else
    if size(A, 1) < options_.simul.block_diagonal_lu_maxlu
        error('block_diagonal_lu preconditioner: too small problem; either decrease block_diagonal_lu_maxlu option, or use block_diagonal_lu_nperiods or block_diagonal_lu_nlu options')
    end
    if ny > options_.simul.block_diagonal_lu_maxlu
        error('block_diagonal_lu preconditioner: too large problem; either increase block_diagonal_lu_maxlu option, or use block_diagonal_lu_nperiods or block_diagonal_lu_nlu options')
    end
    lusize = floor(options_.simul.block_diagonal_lu_maxlu / ny) * ny;
end

if options_.simul.block_diagonal_lu_relu < 0 || options_.simul.block_diagonal_lu_relu > 1
    error('block_diagonal_lu preconditioner: block_diagonal_lu_relu option must be between 0 and 1')
end

luidx = floor((floor(size(A, 1) / lusize) - 1) * options_.simul.block_diagonal_lu_relu) * lusize + (1:lusize);
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
    fprintf('block_diagonal_lu preconditioner: main LU size = %d, remainder LU size = %d\n', lusize, r)
end


function iter_solver_error_flag(flag)
% Interprets error flag from gmres or bicgstab functions

switch flag
    case 1
        error('Maximum number of iterations exceeded in GMRES/BiCGStab. You may want to increase the iter_maxit option or the gmres_restart option (if using stack_solve_algo=2 for the latter), or decrease iter_tol.')
    case 2
        error('The preconditioner matrix is ill conditioned in GMRES/BiCGStab. You should try another preconditioner.')
    case 3
        error('No progress between two successive iterations in GMRES/BiCGStab. You may want to increase the iter_tol option.')
    case 4
        error('One of the scalar quantities calculated by GMRES/BiCGStab became too small or too large. Try another preconditioner or algorithm.')
end
