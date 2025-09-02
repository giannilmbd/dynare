function [x, umfiter_precond] = lin_solve(A, b, options_, umfiter_precond, force_lu)
% Solves the linear system A·x=b. Used at the heart of the perfect foresight solver when
% stack_solve_algo equals 0 (LU), 2 (GMRES) or 3 (BiCGStab).
%
% If force_lu is true, then the value of options_.stack_solve_algo is ignored and a LU is used.
%
% umfiter_precond corresponds to the preconditioner used when precondioner=umfiter.
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
    umfiter_precond = [];
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
    if strcmp(options_.simul.preconditioner, 'umfiter')
        if isempty(umfiter_precond)
            [L, U, P, Q] = lu(A);
        else
            L = umfiter_precond.L;
            U = umfiter_precond.U;
            P = umfiter_precond.P;
            Q = umfiter_precond.Q;
        end
    elseif strcmp(options_.simul.preconditioner, 'iterstack')
        [L, U, P, Q] = iterstack_preconditioner(A, options_);
    elseif strcmp(options_.simul.preconditioner, 'ilu')
        [L, U, P] = ilu(A, options_.simul.ilu);
        Q = speye(size(A));
    end

    if strcmp(options_.simul.preconditioner, 'umfiter') && isempty(umfiter_precond)
        z = L\(P*b);
        y = U\z;
        umfiter_precond = struct('L', L, 'U', U, 'P', P, 'Q', Q);
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


function [L, U, P, Q] = iterstack_preconditioner(A, options_)

ny = size(A, 1) / options_.periods;

if options_.simul.iterstack_nperiods > 0
    if options_.simul.iterstack_nperiods > options_.periods
        error('iterstack preconditioner: iterstack_nperiods option is greater than the periods options')
    end
    lusize = ny * options_.simul.iterstack_nperiods;
elseif options_.simul.iterstack_nlu ~= 0
    if options_.simul.iterstack_nlu < 0 % To avoid misleading TROLL users
        error('iterstack preconditioner: negative value of iterstack_nlu option is not supported')
    end
    lusize = floor(floor(size(A, 1) / options_.simul.iterstack_nlu) / ny) * ny;
    if lusize == 0
        error('iterstack preconditioner: iterstack_nlu option is too large')
    end
else
    if size(A, 1) < options_.simul.iterstack_maxlu
        error('iterstack preconditioner: too small problem; either decrease iterstack_maxlu option, or use iterstack_nperiods or iterstack_nlu options')
    end
    if ny > options_.simul.iterstack_maxlu
        error('iterstack preconditioner: too large problem; either increase iterstack_maxlu option, or use iterstack_nperiods or iterstack_nlu options')
    end
    lusize = floor(options_.simul.iterstack_maxlu / ny) * ny;
end

if options_.simul.iterstack_relu < 0 || options_.simul.iterstack_relu > 1
    error('iterstack preconditioner: iterstack_relu option must be between 0 and 1')
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
    fprintf('iterstack preconditioner: main LU size = %d, remainder LU size = %d\n', lusize, r)
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
