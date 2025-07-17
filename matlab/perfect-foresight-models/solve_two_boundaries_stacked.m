function [y, T, success, max_res, iter] = solve_two_boundaries_stacked(fh, y, x, steady_state, T, Block_Num, cutoff, options_, M_)
% Computes the deterministic simulation of a block of equations containing
% both lead and lag variables, using a Newton method over the stacked Jacobian
% (in particular, this excludes LBJ).
%
% INPUTS
%   fh                  [handle]        function handle to the dynamic file for the block
%   y                   [matrix]        All the endogenous variables of the model
%   x                   [matrix]        All the exogenous variables of the model
%   steady_state        [vector]        steady state of the model
%   T                   [matrix]        Temporary terms
%   Block_Num           [integer]       block number
%   cutoff              [double]        cutoff to correct the direction in Newton in case
%                                       of singular Jacobian matrix
%   options_             [structure]     storing the options
%   M_                   [structure]     Model description
%
% OUTPUTS
%   y                   [matrix]        All endogenous variables of the model
%   T                   [matrix]        Temporary terms
%   success             [logical]       Whether a solution was found
%   max_res             [double]        ∞-norm of the residual
%   iter                [integer]       Number of iterations
%
% ALGORITHM
%   Newton with LU or GMRES or BiCGStab

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

Blck_size = M_.block_structure.block(Block_Num).mfs;
y_index = M_.block_structure.block(Block_Num).variable(end-Blck_size+1:end);
periods = get_simulation_periods(options_);
y_kmin = M_.maximum_lag;
stack_solve_algo = options_.stack_solve_algo;

if ~ismember(stack_solve_algo, [0 2 3 4])
    error('Unsupported stack_solve_algo value')
end

verbose = options_.verbosity;

cvg=false;
iter=0;
correcting_factor=0.01;
max_resa=1e100;
lambda = 1; % Length of Newton step (unused for stack_solve_algo=4)
while ~(cvg || iter > options_.simul.maxit)
    r = NaN(Blck_size, periods);
    g1a = spalloc(Blck_size*periods, Blck_size*periods, M_.block_structure.block(Block_Num).NNZDerivatives*periods);
    for it_ = y_kmin+(1:periods)
        [yy, T(:, it_), r(:, it_-y_kmin), g1]=fh(dynendo(y, it_, M_), x(it_, :), M_.params, steady_state, ...
                                                 M_.block_structure.block(Block_Num).g1_sparse_rowval, ...
                                                 M_.block_structure.block(Block_Num).g1_sparse_colval, ...
                                                 M_.block_structure.block(Block_Num).g1_sparse_colptr, T(:, it_));
        y(:, it_) = yy(M_.endo_nbr+(1:M_.endo_nbr));
        if periods == 1
            g1a = g1(:, Blck_size+(1:Blck_size));
        elseif it_ == y_kmin+1
            g1a(1:Blck_size, 1:Blck_size*2) = g1(:, Blck_size+1:end);
        elseif it_ == y_kmin+periods
            g1a((periods-1)*Blck_size+1:end, (periods-2)*Blck_size+1:end) = g1(:, 1:2*Blck_size);
        else
            g1a((it_-y_kmin-1)*Blck_size+(1:Blck_size), (it_-y_kmin-2)*Blck_size+(1:3*Blck_size)) = g1;
        end
    end
    ya = reshape(y(y_index, y_kmin+(1:periods)), 1, periods*Blck_size)';
    ra = reshape(r, periods*Blck_size, 1);
    b=-ra+g1a*ya;
    [max_res, max_indx]=max(max(abs(r')));
    if ~isreal(r)
        max_res = (-max_res^2)^0.5;
    end
    if ~isreal(max_res) || isnan(max_res)
        cvg = false;
    elseif M_.block_structure.block(Block_Num).is_linear && iter>0
        cvg = true;
    else
        cvg = (max_res < options_.dynatol.f);
    end
    if ~cvg
        if iter>0
            if ~isreal(max_res) || isnan(max_res) %|| (max_resa<max_res && iter>1)
                if verbose && ~isreal(max_res)
                    disp(['Variable ' M_.endo_names{max_indx} ' (' int2str(max_indx) ') returns an undefined value']);
                end
                if isnan(max_res)
                    detJ=det(g1aa);
                    if abs(detJ)<1e-7
                        max_factor=max(max(abs(g1aa)));
                        ze_elem=sum(diag(g1aa)<cutoff);
                        if verbose
                            disp([num2str(full(ze_elem),'%d') ' elements on the Jacobian diagonal are below the cutoff (' num2str(cutoff,'%f') ')']);
                        end
                        if correcting_factor<max_factor
                            correcting_factor=correcting_factor*4;
                            if verbose
                                disp(['The Jacobian matrix is singular, det(Jacobian)=' num2str(detJ,'%f') '.']);
                                disp('    trying to correct the Jacobian matrix:');
                                disp(['    correcting_factor=' num2str(correcting_factor,'%f') ' max(Jacobian)=' num2str(full(max_factor),'%f')]);
                            end
                            dx = (g1aa+correcting_factor*speye(periods*Blck_size))\ba- ya_save;
                            y(y_index, y_kmin+(1:periods))=reshape((ya_save+lambda*dx)',length(y_index),periods);
                            continue
                        else
                            disp('The singularity of the Jacobian matrix could not be corrected');
                            success = false;
                            return
                        end
                    end
                elseif lambda>1e-8 && stack_solve_algo ~= 4
                    lambda=lambda/2;
                    if verbose
                        disp(['reducing the path length: lambda=' num2str(lambda,'%f')]);
                    end
                    y(y_index, y_kmin+(1:periods))=reshape((ya_save+lambda*dx)',length(y_index),periods);
                    continue
                else
                    if verbose
                        if cutoff==0
                            fprintf('Convergence not achieved in block %d, after %d iterations.\n Increase "maxit".\n',Block_Num, iter);
                        else
                            fprintf('Convergence not achieved in block %d, after %d iterations.\n Increase "maxit" or set "cutoff=0" in model options.\n',Block_Num, iter);
                        end
                    end
                    success = false;
                    return
                end
            else
                if lambda<1 && stack_solve_algo ~= 4
                    lambda=max(lambda*2, 1);
                end
            end
        end
        ya_save=ya;
        g1aa=g1a;
        ba=b;
        max_resa=max_res;
        if stack_solve_algo==0 || (ismember(stack_solve_algo, [2 3]) && strcmp(options_.simul.preconditioner, 'iterstack') ...
                                   && options_.simul.iterstack_nperiods == 0 && options_.simul.iterstack_nlu == 0 && size(g1a, 1) < options_.simul.iterstack_maxlu) % Fallback to LU if block too small for iterstack
            dx = g1a\b- ya;
            ya = ya + lambda*dx;
            y(y_index, y_kmin+(1:periods))=reshape(ya',length(y_index),periods);
        elseif ismember(stack_solve_algo, [2, 3])
            if strcmp(options_.simul.preconditioner, 'umfiter') && iter == 0
                [L, U, P, Q] = lu(g1a);
                zc = L\(P*b);
                zb = U\zc;
            elseif strcmp(options_.simul.preconditioner, 'iterstack')
                [L, U, P, Q] = iterstack_preconditioner(g1a, options_);
            elseif strcmp(options_.simul.preconditioner, 'ilu')
                [L, U, P] = ilu(g1a, options_.simul.ilu);
                Q = speye(size(g1a));
            end
            if ~(strcmp(options_.simul.preconditioner, 'umfiter') && iter == 0)
                [iter_tol, iter_maxit, gmres_restart] = iter_solver_params(options_, g1a, b);
                if stack_solve_algo == 2
                    [zb, flag] = gmres(P*g1a*Q, P*b, gmres_restart, iter_tol, iter_maxit, L, U);
                else
                    [zb, flag] = bicgstab(P*g1a*Q, P*b, iter_tol, iter_maxit, L, U);
                end
                iter_solver_error_flag(flag)
            end
            dx = Q*zb - ya;
            ya = ya + lambda*dx;
            y(y_index, y_kmin+(1:periods))=reshape(ya',length(y_index),periods);
        elseif stack_solve_algo==4
            stpmx = 100 ;
            stpmax = stpmx*max([sqrt(ya'*ya);size(y_index,2)]);
            nn=1:size(ra,1);
            g = (ra'*g1a)';
            f = 0.5*ra'*ra;
            p = -g1a\ra;
            yn = lnsrch1(ya,f,g,p,stpmax,@lnsrch1_wrapper_two_boundaries,nn,nn, options_.solve_tolx, fh, Block_Num, y, y_index,x, M_.params, steady_state, T, periods, Blck_size, M_);
            dx = ya - yn;
            y(y_index, y_kmin+(1:periods))=reshape(yn',length(y_index),periods);
        end
    end
    iter=iter+1;
    if verbose
        disp(['iteration: ' num2str(iter,'%d') ' error: ' num2str(max_res,'%e')]);
    end
end

if iter > options_.simul.maxit
    if verbose
        printline(41)
        %disp(['No convergence after ' num2str(iter,'%4d') ' iterations in Block ' num2str(Block_Num,'%d')])
    end
    success = false;
    return
end

success = true;


function y3n = dynendo(y, it_, M_)
    y3n = reshape(y(:, it_+(-1:1)), 3*M_.endo_nbr, 1);

function ra = lnsrch1_wrapper_two_boundaries(ya, fh, Block_Num, y, y_index, x, ...
                                             params, steady_state, T, periods, ...
                                             y_size, M_)
    y(y_index, M_.maximum_lag+(1:periods)) = reshape(ya',length(y_index),periods);
    ra = NaN(periods*y_size, 1);
    for it_ = M_.maximum_lag+(1:periods)
        [~, ~, ra((it_-M_.maximum_lag-1)*y_size+(1:y_size)), g1] = fh(dynendo(y, it_, M_), x(it_, :), params, steady_state, M_.block_structure.block(Block_Num).g1_sparse_rowval, M_.block_structure.block(Block_Num).g1_sparse_colval, M_.block_structure.block(Block_Num).g1_sparse_colptr, T(:, it_));
    end
