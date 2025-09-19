function [res, A, info] = ep_problem_1(y, x, pfm)

% Evaluate the residuals and stacked Jacobian of a stochastic perfect
% foresight, considering sequences of future innovations in a sparse tree.
%
% INPUTS:
% - y      [double]   m×1 vector (endogenous variables in all periods and future worlds).
% - x      [double]   q×1 vector of exogenous variables.
% - pfm    [struct]   Definition of the perfect foresight model to be solved.
%
% OUTPUTS:
% - res    [double]   m×1 vector, residuals of the stacked equations.
% - A      [double]   m×m sparse matrix, Jacobian of the stacked equations.
% - info   [logical]  scalar
%
% REMARKS:
% [1] The structure pfm holds the given initial condition for the states (pfm.y0) and the terminal condition

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

info = false;
A = [];

dynamic_resid = pfm.dynamic_resid;
dynamic_g1 = pfm.dynamic_g1;
sparse_rowval = pfm.sparse_rowval;
sparse_colval = pfm.sparse_colval;
sparse_colptr = pfm.sparse_colptr;
ny = pfm.ny;
params = pfm.params;
steady_state = pfm.steady_state;
order = pfm.stochastic_order;
nodes = pfm.nodes;
nnodes = pfm.nnodes;
weights = pfm.weights;
h_correction = pfm.h_correction;
dimension = pfm.dimension;
world_nbr = pfm.world_nbr;
nnzA = pfm.nnzA;
periods = pfm.periods;
i_rows = pfm.i_rows';
i_cols = pfm.i_cols;
nyp = pfm.nyp;
nyf = pfm.nyf;
hybrid_order = pfm.hybrid_order;
i_cols_1 = pfm.i_cols_1;
i_cols_j = pfm.i_cols_j;
icA = pfm.icA;
i_cols_T = pfm.i_cols_T;
eq_index = pfm.eq_index;



i_cols_p = 1:ny;
i_cols_s = ny + (1:ny);
i_cols_f = 2*ny + (1:ny);
i_cols_Ap0 = i_cols(1:nyp);
i_cols_As = i_cols(nyp+(1:ny));
i_cols_Af0 = i_cols(nyp+ny+(1:nyf)) - ny;
i_hc = 1:ny;

nzA = cell(periods,world_nbr);
res = zeros(ny,periods,world_nbr);
Y = zeros(ny*(periods+2),world_nbr);
Y(1:ny,1) = pfm.y0;
Y(end-ny+1:end,:) = repmat(steady_state,1,world_nbr);
Y(pfm.i_upd_y) = y;

offset_r0 = 0;
for i = 1:order+1
    i_w_p = 1;
    for j = 1:(1+(nnodes-1)*(i-1))
        innovation = x;
        if i <= order && j == 1
            % first world, integrating future shocks
            if nargout > 1
                A1 = sparse([],[],[],i*(1+(nnodes-1)*(i-1))*ny,dimension,nnzA*world_nbr);
            end
            for k=1:nnodes
                if nargout > 1
                    if i == 2
                        i_cols_Ap = i_cols_Ap0;
                    elseif i > 2
                        i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes- ...
                                                          1)*(i-2)*(i-3)/2);
                    end
                    if k == 1
                        i_cols_Af = i_cols_Af0 + ny*(i-1+(nnodes-1)*i*(i-1)/ ...
                                                     2);
                    else
                        i_cols_Af = i_cols_Af0 + ny*(i-1+(nnodes-1)*i*(i-1)/ ...
                                                     2+(i-1)*(nnodes-1) ...
                                                     + k - 1);
                    end
                end
                if i > 1
                    innovation(i+1,:) = nodes(k,:);
                end
                if k == 1
                    k1 = 1;
                else
                    k1 = (nnodes-1)*(i-1)+k;
                end
                if hybrid_order && (k > 1 || i == order)
                    z = [Y(i_cols_p,1);
                         Y(i_cols_s,1);
                             Y(i_cols_f,k1)+h_correction(i_hc)];
                else
                    z = [Y(i_cols_p,1);
                         Y(i_cols_s,1);
                         Y(i_cols_f,k1)];
                end
                if nargout > 1
                    [d1, T_order, T] = dynamic_resid(z, innovation(i+1,:), params, steady_state);
                    jacobian = dynamic_g1(z, innovation(i+1,:), params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
                    if i == 1
                        % in first period we don't keep track of
                        % predetermined variables
                        i_cols_A = [i_cols_As - ny; i_cols_Af];
                        A1(i_rows,i_cols_A) = A1(i_rows,i_cols_A) + weights(k)*jacobian(eq_index,i_cols_1);
                    else
                        i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                        A1(i_rows,i_cols_A) = A1(i_rows,i_cols_A) + weights(k)*jacobian(eq_index,i_cols_j);
                    end
                else
                    d1 = dynamic_resid(z, innovation(i+1,:), params, steady_state);
                end
                res(:,i,1) = res(:,i,1)+weights(k)*d1(eq_index);
            end
            if nargout > 1
                [ir,ic,v] = find(A1);
                nzA{i,j} = [ir,ic,v]';
            end
        elseif j > 1 + (nnodes-1)*(i-2)
            % new world, using previous state of world 1 and hit
            % by shocks from nodes
            if nargout > 1
                i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes-1)*(i-2)*(i-3)/2);
                i_cols_Af = i_cols_Af0 + ny*(i+(nnodes-1)*i*(i-1)/2+j-2);
            end
            k = j - (nnodes-1)*(i-2);
            innovation(i+1,:) = nodes(k,:);
            z = [Y(i_cols_p,1);
                 Y(i_cols_s,j);
                 Y(i_cols_f,j)];
            if nargout > 1
                [d1, T_order, T] = dynamic_resid(z, innovation(i+1,:), params, steady_state);
                jacobian = dynamic_g1(z, innovation(i+1,:), params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
                i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                [ir,ic,v] = find(jacobian(eq_index,i_cols_j));
                nzA{i,j} = [i_rows(ir),i_cols_A(ic), v]';
            else
                d1 = dynamic_resid(z, innovation(i+1,:), params, steady_state);
            end
            res(:,i,j) = d1(eq_index);
            if nargout > 1
                i_cols_Af = i_cols_Af + ny;
            end
        else
            % existing worlds other than 1
            if nargout > 1
                i_cols_Ap = i_cols_Ap0 + ny*(i-2+(nnodes-1)*(i-2)*(i-3)/2+j-1);
                i_cols_Af = i_cols_Af0 + ny*(i+(nnodes-1)*i*(i-1)/2+j-2);
            end
            z = [Y(i_cols_p,j);
                 Y(i_cols_s,j);
                 Y(i_cols_f,j)];
            if nargout > 1
                [d1, T_order, T] = dynamic_resid(z, innovation(i+1,:), params, steady_state);
                jacobian = dynamic_g1(z, innovation(i+1,:), params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
                i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                [ir,ic,v] = find(jacobian(eq_index,i_cols_j));
                nzA{i,j} = [i_rows(ir),i_cols_A(ic),v]';
                i_cols_Af = i_cols_Af + ny;
            else
                d1 = dynamic_resid(z, innovation(i+1,:), params, steady_state);
            end
            res(:,i,j) = d1(eq_index);
        end
        i_rows = i_rows + ny;
        if mod(j,nnodes) == 0
            i_w_p = i_w_p + 1;
        end
        if nargout > 1 && i > 1
            i_cols_As = i_cols_As + ny;
        end
        offset_r0 = offset_r0 + ny;
    end
    i_cols_p = i_cols_p + ny;
    i_cols_s = i_cols_s + ny;
    i_cols_f = i_cols_f + ny;
end
for j=1:world_nbr
    i_rows_y = (1:3*ny) + (order+1)*ny;
    offset_c = ny*(order+(nnodes-1)*(order-1)*order/2+j-1);
    offset_r = offset_r0+(j-1)*ny;
    for i=order+2:periods
        if nargout > 1
            [d1, T_order, T] = dynamic_resid(Y(i_rows_y,j), x(i+1,:), params, steady_state);
            jacobian = dynamic_g1(Y(i_rows_y,j), x(i+1,:), params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
            if i < periods
                [ir,ic,v] = find(jacobian(eq_index,i_cols_j));
            else
                [ir,ic,v] = find(jacobian(eq_index,i_cols_T));
            end
            nzA{i,j} = [offset_r+ir,offset_c+icA(ic), v]';
        else
            d1 = dynamic_resid(Y(i_rows_y,j), x(i+1,:), params, steady_state);
        end
        res(:,i,j) = d1(eq_index);
        i_rows_y = i_rows_y + ny;
        offset_c = offset_c + world_nbr*ny;
        offset_r = offset_r + world_nbr*ny;
    end
end
if nargout > 1
    iA = [nzA{:}]';
    A = sparse(iA(:,1),iA(:,2),iA(:,3),dimension,dimension);
end
res = res(pfm.i_upd_r);
