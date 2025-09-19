function [res, A, info] = ep_problem_0(y, x, pfm)

% Evaluate the residuals and stacked Jacobian of a stochastic perfect
% foresight, considering sequences of future innovations in a perfect n-ary tree.
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

params = pfm.params;
steady_state = pfm.steady_state;
ny = pfm.ny;
periods = pfm.periods;
dynamic_resid = pfm.dynamic_resid;
dynamic_g1 = pfm.dynamic_g1;
sparse_rowval = pfm.sparse_rowval;
sparse_colval = pfm.sparse_colval;
sparse_colptr = pfm.sparse_colptr;
lead_lag_incidence = pfm.lead_lag_incidence;
i_cols_1 = pfm.i_cols_1;
i_cols_j = pfm.i_cols_j;
i_cols_T = pfm.i_cols_T;
order = pfm.order;
hybrid_order = pfm.hybrid_order;
h_correction = pfm.h_correction;
nodes = pfm.nodes;
weights = pfm.weights;
nnodes = pfm.nnodes;

i_cols_p = pfm.i_cols_p;
i_cols_s = pfm.i_cols_s;
i_cols_f = pfm.i_cols_f;
i_rows = pfm.i_rows;

i_cols_Ap = pfm.i_cols_Ap;
i_cols_As = pfm.i_cols_As;
i_cols_Af = pfm.i_cols_Af;
i_hc = pfm.i_hc;

Y = pfm.Y;
Y(pfm.i_upd_y) = y;

A1 = pfm.A1;
res = pfm.res;

for i = 1:order+1
    i_w_p = 1;
    for j = 1:nnodes^(i-1)
        innovation = x;
        if i > 1
            innovation(i+1,:) = nodes(mod(j-1,nnodes)+1,:);
        end
        if i <= order
            for k=1:nnodes
                if hybrid_order && i==order
                    z = [Y(i_cols_p,i_w_p);
                         Y(i_cols_s,j);
                         Y(i_cols_f,(j-1)*nnodes+k)+h_correction(i_hc)];
                else
                    z = [Y(i_cols_p,i_w_p);
                         Y(i_cols_s,j);
                         Y(i_cols_f,(j-1)*nnodes+k)];
                end

                [d1, T_order, T] = dynamic_resid(z, innovation(i+1,:), params, steady_state);
                jacobian = dynamic_g1(z, innovation(i+1,:), params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
                if i == 1
                    % in first period we don't keep track of
                    % predetermined variables
                    i_cols_A = [i_cols_As - ny; i_cols_Af];
                    A1(i_rows,i_cols_A) = A1(i_rows,i_cols_A) + weights(k)*jacobian(:,i_cols_1);
                else
                    i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                    A1(i_rows,i_cols_A) = A1(i_rows,i_cols_A) + weights(k)*jacobian(:,i_cols_j);
                end
                res(:,i,j) = res(:,i,j)+weights(k)*d1;
                i_cols_Af = i_cols_Af + ny;
            end
        else
            z = [Y(i_cols_p,i_w_p);
                 Y(i_cols_s,j);
                 Y(i_cols_f,j)];
            [d1, T_order, T] = dynamic_resid(z, innovation(i+1,:), params, steady_state);
            jacobian = dynamic_g1(z, innovation(i+1,:), params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
            if i == 1
                % in first period we don't keep track of
                % predetermined variables
                i_cols_A = [i_cols_As - ny; i_cols_Af];
                A1(i_rows,i_cols_A) = jacobian(:,i_cols_1);
            else
                i_cols_A = [i_cols_Ap; i_cols_As; i_cols_Af];
                A1(i_rows,i_cols_A) = jacobian(:,i_cols_j);
            end
            res(:,i,j) = d1;
            i_cols_Af = i_cols_Af + ny;
        end
        i_rows = i_rows + ny;
        if mod(j,nnodes) == 0
            i_w_p = i_w_p + 1;
        end
        if i > 1
            if mod(j,nnodes) == 0
                i_cols_Ap = i_cols_Ap + ny;
            end
            i_cols_As = i_cols_As + ny;
        end
    end
    i_cols_p = i_cols_p + ny;
    i_cols_s = i_cols_s + ny;
    i_cols_f = i_cols_f + ny;
end
nzA = cell(periods,pfm.world_nbr);
for j=1:pfm.world_nbr
    i_rows_y = (1:3*ny) + (order+1)*ny;
    offset_c = ny*(sum(nnodes.^(0:order-1),2)+j-1);
    offset_r = (j-1)*ny;
    for i=order+2:periods
        [d1, T_order, T] = dynamic_resid(Y(i_rows_y,j), x(i+1,:), params, steady_state);
        jacobian = dynamic_g1(Y(i_rows_y,j), x(i+1,:), params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T_order, T);
        if i == periods
            [ir,ic,v] = find(jacobian(:,i_cols_T));
        else
            [ir,ic,v] = find(jacobian(:,i_cols_j));
        end
        nzA{i,j} = [offset_r+ir,offset_c+pfm.icA(ic), v]';
        res(:,i,j) = d1;
        i_rows_y = i_rows_y + ny;
        offset_c = offset_c + pfm.world_nbr*ny;
        offset_r = offset_r + pfm.world_nbr*ny;
    end
end
A2 = [nzA{:}]';
A = [A1; sparse(A2(:,1),A2(:,2),A2(:,3),ny*(periods-order-1)*pfm.world_nbr,pfm.dimension)];
res = res(pfm.i_upd_r);
