function [residuals,JJacobian] = perfect_foresight_mcp_problem(y, dynamic_resid_function, ...
                                                               dynamic_g1_function, Y0, YT, ...
                                                               exo_simul, params, steady_state, ...
                                                               dynamic_g1_sparse_rowval, ...
                                                               dynamic_g1_sparse_colval, ...
                                                               dynamic_g1_sparse_colptr, ...
                                                               maximum_lag, T, ny, eq_index)
% function [residuals,JJacobian] = perfect_foresight_mcp_problem(y, dynamic_resid_function, ...
%                                                                dynamic_g1_function, Y0, YT, ...
%                                                                exo_simul, params, steady_state, ...
%                                                                dynamic_g1_sparse_rowval, ...
%                                                                dynamic_g1_sparse_colval, ...
%                                                                dynamic_g1_sparse_colptr, ...
%                                                                maximum_lag, T, ny, eq_index)
% Computes the residuals and the Jacobian matrix for a perfect foresight problem over T periods
% in a mixed complementarity problem context
%
% INPUTS
%   y                        [double] (ny*T)*1 array, path for the endogenous variables
%   dynamic_resid_function   [handle] function handle to dynamic residuals
%   dynamic_g1_function      [handle] function handle to dynamic Jacobian
%   Y0                       [double] ny*1 array, initial conditions for the endogenous variables
%   YT                       [double] ny*1 array, terminal conditions for the endogenous variables
%   exo_simul                [double] (max_lag+nperiods+max_lead)*M_.exo_nbr matrix of exogenous
%                                     variables (in declaration order) for all simulation periods
%   params                   [double] nparams*1 array, parameter values
%   steady_state             [double] ny*1 vector of steady state values
%   dynamic_g1_sparse_rowval [integer vector] eponymous field in M_
%   dynamic_g1_sparse_colval [integer vector] eponymous field in M_
%   dynamic_g1_sparse_colptr [integer vector] eponymous field in M_
%   maximum_lag              [scalar] maximum lag present in the model
%   T                        [scalar] number of simulation periods
%   ny                       [scalar] number of endogenous variables
%   eq_index                 [double] ny*1 array, index vector describing residual mapping resulting
%                                     from complementarity setup
% OUTPUTS
%   residuals                [double] (ny*T)*1 array, residuals of the stacked problem
%   JJacobian                [double] (ny*T)*(ny*T) array, Jacobian of the stacked problem

% Copyright © 1996-2024 Dynare Team
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

YY = reshape([Y0; y; YT], ny, T+2);

residuals = NaN(T*ny,1);
if nargout == 2
    iJacobian = cell(T,1);
end

offset = 0;

for t = 1:T
    y3n = YY(:, t+(0:2));
    x = exo_simul(t+maximum_lag, :);

    [res, TT_order, TT] = dynamic_resid_function(y3n, x, params, steady_state);
    residuals(offset + (1:ny)) = res(eq_index);

    if nargout == 2
        jacobian = dynamic_g1_function(y3n, x, params, steady_state, dynamic_g1_sparse_rowval, dynamic_g1_sparse_colval, dynamic_g1_sparse_colptr, TT_order, TT);
        if T == 1 % static model
            [rows, cols, vals] = find(jacobian(eq_index, ny+(1:ny)));
        elseif t == 1
            [rows, cols, vals] = find(jacobian(eq_index, ny+(1:2*ny)));
        elseif t == T
            [rows, cols, vals] = find(jacobian(eq_index, 1:2*ny));
            cols = cols - ny;
        else
            [rows, cols, vals] = find(jacobian(eq_index, 1:3*ny));
            cols = cols - ny;
        end
        if numel(eq_index) == 1 % find() will return row vectors in this case
            rows = rows';
            cols = cols';
            vals = vals';
        end
        iJacobian{t} = [offset+rows, offset+cols, vals];
    end
    offset = offset + ny;
end

if nargout == 2
    iJacobian = cat(1,iJacobian{:});
    JJacobian = sparse(iJacobian(:,1),iJacobian(:,2),iJacobian(:,3),T*ny,T*ny);
end
