function [residuals,JJacobian] = linear_perfect_foresight_problem(y, dynamicjacobian, Y0, YT, ...
                                                                  exo_simul, params, steady_state, ...
                                                                  maximum_lag, T, ny)

% Computes the residuals and the Jacobian matrix for a linear perfect foresight problem over T periods.

% Copyright © 2015-2024 Dynare Team
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
    z = [ vec(YY(:, t+(0:2))); transpose(exo_simul(t+maximum_lag, :)) ];
    residuals(offset + (1:ny)) = dynamicjacobian*z;

    if nargout == 2
        if T == 1 % static model
            [rows, cols, vals] = find(dynamicjacobian(:, ny+(1:ny)));
        elseif t == 1
            [rows, cols, vals] = find(dynamicjacobian(:, ny+(1:2*ny)));
        elseif t == T
            [rows, cols, vals] = find(dynamicjacobian(:, 1:2*ny));
            cols = cols - ny;
        else
            [rows, cols, vals] = find(dynamicjacobian(:, 1:3*ny));
            cols = cols - ny;
        end
        if size(dynamicjacobian, 1) == 1 % find() will return row vectors in this case
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
    JJacobian = sparse(iJacobian(:,1), iJacobian(:,2), iJacobian(:,3), T*ny, T*ny);
end
