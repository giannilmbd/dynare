function maxerror = compute_maxerror(endo_simul, exo_simul, steady_state, M_, options_)
% Computes ∞-norm of residuals for a given path of endogenous,
% given the exogenous path, steady state and parameters in M_

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

periods = size(endo_simul, 2) - M_.maximum_lag - M_.maximum_lead;

if options_.bytecode
    residuals = bytecode('dynamic', 'evaluate', M_, options_, endo_simul, exo_simul, M_.params, steady_state, periods);
else
    ny = size(endo_simul, 1);
    if M_.maximum_lag > 0
        y0 = endo_simul(:, M_.maximum_lag);
    else
        y0 = NaN(ny, 1);
    end
    if M_.maximum_lead > 0
        yT = endo_simul(:, M_.maximum_lag+periods+1);
    else
        yT = NaN(ny, 1);
    end
    yy = endo_simul(:,M_.maximum_lag+(1:periods));
    residuals = perfect_foresight_problem(yy(:), y0, yT, exo_simul, M_.params, steady_state, periods, M_, options_);
end

maxerror = norm(vec(residuals), 'Inf'); % Do not use max(max(abs(…))) because it omits NaN
