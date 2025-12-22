function a = get_init_state(a,xparam1,estim_params_,dr,M_,options_)
% a = get_init_state(a,xparam1,estim_params_,dr,M_,options_)
% Builds the Kalman/filter initial state vector from model values.
%
% Given `M_.endo_initial_state.values` and the steady state `dr.ys`, this
% routine fills the initial state vector `a` (decision-rule order) with
% deviations from steady state. In log-linear mode it uses log deviations;
% otherwise level deviations. The function updates parameters from
% `xparam1` before computing the state.
%
% INPUTS
% - a                   [double]     initial state vector placeholder (size M_.endo_nbr).
% - xparam1             [double]     parameter vector used to set `M_.params`.
% - estim_params_       [structure]  parameters-to-estimate metadata.
% - dr                  [structure]  decision rule structure; uses `dr.ys` and `dr.order_var`.
% - M_                  [structure]  model structure; uses `endo_initial_state.values` and `state_var`.
% - options_            [structure]  options; uses `loglinear` and `logged_steady_state`.
%
% OUTPUTS
% - a                   [double]     updated initial state vector in decision-rule order.
%
% NOTES
% - In log-linear mode (`options_.loglinear==true`), computes
%   log(M_.endo_initial_state.values) - log(dr.ys) on state entries.
% - In level mode, computes M_.endo_initial_state.values - dr.ys on state entries.
% - Errors if `options_.logged_steady_state==true` (unsupported combination).

% Copyright © 2024-2026 Dynare Team
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

M_ = set_all_parameters(xparam1,estim_params_,M_);

if options_.loglinear && ~options_.logged_steady_state
    a(M_.state_var) = log(M_.endo_initial_state.values(M_.state_var)) - log(dr.ys(M_.state_var));
elseif ~options_.loglinear && ~options_.logged_steady_state
    a(M_.state_var) = M_.endo_initial_state.values(M_.state_var) - dr.ys(M_.state_var);
else
    error('The steady state is logged. This should not happen. Please contact the developers')
end

a=a(dr.order_var);
