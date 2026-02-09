function [dr, state_var]=set_state_space(dr,M_)
% [dr, state_var]=set_state_space(dr,M_)
% Computes the DR ordering and inverse ordering.
%
% INPUTS
% - dr            [struct]    MATLAB's structure describing decision and transition rules.
% - M_            [struct]    MATLAB's structure describing the model.
%
% OUTPUTS
% - dr            [struct]    MATLAB's structure describing decision and transition rules.
% - state_var     [vector]    (optional) Indices of state variables (variables with lags) in
%                             declaration order. See also dyn_first_order_solver.m where state_var
%                             is computed as [no_both_lag_id, both_id].
%
% CALLED BY
%   check, cli/prior, discretionary_policy/discretionary_policy_1, ep/extended_path_initialization,
%   estimation/dynare_estimation_init, estimation/optimize_prior, estimation/prior_sampler,
%   perfect-foresight-models, stochastic_solver/stochastic_solvers,
%   stochastic_solver/stoch_simul, +osr/run, preprocessor-generated driver.m files

% Copyright © 1996-2026 Dynare Team
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

fwrd_var = find(M_.lead_lag_incidence(M_.maximum_endo_lag+2:end,:))';
if M_.maximum_endo_lag > 0
    pred_var = find(M_.lead_lag_incidence(1,:))';
    both_var = intersect(pred_var,fwrd_var);
    pred_var = setdiff(pred_var,both_var);
    fwrd_var = setdiff(fwrd_var,both_var);
    static_var = setdiff([1:M_.endo_nbr]',union(union(pred_var,both_var),fwrd_var));  % static variables
else
    pred_var = [];
    both_var = [];
    static_var = setdiff([1:M_.endo_nbr]',fwrd_var);
end

dr.order_var = [static_var(:); pred_var(:); both_var(:); fwrd_var(:)];
dr.inv_order_var(dr.order_var) = 1:M_.endo_nbr;

if nargout > 1
    % State variables: variables that appear at t-1 (predetermined)
    % See also: dyn_first_order_solver.m where state_var is computed as [no_both_lag_id, both_id]
    state_var = [pred_var(:); both_var(:)]';
end
