function [UP, XP] = get_init_state_prior(xparam1, options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [UP, XP] = get_init_state_prior(xparam1, options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% Derives prior constraints for endogenous initial states from Pstar.
%
% This routine temporarily disables endogenous-initial-state priors and
% computes the initial covariance `Pstar` (with `lik_init=1`). It then
% performs an SVD of the restricted `Pstar` block over states to obtain
% eigenvectors/eigenvalues above `options_.kalman_tol`. The returned
% matrices `UP` and `XP` describe the subspace used to enforce prior
% consistency of initial states (e.g., null/positive variance directions).
%
% INPUTS
% - xparam1             [double]     parameter vector used to set model parameters.
% - options_            [structure]  options; internally sets `lik_init=1` and disables init-state prior.
% - M_                  [structure]  model structure; uses `endo_initial_state.status`.
% - estim_params_       [structure]  parameters to be estimated.
% - bayestopt_          [structure]  priors/state indexing; uses `mf0`.
% - BoundsInfo          [structure]  prior bounds and definiteness checks.
% - dr                  [structure]  decision rules; used for model resolution.
% - endo_steady_state   [vector]     steady state for endogenous variables.
% - exo_steady_state    [vector]     steady state for exogenous variables.
% - exo_det_steady_state[vector]     steady state for deterministic exogenous variables.
%
% OUTPUTS
% - UP                  [double]     eigenvectors (columns) of the restricted `Pstar` subspace retained.
% - XP                  [double]     diagonal matrix of retained singular values/eigenvalues.

% Copyright © 2024-2025 Dynare Team
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

M_.endo_initial_state.status=false;
options_.lik_init=1;
options_.init_state_endogenous_prior=false;
[Pstar, info] = get_pstar(xparam1,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
if info(1)
    return
end

[UP,XP] = svd(0.5*(Pstar(bayestopt_.mf0,bayestopt_.mf0)+Pstar(bayestopt_.mf0,bayestopt_.mf0)'));
isp = find(diag(XP)>options_.kalman_tol);
UP = UP(:,isp);
XP = XP(isp,isp);
