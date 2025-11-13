function [UP, XP] = get_init_state_prior(xparam1, options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [UP, XP] = get_init_state_prior(xparam1, options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% Computes the endogenous log prior addition to the initial prior
%
% INPUTS
%    xparam1            [double]     n vector of estimated params
%    bayestopt_         [structure]
%    dr                 [structure]
%    M_                 [structure]
%
% OUTPUTS
%    xparam1            [double]     n vector of estimated params
%    icheck             [logical]    flag for the need to update xparam1

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

filter_initial_state=M_.filter_initial_state;
M_.filter_initial_state=[];
options_.lik_init=1;
options_.init_state_endogenous_prior=false;
[Pstar, info] = get_pstar(xparam1,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
if info(1)
    return
end
M_.filter_initial_state=filter_initial_state;

[UP,XP] = svd(0.5*(Pstar(bayestopt_.mf0,bayestopt_.mf0)+Pstar(bayestopt_.mf0,bayestopt_.mf0)'));
isp = find(diag(XP)>options_.kalman_tol);
UP = UP(:,isp);
XP = XP(isp,isp);
