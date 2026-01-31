function [xparam1, icheck] = set_init_state(xparam1, ys0, options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [xparam1, icheck] = set_init_state(xparam1, options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% Projects endogenous initial-state parameters to satisfy Pstar null-space constraints.
%
% This routine temporarily disables endogenous-initial-state priors and
% computes the initial covariance `Pstar` (with `lik_init=1`). If `Pstar`
% is singular, it enforces compatibility of the initial state vector with
% the null space implied by `Pstar` by adjusting the corresponding entries
% in `xparam1` (those named 'init ...'). It returns `icheck=true` only if
% `xparam1` was modified to meet these constraints.
%
% INPUTS
% - xparam1             [double]        parameter vector (includes 'init ' state parameters).
% - ys0                 [double]        original steady state
% - options_            [structure]     options; internally sets `lik_init=1` for checks.
% - M_                  [structure]     model structure (updates `endo_initial_state`).
% - estim_params_       [structure]     parameters to be estimated.
% - bayestopt_          [structure]     priors and state/measurement indices (uses mf0).
% - BoundsInfo          [structure]     prior bounds info.
% - dr                  [structure]     reduced-form decision rules.
% - endo_steady_state   [vector]        steady state for endogenous variables.
% - exo_steady_state    [vector]        steady state for exogenous variables.
% - exo_det_steady_state [vector]       steady state for deterministic exogenous variables.
%
% OUTPUTS
% - xparam1             [double]        possibly updated parameter vector with consistent
%                                       initial-state entries.
% - icheck              [logical]       true if `xparam1` was modified, false otherwise.
%
% SEE ALSO
%   get_init_state_prior - computes the Pstar-based subspace (UP, XP) used here.

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

icheck=false;

M_.endo_initial_state.status=false;
options_.lik_init=1;
options_.estimate_initial_states_endogenous_prior=false;
M_ = set_all_parameters(xparam1,estim_params_,M_);
[Pstar, Q, info] = get_pstar(xparam1,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
if info(1)
    return %JP: error code is neither returned nor handled
end
dr.ys = evaluate_steady_state(endo_steady_state,[exo_steady_state; exo_det_steady_state],M_,options_,true);
M_.endo_initial_state.status=true;

[UP,XP] = svd(0.5*(Pstar(bayestopt_.mf0,bayestopt_.mf0)+Pstar(bayestopt_.mf0,bayestopt_.mf0)'));
% handle zero singular values
isp = find(diag(XP)>options_.kalman_tol);
if length(isp)>rank(Q)
    isp = isp(1:rank(Q));
end
isn = length(isp)+1:length(XP); 

UPN = UP(:,isn);
[~,im]=max(abs(UPN));
if length(unique(im))<length(im)
    % make sure we identify as many states as the null singular values
    in=im*nan;
    for k=1:length(isn)
        [~, tmp] = max(abs(UPN(:,k)));
        if not(ismember(tmp,in))
            in(k) = tmp;
        else
            [~, tmp] = sort(abs(UPN(:,k)));
            not_found = true;
            offset=0;
            while not_found
                offset=offset+1;
                if not(ismember(tmp(end-offset),in))
                    in(k) = tmp(end-offset);
                    not_found = false;
                end
            end
        end
    end
    im=in;
end

a = get_init_state(zeros(M_.endo_nbr,1),xparam1,estim_params_,dr,M_,options_);
if options_.loglinear
    a = a + log(dr.ys(dr.order_var))-log(ys0(dr.order_var));
else
    a = a + dr.ys(dr.order_var)-ys0(dr.order_var);
end

mycheck =max(abs(UPN'*a(dr.restrict_var_list(bayestopt_.mf0))));
if mycheck>options_.kalman_tol
    icheck=true;
    nstates = length(dr.state_var);
    % check that state draws are in the null space of prior states 0|0
    in=setdiff(1:nstates,im);
    a1=a;
    % express [im] states as function of [in] states according to null
    % space structure
    a1(dr.restrict_var_list(bayestopt_.mf0(im)))=-UPN(im,:)'\(UPN(in,:)'*a(dr.restrict_var_list(bayestopt_.mf0(in))));
    alphahat01 = a1(dr.inv_order_var); %declaration order

    % map params associated to init states

    IB = startsWith(bayestopt_.name, 'init ');
    if options_.loglinear
        xparam1(IB) = exp(alphahat01(dr.state_var)).*dr.ys(dr.state_var);
    else
        xparam1(IB) = alphahat01(dr.state_var)+dr.ys(dr.state_var);
    end
end
