function [lnpriorendoinitstate, lnpriorinitstate] = init_state_endogenous_prior(a_0_given_tm0,T,R,Q,xparam1,bayestopt_,options_)
% [lnpriorendoinitstate, lnpriorinitstate] = init_state_endogenous_prior(a_0_given_tm0,T,R,Q,xparam1,bayestopt_,options_)
% Computes log-prior densities for endogenous initial states.
%
% This routine evaluates the endogenous prior for initial states by
% comparing `a_0_given_tm0` against the model-implied unconditional
% covariance `Pstar` (computed via Lyapunov solution). It returns the
% Gaussian log-density under the endogenous prior, and optionally the
% standard user-declared log-prior for the 'init ' parameters.
%
% INPUTS
% - a_0_given_tm0       [double]     initial state vector at time t=0|t=-1.
% - T                   [double]     state transition matrix from decision rules.
% - R                   [double]     shock impact matrix from decision rules.
% - Q                   [double]     shock covariance matrix.
% - xparam1             [double]     parameter vector (includes 'init ' entries).
% - bayestopt_          [structure]  priors and state indexing (uses mf0, pshape, etc.).
% - M_                  [structure]  model structure (not actively used here).
% - options_            [structure]  options; uses kalman_tol for SVD threshold.
%
% OUTPUTS
% - lnpriorendoinitstate[double]     log-density under the endogenous (Pstar-based) prior.
% - lnpriorinitstate    [double]     (optional) log-density under user-declared priors for 'init ' parameters.
%
% NOTES
% - The endogenous prior is a multivariate normal with covariance Pstar
%   restricted to state variables (mf0), retaining only directions above
%   `options_.kalman_tol` via SVD.
% - If second output is requested, extracts 'init ' parameters and computes
%   their prior density using the declared prior shapes.
%
% SEE ALSO
%   get_init_state_prior - derives the Pstar subspace (UP, XP) used similarly.

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

Pstar=lyapunov_solver(T,R,Q,options_);

[UP,XP] = svd(0.5*(Pstar(bayestopt_.mf0,bayestopt_.mf0)+Pstar(bayestopt_.mf0,bayestopt_.mf0)'));
isp = find(diag(XP)>options_.kalman_tol);
dd=diag(XP(isp,isp));
isp = isp(1:sum((dd./dd(1))>eps));
ns = size(XP,1);
isn=length(isp)+1:ns;
UPN = UP(:,isn);
mycheck =max(abs(UPN'*a_0_given_tm0(bayestopt_.mf0)));
UP = UP(:,isp);
S = XP(isp,isp);
log_dS = log(det(S));

vv = UP'*(a_0_given_tm0(bayestopt_.mf0));
lnpriorendoinitstate = -(log_dS + transpose(vv)/S*vv + ns*log(2*pi))/2;
if mycheck>options_.kalman_tol
    do_nothing=true;
    % lnpriorendoinitstate=-inf;
end
if nargout>1
    % now I remove original state prior declared
    IB = startsWith(bayestopt_.name, 'init ');
    lnpriorinitstate = priordens(xparam1(IB),bayestopt_.pshape(IB),bayestopt_.p6(IB),bayestopt_.p7(IB),bayestopt_.p3(IB),bayestopt_.p4(IB));
end


