function [lnpriorendoinitstate, lnpriorinitstate] = init_state_endogenous_prior(a_0_given_tm0,T,R,Q,xparam1,bayestopt_,M_,options_)
% [lnpriorendoinitstate, lnpriorinitstate] = init_state_endogenous_prior(a_0_given_tm0,T,R,Q,xparam1,bayestopt_,M_,options_)
% Computes the endogenous log prior addition to the initial prior
%
% INPUTS
%    a_0_given_tm0      [double]     k vector of initial states
%    T                  [double]     k*k state matrix
%    R                  [double]     k*n impact shock matric
%    Q                  [double]     k*n impact shock matric
%    xparam1            [double]     n vector of estimated params
%    bayestopt_         [structure]  describing the priors
%    M_                 [structure]  Matlab's structure describing the model
%    options_           [structure]  Matlab's structure describing the options
%
% OUTPUTS
%    lnpriorendoinitstate [double]     scalar of log init state endogenous prior value
%    lnpriorinitstate    [double]     scalar of log init state prior value

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

Pstar=lyapunov_solver(T,R,Q,options_);

[UP,XP] = svd(0.5*(Pstar(bayestopt_.mf0,bayestopt_.mf0)+Pstar(bayestopt_.mf0,bayestopt_.mf0)'));
isp = find(diag(XP)>options_.kalman_tol);
ns = size(XP,1);
UP = UP(:,isp);
S = XP(isp,isp);
log_dS = log(det(S));

vv = UP'*(a_0_given_tm0(bayestopt_.mf0));
lnpriorendoinitstate = -(log_dS + transpose(vv)/S*vv + ns*log(2*pi))/2;

if nargout>1
    % now I remove original state prior declared
    pvec=[];
    for ii=1:size(M_.filter_initial_state,1)
        if ~isempty(M_.filter_initial_state{ii,1})
            tmp1 = strrep(M_.filter_initial_state{ii,2},');','');
            tmp1 = strrep(tmp1,'M_.params(','');
            pvec = [pvec eval(tmp1)];
        end
    end
    [~,~,IB] = intersect(M_.param_names(pvec),bayestopt_.name,'stable');

    lnpriorinitstate = priordens(xparam1(IB),bayestopt_.pshape(IB),bayestopt_.p6(IB),bayestopt_.p7(IB),bayestopt_.p3(IB),bayestopt_.p4(IB));
end


