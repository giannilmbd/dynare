function [LIK, likk, a] = kalman_filter_ss(Y,start,last,a,T,K,iF,log_dF,Z,pp,Zflag,analytic_derivation,Da,DT,DYss,analytic_Hessian,DK,DF,D2a,D2T,D2Yss,D2K,D2F)
% Computes the log-likelihood of a stationary state space model (steady-state Kalman filter).

%
% INPUTS
% - Y                       [matrix]      [pp x T] matrix of observed data
% - start                   [integer]     index of the first period processed in Y
% - last                    [integer]     index of the last period processed in Y
% - a                       [vector]      [mm x 1] predicted initial state vector, E_0(alpha_1)
% - T                       [matrix]      [mm x mm] transition matrix of the state equation
% - K                       [matrix]      [mm x pp] steady-state Kalman gain
% - iF                      [matrix]      [pp x pp] inverse of steady-state covariance matrix of prediction errors
% - log_dF                  [double]      log determinant of steady-state covariance matrix of prediction errors
% - Z                       [matrix]      [pp x mm] measurement matrix, or index vector when Zflag=0
% - pp                      [integer]     number of observed variables
% - Zflag                   [integer]     0 if Z is an index vector; 1 if Z is a [pp x mm] matrix
% - analytic_derivation     [logical]     whether to compute analytic derivatives (true/false)
% - Da, DT, DYss            [array]       first-derivative objects used when analytic_derivation is true
% - analytic_Hessian        [string]      Hessian mode: '' (none), 'full' (analytic), 'opg' (outer product), 'asymptotic'
% - DK                      [array]       [mm x pp x k] cached Kalman gain derivatives
% - DF                      [array]       [pp x pp x k] cached forecast error variance derivatives
% - D2a, D2T, D2Yss         [array]       second-derivative objects used when analytic_Hessian=='full'
% - D2K                     [array]       [mm x pp x k x k] cached second Kalman gain derivatives (full Hessian only)
% - D2F                     [array]       [pp x pp x k x k] cached second forecast error variance derivatives (full Hessian only)
%
% OUTPUTS
% - LIK                     [double|cell] minus log-likelihood; if analytic_derivation is true, returns cell array {LIK,DLIK[,Hess]}
% - likk                    [vector|cell] [smpl x 1] period-wise log-likelihood contributions; if analytic_derivation is true, returns {likk,dlikk}
% - a                       [vector]      [mm x 1] filtered state estimate for next period, E_T(alpha_{T+1})
%
% This function is called by: kalman_filter
% This function calls: none

% Copyright © 2011-2026 Dynare Team
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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
t    = start;              % Initialization of the time index.
likk = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
notsteady = 0;
if nargin<12
    analytic_derivation = false;
end
if nargin<16
    analytic_Hessian = '';
end

% Extract Hessian mode flags from string
full_Hess = strcmp(analytic_Hessian, 'full');
asy_hess = strcmp(analytic_Hessian, 'asymptotic');

if ~analytic_derivation
    DLIK=[];
    Hess=[];
else
    k = size(DT,3);                                 % number of structural parameters
    DLIK  = zeros(k,1);                             % Initialization of the score.
    dlikk = zeros(smpl,k);
    if full_Hess || asy_hess
        Hess  = zeros(k,k);                             % Initialization of the Hessian
    else
        Hess=[];
    end
end

while t <= last
    if Zflag
        v = Y(:,t)-Z*a;
    else
        v = Y(:,t)-a(Z);
    end
    tmp = (a+K*v);
    if analytic_derivation
        if full_Hess
            [Da,~,DLIKt,Hesst,DK,DF,D2a,~,D2K,D2F] = computeDLIK(k,tmp,Z,Zflag,v,T,K,[],iF,Da,DYss,DT,[],[],[],notsteady,true,DK,DF,D2a,D2Yss,D2T,[],[],[],D2K,D2F);
        else
            [Da,~,DLIKt,Hesst,DK,DF] = computeDLIK(k,tmp,Z,Zflag,v,T,K,[],iF,Da,DYss,DT,[],[],[],notsteady,false,DK,DF);
        end
        DLIK = DLIK + DLIKt;
        if full_Hess || asy_hess
            Hess = Hess + Hesst;
        end
        dlikk(t-start+1,:)=DLIKt;
    end
    a = T*tmp;
    likk(t-start+1) = transpose(v)*iF*v;
    t = t+1;
end

% Adding constant determinant of F (prediction error covariance matrix)
likk = likk + log_dF;

% Add log-likelihood constants and divide by two
likk = .5*(likk + pp*log(2*pi));

% Sum the observation's densities (minus the likelihood)
LIK = sum(likk);
if analytic_derivation
    dlikk = dlikk/2;
    DLIK = DLIK/2;
    likk = {likk, dlikk};
end
if full_Hess || asy_hess
    if ~asy_hess
        Hess = Hess + tril(Hess,-1)';
    end
    Hess = -Hess/2;
    LIK={LIK,DLIK,Hess};
elseif analytic_derivation
    LIK={LIK,DLIK};
end
