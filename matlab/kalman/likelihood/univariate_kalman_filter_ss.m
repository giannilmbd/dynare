function [LIK,likk,a] = univariate_kalman_filter_ss(Y,start,last,a,P,kalman_tol,T,H,Z,pp,Zflag,analytic_derivation,Da,DT,DYss,DP,DH,analytic_Hessian,D2a,D2T,D2Yss,D2H,D2P)
% Computes the log-likelihood of a stationary state space model (steady-state univariate Kalman filter).

%
% INPUTS
% - Y                       [matrix]      [pp x T] matrix of observed data
% - start                   [integer]     index of the first period processed in Y
% - last                    [integer]     index of the last period processed in Y
% - a                       [vector]      [mm x 1] predicted initial state mean, E_0(alpha_1)
% - P                       [matrix]      [mm x mm] steady-state covariance matrix of the state vector
% - kalman_tol              [double]      tolerance parameter for scalar forecast-error variances
% - T                       [matrix]      [mm x mm] transition matrix of the state equation
% - H                       [vector]      [pp x 1] diagonal of covariance matrix of measurement errors
% - Z                       [matrix]      [pp x mm] measurement matrix, or index vector when Zflag=0
% - pp                      [integer]     number of observed variables
% - Zflag                   [integer]     0 if Z is an index vector; 1 if Z is a [pp x mm] matrix
% - analytic_derivation     [logical]     whether to compute analytic derivatives (true/false)
% - Da, DT, DYss, DP, DH    [array]       first-derivative objects used when analytic_derivation is true
% - analytic_Hessian        [string]      Hessian mode: '' (none), 'full' (analytic), 'opg' (outer product), 'asymptotic'
% - D2a, D2T, D2Yss, D2H, D2P    [array]       second-derivative objects used when analytic_Hessian=='full'
%
% OUTPUTS
% - LIK                     [double|cell] minus log-likelihood; if analytic_derivation is true, returns cell array {LIK,DLIK[,Hess]}
% - likk                    [matrix|cell] [smpl x pp] period/observable log-likelihood contributions; if analytic_derivation is true, returns {likk,dlikk}
% - a                       [vector]      [mm x 1] filtered state mean at the end of the processed sample
%
% This function is called by: univariate_kalman_filter
% This function calls: none
%
% Algorithm: See univariate_kalman_filter.m

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
likk = zeros(smpl,pp);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
l2pi = log(2*pi);

if nargin<12
    analytic_derivation = false;
end
if nargin<18
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

% Steady state kalman filter.
while t<=last
    s  = t-start+1;
    PP = P;
    if analytic_derivation
        DPP = DP;
        if full_Hess
            D2PP = D2P;
        end
    end
    for i=1:pp
        if Zflag
            prediction_error = Y(i,t) - Z(i,:)*a;
            PPZ = PP*Z(i,:)';
            Fi = Z(i,:)*PPZ + H(i);
        else
            prediction_error = Y(i,t) - a(Z(i));
            PPZ = PP(:,Z(i));
            Fi = PPZ(Z(i)) + H(i);
        end
        if Fi>kalman_tol
            Ki = PPZ/Fi;
            a  = a + Ki*prediction_error;
            PP = PP - PPZ*Ki';
            likk(s,i) = log(Fi) + prediction_error*prediction_error/Fi + l2pi;
            if analytic_derivation
                if full_Hess
                    [Da,DPP,DLIKt,D2a,D2PP, Hesst] = univariate_computeDLIK(k,i,Z(i,:),Zflag,prediction_error,Ki,PPZ,Fi,Da,DYss,DPP,DH(i,:),0,D2a,D2Yss,squeeze(D2H(i,:,:)),D2PP);
                else
                    [Da,DPP,DLIKt,Hesst] = univariate_computeDLIK(k,i,Z(i,:),Zflag,prediction_error,Ki,PPZ,Fi,Da,DYss,DPP,DH(i,:),0);
                end
                DLIK = DLIK + DLIKt;
                if full_Hess || asy_hess
                    Hess = Hess + Hesst;
                end
                dlikk(s,:)=dlikk(s,:)+DLIKt';
            end
        else
            % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
            % p. 157, DK (2012)
        end
    end
    if analytic_derivation
        if full_Hess
            [Da,~,D2a] = univariate_computeDstate(k,a,P,T,Da,DP,DT,[],0,D2a,D2P,D2T);
        else
            Da = univariate_computeDstate(k,a,P,T,Da,DP,DT,[],0);
        end
    end
    a = T*a;
    t = t+1;
end

likk = .5*likk;

LIK = sum(sum(likk));
if analytic_derivation
    dlikk = dlikk/2;
    DLIK = DLIK/2;
    likk = {likk, dlikk};
end
if full_Hess || asy_hess
    %     Hess = (Hess + Hess')/2;
    Hess = -Hess/2;
    LIK={LIK,DLIK,Hess};
elseif analytic_derivation
    LIK={LIK,DLIK};
end