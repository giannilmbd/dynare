function [LIK, LIKK, a, P] = kalman_filter(Y,start,last,a,P,kalman_tol,riccati_tol,rescale_prediction_error_covariance,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods,analytic_derivation,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P)
% [LIK, LIKK, a, P] = kalman_filter(Y,start,last,a,P,kalman_tol,riccati_tol,rescale_prediction_error_covariance,presample,T,Q,R,H,Z,mm,pp,rr,Zflag,diffuse_periods,analytic_derivation,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P)
% Computes the log-likelihood of a stationary state space model.

%
% INPUTS
% - Y                       [matrix]      [pp x T] matrix of observed data
% - start                   [integer]     index of the first period processed in Y
% - last                    [integer]     index of the last period processed in Y
% - a                       [vector]      [mm x 1] initial mean of the state vector, E_0(alpha_1)
% - P                       [matrix]      [mm x mm] initial covariance matrix of the state vector, Var_0(alpha_1)
% - kalman_tol              [double]      tolerance parameter for inversion/conditioning of prediction-error covariance matrices
% - riccati_tol             [double]      tolerance parameter for Riccati fixed-point convergence
% - rescale_prediction_error_covariance
%                           [logical]     if true, rescales covariance matrix before inversion to improve numerical stability
% - presample               [integer]     number of initial iterations discarded when evaluating the likelihood
% - T                       [matrix]      [mm x mm] transition matrix of the state equation
% - Q                       [matrix]      [rr x rr] covariance matrix of structural innovations, or 3D array for time-varying Q
% - R                       [matrix]      [mm x rr] mapping from structural innovations to state innovations
% - H                       [matrix]      [pp x pp] covariance matrix of measurement errors
% - Z                       [matrix]      [pp x mm] measurement matrix, or index vector when Zflag=0
% - mm                      [integer]     number of state variables
% - pp                      [integer]     number of observed variables
% - rr                      [integer]     number of structural innovations
% - Zflag                   [integer]     0 if Z is an index vector; 1 if Z is a [pp x mm] matrix
% - diffuse_periods         [integer]     number of diffuse-filter periods already consumed during initialization
% - analytic_derivation     [integer]     derivative mode: 0 (none), 1 (score), 2 (score and Hessian), or asymptotic-Hessian mode
% - DT, DYss, DOm, DH, DP   [array]       first-derivative objects used when analytic_derivation > 0
% - D2T, D2Yss, D2Om, D2H, D2P
%                           [array]       second-derivative objects used when analytic_derivation == 2
%
% OUTPUTS
% - LIK                     [double|cell] minus log-likelihood; if analytic_derivation>0, returns cell array {LIK,DLIK[,Hess]}
% - LIKK                    [vector|cell] [smpl x 1] period-wise log-likelihood contributions; if analytic_derivation>0, returns {LIKK,dlik}
% - a                       [vector]      [mm x 1] filtered state mean at the end of the processed sample
% - P                       [matrix]      [mm x mm] filtered state covariance at the end of the processed sample
%
% This function is called by: dsge_likelihood
% This function calls: kalman_filter_ss

% Copyright © 2004-2026 Dynare Team
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

% Set defaults.
if nargin<17
    Zflag = 0;
end

if nargin<18
    diffuse_periods = 0;
end

if nargin<19
    analytic_derivation = 0;
end

if isempty(Zflag)
    Zflag = 0;
end

if isempty(diffuse_periods)
    diffuse_periods = 0;
end

% Get sample size.
smpl = last-start+1;

% Initialize some variables.
isqvec = false;
if ndims(Q)>2
    Qvec = Q;
    Q=Q(:,:,1);
    isqvec = true;
end
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
likk = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
oldK = Inf;
notsteady   = 1;
F_singular  = true;
asy_hess=0;

if  analytic_derivation == 0
    DLIK=[];
    Hess=[];
    LIKK=[];
else
    k = size(DT,3);                                 % number of structural parameters
    DLIK  = zeros(k,1);                             % Initialization of the score.
    Da    = zeros(mm,k);                            % Derivative State vector.
    dlikk = zeros(smpl,k);

    if analytic_derivation==2
        Hess  = zeros(k,k);                             % Initialization of the Hessian
        D2a    = zeros(mm,k,k);                             % State vector.
    else
        asy_hess=D2T;
        Hess=[];
        D2a=[];
        D2T=[];
        D2Yss=[];
    end
    if asy_hess
        Hess  = zeros(k,k);                             % Initialization of the Hessian
    end
    LIK={inf,DLIK,Hess};
    LIKK={likk,dlikk};
end

rescale_prediction_error_covariance0=rescale_prediction_error_covariance;
while notsteady && t<=last
    s = t-start+1;
    if isqvec
        QQ = R*Qvec(:,:,t+1)*transpose(R);
    end
    if Zflag
        v  = Y(:,t)-Z*a;
        F  = Z*P*Z' + H;
    else
        v  = Y(:,t)-a(Z);
        F  = P(Z,Z) + H;
    end
    badly_conditioned_F = false;
    if rescale_prediction_error_covariance
        sig=sqrt(diag(F));
        if any(diag(F)<kalman_tol) || rcond(F./(sig*sig'))<kalman_tol
            badly_conditioned_F = true;
        end
    else
        if rcond(F)<kalman_tol
            sig=sqrt(diag(F));
            if any(diag(F)<kalman_tol) || rcond(F./(sig*sig'))<kalman_tol
                badly_conditioned_F = true;
            else
                rescale_prediction_error_covariance=1;
            end
        end
    end
    if badly_conditioned_F
        % if ~all(abs(F(:))<kalman_tol), then use univariate filter (will remove
        % observations with zero variance prediction error), otherwise this is a
        % pathological case and the draw is discarded
        return
    else
        F_singular = false;
        if rescale_prediction_error_covariance
            log_dF = log(det(F./(sig*sig')))+2*sum(log(sig));
            iF = inv(F./(sig*sig'))./(sig*sig');
            rescale_prediction_error_covariance=rescale_prediction_error_covariance0;
        else
            log_dF = log(det(F));
            iF = inv(F);
        end
        likk(s) = log_dF+transpose(v)*iF*v;
        if Zflag
            K = P*Z'*iF;
            Ptmp = T*(P-K*Z*P)*transpose(T)+QQ;
        else
            K = P(:,Z)*iF;
            Ptmp = T*(P-K*P(Z,:))*transpose(T)+QQ;
        end
        tmp = (a+K*v);
        if analytic_derivation
            if analytic_derivation==2
                [Da,DP,DLIKt,D2a,D2P, Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady,D2a,D2Yss,D2T,D2Om,D2P);
            else
                [Da,DP,DLIKt,Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady);
            end
            if t>presample
                DLIK = DLIK + DLIKt;
                if analytic_derivation==2 || asy_hess
                    Hess = Hess + Hesst;
                end
            end
            dlikk(s,:)=DLIKt;
        end
        a = T*tmp;
        P = Ptmp;
        if ~isqvec
            notsteady = max(abs(K(:)-oldK))>riccati_tol;
        end
        oldK = K(:);
    end
    t = t+1;
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% Add observation's densities constants and divide by two.
likk(1:s) = .5*(likk(1:s) + pp*log(2*pi));
if analytic_derivation
    DLIK = DLIK/2;
    dlikk = dlikk/2;
    if analytic_derivation==2 || asy_hess
        if asy_hess==0
            Hess = Hess + tril(Hess,-1)';
        end
        Hess = -Hess/2;
    end
end

% Call steady state Kalman filter if needed.
if t <= last
    if analytic_derivation
        if analytic_derivation==2
            [tmp, tmp2] = kalman_filter_ss(Y, t, last, a, T, K, iF, log_dF, Z, pp, Zflag, analytic_derivation, Da, DT, DYss, D2a, D2T, D2Yss);
        else
            [tmp, tmp2] = kalman_filter_ss(Y, t, last, a, T, K, iF, log_dF, Z, pp, Zflag, analytic_derivation, Da, DT, DYss, asy_hess);
        end
        likk(s+1:end) = tmp2{1};
        dlikk(s+1:end,:) = tmp2{2};
        DLIK = DLIK + tmp{2};
        if analytic_derivation==2 || asy_hess
            Hess = Hess + tmp{3};
        end
    else
        [~, likk(s+1:end)] = kalman_filter_ss(Y, t, last, a, T, K, iF, log_dF, Z, pp, Zflag);
    end
end

% Compute minus the log-likelihood.
if presample>diffuse_periods
    LIK = sum(likk(1+(presample-diffuse_periods):end));
else
    LIK = sum(likk);
end

if analytic_derivation
    if analytic_derivation==2 || asy_hess
        LIK={LIK, DLIK, Hess};
    else
        LIK={LIK, DLIK};
    end
    LIKK={likk, dlikk};
else
    LIKK=likk;
end
