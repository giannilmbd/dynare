function [LIK, lik,a,P] = univariate_kalman_filter(data_index,no_more_missing_observations,Y,start,last,a,P,kalman_tol,riccati_tol,presample,T,Q,R,H,Z,mm,pp,Zflag,diffuse_periods,analytic_derivation,DT,DYss,DOm,DH,DP,analytic_Hessian,D2T,D2Yss,D2Om,D2H,D2P)
% Computes the log-likelihood of a stationary state space model (univariate approach).
%
% The function implements the univariate Kalman filter, processing one
% observable at a time within each period. This avoids inversion of the
% (pp x pp) multivariate forecast error variance matrix and improves
% numerical robustness in the presence of near-singular prediction-error
% variances.
%
% The state space model is given by:
%   y_t = Z * alpha_t + epsilon_t,           epsilon_t ~ N(0, H)
%   alpha_{t+1} = T * alpha_t + R * eta_t,   eta_t ~ N(0, Q)
%
% where H is assumed diagonal (required by the univariate approach). Each
% scalar observable y_{t,i} is processed sequentially with scalar forecast
% error variance:
%   F_{t,i} = Z_i * P_{t,i} * Z_i' + H_i
%
% For each scalar observable, two cases are distinguished:
%   (i)  F_{t,i} > kalman_tol: standard univariate update with contribution
%        log(F_{t,i}) + v_{t,i}^2/F_{t,i} + log(2*pi)
%   (ii) F_{t,i} <= kalman_tol: observable treated as uninformative
%        (state mean and covariance unchanged)
%
% Once missing data are no longer present (t >= no_more_missing_observations),
% the function checks convergence of Kalman gains using riccati_tol and,
% if converged, switches to the steady-state univariate Kalman filter.

%
% INPUTS
% - data_index              [cell]      1*T cell of column vectors of indices (in the vector of observed variables)
% - no_more_missing_observations
%                           [integer]   date after which there are no missing observations
% - Y                       [matrix]    [pp x T] matrix of observed data
% - start                   [integer]   index of the first period processed in Y
% - last                    [integer]   index of the last period processed in Y
% - a                       [vector]    [mm x 1] initial mean of the state vector, E_0(alpha_1)
% - P                       [matrix]    [mm x mm] initial covariance matrix of the state vector, Var_0(alpha_1)
% - kalman_tol              [double]    tolerance parameter for scalar forecast-error variances
% - riccati_tol             [double]    tolerance parameter for convergence of Kalman gains
% - presample               [integer]   number of initial iterations discarded when evaluating the likelihood
% - T                       [matrix]    [mm x mm] transition matrix of the state equation
% - Q                       [matrix]    [rr x rr] covariance matrix of structural innovations, or 3D array for time-varying Q
% - R                       [matrix]    [mm x rr] mapping from structural innovations to state innovations
% - H                       [vector]    [pp x 1] diagonal of covariance matrix of measurement errors
% - Z                       [matrix]    [pp x mm] measurement matrix, or index vector when Zflag=0
% - mm                      [integer]   number of state variables
% - pp                      [integer]   number of observed variables
% - Zflag                   [integer]   0 if Z is an index vector; 1 if Z is a [pp x mm] matrix
% - diffuse_periods         [integer]   number of diffuse-filter periods already consumed during initialization
% - analytic_derivation     [logical]   whether to compute analytic derivatives (true/false)
% - DT, DYss, DOm, DH, DP   [array]     first-derivative objects used when analytic_derivation is true
% - analytic_Hessian        [string]    Hessian mode: '' (none), 'full' (analytic), 'opg' (outer product), 'asymptotic'
% - D2T, D2Yss, D2Om, D2H, D2P
%                           [array]     second-derivative objects used when analytic_Hessian=='full'
%
% OUTPUTS
% - LIK                     [double|cell] minus log-likelihood; if analytic_derivation is true, returns cell array {LIK,DLIK[,Hess]}
% - lik                     [matrix|cell] [smpl x pp] observable-level log-likelihood contributions; if analytic_derivation is true, returns {lik,dlik}
% - a                       [vector]      [mm x 1] filtered state mean at the end of the processed sample
% - P                       [matrix]      [mm x mm] filtered state covariance at the end of the processed sample
%
% This function is called by: dsge_likelihood
% This function calls: univariate_kalman_filter_ss
%
% Algorithm:
%
%   Uses the univariate filter as described in Durbin/Koopman (2012): "Time
%   Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 6.4 and 7.2.5


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

% AUTHOR(S) stephane DOT adjemian AT univ DASH lemans DOT fr

if nargin<18 || isempty(Zflag)% Set default value for Zflag ==> Z is a vector of indices.
    Zflag = 0;
    diffuse_periods = 0;
end

if nargin<19
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
lik  = zeros(smpl,pp);     % Initialization of the matrix gathering the densities at each time and each observable
LIK  = Inf;                % Default value of the log likelihood.
l2pi = log(2*pi);
notsteady = 1;

oldK = Inf;
K = NaN(mm,pp);

% Extract Hessian mode flags from string
full_Hess = strcmp(analytic_Hessian, 'full');
asy_hess = strcmp(analytic_Hessian, 'asymptotic');

if ~analytic_derivation
    DLIK=[];
    Hess=[];
else
    k = size(DT,3);                                 % number of structural parameters
    DLIK  = zeros(k,1);                             % Initialization of the score.
    Da    = zeros(mm,k);                            % Derivative State vector.
    dlik  = zeros(smpl,k);

    if Zflag==0
        C = zeros(pp,mm);
        for ii=1:pp, C(ii,Z(ii))=1; end         % SELECTION MATRIX IN MEASUREMENT EQ. (FOR WHEN IT IS NOT CONSTANT)
    else
        C=Z;
    end
    if full_Hess
        Hess  = zeros(k,k);                             % Initialization of the Hessian
        D2a    = zeros(mm,k,k);                             % State vector.
    else
        Hess=[];
        D2a=[];
        D2T=[];
        D2Yss=[];
    end
    if asy_hess
        Hess  = zeros(k,k);                             % Initialization of the Hessian
    end
    LIK={inf,DLIK,Hess};
end

while notsteady && t<=last %loop over t
    s = t-start+1;
    d_index = data_index{t};
    if isqvec
        QQ = R*Qvec(:,:,t+1)*transpose(R);
    end
    if Zflag
        z = Z(d_index,:);
    else
        z = Z(d_index);
    end
    for i=1:rows(z) %loop over i
        if Zflag
            prediction_error = Y(d_index(i),t) - z(i,:)*a;  % nu_{t,i} in 6.13 in DK (2012)
            PZ = P*z(i,:)';                                 % Z_{t,i}*P_{t,i}*Z_{t,i}'
            Fi = z(i,:)*PZ + H(d_index(i));                 % F_{t,i} in 6.13 in DK (2012), relies on H being diagonal
        else
            prediction_error = Y(d_index(i),t) - a(z(i));   % nu_{t,i} in 6.13 in DK (2012)
            PZ = P(:,z(i));                                 % Z_{t,i}*P_{t,i}*Z_{t,i}'
            Fi = PZ(z(i)) + H(d_index(i));                  % F_{t,i} in 6.13 in DK (2012), relies on H being diagonal
        end
        if Fi>kalman_tol
            Ki =  PZ/Fi; %K_{t,i} in 6.13 in DK (2012)
            if t>=no_more_missing_observations
                K(:,i) = Ki;
            end
            lik(s,i) = log(Fi) + (prediction_error*prediction_error)/Fi + l2pi; %Top equation p. 175 in DK (2012)
            if analytic_derivation
                if full_Hess
                    [Da,DP,DLIKt,D2a,D2P, Hesst] = univariate_computeDLIK(k,i,z(i,:),Zflag,prediction_error,Ki,PZ,Fi,Da,DYss,DP,DH(d_index(i),:),notsteady,D2a,D2Yss,squeeze(D2H(d_index(i),:,:)),D2P);
                else
                    [Da,DP,DLIKt,Hesst] = univariate_computeDLIK(k,i,z(i,:),Zflag,prediction_error,Ki,PZ,Fi,Da,DYss,DP,DH(d_index(i),:),notsteady);
                end
                if t>presample
                    DLIK = DLIK + DLIKt;
                    if full_Hess || asy_hess
                        Hess = Hess + Hesst;
                    end
                end
                dlik(s,:)=dlik(s,:)+DLIKt';
            end
            a = a + Ki*prediction_error; %filtering according to (6.13) in DK (2012)
            P = P - PZ*Ki';              %filtering according to (6.13) in DK (2012)
        else
            if Fi<0
                %pathological numerical case where variance is negative
                if analytic_derivation
                    LIK={NaN,DLIK,Hess};
                else
                    LIK = NaN;
                end
                return
            else
                % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
                % p. 157, DK (2012)
            end
        end
    end
    if analytic_derivation
        if full_Hess
            [Da,DP,D2a,D2P] = univariate_computeDstate(k,a,P,T,Da,DP,DT,DOm,notsteady,D2a,D2P,D2T,D2Om);
        else
            [Da,DP] = univariate_computeDstate(k,a,P,T,Da,DP,DT,DOm,notsteady);
        end
    end
    a = T*a;            %transition according to (6.14) in DK (2012)
    P = T*P*T' + QQ;    %transition according to (6.14) in DK (2012)
    if t>=no_more_missing_observations && ~isqvec
        notsteady = max(abs(K(:)-oldK))>riccati_tol;
        oldK = K(:);
    end
    t = t+1;
end

% Divide by two.
lik(1:s,:) = .5*lik(1:s,:);
if analytic_derivation
    DLIK = DLIK/2;
    dlik = dlik/2;
    if full_Hess || asy_hess
        %         Hess = (Hess + Hess')/2;
        Hess = -Hess/2;
    end
end

% Call steady state univariate Kalman filter if needed.
if t <= last
    if analytic_derivation
        if full_Hess
            [tmp, tmp2] = univariate_kalman_filter_ss(Y,t,last,a,P,kalman_tol,T,H,Z,pp,Zflag, ...
                                                      analytic_derivation,Da,DT,DYss,DP,DH,analytic_Hessian,D2a,D2T,D2Yss,D2H,D2P);
        else
            [tmp, tmp2] = univariate_kalman_filter_ss(Y,t,last,a,P,kalman_tol,T,H,Z,pp,Zflag, ...
                                                      analytic_derivation,Da,DT,DYss,DP,DH,analytic_Hessian);
        end
        lik(s+1:end,:)=tmp2{1};
        dlik(s+1:end,:)=tmp2{2};
        DLIK = DLIK + tmp{2};
        if full_Hess || asy_hess
            Hess = Hess + tmp{3};
        end
    else
        [~, lik(s+1:end,:)] = univariate_kalman_filter_ss(Y,t,last,a,P,kalman_tol,T,H,Z,pp,Zflag);
    end
end

% Compute minus the log-likelihood.
if presample > diffuse_periods
    LIK = sum(sum(lik(1+presample-diffuse_periods:end,:)));
else
    LIK = sum(sum(lik));
end

if analytic_derivation
    if full_Hess || asy_hess
        LIK={LIK, DLIK, Hess};
    else
        LIK={LIK, DLIK};
    end
    lik={lik, dlik};
end
