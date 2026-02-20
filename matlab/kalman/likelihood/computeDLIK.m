function [Da,DP,DLIK,D2a,D2P,Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady,D2a,D2Yss,D2T,D2Om,D2H,D2P)
% [Da,DP,DLIK,D2a,D2P,Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady,D2a,D2Yss,D2T,D2Om,D2H,D2P)
% Compute first and second derivatives of the Kalman filter log-likelihood
%
% USAGE:
%   [Da,DP,DLIK,D2a,D2P,Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady,D2a,D2Yss,D2T,D2Om,D2H,D2P)
%
% INPUTS:
%   k           - number of parameters
%   tmp         - temporary storage for state transition computation
%   Z           - observation matrix (or index if Zflag=0)
%   Zflag       - flag indicating Z format (1: matrix, 0: index vector)
%   v           - innovation (measurement residuals)
%   T           - state transition matrix
%   K           - Kalman gain
%   P           - state prediction error covariance
%   iF          - inverse of innovation forecast error variance
%   Da          - derivatives of state prediction (input/output)
%   DYss        - derivatives of steady state
%   DT          - derivatives of state transition matrix
%   DOm         - derivatives of covariance matrix
%   DP          - derivatives of state covariance (input/output)
%   DH          - derivatives of measurement error variance
%   notsteady   - flag indicating if not at steady state
%   D2a         - second derivatives of state prediction (input/output)
%   D2Yss       - second derivatives of steady state
%   D2T         - second derivatives of state transition matrix
%   D2Om        - second derivatives of covariance matrix
%   D2H         - second derivatives of measurement error variance
%   D2P         - second derivatives of state covariance (input/output)
%
% OUTPUTS:
%   Da          - updated derivatives of state prediction
%   DP          - updated derivatives of state covariance
%   DLIK        - (k x 1) first derivatives of log-likelihood
%   D2a         - second derivatives of state prediction (if nargout > 4)
%   D2P         - second derivatives of state covariance (if nargout > 4)
%   Hesst       - Hessian of log-likelihood w.r.t. parameters (if nargout == 6)
%
% APPROACH:
%   The function computes the gradient and Hessian of the Kalman filter log-likelihood
%   by differentiating the likelihood contributions from each observation:
%   1. Compute derivatives of Kalman matrices (K, F) based on model derivatives
%   2. Compute derivatives of innovations (Dv) accounting for observation equation changes
%   3. For each parameter, update state derivatives and accumulate likelihood derivatives
%   4. Optionally compute second derivatives for the Hessian matrix

% Copyright © 2012-2026 Dynare Team
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

% AUTHOR(S) marco.ratto@jrc.ec.europa.eu

persistent DK DF D2K D2F

%% Compute derivatives of Kalman filter matrices (K and F) if not at steady state
if notsteady
    if Zflag
        [DK,DF,DP1] = computeDKalmanZ(T,DT,DOm,P,DP,DH,Z,iF,K);
        if nargout>4
            [D2K,D2F,D2P] = computeD2KalmanZ(T,DT,D2T,D2Om,P,DP,D2P,DH,D2H,Z,iF,K,DK);
        end
    else
        [DK,DF,DP1] = computeDKalman(T,DT,DOm,P,DP,DH,Z,iF,K);
        if nargout>4
            [D2K,D2F,D2P] = computeD2Kalman(T,DT,D2T,D2Om,P,DP,D2P,DH,D2H,Z,iF,K,DK);
        end
    end
    DP=DP1;
    clear DP1
end

%% Compute derivatives of innovations (Dv) w.r.t. parameters
% The innovation is v = y - Z*a where a is the state prediction
% Its derivative accounts for measurement matrix and state prediction changes
Dv=zeros(length(v),k);
for ii = 1:k
    if Zflag
        % When Z is a matrix: Dv = -Z*Da - Z*DYss
        Dv(:,ii) = -Z*Da(:,ii) - Z*DYss(:,ii);
    else
        % When Z is an index vector: extract relevant rows
        Dv(:,ii) = -Da(Z,ii) - DYss(Z,ii);
    end
end

%% Compute first derivatives of log-likelihood and first-order state updates
% For each parameter, update the state derivative and accumulate the likelihood gradient contribution
Hesst = zeros(k,k);
DLIK = zeros(k,1);
dtmp = zeros(size(Da));
jcount = 0;
for ii = 1:k
    % Compute state update: da(t) = DT*tmp + T*dtmp
    % where dtmp accounts for derivatives of Kalman gain and innovations
    dKi  = DK(:,:,ii);
    dtmp(:,ii) = Da(:,ii)+dKi*v+K*Dv(:,ii);

    if nargout>4
        % Compute second derivatives of log-likelihood (Hessian)
        diFi = -iF*DF(:,:,ii)*iF;
        for jj = 1:ii
            jcount=jcount+1;
            dFj    = DF(:,:,jj);
            diFj   = -iF*DF(:,:,jj)*iF;
            dKj  = DK(:,:,jj);
            d2Kij  = D2K(:,:,jj,ii);
            d2Fij  = D2F(:,:,jj,ii);
            d2iFij = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;

            % Compute second derivative of innovations
            if Zflag
                d2vij  = -Z*D2Yss(:,jj,ii)  - Z*D2a(:,jj,ii);
            else
                d2vij  = -D2Yss(Z,jj,ii)  - D2a(Z,jj,ii);
            end

            % Compute second derivative of state update
            d2tmpij = D2a(:,jj,ii) + d2Kij*v + dKj*Dv(:,ii) + dKi*Dv(:,jj) + K*d2vij;
            D2a(:,jj,ii) = reshape(D2T(:,jcount),size(T))*tmp + DT(:,:,jj)*dtmp(:,ii) + DT(:,:,ii)*dtmp(:,jj) + T*d2tmpij;
            D2a(:,ii,jj) = D2a(:,jj,ii);

            if nargout==6
                Hesst(ii,jj) = getHesst_ij(v,Dv(:,ii),Dv(:,jj),d2vij,iF,diFi,diFj,d2iFij,dFj,d2Fij);
            end
        end
    end

    % Update state derivative: Da(t+1) = DT*tmp + T*dtmp
    Da(:,ii)   = DT(:,:,ii)*tmp + T*dtmp(:,ii);

    % Accumulate first derivative of log-likelihood
    % DLIK = trace(iF*DF) + 2*Dv'*iF*v - v'*(iF*DF*iF)*v
    DLIK(ii,1)  = trace( iF*DF(:,:,ii) ) + 2*Dv(:,ii)'*iF*v - v'*(iF*DF(:,:,ii)*iF)*v;
end

%% Compute approximate Hessian when only 4 outputs requested
% Uses Kronecker product approximation to estimate second derivatives
if nargout==4
    vecDPmf = reshape(DF,[],k);
    D2a = 2*Dv'*iF*Dv + (vecDPmf' * kron(iF,iF) * vecDPmf);
end

%% Helper function: Compute (i,j) term of the Hessian
function Hesst_ij = getHesst_ij(e,dei,dej,d2eij,iS,diSi,diSj,d2iSij,dSj,d2Sij)
% GETHESST_IJ Compute the (i,j)-th element of the Hessian of log-likelihood
%
% INPUTS:
%   e       - innovations/residuals
%   dei     - derivative of innovations w.r.t. i-th parameter
%   dej     - derivative of innovations w.r.t. j-th parameter
%   d2eij   - second derivative of innovations w.r.t. i-th and j-th parameters
%   iS      - inverse of forecast error covariance (iF)
%   diSi    - derivative of inverse covariance w.r.t. i-th parameter
%   diSj    - derivative of inverse covariance w.r.t. j-th parameter
%   d2iSij  - second derivative of inverse covariance
%   dSj     - derivative of forecast error covariance w.r.t. j-th parameter
%   d2Sij   - second derivative of forecast error covariance
%
% OUTPUT:
%   Hesst_ij - Hessian element H(i,j)

Hesst_ij = trace(diSi*dSj + iS*d2Sij) + e'*d2iSij*e + 2*(dei'*diSj*e + dei'*iS*dej + e'*diSi*dej + e'*iS*d2eij);

%% Nested function: Compute first derivatives of Kalman matrices (with Z as index)
function [DK,DF,DP1] = computeDKalman(T,DT,DOm,P,DP1,DH,Z,iF,K)
% COMPUTEDKALMAN Compute first derivatives of Kalman filter matrices
%   Used when observation matrix Z is represented as an index vector
%
% Computes derivatives of innovation variance (DF), Kalman gain (DK),
% and filtered state covariance (DP1) with respect to model parameters

k      = size(DT,3);
tmp    = P-K*P(Z,:);
DF = zeros([size(iF),k]);
DK = zeros([size(K),k]);

for ii = 1:k
    % Derivative of innovation variance: DF = DP(Z,Z) + DH
    DF(:,:,ii)  = DP1(Z,Z,ii) + DH(:,:,ii);
    DiF = -iF*DF(:,:,ii)*iF;

    % Derivative of Kalman gain: DK = DP(:,Z)*iF + P(:,Z)*DiF
    DK(:,:,ii)  = DP1(:,Z,ii)*iF + P(:,Z)*DiF;

    % Derivative of filtered covariance: DP(t+1) = DT*tmp*T' + T*Dtmp*T' + T*tmp*DT' + DOm
    Dtmp        = DP1(:,:,ii) - DK(:,:,ii)*P(Z,:) - K*DP1(Z,:,ii);
    DP1(:,:,ii) = DT(:,:,ii)*tmp*T' + T*Dtmp*T' + T*tmp*DT(:,:,ii)' + DOm(:,:,ii);
end

%% Nested function: Compute first derivatives of Kalman matrices (with Z as matrix)
function [DK,DF,DP1] = computeDKalmanZ(T,DT,DOm,P,DP,DH,Z,iF,K)
% COMPUTEDKALMAN_Z Compute first derivatives of Kalman filter matrices
%   Used when observation matrix Z is represented as a matrix
%
% Computes derivatives of innovation variance (DF), Kalman gain (DK),
% and filtered state covariance (DP1) with respect to model parameters

k      = size(DT,3);
tmp    = P-K*Z*P;
DF = zeros([size(iF),k]);
DK = zeros([size(K),k]);
DP1 = zeros([size(P),k]);

for ii = 1:k
    % Derivative of innovation variance: DF = Z*DP*Z' + DH
    DF(:,:,ii)  = Z*DP(:,:,ii)*Z + DH(:,:,ii);
    DiF = -iF*DF(:,:,ii)*iF;

    % Derivative of Kalman gain: DK = DP*Z'*iF + P*Z'*DiF
    DK(:,:,ii)  = DP(:,:,ii)*Z*iF + P(:,:)*Z*DiF;

    % Derivative of filtered covariance: DP(t+1) = DT*tmp*T' + T*Dtmp*T' + T*tmp*DT' + DOm
    Dtmp        = DP(:,:,ii) - DK(:,:,ii)*Z*P(:,:) - K*Z*DP(:,:,ii);
    DP1(:,:,ii) = DT(:,:,ii)*tmp*T' + T*Dtmp*T' + T*tmp*DT(:,:,ii)' + DOm(:,:,ii);
end

%% Nested function: Compute second derivatives of Kalman matrices (with Z as index)
function [d2K,d2S,d2P1] = computeD2Kalman(A,dA,d2A,d2Om,P0,dP0,d2P1,DH,D2H,Z,iF,K0,dK0)
% COMPUTED2KALMAN Compute second derivatives of Kalman filter matrices
%   Used when observation matrix Z is represented as an index vector
%   where A = T (state transition matrix) in the main function
%
% Computes second derivatives of Kalman gain (d2K), innovation variance (d2S),
% and filtered state covariance (d2P1) with respect to model parameters

k      = size(dA,3);
tmp    = P0-K0*P0(Z,:);
[ns,no] = size(K0);

d2K  = zeros(ns,no,k,k);
d2S  = zeros(no,no,k,k);

jcount=0;
for ii = 1:k
    dAi = dA(:,:,ii);
    dFi = dP0(Z,Z,ii) + DH(:,:,ii);
    diFi = -iF*dFi*iF;
    dKi = dK0(:,:,ii);
    for jj = 1:ii
        jcount=jcount+1;
        dAj = dA(:,:,jj);
        dFj = dP0(Z,Z,jj) + DH(:,:,jj);
        diFj = -iF*dFj*iF;
        dKj = dK0(:,:,jj);

        % Unvectorize second derivatives
        d2Aij = reshape(d2A(:,jcount),[ns ns]);
        d2Pij = dyn_unvech(d2P1(:,jcount));
        d2Omij = dyn_unvech(d2Om(:,jcount));

        % Second derivative of innovation variance
        d2Fij = d2Pij(Z,Z) + D2H(:,:,jj,ii);

        % Second derivative of inverse covariance
        d2iF = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;

        % Second derivative of Kalman gain
        d2Kij= d2Pij(:,Z)*iF + P0(:,Z)*d2iF + dP0(:,Z,jj)*diFi + dP0(:,Z,ii)*diFj;

        % Accumulated Kalman gain effect on covariance
        d2KCP = d2Kij*P0(Z,:) + K0*d2Pij(Z,:) + dKi*dP0(Z,:,jj) + dKj*dP0(Z,:,ii);

        % First-order derivatives of temporary covariance term
        dtmpi        = dP0(:,:,ii) - dK0(:,:,ii)*P0(Z,:) - K0*dP0(Z,:,ii);
        dtmpj        = dP0(:,:,jj) - dK0(:,:,jj)*P0(Z,:) - K0*dP0(Z,:,jj);
        d2tmp = d2Pij - d2KCP;

        % Accumulate second derivatives of state transition equation
        d2AtmpA = d2Aij*tmp*A' + A*d2tmp*A' + A*tmp*d2Aij' + dAi*dtmpj*A' + dAj*dtmpi*A' + A*dtmpj*dAi' + A*dtmpi*dAj' + dAi*tmp*dAj' + dAj*tmp*dAi';

        % Store symmetric pairs
        d2K(:,:,ii,jj)  = d2Kij;
        d2P1(:,jcount) = dyn_vech(d2AtmpA  + d2Omij);
        d2S(:,:,ii,jj)  = d2Fij;
        d2K(:,:,jj,ii)  = d2Kij;
        d2S(:,:,jj,ii)  = d2Fij;
    end
end

%% Nested function: Compute second derivatives of Kalman matrices (with Z as matrix)
function [d2K,d2S,d2P1] = computeD2KalmanZ(A,dA,d2A,d2Om,P0,dP0,d2P1,DH,D2H,Z,iF,K0,dK0)
% COMPUTED2KALMAN_Z Compute second derivatives of Kalman filter matrices
%   Used when observation matrix Z is represented as a matrix
%   where A = T (state transition matrix) in the main function
%
% Computes second derivatives of Kalman gain (d2K), innovation variance (d2S),
% and filtered state covariance (d2P1) with respect to model parameters

k      = size(dA,3);
tmp    = P0-K0*Z*P0(:,:);
[ns,no] = size(K0);

d2K  = zeros(ns,no,k,k);
d2S  = zeros(no,no,k,k);

jcount=0;
for ii = 1:k
    dAi = dA(:,:,ii);
    dFi = Z*dP0(:,:,ii)*Z + DH(:,:,ii);
    diFi = -iF*dFi*iF;
    dKi = dK0(:,:,ii);
    for jj = 1:ii
        jcount=jcount+1;
        dAj = dA(:,:,jj);
        dFj = Z*dP0(:,:,jj)*Z + DH(:,:,jj);
        diFj = -iF*dFj*iF;
        dKj = dK0(:,:,jj);

        % Unvectorize second derivatives
        d2Aij = reshape(d2A(:,jcount),[ns ns]);
        d2Pij = dyn_unvech(d2P1(:,jcount));
        d2Omij = dyn_unvech(d2Om(:,jcount));

        % Second derivative of innovation variance
        d2Fij = Z*d2Pij(:,:)*Z + D2H(:,:,jj,ii);

        % Second derivative of inverse covariance
        d2iF = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;

        % Second derivative of Kalman gain
        d2Kij= d2Pij(:,:)*Z*iF + P0(:,:)*Z*d2iF + dP0(:,:,jj)*Z*diFi + dP0(:,:,ii)*Z*diFj;

        % Accumulated Kalman gain effect on covariance
        d2KCP = d2Kij*Z*P0(:,:) + K0*Z*d2Pij(:,:) + dKi*Z*dP0(:,:,jj) + dKj*Z*dP0(:,:,ii);

        % First-order derivatives of temporary covariance term
        dtmpi        = dP0(:,:,ii) - dK0(:,:,ii)*Z*P0(:,:) - K0*Z*dP0(:,:,ii);
        dtmpj        = dP0(:,:,jj) - dK0(:,:,jj)*Z*P0(:,:) - K0*Z*dP0(:,:,jj);
        d2tmp = d2Pij - d2KCP;

        % Accumulate second derivatives of state transition equation
        d2AtmpA = d2Aij*tmp*A' + A*d2tmp*A' + A*tmp*d2Aij' + dAi*dtmpj*A' + dAj*dtmpi*A' + A*dtmpj*dAi' + A*dtmpi*dAj' + dAi*tmp*dAj' + dAj*tmp*dAi';

        % Store symmetric pairs
        d2K(:,:,ii,jj)  = d2Kij;
        d2P1(:,jcount) = dyn_vech(d2AtmpA  + d2Omij);
        d2S(:,:,ii,jj)  = d2Fij;
        d2K(:,:,jj,ii)  = d2Kij;
        d2S(:,:,jj,ii)  = d2Fij;
    end
end
