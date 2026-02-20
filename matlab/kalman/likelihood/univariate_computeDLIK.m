function [Da,DP1,DLIK,D2a,D2P,Hesst] = univariate_computeDLIK(k,indx,Z,Zflag,v,K,PZ,F,Da,DYss,DP,DH,notsteady,D2a,D2Yss,D2H,D2P)
% [Da,DP1,DLIK,D2a,D2P,Hesst] = univariate_computeDLIK(k,indx,Z,Zflag,v,K,PZ,F,Da,DYss,DP,DH,notsteady,D2a,D2Yss,D2H,D2P)
% Compute first and second derivatives of the univariate Kalman filter log-likelihood
%
% USAGE:
%   [Da,DP1,DLIK,D2a,D2P,Hesst] = univariate_computeDLIK(k,indx,Z,Zflag,v,K,PZ,F,Da,DYss,DP,DH,notsteady,D2a,D2Yss,D2H,D2P)
%
% INPUTS:
%   k           - number of parameters
%   indx        - index of current observation in univariate filter sequence
%   Z           - observation matrix (or index if Zflag=0)
%   Zflag       - flag indicating Z format (1: matrix, 0: index vector)
%   v           - innovation (measurement residual)
%   K           - Kalman gain
%   PZ          - P*Z' (state covariance times observation matrix transpose)
%   F           - innovation forecast error variance (scalar)
%   Da          - derivatives of state prediction (input/output)
%   DYss        - derivatives of steady state
%   DP          - derivatives of state covariance (input/output)
%   DH          - derivatives of measurement error variance
%   notsteady   - flag indicating if not at steady state
%   D2a         - second derivatives of state prediction (input/output)
%   D2Yss       - second derivatives of steady state
%   D2H         - second derivatives of measurement error variance
%   D2P         - second derivatives of state covariance (input/output)
%
% OUTPUTS:
%   Da          - updated derivatives of state prediction
%   DP1         - updated derivatives of filtered state covariance
%   DLIK        - (k x 1) first derivatives of log-likelihood
%   D2a         - second derivatives of state prediction (if nargout > 4)
%   D2P         - second derivatives of state covariance (if nargout > 4)
%   Hesst       - Hessian of log-likelihood w.r.t. parameters (if nargout == 6)
%
% APPROACH:
%   The function computes the gradient and Hessian of the univariate Kalman filter
%   log-likelihood by processing observations one at a time. For each observation:
%   1. Compute derivatives of Kalman gain (DK) and innovation variance (DF)
%   2. Compute derivatives of innovations (Dv) from observation equation
%   3. Update state derivatives and accumulate likelihood derivatives
%   4. Optionally compute second derivatives for the Hessian matrix
%   5. Update filtered covariance derivatives for the next iteration

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

%% Persistent storage for derivatives across univariate filter recursions
persistent DDK DDF DD2K DD2F

%% Compute derivatives of Kalman filter quantities if not at steady state
if notsteady
    if Zflag
        % Case 1: Z is a matrix
        % Compute first derivative of innovation: Dv = -Z*Da - Z*DYss
        Dv   = -Z*Da(:,:) - Z*DYss(:,:);
        DF = zeros(k,1);
        DK = zeros([rows(K),k]);

        % Compute first derivatives of forecast error variance (DF) and Kalman gain (DK)
        for j=1:k
            % DF = Z*DP*Z' + DH
            DF(j)=Z*DP(:,:,j)*Z'+DH(j);
            % DK = DP*Z'/F - PZ*DF/F^2
            DK(:,j) = (DP(:,:,j)*Z')/F-PZ*DF(j)/F^2;
        end

        % Compute second derivatives if requested
        if nargout>4
            D2F = zeros(k,k);
            D2v = zeros(k,k);
            D2K = zeros(rows(K),k,k);
            jcount=0;
            for j=1:k
                % Second derivative of innovations
                D2v(:,j)   = -Z*D2a(:,:,j) - Z*D2Yss(:,:,j);
                for i=1:j
                    jcount=jcount+1;
                    % Second derivative of forecast error variance
                    D2F(j,i)=Z*dyn_unvech(D2P(:,jcount))*Z'+D2H(j,i);
                    D2F(i,j)=D2F(j,i);
                    % Second derivative of Kalman gain
                    D2K(:,j,i) = (dyn_unvech(D2P(:,jcount))*Z')/F-(DP(:,:,j)*Z')*DF(i)/F^2-(DP(:,:,i)*Z')*DF(j)/F^2- ...
                        PZ*D2F(j,i)/F^2 + 2*PZ/F^3*DF(i)*DF(j);
                    D2K(:,i,j) = D2K(:,j,i);
                end
            end
        end

    else
        % Case 2: Z is an index vector
        % Compute first derivative of innovation
        Dv   = -Da(Z,:) - DYss(Z,:);
        % Extract derivatives from indexed rows/columns
        DF = squeeze(DP(Z,Z,:))+DH';
        DK = squeeze(DP(:,Z,:))/F-PZ*transpose(DF)/F^2;

        % Compute second derivatives if requested
        if nargout>4
            D2F = zeros(k,k);
            D2K = zeros(rows(K),k,k);
            D2v   = squeeze(-D2a(Z,:,:) - D2Yss(Z,:,:));
            jcount=0;
            for j=1:k
                for i=1:j
                    jcount=jcount+1;
                    tmp = dyn_unvech(D2P(:,jcount));
                    % Second derivative of forecast error variance
                    D2F(j,i) = tmp(Z,Z)+D2H(j,i);
                    D2F(i,j)=D2F(j,i);
                    % Second derivative of Kalman gain
                    D2K(:,i,j) = tmp(:,Z)/F;
                    D2K(:,i,j) = D2K(:,i,j) -PZ*D2F(j,i)/F^2 - squeeze(DP(:,Z,i))*DF(j)/F^2 - ...
                        squeeze(DP(:,Z,j))*DF(i)'/F^2 + 2/F^3*PZ*DF(i)'*DF(j);
                    D2K(:,j,i) = D2K(:,i,j);
                end
            end
        end
    end

    % Store computed derivatives in persistent storage for this observation index
    if nargout>4
        DD2K(:,indx,:,:)=D2K;
        DD2F(indx,:,:)=D2F;
    end
    DDK(:,indx,:)=DK;
    DDF(indx,:)=DF;
else
    % At steady state: retrieve pre-computed derivatives from persistent storage
    DK = squeeze(DDK(:,indx,:));
    DF = DDF(indx,:)';
    if nargout>4
        D2K = squeeze(DD2K(:,indx,:,:));
        D2F = squeeze(DD2F(indx,:,:));
    end

    % Compute innovation derivatives (state-dependent even at steady state)
    if Zflag
        Dv   = -Z*Da(:,:) - Z*DYss(:,:);
        if nargout>4
            D2v = zeros(k,k);
            for j=1:k
                D2v(:,j)   = -Z*D2a(:,:,j) - Z*D2Yss(:,:,j);
            end
        end
    else
        Dv   = -Da(Z,:) - DYss(Z,:);
        if nargout>4
            D2v   = squeeze(-D2a(Z,:,:) - D2Yss(Z,:,:));
        end
    end
end

%% Compute first derivative of log-likelihood contribution
% DLIK = DF/F + 2*Dv'*v/F - v^2*DF/F^2
DLIK = DF/F + 2*Dv'/F*v - v^2/F^2*DF;

%% Compute second derivative of log-likelihood (Hessian) or approximate Hessian
if nargout==6
    % Full Hessian computation
    Hesst = D2F/F-1/F^2*(DF*DF') + 2*D2v/F*v + 2*(Dv'*Dv)/F - 2*(DF*Dv)*v/F^2 ...
            - v^2/F^2*D2F - 2*v/F^2*(Dv'*DF') + 2*v^2/F^3*(DF*DF');
elseif nargout==4
    % Approximate Hessian using outer product of scores
    D2a = 1/F^2*(DF*DF') + 2*(Dv'*Dv)/F ;
end

%% Update state prediction derivatives
% Da(t+1) = Da(t) + DK*v + K*Dv
Da = Da + DK*v+K*Dv;

%% Update second derivatives of state prediction if requested
if nargout>4
    D2a = D2a + D2K*v;
    for j=1:k
        for i=1:j
            % D2a = D2a + DK_i*Dv_j + DK_j*Dv_i + K*D2v
            D2a(:,j,i) = D2a(:,j,i) + DK(:,i)*Dv(j) + DK(:,j)*Dv(i) + K*D2v(j,i);
            D2a(:,i,j) = D2a(:,j,i);
        end
    end
end

%% Update filtered covariance derivatives
if notsteady
    DP1 = DP*0;
    if Zflag
        % Case 1: Z is a matrix
        % DP1 = DP - DP*Z'*K' - P*Z'*DK'
        for j=1:k
            DP1(:,:,j)=DP(:,:,j) - (DP(:,:,j)*Z')*K'-PZ*DK(:,j)';
        end
    else
        % Case 2: Z is an index vector
        for j=1:k
            DP1(:,:,j)=DP(:,:,j) - (DP(:,Z,j))*K'-PZ*DK(:,j)';
        end
    end

    % Update second derivatives of filtered covariance if requested
    if nargout>4
        if Zflag
            jcount = 0;
            for j=1:k
                for i=1:j
                    jcount = jcount+1;
                    tmp = dyn_unvech(D2P(:,jcount));
                    % D2P1 = D2P - D2P*Z'*K' - DP*Z'*DK' - DP*Z'*DK' - P*Z'*D2K'
                    tmp = tmp - (tmp*Z')*K' - (DP(:,:,j)*Z')*DK(:,i)' ...
                          - (DP(:,:,i)*Z')*DK(:,j)' -PZ*D2K(:,j,i)';
                    D2P(:,jcount) = dyn_vech(tmp);
                end
            end
        else
            DPZ = squeeze(DP(:,Z,:));
            jcount = 0;
            for j=1:k
                for i=1:j
                    jcount = jcount+1;
                    tmp = dyn_unvech(D2P(:,jcount));
                    D2PZ = tmp(:,Z);
                    tmp = tmp - D2PZ*K' - DPZ(:,j)*DK(:,i)'- DPZ(:,i)*DK(:,j)' - PZ*squeeze(D2K(:,j,i))';
                    D2P(:,jcount) = dyn_vech(tmp);
                end
            end
        end
    end
else
    % At steady state: filtered covariance derivatives equal predicted covariance derivatives
    DP1=DP;
end
