function [dLIK,dlik,a,Pstar] = kalman_filter_d(Y, start, last, a, Pinf, Pstar, kalman_tol, diffuse_kalman_tol, presample, T, R, Q, H, Z, pp)
% [dLIK,dlik,a,Pstar] = kalman_filter_d(Y, start, last, a, Pinf, Pstar, kalman_tol, diffuse_kalman_tol, presample, T, R, Q, H, Z, pp)
% Computes the diffuse log-likelihood of a state space model.
%
% The function implements the exact diffuse Kalman filter for the initialization
% of the state vector when some or all initial conditions are unknown (diffuse).
% It iterates until the diffuse part of the state covariance matrix has
% converged to zero (i.e., Pinf becomes numerically zero), at which point
% the standard Kalman filter can take over.
%
% The state space model is given by:
%   y_t = Z * alpha_t + epsilon_t,   epsilon_t ~ N(0, H)
%   alpha_{t+1} = T * alpha_t + R * eta_t,   eta_t ~ N(0, Q)
%
% with the diffuse initialization:
%   alpha_1 ~ N(a, Pinf * kappa + Pstar)  as kappa -> infinity
%
% This implies a decomposition of the forecast error variance:
%   F_t = F_{inf,t} * kappa + F_{*,t}    as kappa -> infinity
% where:
%   F_{inf,t} = Z * P_{inf,t} * Z'       diffuse forecast error variance
%   F_{*,t}   = Z * P_{*,t} * Z' + H     stationary forecast error variance
%
% The filter distinguishes two cases at each period t:
%   (i)  F_{inf,t} full rank: the observation resolves diffuse uncertainty.
%        The likelihood contribution is log|F_{inf,t}| (no quadratic term,
%        as the forecast error variance diverges with kappa).
%   (ii) F_{inf,t} = 0: diffuse uncertainty does not affect this observable.
%        The standard likelihood contribution applies:
%        log|F_{*,t}| + v_t' * F_{*,t}^{-1} * v_t
%        where v_t = y_t - Z * a_t is the prediction error.
%
% The diffuse phase ends once P_{inf,t} (and hence F_{inf,t}) has
% converged to zero for all observables, after which the standard
% Kalman filter takes over.
%
% INPUTS
% - Y                       [matrix]    [pp x smpl] matrix of (detrended) data
% - start                   [integer]   index of the first observation to process in Y
% - last                    [integer]   index of the last observation to process in Y
% - a                       [vector]    initial state vector (E_0(alpha_1))
%   Pinf                    [double]    (m x m) diffuse part of the initial state covariance matrix;
%                                       reflects prior uncertainty about nonstationary components
%   Pstar                   [double]    (m x m) stationary part of the initial state covariance matrix;
%                                       reflects prior uncertainty about stationary components
% - kalman_tol              [double]    tolerance parameter (rcond) of F_star
% - diffuse_kalman_tol      [double]    tolerance parameter (rcond) of Pinf
% - presample               [integer]   number of initial iterations to be discarded when evaluating the likelihood
% - T                       [matrix]    transition matrix in the state equations
% - R                       [matrix]    matrix relating the structural innovations to the state vector
% - Q                       [matrix]    covariance matrix of the structural innovations
% - H                       [matrix]    covariance matrix of the measurement errors
% - Z                       [matrix]    matrix relating states to the observed variables
% - pp                      [integer]   number of observed variables
%
% OUTPUTS
% - dLIK                    [double]    minus loglikelihood
% - dlik                    [vector]    contributions to the log-likelihood at each period of the diffuse phase, where s is the number of
%                                       diffuse iterations; each element equals 0.5*(w_t + pp*log(2*pi))
%                                       with w_t as defined in Durbin and Koopman (2012, p. 172)
% - a                       [vector]    estimated state vector at the end of the diffuse phase,
%                                       E_{t_d}(alpha_{t_d+1}), to be used as initial condition for
%                                       the standard Kalman filter
% - Pstar                   [matrix]    state covariance matrix at the end of the diffuse phase,
%                                       Var_{t_d}(alpha_{t_d+1}), to be used as initial condition for
%                                       the standard Kalman filter
%
% This function is called by: dsge_likelihood
% This function calls: none
%
% References:
%   Koopman, S.J. and Durbin, J. (2003), "Filtering and Smoothing of State
%   Vector for Diffuse State Space Models", Journal of Time Series Analysis,
%   vol. 24(1), pp. 85-98.
%
%   Durbin, J. and Koopman, S.J. (2012), "Time Series Analysis by State Space
%   Methods", Oxford University Press, Second Edition, Ch. 5 and 7.2

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
dlik = zeros(smpl,1);      % Initialization of the vector gathering the densities.
dLIK = Inf;                % Default value of the log likelihood.
s    = 0;

while rank(Z*Pinf*Z',diffuse_kalman_tol) && (t<=last)
    s = t-start+1;
    v = Y(:,t)-Z*a;                                                     %get prediction error v^(0) in (5.13) DK (2012)
    Finf  = Z*Pinf*Z';                                                  % (5.7) in DK (2012)
    if isqvec
        QQ = R*Qvec(:,:,t+1)*transpose(R);
    end
                                                                        %do case distinction based on whether F_{\infty,t} has full rank or 0 rank
    if rcond(Finf) < diffuse_kalman_tol                                 %F_{\infty,t} = 0
        if ~all(abs(Finf(:)) < diffuse_kalman_tol)                      %rank-deficient but not rank 0
                                                                        % The univariate diffuse Kalman filter should be used instead.
            return
        else                                                            %rank of F_{\infty,t} is 0
            Fstar  = Z*Pstar*Z' + H;                                    % (5.7) in DK (2012)
            if rcond(Fstar) < kalman_tol                                %F_{*} is singular
                % if ~all(abs(Fstar(:))<kalman_tol), then use univariate diffuse filter,
                % otherwise this is a pathological case and the draw is discarded
                return
            else
                iFstar = inv(Fstar);
                dFstar = det(Fstar);
                Kstar  = Pstar*Z'*iFstar;                               %(5.15) of DK (2012) with Kstar=T^{-1}*K^(0)
                dlik(s)= log(dFstar) + v'*iFstar*v;                     %set w_t to bottom case in bottom equation page 172, DK (2012)
                Pinf   = T*Pinf*transpose(T);                           % (5.16) DK (2012)
                Pstar  = T*(Pstar-Pstar*Z'*Kstar')*T'+QQ;               % (5.17) DK (2012)
                a      = T*(a+Kstar*v);                                 % (5.13) DK (2012)
            end
        end
    else                                                                %F_{\infty,t} positive definite
                                                                        %To compare to DK (2012), this block makes use of the following transformation
                                                                        %Kstar=T^{-1}*K^{(1)}=M_{*}*F^{(1)}+M_{\infty}*F^{(2)}
                                                                        %     =P_{*}*Z'*F^{(1)}+P_{\infty}*Z'*((-1)*(-F_{\infty}^{-1})*F_{*}*(F_{\infty}^{-1}))
                                                                        %     =[P_{*}*Z'-Kinf*F_{*})]*F^{(1)}
                                                                        %Make use of L^{0}'=(T-K^{(0)}*Z)'=(T-T*M_{\infty}*F^{(1)}*Z)'
                                                                        %                  =(T-T*P_{\infty*Z'*F^{(1)}*Z)'=(T-T*Kinf*Z)'
                                                                        %                  = (T*(I-*Kinf*Z))'=(I-Z'*Kinf')*T'
                                                                        %P_{*}=T*P_{\infty}*L^{(1)}+T*P_{*}*L^{(0)}+RQR
                                                                        %     =T*[(P_{\infty}*(-K^{(1)*Z}))+P_{*}*(I-Z'*Kinf')*T'+RQR]
        dlik(s)= log(det(Finf));                                        %set w_t to top case in bottom equation page 172, DK (2012)
        iFinf  = inv(Finf);
        Kinf   = Pinf*Z'*iFinf;                                         %define Kinf=T^{-1}*K_0 with M_{\infty}=Pinf*Z'
        Fstar  = Z*Pstar*Z' + H;                                        %(5.7) DK(2012)
        Kstar  = (Pstar*Z'-Kinf*Fstar)*iFinf;                           %(5.12) DK(2012); note that there is a typo in DK (2003) with "+ Kinf" instead of "- Kinf", but it is correct in their appendix
        Pstar  = T*(Pstar-Pstar*Z'*Kinf'-Pinf*Z'*Kstar')*T'+QQ;         %(5.14) DK(2012)
        Pinf   = T*(Pinf-Pinf*Z'*Kinf')*T';                             %(5.14) DK(2012)
        a      = T*(a+Kinf*v);                                          %(5.13) DK(2012)
    end
    t = t+1;
end

if t>last
    warning('kalman_filter_d: There isn''t enough information to estimate the initial conditions of the nonstationary variables. The diffuse Kalman filter never left the diffuse stage.');
    dLIK = NaN;
    return
end

dlik = dlik(1:s);
dlik = .5*(dlik + pp*log(2*pi));

dLIK = sum(dlik(1+presample:end));
