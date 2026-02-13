function [dLIK, dlikk, a, Pstar, llik] = univariate_kalman_filter_d(data_index, Y, start, last, a, Pinf, Pstar, kalman_tol, diffuse_kalman_tol, presample, T, R, Q, H, Z, pp)
% [dLIK, dlikk, a, Pstar, llik] = univariate_kalman_filter_d(data_index, Y, start, last, a, Pinf, Pstar, kalman_tol, diffuse_kalman_tol, presample, T, R, Q, H, Z, pp)
% Computes the diffuse log-likelihood of a state space model using the
% univariate filter approach.
%
% The function implements the exact diffuse Kalman filter in its univariate
% formulation, processing observations one at a time within each period t.
% This avoids the inversion of the (pp x pp) forecast error variance matrix
% and is numerically more robust than the multivariate version (kalman_filter_d),
% in particular when F_{inf,t} is rank-deficient but not zero (a case where
% the multivariate filter cannot proceed). The univariate approach also
% naturally handles missing observations via data_index.
%
% The state space model is given by:
%   y_t = Z * alpha_t + epsilon_t,           epsilon_t ~ N(0, H)
%   alpha_{t+1} = T * alpha_t + R * eta_t,   eta_t ~ N(0, Q)
%
% where H is assumed to be diagonal (a requirement for the univariate
% approach). Each scalar observation y_{t,i} is processed sequentially
% with scalar forecast error variances:
%   F_{inf,t,i} = Z_i * P_{inf,t,i} * Z_i'
%   F_{*,t,i}   = Z_i * P_{*,t,i} * Z_i' + H_i
%
% where Z_i denotes the i-th row of Z and H_i the i-th diagonal element
% of H. The state vector and covariance matrices are updated after each
% scalar observation (indexed by i) within a period, with the prediction
% step (involving T, R, Q) applied only at the end of each period t.
%
% with the diffuse initialization:
%   alpha_1 ~ N(a, Pinf * kappa + Pstar)  as kappa -> infinity
%
% The filter distinguishes three cases for each scalar observation:
%   (i)   F_{inf,t,i} > 0: the observation resolves diffuse uncertainty.
%         The likelihood contribution is log(F_{inf,t,i}) + log(2*pi)
%         (no quadratic term). See upper case on p. 175 of DK (2012).
%   (ii)  F_{inf,t,i} = 0 and F_{*,t,i} > 0: diffuse uncertainty does not
%         affect this observable. The standard likelihood contribution applies:
%         log(F_{*,t,i}) + v_{t,i}^2 / F_{*,t,i} + log(2*pi)
%   (iii) Both F_{inf,t,i} and F_{*,t,i} are zero (or below tolerance):
%         the observation is uninformative and skipped (a_{t,i+1} = a_{t,i},
%         P_{t,i+1} = P_{t,i}), see p. 157 of DK (2012).
%
% The diffuse phase ends once P_{inf,t} (and hence Z*Pinf*Z') has
% converged to zero, after which the standard Kalman filter takes over.
%
% INPUTS
% - data_index              [cell]      1*T cell of column vectors of indices (in the vector of observed variables)
% - Y                       [matrix]    pp*T matrix of doubles, data
% - start                   [integer]   index of the first observation to process in Y
% - last                    [integer]   index of the last observation to process in Y
% - a                       [vector]    mm*1 vector of doubles, initial mean of the state vector, E_0(alpha_1),
%   Pinf                    [double]    (m*m) diffuse part of the initial state covariance matrix;
%                                       reflects prior uncertainty about nonstationary components
%   Pstar                   [double]    (m*m) stationary part of the initial state covariance matrix;
%                                       reflects prior uncertainty about stationary components
% - kalman_tol              [double]    tolerance parameter (rcond, invertibility of the covariance matrix of the prediction errors)
% - diffuse_kalman_tol      [double]    tolerance parameter for diffuse filter
% - presample               [integer]   number of initial iterations to be discarded when evaluating the likelihood
% - T                       [matrix]    transition matrix of the state equation
% - R                       [matrix]    matrix relating the structural innovations to the state variables
% - Q                       [matrix]    covariance matrix of the structural innovations
% - H                       [vector]    diagonal of covariance matrix of the measurement errors
% - Z                       [matrix]    matrix relating the states to the observed variables
% - pp                      [integer]   number of observed variables
%
% OUTPUTS
%   dLIK                    [double]    scalar, minus the diffuse log-likelihood (up to a constant)
%   dlikk                   [double]    (s x pp) matrix of log-likelihood contributions by period
%                                       and observable during the diffuse phase, where s is the
%                                       number of diffuse iterations; each element equals
%                                       0.5*(w_{t,i}) with w_{t,i} as defined on p. 175 of
%                                       DK (2012). Note: on output, dlikk contains the
%                                       observation-level contributions (same as llik)
%   a                       [double]    (m*1) estimated state vector at the end of the diffuse
%                                       phase, E_{t_d}(alpha_{t_d+1}), to be used as initial
%                                       condition for the standard Kalman filter
%   Pstar                   [double]    (m*m) state covariance matrix at the end of the diffuse
%                                       phase, Var_{t_d}(alpha_{t_d+1}), to be used as initial
%                                       condition for the standard Kalman filter
%   llik                    [double]    (s*pp) matrix of log-likelihood contributions by period
%                                       (row) and observable (column); non-observed entries are
%                                       zero. Each non-zero element equals
%                                       0.5*(w_{t,i}) with w_{t,i} as on p. 175 of DK (2012)
%
% This function is called by: dsge_likelihood
%
% Algorithm:
%   Uses the diffuse univariate filter as described in Durbin/Koopman (2012): "Time
%   Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 5, 6.4 + 7.2.5

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
dlikk= zeros(smpl,1);      % Initialization of the vector gathering the densities.
llik = zeros(smpl,pp);

newRank = rank(Pinf,diffuse_kalman_tol);
l2pi = log(2*pi);
s=0;

while newRank && (t<=last)
    s = t-start+1;
    d_index = data_index{t};
    for i=1:length(d_index)
        Zi = Z(d_index(i),:);
        prediction_error = Y(d_index(i),t) - Zi*a;      % nu_{t,i} in 6.13 in DK (2012)
        Fstar = Zi*Pstar*Zi' + H(d_index(i));           % F_{*,t} in 5.7 in DK (2012), relies on H being diagonal
        Finf  = Zi*Pinf*Zi';                            % F_{\infty,t} in 5.7 in DK (2012), relies on H being diagonal
        Kstar = Pstar*Zi';
        % Conduct check of rank
        % Pinf and Finf are always scaled such that their norm=1: Fstar/Pstar, instead,
        % depends on the actual values of std errors in the model and can be badly scaled.
        % experience is that diffuse_kalman_tol has to be bigger than kalman_tol, to ensure
        % exiting the diffuse filter properly, avoiding tests that provide false non-zero rank for Pinf.
        % Also the test for singularity is better set coarser for Finf than for Fstar for the same reason
        if Finf>diffuse_kalman_tol && newRank           % F_{\infty,t,i} = 0, use upper part of bracket on p. 175 DK (2012) for w_{t,i}
            Kinf   = Pinf*Zi';
            Kinf_Finf = Kinf/Finf;
            a         = a + Kinf_Finf*prediction_error;
            Pstar     = Pstar + Kinf*(Kinf_Finf'*(Fstar/Finf)) - Kstar*Kinf_Finf' - Kinf_Finf*Kstar';
            Pinf      = Pinf - Kinf*Kinf_Finf';
            llik(s,d_index(i)) = log(Finf) + l2pi;
            dlikk(s) = dlikk(s) + llik(s,d_index(i));
        elseif Fstar>kalman_tol
            llik(s,d_index(i)) = log(Fstar) + (prediction_error*prediction_error/Fstar) + l2pi;
            dlikk(s) = dlikk(s) + llik(s,d_index(i));
            a = a+Kstar*(prediction_error/Fstar);
            Pstar = Pstar-Kstar*(Kstar'/Fstar);
        else
            if Fstar<0 || Finf<0
                %pathological numerical case where variance is negative
                dLIK = NaN;
                return
            else
                % do nothing as a_{t,i+1}=a_{t,i} and P_{t,i+1}=P_{t,i}, see
                % p. 157, DK (2012)
            end
        end
    end
    if newRank
        oldRank = rank(Z*Pinf*Z',diffuse_kalman_tol);
    else
        oldRank = 0;
    end
    a     = T*a;
    if isqvec
        QQ = R*Qvec(:,:,t+1)*transpose(R);
    end
    Pstar = T*Pstar*T'+QQ;
    Pinf  = T*Pinf*T';
    if newRank
        newRank = rank(Z*Pinf*Z',diffuse_kalman_tol);
    end
    if oldRank ~= newRank
        disp('univariate_diffuse_kalman_filter:: T does influence the rank of Pinf!')
        disp('This may happen for models with order of integration >1.')
    end
    t = t+1;
end

if (t>last)
    warning('univariate_diffuse_kalman_filter:: There isn''t enough information to estimate the initial conditions of the nonstationary variables');
    dLIK = NaN;
    return
end

% Divide by two.
dlikk = .5*dlikk(1:s);
llik  = .5*llik(1:s,:);

dLIK = sum(dlikk(1+presample:end));
dlikk = llik;