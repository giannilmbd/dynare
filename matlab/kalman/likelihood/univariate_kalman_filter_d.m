function [dLIK, dlikk, a, Pstar, llik] = univariate_kalman_filter_d(data_index, Y, start, last, a, Pinf, Pstar, kalman_tol, diffuse_kalman_tol, presample, T, R, Q, H, Z, pp)
% [dLIK, dlikk, a, Pstar, llik] = univariate_kalman_filter_d(data_index, Y, start, last, a, Pinf, Pstar, kalman_tol, diffuse_kalman_tol, presample, T, R, Q, H, Z, pp)
% Computes the likelihood of a state space model (initialization with diffuse steps, univariate approach).
%
% INPUTS
% - data_index              [cell]      1*T cell of column vectors of indices (in the vector of observed variables)
% - Y                       [matrix]    pp*T matrix of doubles, data
% - start                   [integer]   first period
% - last                    [integer]   last period
% - a                       [vector]    mm*1 vector of doubles, initial mean of the state vector
% - Pinf                    [matrix]    initial covariance matrix of the state vector (non stationary part)
% - Pstar                   [matrix]    initial covariance matrix of the state vector (stationary part)
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
% - dLIK                    [double]    value of (minus) the likelihood for the first d observations
% - dlikk                   [vector]    d*pp matrix of doubles, log-likelihood values
% - a                       [vector]    mm*1 vector of doubles, mean of the state vector at the end of the (sub)sample
% - Pstar                   [matrix]    mm*mm matrix of doubles, covariance of the state vector at the end of the subsample
% - llik                    [matrix]    d*pp matrix of doubles, log-likelihood contributions by observation
%
% This function is called by: dsge_likelihood
%
% Algorithm:
%   Uses the diffuse univariate filter as described in Durbin/Koopman (2012): "Time
%   Series Analysis by State Space Methods", Oxford University Press,
%   Second Edition, Ch. 5, 6.4 + 7.2.5

% Copyright © 2004-2024 Dynare Team
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
dF   = 1;
isqvec = false;
if ndims(Q)>2
    Qvec = Q;
    Q=Q(:,:,1);
    isqvec = true;
end
QQ   = R*Q*transpose(R);   % Variance of R times the vector of structural innovations.
t    = start;              % Initialization of the time index.
dlikk= zeros(smpl,1);      % Initialization of the vector gathering the densities.
dLIK = Inf;                % Default value of the log likelihood.
oldK = Inf;
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