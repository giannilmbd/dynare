function [LIK, LIKK, a, P] = kalman_filter_fast(Y,start,last,a,P,kalman_tol,presample,T,H,Z,pp,Zflag,diffuse_periods)
% [LIK, LIKK, a, P] = kalman_filter_fast(Y,start,last,a,P,kalman_tol,presample,T,H,Z,pp,Zflag,diffuse_periods)
% computes the likelihood of a stationary state space model using the fast
% implementation of the Kalman filter of Edward Herbst, (2015). Using the 'Chandrasekhar Recursions'
% for likelihood evaluation of DSGE models. Computational Economics, 45(4):693–705
% Useful for models where the number of states is much greater than the
% number of observables.
%
% Inputs:
% - Y           [double]        pp*T matrix of data
% - scalar      [integer]       first period.
% - last        [integer]       last period, last- first has to be smaller than T
% - a           [double]        mm*1 vector of predicted initial state variables E_{0}(alpha_1)
% - P           [double]        mm*mm initial covariance matrix of the state vector
% - kalman_tol  [double]        tolerance parameter (rcond, invertibility of the covariance matrix of the prediction errors).
% - presample   [integer]       number of initial iterations to be discarded when evaluating the likelihood
% - T           [double]        mm*mm state transition matrix
% - H           [double]        pp*pp covariance matrix of the measurement errors (if no measurement errors set H as a zero scalar).
% - Z           [double]        pp*mm matrix matrix relating states to  observed variables or vector of indices (depending on the value of Zflag).
% - pp          [integer]       number of observed variables
% - Zflag       [boolean]       equal to 0 if Z is a vector of indices targeting the observed variables in the state vector,
%                               equal to 1 if Z is a pp*mm} matrix
% - diffuse_periods [integer]       number of diffuse filter periods in the initialization step
%
% Outputs
% - LIK         [double]        scalar value of (minus) the likelihood
% - likk        [double]        column vector of density of each observation
% - a           [double]        mm*1 mean of the state vector at the end of the (sub)sample (E_{T}(alpha_{T+1})).
% - P           [double]        mm*mm covariance of the state vector at the end of the (sub)sample
%
% This function is called by: dsge_likelihood

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
if nargin<12
    Zflag = 0;
end

if nargin<13
    diffuse_periods = 0;
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
t    = start;              % Initialization of the time index.
likk = zeros(smpl,1);      % Initialization of the vector gathering the densities.
LIK  = Inf;                % Default value of the log likelihood.
F_singular  = 1;
LIKK=[];

if Zflag
    K = T*P*Z';
    F = Z*P*Z' + H;
else
    K = T*P(:,Z);
    F = P(Z,Z) + H;
end
W = K;
iF = inv(F);
Kg = K*iF;
M = -iF;

while t<=last
    s = t-start+1;
    if Zflag
        v  = Y(:,t)-Z*a;
    else
        v  = Y(:,t)-a(Z);
    end
    if rcond(F) < kalman_tol
        % if ~all(abs(F(:))<kalman_tol), then use univariate diffuse filter, otherwise
        % this is a pathological case and the draw is discarded
        return
    else
        F_singular = 0;
        dF      = det(F);
        likk(s) = log(dF)+transpose(v)*iF*v;
        a = T*a+Kg*v;
        if Zflag
            ZWM = Z*W*M;
            ZWMWp = ZWM*W';
            M = M + ZWM'*iF*ZWM;
            F  = F + ZWMWp*Z';
            iF      = inv(F);
            K = K + T*ZWMWp';
            Kg = K*iF;
            W = (T - Kg*Z)*W;
        else
            ZWM = W(Z,:)*M;
            ZWMWp = ZWM*W';
            M = M + ZWM'*iF*ZWM;
            F  = F + ZWMWp(:,Z);
            iF      = inv(F);
            K = K + T*ZWMWp';
            Kg = K*iF;
            W = T*W - Kg*W(Z,:);
        end
    end
    t = t+1;
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end

% Add observation's densities constants and divide by two.
likk(1:s) = .5*(likk(1:s) + pp*log(2*pi));

% Compute minus the log-likelihood.
if presample>diffuse_periods
    LIK = sum(likk(1+(presample-diffuse_periods):end));
else
    LIK = sum(likk);
end
