function [p, err, funevals] = mvtcdfqmc(a, b, Rho, nu, tol, maxfunevals, verbose)
% [p, err, funevals] = mvtcdfqmc(a, b, Rho, nu)
% [p, err, funevals] = mvtcdfqmc(a, b, Rho, nu, tol)
% [p, err, funevals] = mvtcdfqmc(a, b, Rho, nu, tol, maxfunevals)
% [p, err, funevals] = mvtcdfqmc(a, b, Rho, nu, tol, maxfunevals, verbose)
% -------------------------------------------------------------------------
% Quasi-Monte Carlo computation of multivariate Student's t CDF.
% -------------------------------------------------------------------------
% p = mvtcdfqmc(a, b, Rho, nu) computes the cumulative probability of the
% multivariate Student's t distribution with correlation matrix Rho and
% degrees of freedom nu, evaluated over the hyper-rectangle with lower
% limits a and upper limits b.
%
% For multivariate normal (MVN) integration, set nu = Inf.
%
% p = mvtcdfqmc(..., tol) specifies the absolute error tolerance (default 1e-4).
%
% p = mvtcdfqmc(..., tol, maxfunevals) specifies the maximum number of
% integrand evaluations (default 1e7).
%
% p = mvtcdfqmc(..., tol, maxfunevals, verbose) controls display output:
%   0 = 'off'   - no output (default)
%   1 = 'final' - display the final probability and related error after the
%                 integrand has converged successfully
%   2 = 'iter'  - display the probability and estimated error at each repetition
%
% [p, err, funevals] = mvtcdfqmc(...) also returns an estimate of the
% error and the number of function evaluations used.
% -------------------------------------------------------------------------
% Inputs:
% a           [D-by-1 vector]   Lower integration limits
% b           [D-by-1 vector]   Upper integration limits
% Rho         [D-by-D matrix]   Correlation matrix (positive definite)
% nu          [scalar]          Degrees of freedom (Inf for MVN)
% tol         [scalar]          Absolute error tolerance (optional, default 1e-4)
% maxfunevals [scalar]          Maximum function evaluations (optional, default 1e7)
% verbose     [scalar]          Display level: 0=off, 1=final, 2=iter (optional, default 0)
% -------------------------------------------------------------------------
% Outputs:
% p        [scalar]   Probability estimate in [0, 1]
% err      [scalar]   Error estimate
% funevals [scalar]   Number of function evaluations
% -------------------------------------------------------------------------
% Example:
% % 4-dimensional MVN probability
% Rho = [1.0 0.3 0.2 0.1; 0.3 1.0 0.4 0.2; 0.2 0.4 1.0 0.3; 0.1 0.2 0.3 1.0];
% a = [-Inf; -Inf; -Inf; -Inf];
% b = [1; 2; 3; 4];
% [p, err] = mvtcdfqmc(a, b, Rho, Inf)
% -------------------------------------------------------------------------
% References:
% Genz, A. and F. Bretz (1999) "Numerical Computation of Multivariate
%   t Probabilities with Application to Power Calculation of Multiple Contrasts."
%   Journal of Statistical Computation and Simulation 63(4), pages 361-378.
%   doi: 10.1080/00949659908811962
% Genz, A. and F. Bretz (2002) "Comparison of Methods for the Computation
%   of Multivariate t Probabilities."
%   Journal of Computational and Graphical Statistics 11(4), pages 950-971.
%   doi: 10.1198/106186002394
%
% Code adapted from the GNU Octave statistics package (version 1.7.7):
% https://github.com/gnu-octave/statistics/blob/main/inst/dist_fun/mvtcdfqmc.m
%
% See also mvncdf.

% Copyright © 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
% Copyright © 2025 Dynare Team
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

% input validation
if nargin < 4
    error('mvtcdfqmc: requires at least 4 input arguments (a, b, Rho, nu)');
end

% set defaults
if nargin < 5 || isempty(tol)
    tol = 1e-4;
end
if nargin < 6 || isempty(maxfunevals)
    maxfunevals = 1e7;
end
if nargin < 7 || isempty(verbose)
    verbose = 0;
end

% validate optional arguments
if ~isscalar(tol) || ~isreal(tol) || tol <= 0
    error('mvtcdfqmc: tol must be a positive scalar');
end
if ~isscalar(maxfunevals) || ~isreal(maxfunevals) || maxfunevals <= 0
    error('mvtcdfqmc: maxfunevals must be a positive scalar');
end
maxfunevals = floor(maxfunevals);
if ~isscalar(verbose) || ~ismember(verbose, [0 1 2])
    error('mvtcdfqmc: verbose must be 0, 1, or 2');
end

% ensure column vectors
a = a(:);
b = b(:);

% check for consistent limits
if ~all(a < b)
    if any(a > b)
        error('mvtcdfqmc: lower limits a must be less than upper limits b');
    elseif any(isnan(a)) || any(isnan(b))
        warning('mvtcdfqmc: NaN values in integration limits');
        p = NaN;
        err = NaN;
        funevals = 0;
        return
    else
        % a == b for some dimensions (degenerate case)
        warning('mvtcdfqmc: zero distance between lower upper limits (degenerate case)')
        p = 0;
        err = 0;
        funevals = 0;
        return
    end
end

% dimensions with limits (-Inf, Inf) can be ignored
dblInfLims = (a == -Inf) & (b == Inf);
if any(dblInfLims)
    if all(dblInfLims)
        % all dimensions have infinite limits: probability is 1
        warning('mvtcdfqmc: infinite distance between lower upper limits')
        p = 1;
        err = 0;
        funevals = 0;
        return
    end
    % remove dimensions with infinite limits
    a(dblInfLims) = [];
    b(dblInfLims) = [];
    Rho(:,dblInfLims) = [];
    Rho(dblInfLims,:) = [];
end

% get dimension
m = size(Rho, 1);

% sort integration order by increasing interval length (improves convergence)
[~, ord] = sort(b - a);
a = a(ord);
b = b(ord);
Rho = Rho(ord, ord);

% warn about highly correlated matrices
if any(any(abs(tril(Rho, -1)) > 0.999))
    warning('mvtcdfqmc: correlation matrix has correlations > 0.999');
end

% Cholesky factorization and scaling
C = chol(Rho);
c = diag(C);
a = a ./ c;
b = b ./ c;
C = C ./ repmat(c', m, 1);

% Monte Carlo settings
MCreps = 25;  % number of Monte Carlo replications for variance estimation
MCdims = m - isinf(nu);  % for MVN (nu=Inf), one fewer dimension needed

% initialize outputs
p = 0;
sigsq = Inf;
funevals = 0;
err = NaN;

% prime sequence for adaptive refinement
P = [31 47 73 113 173 263 397 593 907 1361 2053 3079 4621 6947 10427 15641 23473 ...
     35221 52837 79259 118891 178349 267523 401287 601942 902933 1354471 2031713];

% display header for verbose mode
if verbose > 1
    fprintf('mvtcdfqmc: Probability    Error       FunEvals\n');
    fprintf('mvtcdfqmc: -----------    ---------   --------\n');
end

% adaptive quasi-Monte Carlo integration
for i = 5:length(P)
    % check if we would exceed maxfunevals
    if (funevals + 2*MCreps*P(i)) > maxfunevals
        break
    end

    % Niederreiter point set generator
    q = 2.^((1:MCdims) / (MCdims + 1));

    % compute randomized quasi-Monte Carlo estimate with P points
    [THat, sigsqTHat] = estimate_mvtqmc(MCreps, P(i), q, C, nu, a, b);
    funevals = funevals + 2*MCreps*P(i);

    % recursively update estimate and error estimate
    p = p + (THat - p) / (1 + sigsqTHat / sigsq);
    sigsq = sigsqTHat / (1 + sigsqTHat / sigsq);

    % conservative error estimate: 3.5 times the MC standard error
    err = 3.5 * sqrt(sigsq);

    % display iteration info
    if verbose > 1
        fprintf('mvtcdfqmc: %.5g       %.5e   %d\n', p, err, funevals);
    end

    % check convergence
    if err < tol
        if verbose > 0
            fprintf('mvtcdfqmc: Converged. Probability = %.6g, Error = %.4e, FunEvals = %d\n', ...
                    p, err, funevals);
        end
        return
    end
end

% did not converge
warning('mvtcdfqmc: did not converge to tolerance %.4e within %d function evaluations', ...
        tol, maxfunevals);


return % --*-- Unit tests --*--

%@test:1
% test MVN at origin: P(X <= 0) = 0.5^D for independent variables
try
    Rho = eye(4);
    a = -Inf(4, 1);
    b = zeros(4, 1);
    [p, err] = mvtcdfqmc(a, b, Rho, Inf);
    t(1) = abs(p - 0.0625) < 1e-3; % 0.5^4 = 0.0625
    t(2) = err < 1e-3;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:1

%@test:2
% test MVN with correlation
try
    Rho = [1 0.5; 0.5 1];
    a = [-Inf; -Inf];
    b = [0; 0];
    [p, err] = mvtcdfqmc(a, b, Rho, Inf);
    p_expected = 0.25 + asin(0.5) / (2*pi); % analytical result
    t(1) = abs(p - p_expected) < 1e-3;
catch
    t = false;
end
T = all(t);
%@eof:2

%@test:3
% test degenerate case: a == b should give p = 0
try
    Rho = eye(3);
    a = [0; 0; 0];
    b = [0; 0; 0];
    [p, err, funevals] = mvtcdfqmc(a, b, Rho, Inf);
    t(1) = p == 0;
    t(2) = err == 0;
    t(3) = funevals == 0;
catch
    t = false(3, 1);
end
T = all(t);
%@eof:3

%@test:4
% test infinite limits: all dimensions (-Inf, Inf) should give p = 1
try
    Rho = [1 0.3; 0.3 1];
    a = [-Inf; -Inf];
    b = [Inf; Inf];
    [p, err, funevals] = mvtcdfqmc(a, b, Rho, Inf);
    t(1) = p == 1;
    t(2) = err == 0;
    t(3) = funevals == 0;
catch
    t = false(3, 1);
end
T = all(t);
%@eof:4

%@test:5
% test MVT with finite degrees of freedom
try
    Rho = eye(4);
    a = -Inf(4, 1);
    b = zeros(4, 1);
    nu = 5; % degrees of freedom
    [p, err] = mvtcdfqmc(a, b, Rho, nu);
    % for independent t variables: P(T1<=0, ..., T4<=0) = 0.5^4 = 0.0625
    t(1) = abs(p - 0.0625) < 1e-3;
    t(2) = err < 1e-3;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:5

%@test:6
% test error is returned
try
    Rho = [1.0 0.3 0.2; 0.3 1.0 0.4; 0.2 0.4 1.0];
    a = [-Inf; -Inf; -Inf];
    b = [1; 2; 3];
    [p, err] = mvtcdfqmc(a, b, Rho, Inf);
    t(1) = ~isnan(err);
    t(2) = err >= 0;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:6

%@test:7
% test funevals is returned
try
    Rho = eye(5);
    a = -Inf(5, 1);
    b = ones(5, 1);
    [~, ~, funevals] = mvtcdfqmc(a, b, Rho, Inf);
    t(1) = funevals > 0;
catch
    t = false;
end
T = all(t);
%@eof:7

%@test:8
% test custom tolerance
try
    Rho = [1 0.5; 0.5 1];
    a = [-Inf; -Inf];
    b = [1; 1];
    [p1, err1] = mvtcdfqmc(a, b, Rho, Inf, 1e-3);
    [p2, err2] = mvtcdfqmc(a, b, Rho, Inf, 1e-5);
    % both should give similar results
    t(1) = abs(p1 - p2) < 1e-2;
    % tighter tolerance should have smaller error
    t(2) = err2 <= err1 || err2 < 1e-4;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:8

%@test:9
% test probability is in [0, 1]
try
    Rho = [1.0 0.8 0.6; 0.8 1.0 0.7; 0.6 0.7 1.0];
    a = [-2; -1; 0];
    b = [1; 2; 3];
    [p, ~] = mvtcdfqmc(a, b, Rho, Inf);
    t(1) = p >= 0 && p <= 1;
catch
    t = false;
end
T = all(t);
%@eof:9

%@test:10
% test error: a > b should error
try
    Rho = eye(2);
    a = [1; 0];
    b = [0; 1]; % a(1) > b(1)
    mvtcdfqmc(a, b, Rho, Inf);
    t(1) = false; % should not reach here
catch
    t(1) = true;
end
T = all(t);
%@eof:10

end % mvtcdfqmc


function [THat, sigsqTHat] = estimate_mvtqmc(MCreps, P, q, C, nu, a, b)
% Randomized quasi-Monte Carlo estimate of the integral

qq = (1:P)' * q;
THat = zeros(MCreps, 1);

for rep = 1:MCreps
    % generate random lattice shift
    % for MVT this is the m-dimensional unit hypercube
    % for MVN this is the (m-1)-dimensional unit hypercube
    w = abs(2 * mod(qq + repmat(rand(size(q)), P, 1), 1) - 1);

    % compute mean of integrand over all P points and antithetic points
    THat(rep) = (F_qrsvn(a, b, C, nu, w) + F_qrsvn(a, b, C, nu, 1 - w)) / 2;
end

% return MC mean and variance
sigsqTHat = var(THat) / MCreps;
THat = mean(THat);

end % estimate_mvtqmc


function TBar = F_qrsvn(a, b, C, nu, w)
% Integrand for computation of MVT/MVN probabilities
%
% Given box bounds a and b (may be infinite), Cholesky factor C of
% correlation matrix, and degrees of freedom nu, compute the transformed
% integrand at each row of quasi-random lattice w.
% For MVN (nu=Inf), w is in (m-1)-dimensional unit hypercube.
% For MVT, w is in m-dimensional unit hypercube.

N = size(w, 1);  % number of quasi-random points
m = length(a);   % number of dimensions

% compute scaling factor from chi distribution (for MVT)
if isinf(nu)
    snu = 1;
else
    snu = chi_inv(w(:, m), nu) / sqrt(nu);
end

% first dimension
d = norm_cdf(snu .* a(1));        % a is already scaled by diag(C)
emd = norm_cdf(snu .* b(1)) - d;  % b is already scaled by diag(C)
T = emd;
y = zeros(N, m);

% remaining dimensions
for i = 2:m
    % clamp z to avoid infinite values
    z = min(max(d + emd .* w(:, i-1), eps/2), 1 - eps/2);
    y(:, i-1) = norm_inv(z);
    ysum = y * C(:, i);
    d = norm_cdf(snu .* a(i) - ysum);        % a is already scaled by diag(C)
    emd = norm_cdf(snu .* b(i) - ysum) - d;  % b is already scaled by diag(C)
    T = T .* emd;
end

TBar = sum(T) / length(T);

end % F_qrsvn


function p = norm_cdf(z)
% Normal cumulative distribution function (CDF)
p = 0.5 * erfc(-z / sqrt(2));
end % norm_cdf


function z = norm_inv(p)
% Inverse of normal cumulative distribution function (CDF)
z = -sqrt(2) * erfcinv(2 * p);
end % norm_inv


function x = chi_inv(p, nu)
% Inverse chi cumulative distribution function CDF (not chi-squared)
x = sqrt(gammaincinv(p, nu / 2) * 2);
end % chi_inv

