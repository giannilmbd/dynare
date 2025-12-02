function [p, err] = mvncdf(varargin)
% [p, err] = mvncdf(XU)
% [p, err] = mvncdf(XU, MU, SIGMA)
% [p, err] = mvncdf(XL, XU, MU, SIGMA)
% [p, err] = mvncdf(..., tol)
% -------------------------------------------------------------------------
% Multivariate normal cumulative distribution function (cdf).
% -------------------------------------------------------------------------
% p = mvncdf(XU) returns the cumulative probability of the multivariate
% normal distribution with zero mean and identity covariance matrix,
% evaluated at each row of XU. Rows of the N-by-D matrix XU correspond to
% observations or points, and columns correspond to variables or coordinates.
% p is an N-by-1 vector.
%
% p = mvncdf(XU, MU, SIGMA) returns the cumulative probability of the
% multivariate normal distribution with mean MU and covariance SIGMA,
% evaluated at each row of XU. MU is a 1-by-D vector, and SIGMA is a D-by-D
% symmetric, positive definite matrix. MU can also be a scalar value, which
% is replicated to match the size of XU. Pass in the empty matrix for MU to
% use its default value (zero mean) when you want to only specify SIGMA.
% If the covariance matrix is diagonal, SIGMA may also be specified as a
% 1-by-D vector containing just the diagonal (variances).
% This computes Pr{X(1)<=XU(1), X(2)<=XU(2), ..., X(D)<=XU(D)} for each row.
%
% p = mvncdf(XL, XU, MU, SIGMA) returns the multivariate normal cumulative
% probability evaluated over the rectangle with lower limit XL and upper
% limit XU. XL and XU must have the same size (N-by-D).
% This computes Pr{XL(1)<=X(1)<=XU(1), ..., XL(D)<=X(D)<=XU(D)} for each row.
%
% p = mvncdf(..., tol) specifies the absolute error tolerance for the
% quasi-Monte Carlo integration (default 1e-4). Only used for D >= 4.
%
% [p, err] = mvncdf(...) returns an estimate of the error in p. For
% univariate distributions, diagonal covariance matrices, and bivariate
% distributions, err is NaN. For trivariate distributions (D=3), err is the
% quadrature tolerance (default 1e-8). For higher dimensional cases (D >= 4),
% err is the estimated error from the adaptive quasi-Monte Carlo integration.
% -------------------------------------------------------------------------
% Inputs:
% XU      [N-by-D matrix]   Upper integration limits (rows=observations, cols=variables)
% XL      [N-by-D matrix]   Lower integration limits (optional, default -Inf)
% MU      [1-by-D vector]   Mean vector (or D-by-1, or scalar, or empty for zero mean)
% SIGMA   [D-by-D matrix]   Covariance matrix (positive semi-definite), or 1-by-D diagonal
% tol     [scalar]          Absolute error tolerance for QMC integration (optional, default 1e-4)
% -------------------------------------------------------------------------
% Outputs:
% p     [N-by-1 vector]   Probability values in [0, 1]
% err   [N-by-1 vector]   Error estimates (NaN for analytical cases)
% -------------------------------------------------------------------------
% Examples:
% % Standard normal (zero mean, identity covariance)
% XU = [0.5, 0.3; 1.0, 1.0; -0.5, 0.2];
% p = mvncdf(XU)
%
% % Single observation (column vector input)
% mu = [1; -1];
% Sigma = [0.9 0.4; 0.4 0.3];
% xu = [2; 0];
% p = mvncdf(xu, mu, Sigma)
%
% % Multiple observations (matrix input, each row is an observation)
% mu = [1, -1];
% Sigma = [0.9 0.4; 0.4 0.3];
% [X1, X2] = meshgrid(linspace(-1, 3, 25)', linspace(-3, 1, 25)');
% XU = [X1(:), X2(:)];
% p = mvncdf(XU, mu, Sigma);
% surf(X1,X2,reshape(p,25,25));
%
% % Using empty MU for zero mean and diagonal SIGMA as vector
% XU = [1, 2, 3];
% p = mvncdf(XU, [], [1, 4, 9])  % variances: 1, 4, 9
% -------------------------------------------------------------------------
% References:
% Drezner, Z. (1994) "Computation of the Trivariate Normal Integral",
%   Mathematics of Computation 62(205), pages 289-294. doi:10.2307/2153409
% Genz, A. (1992) "Numerical Computation of Multivariate Normal Probabilities."
%   Journal of Computational and Graphical Statistics 1(2), pages 141-149.
%   doi: 10.1080/10618600.1992.10477010
% Genz, A. (2004) "Numerical Computation of Rectangular Bivariate and Trivariate
%   Normal and t Probabilities", Statistics and Computing 14, pages 251-260.
%   doi: 10.1023/B:STCO.0000035304.20635.31
% Genz, A. and F. Bretz (1999) "Numerical Computation of Multivariate
%   t-Probabilities with Application to Power Calculation of Multiple Contrasts."
%   Journal of Statistical Computation and Simulation 63(4), pages 103-117.
%   doi: 10.1080/00949659908811962
% Genz, A. and F. Bretz (2002) "Comparison of Methods for the Computation
%   of Multivariate t Probabilities."
%   Journal of Computational and Graphical Statistics 11(4), pages 950-971.
%   doi: 10.1198/106186002394
%
% Code adapted from the GNU Octave statistics package (version 1.7.7):
% https://github.com/gnu-octave/statistics/blob/main/inst/dist_fun/mvncdf.m
%
% See also normcdf, bvncdf, tvncdf, mvtcdfqmc.

% Copyright © 2008 Arno Onken <asnelt@asnelt.org>
% Copyright © 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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

% default tolerance for quasi-Monte Carlo integration (D >= 4)
tol_qmc = 1e-4;
% parse inputs
if nargin < 1
    error('mvncdf: Expected at least 1 input argument');
elseif nargin == 1
    % mvncdf(XU) - zero mean, identity covariance
    upperLimitOnly = true;
    XU = varargin{1};
    mu = [];
    Sigma = [];
    XL = -inf(size(XU));
elseif nargin == 2 && isscalar(varargin{2})
    % mvncdf(XU, tol) - zero mean, identity covariance, with tol specified
    upperLimitOnly = true;
    XU = varargin{1};
    mu = [];
    Sigma = [];
    XL = -inf(size(XU));
    tol_qmc = varargin{2};
elseif nargin == 3 || (nargin == 4 && isscalar(varargin{4}))
    % mvncdf(XU, MU, SIGMA) or mvncdf(XU, MU, SIGMA, tol)
    upperLimitOnly = true;
    XU = varargin{1};
    mu = varargin{2};
    Sigma = varargin{3};
    XL = -inf(size(XU));
    if nargin == 4
        tol_qmc = varargin{4};
    end
elseif nargin == 4 || nargin == 5
    % mvncdf(XL, XU, MU, SIGMA) or mvncdf(XL, XU, MU, SIGMA, tol)
    upperLimitOnly = false;
    XL = varargin{1};
    XU = varargin{2};
    mu = varargin{3};
    Sigma = varargin{4};
    if nargin == 5
        tol_qmc = varargin{5};
    end
else
    error('mvncdf: Expected 1 to 5 input arguments');
end

% get size of data: rows = observations, columns = variables
[n, d] = size(XU);

% handle edge cases for empty inputs
if d == 0 % no dimensions: error
    error('mvncdf: X must have at least one column');
end
if n == 0
    % no observations but dimensions exist: return 1 (empty product over observations)
    p = 1;
    if nargout > 1
        err = NaN;
    end
    return;
end

% handle column vector input: interpret as single D-dimensional observation
if d == 1 && n > 1
    % special case: if Sigma is provided, use it to determine if XU should be transposed
    if ~isempty(Sigma)
        sz = size(Sigma);
        if (sz(1) == 1 && sz(2) == n) || (sz(1) == n && sz(2) == n)
            % Sigma suggests XU is a single observation with n dimensions
            XU = XU';
            XL = XL';
            [n, d] = size(XU);
        end
    else
        % no Sigma, so transpose column vector to row vector (single observation)
        XU = XU';
        XL = XL';
        [n, d] = size(XU);
    end
end

% handle empty MU: default to zero mean
if isempty(mu)
    mu = zeros(1, d);
elseif isscalar(mu)
    % scalar mu: replicate to all dimensions
    mu = repmat(mu, 1, d);
else
    % ensure mu is a row vector
    mu = mu(:)';
end

% handle empty Sigma: default to identity covariance
if isempty(Sigma)
    Sigma = eye(d);
end

% handle SIGMA as 1-by-D diagonal vector
sigmaIsDiagVector = false;
if size(Sigma, 1) == 1 && size(Sigma, 2) == d && d > 1
    % Sigma is provided as a row vector of variances, expand to full diagonal matrix
    sigmaIsDiagVector = true;
    Sigma = diag(Sigma);
end

% input validation
if size(XL, 1) ~= n || size(XL, 2) ~= d
    error('mvncdf: XL and XU must have the same size');
end
if length(mu) ~= d
    error('mvncdf: the length of MU must match the number of columns in XU');
end
if ~isequal(size(Sigma), [d d])
    error('mvncdf: the covariance matrix SIGMA must be a D-by-D matrix where D is the number of columns in XU');
end
if any(XL > XU, 'all')
    error('mvncdf: the lower integration limit XL must be less than or equal to the upper integration limit XU element-wise');
end

% dimension limit check
if d > 25
    error('mvncdf: number of dimensions D must not exceed 25');
end

% check if covariance matrix is diagonal
isDiagonal = sigmaIsDiagVector || isequal(Sigma, diag(diag(Sigma)));

% check for valid covariance matrix
if isDiagonal
    if any(diag(Sigma) <= 0)
        error('mvncdf: diagonal elements of covariance matrix SIGMA must be positive');
    end
else
    [~, cholErr] = chol(Sigma);
    if cholErr ~= 0
        error('mvncdf: the covariance matrix SIGMA must be positive semi-definite');
    end
end

% center the integration limits by subtracting mu (broadcast over rows)
XL0 = XL - mu;
XU0 = XU - mu;

% standardize: scale limits by standard deviations and convert covariance to correlation
s = sqrt(diag(Sigma))';  % row vector of standard deviations
XL0 = XL0 ./ s;
XU0 = XU0 ./ s;
if ~isDiagonal
    Rho = Sigma ./ (s' * s); % correlation matrix
end

if d == 1 % univariate case: use normcdf directly
    p = normcdf(XU0, 0, 1) - normcdf(XL0, 0, 1);
    if nargout > 1
        err = NaN(n, 1);
    end

elseif isDiagonal % diagonal covariance matrix: product of independent univariate CDFs
    p = prod(normcdf(XU0, 0, 1) - normcdf(XL0, 0, 1), 2);
    if nargout > 1
        err = NaN(n, 1);
    end

elseif d == 2 % bivariate case: use bvncdf (Gauss-Legendre quadrature)
    if upperLimitOnly
        % upper limit only: direct call to bvncdf (handles N observations)
        p = bvncdf(XU, mu, Sigma);
    else
        % both lower and upper limits: use inclusion-exclusion principle
        % P(XL <= X <= XU) = sum over subsets of (-1)^|S| * P(X <= x_S)
        % where x_S has XL for indices in S and XU otherwise
        
        % handle degenerate rectangles (XL == XU) by setting both to -Inf
        equalLimits = (XL == XU);
        XU(equalLimits) = -Inf;
        XL(equalLimits) = -Inf;
        
        p = zeros(n, 1);
        for i = 0:d
            k = nchoosek(1:d, i);
            for j = 1:size(k, 1)
                x = XU;
                x(:, k(j,:)) = XL(:, k(j,:));
                p = p + (-1)^i * bvncdf(x, mu, Sigma);
            end
        end
    end
    if nargout > 1
        err = NaN(n, 1); % bvncdf uses analytical quadrature, no error estimate
    end

elseif d == 3 % trivariate case: use tvncdf (adaptive quadrature)
    % extract correlation coefficients [rho_21, rho_31, rho_32] from correlation matrix
    % Rho([2 3 6]) in column-major order gives these elements
    rho_vec = Rho([2 3 6]);
    tol = 1e-8; % default tolerance for tvncdf
    
    if upperLimitOnly
        % upper limit only: direct call to tvncdf with standardized limits (handles N observations)
        p = tvncdf(XU0, rho_vec, tol);
    else
        % both lower and upper limits: use inclusion-exclusion principle
        
        % handle degenerate rectangles (XL0 == XU0) by setting both to -Inf
        equalLimits = (XL0 == XU0);
        XU0(equalLimits) = -Inf;
        XL0(equalLimits) = -Inf;
        
        p = zeros(n, 1);
        for i = 0:d
            k = nchoosek(1:d, i);
            for j = 1:size(k, 1)
                X = XU0;
                X(:, k(j,:)) = XL0(:, k(j,:));
                p = p + (-1)^i * tvncdf(X, rho_vec, tol/8);
            end
        end
    end
    if nargout > 1
        err = repmat(tol, n, 1); % return tolerance as error estimate
    end

else % general multivariate case (d >= 4): use adaptive quasi-Monte Carlo integration
    % call mvtcdfqmc with nu=Inf for multivariate normal
    % mvtcdfqmc(a, b, Rho, nu, tol) computes P(a <= X <= b) for MVT/MVN
    p = zeros(n, 1);
    err = zeros(n, 1);
    for i = 1:n
        [p(i), err(i)] = mvtcdfqmc(XL0(i,:), XU0(i,:), Rho, Inf, tol_qmc);
    end
end

% clamp probability to [0, 1] to handle numerical errors
p(p < 0) = 0;
p(p > 1) = 1;


return % --*-- Unit tests --*--

%@test:1
% test univariate case: should match normcdf exactly
try
    xu = [0; 1; -1; 2];
    mu = 0.5;
    Sigma = 4;
    p = mvncdf(xu, mu, Sigma);
    p_expected = normcdf(xu, mu, sqrt(Sigma));
    t(1) = isequal(p, p_expected);
catch
    t = false;
end
T = all(t);
%@eof:1

%@test:2
% test standard normal at origin: P(X <= 0) = 0.5^D for independent variables
try
    % D = 2 (uses bvncdf)
    p2 = mvncdf([0, 0]);
    t(1) = isequal(p2, 0.25);
    % D = 3 (uses tvncdf)
    p3 = mvncdf([0, 0, 0]);
    t(2) = isequal(p3, 0.125);
    % D = 4 (uses mvtcdfqmc)
    p4 = mvncdf([0, 0, 0, 0]);
    t(3) = isequal(p4, 0.0625); % QMC has lower precision
catch
    t = false(3, 1);
end
T = all(t);
%@eof:2

%@test:3
% test diagonal covariance: product of independent univariate CDFs
try
    xu = [1, 2, 3];
    mu = [0, 0, 0];
    Sigma = diag([1, 4, 9]);
    p = mvncdf(xu, mu, Sigma);
    p_expected = normcdf(1, 0, 1) * normcdf(2, 0, 2) * normcdf(3, 0, 3);
    t(1) = isequal(p, p_expected);
catch
    t = false;
end
T = all(t);
%@eof:3

%@test:4
% test diagonal SIGMA as 1-by-D vector
try
    xu = [1, 2, 3];
    mu = [0, 0, 0];
    Sigma_vec = [1, 4, 9]; % variances as vector
    Sigma_mat = diag([1, 4, 9]); % variances as matrix
    p1 = mvncdf(xu, mu, Sigma_vec);
    p2 = mvncdf(xu, mu, Sigma_mat);
    t(1) = isequal(p1, p2);
catch
    t = false;
end
T = all(t);
%@eof:4

%@test:5
% test bivariate at origin with correlation: P(X1<=0, X2<=0) = 1/4 + arcsin(rho)/(2*pi)
try
    rho = 0.5;
    xu = [0, 0];
    Sigma = [1 rho; rho 1];
    p = mvncdf(xu, [], Sigma);
    p_expected = 0.25 + asin(rho) / (2*pi);
    t(1) = isequal(p, p_expected);
catch
    t = false;
end
T = all(t);
%@eof:5

%@test:6
% test multiple observations (bivariate, Octave reference values - first 10 points)
try
    mu = [1, -1];
    Sigma = [0.9, 0.4; 0.4, 0.3];
    [X1, X2] = meshgrid(linspace(-1, 3, 25)', linspace(-3, 1, 25)');
    XU = [X1(:), X2(:)];
    p = mvncdf(XU, mu, Sigma);
    p_expected = [0.00011878988774500, 0.00034404112322371, ...
                  0.00087682502191813, 0.00195221905058185, ...
                  0.00378235566873474, 0.00638175749734415, ...
                  0.00943764224329656, 0.01239164888125426, ...
                  0.01472750274376648, 0.01623228313374828]';
    t(1) = max(abs(p(1:10) - p_expected)) < 1e-16;
catch
    t = false;
end
T = all(t);
%@eof:6

%@test:7
% test multiple observations (bivariate, Octave reference values - last 10 points)
try
    mu = [1, -1];
    Sigma = [0.9, 0.4; 0.4, 0.3];
    [X1, X2] = meshgrid(linspace(-1, 3, 25)', linspace(-3, 1, 25)');
    XU = [X1(:), X2(:)];
    p = mvncdf(XU, mu, Sigma);
    p_expected = [0.8180695783608276, 0.8854485749482751, ...
                  0.9308108777385832, 0.9579855743025508, ...
                  0.9722897881414742, 0.9788150170059926, ...
                  0.9813597788804785, 0.9821977956568989, ...
                  0.9824283794464095, 0.9824809345614861]';
    t(1) = max(abs(p(616:625) - p_expected)) < 1e-16;
catch
    t = false;
end
T = all(t);
%@eof:7

%@test:8
% test rectangle integration (bivariate, Octave reference)
try
    mu = [0, 0];
    Sigma = [0.25, 0.3; 0.3, 1];
    [p, err] = mvncdf([0, 0], [1, 1], mu, Sigma);
    t(1) = abs(p - 0.2097424404755626) < 1e-16;
    t(2) = isnan(err); % bivariate returns NaN for error
catch
    t = false(2, 1);
end
T = all(t);
%@eof:8

%@test:9
% test bivariate (Octave reference)
try
    x = [1 2];
    mu = [0.5 1.5];
    Sigma = [1.0, 0.5; 0.5, 1.0];
    p = mvncdf(x, mu, Sigma);
    t(1) = abs(p - 0.546244443857090) < 1e-15;
catch
    t = false;
end
T = all(t);
%@eof:9

%@test:10
% test rectangle integration bivariate (Octave reference)
try
    xu = [1 2];
    mu = [0.5 1.5];
    Sigma = [1.0, 0.5; 0.5, 1.0];
    xl = [-inf 0];
    p = mvncdf(xl, xu, mu, Sigma);
    t(1) = abs(p - 0.482672935215631) < 1e-15;
catch
    t = false;
end
T = all(t);
%@eof:10

%@test:11
% test trivariate case
try
    xu = [0.5, 0.3, 0.2];
    mu = [0, 0, 0];
    Sigma = [1.0 0.3 0.2; 0.3 1.0 0.1; 0.2 0.1 1.0];
    [p, err] = mvncdf(xu, mu, Sigma);
    % result should be valid probability
    t(1) = p >= 0 && p <= 1;
    % error should be the tolerance
    t(2) = isequal(err, 1e-8);
catch
    t = false(2, 1);
end
T = all(t);
%@eof:11

%@test:12
% test 4-dimensional case with QMC integration (Octave reference)
try
    fD = (-2:2)';
    XU = repmat(fD, 1, 4);
    p = mvncdf(XU);
    p_expected = [0; 0.0006; 0.0625; 0.5011; 0.9121];
    t(1) = max(abs(p - p_expected)) < 1e-3; % QMC has ~1e-4 precision
catch
    t = false;
end
T = all(t);
%@eof:12

%@test:13
% test empty input handling
try
    % 0 observations with 1 dimension: returns 1
    p = mvncdf(zeros(0,1), [], []);
    t(1) = isequal(p, 1);
    % also works with double.empty
    p2 = mvncdf(double.empty(0,1), double.empty(0,1), []);
    t(2) = isequal(p2, 1);
catch
    t = false(2, 1);
end
T = all(t);
%@eof:13

%@test:14
% test that d=0 (no columns) throws error
try
    mvncdf(zeros(1,0), [], []);
    t(1) = false;  % should not reach here
catch
    t(1) = true;   % error expected
end
try
    mvncdf([], [], []);
    t(2) = false;  % should not reach here
catch
    t(2) = true;   % error expected
end
T = all(t);
%@eof:14


end % mvncdf