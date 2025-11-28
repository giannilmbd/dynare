function k = kurtosis(x, flag, dim)
% k = kurtosis(x)
% k = kurtosis(x, flag)
% k = kurtosis(x, flag, dim)
% k = kurtosis(x, flag, 'all')
% k = kurtosis(x, flag, vecdim)
% -------------------------------------------------------------------------
% Compute the sample kurtosis of the elements of x.
% -------------------------------------------------------------------------
% K = KURTOSIS(X) returns the sample kurtosis of the values in X.
% - For a vector input, K is the fourth central moment of X, divided by the
%   fourth power of its standard deviation.
% - For a matrix input, K is a row vector containing the sample kurtosis of
%   each column of X.
% - For N-D arrays, KURTOSIS operates along the first non-singleton dimension.
%
% The sample kurtosis is defined as:
%
%                  mean((x - mean(x)).^4)
%   kurtosis(x) = ------------------------
%                       std(x).^4
%
% where std(x) uses N in the denominator (population standard deviation).
%
% Note: This returns the kurtosis (not excess kurtosis). A normal
% distribution has kurtosis = 3. Some references define "excess kurtosis"
% as kurtosis - 3, so a normal distribution has excess kurtosis = 0.
%
% KURTOSIS(X, 0) adjusts the kurtosis for bias.
% KURTOSIS(X, 1) is the same as KURTOSIS(X), and does not adjust for bias.
%
% The adjusted (bias-corrected) kurtosis coefficient is:
%
%                              N - 1
%   kurtosis(x, 0) = 3 + -------------- * ((N + 1) * kurtosis(x, 1) - 3 * (N - 1))
%                        (N - 2)(N - 3)
%
% The bias-corrected kurtosis is obtained by replacing the
% sample second and fourth central moments by their unbiased versions.
% It is an unbiased estimate of the population kurtosis for normal populations.
%
% KURTOSIS(X, FLAG, 'all') is the kurtosis of all the elements of X.
%
% KURTOSIS(X, FLAG, DIM) takes the kurtosis along dimension DIM of X.
%
% KURTOSIS(X, FLAG, VECDIM) finds the kurtosis of the elements of X based
% on the dimensions specified in the vector VECDIM.
%
% KURTOSIS treats NaNs as missing values, and removes them.
% -------------------------------------------------------------------------
% Inputs:
%   x      [numeric array]         Input data (vector, matrix, or N-D array)
%   flag   [scalar]                0 = bias-corrected, 1 = not corrected (default)
%   dim    [scalar/vector/'all']   Dimension(s) to operate along
% -------------------------------------------------------------------------
% Outputs:
%   k      [numeric array]         Sample kurtosis values
% -------------------------------------------------------------------------
% Examples:
%   % Kurtosis of a normal distribution (should be ~3)
%   x = randn(10000, 1);
%   k = kurtosis(x)  % returns approximately 3
%
%   % Uniform distribution has lower kurtosis (platykurtic)
%   x = rand(10000, 1);
%   k = kurtosis(x)  % returns approximately 1.8
%
%   % Laplace/double-exponential has higher kurtosis (leptokurtic)
%   x = log(rand(10000, 1) ./ rand(10000, 1));
%   k = kurtosis(x)  % returns approximately 6
%
%   % Bias-corrected kurtosis
%   k = kurtosis(x, 0)
%
%   % Kurtosis along specific dimension
%   X = randn(10, 5);
%   k = kurtosis(X, 1, 2)  % kurtosis along rows
% -------------------------------------------------------------------------
% References:
% [1] Joanes, D.N. and C.A. Gill (1998), "Comparing measures of sample
%     skewness and kurtosis", The Statistician, 47(1), pages 183-189.
%     doi: 10.1111/1467-9884.00122
%
% Algorithm adapted from GNU Octave (version 9.x):
% https://github.com/gnu-octave/octave/blob/default/scripts/statistics/kurtosis.m

% Copyright © 1996-2025 The Octave Project Developers
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
if ~(isnumeric(x) || islogical(x))
    error('kurtosis: X must be a numeric or logical array');
end

% handle empty input
if isempty(x)
    k = NaN;
    return;
end

% convert logical to double for calculations
if islogical(x)
    x = double(x);
end

% handle flag argument
if nargin < 2 || isempty(flag)
    flag = 1;  % default: do not use bias correction
else
    if ~(isscalar(flag) && (flag == 0 || flag == 1))
        error('kurtosis: FLAG must be 0 or 1');
    end
end

% handle dimension argument
if nargin < 3 || isempty(dim)
    % find the first non-singleton dimension
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
    use_all = false;
    use_vecdim = false;
elseif ischar(dim) && strcmpi(dim, 'all')
    use_all = true;
    use_vecdim = false;
elseif isnumeric(dim) && isvector(dim) && all(dim == fix(dim)) && all(dim > 0)
    if isscalar(dim)
        use_all = false;
        use_vecdim = false;
    else
        use_all = false;
        use_vecdim = true;
    end
else
    error('kurtosis: DIM must be a positive integer scalar, vector, or ''all''');
end

% compute kurtosis
if use_all
    % operate on all elements
    x = x(:);
    n = sum(~isnan(x));
    x0 = x - nanmean(x);
    s2 = nanmean(x0.^2);  % biased variance (1/N)
    m4 = nanmean(x0.^4);
    if s2 == 0
        k = NaN(class(x));
    else
        k = m4 / (s2^2);
    end
elseif use_vecdim
    % operate along multiple dimensions
    % permute and reshape to collapse vecdim into one dimension
    sz = size(x);
    dim = sort(dim);
    % check valid dimensions
    if any(dim > length(sz))
        dim(dim > length(sz)) = [];  % remove dimensions beyond array
    end
    if isempty(dim)
        k = x;
        return;
    end

    % calculate kurtosis over specified dimensions
    n = sum(~isnan(x), dim);
    x0 = x - nanmean(x, dim);
    s2 = nanmean(x0.^2, dim);  % biased variance
    m4 = nanmean(x0.^4, dim);
    k = m4 ./ (s2.^2);
    k(s2 == 0) = NaN;
else
    % operate along single dimension
    n = size(x, dim);
    x0 = x - nanmean(x, dim); % center the data (subtract mean), handling NaN
    s2 = nanmean(x0.^2, dim); % compute biased variance (divide by N, not N-1)
    m4 = nanmean(x0.^4, dim); % compute fourth central moment
    k = m4 ./ (s2.^2);        % compute kurtosis
    k(s2 == 0) = NaN;         % handle zero variance case
    % for bias correction, need to count non-NaN values
    if flag == 0
        n = sum(~isnan(x), dim);
    end
end

% apply bias correction
if flag == 0
    % bias-corrected kurtosis: k0 = 3 + (N-1)/((N-2)*(N-3)) * ((N+1)*k1 - 3*(N-1)) requires n > 3
    % note: n can be an array (not scalar) due to NaN handling, so we use
    % element-wise k(n < 4) = NaN instead of a simple if-statement
    C = (n - 1) ./ ((n - 2) .* (n - 3));
    k = 3 + C .* ((n + 1) .* k - 3 .* (n - 1));
    k(n < 4) = NaN;
end


return % --*-- Unit tests --*--

%@test:1
% test basic kurtosis calculation
try
    x = [-1; 0; 0; 0; 1];
    y = [x, 2*x];
    k = kurtosis(y);
    expected = [2.5, 2.5];
    t(1) = max(abs(k - expected)) < sqrt(eps);
catch
    t = false;
end
T = all(t);
%@eof:1

%@test:2
% test antisymmetry property: kurtosis(-x) = kurtosis(x)
try
    x = [-3, 0, 1];
    t(1) = abs(kurtosis(x) - kurtosis(-x)) < eps;
catch
    t = false;
end
T = all(t);
%@eof:2

%@test:3
% test constant input: kurtosis should be NaN (zero variance)
try
    t(1) = all(isnan(kurtosis(ones(3, 5))));
    t(2) = isnan(kurtosis(5));
    t(3) = isnan(kurtosis(1, [], 3));
catch
    t = false(3, 1);
end
T = all(t);
%@eof:3

%@test:4
% test dimension argument
try
    x = [1:5, 10; 1:5, 10];
    k0 = kurtosis(x, 0, 2);  % bias-corrected, along rows
    k1 = kurtosis(x, 1, 2);  % not corrected, along rows
    expected0 = 5.4377317925288901 * [1; 1];
    expected1 = 2.9786509002956195 * [1; 1];
    t(1) = max(abs(k0 - expected0)) < 1e-13;
    t(2) = max(abs(k1 - expected1)) < 1e-13;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:4

%@test:5
% test empty flag defaults to 1
try
    x = [1:5, 10; 1:5, 10];
    k1 = kurtosis(x, [], 2);  % empty flag, should be same as flag=1
    k2 = kurtosis(x, 1, 2);
    t(1) = isequal(k1, k2);
catch
    t = false;
end
T = all(t);
%@eof:5

%@test:6
% test bias correction requires n >= 4
try
    x = [1, 2];
    k = kurtosis(x, 0);
    t(1) = isnan(k);
    x = [1, 2, 3];
    k = kurtosis(x, 0);
    t(2) = isnan(k);  % n=3 also requires NaN
catch
    t = false(2, 1);
end
T = all(t);
%@eof:6

%@test:7
% test NaN handling
try
    x = [1, 2, NaN, 4, 5, 6];
    k = kurtosis(x);
    % should ignore NaN and compute kurtosis of [1, 2, 4, 5, 6]
    k_expected = kurtosis([1, 2, 4, 5, 6]);
    t(1) = abs(k - k_expected) < eps;
catch
    t = false;
end
T = all(t);
%@eof:7

%@test:8
% test 'all' dimension option
try
    x = reshape(1:12, [3, 4]);
    k = kurtosis(x, 1, 'all');
    k_expected = kurtosis(x(:));
    t(1) = abs(k - k_expected) < eps;
catch
    t = false;
end
T = all(t);
%@eof:8

%@test:9
% test analytical: normal distribution has kurtosis = 3
try
    rng(42);  % for reproducibility
    x = randn(100000, 1);
    k = kurtosis(x);
    t(1) = abs(k - 3) < 0.1;  % within 0.1 of theoretical value
catch
    t = false;
end
T = all(t);
%@eof:9

%@test:10
% test analytical: uniform distribution has kurtosis = 9/5 = 1.8
try
    rng(123);
    x = rand(100000, 1);
    k = kurtosis(x);
    t(1) = abs(k - 1.8) < 0.05;
catch
    t = false;
end
T = all(t);
%@eof:10

%@test:11
% test analytical: exponential distribution has kurtosis = 9
try
    rng(456);
    x = -log(rand(100000, 1));  % exponential(1)
    k = kurtosis(x);
    t(1) = abs(k - 9) < 0.2;
catch
    t = false;
end
T = all(t);
%@eof:11

%@test:12
% test analytical: Laplace/double-exponential has kurtosis = 6
try
    rng(789);
    % Laplace = difference of two exponentials
    x = -log(rand(100000, 1)) + log(rand(100000, 1));
    k = kurtosis(x);
    t(1) = abs(k - 6) < 0.2;
catch
    t = false;
end
T = all(t);
%@eof:12

%@test:13
% test logical input
try
    x = logical([0, 0, 0, 1, 1]);
    k = kurtosis(x);
    k_expected = kurtosis(double(x));
    t(1) = abs(k - k_expected) < eps;
catch
    t = false;
end
T = all(t);
%@eof:13

%@test:14
% test higher dimension (dim=3)
try
    x = randn(3, 4, 5);
    k = kurtosis(x, 1, 3);
    t(1) = isequal(size(k), [3, 4]);
catch
    t = false;
end
T = all(t);
%@eof:14

%@test:15
% test empty input
try
    k = kurtosis([]);
    t(1) = isnan(k);
catch
    t = false;
end
T = all(t);
%@eof:15

%@test:16
% test that scaling does not affect kurtosis
% kurtosis(a*x + b) = kurtosis(x) for any constants a != 0, b
try
    x = [1, 2, 3, 4, 5, 10];
    k1 = kurtosis(x);
    k2 = kurtosis(2*x + 5);
    k3 = kurtosis(-3*x + 100);
    t(1) = abs(k1 - k2) < 10*eps;
    t(2) = abs(k1 - k3) < 10*eps;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:16

%@test:17
% test relationship between biased and unbiased estimators
% verify the correction formula for a specific case
try
    x = [1, 2, 3, 4, 5];
    n = 5;
    k1 = kurtosis(x, 1);  % biased
    k0 = kurtosis(x, 0);  % unbiased
    % k0 = 3 + (n-1)/((n-2)*(n-3)) * ((n+1)*k1 - 3*(n-1))
    k0_expected = 3 + (n-1)/((n-2)*(n-3)) * ((n+1)*k1 - 3*(n-1));
    t(1) = abs(k0 - k0_expected) < eps;
catch
    t = false;
end
T = all(t);
%@eof:17


end % kurtosis

