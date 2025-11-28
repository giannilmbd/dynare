function s = skewness(x, flag, dim)
% s = skewness(x)
% s = skewness(x, flag)
% s = skewness(x, flag, dim)
% s = skewness(x, flag, 'all')
% s = skewness(x, flag, vecdim)
% -------------------------------------------------------------------------
% Compute the sample skewness of the elements of x.
% -------------------------------------------------------------------------
% S = SKEWNESS(X) returns the sample skewness of the values in X.
% - For a vector input, S is the third central moment of X, divided by the
%   cube of its standard deviation.
% - For a matrix input, S is a row vector containing the sample skewness of
%   each column of X.
% - For N-D arrays, SKEWNESS operates along the first non-singleton dimension.
%
% The sample skewness is defined as:
%
%                  mean((x - mean(x)).^3)
%   skewness(x) = ------------------------
%                       std(x).^3
%
% where std(x) uses N in the denominator (population standard deviation).
%
% SKEWNESS(X, 0) adjusts the skewness for bias.
% SKEWNESS(X, 1) is the same as SKEWNESS(X), and does not adjust for bias.
%
% The adjusted (bias-corrected) skewness coefficient is:
%
%                      sqrt(N*(N-1))
%   skewness(x, 0) = --------------- * skewness(x, 1)
%                         N - 2
%
% The bias-corrected skewness coefficient is obtained by replacing the
% sample second and third central moments by their unbiased versions.
% It is an unbiased estimate of the population skewness for normal populations.
%
% SKEWNESS(X, FLAG, 'all') is the skewness of all the elements of X.
%
% SKEWNESS(X, FLAG, DIM) takes the skewness along dimension DIM of X.
%
% SKEWNESS(X, FLAG, VECDIM) finds the skewness of the elements of X based
% on the dimensions specified in the vector VECDIM.
%
% SKEWNESS treats NaNs as missing values, and removes them.
% -------------------------------------------------------------------------
% Inputs:
%   x      [numeric array]         Input data (vector, matrix, or N-D array)
%   flag   [scalar]                0 = bias-corrected, 1 = not corrected (default)
%   dim    [scalar/vector/'all']   Dimension(s) to operate along
% -------------------------------------------------------------------------
% Outputs:
%   s      [numeric array]         Sample skewness values
% -------------------------------------------------------------------------
% Examples:
%   % Skewness of a symmetric distribution (should be ~0)
%   x = [-1, 0, 1];
%   s = skewness(x)  % returns 0
%
%   % Right-skewed distribution (positive skewness)
%   x = [1, 2, 3, 4, 10];
%   s = skewness(x)  % returns positive value (1.1384)
%
%   % Left-skewed distribution (negative skewness)
%   x = [1, 7, 8, 9, 10];
%   s = skewness(x)  % returns negative value (-1.1384)
%
%   % Bias-corrected skewness
%   x = [1, 7, 8, 9, 10];
%   s = skewness(x, 0) % returns -1.6971
%
%   % Skewness along specific dimension
%   X = randn(10, 5);
%   s = skewness(X, 1, 1)  % skewness along columns
%   s = skewness(X, 1, 2)  % skewness along rows
% -------------------------------------------------------------------------
% References:
% [1] Joanes, D.N. and C.A. Gill (1998), "Comparing measures of sample
%     skewness and kurtosis", The Statistician, 47(1), pages 183-189.
%     doi: 10.1111/1467-9884.00122
%
% Algorithm adapted from GNU Octave (version 9.x):
% https://github.com/gnu-octave/octave/blob/default/scripts/statistics/skewness.m

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
    error('skewness: X must be a numeric or logical array');
end

% handle empty input
if isempty(x)
    s = NaN;
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
        error('skewness: FLAG must be 0 or 1');
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
    error('skewness: DIM must be a positive integer scalar, vector, or ''all''');
end

% compute skewness
if use_all
    % operate on all elements
    x = x(:);
    n = sum(~isnan(x));
    x0 = x - nanmean(x);
    s2 = nanmean(x0.^2);  % biased variance (1/N)
    m3 = nanmean(x0.^3);
    if s2 == 0
        s = NaN(class(x));
    else
        s = m3 / (s2^1.5);
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
        s = x;
        return;
    end

    % calculate skewness over specified dimensions
    n = sum(~isnan(x), dim);
    x0 = x - nanmean(x, dim);
    s2 = nanmean(x0.^2, dim);  % biased variance
    m3 = nanmean(x0.^3, dim);
    s = m3 ./ (s2.^1.5);
    s(s2 == 0) = NaN;
else
    % operate along single dimension
    n = size(x, dim);
    x0 = x - nanmean(x, dim); % center the data (subtract mean), handling NaN
    s2 = nanmean(x0.^2, dim); % compute biased variance (divide by N, not N-1)
    m3 = nanmean(x0.^3, dim); % compute third central moment
    s = m3 ./ (s2.^1.5);      % compute skewness
    s(s2 == 0) = NaN;         % handle zero variance case
    % for bias correction, need to count non-NaN values
    if flag == 0
        n = sum(~isnan(x), dim);
    end
end

% apply bias correction
if flag == 0
    % bias-corrected skewness: s0 = sqrt(N*(N-1))/(N-2) * s1 requires n > 2
    % note: n can be an array (not scalar) due to NaN handling, so we use
    % element-wise s(n < 3) = NaN instead of a simple if-statement
    C = sqrt(n .* (n - 1)) ./ (n - 2);
    s = C .* s;
    s(n < 3) = NaN;
end


return % --*-- Unit tests --*--

%@test:1
% test symmetric distribution: skewness should be zero
try
    x = [-1, 0, 1];
    s = skewness(x);
    t(1) = isequal(s,0);
catch
    t = false;
end
T = all(t);
%@eof:1

%@test:2
% test that left and right skew have opposite signs
try
    s1 = skewness([-2, 0, 1]);
    s2 = skewness([-1, 0, 2]);
    t(1) = s1 < 0;  % left-skewed (tail on left)
    t(2) = s2 > 0;  % right-skewed (tail on right)
catch
    t = false(2, 1);
end
T = all(t);
%@eof:2

%@test:3
% test antisymmetry: skewness(-x) = -skewness(x)
try
    x = [-3, 0, 1];
    t(1) = isequal(skewness(x) + skewness(-x), 0);
    t(2) = isequal(skewness([-3, 0, 1]) + skewness([-1, 0, 3]), 0);
catch
    t = false(2, 1);
end
T = all(t);
%@eof:3

%@test:4
% test constant input: skewness should be NaN (zero variance)
try
    t(1) = all(isnan(skewness(ones(3, 5))));
    t(2) = isnan(skewness(5));
catch
    t = false(2, 1);
end
T = all(t);
%@eof:4

%@test:5
% test matrix input: operates along first non-singleton dimension
try
    x = [0; 0; 0; 1];
    y = [x, 2*x];
    s = skewness(y);
    % expected value: 1.154700538379251 for both columns
    expected = 1.154700538379251 * [1, 1];
    t(1) = max(abs(s - expected)) < 1e-15;
catch
    t = false;
end
T = all(t);
%@eof:5

%@test:6
% test dimension argument
try
    x = [1:5, 10; 1:5, 10];
    s1 = skewness(x, 0, 2);  % bias-corrected, along rows
    s2 = skewness(x, 1, 2);  % not corrected, along rows
    expected1 = 1.439590274527954 * [1; 1];
    expected2 = 1.051328089232020 * [1; 1];
    t(1) = max(abs(s1 - expected1)) < 1e-15;
    t(2) = max(abs(s2 - expected2)) < 1e-15;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:6

%@test:7
% test empty flag defaults to 1
try
    x = [1:5, 10; 1:5, 10];
    s1 = skewness(x, [], 2);  % empty flag, should be same as flag=1
    s2 = skewness(x, 1, 2);
    t(1) = isequal(s1, s2);
catch
    t = false;
end
T = all(t);
%@eof:7

%@test:8
% test bias correction requires n >= 3
try
    x = [1, 2];
    s = skewness(x, 0);
    t(1) = isnan(s);
catch
    t = false;
end
T = all(t);
%@eof:8

%@test:9
% test NaN handling
try
    x = [1, 2, NaN, 4, 5, 6];
    s = skewness(x);
    % should ignore NaN and compute skewness of [1, 2, 4, 5]
    s_expected = skewness([1, 2, 4, 5, 6]);
    t(1) = isequal(s, s_expected);
catch
    t = false;
end
T = all(t);
%@eof:9

%@test:10
% test 'all' dimension option
try
    x = randn(3,4);
    s = skewness(x, 1, 'all');
    s_expected = skewness(x(:));
    t(1) = isequal(s, s_expected);
catch
    t = false;
end
T = all(t);
%@eof:10

%@test:11
% test analytical: exponential distribution has skewness = 2
% For large n, sample skewness approaches theoretical value
try
    rng(42);  % for reproducibility
    x = -log(rand(100000, 1));  % exponential(1) distribution
    s = skewness(x);
    t(1) = abs(s - 2) < 0.05;  % within 5% of theoretical value
catch
    t = false;
end
T = all(t);
%@eof:11

%@test:12
% test analytical: chi-squared(k) has skewness = sqrt(8/k)
% For k=4: skewness = sqrt(2) ≈ 1.414
try
    rng(123);
    k = 4;  % degrees of freedom
    x = 2 * sum(-log(rand(100000, k/2)), 2);  % chi-squared(k) = 2 * Gamma(k/2, 1)
    s = skewness(x);
    expected = sqrt(8/k);  % = sqrt(2) ≈ 1.414
    t(1) = abs(s - expected) < 0.05;  % within 5% of theoretical value
catch
    t = false;
end
T = all(t);
%@eof:12

%@test:13
% test logical input
try
    x = logical([0, 0, 0, 1, 1]);
    s = skewness(x);
    s_expected = skewness(double(x));
    t(1) = abs(s - s_expected) < eps;
catch
    t = false;
end
T = all(t);
%@eof:13

%@test:14
% test higher dimension (dim=3)
try
    x = randn(3, 4, 5); x(1,1,1) = 10;
    s = skewness(x, 1, 3);
    t(1) = isequal(size(s), [3, 4]);
catch
    t = false;
end
T = all(t);
%@eof:14

%@test:15
% test empty input
try
    s = skewness([]);
    t(1) = isnan(s);
catch
    t = false;
end
T = all(t);
%@eof:15


end % skewness

