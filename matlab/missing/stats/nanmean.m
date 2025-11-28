function y = nanmean(x, dim)
% y = nanmean(x)
% y = nanmean(x, dim)
% -------------------------------------------------------------------------
% Mean value, ignoring NaNs.
% -------------------------------------------------------------------------
% M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing values.
% - For vector input, M is the mean value of the non-NaN elements in X.
% - For matrix input, M is a row vector containing the mean value of non-NaN
%   elements in each column.
% - For N-D arrays, NANMEAN operates along the first non-singleton dimension.
%
% NANMEAN(X, DIM) takes the mean along dimension DIM of X.
% -------------------------------------------------------------------------
% Inputs:
%   x     [numeric array]   Input data (vector, matrix, or N-D array)
%   dim   [integer]         Dimension along which to compute the mean
% -------------------------------------------------------------------------
% Outputs:
%   y     [numeric array]   Mean values with NaNs treated as missing
% -------------------------------------------------------------------------
% Examples:
%   x = [1 2 NaN; 4 NaN 6];
%   nanmean(x)      % returns [2.5, 2, 6]
%   nanmean(x, 2)   % returns [1.5; 5]

% Copyright © 2011-2025 Dynare Team
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

if nargin < 2 || isempty(dim)
    % find the first non-singleton dimension
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end

notnan = ~isnan(x);   % create a mask for non-NaN values
x(isnan(x)) = 0;      % replace NaN with 0 for summing
s = sum(x, dim);      % compute sum of non-NaN values
n = sum(notnan, dim); % count non-NaN values along the dimension
y = s ./ n;           % compute mean (0/0 = NaN is handled automatically)


return % --*-- Unit tests --*--

%@test:1
% test basic vector with NaN
try
    x = [1, 2, NaN, 4];
    y = nanmean(x);
    t(1) = isequal(y, (1+2+4)/3);
catch
    t = false;
end
T = all(t);
%@eof:1

%@test:2
% test matrix along dim=1 (columns)
try
    x = [1 2 NaN; 4 NaN 6];
    y = nanmean(x, 1);
    expected = [2.5, 2, 6];
    t(1) = isequal(y, expected);
catch
    t = false;
end
T = all(t);
%@eof:2

%@test:3
% test matrix along dim=2 (rows)
try
    x = [1 2 NaN; 4 NaN 6];
    y = nanmean(x, 2);
    expected = [1.5; 5];
    t(1) = isequal(y, expected);
catch
    t = false;
end
T = all(t);
%@eof:3

%@test:4
% test all NaN returns NaN
try
    x = [NaN, NaN, NaN];
    t(1) = isnan(nanmean(x));
catch
    t = false;
end
T = all(t);
%@eof:4

%@test:5
% test no NaN matches regular mean
try
    x = [1, 2, 3, 4, 5];
    t(1) = isequal(nanmean(x), mean(x));
    X = [1 2 3; 4 5 6];
    t(2) = isequal(nanmean(X, 1), mean(X, 1));
    t(3) = isequal(nanmean(X, 2), mean(X, 2));
catch
    t = false(3, 1);
end
T = all(t);
%@eof:5

%@test:6
% test 3-D array along dim=1
try
    x = randn(3, 4, 5);
    x(1, 2, 3) = NaN;
    x(2, 1, 1) = NaN;
    y = nanmean(x, 1);
    t(1) = isequal(size(y), [1, 4, 5]);
    % check specific value where NaN was removed
    t(2) = isequal(y(1, 2, 3), mean([x(2,2,3), x(3,2,3)]));
catch
    t = false(2, 1);
end
T = all(t);
%@eof:6

%@test:7
% test 3-D array along dim=2
try
    x = ones(2, 3, 4);
    x(1, 2, :) = NaN;
    y = nanmean(x, 2);
    t(1) = isequal(size(y), [2, 1, 4]);
    t(2) = all(y(1, 1, :) == 1);  % mean of [1, NaN, 1] = 1
    t(3) = all(y(2, 1, :) == 1);  % mean of [1, 1, 1] = 1
catch
    t = false(3, 1);
end
T = all(t);
%@eof:7

%@test:8
% test 3-D array along dim=3
try
    x = randn(2, 3, 5);
    x(1, 1, 3) = NaN;
    y = nanmean(x, 3);
    t(1) = isequal(size(y), [2, 3]);
    % verify the element with NaN is computed from 4 values instead of 5
    x_subset = squeeze(x(1, 1, :));
    t(2) = isequal(y(1, 1), mean(x_subset(~isnan(x_subset))));
catch
    t = false(2, 1);
end
T = all(t);
%@eof:8

%@test:9
% test default dimension (first non-singleton)
try
    x = [1 NaN 3; 4 5 NaN];
    y1 = nanmean(x);      % should default to dim=1
    y2 = nanmean(x, 1);
    t(1) = isequal(y1, y2);
catch
    t = false;
end
T = all(t);
%@eof:9

%@test:10
% test column vector defaults to dim=1
try
    x = [1; NaN; 3; 4];
    y = nanmean(x);
    t(1) = isequal(y, (1 + 3 + 4) / 3);
catch
    t = false;
end
T = all(t);
%@eof:10

%@test:11
% test row vector defaults to dim=2
try
    x = [1, NaN, 3, 4];
    y = nanmean(x);
    t(1) = isequal(y, (1 + 3 + 4) / 3);
catch
    t = false;
end
T = all(t);
%@eof:11

%@test:12
% test scalar input
try
    t(1) = nanmean(5) == 5;
    t(2) = isnan(nanmean(NaN));
catch
    t = false(2, 1);
end
T = all(t);
%@eof:12
