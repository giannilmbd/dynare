function inv = betainv (x, a, b)
% BETAINV  Quantile function of the Beta distribution
%  INV = betainv(X, A, B) computes, for each element of X, the
%  quantile (the inverse of the CDF) at X of the Beta distribution
%  with parameters A and B (i.e. mean of the distribution is
%  A/(A+B) and variance is A*B/(A+B)^2/(A+B+1) ).

% Adapted for MATLAB (R) from GNU Octave 3.0.1
% Original file: statistics/distributions/betainv.m
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>

% Copyright © 2012 Rik Wehbring
% Copyright © 1995-2016 Kurt Hornik
% Copyright © 2022-2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
% Copyright © 2008-2026 Dynare Team
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

if (nargin ~= 3)
    error ('betainv: you must give three arguments');
end

if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
        error ('betainv: x, a and b must be of common size or scalars');
    end
end

if (~isreal (x) || ~isreal (a) || ~isreal (b))
    error ('betainv: x, a and b must not be complex');
end

if (isa (x, 'single') || isa (a, 'single') || isa (b, 'single'))
    inv = zeros (size (x), 'single');
    myeps = eps ('single');
else
    inv = zeros (size (x));
    myeps = eps;
end

k = find ((x < 0) | (x > 1) | ~(a > 0) | ~(b > 0) | isnan (x));
if (any (k))
    inv (k) = NaN;
end

k = find ((x == 1) & (a > 0) & (b > 0));
if (any (k))
    inv (k) = 1;
end

k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
if (any (k))
    if (~isscalar(a) || ~isscalar(b))
        a = a (k);
        b = b (k);
        y = a ./ (a + b);
    else
        y = a / (a + b) * ones (size (k));
    end
    x = x (k);
    l = find (y < myeps);
    if (any (l))
        y(l) = sqrt (myeps) * ones (length (l), 1);
    end
    l = find (y > 1 - myeps);
    if (any (l))
        y(l) = 1 - sqrt (myeps) * ones (length (l), 1);
    end

    y_new = y;
    loopcnt = 0;
    h = inf (size (y_new));
    while (max (abs (h)) >= sqrt (myeps) && loopcnt < 40)
        y_old = y_new;
        h     = (betacdf (y_old, a, b) - x) ./ betapdf (y_old, a, b);
        y_new = y_old - h;
        ind   = find (y_new <= myeps);
        if (any (ind))
            y_new (ind) = y_old (ind) / 10;
        end
        ind = find (y_new >= 1 - myeps);
        if (any (ind))
            y_new (ind) = 1 - (1 - y_old (ind)) / 10;
        end
        h = y_old - y_new;
        loopcnt = loopcnt + 1;
    end
    if loopcnt == 40
        warning ('betainv: calculation failed to converge for some values');
    end

    inv (k) = y_new;
end

return % --*-- Unit tests --*--

%@test:1
% Basic output: scalar parameters, vector x
p = [-1 0 0.75 1 2];
try
    result = betainv (p, ones (1, 5), 2*ones (1, 5));
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isequaln (isnan (result), [true false false false true]);
    t(3) = (result(2) == 0);
    t(4) = (abs (result(3) - 0.5) < eps);
    t(5) = (result(4) == 1);
end
T = all(t);
%@eof:1

%@test:2
% Scalar a, vector x and b
p = [-1 0 0.75 1 2];
try
    result = betainv (p, 1, 2*ones (1, 5));
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isequaln (isnan (result), [true false false false true]);
    t(3) = (result(2) == 0);
    t(4) = (abs (result(3) - 0.5) < eps);
    t(5) = (result(4) == 1);
end
T = all(t);
%@eof:2

%@test:3
% Scalar b, vector x and a
p = [-1 0 0.75 1 2];
try
    result = betainv (p, ones (1, 5), 2);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isequaln (isnan (result), [true false false false true]);
    t(3) = (result(2) == 0);
    t(4) = (abs (result(3) - 0.5) < eps);
    t(5) = (result(4) == 1);
end
T = all(t);
%@eof:3

%@test:4
% Invalid a values (zero, NaN) produce NaN; valid element returns correct value
p = [-1 0 0.75 1 2];
try
    result = betainv (p, [1 0 NaN 1 1], 2);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isequaln (result, [NaN NaN NaN 1 NaN]);
end
T = all(t);
%@eof:4

%@test:5
% Invalid b values (zero, NaN) produce NaN; valid element returns correct value
p = [-1 0 0.75 1 2];
try
    result = betainv (p, 1, 2*[1 0 NaN 1 1]);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isequaln (result, [NaN NaN NaN 1 NaN]);
end
T = all(t);
%@eof:5

%@test:6
% NaN in x propagates to output
p = [-1 0 0.75 1 2];
try
    result = betainv ([p(1:2) NaN p(4:5)], 1, 2);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isequaln (result, [NaN 0 NaN 1 NaN]);
end
T = all(t);
%@eof:6

%@test:7
% Single precision input x yields single output
p = [-1 0 0.75 1 2];
try
    result = betainv (single ([p, NaN]), 1, 2);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isa (result, 'single');
    t(3) = isequaln (isnan (result), [true false false false true true]);
    t(4) = (result(2) == 0);
    t(5) = (abs (result(3) - single (0.5)) <= eps ('single'));
    t(6) = (result(4) == 1);
end
T = all(t);
%@eof:7

%@test:8
% Single precision a yields single output
p = [-1 0 0.75 1 2];
try
    result = betainv ([p, NaN], single (1), 2);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isa (result, 'single');
    t(3) = isequaln (isnan (result), [true false false false true true]);
end
T = all(t);
%@eof:8

%@test:9
% Single precision b yields single output
p = [-1 0 0.75 1 2];
try
    result = betainv ([p, NaN], 1, single (2));
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    t(2) = isa (result, 'single');
    t(3) = isequaln (isnan (result), [true false false false true true]);
end
T = all(t);
%@eof:9

%@test:10
% Error: too few arguments
try
    betainv (1, 2);
    t(1) = false;
catch
    t(1) = true;
end
T = all(t);
%@eof:10

%@test:11
% Error: non-conformant array sizes
try
    betainv (ones (3), ones (2), ones (2));
    t(1) = false;
catch
    t(1) = true;
end
T = all(t);
%@eof:11

%@test:12
% Error: complex input
try
    betainv (1i, 2, 2);
    t(1) = false;
catch
    t(1) = true;
end
T = all(t);
%@eof:12

%@test:13
% Consistency check: betainv is the inverse of betacdf
q = 0.1:0.1:0.9;
a = 2; b = 5;
try
    x = betainv (q, a, b);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    for i = 1:length(q)
        t(i+1) = (abs (betacdf (x(i), a, b) - q(i)) < sqrt (eps));
    end
end
T = all(t);
%@eof:13
