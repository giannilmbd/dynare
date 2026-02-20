function d=hess_element(func,element1,element2,args)
% function d=hess_element(func,element1,element2,args)
% returns an entry of the finite differences approximation to the Hessian of func
%
% INPUTS
%    func       [function name]    string with name of the function
%    element1   [int]              the indices showing the element within the Hessian that should be returned
%    element2   [int]
%    args       [cell array]       arguments provided to func
%
% OUTPUTS
%    d          [double]           the (element1,element2) entry of the Hessian
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2010-2025 Dynare Team
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

assert(element1 <= length(args) && element2 <= length(args));

func = str2func(func);

h=1e-6;
p10 = args;
p01 = args;
m10 = args;
m01 = args;
p11 = args;
m11 = args;

p10{element1} = p10{element1} + h;
m10{element1} = m10{element1} - h;

p11{element1} = p11{element1} + h;
m11{element1} = m11{element1} - h;

p01{element2} = p01{element2} + h;
m01{element2} = m01{element2} - h;

p11{element2} = p11{element2} + h;
m11{element2} = m11{element2} - h;

% From Abramowitz and Stegun. Handbook of Mathematical Functions (1965)
% formulas 25.3.24 and 25.3.27 p. 884
if element1==element2
    d = (16*func(p10{:})...
         +16*func(m10{:})...
         -30*func(args{:})...
         -func(p11{:})...
         -func(m11{:}))/(12*h^2);
else
    d = (func(p10{:})...
         +func(m10{:})...
         +func(p01{:})...
         +func(m01{:})...
         -2*func(args{:})...
         -func(p11{:})...
         -func(m11{:}))/(-2*h^2);
end

return % --*-- Unit tests --*--

%@test:1
% Test polynomial function: f(x,y) = x^3 + 2*x^2*y + 3*y^2 + 4*x + 5*y + 6
% Analytical derivatives:
%   df/dx = 3*x^2 + 4*x*y + 4
%   df/dy = 2*x^2 + 6*y + 5
%   d2f/dx2 = 6*x + 4*y
%   d2f/dy2 = 6
%   d2f/dxdy = d2f/dydx = 4*x

% Test at point (x,y) = (2, 3)
x0 = 2;
y0 = 3;

% Analytical Hessian at (2,3):
% H = [ 6*2 + 4*3,  4*2 ]  = [ 24,  8 ]
%     [ 4*2,        6   ]    [  8,  6 ]

try
    h11 = hess_element('test_poly_2vars', 1, 1, {x0, y0});
    h22 = hess_element('test_poly_2vars', 2, 2, {x0, y0});
    h12 = hess_element('test_poly_2vars', 1, 2, {x0, y0});
    h21 = hess_element('test_poly_2vars', 2, 1, {x0, y0});
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = abs(h11 - 24) < 5e-2;
    t(3) = abs(h22 - 6) < 5e-2;
    t(4) = abs(h12 - 8) < 5e-2;
    t(5) = abs(h21 - 8) < 5e-2;
    t(6) = abs(h12 - h21) < 1e-10; % Verify symmetry
end
T = all(t);
%@eof:1

%@test:2
% Test at origin (0,0) for simpler verification
x0 = 0;
y0 = 0;

% Analytical Hessian at (0,0):
% H = [ 0,  0 ]
%     [ 0,  6 ]
t=false(5,1);
try
    h11 = hess_element('test_poly_2vars', 1, 1, {x0, y0});
    h12 = hess_element('test_poly_2vars', 1, 2, {x0, y0});
    h21 = hess_element('test_poly_2vars', 2, 1, {x0, y0});
    h22 = hess_element('test_poly_2vars', 2, 2, {x0, y0});
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = abs(h11 - 0) < 1e-4;
    t(3) = abs(h12 - 0) < 1e-3;
    t(4) = abs(h21 - 0) < 1e-3;
    t(5) = abs(h22 - 6) < 1e-3;
end
T = all(t);
%@eof:2

%@test:3
% Test negative values
x0 = -1;
y0 = -2;

% Analytical Hessian at (-1,-2):
% H = [ 6*(-1) + 4*(-2),  4*(-1) ]  = [ -14, -4 ]
%     [ 4*(-1),           6      ]    [ -4,   6 ]
t=false(5,1);
try
    h11 = hess_element('test_poly_2vars', 1, 1, {x0, y0});
    h12 = hess_element('test_poly_2vars', 1, 2, {x0, y0});
    h21 = hess_element('test_poly_2vars', 2, 1, {x0, y0});
    h22 = hess_element('test_poly_2vars', 2, 2, {x0, y0});
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = abs(h11 - (-14)) < 5e-3;
    t(3) = abs(h12 - (-4)) < 1e-3;
    t(4) = abs(h21 - (-4)) < 1e-3;
    t(5) = abs(h22 - 6) < 1e-3;
end
T = all(t);
%@eof:3

function f = test_poly_2vars(x, y)
    % Test polynomial: f(x,y) = x^3 + 2*x^2*y + 3*y^2 + 4*x + 5*y + 6
    f = x^3 + 2*x^2*y + 3*y^2 + 4*x + 5*y + 6;
