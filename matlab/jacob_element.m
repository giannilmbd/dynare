function d=jacob_element(func,element,args)
% function d=jacob_element(func,element,args)
% returns an entry of the finite differences approximation to the Jacobian of func
%
% INPUTS
%    func       [function name]    string with name of the function
%    element    [int]              the index showing the element within the Jacobian that should be returned
%    args       [cell array]       arguments provided to func
%
% OUTPUTS
%    d          [double]           Jacobian [element]
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

assert(element <= length(args));

func = str2func(func);

h=1e-6;
margs=args;

args{element} = args{element} + h;
margs{element} = margs{element} - h;

d=(func(args{:})-func(margs{:}))/(2*h);

return % --*-- Unit tests --*--

%@test:1
% Test polynomial function: f(x,y) = x^3 + 2*x^2*y + 3*y^2 + 4*x + 5*y + 6
% Analytical derivatives:
%   df/dx = 3*x^2 + 4*x*y + 4
%   df/dy = 2*x^2 + 6*y + 5

% Test at point (x,y) = (2, 3)
x0 = 2;
y0 = 3;

% Analytical Jacobian at (2,3):
% df/dx = 3*4 + 4*2*3 + 4 = 12 + 24 + 4 = 40
% df/dy = 2*4 + 6*3 + 5 = 8 + 18 + 5 = 31

try
    jac_x = jacob_element('test_poly_2vars', 1, {x0, y0});
    jac_y = jacob_element('test_poly_2vars', 2, {x0, y0});
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = abs(jac_x - 40) < 1e-6;
    t(3) = abs(jac_y - 31) < 1e-6;
end
T = all(t);
%@eof:1

%@test:2
% Test at origin (0,0) for simpler verification
x0 = 0;
y0 = 0;

% Analytical Jacobian at (0,0):
% df/dx = 0 + 0 + 4 = 4
% df/dy = 0 + 0 + 5 = 5

try
    jac_x = jacob_element('test_poly_2vars', 1, {x0, y0});
    jac_y = jacob_element('test_poly_2vars', 2, {x0, y0});
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = abs(jac_x - 4) < 1e-6;
    t(3) = abs(jac_y - 5) < 1e-6;
end
T = all(t);
%@eof:2

%@test:3
% Test with negative values
x0 = -1;
y0 = -2;

% Analytical Jacobian at (-1,-2):
% df/dx = 3*1 + 4*(-1)*(-2) + 4 = 3 + 8 + 4 = 15
% df/dy = 2*1 + 6*(-2) + 5 = 2 - 12 + 5 = -5

try
    jac_x = jacob_element('test_poly_2vars', 1, {x0, y0});
    jac_y = jacob_element('test_poly_2vars', 2, {x0, y0});
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = abs(jac_x - 15) < 1e-6;
    t(3) = abs(jac_y - (-5)) < 1e-6;
end
T = all(t);
%@eof:3

%@test:4
% Test with large values to check numerical stability
x0 = 10;
y0 = 20;

% Analytical Jacobian at (10,20):
% df/dx = 3*100 + 4*10*20 + 4 = 300 + 800 + 4 = 1104
% df/dy = 2*100 + 6*20 + 5 = 200 + 120 + 5 = 325

try
    jac_x = jacob_element('test_poly_2vars', 1, {x0, y0});
    jac_y = jacob_element('test_poly_2vars', 2, {x0, y0});
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = abs(jac_x - 1104) < 1e-4;
    t(3) = abs(jac_y - 325) < 1e-6;
end
T = all(t);
%@eof:4

function f = test_poly_2vars(x, y)
    % Test polynomial: f(x,y) = x^3 + 2*x^2*y + 3*y^2 + 4*x + 5*y + 6
    f = x^3 + 2*x^2*y + 3*y^2 + 4*x + 5*y + 6;
