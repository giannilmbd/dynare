function t = wblinv(proba, scale, shape)

% Inverse cumulative distribution function.
%
% INPUTS
% - proba [double] Probability, scalar or vector with values in [0,1].
% - scale [double] Positive hyperparameter (scalar).
% - shape [double] Positive hyperparameter (scalar).
%
% OUTPUTS
% - t     [double] Quantile(s), same size as proba.

% Copyright © 2015-2023 Dynare Team
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

% Check input arguments.

if nargin<3
    error('Three input arguments required!')
end

if ~isnumeric(proba) || ~isreal(proba) || any(proba(:)<0) || any(proba(:)>1)
    error('First input argument must be a real scalar or vector with values in [0,1] (probability)!')
end

if ~isnumeric(scale) || ~isscalar(scale) || ~isreal(scale) || scale<=0
    error('Second input argument must be a real positive scalar (scale parameter of the Weibull distribution)!')
end

if ~isnumeric(shape) || ~isscalar(shape) || ~isreal(shape) || shape<=0
    error('Third input argument must be a real positive scalar (shape parameter of the Weibull distribution)!')
end

% Vectorized closed-form inverse CDF.
t = zeros(size(proba));
t(proba > 1-2*eps()) = Inf;
k = proba >= 2*eps() & proba <= 1-2*eps();
if any(k(:))
    t(k) = exp(log(scale) + log(-log(1-proba(k)))./shape);
end

return % --*-- Unit tests --*--

%@test:1
try
   x = wblinv(0, 1, 2);
   t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = isequal(x, 0);
end
T = all(t);
%@eof:1

%@test:2
try
   x = wblinv(1, 1, 2);
   t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = isinf(x);
end
T = all(t);
%@eof:2

%@test:3
scales = [.5, 1, 5];
shapes = [.1, 1, 2];
x = NaN(9,1);

try
   k = 0;
   for i=1:3
      for j=1:3
          k = k+1;
          x(k) = wblinv(.5, scales(i), shapes(j));
      end
   end
   t(1) = true;
catch
   t(1) = false;
end

if t(1)
    k = 1;
    for i=1:3
        for j=1:3
            k = k+1;
            t(k) = abs(x(k-1)-scales(i)*log(2)^(1/shapes(j)))<1e-12;
        end
    end
end
T = all(t);
%@eof:3

%@test:4
debug = false;
scales = [ .5, 1, 5];
shapes = [ 1, 2, 3];
x = NaN(9,1);
p = 1e-1;

try
   k = 0;
   for i=1:3
      for j=1:3
          k = k+1;
          x(k) = wblinv(p, scales(i), shapes(j));
      end
   end
   t(1) = true;
catch
   t(1) = false;
end

if t(1)
   k = 1;
   for i=1:3
      for j=1:3
          k = k+1;
          shape = shapes(j);
          scale = scales(i);
          density = @(z) exp(lpdfgweibull(z,shape,scale));
          if debug
              [shape, scale, x(k-1)]
          end
          if isoctave
              s = quadv(density, 0, x(k-1),1e-10);
          else
              s = integral(density, 0, x(k-1));
          end
          if debug
              [s, abs(p-s)]
          end
          if isoctave
              t(k) = abs(p-s)<1e-9;
          else
              t(k) = abs(p-s)<1e-12;
          end
      end
   end
end
T = all(t);
%@eof:4

%@test:5
% Test with vector-valued proba input.

% Set the hyperparameters of a Weibull distribution.
scale = .5;
shape = 1.5;

% Build a probability vector including boundary cases.
q = [0, 0.25, 0.5, 0.75, 1];

try
    x = wblinv(q, scale, shape);
    t(1) = true;
catch
    t(1) = false;
end

% Check the results.
if t(1)
    t(2) = isequal(size(x), size(q));
    t(3) = isequal(x(1), 0);     % p = 0 → x = 0
    t(4) = isinf(x(end));         % p = 1 → x = Inf
    % Check round-trip consistency with wblcdf for interior points.
    for i = 2:length(q)-1
        t(i+3) = abs(wblcdf(x(i), scale, shape) - q(i)) < 1e-12;
    end
end
T = all(t);
%@eof:5
