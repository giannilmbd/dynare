function p = tvncdf(x, rho, tol)
% p = tvncdf(x, rho, tol)
% -------------------------------------------------------------------------
% Trivariate normal cumulative distribution function (CDF).
% -------------------------------------------------------------------------
% p = tvncdf(x, rho, tol) computes the trivariate normal cumulative
% distribution function for standardized variables (zero mean, unit variance)
% with correlation coefficients specified in rho.
%
% This implements equation (14) in Section 3.2 of Genz (2004), integrating
% each term separately using adaptive quadrature.
% -------------------------------------------------------------------------
% Inputs:
% x       [N-by-3 matrix]   Each row is a point, columns are variables (standardized)
% rho     [3-by-1 vector]   Correlation coefficients: [rho_21; rho_31; rho_32]
%                           where rho_ij is the correlation between variables i and j
% tol     [scalar]          Absolute error tolerance for quadrature (default 1e-8)
% -------------------------------------------------------------------------
% Outputs:
% p       [N-by-1 vector]   Probability values in [0, 1]
% -------------------------------------------------------------------------
% Example:
% x = [0.5, 0.3, -0.2];
% rho = [0.3; 0.1; 0.2];  % correlations: rho_21=0.3, rho_31=0.1, rho_32=0.2
% p = tvncdf(x, rho)
% -------------------------------------------------------------------------
% References:
% Drezner, Z. (1994) "Computation of the Trivariate Normal Integral",
%   Mathematics of Computation 62(205), pages 289-294. doi:10.2307/2153409
% Genz, A. (2004) "Numerical Computation of Rectangular Bivariate and
%   Trivariate Normal and t Probabilities", Statistics and Computing 14, pages 251-260.
%   doi: 10.1023/B:STCO.0000035304.20635.31
%
% Code adapted from the GNU Octave statistics package (version 1.7.7):
% https://github.com/gnu-octave/statistics/blob/main/inst/dist_fun/mvncdf.m
%
% See also mvncdf, bvncdf, normcdf.

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

% default tolerance
if nargin < 3 || isempty(tol)
    tol = 1e-8;
end

% input validation
if size(x, 2) ~= 3
    error('tvncdf: x must be an N-by-3 matrix with each variable as a column vector');
end
if numel(rho) ~= 3
    error('tvncdf: rho must be a 3-element vector [rho_21; rho_31; rho_32]');
end
rho = rho(:); % ensure column vector

% get size of data
n = size(x, 1);

% find a permutation that makes rho_32 == max(rho) to improve numerical stability
[~, imax] = max(abs(rho));
if imax == 1     % swap variables 1 and 3
    rho_21 = rho(3); rho_31 = rho(2); rho_32 = rho(1);
    x = x(:, [3 2 1]);
elseif imax == 2 % swap variables 1 and 2
    rho_21 = rho(1); rho_31 = rho(3); rho_32 = rho(2);
    x = x(:, [2 1 3]);
else             % imax == 3, no swap needed
    rho_21 = rho(1); rho_31 = rho(2); rho_32 = rho(3);
end

% first term: Phi(x1) * Phi_2(x2, x3; rho_32)
phi_x1 = 0.5 * erfc(-x(:,1) / sqrt(2));
p1 = phi_x1 .* bvncdf(x(:, 2:3), [], rho_32);

% second term: integral involving rho_21
if abs(rho_21) > 0
    loLimit = 0;
    hiLimit = asin(rho_21);
    p2 = zeros(n, 1);
    for i = 1:n
        b1 = x(i, 1); bj = x(i, 2); bk = x(i, 3);
        if isfinite(b1) && isfinite(bj) && ~isnan(bk)
            p2(i) = quadgk(@(theta) tvnIntegrand(theta, b1, bj, bk, rho_21, rho_31, rho_32), ...
                           loLimit, hiLimit, 'AbsTol', tol/3, 'RelTol', 0);
        end
    end
else
    p2 = zeros(n, 1);
end

% third term: integral involving rho_31
if abs(rho_31) > 0
    loLimit = 0;
    hiLimit = asin(rho_31);
    p3 = zeros(n, 1);
    for i = 1:n
        b1 = x(i, 1); bj = x(i, 3); bk = x(i, 2);
        if isfinite(b1) && isfinite(bj) && ~isnan(bk)
            p3(i) = quadgk(@(theta) tvnIntegrand(theta, b1, bj, bk, rho_31, rho_21, rho_32), ...
                           loLimit, hiLimit, 'AbsTol', tol/3, 'RelTol', 0);
        end
    end
else
    p3 = zeros(n, 1);
end

% combine terms and clamp to [0, 1] to ensure the result is a valid probability
p = p1 + (p2 + p3) ./ (2 * pi);
p = max(0, min(1, p));


return % --*-- Unit tests --*--

%@test:1
% test with independent variables (all correlations zero)
% P(X1 <= x1, X2 <= x2, X3 <= x3) = Phi(x1) * Phi(x2) * Phi(x3)
try
    x = [0.5, -0.3, 0.8];
    rho = [0; 0; 0];
    p = tvncdf(x, rho);
    p_expected = normcdf(0.5) * normcdf(-0.3) * normcdf(0.8);
    t(1) = isequal(p, p_expected);
catch
    t = false;
end
T = all(t);
%@eof:1

%@test:2
% test at origin with positive correlations
try
    x = [0, 0, 0];
    rho = [0.3; 0.2; 0.1];
    p = tvncdf(x, rho);
    % result should be around 0.125 (1/8) for uncorrelated, slightly higher for positive correlations
    t(1) = ( abs(p - 0.173241295324728) > 1e-16 ) ;
catch
    t = false;
end
T = all(t);
%@eof:2

%@test:3
% test symmetry: swapping variables with same correlation should give same result
try
    x1 = [0.5, 0.3, 0.2];
    x2 = [0.3, 0.5, 0.2];  % swap first two
    rho = [0.4; 0.2; 0.3]; % rho_21=0.4, rho_31=0.2, rho_32=0.3
    % after swap: rho_21=0.4 (same), rho_31=0.3, rho_32=0.2
    rho_swapped = [0.4; 0.3; 0.2];
    p1 = tvncdf(x1, rho);
    p2 = tvncdf(x2, rho_swapped);
    t(1) = isequal(p1, p2);
catch
    t = false;
end
T = all(t);
%@eof:3

%@test:4
% test with high correlations
try
    x = [1, 1, 1];
    rho = [0.8; 0.7; 0.9];
    p = tvncdf(x, rho);
    % result should be valid probability
    t(1) = p >= 0 && p <= 1;
    % should be close to univariate Phi(1) for high positive correlations
    t(2) = p > 0.7 && p < normcdf(1);
catch
    t = false(2, 1);
end
T = all(t);
%@eof:4

%@test:5
% test multiple points at once
try
    x = [0, 0, 0; 1, 1, 1; -1, -1, -1];
    rho = [0.3; 0.2; 0.1];
    p = tvncdf(x, rho);
    t(1) = length(p) == 3;
    t(2) = all(p >= 0 & p <= 1);
    t(3) = p(2) > p(1) && p(1) > p(3); % monotonicity
catch
    t = false(3, 1);
end
T = all(t);
%@eof:5

%@test:6
% test with negative correlations
try
    x = [0, 0, 0];
    rho = [-0.3; -0.2; 0.1];
    p = tvncdf(x, rho);
    % result should be valid probability, less than 0.125 due to negative correlations
    t(1) = p >= 0 && p <= 1;
    t(2) = p < 0.093;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:6


end % tvncdf


function integrand = tvnIntegrand(theta, b1, bj, bk, rho_j1, rho_k1, rho_32)
% Integrand for trivariate normal CDF computation
% Integrand is exp(-(b1^2 + bj^2 - 2*b1*bj*sin(theta))/(2*cos(theta)^2)) * Phi(...)

sintheta = sin(theta);
cossqtheta = cos(theta) .^ 2;
expon = ((b1 * sintheta - bj) .^ 2 ./ cossqtheta + b1 .^ 2) / 2;

sinphi = sintheta .* rho_k1 ./ rho_j1;
numeru = bk .* cossqtheta - b1 .* (sinphi - rho_32 .* sintheta) ...
                          - bj .* (rho_32 - sintheta .* sinphi);
denomu = sqrt(cossqtheta .* (cossqtheta - sinphi .* sinphi ...
                           - rho_32 .* (rho_32 - 2 .* sintheta .* sinphi)));
phi = 0.5 * erfc(-(numeru ./ denomu) / sqrt(2));
integrand = exp(-expon) .* phi;

end % tvnIntegrand

