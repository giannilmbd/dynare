function p = bvncdf(x, mu, Sigma)
% p = bvncdf(x, mu, Sigma)
% -------------------------------------------------------------------------
% Bivariate normal cumulative distribution function (CDF).
% -------------------------------------------------------------------------
% p = bvncdf(x, mu, Sigma) computes the bivariate normal cumulative
% distribution function of x given mean parameter mu and covariance matrix Sigma.
% -------------------------------------------------------------------------
% Inputs:
% x       [N-by-2 matrix]   Each row is a point, columns are variables
% mu      [2-by-1 vector]   Mean vector (or scalar for common mean, or empty for zero mean)
% Sigma   [2-by-2 matrix]   Covariance matrix (positive definite)
% -------------------------------------------------------------------------
% Outputs:
% p       [N-by-1 vector]   Probability values in [0, 1]
% -------------------------------------------------------------------------
% Example:
% mu = [1; -1];
% Sigma = [0.9 0.4; 0.4 0.3];
% [X1, X2] = meshgrid(linspace(-1, 3, 25)', linspace(-3, 1, 25)');
% x = [X1(:), X2(:)];
% p = bvncdf(x, mu, Sigma);
% Z = reshape(p, 25, 25);
% surf(X1, X2, Z);
% -------------------------------------------------------------------------
% References:
% Drezner, Z. and G.O. Wesolowsky (1989) "On the Computation of the Bivariate Normal Integral".
%   Journal of Statistical Computation and Simulation 35(1-2), pages 101-107.
%   doi: 10.1080/00949659008811236
% Code adapted from the GNU Octave statistics package (version 1.7.7) by Andreas Bertsatos:
% https://github.com/gnu-octave/statistics/blob/main/inst/dist_fun/bvncdf.m
%
% See also mvncdf, normcdf, tvncdf.

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

% input validation
if size(x, 2) ~= 2
    error('bvncdf: x must be an N-by-2 matrix with each variable as a column vector');
end

% handle mu
if isempty(mu)
    mu = [0; 0];
elseif isscalar(mu)
    mu = [mu; mu];
elseif numel(mu) == 2
    mu = mu(:);
else
    error('bvncdf: mu must be a scalar or a two-element vector');
end

% handle Sigma
if isscalar(Sigma)
    % scalar Sigma means correlation coefficient
    Sigma = [1, Sigma; Sigma, 1];
elseif ~isequal(size(Sigma), [2, 2])
    error('bvncdf: Sigma must be either a scalar (i.e. correlation) or a 2-by-2 covariance matrix');
end

% check for symmetric positive definite covariance matrix
[~, err] = chol(Sigma);
if err ~= 0
    error('bvncdf: the covariance matrix Sigma is not positive definite and/or symmetric');
end

% standardize: center by mean and scale by standard deviation
dh = (x(:,1) - mu(1)) / sqrt(Sigma(1,1));
dk = (x(:,2) - mu(2)) / sqrt(Sigma(2,2));

% compute correlation coefficient
r = Sigma(1,2) / sqrt(Sigma(1,1) * Sigma(2,2));

% initialize output
p = NaN(size(dh));

% handle special cases for infinite integration limits
p(dh == Inf & dk == Inf) = 1; % both limits are infinite: P(X1 ≤ ∞, X2 ≤ ∞) = 1
p(dk == Inf & dh ~= Inf) = 0.5 * erfc(-dh(dk == Inf & dh ~= Inf) / sqrt(2)); % x2 → +∞ (finite x1): P(X1 ≤ x1, X2 ≤ ∞) = P(X1 ≤ x1) = Φ(dh)
p(dh == Inf & dk ~= Inf) = 0.5 * erfc(-dk(dh == Inf & dk ~= Inf) / sqrt(2)); % x1 → +∞ (finite x2): P(X1 ≤ ∞, X2 ≤ x2) = P(X2 ≤ x2) = Φ(dk)
p(dh == -Inf | dk == -Inf) = 0; % x1 → -∞ or x2 → -∞: P(X1 ≤ -∞, X2 ≤ x2) = P(X1 ≤ x1, X2 ≤ -∞) = 0

% compute for finite values
ind = (dh > -Inf & dh < Inf & dk > -Inf & dk < Inf);
if sum(ind) > 0
    p(ind) = calculate_bvncdf(-dh(ind), -dk(ind), r);
end


return % --*-- Unit tests --*--

%@test:1
% test output against Octave reference values (first 10 points)
try
    mu = [1, -1];
    Sigma = [0.9, 0.4; 0.4, 0.3];
    [X1, X2] = meshgrid(linspace(-1, 3, 25)', linspace(-3, 1, 25)');
    x = [X1(:), X2(:)];
    p = bvncdf(x, mu, Sigma);
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
%@eof:1

%@test:2
% test output against Octave reference values (last 10 points)
try
    mu = [1, -1];
    Sigma = [0.9, 0.4; 0.4, 0.3];
    [X1, X2] = meshgrid(linspace(-1, 3, 25)', linspace(-3, 1, 25)');
    x = [X1(:), X2(:)];
    p = bvncdf(x, mu, Sigma);
    p_expected = [0.8180695783608276, 0.8854485749482751, ...
                  0.9308108777385832, 0.9579855743025508, ...
                  0.9722897881414742, 0.9788150170059926, ...
                  0.9813597788804785, 0.9821977956568989, ...
                  0.9824283794464095, 0.9824809345614861]';
    t(1) = max(abs(p(616:625) - p_expected)) < 3e-16;
catch
    t = false;
end
T = all(t);
%@eof:2

%@test:3
% test with independent variables (identity covariance)
% for independent: P(X1 <= x1, X2 <= x2) = P(X1 <= x1) * P(X2 <= x2)
try
    x = [0.5 -0.3];
    mu = [0; 0];
    Sigma = eye(2);
    p = bvncdf(x, mu, Sigma);
    p_expected = normcdf(0.5) * normcdf(-0.3);
    t(1) = isequal(p, p_expected);
catch
    t = false;
end
T = all(t);
%@eof:3

%@test:4
% test at z = [0, 0] with exact formula: P(X1 <= 0, X2 <= 0) = 1/4 + arcsin(rho)/(2*pi)
try
    rho = 0.5;
    x = [0, 0];
    mu = [];
    Sigma = [1 rho; rho 1];
    p = bvncdf(x, mu, Sigma);
    p_expected = 0.25 + asin(rho) / (2*pi);
    t(1) = isequal(p, p_expected);
catch
    t = false;
end
T = all(t);
%@eof:4

%@test:5
% test with empty mu (should default to zero mean)
try
    x = [1, 1];
    Sigma = [1 0.3; 0.3 1];
    p1 = bvncdf(x, [], Sigma);
    p2 = bvncdf(x, [0; 0], Sigma);
    t(1) = isequal(p1, p2);
catch
    t = false;
end
T = all(t);
%@eof:5

%@test:6
% test with scalar mu (should replicate to both dimensions)
try
    x = [1, 1];
    Sigma = [1 0.3; 0.3 1];
    p1 = bvncdf(x, 0.5, Sigma);
    p2 = bvncdf(x, [0.5; 0.5], Sigma);
    t(1) = isequal(p1, p2);
catch
    t = false;
end
T = all(t);
%@eof:6

%@test:7
% test with scalar Sigma (correlation coefficient)
try
    x = [0.5, 0.25];
    mu = [0; 2];
    rho = 0.6;
    p1 = bvncdf(x, mu, rho);
    p2 = bvncdf(x, mu, [1 rho; rho 1]);
    t(1) = isequal(p1, p2);
catch
    t = false;
end
T = all(t);
%@eof:7

%@test:8
% test infinite limits:
% (1): x1 = Inf, x2 = Inf should give p = 1
% (2): x1 = -Inf should give p = 0
% (3): x2 = -Inf should give p = 0
% (4): x2 = Inf should give marginal CDF of x1
% (5): x1 = Inf should give marginal CDF of x2
try
    mu = [0; 0];
    Sigma = [1 0.5; 0.5 1];
    x = [Inf, Inf];
    t(1) = isequal(bvncdf(x, mu, Sigma), 1);
    x = [-Inf, 2];
    t(2) = isequal(bvncdf(x, mu, Sigma), 0);
    x = [1, -Inf];
    t(3) = isequal(bvncdf(x, mu, Sigma), 0);
    x = [0.5, Inf];
    t(4) = isequal(bvncdf(x, mu, Sigma), normcdf(0.5));
    x = [Inf, 0.5];
    t(5) = isequal(bvncdf(x, mu, Sigma), normcdf(0.5));
catch
    t = false;
end
T = all(t);
%@eof:8

%@test:9
% test high correlation (|r| > 0.925) which uses different algorithm branch
try
    rho = 0.98;
    x = [0.5, 0.5];
    mu = [0; 0];
    Sigma = [1 rho; rho 1];
    p = bvncdf(x, mu, Sigma);
    % result should be valid probability
    t(1) = ( p >= 0 && p <= 1 );
    % for high positive correlation, P(X1 <= 0.5, X2 <= 0.5) should be close to P(X1 <= 0.5)
    t(2) = ( p > 0.5 && p < normcdf(0.5) );
catch
    t = false(2, 1);
end
T = all(t);
%@eof:9

%@test:10
% test negative correlation
try
    rho = -0.7;
    x = [0, 0];
    mu = [];
    Sigma = [1 rho; rho 1];
    p = bvncdf(x, mu, Sigma);
    p_expected = 0.25 + asin(rho) / (2*pi);
    t(1) = isequal(p, p_expected);
catch
    t = false;
end
T = all(t);
%@eof:10

%@test:11
% test multiple points at once
try
    x = [0, 0; 1, 1; -1, -1];
    mu = [0; 0];
    Sigma = [1 0.3; 0.3 1];
    p = bvncdf(x, mu, Sigma);
    t(1) = ( length(p) == 3 );
    t(2) = ( all(p >= 0 & p <= 1) );
    t(3) = ( p(2) > p(1) && p(1) > p(3) ); % monotonicity
catch
    t = false(3, 1);
end
T = all(t);
%@eof:11

end % bvncdf


function p = calculate_bvncdf(dh, dk, r)
% core computation of bivariate normal CDF using Gauss-Legendre quadrature

% select Gauss-Legendre points and weights based on correlation magnitude
if abs(r) < 0.3
    lg = 3;
    % Gauss Legendre points and weights, n = 6
    w = [0.1713244923791705, 0.3607615730481384, 0.4679139345726904];
    x = [0.9324695142031522, 0.6612093864662647, 0.2386191860831970];
elseif abs(r) < 0.75
    lg = 6;
    % Gauss Legendre points and weights, n = 12
    w = [0.04717533638651177, 0.1069393259953183, 0.1600783285433464, ...
         0.2031674267230659, 0.2334925365383547, 0.2491470458134029];
    x = [0.9815606342467191, 0.9041172563704750, 0.7699026741943050, ...
         0.5873179542866171, 0.3678314989981802, 0.1252334085114692];
else
    lg = 10;
    % Gauss Legendre points and weights, n = 20
    w = [0.01761400713915212, 0.04060142980038694, 0.06267204833410906, ...
         0.08327674157670475, 0.1019301198172404, 0.1181945319615184, ...
         0.1316886384491766, 0.1420961093183821, 0.1491729864726037, ...
         0.1527533871307259];
    x = [0.9931285991850949, 0.9639719272779138, 0.9122344282513259, ...
         0.8391169718222188, 0.7463319064601508, 0.6360536807265150, ...
         0.5108670019508271, 0.3737060887154196, 0.2277858511416451, ...
         0.07652652113349733];
end

dim1 = ones(size(dh, 1), 1);
dim2 = ones(1, lg);
hk = dh .* dk;
bvn = zeros(size(dh));
phi_dh = 0.5 * erfc(dh / sqrt(2));
phi_dk = 0.5 * erfc(dk / sqrt(2));

if abs(r) < 0.925
    hs = (dh .* dh + dk .* dk) / 2;
    asr = asin(r);
    sn1 = sin(asr * (1 - x) / 2);
    sn2 = sin(asr * (1 + x) / 2);
    bvn = sum((dim1 * w) .* exp(((dim1 * sn1) .* (hk * dim2) - ...
              hs * dim2) ./ (1 - dim1 * (sn1 .^ 2))) + ...
              (dim1 * w) .* exp(((dim1 * sn2) .* (hk * dim2) - ...
              hs * dim2) ./ (1 - dim1 * (sn2 .^ 2))), 2) * ...
              asr / (4 * pi) + phi_dh .* phi_dk;
else
    twopi = 2 * pi;
    if r < 0
        dk = -dk;
        hk = -hk;
    end
    if abs(r) < 1
        as = (1 - r) * (1 + r);
        a = sqrt(as);
        bs = (dh - dk) .^ 2;
        c = (4 - hk) / 8;
        d = (12 - hk) / 16;
        asr = -(bs ./ as + hk) / 2;
        ind = asr > -100;
        bvn(ind) = a * exp(asr(ind)) .* (1 - (c(ind) .* (bs(ind) - as)) ...
                    .* (1 - d(ind) .* bs(ind) / 5) / 3 ...
                    + (c(ind) .* d(ind)) .* as .^ 2 / 5);
        ind = hk > -100;
        b = sqrt(bs);
        phi_ba = 0.5 * erfc((b/a) / sqrt(2));
        sp = sqrt(twopi) * phi_ba;
        bvn(ind) = bvn(ind) - (exp(-hk(ind) / 2) .* sp(ind)) ...
                   .* b(ind) .* (1 - c(ind) .* bs(ind) ...
                   .* (1 - d(ind) .* bs(ind) / 5) / 3);
        a = a / 2;
        for is = -1:2:1
            xs = (a + a * is * x) .^ 2;
            rs = sqrt(1 - xs);
            asr1 = -((bs * dim2) ./ (dim1 * xs) + hk * dim2) / 2;
            ind1 = (asr1 > -100);
            sp1 = (1 + (c * dim2) .* (dim1 * xs) .* ...
                  (1 + (d * dim2) .* (dim1 * xs)));
            ep1 = exp(-(hk * dim2) .* (1 - dim1 * rs) ./ ...
                      (2 * (1 + dim1 * rs))) ./ (dim1 * rs);
            bvn = bvn + sum(a .* (dim1 * w) .* exp(asr1 .* ind1) ...
                               .* (ep1 .* ind1 - sp1 .* ind1), 2);
        end
        bvn = -bvn / twopi;
    end
    if r > 0
        tmp = max(dh, dk);
        bvn = bvn + 0.5 * erfc(tmp / sqrt(2));
    elseif r < 0
        phi_dh = 0.5 * erfc(dh / sqrt(2));
        phi_dk = 0.5 * erfc(dk / sqrt(2));
        bvn = -bvn + max(0, phi_dh - phi_dk);
    end
end

p = max(0, min(1, bvn));

end % calculate_bvncdf
