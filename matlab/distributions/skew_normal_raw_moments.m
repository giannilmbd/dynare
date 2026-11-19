function m = skew_normal_raw_moments(variance, skewness, D)
% m = skew_normal_raw_moments(variance, skewness, D)
% Raw moments E[e^d], d = 1..D, of a mean-zero skew normal shock with the given variance and standardized skewness,
% the distribution of a shock declared skewed in M_.Skew_e (csn_update_specification.m). Row vector, suitable for one
% row of M_.Raw_moments_e; the k-order MEX computes the same moments (skewNormalRawMoments, libkorder/tl).
%
% INPUTS
% - variance  [double] variance of e
% - skewness  [double] standardized skewness, |skewness| < 0.9953
% - D         [integer] highest moment order
%
% OUTPUTS
% - m         [double] 1 x D vector, m(d) = E[e^d]; m(1) = 0, m(2) = variance
%
% e = omega (delta (|U0| - sqrt(2/pi)) + sqrt(1 - delta^2) U1), U0, U1 iid N(0,1) (Azzalini and Capitanio 2014,
% Sect. 2.1), with delta and omega from (variance, skewness) as in csn_update_specification.m;
% E|U0|^i = 2^(i/2) Gamma((i+1)/2) / sqrt(pi).

% Copyright © 2026 Dynare Team
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

if variance < 0
    error('skew_normal_raw_moments: the variance must be non-negative')
end
delta = 0;
if skewness ~= 0
    s23 = abs(skewness)^(2/3);
    delta = sign(skewness) * sqrt(pi/2 * s23 / (s23 + ((4 - pi)/2)^(2/3)));
end
if abs(delta) >= 1
    error('skew_normal_raw_moments: the skewness must be smaller than 0.9953 in absolute value')
end
omega = sqrt(variance / (1 - 2/pi * delta^2));
c = sqrt(2/pi);
absmom = 2.^((0:D)/2) .* gamma(((0:D) + 1)/2) / sqrt(pi);   % absmom(i+1) = E|U0|^i
mb = zeros(1, D + 1); mb(1) = 1;                              % mb(n+1) = E[(sqrt(1-delta^2) U1)^n]
for n = 2:2:D
    mb(n + 1) = mb(n - 1) * (n - 1) * (1 - delta^2);
end
ma = zeros(1, D + 1);                                         % ma(j+1) = E[(delta (|U0| - c))^j]
for j = 0:D
    i = 0:j;
    ma(j + 1) = delta^j * sum(arrayfun(@(k) nchoosek(j, k), i) .* absmom(i + 1) .* (-c).^(j - i));
end
m = zeros(1, D);
for d = 1:D
    j = 0:d;
    m(d) = omega^d * sum(arrayfun(@(k) nchoosek(d, k), j) .* ma(j + 1) .* mb(d - j + 1));
end
m(1) = 0;
