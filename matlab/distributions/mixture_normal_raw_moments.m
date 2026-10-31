function m = mixture_normal_raw_moments(w, mu, sd, D)
% m = mixture_normal_raw_moments(w, mu, sd, D)
% Raw moments E[e^d], d = 1..D, of the mean-zero shock e = x - E[x], with x a finite mixture of normals,
% x ~ sum_j w_j N(mu_j, sd_j^2). Row vector, suitable for one row of M_.Raw_moments_e.
%
% INPUTS
% - w     [double] vector of mixture weights (non-negative, summing to 1)
% - mu    [double] vector of component means
% - sd    [double] vector of component standard deviations
% - D     [integer] highest moment order
%
% OUTPUTS
% - m     [double] 1 x D vector, m(d) = E[e^d]; m(1) = 0, m(2) = variance of x
%
% E[(c + s Z)^d] = sum_j C(d,j) c^(d-j) s^j E[Z^j], Z ~ N(0,1), E[Z^j] = (j-1)!! for even j, 0 for odd j.

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

w = w(:); mu = mu(:); sd = sd(:);
if any(w < 0) || abs(sum(w) - 1) > 1e-12
    error('mixture_normal_raw_moments: weights must be non-negative and sum to 1')
end
if any(sd < 0)
    error('mixture_normal_raw_moments: standard deviations must be non-negative')
end
mbar = w'*mu;
zmom = zeros(1, D + 1);                     % zmom(j+1) = E[Z^j]
zmom(1) = 1;
for j = 2:2:D
    zmom(j + 1) = zmom(j - 1)*(j - 1);
end
m = zeros(1, D);
for d = 1:D
    for c = 1:numel(w)
        s = 0;
        for j = 0:d
            s = s + nchoosek(d, j)*(mu(c) - mbar)^(d - j)*sd(c)^j*zmom(j + 1);
        end
        m(d) = m(d) + w(c)*s;
    end
end
m(1) = 0;
