function [csn_draws] = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta)
% [csn_draws] = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta)
% -------------------------------------------------------------------------
% Draws random numbers from multivariate closed skew normal distribution
% CSN(mu,Sigma,Gamma,nu,Delta) using accept-reject method
% -------------------------------------------------------------------------
% INPUTS
% - n          [scalar]   number of random draws
% - mu         [p by 1]   location parameter of the CSN distribution (does not equal expectation vector unless Gamma=0)
% - Sigma      [p by p]   scale parameter of the CSN distribution (does not equal covariance matrix unless Gamma=0)
% - Gamma      [q by p]   skewness shape parameter of the CSN distribution (if 0 then CSN reduces to Gaussian)
% - nu         [q by 1]   skewness conditioning parameter of the CSN distribution (enables closure of CSN distribution under conditioning, irrelevant if Gamma=0)
% - Delta      [q by q]   marginalization parameter of the CSN distribution (enables closure of CSN distribution under marginalization, irrelevant if Gamma=0)
% -------------------------------------------------------------------------
% OUTPUTS
% - csn_draws  [p by n]   random draws of the CSN(mu,Sigma,Gamma,nu,Delta) distribution

% Copyright © 2015 GPL v2 Dmitry Pavlyuk, Eugene Girtcius (rcsn.R function in the "csn" R package)
% Copyright © 2022-2023 Gaygysyz Guljanov
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

q      = length(nu);
v      = nan(q, n);
Sigma2 = Delta + Gamma * Sigma * Gamma';
L      = chol(Sigma2, 'lower');
for i = 1:n
    draws = -1;
    while sum(draws >= 0) ~= q
        draws = -nu + L * randn(q, 1);
    end
    v(:, i) = draws;
end

tmp = Sigma * Gamma' / Sigma2;
Exp = mu + tmp * (v + nu);
Var = Sigma - tmp * Gamma * Sigma;
L   = chol(Var, 'lower');
csn_draws = Exp + L * randn(length(mu), n);


return % --*-- Unit tests --*--

%@test:1
% verify output dimensions
try
    mu = [1; 2; 3];
    Sigma = eye(3);
    Gamma = [0.5 0 0; 0 0.5 0];
    nu = [0; 0];
    Delta = eye(2);
    n = 100;
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);
    t(1) = (size(X, 1) == length(mu));
    t(2) = (size(X, 2) == n);
catch
    t = false(2, 1);
end
T = all(t);
%@eof:1

%@test:2
% reproducibility: same seed should give same results
try
    mu = [0; 1];
    Sigma = [1 0.3; 0.3 1];
    Gamma = [1 0.5];
    nu = 0;
    Delta = 1;
    n = 50;

    rng(12345);
    X1 = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);
    rng(12345);
    X2 = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);

    t(1) = isequal(X1, X2);
catch
    t = false;
end
T = all(t);
%@eof:2

%@test:3
% special case: Gamma=0 should reduce to multivariate normal with mean mu and covariance Sigma
try
    mu = [2; -1];
    Sigma = [1 0.4; 0.4 1.5];
    Gamma = [0 0];
    nu = 0;
    Delta = 1;

    rng(111);
    n = 1e5;
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);

    % for Gamma=0, mean should be mu, covariance should be Sigma,
    % skewness should be 0, kurtosis should be 3
    empirical_mean = mean(X, 2);
    empirical_cov = cov(X');
    empirical_skew = skewness(X');
    empirical_kurt = kurtosis(X');

    t(1) = norm(empirical_mean - mu, 'Inf') < 0.006;
    t(2) = norm(empirical_cov - Sigma, 'Inf') < 0.02;
    t(3) = norm(empirical_skew, 'Inf') < 0.01;
    t(4) = norm(empirical_kurt-3, 'Inf') < 0.03;
catch
    t = false(4, 1);
end
T = all(t);
%@eof:3

%@test:4
% verify skewness: for positive Gamma, distribution should be right-skewed (positive skewness)
try
    mu = 0; Sigma = 1; nu = 0; Delta = 1;
    n = 1e5;
    Gamma = 5;
    rng(444);
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);
    t(1) = skewness(X) > 0;
    Gamma = -5;
    rng(555);
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);
    t(2) = skewness(X) < 0;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:4

%@test:5
% verify samples satisfy the skewness constraint: should have Gamma*x + nu < 0 rejected
% i.e., all accepted draws should come from the truncated region
try
    mu = [0; 0];
    Sigma = eye(2);
    Gamma = [1 1]; % constraint: x1 + x2 + nu < 0
    nu = 0;
    Delta = 1;

    rng(666);
    n = 1000;
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);

    % the CSN distribution is constructed such that the density is
    % proportional to phi(x) * Phi(Gamma*x + nu), so samples should
    % have a tendency towards Gamma*x + nu > 0 (i.e., x1 + x2 > 0)
    % Check that mean of Gamma*X + nu is positive
    constraint_vals = Gamma * X + nu;
    t(1) = mean(constraint_vals) > 0;
catch
    t = false;
end
T = all(t);
%@eof:5

%@test:6
% different dimensions: p=4, q=3
try
    p = 4; q = 3;
    mu = randn(p, 1);
    Sigma = eye(p) + 0.1 * ones(p);
    Gamma = randn(q, p);
    nu = zeros(q, 1);
    Delta = eye(q) + 0.05 * ones(q);
    rng(777);
    n = 1e6;
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);
    E_theoretical = csn_mean(mu, Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston');
    V_theoretical = csn_variance(Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston');
    E_empirical = mean(X, 2);
    V_empirical = cov(X');
    t(1) = norm(E_theoretical - E_empirical, 'Inf') < 0.02;
    t(2) = norm(dyn_vech(V_theoretical) - dyn_vech(V_empirical), 'Inf') < 0.05;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:6