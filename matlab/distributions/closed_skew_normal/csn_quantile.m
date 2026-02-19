function [q_alph] = csn_quantile(alph, mu, Sigma, Gamma, nu, Delta, mvnlogcdf, optim_fct)
% [q_alph] = csn_quantile(alph, mu, Sigma, Gamma, nu, Delta, mvnlogcdf, optim_fct)
% -------------------------------------------------------------------------
% Evaluates the quantile of a CSN(mu,Sigma,Gamma,nu,Delta) distributed random variable for the cumulative probability alph in the interval [0,1]
% Idea:
% - q_alph = min{x: F(x) >= alph}
% - Lemma 2.2.1 of chapter 2 of Genton (2004) - "Skew-elliptical distributions and their applications: a journey beyond normality"
%   gives the cdf, F_X(x) of a CSN distributed random variable X
% - We need to solve F_X(q_alpha) = alpha, to find alpha quantile q_alpha
% - Minimization problem: q_alpha = argmin(log(F_X(q_alpha)) - log(alpha))
% -------------------------------------------------------------------------
% INPUTS
% - alph        [0<alph<1]   cumulative probability alph in the interval [0,1]
% - mu          [p by 1]     location parameter of the CSN distribution (does not equal expectation vector unless Gamma=0)
% - Sigma       [p by p]     scale parameter of the CSN distribution (does not equal covariance matrix unless Gamma=0)
% - Gamma       [q by p]     skewness shape parameter of the CSN distribution (if 0 then CSN reduces to Gaussian)
% - nu          [q by 1]     skewness conditioning parameter of the CSN distribution (enables closure of CSN distribution under conditioning, irrelevant if Gamma=0)
% - Delta       [q by q]     marginalization parameter of the CSN distribution (enables closure of CSN distribution under marginalization, irrelevant if Gamma=0)
% - mvnlogcdf   [string]     name of function to compute log Gaussian cdf, possible values: 'gaussian_log_mvncdf_mendell_elston', 'mvncdf'
% - optim_fct   [string]     name of minimization function, possible values: 'fminunc', 'lsqnonlin'
% -------------------------------------------------------------------------
% OUTPUTS
% - q_alph      [p by 1]    alpha quantile vector of the CSN(mu,Sigma,Gamma,nu,Delta) distribution

% Copyright © 2022-2023 Gaygysyz Guljanov, Willi Mutschler
% Copyright © 2025-2026 Dynare Team
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

if nargin < 8
    optim_fct = 'lsqnonlin';
end
p = size(Sigma,1);
q = size(Delta, 1);
if p>1 || q>1
    warning('csn_quantile:multivariate', '''csn_quantile'' works reliably only in the univarite case,as there is no consensus on multivariate extensions of quantiles.\n          Results for quantiles are most likely wrong!')
end
Cov_mat = -Sigma * Gamma';

% evaluate the first log-CDF
Var_Cov1 = Delta - Gamma * Cov_mat;
if strcmp(mvnlogcdf,"gaussian_log_mvncdf_mendell_elston")
    normalization1 = diag(1 ./ sqrt(diag(Var_Cov1)));
    eval_point1    = normalization1 * (zeros(q, 1) - nu);
    Corr_mat1      = normalization1 * Var_Cov1 * normalization1;
    Corr_mat1      = 0.5*(Corr_mat1 + Corr_mat1');
    cdf1           = gaussian_log_mvncdf_mendell_elston(eval_point1, Corr_mat1);
elseif strcmp(mvnlogcdf,"mvncdf")
    cdf1 = log(mvncdf(zeros(1, q), nu', Var_Cov1));
end

% prepare helper function to get quantile
Var_Cov2 = [Sigma, Cov_mat; Cov_mat', Var_Cov1];

if strcmp(mvnlogcdf,"gaussian_log_mvncdf_mendell_elston")
    % requires zero mean and correlation matrix as inputs
    normalization2 = diag(1 ./ sqrt(diag(Var_Cov2)));
    Corr_mat2      = normalization2 * Var_Cov2 * normalization2;
    Corr_mat2      = 0.5 * (Corr_mat2 + Corr_mat2');
    if strcmp(optim_fct,'lsqnonlin')
        fun = @(x) gaussian_log_mvncdf_mendell_elston( (normalization2 * ([x; zeros(q, 1)] - [mu; nu])), Corr_mat2 ) - cdf1 - log(alph);
    elseif strcmp(optim_fct,'fminunc')
        fun = @(x) ( gaussian_log_mvncdf_mendell_elston( (normalization2 * ([x; zeros(q, 1)] - [mu; nu])), Corr_mat2 ) - cdf1 - log(alph))^2;
    elseif strcmp(optim_fct,'fminsearch')
        fun = @(x) ( gaussian_log_mvncdf_mendell_elston( (normalization2 * ([x; zeros(q, 1)] - [mu; nu])), Corr_mat2 ) - cdf1 - log(alph))^2;
    end
elseif strcmp(mvnlogcdf,"mvncdf")
    if strcmp(optim_fct,'lsqnonlin')
        fun = @(x) log(mvncdf([x', zeros(1, q)], [mu', nu'], Var_Cov2)) - cdf1 - log(alph);
    elseif strcmp(optim_fct,'fminunc')
        fun = @(x) (log(mvncdf([x', zeros(1, q)], [mu', nu'], Var_Cov2)) - cdf1 - log(alph))^2;
    elseif strcmp(optim_fct,'fminsearch')
        fun = @(x) (log(mvncdf([x', zeros(1, q)], [mu', nu'], Var_Cov2)) - cdf1 - log(alph))^2;
    end
end

% run minimization
optim_opt = optimset('MaxIter', 10000, 'MaxFunEvals', 10000, 'Display', 'off');
if strcmp(optim_fct,'lsqnonlin')
    q_alph = lsqnonlin(fun, mu, [], [], optim_opt);
elseif strcmp(optim_fct,'fminunc')
    q_alph = fminunc(fun, mu, optim_opt);
elseif strcmp(optim_fct,'fminsearch')
    q_alph = fminsearch(fun, mu, optim_opt);
end


return % --*-- Unit tests --*--


%@test:1
% univariate test: Gamma=0 reduces to standard normal, median should be at mu
try
    mu = 2.5; Sigma = 1.0; Gamma = 0; nu = 0; Delta = 1;
    q_05_mvncdf = csn_quantile(0.5, mu, Sigma, Gamma, nu, Delta, 'mvncdf', 'lsqnonlin');
    q_05_mendell = csn_quantile(0.5, mu, Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston', 'lsqnonlin');
    % for symmetric distribution (Gamma=0), the median should equal mu
    t(1) = abs(q_05_mvncdf - mu) < 1e-6;
    t(2) = abs(q_05_mendell - mu) < 1e-6;
catch
    t = [false, false];
end
T = all(t);
%@eof:1

%@test:2
% univariate test: for standard normal, check multiple quantiles
try
    mu = 0; Sigma = 1; Gamma = 0; nu = 2; Delta = 5; % nu and Delta are irrelevant
    % check 0.025 quantile (should be approximately -1.96)
    q_025 = csn_quantile(0.025, mu, Sigma, Gamma, nu, Delta, 'mvncdf', 'lsqnonlin');
    % check 0.975 quantile (should be approximately 1.96)
    q_975 = csn_quantile(0.975, mu, Sigma, Gamma, nu, Delta, 'mvncdf', 'lsqnonlin');
    t(1) = abs(q_025 - norminv(0.025, 0, 1)) < 1e-4;
    t(2) = abs(q_975 - norminv(0.975, 0, 1)) < 1e-4;
catch
    t = [false, false];
end
T = all(t);
%@eof:2

%@test:3
% univariate skew normal test: verify consistency of quantile with CDF
try
    xi = 2; omega = 0.1; alpha = -4;
    mu = xi; Sigma = omega^2; Gamma = alpha/omega; nu = 0; Delta = 1;
    alph = 0.25;
    % compute quantile
    q_alph = csn_quantile(alph, mu, Sigma, Gamma, nu, Delta, 'mvncdf', 'lsqnonlin');
    % verify by computing CDF at the quantile (should equal alph)
    Cov_mat = -Sigma * Gamma';
    Var_Cov1 = Delta - Gamma * Cov_mat;
    Var_Cov2 = [Sigma, Cov_mat; Cov_mat', Var_Cov1];
    cdf1 = mvncdf(0, nu', Var_Cov1);
    cdf_at_q = mvncdf([q_alph', 0], [mu', nu'], Var_Cov2) / cdf1;
    t(1) = abs(cdf_at_q - alph) < 1e-6;
catch
    t = false;
end
T = all(t);
%@eof:3

%@test:4
% simulation-based test: compare theoretical quantiles with empirical quantiles from random samples (univariate)
try
    xi = 2; omega = 0.1; alpha = -4;
    mu = xi; Sigma = omega^2; Gamma = alpha/omega; nu = 0; Delta = 1;
    rng(789);
    n = 1e6; tol = 5e-3; % tolerance for empirical vs theoretical quantiles
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);
    % test multiple quantiles
    alphas = [0.25, 0.50, 0.75];
    test_results = true(length(alphas), 1);
    for i = 1:length(alphas)
        q_theoretical = csn_quantile(alphas(i), mu, Sigma, Gamma, nu, Delta, 'mvncdf', 'lsqnonlin');
        q_empirical = quantile(X, alphas(i));
        test_results(i) = abs(q_theoretical - q_empirical) < tol;
    end
    t(1) = all(test_results);
catch
    t = false;
end
T = all(t);
%@eof:4
