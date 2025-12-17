function [term1, term2, evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr] = gaussian_cdf_derivative(mu, Sigma, which_element, mvnlogcdf)
% [term1, term2, evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr] = gaussian_cdf_derivative(mu, Sigma, which_element, mvnlogcdf)
% -------------------------------------------------------------------------
% Differentiates the multivariate normal cdf according to steps given in
% Rodenburger (2015, p.28-29) - On Approximations of Normal Integrals in the Context of Probit based Discrete Choice Modeling (2015)
% - First derivative is equal to term1 * term2
% - Output variables "evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr"
%   are used for differentiating the first derivative for the second time,
%   for details look at the double loop in the code of csn_variance.m
% -------------------------------------------------------------------------
% INPUTS
% - mu             [p by 1]   mean of normal distribution
% - Sigma          [p by p]   covariance matrix of normal distribution
% - which_element  [scalar]   indicator (1,...,p) for which variable in multivariate normal cdf the derivative is computed for
% - mvnlogcdf      [string]   name of function to compute log Gaussian cdf, possible values: 'gaussian_log_mvncdf_mendell_elston', 'mvncdf'
% -------------------------------------------------------------------------
% OUTPUTS
% - term1          [double]           evaluated normal PDF, see the first part of the derivative
% - term2          [double]           evaluated normal CDF, see the second part of the derivative
% - evalp_t1       [1 by 1]           evaluation point of normal PDF in term1
% - covar_t1       [1 by 1]           variance term of normal PDF in term1
% - evalp_t2       [(p-1) by 1]       evaluation point of normal CDF in term2
% - covar_t2       [(p-1) by (p-1)]   covariance matrix of normal CDF in term2
% - mult_matr      [(p-1) by 1]       auxiliary expression Sigma_12/Sigma_22, this expression shows up when evaluating conditional mean and conditional variance out of partitions of mean and variance of joint distribution

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

p = length(mu);
logic_vec = true(p, 1);
logic_vec(which_element, 1) = false;

% permute the matrices
mu2        = zeros(p, 1);
mu2(1:p-1) = mu(logic_vec);
mu2(p)     = mu(which_element);

Sigma2           = zeros(size(Sigma));
Sigma2(1:p-1, :) = Sigma(logic_vec, :);
Sigma2(p, :)     = Sigma(which_element, :);
Sigma2           = [Sigma2(:, logic_vec), Sigma2(:, which_element)];

% Conditional mean and conditional var-cov matrix
eval_point = -mu2;
mult_matr  = Sigma2(1:p-1, p) / Sigma2(p, p);
cond_mean  = mult_matr * eval_point(p);
cond_var   = Sigma2(1:p-1, 1:p-1) - mult_matr * Sigma2(p, 1:p-1);

evalp_t1 = eval_point(p);
covar_t1 = Sigma2(p, p);

evalp_t2 = eval_point(1:p-1) - cond_mean;
covar_t2 = cond_var;

% Evaluate the first term in the derivative
term1 = normpdf(evalp_t1, 0, sqrt(covar_t1));

% Evaluate the second term in the derivative
if isempty(cond_var)
    term2 = 1;
else
    if strcmp(mvnlogcdf,"gaussian_log_mvncdf_mendell_elston")
        stdnrd   = diag(1./sqrt(diag(cond_var)));
        haspl    = stdnrd * evalp_t2;
        Uytgesme = stdnrd * cond_var * stdnrd;
        Uytgesme = 0.5 * (Uytgesme + Uytgesme');
        term2    = exp(gaussian_log_mvncdf_mendell_elston(haspl, Uytgesme));
    elseif strcmp(mvnlogcdf,"mvncdf")
        term2 = mvncdf(eval_point(1:p-1)', cond_mean', cond_var);
    end
end


return % --*-- Unit tests --*--


%@test:1
% univariate case (p=1): term2 should be 1, derivative equals normpdf
try
    mu = 0.5; Sigma = 2.0;
    [term1, term2, evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr] = gaussian_cdf_derivative(mu, Sigma, 1, 'mvncdf');
    % for univariate case, term2 = 1
    t(1) = abs(term2 - 1) < 1e-14;
    % term1 should be normpdf(-mu, 0, sqrt(Sigma))
    t(2) = abs(term1 - normpdf(-mu, 0, sqrt(Sigma))) < 1e-14;
    % evalp_t2 and covar_t2 should be empty
    t(3) = isempty(evalp_t2) && isempty(covar_t2);
    % mult_matr should be empty
    t(4) = isempty(mult_matr);
    % covar_t1 should equal Sigma
    t(5) = abs(covar_t1 - Sigma) < 1e-14;
catch
    t = false(5, 1);
end
T = all(t);
%@eof:1

%@test:2
% bivariate case: numerical verification using finite differences
% note: the function computes d/d(-mu) Phi(-mu), so analytical_deriv = -numerical_deriv
try
    mu = [0.3; -0.7]; Sigma = [1.0 0.4; 0.4 1.5];
    h = 1e-6; % step size for finite difference
    
    % test derivative with respect to first element
    mu_plus = mu; mu_plus(1) = mu(1) + h;
    mu_minus = mu; mu_minus(1) = mu(1) - h;
    cdf_plus = mvncdf(zeros(1,2), mu_plus', Sigma);
    cdf_minus = mvncdf(zeros(1,2), mu_minus', Sigma);
    numerical_deriv1 = (cdf_plus - cdf_minus) / (2*h);
    [term1, term2] = gaussian_cdf_derivative(mu, Sigma, 1, 'mvncdf');
    analytical_deriv1 = term1 * term2;
    % function returns magnitude, numerical has opposite sign due to chain rule
    t(1) = abs(-numerical_deriv1 - analytical_deriv1) < 1e-5;
    
    % test derivative with respect to second element
    mu_plus = mu; mu_plus(2) = mu(2) + h;
    mu_minus = mu; mu_minus(2) = mu(2) - h;
    cdf_plus = mvncdf(zeros(1,2), mu_plus', Sigma);
    cdf_minus = mvncdf(zeros(1,2), mu_minus', Sigma);
    numerical_deriv2 = (cdf_plus - cdf_minus) / (2*h);
    [term1, term2] = gaussian_cdf_derivative(mu, Sigma, 2, 'mvncdf');
    analytical_deriv2 = term1 * term2;
    t(2) = abs(-numerical_deriv2 - analytical_deriv2) < 1e-5;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:2

%@test:3
% trivariate case: numerical verification and compare both mvnlogcdf methods
% note: the function computes d/d(-mu) Phi(-mu), so analytical_deriv = -numerical_deriv
try
    mu = [0.2; -0.5; 0.8]; Sigma = [1.0 0.3 0.1; 0.3 1.2 -0.2; 0.1 -0.2 0.9];
    h = 1e-6;
    
    test_results = true(3, 2);
    for which_el = 1:3
        % numerical derivative of d/d(mu) Phi(0 | mu, Sigma)
        mu_plus = mu; mu_plus(which_el) = mu(which_el) + h;
        mu_minus = mu; mu_minus(which_el) = mu(which_el) - h;
        cdf_plus = mvncdf(zeros(1,3), mu_plus', Sigma);
        cdf_minus = mvncdf(zeros(1,3), mu_minus', Sigma);
        numerical_deriv = (cdf_plus - cdf_minus) / (2*h);
        
        % analytical derivative using mvncdf (returns -d/d(mu))
        [term1_mv, term2_mv] = gaussian_cdf_derivative(mu, Sigma, which_el, 'mvncdf');
        analytical_deriv_mv = term1_mv * term2_mv;
        test_results(which_el, 1) = abs(-numerical_deriv - analytical_deriv_mv) < 1e-4;
        
        % analytical derivative using gaussian_log_mvncdf_mendell_elston
        [term1_me, term2_me] = gaussian_cdf_derivative(mu, Sigma, which_el, 'gaussian_log_mvncdf_mendell_elston');
        analytical_deriv_me = term1_me * term2_me;
        test_results(which_el, 2) = abs(-numerical_deriv - analytical_deriv_me) < 1e-4;
    end
    t(1) = all(test_results(:));
catch
    t = false;
end
T = all(t);
%@eof:3

%@test:4
% verify conditional distribution properties: covar_t2 should be conditional covariance
try
    mu = [1; 2; 3]; 
    Sigma = [4 1 0.5; 1 3 0.8; 0.5 0.8 2];
    which_el = 2;
    
    [~, ~, evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr] = gaussian_cdf_derivative(mu, Sigma, which_el, 'mvncdf');
    
    % after permutation, element 2 moves to last position
    % remaining elements are [1, 3] in positions [1, 2]
    idx_remain = [1, 3];
    Sigma_11 = Sigma(idx_remain, idx_remain);
    Sigma_12 = Sigma(idx_remain, which_el);
    Sigma_22 = Sigma(which_el, which_el);
    
    % conditional covariance: Sigma_11 - Sigma_12 * Sigma_22^(-1) * Sigma_12'
    expected_cond_var = Sigma_11 - (Sigma_12 / Sigma_22) * Sigma_12';
    
    t(1) = norm(covar_t2 - expected_cond_var, 'Inf') < 1e-14;
    t(2) = abs(covar_t1 - Sigma_22) < 1e-14;
    t(3) = norm(mult_matr - Sigma_12 / Sigma_22, 'Inf') < 1e-14;
catch
    t = false(3, 1);
end
T = all(t);
%@eof:4

%@test:5
% standard normal case: verify known values at mu=0
try
    mu = [0; 0]; Sigma = eye(2);
    
    [term1, term2] = gaussian_cdf_derivative(mu, Sigma, 1, 'mvncdf');
    % for standard normal at mu=0: term1 = normpdf(0) = 1/sqrt(2*pi)
    t(1) = abs(term1 - normpdf(0)) < 1e-14;
    % term2 = normcdf(0) = 0.5 (since independent and conditional mean is 0)
    t(2) = abs(term2 - 0.5) < 1e-10;
    
    % derivative should be normpdf(0) * normcdf(0) = 0.5/sqrt(2*pi)
    expected_deriv = normpdf(0) * normcdf(0);
    t(3) = abs(term1 * term2 - expected_deriv) < 1e-10;
catch
    t = false(3, 1);
end
T = all(t);
%@eof:5

%@test:6
% consistency: both mvnlogcdf methods should give same results
try
    mu = [-0.4; 0.6; 1.2; -0.3];
    Sigma = [2.0  0.5  0.2  0.1;
             0.5  1.5 -0.3  0.4;
             0.2 -0.3  1.8  0.2;
             0.1  0.4  0.2  1.2];
    
    test_results = true(4, 1);
    for which_el = 1:4
        [term1_mv, term2_mv] = gaussian_cdf_derivative(mu, Sigma, which_el, 'mvncdf');
        [term1_me, term2_me] = gaussian_cdf_derivative(mu, Sigma, which_el, 'gaussian_log_mvncdf_mendell_elston');
        % term1 should be identical (both use normpdf)
        test_results(which_el) = abs(term1_mv - term1_me) < 1e-14 && abs(term2_mv - term2_me) < 1e-3;
    end
    t(1) = all(test_results);
catch
    t = false;
end
T = all(t);
%@eof:6