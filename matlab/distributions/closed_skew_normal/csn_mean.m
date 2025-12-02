function [E_csn] = csn_mean(mu, Sigma, Gamma, nu, Delta, mvnlogcdf)
% [E_csn] = csn_mean(mu, Sigma, Gamma, nu, Delta, mvnlogcdf)
% -------------------------------------------------------------------------
% Evaluates the unconditional expectation vector of a CSN(mu,Sigma,Gamma,nu,Delta) distributed random variable based on expressions derived in
% - Dominguez-Molina, Gonzalez-Farias, Gupta (2003, sec. 2.4) - The multivariate closed skew normal distribution
% - Rodenburger (2015) - On Approximations of Normal Integrals in the Context of Probit based Discrete Choice Modeling (2015)
% -------------------------------------------------------------------------
% INPUTS
% - mu          [p by 1]   location parameter of the CSN distribution (does not equal expectation vector unless Gamma=0)
% - Sigma       [p by p]   scale parameter of the CSN distribution (does not equal covariance matrix unless Gamma=0)
% - Gamma       [q by p]   skewness shape parameter of the CSN distribution (if 0 then CSN reduces to Gaussian)
% - nu          [q by 1]   skewness conditioning parameter of the CSN distribution (enables closure of CSN distribution under conditioning, irrelevant if Gamma=0)
% - Delta       [q by q]   marginalization parameter of the CSN distribution (enables closure of CSN distribution under marginalization, irrelevant if Gamma=0)
% - mvnlogcdf   [string]   name of function to compute log Gaussian cdf, possible values: 'gaussian_log_mvncdf_mendell_elston', 'mvncdf'
% -------------------------------------------------------------------------
% OUTPUTS
% - E_csn         [p by 1]   expectation vector of the CSN(mu,Sigma,Gamma,nu,Delta) distribution

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

q       = size(Delta, 1); % skewness dimension
derVect = nan(q, 1);      % initialize gradient vector
Delta1  = Delta + Gamma * Sigma * Gamma'; % var-cov matrix inside the expression "psi" of equation (1) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10)

% differentiate the normal cdf, steps are given in Rodenburger (2015, p.28-29)
for jj = 1:q
    % permutation of var-cov matrix
    logi   = true(q, 1); logi(jj, 1)=false;
    nu2    = nan(q, 1); nu2(1:q-1)=nu(logi); nu2(q)=nu(jj);
    Delta2 = nan(size(Delta1));
    Delta2(1:q-1, :) = Delta1(logi, :);
    Delta2(q, :) = Delta1(jj, :);
    Delta2 = [Delta2(:, logi), Delta2(:, jj)];

    % conditional mean and var-cov matrix
    nu2      = -nu2;
    condMean = (Delta2(1:q-1, q)*nu2(q))/Delta2(q, q);
    condVar  = Delta2(1:q-1, 1:q-1) - (Delta2(1:q-1, q)/Delta2(q, q))*Delta2(q, 1:q-1);

    % evaluate the gradient vector
    term1 = normpdf(nu2(q), 0, sqrt(Delta2(q, q)));
    if isempty(condVar)
        term2 = 1;
    else
        % different functions to evaluate log Gaussian cdf
        if mvnlogcdf == "gaussian_log_mvncdf_mendell_elston"
            % requires zero mean and correlation matrix as inputs
            normalization2 = diag(1./sqrt(diag(condVar)));
            eval_point2 = normalization2*(nu2(1:q-1) - condMean);
            Corr_mat2 = normalization2*condVar*normalization2; 
            Corr_mat2 = 0.5*(Corr_mat2 + Corr_mat2');
            term2 = exp(gaussian_log_mvncdf_mendell_elston(eval_point2, Corr_mat2));
        elseif mvnlogcdf == "mvncdf"
            term2 = mvncdf(nu2(1:q-1)', condMean', condVar);
        end
    end
    derVect(jj) = term1*term2;
end

% evaluate end result using equation (1) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10)
Delta1 = 0.5*(Delta1+Delta1');
if isempty(Gamma) || isempty(nu) || isempty(Delta1)
    if isempty(Gamma) && isempty(nu) && isempty(Delta1)
        term3 = 1;
    else
        error("Problem with Gamma, nu, Delta being empty / not empty")
    end
else
    if mvnlogcdf == "gaussian_log_mvncdf_mendell_elston"
        % requires zero mean and correlation matrix as inputs
        normalization3 = diag(1./sqrt(diag(Delta1)));
        eval_point3 = normalization3*(zeros(q, 1) - nu);
        Corr_mat3 = normalization3*Delta1*normalization3;
        Corr_mat3 = 0.5*(Corr_mat3 + Corr_mat3');
        term3 = exp(gaussian_log_mvncdf_mendell_elston(eval_point3, Corr_mat3));
    elseif mvnlogcdf == "mvncdf"
        term3 = mvncdf(zeros(1, q), nu', Delta1);
    end
end
E_csn = mu + Sigma * Gamma' * (derVect/term3); % equation (1) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10)


return % --*-- Unit tests --*--


%@test:1
% univariate test: univariate skew normal SN(xi, omega, alpha) is special case of closed skew normal CSN_{1,1}(xi, omega^2, alpha/omega, 0, 1)
% use closed-form formulas for skew normal distribution
try
    xi = 2; omega = 0.1; alpha = -4;
    delta = alpha / sqrt(1+alpha^2);
    E = xi + omega * delta * sqrt(2/pi);
    E_csn_mvncdf = csn_mean(xi, omega^2, alpha/omega, 0, 1, 'mvncdf');
    E_csn_mendell = csn_mean(xi, omega^2, alpha/omega, 0, 1, 'gaussian_log_mvncdf_mendell_elston');
    t(1) = norm(E - E_csn_mvncdf, 'Inf') < 1e-20;
    t(2) = norm(E - E_csn_mendell, 'Inf') < 1e-20;
catch
    t = [false, false];
end
T = all(t);
%@eof:1

%@test:2
% univariate test: use auxiliary parameter reparametrization such that closed-form formulas are available
try
    mu = 0.3; Sigma = 0.64; nu = 0; lambda = -0.89;
    Gamma = lambda * Sigma^(-1/2);
    Delta = (1-lambda^2) * eye(size(mu,1));
    E = mu + (sqrt(2/pi) * lambda * Sigma^(1/2)) * ones(size(mu,1),1);
    E_csn_mvncdf = csn_mean(mu, Sigma, Gamma, nu, Delta, 'mvncdf');
    E_csn_mendell = csn_mean(mu, Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston');
    t(1) = norm(E - E_csn_mvncdf, 'Inf') < 1e-20;
    t(2) = norm(E - E_csn_mendell, 'Inf') < 1e-20;
catch
    t = [false, false];
end
T = all(t);
%@eof:2

%@test:3
% multivariate test: use auxiliary parameter reparametrization such that closed-form formulas are available
try
    mu = [0.3455; -1.8613; 0.7765; -0.5964];
    Sigma = [ 0.0013 -0.0111  0.0116 -0.0089;
             -0.0111  0.1009 -0.2301  0.1014;
              0.0116 -0.2301  3.3198 -1.0618;
             -0.0089  0.1014 -1.0618  1.0830];
    nu = [0; 0; 0; 0];
    lambda = 0.89;
    Gamma = lambda * Sigma^(-1/2);
    Delta = (1-lambda^2) * eye(size(mu,1));
    E = mu + (sqrt(2/pi) * lambda * Sigma^(1/2)) * ones(size(mu,1),1);
    E_csn_mvncdf = csn_mean(mu, Sigma, Gamma, nu, Delta, 'mvncdf');
    E_csn_mendell = csn_mean(mu, Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston');
    t(1) = norm(E - E_csn_mvncdf, 'Inf') < 1e-14;
    t(2) = norm(E - E_csn_mendell, 'Inf') < 1e-14;
catch
    t = [false, false];
end
T = all(t);
%@eof:3

%@test:4
% simulation-based test: compare theoretical mean with empirical mean from random samples
try
    mu = [2; -1; 0];
    Sigma = [0.16  0.10  0.00;
             0.10  0.64 -0.05;
             0.00 -0.05  0.25];
    Gamma = [2  0   0.7;
             1 -0.3 0  ];
    nu = [2; -0.5];
    Delta = [1 0.12; 0.12 1];
    
    rng(123);
    n = 1e6; tol = 1e-3; % tolerance for empirical vs theoretical
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);
    E_theoretical_mvncdf = csn_mean(mu, Sigma, Gamma, nu, Delta, 'mvncdf');
    E_theoretical_mendell = csn_mean(mu, Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston');
    E_empirical = mean(X, 2);
    t(1) = all( abs(E_theoretical_mvncdf - E_empirical) < tol );
    t(2) = all( abs(E_theoretical_mendell - E_empirical) < tol );
catch
    t = false;
end
T = all(t);
%@eof:4