function V_csn = csn_variance(Sigma, Gamma, nu, Delta, mvnlogcdf)
% V_csn = csn_variance(Sigma, Gamma, nu, Delta, mvnlogcdf)
% -------------------------------------------------------------------------
% Evaluates the unconditional covariance matrix of a CSN(mu,Sigma,Gamma,nu,Delta)
% distributed random variable based on expressions derived in
% Dominguez-Molina, Gonzalez-Farias, Gupta (2003, sec. 2.4) - The multivariate closed skew normal distribution.
% -------------------------------------------------------------------------
% INPUTS
% - Sigma       [p by p]   scale parameter of the CSN distribution (does not equal covariance matrix unless Gamma=0)
% - Gamma       [q by p]   skewness shape parameter of the CSN distribution (if 0 then CSN reduces to Gaussian)
% - nu          [q by 1]   skewness conditioning parameter of the CSN distribution (enables closure of CSN distribution under conditioning, irrelevant if Gamma=0)
% - Delta       [q by q]   marginalization parameter of the CSN distribution (enables closure of CSN distribution under marginalization, irrelevant if Gamma=0)
% - mvnlogcdf   [string]   name of function to compute log Gaussian cdf, possible values: 'gaussian_log_mvncdf_mendell_elston', 'mvncdf'
% -------------------------------------------------------------------------
% OUTPUTS
% - V_csn       [p by p]   covariance matrix of the CSN(mu,Sigma,Gamma,nu,Delta) distribution

% Copyright © 2022-2023 Gaygysyz Guljanov
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

q = size(Delta, 1); % skewness dimension

% Eq. (1) and (3) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10-11)
Delta2 = Delta + Gamma * Sigma * Gamma';
Delta2 = 0.5 * (Delta2 + Delta2'); % ensure symmetry
% Evaluate first and second derivatives
first_der  = zeros(q, 1); % first derivatives
second_der = zeros(q, q); % second derivatives
for ii = 1:q
    [term1, term2, evalp_t1, covar_t1, evalp_t2, covar_t2, mult_matr] = gaussian_cdf_derivative(nu, Delta2, ii, mvnlogcdf);
    first_der(ii) = term1 * term2;
    for jj = 1:q
        if ii ~= jj
            which_element = jj - (ii < jj);
            [expr1, expr2] = gaussian_cdf_derivative(-evalp_t2, covar_t2, which_element, mvnlogcdf);
            second_der(ii, jj) = term1 * expr1 * expr2;
        else
            expr3 = -evalp_t1 / covar_t1 * term1 * term2;
            expr4 = zeros(1, q-1);
            for kk = 1:q-1
                [first_term, second_term] = gaussian_cdf_derivative(-evalp_t2, covar_t2, kk, mvnlogcdf);
                expr4(kk) = first_term * second_term;
            end
            expr5 = term1 * expr4 * (-mult_matr);
            second_der(jj, jj) = expr3 + expr5;
        end
    end
end

% denominator in psi and Lambda of eq. (1) and (3) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10-11)
if strcmp(mvnlogcdf,"gaussian_log_mvncdf_mendell_elston")
    % requires zero mean and correlation matrix as inputs
    normalization3 = diag(1 ./ sqrt(diag(Delta2)));
    eval_point3 = normalization3 * (zeros(q, 1) - nu);
    Corr_mat3 = normalization3 * Delta2 * normalization3;
    Corr_mat3 = 0.5 * (Corr_mat3 + Corr_mat3');
    term3 = exp(gaussian_log_mvncdf_mendell_elston(eval_point3, Corr_mat3));
elseif strcmp(mvnlogcdf,"mvncdf")
    term3 = mvncdf(zeros(1, q), nu', Delta2);
end

% definition of psi in eq (1) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 10)
psi = first_der / term3;
% definition of Lambda in eq (3) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 11)
Lambda = second_der / term3; Lambda = 0.5 * (Lambda + Lambda');
% eq (3) of Dominguez-Molina, Gonzalez-Farias, Gupta (2003, p. 11)
V_csn = Sigma + Sigma * Gamma' * Lambda * Gamma * Sigma - Sigma * Gamma' * (psi * psi') * Gamma * Sigma;


return % --*-- Unit tests --*--


%@test:1
% univariate test: univariate skew normal SN(xi, omega, alpha) is special case of closed skew normal CSN_{1,1}(xi, omega^2, alpha/omega, 0, 1)
% use closed-form formulas for skew normal distribution
try
    xi = 2; omega = 0.1; alpha = -4;
    delta = alpha / sqrt(1+alpha^2);
    V = omega^2 * (1 - 2/pi * delta^2);
    V_csn_mvncdf = csn_variance(omega^2, alpha/omega, 0, 1, 'mvncdf');
    V_csn_mendell = csn_variance(omega^2, alpha/omega, 0, 1, 'gaussian_log_mvncdf_mendell_elston');
    t(1) = norm(V - V_csn_mvncdf, 'Inf') < 1e-18;
    t(2) = norm(V - V_csn_mendell, 'Inf') < 1e-18;
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
    V = Sigma * (1 - 2/pi * lambda^2);
    V_csn_mvncdf = csn_variance(Sigma, Gamma, nu, Delta, 'mvncdf');
    V_csn_mendell = csn_variance(Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston');
    t(1) = norm(V - V_csn_mvncdf, 'Inf') < 1e-16;
    t(2) = norm(V - V_csn_mendell, 'Inf') < 1e-16;
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
    V = Sigma * (1 - 2/pi * lambda^2);
    V_csn_mvncdf = csn_variance(Sigma, Gamma, nu, Delta, 'mvncdf');
    V_csn_mendell = csn_variance(Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston');
    t(1) = norm(V - V_csn_mvncdf, 'Inf') < 1e-14;
    t(2) = norm(V - V_csn_mendell, 'Inf') < 1e-14;
catch
    t = [false, false];
end
T = all(t);
%@eof:3

%@test:4
% simulation-based test: compare theoretical variance with empirical variance from random samples
try
    mu = [2; -1; 0];
    Sigma = [0.16  0.10  0.00;
             0.10  0.64 -0.05;
             0.00 -0.05  0.25];
    Gamma = [2  0   0.7;
             1 -0.3 0  ];
    nu = [2; -0.5];
    Delta = [1 0.12; 0.12 1];

    rng(456);
    n = 1e6; tol = 1e-3; % tolerance for empirical vs theoretical
    X = rand_multivariate_csn(n, mu, Sigma, Gamma, nu, Delta);
    V_theoretical_mvncdf = dyn_vech(csn_variance(Sigma, Gamma, nu, Delta, 'mvncdf'));
    V_theoretical_mendell = dyn_vech(csn_variance(Sigma, Gamma, nu, Delta, 'gaussian_log_mvncdf_mendell_elston'));
    V_empirical = dyn_vech(cov(X'));
    t(1) = all( abs(V_theoretical_mvncdf - V_empirical) < tol );
    t(2) = all( abs(V_theoretical_mendell - V_empirical) < tol );
catch
    t = false;
end
T = all(t);
%@eof:4
