function [Sigma, Gamma, nu, Delta] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol)
% [Sigma, Gamma, nu, Delta] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol)
% -------------------------------------------------------------------------
% Reduces the dimension of a CSN distributed random variable
% Idea:
% (1) Representation of CSN distribution:
%     X ~ CSN(mu,Sigma,Gamma,nu,Delta) is equivalent to X = W|Z>=0 where
%     [W;Z] ~ N([mu;-nu], [Sigma, Sigma*Gamma'; Gamma*Sigma, Delta+Gamma*Sigma*Gamma']
% (2) Correlation introduces skewness:
%     Skewness in CSN is based on the correlation between W and Z:
%     - CSN and N are the same distribution if Gamma=0
%     - if correlation is small, then CSN very close to N
% (3) Prune dimensions if correlation in absolute value is below tol.
%     This reduces the overall skewness dimension to qq < q.
% -------------------------------------------------------------------------
% INPUTS
% - Sigma         [p by p]   scale parameter of the CSN distribution (does not equal covariance matrix unless Gamma=0)
% - Gamma         [q by p]   skewness shape parameter of the CSN distribution (if 0 then CSN reduces to Gaussian)
% - nu            [q by 1]   skewness conditioning parameter of the CSN distribution (enables closure of CSN distribution under conditioning, irrelevant if Gamma=0)
% - Delta         [q by q]   marginalization parameter of the CSN distribution (enables closure of CSN distribution under marginalization, irrelevant if Gamma=0)
% - tol           [double]    threshold value for correlation (in absolute value) below which correlated variables are pruned
% where p is the normal dimension and q the skewness dimension
% -------------------------------------------------------------------------
% OUTPUTS
% - Sigma         [p by p]    skewness shape parameter of the pruned CSN distribution
% - Gamma         [qq by p]   skewness shape parameter of the pruned CSN distribution
% - nu            [qq by 1]   skewness conditioning parameter of the pruned CSN distribution
% - Delta         [qq by q]   skewness marginalization parameter of the pruned CSN distribution
% where qq < q is the skewness dimension of the pruned distribution

% Copyright © 2022-2023 Gaygysyz Guljanov, Willi Mutschler
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

p = length(Sigma); % normal dimension

% create correlation matrix of conditional definition of CSN distributed variable
P = [Sigma, Sigma*Gamma'; Gamma*Sigma, Delta+Gamma*Sigma*Gamma']; P = 0.5*(P + P'); % covariance matrix
n = length(P);
normalization = diag(1./sqrt(diag(P)));
R = abs(normalization*P*normalization);

% prune dimensions in R if they are in absolute value lower than tol
R = R-diag(repelem(Inf, n));
logi2 = max(R(p+1:end, 1:p), [], 2);
logi2 = (logi2 < tol);

% prune dimensions in CSN parameters
Gamma(logi2, :) = [];
nu(logi2)       = [];
Delta(logi2, :) = [];
Delta(:, logi2) = [];


return % --*-- Unit tests --*--


%@test:1
% no pruning: high correlations should not be pruned with low tolerance
try
    Sigma = [1 0.5; 0.5 1];
    Gamma = [1 0; 0 1]; % creates high correlation
    nu = [0; 0];
    Delta = eye(2);
    tol = 0.01; % very low tolerance
    [Sigma_out, Gamma_out, nu_out, Delta_out] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol);
    % nothing should be pruned
    t(1) = isequal(size(Gamma_out), size(Gamma));
    t(2) = isequal(size(nu_out), size(nu));
    t(3) = isequal(size(Delta_out), size(Delta));
    t(4) = isequal(Sigma_out, Sigma);
catch
    t = false(4, 1);
end
T = all(t);
%@eof:1

%@test:2
% full pruning: very small Gamma should be fully pruned with reasonable tolerance
try
    Sigma = eye(3);
    Gamma = 1e-10 * ones(2, 3); % very small, creates negligible correlation
    nu = [0; 0];
    Delta = eye(2);
    tol = 0.1;
    [Sigma_out, Gamma_out, nu_out, Delta_out] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol);
    % everything should be pruned
    t(1) = isempty(Gamma_out);
    t(2) = isempty(nu_out);
    t(3) = isempty(Delta_out);
    t(4) = isequal(Sigma_out, Sigma);
catch
    t = false(4, 1);
end
T = all(t);
%@eof:2

%@test:3
% partial pruning: one dimension has high correlation, one has low
try
    Sigma = eye(2);
    Gamma = [1 0;      % high correlation with x1
             0.001 0]; % very low correlation
    nu = [0; 0];
    Delta = eye(2);
    tol = 0.1;
    [~, Gamma_out, nu_out, Delta_out] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol);
    % second row should be pruned, first kept
    t(1) = size(Gamma_out, 1) == 1;
    t(2) = isscalar(nu_out);
    t(3) = isequal(size(Delta_out), [1, 1]);
catch
    t = false(3, 1);
end
T = all(t);
%@eof:3

%@test:4
% verify output dimensions are consistent
try
    p = 4; q = 3;
    Sigma = eye(p) + 0.1 * ones(p);
    Gamma = randn(q, p);
    nu = zeros(q, 1);
    Delta = eye(q);
    tol = 0.5;
    [Sigma_out, Gamma_out, nu_out, Delta_out] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol);
    qq = size(Gamma_out, 1); % new skewness dimension
    t(1) = ( size(Sigma_out, 1) == p && size(Sigma_out, 2) == p );
    t(2) = ( size(Gamma_out, 2) == p );
    t(3) = ( length(nu_out) == qq );
    t(4) = ( size(Delta_out, 1) == qq && size(Delta_out, 2) == qq );
catch
    t = false(4, 1);
end
T = all(t);
%@eof:4

%@test:5
% univariate case: p=1, q=1
try
    Sigma = 2;
    Gamma = 0.5;
    nu = 0;
    Delta = 1;
    tol = 0.01; % low tolerance, should not prune
    [Sigma_out, Gamma_out, nu_out, Delta_out] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol);
    t(1) = (Sigma_out == Sigma);
    t(2) = (Gamma_out == Gamma);
    t(3) = (nu_out == nu);
    t(4) = (Delta_out == Delta);
catch
    t = false(4, 1);
end
T = all(t);
%@eof:5

%@test:6
% tolerance of 1 should prune everything (correlation always <= 1)
try
    Sigma = [1 0.3; 0.3 1];
    Gamma = [1 0.5; 0.2 1];
    nu = [0; 0];
    Delta = eye(2);
    tol = 1.0; % maximum possible correlation
    [Sigma_out, Gamma_out, nu_out, Delta_out] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol);
    % everything should be pruned
    t(1) = isempty(Gamma_out);
    t(2) = isempty(nu_out);
    t(3) = isempty(Delta_out);
catch
    t = false(3, 1);
end
T = all(t);
%@eof:6

%@test:7
% verify moments are approximately preserved after appropriate pruning
try
    mu = [0; 0];
    Sigma = [1 0.3; 0.3 1];
    Gamma = [2 0;        % high correlation
             0.001 0];   % negligible correlation
    nu = [0; 0];
    Delta = eye(2);
    tol = 0.1;
    [Sigma_out, Gamma_out, nu_out, Delta_out] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol);
    E_before = csn_mean(mu, Sigma, Gamma, nu, Delta, 'mvncdf');
    V_before = csn_variance(Sigma, Gamma, nu, Delta, 'mvncdf');
    E_after = csn_mean(mu, Sigma_out, Gamma_out, nu_out, Delta_out, 'mvncdf');
    V_after = csn_variance(Sigma_out, Gamma_out, nu_out, Delta_out, 'mvncdf');
    % moments should be very close since pruned dimension had negligible effect
    t(1) = norm(E_before - E_after, 'Inf') < 1e-3;
    t(2) = norm(V_before - V_after, 'Inf') < 1e-3;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:7

%@test:8
% zero tolerance should never prune (unless correlation is exactly 0)
try
    Sigma = eye(2);
    Gamma = [0.1 0; 0 0.1]; % small but non-zero
    nu = [0; 0];
    Delta = eye(2);
    tol = 0; % zero tolerance
    [~, Gamma_out, nu_out, Delta_out] = csn_prune_distribution(Sigma, Gamma, nu, Delta, tol);
    % nothing should be pruned (correlations are small but > 0)
    t(1) = isequal(size(Gamma_out), size(Gamma));
    t(2) = isequal(size(nu_out), size(nu));
    t(3) = isequal(size(Delta_out), size(Delta));
catch
    t = false(3, 1);
end
T = all(t);
%@eof:8