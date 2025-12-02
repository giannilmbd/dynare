% Various tests within Maximum Likelihood estimation
% with skew normally distributed shocks
% -------------------------------------------------------------------------

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

@#include "ireland2004_model.inc"

estimated_params;
stderr eta_a, 3.01673,  0,    10;
stderr eta_e, 0.02476,  0,    10;
stderr eta_z, 0.88648,  0,    10;
stderr eta_r, 0.27903,  0,    10;
OMEGA,        0.05809,  0,     1;
RHO_PI,       0.38648,  0,     1;
RHO_G,        0.39601,  0,     1;
RHO_X,        0.16539,  0,     1;
RHO_A,        0.90480,  0,     1;
RHO_E,        0.99067,  0,     1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1: Gaussian and pruned skewed kalman filter should run without error %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use Gaussian Kalman filter
estimation(datafile = 'ireland2004_post1980.csv'
          ,mode_compute = 8
          ,silent_optimizer
          ,lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
          ,kalman_algo = 1  % use Gaussian Kalman filter
          );
struct0 = cell2struct([struct2cell(oo_.mle_mode.parameters); struct2cell(oo_.mle_mode.shocks_std); oo_.posterior.optimization.log_density], ...
                      [fieldnames(oo_.mle_mode.parameters); ["stderr_" + fieldnames(oo_.mle_mode.shocks_std); "loglik"]]);
smoothed_shocks_0 = oo_.SmoothedShocks;

% use Pruned Skewed Kalman filter as Gaussian distribution is special case of CSN distribution
estimation(datafile = 'ireland2004_post1980.csv'
          ,mode_compute = 8
          ,silent_optimizer
          ,lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
          ,kalman_algo = 5  % use pruned skewed Kalman filter routine even in Gaussian case for comparability;
                            % note that Dynare's implementation is much faster because it switches to the steady-state Kalman filter which we have not implemented yet for the PSKF
          ,skewed_kalman_mvnlogcdf = 'mvncdf' % test whether the option is passed correctly, in Gaussian case there is no curse of increasing skewness dimension (as there is no need to compute Gaussian cdfs)
          );
struct5 = cell2struct([struct2cell(oo_.mle_mode.parameters); struct2cell(oo_.mle_mode.shocks_std); oo_.posterior.optimization.log_density], ...
                      [fieldnames(oo_.mle_mode.parameters); ["stderr_" + fieldnames(oo_.mle_mode.shocks_std); "loglik"]]);
smoothed_shocks_5 = oo_.SmoothedShocks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2: skewed_kalman_mvnlogcdf option should have been set to 'mvncdf' %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(options_.skewed_kalman.mvnlogcdf,'mvncdf')
    error('option ''skewed_kalman_mvnlogcdf'' was not set correctly')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 3: Estimation of Gaussian model with either Gaussian    %
%         or Pruned Skewed Kalman filter should be numerically %
%         equal to each other                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Comparing estimation results:\n');
fprintf('%-15s %15s %15s %15s %10s\n', 'Parameter', 'Gaussian KF', 'Skewed KF', 'Max Abs Dev', 'Status');
fprintf('%s\n', repmat('-', 1, 85));
status = nan(length(fieldnames(struct0)),1);
for i = 1:length(fieldnames(struct0))
    field_name = fieldnames(struct0);
    max_abs_dev = max(abs(struct0.(field_name{i}) - struct5.(field_name{i})));
    status(i) = (max_abs_dev > 1e-4);
    if status(i); status_str = 'FAIL'; else status_str = 'PASS'; end
    fprintf('%-15s %15f %15f %15.6e %10s\n', field_name{i}, struct0.(field_name{i}), struct5.(field_name{i}), max_abs_dev, status_str);
end
if any(status)
    error('Estimates with Gaussian KF (kalman_algo=0) are not numerically equal to Pruned Skewed KF (kalman_algo=5), something is wrong!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4: Smoothed shock values of Gaussian model with either Gaussian %
%         or Pruned Skewed Kalman smoother should be numerically       %
%         equal to each other                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n<strong>smoothed eta_r</strong>\n');
disp(array2table([smoothed_shocks_0.eta_r smoothed_shocks_5.eta_r (smoothed_shocks_0.eta_r-smoothed_shocks_5.eta_r)], ...
     'RowNames', "t="+string(1:options_.nobs), ...
     'VariableNames', ["kalman_algo=0", "kalman_algo=5", "dev"]));
fprintf('\n<strong>smoothed eta_z</strong>\n');
disp(array2table([smoothed_shocks_0.eta_z smoothed_shocks_5.eta_z (smoothed_shocks_0.eta_z-smoothed_shocks_5.eta_z)], ...
     'RowNames', "t="+string(1:options_.nobs), ...
     'VariableNames', ["kalman_algo=0", "kalman_algo=5", "dev"]));
fprintf('\n<strong>smoothed eta_e</strong>\n');
disp(array2table([smoothed_shocks_0.eta_e smoothed_shocks_5.eta_e (smoothed_shocks_0.eta_e-smoothed_shocks_5.eta_e)], ...
     'RowNames', "t="+string(1:options_.nobs), ...
     'VariableNames', ["kalman_algo=0", "kalman_algo=5", "dev"]));
fprintf('\n<strong>smoothed eta_a</strong>\n');
disp(array2table([smoothed_shocks_0.eta_a smoothed_shocks_5.eta_a (smoothed_shocks_0.eta_a-smoothed_shocks_5.eta_a)], ...
     'RowNames', "t="+string(1:options_.nobs), ...
     'VariableNames', ["kalman_algo=0", "kalman_algo=5", "dev"]));

if any(abs(smoothed_shocks_0.eta_r - smoothed_shocks_5.eta_r) > 1e-3) || ...
     any(abs(smoothed_shocks_0.eta_z - smoothed_shocks_5.eta_z) > 1e-3) || ...
       any(abs(smoothed_shocks_0.eta_e - smoothed_shocks_5.eta_e) > 1e-3) || ...
         any(abs(smoothed_shocks_0.eta_a - smoothed_shocks_5.eta_a) > 1e-3)
    error('Smoothed values should be numerically equal to each other')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 3: Estimation of skewness coefficients can be successfully done with %
%         Pruned Skewed Kalman filter, estimation should run without error  %
%         and estimation options should have be correctly changed           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OMEGA  = 0.15356;
RHO_PI = 0.28761;
RHO_G  = 0.33795;
RHO_X  = 0.28325;
RHO_A  = 0.91704;
RHO_E  = 0.98009;

shocks;
var eta_a; stderr 2.5232;
var eta_e; stderr 0.021228;
var eta_z; stderr 0.79002;
var eta_r; stderr 0.28384;
skew eta_z = -0.9950;
end;

estimated_params(overwrite);
OMEGA, 0.15356;
stderr eta_e, 0.021228;
skew eta_a, -0.1948;
skew eta_e, -0.2140;
skew eta_r,  0.8128;
end;

% should run without error
estimation(datafile = 'ireland2004_post1980.csv'
          ,mode_compute = 5
          ,silent_optimizer
          ,kalman_algo = 5  % use pruned skewed Kalman filter
          ,lik_init = 1     % initialize Kalman filter at Gaussian steady-state distribution
          ,skewed_kalman_prune_tol = 0.02
          ,skewed_kalman_rank_deficiency_transform
          ,skewed_kalman_mvnlogcdf = 'gaussian_log_mvncdf_mendell_elston'
          ,skewed_kalman_smoother_skip
          ,tex
          );

% check whether options have been correctly passed by the preprocessor
if (options_.skewed_kalman.prune_tol ~= 0.02) || ...
     (options_.skewed_kalman.rank_deficiency_transform ~= true) || ...
       ~strcmp(options_.skewed_kalman.mvnlogcdf, 'gaussian_log_mvncdf_mendell_elston') || ...
       (options_.skewed_kalman.skip_smoother ~= 1)
    error('skewed kalman options were not passed correctly')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4: Likelihood with skew normal model should be higher than with Gaussian model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loglik_csn = oo_.posterior.optimization.log_density; % expected to be approx 1215.68
if (loglik_csn-1215.6815424) > 1e-5
    error('Something wrong with CSN likelihood')
end
loglik_gauss = struct5.loglik; % expected to be approx 1207.56
if (loglik_gauss-1207.5618755) > 1e-5
    error('Something wrong with Gaussian likelihood')
end
lr_stat = 2 * (loglik_csn - loglik_gauss);
lr_pval = 1 - chi2cdf(lr_stat, 4); % 4 additional shock parameters with CSN (only if one would have estimated all parameters, here we don't do this to speed up testsuite)
lr_str = sprintf('test statistic of $%.2f$ and a p-value of $%.4f$.', lr_stat, lr_pval);
fprintf('\n%s\n* LIKELIHOOD RATIO TEST *\n%s\n', repmat('*',1,25), repmat('*',1,25));
fprintf('- unconstrained log-likelihood (CSN): %.4f\n', loglik_csn);
fprintf('- constrained log-likelihood (Gaussian): %.4f\n', loglik_gauss);
fprintf('- %s\n\n',lr_str);