% Various tests for simulations with skew normally distributed shocks
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
var rAnnualized piAnnualized;
model;
rAnnualized = 4*rhat;
piAnnualized = 4*pihat;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1: check user-facing set_*_value functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set_param_value('OMEGA', 0.15356);
set_param_value('RHO_PI', 0.28761);
set_param_value('RHO_G', 0.33795);
set_param_value('RHO_X', 0.28325);
set_param_value('RHO_A', 0.91704);
set_param_value('RHO_E', 0.98009);
set_shock_stderr_value('eta_a',2.5232);
set_shock_stderr_value('eta_e',0.021228);
set_shock_stderr_value('eta_z',0.79002);
set_shock_stderr_value('eta_r',0.28384);
set_shock_skew_value('eta_a',-0.1948); % single shock input
set_shock_skew_value('eta_e','eta_e','eta_e',-0.21401); % triple shock input
set_shock_skew_value('eta_z',-0.99527);
set_shock_skew_value('eta_r',0.81275);

if ~isequal(M_.params, [0.99; 0.1; 0; 0; 0.91704; 0.98009; 0.15356; 0.28761; 0.33795; 0.28325; 1])
    error('set_param_value did not set values correctly')
end
if ~isequal(M_.Sigma_e,diag([2.5232^2; 0.021228^2; 0.79002^2; 0.28384^2]))
    error('set_shock_stderr_value did not set values correctly')
end
if ~isequal(M_.Skew_e(1,1,1),-0.1948) || ~isequal(M_.Skew_e(2,2,2),-0.21401) || ~isequal(M_.Skew_e(3,3,3),-0.99527) || ~isequal(M_.Skew_e(4,4,4),0.81275)
    error('set_shock_skew_value did not set values correctly')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2: compute quantiles of CSN distribution and use for IRFs %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% note that skew normal is special case of CSN
csn = csn_update_specification(M_.Sigma_e, M_.Skew_e);
csn_q16 = csn_quantile(0.16, csn.mu_e(4,1), csn.Sigma_e(4,4), csn.Gamma_e(4,4), 0, 1, 'mvncdf');
csn_q84 = csn_quantile(0.84, csn.mu_e(4,1), csn.Sigma_e(4,4), csn.Gamma_e(4,4), 0, 1, 'mvncdf');
if abs(csn_q16 - (-0.273121782371810)) > 1e-15 || abs(csn_q84 - (0.283270501374202)) > 1e-15
    error('CSN quantiles have not been computed correctly');
end

% Note on IRF computations:
% Dynare always computes irfs with respect to a one standard deviation shock,
% which is provided by the user in the shocks block.
% Therefore, we use the shocks block to set the stderr to the quantiles
% and adjust the auxiliary parameter SIGN_SHOCKS in the mod file to get the
% correct sign on shocks (whether positive or negative).

% showcase how to simulate negative shocks, i.e. 16th quantiles
SIGN_SHOCKS = -1; % flip sign of shocks for 16th quantiles as they are negative
shocks;
var eta_r; stderr (csn_q16);
end;
stoch_simul(order=1, periods=0, irf=15, nodecomposition, nomoments, nocorr, nofunctions, nograph) rAnnualized piAnnualized xhat;
irfs_csn_neg = oo_.irfs;

% showcase how to simulate positive shocks, i.e. 84th quantile 
SIGN_SHOCKS = 1; % don't flip sign of shocks for 84th quantiles as they are positive
shocks;
var eta_r; stderr (csn_q84);
end;
stoch_simul(order=1, periods=0, irf=15, nodecomposition, nomoments, nocorr, nofunctions, nograph) rAnnualized piAnnualized xhat;
irfs_csn_pos = oo_.irfs;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 3: simulate data with skew normally distributed shocks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#for DISTRIB in ["N", "SN"] // N: Gaussian, SN: skew normal distribution
shocks;
var eta_a = 0.1^2;
var eta_e = 0.2^2;
var eta_z = 0.3^2;
var eta_r = 0.4^2;
  @#if DISTRIB == "N"
skew eta_a = 0;
skew eta_e = 0;
skew eta_z = 0;
skew eta_r = 0;
  @#elseif DISTRIB == "SN"
skew eta_a = -0.2;
skew eta_e = -0.6;
skew eta_z = -0.9;
skew eta_r = +0.8;
  @#endif
end;

set_dynare_seed(132);
stoch_simul(order=1,periods=500000,irf=0,nodecomposition,nomoments,nocorr,nofunctions,nograph);
send_exogenous_variables_to_workspace;
exo_@{DISTRIB} = oo_.exo_simul;
hh_fig = dyn_figure(options_.nodisplay,'Name','Histogram of @{DISTRIB} shocks');
subplot(2,2,1); histogram(eta_a,'normalization','pdf'); title('\eta_a');
subplot(2,2,2); histogram(eta_e,'normalization','pdf'); title('\eta_e');
subplot(2,2,3); histogram(eta_z,'normalization','pdf'); title('\eta_z');
subplot(2,2,4); histogram(eta_r,'normalization','pdf'); title('\eta_r');
sgtitle(hh_fig.Name);

@#endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4: compare empirical statistics to theoretical moments for exogenous shocks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
statNames = {'Mean', 'Median', 'StdDev', 'Variance', 'Skewness', 'Kurtosis', 'Min', 'Max'};
stats_N = nan(M_.exo_nbr, length(statNames));
stats_SN = nan(M_.exo_nbr, length(statNames));
for i = 1:M_.exo_nbr
    @#for DISTRIB in ["N", "SN"]
    stats_@{DISTRIB}(i, 1) = mean(exo_@{DISTRIB}(:, i)); 
    stats_@{DISTRIB}(i, 2) = median(exo_@{DISTRIB}(:, i));
    stats_@{DISTRIB}(i, 3) = std(exo_@{DISTRIB}(:, i));
    stats_@{DISTRIB}(i, 4) = var(exo_@{DISTRIB}(:, i));
    stats_@{DISTRIB}(i, 5) = skewness(exo_@{DISTRIB}(:, i));
    stats_@{DISTRIB}(i, 6) = kurtosis(exo_@{DISTRIB}(:, i));
    stats_@{DISTRIB}(i, 7) = min(exo_@{DISTRIB}(:, i));
    stats_@{DISTRIB}(i, 8) = max(exo_@{DISTRIB}(:, i));
    @#endfor
end
fprintf('<strong>DESCRIPTIVE STATISTICS FOR NORMAL SHOCKS</strong>\n');
disp(array2table(stats_N, 'VariableNames', statNames, 'RowNames', M_.exo_names));
fprintf('<strong>DESCRIPTIVE STATISTICS FOR SKEW NORMAL SHOCKS</strong>\n');
disp(array2table(stats_SN, 'VariableNames', statNames, 'RowNames', M_.exo_names));

% error out if something is not close to theoretical value
if any((abs(mean(exo_N)) > 1e-3))
    error('Mean of normal shocks should be close to zero!')
end
if any(abs(mean(exo_SN)) > 1e-3)
    error('Mean of skew normal shocks should be close to zero!')
end
if any(abs(std(exo_N) - transpose(sqrt(diag(M_.Sigma_e)))) > 1e-3)
    error('Standard deviation of normal shocks should be close to theoretical values in M_.Sigme_e!')
end
if any(abs(std(exo_SN) - transpose(sqrt(diag(M_.Sigma_e)))) > 1e-3)
    error('Standard deviation of skew normal shocks should be close to theoretical values in M_.Sigme_e!')
end
if any(abs(skewness(exo_N)) > 1e-2)
    error('Skewness coefficient of normal shocks should be close to 0!')
end
if any(abs(skewness(exo_SN)-transpose(nonzeros(M_.Skew_e(:)))) > 1e-2)
    error('Skewness coefficient of skew normal shocks should be close to theoretical values in M_.Skew_e!')
end