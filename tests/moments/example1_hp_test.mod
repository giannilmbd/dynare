@#include "example1_common.inc"

// Parameter values
alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

steady(solve_algo=4,maxit=1000);

stoch_simul(order=1,nofunctions,hp_filter=1600,irf=0,filtered_theoretical_moments_grid=8192);

total_var_filtered=diag(oo_.var);
oo_filtered_all_shocks=oo_;

stoch_simul(order=1,nofunctions,hp_filter=0,periods=2500000,nomoments);
send_endogenous_variables_to_workspace;
options_.nomoments=0;
oo_unfiltered_all_shocks=oo_;

[junk, y_filtered]=sample_hp_filter(y,1600);
[junk, c_filtered]=sample_hp_filter(c,1600);
[junk, k_filtered]=sample_hp_filter(k,1600);
[junk, a_filtered]=sample_hp_filter(a,1600);
[junk, h_filtered]=sample_hp_filter(h,1600);
[junk, b_filtered]=sample_hp_filter(b,1600);

verbatim;
total_std_all_shocks_filtered_sim=std([y_filtered c_filtered k_filtered a_filtered h_filtered b_filtered]);
cov_filtered_all_shocks=cov([y_filtered c_filtered k_filtered a_filtered h_filtered b_filtered]);
acf = zeros(6);
[junk, acf(:,1)] = sample_autocovariance([y_filtered ],5);
[junk, acf(:,2)] = sample_autocovariance([c_filtered ],5);
[junk, acf(:,3)] = sample_autocovariance([k_filtered ],5);
[junk, acf(:,4)] = sample_autocovariance([a_filtered ],5);
[junk, acf(:,5)] = sample_autocovariance([h_filtered ],5);
[junk, acf(:,6)] = sample_autocovariance([b_filtered ],5);
autocorr_filtered_all_shocks=acf(2:end,:)';
end;

shocks;
var e; stderr 0;
var u; stderr 0.009;
var e, u = phi*0.009*0;
end;

stoch_simul(order=1,nofunctions,hp_filter=1600,irf=0,periods=0);

total_var_filtered_one_shock=diag(oo_.var);
oo_filtered_one_shock=oo_;

stoch_simul(order=1,nofunctions,hp_filter=0,periods=2500000,nomoments);
send_endogenous_variables_to_workspace;
oo_unfiltered_one_shock=oo_;

[junk, y_filtered]=sample_hp_filter(y,1600);
[junk, c_filtered]=sample_hp_filter(c,1600);
[junk, k_filtered]=sample_hp_filter(k,1600);
[junk, a_filtered]=sample_hp_filter(a,1600);
[junk, h_filtered]=sample_hp_filter(h,1600);
[junk, b_filtered]=sample_hp_filter(b,1600);

verbatim;
total_std_one_shock_filtered_sim=std([y_filtered c_filtered k_filtered a_filtered h_filtered b_filtered]);
cov_filtered_one_shock=cov([y_filtered c_filtered k_filtered a_filtered h_filtered b_filtered]);
acf = zeros(6);
[junk, acf(:,1)] = sample_autocovariance([y_filtered ],5);
[junk, acf(:,2)] = sample_autocovariance([c_filtered ],5);
[junk, acf(:,3)] = sample_autocovariance([k_filtered ],5);
[junk, acf(:,4)] = sample_autocovariance([a_filtered ],5);
[junk, acf(:,5)] = sample_autocovariance([h_filtered ],5);
[junk, acf(:,6)] = sample_autocovariance([b_filtered ],5);
autocorr_filtered_one_shock=acf(2:end,:)';
end;

if max(abs((1-(total_std_one_shock_filtered_sim.^2)./(total_std_all_shocks_filtered_sim.^2))*100-oo_filtered_all_shocks.variance_decomposition(:,1)'))>2
    error('Variance Decomposition wrong')
end

if max(max(abs(oo_filtered_all_shocks.var-cov_filtered_all_shocks)))>1e-4;
    error('Covariance wrong')
end

if max(max(abs(oo_filtered_one_shock.var-cov_filtered_one_shock)))>5e-5;
    error('Covariance wrong')
end

for ii=1:options_.ar
    autocorr_model_all_shocks(:,ii)=diag(oo_filtered_all_shocks.autocorr{ii});
    autocorr_model_one_shock(:,ii)=diag(oo_filtered_one_shock.autocorr{ii});
end

if max(max(abs(autocorr_model_all_shocks-autocorr_filtered_all_shocks)))>1e-2;
    error('Covariance wrong')
end

if max(max(abs(autocorr_model_one_shock-autocorr_filtered_one_shock)))>1e-2;
    error('Covariance wrong')
end
