@#include "example1_common.inc"

// Parameter values
alpha = 0.36;
rho   = 0.5;
tau   = 0.025;
beta  = 0.98;
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

stoch_simul(order=1,nofunctions,irf=0,bandpass_filter=[6 32],filtered_theoretical_moments_grid=8192);
oo_filtered_all_shocks_theoretical=oo_;

stoch_simul(order=1,nofunctions,periods=1000000);
oo_filtered_all_shocks_simulated=oo_;


shocks;
var e; stderr 0;
var u; stderr 0.009;
var e, u = phi*0.009*0;
end;

stoch_simul(order=1,nofunctions,irf=0,periods=0);

oo_filtered_one_shock_theoretical=oo_;

stoch_simul(order=1,nofunctions,periods=1000000);
oo_filtered_one_shock_simulated=oo_;


if max(abs((1-diag(oo_filtered_one_shock_simulated.var)./(diag(oo_filtered_all_shocks_simulated.var)))*100-oo_filtered_all_shocks_theoretical.variance_decomposition(:,1)))>2
    error('Variance Decomposition wrong')
end

if max(max(abs(oo_filtered_all_shocks_simulated.var-oo_filtered_all_shocks_theoretical.var)))>5e-4;
    error('Covariance wrong')
end

if max(max(abs(oo_filtered_one_shock_simulated.var-oo_filtered_one_shock_theoretical.var)))>5e-4;
    error('Covariance wrong')
end

for ii=1:options_.ar
    autocorr_model_all_shocks_simulated(:,ii)=diag(oo_filtered_all_shocks_simulated.autocorr{ii});
    autocorr_model_all_shocks_theoretical(:,ii)=diag(oo_filtered_all_shocks_theoretical.autocorr{ii});
    autocorr_model_one_shock_simulated(:,ii)=diag(oo_filtered_one_shock_simulated.autocorr{ii});
    autocorr_model_one_shock_theoretical(:,ii)=diag(oo_filtered_one_shock_theoretical.autocorr{ii});
end

if max(max(abs(autocorr_model_all_shocks_simulated-autocorr_model_all_shocks_theoretical)))>0.05;
    error('Correlation wrong')
end

if max(max(abs(autocorr_model_one_shock_simulated-autocorr_model_one_shock_theoretical)))>0.05;
    error('Correlation wrong')
end
