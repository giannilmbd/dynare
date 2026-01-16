@#include "fs2000.inc"

estimated_params;
alp, beta_pdf, 0.356, 0.02;
bet, beta_pdf, 0.993, 0.002;
gam, normal_pdf, 0.0085, 0.003;
mst, normal_pdf, 1.0002, 0.007;
rho, beta_pdf, 0.129, 0.05;
psi, beta_pdf, 0.65, 0.05;
del, beta_pdf, 0.01, 0.005;
stderr e_a, inv_gamma_pdf, 0.035449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
end;

varobs gp_obs gy_obs;

stoch_simul(order=1,periods=200);
forecast(periods=50);
options_.TeX=true;
shock_decomposition(forecast_type=unconditional,datafile=fsdat_simul,parameter_set=calibration) gy_obs gp_obs;

conditional_forecast_paths;
var gy_obs;
periods  1  2  3:5;
values   (gy_obs+0.01) (gy_obs-0.02) (gy_obs);
var gp_obs;
periods 1:5;
values  (gp_obs+0.05);
end;

conditional_forecast(periods=50, parameter_set=calibration, controlled_varexo=(e_a,e_m));
plot_conditional_forecast(periods=50) gy_obs gp_obs;
shock_decomposition(forecast_type=conditional,parameter_set=calibration) gy_obs gp_obs;


// Metropolis replications are too few, this is only for testing purpose
estimation(order=1,datafile=fsdat_simul,nobs=192,mh_replic=2000,mh_nblocks=1,mh_jscale=0.8,forecast=50,silent_optimizer,plot_priors=0) gy_obs gp_obs;
shock_decomposition(forecast_type=unconditional) gy_obs gp_obs;

conditional_forecast_paths;
var gy_obs;
periods  1  2  3:5;
values   (gy_obs+0.01) (gy_obs-0.02) (gy_obs);
var gp_obs;
periods 1:5;
values  (gp_obs+0.05);
end;

conditional_forecast(periods=50, parameter_set=posterior_mean, controlled_varexo=(e_a,e_m));
plot_conditional_forecast(periods=50) gy_obs gp_obs;
shock_decomposition(forecast_type=conditional) gy_obs gp_obs;
collect_latex_files;

[status, cmdout]=system(['pdflatex -halt-on-error -interaction=nonstopmode ' M_.fname '_TeX_binder.tex']);
if status
    cmdout
    error('TeX-File did not compile.')
end
