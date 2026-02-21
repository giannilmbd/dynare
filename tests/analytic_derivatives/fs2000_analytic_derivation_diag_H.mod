// Tests the analytic_derivation option

@#include "fs2000.inc"

shocks;
var e_a; stderr 0.02;
var e_m; stderr 0.02;
var gy_obs= 0.005^2;
var gp_obs= 0.005^2;
% corr gp_obs, gy_obs= 0.05;
end;

stoch_simul(order=1,periods=500, irf=0,nomoments,noprint);

send_endogenous_variables_to_workspace;
% stoch_simul does not add measurement error; draw it from M_.H
% varobs order: gp_obs (1), gy_obs (2)
verbatim;
T = length(gy_obs);
me = chol(M_.H, 'lower') * randn(size(M_.H, 1), T);
gp_obs = gp_obs + me(1,:)';
gy_obs = gy_obs + me(2,:)';
end;
save('my_data.mat','gp_obs','gy_obs');


estimated_params;
alp, 0.33;
stderr e_a, 0.02, 0, inf;
stderr e_m, 0.02, 0, inf;
corr e_a, e_m, 0 ,-1, 1;
stderr gy_obs, 0.005, 0, inf;
stderr gp_obs, 0.005, 0, inf;
% corr gp_obs, gy_obs, 0.05, -1, 1;
end;


varobs gp_obs gy_obs;

options_.solve_tolf = 1e-12;

estimation(order=1,mode_compute=4,silent_optimizer,analytic_derivation,kalman_algo=1,datafile=my_data,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,prior_trunc=0,frequentist_smoother=false);
if (isoctave && user_has_octave_package('optim', '1.6')) || (~isoctave && user_has_matlab_license('optimization_toolbox'))
    estimation(order=1,mode_compute=1,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false
    %,optim = ('CheckGradients', true,'FiniteDifferenceType','central')
    );
    estimation(order=1,mode_compute=3,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
    estimation(order=1,mode_compute=101,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
end
if ~isoctave % This estimation randomly fails on Octave
estimation(order=1,mode_compute=5,silent_optimizer,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=2,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
end
estimation(order=1,mode_compute=4,silent_optimizer,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
estimation(order=1,mode_compute=4,silent_optimizer,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=2,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
options_.debug=1;
estimation(order=1,mode_compute=0,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,plot_priors=0,frequentist_smoother=false);
fval_ML_1=oo_.likelihood_at_initial_parameters;
estimation(order=1,mode_compute=0,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=2,datafile=my_data,nobs=500,mh_replic=0,plot_priors=0,frequentist_smoother=false);
fval_ML_2=oo_.likelihood_at_initial_parameters;
options_.analytic_derivation=false;
estimation(order=1,mode_compute=0,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,plot_priors=0,frequentist_smoother=false);
fval_ML_3=oo_.likelihood_at_initial_parameters;

if abs(fval_ML_1-fval_ML_2)>1e-5 || abs(fval_ML_1-fval_ML_3)>1e-5
    error('Likelihood does not match')
end
options_.debug=0;

estimated_params(overwrite);
alp, beta_pdf, 0.356, 0.02;
stderr e_a, inv_gamma_pdf, 0.035449, inf;
stderr e_m, inv_gamma_pdf, 0.008862, inf;
corr e_a, e_m, uniform_pdf, , ,-1, 1;
stderr gy_obs, inv_gamma_pdf, 0.005, inf;
stderr gp_obs, inv_gamma_pdf, 0.005, inf;
% corr gp_obs, gy_obs, uniform_pdf, , ,-1, 1;
end;

estimation(order=1,mode_compute=9,silent_optimizer,analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,prior_trunc=0,frequentist_smoother=false);
if (isoctave && user_has_octave_package('optim', '1.6')) || (~isoctave && user_has_matlab_license('optimization_toolbox'))
    estimation(order=1,mode_compute=1,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false
    %,optim = ('DerivativeCheck', 'on','FiniteDifferenceType','central')
    );
    estimation(order=1,mode_compute=3,silent_optimizer,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
    estimation(order=1,mode_compute=101,silent_optimizer,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
end
if ~isoctave % This estimation randomly fails on Octave
estimation(order=1,mode_compute=5,silent_optimizer,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=2,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
end
estimation(order=1,mode_compute=4,silent_optimizer,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
estimation(order=1,mode_compute=4,silent_optimizer,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=2,datafile=my_data,nobs=500,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0,frequentist_smoother=false);
options_.debug=1;
estimation(order=1,mode_compute=0,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,plot_priors=0,frequentist_smoother=false);
fval_Bayes_1=oo_.likelihood_at_initial_parameters;
estimation(order=1,mode_compute=0,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',analytic_derivation,kalman_algo=2,datafile=my_data,nobs=500,mh_replic=0,plot_priors=0,frequentist_smoother=false);
fval_Bayes_2=oo_.likelihood_at_initial_parameters;
options_.analytic_derivation=false;
estimation(order=1,mode_compute=0,mode_file='fs2000_analytic_derivation_diag_H/Output/fs2000_analytic_derivation_diag_H_mode',kalman_algo=1,datafile=my_data,nobs=500,mh_replic=0,plot_priors=0,frequentist_smoother=false);
fval_Bayes_3=oo_.likelihood_at_initial_parameters;

if abs(fval_Bayes_1-fval_Bayes_2)>1e-5 || abs(fval_Bayes_1-fval_Bayes_3)>1e-5
    error('Likelihood does not match')
end
