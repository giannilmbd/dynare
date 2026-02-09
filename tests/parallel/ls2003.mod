
@#includepath "folder"
@#include "endovars.mod"
@#include "_ls2003.inc"

estimation(datafile=data_ca1,first_obs=8,nobs=79,silent_optimizer,mh_replic=0,nodisplay);
if options_.parallel ~=0
estimation(datafile=data_ca1,first_obs=8,nobs=79,mode_compute=0,nodisplay, mode_file='ls2003/Output/ls2003_mode', mh_nblocks=4, prefilter=1, mh_jscale=0.5, mh_replic=2000);
estimation(datafile=data_ca1,first_obs=8,nobs=79,mode_compute=0,nodisplay, mode_file='ls2003/Output/ls2003_mode', mh_nblocks=4,prefilter=1,mh_jscale=0.5,mh_replic=2000,bayesian_irf,load_mh_file,smoother,forecast=12, filtered_vars, filter_step_ahead=[1 2 3 4]) y y_s R pie dq pie_s de A y_obs pie_obs R_obs;
end