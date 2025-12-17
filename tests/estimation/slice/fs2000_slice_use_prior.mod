@#include "fs2000.inc"

//options_.posterior_sampling_method = 'slice';
estimation(order=1,datafile='../fsdat_simul',nobs=192,silent_optimizer,loglinear,mh_replic=50,mh_nblocks=2,mh_drop=0.2, //mode_compute=0,cova_compute=0,
posterior_sampling_method='slice',
posterior_sampler_options=('use_prior_draws',1)
);

options_.TeX=1;
generate_trace_plots(1:2);
