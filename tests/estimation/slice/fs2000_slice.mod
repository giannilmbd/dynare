@#include "fs2000.inc"

//options_.posterior_sampling_method = 'slice';
estimation(order=1,datafile='../fsdat_simul',nobs=192,silent_optimizer,loglinear,mh_replic=50,mh_nblocks=2,mh_drop=0.2, //mode_compute=0,cova_compute=0,
posterior_sampling_method='slice'
);
// continue with rotated slice
estimation(order=1,datafile='../fsdat_simul',silent_optimizer,nobs=192,loglinear,mh_replic=100,mh_nblocks=2,mh_drop=0.5,load_mh_file,//mode_compute=0,
posterior_sampling_method='slice',
posterior_sampler_options=('rotated',1,'use_mh_covariance_matrix',1)
);
// continue with optimizer within slice
options_.gradient_epsilon = 1e-5;
estimation(order=1,datafile='../fsdat_simul',silent_optimizer,nobs=192,loglinear,mh_replic=1,mh_nblocks=2,mh_drop=0.5,load_mh_file,//mode_compute=0,
posterior_sampling_method='slice',
posterior_sampler_options=('maximize',1,'maximize_using_mh_bounds',1),
optim = ('MaxIter', 10, 'TolFun',1.e-4, 'robust',1));

options_.TeX=1;
generate_trace_plots(1:2);