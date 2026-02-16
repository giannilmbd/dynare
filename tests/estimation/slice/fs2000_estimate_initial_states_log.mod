// Test estimate_initial_states_diffuse_prior option with slice sampler
// Based on fs2000 model

@#include "fs2000.inc"


estimation(order=1,datafile='../fsdat_simul',nobs=192,loglinear,
    mh_replic=100,mh_nblocks=2,mh_drop=0.2,
    mode_compute=0,
    posterior_sampling_method='slice',
    estimate_initial_states_endogenous_prior,
    smoother, filtered_vars, consider_all_endogenous
);

generate_trace_plots(1:2);