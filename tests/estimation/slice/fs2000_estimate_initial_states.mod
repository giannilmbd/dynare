// Test estimate_initial_states_diffuse_prior option with slice sampler
// Based on fs2000 model

@#include "fs2000.inc"

stoch_simul(order=1,periods=200,nomoments,nofunctions,irf=0);
send_endogenous_variables_to_workspace;
datatomfile('fsdat_simul',{'gp_obs','gy_obs'})

estimation(order=1,datafile=fsdat_simul,nobs=192,
    mh_replic=50,mh_nblocks=1,mh_drop=0.2,
    mode_compute=0,
    posterior_sampling_method='slice',
    estimate_initial_states_endogenous_prior,
    smoother, filtered_vars, consider_all_endogenous
);
