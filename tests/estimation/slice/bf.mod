// test slice sampler with conditional likelihood initial state estimation of OccBin
var q r rlag rnot u;    
 
varexo  epsu
epsr
;
parameters betap phip rhop rhor rhou rlb shock_scale_zlb sigmap;
betap   = 0.99;
phip    = 0.2;
qinit   = 0;
rinit   = 0;
rhop    = 0.5;
rhor    = 0.8; 
rhou    = 0; 
rlb     = -(1/betap-1);
shock_scale_zlb =1;
sigmap  = 5;
model;
    [name = 'Asset price']
    q = betap*(1-rhop)*q(1)+rhop*q(-1)-sigmap*r+u;
    
    [name = 'Shock process u']
    u = rhou*u(-1)+epsu;
    
    [name = 'Notional rate']
    rnot = (1-rhor)*phip*q+rhor*rlag(-1);
    
    [name = 'Observed interest rate',relax='zlb']
    r = rnot + epsr;
        
    [name = 'Observed interest rate',bind='zlb']
    r = rlb + shock_scale_zlb*epsr;
    
    [name = 'Lag term TR',bind='zlb']
    rlag = rnot;
    
    [name = 'Lag term TR',relax='zlb']
    rlag = r;
end;
steady;

shocks;
    var epsr;
    stderr 0.001;
    var epsu;
    stderr 0.15;
end;

occbin_constraints;
    name 'zlb'; bind rnot<=rlb;
end;
// ------------------------------------------------------------------------
//.......................... COMPUTATIONS .................................
// ------------------------------------------------------------------------
// Prepare data and observations
// ----------------------------
varobs q r;
// Priors
// -----------------------------
estimated_params;
        phip, 0.2, 0.01, 2, NORMAL_PDF, 0.2, 0.05;
        rhop, 0.5, 0.01000, 0.9999, BETA_PDF, 0.5, 0.2;
        rhor, 0.8, 0.01000, 0.9999, BETA_PDF, 0.8, 0.1;
		stderr epsu, GAMMA_PDF, 0.15, 0.015, 0, inf; 
        stderr epsr, GAMMA_PDF, 0.001, 0.0003, 0, inf;  
end;       

//occbin options
options_.occbin.likelihood.max_number_of_iterations = 30;
options_.occbin.simul.periodic_solution=true;
options_.occbin.smoother.periodic_solution=true;
options_.occbin.likelihood.periodic_solution=true;
// end occbin options

steady;
check;
// Estimation
// -----------------------------

estimation(estimate_initial_states_endogenous_prior,
    datafile='datafile', 
    order=1,
	use_univariate_filters_if_singularity_is_detected=0,
	mh_replic=20, 
	mode_compute=0, 
	posterior_sampling_method='slice', 
    posterior_sampler_options = ('save_iter_info_file', 0, 
                                 'draw_init_state_from_smoother',1,
                                 'draw_init_state_with_rotated_slice',0,
                                 'fast_likelihood_evaluation_for_rejection',1),
	mh_nblocks=1,
	mh_drop=0.5,
	sub_draws = 10,
	//bayesian_irf,
    filter_covariance, smoothed_state_uncertainty,
	filtered_vars, smoother,consider_all_endogenous);

      
