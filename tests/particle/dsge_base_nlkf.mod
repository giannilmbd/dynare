@#include "dsge_base.inc"

estimation(order=2,nograph,filter_algorithm=nlkf,number_of_particles=1000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation);
estimation(order=2,nograph,filter_algorithm=nlkf,number_of_particles=1000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,proposal_approximation=cubature);
estimation(order=2,nograph,filter_algorithm=nlkf,number_of_particles=1000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,proposal_approximation=montecarlo);
estimation(order=2,nograph,filter_algorithm=nlkf,number_of_particles=1000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,proposal_approximation=unscented,distribution_approximation=cubature);
estimation(order=2,nograph,filter_algorithm=nlkf,number_of_particles=5000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,proposal_approximation=unscented,distribution_approximation=montecarlo);
