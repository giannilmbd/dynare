@#include "dsge_base.inc"

estimation(order=2,nograph,filter_algorithm=sis,number_of_particles=2000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation);
estimation(order=2,nograph,filter_algorithm=sis,number_of_particles=2000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,resampling_method=kitagawa);
estimation(order=2,nograph,filter_algorithm=sis,number_of_particles=2000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,resampling_method=stratified);
estimation(order=2,nograph,filter_algorithm=sis,number_of_particles=2000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,resampling_method=smooth);
estimation(order=2,nograph,filter_algorithm=sis,number_of_particles=5000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,resampling=generic);

estimation(order=2,nograph,filter_algorithm=sis,number_of_particles=2000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,particle_filter_options=('pruning',true),pruning);