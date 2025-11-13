@#include "dsge_base.inc"

estimation(order=2,nograph,filter_algorithm=cpf,number_of_particles=1000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation);
//estimation(order=2,nograph,filter_algorithm=cpf,number_of_particles=1000,mh_replic=0,mode_compute=0,mh_posterior_mode_estimation,particle_filter_options=('pruning',true),pruning);
