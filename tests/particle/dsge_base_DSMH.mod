@#include "dsge_base.inc"

set_dynare_seed(10);
estimation(order=1,posterior_sampling_method='dsmh',posterior_sampler_options=('particles',500));
%  estimation(order=2,posterior_sampling_method='dsmh',posterior_sampler_options=('particles',1000));
%  estimation(order=3,posterior_sampling_method='dsmh',filter_algorithm=nlkf,proposal_approximation=montecarlo,number_of_particles=500,posterior_sampler_options=('particles',500));
