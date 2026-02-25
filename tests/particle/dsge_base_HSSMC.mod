@#include "dsge_base.inc"
warning('off'); % shut off warnings to not clutter log and speed up test
estimation(order=1,posterior_sampling_method='hssmc',posterior_sampler_options=('particles',1000));
%  estimation(order=2,posterior_sampling_method='hssmc',posterior_sampler_options=('particles',1000));
%  estimation(order=3,posterior_sampling_method='hssmc',filter_algorithm=nlkf,proposal_approximation=montecarlo,number_of_particles=500,posterior_sampler_options=('particles',500));
