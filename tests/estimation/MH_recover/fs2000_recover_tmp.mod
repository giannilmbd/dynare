//Test mh_recover function for RW-MH

@#include "fs2000.common.inc"

estimation(order=1, datafile='../fsdat_simul',nobs=192, silent_optimizer,loglinear, mh_replic=1000, mh_nblocks=1, mh_jscale=0.8,posterior_sampler_options=('save_tmp_file',1));
copyfile([M_.dname filesep 'metropolis' filesep M_.dname '_mh1_blck1.mat'],[M_.dname '_mh1_blck1.mat'])
delete([M_.dname filesep 'metropolis' filesep M_.dname '_mh1_blck1.mat'])
copyfile([M_.fname '_mh_tmp_blck1.mat'],[M_.dname filesep 'metropolis' filesep M_.dname '_mh_tmp_blck1.mat'])

estimation(order=1, datafile='../fsdat_simul',mode_compute=0,mode_file='fs2000_recover_tmp/Output/fs2000_recover_tmp_mode', nobs=192, loglinear, mh_replic=1000, mh_nblocks=1, mh_jscale=0.8,mh_recover);

%check first unaffected chain
temp1=load([M_.dname '_mh1_blck1.mat']);
temp2=load([M_.dname filesep 'metropolis' filesep M_.dname '_mh1_blck1.mat']);

if max(max(abs(temp1.x2./temp2.x2-1)))>5e-3
    max(max(abs(temp1.x2./temp2.x2-1)))
    format long
    temp1.x2
    temp2.x2
    temp1.x2-temp2.x2
    error('Draws of unaffected chain are not the same')
end

if max(max(abs(temp1.feval_this_chain-temp2.feval_this_chain)))>1e-6
    temp1.feval_this_chain
    temp2.feval_this_chain
    error('feval_this_chain are not the same')
end

if max(max(abs(temp1.logpo2-temp2.logpo2)))>1e-2
    max(max(abs(temp1.logpo2-temp2.logpo2)))
    temp1.logpo2
    temp2.logpo2
    error('logpo2 is not the same')
end

if max(max(abs(temp1.accepted_draws_this_file-temp2.accepted_draws_this_file)))>1e-6
    temp1.accepted_draws_this_file
    temp2.accepted_draws_this_file
    error('accepted_draws_this_file are not the same')
end

if max(max(abs(temp1.accepted_draws_this_chain-temp2.accepted_draws_this_chain)))>1e-6
    temp1.accepted_draws_this_chain
    temp2.accepted_draws_this_chain
    error('accepted_draws_this_chain are not the same')
end