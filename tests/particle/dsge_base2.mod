// ALGORITHM
@#ifndef LINEAR_KALMAN
    @#define LINEAR_KALMAN = 0
@#endif

@#ifndef NONLINEAR_KALMAN
    @#define NONLINEAR_KALMAN = 0
@#endif

@#ifndef ALGO_SIR
    @#define ALGO_SIR = 1
@#endif

@#ifndef ALGO_APF
    @#define ALGO_APF = 0
@#endif

@#ifndef ALGO_GCF
    @#define ALGO_GCF = 0
@#endif

@#ifndef ALGO_GUF
    @#define ALGO_GUF = 0
@#endif

@#ifndef ALGO_ONLINE
    @#define ALGO_ONLINE = 0
@#endif

@#ifndef MCMC
    @#define MCMC = 0
@#endif

@#ifndef HSSMC
    @#define HSSMC = 0
@#endif

@#ifndef DSMH
    @#define DSMH = 0
@#endif

@#include "dsge_base.inc"

@#if LINEAR_KALMAN
	estimation(nograph,order=1,mode_compute=8,silent_optimizer,mh_replic=0,additional_optimizer_steps=[8 4],mode_check);
@#endif

@#if NONLINEAR_KALMAN
    estimation(nograph,order=2,filter_algorithm=nlkf,number_of_particles=10000,silent_optimizer,proposal_approximation=montecarlo,mode_compute=8,additional_optimizer_steps=[8,4],mh_replic=0);
%    estimation(nograph,order=3,filter_algorithm=nlkf,number_of_particles=10000,silent_optimizer,proposal_approximation=montecarlo,mode_compute=8,additional_optimizer_steps=[8,4],mh_replic=0);
	estimation(nograph,order=2,filter_algorithm=nlkf,proposal_approximation=cubature,silent_optimizer,mode_compute=8,additional_optimizer_steps=[4,8],mh_replic=0,mode_check);
%	estimation(nograph,order=3,filter_algorithm=nlkf,proposal_approximation=cubature,silent_optimizer,mode_compute=8,additional_optimizer_steps=[4,8],mh_replic=0,mode_check);
@#endif

@#if ALGO_SIR
	estimation(order=2,nograph,number_of_particles=1000,mh_replic=0,mode_compute=8,cova_compute=0);
%	estimation(order=2,nograph,number_of_particles=1000,mh_replic=0,silent_optimizer,mode_compute=8,additional_optimizer_steps=[8 8],cova_compute=0);
%	estimation(order=3,nograph,number_of_particles=100,mh_replic=0,mode_file='dsge_base2/Output/dsge_base2_mode.mat',mode_compute=8,cova_compute=0);
@#endif

@#if ALGO_APF
    estimation(order=2,nograph,filter_algorithm=apf,number_of_particles=10000,resampling=none,mh_replic=0,mode_compute=8,cova_compute=0);
    estimation(order=2,nograph,filter_algorithm=apf,number_of_particles=10000,resampling=none,mh_replic=0,mode_compute=8,additional_optimizer_steps=[8,8],silent_optimizer,cova_compute=0);
    estimation(order=3,nograph,filter_algorithm=apf,number_of_particles=10000,resampling=none,mh_replic=0,silent_optimizer,mode_compute=8,additional_optimizer_steps=[8,8],cova_compute=0);
@#endif

@#if ALGO_GCF
estimation(order=2,nograph,filter_algorithm=gf,proposal_approximation=montecarlo,distribution_approximation=montecarlo,number_of_particles=10000,mh_replic=0,mode_compute=8,additional_optimizer_steps=[8,8],silent_optimizer);
estimation(order=3,nograph,filter_algorithm=gf,proposal_approximation=montecarlo,distribution_approximation=montecarlo,number_of_particles=10000,mh_replic=0,mode_compute=8,additional_optimizer_steps=[8,8],silent_optimizer);

%estimation(order=2,nograph,filter_algorithm=gf,proposal_approximation=cubature,distribution_approximation=montecarlo,mh_replic=0,mode_compute=8,additional_optimizer_steps=[4,8],silent_optimizer);
%estimation(order=3,nograph,filter_algorithm=gf,proposal_approximation=cubature,distribution_approximation=montecarlo,mh_replic=0,mode_compute=8,additional_optimizer_steps=[4,8],silent_optimizer);

%estimation(order=2,nograph,filter_algorithm=gf,proposal_approximation=cubature,distribution_approximation=cubature,mh_replic=0,mode_compute=8,additional_optimizer_steps=[4,8],silent_optimizer);
%estimation(order=3,nograph,filter_algorithm=gf,proposal_approximation=cubature,distribution_approximation=cubature,mh_replic=0,mode_compute=8,additional_optimizer_steps=[4,8],silent_optimizer);

%estimation(order=2,nograph,filter_algorithm=gf,proposal_approximation=montecarlo,distribution_approximation=montecarlo,number_of_particles=100000,resampling=systematic,mh_replic=0,mode_compute=8,additional_optimizer_steps=[4,8],silent_optimizer);
%estimation(order=3,nograph,filter_algorithm=gf,proposal_approximation=montecarlo,distribution_approximation=montecarlo,number_of_particles=100000,resampling=systematic,mh_replic=0,mode_compute=8,additional_optimizer_steps=[4,8],silent_optimizer);
@#endif

@#if ALGO_GUF
  estimation(order=2,nograph,filter_algorithm=gf,mh_replic=0,silent_optimizer,mode_compute=8);
  estimation(order=3,nograph,filter_algorithm=gf,mh_replic=0,mode_compute=8,silent_optimizer,mode_check);
@#endif

@#if ALGO_ONLINE
%  estimation(order=1,nograph,mode_compute=11,mh_replic=0,particle_filter_options=('liu_west_delta',0.9));
%  estimation(order=2,nograph,number_of_particles=10000,mode_compute=11,mh_replic=0,particle_filter_options=('liu_west_delta',0.9));
%  estimation(order=3,nograph,number_of_particles=10000,mode_compute=11,mh_replic=0,particle_filter_options=('liu_west_delta',0.9));
  estimation(order=1,posterior_sampling_method='online',posterior_sampler_options=('particles',1000));
  estimation(order=2,posterior_sampling_method='online',posterior_sampler_options=('particles',1000));
  estimation(order=3,posterior_sampling_method='online',filter_algorithm=nlkf,proposal_approximation=montecarlo,number_of_particles=100,posterior_sampler_options=('particles',100));
@#endif

@#if MCMC

  estimated_params_init;
      alp, 0.3980;
      bet, 0.9907;
      tet, 0.3565;
      tau, 48.054;
      delt, 0.0197;
      rho, 0.9511;
      stderr e_a, .0346;
      stderr y, 0.0004;
      stderr l, 0.0001;
      stderr i, 0.0002;
  end;
%  estimation(order=3,number_of_particles=10000,silent_optimizer,mode_compute=8,cova_compute=0,MCMC_jumping_covariance=prior_variance,mh_init_scale_factor=0.01);

  estimated_params_init;
      alp, 0.3982;
      bet, 0.9910;
      tet, 0.3566;
      tau, 46.883;
      delt, 0.0197;
      rho, 0.9511;
      stderr e_a, .0346;
      stderr y, 0.0004;
      stderr l, 0.0001;
      stderr i, 0.0001;
  end;
estimation(order=3,filter_algorithm=nlkf,number_of_particles=10000,proposal_approximation=montecarlo,resampling=none,silent_optimizer,mode_compute=0,cova_compute=0,MCMC_jumping_covariance=prior_variance,mh_init_scale_factor=0.01);
@#endif
