// ALGORITHM
@#ifndef LINEAR_KALMAN
    @#define LINEAR_KALMAN = 0
@#endif

@#ifndef NONLINEAR_KALMAN
    @#define NONLINEAR_KALMAN = 0
@#endif

@#ifndef ALGO_SIR
    @#define ALGO_SIR = 0
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
    @#define ALGO_ONLINE = 1
@#endif

@#ifndef MCMC
    @#define MCMC = 0
@#endif

@#ifndef SMC
    @#define SMC = 0
@#endif

%%

var k A c l i y;
varexo e_a;

parameters alp bet tet tau delt rho ;
alp = 0.4;
bet = 0.99;
tet = 0.357 ;
tau =  50 ;
delt = 0.02;
rho = 0.95;

model;
	c = ((1 - alp)*tet/(1-tet))*A*(1-l)*((k(-1)/l)^alp) ;
	y = A*(k(-1)^alp)*(l^(1-alp)) ;
	i = y-c ;
	k = (1-delt)*k(-1) + i ;
	log(A) = rho*log(A(-1)) + e_a ;
	(((c^(tet))*((1-l)^(1-tet)))^(1-tau))/c - bet*((((c(+1)^(tet))*((1-l(+1))^(1-tet)))^(1-tau))/c(+1))*(1 -delt+alp*(A(1)*(k^alp)*(l(1)^(1-alp)))/k)=0 ;
end;

steady_state_model;
	k = -(alp-1)*(alp^(1/(1-alp)))*(bet^(1/(1-alp)))*((bet*(delt-1)+1)^(alp/(alp-1)))*tet/(-alp*delt*bet+delt*bet+alp*tet*bet-bet-alp*tet+1);
	l = (alp-1)*(bet*(delt-1)+1)*tet/(alp*tet+bet*((alp-1)*delt-alp*tet+1)-1) ;
	y = (k^alp)*(l^(1-alp)) ;
	i = delt*k ;
	c = y - i ;
	A = 1;
end;

shocks;
    var e_a; stderr 0.035;
end;

steady;

estimated_params;
    alp, uniform_pdf,,, 0.0001, 0.99;
    bet, uniform_pdf,,, 0.0001, 0.99999;
    tet, uniform_pdf,,, 0.0001, .999;
    tau, uniform_pdf,,, 0.0001, 100;
    delt, uniform_pdf,,, 0.0001, 0.05;
    rho, uniform_pdf,,, 0.0001, 0.9999;
    stderr e_a, uniform_pdf,,, 0.00001, 0.1;
    stderr y, uniform_pdf,,, 0.00001, 0.1;
    stderr l, uniform_pdf,,, 0.00001, 0.1;
    stderr i, uniform_pdf,,, 0.00001, 0.1;
end;

estimated_params_init;
	alp, 0.4;
	bet, 0.99;
	tet, 0.357;
	tau, 50;
	delt, 0.02;
	rho, 0.95;
	stderr e_a, .035;
	stderr y, .001;
	stderr l, .001;
	stderr i, .001;
end;

varobs y l i ;

%data(file='./risky.m');
%stoch_simul(order=5,periods=1000);
%datatomfile('mysample')
%return;

data(file='./mysample.m',first_obs=801Y,nobs=200); %no measurement errors added in the simulated data

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
	estimation(order=2,nograph,number_of_particles=10000,mh_replic=0,mode_compute=8,cova_compute=0);
	estimation(order=2,nograph,number_of_particles=10000,mh_replic=0,silent_optimizer,mode_compute=8,additional_optimizer_steps=[8 8],cova_compute=0);
	estimation(order=3,nograph,number_of_particles=10000,mh_replic=0,silent_optimizer,mode_compute=8,additional_optimizer_steps=[8 8],cova_compute=0);
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
  estimation(order=3,posterior_sampling_method='online',filter_algorithm=nlkf,proposal_approximation=montecarlo,number_of_particles=500,posterior_sampler_options=('particles',500));
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

@#if SMC
  estimation(order=1,posterior_sampling_method='hssmc',posterior_sampler_options=('particles',1000));
  estimation(order=2,posterior_sampling_method='hssmc',posterior_sampler_options=('particles',1000));
  estimation(order=3,posterior_sampling_method='hssmc',filter_algorithm=nlkf,proposal_approximation=montecarlo,number_of_particles=500,posterior_sampler_options=('particles',500));
@#endif
