// Tests the analytic_derivation option

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

alp = 0.33;
bet = 0.99;
gam = 0.003;
mst = 1.011;
rho = 0.7;
psi = 0.787;
del = 0.02;

model;
dA = exp(gam+e_a);
log(m) = (1-rho)*log(mst) + rho*log(m(-1))+e_m;
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
W = l/n;
-(psi/(1-psi))*(c*P/(1-n))+l/n = 0;
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
P*c = m;
m-1+d = l;
e = exp(e_a);
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
gy_obs = dA*y/y(-1);
gp_obs = (P/P(-1))*m(-1)/dA;
end;

steady_state_model;
  dA = exp(gam);
  gst = 1/dA;
  m = mst;
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/mst )^(-1);
  nust = psi*mst^2/( (1-alp)*(1-psi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  k  = khst*n;

  l  = psi*mst*n/( (1-psi)*(1-n) );
  c  = mst/P;
  d  = l - mst + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  R  = mst/bet;
  W  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;
  
  gp_obs = m/dA;
  gy_obs = dA;
end;

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
var gy_obs= 0.005^2;
var gp_obs= 0.005^2;
corr gp_obs, gy_obs= 0.05;
end;

stoch_simul(order=1,periods=1000, irf=0,nomoments,noprint);

send_endogenous_variables_to_workspace;
% stoch_simul does not add measurement error; draw it from M_.H
% varobs order: gp_obs (1), gy_obs (2)
verbatim;
T = length(gy_obs);
me = chol(M_.H, 'lower') * randn(size(M_.H, 1), T);
gp_obs = gp_obs + me(1,:)';
gy_obs = gy_obs + me(2,:)';
end;
save('my_data.mat','gp_obs','gy_obs');


estimated_params;
alp, 0.33;
stderr e_a, 0.02, 0, inf;
stderr e_m, 0.02, 0, inf;
corr e_a, e_m, 0 ,-1, 1;
stderr gy_obs, 0.005, 0, inf;
stderr gp_obs, 0.005, 0, inf;
corr gp_obs, gy_obs, 0.05, -1, 1;
end;

varobs gp_obs gy_obs;

options_.solve_tolf = 1e-12;

estimation(order=1,mode_compute=5,analytic_derivation,kalman_algo=2,datafile=my_data,nobs=192,mh_replic=0,mh_nblocks=2,mh_jscale=0.8,plot_priors=0);
