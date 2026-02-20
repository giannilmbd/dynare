// Test that histval_file works with a file created by smoother2histval (after a Metropolis)

var m P c e W R k d n l gy_obs gp_obs y dA;
varexo e_a e_m;

parameters alp bet gam mst rho psi del;

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
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a(-1)));
gy_obs = dA*y/y(-2);
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

% Use estimated parameters (posterior mean) as calibration
results_estimation=load(['fs2000_smooth' filesep 'Output' filesep 'fs2000_smooth_results']);
M_.params=results_estimation.M_.params;

steady;

histval_file(datafile = 'fs2000_histval.mat');

perfect_foresight_setup(periods = 100);
perfect_foresight_solver;


hf = dseries('fs2000_histval.mat');

% Check that y(0) is correct
assert(oo_.endo_simul(13,1) == hf.y('2Y').data)

% Check that y(-1) is correct
% First find the corresponding auxiliary variable
% NB: we must use cellfun and brackets become some orig_index and orig_lead_lag fields are empty
auxid = find([M_.aux_vars.type] == 1 & cellfun(@(x) isequal(x, 13), {M_.aux_vars.orig_index}) & cellfun(@(x) isequal(x, -1), {M_.aux_vars.orig_lead_lag}));
endoid = M_.aux_vars(auxid).endo_index;
assert(oo_.endo_simul(endoid,1) == hf.y('1Y').data);
