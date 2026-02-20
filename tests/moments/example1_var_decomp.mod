@#include "example1_common.inc"

// Parameter values
alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul(relative_irf,order=1,periods=0);
oo_1_theoretic=oo_;
stoch_simul(relative_irf,order=1,periods=100000);
oo_1_simul=oo_;
stoch_simul(relative_irf,order=2,periods=0);
oo_2_theoretic=oo_;
set_dynare_seed('default');
stoch_simul(relative_irf,order=2,periods=100000);
oo_2_simul=oo_;

if max(max(abs(oo_1_theoretic.variance_decomposition-oo_1_simul.variance_decomposition)))>2
   error('Variance Decomposition wrong')
end

if max(max(abs(oo_1_theoretic.variance_decomposition-oo_2_theoretic.variance_decomposition)))>1e-10
   error('Theoretical Variance Decomposition wrong')
end

if max(max(abs(oo_2_theoretic.variance_decomposition-oo_2_simul.variance_decomposition)))>3
   error('Variance Decomposition wrong')
end

