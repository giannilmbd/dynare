@#include "example1_common.inc"

// Parameter values
alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

steady(solve_algo=4,maxit=1000);

stoch_simul(order=1,nofunctions,one_sided_hp_filter=1600,irf=0,periods=5000,filtered_theoretical_moments_grid=8192);

varobs k c y;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
var c; stderr 0.01;
end;

stoch_simul(order=1,nofunctions,one_sided_hp_filter=1600,irf=0,periods=5000,filtered_theoretical_moments_grid=8192);
