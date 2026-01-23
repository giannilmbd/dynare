// Test against regression in ordering filtered autocorrelations (were transposed)
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

stoch_simul(order=1,nofunctions,irf=0,bandpass_filter) y c;
out1=oo_.autocorr{1};

stoch_simul(order=1,nofunctions,irf=0,periods=10000000,bandpass_filter) y c ;
out2=oo_.autocorr{1};
if max(max(abs(out1-out2)))>0.03
    error('The autocorrelations are wrong.')
end

