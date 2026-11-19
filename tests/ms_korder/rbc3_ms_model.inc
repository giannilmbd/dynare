// RBC with three states (capital k, AR(1) TFP z, AR(1) depreciation d), two correlated innovations, and switching
// TFP level and depreciation rate (aA, aD: [ms] variables, steady state 0).
var k c z d aA aD;
varexo e1 e2;
parameters ms_zero beta zeta gam Abar dbar rhoz rhod sigz sigd;
beta = 0.99; zeta = 0.36; gam = 2; Abar = 1; dbar = 0.025; rhoz = 0.9; rhod = 0.8; sigz = 0.01; sigd = 0.05;
ms_zero = 0;
model;
c + k = (Abar + aA)*exp(z)*k(-1)^zeta + (1 - dbar*exp(d) - aD)*k(-1);
c^(-gam) = beta*c(+1)^(-gam)*(1 - dbar*exp(d(+1)) - aD(+1) + zeta*(Abar + aA(+1))*exp(z(+1))*k^(zeta-1));
z = rhoz*z(-1) + sigz*e1;
d = rhod*d(-1) + sigd*e2;
[ms='aA'] aA = ms_zero;
[ms='aD'] aD = ms_zero;
end;
steady_state_model;
k = (zeta*Abar/(1/beta - 1 + dbar))^(1/(1 - zeta));
c = Abar*k^zeta - dbar*k;
z = 0; d = 0; aA = 0; aD = 0;
end;
shocks; var e1; stderr 1; var e2; stderr 1; corr e1, e2 = 0.3; end;
steady;
