// Simple test of the exogenize stanza (controlled paths) within shock_paths block

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
end;

initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

steady;

shock_paths;
  exogenize c;
  periods 2019Y, 4:5;
  values 1.6, 1.7;
  endogenize x;
end;

perfect_foresight_setup(periods=200, first_simulation_period = 2018Y);
perfect_foresight_solver;

assert(oo_.endo_simul(1,3) == 1.6)
assert(all(oo_.endo_simul(1,5:6) == 1.7))
assert(isequal(find(oo_.exo_simul ~= 1), [3; 5; 6]))
