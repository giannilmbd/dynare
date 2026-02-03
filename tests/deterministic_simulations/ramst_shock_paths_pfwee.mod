// Simple test of shock_paths with perfect_foresight_with_expectation_errors_{setup,solver}

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

check;

shock_paths;
  var x;
  periods 2018Y, 2019Y, 2020Y:end;
  values 1.2, self.x(-1) + 0.01, initval.x + 0.02;
end;

shock_paths(learnt_in = 2);
  var x;
  periods 2019Y;
  values prev.x + 0.03;
end;

shock_paths(learnt_in = 2024Y);
  var x;
  periods 2026Y, 2027Y:end;
  values learnt_in(2019Y).x(+1) + 0.04, learnt_in(2018Y).x(-8);
end;

perfect_foresight_with_expectation_errors_setup(periods=200, first_simulation_period = 2018Y);
perfect_foresight_with_expectation_errors_solver;

assert(oo_.exo_simul(3) == 1.2 + 0.01 + 0.03)
assert(all(oo_.exo_simul(4:9) == 1 + 0.02))
assert(oo_.exo_simul(10) == 1 + 0.02 + 0.04)
assert(oo_.exo_simul(11) == 1.2 + 0.01)
assert(all(oo_.exo_simul(12:202) == 1 + 0.02))
