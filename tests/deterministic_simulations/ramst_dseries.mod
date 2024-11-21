// Test dates in shocks block and in perfect_foresight_setup

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

shocks;
  var x;
  periods 2000Q2:2000Q4;
  values 1.2;
end;

perfect_foresight_setup(first_simulation_period = 2000Q1, last_simulation_period = 2049Q4);
perfect_foresight_solver;

if ~isequal(find(oo_.exo_simul == 1.2), (3:5)')
  error('Incorrect handling of dates in the shocks block')
end

if ~isequal(Simulated_time_series.dates, dates('1999Q4'):dates('2050Q1'))
    error('Incorrect Simulated_time_series object')
end
