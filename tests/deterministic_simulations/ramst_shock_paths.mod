// Simple test of shock_paths with perfect_foresight_{setup,solver}

var c k;
varexo x;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

verbatim;
mydb = dseries([linspace(0, 1, 210); linspace(-2, 0, 210)]','2017Y',{'A1'; 'A2'},{'A_1'; 'A_2'});
mydb2 = table(linspace(0, 1, 210)', linspace(-2, 0, 210)', 'VariableNames', {'A1', 'A2' });
mydb3 = table(1.04, 2, 'VariableNames', {'Y1', 'Y2' }); % Special case of a single-row table
end;

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

database mydb mydb2 mydb3;

shock_paths;
  var x;
  periods 2018Y, 2019Y, 2020Y, 2021Y, 2022Y, 2023Y:end;
  values 1.2, mydb.A1(+1), mydb2.A2(-1), mydb3.Y1, self.x(-1) + 0.01, initval.x + log1p(0.01);
end;

perfect_foresight_setup(periods=200, first_simulation_period = 2018Y);
perfect_foresight_solver;

assert(oo_.exo_simul(3) == mydb.A1('2020Y').data)
assert(oo_.exo_simul(4) == mydb2.A2(2))
assert(oo_.exo_simul(5) == mydb3.Y1)
assert(oo_.exo_simul(6) == oo_.exo_simul(5) + 0.01)
assert(all(oo_.exo_simul(7:202) == 1 + log1p(0.01)))
