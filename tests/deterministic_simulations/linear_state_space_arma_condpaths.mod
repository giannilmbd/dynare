% Tests perfect simulation of linear model (sim1_linear.m) with controlled paths
%
% In particular, the endo and exo flipped (z vs u) have a nonzero steady state,
% to verify that this case is correctly handled.

var x
    y
    z;

varexo u
       v;

parameters a1 a2 a3 a4
	   b1 b2 b3
	   c1;

a1 =  .50;
a2 =  .00;
a3 =  .70;
a4 =  .40;
b1 =  .90;
b2 =  .00;
b3 =  .80;
c1 =  .95;


model(linear);
  y = a1*x(-1) + a2*x(+1) + a3*z + a4*y(-1);
  z = 0.5 + b1*z(-1) + b2*z(+1) + b3*x + u;
  x = c1*x(-1) + v +v(-1)+v(+1);
end;

initval;
  y=-1;
  x=-1;
  z=-1;
end;

endval;
  y=0;
  x=0;
  z=0;
  u=0.1;
end;

steady;

perfect_foresight_controlled_paths;
  exogenize z;
  periods 50;
  values 0.5;
  endogenize u;
end;

perfect_foresight_setup(periods=100);
perfect_foresight_solver;

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if ~(oo_.endo_simul(3, 51) == 0.5 && isequal(find(oo_.exo_simul(:,1) ~= 0.1), [1; 51]))
   error('Linear perfect foresight simulation with controlled paths failed')
end
