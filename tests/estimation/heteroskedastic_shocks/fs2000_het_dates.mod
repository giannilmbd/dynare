// Test heteroskedastic filter/smoother with dates syntax

@#include "fs2000_het_model.inc"

shocks;
var e_a; stderr 0.014;
var e_m; stderr 0.005;
end;

steady;

heteroskedastic_shocks;
  var e_b;
  periods 2000Y:2020Y;
  values 0.01;

  var e_a;
  periods 2000Y:2020Y;
  scales 0;
end;

ts = dseries('../fsdat_simul.m', 1901Y);

data(series = ts, nobs = 192);

estimation(order=1,silent_optimizer,mode_compute=5,loglinear,mh_replic=0,smoother,filtered_vars,forecast=8,filter_step_ahead=[1:8],consider_all_endogenous,heteroskedastic_filter);
