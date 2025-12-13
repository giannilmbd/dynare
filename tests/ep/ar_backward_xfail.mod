// test purely backward model
var Efficiency, efficiency;

varexo EfficiencyInnovation;

parameters rho, effstar, sigma;

/*
** Calibration
*/


rho     =  0.950;
effstar =  1.000;
sigma   =  0.0001;

external_function(name=mean_preserving_spread,nargs=2);

model;

  // Eq. n°1:
  efficiency = rho*efficiency(-1) + sigma*EfficiencyInnovation;

  // Eq. n°2:
  Efficiency = effstar*exp(efficiency-mean_preserving_spread(rho,sigma));

end;

shocks;
var EfficiencyInnovation = 1;
end;

steady_state_model;
efficiency=0;
Efficiency=effstar;
end;

steady;

extended_path(order=0,periods=100);
ts=Simulated_time_series;

if ~oo_.extended_path.status
    error('Extended path did not find solution in ar.mod')
end

extended_path(order=1,periods=100);
sts=Simulated_time_series;

if ~oo_.extended_path.status
    error('Extended path did not find solution in ar.mod')
end

// The model is backward, we do not care about future uncertainty, extended path and stochastic extended path
// should return the same results.
if max(max(abs(ts.data-sts.data)))>pi*options_.dynatol.x
   error('Stochastic Extended Path:: Something is wrong here (potential bug in extended_path.m)!!!')
end

extended_path(order=1,periods=100,tree=sparse);
sts_sparse=Simulated_time_series;

if ~oo_.extended_path.status
    error('Extended path did not find solution in ar.mod')
end

// The model is backward, we do not care about future uncertainty, extended path and stochastic extended path
// should return the same results.
if max(max(abs(ts.data-sts_sparse.data)))>pi*options_.dynatol.x
   error('Stochastic Extended Path:: Something is wrong here (potential bug in extended_path.m)!!!')
end
