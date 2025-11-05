%test whether forecast for AR(1) is correct
// Endogenous variables: consumption and capital
var junk A junk2;

// Exogenous variable: technology level
varexo epsilon;

// Parameters declaration and calibration
parameters rho;

rho=0.9;

// Equilibrium conditions
model;
  A= rho*A(-1)+epsilon;
  junk=junk(+1)^2;
  junk2=0.9*junk2(-1);
end;

steady_state_model;
  A = 0;
  junk=0;
  junk2=0;
end;
steady;

// Declare a positive technological shock in period 1
shocks;
  var epsilon=0.01^2;
end;

stoch_simul(order=2,irf=0);

options_.forecast_replic=50000;
forecast(periods=100,conf_sig=0.95) A;
%% get theoretical forecast error variance for AR(1) process:
var_analytical=NaN(1+options_.forecast,1);
var_analytical(1)=0; %zero variance in initial period
for ii=1:100
    var_analytical(ii+1,1)=var_analytical(ii,1)+(rho^(ii-1))^2*0.01^2;
end
std_analytical=sqrt(var_analytical);

if norm(oo_.forecast.Mean.A-0,Inf)>1e-9
    error('Mean forecast does not match theoretical moments')
end
if norm(oo_.forecast.HPDsup.A-norminv(0.975)*std_analytical(2:end),Inf)>3e-3
    error('Forecast variance does not match theoretical moments')
end

varobs A junk2;

shocks;
  var epsilon=0.01^2;
  var A=0.01^2;
  var junk2=0.1^2;
end;

forecast(periods=100,conf_sig=0.95) A;
std_analytical_ME=sqrt(var_analytical+0.01^2);
if norm(oo_.forecast.Mean.A-0,Inf)>1e-9
    error('Mean forecast with ME does not match theoretical moments')
end
if norm(oo_.forecast.HPDsup.A-norminv(0.975)*std_analytical(2:end),Inf)>1e-3
    error('Forecast variance with ME does not match theoretical moments')
end
if norm(oo_.forecast.HPDsup_ME.A-norminv(0.975)*std_analytical_ME(2:end),Inf)>1e-3
    error('Forecast variance with ME does not match theoretical moments')
end
