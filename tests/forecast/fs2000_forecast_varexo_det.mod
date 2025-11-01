/*
 * Copyright © 2004-2025 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

//checks correctness of forecasts of AR1 within larger model such that variable ordering can be tested

var m       ${m}$           (long_name='money growth')
    P       ${P}$           (long_name='Price level')
    c       ${c}$           (long_name='consumption')
    e       ${e}$           (long_name='capital stock')
    W       ${W}$           (long_name='Wage rate')
    R       ${R}$           (long_name='interest rate')
    k       ${k}$           (long_name='capital stock')
    d       ${d}$           (long_name='dividends')
    n       ${n}$           (long_name='labor')
    l       ${l}$           (long_name='loans')
    A
    gy_obs  ${\Delta \ln GDP}$  (long_name='detrended capital stock')
    gp_obs  ${\Delta \ln P}$    (long_name='detrended capital stock')
    y       ${y}$           (long_name='detrended output')
    dA      ${\Delta A}$    (long_name='TFP growth')
    ;
varexo e_a  ${\epsilon_A}$      (long_name='TFP shock')
    e_m     ${\epsilon_M}$      (long_name='Money growth shock')
    ;
varexo_det eps_det;

parameters alp  ${\alpha}$       (long_name='capital share')
    bet         ${\beta}$        (long_name='discount factor')
    gam         ${\gamma}$       (long_name='long-run TFP growth')
    logmst         ${\log(m^*)}$ (long_name='long-run money growth')
    rho         ${\rho}$         (long_name='autocorrelation money growth')
    phi         ${\phi}$         (long_name='labor weight in consumption')
    del         ${\delta}$       (long_name='depreciation rate')
    ;

% roughly picked values to allow simulating the model before estimation
alp = 0.33;
bet = 0.99;
gam = 0.003;
logmst = log(1.011);
rho = 0.7;
phi = 0.787;
del = 0.02;

model;
[name='NC before eq. (1), TFP growth equation']
dA = exp(gam+e_a);
[name='NC eq. (2), money growth rate']
log(m) = (1-rho)*logmst + rho*log(m(-1))+e_m;
[name='NC eq. (A1), Euler equation']
-P/(c(+1)*P(+1)*m)+bet*P(+1)*(alp*exp(-alp*(gam+log(e(+1))))*k^(alp-1)*n(+1)^(1-alp)+(1-del)*exp(-(gam+log(e(+1)))))/(c(+2)*P(+2)*m(+1))=0;
[name='NC below eq. (A1), firm borrowing constraint']
W = l/n;
[name='NC eq. (A2), intratemporal labour market condition']
-(phi/(1-phi))*(c*P/(1-n))+l/n = 0;
[name='NC below eq. (A2), credit market clearing']
R = P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(-alp)/W;
[name='NC eq. (A3), credit market optimality']
1/(c*P)-bet*P*(1-alp)*exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)/(m*l*c(+1)*P(+1)) = 0;
[name='NC eq. (18), aggregate resource constraint']
c+k = exp(-alp*(gam+e_a))*k(-1)^alp*n^(1-alp)+(1-del)*exp(-(gam+e_a))*k(-1);
[name='NC eq. (19), money market condition']
P*c = m;
[name='NC eq. (20), credit market equilibrium condition']
m-1+d = l;
[name='Definition TFP shock']
e = exp(e_a);
[name='Implied by NC eq. (18), production function']
y = k(-1)^alp*n^(1-alp)*exp(-alp*(gam+e_a));
[name='Observation equation GDP growth']
gy_obs = dA*y/y(-1);
[name='Observation equation price level']
gp_obs = (P/P(-1))*m(-1)/dA;
A= rho*A(-1)+e_a + eps_det;
end;

shocks;
var e_a; stderr 0.01;
var e_m; stderr 0.005;
end;

steady_state_model;
  dA = exp(gam);
  gst = 1/dA;
  m = exp(logmst);
  khst = ( (1-gst*bet*(1-del)) / (alp*gst^alp*bet) )^(1/(alp-1));
  xist = ( ((khst*gst)^alp - (1-gst*(1-del))*khst)/m )^(-1);
  nust = phi*m^2/( (1-alp)*(1-phi)*bet*gst^alp*khst^alp );
  n  = xist/(nust+xist);
  P  = xist + nust;
  k  = khst*n;

  l  = phi*m*n/( (1-phi)*(1-n) );
  c  = m/P;
  d  = l - m + 1;
  y  = k^alp*n^(1-alp)*gst^alp;
  R  = m/bet;
  W  = l/n;
  ist  = y-c;
  q  = 1 - d;

  e = 1;

  gp_obs = m/dA;
  gy_obs = dA;
end;


stoch_simul(order=1,irf=0);

forecast(periods=500,conf_sig=0.95) A;
%% get theoretical forecast error variance for AR(1) process:
var_analytical=NaN(1+500,1);
var_analytical(1)=0; %zero variance in initial period
for ii=1:500
    var_analytical(ii+1)=var_analytical(ii)+(rho^(ii-1))^2*0.01^2;
end
std_analytical=sqrt(var_analytical);

if norm(oo_.forecast.Mean.A-0,Inf)>1e-9
    error('Mean forecast does not match theoretical moments')
end
if norm(oo_.forecast.HPDsup.A-norminv(0.975)*std_analytical(2:end),Inf)>1e-9
    error('Forecast variance does not match theoretical moments')
end

varobs A;

shocks;
  var e_a=0.01^2;
  var A=0.01^2;
end;

forecast(periods=500,conf_sig=0.95) A;
std_analytical_ME=sqrt(var_analytical+0.01^2);
if norm(oo_.forecast.Mean.A-0,Inf)>1e-9
    error('Mean forecast with ME does not match theoretical moments')
end
if norm(oo_.forecast.HPDsup.A-norminv(0.975)*std_analytical(2:end),Inf)>1e-9
    error('Forecast variance with ME does not match theoretical moments')
end
if norm(oo_.forecast.HPDsup_ME.A-norminv(0.975)*std_analytical_ME(2:end),Inf)>1e-9
    error('Forecast variance with ME does not match theoretical moments')
end


forecast(periods=500,conf_sig=0.95);
std_analytical_ME=sqrt(var_analytical+0.01^2);
if norm(oo_.forecast.Mean.A-0,Inf)>1e-9
    error('Mean forecast with ME does not match theoretical moments')
end
if norm(oo_.forecast.HPDsup.A-norminv(0.975)*std_analytical(2:end),Inf)>1e-9
    error('Forecast variance with ME does not match theoretical moments')
end
if norm(oo_.forecast.HPDsup_ME.A-norminv(0.975)*std_analytical_ME(2:end),Inf)>1e-9
    error('Forecast variance with ME does not match theoretical moments')
end
