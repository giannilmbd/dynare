var k A c l i y;
varexo e_a;

parameters alp bet tet tau delt rho ;
alp = 0.4;
bet = 0.99;
tet = 0.357 ;
tau =  50 ;
delt = 0.02;
rho = 0.95;

model;
	c = ((1 - alp)*tet/(1-tet))*A*(1-l)*((k(-1)/l)^alp) ;
	y = A*(k(-1)^alp)*(l^(1-alp)) ;
	i = y-c ;
	k = (1-delt)*k(-1) + i ;
	log(A) = rho*log(A(-1)) + e_a ;
	(((c^(tet))*((1-l)^(1-tet)))^(1-tau))/c - bet*((((c(+1)^(tet))*((1-l(+1))^(1-tet)))^(1-tau))/c(+1))*(1 -delt+alp*(A(1)*(k^alp)*(l(1)^(1-alp)))/k)=0 ;
end;

steady_state_model;
	k = -(alp-1)*(alp^(1/(1-alp)))*(bet^(1/(1-alp)))*((bet*(delt-1)+1)^(alp/(alp-1)))*tet/(-alp*delt*bet+delt*bet+alp*tet*bet-bet-alp*tet+1);
	l = (alp-1)*(bet*(delt-1)+1)*tet/(alp*tet+bet*((alp-1)*delt-alp*tet+1)-1) ;
	y = (k^alp)*(l^(1-alp)) ;
	i = delt*k ;
	c = y - i ;
	A = 1;
end;

shocks;
    var e_a; stderr 0.035;
end;

steady;

estimated_params;
    alp, uniform_pdf,,, 0.0001, 0.99;
    bet, uniform_pdf,,, 0.0001, 0.99999;
    tet, uniform_pdf,,, 0.0001, .999;
    tau, uniform_pdf,,, 0.0001, 100;
    delt, uniform_pdf,,, 0.0001, 0.05;
    rho, uniform_pdf,,, 0.0001, 0.9999;
    stderr e_a, uniform_pdf,,, 0.00001, 0.1;
    stderr y, uniform_pdf,,, 0.00001, 0.1;
    stderr l, uniform_pdf,,, 0.00001, 0.1;
    stderr i, uniform_pdf,,, 0.00001, 0.1;
end;

estimated_params_init;
	alp, 0.4;
	bet, 0.99;
	tet, 0.357;
	tau, 50;
	delt, 0.02;
	rho, 0.95;
	stderr e_a, .035;
	stderr y, .001;
	stderr l, .001;
	stderr i, .001;
end;

varobs y l i ;

%data(file='./risky.m');
%stoch_simul(order=5,periods=1000);
%datatomfile('mysample')
%return;

data(file='./mysample.m',first_obs=801Y,nobs=50); %no measurement errors added in the simulated data