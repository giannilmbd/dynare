// (Stochastic) Extended path cannot simulate purely forward models.

var a;
varexo epsil ;
parameters betta;
betta = 0.97;

model;
a = betta*a(1)+epsil;
end;

initval;
a=0;
end;
steady;

shocks;
var epsil; stderr .2;
end;

steady;

extended_path(order=0, periods=100, solver_periods=200);
