var pi x R;

varexo epsilon;

parameters  beta kappa lambdax lambdaR sigma;

sigma = 1;
beta = 0.99;
kappa = 0.102;
lambdax = 0.1;

model(linear);
pi = beta*pi(+1)+kappa*x+epsilon;
x = x(+1)-(1/sigma)*(R-pi(+1));
end;

shocks;
var epsilon; stderr 0.01;
end;

planner_objective (pi^2 +lambdax*x^2);

discretionary_policy(instruments=(R),irf=20,planner_discount= beta);




