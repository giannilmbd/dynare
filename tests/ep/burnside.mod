var y x;

varexo e;

parameters beta theta rho xbar;
xbar = 0.0179;
rho =  -0.139;
theta = -1.5;
beta = 0.95;

model;
1 = beta*exp(theta*x(+1))*(1+y(+1))/y;
x = (1-rho)*xbar + rho*x(-1)+e;
end;

shocks;
var e; stderr 0.0348;
end;

initval;
x = xbar;
y = beta*exp(theta*xbar)/(1-beta*exp(theta*xbar));
end;

resid;

steady;

check;

if beta*exp(theta*xbar+.5*theta^2*M_.Sigma_e/(1-rho)^2)>1-eps
   disp('The model doesn''t have a solution!')
   return
end

seed = 31415;

set_dynare_seed(seed);
stoch_simul(order=1,irf=0,periods=200);
e_1 = oo_.exo_simul;

set_dynare_seed(seed);
stoch_simul(order=2,irf=0,periods=200);
e_2 = oo_.exo_simul;

set_dynare_seed(seed);
stoch_simul(order=4,irf=0,periods=200);
e_4 = oo_.exo_simul;


T = 100;

options_.ep.stochastic.algo=1; // Default is to use a sparse tree

tic

options_.ep.stochastic.order = 0;
ts0 = extended_path([], T, e_1, options_, M_, oo_);


options_.ep.stochastic.order = 1;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
ts1 = extended_path([], T, e_1, options_, M_, oo_);

options_.ep.stochastic.order = 1;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = 2;
ts1h = extended_path([], T, e_1, options_, M_, oo_);
options_.ep.stochastic.hybrid_order = 0;


options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
ts2 = extended_path([], T, e_1, options_, M_, oo_);


options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.algo=0; // Full tree of future innovations
ts2__ = extended_path([], T, e_1, options_, M_, oo_);
options_.ep.stochastic.algo=1;


options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = 2;
ts2h = extended_path([], T, e_1, options_, M_, oo_);

options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = 4;
ts2h = extended_path([], T, e_1, options_, M_, oo_);

options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = 2;
options_.ep.stochastic.algo = 0;
ts2h__ = extended_path([], T, e_1, options_, M_, oo_);
options_.ep.stochastic.algo = 1;

toc
