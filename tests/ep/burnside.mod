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
[ts0, o0, errorflag] = extended_path([], T, e_1, options_, M_, oo_);

if errorflag, error('EP failed'), end

options_.ep.stochastic.order = 1;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
[ts1, o1, errorflag] = extended_path([], T, e_1, options_, M_, oo_);

if errorflag, error('SEP(1) failed'), end

options_.ep.stochastic.order = 1;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = 2;
[ts1h, o1h, errorflag] = extended_path([], T, e_1, options_, M_, oo_);
options_.ep.stochastic.hybrid_order = 0;

if errorflag, error('Hybrid SEP(1) failed'), end

options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
[ts2, o2, errorflag] = extended_path([], T, e_1, options_, M_, oo_);

if errorflag, error('SEP(2) failed'), end

options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.algo=0; // Full tree of future innovations
[ts2__, o2full, errorflag] = extended_path([], T, e_1, options_, M_, oo_);
options_.ep.stochastic.algo=1;

if errorflag, error('SEP(2) with perfect tree failed'), end

options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = 2;
[ts2h, o2h, errorflag] = extended_path([], T, e_1, options_, M_, oo_);

if errorflag, error('Hybrid (order 2) SEP(2) failed'), end

options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = 4;
[ts2hh, o2hh, errorflag] = extended_path([], T, e_1, options_, M_, oo_);

if errorflag, error('Hybrid (order 4) SEP(2) failed'), end

options_.ep.stochastic.order = 2;
options_.ep.stochastic.IntegrationAlgorithm='Tensor-Gaussian-Quadrature';
options_.ep.stochastic.quadrature.nodes = 3;
options_.ep.stochastic.hybrid_order = 2;
options_.ep.stochastic.algo = 0;
[ts2h__, o2hfull, errorflag] = extended_path([], T, e_1, options_, M_, oo_);
options_.ep.stochastic.algo = 1;

if errorflag, error('Hybrid (order 2) SEP(2) with perfect tree failed'), end

toc

testresults = false;

if testresults

   saveresults = false;

   if saveresults
   save('burnside-simulations.mat', 'o0', 'o1', 'o1h', 'o2', 'o2full', 'o2h', 'o2hh', 'o2hfull')
   end

   if isfile('burnside-simulations.mat')
       simulations = load('burnside-simulations.mat');
       errorflag = false;
       if max(abs(o0.endo_simul(:)-simulations.o0.endo_simul(:)))>1e-16
           errorflag = true;
           dprintf('EP simulations are wrong.')
       end
       if max(abs(o1.endo_simul(:)-simulations.o1.endo_simul(:)))>1e-16
           errorflag = true;
           dprintf('SEP(1) simulations are wrong.')
       end
       if max(abs(o1h.endo_simul(:)-simulations.o1h.endo_simul(:)))>1e-16
           errorflag = true;
           dprintf('Hybrid SEP(1) simulations are wrong.')
       end
       if max(abs(o2.endo_simul(:)-simulations.o2.endo_simul(:)))>1e-16
           errorflag = true;
           dprintf('SEP(2) simulations are wrong.')
       end
       if max(abs(o2full.endo_simul(:)-simulations.o2full.endo_simul(:)))>1e-16
           errorflag = true;
           dprintf('SEP(2) with perfect tree simulations are wrong.')
       end
       if max(abs(o2h.endo_simul(:)-simulations.o2h.endo_simul(:)))>1e-16
           errorflag = true;
           dprintf('Hybrid SEP(2) simulations are wrong.')
       end
       if max(abs(o2hh.endo_simul(:)-simulations.o2hh.endo_simul(:)))>1e-16
           errorflag = true;
           dprintf('Hybrid (fourth order) SEP(2) simulations are wrong.')
       end
       if max(abs(o2hfull.endo_simul(:)-simulations.o2hfull.endo_simul(:)))>1e-16
           errorflag = true;
           dprintf('Hybrid SEP(2) with perfect tree simulations are wrong.')
       end
       if errorflag
           error('Some simulations do not match the expected results.')
       end
    end
end
