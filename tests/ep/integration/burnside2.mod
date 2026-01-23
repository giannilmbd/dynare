var y x z u;

varexo e f g;

parameters beta theta rho xbar;
xbar = 0.0179;
rho =  -0.139;
theta = -1.5;
beta = 0.95;

model;
  1 = beta*exp(theta*u(+1))*(1+y(+1))/y;
  x = (1-rho)*xbar + rho*x(-1) + e;
  z = (1-rho)*xbar + rho*z(-1) + .3*x(-1) + f;
  u = (x+z)/2 + g;
end;

shocks;
  var e; stderr 0.0348;
  var f; stderr 0.0100;
  var g; stderr 0.0100;
end;

initval;
  x = xbar;
  z = xbar;
  u = xbar;
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

tic
extended_path(periods=10, order=2, tree=perfect, integration=unscented);
dprintf('SEP(2) with perfect tree and unscented integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside2.mod (tree=perfect, order=2, integration=unscented).')
end

seed = 31415;

tic
extended_path(periods=10, order=2, tree=perfect, number_of_quadrature_nodes=3, integration=tensor_gaussian_quadrature);
dprintf('SEP(2) with perfect tree and gaussian quadrature integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside2.mod (tree=perfect, order=2, number_of_quadrature_nodes=3, integration=tensor_gaussian_quadrature).')
end

seed = 31415;

tic
extended_path(periods=10, order=2, tree=perfect, integration=stroud);
dprintf('SEP(2) with perfect tree and cubature integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside2.mod (tree=perfect, order=2, integration=stroud).')
end

seed = 31415;

tic
extended_path(periods=10, order=2, tree=sparse, integration=unscented);
dprintf('SEP(2) with sparse tree and unscented integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside2.mod (tree=sparse, order=2, integration=unscented).')
end

seed = 31415;

tic
extended_path(periods=10, order=2, tree=sparse, number_of_quadrature_nodes=3, integration=tensor_gaussian_quadrature);
dprintf('SEP(2) with sparse tree and gaussian quadrature integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside2.mod (tree=sparse, order=2, number_of_quadrature_nodes=3, integration=tensor_gaussian_quadrature).')
end

seed = 31415;

tic
extended_path(periods=10, order=2, tree=sparse, integration=stroud);
dprintf('SEP(2) with sparse tree and cubature integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside2.mod (tree=sparse, order=2, integration=stroud).')
end
