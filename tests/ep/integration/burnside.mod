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

tic
extended_path(periods=10, order=4, tree=perfect, integration=unscented);
dprintf('SEP(4) with perfect tree and unscented integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside.mod (tree=perfect, order=4, integration=unscented).')
end

seed = 31415;

tic
extended_path(periods=10, order=4, tree=perfect, number_of_quadrature_nodes=3, integration=tensor_gaussian_quadrature);
dprintf('SEP(4) with perfect tree and gaussian quadrature integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside.mod (tree=perfect, order=4, number_of_quadrature_nodes=3, integration=tensor_gaussian_quadrature).')
end

seed = 31415;

tic
extended_path(periods=10, order=4, tree=perfect, integration=stroud);
dprintf('SEP(4) with perfect tree and cubature integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside.mod (tree=perfect, order=4, integration=stroud).')
end

seed = 31415;

tic
extended_path(periods=10, order=4, tree=sparse, integration=unscented);
dprintf('SEP(4) with sparse tree and unscented integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside.mod (tree=sparse, order=4, integration=unscented).')
end

seed = 31415;

tic
extended_path(periods=10, order=4, tree=sparse, number_of_quadrature_nodes=3, integration=tensor_gaussian_quadrature);
dprintf('SEP(4) with sparse tree and gaussian quadrature integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside.mod (tree=sparse, order=4, number_of_quadrature_nodes=3, integration=tensor_gaussian_quadrature).')
end

seed = 31415;

tic
extended_path(periods=10, order=4, tree=sparse, integration=stroud);
dprintf('SEP(4) with sparse tree and cubature integration took %6.4f seconds.\n', toc);

if ~oo_.extended_path.status
    error('Stochastic extended path did not find solution in integration/burnside.mod (tree=sparse, order=4, integration=stroud).')
end
