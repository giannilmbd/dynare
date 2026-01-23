// --+ options: stochastic

var y pie r;
varexo e_y e_pie;

parameters delta sigma alpha kappa gamma1 gamma2;

delta =  0.44;
kappa =  0.18;
alpha =  0.48;
sigma = -0.06;

gamma1 = 1.5;
gamma2 = 0.5;

model(use_dll);
y  = delta * y(-1)  + (1-delta)*y(+1)+sigma *(r - pie(+1)) + e_y;
pie  =   alpha * pie(-1) + (1-alpha) * pie(+1) + kappa*y + e_pie;
r = gamma1*pie+gamma2*y;
end;

shocks;
var e_y;
stderr 0.63;
var e_pie;
stderr 0.4;
end;

steady;


// Extended path simulation
options_.ep.order = 0;
[ts, oo_] = extended_path([], 100, [], options_, M_, oo_);
if ~oo_.extended_path.status
   error('Extended path algorithm failed in ep/linearmodel1.mod')
end

// Stochastic extended path simulation (should be identical to extended path because the model is linear)
options_.ep.order = 1;
[sts, oo_] = extended_path([], 100, [], options_, M_, oo_);
if ~oo_.extended_path.status
   error('Stochastic extended path simulation failed in ep/linearmodel1.mod')
end

// The generated paths should be identical (because the model is linear)
if max(max(abs(ts.data-sts.data))) > 1e-12
   error('Stochastic extended path simulation are wrong in ep/linearmodel1.mod')
end

// Stochastic extended path simulation (should be identical to extended path because the model is linear, hybrid correction )
options_.ep.hybrid = 2;
[sts, oo_] = extended_path([], 100, [], options_, M_, oo_);
if ~oo_.extended_path.status
   error('Stochastic extended path (hybrid) simulation failed in ep/linearmodel1.mod')
end

// The generated paths should be identical (because the model is linear)
if max(max(abs(ts.data-sts.data))) > 1e-12
   error('Stochastic extended path (hybrid) simulation are wrong in ep/linearmodel1.mod')
end
