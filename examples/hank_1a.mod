// One-Asset HANK Model (Heterogeneous Agent New Keynesian)
// Example: Stochastic Simulation
//
// This example demonstrates:
// 1. Loading a pre-computed steady state
// 2. Solving a HANK model with one liquid asset
// 3. Computing impulse response functions (IRFs)
// 4. Running stochastic simulations to generate time series

// Declare heterogeneity dimension
heterogeneity_dimension households;

// Household-level variables
var(heterogeneity=households)
   c      (long_name = 'consumption')
   n      (long_name = 'labor supply')
   ns     (long_name = 'effective labor supply')
   a      (long_name = 'assets')
   Va     (long_name = 'derivative of value function w.r.t. a')
;

// Household-level shocks
varexo(heterogeneity=households)
   e      (long_name = 'idiosyncratic efficiency')
;

// Aggregate variables
var
   Y    (long_name = 'Aggregate output')
   L    (long_name = 'Aggregate labor')
   w    (long_name = 'Real wage')
   pi   (long_name = 'Inflation')
   Div  (long_name = 'Dividends')
   Tax  (long_name = 'Taxes')
   r    (long_name = 'Real interest rate')
;

// Aggregate shocks
varexo G markup rstar;

// Aggregate shock specification
shocks;
    var G; stderr 0.01;
    var markup; stderr 0.01;
    var rstar; stderr 0.01;
end;

// Parameters
parameters
   beta vphi
   eis frisch
   rho_e sig_e
   rho_Z sig_Z
   mu kappa phi
   Z B r_ss
;

// Household optimization problem
model(heterogeneity=households);
   [name='Euler equation with borrowing constraint']
   beta * Va(+1) - c^(-1/eis) = 0 ⟂ a >= 0;

   [name='Budget constraint']
   (1 + r) * a(-1) + w * n * e + (Div-Tax) * e - c - a;

   [name='Envelope condition']
   Va = (1 + r) * c^(-1/eis);

   [name='Labor supply']
   vphi*n^(1/frisch) - w*e*c^(-1/eis);

   [name='Effective labor supply']
   ns = n * e;
end;

// Aggregate equilibrium conditions
model;
   [name='Labor demand']
   L - Y / Z;

   [name='Dividends']
   Div - (Y - w * L - mu / (mu - 1) / (2 * kappa) * log(1 + pi)^2 * Y);

   [name='Taylor rule']
   (1 + r_ss + rstar(-1) + phi * pi(-1)) / (1 + pi) - 1 - r;

   [name='Government budget constraint']
   Tax - (r * B) - G;

   [name='New Keynesian Phillips curve']
   kappa * (w / Z - 1 / mu)
   + Y(+1)/Y * log(1 + pi(+1)) / (1 + r(+1))
   + markup
   - log(1 + pi);

   [name='Asset market clearing']
   SUM(a) - B;

   [name='Labor market clearing']
   sum(ns) - L;
end;

B = 5.6;
Z = 1;
beta = 0.9822435537831447;
eis = 0.5;
frisch = 0.5;
kappa = 0.1;
mu = 1.2;
phi = 1.5;
r_ss = 0.005;
vphi = 0.7864334221640324;

//==========================================================================
// STEP 1: Load pre-computed steady state
//==========================================================================
heterogeneity_load_steady_state(filename = hank_1a);

//==========================================================================
// STEP 2: Solve the model
//==========================================================================
// Compute the linearized solution using sequence-space Jacobians
heterogeneity_solve(truncation_horizon = 300);

//==========================================================================
// STEP 3: Run stochastic simulation
//==========================================================================
// Simulate 1000 periods, drop first 100 as burn-in
heterogeneity_simulate(periods = 1000);

// Alternative examples:
// Compute only IRFs (no stochastic simulation):
// heterogeneity_simulate(irf = 80);
//
// Compute both IRFs and stochastic simulation together:
// heterogeneity_simulate(periods = 1000, irf = 80);
