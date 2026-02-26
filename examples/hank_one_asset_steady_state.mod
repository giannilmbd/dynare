/*
 * One-Asset HANK Model (Heterogeneous Agent New Keynesian)
 * Example: Computing the Steady State with Parameter Calibration
 *
 * This example demonstrates:
 * 1. Computing the steady state numerically from an initial guess
 * 2. Calibrating a free parameter (beta) to match a market-clearing condition
 * 3. Solving a HANK model with one liquid asset
 * 4. Running stochastic simulations to generate time series
 */

/*
 * Copyright © 2026 Dynare Team
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

// Declare heterogeneity dimension
heterogeneity_dimension households;

// Household-level variables
var(heterogeneity=households)
   c      (long_name = 'consumption')
   n      (long_name = 'labor supply')
   ns     (long_name = 'effective labor supply')
   a      (long_name = 'assets')
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

// Parameters
parameters
   beta vphi
   eis frisch
   rho_e sig_e
   rho_Z sig_Z
   mu kappa phi
   Z B r_ss
;

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

// Household optimization problem
model(heterogeneity=households);
   [name='Euler equation with borrowing constraint']
   beta * (1 + r(+1)) * c(+1)^(-1/eis) - c^(-1/eis) = 0 ⟂ a >= 0;

   [name='Budget constraint']
   (1 + r) * a(-1) + w * n * e + (Div-Tax) * e - c - a;

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

// Aggregate shock specification
shocks;
    var G; stderr 0.01;
    var markup; stderr 0.01;
    var rstar; stderr 0.01;
end;

//==========================================================================
// STEP 1: Compute steady state with parameter calibration
//==========================================================================
// The initial guess is loaded from hank_one_asset_sp.mat. The discount factor
// beta is calibrated so that the asset market clearing condition holds.
// The .mat file contains: steady_state.free_parameters.beta.initial_guess
heterogeneity_compute_steady_state(filename = hank_one_asset_sp,
    calibration_target_equations=['Asset market clearing'],
    time_iteration_solver_stop_on_error);

//==========================================================================
// STEP 2: Solve the model
//==========================================================================
// Compute the linearized solution using sequence-space Jacobians
heterogeneity_solve(truncation_horizon = 300);

//==========================================================================
// STEP 3: Run stochastic simulation
//==========================================================================
// Simulate 1000 periods
heterogeneity_simulate(periods = 1000);
