/*
 * Krusell and Smith (1998)
 * Example: Computing the Steady State
 *
 * This example demonstrates:
 * 1. Computing the steady state numerically from an initial guess
 * 2. Solving the model
 * 3. Computing impulse response functions to aggregate shocks
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
    c  (long_name = 'Consumption')
    a  (long_name = 'Assets')
;

// Household-level shocks
varexo(heterogeneity=households)
    e  (long_name = 'Idiosyncratic productivity shock')
;

// Aggregate variables
var
    Y (long_name = 'Aggregate output')
    r (long_name = 'Rate of return on capital net of depreciation')
    w (long_name = 'Wage rate')
    K (long_name = 'Aggregate capital')
;

// Aggregate shock
varexo Z (long_name = 'Aggregate productivity shock');

// Parameters
parameters
    L     (long_name = 'Labor')
    alpha (long_name = 'Share of capital in production function')
    beta  (long_name = 'Subjective discount rate of households')
    delta (long_name = 'Capital depreciation rate')
    eis   (long_name = 'Elasticity of intertemporal substitution')
    Z_ss  (long_name = 'Aggregate TFP shock average value')
;

L = 1;
Z_ss = 0.8816460975214567;
alpha = 0.11;
beta = 0.9819527880123727;
delta = 0.025;
eis = 1;

// Household optimization problem
model(heterogeneity=households);
    [name='Euler equation with borrowing constraint']
    beta*(1+r(+1))*c(+1)^(-1/eis)-c^(-1/eis)=0 ⟂ a>=0;

    [name='Budget constraint']
    (1+r)*a(-1)+w*e-c-a;
end;

// Aggregate equilibrium conditions
model;
    [name='Production function']
    Y = (Z_ss+Z) * K(-1)^alpha * L^(1 - alpha);

    [name='Capital rental rate']
    r = alpha * (Z_ss+Z) * (K(-1) / L)^(alpha - 1) - delta;

    [name='Wage rate']
    w = (1 - alpha) * (Z_ss+Z) * (K(-1) / L)^alpha;

    [name='Capital market clearing']
    K - SUM(a);
end;

// Aggregate shock specification
shocks;
    var Z; stderr 0.01;
end;

//==========================================================================
// STEP 1: Compute steady state from initial guess
//==========================================================================
// The initial guess (policy functions, grids, shock discretization)
// is loaded from krusell_smith_1998_sp.mat. No free parameters: the command
// iterates on policy functions and the distribution, then reports aggregate
// residuals.
heterogeneity_compute_steady_state(filename = krusell_smith_1998_sp);

//==========================================================================
// STEP 2: Solve the model
//==========================================================================
// Compute the linearized solution using sequence-space Jacobians
// truncation_horizon = 400 specifies the time horizon for fake news algorithm
heterogeneity_solve(truncation_horizon = 400);

//==========================================================================
// STEP 3: Compute impulse response functions
//==========================================================================
// Compute IRFs with 80 periods horizon
heterogeneity_simulate(irf = 80);
