/*
 * An elementary RBC model (same as ramst.mod), simulated in perfect foresight
 * with expectation errors: agents behave as under perfect foresight, but they
 * can still be surprised by unexpected shocks, and thus recompute their
 * optimal plans when such an unexpected shock happens.
 */

/*
 * Copyright © 2001-2025 Dynare Team
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


// Endogenous variables: consumption and capital
var c k;

// Exogenous variable: technology level
varexo x;

// Parameters declaration and calibration
parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

// Equilibrium conditions
model;
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1); // Resource constraint
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam); // Euler equation
end;

// Steady state (analytically solved)
initval;
x = 1;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
end;

// Check that this is indeed the steady state
steady;

// Declare a positive technological shock in period 1
shocks;
  var x;
  periods 1;
  values 1.2;
end;

// Declare a positive technological shock in period 3, that agents only learn about in period 2
shocks(learnt_in = 2);
  var x;
  periods 3;
  values 1.1;
end;

// Declare a permanent positive technological shock, that agents only learn about in period 5.
// Note that this will automatically trigger the recomputation of the terminal steady state.
endval(learnt_in = 5);
  x = 1.05;
end;

// Prepare the deterministic simulation of the model over 200 periods
perfect_foresight_with_expectation_errors_setup(periods=200);

// Perform the simulation
perfect_foresight_with_expectation_errors_solver;

// Display the path of consumption and capital
rplot c;
rplot k;
