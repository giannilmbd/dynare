// Performs a pf_with_expectation_errors simulation with a perfect_foresight_controlled_paths block
// In particular, tests that it’s possible to change the expected value of a given variable (c and y in period 4)

@#define shocks = true
@#include "ramst_pf_controlled_paths_common.inc"

perfect_foresight_controlled_paths(learnt_in=2004Y);
  exogenize c;
  periods 4:6;
  values 1.8;
  endogenize x;

  exogenize y;
  periods 2004Y:2006Y;
  values 0.6;
  endogenize e_y;
end;

perfect_foresight_with_expectation_errors_setup(periods=30, first_simulation_period = 2001Y);
perfect_foresight_with_expectation_errors_solver;

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if oo_.endo_simul(1, 3) ~= 1.6 || any(oo_.endo_simul(1, 5:7) ~= 1.8) || any(oo_.endo_simul(2, 8:10) ~= 13) ...
   || oo_.endo_simul(3, 4) ~= 0.5 || any(oo_.endo_simul(3, 5:7) ~= 0.6) || oo_.exo_simul(2, 1) ~= 1.2
   error('Perfect foresight with expectation errors and controlled paths failed')
end
