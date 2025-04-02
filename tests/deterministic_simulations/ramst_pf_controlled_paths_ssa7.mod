// Performs a simulation with a perfect_foresight_controlled_paths block
// with the stack_solve_algo=7

@#define shocks = true
@#include "ramst_pf_controlled_paths_common.inc"

perfect_foresight_setup(periods=30, first_simulation_period = 2001Y);
perfect_foresight_solver(stack_solve_algo = 7);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if oo_.endo_simul(1, 3) ~= 1.6 || any(oo_.endo_simul(1, 5:6) ~= 1.7) || any(oo_.endo_simul(2, 8:10) ~= 13) ...
   || any(oo_.endo_simul(3, 4:6) ~= 0.5) || oo_.exo_simul(2, 1) ~= 1.2
   error('Perfect foresight with controlled paths failed')
end
