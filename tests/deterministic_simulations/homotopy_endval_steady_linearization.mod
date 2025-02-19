// Example that triggers homotopy in perfect foresight simulation.
// Tests the endval_steady, homotopy_linearization_fallback options, and a few more.

@#include "homotopy_linearization.inc"

perfect_foresight_setup(periods=200, endval_steady);

perfect_foresight_solver(homotopy_initial_step_size = 0.5,
                         homotopy_max_completion_share = 0.6,
                         homotopy_min_step_size = 0.0001,
                         homotopy_linearization_fallback,
                         steady_solve_algo = 13, steady_tolf=1e-7);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end
