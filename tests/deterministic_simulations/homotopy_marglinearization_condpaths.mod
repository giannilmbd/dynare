// Example that triggers homotopy in perfect foresight simulation.
// Tests homotopy_marginal_linearization_fallback with perfect_foresight_controlled_paths

@#include "homotopy_linearization.inc"

perfect_foresight_controlled_paths;
  exogenize LoggedProductivity;
  periods 8:9;
  values 5;
  endogenize LoggedProductivityInnovation;
end;

perfect_foresight_setup(periods=200, endval_steady);

perfect_foresight_solver(homotopy_max_completion_share = 0.7,
                         homotopy_marginal_linearization_fallback,
                         steady_solve_algo = 13);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if ~(all(abs(oo_.endo_simul(3, 9:10) - 5) < 1e-14) && abs(oo_.exo_simul(10) - 0.5) < 1e-13)
   error('Homotopy with marginal linearization and controlled paths failed')
end
