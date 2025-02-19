/* Example that triggers homotopy in perfect foresight simulation with
   expectation errors + controlled paths, and tests marginal linearization. */

@#include "pfwee_homotopy.inc"

perfect_foresight_controlled_paths(learnt_in = 7);
  exogenize LoggedProductivity;
  periods 8:9;
  values 5;
  endogenize LoggedProductivityInnovation;
end;

perfect_foresight_with_expectation_errors_setup(periods=200);
perfect_foresight_with_expectation_errors_solver(homotopy_max_completion_share = 0.8, homotopy_marginal_linearization_fallback, steady_solve_algo = 13);

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

if ~(all(abs(oo_.endo_simul(3, 9:10) - 5) < 1e-14) && abs(oo_.exo_simul(10) - 0.5) < 1e-14)
   error('Homotopy with marginal linearization and controlled paths failed')
end
