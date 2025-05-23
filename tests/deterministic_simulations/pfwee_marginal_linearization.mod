/* Test that perfect foresight simulation with expectation errors gives correct
   results with marginal linearization. Homotopy is not really triggered, but
   the use of homotopy_max_completion_share is enough to test the computations.

   This test is similar to pfwee_homotopy_marginal_linearization.mod, except that
   the numerical values are different (the shock is smaller to facilitate computations).
   */

var Consumption, Capital, LoggedProductivity;

varexo LoggedProductivityInnovation;

parameters beta, alpha, delta, rho;

beta = .985;
alpha = 1/3;
delta = alpha/10;
rho = .4;

model;
  1/Consumption = beta/Consumption(1)*(alpha*exp(LoggedProductivity(1))*Capital^(alpha-1)+1-delta);
  Capital = exp(LoggedProductivity)*Capital(-1)^alpha+(1-delta)*Capital(-1)-Consumption;
  LoggedProductivity = rho*LoggedProductivity(-1)+LoggedProductivityInnovation;
end;

initval;
  LoggedProductivityInnovation = 0.1;
end;

steady;

endval;
  LoggedProductivityInnovation = 0.3;
end;

shocks(learnt_in = 5);
  var LoggedProductivityInnovation;
  periods 7;
  values -0.5;
end;

endval(learnt_in = 10);
  LoggedProductivityInnovation = 0.6;
end;

// Save initial steady state (it will be modified by pfwee)
orig_steady_state = oo_.steady_state;
orig_exo_steady_state = oo_.exo_steady_state;

perfect_foresight_with_expectation_errors_setup(periods = 50);
perfect_foresight_with_expectation_errors_solver(homotopy_max_completion_share = 0.8, homotopy_marginal_linearization_fallback, steady_solve_algo = 13);

verbatim;

if ~oo_.deterministic_simulation.status
   error('Perfect foresight simulation failed')
end

pfwee_endo_simul = oo_.endo_simul;
pfwee_exo_simul = oo_.exo_simul;

oo_.steady_state = orig_steady_state;
oo_.exo_steady_state = orig_exo_steady_state;

oo_ = make_ex_(M_,options_,oo_);
oo_ = make_y_(M_,options_,oo_);

options_.simul.homotopy_marginal_linearization_fallback = 0;
options_.simul.homotopy_max_completion_share = 1;
options_.solve_algo = 13;

%% Information arriving in period 1 (permanent shock now)
% Simulation at 80%
disp('*** info_period=1, share=80%')
oo_.exo_steady_state = 0.1+0.2*0.8;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:,end) = oo_.steady_state;
oo_.exo_simul(2:end,:) = oo_.exo_steady_state;
sim1_1 = perfect_foresight_solver(M_, options_, oo_);
% Simulation at 79%
disp('*** info_period=1, share=79%')
oo_.exo_steady_state = 0.1+0.2*0.79;
oo_.steady_state = evaluate_steady_state(oo_.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:,end) = oo_.steady_state;
oo_.exo_simul(2:end,:) = oo_.exo_steady_state;
sim2_1 = perfect_foresight_solver(M_, options_, oo_);

%% Information arriving in period 5 (temporary shock in period 7)
options_.periods = 46;
% Simulation at 80%
disp('*** info_period=5, share=80%')
oo_.endo_simul = sim1_1.endo_simul(:,5:end);
oo_.exo_simul = sim1_1.exo_simul(5:end,:);
oo_.steady_state = sim1_1.steady_state;
oo_.exo_steady_state = sim1_1.exo_steady_state;
oo_.exo_simul(4,1) = 0.1-0.6*0.8;
sim1_5 = perfect_foresight_solver(M_, options_, oo_);
% Simulation at 79%
disp('*** info_period=5, share=79%')
oo_.endo_simul = sim2_1.endo_simul(:,5:end);
oo_.exo_simul = sim2_1.exo_simul(5:end,:);
oo_.steady_state = sim2_1.steady_state;
oo_.exo_steady_state = sim2_1.exo_steady_state;
oo_.exo_simul(4,1) = 0.1-0.6*0.79;
sim2_5 = perfect_foresight_solver(M_, options_, oo_);

%% Information arriving in period 10 (permanent shock now)
options_.periods = 41;
% Simulation at 80%
disp('*** info_period=10, share=80%')
oo_.endo_simul = sim1_5.endo_simul(:,6:end);
oo_.exo_simul = sim1_5.exo_simul(6:end,:);
oo_.exo_steady_state = 0.1+0.5*0.8;
oo_.steady_state = evaluate_steady_state(sim1_5.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:,end) = oo_.steady_state;
oo_.exo_simul(2:end,:) = oo_.exo_steady_state;
sim1_10 = perfect_foresight_solver(M_, options_, oo_);
% Simulation at 79%
disp('*** info_period=10, share=79%')
oo_.endo_simul = sim2_5.endo_simul(:,6:end);
oo_.exo_simul = sim2_5.exo_simul(6:end,:);
oo_.exo_steady_state = 0.1+0.5*0.79;
oo_.steady_state = evaluate_steady_state(sim2_5.steady_state, oo_.exo_steady_state, M_, options_, true);
oo_.endo_simul(:,end) = oo_.steady_state;
oo_.exo_simul(2:end,:) = oo_.exo_steady_state;
sim2_10 = perfect_foresight_solver(M_, options_, oo_);

%% Reconstruct the two simulations
sim1_endo_simul = [ sim1_1.endo_simul(:,1:5), sim1_5.endo_simul(:,2:6), sim1_10.endo_simul(:,2:end);];
sim2_endo_simul = [ sim2_1.endo_simul(:,1:5), sim2_5.endo_simul(:,2:6), sim2_10.endo_simul(:,2:end);];
sim1_exo_simul = [ sim1_1.exo_simul(1:5,:); sim1_5.exo_simul(2:6,:); sim1_10.exo_simul(2:end,:);];
sim2_exo_simul = [ sim2_1.exo_simul(1:5,:); sim2_5.exo_simul(2:6,:); sim2_10.exo_simul(2:end,:);];

%% Compute marginal linearization
endo_simul = sim1_endo_simul + 20*(sim1_endo_simul - sim2_endo_simul);
exo_simul = sim1_exo_simul + 20*(sim1_exo_simul - sim2_exo_simul);

endo_abs_diff = max(max(abs(pfwee_endo_simul-endo_simul)));
endo_rel_diff = max(max(abs(endo_simul./pfwee_endo_simul-1)));

fprintf('Max absolute difference on endogenous: %g\nMax relative difference on endogenous: %g%%\n', endo_abs_diff, 100*endo_rel_diff)

exo_abs_diff = max(max(abs(pfwee_exo_simul-exo_simul)));
exo_rel_diff = max(max(abs(exo_simul./pfwee_exo_simul-1)));

fprintf('Max absolute difference on exogenous: %g\nMax relative difference on exogenous: %g%%\n', exo_abs_diff, 100*exo_rel_diff)

if endo_abs_diff > 1e-3 || endo_rel_diff > 3e-3 || exo_abs_diff > 1e-15 || exo_rel_diff > 1e-13
  error('perfect foresight with expectation errors + marginal linearization gives wrong result')
end

end; // verbatim
