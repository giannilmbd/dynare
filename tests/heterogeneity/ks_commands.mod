// Comprehensive test for all heterogeneity command variations
// Tests all options and syntax variations
//
// USAGE: This file tests two distinct simulation modes:
//
// 1. STOCHASTIC SIMULATION MODE (default, no flags):
//    dynare-preprocessor ks_commands.mod
//    Tests: heterogeneity_simulate with various options:
//      - IRF options: irf, relative_irf, irf_shocks
//      - Simulation options: periods, drop
//      - Combined: Both IRFs and stochastic simulation together
//    Output: oo_.irfs (always computed) and oo_.endo_simul (when periods > 0)
//
// 2. NEWS SHOCK SEQUENCE MODE (with -DDET_SHOCKS=1):
//    dynare-preprocessor ks_commands.mod -DDET_SHOCKS=1
//    Tests: heterogeneity_simulate with anticipated shock sequences (news shocks)
//    When shocks are defined with periods/values in shocks block, agents learn
//    at t=0 about shocks at various future dates (Auclert et al. 2021)
//    Output: oo_.endo_simul with perfect foresight paths
//
// All modes test heterogeneity_load_steady_state and heterogeneity_solve commands

heterogeneity_dimension households;

var(heterogeneity=households)
    c  // Consumption
    a  // Assets
    Va // Derivative of the value function w.r.t assets
;

varexo(heterogeneity=households)
    e  // Idiosyncratic productivity shock
;

var
    Y // Aggregate output
    r // Rate of return on capital net of depreciation
    w // Wage rate
    K // Aggregate capital
;

varexo Z; // Aggregate productivity shock

parameters
    L     // Labor
    alpha // Share of capital in production function
    beta  // Subjective discount rate of households
    delta // Capital depreciation rate
    eis   // Elasticity of intertemporal substitution
    rho_Z // Aggregate TFP shock persistence
    sig_Z // Aggregate TFP shock innovation std err
    Z_ss  // Aggregate TFP shock average value
;

model(heterogeneity=households);
    beta*Va(+1)-c^(-1/eis)=0 ⟂ a>=0;
    (1+r)*a(-1)+w*e-c-a;
    Va = (1+r)*c^(-1/eis);
end;

model;
    (Z_ss+Z) * K(-1)^alpha * L^(1 - alpha) - Y;
    alpha * (Z_ss+Z) * (K(-1) / L)^(alpha - 1) - delta - r;
    (1 - alpha) * (Z_ss+Z) * (K(-1) / L)^alpha - w;
    K - SUM(a);
end;

@#ifdef DET_SHOCKS
// News shock sequence (anticipated shocks known at t=0)
shocks;
    var Z;
    periods 1:10;
    values 0.01;
end;
@#else
// Stochastic shocks (for IRFs and stochastic simulation)
shocks;
    var Z; stderr 0.01;
end;
@#endif

//==========================================================================
// TEST 1: heterogeneity_load_steady_state with both options
//==========================================================================
heterogeneity_load_steady_state(filename = ks);

verbatim;
    disp('TEST 1: heterogeneity_load_steady_state(filename=ks)');
    assert(isfield(oo_, 'heterogeneity'), 'TEST 1 FAILED: Missing oo_.heterogeneity');
    assert(isfield(oo_.heterogeneity, 'steady_state'), 'TEST 1 FAILED: Missing steady state');
    assert(isfield(oo_.heterogeneity.steady_state, 'agg'), 'TEST 1 FAILED: Missing aggregate variables');
    assert(isfield(oo_.heterogeneity.steady_state, 'pol'), 'TEST 1 FAILED: Missing policy functions');
    disp('  ✓ Steady state loaded successfully');
    disp(' ');
end;

//==========================================================================
// TEST 2: heterogeneity_solve with custom time horizon
//==========================================================================
heterogeneity_solve(truncation_horizon = 400);

verbatim;
    disp('TEST 2: heterogeneity_solve(truncation_horizon=400)');
    assert(isfield(options_.heterogeneity.solve, 'truncation_horizon'), 'TEST 2 FAILED: truncation_horizon option not set');
    assert(options_.heterogeneity.solve.truncation_horizon == 400, 'TEST 2 FAILED: truncation_horizon option not set to 400');
    assert(isfield(oo_.heterogeneity, 'dr'), 'TEST 2 FAILED: Missing solution');
    assert(isfield(oo_.heterogeneity.dr, 'G'), 'TEST 2 FAILED: Missing Jacobians');
    disp('  ✓ Solution computed with truncation_horizon=400');
    disp(' ');
end;

@#ifndef DET_SHOCKS
//==========================================================================
// STOCHASTIC SIMULATION MODE TESTS (default mode)
// Tests IRF computation and stochastic simulation options
//==========================================================================

//==========================================================================
// TEST 3: heterogeneity_simulate - basic IRFs
//==========================================================================
heterogeneity_simulate;

verbatim;
    disp('TEST 3: heterogeneity_simulate (basic IRFs)');
    assert(isfield(oo_, 'irfs'), 'TEST 3 FAILED: Missing oo_.irfs');
    assert(~isempty(fieldnames(oo_.irfs)), 'TEST 3 FAILED: IRFs are empty');
    disp('  ✓ IRFs computed successfully');
    disp(' ');
end;

//==========================================================================
// TEST 4: heterogeneity_simulate with irf option
//==========================================================================
heterogeneity_simulate(irf = 80);

verbatim;
    disp('TEST 4: heterogeneity_simulate(irf=80)');
    assert(options_.irf == 80, 'TEST 4 FAILED: irf option not set to 80');
    disp('  ✓ irf option correctly set to 80');
    disp(' ');
end;

//==========================================================================
// TEST 5: heterogeneity_simulate with relative_irf
//==========================================================================
heterogeneity_simulate(relative_irf);

verbatim;
    disp('TEST 5: heterogeneity_simulate(relative_irf)');
    assert(isfield(options_, 'relative_irf'), 'TEST 5 FAILED: relative_irf option not set');
    assert(options_.relative_irf == true, 'TEST 5 FAILED: relative_irf not true');
    disp('  ✓ relative_irf option correctly set');
    disp(' ');
end;

//==========================================================================
// TEST 6: heterogeneity_simulate with multiple options
//==========================================================================
heterogeneity_simulate(irf = 60, relative_irf);

verbatim;
    disp('TEST 6: heterogeneity_simulate(irf=60, relative_irf)');
    assert(options_.irf == 60, 'TEST 6 FAILED: irf option not set to 60');
    assert(options_.relative_irf == true, 'TEST 6 FAILED: relative_irf not set');
    disp('  ✓ Both options correctly set');
    disp(' ');
end;

//==========================================================================
// TEST 7: heterogeneity_simulate with irf_shocks option
//==========================================================================
heterogeneity_simulate(irf_shocks = (Z));

verbatim;
    disp('TEST 7: heterogeneity_simulate(irf_shocks=(Z))');
    assert(isfield(options_, 'irf_shocks'), 'TEST 7 FAILED: irf_shocks option not set');
    assert(isequal(options_.irf_shocks, {'Z'}), 'TEST 7 FAILED: irf_shocks not set to {''Z''}');
    disp('  ✓ irf_shocks correctly set to {''Z''}');
    disp(' ');
end;

//==========================================================================
// TEST 8: heterogeneity_simulate with variable list
//==========================================================================
heterogeneity_simulate Y K;

verbatim;
    disp('TEST 8: heterogeneity_simulate Y K');
    % Variable list is passed to simulate function, no option to check
    disp('  ✓ Variable list syntax accepted');
    disp(' ');
end;

//==========================================================================
// TEST 9: heterogeneity_simulate with options AND variable list
//==========================================================================
heterogeneity_simulate(irf = 50) w r;

verbatim;
    disp('TEST 9: heterogeneity_simulate(irf=50) w r');
    assert(options_.irf == 50, 'TEST 9 FAILED: irf option not set to 50');
    disp('  ✓ Options and variable list both accepted');
    disp(' ');
end;

//==========================================================================
// TEST 10: heterogeneity_simulate with many options
//==========================================================================
heterogeneity_simulate(irf = 70, irf_shocks = (Z), relative_irf);

verbatim;
    disp('TEST 10: heterogeneity_simulate(irf=70, irf_shocks=(Z), relative_irf)');
    assert(options_.irf == 70, 'TEST 10 FAILED: irf option not set to 70');
    assert(isequal(options_.irf_shocks, {'Z'}), 'TEST 10 FAILED: irf_shocks not correct');
    assert(options_.relative_irf == true, 'TEST 10 FAILED: relative_irf not set');
    disp('  ✓ All three options correctly set');
    disp(' ');
end;

//==========================================================================
// TEST 13: heterogeneity_simulate with periods option
//==========================================================================
heterogeneity_simulate(periods = 200);

verbatim;
    disp('TEST 13: heterogeneity_simulate(periods=200)');
    assert(isfield(options_, 'periods'), 'TEST 13 FAILED: periods option not set');
    assert(options_.periods == 200, 'TEST 13 FAILED: periods not set to 200');
    assert(isfield(oo_, 'endo_simul'), 'TEST 13 FAILED: Missing oo_.endo_simul');
    disp('  ✓ Stochastic simulation with periods=200');
    disp(' ');
end;

//==========================================================================
// TEST 14: heterogeneity_simulate with periods and drop options
//==========================================================================
heterogeneity_simulate(periods = 300, drop = 50);

verbatim;
    disp('TEST 14: heterogeneity_simulate(periods=300, drop=50)');
    assert(options_.periods == 300, 'TEST 14 FAILED: periods not set to 300');
    assert(isfield(options_, 'drop'), 'TEST 14 FAILED: drop option not set');
    assert(options_.drop == 50, 'TEST 14 FAILED: drop not set to 50');
    disp('  ✓ Both periods and drop options correctly set');
    disp(' ');
end;

//==========================================================================
// TEST 15: heterogeneity_simulate with periods, drop, and variable list
//==========================================================================
heterogeneity_simulate(periods = 250, drop = 25) Y K w r;

verbatim;
    disp('TEST 15: heterogeneity_simulate(periods=250, drop=25) Y K w r');
    assert(options_.periods == 250, 'TEST 15 FAILED: periods not set to 250');
    assert(options_.drop == 25, 'TEST 15 FAILED: drop not set to 25');
    disp('  ✓ Options and variable list both accepted');
    disp(' ');
end;

//==========================================================================
// TEST 16: heterogeneity_simulate with combined IRFs and stochastic simulation
//==========================================================================
heterogeneity_simulate(irf = 50, periods = 200, drop = 50);

verbatim;
    disp('TEST 16: heterogeneity_simulate(irf=50, periods=200, drop=50)');
    assert(options_.irf == 50, 'TEST 16 FAILED: irf option not set to 50');
    assert(options_.periods == 200, 'TEST 16 FAILED: periods not set to 200');
    assert(options_.drop == 50, 'TEST 16 FAILED: drop not set to 50');
    assert(isfield(oo_, 'irfs'), 'TEST 16 FAILED: Missing oo_.irfs');
    assert(isfield(oo_, 'endo_simul'), 'TEST 16 FAILED: Missing oo_.endo_simul');
    disp('  ✓ Combined IRFs and stochastic simulation work together');
    disp(' ');
end;

@#else
//==========================================================================
// NEWS SHOCK SEQUENCE MODE TESTS (when DET_SHOCKS is defined)
// Tests simulation with anticipated shocks known at t=0
//==========================================================================

//==========================================================================
// TEST 11: heterogeneity_simulate with news shocks
//==========================================================================
// When news shocks are defined in the shocks block (with periods/values),
// heterogeneity_simulate automatically uses news shock mode
heterogeneity_simulate;

verbatim;
    disp('TEST 11: heterogeneity_simulate (news shocks)');
    assert(isfield(oo_, 'endo_simul'), 'TEST 11 FAILED: Missing oo_.endo_simul');
    assert(~isempty(M_.det_shocks), 'TEST 11 FAILED: News shocks not configured');
    disp('  ✓ News shock simulation computed successfully');
    disp(' ');
end;

//==========================================================================
// TEST 12: heterogeneity_simulate with news shocks and variable list
//==========================================================================
heterogeneity_simulate Y K w;

verbatim;
    disp('TEST 12: heterogeneity_simulate Y K w (news shocks with var list)');
    % Variable list is passed to simulate function, no option to check
    disp('  ✓ Variable list syntax accepted in news shock mode');
    disp(' ');
end;

@#endif

//==========================================================================
// Summary
//==========================================================================
verbatim;
    disp('=================================================================');
@#ifndef DET_SHOCKS
    disp('ALL 16 TESTS PASSED (STOCHASTIC SIMULATION MODE)! ✓✓✓');
    disp('  - Tests 1-2: Load and solve');
    disp('  - Tests 3-10: IRF options');
    disp('  - Tests 13-16: Stochastic simulation and combined mode');
@#else
    disp('ALL 4 TESTS PASSED (NEWS SHOCK SEQUENCE MODE)! ✓✓✓');
    disp('  - Tests 1-2: Load and solve');
    disp('  - Tests 11-12: News shock simulation');
@#endif
    disp('=================================================================');
end;
