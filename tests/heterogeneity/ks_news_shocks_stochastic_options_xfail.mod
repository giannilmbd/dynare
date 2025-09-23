// EXPECTEDFAIL Test: News shocks + heterogeneity_simulate stochastic options
//
// This test verifies that the preprocessor correctly rejects .mod files that
// combine news shocks (anticipated shocks known at t=0) with stochastic
// simulation options in heterogeneity_simulate.
//
// NEWS SHOCKS are triggered by shocks block with 'periods' and 'values' keywords
// STOCHASTIC OPTIONS include: irf, periods, drop, irf_shocks, relative_irf
//
// These modes are mutually exclusive (Auclert et al. 2021 terminology).
//
// USAGE: Set TEST_OPTION macro to select which option to test:
//   -DTEST_OPTION=periods      : Test periods option
//   -DTEST_OPTION=irf          : Test irf option
//   -DTEST_OPTION=irf_shocks   : Test irf_shocks option
//   -DTEST_OPTION=relative_irf : Test relative_irf option
//   -DTEST_OPTION=combined     : Test combined irf+periods
//
// Expected result: Preprocessor error during dynare command execution

@#echomacrovars

// Minimal Krusell-Smith model definition
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
    Z_ss  // Aggregate TFP shock average value
;

// Household optimization problem
model(heterogeneity=households);
    beta*Va(+1)-c^(-1/eis)=0 ⟂ a>=0;
    (1+r)*a(-1)+w*e-c-a;
    Va = (1+r)*c^(-1/eis);
end;

// Aggregate equilibrium conditions
model;
    (Z_ss+Z) * K(-1)^alpha * L^(1 - alpha) - Y;
    alpha * (Z_ss+Z) * (K(-1) / L)^(alpha - 1) - delta - r;
    (1 - alpha) * (Z_ss+Z) * (K(-1) / L)^alpha - w;
    K - SUM(a);
end;

//==========================================================================
// NEWS SHOCKS DEFINITION
// Using 'periods' and 'values' keywords triggers news shock sequence mode
//==========================================================================
shocks;
    var Z;
    periods 1:10;
    values 0.01;
end;

heterogeneity_load_steady_state(filename=ks);

heterogeneity_solve;

//==========================================================================
// CONFLICTING HETEROGENEITY_SIMULATE COMMAND
// This should trigger preprocessor validation error
//==========================================================================
@#if TEST_OPTION == "periods"
    // ERROR: 'periods' option conflicts with news shock mode
    heterogeneity_simulate(periods = 100);

@#elseif TEST_OPTION == "irf"
    // ERROR: 'irf' option conflicts with news shock mode
    heterogeneity_simulate(irf = 50);

@#elseif TEST_OPTION == "irf_shocks"
    // ERROR: 'irf_shocks' option conflicts with news shock mode
    heterogeneity_simulate(irf_shocks = (Z));

@#elseif TEST_OPTION == "relative_irf"
    // ERROR: 'relative_irf' option conflicts with news shock mode
    heterogeneity_simulate(relative_irf);

@#elseif TEST_OPTION == "combined"
    // ERROR: Combined options conflict with news shock mode
    heterogeneity_simulate(irf = 50, periods = 100);

@#else
@#error "Invalid TEST_OPTION: must be periods, irf, irf_shocks, relative_irf, or combined"
@#endif
