heterogeneity_dimension households;

var(heterogeneity=households)
   c      // consumption
   n      // labor supply
   ns     // effective labor supply
   a      // assets
   Va     // derivative of value function w.r.t. a
;

varexo(heterogeneity=households)
   e      // idiosyncratic efficiency
;

var
   Y L w pi Div Tax r
;

varexo G markup rstar;

@#ifdef DET_SHOCKS
seq = 0.01*0.8.^(0:49);
@#endif

shocks;
@#ifndef DET_SHOCKS
    var G; stderr 0.01;
    var markup; stderr 0.01;
    var rstar; stderr 0.01;
@#else
    var G;
    periods 1;
    values 0.01;
    var markup;
    periods 4:5;
    values 0.01;
    var rstar;
    periods 1:50;
    values (seq);
@#endif
end;

parameters
   beta vphi
   eis frisch
   rho_e sig_e
   rho_Z sig_Z
   mu kappa phi
   Z B r_ss
;

model(heterogeneity=households);
   // Euler equation with borrowing constraint
   beta * Va(+1) - c^(-1/eis) = 0 ⟂ a >= 0;

   // Budget constraint
   (1 + r) * a(-1) + w * n * e + (Div-Tax) * e - c - a;

   // Envelope condition
   Va = (1 + r) * c^(-1/eis);

   // Intratemporal FOC for labor supply
   vphi*n^(1/frisch) - w*e*c^(-1/eis);

   // Definition of the effective labor supply 
   ns = n * e;
end;

model;
   // Firm
   L - Y / Z;
   Div - (Y - w * L - mu / (mu - 1) / (2 * kappa) * log(1 + pi)^2 * Y);

   // Monetary policy (Taylor rule)
   (1 + r_ss + rstar(-1) + phi * pi(-1)) / (1 + pi) - 1 - r;

   // Fiscal policy
   Tax - (r * B) - G;

   // NKPC
   kappa * (w / Z - 1 / mu)
   + Y(+1)/Y * log(1 + pi(+1)) / (1 + r(+1))
   + markup
   - log(1 + pi);

   // Market clearing
   SUM(a) - B;  // asset market
   sum(ns) - L; // labor market
end;

load 'hank_1a.mat';

param_names = fieldnames(steady_state.params);
for i=1:numel(param_names)
    param = param_names{i};
    set_param_value(param, steady_state.params.(param));
end

% Initialize steady state once for all tests
heterogeneity_load_steady_state(filename = hank_1a);
% Check that required fields exist
assert(isfield(oo_, 'heterogeneity'), 'oo_.heterogeneity field missing');
assert(isfield(oo_.heterogeneity, 'ss'), 'oo_.heterogeneity.ss field missing');
assert(isfield(oo_.heterogeneity, 'sizes'), 'oo_.heterogeneity.sizes field missing');
assert(isfield(oo_.heterogeneity, 'mat'), 'oo_.heterogeneity.mat field missing');
assert(isfield(oo_.heterogeneity, 'indices'), 'oo_.heterogeneity.indices field missing');

% Solve once for all tests
heterogeneity_solve;

verbatim;
    testFailed = 0;
    testResults = [];  % Array to collect all test results

% Test heterogeneity.simulate function (reuse initialized and solved model)
@#ifdef DET_SHOCKS
    [testFailed, testResults] = test_simulate_news_shocks(M_, options_, oo_, steady_state, sim_result, testFailed, testResults);
@#else
    [testFailed, testResults] = test_simulate_stochastic(M_, options_, oo_, steady_state, testFailed, testResults);
@#endif

    % Print test summary
    print_test_summary(testResults);

    if testFailed > 0
        error('Some unit tests failed!');
    end
end;