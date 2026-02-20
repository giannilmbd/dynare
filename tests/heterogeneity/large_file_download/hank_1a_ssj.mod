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

shocks;
    var G; stderr 0.01;
    var markup; stderr 0.01;
    var rstar; stderr 0.01;
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

mat_file = 'hank_1a_ssj.mat';
if ~exist(mat_file, 'file')
    url = sprintf('https://dynare.org/test-data/%s', mat_file);
    websave(mat_file, url);
end
load(mat_file);

param_names = fieldnames(steady_state.params);
for i=1:numel(param_names)
    param = param_names{i};
    set_param_value(param, steady_state.params.(param));
end

% Initialize steady state once for all tests
heterogeneity_load_steady_state(filename = hank_1a_ssj);

% Solve once for all tests
heterogeneity_solve;

verbatim;
    testFailed = 0;
    testResults = [];  % Array to collect all test results

    % Test the solve routine and SSJ comparisons (validation only)
    [testFailed, result] = run_test('heterogeneity.solve functionality and SSJ validation', @() test_solve_functionality(M_, oo_, G, true), testFailed);
    testResults = [testResults; result];

    % Test the solve routine with custom impulse responses (calls solve with custom impulses)
    [testFailed, result] = run_test('heterogeneity.solve with custom impulse responses', @() test_solve_with_custom_impulse_responses(M_, options_, oo_, x_hat_dash, G, J, F, curlyDs, curlyYs, a_i, a_pi), testFailed);
    testResults = [testResults; result];

    % Test heterogeneity.simulate function (reuse initialized and solved model)
    [testFailed, testResults] = test_simulate_stochastic(M_, options_, oo_, steady_state, testFailed, testResults);

    % Print test summary
    print_test_summary(testResults);

    if testFailed > 0
        error('Some unit tests failed!');
    end
end;