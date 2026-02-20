heterogeneity_dimension households;

var(heterogeneity=households)
    c  // Consumption
    a  // Assets
    Va // Derivative of the value function w.r.t assets
;

varexo(heterogeneity=households)
    e
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
    alpha // Share of capital in production fuction
    beta  // Subjective discount rate of houselholds
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

shocks;
    var Z; stderr 0.01;
end;

verbatim;
    mat_file = 'ks_ssj.mat';
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
end;

% Initialize steady state once for all tests
heterogeneity_load_steady_state(filename = ks_ssj);

verbatim;
    testFailed = 0;
    testResults = [];  % Array to collect all test results

    % Test the solve routine with custom impulse responses (calls solve with custom impulses)
    [testFailed, result] = run_test('heterogeneity.solve with custom impulse responses', @() test_solve_with_custom_impulse_responses(M_, options_, oo_, x_hat_dash, G, J, F, curlyDs, curlyYs, interpolation_indices, interpolation_weights), testFailed);
    testResults = [testResults; result];

    % Print test summary
    print_test_summary(testResults);

    if testFailed > 0
        error('Some unit tests failed!');
    end
end;