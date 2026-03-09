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

load 'ks.mat';

param_names = fieldnames(steady_state.params);
for i=1:numel(param_names)
    param = param_names{i};
    set_param_value(param, steady_state.params.(param));
end

verbatim;
    skipline()
    disp('*** TESTING: heterogeneity.check_steady_state_input.m ***');
    options_.heterogeneity.check.no_warning_redundant = false;
    options_.heterogeneity.check.no_warning_d_grids = false;
    testFailed = 0;
    testResults = [];  % Array to collect all test results

    % Test heterogeneity.check_steady_state_input function
    [testFailed, testResults] = test_check_steady_state(M_, options_, steady_state, testFailed, testResults);
end;

%% === FIELD VALIDATION TESTS ===
verbatim;
    % Field validation tests
    verbose = true;
    [testFailed, testResults] = test_field_validation(M_, options_, steady_state, testFailed, testResults);
end;

%% === FILE LOADING TESTS ===
% Test file loading with various path formats, custom variables, and error handling

% Test 1: Load with explicit .mat extension + validate residuals
heterogeneity_load_steady_state(filename = 'ks.mat');
assert(isfield(oo_, 'heterogeneity'), 'oo_.heterogeneity field missing');
assert(isfield(oo_.heterogeneity, 'steady_state'), 'oo_.heterogeneity.ss field missing');
assert(isfield(oo_.heterogeneity, 'sizes'), 'oo_.heterogeneity.sizes field missing');
assert(isfield(oo_.heterogeneity, 'mat'), 'oo_.heterogeneity.mat field missing');
assert(isfield(oo_.heterogeneity, 'indices'), 'oo_.heterogeneity.indices field missing');
verbatim;
    % Validate steady-state residuals
    agg_r = max(abs(oo_.heterogeneity.mat.G));
    assert(agg_r <= 1e-9, sprintf('Aggregate residuals too large: %.2e', agg_r));

    % Check heterogeneous residuals with Euler equation scaling
    het_r = oo_.heterogeneity.mat.F;
    euler_idx = 1;  % First equation is Euler equation in ks.mod
    scaling = reshape(steady_state.pol.values.c .^ (-1/eis), 1, []);

    if isoctave()
        het_r_euler = max(max(abs(het_r(euler_idx,:) ./ scaling)));
        het_r_other = max(max(abs(het_r(setdiff(1:size(het_r,1), euler_idx),:))));
    else
        het_r_euler = max(abs(het_r(euler_idx,:) ./ scaling), [], "all");
        het_r_other = max(abs(het_r(setdiff(1:size(het_r,1), euler_idx),:)), [], "all");
    end

    assert(het_r_euler <= 1e-2, sprintf('Euler residuals (scaled) too large: %.2e', het_r_euler));
    assert(het_r_other <= 5e-4, sprintf('Other heterogeneous residuals too large: %.2e', het_r_other));

    fprintf('✔ Test 1: Explicit .mat extension + residuals validated\n');
    testResults = [testResults; struct('name', 'Load: explicit .mat + residuals', 'passed', true, 'message', '')];
end;

% Test 2: Load without extension (auto-detect)
heterogeneity_load_steady_state(filename = ks);
assert(isfield(oo_, 'heterogeneity'), 'oo_.heterogeneity field missing');
assert(isfield(oo_.heterogeneity, 'steady_state'), 'oo_.heterogeneity.ss field missing');
assert(isfield(oo_.heterogeneity, 'sizes'), 'oo_.heterogeneity.sizes field missing');
assert(isfield(oo_.heterogeneity, 'mat'), 'oo_.heterogeneity.mat field missing');
assert(isfield(oo_.heterogeneity, 'indices'), 'oo_.heterogeneity.indices field missing')
verbatim;
    fprintf('✔ Test 2: Auto-detect .mat extension works\n');
    testResults = [testResults; struct('name', 'Path: auto-detect extension', 'passed', true, 'message', '')];
end;

% Test 3: Load with relative path to current directory
heterogeneity_load_steady_state(filename = './ks.mat');
assert(isfield(oo_, 'heterogeneity'), 'oo_.heterogeneity field missing');
assert(isfield(oo_.heterogeneity, 'steady_state'), 'oo_.heterogeneity.ss field missing');
assert(isfield(oo_.heterogeneity, 'sizes'), 'oo_.heterogeneity.sizes field missing');
assert(isfield(oo_.heterogeneity, 'mat'), 'oo_.heterogeneity.mat field missing');
assert(isfield(oo_.heterogeneity, 'indices'), 'oo_.heterogeneity.indices field missing')
verbatim;
    fprintf('✔ Test 3: Relative path (./) works\n');
    testResults = [testResults; struct('name', 'Path: relative path ./', 'passed', true, 'message', '')];
end;

% Test 4: Load with custom variable name
verbatim;
    % Create temp file with custom variable name
    temp_file = 'ks_test_temp.mat';
    steady_state_temp = steady_state;
    save(temp_file, 'steady_state_temp', '-v7');

    % Set options for next command
    options_.heterogeneity.steady_state_file_name = temp_file;
    options_.heterogeneity.steady_state_variable_name = 'steady_state_temp';
end;

heterogeneity_load_steady_state;

verbatim;
    assert(isfield(oo_.heterogeneity, 'steady_state'), 'Failed to load custom variable');
    delete(temp_file);  % Clean up

    fprintf('✔ Test 4: Custom variable name works\n');
    testResults = [testResults; struct('name', 'Load: custom variable name', 'passed', true, 'message', '')];

    % Reset to defaults
    options_.heterogeneity.steady_state_file_name = '';
    options_.heterogeneity.steady_state_variable_name = 'steady_state';
end;

% Test 5: Error on missing variable
verbatim;
    error_caught = false;
    try
        options_.heterogeneity.steady_state_file_name = 'ks.mat';
        options_.heterogeneity.steady_state_variable_name = 'nonexistent_field';
        oo_test = heterogeneity.load_steady_state(M_, options_.heterogeneity, oo_.heterogeneity);
    catch ME
        error_caught = true;
        assert(contains(ME.message, 'Variable ''nonexistent_field'' not found'), ...
            sprintf('Wrong error message: %s', ME.message));
    end
    assert(error_caught, 'Expected error for missing variable not raised');

    fprintf('✔ Test 5: Error on missing variable works\n');
    testResults = [testResults; struct('name', 'Error: missing variable', 'passed', true, 'message', '')];

    % Reset
    options_.heterogeneity.steady_state_file_name = '';
    options_.heterogeneity.steady_state_variable_name = 'steady_state';
end;

% Test 6: Error on missing file
verbatim;
    error_caught = false;
    try
        options_.heterogeneity.steady_state_file_name = 'nonexistent_file.mat';
        oo_test = heterogeneity.load_steady_state(M_, options_.heterogeneity, oo_.heterogeneity);
    catch ME
        error_caught = true;
        assert(contains(ME.message, 'nonexistent_file.mat') || contains(lower(ME.message), 'not found'), ...
            sprintf('Wrong error message: %s', ME.message));
    end
    assert(error_caught, 'Expected error for missing file not raised');

    fprintf('✔ Test 6: Error on missing file works\n');
    testResults = [testResults; struct('name', 'Error: missing file', 'passed', true, 'message', '')];

    % Reset
    options_.heterogeneity.steady_state_file_name = '';
end;

% Test 7: Error when filename is empty and workspace variable does not exist
verbatim;
    error_caught = false;
    try
        options_.heterogeneity.steady_state_file_name = '';
        options_.heterogeneity.steady_state_variable_name = 'nonexistent_variable';
        oo_test = heterogeneity.load_steady_state(M_, options_.heterogeneity, oo_.heterogeneity);
    catch ME
        error_caught = true;
        assert(contains(ME.message, 'not found in workspace'), ...
            sprintf('Wrong error message: %s', ME.message));
    end
    assert(error_caught, 'Expected error for missing workspace variable not raised');

    % Reset
    options_.heterogeneity.steady_state_variable_name = 'steady_state';

    fprintf('✔ Test 7: Error on missing workspace variable works\n');
    testResults = [testResults; struct('name', 'Error: missing workspace variable', 'passed', true, 'message', '')];

    fprintf('\n');
end;

% Initialize steady state once for all tests
heterogeneity_load_steady_state(filename = ks);

% Solve once for all tests
heterogeneity_solve;

verbatim;
    % Test the solve routine and SSJ comparisons (validation only)
    [testFailed, result] = run_test('heterogeneity.solve functionality and SSJ validation', @() test_solve_functionality(M_, oo_, G_eq, true), testFailed);
    testResults = [testResults; result];

    % Test heterogeneity.simulate function (reuse initialized and solved model)
    [testFailed, testResults] = test_simulate_stochastic(M_, options_, oo_, steady_state, testFailed, testResults);

    % Test permutation handling
    [testFailed, testResults] = test_permutation(M_, options_, oo_, steady_state, testFailed, testResults);

    % Print test summary
    print_test_summary(testResults);

    if testFailed > 0
        error('Some unit tests failed!');
    end
end;