function [testFailed, testResults] = test_simulate_news_shocks(M_, options_, oo_, ss, sim_result, testFailed, testResults)
%test_simulate_news_shocks Tests for heterogeneity.simulate with news shock sequences
%
% Tests news shock sequence mode using M_.det_shocks populated from shocks block
% with anticipated shocks (periods/values keywords) known at t=0.
% Initializes and solves the model before running simulation tests.
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure
%   oo_          [struct]  Dynare results structure
%   ss           [struct]  Steady state structure
%   sim_result   [struct]  Python simulation results (can be empty)
%   testFailed   [scalar]  Number of tests failed so far
%   testResults  [array]   Array of test result structures
%
% OUTPUTS:
%   testFailed   [scalar]  Updated number of tests failed
%   testResults  [array]   Updated array of test results

fprintf('\n');
fprintf('Testing heterogeneity.simulate with news shock sequences\n');
fprintf('=========================================================\n');

% Verify that model was already initialized and solved
try
    assert(isfield(oo_, 'heterogeneity') && isfield(oo_.heterogeneity, 'dr') && isfield(oo_.heterogeneity.dr, 'G'), ...
        'Solution not available (model not initialized/solved)');
    fprintf('✔ Using pre-initialized and pre-solved model\n');
catch ME
    fprintf('❌ Model not properly initialized: %s\n', ME.message);
    testResults = [testResults; struct('name', 'Model verification (news shocks)', 'passed', false, 'message', ME.message)];
    testFailed = testFailed + 1;
    return;
end

%% Test 1: News shock simulation with M_.det_shocks
try
    % Verify M_.det_shocks exists (should be populated by shocks block with periods/values)
    assert(isfield(M_, 'det_shocks') && ~isempty(M_.det_shocks), ...
        'M_.det_shocks must exist for news shock test');

    % Run news shock simulation (no 'periods' parameter needed - inferred from options_.periods)
    options_.nograph = true;
    oo_det = heterogeneity.simulate(M_, options_, oo_, {});

    % Verify outputs
    assert(isfield(oo_det, 'endo_simul'), 'oo_.endo_simul must exist for news shock simulation');
    assert(isfield(oo_det, 'exo_simul'), 'oo_.exo_simul must exist (populated by make_ex_)');
    assert(size(oo_det.endo_simul, 1) == M_.endo_nbr, 'oo_.endo_simul should have n_endo rows');

    fprintf('  News shock simulation completed with %d periods\n', size(oo_det.endo_simul, 2) - 1);
    testResults = [testResults; struct('name', 'News shock simulation basic', 'passed', true, 'message', '')];
    fprintf('✔ News shock simulation basic\n');
catch ME
    testResults = [testResults; struct('name', 'News shock simulation basic', 'passed', false, 'message', ME.message)];
    fprintf('❌ News shock simulation basic: %s\n', ME.message);
    testFailed = testFailed + 1;
    return;  % Cannot proceed with other tests
end

%% Test 2: Verify shock values in oo_.exo_simul
try
    % Check that oo_.exo_simul contains the shock values from M_.det_shocks
    for i = 1:length(M_.det_shocks)
        shock_entry = M_.det_shocks(i);
        if shock_entry.exo_det
            continue;  % Skip deterministic exogenous variables
        end

        shock_idx = shock_entry.exo_id;
        shock_name = M_.exo_names{shock_idx};
        periods = shock_entry.periods;
        expected_value = shock_entry.value;

        % Periods in exo_simul are offset by M_.maximum_lag
        exo_periods = M_.maximum_lag + periods;

        for p_idx = 1:length(exo_periods)
            p = exo_periods(p_idx);
            if p <= size(oo_det.exo_simul, 1)
                actual_value = oo_det.exo_simul(p, shock_idx);
                if abs(actual_value - expected_value) >= 1e-10
                    error('Shock %s at period %d incorrect: expected %.6f, got %.6f', ...
                            shock_name, periods(p_idx), expected_value, actual_value);
                end
            end
        end
    end

    fprintf('  Shock values verified in oo_.exo_simul\n');
    testResults = [testResults; struct('name', 'Shock values verification', 'passed', true, 'message', '')];
    fprintf('✔ Shock values verification\n');
catch ME
    testResults = [testResults; struct('name', 'Shock values verification', 'passed', false, 'message', ME.message)];
    fprintf('❌ Shock values verification: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test 3: Verify endogenous response is non-zero

% Build the vector of steady-state aggregate endogenous variables
steady_state = zeros(M_.endo_nbr, 1);
agg_vars = fieldnames(oo_.heterogeneity.ss.agg);
for i = 1:length(agg_vars)
    var_name = agg_vars{i};
    var_idx = find(strcmp(var_name, M_.endo_names));
    if ~isempty(var_idx)
        steady_state(var_idx) = oo_.heterogeneity.ss.agg.(var_name);
    end
end

try
    % Check that at least one endogenous variable responds to news shocks
    Y_idx = find(strcmp('Y', M_.endo_names));
    if ~isempty(Y_idx)
        Y_path = oo_det.endo_simul(Y_idx, :);
        Y_deviation = Y_path - steady_state(Y_idx);
        max_deviation = max(abs(Y_deviation));

        assert(max_deviation > 1e-10, ...
            'Output Y should respond to news shocks');

        fprintf('  Output Y max deviation from steady state: %.6f\n', max_deviation);
        testResults = [testResults; struct('name', 'Endogenous response verification', 'passed', true, 'message', '')];
        fprintf('✔ Endogenous response verification\n');
    else
        fprintf('  Skipped: Variable Y not found in model\n');
        testResults = [testResults; struct('name', 'Endogenous response verification', 'passed', true, 'message', 'Skipped (Y not found)')];
    end
catch ME
    testResults = [testResults; struct('name', 'Endogenous response verification', 'passed', false, 'message', ME.message)];
    fprintf('❌ Endogenous response verification: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test 4: Compare MATLAB simulation with Python sim_result
try
    % Check if Python sim_result was passed as input
    if ~isempty(sim_result) && isstruct(sim_result)
        fprintf('  Comparing MATLAB simulation with Python sim_result...\n');

        max_diff = 0;
        max_diff_var = '';
        n_compared = 0;

        % Compare each variable that exists in both MATLAB and Python results
        var_names = fieldnames(sim_result);
        for i = 1:length(var_names)
            var_name = var_names{i};

            % Find variable index in MATLAB output
            var_idx = find(strcmp(var_name, M_.endo_names));
            if isempty(var_idx)
                continue;
            end

            % Extract MATLAB simulation (deviations from steady state)
            matlab_sim = oo_det.endo_simul(var_idx, :)' - steady_state(var_idx);

            % Extract Python simulation
            python_sim = sim_result.(var_name).';

            % Compare only over common length
            T_compare = min(length(matlab_sim), length(python_sim));
            diff = abs(matlab_sim(1:T_compare) - python_sim(1:T_compare));
            var_max_diff = max(diff);

            if var_max_diff > max_diff
                max_diff = var_max_diff;
                max_diff_var = var_name;
            end

            n_compared = n_compared + 1;
        end

        % Tolerance for numerical differences
        tol = 1e-5;

        if n_compared > 0
            fprintf('  Compared %d variables\n', n_compared);
            fprintf('  Maximum difference: %.2e (variable: %s)\n', max_diff, max_diff_var);

            assert(max_diff < tol, ...
                sprintf('MATLAB-Python mismatch exceeds tolerance: %.2e > %.2e for variable %s', ...
                max_diff, tol, max_diff_var));

            testResults = [testResults; struct('name', 'MATLAB vs Python simulation consistency', 'passed', true, 'message', '')];
            fprintf('✔ MATLAB vs Python simulation consistency\n');
        else
            fprintf('  No variables found for comparison\n');
            testResults = [testResults; struct('name', 'MATLAB vs Python simulation consistency', 'passed', true, 'message', 'Skipped (no common variables)')];
        end
    else
        fprintf('  Python sim_result not provided, skipping comparison\n');
        testResults = [testResults; struct('name', 'MATLAB vs Python simulation consistency', 'passed', true, 'message', 'Skipped (no Python data)')];
    end
catch ME
    testResults = [testResults; struct('name', 'MATLAB vs Python simulation consistency', 'passed', false, 'message', ME.message)];
    fprintf('❌ MATLAB vs Python simulation consistency: %s\n', ME.message);
    testFailed = testFailed + 1;
end

end
