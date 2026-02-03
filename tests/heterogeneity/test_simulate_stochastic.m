function [testFailed, testResults] = test_simulate_stochastic(M_, options_, oo_, steady_state, testFailed, testResults)
%test_simulate_stochastic Tests for heterogeneity.simulate in stochastic mode
%
% Tests core functionality for IRF computation and stochastic simulations
% including combined IRF + simulation mode
% Initializes and solves the model before running simulation tests
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure
%   oo_          [struct]  Dynare results structure
%   steady_state [struct]  Steady state structure
%   testFailed   [scalar]  Number of tests failed so far
%   testResults  [array]   Array of test result structures
%
% OUTPUTS:
%   testFailed   [scalar]  Updated number of tests failed
%   testResults  [array]   Updated array of test results

fprintf('\n');
fprintf('Testing heterogeneity.simulate function\n');
fprintf('===============================\n');

% Verify that model was already initialized and solved
try
    assert(isfield(oo_, 'heterogeneity') && isfield(oo_.heterogeneity, 'dr') && isfield(oo_.heterogeneity.dr, 'G'), ...
        'Solution not available (model not initialized/solved)');
    fprintf('✔ Using pre-initialized and pre-solved model\n');
catch ME
    fprintf('❌ Model not properly initialized: %s\n', ME.message);
    testResults = [testResults; struct('name', 'Model verification', 'passed', false, 'message', ME.message)];
    testFailed = testFailed + 1;
    return;
end

%% Test 1: Default IRFs for all shocks
try
    options_.nograph = true;
    oo_ = heterogeneity.simulate(M_, options_, oo_, {});
    assert(isfield(oo_, 'irfs'), 'oo_.irfs must exist');

    % Check that IRFs exist for all shock-variable pairs
    n_irfs = 0;
    for i_shock = 1:M_.exo_nbr
        shock_name = M_.exo_names{i_shock};
        for i_var = 1:M_.endo_nbr
            var_name = M_.endo_names{i_var};
            field_name = [var_name '_' shock_name];
            if isfield(oo_.irfs, field_name)
                n_irfs = n_irfs + 1;
                assert(isvector(oo_.irfs.(field_name)), sprintf('IRF for %s must be a vector', field_name));
            end
        end
    end
    fprintf('  oo_.irfs contains %d IRFs for %d shocks\n', n_irfs, M_.exo_nbr);
    testResults = [testResults; struct('name', 'Default IRFs for all shocks', 'passed', true, 'message', '')];
    fprintf('✔ Default IRFs for all shocks\n');
catch ME
    testResults = [testResults; struct('name', 'Default IRFs for all shocks', 'passed', false, 'message', ME.message)];
    fprintf('❌ Default IRFs for all shocks: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test 2: IRF for subset of shocks using irf_shocks option
try
    shock_name = M_.exo_names{1};
    options_.irf_shocks = {shock_name};
    options_.nograph = true;
    oo_ = heterogeneity.simulate(M_, options_, oo_, {});
    assert(isfield(oo_, 'irfs'), 'oo_.irfs must exist');

    % Check that IRFs exist only for the specified shock
    for i_var = 1:M_.endo_nbr
        var_name = M_.endo_names{i_var};
        field_name = [var_name '_' shock_name];
        if isfield(oo_.heterogeneity.dr.G, var_name) && isfield(oo_.heterogeneity.dr.G.(var_name), shock_name)
            assert(isfield(oo_.irfs, field_name), sprintf('Variable %s_%s missing from oo_.irfs', var_name, shock_name));
        end
    end
    fprintf('  IRFs computed for shock "%s" only\n', shock_name);
    testResults = [testResults; struct('name', 'IRF for subset of shocks', 'passed', true, 'message', '')];
    fprintf('✔ IRF for subset of shocks\n');
catch ME
    testResults = [testResults; struct('name', 'IRF for subset of shocks', 'passed', false, 'message', ME.message)];
    fprintf('❌ IRF for subset of shocks: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test 3: Variable subset using positional arguments
try
    shock_name = M_.exo_names{1};
    if M_.endo_nbr >= 3
        var_subset = M_.endo_names(1:3);
    else
        var_subset = M_.endo_names;
    end
    options_.irf_shocks = {shock_name};
    options_.nograph = true;
    oo_ = heterogeneity.simulate(M_, options_, oo_, var_subset);

    % Only the requested variables should have IRFs stored
    for i = 1:length(var_subset)
        field_name = [var_subset{i} '_' shock_name];
        if isfield(oo_.heterogeneity.dr.G, var_subset{i}) && isfield(oo_.heterogeneity.dr.G.(var_subset{i}), shock_name)
            assert(isfield(oo_.irfs, field_name), ...
                sprintf('Variable %s_%s missing from oo_.irfs', var_subset{i}, shock_name));
        end
    end
    fprintf('  Subset of %d variables computed correctly\n', length(var_subset));
    testResults = [testResults; struct('name', 'Variable subset', 'passed', true, 'message', '')];
    fprintf('✔ Variable subset\n');
catch ME
    testResults = [testResults; struct('name', 'Variable subset', 'passed', false, 'message', ME.message)];
    fprintf('❌ Variable subset: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test 4: Custom IRF horizon
try
    shock_name = M_.exo_names{1};
    T_custom = 50;
    options_.irf = T_custom;
    options_.irf_shocks = {shock_name};
    options_.nograph = true;
    oo_ = heterogeneity.simulate(M_, options_, oo_, {});
    var_name = M_.endo_names{1};
    field_name = [var_name '_' shock_name];
    assert(length(oo_.irfs.(field_name)) == T_custom, ...
        sprintf('Expected horizon %d, got %d', T_custom, length(oo_.irfs.(field_name))));
    fprintf('  Custom horizon irf=%d verified\n', T_custom);
    testResults = [testResults; struct('name', 'Custom IRF horizon', 'passed', true, 'message', '')];
    fprintf('✔ Custom IRF horizon\n');
catch ME
    testResults = [testResults; struct('name', 'Custom IRF horizon', 'passed', false, 'message', ME.message)];
    fprintf('❌ Custom IRF horizon: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test 5: Stochastic simulation with periods option
try
    T_sim = 100;
    T_drop = 20;

    % Check if all shocks have non-zero variance
    if isfield(M_, 'Sigma_e') && ~isempty(M_.Sigma_e) && all(diag(M_.Sigma_e) > 0)
        options_.periods = T_sim;
        options_.drop = T_drop;
        options_.nograph = true;
        oo_ = heterogeneity.simulate(M_, options_, oo_, {});
        assert(isfield(oo_, 'endo_simul'), 'oo_.endo_simul must exist for stochastic simulation');
        assert(size(oo_.endo_simul, 1) == M_.endo_nbr, 'oo_.endo_simul should have n_endo rows');
        assert(size(oo_.endo_simul, 2) == T_sim, sprintf('Expected %d periods, got %d', T_sim, size(oo_.endo_simul, 2)));
        fprintf('  Stochastic simulation: %d periods, %d drop\n', T_sim, T_drop);
        testResults = [testResults; struct('name', 'Stochastic simulation with periods', 'passed', true, 'message', '')];
        fprintf('✔ Stochastic simulation with periods\n');
    else
        fprintf('  Skipped: Some shocks have zero variance\n');
        testResults = [testResults; struct('name', 'Stochastic simulation with periods', 'passed', true, 'message', 'Skipped (zero variance shocks)')];
        fprintf('✔ Stochastic simulation with periods (skipped)\n');
    end
catch ME
    testResults = [testResults; struct('name', 'Stochastic simulation with periods', 'passed', false, 'message', ME.message)];
    fprintf('❌ Stochastic simulation with periods: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test 6: Error handling - unknown shock in irf_shocks
try
    try
        options_.irf_shocks = {'nonexistent_shock'};
        options_.nograph = true;
        oo_ = heterogeneity.simulate(M_, options_, oo_, {});
        error('Should have thrown an error for unknown shock');
    catch ME
        assert(strcmp(ME.identifier, 'heterogeneity:simulate:UnknownShock'), ...
            'Wrong error identifier');
    end
    fprintf('  Unknown shock correctly rejected\n');
    testResults = [testResults; struct('name', 'Error handling: unknown shock', 'passed', true, 'message', '')];
    fprintf('✔ Error handling: unknown shock\n');
catch ME
    testResults = [testResults; struct('name', 'Error handling: unknown shock', 'passed', false, 'message', ME.message)];
    fprintf('❌ Error handling: unknown shock: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test 7: Plot generation (conditional - only if graphics enabled)
if ~options_.nograph
    try
        shock_name = M_.exo_names{1};
        options_.irf_shocks = {shock_name};
        options_.nograph = false;
        oo_ = heterogeneity.simulate(M_, options_, oo_, {});
        assert(isfield(oo_, 'irfs'), 'oo_.irfs must exist');
        assert(~isempty(fieldnames(oo_.irfs)), 'oo_.irfs should not be empty');

        % Check if any IRFs exceed threshold (only then should plots be created)
        threshold = 1e-10;
        if isfield(options_, 'impulse_responses') && isfield(options_.impulse_responses, 'plot_threshold')
            threshold = options_.impulse_responses.plot_threshold;
        end

        irf_fields = fieldnames(oo_.irfs);
        max_irf = 0;
        for i = 1:length(irf_fields)
            max_irf = max(max_irf, max(abs(oo_.irfs.(irf_fields{i}))));
        end

        if max_irf >= threshold
            % IRFs exceed threshold, plots should be created
            assert(isfolder([M_.dname '/graphs']), 'Graphs directory should exist');
            graph_files = dir([M_.dname '/graphs/' M_.fname '_IRF_' shock_name '*']);
            assert(~isempty(graph_files), 'Plot file should have been created');
            fprintf('  Created %d graph file(s)\n', length(graph_files));
        else
            % IRFs below threshold, no plots should be created
            fprintf('  IRFs below threshold, no plots created (as expected)\n');
        end
        testResults = [testResults; struct('name', 'Plot generation (default)', 'passed', true, 'message', '')];
        fprintf('✔ Plot generation (default)\n');
    catch ME
        testResults = [testResults; struct('name', 'Plot generation (default)', 'passed', false, 'message', ME.message)];
        fprintf('❌ Plot generation (default): %s\n', ME.message);
        testFailed = testFailed + 1;
    end

    %% Test 8: No plot with nograph option
    try
        % Use second shock if available
        shock_idx = min(2, M_.exo_nbr);
        shock_name = M_.exo_names{shock_idx};
        if isfolder([M_.dname '/graphs'])
            existing = dir([M_.dname '/graphs/' M_.fname '_IRF_' shock_name '*']);
            for i = 1:length(existing)
                delete([M_.dname '/graphs/' existing(i).name]);
            end
        end
        options_.irf_shocks = {shock_name};
        options_.nograph = true;
        oo_ = heterogeneity.simulate(M_, options_, oo_, {});
        assert(isfield(oo_, 'irfs'), 'oo_.irfs must exist');
        graph_files = dir([M_.dname '/graphs/' M_.fname '_IRF_' shock_name '*']);
        assert(isempty(graph_files), 'No plot files should have been created');
        fprintf('  Correctly skipped plotting\n');
        testResults = [testResults; struct('name', 'No plot with nograph', 'passed', true, 'message', '')];
        fprintf('✔ No plot with nograph\n');
    catch ME
        testResults = [testResults; struct('name', 'No plot with nograph', 'passed', false, 'message', ME.message)];
        fprintf('❌ No plot with nograph: %s\n', ME.message);
        testFailed = testFailed + 1;
    end
end

%% Test 9: relative_irf option
try
    % Find a shock with non-zero variance
    shock_idx = 1;
    for i = 1:M_.exo_nbr
        if isfield(M_, 'Sigma_e') && ~isempty(M_.Sigma_e) && size(M_.Sigma_e, 1) >= i
            if M_.Sigma_e(i, i) > 0
                shock_idx = i;
                break;
            end
        end
    end
    shock_name = M_.exo_names{shock_idx};

    % Compute IRFs with and without relative_irf
    options_.periods = 0;  % Reset to default (IRFs only, no simulation paths)
    options_.irf_shocks = {shock_name};
    options_.relative_irf = false;
    options_.nograph = true;
    oo_std = heterogeneity.simulate(M_, options_, oo_, {});

    options_.relative_irf = true;
    oo_rel = heterogeneity.simulate(M_, options_, oo_, {});

    % Get shock standard deviation
    shock_std = sqrt(M_.Sigma_e(shock_idx, shock_idx));

    % Check relationship: relative_irf should give (1/shock_std) * 100 times the standard IRF
    var_name = M_.endo_names{1};
    field_name = [var_name '_' shock_name];

    if isfield(oo_std.irfs, field_name) && isfield(oo_rel.irfs, field_name)
        % Check if IRF values are non-zero (can be zero if shock has zero variance)
        irf_std_val = oo_std.irfs.(field_name)(1);
        irf_rel_val = oo_rel.irfs.(field_name)(1);

        if abs(irf_std_val) > 1e-14 && abs(irf_rel_val) > 1e-14
            % relative_irf uses unit shock * 100, standard uses shock_std
            % So: oo_rel = (1/shock_std) * 100 * oo_std
            expected_ratio = 100 / shock_std;
            actual_ratio = irf_rel_val / irf_std_val;

            assert(abs(actual_ratio - expected_ratio) < 1e-10, ...
                sprintf('relative_irf ratio incorrect: expected %.6f, got %.6f', expected_ratio, actual_ratio));

            fprintf('  relative_irf scaling verified (ratio: %.2f)\n', actual_ratio);
            testResults = [testResults; struct('name', 'relative_irf option', 'passed', true, 'message', '')];
            fprintf('✔ relative_irf option\n');
        else
            fprintf('  Skipped: IRF values too small (shock has zero/near-zero variance)\n');
            testResults = [testResults; struct('name', 'relative_irf option', 'passed', true, 'message', 'Skipped (zero variance)')];
        end
    else
        fprintf('  Skipped: IRF field not found (likely below threshold)\n');
        testResults = [testResults; struct('name', 'relative_irf option', 'passed', true, 'message', 'Skipped (IRF below threshold)')];
    end
catch ME
    testResults = [testResults; struct('name', 'relative_irf option', 'passed', false, 'message', ME.message)];
    fprintf('❌ relative_irf option: %s\n', ME.message);
    testFailed = testFailed + 1;
end

%% Test: Combined IRFs and stochastic simulation
fprintf('\n--- Testing combined IRFs + stochastic simulation ---\n');
try
    options_.periods = 100;
    options_.drop = 50;
    options_.irf = 50;
    options_.irf_shocks = {};  % All shocks
    options_.nograph = true;

    oo_ = heterogeneity.simulate(M_, options_, oo_, {});

    % Verify both outputs exist
    assert(isfield(oo_, 'irfs'), 'IRFs should be computed when periods > 0');
    assert(isfield(oo_, 'endo_simul'), 'Simulation should be computed when periods > 0');

    % Check simulation dimensions
    expected_cols = options_.periods;  % No initial state column (matches stoch_simul)
    actual_cols = size(oo_.endo_simul, 2);
    assert(actual_cols == expected_cols, ...
           sprintf('Simulation length mismatch: expected %d, got %d', expected_cols, actual_cols));

    % Check IRFs exist for at least one variable-shock pair
    irf_fields = fieldnames(oo_.irfs);
    assert(~isempty(irf_fields), 'At least one IRF should be computed');

    % Check IRF length
    first_irf = oo_.irfs.(irf_fields{1});
    assert(length(first_irf) == options_.irf, ...
           sprintf('IRF length mismatch: expected %d, got %d', options_.irf, length(first_irf)));

    fprintf('✓ Combined IRFs + simulation mode works correctly\n');
    testResults = [testResults; struct('name', 'combined IRFs + simulation', 'passed', true, 'message', '')];
catch ME
    testResults = [testResults; struct('name', 'combined IRFs + simulation', 'passed', false, 'message', ME.message)];
    fprintf('❌ Combined IRFs + simulation: %s\n', ME.message);
    testFailed = testFailed + 1;
end

end
