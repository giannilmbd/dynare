function [testFailed, testResults] = test_permutation(M_, options_, oo_, steady_state, testFailed, testResults)
%test_permutation Tests that permutation logic works correctly for pol.order and d.order
%
% Verifies that:
% 1. Valid alternative orderings are accepted by check_steady_state_input
% 2. Output arrays are correctly permuted to canonical order
% 3. heterogeneity_solve produces identical G matrices regardless of input order
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure
%   oo_          [struct]  Dynare output structure (must have heterogeneity.mat.G)
%   steady_state [struct]  Steady state structure with canonical order
%   testFailed   [scalar]  Number of tests failed so far
%   testResults  [array]   Array of test result structures
%
% OUTPUTS:
%   testFailed   [scalar]  Updated number of tests failed
%   testResults  [array]   Updated array of test results

n_dims = numel(steady_state.pol.order);
if n_dims == 2
    fprintf('\nTesting permutation handling (2D)\n');
    fprintf('=================================\n');
else
    fprintf('\nTesting permutation handling (%dD)\n', n_dims);
    fprintf('==================================\n');
end

% Store baseline G from canonical order (already solved)
G_baseline = oo_.heterogeneity.mat.G;

% Get canonical order and sizes
canonical_order = steady_state.pol.order;
sizes = zeros(1, n_dims);
for i = 1:n_dims
    var = canonical_order{i};
    if isfield(steady_state.shocks.grids, var)
        sizes(i) = numel(steady_state.shocks.grids.(var));
    else
        sizes(i) = numel(steady_state.pol.grids.(var));
    end
end

% Create a non-trivial permutation (circular shift: [2,3,...,n,1])
perm_vec = [2:n_dims, 1];
new_order = canonical_order(perm_vec);

% Create permuted steady state
ss_perm = steady_state;
ss_perm.pol.order = new_order;
pol_fields = fieldnames(ss_perm.pol.values);
for i = 1:numel(pol_fields)
    ss_perm.pol.values.(pol_fields{i}) = permute(steady_state.pol.values.(pol_fields{i}), perm_vec);
end
ss_perm.d.order = new_order;
ss_perm.d.hist = permute(steady_state.d.hist, perm_vec);

%% Unit test: check_steady_state_input restores canonical order
try
    [out_ss, ~, ~] = heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_perm);

    % Verify dimensions are canonical
    assert(isequal(size(out_ss.pol.values.(canonical_order{end})), sizes), ...
        sprintf('pol.values.%s not permuted to canonical dimensions', canonical_order{end}));
    assert(isequal(size(out_ss.d.hist), sizes), 'd.hist not permuted to canonical dimensions');

    % Verify values match original (numerical correctness)
    ref_var = canonical_order{end};  % Use last state variable for comparison
    max_diff_pol = max(abs(out_ss.pol.values.(ref_var)(:) - steady_state.pol.values.(ref_var)(:)));
    assert(max_diff_pol < 1e-14, sprintf('pol.values.%s values differ after permutation: %.2e', ref_var, max_diff_pol));
    max_diff_hist = max(abs(out_ss.d.hist(:) - steady_state.d.hist(:)));
    assert(max_diff_hist < 1e-14, sprintf('d.hist values differ after permutation: %.2e', max_diff_hist));

    fprintf('✔ check_steady_state_input correctly permutes arrays to canonical order\n');
    testResults = [testResults; struct('name', sprintf('Permutation %dD: unit test', n_dims), 'passed', true, 'message', '')];
catch ME
    testFailed = testFailed + 1;
    fprintf('✘ check_steady_state_input permutation test failed: %s\n', ME.message);
    testResults = [testResults; struct('name', sprintf('Permutation %dD: unit test', n_dims), 'passed', false, 'message', ME.message)];
    return;  % Skip integration test if unit test fails
end

%% Integration test: solve with permuted input, compare G matrices
try
    % Load and solve with permuted data directly (no temp file needed)
    oo_.heterogeneity = heterogeneity.load_steady_state(M_, options_.heterogeneity, oo_.heterogeneity, ss_perm);
    oo_.heterogeneity.dr = heterogeneity.solve(M_, options_.heterogeneity.solve, oo_.heterogeneity);

    G_permuted = oo_.heterogeneity.mat.G;

    % Compare G matrices - should be identical
    max_diff_G = max(abs(G_baseline(:) - G_permuted(:)));
    assert(max_diff_G < 1e-10, sprintf('G matrices differ: max diff = %.2e', max_diff_G));

    fprintf('✔ heterogeneity_solve produces identical G with permuted input (diff: %.2e)\n', max_diff_G);
    testResults = [testResults; struct('name', sprintf('Permutation %dD: integration test', n_dims), 'passed', true, 'message', '')];
catch ME
    testFailed = testFailed + 1;
    fprintf('✘ heterogeneity_solve permutation test failed: %s\n', ME.message);
    testResults = [testResults; struct('name', sprintf('Permutation %dD: integration test', n_dims), 'passed', false, 'message', ME.message)];
end

end
