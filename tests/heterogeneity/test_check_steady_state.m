function [testFailed, testResults] = test_check_steady_state(M_, options_, steady_state, testFailed, testResults)
%test_check_steady_state Tests for heterogeneity.check_steady_state_input
%
% Comprehensive tests for the steady-state input validation function
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure
%   steady_state           [struct]  Steady state structure
%   testFailed   [scalar]  Number of tests failed so far
%   testResults  [array]   Array of test result structures
%
% OUTPUTS:
%   testFailed   [scalar]  Updated number of tests failed
%   testResults  [array]   Updated array of test results

fprintf('\n');
fprintf('Testing heterogeneity.check_steady_state_input\n');
fprintf('=============================================\n');

verbose = true;

%% === OUTPUTS ===
[testFailed, result] = run_test('sizes output validation', @() test_sizes_output(M_, options_, steady_state), testFailed);
testResults = [testResults; result];

[testFailed, result] = run_test('out_ss output validation', @() test_out_ss_output(M_, options_, steady_state), testFailed);
testResults = [testResults; result];

function test_sizes_output(M_, options_, steady_state)
    [out_ss, sizes] = heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, steady_state);
    assert(sizes.n_e == numel(fieldnames(steady_state.shocks.grids)), 'sizes.n_e incorrect');
    assert(sizes.n_a == numel(fieldnames(steady_state.pol.grids)), 'sizes.n_a incorrect');
    assert(sizes.N_e == numel(steady_state.shocks.grids.e), 'sizes.N_e incorrect');
    assert(sizes.n_pol == M_.heterogeneity(1).endo_nbr, 'sizes.n_pol incorrect');
    assert(sizes.agg == numel(fieldnames(steady_state.agg)), 'sizes.agg incorrect');
    assert(sizes.shocks.e == numel(steady_state.shocks.grids.e), 'sizes.shocks.e incorrect');
    assert(sizes.pol.N_a == numel(steady_state.pol.grids.a), 'sizes.pol.N_a incorrect');
    assert(sizes.pol.states.a == sizes.pol.N_a, 'sizes.pol.states.a incorrect');
    assert(sizes.d.N_a == sizes.pol.N_a, 'sizes.d.N_a incorrect');
    assert(sizes.d.states.a == sizes.pol.states.a, 'sizes.d.states.a incorrect');
end

function test_out_ss_output(M_, options_, steady_state)
    [out_ss, sizes] = heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, steady_state);
    assert(iscolumn(out_ss.shocks.grids.e), 'out_ss.shocks.grids.e not column');
    assert(iscolumn(out_ss.pol.grids.a), 'out_ss.pol.grids.a not column');
    assert(iscolumn(out_ss.d.grids.a), 'out_ss.d.grids.a not column');
    assert(size(out_ss.shocks.grids.e,1) == sizes.shocks.e, 'out_ss.shocks.grids.e size mismatch');
    assert(size(out_ss.pol.grids.a,1) == sizes.pol.states.a, 'out_ss.pol.grids.a size mismatch');
    assert(size(out_ss.d.grids.a,1) == sizes.d.states.a, 'out_ss.d.grids.a size mismatch');
    assert(numel(fieldnames(out_ss.agg)) == sizes.agg, 'out_ss.agg field count mismatch');
end

%% === NON-STRUCT FIELDS ===

% steady_state is not a struct
ss_test = 123;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), 'steady_state is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.shocks is not a struct
ss_test = steady_state;
ss_test.shocks = 1;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), 'steady_state.shocks is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.shocks.grids is not a struct
ss_test = steady_state;
ss_test.shocks.grids = 42;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), ...
    'steady_state.shocks.grids is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.shocks.Pi is not a struct
ss_test = steady_state;
ss_test.shocks.Pi = "not_a_struct";
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), ...
    'steady_state.shocks.Pi is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.pol is not a struct
ss_test = steady_state;
ss_test.pol = pi;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), ...
    'steady_state.pol is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.pol.grids is not a struct
ss_test = steady_state;
ss_test.pol.grids = "grid";
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), ...
    'steady_state.pol.grids is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.pol.values is not a struct
ss_test = steady_state;
ss_test.pol.values = "values";
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), ...
    'steady_state.pol.values is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.d is not a struct
ss_test = steady_state;
ss_test.d = 1;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), ...
    'steady_state.d is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.d.grids is not a struct
ss_test = steady_state;
ss_test.d.grids = "grid";
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), ...
    'steady_state.d.grids is not a struct', testFailed, verbose);
testResults = [testResults; result];

% steady_state.agg is not a struct
ss_test = steady_state;
ss_test.agg = 0;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity.check, ss_test), ...
    'steady_state.agg is not a struct', testFailed, verbose);
testResults = [testResults; result];

end
