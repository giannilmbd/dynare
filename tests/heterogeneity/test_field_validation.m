function [testFailed, testResults] = test_field_validation(M_, options_, ss, testFailed, testResults)
%test_field_validation Tests for field validation in heterogeneity.check_steady_state_input
%
% Comprehensive tests for field existence, types, and sizes in steady-state structure
%
% INPUTS:
%   M_           [struct]  Dynare model structure
%   options_     [struct]  Dynare options structure
%   ss           [struct]  Steady state structure
%   testFailed   [scalar]  Number of tests failed so far
%   testResults  [array]   Array of test result structures
%
% OUTPUTS:
%   testFailed   [scalar]  Updated number of tests failed
%   testResults  [array]   Updated array of test results

fprintf('\n');
fprintf('Testing field validation\n');
fprintf('========================\n');

verbose = true;

%% === SHOCK FIELD TESTS ===

% Missing `shocks` field
ss_test = ss;
ss_test = rmfield(ss_test, 'shocks');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing `shocks` field', testFailed, verbose);
testResults = [testResults; result];

% Wrong type for `shocks.grids.e`
ss_test = ss;
ss_test.shocks.grids.e = ones(3,3);  % invalid: matrix instead of vector
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.grids.e are not dense real vectors', testFailed, verbose);
testResults = [testResults; result];

% Wrong type for `shocks.Pi.e` (sparse instead of dense)
ss_test = ss;
ss_test.shocks.Pi.e = speye(numel(ss.shocks.grids.e));
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.Pi.e is not a dense real matrix', testFailed, verbose);
testResults = [testResults; result];

% Markov matrix wrong row number
ss_test = ss;
ss_test.shocks.Pi.e = rand(3, 7);  % Should be square matching grid
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.Pi.e row count does not match', testFailed, verbose);
testResults = [testResults; result];

% Markov matrix wrong column number
ss_test = ss;
ss_test.shocks.Pi.e = rand(7, 3);  % Should be square matching grid
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.Pi.e column count does not match', testFailed, verbose);
testResults = [testResults; result];

% Markov matrix has negative elements
ss_test = ss;
ss_test.shocks.Pi.e(1,1) = -0.1;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.Pi.e contains negative elements', testFailed, verbose);
testResults = [testResults; result];

% Markov matrix does not sum to 1
ss_test = ss;
ss_test.shocks.Pi.e = ss.shocks.Pi.e * 2;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.Pi.e has row sums different from 1', testFailed, verbose);
testResults = [testResults; result];

%% === POLICY TESTS ===

% Missing `pol` field
ss_test = ss;
ss_test = rmfield(ss_test, 'pol');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing `pol` field', testFailed, verbose);
testResults = [testResults; result];

% Missing `pol.grids`
ss_test = ss;
ss_test.pol = rmfield(ss.pol, 'grids');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing `pol.grids`', testFailed, verbose);
testResults = [testResults; result];

% Missing one state grid
ss_test = ss;
ss_test.pol.grids = rmfield(ss.pol.grids, 'a');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing state grid for a', testFailed, verbose);
testResults = [testResults; result];

% Wrong type for `pol.grids.a`
ss_test = ss;
ss_test.pol.grids.a = ones(3, 3);  % should be a vector
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.grids.a should be a vector', testFailed, verbose);
testResults = [testResults; result];

% Missing `pol.values`
ss_test = ss;
ss_test.pol = rmfield(ss.pol, 'values');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing `pol.values`', testFailed, verbose);
testResults = [testResults; result];

% Missing policy value for declared variable
ss_test = ss;
ss_test.pol.values = rmfield(ss.pol.values, 'a');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing policy value for a', testFailed, verbose);
testResults = [testResults; result];

% Wrong type for `pol.values.a`
ss_test = ss;
ss_test.pol.values.a = rand(2, 2, 2);
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.values.a should be a matrix', testFailed, verbose);
testResults = [testResults; result];

% Incompatible policy matrix row size (shocks × states)
ss_test = ss;
ss_test.pol.values.a = rand(99, size(ss.pol.values.a, 2));
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Wrong number of rows in policy matrix', testFailed, verbose);
testResults = [testResults; result];

% Incompatible policy matrix column size
ss_test = ss;
ss_test.pol.values.a = rand(size(ss.pol.values.a, 1), 99);
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Wrong number of columns in policy matrix', testFailed, verbose);
testResults = [testResults; result];

%% === DISTRIBUTION TESTS ===

% Missing `d` field
ss_test = ss;
ss_test = rmfield(ss_test, 'd');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing distribution structure `d`', testFailed, verbose);
testResults = [testResults; result];

% Wrong type for `ss.d.grids.a`
ss_test = ss;
ss_test.d.grids.a = rand(3, 3);  % not a vector
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'd.grids.a should be a vector', testFailed, verbose);
testResults = [testResults; result];

% Missing `d.hist`
ss_test = ss;
ss_test.d = rmfield(ss.d, 'hist');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing d.hist', testFailed, verbose);
testResults = [testResults; result];

% Wrong type for `ss.d.hist`
ss_test = ss;
ss_test.d.hist = rand(2,2,2);
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'd.hist should be a matrix', testFailed, verbose);
testResults = [testResults; result];

% Wrong row size in `d.hist`
ss_test = ss;
ss_test.d.hist = rand(99, size(ss.d.hist,2));
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Wrong number of rows in d.hist', testFailed, verbose);
testResults = [testResults; result];

% Wrong column size in `d.hist`
ss_test = ss;
ss_test.d.hist = rand(size(ss.d.hist,1), 99);
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Wrong number of columns in d.hist', testFailed, verbose);
testResults = [testResults; result];

%% === AGGREGATES TESTS ===

% Missing `agg` field
ss_test = ss;
ss_test = rmfield(ss_test, 'agg');
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing `agg` field', testFailed, verbose);
testResults = [testResults; result];

% Remove one required aggregate variable
ss_test = ss;
ss_test.agg = rmfield(ss.agg, 'r');  % remove interest rate
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'Missing aggregate variable r', testFailed, verbose);
testResults = [testResults; result];

% Wrong type for `ss.agg.r`
ss_test = ss;
ss_test.agg.r = 'non-numeric';
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'agg.r should be numeric', testFailed, verbose);
testResults = [testResults; result];

%% === PERMUTATION TESTS ===

% pol.order is not a string or cellstr
ss_test = ss;
ss_test.pol.order = 123;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.order is not a string array or cellstr', testFailed, verbose);
testResults = [testResults; result];

% pol.order has wrong names
ss_test = ss;
ss_test.pol.order = {'wrong_name'};
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.order fieldnames is not consistent', testFailed, verbose);
testResults = [testResults; result];

% d.order is not a string or cellstr
ss_test = ss;
ss_test.d.order = 3.14;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'd.order is not a string array or cellstr', testFailed, verbose);
testResults = [testResults; result];

% d.order inconsistent with pol.order
ss_test = ss;
ss_test.d.order = {'ghost_shock'};
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'd.order not matching pol.order', testFailed, verbose);
testResults = [testResults; result];

%% === NaN TESTS ===

% NaN in shocks.grids.e
ss_test = ss;
ss_test.shocks.grids.e(3) = NaN;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.grids.e contain NaN elements', testFailed, verbose);
testResults = [testResults; result];

% NaN in shocks.Pi.e
ss_test = ss;
ss_test.shocks.Pi.e(2,3) = NaN;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.Pi.e contains NaN elements', testFailed, verbose);
testResults = [testResults; result];

% NaN in pol.grids.a
ss_test = ss;
ss_test.pol.grids.a(10) = NaN;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.grids.a contain NaN elements', testFailed, verbose);
testResults = [testResults; result];

% NaN in pol.values.a
ss_test = ss;
ss_test.pol.values.a(5,10) = NaN;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.values.a contain NaN elements', testFailed, verbose);
testResults = [testResults; result];

% NaN in d.grids.a
ss_test = ss;
ss_test.d.grids.a(20) = NaN;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'd.grids.a contain NaN elements', testFailed, verbose);
testResults = [testResults; result];

% NaN in d.hist
ss_test = ss;
ss_test.d.hist(3,15) = NaN;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'd.hist contains NaN elements', testFailed, verbose);
testResults = [testResults; result];

% NaN in agg.r
ss_test = ss;
ss_test.agg.r = NaN;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'agg.r contain NaN values', testFailed, verbose);
testResults = [testResults; result];

%% === MONOTONICITY TESTS ===

% Non-monotonic shocks.grids.e (decreasing)
ss_test = ss;
ss_test.shocks.grids.e = [1; 2; 3; 2.5; 4; 5; 6];  % element 4 < element 3
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.grids.e are not strictly increasing', testFailed, verbose);
testResults = [testResults; result];

% Non-monotonic shocks.grids.e (equal consecutive values)
ss_test = ss;
ss_test.shocks.grids.e = [1; 2; 3; 3; 4; 5; 6];  % element 3 == element 4
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'shocks.grids.e are not strictly increasing (equal values)', testFailed, verbose);
testResults = [testResults; result];

% Non-monotonic pol.grids.a (decreasing)
ss_test = ss;
ss_test.pol.grids.a(50) = ss.pol.grids.a(49) - 0.1;  % Make element 50 < element 49
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.grids.a are not strictly increasing', testFailed, verbose);
testResults = [testResults; result];

% Non-monotonic pol.grids.a (equal consecutive values)
ss_test = ss;
ss_test.pol.grids.a(25) = ss.pol.grids.a(24);  % Make element 25 == element 24
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.grids.a are not strictly increasing (equal values)', testFailed, verbose);
testResults = [testResults; result];

%% === GRID BOUNDS TESTS ===

% Policy function below pol.grids minimum
ss_test = ss;
ss_test.pol.values.a(1,1) = min(ss.pol.grids.a) - 1.0;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.values.a below pol.grids minimum', testFailed, verbose);
testResults = [testResults; result];

% Policy function above pol.grids maximum
ss_test = ss;
ss_test.pol.values.a(end,end) = max(ss.pol.grids.a) + 1.0;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.values.a exceeds pol.grids maximum', testFailed, verbose);
testResults = [testResults; result];

% Policy function below d.grids minimum (when d.grids is specified)
ss_test = ss;
ss_test.d.grids.a = ss.pol.grids.a;
ss_test.pol.values.a(1,1) = min(ss.pol.grids.a) - 0.5;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.values.a below d.grids minimum', testFailed, verbose);
testResults = [testResults; result];

% Policy function above d.grids maximum (when d.grids is specified)
ss_test = ss;
ss_test.d.grids.a = ss.pol.grids.a;
ss_test.pol.values.a(end,end) = max(ss.pol.grids.a) + 0.5;
[testFailed, result] = expect_error(@() heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test), ...
    'pol.values.a exceeds d.grids maximum', testFailed, verbose);
testResults = [testResults; result];

%% === REDUNDANCY TESTS ===

% shocks.Pi contains redundant entry (not in shocks.grids)
ss_test = ss;
ss_test.shocks.Pi.extra = eye(length(ss.shocks.grids.e));  % not listed in shocks.grids
try
    [out_ss, sizes] = heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test);
    fprintf('✔ Redundant shocks.Pi entry accepted (warning expected)\n');
    testResults = [testResults; struct('name', 'Redundant shocks.Pi entry', 'passed', true, 'message', '')];
catch ME
    testFailed = testFailed+1;
    fprintf('❌ Unexpected error from redundant shocks.Pi entry: %s\n', ME.message);
    testResults = [testResults; struct('name', 'Redundant shocks.Pi entry', 'passed', false, 'message', ME.message)];
end

% pol.values contains redundant field (not in model's pol_symbs)
ss_test = ss;
ss_test.pol.values.redundant_var = ones(size(ss.pol.values.a));
try
    [out_ss, sizes] = heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test);
    fprintf('✔ Redundant pol.values field accepted (warning expected)\n');
    testResults = [testResults; struct('name', 'Redundant pol.values field', 'passed', true, 'message', '')];
catch ME
    testFailed = testFailed+1;
    fprintf('❌ Unexpected error from redundant pol.values field: %s\n', ME.message);
    testResults = [testResults; struct('name', 'Redundant pol.values field', 'passed', false, 'message', ME.message)];
end

% d.grids contains redundant field (not a declared state)
ss_test = ss;
ss_test.d.grids.not_a_state = [0;1;2];  % not in M_.heterogeneity.state_var
try
    [out_ss, sizes] = heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test);
    fprintf('✔ Redundant d.grids field accepted (warning expected)\n');
    testResults = [testResults; struct('name', 'Redundant d.grids field', 'passed', true, 'message', '')];
catch ME
    testFailed = testFailed+1;
    fprintf('❌ Unexpected error from redundant d.grids field: %s\n', ME.message);
    testResults = [testResults; struct('name', 'Redundant d.grids field', 'passed', false, 'message', ME.message)];
end

% agg contains extra field (not in M_.endo_names)
ss_test = ss;
ss_test.agg.extra_agg = 999;
try
    [out_ss, sizes] = heterogeneity.check_steady_state_input(M_, options_.heterogeneity, ss_test);
    fprintf('✔ Redundant agg field accepted (warning expected)\n');
    testResults = [testResults; struct('name', 'Redundant agg field', 'passed', true, 'message', '')];
catch ME
    testFailed = testFailed+1;
    fprintf('❌ Unexpected error from redundant agg field: %s\n', ME.message);
    testResults = [testResults; struct('name', 'Redundant agg field', 'passed', false, 'message', ME.message)];
end

end
