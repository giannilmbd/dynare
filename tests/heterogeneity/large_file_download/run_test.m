function [testFailedOut, testResult] = run_test(testName, testFn, testFailedIn)
    % Safe test runner
    % Usage: [testFailed, result] = run_test('Test Name', @() your_test_code(), testFailed);

    try
        testFn();
        testResult = struct('name', testName, 'passed', true, 'message', '');
        fprintf('✔ %s\n', testName);
        testFailedOut = testFailedIn;
    catch ME
        testResult = struct('name', testName, 'passed', false, 'message', ME.message);
        fprintf('❌ %s: %s\n', testName, ME.message);
        testFailedOut = testFailedIn + 1;
    end
end