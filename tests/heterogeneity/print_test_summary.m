function print_test_summary(testResults)
    % Print a comprehensive test summary

    if isempty(testResults)
        fprintf('\nNo test results to display.\n');
        return;
    end

    totalTests = length(testResults);
    passedTests = sum([testResults.passed]);
    failedTests = totalTests - passedTests;

    fprintf('\n%s\n', repmat('=', 1, 50));
    fprintf('TEST SUMMARY\n');
    fprintf('%s\n', repmat('=', 1, 50));
    fprintf('Total tests: %d\n', totalTests);
    fprintf('Passed: %d (%.1f%%)\n', passedTests, 100*passedTests/totalTests);
    fprintf('Failed: %d (%.1f%%)\n', failedTests, 100*failedTests/totalTests);

    if failedTests > 0
        fprintf('\nFAILED TESTS:\n');
        fprintf('%s\n', repmat('-', 1, 30));
        for i = 1:totalTests
            if ~testResults(i).passed
                fprintf('• %s\n', testResults(i).name);
                if ~isempty(testResults(i).message)
                    fprintf('  ↪ %s\n', testResults(i).message);
                end
            end
        end
    end

    fprintf('%s\n\n', repmat('=', 1, 50));
end