function [incTestFailed, testResult] = expect_error(fn, description, testFailed, show_message)
    if nargin < 4
        show_message = false;
    end
    try
        fn();
        incTestFailed = testFailed+1;
        fprintf('❌ Expected error for: %s\n', description);
        testResult = struct('name', description, 'passed', false, 'message', 'Expected error but none was thrown');
    catch ME
        fprintf('✔ Correctly threw error for: %s\n', description);
        if show_message
            fprintf('   ↪ Error message: %s\n', ME.message);
        end
        incTestFailed = testFailed;
        testResult = struct('name', description, 'passed', true, 'message', '');
    end
end