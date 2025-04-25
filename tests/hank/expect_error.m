function [incTestFailed] = expect_error(fn, description, testFailed, show_message)
    if nargin < 3
        show_message = false;
    end
    try
        fn();
        incTestFailed = testFailed+1;
        fprintf('❌ Expected error for: %s\n', description);
    catch ME
        fprintf('✔ Correctly threw error for: %s\n', description);
        if show_message
            fprintf('   ↪ Error message: %s\n', ME.message);
        end
        incTestFailed = testFailed;
    end
end