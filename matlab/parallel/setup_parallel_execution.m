function [run_with_pct, cleanup_obj] = setup_parallel_execution(use_pct, caller_name)
% [run_with_pct, cleanup_obj] = setup_parallel_execution(use_pct, caller_name)
% ------------------------------------------------------------------------------
% This function handles parallel pool management for Dynare commands that
% can optionally use the Parallel Computing Toolbox:
%   - Checks if PCT is available (requires MATLAB R2024b+ and PCT installed, not Octave)
%   - Respects user preferences (use_pct option)
%   - If user wants serial but pool is open, closes pool and sets up restoration
%   - Returns appropriate settings for the calling function
%
% INPUTS
%   use_pct      [boolean or empty]  User's preference for parallel execution:
%                                      true:  explicitly request parallel execution
%                                      false: explicitly request serial execution
%                                      []:    use default (parallel if available)
%   caller_name  [string]            Name of calling function (for warnings)
%
% OUTPUTS
%   run_with_pct [boolean]           true if parallel execution should be used
%   cleanup_obj  [onCleanup or []]   Cleanup object that restores pool state
%                                    when done; keep in scope until finished.
%                                    Closes pool if created, restores pool
%                                    if closed for serial execution.
%
% CASE MATRIX (see unit tests at the end of this file)
%   has_pct=true (PCT available):
%     Case 1: use_pct=[], pool closed  -> creates pool, run_with_pct=true, cleanup closes pool
%     Case 2: use_pct=[], pool open    -> uses pool, run_with_pct=true
%     Case 3: use_pct=true, pool closed -> creates pool, run_with_pct=true, cleanup closes pool
%     Case 4: use_pct=true, pool open  -> uses pool, run_with_pct=true
%     Case 5: use_pct=false, pool closed -> run_with_pct=false, no cleanup
%     Case 6: use_pct=false, pool open -> run_with_pct=false, closes pool, cleanup restores it
%   has_pct=false (PCT not available, e.g. Octave or MATLAB<R2024b):
%     Case 7: use_pct=[]    -> run_with_pct=false, no warning
%     Case 8: use_pct=true  -> run_with_pct=false, warning issued (specific to reason)
%     Case 9: use_pct=false -> run_with_pct=false, no warning
%
% EXAMPLES
%   % in posterior_sampler_core.m:
%   [run_with_pct, cleanup] = setup_parallel_execution(options_.parallel_info.use_pct.estimation.sampler, 'posterior_sampler_core');
%   if run_with_pct
%       parfor i = 1:n
%           ...
%       end
%   else
%       for i = 1:n
%           ...
%       end
%   end
%   % Pool is restored automatically when cleanup goes out of scope (including on error)

% Copyright © 2026 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

% Check if Parallel Computing Toolbox is available and supported
has_pct = ~isoctave && ~matlab_ver_less_than('24.2') && matlab.internal.parallel.isPCTInstalled && license('test', 'Distrib_Computing_Toolbox');

% Initialize outputs
run_with_pct = false;
cleanup_obj = [];

% Check if a parallel pool is currently open
pool_open = false;
if has_pct
    pool_open = ~isempty(gcp('nocreate'));
end

% Determine desired behavior based on use_pct setting:
%   1. use_pct is empty []: default behavior - use parallel if available
%   2. use_pct is true:    user explicitly wants parallel
%   3. use_pct is false:   user explicitly wants serial
if isempty(use_pct)
    % User didn't specify; default to using PCT if available
    want_parallel = has_pct;
elseif use_pct
    % User explicitly wants parallel
    want_parallel = true;
else
    % User explicitly wants serial (use_pct = false)
    want_parallel = false;
end

% Configure parallel execution based on desire and availability
if want_parallel
    if has_pct
        try
            pool = gcp; % Get existing pool or create new one
            if ~isempty(pool)
                run_with_pct = true;
                if ~pool_open
                    % Pool was not open before; set up cleanup to close it
                    % when the cleanup object goes out of scope
                    cleanup_obj = onCleanup(@() delete(gcp('nocreate')));
                end
            end
        catch ME
            warning_id = sprintf('%s:pct_failure', caller_name);
            warning(warning_id, ...
                'Could not access parallel pool despite ''use_pct=true''. Falling back to serial. Error: %s', ...
                ME.message);
            run_with_pct = false;
        end
    else
        % User wants parallel but PCT not available
        if ~isempty(use_pct) && use_pct
            % Only warn if user explicitly requested parallel;
            % provide specific message based on reason
            if isoctave
                warning_id = sprintf('%s:pct_not_available', caller_name);
                warning(warning_id, 'Parallel Computing Toolbox is not supported in Octave. Falling back to serial.');
            elseif matlab_ver_less_than('24.2') && license('test', 'Distrib_Computing_Toolbox')
                warning_id = sprintf('%s:matlab_version_too_old', caller_name);
                warning(warning_id, 'Parallel Computing Toolbox is installed but requires MATLAB R2024b or later. Falling back to serial.');
            else
                warning_id = sprintf('%s:pct_not_installed', caller_name);
                warning(warning_id, 'Parallel Computing Toolbox is not installed. Falling back to serial.');
            end
        end
        run_with_pct = false;
    end
else
    % User wants serial execution
    run_with_pct = false;
    if has_pct && pool_open
        % Pool is open but user wants serial; close it and set up restoration
        pool_obj = gcp('nocreate');
        pool_profile = pool_obj.Cluster.Profile;
        pool_workers = pool_obj.NumWorkers;
        delete(pool_obj);
        % Create cleanup object to restore pool when it goes out of scope
        cleanup_obj = onCleanup(@() restore_parallel_pool(pool_profile, pool_workers, caller_name));
    end
end

function restore_parallel_pool(profile, num_workers, caller_name)
% Restore parallel pool with original settings
try
    parpool(profile, num_workers);
catch ME
    warning_id = sprintf('%s:pool_restore_failed', caller_name);
    warning(warning_id, 'Failed to restore parallel pool: %s', ME.message);
end
end


return % --*-- Unit tests --*--

% Unit tests for setup_parallel_execution
% Tests are split based on whether PCT is available or not.
% When PCT is available (MATLAB R2024b+ with PCT installed): Cases 1-6
% When PCT is NOT available (Octave or old MATLAB or no PCT): Cases 7-9

%@test:1
% Case 1: use_pct = [] (default), pool NOT open
% Expected: creates pool, run_with_pct = true, cleanup closes pool
% Only runs if PCT is available
if isoctave || matlab_ver_less_than('24.2') || ~matlab.internal.parallel.isPCTInstalled %#ok<UNRCH>
    t(1) = true; % skip test
else
    try
        delete(gcp('nocreate')); % ensure no pool
        [run_with_pct, cleanup] = setup_parallel_execution([], 'unittest1');
        t(1) = run_with_pct == true;
        t(2) = ~isempty(cleanup); % cleanup should exist to close pool
        t(3) = ~isempty(gcp('nocreate')); % pool should exist
        % Trigger cleanup (closes pool) by clearing the variable
        clear cleanup;
        t(4) = isempty(gcp('nocreate')); % pool should be closed
    catch
        t = false(4,1);
    end
end
T = all(t);
%@eof:1

%@test:2
% Case 2: use_pct = [] (default), pool already open
% Expected: uses existing pool, run_with_pct = true
% Only runs if PCT is available
if isoctave || matlab_ver_less_than('24.2') || ~matlab.internal.parallel.isPCTInstalled
    t(1) = true; % skip test
else
    try
        delete(gcp('nocreate'));
        parpool('local', 2);
        pool_before_NumWorkers = gcp('nocreate').NumWorkers;
        [run_with_pct, cleanup] = setup_parallel_execution([], 'unittest2');
        t(1) = run_with_pct == true;
        t(2) = isempty(cleanup); % no cleanup needed
        t(3) = ~isempty(gcp('nocreate')); % pool should still exist
        t(4) = gcp('nocreate').NumWorkers == pool_before_NumWorkers; % same pool
        delete(gcp('nocreate')); % cleanup
    catch
        t = false(4,1);
    end
end
T = all(t);
%@eof:2

%@test:3
% Case 3: use_pct = true, pool NOT open
% Expected: creates pool, run_with_pct = true, cleanup closes pool
% Only runs if PCT is available
if isoctave || matlab_ver_less_than('24.2') || ~matlab.internal.parallel.isPCTInstalled
    t(1) = true; % skip test
else
    try
        delete(gcp('nocreate')); % ensure no pool
        [run_with_pct, cleanup] = setup_parallel_execution(true, 'unittest3');
        t(1) = run_with_pct == true;
        t(2) = ~isempty(cleanup); % cleanup should exist to close pool
        t(3) = ~isempty(gcp('nocreate')); % pool should exist
        % Trigger cleanup (closes pool) by clearing the variable
        clear cleanup;
        t(4) = isempty(gcp('nocreate')); % pool should be closed
    catch
        t = false(4,1);
    end
end
T = all(t);
%@eof:3

%@test:4
% Case 4: use_pct = true, pool already open
% Expected: uses existing pool, run_with_pct = true
% Only runs if PCT is available
if isoctave || matlab_ver_less_than('24.2') || ~matlab.internal.parallel.isPCTInstalled
    t(1) = true; % skip test
else
    try
        delete(gcp('nocreate'));
        parpool('local', 2);
        pool_before_NumWorkers = gcp('nocreate').NumWorkers;
        [run_with_pct, cleanup] = setup_parallel_execution(true, 'unittest4');
        t(1) = run_with_pct == true;
        t(2) = isempty(cleanup); % no cleanup needed
        t(3) = ~isempty(gcp('nocreate')); % pool should still exist
        t(4) = gcp('nocreate').NumWorkers == pool_before_NumWorkers; % same pool
        delete(gcp('nocreate')); % cleanup
    catch
        t = false(4,1);
    end
end
T = all(t);
%@eof:4

%@test:5
% Case 5: use_pct = false, pool NOT open
% Expected: run_with_pct = false, no cleanup needed
% Only runs if PCT is available
if isoctave || matlab_ver_less_than('24.2') || ~matlab.internal.parallel.isPCTInstalled
    t(1) = true; % skip test
else
    try
        delete(gcp('nocreate')); % ensure no pool
        [run_with_pct, cleanup] = setup_parallel_execution(false, 'unittest5');
        t(1) = run_with_pct == false;
        t(2) = isempty(cleanup); % no cleanup needed
        t(3) = isempty(gcp('nocreate')); % no pool should exist
    catch
        t = false(3,1);
    end
end
T = all(t);
%@eof:5

%@test:6
% Case 6: use_pct = false, pool open
% Expected: run_with_pct = false, pool closed, cleanup restores it
% Only runs if PCT is available
if isoctave || matlab_ver_less_than('24.2') || ~matlab.internal.parallel.isPCTInstalled
    t(1) = true; % skip test
else
    try
        delete(gcp('nocreate'));
        parpool('local', 2);
        pool_before_NumWorkers = gcp('nocreate').NumWorkers;
        pool_before_Profile = gcp('nocreate').Cluster.Profile;
        [run_with_pct, cleanup] = setup_parallel_execution(false, 'unittest6');
        t(1) = run_with_pct == false;
        t(2) = ~isempty(cleanup); % cleanup object should exist
        t(3) = isempty(gcp('nocreate')); % pool should be closed
        % Trigger cleanup (i.e. restoring pool) by clearing the variable
        clear cleanup;
        % Verify pool was restored and is the same to before
        pool_after = gcp('nocreate');
        t(4) = ~isempty(pool_after);
        t(5) = pool_after.NumWorkers == pool_before_NumWorkers;
        t(6) = strcmp(pool_after.Cluster.Profile, pool_before_Profile);
        delete(gcp('nocreate')); % final cleanup
    catch
        t = false(6,1);
    end
end
T = all(t);
%@eof:6

%@test:7
% Case 7: use_pct = [] (default), PCT not available
% Expected: run_with_pct = false, no warning
% Only runs if PCT is NOT available
if ~isoctave && ~matlab_ver_less_than('24.2') && matlab.internal.parallel.isPCTInstalled
    t(1) = true; % skip test - PCT is available
else
    try
        lastwarn(''); % clear last warning
        [run_with_pct, cleanup] = setup_parallel_execution([], 'unittest7');
        t(1) = run_with_pct == false;
        t(2) = isempty(cleanup); % no cleanup needed
        [warnMsg, ~] = lastwarn();
        t(3) = isempty(warnMsg); % no warning expected
    catch
        t = false(3,1);
    end
end
T = all(t);
%@eof:7

%@test:8
% Case 8: use_pct = true, PCT not available
% Expected: run_with_pct = false, warning issued
% Only runs if PCT is NOT available
if ~isoctave && ~matlab_ver_less_than('24.2') && matlab.internal.parallel.isPCTInstalled
    t(1) = true; % skip test - PCT is available
else
    try
        lastwarn(''); % clear last warning
        [run_with_pct, cleanup] = setup_parallel_execution(true, 'unittest8');
        t(1) = run_with_pct == false;
        t(2) = isempty(cleanup); % no cleanup needed
        [warnMsg, warnId] = lastwarn();
        t(3) = ~isempty(warnMsg); % warning expected
        t(4) = strcmp(warnId, 'unittest8:pct_not_available') || ...
               strcmp(warnId, 'unittest8:pct_not_installed') || ...
               strcmp(warnId, 'unittest8:matlab_version_too_old');
    catch
        t = false(4,1);
    end
end
T = all(t);
%@eof:8

%@test:9
% Case 9: use_pct = false, PCT not available
% Expected: run_with_pct = false, no warning
% Only runs if PCT is NOT available
if ~isoctave && ~matlab_ver_less_than('24.2') && matlab.internal.parallel.isPCTInstalled
    t(1) = true; % skip test - PCT is available
else
    try
        lastwarn(''); % clear last warning
        [run_with_pct, cleanup] = setup_parallel_execution(false, 'unittest9');
        t(1) = run_with_pct == false;
        t(2) = isempty(cleanup); % no cleanup needed
        [warnMsg, ~] = lastwarn();
        t(3) = isempty(warnMsg); % no warning expected
    catch
        t = false(3,1);
    end
end
T = all(t);
%@eof:9

end
