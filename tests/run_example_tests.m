% This script runs all .mod files found in the examples directory.
% It dynamically discovers mod files, making it flexible to changes
% in the examples folder during development.

% Copyright © 2025 Dynare Team
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

source_dir = getenv('source_root');
addpath([source_dir filesep 'tests' filesep 'utils']);
addpath([source_dir filesep 'matlab']);

if isoctave
    load_octave_packages
end

dynare_config();

% Get the examples directory
examples_dir = [source_dir filesep 'examples'];

% Find all .mod files in the examples directory
mod_files = dir([examples_dir filesep '*.mod']);

if isempty(mod_files)
    error('No .mod files found in examples directory: %s', examples_dir);
end

failedtests = {};

fprintf('\n*** Running example mod files tests ***\n');
fprintf('*** Found %d mod files in %s ***\n\n', length(mod_files), examples_dir);

total_time = tic;

for i = 1:length(mod_files)
    mod_file = mod_files(i).name;
    
    fprintf('\n==========================================================\n');
    fprintf('*** TESTING (%d/%d): %s ***\n', i, length(mod_files), mod_file);
    fprintf('==========================================================\n\n');

    test_time = tic;

    % Change to examples directory so that helper files can be found
    old_dir = pwd;
    cd(examples_dir);

    try
        old_path = path;
        % Save workspace variables that need to persist across dynare calls
        save('wsMat.mat', 'mod_file', 'failedtests', 'mod_files', ...
             'i', 'examples_dir', 'old_dir', 'old_path', 'source_dir', 'total_time', 'test_time');

        dynare(mod_file, 'console');
        close all;

        load('wsMat.mat');
        path(old_path);
        cd(old_dir);

        fprintf('\n*** %s: PASSED (%.1f seconds) ***\n', mod_file, toc(test_time));
    catch exception
        load('wsMat.mat');
        path(old_path);
        cd(old_dir);

        failedtests{end+1} = mod_file;
        printTestError(mod_file, exception);
        clear exception;
    end
end

% Clean up temporary file
ws_mat_file = [examples_dir filesep 'wsMat.mat'];
if exist(ws_mat_file, 'file')
    delete(ws_mat_file);
end

fprintf('\n==========================================================\n');
fprintf('*** Example tests summary ***\n');
fprintf('==========================================================\n');
fprintf('Total tests: %d\n', length(mod_files));
fprintf('Passed: %d\n', length(mod_files) - length(failedtests));
fprintf('Failed: %d\n', length(failedtests));
fprintf('Total elapsed time (in seconds): %.1f\n', toc(total_time));

if ~isempty(failedtests)
    fprintf('\n*** Failed tests:\n');
    for i = 1:length(failedtests)
        fprintf('    - %s\n', failedtests{i});
    end
end

fprintf('\n');

quit(~isempty(failedtests))

