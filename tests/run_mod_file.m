% Copyright © 2011-2026 Dynare Team
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

source_root = getenv('source_root');
addpath([source_root filesep 'tests' filesep 'utils']);
addpath([source_root filesep 'matlab']);
profile_suite_env = getenv('DYNARE_PROFILE_SUITE');
profile_suite = strcmpi(profile_suite_env, 'true') || strcmpi(profile_suite_env, '1') || strcmpi(profile_suite_env, 'yes');

if isoctave
    load_octave_packages
end

fprintf('\n*** TESTING: %s ***\n\n', getenv('mod_file'));

tic;

% NB: all variables will be cleared by the call to Dynare
try
    if ~isoctave && profile_suite
        profile clear;
        profile on;
    end
    % Read arguments from individual environment variables
    dynare_arg_count = str2double(getenv('dynare_arg_count'));
    if isnan(dynare_arg_count) || dynare_arg_count == 0
        dynare(getenv('mod_file'), 'console')
    else
        args_cell = cell(1, dynare_arg_count);
        for i = 1:dynare_arg_count
            args_cell{i} = getenv(sprintf('dynare_arg_%d', i-1));
        end
        dynare(getenv('mod_file'), 'console', args_cell{:})
    end
    testFailed = false;
catch exception
    printTestError(getenv('mod_file'), exception);
    testFailed = true;
end

if ~isoctave && profile_suite
    profile off;
    pinfo = profile('info');
    fprintf('\n*** Profiling summary (top 20 by total time) ***\n');
    if ~isempty(pinfo.FunctionTable)
        [~, idx] = sort([pinfo.FunctionTable.TotalTime], 'descend');
        topN = min(20, numel(idx));
        for k = 1:topN
            ft = pinfo.FunctionTable(idx(k));
            child_time = 0;
            if isfield(ft, 'Children') && ~isempty(ft.Children)
                child_time = sum([ft.Children.TotalTime]);
            end
            self_time = ft.TotalTime - child_time;
            if isfield(ft, 'NumCalls')
                calls = ft.NumCalls;
            elseif isfield(ft, 'TotalCalls')
                calls = ft.TotalCalls;
            else
                calls = NaN;
            end
            fprintf('%2d) %-60s total %8.3f s | self %8.3f s (%d calls)\n', k, ft.FunctionName, ft.TotalTime, self_time, calls);
        end
    else
        fprintf('(no profile data)\n');
    end
end

fprintf('\n*** Elapsed time (in seconds): %.1f\n\n', toc);

% Ensure proper termination and exit code
quit(testFailed)
