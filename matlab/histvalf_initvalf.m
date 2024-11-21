function [series, p] = histvalf_initvalf(caller, M_, options)
% [series, p] = histvalf_initvalf(caller, M_, options)
% Handles options for histval_file and initval_file
%
% INPUTS
% - caller           [char]      row array, name of calling function
% - M_               [struct]    model description, a.k.a
% - options          [struct]    options specific to histval_file and initval_file
%
% OUTPUTS
% - series           [dseries]   selected data from a file or a dseries
% - p                [integer]   number of periods (excluding the initial and terminal conditions)

% Copyright © 2003-2024 Dynare Team
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


% dseries
if isfield(options, 'series')
    series = evalin('base', options.series);
    dseries_ispresent = true;
else
    dseries_ispresent = false;
end

% file
datafile = '';
if isfield(options, 'filename')
    warning('%s_FILE: option FILENAME is deprecated, please use option DATAFILE', caller)
    if dseries_ispresent
        error('%s_FILE: you can''t use option FILENAME and option SERIES at the same time', caller)
    end
    if isfield(options, 'datafile')
        error('%s_FILE: you can''t use option DATAFILE and option FILENAME at the same time', caller)
    end
    datafile = options.filename;
end

if isfield(options, 'datafile')
    if dseries_ispresent
        error('%s_FILE: you can''t use option DATAFILE and option SERIES at the same time', caller)
    end
    datafile = options.datafile;
end

if datafile
    [~,basename,extension] = fileparts(datafile);
    % Auto-detect extension if not provided
    if isempty(extension)
        if isfile([basename '.m'])
            extension = '.m';
        elseif isfile([basename '.mat'])
            extension = '.mat';
        elseif isfile([basename '.xls'])
            extension = '.xls';
        elseif isfile([basename '.xlsx'])
            extension = '.xlsx';
        else
            error('%s_FILE: Can''t find datafile: %s.{m,mat,xls,xlsx}', caller, basename);
        end
    end
    fullname = [basename extension];
    series = dseries(fullname);
end

% checking that all variable are present
error_flag = false;
for i = 1:M_.orig_endo_nbr
    if ~series.exist(M_.endo_names{i})
        dprintf('%s_FILE: endogenous variable %s is missing', caller, M_.endo_names{i})
        error_flag = true;
    end
end

for i = 1:M_.exo_nbr
    if ~series.exist(M_.exo_names{i})
        dprintf('%s_FILE: exogenous variable %s is missing', caller, M_.exo_names{i})
        error_flag = true;
    end
end

for i = 1:M_.exo_det_nbr
    if ~series.exist(M_.exo_det_names{i})
        dprintf('%s_FILE: exo_det variable %s is missing', caller, M_.exo_det_names{i})
        error_flag = true;
    end
end

if error_flag
    error('%s_FILE: some variables are missing', caller)
end

if isfile(sprintf('+%s/dynamic_set_auxiliary_series.m', M_.fname))
    series = feval(sprintf('%s.dynamic_set_auxiliary_series', M_.fname), series, M_.params);
end

% selecting observations
if isfield(options, 'nobs')
    nobs = options.nobs;
else
    nobs = 0;
end

periods = series.dates;
nobs0 = series.nobs;

first_obs_ispresent = false;
last_obs_ispresent = false;

if isfield(options, 'first_obs') && ~isempty(options.first_obs)
    if isa(options.first_obs, 'numeric')
        if options.first_obs < 1
            error('%s_FILE: first_obs must be a positive number', caller)
        elseif options.first_obs > nobs0
            error('%s_FILE: first_obs = %d is larger than the number of observations in the data file (%d)', ...
                  caller, options.first_obs, nobs0)
        else
            first_obs = periods(options.first_obs);
        end
    elseif isa(options.first_obs, 'dates')
        first_obs = options.first_obs;
    else
        error('Incorrect class for the value of first_obs option')
    end
    if isfield(options, 'first_simulation_period')
        if isa(options.first_simulation_period, 'numeric')
            if first_obs ~= periods(options.first_simulation_period) - M_.orig_maximum_lag
                error('%s_FILE: first_obs = %s and first_simulation_period = %d have values inconsistent with a maximum lag of %d periods', ...
                      caller, first_obs, options.first_simulation_period, M_.orig_maximum_lag)
            end
        elseif isa(options.first_simulation_period, 'dates')
            if first_obs ~= options.first_simulation_period - M_.orig_maximum_lag
                error('%s_FILE: first_obs = %s and first_simulation_period = %s have values inconsistent with a maximum lag of %d periods', ...
                      caller, first_obs, options.first_simulation_period, M_.orig_maximum_lag)
            end
        else
            error('Incorrect class for the value of first_simulation_period option')
        end
    end
    first_obs_ispresent = true;
else
    first_obs = periods(1);
end

if ~first_obs_ispresent && isfield(options, 'first_simulation_period')
    if isa(options.first_simulation_period, 'numeric')
        if options.first_simulation_period < M_.orig_maximum_lag
            error('%s_FILE: first_simulation_period = %d must be larger than the maximum lag (%d)', ...
                  caller, options.first_simulation_period, M_.orig_maximum_lag)
        elseif options.first_simulation_period > nobs0
            error('%s_FILE: first_simulations_period = %d is larger than the number of observations in the data file (%d)', ...
                  caller, options.first_obs, nobs0)
        end
        first_obs = periods(options.first_simulation_period) - M_.orig_maximum_lag;
        first_obs_ispresent = true;
    elseif isa(options.first_simulation_period, 'dates')
        first_obs = options.first_simulation_period - M_.orig_maximum_lag;
        first_obs_ispresent = true;
    else
        error('Incorrect class for the value of first_simulation_period option')
    end
end

if isfield(options, 'last_obs')
    if isa(options.last_obs, 'numeric')
        if options.last_obs > nobs0
            error('%s_FILE: last_obs = %d is larger than the number of observations in the dataset (%d)', caller, options.last_obs, nobs0)
        else
            last_obs = periods(options.last_obs);
        end
    elseif isa(options.last_obs, 'dates')
        if options.last_obs > series.last
            error('%s_FILE: last_obs = %s is larger than the number of observations in the dataset (%s)', caller, options.last_obs, series.last)
        else
            last_obs = options.last_obs;
        end
    end

    if first_obs_ispresent
        if nobs > 0 && (last_obs ~= first_obs + nobs - 1)
            error('%s_FILE: FIRST_OBS, LAST_OBS and NOBS contain inconsistent information. Use only two of these options.', caller)
        end
    else
        if nobs > 0
            first_obs = last_obs - nobs + 1;
        end
    end
elseif nobs > 0
    last_obs = first_obs + nobs - 1;
else
    last_obs = series.last;
end

if isfield(options, 'last_simulation_period')
    if isa(options.last_simulation_period, 'numeric')
        lastsimulationperiod = periods(options.last_simulation_period);
    elseif isa(options.last_simulation_period, 'dates')
        lastsimulationperiod = options.last_simulation_period;
    else
        error('Incorrect class for the value of last_simulation_period option')
    end
    if lastsimulationperiod <= last_obs-M_.orig_maximum_lead
        last_obs = lastsimulationperiod+M_.orig_maximum_lead;
    else
        error('%s_FILE: LAST_SIMULATION_PERIOD is too large compared to the available data.', caller)
    end
    p = lastsimulationperiod-(first_obs+M_.orig_maximum_lag)+1;
else
    p = (last_obs-M_.orig_maximum_lead)-(first_obs+M_.orig_maximum_lag)+1;
end

series = series(first_obs:last_obs);
