function [xparam1, hh] = check_mode_file(xparam1, hh, options_, bayestopt_)
% function [xparam1, hh] = check_mode_file(xparam1, hh, options_, bayestopt_)
% -------------------------------------------------------------------------
% Check that the provided mode_file is compatible with the current estimation settings.
% -------------------------------------------------------------------------
% INPUTS
%  o xparam1:                [vector] current vector of parameter values at the mode
%  o hh:                     [matrix] current 'Hessian matrix at the mode
%  o options_:               [structure] information about options
%  o bayestopt_:             [structure] information about priors
% -------------------------------------------------------------------------
% OUTPUTS
%  o xparam1:                [vector] updated vector of parameter values at the mode
%  o hh:                     [matrix] updated Hessian matrix at the mode
% -------------------------------------------------------------------------
% This function is called by
%  o dynare_estimation_init.m
% -------------------------------------------------------------------------

% Copyright © 2023-2025 Dynare Team
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

number_of_estimated_parameters = length(xparam1);
mode_file = load(options_.mode_file);
if number_of_estimated_parameters>length(mode_file.xparam1)
    % More estimated parameters than parameters in the mode_file.
    skipline()
    disp(['The ''mode_file'' ' options_.mode_file ' has been generated using another specification of the model or another model!'])
    disp(['Your file contains estimates for ' int2str(length(mode_file.xparam1)) ' parameters, while you are attempting to estimate ' int2str(number_of_estimated_parameters) ' parameters:'])
    md = []; xd = [];
    for i=1:number_of_estimated_parameters
        id = match_parameter_name_with_stderr(bayestopt_.name{i,:},mode_file.parameter_names);
        if isempty(id)
            disp(['--> Estimated parameter ' bayestopt_.name{i} ' is not present in the loaded ''mode_file'' (prior mean or initialized values will be used, if possible).'])
        else
            xd = [xd; i];
            md = [md; id];
        end
    end
    for i=1:length(mode_file.xparam1)
        id = match_parameter_name_with_stderr(mode_file.parameter_names{i},bayestopt_.name);
        if isempty(id)
            disp(['--> Parameter ' mode_file.parameter_names{i} ' is not estimated according to the current mod file.'])
        end
    end
    if ~options_.mode_compute
        % The mode is not estimated.
        error('Please change the ''mode_file'' option, the list of estimated parameters or set ''mode_compute''>0.')
    else
        % The mode is estimated, the Hessian evaluated at the mode is not needed so we set values for the parameters missing in the mode file using the prior mean or initialized values.
        if ~isempty(xd)
            xparam1(xd) = mode_file.xparam1(md);
        else
            error('Please remove the ''mode_file'' option.')
        end
    end
elseif number_of_estimated_parameters<length(mode_file.xparam1)
    % Less estimated parameters than parameters in the mode_file.
    skipline()
    disp(['The ''mode_file'' ' options_.mode_file ' has been generated using another specification of the model or another model!'])
    disp(['Your file contains estimates for ' int2str(length(mode_file.xparam1)) ' parameters, while you are attempting to estimate only ' int2str(number_of_estimated_parameters) ' parameters:'])
    md = []; xd = [];
    for i=1:number_of_estimated_parameters
        id = match_parameter_name_with_stderr(bayestopt_.name{i,:},mode_file.parameter_names);
        if isempty(id)
            disp(['--> Estimated parameter ' bayestopt_.name{i,:} ' is not present in the loaded ''mode_file'' (prior mean or initialized values will be used, if possible).'])
        else
            xd = [xd; i];
            md = [md; id];
        end
    end
    for i=1:length(mode_file.xparam1)
        id = match_parameter_name_with_stderr(mode_file.parameter_names{i},bayestopt_.name);
        if isempty(id)
            disp(['--> Parameter ' mode_file.parameter_names{i} ' is not estimated according to the current mod file.'])
        end
    end
    if ~options_.mode_compute
        % The posterior mode is not estimated. If possible, fix the mode_file.
        if isequal(length(xd),number_of_estimated_parameters)
            disp('==> Fix mode_file (remove unused parameters).')
            xparam1 = mode_file.xparam1(md);
            if isfield(mode_file,'hh')
                hh = mode_file.hh(md,md);
            end
        else
            error('Please change the ''mode_file'' option, the list of estimated parameters or set ''mode_compute''>0.')
        end
    else
        % The mode is estimated, the Hessian evaluated at the mode is not needed so we set values for the parameters missing in the mode_file using the prior mean or initialized values.
        if ~isempty(xd)
            xparam1(xd) = mode_file.xparam1(md);
        else
            % None of the estimated parameters are present in the mode_file.
            error('Please remove the ''mode_file'' option.')
        end
    end
else
    % The number of declared estimated parameters match the number of parameters in the mode file.
    % Check that the parameters in the mode file and according to the current mod file are identical.
    if ~isfield(mode_file,'parameter_names')
        disp(['The ''mode_file'' ' options_.mode_file ' has been generated using an older version of Dynare. It cannot be verified if it matches the present model. Proceed at your own risk.'])
        mode_file.parameter_names=deblank(bayestopt_.name); %set names
    end
    if isequal(mode_file.parameter_names, bayestopt_.name)
        xparam1 = mode_file.xparam1;
        if isfield(mode_file,'hh')
            hh = mode_file.hh;
        end
    else
        skipline()
        disp(['The ''mode_file'' ' options_.mode_file ' has been generated using another specification of the model or another model!'])
        % Check if this is only an ordering issue or if the missing parameters can be initialized with the prior mean.
        md = []; xd = [];
        for i=1:number_of_estimated_parameters
            id = match_parameter_name_with_stderr(bayestopt_.name{i,:}, mode_file.parameter_names);
            if isempty(id)
                disp(['--> Estimated parameter ' bayestopt_.name{i} ' is not present in the loaded ''mode_file''.'])
            else
                xd = [xd; i];
                md = [md; id];
            end
        end
        if ~options_.mode_compute
            % The mode is not estimated
            if isequal(length(xd), number_of_estimated_parameters)
                % This is an ordering issue.
                xparam1 = mode_file.xparam1(md);
                if isfield(mode_file,'hh')
                    hh = mode_file.hh(md,md);
                end
            else
                error('Please change the ''mode_file'' option, the list of estimated parameters or set ''mode_compute''>0.')
            end
        else
            % The mode is estimated, the Hessian evaluated at the mode is not needed so we set values for the parameters missing in the mode_file using the prior mean or initialized values.
            if ~isempty(xd)
                xparam1(xd) = mode_file.xparam1(md);
                if isfield(mode_file,'hh')
                    hh(xd,xd) = mode_file.hh(md,md);
                end
            else
                % None of the estimated parameters are present in the mode_file.
                error('Please remove the ''mode_file'' option.')
            end
        end
    end
end
skipline()

%--------------------------------------------------------------------------
function id = match_parameter_name_with_stderr(search_name, search_list)
% MATCH_PARAMETER_NAME_WITH_STDERR Find a parameter name, accounting for stderr prefixes
%
% This function handles backward compatibility with mode files that don't
% have the "stderr " prefix for standard deviation parameters. Old mode files
% store parameters without the "stderr " prefix, while new code adds this prefix.
%
% INPUTS
%   search_name:  [string or char] The parameter name to search for (may have "stderr " prefix)
%   search_list:  [cell or char array] The list of parameter names to search in (no prefix)
%
% OUTPUTS
%   id: [int or empty] Index of the match if found, empty otherwise

% Ensure search_name is a char array 
if isstring(search_name)
    search_name = char(search_name);
end

% Try exact match first
id = strmatch(search_name, search_list, 'exact');
if ~isempty(id)
    return
end

% If no exact match and search_name has "stderr " prefix, try without it
if length(search_name) > 7 && strncmp(search_name, 'stderr ', 7)
    name_without_stderr = deblank(search_name(8:end));
    id = strmatch(name_without_stderr, search_list, 'exact');
end