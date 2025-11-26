function mh_autocorrelation_function(options_,M_,estim_params_,type,blck,name1,name2,name3)
% mh_autocorrelation_function(options_,M_,estim_params_,type,blck,name1,name2,name3)
% This function plots the autocorrelation of the sampled draws in the
% posterior distribution.
%
% INPUTS
%   options_        [structure]    Dynare structure.
%   M_              [structure]    Dynare structure (related to model definition).
%   estim_params_   [structure]    Dynare structure (related to estimation).
%   type            [string]       'DeepParameter', 'MeasurementError' (for measurement equation error) or 'StructuralShock' (for structural shock).
%   blck            [integer]      Number of the mh chain.
%   name1           [string]       Object name.
%   name2           [string]       Object name.
%   name3           [string]       Object name.
%
% OUTPUTS
%   None
%
% SPECIAL REQUIREMENTS

% Copyright © 2003-2025 Dynare Team
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

% Get the column index:
if nargin<7
    column = name2index(M_, estim_params_, type, name1);
elseif nargin<8
    column = name2index(M_, estim_params_, type, name1, name2);
elseif nargin<9
    column = name2index(M_, estim_params_, type, name1, name2, name3);
end

if isempty(column)
    return
end

% Get all the posterior draws:
PosteriorDraws = GetAllPosteriorDraws(options_, M_.dname, M_.fname, column, [],[], blck);

% Compute the autocorrelation function:
[~,autocor] = sample_autocovariance(PosteriorDraws,options_.mh_autocorrelation_function_size);

% Plot the posterior draws:

if strcmpi(type,'DeepParameter')
    TYPE = 'parameter ';
elseif strcmpi(type,'StructuralShock')
    if nargin<7
        TYPE = 'the standard deviation of structural shock ';
    elseif nargin<8
        TYPE = 'the correlation between structural shocks ';
    elseif nargin<9
        if strcmp(name1,name2) && strcmp(name1,name3)
            TYPE = 'the skewness coefficient of structural shock ';
        else
            TYPE = 'the skewness coefficient between structural shocks ';
        end
    end
elseif strcmpi(type,'MeasurementError')
    if nargin<7
        TYPE = 'the standard deviation of measurement error ';
    else
        TYPE = 'the correlation between measurement errors ';
    end
end

if nargin<7
    FigureName = ['Autocorrelogram for ' TYPE name1];
elseif nargin<8
    FigureName = ['Autocorrelogram for ' TYPE name1 ' and ' name2];
elseif nargin<9
    if strcmp(name1,name2) && strcmp(name1,name3)
        FigureName = ['Autocorrelogram for ' TYPE name1];
    else
        FigureName = ['Autocorrelogram for ' TYPE name1 ', ' name2 ', and ' name3];
    end
end

if options_.mh_nblck>1
    FigureName = [ FigureName , ' (block number' int2str(blck)  ').'];
end

hh_fig=dyn_figure(options_.nodisplay,'Name',FigureName);

bar(0:options_.mh_autocorrelation_function_size,autocor,'k');
axis tight
% create subdirectory <dname>/graphs if it doesn't exist
if ~isfolder(M_.dname)
    mkdir('.',M_.dname);
end
if ~isfolder([M_.dname filesep 'graphs'])
    mkdir(M_.dname,'graphs');
end

plot_name=get_the_name(column,0,M_,estim_params_,options_.varobs);
dyn_saveas(hh_fig,[M_.dname, filesep, 'graphs', filesep, 'MH_Autocorrelation_' plot_name],options_.nodisplay,options_.graph_format);
