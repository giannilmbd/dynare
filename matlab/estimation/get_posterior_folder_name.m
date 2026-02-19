function folder_name=get_posterior_folder_name(options_)
% function folder_name=get_posterior_folder_name(options_)
% runs the estimation of the model
%
% INPUTS
%   options_            [structure] describing the options
%
% OUTPUTS
%   folder_name         [char]      output folder for different posterior
%                                   samplers
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2024-2025 Dynare Team
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

if ~issmc(options_)
    folder_name='metropolis';
else
    if ishssmc(options_)
        folder_name='hssmc';
    elseif isdsmh(options_)
        folder_name='dsmh';
    elseif isdime(options_)
        folder_name='dime';
    elseif isonline(options_)
        folder_name='online';
    else
        error('get_posterior_folder_name:: case should not happen. Please contact the developers')
    end
end