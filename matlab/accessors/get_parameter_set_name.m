function parameter_set_string=get_parameter_set_name(parameter_set)
% parameter_set_string=get_parameter_set_name(options_)
% Inputs:
%   parameter_set           [string]    content of options_.parameter_set
% 
% Outputs:
%   parameter_set_string    [string]    associated name

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

switch parameter_set
    case 'calibration'
        parameter_set_string='Calibration';
    case 'prior_mode'
        parameter_set_string='Prior Mode';
    case 'prior_mean'
        parameter_set_string='Prior Mean';
    case 'posterior_mode'
        parameter_set_string='Posterior Mode';
    case 'posterior_mean'
        parameter_set_string='Posterior Mean';
    case 'posterior_median'
        parameter_set_string='Posterior Median';
    case 'mle_mode'
        parameter_set_string='MLE Mode';
    otherwise
        parameter_set_string='';
end