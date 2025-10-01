function y = get_filtered_time_series(y, m, options_)
% y = get_filtered_time_series(y, m, options_)
% Computes filtered time series
% INPUTS
%   y                   [double]       nvar*nperiods vector of simulated variables.
%   m                   [double]       nvar mean to be subtracted in case of no filtering 
%   options_            [structure]    Dynare's options structure
%
% OUTPUTS
%   y                   [double]       nvar*nperiods vector of filtered variables.

% Copyright © 2001-2025 Dynare Team
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


if options_.hp_filter && ~options_.one_sided_hp_filter  && ~options_.bandpass.indicator
    [~,y] = sample_hp_filter(y,options_.hp_filter);
elseif ~options_.hp_filter && options_.one_sided_hp_filter && ~options_.bandpass.indicator
    [~,y] = one_sided_hp_filter(y,options_.one_sided_hp_filter);
elseif ~options_.hp_filter && ~options_.one_sided_hp_filter && options_.bandpass.indicator
    data_temp=dseries(y,'0q1');
    data_temp=baxter_king_filter(data_temp,options_.bandpass.passband(1),options_.bandpass.passband(2),options_.bandpass.K);
    y=[NaN(options_.bandpass.K,size(y,2)); data_temp.data; NaN(options_.bandpass.K,size(y,2))]; %pad cut off observations with NaN to preserve length of time series
elseif ~options_.hp_filter && ~options_.one_sided_hp_filter  && ~options_.bandpass.indicator
    y = bsxfun(@minus, y, m);
else
    error('disp_moments:: You cannot use more than one filter at the same time')
end