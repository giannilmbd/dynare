function forecasts = backward_model_forecast(varargin)
% forecasts = backward_model_forecast(varargin)
% wrapper for unconditional forecasts
%
% INPUTS
% INPUTS
% - varargin            inputs for backward_model_irf.forecast
%
% OUTPUTS
% - forecast            [structure]           each field is a dseries object for the point forecast (i.e.
%                                             forecasts without innovations in the future), the mean forecast,
%                                             the median forecast, the standard deviation of the predictive distribution
%                                             and the lower/upper bounds of the interval containing 95% of the predictive distribution.

% Copyright © 2017-2025 Dynare Team
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

global M_ options_ oo_;
forecasts = backward_model.forecast(M_,options_,oo_,varargin{:});