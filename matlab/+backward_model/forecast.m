function forecasts = forecast(M_,options_,oo_,initialcondition, listofvariables, periods, withuncertainty)
% forecasts = forecast(M_,options_,oo_,initialcondition, listofvariables, periods, withuncertainty)
% Returns unconditional forecasts.
%
% INPUTS
% - M_                  [struct]              Dynare's model structure
% - options_            [struct]              options_
% - oo_                 [struct]              Dynare's results structure
% - initialcondition    [dseries]             Initial conditions for the endogenous variables.
% - periods             [integer]             scalar, the number of (forecast) periods.
% - withuncertainty     [logical]             scalar, returns confidence bands if true.
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

% Check that the model is actually backward
if M_.maximum_lead
    error('backward_model_forecast:: The specified model is not backward looking!')
end

% Initialize returned argument.
forecasts = struct();

% Set defaults.
if nargin<5
    listofvariables = M_.endo_names;
    periods = 8;
    withuncertainty = false;
end

if nargin<6
    periods = 8;
    withuncertainty = false;
end

if nargin<7
    withuncertainty = false;
end

start = initialcondition.dates(end)+1;

% Set default initial conditions for the innovations.
for i=1:M_.exo_nbr
    if ~ismember(M_.exo_names{i}, initialcondition.name)
        initialcondition{M_.exo_names{i}} = dseries(zeros(initialcondition.nobs, 1), initialcondition.dates(1), M_.exo_names{i});
    end
end

% Set up initial conditions
[initialcondition, periods, innovations, options_, M_, oo_, endonames, ~, dynamic_resid, dynamic_g1] = ...
    backward_model.initialize(initialcondition, periods, options_, M_, oo_, zeros(periods, M_.exo_nbr));

% Get vector of indices for the selected endogenous variables.
n = length(listofvariables);
idy = zeros(n,1);
for i=1:n
    j = find(strcmp(listofvariables{i}, endonames));
    if isempty(j)
        error('backward_model_forecast:: Variable %s is unknown!', listofvariables{i})
    else
        idy(i) = j;
    end
end

% Set the number of simulations (if required).
if withuncertainty
    B = 1000;
end

% Get the covariance matrix of the shocks.
if withuncertainty
    sigma = get_lower_cholesky_covariance(M_.Sigma_e,options_.add_tiny_number_to_cholesky);
end

% Compute forecast without shock
if options_.linear
    [ysim__0, ~, ~, errorflag] = backward_model.simul_linear_model(initialcondition, periods, options_, M_, oo_, innovations, dynamic_resid, dynamic_g1);
else
    [ysim__0, ~, ~, errorflag] = backward_model.simul_nonlinear_model(initialcondition, periods, options_, M_, oo_, innovations, dynamic_resid, dynamic_g1);
end

if errorflag
    error('Simulation failed.')
end

forecasts.pointforecast = dseries(transpose(ysim__0(idy,:)), initialcondition.init, listofvariables);

% Set first period of forecast
forecasts.start = start;

if withuncertainty
    % Preallocate an array gathering the simulations.
    ArrayOfForecasts = zeros(n, periods+initialcondition.nobs, B);
    for i=1:B
        innovations = transpose(sigma*randn(M_.exo_nbr, periods));
        if options_.linear
            [ysim__, ~, ~, errorflag] = backward_model.simul_linear_model(initialcondition, periods, options_, M_, oo_, innovations, dynamic_resid, dynamic_g1);
        else
            [ysim__, ~, ~, errorflag] = backward_model.simul_nonlinear_model(initialcondition, periods, options_, M_, oo_, innovations, dynamic_resid, dynamic_g1);
        end
        if errorflag
            error('Simulation failed.')
        end
        ArrayOfForecasts(:,:,i) = ysim__(idy,:);
    end
    % Compute mean (over future uncertainty) forecast.
    forecasts.meanforecast = dseries(transpose(mean(ArrayOfForecasts, 3)), initialcondition.init, listofvariables);
    forecasts.medianforecast = dseries(transpose(median(ArrayOfForecasts, 3)), initialcondition.init, listofvariables);
    forecasts.stdforecast = dseries(transpose(std(ArrayOfForecasts, 1,3)), initialcondition.init, listofvariables);
    % Compute lower and upper 95% confidence bands
    ArrayOfForecasts = sort(ArrayOfForecasts, 3);
    forecasts.lb = dseries(transpose(ArrayOfForecasts(:,:,round(0.025*B))), initialcondition.init, listofvariables);
    forecasts.ub = dseries(transpose(ArrayOfForecasts(:,:,round(0.975*B))), initialcondition.init, listofvariables);
end
