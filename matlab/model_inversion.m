function [endogenousvariables, exogenousvariables] = model_inversion(constraints, ...
                                                  exogenousvariables, ...
                                                  initialconditions, M_, options_, oo_)
% function [endogenousvariables, exogenousvariables] = model_inversion(constraints, ...
%                                                   exogenousvariables, ...
%                                                   initialconditions, M_, options_, oo_)
% INPUTS
% - constraints         [dseries]        with N constrained endogenous variables from t1 to t2.
% - exogenousvariables  [dseries]        with Q exogenous variables.
% - initialconditions   [dseries]        with M endogenous variables starting before t1 (M initialcond must contain at least the state variables).
% - M_                  [struct]         Dynare global structure containing information related to the model.
% - options_            [struct]         Dynare global structure containing all the options.
% - oo_                 [struct]         Dynare global structure containing intermediate results.
%
% OUTPUTS
% - endogenousvariables          [dseries]
% - exogenousvariables           [dseries]
%
% REMARKS

% Copyright © 2018-2026 Dynare Team
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

if ~isequal(nargin, 6)
    error('model_inversion: This routine require six input arguments!')
end

if ~isdseries(constraints)
    error('model_inversion: First input argument must be a dseries object!')
end

if ~isdseries(exogenousvariables)
    error('model_inversion: Second input argument must be a dseries object!')
end

if ~isempty(initialconditions) && ~isdseries(initialconditions)
    error('model_inversion: Third input argument must be a dseries object!')
end

if ~isstruct(M_)
    error('model_inversion: Last input argument must be structures (M_)!')
end

% Set range where the endogenous variables are constrained.
crange = constraints.dates;

% Check that the number of instruments match the number of constrained endogenous variables.
instruments = exogenousvariables(crange);
freeinnovations = instruments.name(find(all(isnan(instruments))));
if ~isequal(length(freeinnovations), constraints.vobs)
    error('The number of instruments must be equal to the number of constrained variables!')
end

% Check if some of the exogenous variables are given.
observed_exogenous_variables_flag = false;
if exogenousvariables.vobs>constraints.vobs
    observed_exogenous_variables_flag = true;
end

if M_.maximum_lag
    % Add auxiliary variables in initialconditions object.
    initialconditions = checkdatabase(initialconditions, M_, true, false);
end

% Get the list of endogenous and exogenous variables.
endo_names = M_.endo_names;
exo_names = M_.exo_names;

exogenousvariables = exogenousvariables{exo_names{:}};

% Use specialized routine if the model is backward looking.
if ~M_.maximum_lead
    if M_.maximum_lag
        [endogenousvariables, exogenousvariables] = ...
            backward_model.inversion(constraints, exogenousvariables, initialconditions, ...
                                     endo_names, exo_names, freeinnovations, ...
                                     M_, options_, oo_);
    else
        [endogenousvariables, exogenousvariables] = ...
            static_model_inversion(constraints, exogenousvariables, ...
                                   endo_names, exo_names, freeinnovations, ...
                                   M_, options_, oo_);
    end
    return
end

% Perform a (constrained) extended path simulation

extra_periods = 25; % Taken from (now removed) det_cond_forecast.m
total_periods = crange.ndat + extra_periods;

endo_simul = repmat(oo_.steady_state, 1, M_.maximum_lag + total_periods + M_.maximum_lead);
exo_simul = repmat(oo_.exo_steady_state', M_.maximum_lag + total_periods + M_.maximum_lead, 1);

crange_with_init = (crange(1)-M_.maximum_lag):crange(end);

initialconditions = initialconditions(crange_with_init);
for i = 1:initialconditions.vobs
    iy = find(strcmp(initialconditions.name{i}, M_.endo_names));
    if ~isempty(iy)
        endo_simul(iy,1:initialconditions.nobs) = initialconditions.data(:,i);
    else
        ix = find(strcmp(initialconditions.name{i}, M_.exo_names));
        if ~isempty(ix) && M_.maximum_lag > 0
            exo_simul(1:M_.maximum_lag,ix) = initialconditions.data(1:M_.maximum_lag,i);
        end
    end
end

% Loop over simulation periods for extended path
for p = 1:crange.ndat
    options_.periods = total_periods - p + 1;

    options_.stack_solve_algo = 0;
    options_.linear_approximation = false;
    options_.endogenous_terminal_period = false;
    options_.bytecode = false;
    options_.block = false;
    options_.verbosity = 0;

    % Set the exogenous observed variables that come as a surprise in this period
    if observed_exogenous_variables_flag
        list_of_observed_exogenous_variables = setdiff(exo_names, freeinnovations);
        observed_exogenous_variables = exogenousvariables{list_of_observed_exogenous_variables{:}};
        for i = 1:length(list_of_observed_exogenous_variables)
            ix = find(strcmp(list_of_observed_exogenous_variables{i}, M_.exo_names));
            exo_simul(p+M_.maximum_lag,ix) = observed_exogenous_variables{list_of_observed_exogenous_variables{i}}(crange(p)).data;
        end
    end

    % Set constrained path for the endogenous variables that come as a surprise in this period
    controlled_paths_by_period = struct('exogenize_id', cell(options_.periods, 1), 'endogenize_id', cell(options_.periods, 1), 'values', cell(options_.periods, 1), 'learnt_in', cell(options_.periods, 1));
    for i = 1:constraints.vobs
        iy = find(strcmp(constraints.name{i}, M_.endo_names));
        ix = find(strcmp(freeinnovations{i}, M_.exo_names));
        controlled_paths_by_period(1).exogenize_id(end+1) = iy;
        controlled_paths_by_period(1).endogenize_id(end+1) = ix;
        controlled_paths_by_period(1).values(end+1) = constraints{constraints.name{i}}(crange(p)).data;
        controlled_paths_by_period(1).learnt_in(end+1) = 1;
    end

    % Perform the deterministic simulation for this informational period
    endo_simul_tmp = endo_simul(:,p:end);
    exo_simul_tmp = exo_simul(p:end,:);

    [endo_simul_tmp, success, ~, ~, ~, exo_simul_tmp] = perfect_foresight_solver_core(endo_simul_tmp, exo_simul_tmp, oo_.steady_state, oo_.exo_steady_state, controlled_paths_by_period, M_, options_);
    if ~success
        error('model_inversion: inversion failure in period %d', p)
    end

    endo_simul(:,p:end) = endo_simul_tmp;
    exo_simul(p:end,:) = exo_simul_tmp;
end

endogenousvariables = dseries(endo_simul(:,1:crange_with_init.ndat)', crange_with_init, M_.endo_names);
exogenousvariables = dseries(exo_simul(1:crange_with_init.ndat,:), crange_with_init, M_.exo_names);
