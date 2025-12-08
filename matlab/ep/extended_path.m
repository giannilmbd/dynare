function [ts,oo_,errorflag] = extended_path(initialconditions, samplesize, exogenousvariables, options_, M_, oo_)

% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
% series of size T  is obtained by solving T perfect foresight models.
%
% INPUTS
% - initialconditions      [double]    m*1 array, where m is the number of endogenous variables in the model.
% - samplesize             [integer]   scalar, size of the sample to be simulated.
% - exogenousvariables     [double]    T*n array, values for the structural innovations.
% - options_               [struct]    options_
% - M_                     [struct]    Dynare's model structure
% - oo_                    [struct]    Dynare's results structure
%
% OUTPUTS
% - ts                     [dseries]   m*samplesize array, the simulations.
% - oo_                    [struct]
% - errorflag              [logical]   scalar, true if the nonlinear solver for the auxiliary model failed in some period.
%
% REMARKS
% If errorflag==true, because the nonlinear solver failed in period T<samplesize, ts holds the simulations for periods 1 to T-1.

% Copyright © 2009-2025 Dynare Team
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

errorflag = false;

if ~isempty(M_.perfect_foresight_controlled_paths)
    error('extended_path command is not compatible with perfect_foresight_controlled_paths block')
end

[initialconditions, innovations, pfm, options_, oo_] = ...
    extended_path_initialization(initialconditions, samplesize, exogenousvariables, options_, M_, oo_);


% Set the initial period.
if isdates(options_.initial_period)
    if ischar(options_.initial_period)
        initial_period = dates(options_.initial_period);
    else
        initial_period = options_.initial_period;
    end
elseif isnan(options_.initial_period)
    initial_period = dates(1,1);
else
    error('Type of option initial_period is wrong.')
end

[shocks, spfm_exo_simul, innovations, oo_] = extended_path_shocks(innovations, exogenousvariables, samplesize, M_, options_, oo_);

% Initialize the matrix for the paths of the endogenous variables.
endogenous_variables_paths = NaN(M_.endo_nbr, samplesize+1);
endogenous_variables_paths(:,1) = initialconditions;

% Set waitbar (graphic or text  mode)
[hh_fig, length_of_old_string] = wait_bar.run(0, [], 'Please wait. Extended Path simulations...', options_.console_mode, 0, 'EP simulations.');

% Initialize while-loop index.
t = 1;

% Main loop.
while (t <= samplesize)
    if ~mod(t,10)
        [~,length_of_old_string]=wait_bar.run(t/samplesize,hh_fig,'Please wait. Extended Path simulations...',options_.console_mode,length_of_old_string);
    end
    % Set period index.
    t = t+1;
    spfm_exo_simul(2,:) = shocks(t-1,:);
    if t>2
        % Set initial guess for the solver (using the solution of the
        % previous period problem).
        initialguess = [endogenousvariablespaths(:, 2:end), oo_.steady_state];
    else
        initialguess = [];
    end
    if t>2
        [endogenous_variables_paths(:,t), info_convergence, endogenousvariablespaths, y] = extended_path_core(innovations.positive_var_indx, ...
                                                                                                           spfm_exo_simul, ...
                                                                                                           endogenous_variables_paths(:,t-1), ...
                                                                                                           pfm, ...
                                                                                                           M_, ...
                                                                                                           options_, ...
                                                                                                           oo_, ...
                                                                                                           initialguess, ...
                                                                                                           y);
    else
        [endogenous_variables_paths(:,t), info_convergence, endogenousvariablespaths, y, pfm, options_] = extended_path_core(innovations.positive_var_indx, ...
                                                                                                                          spfm_exo_simul, ...
                                                                                                                          endogenous_variables_paths(:,t-1), ...
                                                                                                                          pfm, ...
                                                                                                                          M_, ...
                                                                                                                          options_, ...
                                                                                                                          oo_, ...
                                                                                                                          initialguess, ...
                                                                                                                          []);
    end
    if ~info_convergence
        msg = sprintf('No convergence of the (stochastic) perfect foresight solver (in period %s)!', int2str(t));
        warning(msg)
        errorflag = true;
        break
    end
end % (while) loop over t

% Close waitbar.
wait_bar.close(hh_fig,options_.console_mode);

% Return the simulated time series.
if any(isnan(endogenous_variables_paths(:)))
    sl = find(~isnan(endogenous_variables_paths));
    nn = size(endogenous_variables_paths, 1);
    endogenous_variables_paths = reshape(endogenous_variables_paths(sl), nn, length(sl)/nn);
end
ts = dseries(transpose(endogenous_variables_paths), initial_period, M_.endo_names);

oo_.endo_simul = transpose(ts.data);
