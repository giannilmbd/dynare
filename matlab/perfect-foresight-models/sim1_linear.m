function [endogenousvariables, success, ERR, exogenousvariables] = sim1_linear(endogenousvariables, exogenousvariables, steadystate_y, steadystate_x, controlled_paths_by_period, M_, options_)
% [endogenousvariables, success, ERR, exogenousvariables] = sim1_linear(endogenousvariables, exogenousvariables, steadystate_y, steadystate_x, M_, options_)
% Solves a linear approximation of a perfect foresight model using sparse matrix.
%
% INPUTS
% - endogenousvariables [double] N*T array, paths for the endogenous variables (initial condition + initial guess + terminal condition).
% - exogenousvariables  [double] T*M array, paths for the exogenous variables.
% - steadystate_y       [double] N*1 array, steady state for the endogenous variables.
% - steadystate_x       [double] M*1 array, steady state for the exogenous variables.
% - controlled_paths_by_period [struct] data from perfect_foresight_controlled_paths block
% - M_                  [struct] contains a description of the model.
% - options_            [struct] contains various options.
%
% OUTPUTS
% - endogenousvariables [double] N*T array, paths for the endogenous variables (solution of the perfect foresight model).
% - success             [logical] Whether a solution was found
% - ERR                 [double] ∞-norm of the residual
% - exogenousvariables  [double] T*M array, paths for the exogenous variables
%                                (may be modified if perfect_foresight_controlled_paths present)
%
% NOTATIONS
% - N is the number of endogenous variables.
% - M is the number of innovations.
% - T is the number of periods (including initial and/or terminal conditions).
%
% REMARKS
% - The structure `M_` describing the structure of the model, must contain the
% following information:
%  + lead_lag_incidence, incidence matrix (given by the preprocessor).
%  + endo_nbr, number of endogenous variables (including aux. variables).
%  + exo_nbr, number of innovations.
%  + maximum_lag,
%  + maximum_endo_lag,
%  + params, values of model's parameters.
%  + fname, name of the model.
%  + NNZDerivatives, number of non zero elements in the Jacobian of the dynamic model.
% - The structure `options_`, must contain the following options:
%  + verbosity, controls the quantity of information displayed.
%  + periods, the number of periods in the perfect foresight model.
%  + debug.
% - The steady state of the exogenous variables is required because we need
% to center the variables around the deterministic steady state to solve the
% perfect foresight model.

% Copyright © 2015-2025 Dynare Team
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

verbose = options_.verbosity;

ny = M_.endo_nbr;

maximum_lag = M_.maximum_lag;

periods = get_simulation_periods(options_);

params = M_.params;

if ~isempty(controlled_paths_by_period)
    for p = 1:periods
        if isempty(controlled_paths_by_period(p).exogenize_id)
            continue
        end
        endogenousvariables(controlled_paths_by_period(p).exogenize_id, p+M_.maximum_lag) = controlled_paths_by_period(p).values;
    end

    if options_.debug
        error('Debugging not available with perfect_foresight_controlled_paths')
    end
end

if verbose
    skipline()
    printline(80)
    disp('MODEL SIMULATION:')
    skipline()
end

% Evaluate the Jacobian of the dynamic model at the deterministic steady state.
y3n = repmat(steadystate_y, 3, 1);
[d1, TT_order, TT] = feval([M_.fname,'.dynamic_resid'], y3n, steadystate_x', params, ...
                           steadystate_y);
jacobian = feval([M_.fname,'.dynamic_g1'], y3n, steadystate_x', params, steadystate_y, ...
                 M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, ...
                 M_.dynamic_g1_sparse_colptr, TT_order, TT);

if options_.debug
    column=find(all(jacobian(:, 1:3*ny)==0,1));
    if ~isempty(column)
        fprintf('The dynamic Jacobian is singular. The problem derives from:\n')
        for iter=1:length(column)
            var_index = mod(column(iter)-1, ny)+1;
            lead_lag = floor((column(iter)-1)/ny)-1;
            if lead_lag == 0
                fprintf('The derivative with respect to %s being 0 for all equations.\n',M_.endo_names{var_index})
            elseif lead_lag == -1
                fprintf('The derivative with respect to the lag of %s being 0 for all equations.\n',M_.endo_names{var_index})
            elseif lead_lag == 1
                fprintf('The derivative with respect to the lead of %s being 0 for all equations.\n',M_.endo_names{var_index})
            end            
        end
    end   
end

% Check that the dynamic model was evaluated at the steady state.
if ~options_.steadystate.nocheck &&  max(abs(d1))>options_.solve_tolf
    error('Jacobian is not evaluated at the steady state!')
end

% Center the endogenous and exogenous variables around the deterministic steady state.
endogenousvariables = bsxfun(@minus, endogenousvariables, steadystate_y);
exogenousvariables = bsxfun(@minus, exogenousvariables, transpose(steadystate_x));

if M_.maximum_lag > 0
    y0 = endogenousvariables(:, maximum_lag);
else
    y0 = NaN(M_.endo_nbr, 1);
end
if M_.maximum_lead > 0
    yT = endogenousvariables(:, maximum_lag+periods+1);
else
    yT = NaN(M_.endo_nbr, 1);
end
Y = vec(endogenousvariables(:,maximum_lag+(1:periods)));

h2 = clock;

[res, A] = linear_perfect_foresight_problem(Y, jacobian, y0, yT, exogenousvariables, params, steadystate_y, maximum_lag, periods, ny);

if ~isempty(controlled_paths_by_period)
    A = controlled_paths_substitute_stacked_jacobian(A, repmat(steadystate_y, periods, 1), steadystate_y, steadystate_y, repmat(steadystate_x', periods, 1), steadystate_y, controlled_paths_by_period, M_);
end

% Evaluation of the maximum residual at the initial guess (steady state for the endogenous variables).
err = norm(res, 'Inf'); % Do not use max(max(abs(…))) because it omits NaN

if options_.debug
    fprintf('\nLargest absolute residual at iteration %d: %10.3f\n', 1, err);
    if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y))
        fprintf('\nWARNING: NaN or Inf detected in the residuals or endogenous variables.\n');
    end
    if ~isreal(res) || ~isreal(Y)
        fprintf('\nWARNING: Imaginary parts detected in the residuals or endogenous variables.\n');
    end
    skipline()
end

% Try to update the vector of endogenous variables.
try
    dY = -lin_solve(A, res, options_);
    Y = Y + dY;
catch
    % Normally, because the model is linear, the solution of the perfect foresight model should
    % be obtained in one Newton step. This is not the case if the model is singular.
    success = false;
    ERR = [];
    if verbose
        skipline()
        disp('Singularity problem! The jacobian matrix of the stacked model cannot be inverted.')
    end
    return
end

if ~isempty(controlled_paths_by_period)
    for p = 1:periods
        endogenize_id = controlled_paths_by_period(p).endogenize_id;
        exogenize_id = controlled_paths_by_period(p).exogenize_id;
        if isempty(endogenize_id)
            continue
        end
        Y(exogenize_id+(p-1)*M_.endo_nbr) = controlled_paths_by_period(p).values - steadystate_y(exogenize_id);
        exogenousvariables(p+M_.maximum_lag,endogenize_id) = exogenousvariables(p+M_.maximum_lag,endogenize_id) + dY(exogenize_id+(p-1)*M_.endo_nbr)';
    end
end

YY = reshape([y0; Y; yT], ny, periods+2);

i_rows = (1:ny)';
for t = 1:periods
    z = [ vec(YY(:,t+(0:2))); transpose(exogenousvariables(t+maximum_lag, :))];
    res(i_rows) = jacobian*z;
    i_rows = i_rows + ny;
end

ERR = norm(res, 'Inf'); % Do not use max(max(abs(…))) because it omits NaN

if verbose
    fprintf('Iter: %s,\t Initial err. = %s,\t err. = %s,\t time = %s\n', num2str(1), num2str(err), num2str(ERR), num2str(etime(clock,h2)));
    printline(80);
end

if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y)) || ~isreal(res) || ~isreal(Y)
    success = false; % NaN or Inf occurred
    if verbose
        skipline()
        if ~isreal(res) || ~isreal(Y)
            disp('Simulation terminated with imaginary parts in the residuals or endogenous variables.')
        else
            disp('Simulation terminated with NaN or Inf in the residuals or endogenous variables.')
        end
        disp('There is most likely something wrong with your model. Try model_diagnostics or another simulation method.')
    end
else
    success = true; % Convergence obtained.
end

endogenousvariables(:,maximum_lag+(1:periods)) = reshape(Y, ny, periods);
endogenousvariables = bsxfun(@plus, endogenousvariables, steadystate_y);
exogenousvariables = bsxfun(@plus, exogenousvariables, steadystate_x');

if verbose
    skipline();
end
