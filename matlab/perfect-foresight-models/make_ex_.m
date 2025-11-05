function oo_ = make_ex_(M_, options_, oo_)

% Forms oo_.exo_simul and oo_.exo_det_simul
%
% INPUTS
% - M_           [struct]   Dynare model structure
% - options_     [struct]   Dynare options structure
% - oo_          [struct]   Dynare results structure
%
% OUTPUTS
% - oo_          [struct]   Updated Dynare results structure

% Copyright © 1996-2024 Dynare Team
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

try
    [periods, first_simulation_period] = get_simulation_periods(options_);
catch ME
    if strcmp(ME.identifier, 'Dynare:periodsNotSet')
        % This function is called from forecasts.run in some contexts where periods is not set
        periods = 0;
        first_simulation_period = options_.simul.first_simulation_period;
    else
        rethrow(ME);
    end
end

if isempty(oo_.exo_steady_state)
    oo_.exo_steady_state = zeros(M_.exo_nbr,1);
end
if M_.exo_det_nbr > 1 && isempty(oo_.exo_det_steady_state)
    oo_.exo_det_steady_state = zeros(M_.exo_det_nbr,1);
end

% Initialize oo_.exo_simul
if isempty(oo_.initval_series)
    if isempty(M_.exo_histval)
        if isempty(oo_.initial_exo_steady_state)
            oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag+periods+M_.maximum_lead,1);
        else
            oo_.exo_simul = [ repmat(oo_.initial_exo_steady_state',M_.maximum_lag,1) ; repmat(oo_.exo_steady_state',periods+M_.maximum_lead,1) ];
        end
    else
        if isempty(oo_.initial_exo_steady_state)
            oo_.exo_simul = [M_.exo_histval'; repmat(oo_.exo_steady_state',periods+M_.maximum_lead,1)];
        else
            error('histval and endval cannot be used simultaneously')
        end
    end
else
    if M_.exo_nbr > 0
        x = oo_.initval_series{M_.exo_names{:}}.data;
        oo_.exo_simul = x(M_.orig_maximum_lag-M_.maximum_lag+1:M_.orig_maximum_lag + periods + M_.maximum_lead,:);
        if ~isempty(M_.exo_histval)
            oo_.exo_simul(1:M_.maximum_lag, :) ...
                = M_.exo_histval(:, 1:M_.maximum_lag)';
        end
    else
        oo_.exo_simul=zeros(M_.maximum_lag + periods + M_.maximum_lead,M_.exo_nbr);
    end
    if M_.exo_det_nbr > 0
        x_det = oo_.initval_series{M_.exo_det_names{:}}.data;
        oo_.exo_det_simul = x_det(M_.orig_maximum_lag-M_.maximum_lag+1:M_.orig_maximum_lag + periods + M_.maximum_lead,:);
        if ~isempty(M_.exo_det_histval)
            oo_.exo_det_simul(1:M_.maximum_lag, :) ...
                = M_.exo_det_histval(:, 1:M_.maximum_lag)';
        end
    end
end
% Initialize oo_.exo_det_simul
if M_.exo_det_nbr > 0
    if isempty(M_.exo_det_histval)
        oo_.exo_det_simul = repmat(oo_.exo_det_steady_state',M_.maximum_lag+periods+M_.maximum_lead,1);
    else
        oo_.exo_det_simul = [M_.exo_det_histval'; repmat(oo_.exo_det_steady_state',periods+M_.maximum_lead,1)];
    end
end

% Add temporary shocks
if isfield(M_, 'det_shocks')
    for i = 1:length(M_.det_shocks)
        if isa(M_.det_shocks(i).periods, 'numeric')
            k = M_.det_shocks(i).periods + M_.maximum_lag;
            if any(k == 0)
                error('make_ex_: The model does not have lags, so you cannot set values for period 0');
            end
        else % dates
            k = M_.det_shocks(i).periods - first_simulation_period + 1 + M_.maximum_lag;
            if any(k <= 1)
                %% Contrary to the numeric case, we don’t allow k=0 for the dates case, because it
                %% was a design mistake to allow them here (histval should rather be used)
                error('make_ex_: some shocks are set at a date before the first simulation period');
            end
        end
        ivar = M_.det_shocks(i).exo_id;
        v = M_.det_shocks(i).value;
        if ~M_.det_shocks(i).exo_det
            switch M_.det_shocks(i).type
                case 'level'
                    oo_.exo_simul(k,ivar) = v;
                case 'multiply_steady_state'
                    oo_.exo_simul(k,ivar) = oo_.exo_steady_state(ivar) * v;
                case 'multiply_initial_steady_state'
                    if isempty(oo_.initial_exo_steady_state)
                        error('Option relative_to_initval of mshocks block cannot be used without an endval block')
                    end
                    oo_.exo_simul(k,ivar) = oo_.initial_exo_steady_state(ivar) * v;
            end
        else
            switch M_.det_shocks(i).type
                case 'level'
                    oo_.exo_det_simul(k,ivar) = v;
                case 'multiply_steady_state'
                    oo_.exo_det_simul(k,ivar) = oo_.exo_det_steady_state(ivar) * v;
                case 'multiply_initial_steady_state'
                    error('Option relative_to_initval of mshocks block cannot be used with a deterministic exogenous variable')
            end
        end
    end
end
