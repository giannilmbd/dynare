function oo_=perfect_foresight_setup(M_, options_, oo_)
% Prepares a deterministic simulation, by filling oo_.exo_simul and oo_.endo_simul
%
% INPUTS
%   M_                  [structure] describing the model
%   options_            [structure] describing the options
%   oo_                 [structure] storing the results
%
% OUTPUTS
%   oo_                 [structure] storing the results
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS
%   none

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

test_for_deep_parameters_calibration(M_);

if size(M_.lead_lag_incidence,2)-nnz(M_.lead_lag_incidence(M_.maximum_endo_lag+1,:)) > 0
    mess = 'PERFECT_FORESIGHT_SETUP: error in model specification : the variable(s) ';
    var_list = M_.endo_names(M_.lead_lag_incidence(M_.maximum_endo_lag+1,:)==0);
    for i=1:length(var_list)
        if i<length(var_list)
            mess = [mess, var_list{i} ', '];
        else
            mess = [mess, var_list{i} ];
        end
    end
    mess = [mess ' don''t appear as current period variables.'];
    error(mess)
end

[periods, ~, last_simulation_period] = get_simulation_periods(options_);

if ~isempty(M_.det_shocks)
    % Check whether some expected shocks happen after the terminal period.
    mess = '';
    for i=1:length(M_.det_shocks)
        if isa(M_.det_shocks(i).periods, 'dates') && isempty(last_simulation_period)
            error('PERFECT_FORESIGHT_SETUP: temporary shocks are specified using dates but neither first_simulation_period nor last_simulation_period option was passed')
        end
        if (isa(M_.det_shocks(i).periods, 'numeric') && any(M_.det_shocks(i).periods > periods)) ...
               || (isa(M_.det_shocks(i).periods, 'dates') && any(M_.det_shocks(i).periods > last_simulation_period))
            mess = sprintf('%s\n   At least one expected value for %s has been declared after the terminal period.', mess, M_.exo_names{M_.det_shocks(i).exo_id});
        end
    end
    if ~isempty(mess)
        disp(sprintf('\nPERFECT_FORESIGHT_SETUP: Problem with the declaration of the expected shocks:\n%s', mess));
        skipline()
        error('PERFECT_FORESIGHT_SETUP: Please check the declaration of the shocks or increase the number of periods in the simulation.')
    end
end

if options_.simul.endval_steady && M_.maximum_lead == 0
    error('PERFECT_FORESIGHT_SETUP: Option endval_steady cannot be used on a purely backward or static model.')
end

oo_ = make_ex_(M_,options_,oo_);
oo_ = make_y_(M_,options_,oo_);
