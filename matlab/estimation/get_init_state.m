function a = get_init_state(a,xparam1,estim_params_,dr,M_,options_)
% a = get_init_state(a,xparam1,estim_params_,dr,M_,options_)
% Computes the endogenous log prior addition to the initial prior
%
% INPUTS
%    a                  [double]     k vector of initial states
%    xparam1            [double]     n vector of estimated params
%    estim_params_      [structure]  characterizing parameters to be estimated
%    dr                 [structure]  decision rule structure
%    M_                 [structure]  Model description
%    options_           [structure]  MATLAB's structure describing the current options
%
% OUTPUTS
%    a                  [double]     k vector of updated initial states

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

M_ = set_all_parameters(xparam1,estim_params_,M_);

pvec=[];
for ii=1:size(M_.filter_initial_state,1)
    if ~isempty(M_.filter_initial_state{ii,1})
        tmp1 = strrep(M_.filter_initial_state{ii,2},');','');
        tmp1 = strrep(tmp1,'M_.params(','');
        pvec = [pvec eval(tmp1)];
    end
end

for ii=1:size(M_.filter_initial_state,1)
    if ~isempty(M_.filter_initial_state{ii,1})
        if options_.loglinear && ~options_.logged_steady_state
            eval(['a(ii) = log(' strrep(M_.filter_initial_state{ii,2},';','') ') - log(dr.ys(ii));']);
        elseif ~options_.loglinear && ~options_.logged_steady_state
            eval(['a(ii) = ' strrep(M_.filter_initial_state{ii,2},';','') '- dr.ys(ii);'])
        else
            error('The steady state is logged. This should not happen. Please contact the developers')
        end
    end
end

a=a(dr.order_var);
