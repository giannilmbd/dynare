function M_=update_parameters_filter_initial_state(M_,alphahat01,ys,options_)
% Updates M_.params to reflect new parameter values
%
% INPUTS
%    M_                 [structure]  Model description
%    a                  [double]     k vector of initial states
%    ys                 [structure]  steady states 
%    options_           [structure]  MATLAB's structure describing the current options
%
% OUTPUTS
%    M_                 [double]     results structure

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


for ii=1:size(M_.filter_initial_state,1)
    if ~isempty(M_.filter_initial_state{ii,1})
        if options_.loglinear && ~options_.logged_steady_state
            eval([strrep(M_.filter_initial_state{ii,2},';','') ' = exp(log(ys(ii))+alphahat01(ii));']);
        elseif ~options_.loglinear && ~options_.logged_steady_state
            eval([strrep(M_.filter_initial_state{ii,2},';','') '= ys(ii)+alphahat01(ii);'])
        else
            error('The steady state is logged. This should not happen. Please contact the developers')
        end
    end
end
