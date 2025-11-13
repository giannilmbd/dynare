function [index_init_state, IS, index_deep_parameters] = get_init_state_estim_params(M_, bayestopt_, dr)
% [index_init_state, IS, index_deep_parameters] = get_init_state_estim_params(M_, bayestopt_, dr)
% Computes the endogenous log prior addition to the initial prior
%
% INPUTS
% - M_                  [structure]   Matlab's structure describing the model
% - bayestopt_          [structure]   describing the priors
% - dr                  [structure]   describing the decision rules
%
% OUTPUTS
% - index_init_state         [integer]     indices associated to initial states in bayestopt_
% - index_deep_parameters [integer]   indices associated to deep parameters in bayestopt_ 
% 
% Copyright © 2024 Dynare Team
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

% varargin
% 1        2            3        4  5             6          7          8   9                  10                11                   12
% dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info

pvec=[];
for ii=1:size(M_.filter_initial_state,1)
    if ~isempty(M_.filter_initial_state{ii,1})
        tmp1 = strrep(M_.filter_initial_state{ii,2},');','');
        tmp1 = strrep(tmp1,'M_.params(','');
        pvec = [pvec eval(tmp1)];
    end
end
[~,~,index_init_state] = intersect(M_.param_names(pvec),bayestopt_.name,'stable');
index_init_state = sort(index_init_state);
index_deep_parameters = 1:length(bayestopt_.name);
index_deep_parameters = index_deep_parameters(not(ismember(bayestopt_.name,bayestopt_.name(index_init_state))));
% index_deep_parameters containts indices in estim params vector that are NOT init states

nam=M_.endo_names(dr.order_var);
nam=nam(dr.restrict_var_list(bayestopt_.mf0));
for k=1:length(nam)
    nam{k} = [nam{k} 'init']; 
end
[~,~,IS] = intersect(bayestopt_.name(index_init_state),nam,'stable');
% IS contains state indices corresponding estim params