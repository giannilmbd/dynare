function [index_init_state, IS, index_deep_parameters] = get_init_state_estim_params(M_, bayestopt_, dr)
% [index_init_state, IS, index_deep_parameters] = get_init_state_estim_params(M_, bayestopt_, dr)
% Locates initial-state parameters and maps them to state indices.
%
% This helper identifies entries in `bayestopt_.name` that correspond to
% endogenous initial-state parameters (those prefixed with 'init '). It
% returns their positions within the estimated-parameters vector
% (`index_init_state`), the complementary set of deep parameters
% (`index_deep_parameters`), and the mapping `IS` from those init-state
% parameter names to the state indices in decision-rule order (restricted
% state vector governed by `dr.restrict_var_list`).
%
% INPUTS
% - M_                  [structure]   model structure; uses `endo_names`.
% - bayestopt_          [structure]   priors/meta; uses `name`, `mf0`.
% - dr                  [structure]   decision rules; uses `order_var`, `restrict_var_list`.
%
% OUTPUTS
% - index_init_state     [integer]   indices of init-state parameters in `bayestopt_.name`.
% - IS                   [integer]   state indices (decision-rule order) matching init parameters.
% - index_deep_parameters[integer]   indices of non-init (deep) parameters in `bayestopt_.name`.
%
% NOTES
% - The name convention is 'init <state_name>' for init-state parameters.
% - `IS` is computed via intersection of init names and the restricted
%   state-name list constructed in decision-rule order for states.
% 
% Copyright © 2024-2026 Dynare Team
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

tmp=startsWith(bayestopt_.name, 'init ');
index_init_state = find(tmp);
index_deep_parameters = find(~tmp);
% index_deep_parameters contains indices in estim params vector that are NOT init states

nam=M_.endo_names(dr.order_var);
nam=nam(dr.restrict_var_list(bayestopt_.mf0));
for k=1:length(nam)
    nam{k} = ['init ' nam{k}]; 
end
[~,~,IS] = intersect(bayestopt_.name(index_init_state),nam,'stable');
% IS contains state indices in decision rule order corresponding estim params