function e = euler_equation_error(y0,x,M_,options_,oo_,pfm,nodes,weights)
% e = euler_equation_error(y0,x,M_,options_,oo_,pfm,nodes,weights)
%
%  o  y0               [matrix]    path of endogenous, used to construct the guess values (initial condition not used; terminal condition used as guess value iff recompute_final_steady_state=true)
%  o  x                [matrix]    path of exogenous, used to construct the guess values (only if oo_.deterministic_simulation.controlled_paths_by_period is not empty)
%  o  M_               [structure] describing the model
%  o  options_         [structure] describing the options
%  o  oo_              [structure] describing the options
%  o  pfm              [struct]    perfect foresight model description
%  o  nodes            [double]    integration nodes
%  o  weights          [double]    weights of nodes
%
% Outputs: 
%  o  e                [double]    matrix of Euler equation errors  

% Called by ep_accuracy_check.m

% Copyright © 2016-2026 Dynare Team
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

dynamic_resid = str2func([M_.fname '.dynamic_resid']);
[y1, ~, ~, ~, pfm, options_] = extended_path_core(x,oo_.steady_state, pfm, M_, options_, oo_, []);

x1 = [x(2:end,:); zeros(M_.maximum_lead,M_.exo_nbr)];
res=NaN(M_.endo_nbr,size(x1,1));
for i=1:length(nodes)
    x2 = x1;
    x2(2,innovations.positive_var_indx) = x2(2,innovations.positive_var_indx) + nodes(i,:);
    y2 = extended_path_core(x2, y1, pfm, M_, options_, oo_, []);
    res(:,i) = dynamic_resid([y0; y1; y2],x(2,:),M_.params,oo_.steady_state);
end
e = res*weights;
