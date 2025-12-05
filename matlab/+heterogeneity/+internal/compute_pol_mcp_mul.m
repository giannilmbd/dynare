function mult_values = compute_pol_mcp_mul(M_, ss, sizes, mat, indices)
% Compute policy function values for multipliers associated with binding
% complementarity constraints in models featuring rich heterogeneity.
%
% This routine evaluates residuals at grid points where inequality constraints are binding, and
% infers the corresponding values of Lagrange multipliers from the model's dynamic equations. It
% distinguishes lower-bound and upper-bound constraints, and only evaluates residuals where needed.
%
% INPUTS
%   M_     [struct]   : Dynare model structure
%   ss     [struct]   : Steady-state input validated by `check_steady_state_input`
%   sizes  [struct]   : Structure containing grid sizes, as returned by `check_steady_state_input`
%   mat    [struct]   : Structure with interpolation matrices and policy matrices
%
% OUTPUT
%   mult_values [struct] : Structure extending `ss.pol.values` with computed policy function
%                          values for multipliers (e.g. `MULT_L_a`, `MULT_U_c`)
%
% Internally, this function uses `M_.fname.dynamic_het1_resid` to evaluate residuals,
% and relies on complementarity mappings from `M_.fname.dynamic_het1_complementarity_conditions`.

% Copyright © 2025 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify it under the terms of the
% GNU General Public License as published by the Free Software Foundation, either version 3 of
% the License, or (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
% even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with Dynare. If not,
% see <https://www.gnu.org/licenses/>.
%
% Original author: Normann Rion <normann@dynare.org>
   % Get complementarity condition bounds
   [lb, ub] = feval(sprintf('%s.dynamic_het1_complementarity_conditions', M_.fname), M_.params);
   % Get constrained variables indices and names
   H_ = M_.heterogeneity(1);
   lb_constrained_var_names = H_.endo_names(lb ~= -Inf);
   ub_constrained_var_names = H_.endo_names(ub ~= Inf);
   lb = lb(lb ~= -Inf);
   ub = ub(ub ~= Inf);
   % Get the associated multipliers
   if isoctave()
      lb_mult = cellfun(@(x) ["MULT_L_" x], lb_constrained_var_names, 'UniformOutput', false);
      ub_mult = cellfun(@(x) ["MULT_U_" x], ub_constrained_var_names, 'UniformOutput', false);
   else
      lb_mult = cellfun(@(x) "MULT_L_"+x, lb_constrained_var_names, 'UniformOutput', false);
      ub_mult = cellfun(@(x) "MULT_U_"+x, ub_constrained_var_names, 'UniformOutput', false);
   end
   mult_ind = arrayfun(@(x) x.endo_index, H_.aux_vars);
   mult_names = H_.endo_names(mult_ind);
   % Get the associated equations
   eq_nbr = arrayfun(@(x) x.eq_nbr, H_.aux_vars);
   if isoctave()
      eq_nbr_lb = cellfun(@(x) eq_nbr(strcmp(x, mult_names)), lb_mult);
      eq_nbr_ub = cellfun(@(x) eq_nbr(strcmp(x, mult_names)), ub_mult);
   else
      eq_nbr_lb = cellfun(@(x) eq_nbr(x == mult_names), lb_mult);
      eq_nbr_ub = cellfun(@(x) eq_nbr(x == mult_names), ub_mult);
   end
   % Initialize the multipliers policy functions
   mult_values = struct;
   for i=1:numel(lb_mult)
      mult_name = lb_mult{i};
      var = lb_constrained_var_names{i};
      if ~isfield(ss.pol.values, mult_name)
         mult_values.(mult_name) = zeros(size(ss.pol.values.(var)));
      end
   end
   for i=1:numel(ub_mult)
      mult_name = ub_mult{i};
      var = ub_constrained_var_names{i};
      if ~isfield(ss.pol.values, mult_name)
         mult_values.(mult_name) = zeros(size(ss.pol.values.(var)));
      end
   end
   % Get the discretized policy function values of multipliers
   % Extended aggregate endogenous variable vector at the steady state
   y = zeros(3*M_.endo_nbr,1);
   for i=1:M_.orig_endo_nbr
      var = M_.endo_names{i};
      y(i) = ss.agg.(var); 
      y(M_.endo_nbr+i) = ss.agg.(var); 
      y(2*M_.endo_nbr+i) = ss.agg.(var); 
   end
   % Aggregate exogenous variable vector at the steady state
   x = zeros(M_.exo_nbr,1);
   % Heterogeneous endogenous variable vector
   yh = zeros(3*H_.endo_nbr,1);
   % Heterogeneous exogenous variable vector
   xh = zeros(H_.exo_nbr,1);
   % Get indices for shocks - Method 1 AR(1) are in H_.endo_names, Method 2a/2b are in H_.exo_names
   order_shocks_endo = cellfun(@(x) find(strcmp(x,H_.endo_names)), indices.shocks.endo);
   order_shocks_exo = cellfun(@(x) find(strcmp(x,H_.exo_names)), indices.shocks.exo);
   % Get the H_.endo_names indices for the steady-state ordering ss.pol.shocks and ss.pol.states 
   order_states = cellfun(@(x) find(strcmp(x,H_.endo_names)), indices.states);
   % Get the indices of declared ss.pol.values element in H_.endo_names
   all_pol_names = H_.endo_names;
   all_pol_names(ismember(all_pol_names, indices.shocks.all)) = [];
   all_pol_ind = cellfun(@(x) find(strcmp(x,H_.endo_names)), all_pol_names);
   % Get the indices of required declared pol values in H_.endo_names
   min_pol_names = H_.endo_names(1:H_.orig_endo_nbr);
   min_pol_names(ismember(min_pol_names, indices.shocks.all)) = [];
   min_pol_ind = cellfun(@(x) find(strcmp(x,H_.endo_names)), min_pol_names);
   % For each lower-bound constrained variables var, find the indices of values for
   % which ss.pol.values <= lb. For each of these indices, store minus the residual
   % of the equation associated with var. Proceed in a similar way for
   % upper-bound constrained variables, but store plus the residual.
   for i=1:numel(lb_constrained_var_names)
      var = lb_constrained_var_names{i};
      mult_name = lb_mult{i};
      if ~isfield(ss.pol.values, mult_name)
         ind = find(ss.pol.values.(var) <= lb(i));
         eq = eq_nbr_lb(i);
         for j=1:numel(ind)
            ind_sm = ind(j);
            % t-1
            yh(order_states) = mat.pol.sm(sizes.n_e+(1:sizes.n_a), ind_sm);
            % t
            yh(H_.endo_nbr+min_pol_ind) = cellfun(@(x) ss.pol.values.(x)(ind_sm), min_pol_names);
            % Method 1 AR(1) shocks go to yh (endogenous), Method 2a/2b go to xh (exogenous)
            if ~isempty(indices.shocks.endo)
               yh(H_.endo_nbr+order_shocks_endo) = mat.pol.sm(1:numel(indices.shocks.endo),ind_sm);
            end
            if ~isempty(indices.shocks.exo)
               xh(order_shocks_exo) = mat.pol.sm(numel(indices.shocks.endo)+(1:numel(indices.shocks.exo)),ind_sm);
            end
            % t+1
            yh(2*H_.endo_nbr+all_pol_ind) = mat.pol.x_bar_dash*mat.pol.Phi_e(:,ind_sm);
            % Residual function call
            r = feval([M_.fname '.dynamic_het1_resid'], y, x, M_.params, [], yh, xh, []);
            mult_values.(mult_name)(ind_sm) = -r(eq);
         end
      end
   end
   for i=1:numel(ub_constrained_var_names)
      var = ub_constrained_var_names{i};
      mult_name = ub_mult{i};
      if ~isfield(ss.pol.values, mult_name)
         ind = find(ss.pol.values.(var) >= ub(i));
         eq = eq_nbr_ub(i);
         for j=1:numel(ind)
            ind_sm = ind(j);
            % t-1
            yh(order_states) = mat.pol.sm(sizes.n_e+(1:sizes.n_a), ind_sm);
            % t
            yh(H_.endo_nbr+min_pol_ind) = cellfun(@(x) ss.pol.values.(x)(ind_sm), min_pol_names);
            % Method 1 AR(1) shocks go to yh (endogenous), Method 2a/2b go to xh (exogenous)
            if ~isempty(indices.shocks.endo)
               yh(H_.endo_nbr+order_shocks_endo) = mat.pol.sm(1:numel(indices.shocks.endo),ind_sm);
            end
            if ~isempty(indices.shocks.exo)
               xh(order_shocks_exo) = mat.pol.sm(numel(indices.shocks.endo)+(1:numel(indices.shocks.exo)),ind_sm);
            end
            % t+1
            yh(2*H_.endo_nbr+all_pol_ind) = mat.pol.x_bar_dash*mat.pol.Phi_e(:,ind_sm);
            % Residual function call
            r = feval([M_.fname '.dynamic_het1_resid'], y, x, M_.params, [], yh, xh, []);
            mult_values.(mult_name)(ind_sm) = r(eq);
         end
      end
   end
end
