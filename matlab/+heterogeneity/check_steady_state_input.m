function [steady_state, sizes, indices] = check_steady_state_input(M_, options_het_check, steady_state, flag_initial_guess)
% Validate and process steady-state input for heterogeneous-agent models
% Supports AR(1) discretization only. I.i.d shocks are handled as AR(1) with zero persistence.
%
% INPUTS
% - M_ [structure] MATLAB's structure describing the model
% - options_het_check [structure] MATLAB's structure describing the current options
% - steady_state [structure] MATLAB's structure with the model steady-state
%                            information. It contains the following fields:
% - flag_initial_guess [boolean] (optional, default: false) If true, skip validation
%                                of d.hist and d.order (useful when the distribution will
%                                be computed rather than loaded). Validation of d.grids
%                                is still performed if d.grids is provided.
%    - pol [structure] MATLAB's structure containing the policy functions
%                      discretization
%       - pol.grids [structure]: Matlab's structure containing the nodes of the
%                                state grids as column vectors.
%       - pol.values [structure]: MATLAB's structure containing the policy
%                                 function as arrays. Dimensions follow the
%                                 lexicographic order specified in pol.order
%       - pol.order [array]: a list containing the order of variables used for
%                            the dimensions of pol.values.
%    - shocks [structure]: Matlab's structure describing the discretization of
%                          individual shocks or innovations:
%       - shocks.grids [structure]: MATLAB's structure containing the nodes of
%                                   the shock grids
%       - shocks.Pi [structure]: Matlab's structure containing the Markov
%                                matrices for AR(1) processes.
%    - d [structure]: Matlab's structure describing the steady-state
%                     distribution
%       - d.grids [structure]: structure containing the states grids as column
%                              vectors. If one of the states grid is not
%                              specified, it is assumed to be identical to
%                              its policy-function counterpart stored in pol.grids
%       - d.hist [array]: the histogram of the distribution as an array.
%                         Dimension indices follow the lexicographic order
%                         specified in d.order.
%       - d.order [array]: a list containing the order of variables used for the
%                          dimensions of d.hist. If it is not specified, it
%                          falls back to the value induced by pol.order.
%    - agg [structure]: Matlab's tructure containing the steady-state values of
%                       aggregate variables
% OUTPUTS
% - steady_state [structure]: validated and normalized steady-state structure.
%                             It has a similar structure as the input steady_state structure.
% - sizes [structure]: structure containing sizes information
%    - n_e [integer]: number of shocks
%    - N_e [integer]: total number of nodes in the shocks grids
%    - shocks [structure]: number of nodes for each shock grid
%    - n_a [integer]: number of states
%    - pol [structure]:
%       - pol.N_a [integer]: total number of nodes in the policy state grids
%       - pol.states [structure]: number of nodes for each state grid
%    - d [structure]:
%       - d.states [structure]: number of nodes for each distribution state grid
%       - d.N_a [integer]: total number of nodes in the distribution grids
%    - n_pol [integer]: number of policy variables
%    - agg [integer]: number of aggregate variables
% - indices [structure]: structure containing index and ordering information

% Copyright © 2025 Dynare Team
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
%
% Original author: Normann Rion <normann@dynare.org>
   % Handle optional skip_distribution parameter
   if nargin < 4
      flag_initial_guess = false;
   end

   %% Checks
   % Retrieve variables information from M_
   H_ = M_.heterogeneity(1);
   if M_.orig_endo_nbr == 0
      M_.endo_names = {};
   end
   agg_symbs = M_.endo_names(1:M_.orig_endo_nbr);
   orig_pol_symbs = H_.endo_names(1:H_.orig_endo_nbr);
   if ~isstruct(steady_state)
      error('Misspecified steady-state input `steady_state`: the steady-state input `steady_state` is not a structure.');
   end

   % Initialize output variables
   sizes = struct;
   indices = struct;
   indices.states = H_.endo_names(H_.state_var);
   indices.shocks = H_.exo_names;
   indices.order = [indices.shocks; indices.states];

   %% Aggregate variables
   if M_.orig_endo_nbr == 0
      steady_state.agg = struct;
   end
   % Check that the field `steady_state.agg` exists
   check_isfield('agg', steady_state, 'steady_state.agg');
   agg = steady_state.agg;
   % Check that `steady_state.agg` is a structure
   check_isstruct(agg, 'steady_state.agg');
   % Check that the aggregate variables specification and the var statement are
   % compatible
   check_missingredundant(agg, 'steady_state.agg', agg_symbs, options_het_check.no_warning_redundant);
   % Store the number of aggregate variables
   sizes.agg = numel(fieldnames(agg));
   % Check the types of `steady_state.agg` values
   check_fun(agg, 'steady_state.agg', agg_symbs, @(x) isreal(x) && isscalar(x), 'are not real scalars');
   % Check for NaN values
   check_fun(agg, 'steady_state.agg', agg_symbs, @(x) ~isnan(x), 'contain NaN values');

   %% Shocks
   % Check shocks structure
   check_isfield('shocks', steady_state, 'steady_state.shocks')
   shocks = steady_state.shocks;
   check_isstruct(shocks, 'steady_state.shocks')
   check_isfield('Pi', shocks);
   check_isstruct(shocks.Pi, 'steady_state.shocks.Pi');
   check_isfield('grids', shocks);
   check_isstruct(shocks.grids, 'steady_state.shocks.grids')

   for i = 1:H_.exo_nbr
      s = H_.exo_names{i};
      % Check if discretization is provided
      has_grid = isfield(shocks.grids, s);
      has_pi = isfield(shocks.Pi, s);
      if ~has_grid || ~has_pi
         % Both grids and Pi must be provided for all shocks
         error('Misspecified steady-state input `steady_state`: shock %s requires both `steady_state.shocks.grids.%s` and `steady_state.shocks.Pi.%s` to be provided.', s, s, s);
      end
      if has_grid
         % Validate grids
         if ~(isnumeric(shocks.grids.(s)) && isreal(shocks.grids.(s)) && isvector(shocks.grids.(s)) && ~issparse(shocks.grids.(s)))
            error('Misspecified steady-state input `steady_state`: `steady_state.shocks.grids.%s` is not a dense real vector.', s);
         end
         if isrow(shocks.grids.(s))
            steady_state.shocks.grids.(s) = shocks.grids.(s)';
         else
            steady_state.shocks.grids.(s) = shocks.grids.(s);
         end
         % Check for NaN values
         temp_struct.(s) = steady_state.shocks.grids.(s);
         check_fun(temp_struct, 'steady_state.shocks.grids', {s}, @(x) ~any(isnan(x)), 'contain NaN elements');
         % Check for strict monotonicity
         check_fun(temp_struct, 'steady_state.shocks.grids', {s}, @(x) all(diff(x) > 0), 'are not strictly increasing');
      end
      if has_pi
         % Validate Pi matrices
         if ~(ismatrix(shocks.Pi.(s)) && isnumeric(shocks.Pi.(s)) && isreal(shocks.Pi.(s)) && ~issparse(shocks.Pi.(s)))
            error('Misspecified steady-state input `steady_state`: `steady_state.shocks.Pi.%s` is not a dense real matrix.', s);
         end
         if has_grid && size(shocks.Pi.(s), 1) ~= numel(steady_state.shocks.grids.(s))
            error('Misspecified steady-state input `steady_state`: `steady_state.shocks.Pi.%s` row count does not match `steady_state.shocks.grids.%s` size.', s, s);
         end
         if has_grid && size(shocks.Pi.(s), 2) ~= numel(steady_state.shocks.grids.(s))
            error('Misspecified steady-state input `steady_state`: `steady_state.shocks.Pi.%s` column count does not match `steady_state.shocks.grids.%s` size.', s, s);
         end
         % Validate Markov matrix properties
         if (any(any(isnan(shocks.Pi.(s)))))
            error('Misspecified steady-state input `steady_state`: `steady_state.shocks.Pi.%s` contains NaN elements.', s);
         end
         if any(any(shocks.Pi.(s) < 0))
            error('Misspecified steady-state input `steady_state`: `steady_state.shocks.Pi.%s` contains negative elements.', s);
         end
         if any(abs(sum(shocks.Pi.(s), 2) - 1) > options_het_check.tol_check_sum)
            error('Misspecified steady-state input `steady_state`: `steady_state.shocks.Pi.%s` has row sums different from 1.', s);
         end
      end
   end

   % Check for redundant shock entries (warning only)
   if ~options_het_check.no_warning_redundant
      grids_fields = fieldnames(shocks.grids);
      grids_redundant = ~ismember(grids_fields, H_.exo_names);
      if any(grids_redundant)
         warning('Steady-state input `steady_state`. The following fields are redundant in `steady_state.shocks.grids`: %s.', strjoin(grids_fields(grids_redundant)));
      end
      pi_fields = fieldnames(shocks.Pi);
      pi_redundant = ~ismember(pi_fields, H_.exo_names);
      if any(pi_redundant)
         warning('Steady-state input `steady_state`. The following fields are redundant in `steady_state.shocks.Pi`: %s.', strjoin(pi_fields(pi_redundant)));
      end
   end

   % Store shock information
   sizes.n_e = H_.exo_nbr;
   sizes.N_e = 1;
   for i=1:H_.exo_nbr
      s = H_.exo_names{i};
      n = numel(steady_state.shocks.grids.(s));
      sizes.shocks.(s) = n;
      sizes.N_e = n*sizes.N_e;
   end

   %% Policy functions

   % Determine variables used in policy functions and classify shocks
   % First, check that `steady_state.pol` exists to determine what variables are used
   check_isfield('pol', steady_state, 'steady_state.pol');
   pol = steady_state.pol;
   check_isstruct(steady_state.pol, 'steady_state.pol');
   % Check that the field `steady_state.pol.values` exists
   check_isfield('values', pol, 'steady_state.pol.values');
   % Check that `steady_state.pol.values` is a struct
   check_isstruct(pol.values, 'steady_state.pol.values');
   % Check the missing and redundant variables in `steady_state.pol.values`
   check_missingredundant(pol.values, 'steady_state.pol.values', orig_pol_symbs, true);
   pol_values = fieldnames(pol.values);
   % Policy variables that could be present (including multipliers)
   pol_values_in_orig_pol_symbs = ismember(pol_values, orig_pol_symbs);
   if ~options_het_check.no_warning_redundant
      if ~all(pol_values_in_orig_pol_symbs)
         warning('Steady-state input `steady_state`. The following fields are redundant in `%s`: %s.', 'steady_state.pol.values', strjoin(pol_values(~pol_values_in_orig_pol_symbs)));
      end
   end
   % Check that `steady_state.pol.values` values are dense real matrices
   check_fun(pol.values, 'steady_state.pol.values', orig_pol_symbs, @(x) isnumeric(x) && isreal(x) && ~issparse(x), 'are not dense real arrays');
   % Check for NaN values
   check_fun(pol.values, 'steady_state.pol.values', orig_pol_symbs, @(x) ~any(isnan(x(:))), 'contain NaN elements');
   % Check that the `steady_state.pol.grids` field exists
   check_isfield('grids', steady_state.pol, 'steady_state.pol.grids');
   % Check that `steady_state.pol.grids` is a structure
   check_isstruct(steady_state.pol.grids, 'steady_state.pol.grids');
   % Check that the state grids specification and the var(heterogeneity=)
   % statement are compatible
   check_missingredundant(pol.grids, 'steady_state.pol.grids', indices.states, options_het_check.no_warning_redundant);
   % Check that `steady_state.pol.grids` values are dense real vectors
   check_fun(pol.grids, 'steady_state.pol.grids', indices.states, @(x) isnumeric(x) && isreal(x) && isvector(x) && ~issparse(x), 'are not dense real vectors');
   % Check for NaN values
   check_fun(pol.grids, 'steady_state.pol.grids', indices.states, @(x) ~any(isnan(x)), 'contain NaN elements');
   % Check for strict monotonicity
   check_fun(pol.grids, 'steady_state.pol.grids', indices.states, @(x) all(diff(x) > 0), 'are not strictly increasing');
   % Store:
   %    - the number of states
   sizes.n_a = numel(indices.states);
   %    - the size of each state grid and its product
   sizes.pol.N_a = 1;
   for i=1:sizes.n_a
      s = indices.states{i};
      n = numel(pol.grids.(s));
      sizes.pol.states.(s) = n;
      sizes.pol.N_a = n*sizes.pol.N_a;
   end
   sizes.N_sp = sizes.N_e*sizes.pol.N_a;
   % Make `steady_state.pol.grids` values column vector in the output steady-state
   % structure
   for i=1:sizes.n_a
      s = indices.states{i};
      if isrow(steady_state.pol.grids.(s))
         steady_state.pol.grids.(s) = steady_state.pol.grids.(s)';
      else
         steady_state.pol.grids.(s) = steady_state.pol.grids.(s);
      end
   end
   % Check the number of dimensions of `steady_state.pol.values` fields
   sizes.n_pol = H_.endo_nbr;
   check_fun(pol.values, 'steady_state.pol.values', orig_pol_symbs, @(x) ndims(x), 'have a number of dimensions that is not consistent with the number of states and shocks', sizes.n_a+sizes.n_e);
   % Validate that user's pol.order matches our canonical indices.order
   check_permutation(pol, 'order', 'steady_state.pol', indices.order);
   order = cellfun(@(x) find(strcmp(x,steady_state.pol.order)), indices.order);
   ind_shock = cellfun(@(x) find(strcmp(x,steady_state.pol.order)), indices.shocks);
   ind_state = cellfun(@(x) find(strcmp(x,steady_state.pol.order)), indices.states);
   % Check the internal size of dimensions of `steady_state.pol.values` fields
   for i=1:numel(ind_shock)
      d = ind_shock(i);
      var = steady_state.pol.order{d};
      check_fun(pol.values, 'steady_state.pol.values', orig_pol_symbs, @(x) size(x, d), sprintf('have a dimension-%s size that is not compatible with the size steady_state.shocks.grids.%s.', num2str(d), var), sizes.shocks.(var));
   end
   for i=1:numel(ind_state)
      d = ind_state(i);
      var = steady_state.pol.order{d};
      check_fun(pol.values, 'steady_state.pol.values', orig_pol_symbs, @(x) size(x, d), sprintf('have a dimension-%s size that is not compatible with the size steady_state.pol.grids.%s.', num2str(d), var), sizes.pol.states.(var));
   end
   % Order the dimensions of steady_state.pol.values in the indices.order
   for i = 1:numel(pol_values)
      var = pol_values{i};
      steady_state.pol.values.(var) = permute(steady_state.pol.values.(var), order);
   end

   %% Distribution
   if ~flag_initial_guess
      % Check that the field `steady_state.d` exists
      check_isfield('d', steady_state, 'steady_state.d');
      d = steady_state.d;
      % Check that the field `steady_state.d.hist` exists
      check_isfield('hist', steady_state.d, 'steady_state.d.hist');
      % Check the type of `steady_state.d.hist`
      if ~(isnumeric(d.hist) && isreal(d.hist) && ~issparse(d.hist))
         error('Misspecified steady-state input `steady_state`: `steady_state.d.hist` is not a dense real array.');
      end
      % Check for NaN values
      if any(isnan(d.hist(:)))
         error('Misspecified steady-state input `steady_state`: `steady_state.d.hist` contains NaN elements.');
      end
      % Validate d.grids (required when not skipping distribution)
      [steady_state, sizes] = validate_d_grids(steady_state, sizes, indices, pol, options_het_check, true);
      % Check the permutation of variables in the distribution histogram
      if isfield(steady_state.d, 'order')
         check_permutation(d, 'order', 'steady_state.d', indices.order);
      else
         steady_state.d.order = steady_state.pol.order;
      end
      % Check the internal size of dimensions of `steady_state.d.hist`
      ind_shock = cellfun(@(x) find(strcmp(x,steady_state.d.order)), indices.shocks);
      ind_state = cellfun(@(x) find(strcmp(x,steady_state.d.order)), indices.states);
      for i=1:numel(ind_shock)
         dim = ind_shock(i);
         var = steady_state.d.order{dim};
         if size(d.hist,dim) ~= sizes.shocks.(var)
            error('Misspecified steady-state input `steady_state`: dimension %s of `steady_state.d.hist` and `steady_state.shocks.grids.%s` have incompatible sizes.', num2str(dim), var);
         end
      end
      for i=1:numel(ind_state)
         dim = ind_state(i);
         var = steady_state.d.order{dim};
         if size(d.hist,dim) ~= sizes.d.states.(var)
            error('Misspecified steady-state input `steady_state`: dimension %s of `steady_state.d.hist` and `steady_state.pol.grids.%s` have incompatible sizes.', num2str(dim), var);
         end
      end
      % Validate that user's d.order matches our canonical indices.order
      order = cellfun(@(x) find(strcmp(x,steady_state.d.order)), indices.order);
      steady_state.d.hist = permute(steady_state.d.hist, order);
   else
      % When skipping d.hist and d.order validation, still validate d.grids if provided
      [steady_state, sizes] = validate_d_grids(steady_state, sizes, indices, pol, options_het_check, false);
   end

   % Compute sizes.d.N_a and sizes.N_om
   sizes.d.N_a = 1;
   for i = 1:sizes.n_a
      s = indices.states{i};
      sizes.d.N_a = sizes.d.N_a * sizes.d.states.(s);
   end
   sizes.N_om = sizes.d.N_a * sizes.N_e;

   % Store final order for variables
   steady_state.pol.order = indices.order;
   steady_state.d.order = indices.order;
end

function [steady_state, sizes] = validate_d_grids(steady_state, sizes, indices, pol, options_het_check, d_grids_required)
% Validate and process steady_state.d.grids
%
% INPUTS
%   steady_state      [struct] : Steady-state structure (modified in place)
%   sizes             [struct] : Sizes structure (modified in place)
%   indices           [struct] : Indices structure with .states field
%   pol               [struct] : Policy structure with .grids and .values
%   options_het_check [struct] : Options for warning suppression
%   d_grids_required  [bool]   : If true, d.grids field is required to exist
%
% OUTPUTS
%   steady_state [struct] : Updated with validated/defaulted d.grids
%   sizes        [struct] : Updated with d.states and d.N_a

   has_d_grids = isfield(steady_state, 'd') && isfield(steady_state.d, 'grids');

   if ~has_d_grids
      if d_grids_required && ~options_het_check.no_warning_d_grids
         warning('In the steady-state input `steady_state.d.grids`, no distribution-specific grid is set for states %s. The policy grids in `steady_state.pol.grids` shall be used.', strjoin(indices.states));
      end
      % Use pol grids as default
      sizes.d = sizes.pol;
      steady_state.d.grids = steady_state.pol.grids;
      return
   end

   d = steady_state.d;
   % Check that `steady_state.d.grids` is a struct
   check_isstruct(d.grids, 'steady_state.d.grids');

   d_grids_symbs = fieldnames(d.grids);
   if isempty(d_grids_symbs)
      if ~options_het_check.no_warning_redundant
         warning('In the steady-state input `steady_state.d.grids`, no distribution-specific grid is set for states %s. The policy grids in `steady_state.pol.grids` shall be used.', strjoin(indices.states));
      end
      % Use pol grids as default
      sizes.d = sizes.pol;
      steady_state.d.grids = steady_state.pol.grids;
      return
   end

   % Filter to valid state symbols
   d_grids_symbs_in_states = ismember(d_grids_symbs, indices.states);
   if ~all(d_grids_symbs_in_states)
      if ~options_het_check.no_warning_redundant
         warning('In the steady-state input `steady_state.d.grids`, the following specification for the states grids in the distribution structure are not useful: %s.', strjoin(d_grids_symbs(~d_grids_symbs_in_states)));
      end
      d_grids_symbs = d_grids_symbs(d_grids_symbs_in_states);
   end

   % Check the types of `steady_state.d.grids` elements
   check_fun(d.grids, 'steady_state.d.grids', d_grids_symbs, @(x) isnumeric(x) && isreal(x) && isvector(x) && ~issparse(x), 'are not dense real vectors');
   % Check for NaN values
   check_fun(d.grids, 'steady_state.d.grids', d_grids_symbs, @(x) ~any(isnan(x)), 'contain NaN elements');
   % Check for strict monotonicity
   check_fun(d.grids, 'steady_state.d.grids', d_grids_symbs, @(x) all(diff(x) > 0), 'are not strictly increasing');

   % Check that d.grids have the same bounds as pol.grids
   % Since both grids are validated to be strictly increasing, first/last elements are min/max
   for i = 1:numel(d_grids_symbs)
      s = d_grids_symbs{i};
      d_grid = d.grids.(s);
      pol_grid = pol.grids.(s);
      if d_grid(1) ~= pol_grid(1)
         error('Misspecified steady-state input `steady_state`: `steady_state.d.grids.%s` has minimum value %.6f which differs from `steady_state.pol.grids.%s` minimum %.6f.', s, d_grid(1), s, pol_grid(1));
      end
      if d_grid(end) ~= pol_grid(end)
         error('Misspecified steady-state input `steady_state`: `steady_state.d.grids.%s` has maximum value %.6f which differs from `steady_state.pol.grids.%s` maximum %.6f.', s, d_grid(end), s, pol_grid(end));
      end
   end

   % Store the size of the states grids for the distribution and ensure column vectors
   for i = 1:numel(d_grids_symbs)
      s = d_grids_symbs{i};
      sizes.d.states.(s) = numel(d.grids.(s));
      if isrow(steady_state.d.grids.(s))
         steady_state.d.grids.(s) = steady_state.d.grids.(s)';
      end
   end

   % Fill in missing states from pol grids
   states_out_of_d = setdiff(indices.states, d_grids_symbs);
   if ~isempty(states_out_of_d)
      if d_grids_required && ~options_het_check.no_warning_d_grids
         warning('heterogeneity.het_stoch_simul: in the steady-state structure, no distribution-specific grid set for states %s. The policy grids specified in pol.grids shall be used.', strjoin(states_out_of_d));
      end
      for i = 1:numel(states_out_of_d)
         s = states_out_of_d{i};
         sizes.d.states.(s) = sizes.pol.states.(s);
         steady_state.d.grids.(s) = pol.grids.(s);
      end
   end
end

function check_permutation(s, f_name, s_name, symbs)
% Check and returns a permutation order for variables.
%
% INPUTS
% - s [structure]: structure containing the permutation information
% - f_name [char]: name of the field within `s` that should contain the permutation list
% - s_name [char]: full name of the structure `s` for error messages
% - symbs [cell array of char]: list of expected variable names
%
% OUTPUTS
% - out_order [cell array or string array]: permutation order of variables
%
% DESCRIPTION
% Checks that the field `f_name` exists in structure `s` and contains a valid
% permutation (i.e., a cell array of character vectors or a string array).
% If the field is missing, returns the default order given by `symbs`.
% Throws an error if the field `f_name` exists in `s` and specified variables
% are not consistent with `symbs`.
   if isfield(s, f_name)
      f = s.(f_name);
      if ~iscellstr(f) && ~isstring(f)
         error('Misspecified steady-state input `steady_state`: the `%s.%s` field should be a cell array of character vectors or a string array.', s_name, f_name);
      end
      err_var = setdiff(f, symbs);
      if ~isempty(err_var)
         error('Misspecified steady-state input `steady_state`: the set of variables of the `%s.%s` field is not consistent with the information in M_ and the other fields in the steady-state structure `steady_state`. Problematic variables: %s.', s_name, f_name, strjoin(err_var));
      end
   end
end

function check_fun(f, f_name, symbs, fun, text, eq)
% Apply a validation function to a set of fields within a structure.
%
% INPUTS
% - f [structure]: structure containing fields to be validated
% - f_name [char]: name of the structure (for error messages)
% - symbs [cell array of strings]: list of field names to validate
% - fun [function handle]: validation function applied to each field
% - text [char]: descriptive text for the error message in case of validation failure
% - eq [optional, scalar]: if provided, the output of `fun` is compared to `eq`
%
% OUTPUTS
% - (none) This function throws an error if any field fails the validation.
%
% DESCRIPTION
% For each symbol in `symbs`, applies the validation function `fun` to the
% corresponding field in `f`. If an optional `eq` argument is provided, the
% output of `fun` is compared to `eq`. Throws a descriptive error listing all
% fields that do not satisfy the validation criteria.
   elt_check = cellfun(@(s) fun(f.(s)), symbs);
   if nargin == 6
      elt_check = (elt_check == eq);
   end
   if ~all(elt_check)
      error('Misspecified steady-state input `steady_state`: the following fields in `%s` %s: %s.', f_name, text, strjoin(symbs(~elt_check)));
   end
end

function check_isfield(f, s, f_name, details)
% Check the existence of a field in a structure.
%
% INPUTS
% - f [char]: name of the field to check
% - s [structure]: structure in which the field should be present
% - f_name [optional, char]: name to display in the error message (defaults to `f`)
% - details [optional, char]: additional details to append to the error message
%
% OUTPUTS
% - (none) This function throws an error if the specified field is missing.
%
% DESCRIPTION
% Checks whether the field `f` exists in the structure `s`.
% If not, throws an error indicating the missing field.
% If provided, `details` is appended to the error message for additional context.
   if nargin < 4
      details = '';
   end
   if nargin < 3
      f_name = f;
   end
   if ~isfield(s, f)
      error('Misspecified steady-state input `steady_state`: the `%s` field is missing.%s', f_name, details);
   end
end

function check_isstruct(f, f_name)
% Check that an input is a structure.
%
% INPUTS
% - f [any type]: variable to check
% - f_name [char]: name to display in the error message
%
% OUTPUTS
% - (none) This function throws an error if the input is not a structure.
%
% DESCRIPTION
% Checks whether `f` is a structure.
% If not, throws an error indicating that the field `f_name` should be a structure.
   if ~isstruct(f)
      error('Misspecified steady-state input `steady_state`: the `%s` field should be a structure.', f_name);
   end
end

function check_missingredundant(f, f_name, symbs, no_warning_redundant)
% Check for missing and redundant fields in a structure.
%
% INPUTS
% - f [structure]: structure to check
% - f_name [char]: name to display in error/warning messages
% - symbs [cell array of char]: list of expected field names
% - no_warning_redundant [logical]: if true, suppresses warnings about redundant fields
%
% OUTPUTS
% - (none) This function throws an error if any expected field is missing,
%   and issues a warning if redundant fields are found (unless no_warning_redundant is true).
%
% DESCRIPTION
% Verifies that all fields listed in `symbs` are present in the structure `f`.
% Throws an error if any expected field is missing.
% If `no_warning_redundant` is false, checks if `f` contains fields not listed
% in `symbs`, and issues a warning if redundant fields are found.
   fields = fieldnames(f);
   symbs_in_fields = ismember(symbs, fields);
   if ~all(symbs_in_fields)
      error('Misspecified steady-state input `steady_state`. The following fields are missing in `%s`: %s.', f_name, strjoin(symbs(~symbs_in_fields)));
   end
   if ~no_warning_redundant
      fields_in_symbs = ismember(fields, symbs);
      if ~all(fields_in_symbs)
         warning('Steady-state input `steady_state`. The following fields are redundant in `%s`: %s.', f_name, strjoin(fields(~fields_in_symbs)));
      end
   end
end