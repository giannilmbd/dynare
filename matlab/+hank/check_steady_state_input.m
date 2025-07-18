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
% Compute the stochastic simulations of heterogeneus-agent models
%
% INPUTS
% - M_       [structure] MATLAB's structure describing the model
% - options_ [structure] MATLAB's structure describing the current options
% - ss       [structure] MATLAB's structure with the model steady-state
%                       information. It contains the following fields:
%    - pol [structure] MATLAB's structure containing the policy functions
%                      discretization
%       - pol.grids [structure]: MATLAB's structure containing the nodes of the
%                                state grids as column vectors.
%       - pol.values [structure]: MATLAB's structure containing the policy
%                                 function as matrices. Row indices and column
%                                 indices follow the lexicographic order
%                                 specified in pol.shocks and pol.states
%       - pol.shocks [array]: a list containing the order of shock/innovation
%                             variables used for the rows of pol.values. If not
%                             specified, it follows the declaration order in the
%                             var(heterogeneity=) statement (resp. varexo(heterogeneity=)
%                             statement) for discretizes AR(1) processes (resp.
%                             discretized gaussian i.i.d innovations).
%       - pol.states [array]: a list containing the order of state variables
%                             used for the columns of pol.values. If not
%                             specified, it follows the declaration order in the
%                             var(heterogeneity=) statement.
%    - shocks [structure]: MATLAB's structure describing the discretization of
%                          individual shocks or innovations:
%       - shocks.grids [structure]: MATLAB's structure containing the nodes of
%                                   the shock grids
%       - shocks.Pi [structure]: MATLAB's structure containing the Markov
%                                matrices if the shock processes are discretized
%                                AR(1) processes. The field should be absent
%                                otherwise.
%       - shocks.w [structure]: MATLAB's structure containing the Gauss-Hermite
%                               weights if the i.i.d gaussian innovation
%                               processes are discretized
%    - d [structure]: MATLAB's structure describing the steady-state
%                     distribution
%       - d.grids [structure]: structure containing the states grids as column
%                              vectors. If one of the states grid is not
%                              specified, it is assumed to be identical to
%                              its policy-function value stored in pol.grids
%       - d.hist [matrix]: the histogram of the distribution as a matrix with
%                          shocks/innovation indices as rows and state indices
%                          as columns. Row and column indices follow the
%                          lexicographic order specified in d.shocks and
%                          d.states
%       - d.shocks [array]: a list containing the order of shock/innovation
%                           variables used for the rows of d.hist. If it is not
%                           specified, it falls back to the value induced by
%                           pol.shocks.
%       - d.states [array]: a list containing the order of state variables used
%                           for the columns of d.hist. If it is not
%                           specified, it falls back to the value induced by
%                           pol.states.
%    - agg [structure]: MATLAB's structure containing the steady-state values of
%                       aggregate variables
% OUTPUTS
% - out_ss [structure]: validated and normalized steady-state input structure.
%                       It has a similar structure as the ss structure in inputs.
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
function [out_ss, sizes] = check_steady_state_input(M_, options_, ss)
   %% Checks
   % Retrieve variables information from M_
   H_ = M_.heterogeneity(1);
   state_symbs = H_.endo_names(H_.state_var);
   inn_symbs = H_.exo_names;
   agg_symbs = M_.endo_names(1:M_.orig_endo_nbr);
   required_pol_symbs = H_.endo_names(1:H_.orig_endo_nbr);
   mult_symbs = {};
   for i = 1:numel(H_.aux_vars)
      if (H_.aux_vars(i).type == 15)
         mult_symbs{end+1} = H_.endo_names{H_.aux_vars(i).endo_index};
      end
   end
   % Initialize output variables
   sizes = struct;
   out_ss = struct;
   if ~isstruct(ss)
      error('Misspecified steady-state input `ss`: the steady-state input `ss` is not a structure.')
   end
   %% Shocks
   % Check that the field `ss.shocks` exists
   check_isfield('shocks', ss, 'ss.shocks', ' Models without idiosyncratic shocks but with ex-post heterogeneity over individual state variables cannot be solved yet.')
   shocks = ss.shocks;
   % Check that `ss.shocks` is a structure
   check_isstruct(shocks, 'ss.shocks');
   % Check that the field `ss.shocks.grids` exists
   check_isfield('grids', ss.shocks, 'ss.shocks.grids');
   % Check that `ss.shocks.grids` is a structure
   check_isstruct(shocks.grids, 'ss.shocks.grids');
   shock_symbs = fieldnames(shocks.grids);
   % Check the types of `ss.shocks.grids` values
   check_fun(shocks.grids, 'ss.shocks.grids', shock_symbs, @(x) isnumeric(x) && isreal(x) && isvector(x) && ~issparse(x), 'are not dense real vectors');
   % Make `ss.shocks.grids` values column vector in the output steady-state
   % structure
   for i=1:numel(shock_symbs)
      s = shock_symbs{i};
      if isrow(ss.shocks.grids.(s))
         out_ss.shocks.grids.(s) = ss.shocks.grids.(s)';
      else
         out_ss.shocks.grids.(s) = ss.shocks.grids.(s);
      end
   end
   % Check shock discretization
   flag_ar1 = isfield(shocks, 'Pi');
   flag_gh = isfield(shocks, 'w');
   if ~flag_ar1 && ~flag_gh
      error('Misspecified steady-state input `ss`: the `ss.shocks.Pi` and `ss.shocks.w` fields are both non-specified, while exactly one of them should be set.');
   end
   % Check that either individual AR(1) processes or i.i.d gaussion innovation
   % processes are discretized
   if flag_ar1 && flag_gh
      error('Misspecified steady-state input `ss`: the `ss.shocks.Pi` and `ss.shocks.w` fields are both specified, while only one of them should be.');
   end
   % Case of discretized AR(1) exogenous processes
   if flag_ar1
      % Check that the grids specification for AR(1) shocks and the
      % var(heterogeneity=) statement are compatible
      check_consistency(shock_symbs, required_pol_symbs, 'ss.shocks.grids', 'M_.heterogeneity(1).endo_names');
      % Check that shocks.Pi is a structure
      check_isstruct(shocks.Pi, 'ss.shocks.Pi');
      % Check that the grids specification for individual AR(1) shocks, the
      % information in the M_ structure and the Markov transition matrices
      % specification are compatible and display a warning if redundancies are
      % detected
      check_missingredundant(shocks.Pi, 'ss.shocks.Pi', shock_symbs, options_.hank.nowarningredundant);
      % Store:
      %    - the number of shocks
      sizes.n_e = numel(shock_symbs);
      %    - the size of the tensor product of the shock grids
      grid_nb = cellfun(@(s) numel(shocks.grids.(s)), shock_symbs); 
      sizes.N_e = prod(grid_nb);
      %    - the size of each shock grid
      for i=1:numel(shock_symbs)
         s = shock_symbs{i};
         sizes.shocks.(s) = numel(shocks.grids.(s));
      end
      % Check the type of shocks.Pi field values
      check_fun(shocks.Pi, 'ss.shocks.Pi', shock_symbs, @(x) ismatrix(x) && isnumeric(x) && isreal(x) && ~issparse(x), 'are not dense real matrices');
      % Check the internal size compatibility of discretized shocks processes
      check_fun(shocks.Pi, 'ss.shocks.Pi', shock_symbs, @(x) size(x,1), 'have a number of rows that is not consistent with the size of their counterparts in `ss.shocks.grids`', grid_nb);
      check_fun(shocks.Pi, 'ss.shocks.Pi', shock_symbs, @(x) size(x,2), 'have a number of columns that is not consistent with the size of their counterparts in `ss.shocks.grids`', grid_nb);
      % Check that the matrices provided in shocks.Pi are Markov matrices:
      %     - all elements are non-negative
      check_fun(shocks.Pi, 'ss.shocks.Pi', shock_symbs, @(x) all(x >= 0, 'all'), 'contain negative elements')
      %     - the sum of each row equals 1.
      check_fun(shocks.Pi, 'ss.shocks.Pi', shock_symbs, @(x) all(abs(sum(x,2)-1) < options_.hank.tol_check_sum), 'have row sums different from 1.')
      % Remove AR(1) shock processes from state and policy variables
      state_symbs = setdiff(state_symbs, shock_symbs);
      required_pol_symbs = setdiff(required_pol_symbs, shock_symbs);
      % Copy relevant shock-related elements in out_ss 
      out_ss.shocks.Pi = ss.shocks.Pi;
   end
   relevant_pol_symbs = [required_pol_symbs; mult_symbs];
   % Case of discretized i.i.d gaussian innovations
   if flag_gh
      % Check that shocks.w is a structure
      check_isstruct(shocks.w, 'shocks.w');
      % Verify that the grids specification for individual i.i.d innovations and
      % the varexo(heterogeneity=) statement are compatible
      check_consistency(shock_symbs, inn_symbs, 'ss.shocks.grids', 'M_.heterogeneity(1).exo_names');
      % Verify that the grids specification for individual i.i.d innovations and
      % the Gauss-Hermite weights specification are compatible and display a
      % warning if redundancies are detected
      check_missingredundant(shocks.w, 'ss.shocks.w', shock_symbs, options_.hank.nowarningredundant);
      % Check the type of `ss.shocks.w` elements
      check_fun(shocks.w, 'ss.shocks.w', shock_symbs, @(x) isvector(x) && isnumeric(x) && isreal(x) && ~issparse(x), 'are not dense real vectors');
      % Store:
      %    - the number of shocks
      sizes.n_e = numel(shock_symbs);
      %    - the size of the tensor product of the shock grids 
      grid_nb = cellfun(@(s) numel(shocks.grids.(s)), shock_symbs); 
      sizes.N_e = prod(grid_nb);
      %    - the size of each shock grid
      for i=1:numel(shock_symbs)
         s = shock_symbs{i};
         sizes.shocks.(s) = numel(shocks.grids.(s));
      end
      % Check the internal size compatibility of discretized shocks processes
      check_fun(shocks.w, 'ss.shocks.w', shock_symbs, @(x) numel(x), 'have incompatible sizes with those of `ss.shocks.grids`', grid_nb);
      % Check that the arrays provided in shocks.w have the properties of
      % Gauss-Hermite weights:
      %     - all elements are non-negative
      check_fun(shocks.w, 'ss.shocks.w', shock_symbs, @(x) all(x >= 0.), 'contain negative elements');
      %     - the sum of Gauss-Hermite weights is 1
      check_fun(shocks.w, 'ss.shocks.w', shock_symbs, @(x) all(abs(sum(x)-1) < options_.hank.tol_check_sum), 'have sums different from 1.');
      % Make `ss.shocks.w` values column vectors in the output steady-state
      % structure
      for i=1:numel(shock_symbs)
         s = shock_symbs{i};
         if isrow(ss.shocks.w.(s))
            out_ss.shocks.w.(s) = ss.shocks.w.(s)';
         else
            out_ss.shocks.w.(s) = ss.shocks.w.(s);
         end
      end
   end
   %% Policy functions
   % Check that the `ss.pol` field exists
   check_isfield('pol', ss, 'ss.pol');
   pol = ss.pol;
   % Check that `ss.pol` is a structure
   check_isstruct(ss.pol, 'ss.pol');
   % Check that the `ss.pol`.grids field exists
   check_isfield('grids', ss.pol, 'ss.pol.grids');
   % Check that `ss.pol.grids` is a structure
   check_isstruct(ss.pol.grids, 'ss.pol.grids');
   % Check that the state grids specification and the var(heterogeneity=)
   % statement are compatible
   check_missingredundant(pol.grids, 'ss.pol.grids', state_symbs, options_.hank.nowarningredundant);
   % Check that `ss.pol.grids` values are dense real vectors
   check_fun(pol.grids, 'ss.pol.grids', state_symbs, @(x) isnumeric(x) && isreal(x) && isvector(x) && ~issparse(x), 'are not dense real vectors');
   % Store:
   %    - the number of states
   sizes.n_a = numel(state_symbs);
   %    - the size of the tensor product of the states grids 
   grid_nb = cellfun(@(s) numel(pol.grids.(s)), state_symbs);
   sizes.pol.N_a = prod(grid_nb);
   %    - the size of each state grid
   for i=1:numel(state_symbs)
      s = state_symbs{i};
      sizes.pol.states.(s) = numel(pol.grids.(s));
   end
   % Make `ss.shocks.grids` values column vector in the output steady-state
   % structure
   for i=1:numel(state_symbs)
      s = state_symbs{i};
      if isrow(ss.pol.grids.(s))
         out_ss.pol.grids.(s) = ss.pol.grids.(s)';
      else
         out_ss.pol.grids.(s) = ss.pol.grids.(s);
      end
   end
   % Check that the field `ss.pol.values` exists 
   check_isfield('values', pol, 'ss.pol.values');
   % Check that `ss.pol.values` is a struct
   check_isstruct(pol.values, 'ss.pol.values');
   % Check the missing and redundant variables in `ss.pol.values`
   check_missingredundant(pol.values, 'ss.pol.values', required_pol_symbs, true);
   pol_values = fieldnames(pol.values);
   pol_values_in_relevant_pol_symbs = ismember(pol_values, relevant_pol_symbs);
   if ~options_.hank.nowarningredundant
      if ~all(pol_values_in_relevant_pol_symbs)
         warning('Steady-state input `ss`. The following fields are redundant in `%s`: %s.', 'ss.pol.values', strjoin(pol_values(~pol_values_in_relevant_pol_symbs)));
      end
   end
   provided_relevant_pol_symbs = pol_values(pol_values_in_relevant_pol_symbs);
   % Check that `ss.pol.values` values are dense real matrices
   check_fun(pol.values, 'ss.pol.values', provided_relevant_pol_symbs, @(x) isnumeric(x) && isreal(x) && ismatrix(x) && ~issparse(x), 'are not dense real matrices');
   % Check the internal size compatibility of `ss.pol.values`
   sizes.n_pol = H_.endo_nbr;
   check_fun(pol.values, 'ss.pol.values', provided_relevant_pol_symbs, @(x) size(x,1), 'have a number of rows that is not consistent with the sizes of `ss.shocks.grids` elements', sizes.N_e);
   check_fun(pol.values, 'ss.pol.values', provided_relevant_pol_symbs, @(x) size(x,2), 'have a number of columns that is not consistent with the sizes of `ss.pol.grids` elements', sizes.pol.N_a);
   % Copy `ss.pol.values` in `out_ss`
   out_ss.pol.values = ss.pol.values;
   % Check the permutation of state variables for policy functions
   out_ss.pol.states = check_permutation(pol, 'states', 'ss.pol', state_symbs);
   % Check the permutation of shock variables for policy functions
   out_ss.pol.shocks = check_permutation(pol, 'shocks', 'ss.pol', shock_symbs);
   %% Distribution
   % Check that the field `ss.d` exists 
   check_isfield('d', ss, 'ss.d');
   d = ss.d;
   % Check that the field `ss.d.hist` exists 
   check_isfield('hist', ss.d, 'ss.d.hist');
   % Check the type of `ss.d.hist`
   if ~(ismatrix(d.hist) && isnumeric(d.hist) && isreal(d.hist) && ~issparse(d.hist))
      error('Misspecified steady-state input `ss`: `ss.d.hist` is not a dense real matrix.');
   end
   % Check the consistency of `ss.d.grids` 
   if ~isfield(d, 'grids')
      if ~options_.hank.nowarningdgrids
         warning('In the steady-state input `ss.d.grids`, no distribution-specific grid is set for states %s. The policy grids in `ss.pol.grids` shall be used.' , strjoin(state_symbs));
      end
      % Copy the relevant sizes from the pol field
      sizes.d = sizes.pol;
      out_ss.d.grids = out_ss.pol.grids;
   else
      % Check that `ss.d.grids` is a struct
      check_isstruct(d.grids, 'ss.d.grids');
      % Check redundant variables in `ss.d.grids`
      d_grids_symbs = fieldnames(d.grids);
      if isempty(d_grids_symbs)
         if ~options_.hank.nowarningredundant
            warning('In the steady-state input `ss.d.grids`, no distribution-specific grid is set for states %s. The policy grids in `ss.pol.grids` shall be used.' , strjoin(state_symbs));
         end
         % Copy the relevant sizes from the pol field
         sizes.d = sizes.pol;
         out_ss.d.grids = out_ss.pol.grids;
      else
         d_grids_symbs_in_states = ismember(d_grids_symbs, state_symbs);
         if ~all(d_grids_symbs_in_states)
            if ~options_.hank.nowarningredundant
               warning('In the steady-state input `ss.d.states`, the following specification for the states grids in the distribution structure are not useful: %s.', strjoin(d_grids_symbs(~d_grids_symbs_in_states)));
            end
            d_grids_symbs = d_grids_symbs(d_grids_symbs_in_states);
         end
         % Check the types of `ss.d.grids` elements
         check_fun(pol.grids, 'ss.d.grids', d_grids_symbs, @(x) isnumeric(x) && isreal(x) && isvector(x) && ~issparse(x), 'are not dense real vectors');
         % Store the size of the states grids for the distribution
         for i=1:numel(d_grids_symbs)
            s = d_grids_symbs{i};
            sizes.d.states.(s) = numel(d.grids.(s));
            out_ss.d.grids.(s) = d.grids.(s);
         end
         states_out_of_d = setdiff(state_symbs, d_grids_symbs);
         if ~isempty(states_out_of_d)
            if ~options_.hank.nowarningdgrids
               warning('hank.bbeg.het_stoch_simul: in the steady-state structure, no distribution-specific grid set for states %s. The policy grids specified in pol.grids shall be used.', strjoin(states_out_of_d));
            end
            % Store the sizes of the unspecified states grids from the pol field
            for i=1:numel(states_out_of_d)
               s = states_out_of_d{i};
               sizes.d.states.(s) = sizes.pol.states.(s);
               out_ss.d.grids.(s) = pol.grids.(s);
            end
         end
      end
   end
   % Check the internal size compatibility of the distribution histogram
   if size(d.hist,1) ~= sizes.N_e
      error('Misspecified steady-state input `ss`: the number of rows of the histogram matrix `ss.d.hist` is not consistent with the sizes of shocks grids `ss.shocks.grids`.');
   end
   sizes.d.N_a = prod(structfun(@(x) x, sizes.d.states));
   if size(d.hist,2) ~= sizes.d.N_a
      error('Misspecified steady-state input `ss`: the number of columns of the histogram matrix `ss.d.hist` is not consistent with the sizes of states grids `ss.d.grids`/`ss.pol.grids`.');
   end
   % Copy `ss.d.hist` in `out_ss` 
   out_ss.d.hist = ss.d.hist;
   % Check the permutation of state variables in the distribution
   out_ss.d.states = check_permutation(d, 'states', 'ss.d', state_symbs);
   % Check the permutation of shock variables in the distribution
   out_ss.d.shocks = check_permutation(d, 'shocks', 'ss.d', shock_symbs);
   %% Aggregate variables
   % Check that the field `ss.agg` exists
   check_isfield('agg', ss, 'ss.agg');
   agg = ss.agg;
   % Check that `ss.agg` is a structure
   check_isstruct(agg, 'ss.agg');
   % Check that the aggregate variables specification and the var statement are
   % compatible
   check_missingredundant(agg, 'ss.agg', agg_symbs, options_.hank.nowarningredundant);
   % Store the number of aggregate variables
   sizes.agg = numel(fieldnames(agg));
   % Check the types of `ss.agg` values
   check_fun(agg, 'ss.agg', agg_symbs, @(x) isreal(x) && isscalar(x), 'are not real scalars');
   % Copy `ss.agg` into `out_ss`
   out_ss.agg = agg;
end