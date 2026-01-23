function [shocks, spfm_exo_simul, oo_] = extended_path_shocks(pfm, exogenousvariables, sample_size, M_, options_, oo_)
% [shocks, spfm_exo_simul, oo_] = extended_path_shocks(pfm, exogenousvariables, sample_size, M_, options_, oo_)
% INPUTS
%  o pfm                 [struct]    description of pfm, containing information on innovations
%  o exogenousvariables  [matrix]    path of exogenous variables, potentially empty
%  o sample_size         [integer]   sample size
%  o M_                  [structure] describing the model
%  o options_            [structure] describing the options
%  o oo_                 [structure] storing the results
%
% OUTPUTS
%  o shocks             [matrix]    path of exogenous variables
%  o spfm_exo_simul     [matrix]    steady state
%  o oo_                [structure] storing the results
%
% Called by: extended_path.m, extended_path_mc.m

% Copyright © 2016-2025 Dynare Team
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

% Simulate shocks.
if isempty(exogenousvariables)
    switch options_.ep.innovation_distribution
      case 'gaussian'
        shocks = zeros(sample_size, M_.exo_nbr); %non-zero mean steady states are filtered out in extended_path.m
        shocks(:,pfm.positive_var_indx) = transpose(transpose(pfm.Omega)*randn(pfm.effective_number_of_shocks,sample_size)); %Omega is covariance_matrix_upper_cholesky 
      case 'calibrated'
        options = options_;
        options.periods = options.ep.periods;
        oo_local = make_ex_(M_, options, oo_);
        shocks = oo_local.exo_simul(2:end,:);
      otherwise
        error(['extended_path:: ' ep.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
    end
else
    shocks = exogenousvariables;
end

% Copy the shocks in exo_simul
oo_.exo_simul = [repmat(oo_.exo_steady_state',M_.maximum_lag,1); shocks];
spfm_exo_simul = repmat(oo_.exo_steady_state',options_.ep.periods+M_.maximum_lag+M_.maximum_lead,1);
