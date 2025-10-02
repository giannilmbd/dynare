% -------------------------------------------------------------------------
% Functionality testing of Bayesian IRF matching with
% - Random-Walk Metropolis-Hastings
% - whether mh_posterior_mode_estimation option works
% -------------------------------------------------------------------------

% Copyright © 2023-2025 Dynare Team
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
@#define LESSPARAMS=1
@#include "cet_model.inc"

options_.prior_interval= 0.95;

method_of_moments(mom_method = irf_matching
%, add_tiny_number_to_cholesky = 1e-14
, irf_matching_file = cet_irf_matching_file
, mh_conf_sig = 0.90
, mh_replic=100
, mh_posterior_mode_estimation
, mode_file = cet_original_mode
, verbose
);
