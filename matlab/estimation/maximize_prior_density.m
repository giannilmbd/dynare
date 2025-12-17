function [xparams, lpd, hessian_mat] = maximize_prior_density(iparams, options_, M_, bayestopt_, estim_params_, oo_)
% [xparams, lpd, hessian_mat] = maximize_prior_density(iparams, options_, M_, bayestopt_, estim_params_, oo_)
% -------------------------------------------------------------------------
% Maximizes the logged prior density.
%
% INPUTS
% - iparams            [double]   vector of initial parameters.
% - options_           [struct]   structure describing the options
% - bayestopt_         [struct]   structure describing the priors
% - estim_params_      [struct]   structure describing the estimated parameters
% - oo_                [struct]   structure describing the results
%
% OUTPUTS
% - xparams            [double]   vector, prior mode.
% - lpd                [double]   scalar, value of the logged prior density at the mode.
% - hessian_mat        [double]   matrix, Hessian matrix at the prior mode.

% Copyright © 2009-2025 Dynare Team
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

[xparams, lpd, ~, hessian_mat] = dynare_minimize_objective('minus_logged_prior_density', ...
    iparams, ...
    options_.mode_compute, ...
    options_, ...
    [bayestopt_.p3, bayestopt_.p4], ...
    bayestopt_.name, ...
    bayestopt_, ...
    [], ...
    bayestopt_, ...
    options_, ...
    M_, ...
    estim_params_, ...
    oo_);

lpd = -lpd;