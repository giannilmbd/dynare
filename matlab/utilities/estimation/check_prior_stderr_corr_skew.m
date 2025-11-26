function check_prior_stderr_corr_skew(estim_params_,bayestopt_)
% function check_prior_stderr_corr_skew(estim_params_,bayestopt_)
% -------------------------------------------------------------------------
% Issue a warning if the prior implies:
% - negative standard deviations
% - correlations larger than +-1
% - skewness coefficient outside theoretical bound (currently only for skew normal distribution)
% -------------------------------------------------------------------------
% INPUTS
%  o estim_params_:           [struct] information on estimated parameters
%  o bayestopt_:              [struct] information on priors
% -------------------------------------------------------------------------
% OUTPUTS
%  No outputs; issues warnings for invalid prior standard deviation, correlation, or skewness coefficients
% -------------------------------------------------------------------------
% This function is called by
%  o initial_estimation_checks.m
% -------------------------------------------------------------------------

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

% estimated stderr parameters for structural shocks (ordered first in xparam1)
if estim_params_.nvx && any(bayestopt_.p3(1:estim_params_.nvx)<0)
    warning('Your prior allows for negative standard deviations for structural shocks. Due to working with variances, Dynare will be able to continue, but it is recommended to change your prior!');
end
offset = estim_params_.nvx;

 % estimated stderr parameters for measurement errors (ordered second in xparam1)
if estim_params_.nvn && any(bayestopt_.p3(1+offset:offset+estim_params_.nvn)<0)
    warning('Your prior allows for negative standard deviations for measurement error. Due to working with variances, Dynare will be able to continue, but it is recommended to change your prior!');
end
offset = estim_params_.nvx+estim_params_.nvn;

% estimated corr parameters for structural shocks (ordered third in xparam1)
if estim_params_.ncx && (any(bayestopt_.p3(1+offset:offset+estim_params_.ncx)<-1) || any(bayestopt_.p4(1+offset:offset+estim_params_.ncx)>1))
    warning('Your prior allows for correlations between structural shocks larger than +-1 and will not integrate to 1 due to truncation. Please change your prior!');
end
offset = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx;

% estimated corr parameters for measurement errors (ordered fourth in xparam1)
if estim_params_.ncn && (any(bayestopt_.p3(1+offset:offset+estim_params_.ncn)<-1) || any(bayestopt_.p4(1+offset:offset+estim_params_.ncn)>1))
    warning('Your prior allows for correlations between measurement errors larger than +-1 and will not integrate to 1 due to truncation. Please change your prior!');
end
offset = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn;

% estimated skew parameters for structural shocks (ordered fifth in xparam1)
theoretical_bound_sn_skew = abs((sqrt(2)*(pi-4))/(pi-2)^(3/2)); % theoretical bound on skewness coefficient of skew normal distributed shock
if estim_params_.nsx && (any(bayestopt_.p3(1+offset:offset+estim_params_.nsx)<-theoretical_bound_sn_skew) || any(bayestopt_.p4(1+offset:offset+estim_params_.nsx)>theoretical_bound_sn_skew))
    warning('Your prior allows for skewness parameters larger than +-%f, which is outside the theoretical bound of the skew normal distribution, and thus will not integrate to 1 due to truncation. Please change your prior!', theoretical_bound_sn_skew);
end