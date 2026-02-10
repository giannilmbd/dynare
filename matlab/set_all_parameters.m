function M_ = set_all_parameters(xparam1,estim_params_,M_)
% M_ = set_all_parameters(xparam1,estim_params_,M_)
% -------------------------------------------------------------------------
% Update parameter values (deep parameters, covariance and skewness matrices)
% in the Dynare model structure M_ using the estimated parameters xparam1 and
% the estimation parameter structure estim_params_.
%
% Inputs:
%   xparam1        [N x 1 double]   Vector of N estimated parameter values.
%   estim_params_  [struct]         Structure describing the estimated parameters.
%   M_             [struct]         Structure describing the model.
%
% Output:
%   M_             [struct]         Updated model structure with new parameters and covariance matrices.
%
% This function is called by:
%   DsgeSmoother, dynare_estimation_1, gsa.monte_carlo_filtering,
%   identification.analysis, PosteriorFilterSmootherAndForecast,
%   prior_posterior_statistics_core, prior_sampler
% -------------------------------------------------------------------------

% Copyright © 2003-2025 Dynare Team
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

if estim_params_.nvx || estim_params_.ncx
    Sigma_e = M_.Sigma_e;
    Correlation_matrix = M_.Correlation_matrix;
end

H = M_.H;
Correlation_matrix_ME = M_.Correlation_matrix_ME;

% setting shocks variance on the diagonal of Covariance matrix; used later for updating covariances
if estim_params_.nvx % stderr VAREXO are ordered first in xparam1
    var_exo = estim_params_.var_exo;
    for i=1:estim_params_.nvx
        k = var_exo(i,1);
        Sigma_e(k,k) = xparam1(i)^2;
    end
end

% setting measurement error variance; on the diagonal of Covariance matrix; used later for updating covariances
offset = estim_params_.nvx;
if estim_params_.nvn % stderr VAROBS are ordered second in xparam1
    for i=1:estim_params_.nvn
        k = estim_params_.nvn_observable_correspondence(i,1);
        H(k,k) = xparam1(i+offset)^2;
    end
end

% correlations among shocks
offset = estim_params_.nvx + estim_params_.nvn;
if estim_params_.ncx % corr among VAREXO are ordered third in xparam1
    corrx = estim_params_.corrx;
    for i=1:estim_params_.ncx
        k1 = corrx(i,1);
        k2 = corrx(i,2);
        Correlation_matrix(k1,k2) = xparam1(i+offset);
        Correlation_matrix(k2,k1) = Correlation_matrix(k1,k2);
    end
end

% setting measurement error covariances
offset = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx;
if estim_params_.ncn % corr among VAROBS are ordered fourth in xparam1
    corrn_observable_correspondence = estim_params_.corrn_observable_correspondence;
    for i=1:estim_params_.ncn
        k1 = corrn_observable_correspondence(i,1);
        k2 = corrn_observable_correspondence(i,2);
        Correlation_matrix_ME(k1,k2) = xparam1(i+offset);
        Correlation_matrix_ME(k2,k1) = Correlation_matrix_ME(k1,k2);
    end
end

% setting skewness coefficients of structural shocks
offset = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn;
if estim_params_.nsx % skew among VAREXO are ordered fifth in xparam1
    skew_exo = estim_params_.skew_exo;
    for i=1:estim_params_.nsx
        k = skew_exo(i,1);
        % Remove existing entry for (k,k,k) in sparse Skew_e
        idx = M_.Skew_e(:,1)==k & M_.Skew_e(:,2)==k & M_.Skew_e(:,3)==k;
        M_.Skew_e(idx,:) = [];
        % Append new entry if non-zero
        if xparam1(i+offset) ~= 0
            M_.Skew_e = [M_.Skew_e; k, k, k, xparam1(i+offset)];
        end
    end
end

% setting endogenous init states
offset = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn + estim_params_.nsx;
if estim_params_.nendoinit % initial states are ordered sixth in xparam1
    endo_init = estim_params_.endo_init_vals;
    for i=1:estim_params_.nendoinit
        k = endo_init(i,1);
        M_.endo_initial_state.values(k) = xparam1(i+offset);
    end
end

% setting structural parameters
offset = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn + estim_params_.nsx + estim_params_.nendoinit;
if estim_params_.np % structural parameters are ordered last in xparam1
    M_.params(estim_params_.param_vals(:,1)) = xparam1(offset+1:end);
end

% build shock covariance matrix from correlation matrix and variances already on diagonal
if estim_params_.nvx || estim_params_.ncx
    Sigma_e = diag(sqrt(diag(Sigma_e)))*Correlation_matrix*diag(sqrt(diag(Sigma_e)));
    if isfield(estim_params_,'calibrated_covariances') % if calibrated covariances, set them now to their stored value
        Sigma_e(estim_params_.calibrated_covariances.position)=estim_params_.calibrated_covariances.cov_value;
    end
    % updating matrices in M_
    M_.Sigma_e = Sigma_e;
    M_.Correlation_matrix=Correlation_matrix;
end

% build measurement error covariance matrix from correlation matrix and variances already on diagonal
if estim_params_.nvn || estim_params_.ncn
    H = diag(sqrt(diag(H)))*Correlation_matrix_ME*diag(sqrt(diag(H)));
    if isfield(estim_params_,'calibrated_covariances_ME') % if calibrated covariances, set them now to their stored value
        H(estim_params_.calibrated_covariances_ME.position)=estim_params_.calibrated_covariances_ME.cov_value;
    end
    % updating matrices in M_
    M_.H = H;
    M_.Correlation_matrix_ME=Correlation_matrix_ME;
end

% update specification of independent closed skew normal distributed shocks
M_.csn = csn_update_specification(M_.Sigma_e, M_.Skew_e);