function M_ = set_parameters_locally(M_,xparam1)
% M_ = set_parameters_locally(M_,xparam1)
% -------------------------------------------------------------------------
% Sets parameters value (except measurement errors)
% This is called for computations such as IRF and forecast
% when measurement errors aren't taken into account; in contrast to
% set_parameters.m, the global M_-structure is not altered
%
% INPUTS
%    xparam1:   vector of parameters to be estimated (initial values)
%    M_:        Dynare model-structure
%
% OUTPUTS
%    M_:        Dynare model-structure
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2017-2025 Dynare Team
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

global estim_params_

Sigma_e = M_.Sigma_e;
Correlation_matrix = M_.Correlation_matrix;

% setting shocks variance on the diagonal of Covariance matrix; used later for updating covariances
if estim_params_.nvx % stderr VAREXO are ordered first in xparam1
    var_exo = estim_params_.var_exo;
    for i=1:estim_params_.nvx
        k = var_exo(i,1);
        Sigma_e(k,k) = xparam1(i)^2;
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
if estim_params_.nendoinit % structural parameters are ordered last in xparam1
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
Sigma_e = diag(sqrt(diag(Sigma_e)))*Correlation_matrix*diag(sqrt(diag(Sigma_e)));
if isfield(estim_params_,'calibrated_covariances') % if calibrated covariances, set them now to their stored value
    Sigma_e(estim_params_.calibrated_covariances.position)=estim_params_.calibrated_covariances.cov_value;
end
% updating matrices in M_
M_.Sigma_e = Sigma_e;
M_.Correlation_matrix = Correlation_matrix;

% update specification of independent closed skew normal distributed shocks
M_.csn = csn_update_specification(M_.Sigma_e, M_.Skew_e);