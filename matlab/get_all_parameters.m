function xparam1 = get_all_parameters(estim_params_, M_)
% xparam1 = get_all_parameters(estim_params_, M_)
% -------------------------------------------------------------------------
% gets parameters values from M_.params into xparam1 (inverse mapping to set_all_parameters)
% This is called if a model was calibrated before estimation to back out
% parameter values
%
% INPUTS
%    estim_params_:  Dynare structure describing the estimated parameters.
%    M_:             Dynare structure describing the model.
%
% OUTPUTS
%    xparam1:       N*1 double vector of parameters from calibrated model that are to be estimated
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2013-2026 Dynare Team
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

if isempty(estim_params_)
    estim_params_.nvx = 0;
    estim_params_.ncx = 0;
    estim_params_.nvn = 0;
    estim_params_.ncn = 0;
    estim_params_.nsx = 0;
    estim_params_.nendoinit = 0;
    estim_params_.np  = 0;
end

xparam1 = NaN(estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx+estim_params_.nendoinit+estim_params_.np,1);

% standard deviation of exogenous shocks (stderr on varexo, ordered first in xparam1)
if estim_params_.nvx
    var_exo = estim_params_.var_exo;
    for i = 1:estim_params_.nvx
        k = var_exo(i,1);
        xparam1(i) = sqrt(M_.Sigma_e(k,k));
    end
end
offset = estim_params_.nvx;

% standard deviation of measurement errors (stderr on varobs, ordered second in xparam1)
if estim_params_.nvn
    for i = 1:estim_params_.nvn
        k = estim_params_.nvn_observable_correspondence(i,1);
        xparam1(offset+i) = sqrt(M_.H(k,k));
    end
end
offset = estim_params_.nvx+estim_params_.nvn;

% correlations among shocks (corr on varexo, ordered third in xparam1)
if estim_params_.ncx
    corrx = estim_params_.corrx;
    for i = 1:estim_params_.ncx
        k1 = corrx(i,1);
        k2 = corrx(i,2);
        xparam1(i+offset) = M_.Correlation_matrix(k1,k2);
    end
end
offset = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx;

% correlations among measurement errors (corr on varobs, ordered fourth in xparam1)
if estim_params_.ncn
    corrn_observable_correspondence = estim_params_.corrn_observable_correspondence;
    for i=1:estim_params_.ncn
        k1 = corrn_observable_correspondence(i,1);
        k2 = corrn_observable_correspondence(i,2);
        xparam1(i+offset) = M_.Correlation_matrix_ME(k1,k2);
    end
end
offset = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn;

% skewness of exogenous shocks (skew on varexo, ordered fifth in xparam1)
if estim_params_.nsx
    skew_exo = estim_params_.skew_exo;
    for i = 1:estim_params_.nsx
        k = skew_exo(i,1);
        idx = M_.Skew_e(:,1)==k & M_.Skew_e(:,2)==k & M_.Skew_e(:,3)==k; % lookup (k,k,k) in sparse Skew_e using element-wise comparison
        if any(idx)
            xparam1(i+offset) = M_.Skew_e(idx,4);
        else
            xparam1(i+offset) = 0;
        end
    end
end
offset = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx;

% initial state variables (ordered sixth in xparam1)
if estim_params_.nendoinit
    endo_init_vals = estim_params_.endo_init_vals;
    for i = 1:estim_params_.nendoinit
        k = endo_init_vals(i,1);
        xparam1(i+offset) = M_.endo_initial_state.values(k);
    end
end
offset = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx+estim_params_.nendoinit;

% structural parameters (ordered last in xparam1)
if estim_params_.np
    xparam1(offset+1:end) = M_.params(estim_params_.param_vals(:,1));
end