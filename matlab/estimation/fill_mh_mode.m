function oo_ = fill_mh_mode(xparam1, stdh, M_, options_, estim_params_, oo_, field_name)
% oo_ = fill_mh_mode(xparam1, stdh, M_, options_, estim_params_, oo_, field_name)
% -------------------------------------------------------------------------
% Fill oo_.<field_name>.mode and oo_.<field_name>.std_at_mode
%
% INPUTS
% - xparam1         [double]  p×1 vector, estimated posterior mode.
% - stdh            [double]  p×1 vector, estimated posterior standard deviation.
% - M_              [struct]  Description of the model.
% - estim_params_   [struct]  Description of the estimated parameters.
% - options_        [struct]  Dynare's options.
% - oo_             [struct]  Estimation and simulation results.
%
% OUTPUTS
% - oo_                       MATLAB's structure gathering the results
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2005-2025 Dynare Team
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

if estim_params_.np % estimated structural parameters
    ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+1; % offset: structural parameters are ordered last in xparam1
    for i=1:estim_params_.np
        k = estim_params_.param_vals(i,1);
        name = M_.param_names{k};
        oo_.([field_name '_mode']).parameters.(name) = xparam1(ip);
        oo_.([field_name '_std_at_mode']).parameters.(name) = stdh(ip);
        ip = ip+1;
    end
end

if estim_params_.nvx % estimated stderr parameters for structural shocks
    ip = 1; % offset: stderr parameters for structural shocks are ordered first in xparam1
    for i=1:estim_params_.nvx
        k = estim_params_.var_exo(i,1);
        name = M_.exo_names{k};
        oo_.([field_name '_mode']).shocks_std.(name)= xparam1(ip);
        oo_.([field_name '_std_at_mode']).shocks_std.(name) = stdh(ip);
        ip = ip+1;
    end
end

if estim_params_.nvn % estimated stderr parameters for measurement errors
    ip = estim_params_.nvx+1; % offset: stderr parameters for measurement errors are ordered second in xparam1
    for i=1:estim_params_.nvn
        name = options_.varobs{estim_params_.nvn_observable_correspondence(i,1)};
        oo_.([field_name '_mode']).measurement_errors_std.(name) = xparam1(ip);
        oo_.([field_name '_std_at_mode']).measurement_errors_std.(name) = stdh(ip);
        ip = ip+1;
    end
end

if estim_params_.ncx % estimated corr parameters for structural shocks
    ip = estim_params_.nvx+estim_params_.nvn+1; % offset: corr parameters for structural shocks are ordered third in xparam1
    for i=1:estim_params_.ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        name = [M_.exo_names{k1} '_' M_.exo_names{k2}];
        oo_.([field_name '_mode']).shocks_corr.(name) = xparam1(ip);
        oo_.([field_name '_std_at_mode']).shocks_corr.(name) = stdh(ip);
        ip = ip+1;
    end
end

if estim_params_.ncn % estimated corr parameters for measurement errors
    ip = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+1; % offset: corr parameters for measurement errors are ordered fourth in xparam1
    for i=1:estim_params_.ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        name = [M_.endo_names{k1} '_' M_.endo_names{k2}];
        oo_.([field_name '_mode']).measurement_errors_corr.(name) = xparam1(ip);
        oo_.([field_name '_std_at_mode']).measurement_errors_corr.(name) = stdh(ip);
        ip = ip+1;
    end
end