function xparam = get_posterior_parameters(type,M_,estim_params_,oo_,options_,field1)
% xparam = get_posterior_parameters(type,M_,estim_params_,oo_,options_,field1)
% -------------------------------------------------------------------------
% Selects (estimated) parameters (posterior mode or posterior mean).
%
% INPUTS
%   o type              [char]     = 'mode' or 'mean'.
%   o M_:               [structure] Dynare structure describing the model.
%   o estim_params_:    [structure] Dynare structure describing the estimated parameters.
%   o oo_:              [structure] Dynare results structure
%   o options_:         [structure] Dynare options structure
%   o field_1           [char]      optional field like 'mle_'.
%
% OUTPUTS
%   o xparam     vector of estimated parameters
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2006-2025 Dynare Team
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

if nargin<6
    field1='posterior_';
end

xparam = zeros(estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.np,1);
ip = 1; % initialize index in xparam

for i=1:estim_params_.nvx % estimated stderr parameters for structural shocks (ordered first in xparam1)
    k1 = estim_params_.var_exo(i,1);
    name1 = M_.exo_names{k1};
    xparam(ip) = oo_.([field1 type]).shocks_std.(name1);
    ip = ip + 1;
end

for i=1:estim_params_.nvn % estimated stderr parameters for measurement errors (ordered second in xparam1)
    k1 = estim_params_.nvn_observable_correspondence(i,1);
    name1 = options_.varobs{k1};
    xparam(ip) = oo_.([field1 type]).measurement_errors_std.(name1);
    ip = ip + 1;
end

for i=1:estim_params_.ncx % estimated corr parameters for structural shocks (ordered third in xparam1)
    k1 = estim_params_.corrx(i,1);
    k2 = estim_params_.corrx(i,2);
    name1 = M_.exo_names{k1};
    name2 = M_.exo_names{k2};
    xparam(ip) = oo_.([field1 type]).shocks_corr.([name1 '_' name2]);
    ip = ip + 1;
end

for i=1:estim_params_.ncn % estimated corr parameters for measurement errors (ordered fourth in xparam1)
    k1 = estim_params_.corrn_observable_correspondence(i,1);
    k2 = estim_params_.corrn_observable_correspondence(i,2);
    name1 = options_.varobs{k1};
    name2 = options_.varobs{k2};
    xparam(ip) = oo_.([field1 type]).measurement_errors_corr.([name1 '_' name2]);
    ip = ip + 1;
end

for i=1:estim_params_.np % estimated structural parameters (ordered last in xparam1)
    name1 = M_.param_names{estim_params_.param_vals(i,1)};
    xparam(ip) = oo_.([field1 type]).parameters.(name1);
    ip = ip + 1;
end