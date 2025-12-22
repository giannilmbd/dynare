function [nam, texnam] = get_the_name(k, TeX, M_, estim_params_, varobs)
% [nam, texnam] = get_the_name(k, TeX, M_, estim_params_, varobs)
% -------------------------------------------------------------------------
% Returns name of estimated parameter number k, following the internal ordering of 
% the estimated parameters.
% Inputs:
%   - k             [integer]   parameter number.
%   - TeX           [bool]      if false, texnam is not returned (empty matrix)
%   - M_            [structure] model
%   - estim_params_ [structure] describing the estimated parameters
%   - varobs        [cell]      name of observed variables
%
% Outputs
%   - nam       [char]      internal name of the variable
%   - texnam    [char]      TeX name of the same variable (if defined in the mod file)
%
% This function is called by:
% get_prior_info, mcmc_diagnostics, mode_check, PlotPosteriorDistributions, plot_priors
%
% This function calls:
% None.
% 

% Copyright © 2004-2025 Dynare Team
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

texnam = [];

if k <= estim_params_.nvx % estimated stderr parameters for structural shocks (ordered first in xparam1)
    vname = M_.exo_names{estim_params_.var_exo(k,1)};
    nam = sprintf('SE_%s', vname);
    if TeX
        tname  = M_.exo_names_tex{estim_params_.var_exo(k,1)};
        texnam = sprintf('$ \\sigma_{%s} $', tname);
    end
elseif  k <= (estim_params_.nvx+estim_params_.nvn) % estimated stderr parameters for measurement errors (ordered second in xparam1)
    vname = varobs{estim_params_.nvn_observable_correspondence(k-estim_params_.nvx,1)};
    nam = sprintf('SE_EOBS_%s', vname);
    if TeX
        tname  = M_.endo_names_tex{estim_params_.var_endo(k-estim_params_.nvx,1)};
        texnam = sprintf('$ \\sigma^{ME}_{%s} $', tname);
    end
elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx) % estimated corr parameters for structural shocks (ordered third in xparam1)
    jj = k - (estim_params_.nvx+estim_params_.nvn);
    k1 = estim_params_.corrx(jj,1);
    k2 = estim_params_.corrx(jj,2);
    vname = sprintf('%s_%s', M_.exo_names{k1}, M_.exo_names{k2});
    nam = sprintf('CC_%s', vname);
    if TeX
        tname  = sprintf('%s,%s', M_.exo_names_tex{k1}, M_.exo_names_tex{k2});
        texnam = sprintf('$ \\rho_{%s} $', tname);
    end
elseif  k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn) % estimated corr parameters for measurement errors (ordered fourth in xparam1)
    jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx);
    k1 = estim_params_.corrn(jj,1);
    k2 = estim_params_.corrn(jj,2);
    vname = sprintf('%s_%s', M_.endo_names{k1}, M_.endo_names{k2});
    nam = sprintf('CC_EOBS_%s', vname);
    if TeX
        tname  = sprintf('%s,%s', M_.endo_names_tex{k1}, M_.endo_names_tex{k2});
        texnam = sprintf('$ \\rho^{ME}_{%s} $', tname);
    end
elseif k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx) % estimated skew parameters for structural shocks (ordered fifth in xparam1)
    jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn);
    k = estim_params_.skew_exo(jj,1);
    vname = sprintf('%s', M_.exo_names{k});
    nam = sprintf('SKEW_%s', vname);
    if TeX
        tname  = sprintf('%s', M_.exo_names_tex{k});
        texnam = sprintf('$ skew({%s}) $', tname);
    end
elseif k <= (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx+estim_params_.nendoinit) % estimated initial states (ordered sixth in xparam1)
    jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx);
    k = estim_params_.endo_init_vals(jj,1);
    vname = sprintf('%s', M_.endo_names{k});
    nam = sprintf('INIT_%s', vname);
    if TeX
        tname  = sprintf('%s', M_.endo_names_tex{k});
        texnam = sprintf('$ init({%s}) $', tname);
    end
else % estimated structural parameters (ordered last in xparam1)
    jj = k - (estim_params_.nvx+estim_params_.nvn+estim_params_.ncx+estim_params_.ncn+estim_params_.nsx+estim_params_.nendoinit);
    jj1 = estim_params_.param_vals(jj,1);
    nam = M_.param_names{jj1};
    if TeX
        texnam = sprintf('$ %s $', M_.param_names_tex{jj1});
    end
end