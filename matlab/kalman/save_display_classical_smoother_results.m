function oo_=save_display_classical_smoother_results(xparam1,M_,oo_,options_,bayestopt_,dataset_,dataset_info,estim_params_)
% oo_=save_display_classical_smoother_results(xparam1,M_,oo_,options_,bayestopt_,dataset_,dataset_info,estim_params_)
% Inputs:
%   xparam1         [double]        current values for the estimated parameters.
%   M_              [structure]     storing the model information
%   oo_             [structure]     storing the results
%   options_        [structure]     storing the options
%   bayestopt_      [structure]     storing information about priors
%   dataset_        [structure]     storing the dataset
%   dataset_info    [structure]     information about the dataset (descriptive statistics and missing observations)
%   estim_params_   [structure]     storing information about estimated parameters
%
% Outputs:
%   oo_             [structure]     storing the results

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

smoother_error=false;
if options_.occbin.smoother.status && options_.occbin.smoother.inversion_filter
    [~, info, ~, ~, ~, ~, ~, ~, oo_.dr, atT, innov, oo_.occbin.smoother.regime_history] = occbin.IVF_posterior(xparam1,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,prior_bounds(bayestopt_,options_.prior_trunc),oo_.dr, oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state);
    if ismember(info(1),[303,304,306])
        fprintf('\nIVF: smoother did not succeed. No results will be written to oo_.\n')
        smoother_error=true;
    else
        updated_variables = atT*nan;
        measurement_error=[];
        ys = oo_.dr.ys;
        trend_coeff = zeros(length(options_.varobs_id),1);
        bayestopt_.mf = bayestopt_.smoother_var_list(bayestopt_.smoother_mf);
        options_nk=options_.nk;
        options_.nk=[]; %unset options_.nk and reset it later
        [oo_, yf]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff);
        options_.nk=options_nk;
    end
else
    if options_.occbin.smoother.status
        [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,~,~,P,PK,decomp,Trend,state_uncertainty,oo_,bayestopt_.mf,a0T,state_uncertainty0] = occbin.DSGE_smoother(xparam1,M_,oo_,options_,bayestopt_,estim_params_,dataset_,dataset_info,true);
        if oo_.occbin.smoother.error_flag(1)==0
            [oo_, yf] = store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty,a0T,state_uncertainty0);
        else
            smoother_error=true;
            fprintf('\nOccBin: smoother did not succeed. No results will be written to oo_.\n')
        end
    else
        [atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,~,~,P,PK,decomp,Trend,state_uncertainty,oo_.dr,bayestopt_.mf,a0T,state_uncertainty0] = DsgeSmoother(xparam1,dataset_,dataset_info,M_,oo_.dr,oo_.steady_state,oo_.exo_steady_state,oo_.exo_det_steady_state,options_,bayestopt_,estim_params_);
        [oo_,yf]=store_smoother_results(M_,oo_,options_,bayestopt_,dataset_,dataset_info,atT,innov,measurement_error,updated_variables,ys,trend_coeff,aK,P,PK,decomp,Trend,state_uncertainty,a0T,state_uncertainty0);
    end
end
if ~smoother_error
    plot_classical_smoother_results(M_,oo_,options_,dataset_info,dataset_,estim_params_,yf);
end