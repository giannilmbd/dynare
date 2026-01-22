function [lprior, base_draws, info]  = state_priordens(yhat,PriorStateInfo,number_of_particles,kalman_tol)
% [lprior, base_draws, info]  = occbin.ppf.state_priordens(yhat,PriorStateInfo,number_of_particles,kalman_tol)
%
% Draw particles in t-1|t from a mixture of normal distributions [up to as many mixtures as the number of OccBin regimes]
%
% INPUTS
%  - qmc_base               [double]    base normal draws (independent standard normal variates)
%  - PriorStateInfo         [struct]    info of particles distribution t-1|t-1
%  - number_of_particles    [integer]
%  - kalman_tol             [double]
%
% OUTPUTS
%  - lprior                 [double]    minus 2*log prior density of particles
%  - base_draws             [double]    standardised state values along SVD eigenvectors
%
% This function is called by: sequential_importance_particle_filter
% This function calls: residual_resampling, traditional_resampling

% Copyright © 2025-2026 Dynare Team
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

info=zeros(1,number_of_particles);
lprior=zeros(1,number_of_particles);
is_prior_mixture = size(PriorStateInfo.Variance,3)>1;
lprior_tmp = nan(size(PriorStateInfo.Variance,3),1);
base_draws  = nan(size(PriorStateInfo.UPR{1},2),number_of_particles);
for k=1:number_of_particles
    for r=1:size(PriorStateInfo.Variance,3)
        p0=1;
        if is_prior_mixture
            p0 = PriorStateInfo.proba(r);
        end
        if isempty(PriorStateInfo.UPR{r})
            lprior_tmp(r) = 0 - log(p0);
        else
            vv = PriorStateInfo.UPR{r}'*(yhat(:,k)-PriorStateInfo.Mean(:,r));
            lprior_tmp(r) = PriorStateInfo.log_dS(r) + transpose(vv)*PriorStateInfo.iS{r}*vv + length(PriorStateInfo.iS{r})*log(2*pi) - log(p0);
            if nargout==2
                base_draws(:,k) = PriorStateInfo.UPR{r}'*(yhat(:,k)-PriorStateInfo.Mean(:,r)).*sqrt(diag(PriorStateInfo.iS{r}));
            end
        end
        vvNULL = PriorStateInfo.UPNULL{r}'*(yhat(:,k)-PriorStateInfo.Mean(:,r));
        if any(vvNULL.^2>3*kalman_tol)
            lprior_tmp(r) = lprior_tmp(r)+100;
            if not(is_prior_mixture)
                info(k)=1;
            end
        end
    end
    if is_prior_mixture
        lprior(k) = -2*log(sum(exp(-(lprior_tmp-min(lprior_tmp))./2)))+min(lprior_tmp);
    else
        lprior(k) = lprior_tmp;
    end
end

