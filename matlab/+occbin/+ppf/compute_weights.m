function [weights, ESS] = compute_weights(lnw0, smc_scale, number_of_particles, number_of_successful_particles, success)
% [weights, ESS] = occbin.ppf.compute_weights(lnw0, smc_scale, number_of_particles, number_of_successful_particles, success)
% INPUTS
%  - lnw0                   [double]    vector of posterior [1 x number_of_particles]
%  - smc_scale              [double]    scale factor used in SMC
%  - number_of_particles    [integer]
%  - number_of_successful_particles    [integer]
%  - success                [logical]   indicator array of successful particles
%
% OUTPUTS
%  - weights                [double]    vector of normalized weights [1 x number_of_particles]
%  - ESS                    [double]    effective sample size
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

lnw0 = lnw0*smc_scale;
% we need to normalize weights, even if they are already proper likelihood values!
% to avoid that exponential gives ZEROS all over the place!
lnw0 = lnw0-max(lnw0);
weights = zeros(1,number_of_particles);
weights0 = ones(1,number_of_successful_particles)/number_of_successful_particles ;
wtilde0 = weights0.*exp(lnw0);
weights0 = wtilde0/sum(wtilde0);
weights(success) = weights0;
ESS = 1/sum(weights.^2);
