function resampled_output = resample(particles,weights,ParticleOptions)
% resampled_output = resample(particles,weights,ParticleOptions)
% Resamples particles.
% if particles = 0, returns the resampling index (except for smooth resampling)
% Otherwise, returns the resampled particles set.
%
% INPUTS
%  - particles              [double]    n*1 vector of particles
%  - weights                [double]    n*1 vector of particles' weights.
%  - ParticleOptions        [structure] filter options
%
% OUTPUTS
%  - resampled_output       [double]    vector of resample particles
%
% This function is called by: sequential_importance_particle_filter
% This function calls: residual_resampling, traditional_resampling

% Copyright © 2011-2025 Dynare Team
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

defaultmethod = 1; % For residual based method set this variable equal to 0.

if defaultmethod
    if ParticleOptions.resampling.method.kitagawa
        resampled_output = traditional_resampling(particles,weights,rand);
    elseif ParticleOptions.resampling.method.stratified
        resampled_output = traditional_resampling(particles,weights,rand(size(weights)));
    elseif ParticleOptions.resampling.method.smooth
        if particles==0
            error('Particle = 0 is incompatible with the resampling_method=smooth method!')
        end
        resampled_output = multivariate_smooth_resampling(particles,weights);
    else
        error('Unknown sampling method!')
    end
else
    if ParticleOptions.resampling.method.kitagawa
        resampled_output = residual_resampling(particles,weights,rand);
    elseif ParticleOptions.resampling.method.stratified
        resampled_output = residual_resampling(particles,weights,rand(size(weights)));
    else
        error('Unknown sampling method!')
    end
end
