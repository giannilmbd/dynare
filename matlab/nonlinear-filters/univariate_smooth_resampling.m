function new_particles = univariate_smooth_resampling(weights,particles,number_of_new_particles)
% new_particles = univariate_smooth_resampling(weights,particles,number_of_new_particles)
% Smooth Resampling of the particles (univariate version).
%
% INPUTS
%  - weights                    [double]    n*1 vector of particles' weights.
%  - particles                  [double]    n*1 vector of particles
%  - number_of_new_particles    [integer]   number of new particles
%
% OUTPUTS
%  - new_particles              [double]    number_of_new_particles*1 vector of doubles, new particles
%
% This function is called by: sequential_importance_particle_filter

% Copyright © 2012-2025 Dynare Team
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

% AUTHOR(S) frederic DOT karame AT univ DASH lemans DOT fr
%           stephane DOT adjemian AT univ DASH lemans DOT fr

M = length(particles) ;
lambda_tilde = [ (.5*weights(1)) ;
                 (.5*(weights(1:M-1)+weights(2:M))) ;
                 (.5*weights(M)) ] ;
lambda_bar = cumsum(lambda_tilde) ;
u = rand(1,1) ;
new_particles = zeros(number_of_new_particles,1) ;
rj = 0 ;
i = 1 ;
j = 1 ;
while i<=number_of_new_particles
    u_j = ( i-1 + u)/number_of_new_particles ;
    while u_j>lambda_bar(j)
        rj = j ;
        j = j+1 ;
    end
    if rj==0
        new_particles(i) = particles(1) ;
    elseif rj==M
        new_particles(i) = particles(M) ;
    else
        u_star = (u_j - lambda_bar(rj))./lambda_tilde(rj+1) ;
        new_particles(i) = (particles(rj+1) - particles(rj))*u_star + particles(rj) ;
    end
    i = i+1 ;
end