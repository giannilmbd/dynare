function  IncrementalWeights = gaussian_mixture_densities(obs, StateMuPrior, StateSqrtPPrior, StateWeightsPrior, ...
                                                      StateMuPost, StateSqrtPPost, StateWeightsPost, StateParticles, H, ...
                                                      ReducedForm, options_, M_)
% IncrementalWeights = gaussian_mixture_densities(obs, StateMuPrior, StateSqrtPPrior, StateWeightsPrior, ...
%                                                       StateMuPost, StateSqrtPPost, StateWeightsPost, StateParticles, H, ...
%                                                       ReducedForm, options_, M_)
% Elements to calculate the importance sampling ratio
%
% INPUTS
%  - obs                 [double]   current observation
%  - StateMuPrior        [double]   prior mean
%  - StateSqrtPPrior     [double]   prior covariance
%  - StateWeightsPrior   [double]   prior covariance
%  - StateMuPost         [double]   posterior mean
%  - StateSqrtPPost      [double]   posterior covariance
%  - StateWeightsPost    [double]   weights of the particles
%  - StateParticles      [double]   particles
%  - ReducedForm         [structure] decision rules
%  - options_            [structure] describing the options
%  - M_                  [structure] describing the model
%
% OUTPUTS
%    IncrementalWeights  [double]   updated weights

% Copyright © 2009-2026 Dynare Team
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

% Compute the density of particles under the prior distribution
[~, ~, prior] = probability(StateMuPrior, StateSqrtPPrior, StateWeightsPrior, StateParticles);
prior = prior';

% Compute the density of particles under the proposal distribution
[~, ~, proposal] = probability(StateMuPost, StateSqrtPPost, StateWeightsPost, StateParticles);
proposal = proposal';

% Compute the density of the current observation conditionally to each particle
yt_t_1_i = measurement_equations(StateParticles, ReducedForm, options_, M_);

% likelihood
likelihood = probability2(obs, sqrt(H), yt_t_1_i);
IncrementalWeights = likelihood.*prior./proposal;
