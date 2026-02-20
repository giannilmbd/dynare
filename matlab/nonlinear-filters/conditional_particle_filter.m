function [LIK,lik] = conditional_particle_filter(ReducedForm, Y, s, ParticleOptions, ThreadsOptions, options_, M_)
% [LIK,lik] = conditional_particle_filter(ReducedForm, Y, s, ParticleOptions, ThreadsOptions, options_, M_)
% Evaluates the likelihood of a non-linear model with a particle filter
%
% INPUTS
% - ReducedForm        [structure]    MATLAB's structure describing the reduced form model.
% - Y                  [double]       p×T matrix of (detrended) data, where p is the number of observed variables.
% - s                  [integer]      scalar, likelihood evaluation starts at s (has to be smaller than T, the sample length provided in Y).
% - ParticleOptions    [structure]    filter options
% - ThreadsOptions     [structure]    options for threading of mex files
% - options_           [structure]    describing the options
% - M_                 [structure]    describing the model
%
% OUTPUTS
% - LIK                [double]        scalar, likelihood
% - lik                [double]        (T-s+1)×1 vector, density of observations in each period.
%
% REMARKS
% - The proposal is built using the Kalman updating step for each particle.
% - We need draws of  the errors distributions
%       - Whether we use Monte-Carlo draws from a multivariate gaussian distribution
%         as in Amisano & Tristani (JEDC 2010).
%       - Whether we use multidimensional Gaussian sparse grids approximations:
%           - a univariate Kronrod-Paterson Gaussian quadrature combined by the Smolyak
%               operator (ref: Winschel & Kratzig, 2010).
%           - a spherical-radial cubature (ref: Arasaratnam & Haykin, 2009a,2009b).
%           - a scaled unscented transform cubature (ref: Julier & Uhlmann 1997, van der
%               Merwe & Wan 2003).
%
% Pros:
% - Allows using current observable information in the proposal
% - The use of sparse grids Gaussian approximation is much faster than the Monte-Carlo approach
% Cons:
% - The use of the Kalman updating step may bias the proposal distribution since
%   it has been derived in a linear context and is implemented in a nonlinear
%   context. That is why particle resampling is performed.

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

if ParticleOptions.pruning
    error('conditional_particle_filter: pruning is not yet implemented.')
end
% Set default for third input argument.
if isempty(s)
    s = 1;
end

T = size(Y,2);
n = ParticleOptions.number_of_particles;

% Get covariance matrices
Q = ReducedForm.Q;
H = ReducedForm.H;
if isempty(H)
    H = 0;
    H_lower_triangular_cholesky = 0;
else
    H_lower_triangular_cholesky = chol(H)';
end

% Get initial condition for the state vector.
StateVectorMean = ReducedForm.StateVectorMean;
StateVectorVarianceSquareRoot = chol(ReducedForm.StateVectorVariance)';
state_variance_rank = size(StateVectorVarianceSquareRoot, 2);
Q_lower_triangular_cholesky = chol(Q)';

% Set seed for randn().
set_dynare_seed_local_options([],false,'default');

% Initialization of the likelihood.
lik = NaN(T, 1);
ks = 0;
StateParticles = bsxfun(@plus, StateVectorVarianceSquareRoot*randn(state_variance_rank, n), StateVectorMean);
SampleWeights = ones(1, n)/n;

for t=1:T
    flags = false(n, 1);
    for i=1:n
        [StateParticles(:,i), SampleWeights(i), flags(i)] = ...
            conditional_filter_proposal(ReducedForm, Y(:,t), StateParticles(:,i), SampleWeights(i), Q_lower_triangular_cholesky, H_lower_triangular_cholesky, H, ParticleOptions, options_, M_);
    end
    if any(flags)
        LIK = -Inf;
        lik(t) = -Inf;
        return
    end
    SumSampleWeights = sum(SampleWeights);
    lik(t) = log(SumSampleWeights);
    SampleWeights = SampleWeights./SumSampleWeights;
    if (ParticleOptions.resampling.status.generic && neff(SampleWeights)<ParticleOptions.resampling.threshold*T) || ParticleOptions.resampling.status.systematic
        ks = ks + 1;
        StateParticles = resample(StateParticles', SampleWeights', ParticleOptions)';
        SampleWeights = ones(1, n)/n;
    end
end

LIK = -sum(lik(s:end));

function [ProposalStateVector, Weights, flag] = conditional_filter_proposal(ReducedForm, y, StateVectors, SampleWeights, Q_lower_triangular_cholesky, H_lower_triangular_cholesky, ...
                                                  H, ParticleOptions, options_, M_)
% Computes the proposal for each past particle using Gaussian approximations
% for the state errors and the Kalman filter
%
% INPUTS
% - ReducedForm                    [structure]    MATLAB's structure describing the reduced form model.
% - y                              [double]       p×1 vector, current observation (p is the number of observed variables).
% - StateVectors
% - SampleWeights
% - Q_lower_triangular_cholesky
% - H_lower_triangular_cholesky
% - H
% - ParticleOptions
% - options_
% - M_
%
% OUTPUTS
% - ProposalStateVector
% - Weights
% - flag

flag = false;

mf0 = ReducedForm.mf0;
mf1 = ReducedForm.mf1;
number_of_state_variables = length(mf0);
number_of_observed_variables = length(mf1);
number_of_structural_innovations = length(ReducedForm.Q);

if ParticleOptions.proposal_approximation.montecarlo
    nodes = randn(ParticleOptions.number_of_particles/10, number_of_structural_innovations);
    weights = 1/(ParticleOptions.number_of_particles/10);
    weights_c = weights;
elseif ParticleOptions.proposal_approximation.cubature
    [nodes, weights] = spherical_radial_sigma_points(number_of_structural_innovations);
    weights_c = weights;
elseif ParticleOptions.proposal_approximation.unscented
    [nodes, weights, weights_c] = unscented_sigma_points(number_of_structural_innovations, ParticleOptions.unscented);
else
    error('Estimation: This approximation for the proposal is not implemented or unknown!')
end

epsilon = Q_lower_triangular_cholesky*nodes';
tmp=iterate_law_of_motion(StateVectors,epsilon,ReducedForm,M_,options_,ReducedForm.use_k_order_solver,options_.pruning);

PredictedStateMean = tmp(mf0,:)*weights;
PredictedObservedMean = tmp(mf1,:)*weights;

if ParticleOptions.proposal_approximation.cubature || ParticleOptions.proposal_approximation.montecarlo
    PredictedStateMean = sum(PredictedStateMean, 2);
    PredictedObservedMean = sum(PredictedObservedMean, 2);
    dState = bsxfun(@minus, tmp(mf0,:), PredictedStateMean)'.*sqrt(weights);
    dObserved = bsxfun(@minus, tmp(mf1,:), PredictedObservedMean)'.*sqrt(weights);
    PredictedStateVariance = dState*dState';
    big_mat = [dObserved dState; H_lower_triangular_cholesky zeros(number_of_observed_variables,number_of_state_variables)];
    [~, mat] = qr2(big_mat,0);
    mat = mat';
    PredictedObservedVarianceSquareRoot = mat(1:number_of_observed_variables, 1:number_of_observed_variables);
    CovarianceObservedStateSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables),1:number_of_observed_variables);
    StateVectorVarianceSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables),number_of_observed_variables+(1:number_of_state_variables));
    Error = y-PredictedObservedMean;
    StateVectorMean = PredictedStateMean+(CovarianceObservedStateSquareRoot/PredictedObservedVarianceSquareRoot)*Error;
    if ParticleOptions.cpf_weights_method.amisanotristani
        Weights = SampleWeights.*probability2(zeros(number_of_observed_variables,1), PredictedObservedVarianceSquareRoot, Error);
    end
else
    dState = bsxfun(@minus, tmp(mf0,:), PredictedStateMean);
    dObserved = bsxfun(@minus, tmp(mf1,:), PredictedObservedMean);
    PredictedStateVariance = dState*diag(weights_c)*dState';
    PredictedObservedVariance = dObserved*diag(weights_c)*dObserved'+H;
    PredictedStateAndObservedCovariance = dState*diag(weights_c)*dObserved';
    KalmanFilterGain = PredictedStateAndObservedCovariance/PredictedObservedVariance;
    Error = y-PredictedObservedMean;
    StateVectorMean = PredictedStateMean+KalmanFilterGain*Error;
    StateVectorVariance = PredictedStateVariance-KalmanFilterGain*PredictedObservedVariance*KalmanFilterGain';
    StateVectorVariance = 0.5*(StateVectorVariance+StateVectorVariance');
    [StateVectorVarianceSquareRoot, p] = chol(StateVectorVariance, 'lower') ;
    if p
        flag = true;
        ProposalStateVector = zeros(number_of_state_variables, 1);
        Weights = 0.0;
        return
    end
    if ParticleOptions.cpf_weights_method.amisanotristani
        Weights = SampleWeights.*probability2(zeros(number_of_observed_variables, 1), chol(PredictedObservedVariance)', Error);
    end
end

ProposalStateVector = StateVectorVarianceSquareRoot*randn(size(StateVectorVarianceSquareRoot, 2), 1)+StateVectorMean;
if ParticleOptions.cpf_weights_method.murrayjonesparslow
    PredictedStateVariance = 0.5*(PredictedStateVariance+PredictedStateVariance');
    [PredictedStateVarianceSquareRoot, p] = chol(PredictedStateVariance, 'lower');
    if p
        flag = true;
        ProposalStateVector = zeros(number_of_state_variables,1);
        Weights = 0.0;
        return
    end
    Prior = probability2(PredictedStateMean, PredictedStateVarianceSquareRoot, ProposalStateVector);
    Posterior = probability2(StateVectorMean, StateVectorVarianceSquareRoot, ProposalStateVector);
    Likelihood = probability2(y, H_lower_triangular_cholesky, measurement_equations(ProposalStateVector, ReducedForm, options_, M_));
    Weights = SampleWeights.*Likelihood.*(Prior./Posterior);
end
