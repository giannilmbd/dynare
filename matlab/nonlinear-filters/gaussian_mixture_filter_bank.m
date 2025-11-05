function [StateMuPrior,StateSqrtPPrior,StateWeightsPrior,StateMuPost,StateSqrtPPost,StateWeightsPost] =...
    gaussian_mixture_filter_bank(ReducedForm, obs, StateMu, StateSqrtP, StateWeights, ...
                                 StructuralShocksMu, StructuralShocksSqrtP, StructuralShocksWeights, ...
                                 ObservationShocksWeights, H, H_lower_triangular_cholesky, normfactO, ...
                                 ParticleOptions, ThreadsOptions, options_, M_)
% [StateMuPrior,StateSqrtPPrior,StateWeightsPrior,StateMuPost,StateSqrtPPost,StateWeightsPost] =...
%     gaussian_mixture_filter_bank(ReducedForm, obs, StateMu, StateSqrtP, StateWeights, ...
%                                  StructuralShocksMu, StructuralShocksSqrtP, StructuralShocksWeights, ...
%                                  ObservationShocksWeights, H, H_lower_triangular_cholesky, normfactO, ...
%                                  ParticleOptions, ThreadsOptions, options_, M_)
%
% Computes the proposal with a Gaussian approximation for importance
% sampling. This proposal is a gaussian distribution calculated à la Kalman
%
% INPUTS
% - ReducedForm             [structure]     MATLAB's structure describing the reduced form model.
% - obs                     [double]        pp*1 vector of (detrended) data, where pp is the maximum number of observed variables
% - StateMu                 [double]        mean of the states
% - StateSqrtP              [double]        square root of the state covariance matrix
% - StateWeights            [double]        weights of the state particles
% - StructuralShocksMu      [double]        mean of the structural shocks
% - StructuralShocksSqrtP	[double]        square root of covariance matrix of the structural shocks    
% - StructuralShocksWeights	[double]        weights of structural shocks
% - ObservationShocksWeights [double]       weights of measurement errors
% - H                              [double]       Measurement error covariance
% - H_lower_triangular_cholesky    [double]       Cholesky of measurement error covariance
% - normfactO               [double]        normalizing constant in likelihood
% - ParticleOptions         [structure]     filter options
% - ThreadsOptions          [structure]     options for threading of mex files
% - options_                [structure]     describing the options
% - M_                      [structure]     describing the model
%
% OUTPUTS
% - StateMuPrior            [double]        prior mean of the states
% - StateSqrtPPrior         [double]        square root of prior covariance of the states
% - StateWeightsPrior       [double]        prior weight of the states
% - StateMuPost             [double]        posterior mean of the states
% - StateSqrtPPost          [double]        square root of posterior covariance of the states
% - StateWeightsPost        [double]        posterior weight of the states

% Copyright © 2009-2025 Dynare Team
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

mf0 = ReducedForm.mf0;
mf1 = ReducedForm.mf1;
number_of_state_variables = length(mf0);
number_of_observed_variables = length(mf1);
number_of_structural_innovations = length(ReducedForm.Q);

numb = number_of_state_variables+number_of_structural_innovations;

if ParticleOptions.proposal_approximation.cubature
    [nodes3, weights3] = spherical_radial_sigma_points(numb);
    weights_c3 = weights3;
elseif ParticleOptions.proposal_approximation.unscented
    [nodes3, weights3, weights_c3] = unscented_sigma_points(numb, ParticleOptions);
else
    error('This approximation for the proposal is unknown!')
end

epsilon = bsxfun(@plus, StructuralShocksSqrtP*nodes3(:,number_of_state_variables+1:number_of_state_variables+number_of_structural_innovations)', StructuralShocksMu);
StateVectors = bsxfun(@plus, StateSqrtP*nodes3(:,1:number_of_state_variables)', StateMu);
yhat = bsxfun(@minus, StateVectors, ReducedForm.state_variables_steady_state);
if ReducedForm.use_k_order_solver
    tmp = local_state_space_iteration_k(yhat, epsilon, ReducedForm.dr, M_, options_, ReducedForm.udr);
else
    if options_.order == 2
        tmp = local_state_space_iteration_2(yhat, epsilon, ReducedForm.ghx, ReducedForm.ghu, ReducedForm.constant, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ThreadsOptions.local_state_space_iteration_2);
    elseif options_.order == 3
        tmp = local_state_space_iteration_3(yhat, epsilon, ReducedForm.ghx, ReducedForm.ghu, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ReducedForm.ghs2, ReducedForm.ghxxx, ReducedForm.ghuuu, ReducedForm.ghxxu, ReducedForm.ghxuu, ReducedForm.ghxss, ReducedForm.ghuss, ReducedForm.steadystate, ThreadsOptions.local_state_space_iteration_3, false);
    else
        error('Order > 3: use_k_order_solver should be set to true');
    end
end
PredictedStateMean = tmp(mf0,:)*weights3;
PredictedObservedMean = tmp(mf1,:)*weights3;

if ParticleOptions.proposal_approximation.cubature
    PredictedStateMean = sum(PredictedStateMean, 2);
    PredictedObservedMean = sum(PredictedObservedMean, 2);
    dState = (bsxfun(@minus, tmp(mf0,:), PredictedStateMean)').*sqrt(weights3);
    dObserved = (bsxfun(@minus, tmp(mf1,:), PredictedObservedMean)').*sqrt(weights3);
    PredictedStateVariance = dState'*dState;
    big_mat = [dObserved,  dState ; H_lower_triangular_cholesky, zeros(number_of_observed_variables, number_of_state_variables)];
    [~, mat] = qr2(big_mat, 0);
    mat = mat';
    PredictedObservedVarianceSquareRoot = mat(1:number_of_observed_variables, 1:number_of_observed_variables);
    CovarianceObservedStateSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables), 1:number_of_observed_variables);
    StateVectorVarianceSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables), number_of_observed_variables+(1:number_of_state_variables));
    iPredictedObservedVarianceSquareRoot = inv(PredictedObservedVarianceSquareRoot);
    iPredictedObservedVariance = iPredictedObservedVarianceSquareRoot'*iPredictedObservedVarianceSquareRoot;
    sqrdet = 1/sqrt(det(iPredictedObservedVariance));
    PredictionError = obs - PredictedObservedMean;
    StateVectorMean = PredictedStateMean + CovarianceObservedStateSquareRoot*iPredictedObservedVarianceSquareRoot*PredictionError;
else
    dState = bsxfun(@minus, tmp(mf0,:), PredictedStateMean);
    dObserved = bsxfun(@minus, tmp(mf1,:), PredictedObservedMean);
    PredictedStateVariance = dState*diag(weights_c3)*dState';
    PredictedObservedVariance = dObserved*diag(weights_c3)*dObserved' + H;
    PredictedStateAndObservedCovariance = dState*diag(weights_c3)*dObserved';
    sqrdet = sqrt(det(PredictedObservedVariance));
    iPredictedObservedVariance = inv(PredictedObservedVariance);
    PredictionError = obs - PredictedObservedMean;
    KalmanFilterGain = PredictedStateAndObservedCovariance*iPredictedObservedVariance;
    StateVectorMean = PredictedStateMean + KalmanFilterGain*PredictionError;
    StateVectorVariance = PredictedStateVariance - KalmanFilterGain*PredictedObservedVariance*KalmanFilterGain';
    StateVectorVariance = .5*(StateVectorVariance+StateVectorVariance');
    StateVectorVarianceSquareRoot = reduced_rank_cholesky(StateVectorVariance)';
end

data_lik_GM_g = exp(-0.5*PredictionError'*iPredictedObservedVariance*PredictionError)/abs(normfactO*sqrdet) + 1e-99;
StateMuPrior = PredictedStateMean;
StateSqrtPPrior = reduced_rank_cholesky(PredictedStateVariance)';
StateWeightsPrior = StateWeights*StructuralShocksWeights;
StateMuPost = StateVectorMean;
StateSqrtPPost = StateVectorVarianceSquareRoot;
StateWeightsPost = StateWeightsPrior*ObservationShocksWeights*data_lik_GM_g;
