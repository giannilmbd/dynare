function [PredictedStateMean, PredictedStateVarianceSquareRoot, StateVectorMean, StateVectorVarianceSquareRoot] = ...
    gaussian_filter_bank(ReducedForm, obs, StateVectorMean, StateVectorVarianceSquareRoot, Q_lower_triangular_cholesky, H_lower_triangular_cholesky, H, ...
                         ParticleOptions, ThreadsOptions, options_, M_)
% [PredictedStateMean, PredictedStateVarianceSquareRoot, StateVectorMean, StateVectorVarianceSquareRoot] = ...
%     gaussian_filter_bank(ReducedForm, obs, StateVectorMean, StateVectorVarianceSquareRoot, Q_lower_triangular_cholesky, H_lower_triangular_cholesky, H, ...
%                          ParticleOptions, ThreadsOptions, options_, M_)
%
% Computes the proposal with a Gaussian approximation for importance
% sampling. This proposal is a Gaussian distribution calculated à la Kalman
% Inputs
%  - ReducedForm                    [structure]    MATLAB's structure describing the reduced form model.
%  - obs                            [double]       p×1 vector of (detrended) data, where p is the number of observed variables.
%  - StateVectorMean                [double]       mean of the states
%  - StateVectorVarianceSquareRoot  [double]       square root of the state covariance matrix
%  - Q_lower_triangular_cholesky    [double]       Cholesky of shock covariance
%  - H_lower_triangular_cholesky    [double]       Cholesky of measurement error covariance
%  - H                              [double]       Measurement error covariance
%  - ParticleOptions                [structure]    filter options
%  - ThreadsOptions                 [structure]    options for threading of mex files
%  - options_                       [structure]    describing the options
%  - M_                             [structure]    describing the model
%
% Outputs
%  - PredictedStateMean                 [double]   one-step ahead predicted mean of the states    
%  - PredictedStateVarianceSquareRoot   [double]   one-step ahead covariance of the states   
%  - StateVectorMean                    [double]   updated mean of the states   
%  - StateVectorVarianceSquareRoot      [double]   updated covariance of the states   

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

if ParticleOptions.proposal_approximation.montecarlo
    nodes = randn(ParticleOptions.number_of_particles, number_of_state_variables+number_of_structural_innovations) ;
    weights = 1/ParticleOptions.number_of_particles ;
    weights_c = weights ;
elseif ParticleOptions.proposal_approximation.cubature 
    [nodes,weights] = spherical_radial_sigma_points(number_of_state_variables+number_of_structural_innovations) ;
    weights_c = weights ;
elseif ParticleOptions.proposal_approximation.unscented
    [nodes,weights,weights_c] = unscented_sigma_points(number_of_state_variables+number_of_structural_innovations, ParticleOptions);
else
    error('This approximation for the proposal is not implemented or unknown!')
end

xbar = [StateVectorMean ; zeros(number_of_structural_innovations,1)] ;
sqr_Px = [ StateVectorVarianceSquareRoot, zeros(number_of_state_variables, number_of_structural_innovations); 
           zeros(number_of_structural_innovations, number_of_state_variables) Q_lower_triangular_cholesky];
sigma_points = bsxfun(@plus, xbar, sqr_Px*(nodes'));
StateVectors = sigma_points(1:number_of_state_variables,:);
epsilon = sigma_points(number_of_state_variables+1:number_of_state_variables+number_of_structural_innovations,:);
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

PredictedStateMean = tmp(mf0,:)*weights;
PredictedObservedMean = tmp(mf1,:)*weights;

if ParticleOptions.proposal_approximation.cubature || ParticleOptions.proposal_approximation.montecarlo
    PredictedStateMean = sum(PredictedStateMean, 2);
    PredictedObservedMean = sum(PredictedObservedMean, 2);
    dState = bsxfun(@minus,tmp(mf0,:), PredictedStateMean)'.*sqrt(weights);
    dObserved = bsxfun(@minus, tmp(mf1,:), PredictedObservedMean)'.*sqrt(weights);
    PredictedStateVarianceSquareRoot = chol(dState'*dState)';
    big_mat = [dObserved, dState ; H_lower_triangular_cholesky, zeros(number_of_observed_variables,number_of_state_variables)];
    [~, mat] = qr2(big_mat, 0);
    mat = mat';
    PredictedObservedVarianceSquareRoot = mat(1:number_of_observed_variables,1:number_of_observed_variables);
    CovarianceObservedStateSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables),1:number_of_observed_variables);
    StateVectorVarianceSquareRoot = mat(number_of_observed_variables+(1:number_of_state_variables),number_of_observed_variables+(1:number_of_state_variables));
    PredictionError = obs - PredictedObservedMean;
    StateVectorMean = PredictedStateMean + (CovarianceObservedStateSquareRoot/PredictedObservedVarianceSquareRoot)*PredictionError;
else
    dState = bsxfun(@minus, tmp(mf0,:), PredictedStateMean);
    dObserved = bsxfun(@minus, tmp(mf1,:), PredictedObservedMean);
    PredictedStateVariance = dState*diag(weights_c)*dState';
    PredictedObservedVariance = dObserved*diag(weights_c)*dObserved' + H;
    PredictedStateAndObservedCovariance = dState*diag(weights_c)*dObserved';
    PredictedStateVarianceSquareRoot = chol(PredictedStateVariance)';
    PredictionError = obs - PredictedObservedMean;
    KalmanFilterGain = PredictedStateAndObservedCovariance/PredictedObservedVariance;
    StateVectorMean = PredictedStateMean + KalmanFilterGain*PredictionError;
    StateVectorVariance = PredictedStateVariance - KalmanFilterGain*PredictedObservedVariance*KalmanFilterGain';
    StateVectorVariance = .5*(StateVectorVariance+StateVectorVariance');
    StateVectorVarianceSquareRoot = reduced_rank_cholesky(StateVectorVariance)';
end
