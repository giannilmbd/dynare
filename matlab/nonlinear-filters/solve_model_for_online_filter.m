function [info, M_, ReducedForm] = ...
    solve_model_for_online_filter(setinitialcondition, xparam1, dataset_, options_, M_, estim_params_, bayestopt_, bounds, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [info, M_, ReducedForm] = ...
%     solve_model_for_online_filter(setinitialcondition, xparam1, dataset_, options_, M_, estim_params_, bayestopt_, bounds, dr , endo_steady_state, exo_steady_state, exo_det_steady_state)

% Solves the DSGE model for an particular parameter set.
%
% INPUTS
% - setinitialcondition      [logical]    return initial condition if true.
% - xparam1                  [double]     n×1 vector, parameter values.
% - dataset_                 [struct]     Dataset for estimation.
% - options_                 [struct]     Dynare options.
% - M_                       [struct]     Model description.
% - estim_params_            [struct]     Estimated parameters.
% - bayestopt_               [struct]     Prior definition.
% - bounds                   [struct]     Prior bounds.
% - dr                       [structure]  Reduced form model.
% - endo_steady_state        [vector]     steady state value for endogenous variables
% - exo_steady_state         [vector]     steady state value for exogenous variables
% - exo_det_steady_state     [vector]     steady state value for exogenous deterministic variables
%
% OUTPUTS
% - info                     [integer]    scalar, nonzero if any problem occur when computing the reduced form.
% - M_                       [struct]     M_ description.
% - ReducedForm              [struct]     Reduced form model.

% Copyright © 2013-2025 Dynare Team
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

info        = zeros(4,1);

%----------------------------------------------------
% 1. Get the structural parameters & define penalties
%----------------------------------------------------

% Ensure that xparam1 is a column vector.
% (Don't do the transformation if xparam1 is empty, otherwise it would become a
%  0×1 matrix, which create issues with older MATLABs when comparing with [] in
%  check_bounds_and_definiteness_estimation)
if ~isempty(xparam1)
    xparam1 = xparam1(:);
end

M_ = set_all_parameters(xparam1,estim_params_,M_);

[~,info,~,Q,H]=check_bounds_and_definiteness_estimation(xparam1, M_, estim_params_, bounds);
if info(1)    
    return
end

% Update M_.Sigma_e and M_.H.
M_.Sigma_e = Q;
M_.H = H;

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------

warning('off', 'MATLAB:nearlySingularMatrix')
[dr, info, M_.params] = ...
    compute_decision_rules(M_, options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
warning('on', 'MATLAB:nearlySingularMatrix')

if info(1)~=0
    if nargout==3
        ReducedForm = 0;
    end
    return
end

% Set persistent variables (first call).
mf0 = bayestopt_.mf0;
mf1 = bayestopt_.mf1;
restrict_variables_idx  = dr.restrict_var_list;
state_variables_idx = restrict_variables_idx(mf0);
number_of_state_variables = length(mf0);


% Return reduced form model.
if nargout>2
    ReducedForm.ghx = dr.ghx(restrict_variables_idx,:);
    ReducedForm.ghu = dr.ghu(restrict_variables_idx,:);
    ReducedForm.steadystate = dr.ys(dr.order_var(restrict_variables_idx));
    if options_.order==2
        ReducedForm.use_k_order_solver = false;
        ReducedForm.ghxx = dr.ghxx(restrict_variables_idx,:);
        ReducedForm.ghuu = dr.ghuu(restrict_variables_idx,:);
        ReducedForm.ghxu = dr.ghxu(restrict_variables_idx,:);
        ReducedForm.constant = ReducedForm.steadystate + .5*dr.ghs2(restrict_variables_idx);
        ReducedForm.ghs2 = dr.ghs2(restrict_variables_idx,:);
    elseif options_.order>=3
        ReducedForm.use_k_order_solver = true;
        ReducedForm.dr = dr;
        ReducedForm.udr = folded_to_unfolded_dr(dr, M_, options_);
    else
        n_states=size(dr.ghx,2);
        n_shocks=size(dr.ghu,2);
        ReducedForm.use_k_order_solver = false;
        ReducedForm.ghxx = zeros(size(restrict_variables_idx,1),n_states^2);
        ReducedForm.ghuu = zeros(size(restrict_variables_idx,1),n_shocks^2);
        ReducedForm.ghxu = zeros(size(restrict_variables_idx,1),n_states*n_shocks);
        ReducedForm.constant = ReducedForm.steadystate;
    end
    ReducedForm.state_variables_steady_state = dr.ys(dr.order_var(state_variables_idx));
    ReducedForm.Q = Q;
    ReducedForm.H = H;
    ReducedForm.mf0 = mf0;
    ReducedForm.mf1 = mf1;
end

% Set initial condition
if setinitialcondition
    switch options_.particle.initialization
      case 1% Initial state vector covariance is the ergodic variance associated to the first order Taylor-approximation of the model.
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        [A,B] = kalman_transition_matrix(dr,dr.restrict_var_list,dr.restrict_columns);
        StateVectorVariance = lyapunov_symm(A, B*ReducedForm.Q*B', options_.lyapunov_fixed_point_tol, ...
                                        options_.qz_criterium, options_.lyapunov_complex_threshold, [], options_.debug);
        StateVectorVariance = StateVectorVariance(mf0,mf0);
      case 2% Initial state vector covariance is a Monte-Carlo based estimate of the ergodic variance (consistent with a k-order Taylor-approximation of the model).
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        options_.periods = 5000;
        options_.pruning = options_.particle.pruning;
        y_ = simult(dr.ys, dr, M_, options_);
        y_ = y_(dr.order_var(state_variables_idx),2001:options_.periods);
        StateVectorVariance = cov(y_');
      case 3% Initial state vector covariance is a diagonal matrix.
        StateVectorMean = ReducedForm.state_variables_steady_state;%.constant(mf0);
        StateVectorVariance = options_.particle.initial_state_prior_std*eye(number_of_state_variables);
      otherwise
        error('Unknown initialization option!')
    end
    ReducedForm.StateVectorMean = StateVectorMean;
    ReducedForm.StateVectorVariance = StateVectorVariance;
end
