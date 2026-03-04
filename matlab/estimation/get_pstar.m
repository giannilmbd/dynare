function [Pstar, info] = get_pstar(xparam1,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [Pstar, info] = get_pstar(xparam1,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% Computes the initial state covariance Pstar used by the Kalman filter.
%
% Based on the chosen `options_.lik_init` strategy, this routine
% constructs the unconditional/state-initial covariance matrix for the
% state vector implied by the linearized model and shock covariance. It
% updates the model parameters from `xparam1`, resolves the model, and
% returns `Pstar` along with an `info` code (with `info(1)~=0` on error).
%
% INPUTS
% - xparam1             [double]        current values for the estimated parameters.
% - options_            [structure]     options controlling filter init via `lik_init`.
% - M_                  [structure]     model structure (updated via `set_all_parameters`).
% - estim_params_       [structure]     parameters to be estimated.
% - bayestopt_          [structure]     priors and measurement mapping (uses mf1/mf0).
% - BoundsInfo          [structure]     prior bounds and definiteness checks.
% - dr                  [structure]     reduced-form decision rules container.
% - endo_steady_state   [vector]        steady state for endogenous variables.
% - exo_steady_state    [vector]        steady state for exogenous variables.
% - exo_det_steady_state [vector]       steady state for deterministic exogenous variables.
%
% OUTPUTS
% - Pstar               [double]        initial state covariance matrix.
% - info                [double]        info vector; `info(1)~=0` indicates an error.
%
% NOTES ON `options_.lik_init` MODES
% - 0: Start from states in first_period-1 (uses `set_Kalman_starting_state`).
% - 1: Steady-state initialization (Lyapunov solution for covariance).
% - 2: Large diagonal covariance (Harvey scaling) for non-stationary models.
% - 3: Diffuse Kalman filter (Durbin/Koopman) — not used here to start the filter.
% - 4: Riccati steady state solution (falls back to 1 if Riccati fails).
% - 5: Old diffuse approach only for non-stationary variables.
%
% This function is called by: posterior_sampler_iteration

% Copyright © 2024-2026 Dynare Team
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

if ismember(options_.lik_init,[3,4])
    error('get_pstar: lik_init=3 or 4 are not yet supported')
end

% Initialization of the returned variables and others...
%------------------------------------------------------------------------------
% 1. Get the structural parameters & define penalties
%------------------------------------------------------------------------------
M_ = set_all_parameters(xparam1,estim_params_,M_);
Pstar=[];
[~,info,~,Q,H]=check_bounds_and_definiteness_estimation(xparam1, M_, estim_params_, BoundsInfo);
if info(1)
    return
end

%------------------------------------------------------------------------------
% 2. call model setup & reduction program
%------------------------------------------------------------------------------
% Linearize the model around the deterministic steady state and extract the matrices of the state equation (T and R).
[T,R,~,info,~, M_.params] = dynare_resolve(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,'restrict');

% Return if dynare_resolve issues an error code
if info(1)
    return
end

% check endogenous prior restrictions

%% Define a vector of indices for the observed variables. Is this really usefull?...
bayestopt_.mf = bayestopt_.mf1;

switch options_.lik_init
  case 0% start from states in first_period-1
    [~,Pstar]=set_Kalman_starting_state(M_);

  case 1% Standard initialization with the steady state of the state equation.
    Pstar=lyapunov_solver(T,R,Q,options_);

  case 2% Initialization with large numbers on the diagonal of the covariance matrix if the states (for non stationary models).
    Pstar = options_.Harvey_scale_factor*eye(mm);

  case 3% Diffuse Kalman filter (Durbin and Koopman)
        % Use standard kalman filter except if the univariate filter is explicitely choosen.
    [Pstar,~] = compute_Pinf_Pstar(Z,T,R,Q,options_.qz_criterium, dr.restrict_columns);

  case 4% Start from the solution of the Riccati equation.
    try
        if isequal(H,0)
            Pstar = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,mm,length(Z))));
        else
            Pstar = kalman_steady_state(transpose(T),R*Q*transpose(R),transpose(build_selection_matrix(Z,mm,length(Z))),H);
        end
    catch ME
        disp(ME.message)
        disp('dsge_likelihood:: I am not able to solve the Riccati equation, so I switch to lik_init=1!');
        options_.lik_init = 1;
        Pstar=lyapunov_solver(T,R,Q,options_);
    end
  case 5            % Old diffuse Kalman filter only for the non stationary variables
    [eigenvect, eigenv] = eig(T);
    eigenv = diag(eigenv);
    nstable = length(find(abs(abs(eigenv)-1) > 1e-7));
    V = eigenvect(:,abs(abs(eigenv)-1) < 1e-7);
    stable = find(sum(abs(V),2)<1e-5);
    nunit = length(eigenv) - nstable;
    Pstar = options_.Harvey_scale_factor*eye(nunit);
    R_tmp = R(stable, :);
    T_tmp = T(stable,stable);
    Pstar_tmp=lyapunov_solver(T_tmp,R_tmp,Q,options_);
    Pstar(stable, stable) = Pstar_tmp;
  otherwise
    error('dsge_likelihood:: Unknown initialization approach for the Kalman filter!')
end
