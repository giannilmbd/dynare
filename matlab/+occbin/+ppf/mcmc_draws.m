function [x, fsim, feval] = mcmc_draws(x0,niter,StateInfo,PriorStateInfo, ...
    a0, P, P0, P1,my_data_index,Z,ZZ,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
    M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
% [x, fsim, feval] = occbin.ppf.mcmc_draws(x0,niter,StateInfo,PriorStateInfo, ...
%     P, P1,my_data_index,Z,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
%     M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
%
% Performs MCMC sampling of the joint state and data density in the OccBin framework.
%
% INPUTS
% - x0                      [double]    Initial standard normal draw.
% - niter                   [integer]   Number of MCMC iterations.
% - StateInfo               [struct]    Posterior state density information.
% - PriorStateInfo          [struct]    Prior state density information.
% - a0                      [double]    PKF updated state estimate t-1.
% - P                       [double]    Particle state covariance t-1 (=0).
% - P0                      [double]    PKF updated state covariance t-1.
% - P1                      [double]    Particle 1 step ahead state covariance.
% - my_data_index           [cell]      1*2 cell of column vectors of indices.
% - Z                       [double]    Observation matrix.
% - ZZ                      [double]    Observation selection matrix.
% - v                       [double]    Observation innovations.
% - Y                       [double]    Observations t-1:t.
% - H                       [double]    Measurement error variance.
% - QQQ                     [double]    Shocks covariance matrix.
% - T0                      [double]    Base regime state transition matrix.
% - R0                      [double]    Base regime state shock matrix.
% - TT                      [double]    PKF Transition matrices t-1:t given t-1 info.
% - RR                      [double]    PKF Shock matrices t-1:t given t-1 info.
% - CC                      [double]    PKF Constant terms in transition t-1:t given t-1 info.
% - regimesy                [struct]    PKF updated regime info t:t+2.
% - regimes0                [struct]    PKF regime info t:t+1 given t-1 info.
% - base_regime             [integer]   Base regime info.
% - info                    [integer]   PKF error flag.
% - M_                      [struct]    Dynare's model structure.
% - dr                      [struct]    Decision rule structure.
% - endo_steady_state       [double]    Endogenous steady state.
% - exo_steady_state        [double]    Exogenous steady state.
% - exo_det_steady_state    [double]    Exogenous deterministic steady state.
% - options_                [struct]    Dynare options.
% - occbin_options          [struct]    OccBin specific options.
%
% OUTPUTS
% - x                       [double]    Matrix of MCMC draws (ns x niter).
% - fsim                    [double]    Log-density values for each draw.
% - feval                   [integer]   Total number of density function evaluations.

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


priortrunc = 1.e-10;
lb = norminv(priortrunc, 0, 1);
ub = norminv(1-priortrunc, 0, 1);

% set default for slice!
sampler_options.proposal_distribution = '';
sampler_options.rotated=0;
sampler_options.slice_initialize_with_mode=0;
sampler_options.use_mh_covariance_matrix=0;
sampler_options.WR=[];
sampler_options.mode_files=[];
sampler_options.mode=[];
sampler_options.maximize=false;
sampler_options.maximize_using_mh_bounds=false;
sampler_options.mode_compute=5;
sampler_options.initial_step_size=0.8;
sampler_options.use_prior_draws.status=false;
sampler_options.use_prior_draws.mh_blck = [];
sampler_options.save_tmp_file=false;
sampler_options.save_iter_info_file=false;
sampler_options.draw_init_state_from_smoother=false;
sampler_options.draw_init_state_with_rotated_slice=false;
sampler_options.fast_likelihood_evaluation_for_rejection=false;
sampler_options.fast_likelihood_evaluation_for_rejection_penalty=10;

theta = x0;
ns = length(x0);
x = nan(ns,niter);
fsim = nan(1,niter);
feval = 0;
thetaprior = repmat([lb ub],ns,1);
sampler_options.W1=sampler_options.initial_step_size*(thetaprior(:,2)-thetaprior(:,1));
sampler_options.curr_block = 0;

bayestopt_.pshape = ones(ns,1)*3; % Gaussian distributions
bayestopt_.p1 = zeros(ns,1); % prior mean
bayestopt_.p2 = ones(ns,1); % prior standard deviation
bayestopt_.p3 = -inf(ns,1); % prior standard deviation
bayestopt_.p4 = inf(ns,1); % prior standard deviation
bayestopt_.p5 = zeros(ns,1); % prior mode
bayestopt_.p6 = zeros(ns,1); % prior mean
bayestopt_.p7 = ones(ns,1); % prior standard deviation
for p=1:ns
    bayestopt_.name{p} = ['UP' int2str(p)];
end

% varargin{6} must be bayestopt_, varargin{3} must be options_
varargin = {StateInfo,PriorStateInfo, options_,occbin_options, inf, bayestopt_, [], [], 0, {}, ...
    a0, P, P0, P1,my_data_index,Z,ZZ,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
    M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state};


for k=1:niter
    [theta, fxsim, neval, sampler_options] = slice_sampler('occbin.ppf.joint_density',theta,thetaprior,sampler_options,varargin{:});
    x(:,k) = theta;
    fsim(1,k) = fxsim;
    feval = feval+neval;
end