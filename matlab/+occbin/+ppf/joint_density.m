function [ldens, infox] = joint_density(x, StateInfo, PriorStateInfo,options_,occbin_options, likxmode, bayestopt_, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes, ...
    a0, P, P0, P1,my_data_index,Z,ZZ,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
    M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state)
% [ldens, infox] = occbin.ppf.joint_density(x, StateInfo, PriorStateInfo,options_,occbin_options, likxmode, bayestopt_, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes, ...
%     a0, P, P0, P1,my_data_index,Z,ZZ,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
%     M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state);
%
% INPUTS
%  - x                      [double]    standard normal draws
%  - StateInfo              [struct]    info of particles sampling distribution(s) t-1|t
%  - PriorStateInfo         [struct]    info of particles distribution t-1|t-1
%  - options_               [structure] MATLAB's structure containing the options
%  - occbin_options         [structure] options structure for OccBin filter
%  - likxmode               [double]    mode of log-likelihood
%  - bayestopt_             [structure] MATLAB's structure containing Bayesian estimation options
%  - updated_regimes        [struct]    occbin info about updated regimes (particles)
%  - updated_sample         [struct]    occbin info about updated state (particles)
%  - number_of_updated_regimes [integer] number of updated regimes across particles
%  - all_updated_regimes    [cell]      list of updated regimes
%  - a0                     [double]    PKF updated state estimate t-1.
%  - P                      [double]    Particle state covariance t-1 (=0).
%  - P0                     [double]    PKF updated state covariance t-1.
%  - P1                     [double]    Particle 1 step ahead state covariance.
%  - my_data_index          [cell]      1*2 cell of column vectors of indices.
%  - Z                      [double]    Observation matrix.
%  - ZZ                     [double]    Observation selection matrix.
%  - v                      [double]    Observation innovations.
%  - my_data                [double]    observed data t-1:t
%  - H                      [double]    Measurement error variance.
%  - QQQ                    [double]    Shocks covariance matrix.
%  - T0                     [double]    Base regime state transition matrix.
%  - R0                     [double]    Base regime state shock matrix.
%  - TT                     [double]    PKF Transition matrices t-1:t given t-1 info.
%  - RR                     [double]    PKF Shock matrices t-1:t given t-1 info.
%  - CC                     [double]    PKF Constant terms in transition t-1:t given t-1 info.
%  - regimesy               [struct]    PKF updated regime info t:t+2.
%  - regimes0               [struct]    PKF regime info t:t+1 given t-1 info.
%  - base_regime            [integer]   Base regime info.
%  - info                   [integer]   PKF error flag.
%  - M_                     [structure] MATLAB's structure describing the model
%  - dr                     [structure] model information structure
%  - endo_steady_state      [vector]    steady state value for endogenous variables
%  - exo_steady_state       [vector]    steady state value for exogenous variables
%  - exo_det_steady_state   [vector]    steady state value for exogenous deterministic variables
%
% OUTPUTS
%  - ldens                  [double]    -logdensity states and data
%  - infox                  [integer]   indicator of successful computation
%
% This function is called by: ppf.mmcmc_draws
% This function calls: ppf.draw_particles, ppf.conditional_data_density, ppf.state_priordens

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

% we enter with standard normal values
yhat = occbin.ppf.draw_particles(transpose(x),StateInfo,PriorStateInfo,1,options_.kalman_tol);

[likxc, infox] = ...
    occbin.ppf.conditional_data_density(yhat, 1, 2, likxmode, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes, ...
    a0, P, P0, P1,my_data_index,Z,ZZ,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
    M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);

ldens = 1.e8;
if infox==0
    lprior  = occbin.ppf.state_priordens(yhat,PriorStateInfo,1,options_.kalman_tol);
    ldens = +likxc/2+lprior/2;
end