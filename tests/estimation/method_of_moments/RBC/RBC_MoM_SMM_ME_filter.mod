% =========================================================================
% Copyright © 2020-2025 Dynare Team
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
% =========================================================================

% Define testscenario
@#define orderApp = 1
@#define estimParams = 0

% Note that we will set the numerical optimization tolerance levels very large to speed up the testsuite
@#define optimizer = 5

@#include "RBC_MoM_common.inc"

shocks;
var u_a; stderr 0.0072;        
var n; stderr 0.01;
end; 

varobs n c iv;

@#if estimParams == 0
estimated_params;
    DELTA,         0.025;
    BETTA,         0.984;
    B,             0.5;
    %ETAl,          1;
    ETAc,          1;
    ALFA,          0.667;
    RHOA,          0.979;
    stderr u_a,    0.0072;
    %THETA,         3.48;
    stderr n,      0.01;

end;
@#endif

@#if estimParams == 1
estimated_params;
    DELTA,         0.02,        0,           1;
    BETTA,         0.90,        0,           1;
    B,             0.40,        0,           1;
    %ETAl,          1,           0,           10;
    ETAc,          1.80,        0,           10;
    ALFA,          0.60,        0,           1;
    RHOA,          0.90,        0,           1;
    stderr u_a,    0.01,        0,           1;
    stderr n,      0.01,       0,           1;
end;
@#endif

@#if estimParams == 2
estimated_params;
    DELTA,         0.02,         0,           1,  normal_pdf, 0.02, 0.5;
    BETTA,         0.90,         0,           1,  beta_pdf, 0.90, 0.25;
    B,             0.40,         0,           1,  normal_pdf, 0.40, 0.5;
    %ETAl,          1,            0,           10, normal_pdf, 0.25, 0.0.1;
    ETAc,          1.80,         0,           10, normal_pdf, 1.80, 0.5;
    ALFA,          0.60,         0,           1,  normal_pdf, 0.60, 0.5;
    RHOA,          0.90,         0,           1,  normal_pdf, 0.90, 0.5;
    stderr u_a,    0.01,         0,           1,  normal_pdf, 0.01, 0.5;
    stderr n,      0.001,        0,           1,  normal_pdf, 0.01, 0.5;
end;
@#endif

% Simulate data
stoch_simul(order=@{orderApp},pruning,nodisplay,nomoments,periods=250);
send_endogenous_variables_to_workspace;
save('RBC_MoM_data_@{orderApp}.mat', options_.varobs{:} );
pause(1);



%--------------------------------------------------------------------------
% Method of Moments Estimation
%--------------------------------------------------------------------------
% matched_moments blocks : We don't have an interface yet
matched_moments;
c*c;
c*iv;
c*n;
iv*iv;
iv*n;
n*n;

c*c(-1);
n*n(-1);
iv*iv(-1);
end;


method_of_moments(
% Necessary options
      mom_method = SMM           % method of moments method; possible values: GMM|SMM
    , datafile   = 'RBC_MoM_data_@{orderApp}.mat' % name of filename with data
% Options for both GMM and SMM
    , order = @{orderApp}                 % order of Taylor approximation in perturbation
    , pruning                             % use pruned state space system at higher-order
    , weighting_matrix = ['identity_matrix'] % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename. Size of cell determines stages in iterated estimation, e.g. two state with ['DIAGONAL','OPTIMAL']
    , weighting_matrix_scaling_factor=100000  % scaling of weighting matrix in objective function
    , burnin=250                          % number of periods dropped at beginning of simulation

% Optimization options that can be set by the user in the mod file, otherwise default values are provided
    , mode_compute = @{optimizer}         % specifies the optimizer for minimization of moments distance
    , silent_optimizer                  % run minimization of moments distance silently without displaying results or saving files in between
    , hp_filter= 1600
);

options_.hp_filter=0;


method_of_moments(
% Necessary options
      mom_method = SMM           % method of moments method; possible values: GMM|SMM
    , datafile   = 'RBC_MoM_data_@{orderApp}.mat' % name of filename with data
% Options for both GMM and SMM
    , order = @{orderApp}                 % order of Taylor approximation in perturbation
    , pruning                             % use pruned state space system at higher-order
    , weighting_matrix = ['identity_matrix'] % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename. Size of cell determines stages in iterated estimation, e.g. two state with ['DIAGONAL','OPTIMAL']
    , weighting_matrix_scaling_factor=100000  % scaling of weighting matrix in objective function
    , burnin=250                          % number of periods dropped at beginning of simulation

% Optimization options that can be set by the user in the mod file, otherwise default values are provided
    , mode_compute = @{optimizer}         % specifies the optimizer for minimization of moments distance
    , silent_optimizer                  % run minimization of moments distance silently without displaying results or saving files in between
    , one_sided_hp_filter= 1600
);

options_.one_sided_hp_filter=0;

method_of_moments(
% Necessary options
      mom_method = SMM           % method of moments method; possible values: GMM|SMM
    , datafile   = 'RBC_MoM_data_@{orderApp}.mat' % name of filename with data
% Options for both GMM and SMM
    , order = @{orderApp}                 % order of Taylor approximation in perturbation
    , pruning                             % use pruned state space system at higher-order
    , weighting_matrix = ['identity_matrix'] % weighting matrix in moments distance objective function; possible values: OPTIMAL|IDENTITY_MATRIX|DIAGONAL|filename. Size of cell determines stages in iterated estimation, e.g. two state with ['DIAGONAL','OPTIMAL']
    , weighting_matrix_scaling_factor=100000 % scaling of weighting matrix in objective function
    , burnin=250                          % number of periods dropped at beginning of simulation

% Optimization options that can be set by the user in the mod file, otherwise default values are provided
    , mode_compute = @{optimizer}         % specifies the optimizer for minimization of moments distance
    , silent_optimizer                  % run minimization of moments distance silently without displaying results or saving files in between
    , bandpass_filter= [8,20]
);

