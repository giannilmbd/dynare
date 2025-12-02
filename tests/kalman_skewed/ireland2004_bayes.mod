% Various tests within Bayesian MCMC estimation
% with skew normally distributed shocks
% -------------------------------------------------------------------------

% Copyright © 2025 Dynare Team
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

@#include "ireland2004_model.inc"

OMEGA  = 0.1377;
RHO_PI = 0.4522;
RHO_G  = 0.3445;
RHO_X  = 0.2006;
RHO_A  = 0.9336;
RHO_E  = 0.9153;

shocks;
var eta_a; stderr 3.2442;
var eta_e; stderr 0.0573;
var eta_z; stderr 0.7546;
var eta_r; stderr 0.2775;
skew eta_a = -0.1925;
skew eta_e = -0.4014;
skew eta_z = -0.5166;
skew eta_r = +0.6066;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1: declaring skewed shocks in different ways in a Bayesian estimation should work %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params;
stderr eta_a,  ,    0,  10, inv_gamma_pdf, (sqrt(30)),   (sqrt(30)) ;
stderr eta_e,  ,    0,  10, inv_gamma_pdf, (sqrt(0.08)), (sqrt(1))  ;
stderr eta_z,  ,    0,  10, inv_gamma_pdf, (sqrt(5)),    (sqrt(15)) ;
stderr eta_r,  ,    0,  10, inv_gamma_pdf, (sqrt(0.50)), (sqrt(2))  ;
OMEGA,         ,    0,   1, beta_pdf,      0.20,         0.10       ;
RHO_PI,        ,    0,   1, gamma_pdf,     0.30,         0.10       ;
RHO_G,         ,    0,   1, gamma_pdf,     0.30,         0.10       ;
RHO_X,         ,    0,   1, gamma_pdf,     0.25,         0.0625     ;
RHO_A,         ,    0,   1, beta_pdf,      0.85,         0.10       ;
RHO_E,         ,    0,   1, beta_pdf,      0.85,         0.10       ;
skew eta_a,                 beta_pdf,         0,         0.40       , -1, 1;
skew eta_e, 0.4, -0.8, 0.9, beta_pdf,         0,         0.40       , -1, 1;
skew eta_z,    ,     ,    , beta_pdf,         0,         0.40       , -1, 1;
skew eta_r,    , -0.9, 0.8, uniform_pdf,       ,                    , -1, 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2: estimated_params_init, estimated_params_remove blocks should work %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimated_params_init(use_calibration);
end;

estimated_params_remove;
stderr eta_a;
stderr eta_z;
stderr eta_r;
OMEGA;
RHO_PI;
RHO_G;
RHO_A;
RHO_E;
skew eta_a;
skew eta_z;
skew eta_r;
end;

estimated_params_init;
RHO_X, 0.200410549459701;
stderr eta_e, 0.057271868850671;
skew eta_e, -0.399902738406484;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 3: Bayesian estimation with RWMH sampler should work without error %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimation(datafile = 'ireland2004_post1980.csv'
          , mode_compute = 4, silent_optimizer
          , lik_init = 1
          , mh_replic = 250
          , mh_nblocks = 1
          , plot_priors = 0
          , kalman_algo = 5  % use pruned skewed Kalman filter
          , skewed_kalman_prune_tol = 0.2 % set extremely high for speed in testsuite
          , mh_jscale = 10
          , tex
          );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4: trace plots and autocorrelation functions should be computed %
%         and plotted for skewness coefficients of structural shocks   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
generate_trace_plots(1);
trace_plot(options_,M_,estim_params_,'StructuralShock',1,'eta_e','eta_e','eta_e');
mh_autocorrelation_function(options_,M_,estim_params_,'StructuralShock',1,'eta_e','eta_e','eta_e');