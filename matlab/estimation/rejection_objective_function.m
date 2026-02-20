function [fval, exit_flag, arg1, arg2] = rejection_objective_function(fcn, x, threshold_val, varargin)
% [fval, exit_flag, arg1, arg2] = rejection_objective_function(fcn, x, threshold_val, varargin)
% Objective function wrapper that provides posterior kernel for
% MCMC draws and provides a base value in options_ (third element of
% varargin ) after computing lnprior (needs bayestopt_ = 6th element of
% varargin)
%
% INPUTS
% - fcn           [fhandle]   objective function.
% - x             [double]    n*1 vector of parameter values.
% - threshold_val [double]    scalar, threshold value for rejection.
% - varargin      [cell]      additional parameters for fcn.
%
% OUTPUTS
% - fval          [double]    scalar, value of the objective function at x.
% - exit_flag     [integer]   scalar, flag returned by fcn (third output).
% - arg1, arg2                fourth and fifth output arguments of the objective function.

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

    % dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state,derivatives_info)
    % 1        2            3        4  5             6          7          8   9                  10                11                   12


lnprior = priordens(x(:),varargin{6}.pshape,varargin{6}.p6,varargin{6}.p7,varargin{6}.p3,varargin{6}.p4);
varargin{3}.likelihood_base_value    = (threshold_val-lnprior);

[fval, info, exit_flag, arg1, arg2] = fcn(x,varargin{:});

