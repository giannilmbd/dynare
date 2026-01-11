function simulation = simul_static_model(varargin)
% simulation = simul_static_model(varargin)
% Wrapper simulating a stochastic static model (with arbitrary precision).
%
% INPUTS
% - varargin            inputs for backward_model_irf.irf
%
% OUTPUTS
% - simulation          [dseries]     Simulated endogenous and exogenous variables.
%

% Copyright © 2019-2025 Dynare Team
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

global M_ options_ oo_;

[simulation, oo_]= backward_model.simul_static_model(M_,options_,oo_,varargin{:});
