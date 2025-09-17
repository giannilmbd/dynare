function decomposition = shock_decomposition_backward(varargin)
% decomposition = shock_decomposition_backward(varargin)
% Computes and possibly plots the shock decomposition of a backward (possibly nonlinear)
% model simulation.
%
% Inputs:
% - varargin            inputs for backward_model_irf.irf
%
% Output:
% - decomposition        a 3D array of size endo_nbr×nshocks×nperiods where nshocks=length(shocklist)
%                        and nperiods is the number of simulation periods (i.e. excluding the initial
%                        conditions)

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

global M_ options_

decomposition = backward_model.shock_decomposition(M_,options_,varargin{:});