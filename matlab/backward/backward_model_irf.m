function [deviations, baseline, irfs] = backward_model_irf(varargin)
% [deviations, baseline, irfs] = backward_model_irf(varargin)
% Wrapper returning impulse response functions.
%
% INPUTS
% - varargin            inputs for backward_model_irf.irf
%
% OUTPUTS
% - irfs                [struct of dseries]
%
% REMARKS
% - The names of the fields in the returned structure are given by the name
%   of the innovations listed in the second input argument. Each field gather
%   the associated paths for endogenous variables listed in the third input
%   argument.
% - If second argument is not empty, periods must not be greater than innovationbaseline.nobs.

% Copyright © 2017-2025 Dynare Team
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

global M_ options_ oo_

switch nargout
    case 1
        [deviations] = backward_model.irf(M_,options_,oo_,varargin{:});
    case 2
        [deviations, baseline] = backward_model.irf(M_,options_,oo_,varargin{:});
    case 3
        [deviations, baseline, irfs] = backward_model.irf(M_,options_,oo_,varargin{:});
end