function M_ = set_shock_stderr_value_locally(M_,exoname,value)
% M_ = set_shock_stderr_value_locally(M_,exoname,value)
% INPUTS
% - M_                  [structure]     Matlab's structure describing the model
% - exoname             [string]        shock name
% - value               [value]         shock standard error
%
% OUTPUTS
% - M_                  [structure]     Updated Matlab structure describing the model
%
% Remarks: this will only update the diagonal and not take care of
% calibrated correlations/covariances

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

i = strmatch(exoname,M_.exo_names,'exact');

if isempty(i)
    error(['Shock name ' exoname ' doesn''t exist'])
end

M_.Sigma_e(i,i) = value^2;
