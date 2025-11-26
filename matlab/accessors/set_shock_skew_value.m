function set_shock_skew_value(varargin)
% set_shock_skew_value(varargin)
% -------------------------------------------------------------------------
% Update entries of the coskewness tensor M_.Skew_e for shocks identified
% by their names. Supports either the univariate-diagonal case or the full
% trivariate entry with automatic filling of all index permutations.
%
% USAGE
% - set_shock_skew_value(exoname, value)
%     Sets M_.Skew_e(i,i,i) = value for the shock named exoname.
%
% - set_shock_skew_value(exoname1, exoname2, exoname3, value)
%     Sets M_.Skew_e(i1,i2,i3) = value and all permutations of (i1,i2,i3), where
%     i1, i2, i3 correspond to exoname1, exoname2, exoname3.
%
% INPUTS
% - exoname, exoname1, exoname2, exoname3  [char]  shock names as in M_.exo_names
% - value                                  [double] scalar skewness/coskewness value
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

global M_

if ~(nargin == 2 || nargin == 4)
    error('set_shock_skew_value must be called with (exoname, value) or (exoname1, exoname2, exoname3, value)');
end

if nargin == 2
    exoname = varargin{1};
    value = varargin{2};

    i = strmatch(exoname, M_.exo_names, 'exact');
    if isempty(i)
        error(['Shock name ' exoname ' doesn''t exist'])
    end

    M_.Skew_e(i,i,i) = value;
else
    exoname1 = varargin{1};
    exoname2 = varargin{2};
    exoname3 = varargin{3};
    value = varargin{4};

    i1 = strmatch(exoname1, M_.exo_names, 'exact');
    i2 = strmatch(exoname2, M_.exo_names, 'exact');
    i3 = strmatch(exoname3, M_.exo_names, 'exact');

    if isempty(i1)
        error(['Shock name ' exoname1 ' doesn''t exist'])
    end
    if isempty(i2)
        error(['Shock name ' exoname2 ' doesn''t exist'])
    end
    if isempty(i3)
        error(['Shock name ' exoname3 ' doesn''t exist'])
    end

    % Assign value to all permutations of (i1, i2, i3)
    M_.Skew_e(i1, i2, i3) = value;
    M_.Skew_e(i1, i3, i2) = value;
    M_.Skew_e(i2, i1, i3) = value;
    M_.Skew_e(i2, i3, i1) = value;
    M_.Skew_e(i3, i1, i2) = value;
    M_.Skew_e(i3, i2, i1) = value;
end
