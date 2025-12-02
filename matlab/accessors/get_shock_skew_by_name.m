function x = get_shock_skew_by_name(varargin)
% x = get_shock_skew_by_name(exoname)
% x = get_shock_skew_by_name(exoname1, exoname2, exoname3)
% returns the skewness/coskewness value for shocks identified by name(s)
%
% INPUTS:
%   exoname:                     shock name (univariate skewness)
%   exoname1, exoname2, exoname3: shock names (coskewness entry)
%
% OUTPUTS
%   x:      skewness/coskewness value
%
% SPECIAL REQUIREMENTS
%   none

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

if ~(nargin == 1 || nargin == 3)
    error('get_shock_skew_by_name must be called with (exoname) or (exoname1, exoname2, exoname3)');
end

if nargin == 1
    exoname = varargin{1};
    i = find(strcmp(exoname, M_.exo_names));
    if isempty(i)
        error(['Can''t find shock ', exoname])
    end
    x = M_.Skew_e(i,i,i);
else
    exoname1 = varargin{1};
    exoname2 = varargin{2};
    exoname3 = varargin{3};

    i1 = find(strcmp(exoname1, M_.exo_names));
    i2 = find(strcmp(exoname2, M_.exo_names));
    i3 = find(strcmp(exoname3, M_.exo_names));

    if isempty(i1)
        error(['Can''t find shock ', exoname1])
    end
    if isempty(i2)
        error(['Can''t find shock ', exoname2])
    end
    if isempty(i3)
        error(['Can''t find shock ', exoname3])
    end

    x = M_.Skew_e(i1,i2,i3);
end