function w = row_header_width(M_,estim_params_,bayestopt_)
% w = row_header_width(M_,estim_params_,bayestopt_)
% -------------------------------------------------------------------------
% This function computes the width of the row headers for the estimation results
%
% INPUTS
%   estim_params_    [structure]
%   M_               [structure]
%   bayestopt_       [structure]
%
% OUTPUTS
%   w                integer
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright © 2006-2025 Dynare Team
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

w = 0;
if estim_params_.np % estimated structural parameters
    w = cellofchararraymaxlength(bayestopt_.name);
end
if estim_params_.nvx % estimated stderr parameters for structural shocks
    w = max(w, cellofchararraymaxlength(M_.exo_names(estim_params_.var_exo(1:estim_params_.nvx,1))));
end
if estim_params_.nvn % estimated stderr parameters for measurement errors
    w = max(w, cellofchararraymaxlength(M_.endo_names(estim_params_.var_endo(1:estim_params_.nvn,1))));
end
if estim_params_.ncx % estimated corr parameters for structural shocks
    for i=1:estim_params_.ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        w = max(w, length(M_.exo_names{k1})+length(M_.exo_names{k2}));
    end
end
if estim_params_.ncn % estimated corr parameters for measurement errors
    for i=1:estim_params_.ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        w = max(w, length(M_.endo_names{k1})+length(M_.endo_names{k2}));
    end
end
if estim_params_.nsx % estimated skew parameters for structural shocks
    for i=1:estim_params_.nsx
        k = estim_params_.skew_exo(i,1);
        w = max(w, length(M_.exo_names{k}));
    end
end