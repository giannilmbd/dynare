function i = name2index(M_, estim_params_, type, name1, name2, name3 )
% i = name2index(M_, estim_params_, type, name1, name2, name3 )
% -------------------------------------------------------------------------
% Returns the index associated to an estimated object (deep parameter,
% variance of a structural shock or measurement error, covariance between
% two structural shocks, covariance between two measurement errors, or
% skewness of structural shocks).
%
% INPUTS:
%   M_              [structure]    Dynare structure (related to model definition).
%   estim_params_   [structure]    Dynare structure (related to estimation).
%   type            [string]       'DeepParameter', 'MeasurementError' (for measurement equation error) or 'StructuralShock' (for structural shock).
%   name1           [string]
%   name2           [string]
%   name3           [string]
% OUTPUTS
%   i               [integer]      column index (in x2, an array of mh draws) associated to name1[,name2].
%
% SPECIAL REQUIREMENTS
%   none

% Copyright © 2008-2025 Dynare Team
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

i = []; % initialize output

if strcmpi(type,'DeepParameter')
    i = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn + estim_params_.nsx + ...
        strmatch(name1, M_.param_names(estim_params_.param_vals(:,1)), 'exact');
    if nargin>4
        disp('The last input argument(s) are useless!');
    end
    if isempty(i)
        disp([name1 ' is not an estimated deep parameter!']);
    end
    return
end

if strcmpi(type,'StructuralShock')
    if nargin<5% Covariance matrix diagonal term.
        i = strmatch(name1, M_.exo_names(estim_params_.var_exo(:,1)), 'exact');
        if isempty(i)
            disp(['The standard deviation of ' name1  ' is not an estimated parameter!']);
            return
        end
    elseif nargin<6% Covariance matrix off-diagonal term
        offset = estim_params_.nvx+estim_params_.nvn;
        try
            list_of_structural_shocks = [ M_.exo_names(estim_params_.corrx(:,1)) , M_.exo_names(estim_params_.corrx(:,2)) ];
            k1 = strmatch(name1, list_of_structural_shocks(:,1), 'exact');
            k2 = strmatch(name2, list_of_structural_shocks(:,2), 'exact');
            i = offset+intersect(k1,k2);
            if isempty(i)
                k1 = strmatch(name1, list_of_structural_shocks(:,2), 'exact');
                k2 = strmatch(name2, list_of_structural_shocks(:,1), 'exact');
                i = offset+intersect(k1, k2);
            end
            if isempty(i)
                if isempty(i)
                    disp(['The correlation between ' name1 ' and ' name2 ' is not an estimated parameter!']);
                    return
                end
            end
        catch
            disp('Off diagonal terms of the covariance matrix are not estimated (state equation)');
        end
    elseif nargin<7% coskewness matrix
        if ~strcmp(name1, name2) || ~strcmp(name1, name3)
            error('name2index: We currently do not support estimating co-skewness parameters where the variable names differ. Only diagonal/skewness parameters are supported.');
        end
        try
            i = strmatch(name1, M_.exo_names(estim_params_.skew_exo(:,1)), 'exact');
            if isempty(i)
                disp(['The co-skewness between ' name1 ', ' name2 ', and ' name3 ' is not an estimated parameter!']);
                return
            else
                i = estim_params_.nvx + estim_params_.nvn + estim_params_.ncx + estim_params_.ncn + i;
            end
        catch
            disp('Off diagonal terms of the co-skewness matrix are not estimated (state equation)');
        end
    end
end

if strcmpi(type,'MeasurementError')
    if nargin<5% Covariance matrix diagonal term
        i = estim_params_.nvx + strmatch(name1, M_.endo_names(estim_params_.var_endo(:,1)), 'exact');
        if isempty(i)
            disp(['The standard deviation of the measurement error on ' name1  ' is not an estimated parameter!']);
            return
        end
    else% Covariance matrix off-diagonal term
        offset = estim_params_.nvx+estim_params_.nvn+estim_params_.ncx;
        try
            list_of_measurement_errors = [M_.endo_names(estim_params_.corrn(:,1),1) , M_.endo_names(estim_params_.corrn(:,2),1)];
            k1 = strmatch(name1,list_of_measurement_errors(:,1),'exact');
            k2 = strmatch(name2,list_of_measurement_errors(:,2),'exact');
            i = offset+intersect(k1,k2);
            if isempty(i)
                k1 = strmatch(name1,list_of_measurement_errors(:,2),'exact');
                k2 = strmatch(name2,list_of_measurement_errors(:,1),'exact');
                i = offset+intersect(k1,k2);
            end
            if isempty(i)
                if isempty(i)
                    disp(['The correlation between the measurement errors on ' name1 ' and ' name2 ' is not an estimated parameter!']);
                    return
                end
            end
        catch
            disp('Off diagonal terms of the covariance matrix are not estimated (measurement equation)');
        end
    end
end
