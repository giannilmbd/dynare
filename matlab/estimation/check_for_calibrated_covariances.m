function estim_params_=check_for_calibrated_covariances(estim_params_,M_,varobs_id)
% function estim_params_=check_for_calibrated_covariances(estim_params_,M_,varobs_id)
% find calibrated covariances to consider during estimation
% Inputs
%   -estim_params_  [structure] describing parameters to be estimated
%   -M_             [structure] describing the model
%   -varobs_id      [vector]    index of observed variables in declaration order
%                               (options_.varobs_id), used to map corrn entries
%                               (endogenous variable indices) to observable
%                               positions in M_.H
%
% Outputs
%   -estim_params_  [structure] describing parameters to be estimated
%
% Notes: M is local to this function and not updated when calling
% set_all_parameters

% Copyright © 2013-2023 Dynare Team
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

if isfield(estim_params_,'calibrated_covariances')
    estim_params_ = rmfield(estim_params_,'calibrated_covariances'); %remove if already present
end
if isfield(estim_params_,'calibrated_covariances_ME')
    estim_params_ = rmfield(estim_params_,'calibrated_covariances_ME'); %remove if already present
end

[rows_calibrated, columns_calibrated]=ind2sub(size(M_.Sigma_e),find(tril(M_.Sigma_e,-1))); %find linear indices of preset lower triangular covariance entries

if estim_params_.ncx %delete preset entries actually estimated
    for i=1:estim_params_.ncx
        shock_1 = estim_params_.corrx(i,1);
        shock_2 = estim_params_.corrx(i,2);
        estimated_corr_pos=find(rows_calibrated==shock_1 & columns_calibrated==shock_2);
        if ~isempty(estimated_corr_pos)
            rows_calibrated(estimated_corr_pos)=[];
            columns_calibrated(estimated_corr_pos)=[];
        end
        estimated_corr_pos=find(rows_calibrated==shock_2 & columns_calibrated==shock_1);
        if ~isempty(estimated_corr_pos)
            rows_calibrated(estimated_corr_pos)=[];
            columns_calibrated(estimated_corr_pos)=[];
        end
    end
    if any(rows_calibrated)
        estim_params_.calibrated_covariances.position=[sub2ind(size(M_.Sigma_e),rows_calibrated,columns_calibrated);sub2ind(size(M_.Sigma_e),columns_calibrated,rows_calibrated)]; %get linear entries of upper triangular parts
        estim_params_.calibrated_covariances.cov_value=M_.Sigma_e(estim_params_.calibrated_covariances.position);
    end
end

[rows_calibrated, columns_calibrated]=ind2sub(size(M_.H),find(tril(M_.H,-1))); %find linear indices of preset lower triangular covariance entries

if estim_params_.ncn %delete preset entries actually estimated
    for i=1:estim_params_.ncn
        % corrn contains endogenous variable indices (declaration order),
        % but M_.H is indexed by observable position. Map via varobs_id.
        shock_1 = find(varobs_id == estim_params_.corrn(i,1));
        shock_2 = find(varobs_id == estim_params_.corrn(i,2));
        estimated_corr_pos=find(rows_calibrated==shock_1 & columns_calibrated==shock_2);
        if ~isempty(estimated_corr_pos)
            rows_calibrated(estimated_corr_pos)=[];
            columns_calibrated(estimated_corr_pos)=[];
        end
        estimated_corr_pos=find(rows_calibrated==shock_2 & columns_calibrated==shock_1);
        if ~isempty(estimated_corr_pos)
            rows_calibrated(estimated_corr_pos)=[];
            columns_calibrated(estimated_corr_pos)=[];
        end
    end
end
if any(rows_calibrated)
    estim_params_.calibrated_covariances_ME.position=[sub2ind(size(M_.H),rows_calibrated,columns_calibrated);sub2ind(size(M_.H),columns_calibrated,rows_calibrated)]; %get linear entries of upper triangular parts
    estim_params_.calibrated_covariances_ME.cov_value=M_.H(estim_params_.calibrated_covariances_ME.position);
end

return % --*-- Unit tests --*--

%@test:1

M_.Sigma_e=[1 0; 0 1];
M_.H=[1 0; 0 1];
M_.Correlation_matrix= [1 -0.5; -0.5 1];
M_.Correlation_matrix_ME=[1 -0.5; -0.5 1];
estim_params_.ncx=1;
estim_params_.ncn=1;

estim_params_.corrx=[2 1 NaN -1 1 3 0 0.2000 NaN NaN NaN];
estim_params_.corrn=[2 1 NaN -1 1 3 0 0.2000 NaN NaN NaN];

varobs_id=[1 2];
estim_params_=check_for_calibrated_covariances(estim_params_,M_,varobs_id);
if isfield(estim_params_,'calibrated_covariances_ME') || isfield(estim_params_,'calibrated_covariances')
    t(1)=false;
else
    t(1)=true;
end

M_.Sigma_e=[1 -0.1; -0.1 1];
M_.H=[1 -0.1; -0.1 1];
M_.Correlation_matrix= [1 -0.5; -0.5 1];
M_.Correlation_matrix_ME=[1 0; 0 1];
estim_params_.ncx=1;
estim_params_.ncn=0;

estim_params_.corrx=[2 1 NaN -1 1 3 0 0.2000 NaN NaN NaN];
estim_params_.corrn=[];
varobs_id=[1 2];
estim_params_=check_for_calibrated_covariances(estim_params_,M_,varobs_id);
t(2)=isequal(estim_params_.calibrated_covariances_ME.position,[2;3]);
t(3)=isequal(estim_params_.calibrated_covariances_ME.cov_value,[-0.1;-0.1]);

M_.Sigma_e=[1 -0.1; -0.1 1];
M_.H=[1 -0.1; -0.1 1];
M_.Correlation_matrix= [1 -0.5; -0.5 1];
M_.Correlation_matrix_ME=[1 0; 0 1];
estim_params_.ncx=1;
estim_params_.ncn=1;

estim_params_.corrx=[2 1 NaN -1 1 3 0 0.2000 NaN NaN NaN];
estim_params_.corrn=[2 1 NaN -1 1 3 0 0.2000 NaN NaN NaN];
varobs_id=[1 2];
estim_params_=check_for_calibrated_covariances(estim_params_,M_,varobs_id);
if isfield(estim_params_,'calibrated_covariances_ME') || isfield(estim_params_,'calibrated_covariances')
    t(4)=false;
else
    t(4)=true;
end

%% Test 5: Non-trivial varobs mapping
% corrn uses endogenous variable indices [12 11], but H is 2x2 indexed by
% observable position. With varobs_id=[12 11], corrn(12,11) maps to
% observable positions (1,2), which should match the calibrated entry.
M_.Sigma_e=[1 0; 0 1];
M_.H=[1 -0.3; -0.3 1];
estim_params_.ncx=0;
estim_params_.ncn=1;
estim_params_.corrx=[];
estim_params_.corrn=[12 11 NaN -1 1 3 0 0.2000 NaN NaN NaN];
varobs_id=[12 11];
estim_params_=check_for_calibrated_covariances(estim_params_,M_,varobs_id);
% The estimated corrn should remove the calibrated entry from H

t(5)=~isfield(estim_params_,'calibrated_covariances_ME'); % should have been removed since we estimate this corr

%% Test 6: Non-trivial mapping with calibrated entry NOT estimated
M_.Sigma_e=[1 0; 0 1];
M_.H=[1 -0.3; -0.3 1];
estim_params_.ncx=0;
estim_params_.ncn=0;
estim_params_.corrx=[];
estim_params_.corrn=[];
varobs_id=[12 11];
estim_params_=check_for_calibrated_covariances(estim_params_,M_,varobs_id);
t(6)=isequal(estim_params_.calibrated_covariances_ME.position,[2;3]);
t(7)=isequal(estim_params_.calibrated_covariances_ME.cov_value,[-0.3;-0.3]);

T = all(t);
%@eof:1
