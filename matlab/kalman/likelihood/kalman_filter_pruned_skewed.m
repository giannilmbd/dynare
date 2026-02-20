function [LIK, LIKK] = kalman_filter_pruned_skewed( ...
    Y, start, last, presample, ...
    mu_tm1_tm1, Sigma_tm1_tm1, Gamma_tm1_tm1, nu_tm1_tm1, Delta_tm1_tm1, ...
    T, R, Z, mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta, Sigma_eps, ...
    kalman_tol, rescale_prediction_error_covariance, prune_tol, mvnlogcdf, rank_deficiency_transform, verbose)
% [LIK, LIKK] = kalman_filter_pruned_skewed(Y,start,last,presample,mu_tm1_tm1,Sigma_tm1_tm1,Gamma_tm1_tm1,nu_tm1_tm1,Delta_tm1_tm1,T,R,Z,mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta,Sigma_eps,kalman_tol,rescale_prediction_error_covariance,prune_tol,mvnlogcdf,rank_deficiency_transform,verbose)
% -------------------------------------------------------------------------
% Evaluate negative log-likelihood value of linear state space model
% with skew normally distributed innovations and normally distributed noise:
%   α(t) = T*α(t-1) + R*η(t)   [state transition equation]
%   y(t) = Z*x(t)   + ε(t)     [observation equation]
%   η(t) ~ CSN(mu_eta,Sigma_eta,Gamma_eta,nu_eta=0,Delta_eta=I) [innovations, shocks]
%   ε(t) ~ N(mu_eps=0,Sigma_eps)                                [noise, measurement error]
% Dimensions:
%   α(t) is (endo_nbr by 1) state vector
%   y(t) is (varobs_nbr by 1) control vector, i.e. observable variables
%   η(t) is (exo_nbr by 1) vector of shocks
%   ε(t) is (varobs_nbr by 1) vector of measurement errors
% Notes:
% - mu_eta is set such that E[η(t)] = 0
% - This filter works for any nu_eta and Delta_eta, but we currently restrict
%   the shocks to be independent skew normally distributed, which is a special
%   case of the closed skew normal distribution with nu_eta=0 and Delta_eta=I
% -------------------------------------------------------------------------
% INPUTS
% - Y                                     [varobs_nbr by nobs]           matrix with data
% - start                                 [integer scalar]               first observation period to use
% - last                                  [integer scalar]               last observation period to use, (last-first) has to be inferior to obs_nbr
% - presample                             [integer scalar]               number of initial iterations to be discarded when evaluating the likelihood
% - mu_tm1_tm1                            [endo_nbr by 1]                mu_0_0: initial filtered value of location parameter of CSN distributed state vector (does not equal expectation vector unless Gamma_0_0=0)
% - Sigma_tm1_tm1                         [endo_nbr by endo_nbr]         Sigma_0_0: initial filtered value of scale parameter of CSN distributed state vector (does not equal covariance matrix unless Gamma_0_0=0)
% - Gamma_tm1_tm1                         [skew_dim by endo_nbr]         Gamma_0_0: initial filtered value of skewness shape parameter of CSN distributed states vector (if 0 then CSN reduces to Gaussian)
% - nu_tm1_tm1                            [skew_dim by endo_nbr]         nu_0_0: initial filtered value of skewness conditioning parameter of CSN distributed states vector (enables closure of CSN distribution under conditioning, irrelevant if Gamma_0_0=0)
% - Delta_tm1_tm1                         [skew_dim by skew_dim]         Delta_0_0: initial filtered value of skewness marginalization parameter of CSN distributed states vector (enables closure of CSN distribution under marginalization, irrelevant if Gamma_0_0=0)
% - T                                     [endo_nbr by endo_nbr]         state transition matrix mapping previous states to current states
% - R                                     [endo_nbr by exo_nbr]          state transition matrix mapping current innovations to current states
% - Z                                     [varobs_nbr by endo_nbr]       observation equation matrix mapping current states into current observables
% - mu_eta                                [exo_nbr by 1]                 value of location parameter of CSN distributed shocks (does not equal expectation vector unless Gamma_eta=0)
% - Sigma_eta                             [exo_nbr by exo_nbr]           value of scale parameter of CSN distributed shocks (does not equal covariance matrix unless Gamma_eta=0)
% - Gamma_eta                             [skeweta_dim by eta_nbr]       value of skewness shape parameter of CSN distributed shocks (if 0 then CSN reduces to Gaussian)
% - nu_eta                                [skeweta_dim by eta_nbr]       value of skewness conditioning parameter of CSN distributed shocks (enables closure of CSN distribution under conditioning)
% - Delta_eta                             [skeweta_dim by skeweta_dim]   value of skewness marginalization parameter of CSN distributed shocks (enables closure of CSN distribution under marginalization)
% - Sigma_eps                             [varobs_nbr by varobs_nbr]     scale parameter of normally distributed measurement errors (equals covariance matrix), if no measurement errors this is a zero scalar
% - kalman_tol                            [double scalar]                tolerance parameter for invertibility of covariance matrix
% - rescale_prediction_error_covariance   [boolean]                      1: rescales the prediction error covariance (Omega) to avoid badly scaled matrix
% - prune_tol                             [double]                       threshold to prune redundant skewness dimensions, if set to 0 no pruning will be done
% - mvnlogcdf                             [string]                       name of function to compute log Gaussian cdf, possible values: 'gaussian_log_mvncdf_mendell_elston', 'mvncdf'
% - rank_deficiency_transform             [boolean]                      indicator if prediction step is done on joint distribution, [x_t', eta_t']', useful for models where T is singular
% - verbose                               [boolean]                      additional output for debuging
% -------------------------------------------------------------------------
% OUTPUTS
% - LIK                                   [double scalar]                value of negative log likelihood
% - LIKK                                  [(last-start+1) by 1]          vector of densities for each observation

% Copyright © 2024-2025 Gaygysyz Guljanov, Willi Mutschler, Mark Trede
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

%%%%%%%%%%%%%%%%%%%
% INITIALIZATIONS %
%%%%%%%%%%%%%%%%%%%
smpl = last-start+1;         % sample size
t = start;                   % time index
likk = zeros(smpl,1);        % vector gathering the densities
LIK = Inf;                   % default value of the log likelihood
Omega_singular = true;       % indicator whether Omega is singular
LIKK = [];                   % output vector for individal density contributions
rescale_prediction_error_covariance0 = rescale_prediction_error_covariance; % store option
const2pi = -0.5*size(Y,1)*log(2*pi); % constant in Gaussian probability density function (pdf)

if verbose && (rcond(T) < kalman_tol)
    warning('kalman_filter_pruned_skewed: state transition matrix is singular');
end
if rank_deficiency_transform
    T_bar = [T, R]; % multiplying matrix of joint distribution
else
    mu_eta = R*mu_eta;
    Sigma_eta = R*Sigma_eta*R';
    Gamma_eta = Gamma_eta/(R'*R)*R';
    Gamma_eta_X_Sigma_eta = Gamma_eta*Sigma_eta;
    Delta22_common = Delta_eta + Gamma_eta_X_Sigma_eta*Gamma_eta';
end

while t <= last
    s = t-start+1;

    %%%%%%%%%%%%%%%%%%%%
    % STATE PREDICTION %
    %%%%%%%%%%%%%%%%%%%%
    if rank_deficiency_transform
        % parameters of joint distribution, [x_t', eta_t']'
        mu_bar_tm1_tm1 = [mu_tm1_tm1; mu_eta];
        nu_bar_tm1_tm1 = [nu_tm1_tm1; nu_eta];
        Sigma_bar_tm1_tm1 = blkdiag_two(Sigma_tm1_tm1, Sigma_eta);
        Gamma_bar_tm1_tm1 = blkdiag_two(Gamma_tm1_tm1, Gamma_eta);
        Delta_bar_tm1_tm1 = blkdiag_two(Delta_tm1_tm1, Delta_eta);
        % linear transformation of the joint distribution
        [mu_t_tm1, Sigma_t_tm1, Gamma_t_tm1, nu_t_tm1, Delta_t_tm1] = csn_statespace_linear_transform(T_bar, mu_bar_tm1_tm1, Sigma_bar_tm1_tm1, Gamma_bar_tm1_tm1, nu_bar_tm1_tm1, Delta_bar_tm1_tm1);
    else
        % auxiliary matrices
        Gamma_tm1_tm1_X_Sigma_tm1_tm1 = Gamma_tm1_tm1*Sigma_tm1_tm1;
        Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT = Gamma_tm1_tm1_X_Sigma_tm1_tm1*T';
        mu_t_tm1  = T*mu_tm1_tm1 + mu_eta;
        Sigma_t_tm1 = T*Sigma_tm1_tm1*T' + Sigma_eta;
        Sigma_t_tm1 = 0.5*(Sigma_t_tm1 + Sigma_t_tm1'); % ensure symmetry
        invSigma_t_tm1 = pinv(Sigma_t_tm1); % pseudo-inverse is valid for expressions that are derived from the conditional CSN distributions
        Gamma_t_tm1 = [Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT; Gamma_eta_X_Sigma_eta]*invSigma_t_tm1;
        nu_t_tm1 = [nu_tm1_tm1; nu_eta];
        Delta11_t_tm1 = Delta_tm1_tm1 + Gamma_tm1_tm1_X_Sigma_tm1_tm1*Gamma_tm1_tm1' - Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT*invSigma_t_tm1*Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT';
        Delta22_t_tm1 = Delta22_common - Gamma_eta_X_Sigma_eta*invSigma_t_tm1*Gamma_eta_X_Sigma_eta';
        Delta12_t_tm1 = -Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT*invSigma_t_tm1*Gamma_eta_X_Sigma_eta';
        Delta_t_tm1 = [Delta11_t_tm1 , Delta12_t_tm1; Delta12_t_tm1' , Delta22_t_tm1];
        Delta_t_tm1 = 0.5*(Delta_t_tm1 + Delta_t_tm1'); % ensure symmetry
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOG-LIKELIHOOD CONTRIBUTIONS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The conditional distribution of y(t) given y(t-1) is
    % (y(t)|y(t-1)) ~ CSN(mu_y,Sigma_y,Gamma_y,nu_y,Delta_y)
    %               = mvncdf(Gamma_y*(y(t)-mu_y),nu_y,Delta_y) / mvncdf(0,nu_y,Delta_y+Gamma_y*Sigma_y*Gamma_y') * mvnpdf(y(t),mu_y,Sigma_y)
    %               = gaussian_cdf_top / gaussian_cdf_bottom * gaussian_pdf;
    % where:
    %   mu_y    = Z*mu_t_tm1 + mu_eps =: y_predicted
    %   Sigma_y = Z*Sigma_t_tm1*Z' + Sigma_eps =: Omega
    %   K_Gauss := Sigma_t_tm1*Z'*inv(Sigma_y)
    %   Gamma_y = Gamma_t_tm1*Sigma_t_tm1*Z'*inv(Z*Sigma_t_tm1*Z' + Sigma_eps) = Gamma_t_tm1*K_Gauss =: K_Skewed
    %   nu_y    = nu_t_tm1
    %   Delta_y = Delta_t_tm1 + Gamma_t_tm1*Sigma_t_tm1*Gamma_t_tm1'...
    %             - Gamma_t_tm1*Sigma_t_tm1*Z'*inv(Z*Sigma_t_tm1*Z')*Z*Sigma_t_tm1*Gamma_t_tm1'...
    %             + (Gamma_t_tm1*Sigma_t_tm1*Z'*inv(Z*Sigma_t_tm1*Z') - Gamma_t_tm1*Sigma_t_tm1*Z'*inv(Z*Sigma_t_tm1*Z' + Sigma_eps))*Z*Sigma_t_tm1*Gamma_t_tm1';
    %           = Delta_t_tm1 + (Gamma_t_tm1-K_Skewed*Z)*Sigma_t_tm1*Gamma_t_tm1'

    % log-likelihood step 1/6: prune redundant skewness dimensions
    if prune_tol > 0
        [Sigma_t_tm1, Gamma_t_tm1, nu_t_tm1, Delta_t_tm1] = csn_prune_distribution(Sigma_t_tm1,Gamma_t_tm1,nu_t_tm1,Delta_t_tm1,prune_tol);
    end

    % log-likelihood step 2/6: compute Gaussian prediction error and its covariance Omega as well as log(det(Omega)) and inv(Omega) in a numerically stable way
    prediction_error = Y(:,t) - Z*mu_t_tm1;
    Omega = Z*Sigma_t_tm1*Z' + Sigma_eps;
    badly_conditioned_Omega = false;
    if rescale_prediction_error_covariance
        sigOmega = sqrt(diag(Omega)); % standard deviations
        if any(diag(Omega)<kalman_tol) || rcond(Omega./(sigOmega*sigOmega'))<kalman_tol % Omega./(sigOmega*sigOmega') effectively converts Omega to a correlation matrix
            badly_conditioned_Omega = true;
        end
    else
        if rcond(Omega) < kalman_tol
            sigOmega = sqrt(diag(Omega)); % standard deviations
            if any(diag(Omega)<kalman_tol) || rcond(Omega./(sigOmega*sigOmega'))<kalman_tol % Omega./(sigOmega*sigOmega') effectively converts Omega to a correlation matrix
                badly_conditioned_Omega = true;
            else
                rescale_prediction_error_covariance = 1;
            end
        end
    end
    if badly_conditioned_Omega
        % if ~all(abs(Omega(:))<kalman_tol), then use univariate filter (will remove observations with zero variance prediction error), otherwise this is a pathological case and the draw is discarded
        if verbose
            warning('kalman_pruned_skewed: Omega is badly conditioned, discard draw as univariate filter will not be used (this overwrites Dynare''s default behavior)');
        end
        return
    else
        Omega_singular = false;
        if rescale_prediction_error_covariance % Omega needs to be rescaled to avoid numerical instability
            % compute the log-determinant of covariance matrix Omega in a numerically stable way:
            % let: corrOmega = Omega./(sigOmega*sigOmega') be the correlation matrix corresponding to Omega
            % then: Omega = (sigOmega*sigOmega') .* corrOmega = diag(sigOmega) * corrOmega * diag(sigOmega)
            % determinant: det(Omega) = det(diag(sigOmega)) * det(corrOmega) * det(diag(sigOmega)) = det(corrOmega) * prod(sigOmega)^2
            % taking logs: log(det(Omega)) = log(det(corrOmega)) + 2*sum(log(sigOmega))
            log_detOmega = log(det(Omega./(sigOmega*sigOmega'))) + 2*sum(log(sigOmega));
            % compute inv(Omega) in a numerically stable way
            % let: corrOmega = Omega./(sigOmega*sigOmega') be the correlation matrix corresponding to Omega
            % then: Omega = (sigOmega*sigOmega') .* corrOmega = diag(sigOmega) * corrOmega * diag(sigOmega)
            %       inv(Omega) = diag(1./sigOmega) * inv(corrOmega) * diag(1./sigOmega) = inv(corrOmega) ./ (sigOmega*sigOmega')
            invOmega = inv(Omega./(sigOmega*sigOmega'))./(sigOmega*sigOmega');
            rescale_prediction_error_covariance = rescale_prediction_error_covariance0; % reset option as it might have been updated above
        else
            log_detOmega = log(det(Omega)); % compute the log-determinant of Omega directly
            invOmega = inv(Omega); % compute the inverse of Omega directly
        end

        % log-likelihood step 3/6: compute Kalman gains
        K_Gauss = Sigma_t_tm1*Z'*invOmega;
        K_Skewed = Gamma_t_tm1*K_Gauss;

        % log-likelihood step 4/6: evaluate Gaussian cdfs (specific to skewed Kalman filter)
        % bottom one: mvncdf(0,nu_y,Delta_y + Gamma_y*Sigma_y*Gamma_y')
        % top one: mvncdf(Gamma_y*(y(t)-mu_y),nu_y,Delta_y)

        tmp = Gamma_t_tm1*Sigma_t_tm1;
        Delta_y = Delta_t_tm1 + tmp*Gamma_t_tm1' - K_Skewed*Z*tmp';
        Delta_y = 0.5 * (Delta_y + Delta_y'); % ensure symmetry

        cdf_bottom_cov = Delta_y + K_Skewed*Omega*K_Skewed';
        cdf_bottom_cov = 0.5*(cdf_bottom_cov + cdf_bottom_cov'); % ensure symmetry
        if strcmp(mvnlogcdf,'gaussian_log_mvncdf_mendell_elston')
            % requires zero mean and correlation matrix as inputs
            normalization_Delta_y = diag(1./sqrt(diag(cdf_bottom_cov)));
            cdf_bottom_cov = normalization_Delta_y*cdf_bottom_cov*normalization_Delta_y; % this is now a correlation matrix
            cdf_bottom_cov = 0.5*(cdf_bottom_cov + cdf_bottom_cov'); % ensure symmetry
            if ~isempty(cdf_bottom_cov)
                try
                    log_gaussian_cdf_bottom = gaussian_log_mvncdf_mendell_elston(-normalization_Delta_y*nu_t_tm1, cdf_bottom_cov);
                catch
                    message = get_error_message(57);
                    if verbose
                        warning('kalman_filter_pruned_skewed: %s.',message);
                    end
                    return
                end
            else
                log_gaussian_cdf_bottom = 0;
            end

            normalization_Delta_y = diag(1./sqrt(diag(Delta_y)));
            Delta_y = normalization_Delta_y*Delta_y*normalization_Delta_y; % this is now a correlation matrix
            Delta_y = 0.5*(Delta_y + Delta_y'); % ensure symmetry
            if ~isempty(Delta_y)
                try
                    log_gaussian_cdf_top = gaussian_log_mvncdf_mendell_elston(normalization_Delta_y*(K_Skewed*prediction_error - nu_t_tm1), Delta_y);
                catch
                    message = get_error_message(58);
                    if verbose
                        warning('kalman_filter_pruned_skewed: %s.',message);
                    end
                    return
                end
            else
                log_gaussian_cdf_top = 0;
            end
        elseif strcmp(mvnlogcdf,'mvncdf')
            try
                log_gaussian_cdf_bottom = log(mvncdf(zeros(size(nu_t_tm1,1),1), nu_t_tm1, cdf_bottom_cov));
            catch
                message = get_error_message(57);
                if verbose
                    warning('kalman_filter_pruned_skewed: %s.',message);
                end
                return
            end
            try
                log_gaussian_cdf_top = log(mvncdf(K_Skewed*prediction_error, nu_t_tm1, Delta_y));
            catch
                message = get_error_message(58);
                if verbose
                    warning('kalman_filter_pruned_skewed: %s.',message);
                end
                return
            end
        end

        % log-likelihood step 5/6: evaluate Gaussian pdf (common with Gaussian Kalman filter)
        % log_gaussian_pdf = log(mvnpdf(Y(:,t), y_predicted, Omega))
        try
            log_gaussian_pdf = const2pi - 0.5*log_detOmega - 0.5*transpose(prediction_error)*invOmega*prediction_error;
        catch
            message = get_error_message(56);
            if verbose
                warning('kalman_filter_pruned_skewed: %s.',message);
            end
            return
        end

        % log-likelihood step 6/6: collect likelihood contribution
        likk(s) = log_gaussian_cdf_top - log_gaussian_cdf_bottom + log_gaussian_pdf;

        %%%%%%%%%%%%%%%%%%%
        % STATE FILTERING %
        %%%%%%%%%%%%%%%%%%%
        % already assign for next time step
        mu_tm1_tm1 = mu_t_tm1 + K_Gauss*prediction_error;
        Sigma_tm1_tm1 = Sigma_t_tm1 - K_Gauss*Z*Sigma_t_tm1;
        Gamma_tm1_tm1 = Gamma_t_tm1;
        nu_tm1_tm1 = nu_t_tm1 - K_Skewed*prediction_error;
        Delta_tm1_tm1 = Delta_t_tm1;
        Sigma_tm1_tm1 = 0.5*(Sigma_tm1_tm1 + Sigma_tm1_tm1'); % ensure symmetry
        Delta_tm1_tm1 = 0.5*(Delta_tm1_tm1 + Delta_tm1_tm1'); % ensure symmetry

    end
    t = t+1;
    if verbose
        fprintf('Skewness Dimension: %d\n',size(Gamma_t_t,1));
    end
end

if Omega_singular
    message = get_error_message(60);
    warning('kalman_filter_pruned_skewed: %s.',message);
    return
end

% compute minus the log-likelihood
LIK = -1*sum(likk(1+presample:end));
LIKK = -1*likk;


%% auxiliary functions
function res_mat = blkdiag_two(mat1, mat2)
    % Makes a block diagonal matrix out of two matrices
    [nrow_mat1, ncol_mat1] = size(mat1); [nrow_mat2, ncol_mat2] = size(mat2);
    upper_mat = zeros(nrow_mat1, ncol_mat2);
    lower_mat = zeros(nrow_mat2, ncol_mat1);
    res_mat = [mat1, upper_mat; lower_mat, mat2];
end % blkdiag_two


end % kalman_filter_pruned_skewed