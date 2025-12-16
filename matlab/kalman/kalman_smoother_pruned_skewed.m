function [x_t_T, epsilonhat, eta_t_T, x_t_t,atilde,P,aK,PK,decomp,V,aalphahat,eetahat,x_0_T,state_uncertainty0] = kalman_smoother_pruned_skewed( ...
    Y, ...
    mu_tm1_tm1, Sigma_tm1_tm1, Gamma_tm1_tm1, nu_tm1_tm1, Delta_tm1_tm1, ...
    G, R, F, ...
    mu_eta, Sigma_eta, Gamma_eta, nu_eta, Delta_eta, Sigma_eps, ...
    kalman_tol, prune_tol, mvnlogcdf, ...
    console_mode)
% [x_t_T, epsilonhat, eta_t_T, x_t_t,atilde,P,aK,PK,decomp,V,aalphahat,eetahat,x_0_T,state_uncertainty0] = kalman_smoother_pruned_skewed(Y,mu_tm1_tm1,Sigma_tm1_tm1,Gamma_tm1_tm1,nu_tm1_tm1,Delta_tm1_tm1,G,R,F,mu_eta,Sigma_eta,Gamma_eta,nu_eta,Delta_eta,Sigma_eps,kalman_tol,prune_tol,mvnlogcdf,console_mode)
% -------------------------------------------------------------------------
% Skewed Kalman state and shock smoother for linear state space model
% with skew normally distributed innovations and normally distributed noise:
%   α(t) = G*α(t-1) + R*η(t)   [state transition equation]
%   y(t) = F*α(t)   + ε(t)     [observation equation]
%   η(t) ~ CSN(mu_eta,Sigma_eta,Gamma_eta,nu_eta=0,Delta_eta=I) [innovations, shocks]
%   ε(t) ~ N(mu_eps=0,Sigma_eps)                                [noise, measurement error]
% Dimensions:
%   α(t) is (x_nbr by 1) state vector
%   y(t) is (y_nbr by 1) control vector, i.e. observable variables
%   η(t) is (eta_nbr by 1) vector of shocks
%   ε(t) is (y_nbr by 1) vector of measurement errors
% Notes:
% - mu_eta is set such that E[η(t)] = 0
% - This smoother works for any nu_eta and Delta_eta, but we currently restrict
%   the shocks to be independent skew normally distributed, which is a special
%   case of the closed skew normal distribution with nu_eta=0 and Delta_eta=I
% -------------------------------------------------------------------------
% INPUTS
% - Y               [y_nbr by obs_nbr]             matrix with data
% - mu_tm1_tm1      [x_nbr by 1]                   mu_0_0: initial filtered value of location parameter of CSN distributed state vector (does not equal expectation vector unless Gamma_0_0=0)
% - Sigma_tm1_tm1   [x_nbr by x_nbr]               Sigma_0_0: initial filtered value of scale parameter of CSN distributed state vector (does not equal covariance matrix unless Gamma_0_0=0)
% - Gamma_tm1_tm1   [skew_dim by x_nbr]            Gamma_0_0: initial filtered value of skewness shape parameter of CSN distributed states vector (if 0 then CSN reduces to Gaussian)
% - nu_tm1_tm1      [skew_dim by x_nbr]            nu_0_0: initial filtered value of skewness conditioning parameter of CSN distributed states vector (enables closure of CSN distribution under conditioning, irrelevant if Gamma_0_0=0)
% - Delta_tm1_tm1   [skew_dim by skew_dim]         Delta_0_0: initial filtered value of skewness marginalization parameter of CSN distributed states vector (enables closure of CSN distribution under marginalization, irrelevant if Gamma_0_0=0)
% - G               [x_nbr by x_nbr]               state transition matrix mapping previous states to current states
% - R               [x_nbr by eta_nbr]             state transition matrix mapping current innovations to current states
% - F               [y_nbr by x_nbr]               observation equation matrix mapping current states into current observables
% - mu_eta          [eta_nbr by 1]                 value of location parameter of CSN distributed shocks (does not equal expectation vector unless Gamma_eta=0)
% - Sigma_eta       [eta_nbr by eta_nbr]           value of scale parameter of CSN distributed shocks (does not equal covariance matrix unless Gamma_eta=0)
% - Gamma_eta       [skeweta_dim by eta_nbr]       value of skewness shape parameter of CSN distributed shocks (if 0 then CSN reduces to Gaussian)
% - nu_eta          [skeweta_dim by eta_nbr]       value of skewness conditioning parameter of CSN distributed shocks (enables closure of CSN distribution under conditioning)
% - Delta_eta       [skeweta_dim by skeweta_dim]   value of skewness marginalization parameter of CSN distributed shocks (enables closure of CSN distribution under marginalization)
% - Sigma_eps       [y_nbr by y_nbr]               scale parameter of normally distributed measurement errors (equals covariance matrix), if no measurement errors this is a zero scalar
% - kalman_tol      [double scalar]                tolerance parameter for invertibility of covariance matrix
% - prune_tol       [double]                       threshold to prune redundant skewness dimensions, if set to 0 no pruning will be done
% - mvnlogcdf       [string]                       name of function to compute log Gaussian cdf, possible values: 'gaussian_log_mvncdf_mendell_elston', 'mvncdf'
% - console_mode    [boolean]                      display mode for waitbar
% -------------------------------------------------------------------------
% OUTPUTS
% - x_t_T           [x_nbr by obs_nbr]             smoothed state estimates (mean of CSN distribution)
% - epsilonhat      [y_nbr by obs_nbr]             smoothed measurement errors
% - eta_t_T         [eta_nbr by obs_nbr]           smoothed shock estimates (mean of CSN distribution)
% - x_t_t           [x_nbr by obs_nbr]             filtered state estimates
% - atilde          []                             (not used, included for compatibility)
% - P               [x_nbr by x_nbr by obs_nbr+1]  one-step ahead forecast error variance matrices
% - aK              []                             (not used yet, included for future use and for compatibility)
% - PK              []                             (not used yet, included for future use and for compatibility)
% - decomp          []                             (not used yet, included for future use and for compatibility)
% - V               []                             (not used yet, included for future use and for compatibility)
% - aalphahat       []                             (not used yet, included for future use and for compatibility)
% - eetahat         []                             (not used yet, included for future use and for compatibility)
% - x_0_T           [x_nbr by 1]                   initial smoothed state
% - state_uncertainty0 []                           (not used yet, included for future use and for compatibility)

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

% initializations
atilde = [];
P = [];  % 3D array of one-step ahead forecast error variance matrices
aK = []; % 3D array of k step ahead filtered state variables (a_{t+k|t)
PK = []; % 4D array of k-step ahead forecast error variance matrices
decomp = []; % decomposition of the effect of shocks on filtered values
V = []; % 3D array of state uncertainty matrices
aalphahat = [];%     filtered states in t-1|t
eetahat = []; %       updated shocks in t|t
state_uncertainty0 =[] ;
% get dimensions
x_nbr = length(G);
eta_nbr = size(mu_eta, 1);
[~, obs_nbr] = size(Y);
if rcond(G) < kalman_tol
    G_singular = true;
else
    G_singular = false;
end
mu_eta = R*mu_eta;
Sigma_eta = R*Sigma_eta*R';
Gamma_eta = Gamma_eta/(R'*R)*R';
Gamma_eta_X_Sigma_eta = Gamma_eta*Sigma_eta;
Delta22_common = Delta_eta + Gamma_eta_X_Sigma_eta*Gamma_eta';
A = [-inv(R'*R)*R'*G, inv(R'*R)*R'];

mu_0 = mu_tm1_tm1;
Sigma_0 = Sigma_tm1_tm1;
Gamma_0 = Gamma_tm1_tm1;
Delta_0 = Delta_tm1_tm1;

% initialize CSN parameters of predicted and filtered states;
% use cell for skewness parameters as dimensions are time-varying
x_t_t = nan(x_nbr,obs_nbr);
x_t_T = nan(x_nbr,obs_nbr);
mu_t_tm1 = nan(x_nbr, obs_nbr+1);
mu_t_t = nan(x_nbr, obs_nbr);
mu_t_T = nan(x_nbr, obs_nbr);
mu_yt = nan(eta_nbr, obs_nbr);
mu_jt = nan(2*x_nbr, obs_nbr);
Sigma_t_tm1 = nan(x_nbr, x_nbr, obs_nbr+1);
Sigma_t_t = nan(x_nbr, x_nbr, obs_nbr);
Sigma_t_T = nan(x_nbr, x_nbr, obs_nbr);
Sigma_yt = nan(eta_nbr, eta_nbr, obs_nbr);
Sigma_jt = nan(2*x_nbr, 2*x_nbr, obs_nbr);
Gamma_t_tm1 = cell(obs_nbr+1,1);
Gamma_t_t = cell(obs_nbr,1);
Gamma_t_T = cell(obs_nbr,1);
Gamma_yt = cell(obs_nbr,1);
Gamma_jt = cell(obs_nbr,1);
nu_t_tm1 = cell(obs_nbr+1,1);
nu_t_t = cell(obs_nbr,1);
nu_t_T = cell(obs_nbr,1);
Delta_t_tm1 = cell(obs_nbr+1,1);
Delta_t_t = cell(obs_nbr,1);
Delta_yt = cell(obs_nbr,1);
Delta_jt = cell(obs_nbr,1);
Delta_tilde_t = cell(obs_nbr,1);
N_t = nan(size(Gamma_eta,1),x_nbr,obs_nbr);
M_t = Sigma_t_tm1;
L_t = Sigma_t_tm1;
O_t = cell(obs_nbr,1);
eta_t_T = nan(eta_nbr, obs_nbr);

%%%%%%%%%%
% Filter %
%%%%%%%%%%
waitbar_str = sprintf('SMOOTHER FORWARD PASS\nComputing the mean of CSN distributed filtered states takes time\nConsider skipping it with ''skewed𛲖kalman𛲖smoother𛲖skip'' option');
waitbar_title = 'Pruned Skewed Kalman Smoother';
[hh_fig, length_of_old_string] = wait_bar.run(0, [], waitbar_str, console_mode, 0, waitbar_title);

for t = 1:(obs_nbr+1)
    % prediction
    if t==1
        Gamma_tm1_tm1_X_Sigma_tm1_tm1 = Gamma_tm1_tm1*Sigma_tm1_tm1;
        Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT = Gamma_tm1_tm1_X_Sigma_tm1_tm1*G';
        mu_t_tm1(:,t) = G*mu_tm1_tm1 + mu_eta;
        Sigma_t_tm1(:,:,t) = G*Sigma_tm1_tm1*G' + Sigma_eta;
        Sigma_t_tm1(:,:,t) = 0.5*(Sigma_t_tm1(:,:,t) + Sigma_t_tm1(:,:,t)'); % ensure symmetry
        if G_singular
            invSigma_t_tm1 = pinv(Sigma_t_tm1(:,:,t));
        else
            invSigma_t_tm1 = inv(Sigma_t_tm1(:,:,t));
        end
        Delta11_t_tm1 = Delta_tm1_tm1 + Gamma_tm1_tm1_X_Sigma_tm1_tm1*Gamma_tm1_tm1' - Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT*invSigma_t_tm1*Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT';
        Delta22_t_tm1 = Delta22_common - Gamma_eta_X_Sigma_eta*invSigma_t_tm1*Gamma_eta_X_Sigma_eta';
        Delta12_t_tm1 = -Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT*invSigma_t_tm1*Gamma_eta_X_Sigma_eta';
    else
        Gamma_tm1_tm1_X_Sigma_tm1_tm1 = Gamma_t_t{t-1}*Sigma_t_t(:,:,t-1);
        Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT = Gamma_tm1_tm1_X_Sigma_tm1_tm1*G';
        mu_t_tm1(:,t) = G*mu_t_t(:,t-1) + mu_eta;
        Sigma_t_tm1(:,:,t) = G*Sigma_t_t(:,:,t-1)*G' + Sigma_eta;
        Sigma_t_tm1(:,:,t) = 0.5*(Sigma_t_tm1(:,:,t) + Sigma_t_tm1(:,:,t)'); % ensure symmetry
        if G_singular
            invSigma_t_tm1 = pinv(Sigma_t_tm1(:,:,t));
        else
            invSigma_t_tm1 = inv(Sigma_t_tm1(:,:,t));
        end
        Delta11_t_tm1 = Delta_t_t{t-1} + Gamma_tm1_tm1_X_Sigma_tm1_tm1*Gamma_t_t{t-1}' - Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT*invSigma_t_tm1*Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT';
        Delta22_t_tm1 = Delta22_common - Gamma_eta_X_Sigma_eta*invSigma_t_tm1*Gamma_eta_X_Sigma_eta';
        Delta12_t_tm1 = -Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT*invSigma_t_tm1*Gamma_eta_X_Sigma_eta';
    end
    Gamma_t_tm1{t} = [Gamma_tm1_tm1_X_Sigma_tm1_tm1_X_GT; Gamma_eta_X_Sigma_eta]*invSigma_t_tm1;
    nu_t_tm1{t} = [nu_tm1_tm1; nu_eta];
    Delta_t_tm1{t} = [Delta11_t_tm1 , Delta12_t_tm1; Delta12_t_tm1' , Delta22_t_tm1];
    Delta_t_tm1{t} = 0.5*(Delta_t_tm1{t} + Delta_t_tm1{t}'); % ensure symmetry
    % Omega is Gaussian prediction error covariance
    Omega = F*Sigma_t_tm1(:,:,t)*F' + Sigma_eps;
    invOmega = inv(Omega); % compute the inverse of Omega directly
    K_Gauss = Sigma_t_tm1(:,:,t)*F'*invOmega;
    K_Skewed = Gamma_t_tm1{t}*K_Gauss; % Skewed Kalman Gain

    if t <= obs_nbr
        % update
        prediction_error = Y(:,t)-F*mu_t_tm1(:,t);
        mu_t_t(:,t) = mu_t_tm1(:,t) + K_Gauss*prediction_error;
        Sigma_t_t(:,:,t) = Sigma_t_tm1(:,:,t) - K_Gauss*F*Sigma_t_tm1(:,:,t);
        Gamma_t_t{t} = Gamma_t_tm1{t};
        nu_t_t{t} = nu_t_tm1{t} - K_Skewed*prediction_error;
        Delta_t_t{t} = Delta_t_tm1{t};
        Sigma_t_t(:,:,t) = 0.5*(Sigma_t_t(:,:,t) + Sigma_t_t(:,:,t)'); % ensure symmetry
        Delta_t_t{t} = 0.5*(Delta_t_t{t} + Delta_t_t{t}'); % ensure symmetry
        %[Sgm, Gam, nu, Dlt] = csn_prune_distribution(Sigma_t_t(:,:,t),Gamma_t_t{t},nu_t_t{t},Delta_t_t{t},prune_tol);
        %x_t_t(:,t) = csn_mean(mu_t_t(:,t), Sgm, Gam, nu, Dlt, mvnlogcdf);

        % assign for next time step
        mu_tm1_tm1 = mu_t_t(:,t);
        Sigma_tm1_tm1 = Sigma_t_t(:,:,t);
        Gamma_tm1_tm1 = Gamma_t_t{t};
        nu_tm1_tm1 = nu_t_t{t};
        Delta_tm1_tm1 = Delta_t_t{t};
    end

    [~, length_of_old_string] = wait_bar.run(t/(obs_nbr+1), hh_fig, waitbar_str, console_mode, length_of_old_string);

end
wait_bar.close(hh_fig, console_mode);

waitbar_str = sprintf('SMOOTHER BACKWARD PASS\nComputing the mean of CSN distributed smoothed states and shocks takes time\nConsider skipping it with ''skewed𛲖kalman𛲖smoother𛲖skip'' option');
[hh_fig, length_of_old_string] = wait_bar.run(0, [], waitbar_str, console_mode, 0, waitbar_title);

% smoothing step for the last time point, i.e. 'T'
mu_t_T(:,obs_nbr) = mu_t_t(:,obs_nbr);
Sigma_t_T(:,:,obs_nbr) = Sigma_t_t(:,:,obs_nbr);
Gamma_t_T{obs_nbr} = Gamma_t_t{obs_nbr};
nu_t_T{obs_nbr} = nu_t_t{obs_nbr};
Delta_t_T{obs_nbr} = Delta_t_t{obs_nbr};
[Sgm, Gam, nu, Dlt] = csn_prune_distribution(Sigma_t_T(:,:,obs_nbr),Gamma_t_T{obs_nbr},nu_t_T{obs_nbr},Delta_t_T{obs_nbr},prune_tol);
x_t_T(:,obs_nbr) = csn_mean(mu_t_T(:,obs_nbr), Sgm, Gam, nu, Dlt, mvnlogcdf);

% smoothing
for t = (obs_nbr-1):-1:1
    Jt = Sigma_t_t(:,:,t)*G'*pinv(Sigma_t_tm1(:,:,t+1));
    mu_t_T(:,t) = mu_t_t(:,t) + Jt*(mu_t_T(:,t+1)-mu_t_tm1(:,t+1));
    mu_jt(:,t) = [mu_t_t(:,t) + Jt*(mu_t_T(:,t+1)-mu_t_tm1(:,t+1)); mu_t_T(:,t+1)];
    Sigma_t_T(:,:,t) = Sigma_t_t(:,:,t) + Jt*(Sigma_t_T(:,:,t+1)-Sigma_t_tm1(:,:,t+1))*Jt';
    Sigma_jt(:,:,t) = [Sigma_t_t(:,:,t)+Jt*(Sigma_t_T(:,:,t+1)-Sigma_t_tm1(:,:,t+1))*Jt', Jt*Sigma_t_T(:,:,t+1);Sigma_t_T(:,:,t+1)*Jt', Sigma_t_T(:,:,t+1)];
    M_t(:,:,t) = Sigma_t_T(:,:,t+1)*Jt'*pinv(Sigma_t_T(:,:,t));
    N_t(:,:,t) = -Gamma_eta*G + Gamma_eta*M_t(:,:,t);
    if t == (obs_nbr-1)
        O_t{t} = N_t(:,:,t);
    else
        O_t{t} = [N_t(:,:,t); O_t{t+1}*M_t(:,:,t)];
    end
    if t == (obs_nbr-1)
        Delta_tilde_t{t} = Delta_eta + Gamma_eta*(Sigma_t_T(:,:,t+1) - M_t(:,:,t)*Sigma_t_T(:,:,t)*M_t(:,:,t)')*Gamma_eta';
        Gamma_jt{t} = [Gamma_t_t{t}, zeros(size(Gamma_t_t{t},1), size(Gamma_eta, 2)); -Gamma_eta*G, Gamma_eta];
    else
        temp_mat = [Gamma_eta; O_t{t+1}];
        Delta_tilde_t{t} = blkdiag(Delta_eta, Delta_tilde_t{t+1}) + temp_mat*(Sigma_t_T(:,:,t+1) - M_t(:,:,t)*Sigma_t_T(:,:,t)*M_t(:,:,t)')*temp_mat';
        Gamma_jt{t} = [Gamma_t_t{t}, zeros(size(Gamma_t_t{t},1), size(Gamma_eta, 2)); -Gamma_eta*G, Gamma_eta; zeros(size(O_t{t+1},1), size(Gamma_eta, 2)), O_t{t+1}];
    end
    Gamma_t_T{t} = [Gamma_t_t{t}; O_t{t}];
    nu_t_T{t} = nu_t_t{obs_nbr, 1};
    Delta_t_T{t} = blkdiag(Delta_t_t{t}, Delta_tilde_t{t});
    Delta_jt{t} = blkdiag(Delta_t_t{t}, Delta_eta, Delta_tilde_t{t+1});
    mu_yt(:,t) = A*mu_jt(:,t);
    Sigma_yt(:,:,t) = A*Sigma_jt(:,:,t)*A';
    Gamma_X_Sigma_At = Gamma_jt{t}*Sigma_jt(:,:,t)*A';
    Gamma_yt{t} = Gamma_X_Sigma_At * pinv(Sigma_yt(:,:,t));
    Delta_yt{t} = Delta_jt{t}+Gamma_jt{t}*Sigma_jt(:,:,t)*Gamma_jt{t}'-Gamma_X_Sigma_At*pinv(Sigma_yt(:,:,t))*Gamma_X_Sigma_At';
    [Sgm, Gam, nu, Dlt] = csn_prune_distribution(Sigma_t_T(:,:,t),Gamma_t_T{t},nu_t_T{t},Delta_t_T{t},prune_tol);
    x_t_T(:,t) = csn_mean(mu_t_T(:,t), Sgm, Gam, nu, Dlt, mvnlogcdf);
    [Sgm, Gam, nu, Dlt] = csn_prune_distribution(Sigma_yt(:,:,t),Gamma_yt{t},nu_t_T{t},Delta_yt{t},prune_tol);
    eta_t_T(:,t+1) = csn_mean(mu_yt(:,t), Sgm, Gam, nu, Dlt, mvnlogcdf);

    [~, length_of_old_string] = wait_bar.run((obs_nbr-t)/obs_nbr, hh_fig, waitbar_str, console_mode, length_of_old_string);

end

% first period for eta_t_T
J0 = Sigma_0*G'*pinv(Sigma_t_tm1(:,:,1));
mu_0_T = mu_0 + J0*(mu_t_T(:,1)-mu_t_tm1(:,1));
mu_j0 = [mu_0 + J0*(mu_t_T(:,1)-mu_t_tm1(:,1)); mu_t_T(:,1)];
Sigma_0_T = Sigma_0 + J0*(Sigma_t_T(:,:,1)-Sigma_t_tm1(:,:,1))*J0';
Sigma_j0 = [Sigma_0+J0*(Sigma_t_T(:,:,1)-Sigma_t_tm1(:,:,1))*J0', J0*Sigma_t_T(:,:,1);Sigma_t_T(:,:,1)*J0', Sigma_t_T(:,:,1)];
M_0 = Sigma_t_T(:,:,1)*J0'*pinv(Sigma_0_T);
N_0 = -Gamma_eta*G + Gamma_eta*M_0;
O_0 = [N_0; O_t{1}*M_0];
temp_mat = [Gamma_eta; O_t{1}];
Delta_tilde_0 = blkdiag(Delta_eta, Delta_tilde_t{1}) + temp_mat*(Sigma_t_T(:,:,1) - M_0*Sigma_0_T*M_0')*temp_mat';
Gamma_j0 = [Gamma_0, zeros(size(Gamma_0,1), size(Gamma_eta, 2)); -Gamma_eta*G, Gamma_eta; zeros(size(O_t{1},1), size(Gamma_eta, 2)), O_t{1}];
Gamma_0_T = [Gamma_0; O_0];
nu_0_T = nu_t_t{obs_nbr, 1};
Delta_0_T = blkdiag(Delta_0, Delta_tilde_0);
Delta_j0 = blkdiag(Delta_0, Delta_eta, Delta_tilde_t{1});
mu_y0 = A*mu_j0;
Sigma_y0 = A*Sigma_j0*A';
Gamma_X_Sigma_At = Gamma_j0*Sigma_j0*A';
Gamma_y0 = Gamma_X_Sigma_At * pinv(Sigma_y0);
Delta_y0 = Delta_j0+Gamma_j0*Sigma_j0*Gamma_j0'-Gamma_X_Sigma_At*pinv(Sigma_y0)*Gamma_X_Sigma_At';
[Sgm, Gam, nu, Dlt] = csn_prune_distribution(Sigma_0_T,Gamma_0_T,nu_0_T,Delta_0_T,prune_tol);
x_0_T = csn_mean(mu_0_T, Sgm, Gam, nu, Dlt, mvnlogcdf);
[Sgm, Gam, nu, Dlt] = csn_prune_distribution(Sigma_y0,Gamma_y0,nu_0_T,Delta_y0,prune_tol);
eta_t_T(:,1) = csn_mean(mu_y0, Sgm, Gam, nu, Dlt, mvnlogcdf);

wait_bar.close(hh_fig, console_mode);

epsilonhat = Y-F*x_t_T;

end