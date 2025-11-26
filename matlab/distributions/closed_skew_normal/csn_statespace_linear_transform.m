function [mu_tr, Sigma_tr, Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform(G_bar, mu_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar)
% [mu_tr, Sigma_tr, Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform(G_bar, mu_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar)
% -------------------------------------------------------------------------
% Compute CSN parameters for joint distribution of [x_t', eta_t']' from
% state transition equation: x(t) = G*x(t-1) + R*eta(t)
% 1) make proper csn random variable first (in case of singularity)
% 2) do either rank deficient or full rank linear transformation
% -----------------------------------------------------------------------
% INPUTS
%   - G_bar = [G, R] multiplying the joint distribution of x_t and eta_t
%   - mu_bar = [mu_t_t; mu_eta]
%   - nu_bar = [nu_t_t; nu_eta];
%   - Sigma_bar = blkdiag_two(Sigma_t_t, Sigma_eta);
%   - Gamma_bar = blkdiag_two(Gamma_t_t, Gamma_eta);
%   - Delta_bar = blkdiag_two(Delta_t_t, Delta_eta);
% -----------------------------------------------------------------------
% OUTPUT
% transformed parameters of CSN joint distribution of x_t and eta_t:
% mu_tr, Sigma_tr, Gamma_tr, nu_tr, Delta_tr

% Copyright © 2023 Gaygysyz Guljanov
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

tol = 1e-8;

mu_tr = G_bar * mu_bar;
Sigma_tr = G_bar * Sigma_bar * G_bar';
mu_bar = zeros(size(mu_bar)); % make mu_bar a zero vector

%% singular skew normal to closed skew normal
if rcond(Sigma_bar) < 1e-10
    % Sigma_bar is rank deficient
    [S_hat, ~, Sigma_bar, Gamma_bar] = csn_singular_to_csn_proper(mu_bar, Sigma_bar, Gamma_bar, tol, 0);
    G_bar = G_bar * S_hat;
end

%% Linear Transformation
% Check if A_mat is rank deficient
dimensions = size(G_bar);
if dimensions(1) > dimensions(2)
    is_singular = rcond(G_bar' * G_bar) < 1e-8;
else
    is_singular = rcond(Sigma_tr) < 1e-8;
end
% If not full rank, use rank deficient linear transformation
if is_singular
    [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_rank_deficient(G_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar, tol);
else
    [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_full_rank(G_bar, Sigma_bar, Sigma_tr, Gamma_bar, nu_bar, Delta_bar);
end

function [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_rank_deficient(G_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar, tol)
% [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_rank_deficient(G_bar, Sigma_bar, Gamma_bar, nu_bar, Delta_bar, tol)
% -------------------------------------------------------------------------
% Compute CSN parameters for joint distribution of [x_t', eta_t']' from
% state transition equation: x(t) = G*x(t-1) + R*eta(t)
% helper function for rank deficient linear transformation of CSN distributed random variable
% -----------------------------------------------------------------------
% INPUTS
%   - G_bar = [G, R] multiplying the joint distribution of x_t and eta_t
%   - Sigma_bar = blkdiag_two(Sigma_t_t, Sigma_eta);
%   - Gamma_bar = blkdiag_two(Gamma_t_t, Gamma_eta);
%   - nu_bar = [nu_t_t; nu_eta];
%   - Delta_bar = blkdiag_two(Delta_t_t, Delta_eta);
%   - tol: tolerance level showing which values should be assumed to be numerically zero
% -----------------------------------------------------------------------
% OUTPUT
% transformed parameters of rank-deficient CSN joint distribution of x_t and eta_t:
% Gamma_tr, nu_tr, Delta_tr (returns only skewness parameters, as other parameters are easy to obtain)

[nrow, ncol] = size(G_bar);
difference = abs(nrow - ncol);
[S_mat, Lambda_mat, T_mat] = svd(G_bar); % apply singular value decomposition
if nrow >= ncol % more rows
    S_mat = S_mat(:, 1:end-difference);
    Lambda_mat = Lambda_mat(1:end-difference, :);
else
    T_mat = T_mat(:, 1:end-difference);
    Lambda_mat = Lambda_mat(:, 1:end-difference);
end    
hold_vec = abs(diag(Lambda_mat)) > tol;      % which dimensions to delete
S_mat = S_mat(:, hold_vec);                  % delete respective columns
Lambda_mat = Lambda_mat(hold_vec, hold_vec); % delete respective rows and columns
T_mat = T_mat(:, hold_vec);                  % delete respective columns

% Apply the final formulas of theorem
Tmp_mat = T_mat' * Sigma_bar * T_mat;
Tmp_mat2 = Gamma_bar * Sigma_bar;
Tmp_mat3 = Tmp_mat2 * T_mat / Tmp_mat;
Gamma_tr = (Tmp_mat3 / Lambda_mat) / (S_mat' * S_mat) * S_mat';
nu_tr = nu_bar;
Delta_tr = Delta_bar + Tmp_mat2 * Gamma_bar' - Tmp_mat3 * T_mat' * Sigma_bar * Gamma_bar';

end % csn_statespace_linear_transform_rank_deficient


function [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_full_rank(G_bar, Sigma_bar, Sigma_tr, Gamma_bar, nu_bar, Delta_bar)
% [Gamma_tr, nu_tr, Delta_tr] = csn_statespace_linear_transform_full_rank(G_bar, Sigma_bar, Sigma_tr, Gamma_bar, nu_bar, Delta_bar)
% -------------------------------------------------------------------------
% Compute CSN parameters for joint distribution of [x_t', eta_t']' from
% state transition equation: x(t) = G*x(t-1) + R*eta(t)
% helper function for full rank linear transformation of CSN distributed random variable
% -----------------------------------------------------------------------
% INPUTS
%   - G_bar = [G, R] multiplying the joint distribution of x_t and eta_t
%   - Sigma_bar = blkdiag_two(Sigma_t_t, Sigma_eta);
%   - Sigma_tr : transformed (rank deficient) parameter matrix
%   - Gamma_bar = blkdiag_two(Gamma_t_t, Gamma_eta);
%   - nu_bar = [nu_t_t; nu_eta];
%   - Delta_bar = blkdiag_two(Delta_t_t, Delta_eta);
% -----------------------------------------------------------------------
% OUTPUT
% transformed parameters of rank-deficient CSN joint distribution of x_t and eta_t:
% Gamma_tr, nu_tr, Delta_tr (returns only skewness parameters, as other parameters are easy to obtain)

[nrow, ncol] = size(G_bar);
% more rows
if nrow >= ncol
    nu_tr = nu_bar;
    Delta_tr = Delta_bar;
    if nrow == ncol
        Gamma_tr = Gamma_bar / G_bar;
        return
    end
    Gamma_tr = Gamma_bar / (G_bar' * G_bar) * G_bar';
    return
end
% more columns
Tmp_mat = Sigma_bar * G_bar' / Sigma_tr;
Gamma_tr = Gamma_bar * Tmp_mat;
nu_tr = nu_bar;
Delta_tr = Delta_bar + Gamma_bar * Sigma_bar * Gamma_bar' - Gamma_bar * Tmp_mat * G_bar * Sigma_bar * Gamma_bar';

end % csn_statespace_linear_transform_full_rank


end % csn_statespace_linear_transform