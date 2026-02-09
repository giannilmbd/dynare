function csn = csn_update_specification(Var, Skew)
% csn = csn_update_specification(Var, Skew)
% -------------------------------------------------------------------------
% Recover CSN parameters from variance and skewness and set location parameter such that shocks have mean of zero
%
% CURRENT IMPLEMENTATION:
% - Handles independent univariate skew normals (diagonal Var, diagonal Skew tensor)
% - Each shock j has univariate SN(xi_j, omega_j, alpha_j)
% - Equivalent to Gupta et al. MSN with diagonal Lambda
% - Represented as CSN_{exo_nbr,exo_nbr}(mu,Sigma,Gamma,0,I) where:
%   * Sigma = diag(omega_1^2,...,omega_exo_nbr^2)
%   * Gamma = diag(alpha_1/omega_1,...,alpha_exo_nbr/omega_exo_nbr)
%   * mu chosen so E[shocks]=0
% -------------------------------------------------------------------------
% HIERARCHY OF SKEW NORMAL DISTRIBUTIONS:
%
% 1. Azzalini (1985): Univariate skew normal SN(xi, omega, alpha)
%    - Single variable with location xi, scale omega, shape alpha
%    - Special case of CSN_{1,1}(xi, omega^2, alpha/omega, 0, 1)
%
% 2. Azzalini & Dalla Valle (1996): Multivariate skew normal MSN(mu, Omega, alpha)
%    - Uses vector alpha for shape (single latent factor controls all skewness)
%    - LIMITATION: Cannot represent independent skew normals with different skewness
%    - Special case of CSN_{d,1}(mu, Omega, Gamma, 0, 1) where:
%      * Gamma is 1×d row vector: Gamma = alpha'*Omega^(-1/2)
%      * k=1 (single latent factor), nu=0, Delta=1
%
% 3. Gupta, González-Farías & Domínguez-Molina (2004): Generalized MSN(mu, Omega, Lambda)
%    - Uses k×d matrix Lambda for shape (k latent factors)
%    - SPECIAL CASES:
%      * Lambda = alpha*1' (rank-1): reduces to Azzalini-Dalla Valle
%      * Lambda = diag(alpha_1,...,alpha_d): independent skew normals ← OUR CASE
%      * Lambda general: complex skewness patterns
%    - Special case of CSN_{d,k}(mu, Omega, Gamma, 0, I_k) where:
%      * Gamma = Lambda*Omega^(-1/2) is the k×d shape matrix
%      * nu = 0 (no truncation), Delta = I_k (identity)
%
% 4. González-Farías, Domínguez-Molina & Gupta (2004): Closed Skew Normal CSN_{d,k}(mu,Sigma,Gamma,nu,Delta)
%    - Most general formulation
%    - Gupta et al. MSN is CSN_{d,k}(mu,Sigma,Gamma,0,I_k) with Gamma=Lambda*Sigma^(-1/2)
%    - Gaussian is CSN with Gamma=0 (nu and Delta become irrelevant)
%
% PARAMETER RELATIONSHIPS TO CSN_{d,k}(mu,Sigma,Gamma,nu,Delta):
% ┌────────────────────────┬───┬───┬───────────────────-──┬─────┬───────┐
% │ Distribution           │ d │ k │ Gamma                │ nu  │ Delta │
% ├────────────────────────┼───┼───┼──────────────────-───┼─────┼───────┤
% │ Azzalini (1985)        │ 1 │ 1 │ alpha/omega          │ 0   │ 1     │
% │ Azzalini-DV (1996)     │ d │ 1 │ alpha'*Omega^(-1/2)  │ 0   │ 1     │
% │ Gupta et al. (2004)    │ d │ k │ Lambda*Omega^(-1/2)  │ 0   │ I_k   │
% │ - Independent case     │ d │ d │ diag(alpha_i/omega_i)│ 0   │ I_d   │
% │ Gaussian               │ d │ - │ 0                    │ -   │ -     │
% └────────────────────────┴───┴───┴─────────────────────-┴─────┴───────┘
%
% -------------------------------------------------------------------------
% INPUTS
% - Var      [double]   unconditional variance matrix (exo_nbr × exo_nbr)
% - Skew     [double]   sparse skewness representation (N × 4 matrix)
%                       each row is [i, j, k, value] storing a non-zero element
%                       of the coskewness tensor with all index permutations
%   Note: Currently only supports diagonal Var and Skew tensors (no covariance/coskewness)
%         This can be generalized but requires careful design and numerical implementation
% -------------------------------------------------------------------------
% OUTPUTS
% - csn      [struct]   CSN parameters structure with fields:
%   * mu_e    [exo_nbr×1]       location parameter (ensures zero mean shocks)
%   * Sigma_e [exo_nbr×exo_nbr] scale matrix
%   * Gamma_e [exo_nbr×exo_nbr] shape matrix
%   * nu_e    [exo_nbr×1]       truncation parameter (=0 for MSN)
%   * Delta_e [exo_nbr×exo_nbr] correlation matrix (=I for MSN)
% -------------------------------------------------------------------------
% REFERENCES
% - Azzalini, A. (1985). "A class of distributions which includes the normal ones"
%   Scandinavian Journal of Statistics, 12(2), 171-178.
% - Azzalini, A. & Dalla Valle, A. (1996). "The multivariate skew-normal distribution"
%   Biometrika, 83(4), 715-726.
% - Gupta, A.K., González-Farías, G. & Domínguez-Molina, J.A. (2004).
%   "A multivariate skew normal distribution"
%   Journal of Multivariate Analysis, 89(1), 181-190.
% - González-Farías, G., Domínguez-Molina, J.A. & Gupta, A.K. (2004).
%   "Addendum: The closed skew-normal distribution"
%   In: Genton, M.G. (Ed.), Skew-Elliptical Distributions and Their Applications, pp. 25-42.
% - Azzalini, A. & Capitanio, A. (2014). "The Skew-Normal and Related Families"
%   Cambridge University Press

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

exo_nbr = size(Var, 1);
csn.mu_e = zeros(exo_nbr,1);
csn.Sigma_e = Var;
csn.Gamma_e = zeros(exo_nbr,exo_nbr);
csn.nu_e = zeros(exo_nbr,1);
csn.Delta_e = eye(exo_nbr);
if size(Skew, 1) > 0
    if ~isdiag(Var)
        error('csn_update_specification: Skewness is currently only supported for independent skew normally distributed shocks, i.e. they cannot be correlated. Remove the corr parameters.')
    end
    % Check for co-skewness (off-diagonal entries where not all 3 indices are equal)
    if any(Skew(:,1) ~= Skew(:,2) | Skew(:,2) ~= Skew(:,3))
        error('csn_update_specification: Skewness is currently only supported for independent skew normally distributed shocks, i.e. there cannot be co-skewness between shocks.')
    end
    sqrtTwoPi = sqrt(2/pi);
    for jexo = 1:exo_nbr
        idx = Skew(:,1)==jexo & Skew(:,2)==jexo & Skew(:,3)==jexo; % lookup (jexo,jexo,jexo) in sparse Skew using element-wise comparison
        if any(idx)
            skew_val = Skew(idx,4);
            if abs(skew_val) > eps % univariate skew normal distribution, 0 would be Gaussian
                [omega_e, alpha_e] = sn_var_skew_to_scale_shape_univariate(Var(jexo,jexo), skew_val);
                csn.Sigma_e(jexo,jexo) = omega_e^2;
                csn.Gamma_e(jexo,jexo) = alpha_e/omega_e;
                % set location parameter such that E[shocks]=0
                csn.mu_e(jexo) = -omega_e*alpha_e/sqrt(1+alpha_e^2)*sqrtTwoPi;
            end
        end
    end
end

function [omega, alpha] = sn_var_skew_to_scale_shape_univariate(Var, Skew)
% invert univariate moment formulas of skew normal distribution of Azzalini (1985)
% see e.g. equations (2.22) to (2.28) in Azzalini and Capitanio (2014)
% or the discussion in Gupta, González-Farías & Domínguez-Molina (2004)
absSkew2_3 = abs(Skew)^(2/3);
delta = sqrt( pi/2 * absSkew2_3 / ( absSkew2_3 + ((4-pi)/2)^(2/3) ) );
alpha = sign(Skew) * delta/sqrt(1-delta^2); % alpha
omega = sqrt(Var/(1-2/pi*delta^2)); % omega
end % sn_var_skew_to_scale_shape_univariate

end % csn_update_specification