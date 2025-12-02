function log_cdf = gaussian_log_mvncdf_mendell_elston(Zj, R, cutoff)
% log_cdf = gaussian_log_mvncdf_mendell_elston(Zj, R, cutoff)
% ------------------------------------------------------------------------- 
% Approximates Gaussian log(CDF) function according to Mendell and Elston (1974)
% ------------------------------------------------------------------------- 
% INPUTS 
% - Zj       [n by 1]   column vector of points where Gaussian CDF is evaluated at
% - R        [n by n]   correlation matrix
% - cutoff   [2 by 1]   optional threshold points at which values in Zj are too low/high to be evaluated, defaults to [6 38]
% ------------------------------------------------------------------------- 
% OUTPUTS 
% log_cdf    [double]   approximate value of Gaussian log(CDF)

% Copyright © 2015 Dietmar Bauer (original implementation)
% Copyright © 2022-2023 Gaygysyz Guljanov
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

if nargin < 3
    cutoff = [6 38]; % default values
end
% cutoff tails as these are too low to be evaluated
Zj(Zj>cutoff(1))  = cutoff(1);
Zj(Zj<-cutoff(1)) = -cutoff(1);
n = length(Zj); % remaining dimension of Z 

% first element
cdf_val = phid(Zj(1));
pdf_val = phip(Zj(1));
log_cdf = log(cdf_val); % perform all calculations in logs
 
for jj = 1:(n-1)
    ajjm1 = pdf_val / cdf_val;

    % Update Zj and Rij
    tZ = Zj + ajjm1 * R(:, 1); % update Zj
    R_jj = R(:, 1) * R(1, :);
    tRij = R - R_jj * ( ajjm1 + Zj(1) ) * ajjm1; % update Rij

    % Convert Rij (i.e. Covariance matrix) to Correlation matrix
    COV_TO_CORR = sqrt( diag(tRij) );
    Zj = tZ ./ COV_TO_CORR;
    R = tRij ./ (COV_TO_CORR * COV_TO_CORR');

    % Cutoff those dimensions if they are too low to be evaluated
    Zj(Zj>cutoff(2))  = cutoff(2);  % use larger cutoff
    Zj(Zj<-cutoff(2)) = -cutoff(2); % use larger cutoff

    % Evaluate jj's probability
    cdf_val = phid(Zj(2));
    pdf_val = phip(Zj(2));

    % Delete unnecessary parts of updated Zj and Corr_mat
    last_el = n - jj + 1;
    Zj = Zj(2:last_el);
    R = R(2:last_el, 2:last_el);

    % Overall probability
    log_cdf = log_cdf + log(cdf_val);

end


return % --*-- Unit tests --*--


%@test:1
% univariate case: should equal log(normcdf(z))
try
    z_vals = [-2; -1; 0; 1; 2];
    test_results = true(length(z_vals), 1);
    for i = 1:length(z_vals)
        z = z_vals(i);
        R = 1; % correlation matrix is just 1 for univariate
        log_cdf_me = gaussian_log_mvncdf_mendell_elston(z, R);
        log_cdf_exact = log(normcdf(z));
        test_results(i) = abs(log_cdf_me - log_cdf_exact) < 1e-14;
    end
    t(1) = all(test_results);
catch
    t = false;
end
T = all(t);
%@eof:1

%@test:2
% bivariate case: compare with mvncdf for independent variables (R = I)
try
    Zj = [0.5; -0.3];
    R = eye(2);
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R);
    % for independent variables: Phi(z1, z2) = Phi(z1) * Phi(z2)
    log_cdf_exact = log(normcdf(Zj(1))) + log(normcdf(Zj(2)));
    t(1) = abs(log_cdf_me - log_cdf_exact) < 1e-10;
catch
    t = false;
end
T = all(t);
%@eof:2

%@test:3
% bivariate case: compare with mvncdf for correlated variables
try
    Zj = [1.2; 0.8];
    rho = 0.6;
    R = [1 rho; rho 1];
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R);
    % mvncdf computes P(X <= x) for X ~ N(0, R), so we pass Zj directly
    log_cdf_matlab = log(mvncdf(Zj', zeros(1,2), R));
    % Mendell-Elston is an approximation, allow some tolerance
    t(1) = abs(log_cdf_me - log_cdf_matlab) < 1e-3;
catch
    t = false;
end
T = all(t);
%@eof:3

%@test:4
% trivariate case: compare with mvncdf
try
    Zj = [0.5; -0.2; 1.0];
    R = [1.0  0.3  0.1;
         0.3  1.0 -0.2;
         0.1 -0.2  1.0];
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R);
    log_cdf_matlab = log(mvncdf(Zj', zeros(1,3), R));
    % allow tolerance for approximation
    t(1) = abs(log_cdf_me - log_cdf_matlab) < 2e-3;
catch
    t = false;
end
T = all(t);
%@eof:4

%@test:5
% test at z = 0 for standard bivariate normal with positive correlation
% P(X1 <= 0, X2 <= 0) = 1/4 + arcsin(rho)/(2*pi) for bivariate standard normal
try
    rho = 0.5;
    Zj = [0; 0];
    R = [1 rho; rho 1];
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R);
    % exact formula for P(X1 <= 0, X2 <= 0)
    log_cdf_exact = log(0.25 + asin(rho) / (2*pi));
    t(1) = abs(log_cdf_me - log_cdf_exact) < 3e-3;
catch
    t = false;
end
T = all(t);
%@eof:5

%@test:6
% test cutoff behavior: very large positive values should give log_cdf close to 0
try
    Zj = [10; 10];
    R = [1 0.3; 0.3 1];
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R);
    % CDF at very large values should be close to 1, so log(CDF) close to 0
    t(1) = abs(log_cdf_me) < 1e-8;
catch
    t = false;
end
T = all(t);
%@eof:6

%@test:7
% test cutoff behavior: very large negative values should give very negative log_cdf
try
    Zj = [-6; -6];
    R = [1 0.5; 0.5 1];
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R);
    % CDF at very negative values should be very small, so log(CDF) should be very negative
    t(1) = log_cdf_me < -25;
catch
    t = false;
end
T = all(t);
%@eof:7

%@test:8
% four-dimensional case: compare with mvncdf
try
    Zj = [0.8; -0.4; 0.2; 1.1];
    R = [1.0  0.2  0.1  0.05;
         0.2  1.0  0.3 -0.1;
         0.1  0.3  1.0  0.2;
         0.05 -0.1 0.2  1.0];
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R);
    log_cdf_matlab = log(mvncdf(Zj', zeros(1,4), R));
    % allow tolerance for approximation in higher dimensions
    t(1) = abs(log_cdf_me - log_cdf_matlab) < 2e-3;
catch
    t = false;
end
T = all(t);
%@eof:8

%@test:9
% test with custom cutoff parameter
try
    Zj = [0.5; 0.3];
    R = [1 0.4; 0.4 1];
    cutoff = [5 30];
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R, cutoff);
    log_cdf_matlab = log(mvncdf(Zj', zeros(1,2), R));
    t(1) = abs(log_cdf_me - log_cdf_matlab) < 2e-3;
catch
    t = false;
end
T = all(t);
%@eof:9

%@test:10
% verify monotonicity: larger z values should give larger log_cdf
try
    R = [1 0.3; 0.3 1];
    Zj1 = [0; 0];
    Zj2 = [1; 1];
    Zj3 = [2; 2];
    log_cdf1 = gaussian_log_mvncdf_mendell_elston(Zj1, R);
    log_cdf2 = gaussian_log_mvncdf_mendell_elston(Zj2, R);
    log_cdf3 = gaussian_log_mvncdf_mendell_elston(Zj3, R);
    t(1) = log_cdf1 < log_cdf2 && log_cdf2 < log_cdf3;
catch
    t = false;
end
T = all(t);
%@eof:10

%@test:11
% high-dimensional case (n=50): mvncdf fails for dimensions > 25, but Mendell-Elston works
% test with independent variables where analytical solution is known
try
    n = 50;
    Zj = 0.5 * ones(n, 1); % all elements equal to 0.5
    R = eye(n); % independent variables
    log_cdf_me = gaussian_log_mvncdf_mendell_elston(Zj, R);
    % verify that mvncdf fails for high dimensions (this is why Mendell-Elston is needed)
    mvncdf_failed = false;
    try
        cdf_matlab = mvncdf(Zj', zeros(1,n), R);
    catch ME
        mvncdf_failed = true;
    end
    if ~mvncdf_failed
        error('mvncdf unexpectedly worked for n=%d dimensions', n);
    end
    % for independent variables: Phi(z1, ..., zn) = prod(Phi(zi)) = Phi(0.5)^n
    log_cdf_exact = n * log(normcdf(0.5));
    t(1) = abs(log_cdf_me - log_cdf_exact) < 1e-10;
    % verify result is finite and reasonable
    t(2) = isfinite(log_cdf_me) && log_cdf_me < 0;
catch
    t = false(2, 1);
end
T = all(t);
%@eof:11

%@test:12
% high-dimensional case (n=100): test with correlated variables
% since mvncdf doesn't work, verify properties: finite result, monotonicity
try
    n = 100;
    % create a valid correlation matrix using compound symmetry structure
    rho = 0.1;
    R = (1 - rho) * eye(n) + rho * ones(n, n);
    
    % test at different z values to verify monotonicity
    Zj1 = zeros(n, 1);
    Zj2 = 0.5 * ones(n, 1);
    Zj3 = ones(n, 1);
    
    log_cdf1 = gaussian_log_mvncdf_mendell_elston(Zj1, R);
    log_cdf2 = gaussian_log_mvncdf_mendell_elston(Zj2, R);
    log_cdf3 = gaussian_log_mvncdf_mendell_elston(Zj3, R);
    
    % verify results are finite
    t(1) = isfinite(log_cdf1) && isfinite(log_cdf2) && isfinite(log_cdf3);
    % verify monotonicity: larger z values should give larger log_cdf
    t(2) = log_cdf1 < log_cdf2 && log_cdf2 < log_cdf3;
    % verify results are negative (log of probability < 1)
    t(3) = log_cdf1 < 0 && log_cdf2 < 0 && log_cdf3 < 0;
catch
    t = false(3, 1);
end
T = all(t);
%@eof:12

end % gaussian_log_mvncdf_mendell_elston

%%%%%%%%%%%%%%
% Normal PDF %
%%%%%%%%%%%%%%
function p = phip(z)
    p = exp(-z.^2/2)/sqrt(2*pi);
end

%%%%%%%%%%%%%%
% Normal CDF %
%%%%%%%%%%%%%%
function p = phid(z)
    p = erfc( -z/sqrt(2) )/2;
end