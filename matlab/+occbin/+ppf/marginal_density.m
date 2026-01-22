function [ModifiedHarmonicMean, crit_flag] = marginal_density(x2, logpo2, tolerance)  
% [ModifiedHarmonicMean, crit_flag] = marginal_density(x2, logpo2, tolerance)  
% Computes the marginal likelihood using the modified harmonic mean estimator.
%
% INPUTS
% - x2                      [double]    Parameter draws (replicas by row).
% - logpo2                  [double]    Log-posterior density values for each draw.
% - tolerance               [double]    Convergence tolerance for the iterative procedure.
%
% OUTPUTS
% - ModifiedHarmonicMean    [double]    Estimated marginal likelihood.
% - crit_flag               [integer]   Flag indicating successful computation (0 = success).
%
% REMARKS
% The function uses an iterative procedure that increases the variance of the Normal
% approximation to the posterior until the coverage probabilities have stabilized.

% Copyright © 2025-2026 Dynare Team
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

posterior_kernel_at_the_mode = max(logpo2);
posterior_mean = mean(x2)';
posterior_covariance = cov(x2);
MU = transpose(posterior_mean);
SIGMA = posterior_covariance;
lpost_mode = posterior_kernel_at_the_mode;
xparam1 = posterior_mean;
hh = inv(SIGMA);
npar = length(xparam1);
crit_flag=0;
% save the posterior mean and the inverse of the covariance matrix
% (useful if the user wants to perform some computations using
% the posterior mean instead of the posterior mode ==> ).

try 
    % use this robust option to avoid inf/nan
    logdetSIGMA = 2*sum(log(diag(chol(SIGMA)))); 
catch
    % in case SIGMA is not positive definite
    logdetSIGMA = nan;
    fprintf('marginal density: the covariance of MCMC draws is not positive definite. You may have too few MCMC draws.');
end
invSIGMA = hh;
marginal = zeros(9,2);
linee = 0;
check_coverage = 1;
increase = 1;
niter=0;
while check_coverage
    niter=niter+1;
    for p = 0.1:0.1:0.9
        critval = chi2inv(p,npar);
        ifil = 1;
        tmp = 0;
        EndOfFile = size(x2,1);
        for i = ifil:EndOfFile
            deviation  = ((x2(i,:)-MU)*invSIGMA*(x2(i,:)-MU)')/increase;
            if deviation <= critval
                lftheta = -log(p)-(npar*log(2*pi)+(npar*log(increase)+logdetSIGMA)+deviation)/2;
                tmp = tmp + exp(lftheta - logpo2(i) + lpost_mode);
            end
        end
        linee = linee + 1;
        warning_old_state = warning;
        warning off;
        marginal(linee,:) = [p, lpost_mode-log(tmp/length(logpo2))];
        warning(warning_old_state);
    end
    tol_check(niter)=abs((marginal(9,2)-marginal(1,2))/marginal(9,2));
    ModifiedHarmonicMeanVec(niter) = mean(marginal(:,2));
    if tol_check(niter) > tolerance || isinf(marginal(1,2))
        if increase == 1
            increase = 1.2*increase;
            linee    = 0;
        else
            increase = 1.2*increase;
            linee    = 0;
            if increase > 20
                check_coverage = 0;
                crit_flag=1;
                clear invSIGMA detSIGMA increase;
            end
        end
    else
        check_coverage = 0;
        clear invSIGMA detSIGMA increase;
    end
end

[~, im ]= min(tol_check);
ModifiedHarmonicMean = ModifiedHarmonicMeanVec(im);