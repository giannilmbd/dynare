function [yhat, lpost, lprior]  = draw_particles(qmc_base,StateInfo,PriorStateInfo,number_of_particles,kalman_tol)
% [yhat, lpost, lprior]  = draw_particles(qmc_base,StateInfo,number_of_particles,kalman_tol)
%
% Draw particles in t-1|t from a mixture of normal distributions [up to as many mixtures as the number of OccBin regimes]
%
% INPUTS
%  - qmc_base               [double]    base normal draws (independent standard normal variates)
%  - StateInfo              [struct]    info of particles sampling distribution(s) t-1|t
%  - PriorStateInfo         [struct]    info of particles distribution t-1|t-1
%  - number_of_particles    [integer]   number of particles to draw
%  - kalman_tol             [double]    tolerance for Kalman filter computations
%
% OUTPUTS
%  - yhat                   [double]    vector of particles
%  - lpost                  [double]    minus log proposal density of particles
%  - lprior                 [double]    minus log prior density of particles
%
% This function is called by: sequential_importance_particle_filter
% This function calls: residual_resampling, traditional_resampling

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

nmixtures = size(StateInfo.Variance,3);
is_mixture = nmixtures>1;
yhat=[];
log_dS0_ = nan(1,nmixtures);
lpost=zeros(1,number_of_particles);
lpost_tmp=zeros(nmixtures,number_of_particles);
for r=1:nmixtures

    % now do svd and sample from gaussian approx thereof ...
    [U,X] = svd(0.5*(StateInfo.Variance(:,:,r)+StateInfo.Variance(:,:,r)'));
    % P= U*X*U'; % symmetric matrix!
    is = (diag(X)>kalman_tol);
    isnull = (diag(X)<=kalman_tol); % null space
    StateVectorVarianceSquareRoot = chol(X(is,is))';%reduced_rank_cholesky(ReducedForm.StateVectorVariance)';

    % Get the rank of StateVectorVarianceSquareRoot
    UR{r} = U(:,is);
    UNULL{r} = U(:,isnull);
    log_dS0_(r) = log(det(X(is,is)));
    iS0_r{r} = 1./diag(X(is,is));
    if sum(is)>size(qmc_base,2)
        qmc_base = norminv(qmc_scrambled(sum(is),number_of_particles,-1));
    end
    % % Factorize the covariance matrix of the structural innovations
    % Q_lower_triangular_cholesky = chol(Q)';
    if is_mixture
        if r==1
            ndraws = number_of_particles-sum(round(StateInfo.proba(2:end)*number_of_particles));
        else
            ndraws = round(StateInfo.proba(r)*number_of_particles);
        end
        this_qmc = qmc_base(1:ndraws,1:sum(is));
    else
        ndraws = number_of_particles;
        this_qmc = qmc_base(:,1:sum(is));
    end
    % Initialization of the weights across particles.
    if any(is)
        yhat_tmp = bsxfun(@plus,U(:,is)*StateVectorVarianceSquareRoot*transpose(this_qmc),StateInfo.Mean(:,r));
    else
        yhat_tmp = repmat(StateInfo.Mean(:,r),[1 ndraws]);
    end
    yhat = [yhat yhat_tmp];
end
for r=1:size(StateInfo.Variance,3)
    p0=1;
    if is_mixture
        p0 = StateInfo.proba(r);
    end
    if isempty(UR{r})
        lpost_tmp(r,:) = 0 - log(p0);
    else
        vv0 = UR{r}'*bsxfun(@plus,yhat,-StateInfo.Mean(:,r));
        lpost_tmp(r,:) = log_dS0_(r) + iS0_r{r}'*vv0.^2 + length(iS0_r{r})*log(2*pi) - log(p0);
    end
    if ~isempty(UNULL{r})
        vvNULL = UNULL{r}'*bsxfun(@plus,yhat,-StateInfo.Mean(:,r));
        if any(any(vvNULL.^2>3*kalman_tol))
            if size(vvNULL,1)>1
                lpost_tmp(r,any(vvNULL.^2>3*kalman_tol)) = max(lpost_tmp(r,:))+100;
            else
                lpost_tmp(r,vvNULL.^2>3*kalman_tol) = max(lpost_tmp(r,:))+100;
            end
        end
    end
end
if is_mixture
    for p=1:number_of_particles
        likvec = lpost_tmp(:,p);
        lpost(p) = -2*log(sum(exp(-(likvec-min(likvec))./2)))+min(likvec);
    end
else
    lpost = lpost_tmp;
end

lprior=[];
if nargout>2 && not(isempty(PriorStateInfo))
    lprior=zeros(1,number_of_particles);
    is_prior_mixture = size(PriorStateInfo.Variance,3)>1;
    for k=1:number_of_particles
        lprior_tmp=NaN(1,size(PriorStateInfo.Variance,3));
        for r=1:size(PriorStateInfo.Variance,3)
            p0=1;
            if is_prior_mixture
                p0 = PriorStateInfo.proba(r);
            end
            if isempty(PriorStateInfo.UPR{r})
                lprior_tmp(r) = 0 - log(p0);
            else
                vv = PriorStateInfo.UPR{r}'*(yhat(:,k)-PriorStateInfo.Mean(:,r));
                lprior_tmp(r) = PriorStateInfo.log_dS(r) + transpose(vv)*PriorStateInfo.iS{r}*vv + length(PriorStateInfo.iS{r})*log(2*pi) - log(p0);
            end
            vvNULL = PriorStateInfo.UPNULL{r}'*(yhat(:,k)-PriorStateInfo.Mean(:,r));
            if any(vvNULL.^2>3*kalman_tol)
                lprior_tmp(r) = lprior_tmp(r)+100;
            end
        end
        if is_prior_mixture
            lprior(k) = -2*log(sum(exp(-(lprior_tmp-min(lprior_tmp))./2)))+min(lprior_tmp);
        else
            lprior(k) = lprior_tmp;
        end
    end
end

