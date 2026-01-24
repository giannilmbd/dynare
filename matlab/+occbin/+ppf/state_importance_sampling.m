function [StateVector, StateVectorsPKF, StateVectorsPPF, liky, updated_regimes, updated_sample, updated_mode, use_pkf_distribution, error_flag] = state_importance_sampling(StateVector0, a,  P, P1, Py, alphahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
    dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
% [StateVector, StateVectorsPKF, StateVectorsPPF, liky, updated_regimes, updated_sample, updated_mode, use_pkf_distribution, error_flag] = state_importance_sampling(StateVector0, a, P, P1, Py, alphahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
%     dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
%
% Performs importance sampling for state vector updates in the OccBin particle filter framework.
%
% INPUTS
% - StateVector0            [struct]    Prior state vector at t-1 with mean, variance, and regime probabilities.
% - a                       [double]    PKF updated state estimate t-1.
% - P                       [double]    PKF updated state covariance t-1.
% - P1                      [double]    PKF one step ahead state covariance t-1:t given t-1.
% - Py                      [double]    PKF updated state covariance t-1:t.
% - alphahaty               [double]    PKF smoothed state t-1:t|t.
% - V                       [double]    PKF smoothed state covariance t-1:t|t.
% - t                       [integer]   Current time period.
% - data_index              [cell]      1*2 cell of column vectors of indices.
% - Z                       [double]    Observation matrix.
% - v                       [double]    Observation innovations.
% - Y                       [double]    Observations t-1:t.
% - H                       [double]    Measurement error variance.
% - QQQ                     [double]    Shocks covariance matrix.
% - T0                      [double]    Base regime state transition matrix.
% - R0                      [double]    Base regime state shock matrix.
% - TT                      [double]    PKF Transition matrices t-1:t given t-1 info.
% - RR                      [double]    PKF Shock matrices t-1:t given t-1 info.
% - CC                      [double]    PKF Constant terms in transition t-1:t given t-1 info.
% - info                    [integer]   PKF error flag.
% - regimes0                [struct]    PKF regime info t:t+1 given t-1 info.
% - base_regime             [integer]   Base regime info.
% - regimesy                [struct]    PKF updated regime info t:t+2.
% - M_                      [struct]    Dynare's model structure.
% - dr                      [struct]    Decision rule structure.
% - endo_steady_state       [double]    Endogenous steady state.
% - exo_steady_state        [double]    Exogenous steady state.
% - exo_det_steady_state    [double]    Exogenous deterministic steady state.
% - options_                [struct]    Dynare options.
% - occbin_options          [struct]    OccBin specific options.
%
% OUTPUTS
%  - StateVector            [struct]    updated state vector structure
%  - StateVectorsPKF        [double]    state draws from PKF updated state density
%  - StateVectorsPPF        [double]    state draws from PPF (particle) updated state density
%  - liky                   [double]    minus 2*log-likelihood of observations
%  - updated_regimes        [struct]    occbin info about updated regimes
%  - updated_sample         [struct]    occbin info about updated state draw
%  - updated_mode           [struct]    occbin info about log-likelihood mode
%  - use_pkf_distribution   [logical]   indicator for using PKF distribution
%  - error_flag             [integer]   indicator of successful computation

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

% prior distribution in t-1

slice_override_iteration = max(2,options_.occbin.filter.particle.state_importance_sampling.slice_override_iteration);
slice_burnin = options_.occbin.filter.particle.state_importance_sampling.slice_burnin;
for r=1:size(StateVector0.Variance,3)
    [UP,XP] = svd(0.5*(StateVector0.Variance(:,:,r)+StateVector0.Variance(:,:,r)'));
    isp = find(diag(XP)>options_.kalman_tol);
    ns(r) = length(isp);
    if ns(r)==0
        disp_verbose(['ppf: time ' int2str(t) ' no state uncertainty!'],options_.debug);
    end
    UPR{r} = UP(:,isp);
    UPNULL{r} = UP(:,diag(XP)<=options_.kalman_tol);
    S = XP(isp,isp);
    log_dS(r) = log(det(S));
    iS{r} = inv(S);
end
PriorStateInfo.Mean = StateVector0.Mean;
PriorStateInfo.Variance = StateVector0.Variance;
PriorStateInfo.proba = StateVector0.proba;
PriorStateInfo.iS = iS;
PriorStateInfo.log_dS = log_dS;
PriorStateInfo.UPR = UPR;
PriorStateInfo.UPNULL = UPNULL;
ns=max(ns);
if r==1
    is_prior_mixture=false;
else
    is_prior_mixture=true;
end
error_flag = 0;
is_mixture=false; % initialise

% draws states in t-1 from PKF
number_of_particles = options_.occbin.filter.particle.number_of_particles;

% Initialization of the weights across particles.
niter = 0;
start_with_pkf_proposal=true;
qmc_base = norminv(qmc_scrambled(ns,number_of_particles,-1));
StateInfoPKF.Mean = alphahaty(:,1);
StateInfoPKF.Variance = V(:,:,1);
[~,X] = svd(0.5*(StateInfoPKF.Variance+StateInfoPKF.Variance'));
is = (diag(X)>options_.kalman_tol);
if length(is)>ns
    % dimension of posterior state distribution cannot exceed prior
    is = is(1:ns);
end
if not(isempty(is)) && StateVector0.use_pkf_distribution
    [yhat, lpost, lprior]  = occbin.ppf.draw_particles(qmc_base,StateInfoPKF,PriorStateInfo,number_of_particles,options_.kalman_tol);
    StateInfo = StateInfoPKF;
else
    if ns
        % initialize particles from the prior => niter=1 to drop initial step in the algo
        niter = 1;

        StateVectorINIT = StateVector0;

        [yhat, lpost, lprior]  = occbin.ppf.draw_particles(qmc_base,StateVectorINIT,PriorStateInfo,number_of_particles,options_.kalman_tol);
        StateInfo = StateVectorINIT;
        start_with_pkf_proposal=false;
    else
        % state variance is degenerate, all particles are the same
        yhat = repmat(alphahaty(:,1),[1 number_of_particles]);
        lpost = zeros(1,number_of_particles);
    end
end


% store input values
P0=P;
try
    chol(P0);
catch
    P0 = 0.5*(P0+P0');
end
a0=a;
P=P*0;
P1(:,:,1)=P;
QQ = RR(:,:,2)*QQQ(:,:,2)*transpose(RR(:,:,2));
di=data_index{2};
ZZ = Z(di,:);

P1(:,:,2) = QQ;
options_.occbin.filter.state_covariance=false;
my_data_index = data_index;
my_data_index{1} = [];
my_data = Y;
my_data(:,1) = nan;

number_of_particles = size(yhat,2);

is_qmc_sample = true;

ESS = 0;
ESS_max = 0;
smc_scale = 1;
posterior_kernel_at_the_mode = -inf;
use_modified_harmonic_mean = false;
while ESS<0.5*number_of_particles || do_resample
    niter = niter+1;
    number_of_successful_particles = 0;
    all_updated_regimes = {};
    updated_regimes = [];
    updated_sample = [];
    number_of_updated_regimes = 0;
    updated_sample.success = false(number_of_particles,1);
    likxx = nan(1,number_of_particles);
    likxmode = inf;
    this_base = qmc_base;
    for k=1:number_of_particles
        infox=1;
        while infox(1)
        [likxc, infox, likxmode, updated_mode, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes] = ...
            occbin.ppf.conditional_data_density(yhat(:,k), k, t, likxmode, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes, ...
            a0, P, P0, P1,my_data_index,Z,ZZ,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
            M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);

        if infox ==0
            % handle conditional likelihood and updated regime
            number_of_successful_particles = number_of_successful_particles + 1;
            likxx(k) = likxc;
        
        else
            this_base(k,:)=randn(1,ns);
            yhat(:,k) = occbin.ppf.draw_particles(this_base(k,:),StateInfo,PriorStateInfo,1,options_.kalman_tol);
        end
        end
    end

    if number_of_successful_particles==0
        error_flag = 340;
        StateVector.Draws=nan(size(yhat));
        return
    end
    sto.likxx = likxx;
    sto.lpost = lpost;
    sto.lprior = lprior;
    sto.smc_scale = smc_scale;
    sto.yhat = yhat;
    sto.updated_regimes = updated_regimes;
    sto.updated_sample = updated_sample;
    smc_scale_history(niter,1) = smc_scale;
    if is_qmc_sample
        do_resample = false;
        [logpogmax, imax] = max(-likxx(updated_sample.success)/2-lprior(updated_sample.success)/2);
        if logpogmax>posterior_kernel_at_the_mode
            posterior_kernel_at_the_mode = logpogmax;
            ytmp = yhat(:,updated_sample.success);
            ymode = ytmp(:,imax);
        end
        
        lnw0 = -likxx(updated_sample.success)/2-lprior(updated_sample.success)/2+lpost(updated_sample.success)/2; % use posterior kernel to weight t-1 draws based on PKF!
        [weights1, ESS1] = occbin.ppf.compute_weights(lnw0, 1, number_of_particles, number_of_successful_particles, updated_sample.success);
        ess1_history(niter,1) = ESS1;
        if ESS1>ESS_max
            ESS_max = ESS1;
            stomax = sto;
            stomax.lnw0 = lnw0;
            stomax.niter = niter;
            stomax.weights = weights1;
            stomax.ESS = ESS1;
            stomax.StateInfo = StateInfo;
        end
        if ESS1>0.5*number_of_particles
            % the proposal worked
            smc_scale=1;
        end
        [weights, ESS_init] = occbin.ppf.compute_weights(lnw0, smc_scale, number_of_particles, number_of_successful_particles, updated_sample.success);
        ess_init_history(niter,1) = ESS_init;
        ESS=0;
        while ESS<0.5*number_of_particles
            [weights, ESS] = occbin.ppf.compute_weights(lnw0, smc_scale, number_of_particles, number_of_successful_particles, updated_sample.success);

            if ESS<0.5*number_of_particles || (niter==1 && length(updated_regimes)>1)
                % if PKF posterior draws do not deliver a unique regime, then
                % restart from the prior
                do_resample = true;
                if niter==1
                    if  length(updated_regimes)>1
                        ESS=0;
                    end
                    break
                elseif smc_scale==1.e-6
                    break
                else
                    smc_scale = max(1.e-6,smc_scale*0.98);
                end
            end
        end

        if not(do_resample) && smc_scale<1
            if ESS>0.5*number_of_particles && smc_scale<1
                smc_scale = min(1,smc_scale*1.1);
                [weights, ESS] = occbin.ppf.compute_weights(lnw0, smc_scale, number_of_particles, number_of_successful_particles, updated_sample.success);
            end
            while ESS<0.5*number_of_particles
                smc_scale = smc_scale*0.98;
                [weights, ESS] = occbin.ppf.compute_weights(lnw0, smc_scale, number_of_particles, number_of_successful_particles, updated_sample.success);
            end
            do_resample=true;
        end

        ParticleOptions.resampling.method.kitagawa=true;
        indx = resample(0,weights',ParticleOptions);
    else
        stomcmc=sto;
        stomcmc.ljoint = ljoint;
        stomcmc.mcmc_base = mcmc_base;
        indx = find(updated_sample.success);
        ESS = number_of_particles;
        % continue importance sampling if modified harmonic mean had
        % problems
        do_resample=~use_modified_harmonic_mean;
    end
    StateVectorsPPF = updated_sample.y11(:,indx);
    % map updated regimes in t
    indunconstrained = indx(~ismember(indx,find(sum(updated_sample.is_constrained,2))));
    nunconstrained = length(indunconstrained);
    is_mixture = false;
    if M_.occbin.constraint_nbr==1
        indconstrained = indx(ismember(indx,find(updated_sample.is_constrained)));
        indcase = {indunconstrained, indconstrained};
        ncases =[nunconstrained; length(indconstrained)];
    else
        indconstr1 = indx(ismember(indx,find(updated_sample.is_constrained,1)+~updated_sample.is_constrained,2));
        indconstr2 = indx(ismember(indx,find(updated_sample.is_constrained,2)+~updated_sample.is_constrained,1));
        indconstr12 = indx(ismember(indx,find(updated_sample.is_constrained,1).*updated_sample.is_constrained,2));
        indcase = {indunconstrained, indconstr1, indconstr2, indconstr12};
        nconstr1 = length(indconstr1);
        nconstr2 = length(indconstr2);
        nconstr12 = length(indconstr12);
        ncases = [nunconstrained; nconstr1; nconstr2; nconstr12];
    end
    regimevec = zeros(length(ncases),1);
    for r=1:length(ncases)
        regimevec(r,1) = (ncases(r)<min(number_of_successful_particles-ns*(ns+1)/2,0.9*number_of_successful_particles) && ncases(r)>max(ns*(ns+1)/2,0.1*number_of_successful_particles));
    end
    regimevec = logical(regimevec);
    if any(regimevec) && isequal(regimevec, ncases>0) && sum(regimevec)>1
        % we split only if there is enough information to define
        % mixtures
        is_mixture = true;
        rcases = sum(regimevec);
        indcase = indcase(regimevec);
        ncases = ncases(regimevec);
        for r=1:rcases
            StateVectorsPPFcases{r} = updated_sample.y11(:,indcase{r});
        end
    end

    if do_resample
        StateInfo.proba = 1;
        if niter==1
            % start from prior
            StateVectorsPPFMean = StateVector0.Mean;
            StateVectorsPPFVariance = StateVector0.Variance;
        else
            if ~is_qmc_sample
                is_qmc_sample=true;
            end
            if is_mixture
                for r=1:rcases
                    StateVectorsPPFVariance(:,:,r)  = cov(yhat(:,indcase{r})');
                    StateVectorsPPFMean(:,r) = mean(yhat(:,indcase{r}),2);
                end
                StateInfo.proba = ncases./number_of_successful_particles;
                StateInfo.proba = StateInfo.proba./sum(StateInfo.proba); % normalise
            else
                StateVectorsPPFVariance  = cov(yhat(:,indx)');
                StateVectorsPPFMean = mean(yhat(:,indx),2);
            end

        end
        if mod(niter,slice_override_iteration) || is_prior_mixture 
            % slice sampling is not ready for prior mixture state
            % distribution
            StateInfo.Mean=StateVectorsPPFMean;
            StateInfo.Variance=StateVectorsPPFVariance;
            [yhat, lpost, lprior]  = occbin.ppf.draw_particles(qmc_base,StateInfo,PriorStateInfo,number_of_particles,options_.kalman_tol);
        else
            % %             % slice sampling!
            is_qmc_sample=false;
        end
        if ~is_qmc_sample
            [~, init_draw]  = occbin.ppf.state_priordens(ymode,PriorStateInfo,1,options_.kalman_tol);
            % check outcome
            [xsim, fsim] = occbin.ppf.mcmc_draws(init_draw',number_of_particles+slice_burnin,StateVector0,PriorStateInfo, ...
                a0, P, P0, P1,my_data_index,Z,ZZ,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
                M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
            mcmc_base = xsim(:,(1:number_of_particles)+slice_burnin)';
            ljoint = fsim((1:number_of_particles)+slice_burnin)';

            [ModifiedHarmonicMean, crit_flag] = occbin.ppf.marginal_density(mcmc_base, ljoint, 2/number_of_particles);            
            [yhat, lpost, lprior]  = occbin.ppf.draw_particles(mcmc_base,StateVector0,PriorStateInfo,number_of_particles,options_.kalman_tol);
            smc_scale = 1;
            use_modified_harmonic_mean = crit_flag==0;
            if crit_flag
                slice_override_iteration=inf;
            end
        end
    end
end

updated_sample.indx = indx;
updated_sample.ESS = ESS;

% updated state mean and covariance
StateVector.Variance = Py(:,:,1);
[U,X] = svd(0.5*(Py(:,:,1)+Py(:,:,1)'));
% P= U*X*U';
iss = find(diag(X)>1.e-12);

StateVector.Mean = alphahaty(:,2);
if any(iss)
    StateVectorVarianceSquareRoot = chol(X(iss,iss))';
    % Get the rank of StateVectorVarianceSquareRoot
    state_variance_rank = size(StateVectorVarianceSquareRoot,2);
    StateVectorsPKF = bsxfun(@plus,U(:,iss)*StateVectorVarianceSquareRoot*transpose(norminv(qmc_scrambled(state_variance_rank,number_of_particles,-1))),StateVector.Mean);
else
    state_variance_rank = numel(StateVector.Mean );
    StateVectorsPKF = bsxfun(@plus,0*transpose(norminv(qmc_scrambled(state_variance_rank,number_of_particles,-1))),StateVector.Mean);
end
StateVector.graph_info = StateVector0.graph_info;

%% PPF density
if use_modified_harmonic_mean
    liky = -ModifiedHarmonicMean*2;
else
    datalik = mean(exp(-likxx(updated_sample.success)/2-lprior(updated_sample.success)/2+lpost(updated_sample.success)/2));
    liky = -log(datalik)*2;
end
updated_sample.likxx = likxx;

ns = length(iss);
if ns==0
    disp_verbose(['ppf: time ' int2str(t) ' no state uncertainty!'],options_.debug);
    [~,X] = svd(cov(StateVectorsPPF'));
    nsppf = sum(diag(X)>1.e-12);
    pkf_indicator = nsppf==0;

else
    UP = U(:,iss);
    S = X(iss,iss); %UP(:,isp)'*(0.5*(P+P'))*WP(:,isp);
    iS = inv(S);

    chi2=nan(number_of_particles,1);
    for k=1:number_of_particles
        vv = UP'*(StateVectorsPPF(:,k)-StateVector.Mean);
        chi2(k,1) = transpose(vv)*iS*vv;
    end
    pkf_indicator =  max(chi2)<chi2inv(1-1/number_of_particles/2,ns);
end

StateVector.proba = 1;
if not( (isequal(regimesy(1),base_regime) ||  (abs(ESS-number_of_particles)<(1/number_of_particles)^2 && niter==1 && start_with_pkf_proposal) || pkf_indicator) ...
        && any([updated_regimes.is_pkf_regime]) && (sum(ismember(indx,updated_regimes([updated_regimes.is_pkf_regime]).index))/number_of_particles)>=options_.occbin.filter.particle.use_pkf_updated_state_threshold )
    StateVector.Draws = StateVectorsPPF;
    if is_mixture
        for r=1:length(StateVectorsPPFcases)
            StateVector.Variance(:,:,r)  = cov(StateVectorsPPFcases{r}');
            StateVector.Mean(:,r)  = mean(StateVectorsPPFcases{r} ,2);
            [~,X] = svd(StateVector.Variance(:,:,r));
            nsppf(r) = sum(diag(X)>options_.kalman_tol );
        end
        nsppf = max(nsppf);
        StateVector.proba = ncases./number_of_successful_particles;
        StateVector.proba = StateVector.proba./sum(StateVector.proba); % normalise
    else
        StateVector.Variance  = cov(StateVector.Draws');
        StateVector.Mean = mean(StateVector.Draws,2);
        [~,X] = svd(StateVector.Variance);
        nsppf = sum(diag(X)>options_.kalman_tol );
    end
    if ~options_.occbin.filter.particle.draw_states_from_empirical_density
        if nsppf>size(qmc_base,2)
            qmc_base = norminv(qmc_scrambled(nsppf,number_of_particles,-1));
        end
        StateVector.Draws  = occbin.ppf.draw_particles(qmc_base,StateVector,[],number_of_particles,options_.kalman_tol);
    end
    use_pkf_distribution=false;
else
    StateVector.Draws = StateVectorsPKF;
    use_pkf_distribution=true;
end
StateVector.Variance_rank = ns;
StateVector.use_modified_harmonic_mean = use_modified_harmonic_mean;
StateVector.use_pkf_distribution = use_pkf_distribution;

if not(options_.occbin.filter.particle.diagnostics.nograph) && ns && ~use_pkf_distribution

    GraphDirectoryName = CheckPath('occbin_ppf_graphs',M_.dname);
    schi2 = sort(chi2);
    x = zeros(size(schi2));
    for j=1:number_of_particles
        p=0.5*1/number_of_particles+(j-1)/number_of_particles;
        x(j,1) = chi2inv(p,ns);
    end

    % qqplot with chi square distribution
    if isnan(StateVector0.graph_info.hfig(6))
        hfig(6) = figure;
        StateVector.graph_info.hfig(6) = hfig(6);
    else
        hfig(6) = StateVector.graph_info.hfig(6);
        figure(hfig(6))
    end
    clf(hfig(6),'reset')
    set(hfig(6),'name',['QQ plot PPF vs PKF updated state density, t = ' int2str(t)])
    plot(x,x )
    hold on, plot(x,schi2,'o' )
    title(['QQ plot PPF vs PKF updated state density, t = ' int2str(t)])
    dyn_saveas(hfig(6),[GraphDirectoryName, filesep, M_.fname,'_QQ_plot_state_density_t',int2str(t)],false,options_.graph_format);

end

