function [StateVector, StateVectorsPKF, StateVectorsPPF, liky, updated_regimes, updated_sample, updated_mode, use_pkf_distribution, error_flag] = pkf_conditional_density(StateVector0, a, a1, P, P1, Py, alphahaty, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
    dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
% function [StateVector, StateVectorsPKF, StateVectorsPPF, liky, updated_regimes, updated_sample, updated_mode, use_pkf_distribution, error_flag] = pkf_conditional_density(StateVector0, a, a1, P, P1, Py, alphahaty, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
%     dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
%
% INPUTS
% - StateVector0            [struct]    Prior state vector at t-1 with mean, variance, and regime probabilities.
% - a                       [double]    PKF updated state estimate t-1.
%  - a1                     [double]    PKF one step ahead state prediction t-1:t
% - P                       [double]    PKF updated state covariance t-1.
% - P1                      [double]    PKF one step ahead state covariance t-1:t given t-1.
% - Py                      [double]    PKF updated state covariance t-1:t.
% - alphahaty               [double]    PKF smoothed state t-1:t|t.
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

yhat = StateVector0.Draws;
error_flag = 0;

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
number_of_successful_particles = 0;
updated_regimes = [];
updated_sample = [];
number_of_updated_regimes = 0;
updated_sample.success = false(number_of_particles,1);
likxmode = inf;
for k=1:number_of_particles
    % parfor k=1:number_of_particles
    a(:,1) = yhat(:,k);
    a1(:,1) = yhat(:,k);
    a1(:,2) = TT(:,:,2)*yhat(:,k);

    % conditional likelihood
    options_.occbin.filter.state_covariance=true;
    use_the_engine = true;
    if info==0
        % start with PKF regime to speed-up search!
        guess_regime = regimesy(1:2);
        options_.occbin.filter.guess_regime = true;
        [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regimesx, infox, ~, likxc, etahatx,alphahatx, Vxc] = ...
            occbin.kalman_update_algo_1(a,a1,P,P1,my_data_index,Z,v,my_data,H,QQQ,T0,R0,TT,RR,CC,guess_regime,M_, ...
            dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
        options_.occbin.filter.guess_regime = false;
        if infox == 0
            use_the_engine = false;
        end
    end
    if use_the_engine
        [ax, a1x, Px, P1x, vx, Tx, Rx, Cx, regimesx, infox, ~, likxc, etahatx,alphahatx, Vxc] = ...
            occbin.kalman_update_engine(a,a1,P,P1,t,my_data_index,Z,v,my_data,H,QQQ,T0,R0,TT,RR,CC,regimes0,base_regime,my_data_index{2},M_, ...
            dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
    end
    if infox ==0
        % handle conditional likelihood and updated regime
        [~, ~, tmp_str]=occbin.backward_map_regime(regimesx(1));
        if M_.occbin.constraint_nbr==1
            tmp_str = tmp_str(2,:);
        else
            tmp_str = [tmp_str(2,:)  tmp_str(4,:)];
        end
        if number_of_updated_regimes>0
            this_updated_regime = find(ismember(all_updated_regimes,{tmp_str}));
        end
        if number_of_updated_regimes==0 || isempty(this_updated_regime)
            number_of_updated_regimes = number_of_updated_regimes+1;
            this_updated_regime = number_of_updated_regimes;
            updated_regimes(number_of_updated_regimes).index = k;
            if isequal(regimesy(1),regimesx(1))
                updated_regimes(number_of_updated_regimes).is_pkf_regime = true;
            else
                updated_regimes(number_of_updated_regimes).is_pkf_regime = false;
            end
            updated_regimes(number_of_updated_regimes).regime = tmp_str;
            all_updated_regimes = {updated_regimes.regime};
            if ~options_.occbin.filter.particle.likelihood_only
                updated_regimes(number_of_updated_regimes).obsvar = ZZ*(Tx(:,:,1)*P0*Tx(:,:,1)'+Rx(:,:,1)*QQQ(:,:,2)*Rx(:,:,1)')*ZZ' + H(di,di);
                updated_regimes(number_of_updated_regimes).obsmean = ZZ*(Tx(:,:,1)*a0+Cx(:,1));
                updated_regimes(number_of_updated_regimes).ss.C = Cx(:,1);
                updated_regimes(number_of_updated_regimes).ss.R = Rx(:,:,1);
                updated_regimes(number_of_updated_regimes).ss.T = Tx(:,:,1);
                uF = updated_regimes(number_of_updated_regimes).obsvar;
                sig=sqrt(diag(uF));
                vv = Y(di,2) - updated_regimes(number_of_updated_regimes).obsmean;
                if options_.rescale_prediction_error_covariance
                    log_duF = log(det(uF./(sig*sig')))+2*sum(log(sig));
                    iuF = inv(uF./(sig*sig'))./(sig*sig');
                else
                    log_duF = log(det(uF));
                    iuF = inv(uF);
                end
                updated_regimes(number_of_updated_regimes).lik = ...
                    log_duF + transpose(vv)*iuF*vv + length(di)*log(2*pi);
            end
        else
            updated_regimes(this_updated_regime).index = [updated_regimes(this_updated_regime).index k];
        end

        number_of_successful_particles = number_of_successful_particles + 1;
        likxx(k) = likxc;
        if likxc<likxmode
            likxmode = likxc;
            updated_mode.likxc = likxc;
            updated_mode.a = ax; 
            updated_mode.a1 = a1x; 
            updated_mode.P = Px;
            updated_mode.P1 = P1x; 
            updated_mode.v = vx; 
            updated_mode.T = Tx; 
            updated_mode.R = Rx; 
            updated_mode.C = Cx; 
            updated_mode.regimes = regimesx; 
            if ~options_.occbin.filter.particle.likelihood_only
                updated_mode.lik = updated_regimes(this_updated_regime).lik;
            end
            updated_mode.etahat = etahatx;
            updated_mode.alphahat = alphahatx;
            updated_mode.V = Vxc;
        end

        updated_sample.epsilon(:,k)=etahatx;
        updated_sample.regimes{k} = regimesx;
        updated_sample.success(k)=true;
        updated_sample.y01(:,k) = alphahatx(:,1); %smoothed state 0|1!
        updated_sample.y11(:,k) = ax(:,1); %updated state!
        if M_.occbin.constraint_nbr==1
            updated_sample.regime_exit(k,1) = max(regimesx(1).regimestart);
            updated_sample.is_constrained(k,1) = logical(regimesx(1).regime(1));
            updated_sample.is_constrained_in_expectation(k,1) = any(regimesx(1).regime);
        else
            updated_sample.regime_exit(k,1) = max(regimesx(1).regimestart1);
            updated_sample.regime_exit(k,2) = max(regimesx(1).regimestart2);
            updated_sample.is_constrained(k,1) = logical(regimesx(1).regime1(1));
            updated_sample.is_constrained(k,2) = logical(regimesx(1).regime2(1));
            updated_sample.is_constrained_in_expectation(k,1) = any(regimesx(1).regime1);
            updated_sample.is_constrained_in_expectation(k,2) = any(regimesx(1).regime1);
        end
    end
end

if number_of_successful_particles==0
    error_flag = 340;
    StateVector.Draws=nan(size(yhat));
    return
end

lnw0 = -likxx(updated_sample.success)/2; % use conditional likelihood to weight prior draws!
% we need to normalize weights, even if they are already proper likelihood values!
% to avoid that exponential gives ZEROS all over the place!
[weights, updated_sample.ESS] = occbin.ppf.compute_weights(lnw0, 1, number_of_particles, number_of_successful_particles, updated_sample.success);

ParticleOptions.resampling.method.kitagawa=true;
indx = resample(0,weights',ParticleOptions);
StateVectorsPPF = updated_sample.y11(:,indx);
updated_sample.indx = indx;

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
    StateVectorsPKF = bsxfun(@plus,U(:,iss)*StateVectorVarianceSquareRoot*transpose(norminv(qmc_scrambled(state_variance_rank,number_of_particles,1))),StateVector.Mean);
else
    state_variance_rank = numel(StateVector.Mean );
    StateVectorsPKF = bsxfun(@plus,0*transpose(norminv(qmc_scrambled(state_variance_rank,number_of_particles,1))),StateVector.Mean);
end

%% PPF density
datalik = mean(exp(-likxx(updated_sample.success)/2));

liky = -log(datalik)*2;
updated_sample.likxx = likxx;

% test PPF updated state draws against PKF
StateVector.graph_info = StateVector0.graph_info;
ns = length(iss);
if ns==0
    disp_verbose(['ppf: time ' int2str(t) ' no state uncertainty!'],options_.debug);
    pkf_indicator=isequal(regimesy(1),base_regime);
else
    UP = U(:,iss);
    % WP = WP(:,isp);
    S = X(iss,iss); %UP(:,isp)'*(0.5*(P+P'))*WP(:,isp);
    iS = inv(S);

    chi2=nan(number_of_particles,1);
    for k=1:number_of_particles
        vv = UP'*(StateVectorsPPF(:,k)-StateVector.Mean);
        chi2(k,1) = transpose(vv)*iS*vv;
    end

    pkf_indicator =  max(chi2)<chi2inv(1-1/number_of_particles/2,ns);

    if not(options_.occbin.filter.particle.nograph)
        
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
        hold all, plot(x,schi2,'o' )
        title(['QQ plot PPF vs PKF updated state density, t = ' int2str(t)])
        dyn_saveas(hfig(6),[GraphDirectoryName, filesep, M_.fname,'_QQ_plot_state_density_t',int2str(t)],false,options_.graph_format);

    end
end

StateVector.proba=1;
if not( pkf_indicator  && any([updated_regimes.is_pkf_regime]) && (sum(ismember(indx,updated_regimes([updated_regimes.is_pkf_regime]).index))/number_of_particles)>=options_.occbin.filter.particle.use_pkf_updated_state_threshold )
    StateVector.Draws = StateVectorsPPF;
    PS  = cov(StateVector.Draws');
    [U,X] = svd(0.5*(PS+PS'));
    % P= U*X*U';
    iss = find(diag(X)>1.e-12);

    StateVector.Mean = mean(StateVector.Draws,2);
    StateVector.Variance = PS;
    if ~options_.occbin.filter.particle.draw_states_from_empirical_density
        if any(iss)
            StateVectorVarianceSquareRoot = chol(X(iss,iss))';
            % Get the rank of StateVectorVarianceSquareRoot
            state_variance_rank = size(StateVectorVarianceSquareRoot,2);
            StateVector.Draws = bsxfun(@plus,U(:,iss)*StateVectorVarianceSquareRoot*transpose(norminv(qmc_scrambled(state_variance_rank,number_of_particles,1))),StateVector.Mean);
        else
            state_variance_rank = numel(StateVector.Mean );
            StateVector.Draws = bsxfun(@plus,0*transpose(norminv(qmc_scrambled(state_variance_rank,number_of_particles,1))),StateVector.Mean);
        end
    end
    use_pkf_distribution=false;
else
    StateVector.Draws = StateVectorsPKF;
    use_pkf_distribution=true;
end
StateVector.Variance_rank = ns;
StateVector.use_pkf_distribution = use_pkf_distribution;



