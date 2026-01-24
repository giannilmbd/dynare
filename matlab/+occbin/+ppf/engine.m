function [StateVectors,liky, updated_regimes, updated_sample, updated_mode, use_pkf_distribution, error_flag] = engine(StateVector0, likpkf, a, a1, P, P1, a1y, Py, P1y, alphahaty, etahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,isqvec,M_, ...
    dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
% [StateVectors,liky, updated_regimes, updated_sample, updated_mode, use_pkf_distribution, error_flag] = engine(StateVector0, likpkf, a, a1, P, P1, a1y, Py, P1y, alphahaty, etahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,isqvec,M_, ...
%     dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
%
% Main engine for the OccBin particle filter framework, combining PKF and PPF methods.
%
% INPUTS
% - StateVector0            [struct]    Prior state vector (can be empty for initialization).
% - likpkf                  [double]    PKF likelihood value.
% - a                       [double]    PKF updated state estimate t-1.
% - a1                      [double]    PKF one step ahead state prediction t-1:t (given t-1 information).
% - P                       [double]    PKF updated state covariance t-1.
% - P1                      [double]    PKF one step ahead state covariance t-1:t (given t-1 information).
% - a1y                     [double]    PKF one step ahead state prediction t-1:t (given t information).
% - Py                      [double]    PKF updated state covariance t-1:t.
% - P1y                     [double]    PKF one step ahead state covariance t-1:t (given t information).
% - alphahaty               [double]    PKF smoothed states t-1:t|t.
% - etahaty                 [double]    PKF smoothed shocks t-1:t|t.
% - V                       [double]    PKF smoothed states covariance t-1:t|t.
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
% - isqvec                  [double]    Index vector for shocks.
% - M_                      [struct]    Dynare's model structure.
% - dr                      [struct]    Decision rule structure.
% - endo_steady_state       [double]    Endogenous steady state.
% - exo_steady_state        [double]    Exogenous steady state.
% - exo_det_steady_state    [double]    Exogenous deterministic steady state.
% - options_                [struct]    Dynare options.
% - occbin_options          [struct]    OccBin specific options.
%
% OUTPUTS
% - StateVectors            [cell]      Updated state vectors from both methods.
% - liky                    [double]    Log-likelihood contribution.
% - updated_regimes         [struct]    Updated regime indicators.
% - updated_sample          [struct]    Updated particle sample.
% - updated_mode            [double]    Updated mode estimate.
% - use_pkf_distribution    [logical]   Flag indicating PKF distribution usage.
% - error_flag              [logical]   Error flag (true if computation failed).

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

info = 0;
number_of_particles = options_.occbin.filter.particle.number_of_particles;
if isempty(StateVector0)
    [U,X] = svd(0.5*(P+P'));
    % P= U*X*U'; % symmetric matrix!
    is = find(diag(X)>1.e-12);
    StateVectorVarianceSquareRoot = chol(X(is,is))';

    % Get the rank of StateVectorVarianceSquareRoot
    state_variance_rank = size(StateVectorVarianceSquareRoot,2);

    % % Factorize the covariance matrix of the structural innovations
    % Q_lower_triangular_cholesky = chol(Q)';

    % Initialization of the weights across particles.
    if not(isempty(is))
        yhat = U(:,is)*StateVectorVarianceSquareRoot*transpose(norminv(qmc_scrambled(state_variance_rank,number_of_particles,1)))+repmat(a,[1 number_of_particles]);
    else
        yhat = repmat(a,[1 number_of_particles]);
    end
    StateVector0.Draws = yhat;
    StateVector0.Mean = a;
    StateVector0.Variance = P;
    StateVector0.Variance_rank = length(is);
    StateVector0.use_pkf_distribution = true;
    StateVector0.graph_info.hfig = nan(8,1);
    StateVector0.proba = 1;
    StateVector0.stop_particles = false;
end

if options_.occbin.filter.particle.diagnostics.status || options_.occbin.filter.particle.empirical_data_density.status || options_.occbin.filter.particle.ensemble_kalman_filter
    % random shocks
    [US,XS] = svd(QQQ(:,:,2));
    % P= U*X*U';
    ishock = find(sqrt(abs(diag(XS)))>1.e-10);
    ShockVarianceSquareRoot = chol(XS(ishock,ishock))';
    % Get the rank of ShockVarianceSquareRoot
    shock_variance_rank = size(ShockVarianceSquareRoot,2);
    ShockVectors = bsxfun(@plus,US(:,ishock)*ShockVarianceSquareRoot*transpose(norminv(qmc_scrambled(shock_variance_rank,number_of_particles,1))),zeros(length(US),1));
end

if  ~options_.occbin.filter.particle.diagnostics.status || (~isempty(options_.occbin.filter.particle.diagnostics.graph_periods) && ~ismember(t,options_.occbin.filter.particle.diagnostics.graph_periods))
    options_.occbin.filter.particle.diagnostics.nograph=true;
end

if options_.occbin.filter.particle.ensemble_kalman_filter
    [liky, StateVectors, sim_regimes, sim_sample, error_flag] = local_EnKF_iteration(StateVector0,ShockVectors,data_index(2),Z,Y(:,2),H,M_,dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
    use_pkf_distribution = false;
    updated_regimes = sim_regimes; % we use 1|0 for 1|1
    updated_sample = sim_sample; % we use 1|1 for 1|1
    updated_mode = [];
else
    [LastSeeds.Unifor, LastSeeds.Normal] = get_dynare_random_generator_state();
    if options_.occbin.filter.particle.state_importance_sampling.pkf_init && (StateVector0.use_pkf_distribution || ~options_.occbin.filter.particle.draw_states_from_empirical_density)
        [StateVectors, StateVectorsPKF, StateVectorsPPF, liky, updated_regimes, updated_sample, updated_mode, use_pkf_distribution, error_flag] = occbin.ppf.state_importance_sampling(StateVector0, a, P, P1, Py, alphahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
            dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
        if abs(liky-likpkf)>options_.occbin.filter.particle.state_importance_sampling.logpost_crit_threshold
            liky1=liky;
            use_new_state_update = false;
            if ~StateVectors.use_modified_harmonic_mean
                orig_slice_override_iteration = options_.occbin.filter.particle.state_importance_sampling.slice_override_iteration; 
                options_.occbin.filter.particle.state_importance_sampling.slice_override_iteration = 1; 
                set_dynare_random_generator_state(LastSeeds.Unifor, LastSeeds.Normal);
                [StateVectors_tmp, StateVectorsPKF_tmp, StateVectorsPPF_tmp, liky_tmp, updated_regimes_tmp, updated_sample_tmp, updated_mode_tmp, use_pkf_distribution_tmp, error_flag_tmp] = occbin.ppf.state_importance_sampling(StateVector0, a, P, P1, Py, alphahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
                    dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
                options_.occbin.filter.particle.state_importance_sampling.slice_override_iteration = orig_slice_override_iteration; 
                if abs(liky_tmp-likpkf)<options_.occbin.filter.particle.state_importance_sampling.logpost_crit_threshold || isinf(liky)
                    use_new_state_update = true;
                end
            end
            if ~use_new_state_update
                set_dynare_random_generator_state(LastSeeds.Unifor, LastSeeds.Normal);
                [~, ~, ~, liky1, ~, ~, ~, use_pkf_distribution, error_flag] = occbin.pkf_conditional_density(StateVector0, a, a1, P, P1, Py, alphahaty, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
                    dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
            end
            if abs(liky1-likpkf)<options_.occbin.filter.particle.state_importance_sampling.logpost_crit_threshold && ~use_new_state_update
                kalman_tol0 = options_.kalman_tol;
                options_.kalman_tol = kalman_tol0*10;
                set_dynare_random_generator_state(LastSeeds.Unifor, LastSeeds.Normal);
                [StateVectors_tmp, StateVectorsPKF_tmp, StateVectorsPPF_tmp, liky_tmp, updated_regimes_tmp, updated_sample_tmp, updated_mode_tmp, use_pkf_distribution_tmp, error_flag_tmp] = occbin.ppf.state_importance_sampling(StateVector0, a, P, P1, Py, alphahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
                    dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
                if abs(liky_tmp-likpkf)<options_.occbin.filter.particle.state_importance_sampling.logpost_crit_threshold
                    use_new_state_update = true;
                else
                    options_.kalman_tol = kalman_tol0*100;
                    set_dynare_random_generator_state(LastSeeds.Unifor, LastSeeds.Normal);
                    [StateVectors_tmp, StateVectorsPKF_tmp, StateVectorsPPF_tmp, liky_tmp, updated_regimes_tmp, updated_sample_tmp, updated_mode_tmp, use_pkf_distribution_tmp, error_flag_tmp] = occbin.ppf.state_importance_sampling(StateVector0, a, P, P1, Py, alphahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
                        dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
                    if abs(liky_tmp-likpkf)<options_.occbin.filter.particle.state_importance_sampling.logpost_crit_threshold
                        use_new_state_update = true;
                    else
                        options_.kalman_tol = kalman_tol0*0.1;
                        set_dynare_random_generator_state(LastSeeds.Unifor, LastSeeds.Normal);
                        [StateVectors_tmp, StateVectorsPKF_tmp, StateVectorsPPF_tmp, liky_tmp, updated_regimes_tmp, updated_sample_tmp, updated_mode_tmp, use_pkf_distribution_tmp, error_flag_tmp] = occbin.ppf.state_importance_sampling(StateVector0, a, P, P1, Py, alphahaty, V, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
                            dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
                        if abs(liky_tmp-likpkf)<options_.occbin.filter.particle.state_importance_sampling.logpost_crit_threshold
                            use_new_state_update = true;
                        end
                    end
                end
                options_.kalman_tol = kalman_tol0;
            end
            if use_new_state_update
                StateVectors = StateVectors_tmp; 
                StateVectorsPKF = StateVectorsPKF_tmp; 
                StateVectorsPPF = StateVectorsPPF_tmp; 
                liky = liky_tmp; 
                updated_regimes = updated_regimes_tmp; 
                updated_sample = updated_sample_tmp; 
                updated_mode = updated_mode_tmp; 
                use_pkf_distribution = use_pkf_distribution_tmp; 
                error_flag = error_flag_tmp;
            end
        end
    else
        [StateVectors, StateVectorsPKF, StateVectorsPPF, liky, updated_regimes, updated_sample, updated_mode, use_pkf_distribution, error_flag] = occbin.pkf_conditional_density(StateVector0, a, a1, P, P1, Py, alphahaty, t, data_index,Z,v,Y,H,QQQ,T0,R0,TT,RR,CC,info,regimes0,base_regime,regimesy,M_, ...
            dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
    end
end
graph_info = StateVectors.graph_info;
StateVectors.stop_particles = false;

if options_.occbin.filter.particle.diagnostics.status || options_.occbin.filter.particle.empirical_data_density.status
    di=data_index{2};
    ZZ = Z(di,:);
    opts_simul = occbin_options.opts_simul;
    if opts_simul.restrict_state_space
        my_order_var = dr.order_var(dr.restrict_var_list);
    else
        my_order_var = dr.order_var;
    end
    number_of_shocks_per_particle = 1;
    [llik, simulated_regimes, simulated_sample] = ...
        occbin.ppf.simulated_density(number_of_particles, number_of_shocks_per_particle,StateVector0.Draws,ShockVectors, ...
        di, H, my_order_var, QQQ, Y, ZZ, base_regime, regimesy, ...
        M_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, options_, opts_simul);

    if  ~options_.debug && StateVector0.use_pkf_distribution && StateVectors.use_pkf_distribution
        % do not produce plots when state updates that are identical by construction
        % even if nograph==false, unless debug
        options_.occbin.filter.particle.diagnostics.nograph=true;
    end
    % we enter in any case to compute density and density data
    [density, density_data, graph_info] = occbin.ppf.graphs(t, updated_sample.indx, ShockVectors, StateVectorsPKF, StateVectorsPPF, updated_sample.likxx, regimes0, regimesy, ...
        simulated_regimes, simulated_sample, updated_regimes, updated_sample, a1, P1, a1y, P1y, alphahaty, etahaty, ...
        di, H, QQQ, Y, ZZ, dr, M_, options_, graph_info);
    StateVectors.graph_info = graph_info;

    % collect results
    success = simulated_sample.success;
    y100 = simulated_sample.y10;
    y10 = y100(:,success);
    for jk=1:size(y10,1)
        [output.filtered.variables.mean(jk,1), ...
            output.filtered.variables.median(jk,1), ...
            output.filtered.variables.var(jk,1), ...
            output.filtered.variables.hpd_interval(jk,:), ...
            output.filtered.variables.deciles(jk,:), ...
            output.filtered.variables.density(jk,:,:)] = posterior_moments(y10(jk,:)', options_.mh_conf_sig);
    end
    output.filtered.variables.particles = y100;
    output.filtered.regimes.sample = simulated_sample.regime;
    output.filtered.regimes.is_constrained = simulated_sample.is_constrained;
    output.filtered.regimes.is_constrained_in_expectation = simulated_sample.is_constrained_in_expectation;
    [aa,bb]=histcounts(simulated_sample.exit(success),'normalization','pdf','binmethod','integers');
    output.filtered.regimes.exit(bb(1:end-1)+0.5,1) = aa';
    isc = simulated_sample.is_constrained;
    isce = simulated_sample.is_constrained_in_expectation;
    rr=simulated_sample.exit;
    output.filtered.regimes.exit_constrained_share(1,1)=0;
    for ka=2:length(output.filtered.regimes.exit)
        iden = length(find(isce(rr==ka)));
        if iden
            output.filtered.regimes.exit_constrained_share(ka,1) = length(find(isc(rr==ka)))/iden;
        else
            output.filtered.regimes.exit_constrained_share(ka,1) = 0;
        end
    end

    [aa, bb]=histcounts(simulated_sample.is_constrained(success),'normalization','pdf','binmethod','integers');
    if isscalar(aa)
        if bb(1)+0.5==0
            aa=[1 0];
        else
            aa=[0 1];
        end
    end
    output.filtered.regimes.is_constrained = isc;
    output.filtered.regimes.prob.is_constrained = aa(2);
    [aa, bb]=histcounts(simulated_sample.is_constrained_in_expectation(success),'normalization','pdf','binmethod','integers');
    if isscalar(aa)
        if bb(1)+0.5==0
            aa=[1 0];
        else
            aa=[0 1];
        end
    end
    output.filtered.regimes.is_constrained_in_expectation = isce;
    output.filtered.regimes.prob.is_constrained_in_expectation = aa(2);

    output.data.marginal_probability_distribution = density;
    output.data.marginal_probability = density_data;
    output.data.likelihood.empirical_kernel =  llik.kernel;
    output.data.likelihood.empirical_non_parametric = llik.non_parametric;
    output.data.likelihood.empirical_weighted_average = llik.ppf;
    output.data.likelihood.ppf = liky;
    output.data.likelihood.pkf = likpkf;
    output.data.likelihood.enkf = llik.enkf;

    for jk=1:size(StateVectors.Draws,1)
        [output.updated.variables.mean(jk,1), ...
            output.updated.variables.median(jk,1), ...
            output.updated.variables.var(jk,1), ...
            output.updated.variables.hpd_interval(jk,:), ...
            output.updated.variables.deciles(jk,:), ...
            output.updated.variables.density(jk,:,:)] = posterior_moments(StateVectors.Draws(jk,:)', options_.mh_conf_sig);
    end
    output.updated.variables.particles = StateVectors.Draws;
    output.updated.variables.pkf = alphahaty(:,2);
    output.updated.regimes.sample = updated_sample.regimes(updated_sample.indx);
    [aa,bb]=histcounts(updated_sample.regime_exit(updated_sample.indx),'normalization','pdf','binmethod','integers');
    output.updated.regimes.exit(bb(1:end-1)+0.5,1) = aa';
    isc = updated_sample.is_constrained(updated_sample.indx);
    isce = updated_sample.is_constrained_in_expectation(updated_sample.indx);
    rr=updated_sample.regime_exit(updated_sample.indx);
    output.updated.regimes.exit_constrained_share(1,1)=0;
    for ka=2:length(output.updated.regimes.exit)
        iden = length(find(isce(rr==ka)));
        if iden
            output.updated.regimes.exit_constrained_share(ka,1) = length(find(isc(rr==ka)))/iden;
        else
            output.updated.regimes.exit_constrained_share(ka,1) = 0;
        end
    end
    [aa, bb]=histcounts(updated_sample.is_constrained(updated_sample.indx),'normalization','pdf','binmethod','integers');
    if isscalar(aa)
        if bb(1)+0.5==0
            aa=[1 0];
        else
            aa=[0 1];
        end
    end
    output.updated.regimes.is_constrained = isc;
    output.updated.regimes.prob.is_constrained = aa(2);

    [aa, bb]=histcounts(updated_sample.is_constrained_in_expectation(updated_sample.indx),'normalization','pdf','binmethod','integers');
    if isscalar(aa)
        if bb(1)+0.5==0
            aa=[1 0];
        else
            aa=[0 1];
        end
    end
    output.updated.regimes.is_constrained_in_expectation = isce;
    output.updated.regimes.prob.is_constrained_in_expectation = aa(2);

elseif ~options_.occbin.filter.particle.ensemble_kalman_filter
    if StateVectors.Variance_rank==0
        StateVectors.stop_particles = true;
    end
end
end

function [liky, StateVectors, sim_regimes, sim_sample, error_flag] = local_EnKF_iteration(StateVector0,ShockVectors,data_index,Z,Y,H,M_,dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
% local_EnKF_iteration runs the ensemble Kalman filter iteration for the PPF EnKF branch.
%
% INPUTS
% - StateVector0           struct   Current particle state structure with draws and metadata.
% - ShockVectors           matrix   Random shocks drawn per particle.
% - data_index             cell     Indices used to select the observable subset.
% - Z                      matrix   Observation matrix linking states to observables.
% - Y                      matrix   Observed data matrix.
% - H                      matrix   Measurement noise covariance.
% - M_                     struct   Dynare internal model representation.
% - dr                     struct   Decision rule structure.
% - endo_steady_state      vector   Endogenous steady-state levels.
% - exo_steady_state       vector   Exogenous steady-state levels.
% - exo_det_steady_state   vector   Deterministic exogenous steady-state terms.
% - options_               struct   Dynare options structure.
% - occbin_options         struct   OccBin-specific options container.
%
% OUTPUTS
% - liky                   vector   Ensemble Kalman filter log-likelihood contributions per particle.
% - StateVectors           struct   Updated particle draws, mean, and variance.
% - sim_regimes            struct   Regime information collected from successful simulations.
% - sim_sample             struct   Simulation diagnostics for each particle.
% - error_flag             logical  Flag that becomes true when no successful particle is available.

base_regime = struct();
if M_.occbin.constraint_nbr==1
    base_regime.regime = 0;
    base_regime.regimestart = 1;
else
    base_regime.regime1 = 0;
    base_regime.regimestart1 = 1;
    base_regime.regime2 = 0;
    base_regime.regimestart2 = 1;
end

number_of_particles = size(StateVector0.Draws,2);

di=data_index{1};
ZZ = Z(di,:);
opts_simul = occbin_options.opts_simul;
if opts_simul.restrict_state_space
    my_order_var = dr.order_var(dr.restrict_var_list);
else
    my_order_var = dr.order_var;
end

% Ensemble Particle filter step ...
[liky, new_draws, sim_regimes, sim_sample, error_flag] = ...
    ensemble_kalman_filter_step(number_of_particles,StateVector0.Draws,ShockVectors, di, H, my_order_var, Y, ZZ, base_regime, M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state, options_, opts_simul);

StateVectors.Draws = new_draws;
StateVectors.Mean = mean(new_draws,2,'omitnan');
StateVectors.Var = cov(new_draws','partialrows');
end

function [lik_EnKF, y11_EnKF, simulated_regimes, simul_sample, error_flag] = ensemble_kalman_filter_step(number_of_particles, StateVector,ShockVectorsInfo, di, H, my_order_var, Y, ZZ, base_regime, M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state, options_, opts_simul)
% ensemble_kalman_filter_step computes the empirical conditional density from the simulated draws.
%
% INPUTS
% - number_of_particles    integer  Number of particles to process.
% - StateVector            matrix   Current particle state draws.
% - ShockVectorsInfo       matrix   Shock realizations for each particle.
% - di                     vector   Indices of observed variables.
% - H                      matrix   Measurement noise covariance.
% - my_order_var           vector   Ordered variable indices for the simulation output.
% - Y                      matrix   Observed data matching di.
% - ZZ                     matrix   Reduced observation matrix for di.
% - base_regime            struct   Regime identifier used for comparison.
% - M_                     struct   Dynare model structure.
% - dr                     struct   Decision rule structure.
% - endo_steady_state      vector   Endogenous steady-state levels.
% - exo_steady_state       vector   Exogenous steady-state levels.
% - exo_det_steady_state   vector   Deterministic exogenous steady-state terms.
% - options_               struct   Dynare options structure.
% - opts_simul             struct   Options used for simulation (including shocks and initial states).
%
% OUTPUTS
% - lik_EnKF               vector   Log-likelihood contributions for each simulated particle.
% - y11_EnKF               matrix   Filtered particle draws after applying the Kalman gain.
% - simulated_regimes      struct   Distinct regimes encountered during simulation.
% - simul_sample           struct   Simulation diagnostics for each particle.
% - error_flag             logical  Flag raised when no particle successfully simulated.

error_flag=0;
ShockVectors = ShockVectorsInfo;
simul_sample.success = false(number_of_particles,1);
lik_EnKF = zeros(number_of_particles,1);
y11_EnKF = zeros(size(StateVector));
opts_simul.SHOCKS=[];
number_of_simulated_regimes = 0;
for k = 1:number_of_particles
    if opts_simul.restrict_state_space
        tmp=zeros(M_.endo_nbr,1);
        tmp(dr.restrict_var_list,1)=StateVector(:,k);
        opts_simul.endo_init = tmp(dr.inv_order_var,1);
    else
        opts_simul.endo_init = StateVector(dr.inv_order_var,k);
    end
    options_.noprint = true;
    opts_simul.SHOCKS(1,:) = ShockVectors(:,k)';
    options_.occbin.simul=opts_simul;
    options_.occbin.simul.piecewise_only = false;
    [~, out] = occbin.solver(M_,options_,dr ,endo_steady_state, exo_steady_state, exo_det_steady_state);
    if out.error_flag==0
        simul_sample.success(k) = true;
        % store particle
        [~, ~, tmp_str]=occbin.backward_map_regime(out.regime_history(1));
        if M_.occbin.constraint_nbr==1
            tmp_str = tmp_str(2,:);
        else
            tmp_str = [tmp_str(2,:)  tmp_str(4,:)];
        end
        if number_of_simulated_regimes>0
            this_simulated_regime = find(ismember(all_simulated_regimes,{tmp_str}));
        end
        if number_of_simulated_regimes==0 || isempty(this_simulated_regime)
            number_of_simulated_regimes = number_of_simulated_regimes+1;
            this_simulated_regime = number_of_simulated_regimes;
            simulated_regimes(number_of_simulated_regimes).index = k;
            if isequal(base_regime,out.regime_history(1))
                simulated_regimes(number_of_simulated_regimes).is_base_regime = true;
            else
                simulated_regimes(number_of_simulated_regimes).is_base_regime = false;
            end
            simulated_regimes(number_of_simulated_regimes).regime = tmp_str;
            all_simulated_regimes = {simulated_regimes.regime};

        else
            simulated_regimes(this_simulated_regime).index = [simulated_regimes(this_simulated_regime).index k];
        end
        simul_sample.regime{k} = out.regime_history;
        simul_sample.y10(:,k) = out.piecewise(1,my_order_var)-out.ys(my_order_var)';
        simul_sample.y10lin(:,k) = out.linear(1,my_order_var)-out.ys(my_order_var)';

        if M_.occbin.constraint_nbr==1
            simul_sample.exit(k,1) = max(out.regime_history.regimestart);
            simul_sample.is_constrained(k,1) = logical(out.regime_history.regime(1));
            simul_sample.is_constrained_in_expectation(k,1) = any(out.regime_history.regime);
        else
            simul_sample.exit(k,1) = max(out.regime_history.regimestart1);
            simul_sample.exit(k,2) = max(out.regime_history.regimestart2);
            simul_sample.is_constrained(k,1) = logical(out.regime_history.regime1(1));
            simul_sample.is_constrained(k,2) = logical(out.regime_history.regime2(1));
            simul_sample.is_constrained_in_expectation(k,1) = any(out.regime_history.regime1);
            simul_sample.is_constrained_in_expectation(k,2) = any(out.regime_history.regime2);
        end
        % the distribution of ZZ*y10 needs to be confronted with ZZ*a1y, ZZ*P1y*ZZ',
        % i.e. the normal approximation using the regime consistent with the observable!
    end
end

success = simul_sample.success;
number_of_successful_particles = sum(success);

if number_of_successful_particles
    % ensemble KF
    mu=mean(simul_sample.y10(:,success),2,'omitnan');
    C=cov(simul_sample.y10(:,success)','partialrows');
    uF = ZZ*C*ZZ' + H(di,di);
    sig=sqrt(diag(uF));
    vv = Y(di) - ZZ*mu;
    if options_.rescale_prediction_error_covariance
        log_duF = log(det(uF./(sig*sig')))+2*sum(log(sig));
        iuF = inv(uF./(sig*sig'))./(sig*sig');
    else
        log_duF = log(det(uF));
        iuF = inv(uF);
    end
    lik_EnKF = log_duF + transpose(vv)*iuF*vv + length(di)*log(2*pi);
    Kalman_gain = C*ZZ'*iuF;
    y11_EnKF = simul_sample.y10*nan;
    y11_EnKF(:,success) = simul_sample.y10(:,success)+Kalman_gain*(Y(di)-ZZ*simul_sample.y10(:,success));
else
    error_flag=1;
end
end
