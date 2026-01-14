function [llik, simulated_regimes, simul_sample, hp_simul_regimes] = simulated_density(number_of_particles, number_of_shocks_per_particle,StateVector,ShockVectorsInfo, di, H, my_order_var, QQQ, Y, ZZ, base_regime, regimesy, M_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, options_, opts_simul)
% [llik, simulated_regimes, simul_sample, hp_simul_regimes] = simulated_density(number_of_particles, number_of_shocks_per_particle,StateVector,ShockVectorsInfo, di, H, my_order_var, QQQ, Y, ZZ, base_regime, regimesy, M_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, options_, opts_simul)
%
% Computes empirical conditional density for the piecewise linear particle filter (PPF).
%
% INPUTS
% - number_of_particles     [integer]   Number of particles to simulate.
% - number_of_shocks_per_particle [integer] Number of shocks per particle.
% - StateVector             [double]    State vector for particles (particles by column).
% - ShockVectorsInfo        [double|struct] Shock draws or structure with variance decomposition.
% - di                      [integer]   Index for data subset.
% - H                       [double]    Measurement error variance.
% - my_order_var            [integer]   Variable ordering index.
% - QQQ                     [double]    Shocks covariance matrix.
% - Y                       [double]    Observations t-1:t.
% - ZZ                      [double]    Observation selection matrix.
% - base_regime             [struct]    Base regime info.
% - regimesy                [struct]    PKF updated regime info.
% - M_                      [struct]    Dynare's model structure.
% - dr                      [struct]    Decision rule structure.
% - endo_steady_state       [double]    Endogenous steady state.
% - exo_steady_state        [double]    Exogenous steady state.
% - exo_det_steady_state    [double]    Exogenous deterministic steady state.
% - options_                [struct]    Dynare options.
% - opts_simul              [struct]    Simulation options.
%
% OUTPUTS
% - llik                    [struct]    Likelihood contributions (kernel, non_parametric, ppf, enkf).
% - simulated_regimes       [cell]      Regime information for each particle.
% - simul_sample            [cell]      Sample information for each particle.
% - hp_simul_regimes        [struct]    Higher-order particle regime information.
%
% REMARKS
% When number_of_shocks_per_particle>1, ShockVectorsInfo must be a structure with fields:
%     .US                   [double]    Left singular vectors from SVD of shock covariance.
%     .VarianceSquareRoot   [double]    Cholesky factor of shock variance.
%     .VarianceRank         [integer]   Rank of shock covariance matrix.

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

conditional_loop = false;
number_of_iterations = number_of_particles;
if number_of_shocks_per_particle==1
    hp_simul_regimes = [];
    ShockVectors = ShockVectorsInfo;
    simul_sample.success = false(number_of_particles,1);
elseif number_of_particles==1
    number_of_iterations = number_of_shocks_per_particle;
    if ~isstruct(ShockVectorsInfo)
        % for non parametric conditional 
        ShockVectors = ShockVectorsInfo;
        conditional_loop = true;        
    end
end
llik.kernel = zeros(number_of_particles,1);
llik.non_parametric = zeros(number_of_particles,1);
llik.ppf = zeros(number_of_particles,1);
llik.enkf = zeros(number_of_particles,1);
if ~conditional_loop
    y11 = zeros(size(StateVector));
    y11_EnKF = zeros(size(StateVector));
else
    a0 = StateVector;
end
opts_simul.SHOCKS=[];
number_of_simulated_regimes = 0;
for k = 1:number_of_iterations
    if ~conditional_loop
        a0 = StateVector(:,k);
        if opts_simul.restrict_state_space
            tmp=zeros(M_.endo_nbr,1);
            tmp(dr.restrict_var_list,1)=StateVector(:,k);
            opts_simul.endo_init = tmp(dr.inv_order_var,1);
        else
            opts_simul.endo_init = StateVector(dr.inv_order_var,k);
        end
    end
    if number_of_shocks_per_particle>1 && isstruct(ShockVectorsInfo)
        ShockVectorsPerParticle = bsxfun(@plus,ShockVectorsInfo.US*ShockVectorsInfo.VarianceSquareRoot*transpose(norminv(qmc_scrambled(ShockVectorsInfo.VarianceRank,number_of_shocks_per_particle,1))),zeros(ShockVectorsInfo.VarianceRank,1));
        [llik(k), simulated_regimes{k}, simul_sample{k}, hp_simul_regimes(k)] = occbin.ppf.simulated_density(1, number_of_shocks_per_particle,a0,ShockVectorsPerParticle, di, H, my_order_var, QQQ, Y, ZZ, base_regime, regimesy, M_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state, options_, opts_simul);
    else
        options_.noprint = true;
        opts_simul.SHOCKS(1,:) = ShockVectors(:,k)';
        options_.occbin.simul=opts_simul;
        options_.occbin.simul.piecewise_only = false;
        [~, out, ss] = occbin.solver(M_,options_,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
        if out.error_flag==0
            simul_sample.success(k) = true;
            % % store particle
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
                if isequal(regimesy(1),out.regime_history(1))
                    simulated_regimes(number_of_simulated_regimes).is_pkf_regime = true;
                else
                    simulated_regimes(number_of_simulated_regimes).is_pkf_regime = false;
                end
                if isequal(base_regime,out.regime_history(1))
                    simulated_regimes(number_of_simulated_regimes).is_base_regime = true;
                else
                    simulated_regimes(number_of_simulated_regimes).is_base_regime = false;
                end
                simulated_regimes(number_of_simulated_regimes).obsvar = ZZ*(ss.R(my_order_var,:)*QQQ(:,:,2)*ss.R(my_order_var,:)')*ZZ' + H(di,di);
                simulated_regimes(number_of_simulated_regimes).obsmean = ZZ*(ss.T(my_order_var,my_order_var)*a0+ss.C(my_order_var));
                simulated_regimes(number_of_simulated_regimes).regime = tmp_str;
                simulated_regimes(number_of_simulated_regimes).ss.C = ss.C(my_order_var);
                simulated_regimes(number_of_simulated_regimes).ss.R = ss.R(my_order_var,:);
                simulated_regimes(number_of_simulated_regimes).ss.T = ss.T(my_order_var,my_order_var);
                all_simulated_regimes = {simulated_regimes.regime};
                uF = simulated_regimes(number_of_simulated_regimes).obsvar;
                sig=sqrt(diag(uF));
                vv = Y(di,2) - simulated_regimes(number_of_simulated_regimes).obsmean;
                if options_.rescale_prediction_error_covariance
                    log_duF = log(det(uF./(sig*sig')))+2*sum(log(sig));
                    iuF = inv(uF./(sig*sig'))./(sig*sig');
                else
                    log_duF = log(det(uF));
                    iuF = inv(uF);
                end
                simulated_regimes(number_of_simulated_regimes).lik = ...
                    log_duF + transpose(vv)*iuF*vv + length(di)*log(2*pi);
                simulated_regimes(number_of_simulated_regimes).Kalman_gain = (ss.R(my_order_var,:)*QQQ(:,:,2)*ss.R(my_order_var,:)')*ZZ'*iuF;
                if conditional_loop
                    simulated_regimes(number_of_simulated_regimes).update = (ss.T(my_order_var,my_order_var)*a0+ss.C(my_order_var)) + simulated_regimes(number_of_simulated_regimes).Kalman_gain*vv;
                end

            else
                simulated_regimes(this_simulated_regime).index = [simulated_regimes(this_simulated_regime).index k];
            end
            simul_sample.regime{k} = out.regime_history;
            simul_sample.y10(:,k) = out.piecewise(1,my_order_var)-out.ys(my_order_var)';
            simul_sample.y10lin(:,k) = out.linear(1,my_order_var)-out.ys(my_order_var)';
            if conditional_loop
                simul_sample.y11(:,k) = simulated_regimes(this_simulated_regime).update;
            else
                % this is in the spirit of EnKF, but with regime specific
                % Kalman gain
                omean = (simulated_regimes(this_simulated_regime).ss.T*a0+simulated_regimes(this_simulated_regime).ss.C);
                simul_sample.y11(:,k) = omean + simulated_regimes(this_simulated_regime).Kalman_gain*(Y(di,2) - ZZ*omean);
            end
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
end

if number_of_shocks_per_particle==1 || number_of_particles==1
    success = simul_sample.success;
    number_of_successful_particles = sum(success);


    if number_of_successful_particles
        nobs = size(ZZ,1);
        z10 = ZZ*simul_sample.y10(:,success);

        % ensemble KF
        mu=mean(simul_sample.y10(:,success),2);
        C=cov(simul_sample.y10(:,success)');
        uF = ZZ*C*ZZ' + H(di,di);
        sig=sqrt(diag(uF));
        vv = Y(di,2) - ZZ*mu;
        if options_.rescale_prediction_error_covariance
            log_duF = log(det(uF./(sig*sig')))+2*sum(log(sig));
            iuF = inv(uF./(sig*sig'))./(sig*sig');
        else
            log_duF = log(det(uF));
            iuF = inv(uF);
        end
        llik.enkf = log_duF + transpose(vv)*iuF*vv + length(di)*log(2*pi);
        if number_of_shocks_per_particle==1
            Kalman_gain = C*ZZ'*iuF;
            y11_EnKF = simul_sample.y11*nan;
            y11_EnKF(:,success) = simul_sample.y10(:,success)+Kalman_gain*(Y(di,2)-ZZ*simul_sample.y10(:,success));
        end
        % full non parametric density estimate
        dlik = log_duF+diag(transpose(Y(di,2)-z10)*iuF*(Y(di,2)-z10))+ length(di)*log(2*pi);
        [sd, iss] = sort(dlik);
        issux = find(success);
        issux = issux(iss);

        ydist = sqrt(sum((Y(di,2)-z10).^2,1));
        [s, is] = sort(ydist);
        isux = find(success);
        isux = isux(is);
        %bin population
        nbins = ceil(sqrt(size(z10,2)));

        binpop = ceil(number_of_successful_particles/nbins/2)*2; % I center the bin on the observed point

        % volume of n-dimensional sphere
        binvolume = pi^(nobs/2)/gamma(nobs/2+1)*(0.5*(s(binpop)+s(binpop+1)))^nobs;
        if number_of_shocks_per_particle==1
            llik.non_parametric = -2*log(binpop/number_of_successful_particles/binvolume);

            % kernel density
            fobs_kernel = nan;
            if ~isoctave && user_has_matlab_license('statistics_toolbox') % ksdensity is missing
                if nobs>2
                    bw = nan(nobs,1);
                    for jd=1:nobs
                        ss=std(transpose(z10(jd,:)));
                        bw(jd) = ss*(4/(nobs+2)/number_of_successful_particles)^(1/(nobs+4));
                    end
                    fobs_kernel = mvksdensity(z10',Y(di,2)','bandwidth',bw,'Function','pdf');
                else
                    fobs_kernel = ksdensity(z10',Y(di,2)');
                end
            end
            llik.kernel = -2*log(fobs_kernel);
        end

        % check for regime of Y
        is_data_regime = false(number_of_simulated_regimes,1);
        iss_data_regime = false(number_of_simulated_regimes,1);
        proba = ones(number_of_simulated_regimes,1);
        probaX = zeros(number_of_simulated_regimes,1);
        probaY = zeros(number_of_simulated_regimes,1);
        for kr=1:number_of_simulated_regimes
            if options_.occbin.filter.particle.tobit && ~simulated_regimes(kr).is_base_regime
                % base regime is not Tobit
                proba(kr) = length(simulated_regimes(kr).index)/number_of_successful_particles;
            end
            % here we check if the data may be spanned for each regime
            % (given we have a bunch of TRUNCATED normals, Y may be out of
            % truncation ...)
            is_data_regime(kr) = any(ismember(simulated_regimes(kr).index,isux(1:binpop)));
            iss_data_regime(kr) = any(ismember(simulated_regimes(kr).index,issux(1:binpop)));
            if is_data_regime(kr)
                probaY(kr) = length(find(ismember(simulated_regimes(kr).index,isux(1:binpop))))/binpop;
                probaY(kr) = length(find(ismember(simulated_regimes(kr).index,isux(1:binpop))))^2/binpop/length(simulated_regimes(kr).index);
            end
            if iss_data_regime(kr)
                thisbin = find(ismember(issux(1:binpop),simulated_regimes(kr).index));
                probaX(kr) = mean(exp(-sd(thisbin)/2))*length(thisbin)/length(simulated_regimes(kr).index);
            end
        end
        probaY = probaY./sum(probaY);
        if ~any(simul_sample.is_constrained_in_expectation(success))
            % Y is unconstrained
        end

        % use most probable regime to update states
        [~,im]=max((exp(-([simulated_regimes.lik]-min([simulated_regimes.lik]))./2)'.*(proba.*probaY)));
        maxbin = ismember(isux(1:binpop),simulated_regimes(im).index);
        simul_sample.y11_hpregime = mean(simul_sample.y11(:,issux(maxbin)),2);
        hp_simul_regimes = simulated_regimes(im);
        if number_of_shocks_per_particle==1

            simul_sample.y11(:,success) = StateVector(:,success) + simulated_regimes(im).Kalman_gain*(Y(di,2) - ZZ*StateVector(:,success));
            % likelihood is the sum of Gaussians, where constrained regimes are
            % weighted by their probabilities
            llik.ppf = -2*log(exp(-[simulated_regimes.lik]./2)*(proba.*probaY));
        end

        if conditional_loop
            % likelihood is the sum of Gaussians, where constrained regimes are
            % weighted by their probabilities
            llik.ppf = -2*(log(exp(-([simulated_regimes.lik]-min([simulated_regimes.lik]))./2)*(proba.*probaY))-min([simulated_regimes.lik]));
            llik.non_parametric = -2*log(binpop/number_of_successful_particles/binvolume);
        else
            simul_sample.y11_EnKF = y11_EnKF;

        end
    end

end
