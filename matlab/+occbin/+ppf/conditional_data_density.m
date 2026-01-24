function [likxc, infox, likxmode, updated_mode, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes] = ...
    conditional_data_density(yhat, k, t, likxmode, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes, ...
    a0, P, P0, P1,my_data_index,Z,ZZ,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
    M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options)
% [likxc, infox, likxmode, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes] = ...
%     occbin.ppf.conditional_data_density(yhat, k, t, likxmode, updated_regimes, updated_sample, number_of_updated_regimes, all_updated_regimes, ...
%     P, P1,my_data_index,Z,v,my_data,H,QQQ,T0,R0,TT,RR,CC, regimesy,regimes0,base_regime, info, ...
%     M_, dr, endo_steady_state,exo_steady_state,exo_det_steady_state,options_,occbin_options);
%
% INPUTS
%  - yhat                   [double]    state vector in t-1 (particle)
%  - k                      [integer]   particle index 
%  - t                      [integer]   time index 
%  - likxmode               [double]    mode of log-likelihood
%  - updated_regimes        [struct]    occbin info about updated regimes (particles)
%  - updated_sample         [struct]    occbin info about updated state (particles)
%  - number_of_updated_regimes [integer] number of updated regimes across particles
%  - all_updated_regimes    [cell]      list of updated regimes 
%  - a0                     [double]    PKF updated state estimate t-1.
%  - P                      [double]    Particle state covariance t-1 (=0).
%  - P0                     [double]    PKF updated state covariance t-1.
%  - P1                     [double]    Particle 1 step ahead state covariance.
%  - my_data_index          [cell]      1*2 cell of column vectors of indices.
%  - Z                      [double]    Observation matrix.
%  - ZZ                     [double]    Observation selection matrix.
%  - v                      [double]    Observation innovations.
%  - my_data                [double]    observed data t-1:t
%  - H                      [double]    Measurement error variance.
%  - QQQ                    [double]    Shocks covariance matrix.
%  - T0                     [double]    Base regime state transition matrix.
%  - R0                     [double]    Base regime state shock matrix.
%  - TT                     [double]    PKF Transition matrices t-1:t given t-1 info.
%  - RR                     [double]    PKF Shock matrices t-1:t given t-1 info.
%  - CC                     [double]    PKF Constant terms in transition t-1:t given t-1 info.
%  - regimesy               [struct]    PKF updated regime info t:t+2.
%  - regimes0               [struct]    PKF regime info t:t+1 given t-1 info.
%  - base_regime            [integer]   Base regime info.
%  - info                   [integer]   PKF error flag.
%  - M_                     [structure] MATLAB's structure describing the model
%  - dr                     [structure] model information structure
%  - endo_steady_state      [vector]    steady state value for endogenous variables
%  - exo_steady_state       [vector]    steady state value for exogenous variables
%  - exo_det_steady_state   [vector]    steady state value for exogenous deterministic variables                                    
%  - options_               [structure] MATLAB's structure containing the options
%  - occbin_options         [structure] options structure for OccBin filter
%
% OUTPUTS
%  - likx                   [double]    minus 2*log-likelihood density
%  - infox                  [integer]   indicator of successful computation
%  - likxmode               [double]    mode of log-likelihood
%  - updated_mode           [struct]    occbin info about log-likelihood mode
%  - updated_regime         [struct]    occbin info about updated regime
%  - updated_sample         [struct]    occbin info about updated state draw
%  - number_of_updated_regimes [integer] 
%  - all_updated_regimes    [cell]      list of updated regimes 
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

a(:,1) = yhat;
a1(:,1) = yhat;
a1(:,2) = TT(:,:,2)*yhat;
di=my_data_index{2};

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
updated_mode=struct();
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
        if options_.occbin.filter.particle.diagnostics.status
            updated_regimes(number_of_updated_regimes).obsvar = ZZ*(Tx(:,:,1)*P0*Tx(:,:,1)'+Rx(:,:,1)*QQQ(:,:,2)*Rx(:,:,1)')*ZZ' + H(di,di);
            updated_regimes(number_of_updated_regimes).obsmean = ZZ*(Tx(:,:,1)*a0+Cx(:,1));
            updated_regimes(number_of_updated_regimes).ss.C = Cx(:,1);
            updated_regimes(number_of_updated_regimes).ss.R = Rx(:,:,1);
            updated_regimes(number_of_updated_regimes).ss.T = Tx(:,:,1);
            uF = updated_regimes(number_of_updated_regimes).obsvar;
            sig=sqrt(diag(uF));
            vv = my_data(di,2) - updated_regimes(number_of_updated_regimes).obsmean;
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

% % % %     number_of_successful_particles = number_of_successful_particles + 1;
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
        if options_.occbin.filter.particle.diagnostics.status
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
