function [out,info] = get_perturbation_params_derivs_numerical_objective(params, outputflag, estim_params, M_, options_, dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
%function [out,info] = get_perturbation_params_derivs_numerical_objective(params, outputflag, estim_params, M_, options_, dr, steady_state, exo_steady_state, exo_det_steady_state)
% -------------------------------------------------------------------------
% Objective function to compute numerically the Jacobians used for get_perturbation_params_derivs
% =========================================================================
% INPUTS
%   params:         [vector] parameter values at which to evaluate objective function
%                   stderr parameters come first, corr parameters second, model parameters third
%   outputflag:     [string] flag which objective to compute (see below)
%   estim_params:   [structure] storing the estimation information
%   M_:             [structure] storing the model information
%   options_:       [structure] storing the options
%   dr              [structure]     Reduced form model.
%   endo_steady_state       [vector]     steady state value for endogenous variables
%   exo_steady_state        [vector]     steady state value for exogenous variables
%   exo_det_steady_state    [vector]     steady state value for exogenous deterministic variables                                    
% -------------------------------------------------------------------------
%
% OUTPUT 
%   out (dependent on outputflag and order of approximation):
%     - 'perturbation_solution':  out = out1 = [vec(Sigma_e);vec(ghx);vec(ghu)]; (order==1)
%                                 out = out2 = [out1;vec(ghxx);vec(ghxu);vec(ghuu);vec(ghs2)]; (order==2)
%                                 out = out3 = [out1;out2;vec(ghxxx);vec(ghxxu);vec(ghxuu);vec(ghuuu);vec(ghxss);vec(ghuss)]; (order==3)
%     - 'dynamic_model':          out = [Yss; vec(g1); vec(g2); vec(g3)]
%     - 'Kalman_Transition':      out = [Yss; vec(KalmanA); dyn_vech(KalmanB*Sigma_e*KalmanB')];
%     all in DR-order
%   info            [integer] output from resol
% -------------------------------------------------------------------------
% This function is called by
%   * get_perturbation_params_derivs.m (previously getH.m)
% -------------------------------------------------------------------------
% This function calls
%   * [M_.fname,'.dynamic']
%   * resol
%   * dyn_vech

% Copyright © 2019-2025 Dynare Team
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

%% Update stderr, corr and model parameters and compute perturbation approximation and steady state with updated parameters
M_ = set_all_parameters(params,estim_params,M_);
[dr,info,M_.params] = compute_decision_rules(M_,options_,dr,endo_steady_state,exo_steady_state,exo_det_steady_state);
Sigma_e = M_.Sigma_e;

if info(1) > 0
    % there are errors in the solution algorithm
    out = [];
    return
else
    ys = dr.ys; %steady state of model variables in declaration order
    ghx = dr.ghx; ghu = dr.ghu;
    if options_.order > 1
        ghxx = dr.ghxx; ghxu = dr.ghxu; ghuu = dr.ghuu; ghs2 = dr.ghs2;
    end
    if options_.order > 2
        ghxxx = dr.ghxxx; ghxxu = dr.ghxxu; ghxuu = dr.ghxuu; ghxss = dr.ghxss; ghuuu = dr.ghuuu; ghuss = dr.ghuss;
    end
end
Yss = ys(dr.order_var); %steady state of model variables in DR order

%% out = [vec(Sigma_e);vec(ghx);vec(ghu);vec(ghxx);vec(ghxu);vec(ghuu);vec(ghs2);vec(ghxxx);vec(ghxxu);vec(ghxuu);vec(ghuuu);vec(ghxss);vec(ghuss)]
if strcmp(outputflag,'perturbation_solution')
    out = [Sigma_e(:); ghx(:); ghu(:)];
    if options_.order > 1
        out = [out; ghxx(:); ghxu(:); ghuu(:); ghs2(:);];
    end
    if options_.order > 2
        out = [out; ghxxx(:); ghxxu(:); ghxuu(:); ghuuu(:); ghxss(:); ghuss(:)];
    end
end

%% out = [Yss; vec(g1); vec(g2); vec(g3)]; of all endogenous variables, in DR order
if strcmp(outputflag,'dynamic_model')
    y3n = repmat(ys, 3, 1);
    [g1, T, T_order] = feval([M_.fname,'.sparse.dynamic_g1'], y3n, exo_steady_state', M_.params, ys, M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr);
    g1 = identification.legacy_dynamic_g1(g1, M_);
    out = [Yss; g1(:)];
    if options_.order >= 2
        [g2_v, T, T_order] = feval([M_.fname,'.sparse.dynamic_g2'], y3n, exo_steady_state', M_.params, ys, T, T_order);
        g2 = identification.legacy_dynamic_g2(g2_v, M_);
        out = [out; g2(:)];
    end
    if options_.order == 3
        g3_v = feval([M_.fname,'.sparse.dynamic_g3'], y3n, exo_steady_state', M_.params, ys, T, T_order);
        g3 = identification.legacy_dynamic_g3(g3_v, M_);
        g3 = identification.unfold_g3(g3, M_.nspred + M_.endo_nbr + M_.nsfwrd + M_.exo_nbr);
        out = [out; g3(:)];
    end
end

%% out = [Yss; vec(KalmanA); dyn_vech(KalmanB*Sigma_e*KalmanB')]; in DR order, where A and B are Kalman transition matrices
if strcmp(outputflag,'Kalman_Transition')
    if options_.order == 1
        KalmanA = zeros(M_.endo_nbr,M_.endo_nbr);
        KalmanA(:,M_.nstatic+(1:M_.nspred)) = ghx;
        Om = ghu*Sigma_e*transpose(ghu);
        out = [Yss; KalmanA(:); dyn_vech(Om)];
    else
        error('''get_perturbation_params_derivs_numerical_objective.m'': Kalman_Transition works only at order=1');
    end
end
