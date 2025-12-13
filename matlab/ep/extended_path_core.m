function [y1, info_convergence, endo_simul, y, pfm, options_] = ...
    extended_path_core(exo_simul,initial_conditions ,...
                        pfm,M_, options_, oo_, initialguess, y)
% [y1, info_convergence, endo_simul, y, pfm, options_] = ...
%     extended_path_core(exo_simul,initial_conditions ,...
%     pfm,M_, options_, oo_, initialguess, y)
% INPUTS
%  o  exo_simul             [matrix]    path of exogenous, used to construct the guess values (only if oo_.deterministic_simulation.controlled_paths_by_period is not empty)
%  o  initial_conditions    [matrix]    path of endogenous, used to construct the guess values (initial condition not used; terminal condition used as guess value iff recompute_final_steady_state=true)
%  o  pfm                   [struct]    perfect foresight model description
%  o  M_                    [structure] describing the model
%  o  options_              [structure] describing the options
%  o  oo_                   [struct]    Dynare's results structure
%  o  y                     [vector]    initial guess
%
% OUTPUTS
%  o  y                [vector]    solution for current period
%  o  info_convergence [Boolean]   scalar if simulation was successful
%  o  endo_simul       [matrix]    path of endogenous
%  o  errorcode        [integer]   error code
%  o  y                [vector]    solution for current period
%  o  pfm              [struct]    perfect foresight model description
%  o  options_         [structure] describing the options


% Copyright © 2016-2025 Dynare Team
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

ep = options_.ep;

debug = options_.verbosity;
order = options_.ep.stochastic.order;

if options_.ep.use_first_order_solution_as_initial_guess% Compute first order solution (Perturbation)...
    endo_simul = simult_(M_,options_,initial_conditions,oo_.dr,exo_simul(2:end,:),1);
else
    if nargin>6 && ~isempty(initialguess)
        % Note that the first column of initialguess should be equal to initial_conditions.
        endo_simul = initialguess;
    else
        endo_simul = [initial_conditions repmat(oo_.steady_state,1,ep.periods+M_.maximum_lead)];
    end
end

if nargin~=8
    y = [];
end

oo_.endo_simul = endo_simul;

if debug
    save ep_test_1.mat endo_simul exo_simul
end

if options_.bytecode && order > 0
    error('extended path: order > 0 is not compatible with bytecode option.')
end
if options_.block && order > 0
    error('extended path: order > 0 is not compatible with block option.')
end

if order == 0
    % Extended Path
    options_.periods = ep.periods;
    options_.block = pfm.block;
    oo_.endo_simul = endo_simul;
    oo_.exo_simul = exo_simul;
    options_.solve_algo = ep.solve_algo;
    options_.stack_solve_algo = ep.stack_solve_algo;
    [endo_simul, info_convergence] = perfect_foresight_solver_core(oo_.endo_simul, oo_.exo_simul, oo_.steady_state, oo_.exo_steady_state, [], M_, options_);
else
    % Stochastic Extended Path
    switch(ep.stochastic.algo)
      case 0
        % Full tree of future trajectories.
        if nargout>4
            [flag, endo_simul, errorcode, y, pfm, options_] = solve_stochastic_perfect_foresight_model_0(endo_simul, exo_simul, y, options_, M_, pfm);
        else
            [flag, endo_simul, errorcode, y] = solve_stochastic_perfect_foresight_model_0(endo_simul, exo_simul, y, options_, M_, pfm);
        end
      case 1
        % Sparse tree of future histories.
        if nargout>4
            [flag, endo_simul, errorcode, y, pfm, options_] = solve_stochastic_perfect_foresight_model_1(endo_simul, exo_simul, y, options_, M_, pfm);
        else
            [flag, endo_simul, errorcode, y] = solve_stochastic_perfect_foresight_model_1(endo_simul, exo_simul, y, options_, M_, pfm);
        end
    end
    info_convergence = ~flag;
end

if ~info_convergence && ~options_.no_homotopy
    [info_convergence, endo_simul] = extended_path_homotopy(endo_simul, exo_simul, M_, options_, oo_, pfm, ep, order, ep.stochastic.algo, 2, debug);
end

if info_convergence
    y1 = endo_simul(:,2);
else
    y1 = NaN(size(M_.endo_nbr,1));
end
