function measure = measurement_equations(StateVectors,ReducedForm,ThreadsOptions, options_, M_)
% measure = measurement_equations(StateVectors,ReducedForm,ThreadsOptions, options_, M_)
% Get measurement of var
% INPUTS
%  - StateVectors           [double]    value of the state variables
%  - ReducedForm            [structure] MATLAB's structure describing the reduced form model.
%  - ThreadsOptions         [structure] options for threading of mex files
%  - options_               [structure] describing the options
%  - M_                     [structure] describing the model
%
% OUTPUTS
%  - measure                [double]    scalar, likelihood


% Copyright © 2013-2025 Dynare Team
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

number_of_structural_innovations = length(ReducedForm.Q);
yhat = bsxfun(@minus, StateVectors, ReducedForm.state_variables_steady_state);
if ReducedForm.use_k_order_solver
    tmp = local_state_space_iteration_k(yhat, zeros(number_of_structural_innovations, size(yhat,2)), ReducedForm.dr, M_, options_, ReducedForm.udr);
    measure = tmp(ReducedForm.mf1,:);
else
    if options_.order == 2
        measure = local_state_space_iteration_2(yhat, zeros(number_of_structural_innovations, size(yhat,2)), ReducedForm.ghx, ReducedForm.ghu, ReducedForm.constant, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ThreadsOptions.local_state_space_iteration_2);
    elseif options_.order == 3
        measure = local_state_space_iteration_3(yhat, zeros(number_of_structural_innovations, size(yhat,2)), ReducedForm.ghx, ReducedForm.ghu, ReducedForm.ghxx, ReducedForm.ghuu, ReducedForm.ghxu, ReducedForm.ghs2, ReducedForm.ghxxx, ReducedForm.ghuuu, ReducedForm.ghxxu, ReducedForm.ghxuu, ReducedForm.ghxss, ReducedForm.ghuss, ReducedForm.steadystate, ThreadsOptions.local_state_space_iteration_3, false);
    else
        error('Order > 3: use_k_order_solver should be set to true');
    end
end
