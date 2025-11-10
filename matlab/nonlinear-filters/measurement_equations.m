function measure = measurement_equations(StateVectors,ReducedForm,options_, M_)
% measure = measurement_equations(StateVectors,ReducedForm,options_, M_)
% Get measurement of var
% INPUTS
%  - StateVectors           [double]    value of the state variables
%  - ReducedForm            [structure] MATLAB's structure describing the reduced form model.
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
measure=iterate_law_of_motion(StateVectors,zeros(number_of_structural_innovations, size(StateVectors,2)),ReducedForm,M_,options_,ReducedForm.use_k_order_solver,false);
measure = measure(ReducedForm.mf1,:);
