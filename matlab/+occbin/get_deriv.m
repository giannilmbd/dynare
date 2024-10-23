function [h_minus_1, h, h_plus_1, h_exo, resid] = get_deriv(M_, ys_)
% function [h_minus_1, h, h_plus_1, h_exo, resid] = get_deriv(M_, ys_)
% Computes dynamic Jacobian and writes it to conformable matrices
%
% INPUTS
% - M_         [struct]     Definition of the model.
% - ys         vector       steady state
%
% OUTPUTS
% - h_minus_1  [N by N]     derivative matrix with respect to lagged endogenous variables
% - h          [N by N]     derivative matrix with respect to contemporanous endogenous variables
% - h_plus_1   [N by N]     derivative matrix with respect to leaded endogenous variables
% - h_exo      [N by N_exo] derivative matrix with respect to exogenous variables
% - resid      [N by 1]     vector of residuals

% Copyright © 2021-2024 Dynare Team
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

dyn_endo_ss = repmat(ys_, 3, 1);
x = zeros(M_.exo_nbr);

[resid, T_order, T] = feval([M_.fname '.sparse.dynamic_resid'], dyn_endo_ss, x, M_.params, ys_);

g1 = feval([M_.fname '.sparse.dynamic_g1'], dyn_endo_ss, x, M_.params, ys_, ...
           M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, ...
           M_.dynamic_g1_sparse_colptr, T_order, T);

h_minus_1 = full(g1(:,1:M_.endo_nbr));
h = full(g1(:,M_.endo_nbr + (1:M_.endo_nbr)));
h_plus_1 = full(g1(:,2*M_.endo_nbr + (1:M_.endo_nbr)));
h_exo = full(g1(:,3*M_.endo_nbr+1:end));
