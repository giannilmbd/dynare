function [nodes,W_m,W_c] = unscented_sigma_points(n,unscented_options)
% [nodes,W_m,W_c] = unscented_sigma_points(n,unscented_options)
% Computes nodes and weights for a scaled unscented transform cubature,
% i.e. compute the nodes and weights for a second-order accurate propagtion
% of a Gaussian variable through a nonlinear function
%
% INPUTS
%    n                  [integer]   scalar, number of variables.
%    unscented_options  [structure] hyperparameters
%
% OUTPUTS
%    nodes          [double]    nodes of the cubature
%    W_m            [double]    associated weights for the mean
%    W_c            [double]    associated weights for the covariance
%
% REFERENCES
%    Formulas follow the ones in Wan/van der Merwe (2001): "The unscented
%    Kalman filter", in Haykin (editor): Kalman Filtering and Neural
%    Networks, Chapter 7, p. 221-280.
%
% NOTES

% Copyright © 2009-2026 Dynare Team
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

lambda = (unscented_options.alpha^2)*(n+unscented_options.kappa) - n ; %below (7.30)
nodes = [zeros(n,1)  sqrt(n+lambda).*[eye(n) -eye(n)]]' ; %prefactor befor P_x in (7.30)
% Implement weights from (7.34)
W_m = lambda/(n+lambda);
W_c = W_m + (1-unscented_options.alpha^2+unscented_options.beta);
temp = ones(2*n,1)/(2*(n+lambda)) ;
W_m = [W_m ; temp] ;
W_c = [W_c ; temp] ;
