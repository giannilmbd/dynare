function [density] = probability2(mu,S,X)
% [density] = probability2(mu,S,X)
% Multivariate gaussian density
%
% INPUTS
%    mu             [double]   mean of distribution
%    S              [double]   lower Cholesky factor of covariance matrix of distribution
%    X              [double]   random vector to be evaluated
%
% OUTPUTS
%    density        [double]   Gaussian pdf of X

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

dim = size(X,1) ;
normfact = bsxfun(@power,(2*pi),(dim/2)) ;
foo = S\(bsxfun(@minus,X,mu)) ;
density = exp(-0.5*sum(foo.*foo)')./abs((normfact*prod(diag(S)))) + 1e-99 ;

%@test:1
%$ x=[2.3;0.34;0.87];
%$ x=[2.3;0.34;0.87];
%$ mu=zeros(3,1);
%$ Sigma=[1 0.5 0; 0.5 1 0; 0 0 0.05];
%$ L=chol(Sigma,"lower");
%$ pdf1=probability2(mu,L,x);
%$ pdf2=mvnpdf(x,mu,Sigma);
% Check results.
%$ t(1) = dassert(pdf1,pdf2,1e-10);
%$ T = all(t);
%@eof:2
