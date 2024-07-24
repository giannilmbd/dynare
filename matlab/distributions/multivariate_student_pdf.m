function density = multivariate_student_pdf(X,Mean,Sigma_upper_chol,df)
% Evaluates the density of a multivariate student, with expectation Mean,
% variance Sigma_upper_chol'*Sigma_upper_chol and degrees of freedom df, at X.
%
% INPUTS
%
%    X                  [double]    dim*n vector
%    Mean               [double]    1*n vector, expectation of the multivariate random variable.
%    Sigma_upper_chol   [double]    n*n matrix, upper triangular Cholesky decomposition of Sigma (the covariance 
%                                   matrix up to a factor df/(df-2)).
%    df                 [integer]   degrees of freedom.
%
% OUTPUTS
%    density            [double]    density.
%
% SPECIAL REQUIREMENTS

% Copyright © 2003-2024 Dynare Team
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
if df <=0
    error('Degrees of freedom ''df'' must be positive')
end
[~, nn] = size(X);
t1 = gamma( .5*(nn+df) ) / ( gamma( .5*df ) * (df*pi)^(.5*nn) );
t2 = t1 / prod(diag(Sigma_upper_chol)) ;
density = t2 ./ ( 1 + sum(((X-Mean)/Sigma_upper_chol).^2, 2)/df ).^(.5*(nn+df));

return % --*-- Unit tests --*--

%@test:1
% Normal density
try
    m1 = multivariate_student_pdf([1 2],0,chol([1 0.5; 0.5 1]),10)
    t(1) = true;
catch
    t(1) = false;
end
%$
if t(1)
    t(2) = dassert(m1,0.02440738691918476,1e-6);
end
T = all(t);
%@eof:1
