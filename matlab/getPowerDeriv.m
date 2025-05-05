function dxp=getPowerDeriv(x,p,k)
%function dxp=getPowerDeriv(x,p,k)
% The k-th derivative of x^p
%
% INPUTS
%    x: base
%    p: power
%    k: derivative order
%
% OUTPUTS
%    dxp: k-th derivative of x^p
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2011-2012 Dynare Team
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

if (abs(x) < 1e-12) && (p >= 0) && (k > p) && (abs(p - round(p)) < 1e-12)
    dxp = 0;
else
    dxp = x^(p-k);
    for i=0:k-1
        dxp = dxp*p;
        p = p-1;
    end
end

return % --*-- Unit tests --*--

%@test:1
x=getPowerDeriv(2,3,1);
t(1)=all(abs(x-3*4)<1e-10)
x=getPowerDeriv(0,2,2);
t(2)=all(abs(x-2)<1e-10)
x=getPowerDeriv(0,2,3); %special case evaluates to 0
t(3)=all(abs(x-0)<1e-10)
x=getPowerDeriv(1e-13,2,3-1e-13); %0 within tolerance
t(4)=all(abs(x-0)<1e-10)
x=getPowerDeriv(0,0,1);
t(5)=all(abs(x-0)<1e-10)
x=getPowerDeriv(0,0,0);
t(6)=all(abs(x-1)<1e-10);
x=getPowerDeriv(0,1/3,1); %derivative evaluating to Inf due to division by 0
t(7)= isinf(x)
T = all(t);
%@eof:1