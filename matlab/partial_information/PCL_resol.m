function [dr,info]=PCL_resol(M_, options_, oo_)
% function [dr,info]=PCL_resol(M_, options_, oo_)
% Computes first and second order approximations
%
% INPUTS
%    ys:             vector of variables in steady state
%    check_flag=0:   all the approximation is computed
%    check_flag=1:   computes only the eigenvalues
%
% OUTPUTS
%    dr:             structure of decision rules for stochastic simulations
%    info=1:         the model doesn't determine the current variables '...' uniquely
%
% SPECIAL REQUIREMENTS
%    none

% Copyright © 2001-2024 Dynare Team
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

if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

% check if ys is steady state
tempex = oo_.exo_simul;
oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
if M_.exo_det_nbr > 0
    tempexdet = oo_.exo_det_simul;
    oo_.exo_det_simul = repmat(oo_.exo_det_steady_state',M_.maximum_lag+M_.maximum_lead+1,1);
end

[dr.ys,M_.params,info] = evaluate_steady_state(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_,options_,~options_.steadystate.nocheck);
if info(1)
    return
end

[dr,info,oo_] = dr1_PI(dr,M_,options_,oo_);
if info(1)
    return
end

if M_.exo_det_nbr > 0
    oo_.exo_det_simul = tempexdet;
end
oo_.exo_simul = tempex;