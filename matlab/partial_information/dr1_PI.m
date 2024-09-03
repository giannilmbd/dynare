function [dr,info,oo_] = dr1_PI(dr,M_,options_,oo_)
% function [dr,info,oo_] = dr1_PI(dr,M_,options_,oo_)
% Computes the reduced form solution of a rational expectation model first
% order
% approximation of the Partial Information stochastic model solver around the deterministic steady state).
% Prepares System as
%        A0*E_t[y(t+1])+A1*y(t)=A2*y(t-1)+c+psi*eps(t)
% with z an exogenous variable process.
% and calls PI_Gensys.m solver
% based on Pearlman et al 1986 paper and derived from
% C.Sims' gensys linear solver.
% to return solution in format
%       [s(t)' x(t)' E_t x(t+1)']'=G1pi [s(t-1)' x(t-1)' x(t)]'+C+impact*eps(t),
%
% INPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   M_         [matlab structure] Definition of the model.
%   options_   [matlab structure] Global options.
%   oo_        [matlab structure] Results
%
% OUTPUTS
%   dr         [matlab structure] Decision rules for stochastic simulations.
%   info       [integer]          info=1: the model doesn't define current variables uniquely
%                                 info=2: problem in mjdgges.dll info(2) contains error code.
%                                 info=3: BK order condition not satisfied info(2) contains "distance"
%                                         absence of stable trajectory.
%                                 info=4: BK order condition not satisfied info(2) contains "distance"
%                                         indeterminacy.
%                                 info=5: BK rank condition not satisfied.
%   oo_        [matlab structure]
%
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   none.
%

% Copyright © 1996-2024 Dynare Team
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

info = 0;

options_ = set_default_option(options_,'qz_criterium',1.000001);

xlen = M_.maximum_endo_lead + M_.maximum_endo_lag + 1;

if options_.aim_solver
    error('Anderson and Moore AIM solver is not compatible with Partial Information models');
end % end if useAIM and...

if (options_.order > 1)
    warning('You can not use order higher than 1 with Part Info solver. Order 1 assumed');
    options_.order =1;
end

klen = M_.maximum_lag + M_.maximum_lead + 1;
iyv = M_.lead_lag_incidence';

if M_.exo_nbr == 0
    oo_.exo_steady_state = [] ;
end


z = repmat(dr.ys,1,klen);
z = z(find(iyv(:))) ;
[~,jacobia_] = feval([M_.fname '.dynamic'],z,[oo_.exo_simul ...
    oo_.exo_det_simul], M_.params, dr.ys, M_.maximum_lag + 1);

if options_.debug
    save([M_.dname filesep 'Output' filesep M_.fname '_debug.mat'],'jacobia_')
end

% find size xlen of the state vector Y and of A0, A1 and A2 transition matrices:
% it is the sum the all i variables's lag/lead representations,
% for each variable i representation being defined as:
% Max (i_lags-1,0)+ Max (i_leads-1,0)+1
% so that if variable x appears with 2 lags and 1 lead, and z
% with 2 lags and 3 leads, the size of the state space is:
% 1+0+1   +   1+2+1   =6
% e.g. E_t Y(t+1)=
%     E_t x(t)
%     E_t x(t+1)
%     E_t z(t)
%     E_t z(t+1)
%     E_t z(t+2)
%     E_t z(t+3)

% partition jacobian:
jlen=M_.nspred+M_.nsfwrd+M_.endo_nbr+M_.exo_nbr; % length of jacobian
% first transpose M_.lead_lag_incidence';
lead_lag=M_.lead_lag_incidence';
if ( M_.maximum_lag <= 1) && (M_.maximum_lead <= 1)
    xlen=size(jacobia_,1);%nendo;
    AA0=zeros(xlen,xlen);  % empty A0
    AA2=AA0; % empty A2 and A3
    AA3=AA0;
    if xlen==M_.endo_nbr % && M_.maximum_lag <=1 && M_.maximum_lead <=1 % apply a shortcut
        AA1=jacobia_(:,M_.nspred+1:M_.nspred+M_.endo_nbr);
        if M_.maximum_lead ==1
            fnd = find(lead_lag(:,M_.maximum_lag+2));
            AA0(:, fnd)= jacobia_(:,nonzeros(lead_lag(:,M_.maximum_lag+2))); %forwd jacobian
        end
        if M_.nspred>0 && M_.maximum_lag ==1
            fnd = find(lead_lag(:,1));
            AA2(:, fnd)= jacobia_(:,nonzeros(lead_lag(:,1))); %backward
        end
    else
        error('More than one lead or lag in the jabian');
    end
    if M_.orig_endo_nbr<M_.endo_nbr
        % findif there are any expectations at time t
        exp_0= strmatch('AUX_EXPECT_LEAD_0_', M_.endo_names);
        if ~isempty(exp_0) 
    	    error('Partial Information does not support EXPECTATION(0) operators')
        end
    end
end
PSI=-[[zeros(xlen-M_.endo_nbr,M_.exo_nbr)];[jacobia_(:, jlen-M_.exo_nbr+1:end)]]; % exog
cc=0;

try
    % call [G1pi,C,impact,nmat,TT1,TT2,gev,eu]=PI_gensys(a0,a1,a2,c,PSI,NX,NETA,NO_FL_EQS)
    % System given as
    %        a0*E_t[y(t+1])+a1*y(t)=a2*y(t-1)+c+psi*eps(t)
    % with eps an exogenous variable process.
    % Returned system is
    %       [s(t)' x(t)' E_t x(t+1)']'=G1pi [s(t-1)' x(t-1)' x(t)]'+C+impact*eps(t),
    %  and (a) the matrix nmat satisfying   nmat*E_t z(t)+ E_t x(t+1)=0
    %      (b) matrices TT1, TT2  that relate y(t) to these states:
    %      y(t)=[TT1 TT2][s(t)' x(t)']'.

    %% for expectational models when complete
    if any(AA3)
        error('Unsupported case')
    else
        [G1pi,CC,impact,nmat,TT1,TT2,gev,eu, DD, E2,E5, GAMMA, FL_RANK]=PI_gensys(AA0,AA1,-AA2,cc,PSI,M_.exo_nbr);
    end

    % reuse some of the bypassed code and tests that may be needed
    if (eu(1) ~= 1 || eu(2) ~= 1)
        info(1) = abs(eu(1)+eu(2));
        info(2) = 1.0e+8;
        %     return
    end

    dr.PI_ghx=G1pi;
    dr.PI_ghu=impact;
    dr.PI_TT1=TT1;
    dr.PI_TT2=TT2;
    dr.PI_nmat=nmat;
    dr.PI_CC=CC;
    dr.PI_gev=gev;
    dr.PI_eu=eu;
    dr.PI_FL_RANK=FL_RANK;
    dr.ghx=G1pi;
    dr.ghu=impact;
    dr.eigval = eig(G1pi);
    dr.rank=FL_RANK;
    
catch ME
    disp('Problem with using Part Info solver');
    rethrow(ME);
end