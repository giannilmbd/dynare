function [y_,int_width,int_width_ME]=simultxdet(y0,ex,ex_det, iorder,var_list,M_,oo_,options_)
%function [y_,int_width,int_width_ME]=simultxdet(y0,ex,ex_det, iorder,var_list,M_,oo_,options_)
%
% Simulates a stochastic model in the presence of deterministic exogenous shocks
%
% INPUTS:
%    y0:        initial values, of length M_.maximum_lag
%    ex:        matrix of stochastic exogenous shocks, starting at period 1
%    ex_det:    matrix of deterministic exogenous shocks, starting at period 1-M_.maximum_lag
%    iorder:    order of approximation
%    var_list:  list of endogenous variables to simulate
%    M_:        Dynare model structure
%    oo_:       Dynare results structure
%    options_:  Dynare options structure
%
% OUTPUTS:
%   yf:          mean forecast (in var_list order)
%   int_width:   distance between upper bound and
%                mean forecast (in var_list order)
%   int_width_ME:distance between upper bound and
%                mean forecast when considering measurement error (in var_list order)
%   int_width_ME:distance between upper bound and
%                mean forecast when considering measurement error (in var_list order)
%
% The forecast horizon is equal to size(ex, 1).
% The condition size(ex,1)+M_.maximum_lag=size(ex_det,1) must be verified
%  for consistency.

% Copyright © 2008-2025 Dynare Team
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

if options_.order>2 
    error('simultxdet.m: Forecasting with varexo_det does not support order>2.')   
elseif options_.order==2 && options_.pruning
    error('simultxdet.m: Forecasting with varexo_det does not support pruning.')
end
dr = oo_.dr;
ykmin = M_.maximum_lag;
endo_nbr = M_.endo_nbr;
exo_det_steady_state = oo_.exo_det_steady_state;
nstatic = M_.nstatic;
nspred = M_.nspred;
iter = size(ex,1);
if size(ex_det, 1) ~= iter+ykmin
    error('Size mismatch: number of forecasting periods for stochastic exogenous and deterministic exogenous don''t match')
end
y_ = zeros(size(y0,1),iter+ykmin);
y_(:,1:ykmin) = y0;
k1 = ykmin:-1:1;
k2 = nstatic+(1:nspred);

nvar = length(var_list);
if nvar == 0
    nvar = endo_nbr;
    ivar = 1:nvar;
else
    ivar=zeros(nvar,1);
    for i=1:nvar
        i_tmp = strmatch(var_list{i}, M_.endo_names, 'exact');
        if isempty(i_tmp)
            disp(var_list{i})
            error ('One of the variable specified does not exist') ;
        else
            ivar(i) = i_tmp;
        end
    end
end
if nargout==3 % add measurement error variance uncertainty for observables
    [loc_H, loc_varlist] = ismember(options_.varobs, var_list);
    loc_varlist(loc_varlist==0) = []; %location of observables in variable list
end

if iorder == 1
    for i = ykmin+1: iter+ykmin
        tempx1 = y_(dr.order_var,k1);
        tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,ykmin);
        tempx = tempx2(k2);
        y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghx*tempx+dr.ghu*ex(i-ykmin,:)';
        for j=1:min(ykmin+dr.exo_det_length+1-i,dr.exo_det_length)
            y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*(ex_det(i+j-1,:)'-exo_det_steady_state);
        end

        k1 = k1+1;
    end
elseif iorder == 2
    for i = ykmin+1: iter+ykmin
        tempx1 = y_(dr.order_var,k1);
        tempx2 = tempx1-repmat(dr.ys(dr.order_var),1,ykmin);
        tempx = tempx2(k2);
        tempu = ex(i-ykmin,:)';
        tempuu = kron(tempu,tempu);
        tempxx = kron(tempx,tempx);
        tempxu = kron(tempx,tempu);
        y_(dr.order_var,i) = dr.ys(dr.order_var)+dr.ghs2/2+dr.ghx*tempx+ ...
            dr.ghu*tempu+0.5*(dr.ghxx*tempxx+dr.ghuu*tempuu)+dr.ghxu* ...
            tempxu; % in declaration order
        for j=1:min(ykmin+dr.exo_det_length+1-i,dr.exo_det_length)
            tempud = ex_det(i+j-1,:)'-exo_det_steady_state;
            tempudud = kron(tempud,tempud);
            tempxud = kron(tempx,tempud);
            tempuud = kron(tempu,tempud);
            y_(dr.order_var,i) = y_(dr.order_var,i) + dr.ghud{j}*tempud + ...
                dr.ghxud{j}*tempxud + dr.ghuud{j}*tempuud + ...
                0.5*dr.ghudud{j,j}*tempudud;
            for k=1:j-1
                tempudk = ex_det(i+k-1,:)'-exo_det_steady_state;
                tempududk = kron(tempudk,tempud);
                y_(dr.order_var,i) = y_(dr.order_var,i) + ...
                    dr.ghudud{k,j}*tempududk;
            end
        end
        k1 = k1+1;
    end
else
    error('simultxdet.m: order>2 not supported.')
end

% Define union of requested and state variables
[var_pos]=union(dr.inv_order_var(dr.state_var),dr.inv_order_var(ivar),'stable');
% set state_pos to positions of state variables for inner part of
% variance loop
[~,state_pos] = ismember(dr.inv_order_var(dr.state_var),var_pos);
% set var_list_pos to positions of var_list_ for reading out requested
% variances
[~,var_list_pos] = ismember(dr.inv_order_var(ivar),var_pos);
ghx1 = dr.ghx(var_pos,:);
ghu1 = dr.ghu(var_pos,:);


sigma_y1 = zeros(size(ghx1,1),size(ghx1,1)); %no initial uncertainty about the states
var_yf=NaN(iter,nvar); %initialize
for i=1:iter
    sigma_y1 = ghx1*sigma_y1(state_pos,state_pos)*ghx1'+ghu1*M_.Sigma_e*ghu1';%only valid at first order, needs to be fixed, see https://git.dynare.org/Dynare/dynare/-/issues/1940
    var_yf(i,:) = diag(sigma_y1(var_list_pos,var_list_pos))'; % first period is initial condition without uncertainty
end

fact = norminv((1-options_.forecasts.conf_sig)/2,0,1);
if nargout==3
    var_yf_ME=var_yf;
    if ~isempty(loc_varlist)
        var_yf_ME(:,loc_varlist)=var_yf(:,loc_varlist)+repmat(diag(M_.H(loc_H,loc_H))',iter,1);
    end
    int_width_ME = zeros(iter,nvar);
end

int_width = zeros(iter,nvar);
for i=1:nvar
    int_width(:,i) = -fact*sqrt(var_yf(:,i));
    if nargout==3
        int_width_ME(:,i) = -fact*sqrt(var_yf_ME(:,i));
    end
end

y_ = y_(ivar,:); %reorder from declaration to var_list order
