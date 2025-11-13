function [xparam1, icheck] = set_init_state(xparam1, options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% [xparam1, icheck] = set_init_state(xparam1, options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state)
% Computes the endogenous log prior addition to the initial prior
%
% INPUTS
% - xparam1             [double]        n vector of estimated params
% - options_            [structure]     Matlab's structure describing the current options
% - M_                  [structure]     Matlab's structure describing the model
% - estim_params_       [structure]     characterizing parameters to be estimated
% - bayestopt_          [structure]     describing the priors
% - BoundsInfo          [structure]     containing prior bounds
% - dr                  [structure]     Reduced form model.
% - endo_steady_state   [vector]        steady state value for endogenous variables
% - exo_steady_state    [vector]        steady state value for exogenous variables
% - exo_det_steady_state [vector]       steady state value for exogenous deterministic variables
%
% OUTPUTS
%    xparam1            [double]     n vector of estimated params
%    icheck             [logical]    flag for the need to update xparam1

% Copyright © 2024-2025 Dynare Team
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

icheck=false;

filter_initial_state=M_.filter_initial_state;
M_.filter_initial_state=[];
options_.lik_init=1;
options_.init_state_endogenous_prior=false;
M_ = set_all_parameters(xparam1,estim_params_,M_);
[Pstar, info] = get_pstar(xparam1,options_,M_,estim_params_,bayestopt_,BoundsInfo,dr, endo_steady_state, exo_steady_state, exo_det_steady_state);
if info(1)
    return %JP: error code is neither returned nor handled
end
dr.ys = evaluate_steady_state(endo_steady_state,[exo_steady_state; exo_det_steady_state],M_,options_,true);
M_.filter_initial_state=filter_initial_state;

[UP,XP] = svd(0.5*(Pstar(bayestopt_.mf0,bayestopt_.mf0)+Pstar(bayestopt_.mf0,bayestopt_.mf0)'));
isn = find(diag(XP)<=options_.kalman_tol);
UPN = UP(:,isn);
[~,im]=max(abs(UPN));
if length(unique(im))<length(im)
    in=im*nan;
    for k=1:length(isn)
        [~, tmp] = max(abs(UPN(:,k)));
        if not(ismember(tmp,in))
            in(k) = tmp;
        else
            [~, tmp] = sort(abs(UPN(:,k)));
            not_found = true;
            offset=0;
            while not_found
                offset=offset+1;
                if not(ismember(tmp(end-offset),in))
                    in(k) = tmp(end-offset);
                    not_found = false;
                end
            end
        end
    end
    im=in;
end

a = get_init_state(zeros(M_.endo_nbr,1),xparam1,estim_params_,dr,M_,options_);

mycheck =max(abs(UPN'*a(dr.restrict_var_list(bayestopt_.mf0))));
if mycheck>options_.kalman_tol
    icheck=true;
    nstates = length(M_.state_var);
    % check that state draws are in the null space of prior states 0|0
    in=setdiff(1:nstates,im);
    a1=a;
    a1(dr.restrict_var_list(bayestopt_.mf0(im)))=-UPN(im,:)'\(UPN(in,:)'*a(dr.restrict_var_list(bayestopt_.mf0(in))));
    alphahat01 = a1(dr.inv_order_var); %#ok<NASGU> %declaration order

    % map params associated to init states
    pvec=[];
    for ii=1:size(M_.filter_initial_state,1)
        if ~isempty(M_.filter_initial_state{ii,1})
            tmp1 = strrep(M_.filter_initial_state{ii,2},');','');
            tmp1 = strrep(tmp1,'M_.params(','');
            pvec = [pvec eval(tmp1)];
        end
    end
    [~,~,IB] = intersect(M_.param_names(pvec),bayestopt_.name,'stable');

    M_=update_parameters_filter_initial_state(M_,alphahat01,ys,options_);
    xparam1(IB) = M_.params(pvec);
end
