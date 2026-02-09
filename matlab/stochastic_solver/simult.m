function [y_out,exo_simul] =simult(y0, dr,M_,options_)
% Simulate a DSGE model (perturbation approach).

%@info:
%! @deftypefn {Function File} {[@var{y_}, @var{exo_simul}] =} simult (@var{y0},@var{dr},@var{M_},@var{options_})
%! @anchor{simult}
%! @sp 1
%! Simulate a DSGE model (perturbation approach).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item y0
%! Vector of doubles, initial conditions.
%! @item dr
%! MATLAB's structure describing decision and transition rules.
%! @item M_
%! MATLAB's structure describing the model
%! @item options_
%! MATLAB's structure describing the current options
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item y_out
%! Matrix of doubles, simulated time series for all the endogenous variables (one per row).
%! @item exo_simul
%! Matrix of doubles, random exogenous shocks
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{non_linear_dsge_likelihood}, @ref{pea/pea_initialization}, @ref{stoch_simul}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! @ref{simult_}
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! If the routine is called with only one output argument, then field exo_simul (structural innovations) is not updated.
%! @end deftypefn
%@eod:

% Copyright © 2001-2025 Dynare Team
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

order = options_.order;
replic = options_.simul_replic;

if replic > 1
    if ~isfolder([M_.dname '/Output'])
        mkdir(M_.dname,'Output');
    end
    fname = [M_.dname filesep 'Output' filesep M_.fname,'_simul'];
    fh = fopen(fname,'w+');
end

% eliminate shocks with 0 variance
i_exo_var = setdiff(1:M_.exo_nbr,find(diag(M_.Sigma_e) == 0));
nxs = length(i_exo_var);
exo_simul = zeros(options_.periods,M_.exo_nbr);
chol_S = chol(M_.Sigma_e(i_exo_var,i_exo_var));

for i=1:replic
    if ~isempty(M_.Sigma_e)
        % we fill the shocks row wise to have the same values independently of the length of the simulation
        if size(M_.Skew_e, 1) > 0 % draw from skew normal distribution (special case of closed skew normal, see csn_update_specification.m for details)
            exo_simul(:,i_exo_var) = transpose(rand_multivariate_csn(options_.periods, M_.csn.mu_e(i_exo_var,1), M_.csn.Sigma_e(i_exo_var,i_exo_var), M_.csn.Gamma_e(i_exo_var,i_exo_var), M_.csn.nu_e(i_exo_var,1), M_.csn.Delta_e(i_exo_var,i_exo_var)));
        else % draw from Gaussian distribution
            exo_simul(:,i_exo_var) = transpose(randn(nxs,options_.periods))*chol_S;
        end
    end
    y_ = simult_(M_,options_,y0,dr,exo_simul,order);
    % eliminating initial value
    y_ = y_(:,2:end);
    if replic > 1
        fwrite(fh,y_,'float64');
    end
    if i==1
        y_out=y_;
    end
end

if replic > 1
    fclose(fh);
end