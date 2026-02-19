function pfm = setup_stochastic_perfect_foresight_model_solver(M_,options_,oo_)

% Sets up the structure for the stochastic perfect foresight model solver.
%
%
% INPUTS
% - M_                     [struct]    Dynare's model structure
% - options_               [struct]    Dynare's options structure
% - oo_                    [struct]    Dynare's results structure
%
% OUTPUTS
% - pfm                    [struct]    Structure containing model information for the perfect foresight solver:
%                                        * lead_lag_incidence: matrix indicating which variables appear with leads/lags
%                                        * ny: number of endogenous variables
%                                        * number_of_shocks: total number of structural shocks
%                                        * positive_var_indx: indices of shocks with positive variance
%                                        * effective_number_of_shocks: number of shocks with positive variance
%                                        * Sigma: covariance matrix of shocks with positive variance
%                                        * Omega: upper Cholesky factor of Sigma (Sigma = Omega'*Omega)
%                                        * stochastic_order: order of stochastic extended path approximation
%                                        * max_lag: maximum lag in the model
%                                        * nyp: number of predetermined (lagged) variables
%                                        * iyp: indices of predetermined variables
%                                        * ny0: number of contemporaneous variables
%                                        * iy0: indices of contemporaneous variables
%                                        * nyf: number of forward-looking (lead) variables
%                                        * iyf: indices of forward-looking variables
%                                        * periods: number of periods for extended path simulation
%                                        * steady_state: steady state values of endogenous variables
%                                        * params: model parameters
%                                        * i_cols_A1: column indices for first block of Jacobian
%                                        * i_cols_1: shifted column indices for first period
%                                        * i_cols_T: column indices for terminal period
%                                        * i_cols_j: column indices for generic period
%                                        * i_upd: indices of variables to update in simulation
%                                        * dynamic_resid: function handle for dynamic residuals (if not bytecode)
%                                        * dynamic_g1: function handle for dynamic Jacobian (if not bytecode)
%                                        * sparse_rowval: row indices for sparse Jacobian (if not bytecode)
%                                        * sparse_colval: column indices for sparse Jacobian (if not bytecode)
%                                        * sparse_colptr: column pointers for sparse Jacobian (if not bytecode)
%                                        * verbose: verbosity level
%                                        * maxit_: maximum iterations for solver
%                                        * tolerance: convergence tolerance
%
% ALGORITHM
% None.
%
% SPECIAL REQUIREMENTS
% None.
%
% OUTPUTS:
%  o  pfm              [struct]    perfect foresight model description
%
% Called by extended_path_initialization.m


% Copyright © 2013-2026 Dynare Team
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

pfm.lead_lag_incidence = M_.lead_lag_incidence;
pfm.ny = M_.endo_nbr;
pfm.number_of_shocks = M_.exo_nbr;
pfm.positive_var_indx = find(diag(M_.Sigma_e)>0);
pfm.effective_number_of_shocks = length(pfm.positive_var_indx);
if pfm.effective_number_of_shocks>0
    pfm.Sigma = M_.Sigma_e(pfm.positive_var_indx,pfm.positive_var_indx);
    pfm.Omega = chol(pfm.Sigma,'upper'); % Sigma = Omega'*Omega
else
    pfm.Sigma = [];
    pfm.Omega = [];
end
pfm.stochastic_order = options_.ep.stochastic.order;
pfm.max_lag = M_.maximum_endo_lag;
if pfm.max_lag > 0
    pfm.nyp = nnz(pfm.lead_lag_incidence(1,:));
    pfm.iyp = find(pfm.lead_lag_incidence(1,:)>0);
else
    pfm.nyp = 0;
    pfm.iyp = [];
end
pfm.ny0 = nnz(pfm.lead_lag_incidence(pfm.max_lag+1,:));
pfm.iy0 = find(pfm.lead_lag_incidence(pfm.max_lag+1,:)>0);
if M_.maximum_endo_lead
    pfm.nyf = nnz(pfm.lead_lag_incidence(pfm.max_lag+2,:));
    pfm.iyf = find(pfm.lead_lag_incidence(pfm.max_lag+2,:)>0);
else
    pfm.nyf = 0;
    pfm.iyf = [];
end
pfm.periods = options_.ep.periods;
pfm.steady_state = oo_.steady_state;
pfm.params = M_.params;
if M_.maximum_endo_lead
    pfm.i_cols_A1 = find(pfm.lead_lag_incidence(pfm.max_lag+(1:2),:)');
else
    pfm.i_cols_A1 = find(pfm.lead_lag_incidence(pfm.max_lag+1,:)');
end
pfm.i_cols_1 = pfm.i_cols_A1 + pfm.max_lag*pfm.ny;
if pfm.max_lag > 0
    pfm.i_cols_T = find(pfm.lead_lag_incidence(1:2,:)');
else
    pfm.i_cols_T = find(pfm.lead_lag_incidence(1,:)') + pfm.ny;
end
pfm.i_cols_j = find(pfm.lead_lag_incidence');
pfm.i_upd = pfm.ny+(1:pfm.periods*pfm.ny);
if ~options_.bytecode
    pfm.dynamic_resid = str2func([M_.fname, '.dynamic_resid']);
    pfm.dynamic_g1 = str2func([M_.fname, '.dynamic_g1']);
    pfm.sparse_rowval = M_.dynamic_g1_sparse_rowval;
    pfm.sparse_colval = M_.dynamic_g1_sparse_colval;
    pfm.sparse_colptr = M_.dynamic_g1_sparse_colptr;
end
pfm.verbose = options_.ep.verbosity;
pfm.maxit_ = options_.simul.maxit;
pfm.tolerance = options_.dynatol.f;
