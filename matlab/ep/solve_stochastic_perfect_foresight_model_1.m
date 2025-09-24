function [errorflag, endo_simul, errorcode, y, pfm, options_] = solve_stochastic_perfect_foresight_model_1(endo_simul, exo_simul, y, options_, M_, pfm)

% Copyright © 2012-2025 Dynare Team
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

update_pfm_struct = false;
update_options_struct = false;

if nargout>4
    update_pfm_struct = true;
end

if nargout>5
    update_options_struct = true;
end


if update_pfm_struct

    periods = pfm.periods;

    order = pfm.stochastic_order;

    ny = pfm.ny;
    lead_lag_incidence = pfm.lead_lag_incidence;
    i_cols_1 = pfm.i_cols_1;
    i_cols_j = pfm.i_cols_j;
    i_cols_T = pfm.i_cols_T;

    nodes = pfm.nodes;
    weights = pfm.weights;
    nnodes = pfm.nnodes;

    % Make sure that there is a node equal to zero and permute nodes and weights to have zero first
    k = find(sum(abs(nodes),2) < 1e-12);
    if ~isempty(k)
        nodes = [nodes(k,:); nodes(1:k-1,:); nodes(k+1:end,:)];
        weights = [weights(k); weights(1:k-1); weights(k+1:end)];
    else
        error('there is no nodes equal to zero')
    end

    pfm.nodes = nodes;
    pfm.weights = weights;

    if pfm.hybrid_order > 0
        if pfm.hybrid_order == 2
            pfm.h_correction = 0.5*pfm.dr.ghs2(pfm.dr.inv_order_var);
        elseif pfm.hybrid_order>2
            pfm.h_correction = pfm.dr.g_0(pfm.dr.inv_order_var);
        else
            pfm.h_correction = 0;
        end
    else
        pfm.h_correction = 0;
    end

    % Each column of Y represents a different world
    % The upper right cells are unused
    % The first row block is ny x 1
    % The second row block is ny x nnodes
    % The third row block is ny x nnodes^2
    % and so on until size ny x nnodes^order
    world_nbr = pfm.world_nbr;

    % The columns of A map the elements of Y such that
    % each block of Y with ny rows are unfolded column wise
    % number of blocks
    block_nbr = pfm.block_nbr;
    % dimension of the problem
    dimension = ny*block_nbr;
    pfm.dimension = dimension;
    i_upd_r = zeros(dimension, 1);
    i_upd_y = i_upd_r;
    i_upd_r(1:ny) = (1:ny);
    i_upd_y(1:ny) = ny+(1:ny);
    i1 = ny+1;
    i2 = 2*ny;
    n1 = ny+1;
    n2 = 2*ny;
    for i=2:periods
        k = n1:n2;
        for j=1:(1+(nnodes-1)*min(i-1,order))
            i_upd_r(i1:i2) = k+(j-1)*ny*periods;
            i_upd_y(i1:i2) = k+ny+(j-1)*ny*(periods+2);
            i1 = i2+1;
            i2 = i2+ny;
        end
        n1 = n2+1;
        n2 = n2+ny;
    end
    icA = [find(lead_lag_incidence(1,:)) find(lead_lag_incidence(2,:))+world_nbr*ny ...
           find(lead_lag_incidence(3,:))+2*world_nbr*ny]';

    pfm.i_rows = 1:ny;
    pfm.i_cols = find(lead_lag_incidence');
    pfm.i_cols_1 = i_cols_1;
    pfm.i_cols_j = i_cols_j;
    pfm.icA = icA;
    pfm.i_cols_T = i_cols_T;
    pfm.i_upd_r = i_upd_r;
    pfm.i_upd_y = i_upd_y;

end

if isempty(y)
    y = repmat(pfm.steady_state, pfm.block_nbr, 1);
end

if update_options_struct
    % Set algorithm
    options_.solve_algo = options_.ep.solve_algo;
    options_.simul.maxit = options_.ep.maxit;
    [lb, ub] = feval(sprintf('%s.dynamic_complementarity_conditions', M_.fname), pfm.params);
    pfm.eq_index = M_.dynamic_mcp_equations_reordering;
    if options_.ep.solve_algo == 10
        options_.lmmcp.lb = repmat(lb, pfm.block_nbr, 1);
        options_.lmmcp.ub = repmat(ub, pfm.block_nbr, 1);
    elseif options_.ep.solve_algo == 11
        options_.mcppath.lb = repmat(lb, pfm.block_nbr, 1);
        options_.mcppath.ub = repmat(ub, pfm.block_nbr, 1);
    end
end

pfm.y0 = endo_simul(:,1);

[y, errorflag, ~, ~, errorcode] = dynare_solve(@ep_problem_1, y, options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, options_, exo_simul, pfm);
endo_simul(:,2) = y(1:pfm.ny);
