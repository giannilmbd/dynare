function [errorflag, endo_simul, errorcode, y, pfm, options_] = solve_stochastic_perfect_foresight_model_0(endo_simul, exo_simul, y, options_, M_, pfm)

% Copyright © 2025 Dynare Team
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

dynamic_g1 = pfm.dynamic_g1;

periods = pfm.periods;

order = pfm.stochastic_order;

lead_lag_incidence = pfm.lead_lag_incidence;
lead_lag_incidence_t = transpose(lead_lag_incidence);
ny = pfm.ny;
nyp = pfm.nyp;
nyf = pfm.nyf;
i_cols_1 = pfm.i_cols_1;
i_cols_j = pfm.i_cols_j;
i_cols_T = pfm.i_cols_T;

nodes = pfm.nodes;
weights = pfm.weights;
nnodes = pfm.nnodes;

if update_pfm_struct

    % make sure that there is a node equal to zero
    % and permute nodes and weights to have zero first
    k = find(sum(abs(nodes),2) < 1e-12);
    if ~isempty(k)
        nodes = [nodes(k,:); nodes(1:k-1,:); nodes(k+1:end,:)];
        weights = [weights(k); weights(1:k-1); weights(k+1:end)];
    else
        error('there is no nodes equal to zero')
    end

    pfm.nodes = nodes;
    pfm.weights = weights;

    if pfm.hybrid_order>0
        if pfm.hybrid_order==2
            pfm.h_correction = 0.5*pfm.dr.ghs2(pfm.dr.inv_order_var);
        elseif pfm.hybrid_order>2
            pfm.h_correction = pfm.dr.g_0(pfm.dr.inv_order_var);
        else
            pfm.h_correction = 0;
        end
    else
        pfm.h_correction = 0;
    end

    z = endo_simul(1:3*ny);
    jacobian = dynamic_g1(z, exo_simul(2,:), pfm.params, pfm.steady_state, pfm.sparse_rowval, pfm.sparse_colval, pfm.sparse_colptr);

    world_nbr = nnodes^order;

    % The columns of A map the elements of Y such that
    % each block of Y with ny rows are unfolded column wise
    dimension = ny*(sum(nnodes.^(0:order-1),2)+(periods-order)*world_nbr);

    i_upd_r = zeros(dimension,1);
    i_upd_y = i_upd_r;
    i_upd_r(1:ny) = (1:ny);
    i_upd_y(1:ny) = ny+(1:ny);
    i1 = ny+1;
    i2 = 2*ny;
    n1 = ny+1;
    n2 = 2*ny;
    for i=2:periods
        for j=1:nnodes^min(i-1,order)
            i_upd_r(i1:i2) = (n1:n2)+(j-1)*ny*periods;
            i_upd_y(i1:i2) = (n1:n2)+ny+(j-1)*ny*(periods+2);
            i1 = i2+1;
            i2 = i2+ny;
        end
        n1 = n2+1;
        n2 = n2+ny;
    end

    if rows(lead_lag_incidence)>2
        icA = [find(lead_lag_incidence(1,:)) find(lead_lag_incidence(2,:))+world_nbr*ny ...
               find(lead_lag_incidence(3,:))+2*world_nbr*ny]';
    else
        if nyf
            icA = [find(lead_lag_incidence(2,:))+world_nbr*ny find(lead_lag_incidence(3,:))+2*world_nbr*ny ]';
        else
            icA = [find(lead_lag_incidence(1,:)) find(lead_lag_incidence(2,:))+world_nbr*ny ]';
        end
    end

    i_rows = 1:ny;
    i_cols = find(lead_lag_incidence');
    pfm.i_cols_Ap = i_cols(1:nyp);;
    pfm.i_cols_As = i_cols(nyp+(1:ny));
    pfm.i_cols_Af = i_cols(nyp+ny+(1:nyf)) - ny;
    pfm.i_hc = 1:ny;
    pfm.i_cols_p = 1:ny;
    pfm.i_cols_s = ny + (1:ny);
    pfm.i_cols_f = 2*ny + (1:ny);
    pfm.i_rows = i_rows;

    pfm.A1 = sparse([],[],[],ny*(sum(nnodes.^(0:order-1),2)+1),dimension,(order+1)*world_nbr*nnz(jacobian));
    pfm.res = zeros(ny,periods,world_nbr);

    pfm.order = order;
    pfm.world_nbr = world_nbr;

    pfm.i_cols_1 = i_cols_1;
    pfm.i_cols_j = i_cols_j;
    pfm.icA = icA;
    pfm.i_cols_T = i_cols_T;
    pfm.i_upd_r = i_upd_r;
    pfm.i_upd_y = i_upd_y;

    pfm.dimension = dimension;

end

pfm.Y = repmat(endo_simul(:),1,pfm.world_nbr);

if isempty(y)
    y = repmat(pfm.steady_state, pfm.dimension/pfm.ny, 1);
end

if update_options_struct
    % Set algorithm
    options_.solve_algo = options_.ep.solve_algo;
    options_.simul.maxit = options_.ep.maxit;
    [lb, ub] = feval(sprintf('%s.dynamic_complementarity_conditions', M_.fname), pfm.params);
    pfm.eq_index = M_.dynamic_mcp_equations_reordering;
    if options_.ep.solve_algo == 10
        options_.lmmcp.lb = repmat(lb, pfm.dimension/pfm.ny, 1);
        options_.lmmcp.ub = repmat(ub, pfm.dimension/pfm.ny, 1);
    elseif options_.ep.solve_algo == 11
        options_.mcppath.lb = repmat(lb, pfm.dimension/pfm.ny, 1);
        options_.mcppath.ub = repmat(ub, pfm.dimension/pfm.ny, 1);
    end
end

[y, errorflag, ~, ~, errorcode] = dynare_solve(@ep_problem_0, y, options_.simul.maxit, options_.dynatol.f, options_.dynatol.x, options_, exo_simul, pfm);

endo_simul(:,2) = y(1:ny);
