function [nodes, weights, nnodes] = setup_integration_nodes(opt, pfm)

% INPUTS:
% - opt         [struct]   EP options
% - pfm         [struct]   perfect foresight model description
%
% OUTPUTS:
% - nodes       [double]    vector of integration nodes
% - weights     [double]    vector of weights
% - nnodes      [integer]   number of nodes

% Copyright © 2012-2026 Dynare Team
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

nodes=[];
weights=[];
nnodes=[];

if opt.stochastic.order
    % Compute weights and nodes for the stochastic version of the extended path.
    % Use only shocks with positive variance (effective shocks).
    n = pfm.effective_number_of_shocks;
    switch opt.stochastic.IntegrationAlgorithm
      case 'Tensor-Gaussian-Quadrature'
        % Get the nodes and weights from a univariate Gauss-Hermite quadrature.
        [nodes0, weights0] = gauss_hermite_weights_and_nodes(opt.stochastic.quadrature.nodes);
        % Replicate the univariate nodes for each innovation and, if needed, correlate them.
        nodes0 = repmat(nodes0, 1, n)*pfm.Omega;
        % Put the nodes and weights in cells
        for i=1:n
            rr(i) = {nodes0(:,i)};
            ww(i) = {weights0};
        end
        % Build the tensorial grid
        nodes = cartesian_product_of_sets(rr{:});
        weights = prod(cartesian_product_of_sets(ww{:}),2);
        nnodes = length(weights);
      case 'Stroud-Cubature-3'
        [nodes,weights] = cubature_with_gaussian_weight(n, 3, 'Stroud');
        nodes = transpose(pfm.Omega'*nodes);
        weights = weights/sum(weights);
        nnodes = length(weights);
      case 'Stroud-Cubature-5'
        if n==1
            info = warning('query', 'backtrace');
            if strcmp(info.state, 'on')
                warning('off', 'backtrace');
            end
            warning('Stroud-Cubature-5 is not defined for a single shock, falling back to Gaussian quadrature with 3 nodes.')
            skipline()
            warning(info.state, 'backtrace');
            [nodes, weights] = gauss_hermite_weights_and_nodes(3);
            nodes = nodes*pfm.Omega;
            nnodes = length(weights);
        else
            [nodes,weights] = cubature_with_gaussian_weight(n,5,'Stroud');
            nodes = transpose(pfm.Omega'*nodes);
            weights = weights/sum(weights);
            nnodes = length(weights);
        end
      case 'Unscented'
        k = 3;
        C = sqrt(n + k)*pfm.Omega';
        nodes = [zeros(1,n); -C; C];
        weights = [k/(n+k); (1/(2*(n+k)))*ones(2*n,1)];
        nnodes = 2*n+1;
      otherwise
        error('Stochastic extended path:: Unknown integration algorithm %s!',opt.stochastic.IntegrationAlgorithm)
    end
end
