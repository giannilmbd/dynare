function [nodes,weights,nnodes] = setup_integration_nodes(EpOptions,pfm)
% [nodes,weights,nnodes] = setup_integration_nodes(EpOptions,pfm)
% Inputs:
%   ep          [structure] EP options
%   pfm         [structure] perfect foresight model description
%
% Ouputs:
%   nodes       [double]    integration nodes
%   weights     [double]    weights of nodes
%   nnodes      [integer]   number of nodes

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

if EpOptions.stochastic.order
    % Compute weights and nodes for the stochastic version of the extended path.
    % Use only shocks with positive variance (effective shocks).
    n = pfm.effective_number_of_shocks;
    switch EpOptions.stochastic.IntegrationAlgorithm
      case 'Tensor-Gaussian-Quadrature'
        % Get the nodes and weights from a univariate Gauss-Hermite quadrature.
        [nodes0,weights0] = gauss_hermite_weights_and_nodes(EpOptions.stochastic.quadrature.nodes);
        % Replicate the univariate nodes for each innovation and dates, and, if needed, correlate them.
        nodes0 = repmat(nodes0,1,n*pfm.stochastic_order)*kron(eye(pfm.stochastic_order),pfm.Omega);
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
        [nodes,weights] = cubature_with_gaussian_weight(n*pfm.stochastic_order,3,'Stroud');
        nodes = kron(eye(pfm.stochastic_order),transpose(pfm.Omega))*nodes;
        nnodes = length(weights);
      case 'Stroud-Cubature-5'
        [nodes,weights] = cubature_with_gaussian_weight(n*pfm.stochastic_order,5,'Stroud');
        nodes = kron(eye(pfm.stochastic_order),transpose(pfm.Omega))*nodes;
        nnodes = length(weights);
      case 'Unscented'
        k = 3;%EpOptions.ut.k;
        C = sqrt(n + k)*pfm.Omega';
        nodes = [zeros(1,n); -C; C];
        weights = [k/(n+k); (1/(2*(n+k)))*ones(2*n,1)];
        nnodes = 2*n+1;
      otherwise
        error('Stochastic extended path:: Unknown integration algorithm %s!',EpOptions.stochastic.IntegrationAlgorithm)
    end
else
    nodes=[];
    weights=[];
    nnodes=[];
end
