function drs = ms_korder_pruning_tensors(ms)
% drs = ms_korder_pruning_tensors(ms)
% Decision rule of each regime in the format of dr.pruning of k_order_perturbation (input of ms_k_order_simul):
% drs(s).g_d [ny x nv^d], d = 1..K, the unfolded fully symmetric tensor D_d / d! of the derivatives with respect to
% v = (y* - ybar*, u, sigma), nv = n_state + nu + 1, rows in dr order. D_d at an index with i coordinates in y*,
% j in u and k equal to sigma is g_{y^i u^j s^k}(s).
K = ms.order; S = numel(ms.g); ny = numel(ms.order_var); n = numel(ms.state_var);
nu = size(ms.g{1}.y0u1s0, 2); nv = n + nu + 1;
yi = 1:n; ui = n + (1:nu); si = nv;
drs = repmat(struct(), 1, S);
for s = 1:S
    for d = 1:K
        G = zeros([ny, nv * ones(1, d), 1]);
        pm = perms(1:d);
        for i = d:-1:0
            for j = (d - i):-1:0
                k = d - i - j;
                Tm = ms.g{s}.(sprintf('y%du%ds%d', i, j, k));
                if ~any(Tm(:)), continue, end
                % Dynare unfolded columns kron(y*^i, u^j), first factor slowest: modes (u_j, ..., u_1, y_i, ..., y_1)
                E = zeros([ny, nv * ones(1, d), 1]);
                idx = [{':'}, repmat({ui}, 1, j), repmat({yi}, 1, i), repmat({si}, 1, k)];
                E(idx{:}) = reshape(Tm, [ny, nu * ones(1, j), n * ones(1, i), ones(1, k), 1]);
                Es = zeros(size(E));
                for r = 1:size(pm, 1), Es = Es + permute(E, [1, 1 + pm(r, :), d + 2]); end
                G = G + Es / (factorial(i) * factorial(j) * factorial(k));
            end
        end
        drs(s).(sprintf('g_%d', d)) = reshape(G, ny, nv^d) / factorial(d);
    end
end
end
