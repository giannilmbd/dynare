function [y, layers] = ms_korder_simul(ms, s_path, u_path, pruning, y0dev, sigma, engine)
% [y, layers] = ms_korder_simul(ms, s_path, u_path, pruning[, y0dev, sigma, engine])
% Simulates the order-k Markov-switching policy returned by ms_korder_solve along a given regime path.
%
% INPUTS  s_path [1 x T]   regime s_t in 1..S of periods t = 1..T (ms_markov_chain draws one)
%         u_path [nu x T]  innovations u_t
%         pruning  logical: pruned state space (layer rule of Dynare's simulate_pruning, libkordersim/simulation.f08,
%                  with regime-dependent coefficients) or plain iteration of the policy
%         y0dev [n_state x 1] deviation of the state variables (ms.state_var order) at t = 0 (default 0); in the pruned
%                  simulation it is put in layer 1
%         sigma    perturbation parameter (default 1)
%         engine   'auto' (default): the MEX ms_k_order_simul (Fortran, layer rule of libkordersim) for a pruned
%                  simulation at sigma = 1 when only y is requested, this file's loop otherwise; 'matlab' forces the loop
% Policy of regime s:  y_t - ybar = sum_{i+j+k<=K} g_{y^i u^j s^k}(s) kron(y*_{t-1}^{(x)i}, u_t^{(x)j}) sigma^k / (i! j! k!).
% Pruned:  layer q of y_t collects the terms kron(x^(a_1),...,x^(a_i), u^{(x)j}) sigma^k with a_1+...+a_i+j+k = q,
%          x^(a) the layer-a state of t-1, q = 1..K; y_t = ybar + sum_q layer q.
% OUTPUT  y [endo_nbr x T] levels, declaration order; layers [endo_nbr x T x K] (pruned; declaration order), else [].
if nargin < 5 || isempty(y0dev), y0dev = zeros(numel(ms.state_var), 1); end
if nargin < 6 || isempty(sigma), sigma = 1; end
if nargin < 7 || isempty(engine), engine = 'auto'; end
K = ms.order; T = numel(s_path); ny = numel(ms.order_var);
if strcmp(engine, 'auto') && pruning && sigma == 1 && nargout < 2 && exist('ms_k_order_simul', 'file') == 3 ...
        && isfield(ms, 'nstatic')
    ys = ms.ys(ms.order_var); ystart = ys;
    is = ms.nstatic + (1:numel(ms.state_var));
    ystart(is) = ystart(is) + y0dev(:);
    sim = ms_k_order_simul(K, ms.nstatic, ms.npred, ms.nboth, ms.nfwrd, size(u_path, 1), ystart, u_path, ys, ...
                           ms_korder_pruning_tensors(ms), double(s_path(:)'));
    inv_ov(ms.order_var) = 1:ny;
    y = sim(inv_ov, 2:end);
    return
end
nstate = numel(ms.state_var);
srow = find(ismember(ms.order_var, ms.state_var));            % rows of the state variables in dr order
[~, loc] = ismember(ms.state_var, ms.order_var(srow)); srow = srow(loc);
inv_ov(ms.order_var) = 1:ny;                                  % dr row -> declaration index mapping
terms = {};                                                    % {i, j, k, field name, 1/(i!j!k!)}
for d = 1:K
    for i = d:-1:0
        for j = (d - i):-1:0
            k = d - i - j;
            terms(end+1, :) = {i, j, k, sprintf('y%du%ds%d', i, j, k), 1/(factorial(i)*factorial(j)*factorial(k))}; %#ok<AGROW>
        end
    end
end
comps = cell(K, K);                                           % comps{i, n}: compositions of n into i parts >= 1
for i = 1:K
    for n = i:K
        comps{i, n} = compositions(n, i);
    end
end
dev = zeros(ny, T);
if pruning
    L = zeros(ny, T, K);
    xl = zeros(nstate, K); xl(:, 1) = y0dev(:);              % state layers at t-1
    for t = 1:T
        g = ms.g{s_path(t)}; u = u_path(:, t);
        for q = 1:numel(terms(:, 1))
            [i, j, k, nm, c] = terms{q, :};
            Tm = g.(nm);
            if ~any(Tm(:)), continue, end
            uk = 1; for a = 1:j, uk = kron(uk, u); end
            uk = uk * sigma^k * c;
            if i == 0
                L(:, t, j + k) = L(:, t, j + k) + Tm * uk;
            else
                for lay = (i + j + k):K                        % layer receiving the term
                    C = comps{i, lay - j - k};
                    v = zeros(nstate^i, 1);
                    for r = 1:size(C, 1)
                        w = xl(:, C(r, 1));
                        for a = 2:i, w = kron(w, xl(:, C(r, a))); end
                        v = v + w;
                    end
                    L(:, t, lay) = L(:, t, lay) + Tm * kron(v, uk);
                end
            end
        end
        xl = squeeze(L(srow, t, :)); xl = reshape(xl, nstate, K);
    end
    dev = sum(L, 3);
    layers = L(inv_ov, :, :);
else
    x = y0dev(:);
    for t = 1:T
        g = ms.g{s_path(t)}; u = u_path(:, t);
        for q = 1:numel(terms(:, 1))
            [i, j, k, nm, c] = terms{q, :};
            v = 1;
            for a = 1:i, v = kron(v, x); end
            for a = 1:j, v = kron(v, u); end
            dev(:, t) = dev(:, t) + g.(nm) * v * sigma^k * c;
        end
        x = dev(srow, t);
    end
    layers = [];
end
y = ms.ys(:) + dev(inv_ov, :);
end

function C = compositions(n, i)
% all ordered i-tuples of integers >= 1 summing to n (one per row)
if i == 1, C = n; return, end
C = zeros(0, i);
for a = 1:(n - i + 1)
    R = compositions(n - a, i - 1);
    C = [C; a * ones(size(R, 1), 1), R]; %#ok<AGROW>
end
end
