function mom = ms_korder_moments(ms, Sigma_u, second, solver)
% mom = ms_korder_moments(ms, Sigma_u[, second[, solver]])
% Exact unconditional and regime-conditional moments of the PRUNED order-K Markov-switching solution (the process
% simulated by ms_korder_simul(..., pruning = true)), at sigma = 1, innovations u_t ~ iid N(0, Sigma_u)
% (pass M_.Sigma_e), regime chain stationary with transition ms.P.
%
% Method. x^(a)_t (a = 1..K) are the pruned layers of the n states. By the layer rule (Lombardo and Uhlig 2018,
% Theorem 1) layer q is
%     x^(q)_t = sum_tau G_tau(s_t) (x^(tau_1)_{t-1} (x) ... (x) x^(tau_i)_{t-1} (x) u_t^{(x)j}),
% G_tau = g_{y^i u^j s^k}(s_t)/(i! j! k!) (state rows), summed over the ordered compositions tau with
% tau_1 + ... + tau_i + j + k = q. For a sorted layer tuple b = (b_1 <= ... <= b_m) the unknown block is
%     mu_b(s) = E[ F_b (x^(b_1)_t (x) ... (x) x^(b_m)_t) 1{s_t = s} ],
% F_b the folding of each group of equal layers to its symmetric coordinates (sorted multi-indices; the
% generalisation of vech). With nu_c = mu_c P (one-step predicted moments),
%     mu_b(:, s') = F_b sum_choices ( (x)_f G_{tau_f}(s') ) (nu_{tau}(:, s') (x) E[u^{(x)J}]),
% a sum over the terms chosen for each factor; within a group of equal layers only one arrangement per multiset of
% terms is computed, times its multiplicity, and F_b applied to the group-symmetrised sum (F_b Sym = D^{-1} Dup').
% Lagged blocks have weight sum(b) or less, and weight sum(b) only with more factors or (all factors g_y) the block
% itself, so the system is block triangular: blocks are solved by weight ascending and, within a weight, by number of
% factors descending. The diagonal block is  M - Phi^{[b]} M P = R,  Phi^{[b]} = F_b Phi^{(x)m} Dup_b, Phi = g_y
% (regime free); with P V = V diag(lam) (real eigenvalues),  (I - lam_k Phi^{[b]}) Mt(:, k) = (R V)(:, k),
% M = Mt V^{-1}, each solved by Smith doubling with (Phi^{[b]})^(2^j) = (Phi^(2^j))^{[b]} applied as Kronecker-vector
% products. W = K gives the means, W = 2K the second moments. y_t - ys is a sum of terms of the same form.
%
% INPUTS  ms        output of ms_korder_solve
%         Sigma_u   [nu x nu] innovation covariance (M_.Sigma_e)
%         second    logical, also compute second moments (default true)
%         solver    diagonal-block solver: 'doubling' (default), 'gensylv' (Dynare Schur recursion on the unfolded
%                   block), 'lu' (dense LU of the folded block)
% OUTPUT  mom.pi [S x 1] ergodic regime probabilities
%         mom.mean [ny x 1] unconditional mean of y (levels, declaration order)
%         mom.mean_regime [ny x S] E[y_t | s_t = s]
%         mom.var [ny x ny] unconditional covariance; mom.second_regime{s} = E[(y-ys)(y-ys)' | s_t = s]  (if second)
%         mom.blocks (layer tuples, in solution order), mom.n_unknowns (folded unknowns, all regimes),
%         mom.residual (max over blocks of the relative residual of the diagonal-block equation)
if nargin < 3, second = true; end
if nargin < 4, solver = 'doubling'; end
K = ms.order; S = numel(ms.g); P = ms.P;
ny = numel(ms.order_var); n = numel(ms.state_var); nu = size(Sigma_u, 1);
srow = zeros(n, 1);
for i = 1:n, srow(i) = find(ms.order_var == ms.state_var(i)); end
inv_ov(ms.order_var) = 1:ny;

[V, D] = eig(P'); [~, i1] = min(abs(diag(D) - 1)); pi_ = real(V(:, i1)); pi_ = pi_ / sum(pi_);
[Vp, Lp] = eig(P); lam = diag(Lp);
if ~isreal(lam), error('ms_korder_moments: P has complex eigenvalues'), end
if cond(Vp) > 1e10, error('ms_korder_moments: P is not diagonalisable (cond of eigenvectors %.1e)', cond(Vp)), end
Vpi = inv(Vp);

W = K * (1 + logical(second));
[Tst, Ty] = build_terms(ms, srow, n, nu, K, S, ny);
Phi = ms.g{1}.y1u0s0(srow, :);
for s = 2:S
    if ~isequal(ms.g{s}.y1u0s0(srow, :), Phi), error('ms_korder_moments: g_y differs across regimes'), end
end
Mu = gauss_kron(Sigma_u, nu, W);

% ---- blocks, block-triangular order
B = block_list(K, W);
wt = cellfun(@sum, B); ln = cellfun(@numel, B);
[~, o] = sortrows([wt(:), -ln(:)]); B = B(o);

store = containers.Map('KeyType', 'char', 'ValueType', 'any');
symc = cell(1, W);
getnu = @(gam) lagged(gam, store, pi_, P);
res = 0; nunk = 0;
for b = 1:numel(B)
    beta = B{b}; m = numel(beta);
    ua = unique(beta); pa = arrayfun(@(a) sum(beta == a), ua);
    ic = 1; C = 1;
    for g = 1:numel(ua)
        if isempty(symc{pa(g)}), [symc{pa(g)}{1}, symc{pa(g)}{2}] = sym_index(n, pa(g)); end
        ic = reshape(ic, [], 1) + C * (reshape(symc{pa(g)}{1}, 1, []) - 1); ic = ic(:); C = C * symc{pa(g)}{2};
    end
    [~, rep] = unique(ic); cnt = accumarray(ic, 1, [C, 1]);
    nunk = nunk + C * S;
    if any(arrayfun(@(a) isempty(Tst{a}), ua))                % a layer without terms is identically zero
        store(bkey(beta)) = struct('mu', zeros(C, S), 'ic', ic); continue
    end
    % choices: per group a multiset of terms of that layer
    ch = cell(1, numel(ua)); mlt = cell(1, numel(ua));
    for g = 1:numel(ua)
        ch{g} = multisets(numel(Tst{ua(g)}), pa(g));
        mlt{g} = arrayfun(@(r) multinom(ch{g}(r, :)), (1:size(ch{g}, 1))');
    end
    nch = cellfun(@(c) size(c, 1), ch);
    Tsum = zeros(n^m, S);
    for c = 1:prod(nch)
        sub = cell(1, numel(nch)); [sub{:}] = ind2sub([nch, 1], c);
        terms = []; mul = 1;
        for g = 1:numel(ua)
            terms = [terms, Tst{ua(g)}(ch{g}(sub{g}, :))]; %#ok<AGROW>
            mul = mul * mlt{g}(sub{g});
        end
        if all([terms.dg]) || mod(sum([terms.j]), 2) == 1, continue, end
        Tsum = Tsum + mul * contract(terms, getnu, Mu, n, nu, S);
    end
    R = zeros(C, S);
    for s = 1:S, R(:, s) = accumarray(ic, Tsum(:, s), [C, 1]) ./ cnt; end
    Rt = R * Vp; Mt = zeros(C, S);
    for k = 1:S, Mt(:, k) = diag_solve(lam(k), Phi, Rt(:, k), ic, rep, m, n, solver); end
    M = Mt * Vpi;
    MP = M * P; AM = zeros(C, S);
    for s = 1:S, AM(:, s) = apply_sym(Phi, MP(:, s), ic, rep, m, n); end
    sc = max(abs([M(:); R(:)]));
    if sc > 0, res = max(res, max(max(abs(M - AM - R))) / sc); end
    store(bkey(beta)) = struct('mu', M, 'ic', ic);
end
mom.pi = pi_; mom.blocks = B; mom.n_unknowns = nunk; mom.residual = res;

% ---- moments of y
ys = ms.ys(ms.order_var);
Ey = zeros(ny, S);
for t = 1:numel(Ty)
    if mod(Ty(t).j, 2) == 0, Ey = Ey + contract(Ty(t), getnu, Mu, n, nu, S); end
end
mom.mean_regime = ys(inv_ov) + Ey(inv_ov, :) ./ pi_';
mom.mean = ys(inv_ov) + sum(Ey(inv_ov, :), 2);
if second
    E2 = zeros(ny * ny, S);
    for t1 = 1:numel(Ty)
        for t2 = t1:numel(Ty)
            if mod(Ty(t1).j + Ty(t2).j, 2) == 1, continue, end
            T = contract([Ty(t1), Ty(t2)], getnu, Mu, n, nu, S);
            E2 = E2 + T;
            if t2 > t1, E2 = E2 + reshape(permute(reshape(T, ny, ny, S), [2 1 3]), [], S); end
        end
    end
    E2 = reshape(E2, ny, ny, S);
    m1 = sum(Ey, 2);
    V2 = sum(E2, 3) - m1 * m1';
    mom.var = V2(inv_ov, inv_ov);
    for s = 1:S
        mom.second_regime{s} = E2(inv_ov, inv_ov, s) / pi_(s);
    end
end
end

% ======================================================================================================================
function [Tst, Ty] = build_terms(ms, srow, n, nu, K, S, ny)
% Tst{q}: terms of state layer q; Ty: terms of y - ys (all layers, merged over the sigma power).
% Fields: G [rows x n^i nu^j x S] (columns with the first factor fastest), tau (lag layers), j (u power),
% rd (regime dependent), dg (the diagonal term g_y x^(q)).
mk = @() struct('G', {}, 'tau', {}, 'j', {}, 'rd', {}, 'dg', {});
Tst = repmat({mk()}, K, 1); Ty = mk(); ykeys = {};
for d = 1:K
    for i = d:-1:0
        for j = (d - i):-1:0
            k = d - i - j;
            nm = sprintf('y%du%ds%d', i, j, k);
            rad = [n * ones(1, i), nu * ones(1, j)];
            Gall = zeros(ny, prod(rad), S);
            for s = 1:S
                Gall(:, :, s) = first_fastest(ms.g{s}.(nm), rad) / (factorial(i) * factorial(j) * factorial(k));
            end
            if ~any(Gall(:)), continue, end
            Gst = Gall(srow, :, :);
            if i == 0
                taus = {zeros(1, 0)}; qs = j + k;
            else
                taus = {}; qs = [];
                for q = (i + j + k):K
                    Cq = compositions(q - j - k, i);
                    for r = 1:size(Cq, 1), taus{end+1} = Cq(r, :); qs(end+1) = q; end %#ok<AGROW>
                end
            end
            for z = 1:numel(taus)
                if any(Gst(:))
                    Tst{qs(z)}(end+1) = struct('G', Gst, 'tau', taus{z}, 'j', j, 'rd', regdep(Gst), ...
                                               'dg', i == 1 && j == 0 && k == 0);
                end
                key = sprintf('%d,', [taus{z}, 0, j]);
                t = find(strcmp(ykeys, key));
                if isempty(t)
                    ykeys{end+1} = key; %#ok<AGROW>
                    Ty(end+1) = struct('G', Gall, 'tau', taus{z}, 'j', j, 'rd', false, 'dg', false); %#ok<AGROW>
                else
                    Ty(t).G = Ty(t).G + Gall;
                end
            end
        end
    end
end
for t = 1:numel(Ty), Ty(t).rd = regdep(Ty(t).G); end
end

function r = regdep(G)
r = size(G, 3) > 1 && any(reshape(G ~= G(:, :, 1), [], 1));
end

function G = first_fastest(G, rad)
% Dynare unfolded columns (first factor slowest) -> first factor fastest
p = numel(rad);
if p <= 1, return, end
r = size(G, 1);
G = reshape(permute(reshape(G, [r, fliplr(rad)]), [1, p+1:-1:2]), r, []);
end

function C = compositions(n, i)
if i == 1, C = n; return, end
C = zeros(0, i);
for a = 1:(n - i + 1)
    R = compositions(n - a, i - 1);
    C = [C; a * ones(size(R, 1), 1), R]; %#ok<AGROW>
end
end

function B = block_list(K, W)
% all nondecreasing tuples of layers 1..K with sum <= W (nonempty)
B = {};
stack = {zeros(1, 0)};
while ~isempty(stack)
    t = stack{end}; stack(end) = [];
    lo = 1; if ~isempty(t), lo = t(end); end
    for a = lo:K
        if sum(t) + a <= W, B{end+1} = [t, a]; stack{end+1} = [t, a]; end %#ok<AGROW>
    end
end
end

function k = bkey(b)
k = char(b + 48);
end

function X = lagged(gam, store, pi_, P)
% nu_gam = mu_gam P, unfolded (n^L x S, sorted modes, first fastest); the empty tuple gives pi' P = pi'
if isempty(gam), X = pi_'; return, end
e = store(bkey(gam));
X = e.mu * P; X = X(e.ic, :);
end

function [ic, C] = sym_index(n, p)
% unfolded multi-index (first factor fastest) -> symmetric (sorted multi-index) coordinate
idx = cell(1, p);
[idx{:}] = ndgrid(1:n);
Mi = sort(reshape(cat(p + 1, idx{:}), [], p), 2);
[~, ~, ic] = unique(Mi, 'rows'); C = max(ic);
end

function c = multisets(T, p)
c = nchoosek(1:T + p - 1, p) - (0:p - 1);
if T + p - 1 == p, c = c(:)'; end
end

function m = multinom(row)
m = factorial(numel(row)) / prod(factorial(accumarray(row(:), 1)));
end

function Mu = gauss_kron(Sig, nu, W)
% Mu{J+1} = E[u^{(x)J}], u ~ N(0, Sig) (Isserlis: pair mode 1 with each mode k)
Mu = cell(1, W + 1); Mu{1} = 1;
for J = 1:W
    if mod(J, 2) == 1, Mu{J + 1} = zeros(nu^J, 1); continue, end
    base = kron(Mu{J - 1}, Sig(:));                   % modes (1, 2) Sigma, modes 3..J the J-2 moment
    acc = zeros(nu^J, 1);
    for k = 2:J
        perm = zeros(1, J); perm(1) = 1; perm(k) = 2; perm(setdiff(1:J, [1 k])) = 3:J;
        acc = acc + reshape(permute(reshape(base, [nu * ones(1, J), 1]), [perm, J + 1]), [], 1);
    end
    Mu{J + 1} = acc;
end
end

function T = contract(terms, getnu, Mu, n, nu, S)
% E[ ((x)_f G_f(s_t) (z_f (x) u_t^{(x)j_f})) 1{s_t = s'} ] as [(prod rows_f) x S], output factor 1 fastest
m = numel(terms);
Lam = [terms.tau]; L = numel(Lam); J = sum([terms.j]);
[gam, ord] = sort(Lam); invord(ord) = 1:L;
X = getnu(gam);                                       % n^L x S
X = X(:) * Mu{J + 1}(:).';                            % modes [gam, regime, u]
perm = []; lc = 0; uc = 0;
for f = 1:m
    li = numel(terms(f).tau);
    perm = [perm, invord(lc + 1:lc + li), L + 1 + (uc + 1:uc + terms(f).j)]; %#ok<AGROW>
    lc = lc + li; uc = uc + terms(f).j;
end
perm = [perm, L + 1];
X = permute(reshape(X, [n * ones(1, L), S, nu * ones(1, J), 1]), [perm, L + J + 2]);
X = X(:);                                             % modes [f1 lags, f1 u, f2 ..., regime]
d = arrayfun(@(t) size(t.G, 2), terms); r = arrayfun(@(t) size(t.G, 1), terms);
for f = 1:m
    G = terms(f).G;
    if ~terms(f).rd
        Y = G(:, :, 1) * reshape(X, d(f), []);
    else
        A = prod(d(f + 1:m)); Bm = prod(r(1:f - 1));
        Xr = reshape(X, d(f), A, S, Bm); Y = zeros(r(f), A, S, Bm);
        for s = 1:S
            Y(:, :, s, :) = reshape(G(:, :, s) * reshape(Xr(:, :, s, :), d(f), []), r(f), A, 1, Bm);
        end
        Y = reshape(Y, r(f), []);
    end
    X = reshape(Y.', [], 1);                          % contracted block -> output mode last
end
T = reshape(X, S, []).';
end

function y = apply_sym(F, x, ic, rep, m, n)
% F_b (F (x) ... (x) F) Dup_b x  (m factors, cyclic mode products)
y = x(ic);
for a = 1:m
    y = reshape(F * reshape(y, n, []), n, []).'; y = y(:);
end
y = y(rep);
end

function x = diag_solve(lambda, Phi, r, ic, rep, m, n, solver)
% (I - lambda Phi^{[b]}) x = r
switch solver
    case 'doubling'
        x = r; l = lambda; F = Phi;
        for it = 1:60
            dx = l * apply_sym(F, x, ic, rep, m, n); x = x + dx;
            if norm(dx) <= 1e-17 * norm(x), return, end
            l = l^2; F = F * F;
        end
        error('ms_korder_moments: doubling did not converge (spectral radius of lambda Phi^{(x)%d} >= 1?)', m);
    case 'gensylv'
        ru = r(ic);
        if m == 1
            xu = (eye(n) - lambda * Phi) \ ru;
        else
            xu = gensylv(m - 1, eye(n), -lambda * Phi, Phi', reshape(ru, n, n^(m - 1)));
        end
        x = xu(rep);
    case 'lu'
        C = numel(r); A = zeros(C);
        for c = 1:C
            e = zeros(C, 1); e(c) = 1; A(:, c) = apply_sym(Phi, e, ic, rep, m, n);
        end
        x = (eye(C) - lambda * A) \ r;
    otherwise
        error('ms_korder_moments: unknown solver %s', solver);
end
end
