function mom = ms_korder_moments(ms, shocks, order, solver, vars)
% mom = ms_korder_moments(ms, shocks[, order[, solver[, vars]]])
% Exact unconditional and regime-conditional moments of the PRUNED order-K Markov-switching solution (the process
% simulated by ms_korder_simul(..., pruning = true)), at sigma = 1, innovations u_t iid with the moments given by
% shocks (see ms_korder_shock_moments), regime chain stationary with transition ms.P.
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
% (regime free). It is solved by Smith doubling on the matrix equation, M <- M + (Phi^(2^j))^{[b]} M P^(2^j), with
% (Phi^{[b]})^(2^j) = (Phi^(2^j))^{[b]} applied as Kronecker-vector products; no decomposition of P is used, so P may
% have complex eigenvalues or be defective. Moments of order p of y need the blocks of weight up to W = pK.
% y_t - ys is a sum of terms of the same form, and E[(y - ys)^{(x)p} 1{s_t = s}] is the sum over multisets of p terms
% of the contraction of one arrangement times its multiplicity, symmetrised over the p factors.
%
% INPUTS  ms        output of ms_korder_solve
%         shocks    M_ (innovations defined by Sigma_e, and Raw_moments_e and Skew_e for non-Gaussian independent
%                   shocks, as in the k-order solver), or the covariance [nu x nu] of Gaussian innovations; moments of
%                   order p need E[e^d] up to d = pK
%         order     highest moment, 1 to 4 (default 2); true/false are read as 2/1
%         solver    diagonal-block solver: 'doubling' (default), 'gensylv' (Dynare Schur recursion on the unfolded
%                   block, regimes stacked: X - (P' (x) Phi) X Phi'^{(x)(m-1)} = R), 'lu' (dense LU of
%                   I - P' (x) Phi^{[b]} on the folded block)
%         vars      variables (declaration indices) of the moments of order 3 and 4 (default: all)
% OUTPUT  mom.pi [S x 1] ergodic regime probabilities
%         mom.mean [ny x 1] unconditional mean of y (levels, declaration order)
%         mom.mean_regime [ny x S] E[y_t | s_t = s]
%         mom.var [ny x ny] unconditional covariance; mom.second_regime{s} = E[(y-ys)(y-ys)' | s_t = s]  (order >= 2)
%         order >= 3, for the variables v = vars, yh = y - ys:
%           mom.vars; mom.raw_regime{p} [numel(v)^p x S] = E[yh_v^{(x)p} 1{s_t = s}] (symmetric), p = 3..order;
%           mom.skewness [numel(v) x 1], mom.skewness_regime [numel(v) x S] (about E[y | s_t = s], given s_t = s);
%           order 4: mom.kurtosis, mom.kurtosis_regime (not in excess)
%         mom.blocks (layer tuples, in solution order), mom.n_unknowns (folded unknowns, all regimes),
%         mom.residual (max over blocks of the relative residual of the diagonal-block equation)
if nargin < 3, order = 2; end
if nargin < 4 || isempty(solver), solver = 'doubling'; end
if islogical(order), order = 1 + double(order); end
if ~(isscalar(order) && any(order == 1:4)), error('ms_korder_moments: order must be 1, 2, 3 or 4'), end
second = order >= 2;
K = ms.order; S = numel(ms.g); P = ms.P;
ny = numel(ms.order_var); n = numel(ms.state_var);
srow = zeros(n, 1);
for i = 1:n, srow(i) = find(ms.order_var == ms.state_var(i)); end
inv_ov(ms.order_var) = 1:ny;

pi_ = null(P' - eye(S));
if size(pi_, 2) ~= 1, error('ms_korder_moments: the chain P has %d ergodic classes (one required)', size(pi_, 2)), end
pi_ = pi_ / sum(pi_);

W = K * order;
Mu = ms_korder_shock_moments(shocks, W); nu = numel(Mu{2});
nz = cellfun(@(x) any(x(:) ~= 0), Mu);                % nz(J+1): E[u^{(x)J}] is not identically zero
[Tst, Ty] = build_terms(ms, srow, n, nu, K, S, ny);
Phi = ms.g{1}.y1u0s0(srow, :);
for s = 2:S
    if ~isequal(ms.g{s}.y1u0s0(srow, :), Phi), error('ms_korder_moments: g_y differs across regimes'), end
end

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
        if all([terms.dg]) || ~nz(sum([terms.j]) + 1), continue, end
        Tsum = Tsum + mul * contract(terms, getnu, Mu, n, nu, S);
    end
    R = zeros(C, S);
    for s = 1:S, R(:, s) = accumarray(ic, Tsum(:, s), [C, 1]) ./ cnt; end
    M = diag_solve(P, Phi, R, ic, rep, m, n, solver);
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
    if nz(Ty(t).j + 1), Ey = Ey + contract(Ty(t), getnu, Mu, n, nu, S); end
end
mom.mean_regime = ys(inv_ov) + Ey(inv_ov, :) ./ pi_';
mom.mean = ys(inv_ov) + sum(Ey(inv_ov, :), 2);
if second
    E2 = zeros(ny * ny, S);
    for t1 = 1:numel(Ty)
        for t2 = t1:numel(Ty)
            if ~nz(Ty(t1).j + Ty(t2).j + 1), continue, end
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
if order >= 3
    if nargin < 5 || isempty(vars), vars = 1:ny; end
    vars = vars(:)'; r = numel(vars); rs = inv_ov(vars);
    Tyv = Ty;
    for t = 1:numel(Tyv), Tyv(t).G = Tyv(t).G(rs, :, :); end
    mom.vars = vars; mom.raw_regime = cell(1, order);
    for p = 3:order, mom.raw_regime{p} = ymoment(Tyv, p, getnu, Mu, nz, n, nu, S, r); end
    a = cell(1, order);                                % a{p}(i, s) = E[yh_i^p 1{s_t = s}]
    a{1} = Ey(rs, :);
    a{2} = zeros(r, S);
    for s = 1:S, a{2}(:, s) = diag(E2(rs, rs, s)); end
    for p = 3:order, a{p} = mom.raw_regime{p}(1 + (0:r - 1)' * sum(r.^(0:p - 1)), :); end
    A = cellfun(@(x) sum(x, 2), a, 'UniformOutput', false);
    [mom.skewness, mom.kurtosis] = std_moments(A, order);
    c = cellfun(@(x) x ./ pi_', a, 'UniformOutput', false);
    [mom.skewness_regime, mom.kurtosis_regime] = std_moments(c, order);
    if order < 4, mom = rmfield(mom, {'kurtosis', 'kurtosis_regime'}); end
end
end

function [sk, ku] = std_moments(a, order)
% skewness and kurtosis from the raw moments a{p} = E[yh^p] (elementwise)
m = a{1}; v = a{2} - m.^2;
sk = (a{3} - 3 * m .* a{2} + 2 * m.^3) ./ v.^1.5;
ku = [];
if order >= 4, ku = (a{4} - 4 * m .* a{3} + 6 * m.^2 .* a{2} - 3 * m.^4) ./ v.^2; end
end

function E = ymoment(Ty, p, getnu, Mu, nz, n, nu, S, r)
% E[yh^{(x)p} 1{s_t = s}] [r^p x S]: sum over multisets of p terms of one arrangement times its multiplicity,
% symmetrised over the p factors
C = multisets(numel(Ty), p); js = [Ty.j];
E = zeros(r^p, S);
for c = 1:size(C, 1)
    idx = C(c, :);
    if ~nz(sum(js(idx)) + 1), continue, end
    E = E + multinom(idx) * contract(Ty(idx), getnu, Mu, n, nu, S);
end
E = reshape(E, [r * ones(1, p), S]); Es = zeros(size(E));
pm = perms(1:p);
for k = 1:size(pm, 1), Es = Es + permute(E, [pm(k, :), p + 1]); end
E = reshape(Es / size(pm, 1), r^p, S);
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

function M = diag_solve(P, Phi, R, ic, rep, m, n, solver)
% M - Phi^{[b]} M P = R  (M, R: C x S)
[C, S] = size(R);
switch solver
    case 'doubling'
        M = R; F = Phi; Q = P;
        for it = 1:60
            AM = zeros(C, S);
            for s = 1:S, AM(:, s) = apply_sym(F, M(:, s), ic, rep, m, n); end
            dM = AM * Q; M = M + dM;
            if norm(dM, 'fro') <= 1e-17 * norm(M, 'fro'), return, end
            F = F * F; Q = Q * Q;
        end
        error('ms_korder_moments: doubling did not converge (spectral radius of Phi^{(x)%d} >= 1?)', m);
    case 'gensylv'
        % unfolded regime s: X_s (n x n^(m-1)); X_s' - Phi sum_s X_s P(s, s') Phi'^{(x)(m-1)} = R_s'
        Ru = R(ic, :);
        Xb = reshape(permute(reshape(Ru, n, n^(m - 1), S), [1 3 2]), n * S, n^(m - 1));
        if m == 1
            Xb = (eye(n * S) - kron(P', Phi)) \ Xb;
        else
            Xb = gensylv(m - 1, eye(n * S), -kron(P', Phi), Phi', Xb);
        end
        Xu = reshape(permute(reshape(Xb, n, S, n^(m - 1)), [1 3 2]), n^m, S);
        M = Xu(rep, :);
    case 'lu'
        A = zeros(C);
        for c = 1:C
            e = zeros(C, 1); e(c) = 1; A(:, c) = apply_sym(Phi, e, ic, rep, m, n);
        end
        M = reshape((eye(C * S) - kron(P', A)) \ R(:), C, S);
    otherwise
        error('ms_korder_moments: unknown solver %s', solver);
end
end
