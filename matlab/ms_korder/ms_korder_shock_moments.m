function [Mu, ng, raw] = ms_korder_shock_moments(shocks, W)
% [Mu, ng, raw] = ms_korder_shock_moments(shocks, W)
% Moments of the innovations, Mu{J+1} = E[u^{(x)J}] [nu^J x 1], J = 0..W.
%
% INPUTS  shocks  either the covariance [nu x nu] (Gaussian innovations), or a structure with the fields of M_ that
%                 define the innovations: Sigma_e, and optionally Raw_moments_e and Skew_e, read as by the k-order
%                 solver: a shock with a fully specified row of Raw_moments_e (entry (i, d) = E[e_i^d]) or with a
%                 nonzero entry (i, i, i) of Skew_e (independent skew normal) is non-Gaussian and independent of the
%                 other shocks; the remaining shocks are jointly Gaussian with covariance Sigma_e
%         W       highest order
% OUTPUT  Mu      cell 1 x (W+1)
%         ng      indices of the non-Gaussian shocks
%         raw     [numel(ng) x W] their raw moments E[e^d], d = 1..W
%
% Gaussian block: Isserlis, pairing mode 1 with each other mode. With non-Gaussian shocks, E[u_{i_1} ... u_{i_J}]
% depends only on the counts c_i of each shock among (i_1, ..., i_J) and equals
%     E_G[prod_{i in G} u_i^{c_i}] * prod_{i in ng} E[e_i^{c_i}],
% by independence; it is computed once per count vector.
if isnumeric(shocks), shocks = struct('Sigma_e', shocks); end
Sig = shocks.Sigma_e; nu = size(Sig, 1);
rawin = NaN(nu, 0);
if isfield(shocks, 'Raw_moments_e') && ~isempty(shocks.Raw_moments_e), rawin = shocks.Raw_moments_e; end
if size(rawin, 1) ~= nu, error('ms_korder_shock_moments: Raw_moments_e must have %d rows', nu), end
skew = zeros(nu, 1);
if isfield(shocks, 'Skew_e') && ~isempty(shocks.Skew_e)
    for r = 1:size(shocks.Skew_e, 1)
        e = shocks.Skew_e(r, :);
        if e(4) == 0, continue, end
        if e(1) ~= e(2) || e(2) ~= e(3)
            error('ms_korder_shock_moments: co-skewness entries in Skew_e are not supported')
        end
        skew(e(1)) = e(4);
    end
end
ng = []; raw = zeros(0, W);
for i = 1:nu
    row = rawin(i, :);
    given = ~isempty(row) && ~all(isnan(row));
    if ~given && skew(i) == 0, continue, end
    if given
        if any(isnan(row)), error('ms_korder_shock_moments: row %d of Raw_moments_e is partly NaN', i), end
        if skew(i) ~= 0, error('ms_korder_shock_moments: shock %d is both in Raw_moments_e and in Skew_e', i), end
        if numel(row) < W
            error(['ms_korder_shock_moments: these moments need E[e^d] up to d = %d (moment order times ' ...
                   'perturbation order); Raw_moments_e has %d columns'], W, numel(row))
        end
        if abs(row(1)) > 1e-12 * max(1, sqrt(Sig(i, i))), error('ms_korder_shock_moments: shock %d must have mean zero', i), end
        if abs(row(2) - Sig(i, i)) > 1e-10 * max(1, Sig(i, i))
            error('ms_korder_shock_moments: Raw_moments_e(%d, 2) must equal Sigma_e(%d, %d)', i, i, i)
        end
        row = row(1:W);
    else
        row = skew_normal_raw_moments(Sig(i, i), skew(i), W);
    end
    if any(Sig(i, [1:i-1, i+1:nu]) ~= 0)
        error('ms_korder_shock_moments: non-Gaussian shock %d must be uncorrelated with the other shocks', i)
    end
    ng(end + 1) = i; raw(end + 1, :) = row; %#ok<AGROW>
end
G = setdiff(1:nu, ng);
if isempty(ng), Mu = gauss_kron(Sig, nu, W); return, end
nG = numel(G);
if nG > 0, MuG = gauss_kron(Sig(G, G), nG, W); end
Mu = cell(1, W + 1); Mu{1} = 1;
cnt = zeros(1, nu);                                   % count vectors of all multi-indices, first mode fastest
for J = 1:W
    cnt = repmat(cnt, nu, 1) + kron(eye(nu), ones(size(cnt, 1), 1));
    [uc, ~, ic] = unique(cnt, 'rows');
    val = zeros(size(uc, 1), 1);
    for q = 1:size(uc, 1)
        c = uc(q, :);
        v = 1;
        for k = 1:numel(ng)
            if c(ng(k)) > 0, v = v * raw(k, c(ng(k))); end
        end
        cg = c(G); JG = sum(cg);
        if JG > 0
            if mod(JG, 2) == 1, v = 0;
            else
                list = repelem(1:nG, cg);
                v = v * MuG{JG + 1}(1 + sum((list - 1) .* nG.^(0:JG - 1)));
            end
        end
        val(q) = v;
    end
    Mu{J + 1} = val(ic);
end
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
