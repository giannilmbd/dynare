function sm = ms_korder_sim_moments(ms, shocks, vars, T, opts)
% sm = ms_korder_sim_moments(ms, shocks, vars, T[, opts])
% Moments of the pruned Markov-switching solution by simulation, with batch-means standard errors.
%
% INPUTS  ms       output of ms_korder_solve
%         shocks   M_ or the covariance of Gaussian innovations; innovations drawn by ms_korder_draw_shocks
%         vars     variables (declaration indices)
%         T        number of periods kept after the burn-in
%         opts     .burn (default 1000), .batches (default 100; T must be a multiple), .s0 (initial regime, default
%                  1), .seed (default 0: regime path seed 'seed', innovations rng(seed + 1)), .draw (function handle
%                  draw(T) -> [nu x T], see ms_korder_draw_shocks; required for shocks given by Raw_moments_e)
% OUTPUT  sm.vars; sm.T; sm.regime_share [S x 1];
%         sm.mean, sm.std, sm.skewness, sm.kurtosis [numel(vars) x 1] (kurtosis not in excess),
%         sm.mean_regime, sm.std_regime, sm.skewness_regime, sm.kurtosis_regime [numel(vars) x S] (given s_t = s);
%         sm.se.<same fields>: batch-means standard errors (standard deviation of the batch estimates / sqrt(batches))
if nargin < 5, opts = struct(); end
def = struct('burn', 1000, 'batches', 100, 's0', 1, 'seed', 0, 'draw', []);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts, f{i}), opts.(f{i}) = def.(f{i}); end
end
B = opts.batches;
if mod(T, B) ~= 0, error('ms_korder_sim_moments: T must be a multiple of the number of batches (%d)', B), end
S = numel(ms.g); Tb = T / B;
s = ms_markov_chain(ms.P, T + opts.burn, opts.s0, opts.seed);
rng(opts.seed + 1, 'twister');
u = ms_korder_draw_shocks(shocks, T + opts.burn, opts.draw);
y = ms_korder_simul(ms, s, u, true);
y = y(vars(:), opts.burn + 1:end); s = s(opts.burn + 1:end);
stats = @(x) [mean(x, 2), std(x, 1, 2), skw(x), krt(x)];
nv = numel(vars);
full = stats(y);
reg = zeros(nv, 4, S);
for r = 1:S, reg(:, :, r) = stats(y(:, s == r)); end
bf = zeros(nv, 4, B); br = nan(nv, 4, S, B);
for b = 1:B
    j = (b - 1) * Tb + (1:Tb); yb = y(:, j); sb = s(j);
    bf(:, :, b) = stats(yb);
    for r = 1:S
        if nnz(sb == r) > 3, br(:, :, r, b) = stats(yb(:, sb == r)); end
    end
end
sef = std(bf, 0, 3) / sqrt(B);
ser = zeros(nv, 4, S);
for r = 1:S
    x = squeeze(br(:, :, r, :)); x = reshape(x, nv, 4, []);
    nb = sum(~isnan(x(1, 1, :)));
    ser(:, :, r) = std(x, 0, 3, 'omitnan') / sqrt(nb);
end
nm = {'mean', 'std', 'skewness', 'kurtosis'};
sm = struct('vars', vars(:)', 'T', T, 'regime_share', arrayfun(@(r) mean(s == r), (1:S)'));
for q = 1:4
    sm.(nm{q}) = full(:, q);
    sm.([nm{q} '_regime']) = reshape(reg(:, q, :), nv, S);
    sm.se.(nm{q}) = sef(:, q);
    sm.se.([nm{q} '_regime']) = reshape(ser(:, q, :), nv, S);
end
end

function k = skw(x)
d = x - mean(x, 2); k = mean(d.^3, 2) ./ mean(d.^2, 2).^1.5;
end

function k = krt(x)
d = x - mean(x, 2); k = mean(d.^4, 2) ./ mean(d.^2, 2).^2;
end
