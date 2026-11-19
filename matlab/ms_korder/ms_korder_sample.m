function smp = ms_korder_sample(ms, M_, T, opts)
% smp = ms_korder_sample(ms, M_, T[, opts])
% Labelled sample of the ergodic set of the Markov-switching order-K solution, for supervised pre-training.
% A regime path and innovations are drawn, the PRUNED solution is simulated (ms_korder_simul, Fortran engine), and in
% each kept period t the label is the UNPRUNED policy of regime s_t at the simulated state:
%     label_t = ys + g(s_t, x_{t-1}, u_t),   x_{t-1} = y_{t-1}(states) - ys(states)  (pruned state, sum of its layers),
% g as in ms_korder_eval (sigma = 1). The pruned y_t is stored for reference (label_t - y_t is the pruning gap).
%
% INPUTS  ms  output of ms_korder_solve;  M_  the model structure (Sigma_e, Raw_moments_e, Skew_e);  T  kept periods
%         opts.burn (default 1000), opts.thin (default 1: keep every thin-th period after the burn-in),
%         opts.s0 (initial regime, default 1), opts.seed (default 0: regimes seed, innovations rng(seed + 1)),
%         opts.draw (innovation draw function, see ms_korder_draw_shocks; required for shocks given by Raw_moments_e),
%         opts.file (base name: writes <file>.mat and <file>.csv; default none)
% OUTPUT  smp.t [T x 1] period index (after the burn-in), smp.s [T x 1] regime, smp.x [T x n] state deviations
%         x_{t-1} (columns smp.state_names), smp.u [T x nu] innovations u_t (smp.exo_names), smp.label [T x ny] and
%         smp.y_pruned [T x ny] levels (smp.endo_names, declaration order), smp.ys, smp.P.
% CSV columns: t, s, x_<state>, u_<exo>, label_<endo>, pruned_<endo>; written with 17 significant digits.
if nargin < 4, opts = struct(); end
def = struct('burn', 1000, 'thin', 1, 's0', 1, 'seed', 0, 'draw', [], 'file', '');
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts, f{i}), opts.(f{i}) = def.(f{i}); end
end
Ttot = opts.burn + T * opts.thin;
s = ms_markov_chain(ms.P, Ttot, opts.s0, opts.seed);
rng(opts.seed + 1, 'twister');
u = ms_korder_draw_shocks(M_, Ttot, opts.draw);
y = ms_korder_simul(ms, s, u, true);
ny = numel(ms.order_var); inv_ov(ms.order_var) = 1:ny;
ys = ms.ys(:); st = ms.state_var(:);
xall = [zeros(numel(st), 1), y(st, 1:end - 1) - ys(st)];       % x_{t-1}, t = 1..Ttot (x_0 = 0)
keep = opts.burn + (1:T) * opts.thin;
lab = zeros(ny, T);
for q = 1:T
    t = keep(q);
    g = ms_korder_eval(ms, s(t), xall(:, t), u(:, t), 1);
    lab(:, q) = ys + g(inv_ov);
end
smp = struct('t', keep(:) - opts.burn, 's', s(keep)', 'x', xall(:, keep)', 'u', u(:, keep)', 'label', lab', ...
             'y_pruned', y(:, keep)', 'ys', ys, 'P', ms.P, 'endo_names', {M_.endo_names(:)'}, ...
             'exo_names', {M_.exo_names(:)'}, 'state_names', {M_.endo_names(st)'}, 'opts', rmfield(opts, 'draw'));
if ~isempty(opts.file)
    save([opts.file '.mat'], '-struct', 'smp', '-v7');
    hdr = [{'t', 's'}, strcat('x_', smp.state_names), strcat('u_', smp.exo_names), strcat('label_', smp.endo_names), ...
           strcat('pruned_', smp.endo_names)];
    fid = fopen([opts.file '.csv'], 'w');
    fprintf(fid, '%s\n', strjoin(hdr, ','));
    D = [smp.t, smp.s, smp.x, smp.u, smp.label, smp.y_pruned];
    fmt = [strjoin(repmat({'%.17g'}, 1, size(D, 2)), ','), '\n'];
    fprintf(fid, fmt, D');
    fclose(fid);
end
end
