// Markov-switching k-order perturbation with one regime (P = 1, ahat = 0) reproduces stock Dynare:
//  (1) policy derivatives up to order 3 vs oo_.dr (k_order_solver); odd sigma derivatives are 0;
//  (2) pruned order-3 simulation vs simult_;
//  (3) exact pruned moments (ms_korder_moments) vs pruned_SS.pruned_state_space_system.
// Tolerances (fixed before the first run): (1) 1e-12 x max(|stock|, 1), odd sigma <= 1e-14; (2) 1e-12 x path scale
// per variable; (3) 1e-10 x scale (means: largest deviation from steady state; covariance: largest entry).
@#include "rbc3_ms_model.inc"
stoch_simul(order=3, k_order_solver, pruning, irf=0, noprint, nomoments, nocorr, nofunctions);

ms = ms_korder_solve(M_, options_, oo_.dr, 1, zeros(2, 1));
fails = {};
pairs = {'y1u0s0','ghx'; 'y0u1s0','ghu'; 'y2u0s0','ghxx'; 'y1u1s0','ghxu'; 'y0u2s0','ghuu'; 'y0u0s2','ghs2';
         'y3u0s0','ghxxx'; 'y2u1s0','ghxxu'; 'y1u2s0','ghxuu'; 'y0u3s0','ghuuu'; 'y1u0s2','ghxss'; 'y0u1s2','ghuss'};
for q = 1:size(pairs, 1)
    a = ms.g{1}.(pairs{q, 1}); b = oo_.dr.(pairs{q, 2});
    e = max(abs(a(:) - b(:))) / max(max(abs(b(:))), 1);
    fprintf('%-7s vs %-6s rel err %.2e\n', pairs{q, 1}, pairs{q, 2}, e);
    if e > 1e-12, fails{end+1} = pairs{q, 1}; end
end
odd = {'y0u0s1', 'y1u0s1', 'y0u1s1', 'y0u0s3', 'y2u0s1', 'y1u1s1', 'y0u2s1'};
for q = 1:numel(odd)
    e = max(abs(ms.g{1}.(odd{q})(:)));
    if e > 1e-14, fails{end+1} = odd{q}; end
end

rng(1, 'twister'); T = 200; u = chol(M_.Sigma_e, 'lower') * randn(2, T);
ik = strcmp(M_.endo_names, 'k'); iz = strcmp(M_.endo_names, 'z'); idd = strcmp(M_.endo_names, 'd');
y0 = oo_.dr.ys; y0(ik) = 1.02 * y0(ik); y0(iz) = 0.01; y0(idd) = -0.02;
ystock = simult_(M_, options_, y0, oo_.dr, u', 3);
x0 = y0(ms.state_var) - oo_.dr.ys(ms.state_var);
yms = ms_korder_simul(ms, ones(1, T), u, true, x0);
for v = 1:4
    sc = max(abs(ystock(v, 2:end) - oo_.dr.ys(v))); e = max(abs(yms(v, :) - ystock(v, 2:end)));
    fprintf('simulation %-2s max abs err %.2e scale %.2e\n', M_.endo_names{v}, e, sc);
    if e > 1e-12 * sc, fails{end+1} = ['sim ' M_.endo_names{v}]; end
end

mom = ms_korder_moments(ms, M_.Sigma_e, true);
pss = pruned_SS.pruned_state_space_system(M_, options_, oo_.dr, (1:M_.endo_nbr)', 0, false, false);
inv_ov(oo_.dr.order_var) = 1:M_.endo_nbr;
Es = pss.E_y(inv_ov); Vs = pss.Var_y(inv_ov, inv_ov);
em = max(abs(mom.mean - Es)) / max(abs(Es - oo_.dr.ys)); ev = max(abs(mom.var(:) - Vs(:))) / max(abs(Vs(:)));
fprintf('moments: mean rel err %.2e, covariance rel err %.2e\n', em, ev);
if em > 1e-10, fails{end+1} = 'mean'; end
if ev > 1e-10, fails{end+1} = 'covariance'; end
if ~isempty(fails), error('ms_one_regime: failed: %s', strjoin(fails, ', ')), end
disp('ms_one_regime: all checks pass')
