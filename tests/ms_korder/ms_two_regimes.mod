// Markov-switching perturbation of the 3-state RBC (rbc3_ms_model.inc), two regimes, order 3:
//  (1) exact pruned moments with the three diagonal-block solvers (doubling, Dynare gensylv, dense LU) agree;
//  (2) the regime-conditional moments aggregate: mean = sum_s pi(s) mean|s, and
//      var = sum_s pi(s) second|s - (mean - ys)(mean - ys)';
//  (3) ms_korder (solution + moments + simulation) runs and stores its output in oo_.ms.
// Tolerances (fixed before the first run): (1), (2) 1e-10 relative to the largest entry.
@#include "rbc3_ms_model.inc"
stoch_simul(order=3, k_order_solver, irf=0, noprint, nomoments, nocorr, nofunctions);

regimes = struct('P', [0.95 0.05; 0.20 0.80], 'values', struct('aA', [0 -0.10], 'aD', [0 0.03]), ...
                 'names', {{'normal', 'bad'}});
oo_ = ms_korder(M_, options_, oo_, regimes, struct('periods', 100));
ms = oo_.ms.solution; md = oo_.ms.moments; fails = {};
r = @(a, b) max(abs(a(:) - b(:))) / max(abs(b(:)));
for sv = {'gensylv', 'lu'}
    mo = ms_korder_moments(ms, M_.Sigma_e, true, sv{1});
    e = max([r(md.mean - oo_.dr.ys, mo.mean - oo_.dr.ys), r(md.var, mo.var), ...
             r(cat(3, md.second_regime{:}), cat(3, mo.second_regime{:}))]);
    fprintf('doubling vs %-8s rel err %.2e\n', sv{1}, e);
    if e > 1e-10, fails{end+1} = sv{1}; end
end
ys = oo_.dr.ys;
em = r(md.mean_regime * md.pi, md.mean);
SR = zeros(size(md.var));
for s = 1:2, SR = SR + md.pi(s) * md.second_regime{s}; end
ev = r(SR - (md.mean - ys) * (md.mean - ys)', md.var);
fprintf('aggregation: mean %.2e, variance %.2e\n', em, ev);
if em > 1e-10, fails{end+1} = 'mean aggregation'; end
if ev > 1e-10, fails{end+1} = 'variance aggregation'; end
if ~isequal(size(oo_.ms.simulation.y), [M_.endo_nbr 100]), fails{end+1} = 'simulation'; end
if ~isempty(fails), error('ms_two_regimes: failed: %s', strjoin(fails, ', ')), end
disp('ms_two_regimes: all checks pass')
