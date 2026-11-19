// Markov-switching perturbation of a linear model, three regimes, against its exact solution:
//   y_t = rho y_{t-1} + a_t + sig u_t,   p_t = y_t + beta E_t p_{t+1},   a_t = ahat(s_t).
// Exact: p_t = c y_t + d(s_t), c = 1/(1 - beta rho), d = (I - beta P)^{-1} beta P c ahat; the order-1 MS solution is
// exact and all higher derivatives are 0. Exact moments by the Markov jump linear recursion
//   m(s') = E[y_t 1{s_t=s'}] = sum_s P(s,s') rho m(s) + ahat(s') pi(s'),
//   q(s') = E[y_t^2 1{s_t=s'}] = sum_s P(s,s') (rho^2 q(s) + 2 rho ahat(s') m(s)) + (ahat(s')^2 + sig^2) pi(s').
// Tolerances (fixed before the first run): derivatives 1e-12 abs; moments 1e-10 relative.
var y p a;
varexo u;
parameters ms_zero rho beta sig;
rho = 0.8; beta = 0.95; sig = 0.1; ms_zero = 0;
model(linear);
y = rho*y(-1) + a + sig*u;
p = y + beta*p(+1);
[ms='a'] a = ms_zero;
end;
shocks; var u; stderr 1; end;
steady; check;
stoch_simul(order=3, k_order_solver, irf=0, noprint, nomoments, nocorr, nofunctions);

P = [0.90 0.07 0.03; 0.10 0.85 0.05; 0.05 0.15 0.80];
ahat = [0.5 -0.2 -1.0];
regimes = struct('P', P, 'values', struct('a', ahat));
oo_ = ms_korder(M_, options_, oo_, regimes, struct('periods', 50));
ms = oo_.ms.solution; mom = oo_.ms.moments;
S = 3; fails = {};
c = 1 / (1 - beta*rho); dex = (eye(S) - beta*P) \ (beta*P*(c*ahat'));
iy = find(ms.order_var == find(strcmp(M_.endo_names, 'y'))); ip = find(ms.order_var == find(strcmp(M_.endo_names, 'p')));
for s = 1:S
    gs = ms.g{s}.y0u0s1;
    e = max(abs([gs(iy) - ahat(s), gs(ip) - (c*ahat(s) + dex(s))]));
    if e > 1e-12, fails{end+1} = sprintf('g_sigma regime %d (%.1e)', s, e); end
    for f = {'y2u0s0', 'y1u1s0', 'y0u0s2', 'y1u0s1', 'y0u1s2', 'y0u0s3'}
        if max(abs(ms.g{s}.(f{1})(:))) > 1e-12, fails{end+1} = sprintf('%s regime %d', f{1}, s); end
    end
end
[V, D] = eig(P'); [~, i1] = min(abs(diag(D) - 1)); pi_ = V(:, i1) / sum(V(:, i1));
m = (eye(S) - rho*P') \ (ahat' .* pi_);
q = (eye(S) - rho^2*P') \ (2*rho*ahat' .* (P'*m) + (ahat'.^2 + sig^2) .* pi_);
Ey = m ./ pi_; Ep = c*Ey + dex;
Vy = sum(q) - sum(m)^2;
Epp = c^2*q + 2*c*dex.*m + dex.^2.*pi_; Eyp = c*q + dex.*m;
Vp = sum(Epp) - sum(c*m + dex.*pi_)^2; Cyp = sum(Eyp) - sum(m)*sum(c*m + dex.*pi_);
iyd = strcmp(M_.endo_names, 'y'); ipd = strcmp(M_.endo_names, 'p');
chk = {'E[y|s]', mom.mean_regime(iyd, :)', Ey; 'E[p|s]', mom.mean_regime(ipd, :)', Ep; 'Var(y)', mom.var(iyd, iyd), Vy;
       'Var(p)', mom.var(ipd, ipd), Vp; 'Cov(y,p)', mom.var(iyd, ipd), Cyp;
       'E[y^2|s]', cellfun(@(M) M(iyd, iyd), mom.second_regime)', q ./ pi_; 'pi', mom.pi, pi_};
for i = 1:size(chk, 1)
    e = max(abs(chk{i, 2} - chk{i, 3})) / max(abs(chk{i, 3}));
    fprintf('%-9s rel err %.2e\n', chk{i, 1}, e);
    if e > 1e-10, fails{end+1} = chk{i, 1}; end
end
if size(oo_.ms.simulation.y, 2) ~= 50, fails{end+1} = 'simulation length'; end
if ~isempty(fails), error('ms_linear_exact: failed: %s', strjoin(fails, ', ')), end
disp('ms_linear_exact: all checks pass')
