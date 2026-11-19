// Markov-switching perturbation of a linear model, three regimes, against its exact solution, for a transition matrix
// with real eigenvalues and for one with complex eigenvalues (circulant, eigenvalues 1, 0.4 +- 0.26i):
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

Pcases = {[0.90 0.07 0.03; 0.10 0.85 0.05; 0.05 0.15 0.80], [0.60 0.35 0.05; 0.05 0.60 0.35; 0.35 0.05 0.60]};
ahat = [0.5 -0.2 -1.0];
S = 3; fails = {};
c = 1 / (1 - beta*rho);
for ic = 1:numel(Pcases)
P = Pcases{ic};
regimes = struct('P', P, 'values', struct('a', ahat));
oo_ = ms_korder(M_, options_, oo_, regimes, struct('periods', 50));
ms = oo_.ms.solution; mom = oo_.ms.moments;
fprintf('case %d: eigenvalues of P %s\n', ic, mat2str(eig(P).', 4));
dex = (eye(S) - beta*P) \ (beta*P*(c*ahat'));
iy = find(ms.order_var == find(strcmp(M_.endo_names, 'y'))); ip = find(ms.order_var == find(strcmp(M_.endo_names, 'p')));
for s = 1:S
    gs = ms.g{s}.y0u0s1;
    e = max(abs([gs(iy) - ahat(s), gs(ip) - (c*ahat(s) + dex(s))]));
    if e > 1e-12, fails{end+1} = sprintf('case %d g_sigma regime %d (%.1e)', ic, s, e); end
    for f = {'y2u0s0', 'y1u1s0', 'y0u0s2', 'y1u0s1', 'y0u1s2', 'y0u0s3'}
        if max(abs(ms.g{s}.(f{1})(:))) > 1e-12, fails{end+1} = sprintf('case %d %s regime %d', ic, f{1}, s); end
    end
end
pi_ = null(P' - eye(S)); pi_ = pi_ / sum(pi_);
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
    fprintf('case %d %-9s rel err %.2e\n', ic, chk{i, 1}, e);
    if e > 1e-10, fails{end+1} = sprintf('case %d %s', ic, chk{i, 1}); end
end
if size(oo_.ms.simulation.y, 2) ~= 50, fails{end+1} = sprintf('case %d simulation length', ic); end
% Third and fourth moments: mk{k+1}(s) = E[y^k 1{s_t = s}] by the jump-linear recursion
%   mk{k+1} = (I - rho^k P')^{-1} sum_{j<k} C(k,j) rho^j (P' mk{j+1}) .* E[(ahat + sig u)^(k-j)],
% E[p^k 1{s}] = sum_j C(k,j) c^j mk{j+1} .* d^(k-j). Tolerance 1e-10 relative (fixed before the first run).
eu = {ones(S, 1), ahat', ahat'.^2 + sig^2, ahat'.^3 + 3*sig^2*ahat', ahat'.^4 + 6*sig^2*ahat'.^2 + 3*sig^4};
mk = {pi_};
for k = 1:4
    rhs = zeros(S, 1);
    for j = 0:k-1, rhs = rhs + nchoosek(k, j) * rho^j * (P' * mk{j+1}) .* eu{k-j+1}; end
    mk{k+1} = (eye(S) - rho^k * P') \ rhs;
end
pk = @(k) sum(cell2mat(arrayfun(@(j) nchoosek(k, j) * c^j * mk{j+1} .* dex.^(k-j), 0:k, 'UniformOutput', false)), 2);
mo4 = ms_korder_moments(ms, M_.Sigma_e, 4, 'doubling', [find(iyd), find(ipd)]);
r = 2; dg = @(p, i) 1 + (i - 1) * sum(r.^(0:p-1));
sk = @(a1, a2, a3) (a3 - 3*a1.*a2 + 2*a1.^3) ./ (a2 - a1.^2).^1.5;
ku = @(a1, a2, a3, a4) (a4 - 4*a1.*a3 + 6*a1.^2.*a2 - 3*a1.^4) ./ (a2 - a1.^2).^2;
py = cellfun(@(x) sum(x), mk); pp = arrayfun(@(k) sum(pk(k)), 1:4);
chk = {'E[y^3 1{s}]', mo4.raw_regime{3}(dg(3, 1), :)', mk{4}; 'E[y^4 1{s}]', mo4.raw_regime{4}(dg(4, 1), :)', mk{5};
       'E[p^3 1{s}]', mo4.raw_regime{3}(dg(3, 2), :)', pk(3); 'E[p^4 1{s}]', mo4.raw_regime{4}(dg(4, 2), :)', pk(4);
       'E[y^2 p 1{s}]', mo4.raw_regime{3}(1 + r^2, :)', c*mk{4} + dex.*mk{3};
       'skew(y)', mo4.skewness(1), sk(py(2), py(3), py(4)); 'kurt(y)', mo4.kurtosis(1), ku(py(2), py(3), py(4), py(5));
       'skew(p)', mo4.skewness(2), sk(pp(1), pp(2), pp(3)); 'kurt(p)', mo4.kurtosis(2), ku(pp(1), pp(2), pp(3), pp(4));
       'skew(y|s)', mo4.skewness_regime(1, :)', sk(mk{2}./pi_, mk{3}./pi_, mk{4}./pi_);
       'kurt(p|s)', mo4.kurtosis_regime(2, :)', ku(pk(1)./pi_, pk(2)./pi_, pk(3)./pi_, pk(4)./pi_)};
for i = 1:size(chk, 1)
    e = max(abs(chk{i, 2} - chk{i, 3})) / max(abs(chk{i, 3}));
    fprintf('case %d %-13s rel err %.2e\n', ic, chk{i, 1}, e);
    if e > 1e-10, fails{end+1} = sprintf('case %d %s', ic, chk{i, 1}); end
end
end
if ~isempty(fails), error('ms_linear_exact: failed: %s', strjoin(fails, ', ')), end
disp('ms_linear_exact: all checks pass')
