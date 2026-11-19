// Markov-switching perturbation with non-Gaussian innovations, three regimes (P circulant, complex eigenvalues),
// against exact results. Model
//   y_t = rho y_{t-1} + a_t + e_t,  p_t = y_t + beta E_t p_{t+1},  w_t = e_t^2,  v_t = e_t^3,
//   z_t = E_t[v_{t+1} + a_{t+1} w_{t+1}],  a_t = ahat(s_t),
// e_t iid, mean 0, raw moments m_d = E[e^d]: (1) a normal mixture given by M_.Raw_moments_e, (2) a skew normal given by
// M_.Skew_e. Exact: z_t = m_3 + m_2 (P ahat)(s_t), so the order-3 term of z in regime s is m_3 + m_2 (P ahat)(s) and the
// other pure sigma terms of z are 0; p_t = c y_t + d(s_t) as in ms_linear_exact.mod. Moments of y up to order 4 by
// the jump-linear recursion with E[(ahat(s) + e)^k] = sum_i C(k,i) ahat(s)^(k-i) m_i.
// Tolerances (fixed before the first run): skew normal raw moments vs numerical integration 1e-9 (relative to
// max(1, |m_d|)); sigma terms 1e-12 abs; moments 1e-10 relative.
var y p w v z a;
varexo eps;
parameters ms_zero rho beta;
rho = 0.8; beta = 0.95; ms_zero = 0;
model;
y = rho*y(-1) + a + eps;
p = y + beta*p(+1);
w = eps^2;
v = eps^3;
z = v(+1) + a(+1)*w(+1);
[ms='a'] a = ms_zero;
end;
initval; y = 0; p = 0; w = 0; v = 0; z = 0; a = 0; end;
shocks; var eps; stderr 0.3; end;
steady; check;
stoch_simul(order=3, k_order_solver, irf=0, noprint, nomoments, nocorr, nofunctions);

P = [0.60 0.35 0.05; 0.05 0.60 0.35; 0.35 0.05 0.60];
ahat = [0.5 -0.2 -1.0];
S = 3; K = 3; D = 4*K; fails = {};
c = 1 / (1 - beta*rho);
pi_ = null(P' - eye(S)); pi_ = pi_ / sum(pi_);
dex = (eye(S) - beta*P) \ (beta*P*(c*ahat'));
% (2) skew normal raw moments vs numerical integration of the density (Azzalini parameterisation)
sk = -0.6; vr = 0.09;
s23 = abs(sk)^(2/3); dl = sign(sk)*sqrt(pi/2*s23/(s23 + ((4-pi)/2)^(2/3)));
om = sqrt(vr/(1 - 2*dl^2/pi)); al = dl/sqrt(1 - dl^2); m0 = om*dl*sqrt(2/pi);
fsn = @(x) 2/om*normpdf((x + m0)/om).*normcdf(al*(x + m0)/om);
msn = skew_normal_raw_moments(vr, sk, D);
mint = arrayfun(@(d) integral(@(x) x.^d.*fsn(x), -12, 12, 'AbsTol', 1e-16, 'RelTol', 1e-13), 1:D);
e0 = max(abs(msn - mint) ./ max(1, abs(mint)));
fprintf('skew normal raw moments vs integration: max rel err %.2e\n', e0);
if e0 > 1e-9, fails{end+1} = 'skew_normal_raw_moments'; end
cases = {'mixture (Raw_moments_e)', mixture_normal_raw_moments([0.9 0.1], [0 -0.8], [0.2 0.4], D), zeros(0, 4);
         'skew normal (Skew_e)', msn, [1 1 1 sk]};
for ic = 1:2
lbl = cases{ic, 1}; m = cases{ic, 2}; sk_e = cases{ic, 3};
M_.Sigma_e = m(2);
if isempty(sk_e), M_.Raw_moments_e = m; else, M_.Raw_moments_e = []; end
M_.Skew_e = sk_e;
[info, oo_] = stoch_simul(M_, options_, oo_, {});
if info(1), error('stoch_simul failed (%s)', lbl), end
oo_ = ms_korder(M_, options_, oo_, struct('P', P, 'values', struct('a', ahat)), struct('order', 4, 'moment_method', 'exact'));
ms = oo_.ms.solution; mom = oo_.ms.moments;
iz = find(ms.order_var == find(strcmp(M_.endo_names, 'z')));
zex = m(3) + m(2) * (P * ahat');
for s = 1:S
    e = max(abs([ms.g{s}.y0u0s3(iz)/6 - zex(s), ms.g{s}.y0u0s1(iz), ms.g{s}.y0u0s2(iz)]));
    fprintf('%s: sigma terms of z, regime %d, max abs err %.2e\n', lbl, s, e);
    if e > 1e-12, fails{end+1} = sprintf('%s sigma terms of z regime %d', lbl, s); end
end
mm = [1 m];                                          % mm(i+1) = E[e^i]
eu = cell(1, 5);
for k = 0:4, eu{k+1} = sum(cell2mat(arrayfun(@(i) nchoosek(k, i) * ahat'.^(k-i) * mm(i+1), 0:k, 'UniformOutput', false)), 2); end
mk = {pi_};
for k = 1:4
    rhs = zeros(S, 1);
    for j = 0:k-1, rhs = rhs + nchoosek(k, j) * rho^j * (P' * mk{j+1}) .* eu{k-j+1}; end
    mk{k+1} = (eye(S) - rho^k * P') \ rhs;
end
pk = @(k) sum(cell2mat(arrayfun(@(j) nchoosek(k, j) * c^j * mk{j+1} .* dex.^(k-j), 0:k, 'UniformOutput', false)), 2);
iv = cellfun(@(n) find(strcmp(M_.endo_names, n)), {'y', 'p', 'z'});
r = numel(mom.vars); ky = find(mom.vars == iv(1)); kp = find(mom.vars == iv(2));
dg = @(q, i) 1 + (i - 1) * sum(r.^(0:q-1));
skf = @(a1, a2, a3) (a3 - 3*a1.*a2 + 2*a1.^3) ./ (a2 - a1.^2).^1.5;
kuf = @(a1, a2, a3, a4) (a4 - 4*a1.*a3 + 6*a1.^2.*a2 - 3*a1.^4) ./ (a2 - a1.^2).^2;
py = cellfun(@(x) sum(x), mk); pp = arrayfun(@(k) sum(pk(k)), 1:4);
Vy = py(3) - py(2)^2;
chk = {'E[y|s]', mom.mean_regime(iv(1), :)', mk{2} ./ pi_; 'E[z|s]', mom.mean_regime(iv(3), :)', zex;
       'Var(y)', mom.var(iv(1), iv(1)), Vy;
       'E[y^3 1{s}]', mom.raw_regime{3}(dg(3, ky), :)', mk{4}; 'E[y^4 1{s}]', mom.raw_regime{4}(dg(4, ky), :)', mk{5};
       'E[p^3 1{s}]', mom.raw_regime{3}(dg(3, kp), :)', pk(3); 'E[p^4 1{s}]', mom.raw_regime{4}(dg(4, kp), :)', pk(4);
       'skew(y)', mom.skewness(ky), skf(py(2), py(3), py(4)); 'kurt(y)', mom.kurtosis(ky), kuf(py(2), py(3), py(4), py(5));
       'skew(p)', mom.skewness(kp), skf(pp(1), pp(2), pp(3)); 'kurt(p)', mom.kurtosis(kp), kuf(pp(1), pp(2), pp(3), pp(4));
       'skew(y|s)', mom.skewness_regime(ky, :)', skf(mk{2}./pi_, mk{3}./pi_, mk{4}./pi_);
       'kurt(p|s)', mom.kurtosis_regime(kp, :)', kuf(pk(1)./pi_, pk(2)./pi_, pk(3)./pi_, pk(4)./pi_)};
for i = 1:size(chk, 1)
    e = max(abs(chk{i, 2} - chk{i, 3})) / max(abs(chk{i, 3}));
    fprintf('%s: %-12s rel err %.2e\n', lbl, chk{i, 1}, e);
    if e > 1e-10, fails{end+1} = sprintf('%s %s', lbl, chk{i, 1}); end
end
end
M_.Raw_moments_e = []; M_.Skew_e = zeros(0, 4);
if ~isempty(fails), error('ms_nongaussian_exact: failed: %s', strjoin(fails, ', ')), end
disp('ms_nongaussian_exact: all checks pass')
