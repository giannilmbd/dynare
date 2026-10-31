/*
 * Checks user-supplied raw moments of independent shocks (M_.Raw_moments_e) in the k-order solver.
 *
 * Model as in tests/skewed_shocks_korder/skew_korder_moments.mod: with ORDER = K,
 *   y = E_t[e(+1)^2 + … + e(+1)^K + e(+1)^2*z(+1)^2]  (last term only for K ≥ 4)
 *     = E[e^2] + … + E[e^K] + E[e^2]*E[z^2].
 * The K-th order perturbation is exact, so the pure sigma columns of dr.pruning.g_k for y are
 *   g_2 = E[e^2], g_3 = E[e^3], g_4 = E[e^4] + E[e^2]*E[z^2], g_k = E[e^k] for k ≥ 5.
 *
 * Test 1: e is a mean-zero two-component normal mixture (weights 0.98/0.02, means 0/-5, s.d. 1/1.5),
 *         raw moments from mixture_normal_raw_moments.m, checked against numerical integration of the
 *         mixture density; z Gaussian (NaN row).
 * Test 2: e skew normal: Raw_moments_e with the skew normal raw moments reproduces the Skew_e solution.
 * Test 3: Gaussian raw moments for both shocks reproduce the default (Isserlis) solution.
 * Test 4: inconsistent inputs (nonzero mean, variance not equal to Sigma_e) raise errors.
 *
 * Usage: dynare raw_moments_korder -DORDER=4
 */

@#ifndef ORDER
@#define ORDER = 4
@#endif

var x y
@#if ORDER >= 4
  q
@#endif
@#for k in 2:ORDER
  p@{k}
@#endfor
;
varexo e z;

model;
x = 0.5*x(-1) + e + z;
@#for k in 2:ORDER
p@{k} = e^@{k};
@#endfor
@#if ORDER >= 4
q = e^2*z^2;
@#endif
y = 0
@#if ORDER >= 4
  + q(+1)
@#endif
@#for k in 2:ORDER
  + p@{k}(+1)
@#endfor
  ;
end;

steady_state_model;
x = 0; y = 0;
@#if ORDER >= 4
q = 0;
@#endif
@#for k in 2:ORDER
p@{k} = 0;
@#endfor
end;

shocks;
var e = 1;
var z = 2;
end;

stoch_simul(order=@{ORDER}, pruning, k_order_solver, irf=0, periods=0, noprint);

K = @{ORDER};
ie = strmatch('e', M_.exo_names, 'exact');
iz = strmatch('z', M_.exo_names, 'exact');
iy = find(oo_.dr.order_var == strmatch('y', M_.endo_names, 'exact'));
pure = @(dr, k) dr.pruning.(sprintf('g_%d', k))(iy, end);
gdiff = @(a, b) max(arrayfun(@(k) max(max(abs(a.pruning.(sprintf('g_%d', k)) - b.pruning.(sprintf('g_%d', k))))), 1:K));
dr_gauss_default = oo_.dr;

% ---- Test 1: normal mixture ----
w = [0.98 0.02]; mu = [0 -5]; sd = [1 1.5];
mraw = mixture_normal_raw_moments(w, mu, sd, K);
mbar = w*mu';
f = @(x) w(1)*normpdf(x + mbar, mu(1), sd(1)) + w(2)*normpdf(x + mbar, mu(2), sd(2));
mint = arrayfun(@(d) integral(@(x) x.^d.*f(x), -40, 40, 'AbsTol', 1e-14, 'RelTol', 1e-13), 1:K);
err_int = max(abs(mraw - mint)./max(1, abs(mint)));
fprintf('Test 1: mixture raw moments vs numerical integration, max rel. error %.2e\n', err_int);
if err_int > 1e-9
    error('mixture_normal_raw_moments does not match numerical integration')
end
M_.Sigma_e(ie, ie) = mraw(2);
M_.Raw_moments_e = NaN(M_.exo_nbr, K);
M_.Raw_moments_e(ie, :) = mraw;
[info, oo_] = stoch_simul(M_, options_, oo_, {});
if info(1), error('stoch_simul failed (test 1)'), end
target = mraw(2:K);
if K >= 4, target(3) = target(3) + mraw(2)*M_.Sigma_e(iz, iz); end
got = arrayfun(@(k) pure(oo_.dr, k), 2:K);
err1 = max(abs(got - target)./max(1, abs(target)));
fprintf('Test 1: pure sigma^k terms vs exact moments, max rel. error %.2e\n', err1);
disp([2:K; got; target]);
if err1 > 1e-10, error('Test 1 failed'), end

% ---- Test 2: skew normal via Skew_e and via Raw_moments_e ----
M_.Sigma_e(ie, ie) = 1;
M_.Raw_moments_e = [];
M_.Skew_e = [ie ie ie -0.8];
[info, oo_] = stoch_simul(M_, options_, oo_, {});
if info(1), error('stoch_simul failed (test 2, Skew_e)'), end
dr_skew = oo_.dr;
% skew normal raw moments by numerical integration of the density (csn_update_specification.m mapping)
s23 = 0.8^(2/3);
dl = -sqrt(pi/2*s23/(s23 + ((4-pi)/2)^(2/3)));
om = 1/sqrt(1 - 2*dl^2/pi); al = dl/sqrt(1 - dl^2); m0 = om*dl*sqrt(2/pi);
fsn = @(x) 2/om*normpdf((x + m0)/om).*normcdf(al*(x + m0)/om);
msn = arrayfun(@(d) integral(@(x) x.^d.*fsn(x), -40, 40, 'AbsTol', 1e-14, 'RelTol', 1e-13), 1:K);
msn(1) = 0; msn(2) = 1;
M_.Skew_e = zeros(0, 4);
M_.Raw_moments_e = NaN(M_.exo_nbr, K);
M_.Raw_moments_e(ie, :) = msn;
[info, oo_] = stoch_simul(M_, options_, oo_, {});
if info(1), error('stoch_simul failed (test 2, Raw_moments_e)'), end
err2 = gdiff(oo_.dr, dr_skew);
fprintf('Test 2: Raw_moments_e (skew normal, by integration) vs Skew_e, max |g_k difference| %.2e\n', err2);
if err2 > 1e-10, error('Test 2 failed'), end

% ---- Test 3: Gaussian raw moments vs default ----
M_.Raw_moments_e = [mixture_normal_raw_moments(1, 0, 1, K); mixture_normal_raw_moments(1, 0, sqrt(2), K)];
M_.Raw_moments_e([ie iz], :) = M_.Raw_moments_e;
[info, oo_] = stoch_simul(M_, options_, oo_, {});
if info(1), error('stoch_simul failed (test 3)'), end
err3 = gdiff(oo_.dr, dr_gauss_default);
fprintf('Test 3: Gaussian Raw_moments_e vs default, max |g_k difference| %.2e\n', err3);
if err3 > 1e-12, error('Test 3 failed'), end

% ---- Test 4: inconsistent inputs must raise errors ----
bad = {M_.Raw_moments_e, M_.Raw_moments_e};
bad{1}(ie, 1) = 0.1;            % nonzero mean
bad{2}(ie, 2) = 1.5;            % variance different from Sigma_e
for b = 1:2
    M_.Raw_moments_e = bad{b};
    raised = false;
    try
        [info, oo_] = stoch_simul(M_, options_, oo_, {});
        raised = info(1) ~= 0;
    catch
        raised = true;
    end
    fprintf('Test 4.%d: error raised: %d\n', b, raised);
    if ~raised, error('Test 4 failed'), end
end
M_.Raw_moments_e = [];
fprintf('raw_moments_korder ORDER=%d: all tests passed\n', K);
