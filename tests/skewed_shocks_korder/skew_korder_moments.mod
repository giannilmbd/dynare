/*
 * Checks the stochastic constants of a k-th order (pruned) k-order solution
 * when a shock is skew normal (declared with `skew` in the shocks block).
 *
 * With ORDER = K, y = E_t[e(+1)^2 + … + e(+1)^K + e(+1)^2*z(+1)^2]
 * (the last term only for K ≥ 4)
 *                   = E[e^2] + … + E[e^K] + E[e^2]*E[z^2],
 * where e is skew normal and z is an independent Gaussian shock. The K-th
 * order perturbation is exact, and the pure sigma columns of dr.pruning.g_k
 * (which carry the 1/k! factor) must equal, for y:
 *   g_2(sigma^2) = E[e^2], g_3(sigma^3) = E[e^3],
 *   g_4(sigma^4) = E[e^4] + E[e^2]*E[z^2], g_k(sigma^k) = E[e^k] for k ≥ 5.
 * The raw moments of the centred skew normal are obtained by numerical
 * integration of its density, using the (omega, alpha) mapping of
 * csn_update_specification.m.
 * A pruned K-th order simulation must return the constant y in every period.
 *
 * Usage: dynare skew_korder_moments -DSKEW=0.5 -DORDER=4
 */

@#ifndef SKEW
@#define SKEW = 0.5
@#endif
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
@#if SKEW != 0
skew e = @{SKEW};
@#endif
end;

stoch_simul(order=@{ORDER}, pruning, irf=0, periods=0, noprint);

K = @{ORDER};
skew_target = @{SKEW};
ie = strmatch('e', M_.exo_names, 'exact');
iz = strmatch('z', M_.exo_names, 'exact');
var_e = M_.Sigma_e(ie,ie);
var_z = M_.Sigma_e(iz,iz);

%% True raw moments of the centred skew normal e: m(k+1) = E[e^k]
if skew_target == 0
    omega = sqrt(var_e); alpha = 0;
else
    a23 = abs(skew_target)^(2/3);
    delta = sqrt(pi/2*a23/(a23 + ((4-pi)/2)^(2/3)));
    alpha = sign(skew_target)*delta/sqrt(1-delta^2);
    omega = sqrt(var_e/(1-2/pi*delta^2));
end
mu = -omega*alpha/sqrt(1+alpha^2)*sqrt(2/pi);
sn_pdf = @(e) 2/omega*normpdf((e-mu)/omega).*normcdf(alpha*(e-mu)/omega);
m = zeros(K+1,1);
for k = 0:K
    m(k+1) = integral(@(e) e.^k.*sn_pdf(e), -Inf, Inf, 'AbsTol', 1e-14, 'RelTol', 1e-12);
end
fprintf('\nSkew normal e (skewness %g, variance %g), raw moments by numerical integration:\n', skew_target, var_e);
for k = 0:K
    fprintf('  E[e^%d] = %.12f\n', k, m(k+1));
end

%% Targets for the pure sigma terms of y
target = zeros(K,1);
target(2:K) = m(3:K+1);
if K >= 4
    target(4) = target(4) + var_e*var_z;
end

%% Pure sigma columns of the pruning tensors (last unfolded column)
iy = find(oo_.dr.order_var == strmatch('y', M_.endo_names, 'exact'));
g_s = zeros(K,1);
for k = 1:K
    gk = oo_.dr.pruning.(['g_' num2str(k)]);
    g_s(k) = gk(iy, end);
end
fprintf('\nPure sigma terms for y (g_k(sigma^k), 1/k! included):\n');
for k = 1:K
    fprintf('  order %d: %.12f   (target %.12f)\n', k, g_s(k), target(k));
end

%% Pruned K-th order simulation
rng(1);
shocks_sim = randn(50,2)*sqrt(M_.Sigma_e);
y_sim = simult_(M_, options_, oo_.dr.ys, oo_.dr, shocks_sim, K);
iy_decl = strmatch('y', M_.endo_names, 'exact');
fprintf('\nPruned order-%d simulation of y, periods 2..51: min %.12f, max %.12f (target %.12f)\n', ...
        K, min(y_sim(iy_decl,2:end)), max(y_sim(iy_decl,2:end)), sum(target));

%% Verdicts
tol = 1e-8;
fprintf('\nResults:\n');
results = [abs(g_s - target)' < tol, max(abs(y_sim(iy_decl,2:end) - sum(target))) < 1e-6];
names = [arrayfun(@(k) sprintf('order-%d sigma term', k), 1:K, 'UniformOutput', false), ...
         {sprintf('pruned order-%d simulation = constant', K)}];
for i = 1:numel(results)
    if results(i)
        fprintf('  %-40s PASS\n', names{i});
    else
        fprintf('  %-40s FAIL\n', names{i});
    end
end
if ~all(results)
    error('skew_korder_moments: %d of %d checks failed', sum(~results), numel(results));
end
