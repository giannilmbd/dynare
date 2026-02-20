heterogeneity_dimension households;

var(heterogeneity=households)
   b
   a
   c
   Va
   Vb
   u
;

varexo(heterogeneity=households)
   e
;

var
    piw
    psiw
    rb
    ra
    tax
    i
    psip
    I
    Q
    w
    N
    K
    div
    p
    pi
    mc
    r
    Y
;

varexo
    rstar
    markup
    G
    beta
    Z
    rinv_shock
    markup_w
;

shocks;
    var rstar; stderr 0.01;
    var markup; stderr 0.01;
    var G; stderr 0.01;
    var beta; stderr 0.01;
    var Z; stderr 0.01;
    var rinv_shock; stderr 0.01;
    var markup_w; stderr 0.01;
end;

parameters
    kappap
    alpha
    epsI
    muw
    phi
    omega
    Bg
    pshare
    delta
    kappaw
    frisch
    mup
    vphi
    eis
    chi0
    chi1
    chi2
    Z_ss
    beta_ss
    r_ss
    G_ss
;


model(heterogeneity=households);
   // Euler equations
   (beta_ss+beta) * Vb(+1) - c^(-1/eis) = 0 ⟂ b >= 0;
   (beta_ss+beta) * Va(+1) - c^(-1/eis)*(1 + chi1 * sign(a-(1+ra)*a(-1)) * (abs(a-(1+ra)*a(-1))/((1+ra)*a(-1)+chi0))^(chi2-1)) = 0 ⟂ a >= 0;

   // Budget constraint
   (1 + ra) * a(-1) + (1 + rb) * b(-1) - (chi1 / chi2) * abs(a-(1+ra)*a(-1))^chi2 * ((1+ra)*a(-1)+chi0)^(1-chi2) + (1-tax) * w * N * e - c - a - b;

   // Effective labor definition
   u = e * c^(-1/eis);

   // Envelope conditions
   Va = (1 + ra)*(1 - (chi1/chi2) * ( - chi2 * sign(a-(1+ra)*a(-1)) * (abs(a-(1+ra)*a(-1))/((1+ra)*a(-1)+chi0))^(chi2-1) + (1-chi2) * (abs(a-(1+ra)*a(-1))/((1+ra)*a(-1)+chi0))^chi2 )) * c^(-1/eis);
   Vb = (1 + rb) * c^(-1/eis);
end;

model;
   // NKPC
   kappap * (mc - 1 / mup) + Y(+1) / Y * log(1 + pi(+1)) / (1 + r(+1)) + markup - log(1 + pi);

   // Equity price
   div(+1) + p(+1) - p * (1 + r(+1));

   // Production function
   N = (Y / (Z_ss+Z) / K(-1) ^ alpha) ^ (1 / (1 - alpha));

   // Labor demand
   mc = w * N / (1 - alpha) / Y;

   // Tobin's Q
   (K / K(-1) - 1) / (delta * epsI) + 1 - Q;

   // Valuation equation
   alpha * (Z_ss+Z(+1)) * (N(+1) / K) ^ (1 - alpha) * mc(+1) - (K(+1) / K -
   (1 - delta) + (K(+1) / K - 1) ^ 2 / (2 * delta * epsI)) + K(+1) / K * Q(+1) - (1 + r(+1) + rinv_shock) * Q;

   // Price adjustment cost
   mup / (mup - 1) / 2 / kappap * log(1 + pi) ^ 2 * Y - psip;

   // Aggregate investment
   K - (1 - delta) * K(-1) + K(-1) * (K / K(-1) - 1) ^ 2 / (2 * delta * epsI) - I;

   // Resource constraint
   Y - w * N - I - psip - div;

   // Taylor rule
   rstar + r_ss + phi * pi - i;

   // Fiscal policy
   (r * Bg + G_ss + G) / w / N - tax;

   // Returns on liquid asset
   r - omega - rb;

   // Returns on illiquid asset
   pshare * (div + p) / p(-1) + (1 - pshare) * (1 + r) - 1 - ra;

   // Fisher equation
   1 + i(-1) - (1 + r) * (1 + pi);

   // Wage inflation
   (1 + pi) * w / w(-1) - 1 - piw;

   // Wage adjustment cost
   muw / (1 - muw) / 2 / kappaw * log(1 + piw) ^ 2 * N - psiw;

   // Wage NKPC
   kappaw * (vphi * N ^ (1 + 1 / frisch) - (1 - tax) * w * N * SUM(u) / muw) + (beta_ss+beta) * log(1 + piw(+1)) + markup_w - log(1 + piw);

   // Total asset market clearing
   p + Bg - SUM(b) - SUM(a);
end;

load 'hank_2a.mat';

param_names = fieldnames(steady_state.params);
for i=1:numel(param_names)
   param = param_names{i};
   set_param_value(param, steady_state.params.(param));
end

% Initialize steady state once for all tests
heterogeneity_load_steady_state(filename = hank_2a);
assert(isfield(oo_, 'heterogeneity'), 'oo_.heterogeneity field missing');
assert(isfield(oo_.heterogeneity, 'steady_state'), 'oo_.heterogeneity.ss field missing');
assert(isfield(oo_.heterogeneity, 'sizes'), 'oo_.heterogeneity.sizes field missing');
assert(isfield(oo_.heterogeneity, 'mat'), 'oo_.heterogeneity.mat field missing');
assert(isfield(oo_.heterogeneity, 'indices'), 'oo_.heterogeneity.indices field missing');

% Solve once for all tests
heterogeneity_solve;

testFailed = 0;
testResults = [];  % Array to collect all test results

verbatim;
    % Test heterogeneity.simulate function (reuse initialized and solved model)
    [testFailed, testResults] = test_simulate_stochastic(M_, options_, oo_, steady_state, testFailed, testResults);

    % Test permutation handling
    [testFailed, testResults] = test_permutation(M_, options_, oo_, steady_state, testFailed, testResults);

    % Print test summary
    print_test_summary(testResults);

    if testFailed > 0
        error('Some unit tests failed!');
    end
end;