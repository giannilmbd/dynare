# Markov-switching k-order perturbation (`ms_korder`)

Perturbation solution of order k ≥ 2 for models whose parameters switch with an exogenous Markov chain
s_t ∈ {1, …, S}, its pruned simulation, and the exact moments of the pruned solution.

## Model

Each switching parameter θ is written as θ(s) = θ̄ + σ â(s), with σ the perturbation parameter that also scales
the innovations. The policy of regime s is

    y_t − ȳ = Σ_{i+j+k ≤ K} g_{y^i u^j σ^k}(s_t) (y*_{t−1} − ȳ*)^{⊗i} ⊗ u_t^{⊗j} σ^k / (i! j! k!),

where y* are the state variables and ȳ is the steady state at θ̄. Derivatives without σ do not depend on the regime.
Derivatives with σ do, including g_σ(s) at first order.

## Writing the model

Declare each switching parameter as an endogenous variable with steady state 0 and a dummy equation tagged with
its name. The right-hand side of the dummy equation must be a parameter set to 0. A literal `0` is substituted
into the other equations by the preprocessor and removes the variable from the model.

```
var k c aA;
varexo e;
parameters ms_zero Abar ...;
ms_zero = 0;
model;
  c + k = (Abar + aA)*k(-1)^zeta + (1 - delta)*k(-1);
  ...
  [ms='aA'] aA = ms_zero;
end;
steady_state_model; ... aA = 0; end;
stoch_simul(order=3, k_order_solver, irf=0);
```

Then:

```
regimes = struct('P', [0.95 0.05; 0.20 0.80], ...          % P(s, s') = Pr(s_{t+1} = s' | s_t = s)
                 'values', struct('aA', [0 -0.10]), ...    % â(s) of each tagged variable
                 'names', {{'normal', 'bad'}});
oo_ = ms_korder(M_, options_, oo_, regimes, struct('periods', 1000));
```

`ms_korder` stores its output in `oo_.ms`:

| field | content |
|---|---|
| `oo_.ms.solution` | `g{s}.y<i>u<j>s<k>`: unfolded derivatives of regime s. Rows are in `oo_.dr.order_var` order; columns are `kron(y*^{⊗i}, u^{⊗j})` with y* = `solution.state_var`. |
| `oo_.ms.moments` | `pi` (ergodic regime probabilities), `mean`, `mean_regime` (E[y \| s]), `var`, `second_regime{s}` (E[(y − ȳ)(y − ȳ)′ \| s]); declaration order |
| `oo_.ms.simulation` | `y` (levels, declaration order), regime path `s`, innovations `u` |

The building blocks can also be called directly:

| function | purpose |
|---|---|
| `ms_korder_solve(M_, options_, oo_.dr, P, ahat)` | solution; `ahat` is [switching variables × S] in the order of the tags |
| `ms_korder_simul(ms, s_path, u_path, pruning, y0dev)` | simulation along a given regime path, pruned or not |
| `ms_korder_moments(ms, Sigma_u, second, solver)` | exact moments of the pruned solution |
| `ms_korder_eval(ms, s, ystate, u, sigma)` | evaluation of the policy of regime s |
| `ms_markov_chain(P, T, s0, seed)` | draw of a regime path |

## Method

**Solution.** In the stock k-order solver, next period's policy is composed with this period's. With regimes,
the composition depends on the regime pair (s, s′). The MEX `k_order_perturbation_ms` (`libkorder/kord/korder_ms.{hh,cc}`)
keeps one stock `KOrder` object per pair and takes the expectation of the equations of regime s over s′ with
weights P(s, s′). The derivatives with σ satisfy systems coupled across regimes. They are decoupled with the
eigen-decomposition P′W = WΛ and solved with the stock Sylvester solver and the stock σ^k matrix, each scaled
by an eigenvalue λ. The derivation is in `doc/ms_korder_algebra.md`.

**Pruning.** The pruned solution is the series expansion of the solution in σ (Lombardo and Uhlig 2018,
Theorem 1). The term g · x^(a_1) ⊗ … ⊗ x^(a_i) ⊗ u^{⊗j} σ^k goes to layer a_1 + … + a_i + j + k, with g evaluated
at s_t. This is the layer rule of Dynare's pruned simulator.

**Moments.** The unknowns are E[F_b(x^(b_1) ⊗ … ⊗ x^(b_m)) 1{s_t = s}] for sorted layer tuples b with
b_1 + … + b_m ≤ 2K. F_b maps each group of equal layers to its symmetric coordinates, which generalises vech. The
loading matrices on the lagged blocks follow from the layer laws by Kronecker algebra. The system is block
triangular: by weight ascending, and within a weight by number of factors descending. Each diagonal block is
M − Φ^{[b]} M P = R, with Φ = g_y. Using P V = V diag(λ) it splits into (I − λ Φ^{[b]}) x = r for each eigenvalue,
solved by Smith doubling with Kronecker-vector products. The options `'gensylv'` (Dynare's Schur-recursion solver
on the unfolded block) and `'lu'` give independent solutions for checking.

## Requirements and limitations

- All switching parameters are perturbed (θ(s) = θ̄ + σ â(s)), so ȳ, g_y and g_u do not depend on the regime.
  The partition perturbation of Foerster, Rubio-Ramírez, Waggoner and Zha (2016) keeps parameters that do not
  affect the steady state out of the perturbation. That makes first-order coefficients regime dependent and is
  not implemented here.
- P must have real eigenvalues and independent eigenvectors.
- Moments assume Gaussian innovations. The solution reads the same shock moments as `k_order_perturbation`
  (including `M_.Skew_e`), but only Gaussian shocks have been tested with regimes.
- No autocovariances or impulse responses.

## Tests

`tests/ms_korder/`:

- `ms_one_regime.mod`: with one regime, the derivatives, pruned simulation and pruned moments equal stock Dynare.
- `ms_linear_exact.mod`: three regimes, linear model with a forward-looking variable, against its exact solution
  and moments.
- `ms_two_regimes.mod`: three moment solvers agree; regime-conditional moments aggregate to the unconditional
  ones; end-to-end run of `ms_korder`.

## References

- Foerster, A., J. F. Rubio-Ramírez, D. F. Waggoner and T. Zha (2016), Perturbation methods for Markov-switching
  dynamic stochastic general equilibrium models, *Quantitative Economics* 7(2), 637–669.
- Lombardo, G. and H. Uhlig (2018), A theory of pruning, *International Economic Review* 59(4), 1825–1836.
- Costa, O. L. V., M. D. Fragoso and R. P. Marques (2005), *Discrete-Time Markov Jump Linear Systems*, Springer.
