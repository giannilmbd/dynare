# Third Moment Hack

This note documents the local changes made to make Dynare's third-order k-order solver account
for skewed shocks when computing stochastic constants.

## Problem

Dynare already accounts for second moments of innovations in the forward-looking part of the
third-order solution. In the k-order C++ path, this is done through the normal-moment containers and
objects such as `gss`/`ghs2`.

The missing piece was that the code assumed unconditional shock moments only existed at even orders.
That is correct for Gaussian innovations, but it is not correct once `M_.Skew_e` introduces shocks
with nonzero third moments. Consequently, terms involving `E[u^3]` were dropped from the forward
solution. This affected the stochastic constant, especially for forward-looking variables.

## Main Idea

Skewed shocks declared with `skew` in the `shocks` block are independent skew normal variables, as in
`matlab/distributions/closed_skew_normal/csn_update_specification.m`. When `M_.Skew_e` has a
non-zero entry, the C++ MEX bridge builds the moment tensors of every order `d = 2, ..., order` from
the raw moments of each shock:

```text
E[e_{i1} ... e_{id}] = prod_i E[e_i^{n_i}],   n_i = number of occurrences of i in (i1, ..., id)
```

The raw moments of a zero-mean skew normal shock with variance `sigma^2` and standardized skewness
`gamma` come from the representation

```text
e = omega * ( delta * (|U0| - sqrt(2/pi)) + sqrt(1 - delta^2) * U1 ),   U0, U1 iid N(0,1)
```

with `delta` and `omega` mapped from `(sigma^2, gamma)` exactly as in `csn_update_specification.m`.
Gaussian shocks are the case `delta = 0`. A moment tensor is inserted only if it has a non-zero
element, so odd-order tensors exist only when some shock is skewed. When no shock is skewed, the
original Gaussian (Isserlis) construction is used unchanged.

The third moment is therefore `E[e_i^3] = gamma_i * sigma_i^3`, and the fourth moment is the skew
normal one, `E[e_i^4] = sigma_i^4 * (3 + 2(pi-3) (delta sqrt(2/pi))^4 / (1 - 2 delta^2/pi)^2)`, not the
Gaussian `3 sigma_i^4`. Higher orders follow the same construction.

Restrictions enforced by the MEX bridge (as in `csn_update_specification.m`):

- co-skewness entries (`i`, `j`, `k` not all equal) are rejected;
- skewed shocks cannot be correlated (`Sigma_e` must be diagonal);
- `|gamma| < 0.99527...`, the skew normal bound.

## Files Changed

The main C++ moment plumbing is in:

- `mex/sources/libkorder/tl/normal_moments.hh`
- `mex/sources/libkorder/tl/normal_moments.cc`

They define `skewNormalRawMoments(variance, skewness, maxdim)` and a `UNormalMoments` constructor
that builds the moment tensors of independent shocks from per-shock raw moments.

The k-order recursion changes are in:

- `mex/sources/libkorder/kord/korder.hh`
- `mex/sources/libkorder/kord/korder.cc`
- `mex/sources/libkorder/kord/korder_stoch.hh`
- `mex/sources/libkorder/kord/approximation.hh`
- `mex/sources/libkorder/kord/approximation.cc`

The original logic used checks like "is this order even?" to decide whether a shock moment exists.
Those checks were replaced by `hasMoment(...)`, which asks the moment container whether the relevant
unconditional moment is present. This lets the solver keep the Gaussian behavior unchanged while
allowing third moments when skewed shocks are declared.

The MATLAB MEX bridge is in:

- `mex/sources/k_order_perturbation/k_order_perturbation.cc`

It reads `M_.Skew_e`, checks the restrictions above, builds the shock moment container, passes it
into `Approximation`, and exports the pure third stochastic derivative as `gsss`.

The MATLAB decision-rule interface is in:

- `matlab/stochastic_solver/k_order_pert.m`

It maps `dyn_derivs.gsss` to `dr.ghs3`.

Third-order simulation/pruning support was added in:

- `matlab/stochastic_solver/simult_.m`
- `mex/sources/local_state_space_iterations/local_state_space_iteration_3.f08`
- `matlab/nonlinear-filters/iterate_law_of_motion.m`
- `matlab/estimation/non_linear_dsge_likelihood.m`
- `matlab/estimation/online/online_auxiliary_filter.m`

The Fortran local state-space MEX accepts `ghs3` as an optional extra argument. Old 17-argument calls
still work and behave as before with `ghs3 = 0`.

## What `ghs3` Means

`ghs3/6` is the third-order stochastic constant coming from forward expectations and the third
unconditional innovation moment. It is not the same object as the contemporaneous `ghuuu` term.

For example, in the model

```dynare
u = e1^3;
y = u(+1);
```

with `Var(e1)=1` and `E[e1^3]=0.5`, the relevant output in decision-rule order is:

```text
ghs3/6
    0.5000   % y
         0   % dummy state
         0   % u

ghuuu*E[e1^3]/6
         0   % y
         0   % dummy state
    0.5000   % u
```

So the forward-looking variable `y` receives the skewness contribution through `ghs3`, while the
contemporaneous variable `u` receives it through `ghuuu`.

## Validation Trick

To verify that the stochastic constant `g_0` receives both second- and third-moment forward terms,
use:

```dynare
v = e1^2;
u = e1^3;
y = v(+1) + u(+1);
```

With `Var(e1)=1` and `E[e1^3]=0.5`, the expected decomposition is:

```text
ghs2/2   = 1.0000 for y
ghs3/6   = 0.5000 for y
g_0      = 1.5000 for y
```

This is what the patched solver returns.

## Orders above 3

At `order >= 4`, `stoch_simul` simulates through `k_order_simul`; with `pruning` this is the
arbitrary-order pruned state space of `mex/sources/libkordersim/simulation.f08`. It uses the full
decision-rule tensors in `dr.pruning.g_k`, including the pure sigma columns (`g_sss`, `g_ssss`, ...),
so the skewed-shock moments enter without further changes. The pure sigma term of order `k`, with the
`1/k!` factor included, is the last column of `dr.pruning.g_k`.

Pruning at order above 3 is not available for particle filtering, estimation, or theoretical moments.

## Regression test

`tests/skewed_shocks_korder/skew_korder_moments.mod` solves

```dynare
p2 = e^2; ... pK = e^K;  q = e^2*z^2;
y = p2(+1) + ... + pK(+1) + q(+1);
```

with `e` skew normal and `z` Gaussian, for which the order-`K` solution is exact. It compares the pure
sigma terms of `y` with the skew normal raw moments obtained by numerical integration of the density,
and checks that a pruned order-`K` simulation returns that constant. Usage:

```matlab
dynare skew_korder_moments -DSKEW=0.5 -DORDER=4
```

It passes for `SKEW` in {0, 0.5, -0.9} and `ORDER` in {3, 4, 5}; the logs are stored next to it.

## Caveats

Analytical identification derivatives do not include `ghs3` or the non-Gaussian moments. Other
k-order MEX files (e.g. `k_order_welfare`) still use Gaussian moments.

## User-supplied raw moments (branch `raw-moments-shocks`)

`M_.Raw_moments_e` (optional, set in MATLAB before `stoch_simul`) is an `exo_nbr × D` matrix with entry
`(i, d) = E[e_i^d]`, `d = 1..D`, `D >= order`. A row of `NaN` leaves shock `i` to the Gaussian or skew normal
default. A row that is given must satisfy `E[e_i] = 0` and `E[e_i^2] = M_.Sigma_e(i,i)` (relative tolerance
1e-10); the shock must not also appear in `M_.Skew_e`, and it must be uncorrelated with the other shocks. The MEX
bridge passes the rows to the same `UNormalMoments(maxdim, raw_moments)` constructor used for skew normal shocks,
so any independent non-Gaussian shock with finite moments up to the solution order can be used.

`matlab/distributions/mixture_normal_raw_moments.m` returns the raw moments of a mean-zero finite mixture of
normals (one row of `M_.Raw_moments_e`).

Shocks drawn by Dynare's own simulation routines are still Gaussian or skew normal: with user-supplied moments,
simulate by passing the draws to `simult_`.

Test: `tests/raw_moments_korder/raw_moments_korder.mod` (`-DORDER=3,4,5`): mixture moments vs numerical
integration; exact pure-sigma terms of an exactly solvable model; skew normal moments supplied as raw moments
reproduce the `Skew_e` solution (max difference 3.6e-15); Gaussian raw moments reproduce the default solution
exactly; inconsistent inputs raise errors.
