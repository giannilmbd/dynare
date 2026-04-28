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

The fix is to treat the third unconditional innovation moment as another shock moment available to
the k-order recursions.

For a skewed shock entry

```matlab
M_.Skew_e(row,:) = [i, j, k, skewness]
```

the C++ MEX bridge converts the standardized skewness into a raw third moment using

```cpp
raw_moment = skewness * sqrt(Sigma_e(i,i) * Sigma_e(j,j) * Sigma_e(k,k));
```

For diagonal skew normal shocks this reduces to

```text
E[e_i^3] = skewness_i * variance_i^(3/2)
```

The raw third moments are then inserted into the same moment infrastructure that already stores
second and fourth Gaussian moments.

## Files Changed

The main C++ moment plumbing is in:

- `mex/sources/libkorder/tl/normal_moments.hh`
- `mex/sources/libkorder/tl/normal_moments.cc`

These files now define a `ThirdMoment` record and an overload of `UNormalMoments` that can insert a
third-order moment tensor. The tensor is filled symmetrically over all permutations of `(i,j,k)`.

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

It now reads `M_.Skew_e`, builds the third-moment list, passes it into `Approximation`, and exports
the pure third stochastic derivative as `gsss`.

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

## Caveats

This patch targets the third-order k-order path and the associated pruning/local state-space
simulation path.

It does not yet add a formal Dynare regression test. It also does not update analytical
identification derivatives for `ghs3`; those paths may need separate work if full upstream coverage
is required.
