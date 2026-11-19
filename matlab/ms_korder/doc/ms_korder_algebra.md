---
title: "Markov-switching k-order perturbation in Dynare: algebra, pruning and moments"
geometry: margin=2.5cm
---
# 1. Model

$$\mathbb{E}_t\, f\big(y_{t+1},\, y_t,\, y_{t-1},\, u_t,\, a_{t+1},\, a_t\big) = 0 ,$$

- $y_t \in \mathbb{R}^{n_y}$: endogenous variables in deviation from the deterministic steady state $\bar y$.
- $u_t$: iid innovations, $\mathbb{E}u = 0$, $\mathbb{E}\,u\otimes u = \operatorname{vec}\Sigma$.
- $s_t \in \{1,\dots,S\}$: exogenous Markov chain with transition matrix $P$, $P_{ss'} = \Pr(s_{t+1}=s' \mid s_t=s)$.
- $a_t = \sigma\,\hat a(s_t) \in \mathbb{R}^{n_a}$: the switching parameters in deviation from $\bar\theta$, so that
  $\theta(s) = \bar\theta + \sigma \hat a(s)$. The parameter $\sigma$ also scales future innovations,
  $u_{t+1} \to \sigma u_{t+1}$, and the model of interest is $\sigma = 1$.

At $\sigma = 0$ every regime has the same deterministic steady state. The policy is
$y_t = g(y_{t-1}, u_t, \sigma; s_t)$; derivatives at $(0,0,0)$ are written without arguments.

Notation. For a function $x(s)$ on regimes, $(Px)(s) = \sum_{s'} P_{ss'}\, x(s')$, which is
$\mathbb{E}_t\, x(s_{t+1})$ given $s_t = s$. In matrix form, with $X = [x(1)\ \cdots\ x(S)]$, it is $XP^\top$.
First derivatives of $f$ (with respect to $y_{t+1}, y_t, y_{t-1}, u_t, a_{t+1}, a_t$):
$A, B, C, D, E_a, F_a$. Second derivatives: $f_2$, acting on $w \otimes w$ with
$w = (y_{t+1}, y_t, y_{t-1}, u_t, a_{t+1}, a_t)$.

# 2. First order

Guess $y_t = G y_{t-1} + H u_t + \sigma c(s_t)$. Then $\mathbb{E}_t y_{t+1} = G y_t + \sigma (Pc)(s_t)$, and
substituting:

$$\begin{aligned}
y_{t-1}:&\quad (AG+B)G + C = 0,\\
u_t:&\quad (AG+B)H + D = 0,\\
1:&\quad (AG+B)\,c(s) + A\,(Pc)(s) = -\big[E_a (P\hat a)(s) + F_a \hat a(s)\big] \equiv -r(s).
\end{aligned}$$

**(i)** $G$ and $H$ are those of the constant-parameter model: Dynare's `ghx`, `ghu`.

**(ii)** The regime enters through $c(s) = g_\sigma(s)$, the solution of an $S$-regime coupled linear system. In
matrix form, with $M = AG+B$,
$$M\,\mathcal{C} + A\,\mathcal{C}P^\top = -\mathcal{R}, \qquad
(I_S \otimes M + P \otimes A)\operatorname{vec}\mathcal{C} = -\operatorname{vec}\mathcal{R}.$$
The solver uses the stacked form: one LU factorisation of $I_S \otimes M + P \otimes A$, of order $S\,n_y$, valid
for any $P$.

**(iii) Forward solution.** Iterating $c = -M^{-1}r - M^{-1}A\,Pc$ gives, with $F = -M^{-1}A$,
$$c(s_t) = -\sum_{i\ge 0} F^i M^{-1}\Big[E_a\, \mathbb{E}_t \hat a(s_{t+i+1}) + F_a\, \mathbb{E}_t \hat a(s_{t+i})\Big],
\qquad \mathbb{E}_t\, \hat a(s_{t+i}) = (P^i \hat a)(s_t).$$
This is the present discounted sum of expected future exogenous variables. The Markov chain enters at first order
only through $\mathbb{E}_t \hat a(s_{t+i})$.

**(iv) Eigen-components.** If $P$ is diagonalisable, write $P^\top = W\Lambda W^{-1}$ and $\tilde{\mathcal C} = \mathcal{C}W$. Then
$$(M + \lambda_m A)\,\tilde c_m = -(\mathcal{R}W)_m, \qquad m = 1,\dots,S.$$
With an AR(1) process $\mathbb{E}_t a_{t+i} = \rho^i a_t$ and $c = \Gamma a_t$, $(M + \rho A)\Gamma = -(\rho E_a + F_a)$.
The chain is the sum of $S$ such components, one per eigenvector of $P$, with persistence $\lambda_m$.
$\lambda = 1$ gives $M + A = AG + B + A$, which is Dynare's matrix for the $\sigma$-derivatives (`MatrixS`).

**(v) Centring.** If $\bar\theta$ is the ergodic mean, $\sum_s \pi(s)\hat a(s) = 0$ with $\pi^\top P = \pi^\top$. Then
$\sum_s \pi(s) c(s) = 0$: $c$ has no unconditional intercept. Proof: multiply the matrix equation by $\pi$ on the right;
$P^\top\pi = \pi$ gives $(M+A)\,\mathcal{C}\pi = -(E_a + F_a)\hat{\mathcal A}\pi = 0$, and $M + A$ is invertible.


# 3. Second order

$$y_t = G y_{t-1} + H u_t + \sigma c(s_t)
+ \tfrac12\big[G_{yy}(y_{t-1}\otimes y_{t-1}) + 2G_{yu}(y_{t-1}\otimes u_t) + G_{uu}(u_t\otimes u_t)\big]
+ \sigma K(s_t)\,y_{t-1} + \sigma L(s_t)\,u_t + \tfrac12\sigma^2 d(s_t),$$
with $K = g_{y\sigma}$, $L = g_{u\sigma}$, $d = g_{\sigma\sigma}$.

The $\sigma$-free coefficients $G_{yy}, G_{yu}, G_{uu}$ come from equations evaluated at $\sigma = 0$, where no regime
appears. The constant-parameter solution (Dynare's `ghxx`, `ghxu`, `ghuu`) solves them.

Order-$\sigma$ part of the arguments of $f$, given $s_t = s$ and $s_{t+1} = s'$:
$$w_\sigma(s,s') = \big(Gc(s) + c(s') + Hu_{t+1},\ c(s),\ 0,\ 0,\ \hat a(s'),\ \hat a(s)\big).$$
Its parts linear in $y_{t-1}$ and $u_t$ are $W_y = (G^2, G, I, 0, 0, 0)$ and $W_u = (GH, H, 0, I, 0, 0)$.

**$K(s)$: terms in $\sigma\, y_{t-1}$.**
$$(AG+B)\,K(s) + A\,(PK)(s)\,G = -\Big[A\,G_{yy}\big(G\otimes c(s)\big) + f_2\big(W_y \otimes \bar w_\sigma(s)\big)\Big],$$
where $\bar w_\sigma(s) = \sum_{s'}P_{ss'}w_\sigma(s,s')$; the bracket is linear in $w_\sigma$, so only $Pc$ and
$P\hat a$ enter. This is a generalised Sylvester equation coupled across regimes. Stacked over regimes,
$\mathcal K = [K(1); \dots; K(S)]$,
$$(I_S\otimes M)\,\mathcal K + (P\otimes A)\,\mathcal K\,G = \mathcal R,$$
which is the equation of Dynare's `recover_ys` with $n_y$ replaced by $S n_y$. $A$ acts only on the forward-looking
rows, so ordering the unknowns with the forward-looking rows of all regimes last keeps the zero-column structure
that Dynare's `GeneralSylvester` exploits.

**$L(s)$: terms in $\sigma\, u_t$.**
$$(AG+B)\,L(s) = -\Big[A\,(PK)(s)\,H + A\,G_{yy}\big(H\otimes c(s)\big) + f_2\big(W_u\otimes \bar w_\sigma(s)\big)\Big].$$
There is no coupling: $L(s)$ reaches $t+1$ only through today's regime. Each regime is one solve with $M$, which is
Dynare's `recover_yus` (`matA`).

**$d(s)$: terms in $\sigma^2$.**
$$\begin{aligned}
(AG+B)\,d(s) + A\,(Pd)(s) = -\Big[\, &A\,G_{yy}\big(c(s)\otimes c(s)\big) + 2A\,(PK)(s)\,c(s) + A\,G_{uu}\operatorname{vec}\Sigma\\
&+ \textstyle\sum_{s'} P_{ss'}\ \mathbb{E}_{u}\, f_2\big(w_\sigma(s,s')\otimes w_\sigma(s,s')\big)\Big].
\end{aligned}$$
The matrix is the one of first order, $I_S \otimes M + P \otimes A$ (Dynare's `recover_s` with $S n_y$ rows). The last term contains
$$\sum_{s'}P_{ss'}\, c(s')\otimes c(s'),\qquad \sum_{s'}P_{ss'}\, c(s')\otimes \hat a(s'),\qquad
\sum_{s'}P_{ss'}\, \hat a(s')\otimes \hat a(s'),$$
the conditional second moments of the chain one step ahead. They differ from $(Pc)\otimes(Pc)$ etc. by the
conditional variance of the regime transition: the precautionary effect of regime risk. The innovation part gives
$(H\otimes H)\operatorname{vec}\Sigma$ in the $y_{t+1}$ block, as in the constant model.

**Order of solution:** $G_{yy}, G_{yu}, G_{uu}$, then $K$ (needs $c$), then $L$ (needs $K$), then $d$ (needs
$c, K$). This is the sequence of Dynare's `KOrder::performStep`.

# 4. What changes with the order

| | Order 1 | Order 2 | Order $k$ |
|---|---|---|---|
| $\sigma$-free coefficients | Dynare | Dynare | Dynare |
| Coupled equations ($P$ enters the matrix) | $c$ | $K$, $d$ | every $g_{y^i\sigma^j}$, $g_{\sigma^j}$ |
| Uncoupled, per regime | none | $L$ | every $g_{y^i u^l \sigma^j}$, $l\ge1$ |
| Moments of the chain used | $P\hat a$ | $P$ of products of two regime functions | $P$ of products of up to $k$ |

Only one-step conditional moments enter: $\sum_{s'}P_{ss'}$ of products of regime functions ($\hat a$, and the
lower-order coefficients evaluated at $s'$). The forward sums over horizons ((iii) above) are produced by the
coupled solves, as in Dynare. This is what keeps order $k$ tractable: no moments of the chain at horizons beyond one
appear in the solution step. Multi-horizon moments appear only in the unconditional moments of the pruned
solution (section 6).

# 5. Order $k$: the equation for each symmetry

Let $X(s) = g_{y^i u^j \sigma^l}(s)$ with $i + j + l = k$, written as in Dynare: $y$ means $y^*_{t-1}$, $u$ means
$u_t$. Next period, $y_{t+1} = g(g^*(y^*_{t-1}, u_t, \sigma; s), \sigma u_{t+1}, \sigma; s_{t+1})$. Faa di Bruno for
this composition has one outer derivative (of $g(\cdot; s')$) per term, times inner derivatives (of $g^*(\cdot; s)$).
The unknown enters in three places:

1. $f_y\, X(s)$, the current policy;
2. $f_{y^{**}_+}\, g^{**}_{y^*}\, X^*(s)$, the inner function at the current regime (outer at first order, which is
   regime-free);
3. $f_{y^{**}_+} \sum_{s'} P_{ss'} X^{**}(s') \big(g^*_y\big)^{\otimes i}$: the outer function at $s'$ with every inner
   factor at first order. This term exists only if $j = 0$. A current-shock index $u_t$ reaches the outer function
   only through the inner $g^*_u$. The outer derivative is then $g_{y^{i+j}\sigma^l}(s')$, a different symmetry that
   is solved earlier.

Hence, with $M = f_y + f_{y^{**}_+} g^{**}_{y^*}[\text{sel}^*] = $ Dynare's `matA`:

| Symmetry | Equation | Dynare routine, MS version |
|---|---|---|
| $y^i u^j$ ($l = 0$) | regime-free | `recover_y`, `recover_yu`: unchanged |
| $y^i\sigma^l$, $i \ge 1$, $l \ge 1$ | $M X(s) + A\sum_{s'}P_{ss'}X(s')(g^*_y)^{\otimes i} = R(s)$ | `recover_ys`: one Sylvester solve of order $S n_y$ |
| $y^i u^j\sigma^l$, $j \ge 1$, $l \ge 1$ | $M X(s) = R(s)$ | `recover_yus`: per regime, same `matA` |
| $\sigma^l$ | $M X(s) + A\sum_{s'}P_{ss'}X(s') = R(s)$ | `recover_s`: one LU of $I_S \otimes M + P \otimes A$ |

Here $A$ means $f_{y^{**}_+}$ with the selection of forward-looking variables.

**Right-hand sides** contain only symmetries solved earlier in `performStep`'s order. Two kinds of regime dependence
enter them:

- Faa di Bruno terms of the composition have one outer factor at $s'$, so they need $\sum_{s'}P_{ss'}$ of one
  derivative tensor: linear in $P$.
- Higher derivatives of $f$ multiply several next-period arguments ($y_{t+1}$, $a_{t+1}$), all at $s'$. These need
  $\sum_{s'} P_{ss'}$ of products of up to $k$ regime functions, the one-step conditional moments of the chain.

Neither involves horizons beyond one.

**Inner first-order terms in $\sigma$.** In the constant model $g^*_\sigma = 0$, so $\sigma$ reaches the outer
function only through its own $\sigma$ and $u' = \sigma\xi$ arguments. Here $g^*_\sigma(s) = c^*(s) \ne 0$. The
outer derivatives with respect to $y^*$ are then also contracted with $c^*(s)$, which adds terms such as
$g_{y^{i+1}}(c^*(s) \otimes \cdots)$ to the right-hand sides; all are known when needed. Dynare's `fillG` and
`faaDiBrunoG` must therefore keep the $\sigma$-column of the inner stack ($g^*_\sigma$, regime-dependent). Stock
Dynare omits it: "g_σ is zero by construction and we do not insert it to the container ... since g_σ is zero, so is
G_σ" (`libkorder/kord/korder.cc`, KOrder constructor). The MS solver inserts $g_\sigma(s)$ and $G_\sigma(s)$ after
the first-order solve and sets `GContainer::inner_sigma_zero = false` (`tl/stack_container.hh`; the default `true`
keeps the stock behaviour).

**Other departures from the stock step.** The stock step skips symmetries with an odd $\sigma$ order
(`hasMoment(k)` in `recover_ys`, `recover_yus`, `recover_s`, and `hasMoment(k-m)` in `fillG`), because they vanish
without regimes. With regimes they do not vanish, and the MS step solves them. After a $y^i\sigma^j$ or $\sigma^j$
symmetry is solved, $G$ is updated with the new outer tensor only. The stock update with the whole $g^{**}$
container would add a second time the terms $g^{**}_{y^{i+m}\sigma^{j-m}}(c^*)^{\otimes m}$ that the conditional
$G$ already holds. These terms vanish in the stock model because $c^* = 0$.

# 6. Pruning

Write $x^{(a)}_t$, $a = 1,\dots,K$, for the layers of the state variables, and $\hat y_t = y_t - \bar y$. The pruned
solution of order $K$ is
$$x^{(q)}_t = \sum \frac{1}{i!\,j!\,k!}\, g^*_{y^i u^j \sigma^k}(s_t)\,
\big(x^{(a_1)}_{t-1}\otimes\cdots\otimes x^{(a_i)}_{t-1}\otimes u_t^{\otimes j}\big)\sigma^k,$$
where the sum runs over $i + j + k \le K$ and ordered compositions $a_1 + \cdots + a_i = q - j - k$, $a_r \ge 1$.
For $i = 0$ the term belongs to layer $q = j + k$. The same rule with all rows of $g$ gives the layers of $\hat y_t$,
and $\hat y_t = \sum_q \hat y^{(q)}_t$. It is the layer rule of Dynare's pruned simulator
(`libkordersim/simulation.f08`) with coefficients evaluated at $s_t$. By Lombardo and Uhlig (2018, Theorem 1),
$x^{(q)}_t$ is the $q$-th term of the series expansion of the solution in $\sigma$, with the initial state of order
$\sigma$. With regimes, $g_\sigma(s_t)$ lies in layer 1 and $g_{\sigma\sigma}(s_t)/2$ in layer 2.

# 7. Exact moments of the pruned solution

**Unknowns.** For a sorted tuple of layers $b = (b_1 \le \cdots \le b_m)$ with weight $|b| = b_1 + \cdots + b_m$,
$$\mu_b(s) = \mathbb{E}\big[F_b\,(x^{(b_1)}_t\otimes\cdots\otimes x^{(b_m)}_t)\,\mathbf 1\{s_t = s\}\big].$$
$F_b$ maps each group of $p$ equal layers to the $\binom{n+p-1}{p}$ symmetric coordinates of its $p$-fold
Kronecker power (sorted multi-indices; $p = 2$ is vech). $D_b$ (duplication) is its right inverse on symmetric
vectors. Means of $\hat y$ need $|b| \le K$, second moments $|b| \le 2K$. The empty tuple has $\mu_\emptyset = \pi$,
the ergodic distribution of the chain.

**Recursion.** Substitute the layer laws of section 6 into each factor and expand. Each term picks one law term per
factor. With $\nu_c = \mu_c P$ (moments of the lagged block given the current regime) and Gaussian
$\mathbb{E}\,u^{\otimes J}$ (Isserlis),
$$\mu_b(\cdot, s') = F_b \sum_{\text{choices}} \Big(\bigotimes_f G_f(s')\Big)\,\Pi\,
\big(D_c\,\nu_c(\cdot, s') \otimes \mathbb{E}\,u^{\otimes J}\big),$$
where $c$ is the sorted concatenation of the lagged layers of the chosen terms, $J$ the total power of $u$, and $\Pi$
the permutation that puts every factor's lagged and innovation modes next to its coefficient tensor $G_f$. The
choices that differ only by the order of the factors within a group of equal layers give the same result up to a
permutation of the factors. Only one is computed, multiplied by the number of distinct arrangements, and
$F_b$ is applied to the group-symmetrised sum: $F_b\,\mathrm{Sym} = \mathrm{diag}(\text{counts})^{-1} D_b^\top$.

**Block triangular structure.** A lagged block $c$ has $|c| \le |b|$. Equality requires $j = k = 0$ in every
factor. Then each factor term has $i \ge 1$ lagged layers summing to $b_f$, so $c$ has at least $m$ factors, and
exactly $m$ only if every factor takes the term $g_y\,x^{(b_f)}$. Blocks are therefore solved by weight ascending
and, within a weight, by number of factors descending. The diagonal block is
$$M - \Phi^{[b]}\,M P = R, \qquad \Phi^{[b]} = F_b\,\Phi^{\otimes m} D_b,\quad \Phi = g^*_y,$$
with $M = [\mu_b(1)\ \cdots\ \mu_b(S)]$. $\Phi$ is the same in every regime.

**Solving a diagonal block.** The solution is $M = \sum_{i\ge0} (\Phi^{[b]})^i R P^i$. Smith doubling: $M_0 = R$,
$M_{j+1} = M_j + (\Phi^{2^j})^{[b]} M_j P^{2^j}$. It uses $(\Phi^{[b]})^{2^j} = (\Phi^{2^j})^{[b]}$, which holds
because $\Phi^{\otimes m}$ maps symmetric vectors to symmetric vectors. Each product is applied as $D_b$, then $m$
cyclic mode products with $\Phi^{2^j}$, then $F_b$, at cost $O(m\,n^{m+1})$ per regime. No decomposition of $P$ is
used. Since $\rho(P) = 1$, it converges when $\rho(\Phi) < 1$. Two direct solvers are available for checking:
Dynare's `gensylv` on the unfolded block with the regimes stacked,
$\mathcal X - (P^\top\otimes\Phi)\,\mathcal X\,(\Phi^\top)^{\otimes(m-1)} = \mathcal R$, and dense LU of
$I - P^\top\otimes\Phi^{[b]}$. For $n = 6$, $m = 8$ and a non-normal
$\Phi$ with $\rho(\Phi) = 0.98$, doubling took 9 steps and 0.1 s on the unfolded block ($6^8$ unknowns), against
2 s for Dynare's Schur-recursion solver `gensylv`, with solutions agreeing to $10^{-13}$.

**Moments of $\hat y$.** $\hat y_t$ is a sum of terms $G(s_t)(z_{t-1}\otimes u_t^{\otimes j})$, with $z$ a Kronecker
product of lagged layers. Terms with equal $(z, j)$ are merged. $\mathbb{E}[\hat y_t \mathbf 1\{s_t = s\}]$ and
$\mathbb{E}[\hat y_t\hat y_t^\top\mathbf 1\{s_t = s\}]$ follow from one- and two-factor versions of the same
contraction.

# 8. Implementation

| Component | File |
|---|---|
| Solver class `KOrderMS` (pairs of regimes, solves stacked over regimes) | `mex/sources/libkorder/kord/korder_ms.{hh,cc}` |
| MEX | `mex/sources/k_order_perturbation_ms/k_order_perturbation_ms.cc` |
| Stock changes | `korder.hh` (`friend class KOrderMS`); `tl/stack_container.hh` (`inner_sigma_zero`) |
| MATLAB | `matlab/ms_korder/`: `ms_korder` (driver), `ms_korder_solve`, `ms_korder_simul`, `ms_korder_moments`, `ms_korder_eval`, `ms_markov_chain` |
| Tests | `tests/ms_korder/` |

# 9. Validation

- **One regime against stock Dynare:**
  - derivatives up to order 4 agree with `k_order_perturbation` to $10^{-16}$;
  - the pruned and unpruned simulations agree with `simult_` to $10^{-14}$;
  - the pruned moments agree with `pruned_state_space_system` at orders 2 and 3 to $10^{-15}$.
- **Lucas tree with switching dividend growth** (exact closed form):
  - policy derivatives up to order 4 in both regimes agree to $1.3\times10^{-14}$ relative;
  - pruned moments up to order 4 agree to $2.3\times10^{-15}$. The oracle is the Taylor polynomial of the closed
    form, with expectations by Gauss–Hermite quadrature.
- **RBC with two switching parameters, order 4:**
  - against an independent dense prototype solver, homogeneous parts of the policy at 20 points: $6\times10^{-11}$
    of the size of the order-4 part, $3\times10^{-12}$ at orders 2 and 3;
  - pruned simulation against layers computed as $\sigma$-derivatives of the prototype policy by automatic
    differentiation: $1.4\times10^{-11}$ relative;
  - moments against a Monte Carlo of $2\times10^6$ periods: 13 quantities within 1.7 standard errors.
- **Moment engine:**
  - against an earlier scalar-monomial implementation of the same recursion: $6\times10^{-15}$ (3 states,
    2 correlated shocks, orders 2–4);
  - the three diagonal-block solvers (doubling, `gensylv`, dense LU) agree to $7\times10^{-15}$.
- **Linear model with three regimes and a forward-looking variable** (`tests/ms_korder/ms_linear_exact.mod`):
  solution and moments against the exact solution and the Markov jump linear recursion, $3\times10^{-15}$, for a
  transition matrix with real eigenvalues and one with complex eigenvalues.
- **Transition matrix with complex eigenvalues** ($1$, $0.4 \pm 0.26i$), RBC, three regimes, order 4:
  - against the dense prototype solver: $6\times10^{-11}$ of the size of the order-4 part;
  - moments: the three diagonal-block solvers agree to $10^{-14}$; 14 quantities within 0.9 standard errors of a
    Monte Carlo of $2\times10^6$ periods.

# References

- Costa, O. L. V., M. D. Fragoso and R. P. Marques (2005), *Discrete-Time Markov Jump Linear Systems*, Springer.
- Foerster, A., J. F. Rubio-Ramírez, D. F. Waggoner and T. Zha (2016), Perturbation methods for Markov-switching
  dynamic stochastic general equilibrium models, *Quantitative Economics* 7(2), 637–669.
- Lombardo, G. and H. Uhlig (2018), A theory of pruning, *International Economic Review* 59(4), 1825–1836.
