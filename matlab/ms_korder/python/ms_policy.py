"""Evaluator of the Markov-switching order-K policy exported by ms_korder_export.m.

    pol = MSPolicy("policy.mat")                       # numpy
    y = pol(s, x, u)                                   # levels, (B, endo_nbr), declaration order
    pol_j = MSPolicy("policy.mat", xp=jax.numpy)       # same, with jax.numpy (jit/grad-able; enable x64 for float64)

Policy of regime s (deviation from the steady state ys):
    y - ys = sum_{1 <= i+j+k <= K} C[s, i, j, k] @ kron(x^{(x)i}, u^{(x)j}) * sigma^k,
x = deviation of the state variables (pol.state_names) at t-1, u = innovations at t (pol.exo_names), s in 1..S.
kron as in MATLAB: first factor slowest.

The array module xp must provide asarray, ones, arange, reshape, sum and broadcasting (numpy, jax.numpy,
tensorflow.experimental.numpy).
"""
import numpy as np
import scipy.io


def _cellstr(a):
    return [str(v) for v in np.atleast_1d(a)]


class MSPolicy:
    def __init__(self, path, xp=np, dtype="float64"):
        d = scipy.io.loadmat(path, simplify_cells=True)
        self.xp = xp
        self.K = int(d["K"])
        self.S = int(d["S"])
        self.P = np.atleast_2d(d["P"]).astype(float)
        self.pi = np.atleast_1d(d["pi"]).astype(float)
        self.endo_names = _cellstr(d["endo_names"])
        self.exo_names = _cellstr(d["exo_names"])
        self.state_names = _cellstr(d["state_names"])
        self.switch_vars = _cellstr(d["switch_vars"])
        self.ys = np.atleast_1d(d["ys"]).astype(float)
        self.terms = np.atleast_2d(d["terms"]).astype(int)
        self.n = len(self.state_names)
        self.nu = len(self.exo_names)
        self.ny = len(self.endo_names)
        self.dtype = dtype
        # coef[(i, j, k)] = array (S, ny, n^i nu^j)
        coef = d["coef"]
        self.coef_np = {}
        for i, j, k in self.terms:
            mats = []
            for s in range(1, self.S + 1):
                c = np.asarray(coef[f"r{s}_y{i}u{j}s{k}"], dtype=float).reshape(self.ny, self.n ** i * self.nu ** j)
                mats.append(c)
            self.coef_np[(int(i), int(j), int(k))] = np.stack(mats)
        self._sigma = None
        self._collapse(1.0)

    def _collapse(self, sigma):
        """Sum the sigma^k terms: A[(i, j)] (S, ny, n^i nu^j) at this sigma, as xp arrays."""
        if self._sigma == sigma:
            return
        acc = {}
        for (i, j, k), c in self.coef_np.items():
            acc[(i, j)] = acc.get((i, j), 0.0) + c * sigma ** k
        self.A = {key: self.xp.asarray(v, dtype=self.dtype) for key, v in acc.items()}
        self.ys_x = self.xp.asarray(self.ys, dtype=self.dtype)
        self._sigma = sigma

    def deviation(self, s, x, u, sigma=1.0):
        """y - ys, (B, ny). s: (B,) integer regimes in 1..S; x: (B, n) state deviations; u: (B, nu) innovations."""
        xp = self.xp
        self._collapse(float(sigma))
        x = xp.asarray(x, dtype=self.dtype)
        u = xp.asarray(u, dtype=self.dtype)
        B = x.shape[0]
        # features kron(x^{(x)i}, u^{(x)j}) for all (i, j) needed, built by appending fastest factors
        xpow = [xp.ones((B, 1), dtype=self.dtype)]
        for i in range(1, self.K + 1):
            xpow.append(xp.reshape(xpow[-1][:, :, None] * x[:, None, :], (B, -1)))
        s = xp.reshape(xp.asarray(s), (-1, 1))
        onehot = xp.asarray(s == xp.reshape(xp.arange(1, self.S + 1), (1, -1)), dtype=self.dtype)
        y = xp.zeros((B, self.ny), dtype=self.dtype)
        for (i, j), A in self.A.items():
            f = xpow[i]
            for _ in range(j):
                f = xp.reshape(f[:, :, None] * u[:, None, :], (B, -1))
            for r in range(self.S):
                y = y + onehot[:, r:r + 1] * (f @ xp.transpose(A[r]))
        return y

    def __call__(self, s, x, u, sigma=1.0):
        """Levels y, (B, ny), declaration order."""
        return self.ys_x[None, :] + self.deviation(s, x, u, sigma)
