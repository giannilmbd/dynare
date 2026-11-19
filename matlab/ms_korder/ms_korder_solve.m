function ms = ms_korder_solve(M_, options_, dr, P, ahat)
% ms = ms_korder_solve(M_, options_, dr, P, ahat)
% Order-k perturbation with an exogenous Markov chain (MEX k_order_perturbation_ms, libkorder/kord/korder_ms.hh).
%
% Switching parameters are endogenous variables, each with a dummy equation tagged [ms='<its name>']:
%     parameters ms_zero; ms_zero = 0;  ...  [ms='aD'] aD = ms_zero;
% (a literal 'aD = 0' is substituted into the other equations by the preprocessor, which removes aD(+1).)
% The regime value is aD_t = sigma * ahat(r, s_t), r the row of aD in ms.switch_vars (order of the tags in the model).
% Call after stoch_simul(order=K, k_order_solver, ...) (sets dr.ys, dr.order_var, options_.order = K >= 2).
%
% INPUTS   P     [S x S]  transition matrix, P(s, s') = Pr(s_{t+1} = s' | s_t = s); real eigenvalues required
%          ahat  [n_switch x S]  regime values of the switching variables (deviation from their steady state 0)
% OUTPUT   ms.g{s}.y<i>u<j>s<k>: unfolded derivative of the policy of regime s, rows in dr.order_var order, columns
%          kron(ystate^{(x)i}, u^{(x)j}), ystate = dr.order_var(M_.nstatic + (1:M_.nspred)); i + j + k <= order.
%          ms.switch_vars, ms.state_var, ms.order_var, ms.P, ms.ahat, ms.order.

tags = M_.equations_tags;
is_ms = strcmp(tags(:, 2), 'ms');
eqs = cell2mat(tags(is_ms, 1));
names = tags(is_ms, 3);
if isempty(eqs)
    error('ms_korder_solve: no equation tagged [ms=''<var>'']')
end
S = size(P, 1);
if size(P, 2) ~= S || any(abs(sum(P, 2) - 1) > 1e-12) || any(P(:) < 0)
    error('ms_korder_solve: P must be a square stochastic matrix (rows sum to 1)')
end
if ~isequal(size(ahat), [numel(eqs), S])
    error('ms_korder_solve: ahat must be %d x %d (switching variables x regimes)', numel(eqs), S)
end
for r = 1:numel(eqs)
    iv = find(strcmp(M_.endo_names, names{r}));
    if isempty(iv), error('ms_korder_solve: tagged variable %s is not an endogenous variable', names{r}), end
    if abs(dr.ys(iv)) > 1e-12, error('ms_korder_solve: steady state of %s must be 0', names{r}), end
    n = M_.endo_nbr; cp = M_.dynamic_g1_sparse_colptr; rv = M_.dynamic_g1_sparse_rowval;
    rows = [];
    for col = [iv, n + iv, 2*n + iv]
        rows = [rows; reshape(rv(cp(col):cp(col + 1) - 1), [], 1)];
    end
    if ~any(rows ~= eqs(r))
        error(['ms_korder_solve: %s enters no equation other than its dummy equation (write the dummy as ' ...
               '%s = <parameter equal to 0>; a literal 0 is substituted out by the preprocessor)'], names{r}, names{r})
    end
end
forcing = zeros(M_.endo_nbr, S);
forcing(eqs, :) = -ahat;          % residual of the dummy equation: a_t - sigma * ahat(s_t)
[W, L] = eig(P');
lam = diag(L);
if any(abs(imag(lam)) > 1e-12) || any(abs(imag(W(:))) > 1e-12)
    error('ms_korder_solve: P has complex eigenvalues (not implemented: real Schur form needed)')
end
lam = real(lam); W = real(W);
if rcond(W) < 1e-10
    error('ms_korder_solve: eigenvectors of P are (nearly) dependent, rcond %g', rcond(W))
end
msin = struct('P', P, 'lam', lam, 'W', W, 'Winv', inv(W), 'forcing', forcing);
ms.g = k_order_perturbation_ms(dr, M_, options_, msin);
ms.switch_vars = names;
ms.order_var = dr.order_var;
ms.state_var = dr.order_var(M_.nstatic + (1:M_.nspred));
ms.P = P; ms.ahat = ahat; ms.order = options_.order;
ms.ys = dr.ys;
end
