function ex = ms_korder_export(ms, M_, file)
% ex = ms_korder_export(ms, M_[, file])
% Self-describing export of the Markov-switching order-K policy (output of ms_korder_solve), for evaluation outside
% Dynare (matlab/ms_korder/python/ms_policy.py). Written with save -v7 (readable by scipy.io.loadmat).
%
% Policy of regime s (deviation from the steady state, declaration order of the endogenous variables):
%   y_t - ys = sum_{1 <= i+j+k <= K} C{s}.y<i>u<j>s<k> * kron(x_{t-1}^{(x)i}, u_t^{(x)j}) * sigma^k,
%   C = g_{y^i u^j sigma^k}(s) / (i! j! k!)  (Taylor coefficients),
% x_{t-1} = deviation of the state variables (state_names order) from their steady state, u_t the innovations
% (exo_names order); kron as in MATLAB (first factor slowest): column index of x_{a1}...x_{ai} u_{b1}...u_{bj} is
% ((((a1-1) n + a2-1) n ...) nu + b1-1) ... + 1 with n states and nu innovations. sigma = 1 is the model.
%
% INPUTS  ms    output of ms_korder_solve;  M_  the model structure;  file  (optional) .mat file to write
% OUTPUT  ex    the exported structure:
%         K, S, P [S x S] (P(s,s') = Pr(s_{t+1} = s' | s_t = s)), pi [S x 1] ergodic probabilities, ahat, switch_vars,
%         endo_names, exo_names, ys [endo_nbr x 1], state_names, state_index (declaration indices of the states),
%         Sigma_e, Raw_moments_e, Skew_e, terms [nterm x 3] = (i, j, k),
%         coef: struct with fields r<s>_y<i>u<j>s<k> [endo_nbr x n^i nu^j], one per regime s and term (i, j, k),
%         convention (text).
K = ms.order; S = numel(ms.g); ny = numel(ms.order_var);
inv_ov(ms.order_var) = 1:ny;
pi_ = null(ms.P' - eye(S));
if size(pi_, 2) == 1, pi_ = pi_ / sum(pi_); else, pi_ = NaN(S, 1); end
terms = zeros(0, 3);
for d = 1:K
    for i = d:-1:0
        for j = (d - i):-1:0
            terms(end + 1, :) = [i, j, d - i - j]; %#ok<AGROW>
        end
    end
end
coef = struct();
for s = 1:S
    for q = 1:size(terms, 1)
        i = terms(q, 1); j = terms(q, 2); k = terms(q, 3);
        G = ms.g{s}.(sprintf('y%du%ds%d', i, j, k));
        coef.(sprintf('r%d_y%du%ds%d', s, i, j, k)) = G(inv_ov, :) / (factorial(i) * factorial(j) * factorial(k));
    end
end
fld = @(f) ifelse_field(M_, f);
ex = struct('K', K, 'S', S, 'P', ms.P, 'pi', pi_, 'ahat', ms.ahat, 'switch_vars', {ms.switch_vars(:)'}, ...
            'endo_names', {M_.endo_names(:)'}, 'exo_names', {M_.exo_names(:)'}, 'ys', ms.ys(:), ...
            'state_names', {M_.endo_names(ms.state_var)'}, 'state_index', ms.state_var(:), ...
            'Sigma_e', M_.Sigma_e, 'Raw_moments_e', fld('Raw_moments_e'), 'Skew_e', fld('Skew_e'), ...
            'terms', terms, 'coef', coef, ...
            'convention', ['y - ys = sum over regime terms coef.r<s>_y<i>u<j>s<k> * kron(x^(x)i, u^(x)j) * sigma^k; ' ...
                           'rows: endo_names order; x: deviation of state_names from ys; kron first factor slowest']);
if nargin >= 3 && ~isempty(file)
    save(file, '-struct', 'ex', '-v7');
end
end

function v = ifelse_field(M_, f)
if isfield(M_, f), v = M_.(f); else, v = []; end
if isempty(v), v = zeros(0, 0); end
end
