function oo_ = ms_korder(M_, options_, oo_, regimes, opts)
% oo_ = ms_korder(M_, options_, oo_, regimes[, opts])
% Order-k perturbation of a model with parameters that switch with an exogenous Markov chain s_t in {1,...,S}:
% solution, and optionally exact moments and a simulation of the pruned solution. Call after
% stoch_simul(order = K, k_order_solver, ...), K >= 2. See matlab/ms_korder/README.md.
%
% Model convention: each switching parameter is an endogenous variable a with steady state 0 and a dummy equation
% tagged with its name, reading a = <a parameter set to 0>, e.g.
%     var ... aA;  parameters ms_zero ...;  ms_zero = 0;
%     model;  ... (Abar + aA) ... ;  [ms='aA'] aA = ms_zero;  end;
% In regime s the variable equals sigma * ahat_a(s) (sigma = 1 at the solution), so Abar + aA is Abar + ahat_aA(s).
%
% INPUTS  regimes.P       [S x S] transition matrix, P(s, s') = Pr(s_{t+1} = s' | s_t = s)
%         regimes.values  struct: one field per tagged variable, value [1 x S] = ahat(s); untagged fields are an
%                         error, tagged variables without a field are 0 in every regime
%         regimes.names   (optional) cellstr of S regime labels
%         opts.moments    moments of the pruned solution (default true)
%         opts.order      highest moment, 1 to 4 (default 2); orders 3, 4 give skewness and kurtosis of the
%                         non-switching variables (opts.second = true/false is read as order 2/1)
%         opts.moment_method  for orders 3, 4: 'auto' (default: exact when K <= 2, or K = 3 and order 3; simulation
%                         otherwise), 'exact', 'simulation'. Means and variances are always exact.
%         opts.sim_periods  periods of the simulated moments (default 1e6, a multiple of 100)
%         opts.solver     diagonal-block solver of the exact moments: 'doubling' (default), 'gensylv', 'lu'
%         opts.periods    length of a pruned simulation (default 0 = none)
%         opts.s0, opts.seed  initial regime (default 1) and seed (default 0) of the simulations
% OUTPUT  oo_.ms.solution  output of ms_korder_solve (policy derivatives by regime)
%         oo_.ms.moments   output of ms_korder_moments (if opts.moments)
%         oo_.ms.sim_moments  output of ms_korder_sim_moments (if orders 3, 4 are simulated)
%         oo_.ms.simulation.y [endo_nbr x T], .s [1 x T], .u [exo_nbr x T] (if opts.periods > 0)
if nargin < 5, opts = struct(); end
if isfield(opts, 'second') && ~isfield(opts, 'order'), opts.order = 1 + double(logical(opts.second)); end
def = struct('moments', true, 'order', 2, 'moment_method', 'auto', 'sim_periods', 1e6, 'solver', 'doubling', ...
             'periods', 0, 's0', 1, 'seed', 0);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts, f{i}), opts.(f{i}) = def.(f{i}); end
end
if options_.order < 2 || ~options_.k_order_solver
    error('ms_korder: run stoch_simul(order = K, k_order_solver) with K >= 2 first')
end
tags = M_.equations_tags;
names = tags(strcmp(tags(:, 2), 'ms'), 3);
P = regimes.P; S = size(P, 1);
ahat = zeros(numel(names), S);
if isfield(regimes, 'values')
    vf = fieldnames(regimes.values);
    for i = 1:numel(vf)
        r = find(strcmp(names, vf{i}));
        if isempty(r), error('ms_korder: regimes.values.%s: no equation tagged [ms=''%s'']', vf{i}, vf{i}), end
        v = regimes.values.(vf{i});
        if numel(v) ~= S, error('ms_korder: regimes.values.%s must have %d entries (one per regime)', vf{i}, S), end
        ahat(r, :) = v(:)';
    end
end
ms = ms_korder_solve(M_, options_, oo_.dr, P, ahat);
if isfield(regimes, 'names'), ms.regime_names = regimes.names; else, ms.regime_names = compose('regime %d', 1:S); end
oo_.ms = struct('solution', ms);
if opts.moments
    keep = find(~ismember(M_.endo_names, ms.switch_vars))';
    K = ms.order;
    switch opts.moment_method
        case 'auto', exact = opts.order <= 2 || K <= 2 || (K == 3 && opts.order == 3);
        case 'exact', exact = true;
        case 'simulation', exact = opts.order <= 2;
        otherwise, error('ms_korder: opts.moment_method must be ''auto'', ''exact'' or ''simulation''')
    end
    if exact
        mom = ms_korder_moments(ms, M_.Sigma_e, opts.order, opts.solver, keep);
        hm = mom;
    else
        mom = ms_korder_moments(ms, M_.Sigma_e, 2, opts.solver, keep);
        hm = ms_korder_sim_moments(ms, M_.Sigma_e, keep, opts.sim_periods, struct('s0', opts.s0, 'seed', opts.seed));
        oo_.ms.sim_moments = hm;
    end
    oo_.ms.moments = mom;
    if ~options_.noprint, print_moments(M_, ms, mom, hm, opts.order, exact); end
end
if opts.periods > 0
    T = opts.periods;
    s = ms_markov_chain(P, T, opts.s0, opts.seed);
    rng(opts.seed + 1, 'twister');
    [Ve, De] = eig((M_.Sigma_e + M_.Sigma_e') / 2);           % factor valid for a singular covariance
    u = Ve * diag(sqrt(max(diag(De), 0))) * randn(M_.exo_nbr, T);
    oo_.ms.simulation = struct('y', ms_korder_simul(ms, s, u, true), 's', s, 'u', u);
end
end

function print_moments(M_, ms, mom, hm, order, exact)
skip = ismember(M_.endo_names, ms.switch_vars);
second = order >= 2;
fprintf('\nMARKOV-SWITCHING PERTURBATION, order %d, %d regimes (pruned)\n', ms.order, numel(mom.pi));
fprintf('Ergodic regime probabilities: %s\n', sprintf('%.4f ', mom.pi));
if order >= 3 && ~exact
    fprintf(['Means, standard deviations: exact. Skewness, kurtosis: simulation of %d periods, ' ...
             'batch-means standard errors in brackets.\n'], hm.T);
else
    fprintf('Exact moments.\n');
end
w = 12 + 12 * (order >= 3 && ~exact);
hdr = sprintf('%-12s %12s', 'VARIABLE', 'MEAN');
if second, hdr = [hdr, sprintf(' %12s', 'STD. DEV.')]; end
if order >= 3, hdr = [hdr, sprintf(' %*s', w, 'SKEWNESS')]; end
if order >= 4, hdr = [hdr, sprintf(' %*s', w, 'KURTOSIS')]; end
for s = 1:numel(mom.pi), hdr = [hdr, sprintf(' %12s', sprintf('MEAN|s=%d', s))]; end %#ok<AGROW>
fprintf('%s\n', hdr);
for v = find(~skip(:))'
    line = sprintf('%-12s %12.6g', M_.endo_names{v}, mom.mean(v));
    if second, line = [line, sprintf(' %12.6g', sqrt(max(mom.var(v, v), 0)))]; end %#ok<AGROW>
    for f = {'skewness', 'kurtosis'}
        if (strcmp(f{1}, 'skewness') && order >= 3) || (strcmp(f{1}, 'kurtosis') && order >= 4)
            x = hm.(f{1})(hm.vars == v);
            if exact, line = [line, sprintf(' %12.6g', x)]; %#ok<AGROW>
            else, line = [line, sprintf(' %12.6g %11s', x, sprintf('(%.2g)', hm.se.(f{1})(hm.vars == v)))]; end %#ok<AGROW>
        end
    end
    fprintf('%s%s\n', line, sprintf(' %12.6g', mom.mean_regime(v, :)));
end
end
