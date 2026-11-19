function s = ms_markov_chain(P, T, s0, seed)
% s = ms_markov_chain(P, T, s0[, seed]): regime path s(1..T), s(1) drawn from P(s0, :), P(a, b) = Pr(b | a).
if nargin >= 4, rng(seed, 'twister'); end
cP = cumsum(P, 2); s = zeros(1, T); prev = s0;
for t = 1:T
    s(t) = find(rand <= cP(prev, :), 1, 'first');
    if isempty(s(t)), s(t) = size(P, 1); end
    prev = s(t);
end
end
