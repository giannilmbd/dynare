function y = ms_korder_eval(ms, s, ystate, u, sigma, orders)
% y = ms_korder_eval(ms, s, ystate, u, sigma[, orders])
% Policy of regime s in deviation from the steady state, rows in ms.order_var order:
%   sum over i + j + k in orders of g_{y^i u^j sigma^k} kron(ystate^{(x)i}, u^{(x)j}) sigma^k / (i! j! k!).
% ystate: deviation of the state variables (ms.state_var order); u: innovations; orders default 1:ms.order.
if nargin < 6, orders = 1:ms.order; end
ystate = ystate(:); u = u(:);
y = zeros(size(ms.g{s}.y1u0s0, 1), 1);
for d = orders
    for i = d:-1:0
        for j = (d - i):-1:0
            k = d - i - j;
            v = 1;
            for a = 1:i, v = kron(v, ystate); end
            for a = 1:j, v = kron(v, u); end
            T = ms.g{s}.(sprintf('y%du%ds%d', i, j, k));
            y = y + T * v * sigma^k / (factorial(i) * factorial(j) * factorial(k));
        end
    end
end
end
