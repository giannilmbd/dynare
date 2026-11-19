function u = ms_korder_draw_shocks(shocks, T, draw)
% u = ms_korder_draw_shocks(shocks, T[, draw])
% Draw of T innovation vectors [nu x T] from the current random stream.
%
% INPUTS  shocks  M_ (Sigma_e, and optionally Raw_moments_e, Skew_e) or the covariance of Gaussian innovations
%         T       number of periods
%         draw    (optional) function handle, draw(T) returns [nu x T]; required when a shock is defined only by
%                 its raw moments (Raw_moments_e), which do not determine a distribution to draw from
% Without draw: the Gaussian block is drawn from N(0, Sigma_e(G, G)) (eigenvalue factor, valid for a singular
% covariance); a skew normal shock (Skew_e) as omega (delta (|U0| - sqrt(2/pi)) + sqrt(1 - delta^2) U1), with delta and
% omega as in skew_normal_raw_moments.
if isnumeric(shocks), shocks = struct('Sigma_e', shocks); end
Sig = shocks.Sigma_e; nu = size(Sig, 1);
if nargin >= 3 && ~isempty(draw)
    u = draw(T);
    if ~isequal(size(u), [nu, T]), error('ms_korder_draw_shocks: draw(T) must return a %d x T matrix', nu), end
    return
end
if isfield(shocks, 'Raw_moments_e') && ~isempty(shocks.Raw_moments_e) && any(~all(isnan(shocks.Raw_moments_e), 2))
    error(['ms_korder_draw_shocks: shocks defined by Raw_moments_e need a draw function (the raw moments do not ' ...
           'define a distribution)'])
end
skew = zeros(nu, 1);
if isfield(shocks, 'Skew_e') && ~isempty(shocks.Skew_e)
    for r = 1:size(shocks.Skew_e, 1)
        if shocks.Skew_e(r, 4) ~= 0, skew(shocks.Skew_e(r, 1)) = shocks.Skew_e(r, 4); end
    end
end
sk = find(skew ~= 0)'; G = setdiff(1:nu, sk);
u = zeros(nu, T);
if ~isempty(G)
    [V, D] = eig((Sig(G, G) + Sig(G, G)') / 2);
    u(G, :) = V * diag(sqrt(max(diag(D), 0))) * randn(numel(G), T);
end
for i = sk
    if any(Sig(i, [1:i-1, i+1:nu]) ~= 0), error('ms_korder_draw_shocks: skewed shock %d must be uncorrelated', i), end
    s23 = abs(skew(i))^(2/3);
    delta = sign(skew(i)) * sqrt(pi/2 * s23 / (s23 + ((4 - pi)/2)^(2/3)));
    omega = sqrt(Sig(i, i) / (1 - 2/pi * delta^2));
    u(i, :) = omega * (delta * (abs(randn(1, T)) - sqrt(2/pi)) + sqrt(1 - delta^2) * randn(1, T));
end
end
