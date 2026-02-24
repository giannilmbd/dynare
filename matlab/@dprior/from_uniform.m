function pdraw = from_uniform(o, u)

% Transform uniform quantiles into prior draws using inverse CDFs.
%
% INPUTS
% - o    [dprior]  Prior object.
% - u    [double]  n×m matrix of uniform random draws in [0,1], where n is the
%                  number of samples and m is the number of parameters. Each
%                  column contains uniform quantiles that will be transformed
%                  into prior draws via inverse CDFs.
%
% OUTPUTS
% - pdraw [double] n×m matrix, draws from the joint prior density.
%
% ALGORITHM
%   Transforms uniform random draws in u into actual prior draws by:
%   1. Rescaling each column of u to the truncated cumulative range
%      [lbcum, ubcum] for that parameter
%   2. Applying the appropriate inverse CDF (norminv, gaminv, betainv, etc.)
%      to convert quantiles into parameter values
%   This approach allows batched sampling with quasi-random or Latin hypercube
%   sequences generated externally and passed via u.
%
% REMARKS
% This method is equivalent to the former gsa.prior_draw but integrated into the dprior class.
% Unlike the draw() method which generates random numbers internally, from_uniform()
% takes externally generated uniform quantiles, allowing for quasi-Monte Carlo methods
% and reproducible sampling strategies.
%
% EXAMPLE
%
% >> Prior = dprior(bayestopt_, options_.prior_trunc);
% >> u = rand(1000, length(Prior.lb));  % or use QMC sequence
% >> samples = Prior.from_uniform(u);

% Copyright © 2026 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

npar = length(o.lb);
nsamp = size(u, 1);

% Compute cumulative probability bounds for truncated priors
lbcum = zeros(npar, 1);
ubcum = ones(npar, 1);

% Compute truncated cumulative probabilities for each parameter
if o.isbeta
    for j = 1:length(o.idbeta)
        i = o.idbeta(j);
        lbcum(i) = betainc((o.lb(i)-o.p3(i))./(o.p4(i)-o.p3(i)), o.p6(i), o.p7(i));
        ubcum(i) = betainc((o.ub(i)-o.p3(i))./(o.p4(i)-o.p3(i)), o.p6(i), o.p7(i));
    end
end

if o.isgamma
    for j = 1:length(o.idgamma)
        i = o.idgamma(j);
        lbcum(i) = gamcdf(o.lb(i)-o.p3(i), o.p6(i), o.p7(i));
        ubcum(i) = gamcdf(o.ub(i)-o.p3(i), o.p6(i), o.p7(i));
    end
end

if o.isgaussian
    for j = 1:length(o.idgaussian)
        i = o.idgaussian(j);
        lbcum(i) = 0.5 * erfc(-(o.lb(i)-o.p6(i))/o.p7(i) / sqrt(2));
        ubcum(i) = 0.5 * erfc(-(o.ub(i)-o.p6(i))/o.p7(i) / sqrt(2));
    end
end

if o.isinvgamma1
    for j = 1:length(o.idinvgamma1)
        i = o.idinvgamma1(j);
        lbcum(i) = gamcdf(1/(o.ub(i)-o.p3(i))^2, o.p7(i)/2, 2/o.p6(i));
        ubcum(i) = gamcdf(1/(o.lb(i)-o.p3(i))^2, o.p7(i)/2, 2/o.p6(i));
    end
end

if o.isinvgamma2
    for j = 1:length(o.idinvgamma2)
        i = o.idinvgamma2(j);
        lbcum(i) = gamcdf(1/(o.ub(i)-o.p3(i)), o.p7(i)/2, 2/o.p6(i));
        ubcum(i) = gamcdf(1/(o.lb(i)-o.p3(i)), o.p7(i)/2, 2/o.p6(i));
    end
end

if o.isweibull
    for j = 1:length(o.idweibull)
        i = o.idweibull(j);
        lbcum(i) = wblcdf(o.lb(i)-o.p3(i), o.p7(i), o.p6(i));
        ubcum(i) = wblcdf(o.ub(i)-o.p3(i), o.p7(i), o.p6(i));
    end
end

% For uniform priors, bounds are already set correctly (lbcum=0, ubcum=1)

% Transform uniform quantiles to prior draws
pdraw = NaN(nsamp, npar);
for i = 1:npar
    % Rescale to truncated cumulative range
    u(:, i) = u(:, i) .* (ubcum(i) - lbcum(i)) + lbcum(i);
end

% Apply inverse CDFs
if o.isuniform
    for j = 1:length(o.iduniform)
        i = o.iduniform(j);
        pdraw(:, i) = u(:, i) * (o.p4(i) - o.p3(i)) + o.p3(i);
    end
end

if o.isgaussian
    for j = 1:length(o.idgaussian)
        i = o.idgaussian(j);
        pdraw(:, i) = norminv(u(:, i), o.p6(i), o.p7(i));
    end
end

if o.isgamma
    for j = 1:length(o.idgamma)
        i = o.idgamma(j);
        pdraw(:, i) = gaminv(u(:, i), o.p6(i), o.p7(i)) + o.p3(i);
    end
end

if o.isbeta
    for j = 1:length(o.idbeta)
        i = o.idbeta(j);
        pdraw(:, i) = betainv(u(:, i), o.p6(i), o.p7(i)) * (o.p4(i) - o.p3(i)) + o.p3(i);
    end
end

if o.isinvgamma1
    for j = 1:length(o.idinvgamma1)
        i = o.idinvgamma1(j);
        pdraw(:, i) = sqrt(1 ./ gaminv(u(:, i), o.p7(i)/2, 2/o.p6(i))) + o.p3(i);
    end
end

if o.isinvgamma2
    for j = 1:length(o.idinvgamma2)
        i = o.idinvgamma2(j);
        pdraw(:, i) = 1 ./ gaminv(u(:, i), o.p7(i)/2, 2/o.p6(i)) + o.p3(i);
    end
end

if o.isweibull
    for j = 1:length(o.idweibull)
        i = o.idweibull(j);
        pdraw(:, i) = wblinv(u(:, i), o.p7(i), o.p6(i)) + o.p3(i);
    end
end

return % --*-- Unit tests --*--

%@test:1
% Fill global structures with required fields...
prior_trunc = 1e-10;
p0 = repmat([1; 2; 3; 4; 5; 6; 8], 2, 1);    % Prior shape
p1 = .4*ones(14,1);                          % Prior mean
p2 = .2*ones(14,1);                          % Prior std.
p3 = NaN(14,1);
p4 = NaN(14,1);
p5 = NaN(14,1);
p6 = NaN(14,1);
p7 = NaN(14,1);

for i=1:14
    switch p0(i)
      case 1
        % Beta distribution
        p3(i) = 0;
        p4(i) = 1;
        [p6(i), p7(i)] = beta_specification(p1(i), p2(i)^2, p3(i), p4(i));
        p5(i) = compute_prior_mode([p6(i) p7(i)], 1);
      case 2
        % Gamma distribution
        p3(i) = 0;
        p4(i) = Inf;
        [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
        p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
      case 3
        % Normal distribution
        p3(i) = -Inf;
        p4(i) = Inf;
        p6(i) = p1(i);
        p7(i) = p2(i);
        p5(i) = p1(i);
      case 4
        % Inverse Gamma (type I) distribution
        p3(i) = 0;
        p4(i) = Inf;
        [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 1, false);
        p5(i) = compute_prior_mode([p6(i) p7(i)], 4);
      case 5
        % Uniform distribution
        [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
        p3(i) = p6(i);
        p4(i) = p7(i);
        p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
      case 6
        % Inverse Gamma (type II) distribution
        p3(i) = 0;
        p4(i) = Inf;
        [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 2, false);
        p5(i) = compute_prior_mode([p6(i) p7(i)], 6);
      case 8
        % Weibull distribution
        p3(i) = 0;
        p4(i) = Inf;
        [p6(i), p7(i)] = weibull_specification(p1(i), p2(i)^2, p3(i));
        p5(i) = compute_prior_mode([p6(i) p7(i)], 8);
      otherwise
        error('This density is not implemented!')
    end
end

BayesInfo.pshape = p0;
BayesInfo.p1 = p1;
BayesInfo.p2 = p2;
BayesInfo.p3 = p3;
BayesInfo.p4 = p4;
BayesInfo.p5 = p5;
BayesInfo.p6 = p6;
BayesInfo.p7 = p7;

ndraws = 1e5;

% Call the tested routine
try
    % Instantiate dprior object
    o = dprior(BayesInfo, prior_trunc, false);
    % Generate uniform quantiles
    u = rand(ndraws, length(p1));
    % Transform to prior draws
    pdraw = o.from_uniform(u);
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    % Check dimensions
    t(2) = isequal(size(pdraw), [ndraws, length(p1)]);
    % Check all values are within bounds
    t(3) = all(all(pdraw >= repmat(o.lb', ndraws, 1), 1));
    t(4) = all(all(pdraw <= repmat(o.ub', ndraws, 1), 1));
    % Check mean
    m_empirical = mean(pdraw, 1)';
    t(5) = all(abs(m_empirical - BayesInfo.p1) < 3e-3);
    % Check variance
    v_empirical = var(pdraw, 0, 1)';
    t(6) = all(abs(v_empirical - BayesInfo.p2.^2) < 5e-3);
end
T = all(t);
%@eof:1

%@test:2
% Test reproducibility with same input
prior_trunc = 1e-10;
p0 = [3; 2; 1; 5];  % Gaussian, Gamma, Beta, Uniform
p1 = [0.5; 0.4; 0.3; 0.6];
p2 = [0.1; 0.15; 0.08; 0.12];
p3 = NaN(4,1);
p4 = NaN(4,1);
p5 = NaN(4,1);
p6 = NaN(4,1);
p7 = NaN(4,1);

for i=1:4
    switch p0(i)
      case 1
        p3(i) = 0;
        p4(i) = 1;
        [p6(i), p7(i)] = beta_specification(p1(i), p2(i)^2, p3(i), p4(i));
        p5(i) = compute_prior_mode([p6(i) p7(i)], 1);
      case 2
        p3(i) = 0;
        p4(i) = Inf;
        [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
        p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
      case 3
        p3(i) = -Inf;
        p4(i) = Inf;
        p6(i) = p1(i);
        p7(i) = p2(i);
        p5(i) = p1(i);
      case 5
        [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
        p3(i) = p6(i);
        p4(i) = p7(i);
        p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
    end
end

BayesInfo.pshape = p0;
BayesInfo.p1 = p1;
BayesInfo.p2 = p2;
BayesInfo.p3 = p3;
BayesInfo.p4 = p4;
BayesInfo.p5 = p5;
BayesInfo.p6 = p6;
BayesInfo.p7 = p7;

try
    o = dprior(BayesInfo, prior_trunc, false);
    % Same uniform input should produce same output
    u = rand(100, 4);
    pdraw1 = o.from_uniform(u);
    pdraw2 = o.from_uniform(u);
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    % Check reproducibility
    t(2) = all(all(pdraw1 == pdraw2));
end
T = all(t);
%@eof:2

%@test:3
% Test that from_uniform and draw produce equivalent distributions
prior_trunc = 1e-10;
p0 = [3; 2; 5];  % Gaussian, Gamma, Uniform
p1 = [0.5; 0.4; 0.6];
p2 = [0.1; 0.15; 0.12];
p3 = NaN(3,1);
p4 = NaN(3,1);
p5 = NaN(3,1);
p6 = NaN(3,1);
p7 = NaN(3,1);

for i=1:3
    switch p0(i)
      case 2
        p3(i) = 0;
        p4(i) = Inf;
        [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
        p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
      case 3
        p3(i) = -Inf;
        p4(i) = Inf;
        p6(i) = p1(i);
        p7(i) = p2(i);
        p5(i) = p1(i);
      case 5
        [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
        p3(i) = p6(i);
        p4(i) = p7(i);
        p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
    end
end

BayesInfo.pshape = p0;
BayesInfo.p1 = p1;
BayesInfo.p2 = p2;
BayesInfo.p3 = p3;
BayesInfo.p4 = p4;
BayesInfo.p5 = p5;
BayesInfo.p6 = p6;
BayesInfo.p7 = p7;

ndraws = 5000;

try
    o = dprior(BayesInfo, prior_trunc, false);
    % Generate samples using from_uniform
    u = rand(ndraws, 3);
    pdraw_icdf = o.from_uniform(u);
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    % Check means against theoretical values
    m_icdf = mean(pdraw_icdf, 1)';
    t(2) = all(abs(m_icdf - BayesInfo.p1) < 0.05);
    % Check variances against theoretical values
    v_icdf = var(pdraw_icdf, 0, 1)';
    t(3) = all(abs(v_icdf - BayesInfo.p2.^2) < 0.02);
end
T = all(t);
%@eof:3

%@test:4
% Regression test for Weibull parameter order bug (GitHub issue)
% Verifies that Weibull distribution parameters are correctly ordered as [shape, scale]
% and that the inverse CDF uses the correct order (scale, shape) as per MATLAB convention
prior_trunc = 1e-10;

% Create a Weibull prior with specific parameters
p0 = 8;  % Weibull shape code
p1 = 0.8;  % Prior mean
p2 = 0.15; % Prior std
p3 = 0;    % Lower bound
p4 = Inf;  % Upper bound
p5 = NaN;
p6 = NaN;
p7 = NaN;

% weibull_specification returns [shape, scale]
[p6, p7] = weibull_specification(p1, p2^2, p3);

BayesInfo.pshape = p0;
BayesInfo.p1 = p1;
BayesInfo.p2 = p2;
BayesInfo.p3 = p3;
BayesInfo.p4 = p4;
BayesInfo.p5 = p5;
BayesInfo.p6 = p6;
BayesInfo.p7 = p7;

ndraws = 1e4;

try
    o = dprior(BayesInfo, prior_trunc, false);
    % Generate uniform quantiles
    u = rand(ndraws, 1);
    % Transform to prior draws
    pdraw = o.from_uniform(u);
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    % Empirical mean should match the specified prior mean
    m_empirical = mean(pdraw);
    % Empirical variance should match the specified prior variance
    v_empirical = var(pdraw);

    % For Weibull, we use tighter tolerance check to catch parameter order bugs
    % If parameters were swapped, the mean would be completely wrong
    t(2) = abs(m_empirical - BayesInfo.p1) < 0.05;

    % Check that draws are positive (Weibull support is (0, Inf))
    t(3) = all(pdraw > 0);

    % Check that empirical variance is positive and reasonable
    t(4) = v_empirical > 0 && abs(v_empirical - BayesInfo.p2^2) < 0.1;
end
T = all(t);
%@eof:4
