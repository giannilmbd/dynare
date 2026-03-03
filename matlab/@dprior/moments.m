function o = moments(o, name)

% Compute the prior moments.
%
% INPUTS
% - o       [dprior]
%
% OUTPUTS
% - o       [dprior]

% Copyright © 2023 Dynare Team
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

switch name
  case 'mean'
    m = o.p1;
  case 'median'
    m = o.p11;
  case 'std'
    m = o.p2;
  case 'mode'
    m = o.p5;
  otherwise
    error('%s is not an implemented moemnt.', name)
end
id = isnan(m);
if any(id)
    % For some parameters the prior mean is not defined.
    % We compute the first order moment from the
    % hyperparameters, if the hyperparameters are not
    % available an error is thrown.
    if o.isuniform
        jd = intersect(o.iduniform, find(id));
        if ~isempty(jd)
            if any(isnan(o.p3(jd))) || any(isnan(o.p4(jd)))
                error('dprior::mean: Some hyperparameters are missing (uniform distribution).')
            end
            switch name
              case 'mean'
                m(jd) = o.p3(jd) + .5*(o.p4(jd)-o.p3(jd));
              case 'median'
                m(jd) = o.p3(jd) + .5*(o.p4(jd)-o.p3(jd));
              case 'std'
                m(jd) = (o.p4(jd)-o.p3(jd))/sqrt(12);
              case 'mode' % Actually we have a continuum of modes with the uniform distribution.
                m(jd) = o.p3(jd) + .5*(o.p4(jd)-o.p3(jd));
            end
        end
    end
    if o.isgaussian
        jd = intersect(o.idgaussian, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd)))
                error('dprior::mean: Some hyperparameters are missing (gaussian distribution).')
            end
            switch name
              case 'mean'
                m(jd) = o.p6(jd);
              case 'median'
                m(jd) = o.p6(jd);
              case 'std'
                m(jd) = o.p7(jd);
              case 'mode' % Actually we have a continuum of modes with the uniform distribution.
                m(jd) = o.p6(jd);
            end
        end
    end
    if o.isgamma
        jd = intersect(o.idgamma, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd)))
                error('dprior::mean: Some hyperparameters are missing (gamma distribution).')
            end
            % α → o.p6, β → o.p7
            switch name
              case 'mean'
                m(jd) = o.p3(jd) + o.p6(jd).*o.p7(jd);
              case 'median'
                m(jd) = o.p3(jd) + gaminv(.5, o.p6(jd), o.p7(jd));
              case 'std'
                m(jd) = sqrt(o.p6(jd)).*o.p7(jd);
              case 'mode'
                m(jd) = o.p3(jd);
                hd = o.p6(jd)>1;
                m(jd(hd)) = o.p3(jd(hd)) + (o.p6(jd(hd))-1).*o.p7(jd(hd));
            end
        end
    end
    if o.isbeta
        jd = intersect(o.idbeta, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd))) || any(isnan(o.p4(jd)))
                error('dprior::mean: Some hyperparameters are missing (beta distribution).')
            end
            % α → o.p6, β → o.p7
            switch name
              case 'mean'
                m(jd) = o.p3(jd) + (o.p6(jd)./(o.p6(jd)+o.p7(jd))).*(o.p4(jd)-o.p3(jd));
              case 'median'
                m(jd) = o.p3(jd) + betainv(.5, o.p6(jd), o.p7(jd)).*(o.p4(jd)-o.p3(jd));
              case 'std'
                m(jd) = (o.p4(jd)-o.p3(jd)).*sqrt(o.p6(jd).*o.p7(jd)./((o.p6(jd)+o.p7(jd)).^2.*(o.p6(jd)+o.p7(jd)+1)));
              case 'mode'
                h0 = true(length(jd), 1);
                h1 = o.p6(jd)<=1 & o.p7(jd)>1; h0 = h0 & ~h1;
                h2 = o.p7(jd)<=1 & o.p6(jd)>1; h0 = h0 & ~h2;
                h3 = o.p6(jd)<1 & o.p7(jd)<1; h0 = h0 & ~h3;
                h4 = ismembertol(o.p6(jd), 1) & ismembertol(o.p7(jd),1); h0 = h0 & ~h4;
                m(jd(h1)) = o.p3(jd(h1));                                  % Standard β has a mode at 0
                m(jd(h2)) = o.p4(jd(h2));                                  % Standard β has a mode at 1
                m(jd(h3)) = o.p3(jd(h3));                                  % Standard β is bimodal, we pick the lowest mode (0)
                m(jd(h4)) = o.p3(jd(h4)) + .5*(o.p4(jd(h4))-o.p3(jd(h4))); % Standard β is the uniform distribution (continuum of modes), we pick the mean as the mode
                m(jd(h0)) = o.p3(jd(h0))+(o.p4(jd(h0))-o.p3(jd(h0))).*((o.p6(jd(h0))-1)./(o.p6(jd(h0))+o.p7(jd(h0))-2)); % β distribution is concave and has a unique interior mode.
            end
        end
    end
    if o.isinvgamma1
        jd = intersect(o.idinvgamma1, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd)))
                error('dprior::mean: Some hyperparameters are missing (inverse gamma type 1 distribution).')
            end
            % s → o.p6, ν → o.p7
            switch name
              case 'mean'
                m(jd) = o.p3(jd) + sqrt(.5*o.p6(jd)) .*(gamma(.5*(o.p7(jd)-1))./gamma(.5*o.p7(jd)));
              case 'median'
                m(jd) = o.p3(jd) + 1.0/sqrt(gaminv(.5, o.p7(jd)/2.0, 2.0/o.p6(jd)));
              case 'std'
                m(jd) = sqrt( o.p6(jd)./(o.p7(jd)-2)-(.5*o.p6(jd)).*(gamma(.5*(o.p7(jd)-1))./gamma(.5*o.p7(jd))).^2);
              case 'mode'
                m(jd) = o.p3(jd) + sqrt((o.p7(jd)-1)./o.p6(jd));
            end
        end
    end
    if o.isinvgamma2
        jd = intersect(o.idinvgamma2, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd)))
                error('dprior::mean: Some hyperparameters are missing (inverse gamma type 2 distribution).')
            end
            % s → o.p6, ν → o.p7
            switch name
              case 'mean'
                m(jd) =  o.p3(jd) + o.p6(jd)./(o.p7(jd)-2);
              case 'median'
                m(jd) = o.p3(jd) + 1.0/gaminv(.5, o.p7(jd)/2.0, 2.0/o.p6(jd));
              case 'std'
                m(jd) = sqrt(2./(o.p7(jd)-4)).*o.p6(jd)./(o.p7(jd)-2);
              case 'mode'
                m(jd) = o.p3(jd) + o.p6(jd)./(o.p7(jd)+2);
            end
        end
    end
    if o.isweibull
        jd = intersect(o.idweibull, find(id));
        if ~isempty(jd)
            if any(isnan(o.p6(jd))) || any(isnan(o.p7(jd))) || any(isnan(o.p3(jd)))
                error('dprior::mean: Some hyperparameters are missing (weibull distribution).')
            end
            % k → o.p6 (shape parameter), λ → o.p7 (scale parameter)
            % See https://en.wikipedia.org/wiki/Weibull_distribution
            switch name
              case 'mean'
                m(jd) =  o.p3(jd) + o.p7(jd).*gamma(1+1./o.p6(jd));
              case 'median'
                m(jd) = o.p3(jd) + o.p7(jd).*log(2).^(1./o.p6(jd));
              case 'std'
                m(jd) = o.p7(jd).*sqrt(gamma(1+2./o.p6(jd))-gamma(1+1./o.p6(jd)).^2);
              case 'mode'
                m(jd) = 0;
                hd = o.p6(jd)>1;
                m(jd(hd)) = o.p3(jd(hd)) + o.p7(jd(hd)).*((o.p6(jd(hd))-1)./o.p6(jd(hd))).^(1./o.p6(jd(hd)));
            end
        end
    end
    switch name
      case 'mean'
        o.p1 = m;
      case 'median'
        o.p11 = m;
      case 'std'
        o.p2 = m;
      case 'mode'
        o.p5 = m;
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

% Call the tested routine
try
    Prior = dprior(BayesInfo, prior_trunc, false);
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    t(2) = all(Prior.mean()==.4);
    t(3) = all(ismembertol(Prior.mean(true),.4));
    t(4) = all(ismembertol(Prior.variance(),.04));
    t(5) = all(ismembertol(Prior.variance(true),.04));
end
T = all(t);
%@eof:1

%@test:2
%  Verifies median values analytically for every distribution type.
try
    Prior = dprior();
    n = 7;
    %              Gamma  Gamma(shift) Gamma   Weibull Weibull(shift) InvGamma1  Gaussian
    Prior.p3  = [  0;     1.5;         0;      0;      0.5;           0;        -Inf];
    Prior.p6  = [  2.0;   3.0;         1.5;    2.0;    3.0;           2.0;       0.5]; % α / k / s / μ
    Prior.p7  = [  1.0;   2.0;         0.5;    1.0;    0.5;           4.0;       0.2]; % β / λ / ν / σ
    Prior.p11 = NaN(n, 1);
    Prior.idgamma     = [1; 2; 3];
    Prior.isgamma     = true;
    Prior.idweibull   = [4; 5];
    Prior.isweibull   = true;
    Prior.idinvgamma1 = 6;
    Prior.isinvgamma1 = true;
    Prior.idgaussian  = 7;
    Prior.isgaussian  = true;
    m = Prior.median(true);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    expected = [ ...
        gaminv(.5, 2.0, 1.0); ...                        % Gamma, p3=0, α=2, β=1
        1.5 + gaminv(.5, 3.0, 2.0); ...                  % Gamma shifted by 1.5, α=3, β=2
        gaminv(.5, 1.5, 0.5); ...                         % Gamma, p3=0, α=1.5, β=0.5
        log(2).^(1/2.0); ...                              % Weibull, p3=0, k=2, λ=1
        0.5 + 0.5.*log(2).^(1/3.0); ...                  % Weibull shifted by 0.5, k=3, λ=0.5
        1.0/sqrt(gaminv(.5, 4.0/2.0, 2.0/2.0)); ...      % InvGamma1, p3=0, s=2, ν=4
        0.5];                                             % Gaussian, median = μ
    t(2) = all(abs(m - expected) < 1e-10);
end
T = all(t);
%@eof:2

%@test:3
% Integration tests for beta distribution with vector inputs
% Covers: 5 Beta parameters exercising all 5 mode branches (h1–h4, h0) including
% a shifted Beta (p3≠0), plus a Gamma at index 6 so that idbeta and idgamma
% are both active and intersect() receives a proper subset of indices.
try
    Prior = dprior();
    n = 6;
    %            h1-β   h2-β   h3-β   h4-β   h0-β(shift)  Gamma
    Prior.p3 = [  0;     0;     0;     0;     0.2;         0  ];
    Prior.p4 = [  1;     1;     1;     1;     0.8;         Inf];
    Prior.p6 = [  0.5;   2.0;   0.7;   1.0;   3.0;         2.5]; % α (beta) / α (gamma)
    Prior.p7 = [  2.0;   0.5;   0.8;   1.0;   2.0;         1.5]; % β (beta) / β (gamma)
    Prior.p5 = NaN(n, 1);
    Prior.idbeta  = [1; 2; 3; 4; 5];
    Prior.isbeta  = true;
    Prior.idgamma = 6;
    Prior.isgamma = true;
    m = Prior.mode(true);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    expected = zeros(n, 1);
    expected(1) = 0;                                         % h1: α≤1 & β>1  → lower bound
    expected(2) = 1;                                         % h2: β≤1 & α>1  → upper bound
    expected(3) = 0;                                         % h3: α<1 & β<1  → lower bound (bimodal)
    expected(4) = 0 + .5*(1-0);                              % h4: α=β=1      → midpoint (uniform)
    expected(5) = 0.2 + (0.8-0.2)*(3.0-1)/(3.0+2.0-2);     % h0: interior    → scaled formula, shifted
    expected(6) = (2.5-1)*1.5;                               % Gamma mode=(α-1)·β
    t(2) = all(abs(m - expected) < 1e-10);
end
T = all(t);
%@eof:3

%@test:4
% Test that the location parameter p3 shifts mean, median and mode by exactly
% delta for all distribution families where the formulas include p3.

delta = 0.7;
try
    n = 5;
    %           Gamma  Beta   InvGamma1  InvGamma2  Weibull
    p3_base = [ 0;     0;     0;         0;         0    ];
    p4_base = [ Inf;   1;     Inf;       Inf;       Inf  ];
    p6      = [ 3.0;   3.0;   2.0;       2.0;       2.0  ];  % α, α, s, s, k
    p7      = [ 0.5;   2.0;   5.0;       6.0;       1.0  ];  % β, β, ν, ν, λ

    % Unshifted object
    A = dprior();
    A.p3 = p3_base; A.p4 = p4_base; A.p6 = p6; A.p7 = p7;
    A.p1 = NaN(n,1); A.p11 = NaN(n,1); A.p5 = NaN(n,1);
    A.idgamma = 1;   A.isgamma     = true;
    A.idbeta  = 2;   A.isbeta      = true;
    A.idinvgamma1 = 3; A.isinvgamma1 = true;
    A.idinvgamma2 = 4; A.isinvgamma2 = true;
    A.idweibull   = 5; A.isweibull   = true;

    % Shifted object: add delta to p3 for all; add delta to p4 for Beta too
    p4_shifted = p4_base;
    p4_shifted(2) = p4_base(2) + delta;
    B = dprior();
    B.p3 = p3_base + delta; B.p4 = p4_shifted; B.p6 = p6; B.p7 = p7;
    B.p1 = NaN(n,1); B.p11 = NaN(n,1); B.p5 = NaN(n,1);
    B.idgamma = 1;   B.isgamma     = true;
    B.idbeta  = 2;   B.isbeta      = true;
    B.idinvgamma1 = 3; B.isinvgamma1 = true;
    B.idinvgamma2 = 4; B.isinvgamma2 = true;
    B.idweibull   = 5; B.isweibull   = true;

    mean_A   = A.mean(true);   mean_B   = B.mean(true);
    median_A = A.median(true); median_B = B.median(true);
    mode_A   = A.mode(true);   mode_B   = B.mode(true);
    t(1) = true;
catch
    t(1) = false;
end
if t(1)
    % All means shift by delta (covers Gamma, Beta, InvGamma1, InvGamma2, Weibull)
    t(2) = all(abs((mean_B   - mean_A)   - delta) < 1e-10);
    % All medians shift by delta
    t(3) = all(abs((median_B - median_A) - delta) < 1e-10);
    % Mode shifts by delta for all 5 distributions:
    %   Gamma (α=3>1), Beta (h0: α>1,β>1), InvGamma1 (ν=5), InvGamma2 (ν=2), Weibull (k=2>1)
    t(4) = all(abs((mode_B - mode_A) - delta) < 1e-10);
end
T = all(t);
%@eof:4
