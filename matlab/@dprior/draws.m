function P = draws(o, n)

% Return n independent random draws from the prior distribution.
%
% INPUTS
% - o    [dprior]
%
% OUTPUTS
% - P    [double]   m×n matrix, random draw from the prior distribution.
%
% REMARKS
% If the Parallel Computing Toolbox is available, the main loop is run in parallel.
%
% EXAMPLE
%
% >> Prior = dprior(bayestopt_, options_.prior_trunc);
% >> Prior.draws(1e6)

% Copyright © 2023-2026 Dynare Team
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

P = NaN(rows(o.lb), n);
if ~isoctave
    % Derive a deterministic seed from the current global RNG state. This ensures:
    % (1) Reproducibility: same initial RNG state produces identical draws
    % (2) Progression: successive calls yield different draws (via rand(1,n) at the end)
    % (3) Serial/parallel equivalence: both code paths use the same seed and substream scheme
    % The hash uses position-weighted sum to avoid collisions from permuted states.
    global_stream = RandStream.getGlobalStream();
    state = uint64(global_stream.State);
    idx = uint64(1:numel(state))';
    base_seed = double(mod(sum(state.*idx), uint64(2^32)));
    if matlab.internal.parallel.isPCTInstalled
        % Use the same stream/substream scheme as serial to keep results identical across modes.
        sc = parallel.pool.Constant(RandStream('Threefry','Seed',base_seed)); % create a constant stream for the worker pool
        parfor i=1:n
            % RandStream.setGlobalStream must be called inside the parfor body, not before it.
            % Each worker is a separate MATLAB process with its own independent global stream
            % state. Calling setGlobalStream on the client (outside parfor) would only affect
            % the client process and have no effect on the workers, which would then each use
            % their default stream and produce non-reproducible results.
            stream = sc.Value; % use substream per iteration to make parfor scheduling deterministic
            stream.Substream = i; % set the substream for the current iteration
            RandStream.setGlobalStream(stream); % set the stream on the worker process
            P(:,i) = draw(o);
        end
    else
        % In the serial case there is only one process, so it is sufficient to set the global
        % stream once before the loop. Inside the loop only the cheap substream index is
        % advanced, avoiding the overhead of replacing the global stream object on every
        % iteration.
        stream = RandStream('Threefry','Seed',base_seed);
        RandStream.setGlobalStream(stream);
        for i=1:n
            stream.Substream = i;
            P(:,i) = draw(o);
        end
    end
    RandStream.setGlobalStream(global_stream);
    if n > 0
        % Advance the global RNG so successive calls differ even though draws used worker-local streams.
        rand(1, n);
    end
else
    for i=1:n
        P(:,i) = draw(o);
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

% Required options for set_dynare_seed_local_options
options_.DynareRandomStreams.seed = 123;
options_.DynareRandomStreams.algo = 'mt19937ar';
options_.parallel_info.isHybridMatlabOctave = false;

ndraws = 1e5;

% Call the tested routine
try
    % Instantiate dprior object.
    o = dprior(BayesInfo, prior_trunc, false);
    X = o.draws(ndraws);
    m = mean(X, 2);
    v = var(X, 0, 2);
    t(1) = true;
catch
    t(1) = false;
end

if t(1)
    % Check sample mean matches theoretical mean
    t(2) = all(abs(m-BayesInfo.p1)<3e-3);
    % Check sample variance matches theoretical variance
    t(3) = all(all(abs(v-BayesInfo.p2.^2)<5e-3));
    try
        % Different seeds should produce different draws
        set_dynare_seed_local_options(options_.DynareRandomStreams,options_.parallel_info.isHybridMatlabOctave,456);
        X1 = o.draws(20);
        set_dynare_seed_local_options(options_.DynareRandomStreams,options_.parallel_info.isHybridMatlabOctave,457);
        X2 = o.draws(20);
        t(4) = ~isequaln(X1, X2);
    catch
        t(4) = false;
    end
    try
        % Re-seeding should reproduce identical draws
        set_dynare_seed_local_options(options_.DynareRandomStreams,options_.parallel_info.isHybridMatlabOctave,456);
        X3 = o.draws(20);
        t(5) = isequaln(X1, X3);
    catch
        t(5) = false;
    end
    try
        % Successive calls without reseeding should differ
        set_dynare_seed_local_options(options_.DynareRandomStreams,options_.parallel_info.isHybridMatlabOctave,111);
        X4 = o.draws(20);
        X5 = o.draws(20);
        t(6) = ~isequaln(X4, X5);
    catch
        t(6) = false;
    end
    if ~isoctave && matlab.internal.parallel.isPCTInstalled
        try
            % Parallel mode should be reproducible after reseeding
            set_dynare_seed_local_options(options_.DynareRandomStreams,options_.parallel_info.isHybridMatlabOctave,321);
            X7 = o.draws(20);
            set_dynare_seed_local_options(options_.DynareRandomStreams,options_.parallel_info.isHybridMatlabOctave,321);
            X8 = o.draws(20);
            t(7) = isequaln(X7, X8);
        catch
            t(7) = false;
        end
        try
            % Serial and parallel should produce identical results with same seed
            ntest = 20;
            set_dynare_seed_local_options(options_.DynareRandomStreams,options_.parallel_info.isHybridMatlabOctave,789);
            % Compute draws using serial code path (same logic as draws.m but with for loop)
            global_stream_test = RandStream.getGlobalStream();
            state_test = uint64(global_stream_test.State);
            idx_test = uint64(1:numel(state_test))';
            base_seed_test = double(mod(sum(state_test.*idx_test), uint64(2^32)));
            stream_serial = RandStream('Threefry','Seed',base_seed_test);
            X_serial = NaN(rows(o.lb), ntest);
            for ii=1:ntest
                stream_serial.Substream = ii;
                RandStream.setGlobalStream(stream_serial);
                X_serial(:,ii) = o.draw();
            end
            RandStream.setGlobalStream(global_stream_test);
            rand(1, ntest); % advance global RNG as draws() does
            % Now compute draws using parallel code path (via draws method)
            set_dynare_seed_local_options(options_.DynareRandomStreams,options_.parallel_info.isHybridMatlabOctave,789);
            X_parallel = o.draws(ntest);
            t(8) = isequaln(X_serial, X_parallel);
        catch
            t(8) = false;
        end
    else
        t(7) = true; % skip parallel mode test on Octave
        t(8) = true; % skip serial vs parallel equivalence test on Octave
    end
end
T = all(t);
%@eof:1
