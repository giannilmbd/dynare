function Phi_tilde_e = build_Phi_tilde_e(sizes, indices, mat)
% Builds interpolation matrix for heterogeneous agent transitions.
%
% INPUTS
% - sizes [structure]: size information for various dimensions
% - indices [structure]: permutation and indexing information
% - mat [structure]: matrices from steady-state computation
%
% OUTPUTS
% - Phi_tilde_e [sparse matrix]: interpolation matrix (N_sp × N_sp)
%
% DESCRIPTION
% Returns sparse Phi_tilde_e (N_sp × N_sp) for linear interpolation in all states
% with θ mixing by Mu(jθ, ℓθ). Columns j are in θ-fastest order: j = (j_a-1)*N_e + j_e.
% This implementation uses a Gray walk for index updates and a zero-safe incremental
% update for weights: it tracks per row (i) the count of zero factors among selected
% dimensions and (ii) the sum of logs of the non-zero factors. This makes the updates
% robust when some interpolation weights are exactly 0 or 1.

% Copyright © 2025 Dynare Team
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
%
% Original author: Normann Rion <normann@dynare.org>
% Builds interpolation matrix Φ̃ₑ (N_sp × N_sp) with θ mixing by Mu.
% Robust Gray-walk with zero-safe incremental weights using the assumption:
%   - At the "all lowers" corner, no lower weight is zero.
%   - Zeros can only appear when flipping a dim with L==1 (thus U==0).

    states = indices.states;           % fastest-first among state dims
    n  = sizes.n_a;                    % # state dims
    N_e = sizes.N_e;

    % Policy grid sizes
    mk = zeros(1,n);
    for k = 1:n
        mk(k) = sizes.pol.states.(states{k});
    end
    N_a  = prod(mk);
    N_sp = N_e * N_a;

    % ---- Per-state bracketing (flattened with θ fastest) ----
    i0 = zeros(N_sp, n);
    w  = zeros(N_sp, n);
    for k = 1:n
        i0(:,k) = double(mat.pol.ind.(states{k})(:));   % lower index
        w(:,k)  = mat.pol.w.(states{k})(:);             % lower weight in [0,1]
    end

    % Optional: tiny tolerance to classify near-1 as hard; near-0 not expected initially
    L = w;           % lower-side factors
    U = 1 - w;       % upper-side factors
    hard = (L == 1); % dims where upper is exactly zero when flipped
    % soft dims have 0 < L < 1 and 0 < U < 1

    % ---- Strides and base states-only linear index (1..N_a) ----
    stride = ones(n,1);
    for k = 2:n, stride(k) = stride(k-1) * mk(k-1); end
    jidx = 1 + sum( (i0 - 1) .* stride.', 2 );          % N_sp×1

    % θ index per column j (θ fastest)
    j_e = repmat((1:N_e).', N_a, 1);                    % N_sp×1

    % ---- Gray code over n dims ----
    Kcorn   = bitshift(1, n);                           % 2^n
    g       = bitxor(uint32(0:Kcorn-1), bitshift(uint32(0:Kcorn-1), -1));
    flipIdx = zeros(1, Kcorn-1, 'uint8');
    for t = 2:Kcorn
        x = bitxor(g(t), g(t-1));
        for kk = 1:n
            if bitget(x, kk), flipIdx(t-1) = kk; break; end
        end
    end

    % θ mixing rows: for each column j, we need Mu(j_e(j), :)
    Mu_rows = mat.Mu(j_e, :);                           % N_sp×N_e

    % ---- Sparse triplets allocation: each corner adds N_e entries per column ----
    TOT = Kcorn * N_e * N_sp;
    I = zeros(TOT,1);
    J = zeros(TOT,1);
    V = zeros(TOT,1);
    ptr = 0;

    % Emit helper
    function emit_block(jidx_here, beta_here)
        base = (jidx_here - 1) * N_e;                   % N_sp×1
        R = base + (1:N_e);                             % N_sp×N_e
        Jm = (1:N_sp).' * ones(1, N_e);                 % N_sp×N_e
        Vm = Mu_rows .* beta_here;                      % N_sp×N_e
        nb = numel(R);
        I(ptr+1:ptr+nb) = reshape(R.', [], 1);
        J(ptr+1:ptr+nb) = reshape(Jm.', [], 1);
        V(ptr+1:ptr+nb) = reshape(Vm.', [], 1);
        ptr = ptr + nb;
    end

    % ---- Initialize Gray walk state ----
    sk = false(1,n);            % all lowers
    p = prod(L, 2);

    % z = count of active hard dims flipped to upper (i.e., contributing 0)
    z = zeros(N_sp,1); % initially 0 since all lowers

    % beta = p when z==0, else 0
    beta = p;                   % all rows have z==0 initially

    % Emit initial corner (all lowers)
    emit_block(jidx, beta);

    % ---- Walk remaining corners ----
    for t = 1:Kcorn-1
        kf = flipIdx(t);

        % Partition rows by whether this dim is hard or soft
        is_hard = hard(:, kf);
        is_soft = ~is_hard;

        if ~sk(kf)
            % lower -> upper on dim kf
            % Indices
            jidx = jidx + stride(kf);
            sk(kf) = true;

            % --- Hard rows: L==1, U==0 -> these rows become zeroed
            if any(is_hard)
                z(is_hard) = z(is_hard) + 1;
                beta(is_hard) = 0;               % zero out
            end

            % --- Soft rows: safe ratio update (both sides > 0)
            if any(is_soft)
                ratio = U(is_soft, kf) ./ L(is_soft, kf);
                p(is_soft) = p(is_soft) .* ratio;
                % Only apply to beta where z==0
                ok = is_soft & (z == 0);
                if any(ok)
                    beta(ok) = beta(ok) .* ratio(ok(is_soft)); % align indexing
                end
            end

        else
            % upper -> lower on dim kf
            % Indices
            jidx = jidx - stride(kf);
            sk(kf) = false;

            % --- Hard rows: leaving zeroed state if this was the last active zero
            if any(is_hard)
                z(is_hard) = z(is_hard) - 1;
                resurrect = is_hard & (z == 0);
                if any(resurrect)
                    beta(resurrect) = p(resurrect); % resurrect from p
                end
            end

            % --- Soft rows: safe inverse ratio
            if any(is_soft)
                ratio = L(is_soft, kf) ./ U(is_soft, kf);
                p(is_soft) = p(is_soft) .* ratio;
                ok = is_soft & (z == 0);
                if any(ok)
                    beta(ok) = beta(ok) .* ratio(ok(is_soft));
                end
            end
        end

        emit_block(jidx, beta);
    end

    % Assemble sparse Φ̃ₑ
    Phi_tilde_e = sparse(I, J, V, N_sp, N_sp);
end