heterogeneity_dimension households;

var(heterogeneity=households) k y c rk;

varexo(heterogeneity=households) z;

parameters alpha beta theta rho sig
           ss_l ss_k ss_c ss_y ss_rk;

alpha   = 1/3;
beta    = 0.99;
rho     = 0.95;
sig     = 0.01;
ss_l    = 1/3;
ss_k    = ss_l*(beta*alpha)^(1/(1-alpha));
ss_y    = ss_k^alpha*ss_l^(1-alpha);
ss_c    = ss_y-ss_k;
ss_rk   = alpha*ss_y/ss_k;
theta   = (1-alpha)*(1-ss_l)*ss_y/(ss_l*ss_c);

model(heterogeneity=households); 
    1/c - beta*rk(+1)/c(+1) = 0;
    y   = z*k(-1)^alpha*ss_l^(1-alpha);
    rk  = alpha*y/k(-1);
    y   = c+k;
end;

M = 10;
N = 30;
max_iter = options_.heterogeneity.forward.max_iter;
tol = options_.heterogeneity.forward.tol;
[z_grid, p, Q] = rouwenhorst(rho, sig, M, tol, max_iter);
k_grid = linspace(0.8*ss_k, 1.2*ss_k, N);
c_th = zeros(M, N);
k_th = zeros(M, N);
for i_z = 1:M
    Z = z_grid(i_z);
    for i_k = 1:N
        K = k_grid(i_k);
        Y = Z * K^alpha * ss_l^(1-alpha);
        c_th(i_z, i_k) = (1 - alpha*beta) * Y;
        k_th(i_z, i_k) = alpha * beta * Y;
    end
end

steady_state = struct;
steady_state.shocks = struct;
steady_state.shocks.grids.z = z_grid;
steady_state.shocks.Pi.z = Q;
steady_state.pol.grids.k = k_grid
steady_state.pol.values.k = ss_k*ones(M,N);
steady_state.pol.values.y = ss_y*ones(M,N);
steady_state.pol.values.c = ss_c*ones(M,N);
steady_state.pol.values.rk = ss_rk*ones(M,N);
steady_state.pol.values.AUX_HET_EXPECT_15 = beta*ss_rk/ss_c;
steady_state.pol.order = {'z', 'k'};

%% ======== TIME ITERATION CHECK ============
oo_het = struct;
oo_het = heterogeneity.compute_steady_state(M_, options_.heterogeneity, oo_het, steady_state);

c_ti = oo_het.steady_state.pol.values.c;
k_ti = oo_het.steady_state.pol.values.k;

% Compute errors
max_error_c = max(abs(c_th(:) - c_ti(:)) ./ c_th(:));
max_error_k = max(abs(k_th(:) - k_ti(:)) ./ k_th(:));
fprintf('Max error in c(K,Z): %e\n', max_error_c);
fprintf("Max error in k\'(K,Z): %e\n", max_error_k);

tol = 1e-4;
if max_error_c > tol || max_error_k > tol
    error('Time iteration failed accuracy test: errors exceed tolerance %e', tol);
else
    fprintf('Time iteration matches analytical solution within tolerance %e\n', tol);
end

%% ======== FORWARD ITERATION CHECK ============
%% Build sparse transition matrix T
% T(dest, src) = probability of going from src to dest
% State indexing: linear index = (i_k - 1)*M + i_z for (z_i, k_i)
% This matches Fortran column-major with d_mat(1:N_e, 1:N_a) where z varies fastest
n_states = M * N;
rows = [];
cols = [];
vals = [];

for i_z = 1:M
    for i_k = 1:N
        src_idx = (i_k - 1)*M + i_z;
        kp = k_ti(i_z, i_k);  % time-iterated policy (same as Fortran uses)

        % Use find_bracket_linear_weight MEX for bracketing
        % Returns: i_lo (1-based index), w_lo (weight on lower bracket)
        [i_lo, w_lo] = find_bracket_linear_weight(k_grid, kp);

        % Handle boundary: if kp >= k_grid(end), i_lo = N and w_lo = 1
        % We need i_hi = i_lo + 1, but cap at N
        if i_lo >= N
            i_lo = N - 1;
            w_lo = 0.0;  % all weight on upper (index N)
        end
        i_hi = i_lo + 1;
        w_hi = 1 - w_lo;

        % Distribute mass to destination states for each z' realization
        for j_z = 1:M
            prob_z = Q(i_z, j_z);  % shock transition prob P(z'|z)

            % Lower bracket contribution
            if w_lo > 0
                dest_lo = (i_lo - 1)*M + j_z;
                rows = [rows; dest_lo];
                cols = [cols; src_idx];
                vals = [vals; prob_z * w_lo];
            end

            % Upper bracket contribution
            if w_hi > 0
                dest_hi = (i_hi - 1)*M + j_z;
                rows = [rows; dest_hi];
                cols = [cols; src_idx];
                vals = [vals; prob_z * w_hi];
            end
        end
    end
end

T = sparse(rows, int32(cols), vals, n_states, n_states);

%% Forward iteration
tol = 1e-10;
max_iter = 10000;

% Initialize: uniform over k, stationary over z
% D_mat(i_z, i_k) with z varying fastest when flattened
D = zeros(M, N);
for i_z = 1:M
    D(i_z, :) = p(i_z) / N;
end
D = D(:);  % flatten to vector (column-major: z varies fastest, then k)

for iter = 1:max_iter
    D_next = T * D;
    diff_norm = max(abs(D_next - D));
    if diff_norm < tol
        fprintf('Converged in %d iterations (diff = %e)\n', iter, diff_norm);
        break;
    end
    D = D_next;
end

if iter == max_iter
    fprintf('Warning: did not converge (diff = %e)\n', diff_norm);
end

%% Reshape and display results
D_mat = reshape(D, M, N);  % D_mat(i_z, i_k) = distribution mass

fprintf('\n=== Distribution Summary ===\n');
fprintf('Distribution sums to: %f\n', sum(D));
fprintf('Marginal over z: %s\n', mat2str(sum(D_mat, 2)', 4));
fprintf('Expected (stationary p): %s\n', mat2str(p', 4));
% Mean capital: sum over z of (D_mat(i_z, :) * k_grid')
mean_k = sum(D_mat * k_grid');
fprintf('Mean capital: %f (ss_k = %f)\n', mean_k, ss_k);

max_error_d = max(abs(oo_het.steady_state.d.hist(:)-D_mat(:)));
tol = 1e-9;
if (max_error_d > tol)
    error('Forward iteration failed accuracy test: errors exceed tolerance %e', tol);
else
    fprintf('Forward iteration matches the manual computation within tolerance %e\n', tol);
end