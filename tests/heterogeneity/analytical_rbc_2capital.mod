// Two-sector analytical RBC model with 2D endogenous state space
// This model validates time iteration and forward iteration algorithms
// for heterogeneous-agent models with multiple endogenous state variables.
//
// Key features:
// - Two independent sectors with separate consumption goods
// - Full depreciation (delta = 1) allows closed-form policy functions
// - 2D endogenous state space: (k1, k2)
// - 1D exogenous state space: productivity shock z

heterogeneity_dimension households;

var(heterogeneity=households) k1 k2 y1 y2 c1 c2 rk1 rk2;

varexo(heterogeneity=households) z;

parameters alpha1 alpha2 beta l1 l2 rho sig
           ss_k1 ss_k2 ss_y1 ss_y2 ss_c1 ss_c2 ss_rk1 ss_rk2;

// Calibration
alpha1  = 0.33;     // Capital share sector 1
alpha2  = 0.30;     // Capital share sector 2
beta    = 0.96;     // Discount factor
l1      = 0.33;     // Labor supply sector 1 (fixed)
l2      = 0.33;     // Labor supply sector 2 (fixed)
rho     = 0.95;     // TFP persistence
sig     = 0.01;     // TFP shock std dev

// Analytical steady-state values (at z=1)
ss_k1 = (beta * alpha1)^(1/(1-alpha1)) * l1;
ss_k2 = (beta * alpha2)^(1/(1-alpha2)) * l2;
ss_y1 = ss_k1^alpha1 * l1^(1-alpha1);
ss_y2 = ss_k2^alpha2 * l2^(1-alpha2);
ss_c1 = (1 - beta * alpha1) * ss_y1;
ss_c2 = (1 - beta * alpha2) * ss_y2;
ss_rk1 = alpha1 * ss_y1 / ss_k1;
ss_rk2 = alpha2 * ss_y2 / ss_k2;

model(heterogeneity=households);
    // Sector 1
    y1  = z*k1(-1)^alpha1*l1^(1-alpha1);
    rk1 = alpha1*y1/k1(-1);
    y1  = c1 + k1;
    1/c1 - beta*rk1(+1)/c1(+1) = 0;

    // Sector 2
    y2  = z*k2(-1)^alpha2*l2^(1-alpha2);
    rk2 = alpha2*y2/k2(-1);
    y2  = c2 + k2;
    1/c2 - beta*rk2(+1)/c2(+1) = 0;
end;

%% ======== SHOCK DISCRETIZATION ============
M = 5;  // Number of TFP states
max_iter = options_.heterogeneity.forward.max_iter;
tol = options_.heterogeneity.forward.tol;
[z_grid, p, Q] = rouwenhorst(rho, sig, M, tol, max_iter);

%% ======== CAPITAL GRIDS ============
N1 = 20;  // Grid points for k1
N2 = 20;  // Grid points for k2
k1_grid = linspace(0.8*ss_k1, 1.2*ss_k1, N1);
k2_grid = linspace(0.8*ss_k2, 1.2*ss_k2, N2);

%% ======== ANALYTICAL SOLUTION ============
% Compute closed-form policy functions on 3D state space
y1_th  = zeros(M, N1, N2);
y2_th  = zeros(M, N1, N2);
c1_th  = zeros(M, N1, N2);
c2_th  = zeros(M, N1, N2);
k1_th  = zeros(M, N1, N2);
k2_th  = zeros(M, N1, N2);
rk1_th = zeros(M, N1, N2);
rk2_th = zeros(M, N1, N2);

for i_z = 1:M
    Z = z_grid(i_z);
    for i_k1 = 1:N1
        K1_lag = k1_grid(i_k1);  % k1(-1) is the state variable
        for i_k2 = 1:N2
            K2_lag = k2_grid(i_k2);  % k2(-1) is the state variable

            % Sector 1 (independent of sector 2)
            Y1 = Z * K1_lag^alpha1 * l1^(1-alpha1);
            k1_th(i_z, i_k1, i_k2)  = beta * alpha1 * Y1;
            c1_th(i_z, i_k1, i_k2)  = (1 - beta * alpha1) * Y1;
            y1_th(i_z, i_k1, i_k2)  = Y1;
            rk1_th(i_z, i_k1, i_k2) = alpha1 * Y1 / K1_lag;

            % Sector 2 (independent of sector 1)
            Y2 = Z * K2_lag^alpha2 * l2^(1-alpha2);
            k2_th(i_z, i_k1, i_k2)  = beta * alpha2 * Y2;
            c2_th(i_z, i_k1, i_k2)  = (1 - beta * alpha2) * Y2;
            y2_th(i_z, i_k1, i_k2)  = Y2;
            rk2_th(i_z, i_k1, i_k2) = alpha2 * Y2 / K2_lag;
        end
    end
end

%% ======== STEADY-STATE STRUCTURE ============
steady_state = struct;
steady_state.shocks = struct;
steady_state.shocks.grids.z = z_grid;
steady_state.shocks.Pi.z = Q;

steady_state.pol = struct;
steady_state.pol.grids.k1 = k1_grid;
steady_state.pol.grids.k2 = k2_grid;
steady_state.pol.values.k1  = ss_k1 * ones(M, N1, N2);
steady_state.pol.values.k2  = ss_k2 * ones(M, N1, N2);
steady_state.pol.values.y1  = ss_y1 * ones(M, N1, N2);
steady_state.pol.values.y2  = ss_y2 * ones(M, N1, N2);
steady_state.pol.values.c1  = ss_c1 * ones(M, N1, N2);
steady_state.pol.values.c2  = ss_c2 * ones(M, N1, N2);
steady_state.pol.values.rk1 = ss_rk1 * ones(M, N1, N2);
steady_state.pol.values.rk2 = ss_rk2 * ones(M, N1, N2);
steady_state.pol.values.AUX_HET_EXPECT_15 = beta*ss_rk1/ss_c1;  % Auxiliary for sector 1 Euler
steady_state.pol.values.AUX_HET_EXPECT_19 = beta*ss_rk2/ss_c2;  % Auxiliary for sector 2 Euler
steady_state.pol.order = {'z', 'k1', 'k2'};

%% ======== TIME ITERATION CHECK ============
fprintf('\n=== TIME ITERATION VALIDATION ===\n');
oo_het = struct;
oo_het = heterogeneity.compute_steady_state(M_, options_.heterogeneity, oo_het, steady_state);

% Extract time-iterated policies
k1_ti  = oo_het.steady_state.pol.values.k1;
k2_ti  = oo_het.steady_state.pol.values.k2;
c1_ti  = oo_het.steady_state.pol.values.c1;
c2_ti  = oo_het.steady_state.pol.values.c2;
y1_ti  = oo_het.steady_state.pol.values.y1;
y2_ti  = oo_het.steady_state.pol.values.y2;
rk1_ti = oo_het.steady_state.pol.values.rk1;
rk2_ti = oo_het.steady_state.pol.values.rk2;

% Compute relative errors
max_error_k1 = max(abs(k1_th(:) - k1_ti(:)) ./ k1_th(:));
max_error_k2 = max(abs(k2_th(:) - k2_ti(:)) ./ k2_th(:));
max_error_c1 = max(abs(c1_th(:) - c1_ti(:)) ./ c1_th(:));
max_error_c2 = max(abs(c2_th(:) - c2_ti(:)) ./ c2_th(:));

fprintf('Max error in k1''(Z,K1,K2): %e\n', max_error_k1);
fprintf('Max error in k2''(Z,K1,K2): %e\n', max_error_k2);
fprintf('Max error in c1(Z,K1,K2):  %e\n', max_error_c1);
fprintf('Max error in c2(Z,K1,K2):  %e\n', max_error_c2);

% This is expected for coarse grids and validates the algorithm correctly
tol = 1e-3;
if max_error_k1 > tol || max_error_k2 > tol || max_error_c1 > tol || max_error_c2 > tol
    error('Time iteration failed accuracy test: errors exceed tolerance %e', tol);
else
    fprintf('✓ Time iteration matches analytical solution within tolerance %e\n', tol);
end

%% ======== FORWARD ITERATION CHECK ============
fprintf('\n=== FORWARD ITERATION VALIDATION ===\n');

% Build sparse transition matrix T
% T(dest, src) = probability of going from src to dest
% State indexing: linear index = ((i_k2 - 1)*N1 + (i_k1 - 1))*M + i_z
n_states = M * N1 * N2;
rows = [];
cols = [];
vals = [];

for i_z = 1:M
    for i_k1 = 1:N1
        for i_k2 = 1:N2
            src_idx = ((i_k2 - 1)*N1 + (i_k1 - 1))*M + i_z;

            % Next-period capitals from time-iterated policy
            k1p = k1_ti(i_z, i_k1, i_k2);
            k2p = k2_ti(i_z, i_k1, i_k2);

            % Bilinear interpolation: find brackets and weights
            [i1_lo, w1_lo] = find_bracket_linear_weight(k1_grid, k1p);
            [i2_lo, w2_lo] = find_bracket_linear_weight(k2_grid, k2p);

            i1_hi = i1_lo + 1;
            i2_hi = i2_lo + 1;
            w1_hi = 1 - w1_lo;
            w2_hi = 1 - w2_lo;

            % Distribute mass over shock transitions and 4 spatial corners
            for j_z = 1:M
                prob_z = Q(i_z, j_z);

                % Corner 1: (i1_lo, i2_lo)
                dest_idx = ((i2_lo - 1)*N1 + (i1_lo - 1))*M + j_z;
                rows = [rows; dest_idx];
                cols = [cols; src_idx];
                vals = [vals; prob_z * w1_lo * w2_lo];

                % Corner 2: (i1_hi, i2_lo)
                dest_idx = ((i2_lo - 1)*N1 + (i1_hi - 1))*M + j_z;
                rows = [rows; dest_idx];
                cols = [cols; src_idx];
                vals = [vals; prob_z * w1_hi * w2_lo];

                % Corner 3: (i1_lo, i2_hi)
                dest_idx = ((i2_hi - 1)*N1 + (i1_lo - 1))*M + j_z;
                rows = [rows; dest_idx];
                cols = [cols; src_idx];
                vals = [vals; prob_z * w1_lo * w2_hi];

                % Corner 4: (i1_hi, i2_hi)
                dest_idx = ((i2_hi - 1)*N1 + (i1_hi - 1))*M + j_z;
                rows = [rows; dest_idx];
                cols = [cols; src_idx];
                vals = [vals; prob_z * w1_hi * w2_hi];
            end
        end
    end
end

T = sparse(rows, int32(cols), vals, n_states, n_states);

% Forward iteration
tol = 1e-10;
max_iter = 10000;

% Initialize: uniform over (k1,k2), stationary over z
D = zeros(M, N1, N2);
for i_z = 1:M
    D(i_z, :, :) = p(i_z) / (N1 * N2);
end
D = D(:);  % Flatten to vector (z varies fastest, then k1, then k2)

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

% Reshape and display results
D_mat = reshape(D, M, N1, N2);  % D_mat(i_z, i_k1, i_k2) = distribution mass

fprintf('\n=== Distribution Summary ===\n');
fprintf('Distribution sums to: %f\n', sum(D));
fprintf('Marginal over z: %s\n', mat2str(squeeze(sum(sum(D_mat, 2), 3))', 4));
fprintf('Expected (stationary p): %s\n', mat2str(p', 4));

% Mean capitals
mean_k1 = 0;
mean_k2 = 0;
for i_z = 1:M
    for i_k1 = 1:N1
        for i_k2 = 1:N2
            mean_k1 = mean_k1 + D_mat(i_z, i_k1, i_k2) * k1_grid(i_k1);
            mean_k2 = mean_k2 + D_mat(i_z, i_k1, i_k2) * k2_grid(i_k2);
        end
    end
end
fprintf('Mean k1: %f (ss_k1 = %f)\n', mean_k1, ss_k1);
fprintf('Mean k2: %f (ss_k2 = %f)\n', mean_k2, ss_k2);

% Compare to MEX implementation
max_error_d = max(abs(oo_het.steady_state.d.hist(:) - D_mat(:)));
tol = 1e-9;
if max_error_d > tol
    error('Forward iteration failed accuracy test: errors exceed tolerance %e', tol);
else
    fprintf('✓ Forward iteration matches MEX computation within tolerance %e\n', tol);
end

fprintf('\n=== ALL VALIDATION CHECKS PASSED ===\n');
