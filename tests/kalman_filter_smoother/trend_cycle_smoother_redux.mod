// Unobserved-components (trend + cycle) model with unit roots.
// Tests smoother_redux Ptt/V_eta reconstruction for diffuse Kalman filter.
//
// Model:
//   y = yp + z                    (static: observed = trend + cycle)
//   yp = mu + yp(-1) + eyp        (unit root: local level)
//   mu = mu(-1) + emu              (unit root: drift random walk)
//   z  = alpha*z(-1) + ez          (stationary: cycle)
//
// Variables yp, mu are non-stationary => Pinf != 0, d > 0.
// Variable y is static (no lag) => exercises smoother_redux static reconstruction.

var y yp z mu;

varexo ez eyp emu;

parameters alpha;

alpha = 0.889;

model(linear);
  y  = yp + z;
  yp = mu + yp(-1) + eyp;
  mu = mu(-1) + emu;
  z  = alpha*z(-1) + ez;
end;

initval;
  y  = 0;
  z  = 0;
  yp = 0;
  mu = 0;
end;

steady(nocheck);

shocks;
  var ez;  stderr 0.06;
  var eyp; stderr 0.002;
  var emu; stderr 0.002;
end;

varobs y;

options_.qz_criterium = 1+1e-6;

// =========================================================================
// Run 1: full smoother (no smoother_redux) — multivariate baseline
// =========================================================================
calib_smoother(datafile='../filter_step_ahead/trend_cycle_decomposition_data', diffuse_filter, filtered_vars, filter_covariance, smoothed_state_uncertainty) y yp z mu;
oo0 = oo_;

// =========================================================================
// Run 2: smoother_redux with kalman_algo=1 (multivariate diffuse)
// =========================================================================
calib_smoother(datafile='../filter_step_ahead/trend_cycle_decomposition_data', diffuse_filter, filtered_vars, filter_covariance, smoothed_state_uncertainty, smoother_redux) y yp z mu;
oo1 = oo_;

// =========================================================================
// Consistency checks: smoother_redux vs baseline (multivariate)
// =========================================================================

// --- Smoothed shocks ---
for k=1:M_.exo_nbr
    mserr(k) = max(abs(oo0.SmoothedShocks.(M_.exo_names{k}) - oo1.SmoothedShocks.(M_.exo_names{k})));
end
if max(mserr) > 1e-12
    error('smoother_redux does not replicate original smoother for shocks!')
end

// --- Smoothed variables ---
for k=1:M_.endo_nbr
    merr(k) = max(abs(oo0.SmoothedVariables.(M_.endo_names{k}) - oo1.SmoothedVariables.(M_.endo_names{k})));
end
if max(merr) > 1e-12
    error('smoother_redux does not replicate original smoothed variables!')
end

// --- Updated variables ---
for k=1:M_.endo_nbr
    merrU(k) = max(abs(oo0.UpdatedVariables.(M_.endo_names{k}) - oo1.UpdatedVariables.(M_.endo_names{k})));
end
if max(merrU) > 1e-12
    error('smoother_redux does not replicate original updated variables!')
end

// --- Filtered variables ---
for k=1:M_.endo_nbr
    merrF(k) = max(abs(oo0.FilteredVariables.(M_.endo_names{k}) - oo1.FilteredVariables.(M_.endo_names{k})));
end
if max(merrF) > 1e-12
    error('smoother_redux does not replicate original filtered variables!')
end

// --- Filter covariance (P_{t|t-1}) ---
verr = max(max(max(abs(oo0.Smoother.Variance - oo1.Smoother.Variance))));
if verr > 1e-10
    error('smoother_redux does not replicate original filter covariance! max diff = %e', verr)
end

// --- State uncertainty (V_{t|T}): smoothed state covariance ---
verrS = max(max(max(abs(oo0.Smoother.State_uncertainty - oo1.Smoother.State_uncertainty))));
if verrS > 1e-10
    error('smoother_redux does not replicate original state uncertainty! max diff = %e', verrS)
end

fprintf('All smoother_redux (kalman_algo=1) diffuse consistency checks passed.\n');

// =========================================================================
// Run 3: baseline with kalman_algo=4 (univariate diffuse, no smoother_redux)
// =========================================================================
calib_smoother(datafile='../filter_step_ahead/trend_cycle_decomposition_data', diffuse_filter, filtered_vars, filter_covariance, smoothed_state_uncertainty, kalman_algo=4) y yp z mu;
oo2 = oo_;

// =========================================================================
// Run 4: smoother_redux with kalman_algo=4 (univariate diffuse)
// =========================================================================
calib_smoother(datafile='../filter_step_ahead/trend_cycle_decomposition_data', diffuse_filter, filtered_vars, filter_covariance, smoothed_state_uncertainty, kalman_algo=4, smoother_redux) y yp z mu;
oo3 = oo_;

// =========================================================================
// Consistency checks: kalman_algo=4 smoother_redux vs kalman_algo=4 baseline
// =========================================================================

// --- Smoothed shocks ---
for k=1:M_.exo_nbr
    mserr2(k) = max(abs(oo2.SmoothedShocks.(M_.exo_names{k}) - oo3.SmoothedShocks.(M_.exo_names{k})));
end
if max(mserr2) > 1e-12
    error('kalman_algo=4: smoother_redux does not replicate original smoother for shocks!')
end

// --- Smoothed variables ---
for k=1:M_.endo_nbr
    merr2(k) = max(abs(oo2.SmoothedVariables.(M_.endo_names{k}) - oo3.SmoothedVariables.(M_.endo_names{k})));
end
if max(merr2) > 1e-12
    error('kalman_algo=4: smoother_redux does not replicate original smoothed variables!')
end

// --- Updated variables ---
for k=1:M_.endo_nbr
    merrU2(k) = max(abs(oo2.UpdatedVariables.(M_.endo_names{k}) - oo3.UpdatedVariables.(M_.endo_names{k})));
end
if max(merrU2) > 1e-12
    error('kalman_algo=4: smoother_redux does not replicate original updated variables!')
end

// --- Filtered variables ---
for k=1:M_.endo_nbr
    merrF2(k) = max(abs(oo2.FilteredVariables.(M_.endo_names{k}) - oo3.FilteredVariables.(M_.endo_names{k})));
end
if max(merrF2) > 1e-12
    error('kalman_algo=4: smoother_redux does not replicate original filtered variables!')
end

// --- Filter covariance ---
verr2 = max(max(max(abs(oo2.Smoother.Variance - oo3.Smoother.Variance))));
if verr2 > 1e-10
    error('kalman_algo=4: smoother_redux does not replicate original filter covariance! max diff = %e', verr2)
end

// --- State uncertainty ---
verrS2 = max(max(max(abs(oo2.Smoother.State_uncertainty - oo3.Smoother.State_uncertainty))));
if verrS2 > 1e-10
    error('kalman_algo=4: smoother_redux does not replicate original state uncertainty! max diff = %e', verrS2)
end

fprintf('All smoother_redux (kalman_algo=4) diffuse consistency checks passed.\n');

// =========================================================================
// Cross-filter checks: multivariate (oo0) vs univariate (oo2) baseline
// =========================================================================

// --- Smoothed shocks ---
for k=1:M_.exo_nbr
    mserr_cross(k) = max(abs(oo0.SmoothedShocks.(M_.exo_names{k}) - oo2.SmoothedShocks.(M_.exo_names{k})));
end
if max(mserr_cross) > 1e-12
    error('Cross-filter: multivariate vs univariate do not agree on smoothed shocks! max diff = %e', max(mserr_cross))
end

// --- Smoothed variables ---
for k=1:M_.endo_nbr
    merr_cross(k) = max(abs(oo0.SmoothedVariables.(M_.endo_names{k}) - oo2.SmoothedVariables.(M_.endo_names{k})));
end
if max(merr_cross) > 1e-12
    error('Cross-filter: multivariate vs univariate do not agree on smoothed variables! max diff = %e', max(merr_cross))
end

// --- Updated variables ---
for k=1:M_.endo_nbr
    merrU_cross(k) = max(abs(oo0.UpdatedVariables.(M_.endo_names{k}) - oo2.UpdatedVariables.(M_.endo_names{k})));
end
if max(merrU_cross) > 1e-12
    error('Cross-filter: multivariate vs univariate do not agree on updated variables! max diff = %e', max(merrU_cross))
end

// --- Filtered variables ---
for k=1:M_.endo_nbr
    merrF_cross(k) = max(abs(oo0.FilteredVariables.(M_.endo_names{k}) - oo2.FilteredVariables.(M_.endo_names{k})));
end
if max(merrF_cross) > 1e-12
    error('Cross-filter: multivariate vs univariate do not agree on filtered variables! max diff = %e', max(merrF_cross))
end

// --- Filter covariance (P_{t|t-1}) ---
verr_cross = max(max(max(abs(oo0.Smoother.Variance - oo2.Smoother.Variance))));
if verr_cross > 1e-10
    error('Cross-filter: multivariate vs univariate do not agree on filter covariance! max diff = %e', verr_cross)
end

// --- State uncertainty (V_{t|T}) ---
verrS_cross = max(max(max(abs(oo0.Smoother.State_uncertainty - oo2.Smoother.State_uncertainty))));
if verrS_cross > 1e-10
    error('Cross-filter: multivariate vs univariate do not agree on state uncertainty! max diff = %e', verrS_cross)
end

fprintf('All cross-filter (multivariate vs univariate) baseline consistency checks passed.\n');

// =========================================================================
// Cross-filter checks: multivariate (oo1) vs univariate (oo3) smoother_redux
// =========================================================================

// --- Smoothed shocks ---
for k=1:M_.exo_nbr
    mserr_cross_sr(k) = max(abs(oo1.SmoothedShocks.(M_.exo_names{k}) - oo3.SmoothedShocks.(M_.exo_names{k})));
end
if max(mserr_cross_sr) > 1e-12
    error('Cross-filter smoother_redux: multivariate vs univariate do not agree on smoothed shocks! max diff = %e', max(mserr_cross_sr))
end

// --- Smoothed variables ---
for k=1:M_.endo_nbr
    merr_cross_sr(k) = max(abs(oo1.SmoothedVariables.(M_.endo_names{k}) - oo3.SmoothedVariables.(M_.endo_names{k})));
end
if max(merr_cross_sr) > 1e-12
    error('Cross-filter smoother_redux: multivariate vs univariate do not agree on smoothed variables! max diff = %e', max(merr_cross_sr))
end

// --- Updated variables ---
for k=1:M_.endo_nbr
    merrU_cross_sr(k) = max(abs(oo1.UpdatedVariables.(M_.endo_names{k}) - oo3.UpdatedVariables.(M_.endo_names{k})));
end
if max(merrU_cross_sr) > 1e-12
    error('Cross-filter smoother_redux: multivariate vs univariate do not agree on updated variables! max diff = %e', max(merrU_cross_sr))
end

// --- Filtered variables ---
for k=1:M_.endo_nbr
    merrF_cross_sr(k) = max(abs(oo1.FilteredVariables.(M_.endo_names{k}) - oo3.FilteredVariables.(M_.endo_names{k})));
end
if max(merrF_cross_sr) > 1e-12
    error('Cross-filter smoother_redux: multivariate vs univariate do not agree on filtered variables! max diff = %e', max(merrF_cross_sr))
end

// --- Filter covariance ---
verr_cross_sr = max(max(max(abs(oo1.Smoother.Variance - oo3.Smoother.Variance))));
if verr_cross_sr > 1e-10
    error('Cross-filter smoother_redux: multivariate vs univariate do not agree on filter covariance! max diff = %e', verr_cross_sr)
end

// --- State uncertainty ---
verrS_cross_sr = max(max(max(abs(oo1.Smoother.State_uncertainty - oo3.Smoother.State_uncertainty))));
if verrS_cross_sr > 1e-10
    error('Cross-filter smoother_redux: multivariate vs univariate do not agree on state uncertainty! max diff = %e', verrS_cross_sr)
end

fprintf('All cross-filter (multivariate vs univariate) smoother_redux consistency checks passed.\n');
