@#include "fs2000.common.inc"

profile clear; profile on;
estimation(mode_compute=5, silent_optimizer, order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0);
@#include "optimizer_function_count.inc"

profile clear; profile on;
estimation(mode_compute=5, silent_optimizer, order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0, analytic_derivation);
@#include "optimizer_function_count.inc"
options_.analytic_derivation = 0; % reset to numerical derivatives

% note that with optim=('Hessian',2) dynare_estimation_1 calls mr_hessian.m instead of hessian.m, function calls in mr_hessian are however adaptive, so we skip profiling this option
estimation(mode_compute=5, silent_optimizer, order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0, optim=('Hessian',2));

profile clear; profile on;
estimation(mode_compute=5, silent_optimizer, order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0, optim=('Hessian',1,'robust',true));
@#include "optimizer_function_count.inc"