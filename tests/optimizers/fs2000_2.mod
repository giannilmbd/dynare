@#include "fs2000.common.inc"

profile clear; profile on;
estimation(mode_compute=2,silent_optimizer,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0,
optim=(
'MaxIter',20000,
'TolFun',1e-4,
'TolX',1e-4)
);
@#include "optimizer_function_count.inc"