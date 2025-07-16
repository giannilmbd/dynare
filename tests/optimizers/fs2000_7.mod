@#include "fs2000.common.inc"

if exist('fminsearch','file')
  profile clear; profile on;
  estimation(mode_compute=7,silent_optimizer,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0);
  @#include "optimizer_function_count.inc"
end
