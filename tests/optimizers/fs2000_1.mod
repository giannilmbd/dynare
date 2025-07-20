@#include "fs2000.common.inc"

if (isoctave && user_has_octave_forge_package('optim', '1.6')) || (~isoctave && user_has_matlab_license('optimization_toolbox'))
  % numerical derivatives
  profile clear; profile on;
  estimation(mode_compute=1, silent_optimizer, order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0);
  @#include "optimizer_function_count.inc"

  % analytical derivatives
  profile clear; profile on;
  estimation(mode_compute=1, silent_optimizer, order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0, analytic_derivation);
  @#include "optimizer_function_count.inc"
end

