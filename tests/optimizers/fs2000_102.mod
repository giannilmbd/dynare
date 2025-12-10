@#define SPEED_UP_TESTSUITE
@#include "fs2000.common.inc"

estimated_params_init;
stderr e_a, 0.015759023653034;
stderr e_m, 0.007074670436290;
alp, 0.360194076098251;
bet, 0.992147804788292;
gam, 0.003558883548504;
mst, 1.012615945228039;
rho, 0.314245018379123;
psi, 0.667051411177380;
del, 0.007873856590655;
end;

if ~isoctave() && exist('simulannealbnd','file')
  profile clear; profile on;
  estimation(mode_compute=102, silent_optimizer, order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0, mh_nblocks=2, mh_jscale=0.8,
             optim=('TolFun',1e-2)
  );
  @#include "optimizer_function_count.inc"
end
