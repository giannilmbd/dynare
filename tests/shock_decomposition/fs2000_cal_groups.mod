@#include "fs2000.inc"

varobs gy_obs gp_obs;

shock_groups;
supply = e_a ;
demand = e_m ;
end;

mymap = [144/255 212/255 164/255;1 128/255 0; 51/255 51/255 1];

shock_decomposition(datafile=fsdat_simul,parameter_set=calibration,use_shock_groups,colormap=mymap) gp_obs gy_obs;
