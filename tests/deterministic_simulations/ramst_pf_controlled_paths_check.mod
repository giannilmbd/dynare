/* Checks that exogenous computed by ramst_pf_controlled_paths{,lbj,ssa7}.mod actually give back the
   expected endogenous variables */

@#define shocks = false
@#include "ramst_pf_controlled_paths_common.inc"


// stack_solve_algo=0

S0 = load('ramst_pf_controlled_paths/Output/ramst_pf_controlled_paths_results.mat');

perfect_foresight_setup(periods=30);

oo_.exo_simul = S0.oo_.exo_simul;

perfect_foresight_solver;

if max(max(abs(oo_.endo_simul - S0.oo_.endo_simul))) > options_.dynatol.x
    error('perfect_foresight_controlled_paths gives wrong result with stack_solve_algo=0')
end


// stack_solve_algo=1

S1 = load('ramst_pf_controlled_paths_lbj/Output/ramst_pf_controlled_paths_lbj_results.mat');

perfect_foresight_setup(periods=30);

oo_.exo_simul = S1.oo_.exo_simul;

perfect_foresight_solver;

if max(max(abs(oo_.endo_simul - S1.oo_.endo_simul))) > options_.dynatol.x
    error('perfect_foresight_controlled_paths gives wrong result with stack_solve_algo=1')
end


// stack_solve_algo=7

S7 = load('ramst_pf_controlled_paths_ssa7/Output/ramst_pf_controlled_paths_ssa7_results.mat');

perfect_foresight_setup(periods=30);

oo_.exo_simul = S7.oo_.exo_simul;

perfect_foresight_solver;

if max(max(abs(oo_.endo_simul - S7.oo_.endo_simul))) > options_.dynatol.x
    error('perfect_foresight_controlled_paths gives wrong result with stack_solve_algo=7')
end
