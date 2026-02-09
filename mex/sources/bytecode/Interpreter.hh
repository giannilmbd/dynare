/*
 * Copyright © 2007-2026 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INTERPRETER_HH
#define INTERPRETER_HH

#include <cstddef>
#include <fstream>
#include <map>
#include <stack>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "dynmex.h"
#include "dynumfpack.h"

#include "ErrorHandling.hh"
#include "Evaluate.hh"
#include "Mem_Mngr.hh"

using namespace std;

struct t_save_op_s
{
  short int lag, operat;
  int first, second;
};

constexpr int IFLD = 0, IFDIV = 1, IFLESS = 2, IFSUB = 3, IFLDZ = 4, IFMUL = 5, IFSTP = 6,
              IFADD = 7;
constexpr double eps = 1e-15, very_big = 1e24;
constexpr int alt_symbolic_count_max = 1;
constexpr double mem_increasing_factor = 1.1;

constexpr int NO_ERROR_ON_EXIT {0}, ERROR_ON_EXIT {1};

// Stores various options for iterative solvers (GMRES, BiCGStab), from options_.simul
struct iter_solver_opts_t
{
  const char* preconditioner;
  mxArray* iter_tol;
  mxArray* iter_maxit;
  mxArray* gmres_restart;
  mxArray* incomplete_lu;

  void set_static_values(const mxArray* options_);
  void set_dynamic_values(const mxArray* options_);
};

class Interpreter
{
private:
  double g0, gp0, glambda2;
  int try_at_iteration;

  void *Symbolic {nullptr}, *Numeric {nullptr};

  const BasicSymbolTable& symbol_table;
  const bool steady_state; // Whether this is a static or dynamic model

  // Whether to use the block-decomposed version of the bytecode file
  bool block_decomposed;

  Evaluate& evaluator;

  fstream SaveCode;

  Mem_Mngr mem_mngr;
  vector<int> u_liste;
  int *NbNZRow, *NbNZCol;
  NonZeroElem **FNZE_R, **FNZE_C;
  int u_count_init;

  int *pivot, *pivotk, *pivot_save;
  double *pivotv, *pivotva;
  int* b_perm;
  bool* line_done;
  bool symbolic, alt_symbolic;
  int alt_symbolic_count;
  double markowitz_c_s;
  double res1a;
  long int nop1;
  map<tuple<int, int, int>, int> IM_i;
  int u_count_alloc, u_count_alloc_save;
  double slowc, slowc_save, prev_slowc_save, markowitz_c;
  int* index_equa; // Actually unused
  int u_count, tbreak_g;
  int iter;
  int start_compare;
  int restart;
  iter_solver_opts_t iter_solver_opts;

  SuiteSparse_long *Ap_save, *Ai_save;
  double *Ax_save, *b_save;

  int stack_solve_algo, solve_algo;
  int minimal_solving_periods;
  int Per_u_, Per_y_;
  int maxit_;
  double* direction;
  double solve_tolf;
  // 1-norm error, square of 2-norm error, ∞-norm error
  double res1, res2, max_res;
  int max_res_idx;
  int* index_vara;

  double *y, *ya;
  int y_size;
  double* T;
  int nb_row_x;
  int y_kmin, y_kmax, periods;
  double *x, *params;
  double g1; // The derivative of a simple (single-equation) block in simulate mode
  double* u;
  double* steady_y;
  double *r, *res;
  vector<mxArray*> jacobian_block, jacobian_exo_block, jacobian_det_exo_block;
  mxArray* GlobalTemporaryTerms;
  int it_;
  map<int, double> TEF;
  map<pair<int, int>, double> TEFD;
  map<tuple<int, int, int>, double> TEFDD;

  // Information about the current block
  int block_num; // Index of the current block
  int size;      // Size of the current block
  BlockSimulationType type;
  bool is_linear;
  int u_count_int;
  vector<int> variables, equations;

  int verbosity; // Corresponds to options_.verbosity

  bool print; // Whether the “print” command is requested
  int col_x, col_y;
  vector<double> residual;

  void evaluate_over_periods(bool forward);
  void solve_simple_one_periods();
  void solve_simple_over_periods(bool forward);
  void compute_complete_2b();
  void evaluate_a_block(bool initialization, bool single_block, const string& bin_base_name);
  int simulate_a_block(const string& bin_base_name);
  void Simulate_Newton_Two_Boundaries(bool cvg);
  void Simulate_Newton_One_Boundary(bool forward);
  void fixe_u();
  void Read_SparseMatrix(const string& file_name, bool two_boundaries);
  void Singular_display();
  void End_Solver();
  void Init_Gaussian_Elimination();
  void Init_Matlab_Sparse_Two_Boundaries(const mxArray* A_m, const mxArray* b_m,
                                         const mxArray* x0_m) const;
  tuple<SuiteSparse_long*, SuiteSparse_long*, double*, double*>
  Init_UMFPACK_Sparse_Two_Boundaries(const mxArray* x0_m) const;
  bool Init_Matlab_Sparse_One_Boundary(const mxArray* A_m, const mxArray* b_m,
                                       const mxArray* x0_m) const;
  tuple<bool, SuiteSparse_long*, SuiteSparse_long*, double*, double*>
  Init_UMFPACK_Sparse_One_Boundary(const mxArray* x0_m) const;
  bool Simple_Init();
  void End_Gaussian_Elimination();
  tuple<bool, double, double, double, double> mnbrak(double& ax, double& bx);
  pair<bool, double> golden(double ax, double bx, double cx, double tol);
  void Solve_ByteCode_Symbolic_Sparse_GaussianElimination(bool symbolic);
  bool Solve_ByteCode_Sparse_GaussianElimination();
  void Solve_LU_UMFPack_Two_Boundaries(SuiteSparse_long* Ap, SuiteSparse_long* Ai, double* Ax,
                                       double* b);
  void Solve_LU_UMFPack_One_Boundary(SuiteSparse_long* Ap, SuiteSparse_long* Ai, double* Ax,
                                     double* b);

  void Solve_Matlab_GMRES(mxArray* A_m, mxArray* b_m, bool is_two_boundaries, mxArray* x0_m);
  void Solve_Matlab_BiCGStab(mxArray* A_m, mxArray* b_m, bool is_two_boundaries, mxArray* x0_m);
  void Check_and_Correct_Previous_Iteration();
  bool Simulate_One_Boundary();
  bool solve_linear(bool do_check_and_correct);
  void solve_non_linear();
  bool compare(int* save_op, int* save_opa, int* save_opaa, int beg_t, long nop4);
  void Insert(int r, int c, int u_index, int lag_index);
  void Delete(int r, int c);
  pair<int, NonZeroElem*> At_Row(int r) const;
  NonZeroElem* At_Pos(int r, int c) const;
  pair<int, NonZeroElem*> At_Col(int c) const;
  pair<int, NonZeroElem*> At_Col(int c, int lag) const;
  int NRow(int r) const;
  int NCol(int c) const;
  int Get_u();
  void Delete_u(int pos);
  void Clear_u();

  int complete(int beg_t);
  void bksub(int tbreak, int last_period);
  void simple_bksub();

  void compute_block_time(int my_Per_u_, bool evaluate, bool no_derivatives);
  bool compute_complete(bool no_derivatives);

  pair<bool, double> compute_complete(double lambda);

public:
  Interpreter(Evaluate& evaluator_arg, double* params_arg, double* y_arg, double* ya_arg,
              double* x_arg, double* steady_y_arg, double* direction_arg, int y_size_arg,
              int nb_row_x_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg, int maxit_arg_,
              double solve_tolf_arg, double markowitz_c_arg, int minimal_solving_periods_arg,
              int stack_solve_algo_arg, int solve_algo_arg, bool print_arg,
              const mxArray* GlobalTemporaryTerms_arg, bool steady_state_arg,
              bool block_decomposed_arg, int col_x_arg, int col_y_arg,
              const BasicSymbolTable& symbol_table_arg, int verbosity_arg,
              iter_solver_opts_t iter_solver_opts_arg);
  pair<bool, vector<int>> compute_blocks(const string& file_name, bool evaluate, int block);
  void Close_SaveCode();

  inline mxArray*
  get_jacob(int block_num) const
  {
    return jacobian_block[block_num];
  }
  inline mxArray*
  get_jacob_exo(int block_num) const
  {
    return jacobian_exo_block[block_num];
  }
  inline mxArray*
  get_jacob_exo_det(int block_num) const
  {
    return jacobian_det_exo_block[block_num];
  }
  inline vector<double>
  get_residual() const
  {
    return residual;
  }
  inline mxArray*
  get_Temporary_Terms() const
  {
    return GlobalTemporaryTerms;
  }
};

#endif
