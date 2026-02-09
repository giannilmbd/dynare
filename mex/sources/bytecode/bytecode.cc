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

#include <algorithm>
#include <cmath>
#include <cstring>
#include <type_traits>

#include "ErrorHandling.hh"
#include "Interpreter.hh"

string
Get_Argument(const mxArray* prhs)
{
  const mxArray* mxa = prhs;
  auto buflen = mwSize(mxGetM(mxa) * mxGetN(mxa) + 1);
  char* first_argument;
  first_argument = static_cast<char*>(mxCalloc(buflen, sizeof(char)));
  size_t status = mxGetString(mxa, first_argument, buflen);
  if (status != 0)
    mexWarnMsgTxt("Not enough space. The first argument is truncated.");
  string f(first_argument);
  mxFree(first_argument);
  return f;
}

string
deblank(string x)
{
  for (int i = 0; i < static_cast<int>(x.length()); i++)
    if (x[i] == ' ')
      x.erase(i--, 1);
  return x;
}

// NOLINTBEGIN(modernize-avoid-c-arrays)
void
Get_Arguments_and_global_variables(int nrhs, const mxArray* prhs[], double* yd[], size_t& row_y,
                                   size_t& col_y, double* xd[], size_t& row_x, size_t& col_x,
                                   double* params[], double* steady_yd[], size_t& steady_row_y,
                                   size_t& steady_col_y, int& periods, mxArray** block_structur,
                                   bool& steady_state, bool& block_decomposed, bool& evaluate,
                                   int& block, const mxArray** M_, const mxArray** options_,
                                   bool& print, const mxArray** GlobalTemporaryTerms)
// NOLINTEND(modernize-avoid-c-arrays)
{
  int count_array_argument {0};
  for (int i = 0; i < nrhs; i++)
    {
#ifdef DEBUG
      if (mxIsChar(prhs[i]))
        mexPrintf("Arg %d: %s\n", i, Get_Argument(prhs[i]).c_str());
#endif
      if (!mxIsChar(prhs[i]))
        {
          switch (count_array_argument)
            {
            case 0:
              *M_ = prhs[i];
              break;
            case 1:
              *options_ = prhs[i];
              break;
            case 2:
              *yd = mxGetDoubles(prhs[i]);
              row_y = mxGetM(prhs[i]);
              col_y = mxGetN(prhs[i]);
              break;
            case 3:
              *xd = mxGetDoubles(prhs[i]);
              row_x = mxGetM(prhs[i]);
              col_x = mxGetN(prhs[i]);
              break;
            case 4:
              *params = mxGetDoubles(prhs[i]);
              break;
            case 5:
              *steady_yd = mxGetDoubles(prhs[i]);
              steady_row_y = mxGetM(prhs[i]);
              steady_col_y = mxGetN(prhs[i]);
              break;
            case 6:
              periods = static_cast<int>(mxGetScalar(prhs[i]));
              break;
            case 7:
              *block_structur = mxDuplicateArray(prhs[i]);
              break;
            case 8:
              *GlobalTemporaryTerms = prhs[i];
              break;
            default:
              mexPrintf("Unknown argument count_array_argument=%d\n", count_array_argument);
              break;
            }
          count_array_argument++;
        }
      else if (Get_Argument(prhs[i]) == "static")
        steady_state = true;
      else if (Get_Argument(prhs[i]) == "dynamic")
        steady_state = false;
      else if (Get_Argument(prhs[i]) == "block_decomposed")
        block_decomposed = true;
      else if (Get_Argument(prhs[i]) == "evaluate")
        evaluate = true;
      else if (Get_Argument(prhs[i]) == "print")
        print = true;
      else if (Get_Argument(prhs[i]).substr(0, 6) == "block=")
        {
          try
            {
              block = stoi(Get_Argument(prhs[i]).substr(6)) - 1;
            }
          catch (...)
            {
              throw FatalException {"ERROR: incorrect syntax for the 'block=' option"};
            }
        }
      else
        throw FatalException {"In main, unknown argument : " + Get_Argument(prhs[i])};
    }
  if (steady_state)
    {
      if (count_array_argument < 5)
        throw FatalException {"In a static context, the following arguments have to be indicated: "
                              "M_, options_, y, x, params"};
      if (count_array_argument < 7)
        periods = 1;
    }
  else
    {
      if (count_array_argument < 7)
        throw FatalException {"In a dynamic context, the following arguments have to be indicated: "
                              "M_, options_, y, x, params, steady_state, periods"};
    }
}

/* The gateway routine */
void
mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
  const mxArray *M_, *options_;
  const mxArray* GlobalTemporaryTerms {nullptr};
  mxArray* block_structur = nullptr;
  size_t i, row_y = 0, col_y = 0, row_x = 0, col_x = 0;
  size_t steady_row_y, steady_col_y;
  int y_kmin = 0, y_kmax = 0;
  int periods {1};
  double* direction;
  bool steady_state = false;
  bool block_decomposed {false};
  bool evaluate = false;
  int block = -1;
  double* params = nullptr;
  double *yd = nullptr, *xd = nullptr;
  bool print = false; // Whether the “print” command is requested
  int verbosity {1};  // Corresponds to options_.verbosity
  double* steady_yd = nullptr;

#ifdef DEBUG
  mexPrintf("**************************************\n");
  mexPrintf("ENTERING BYTECODE: nargin=%d, nargout=%d\n", nrhs, nlhs);
#endif

  try
    {
      Get_Arguments_and_global_variables(nrhs, prhs, &yd, row_y, col_y, &xd, row_x, col_x, &params,
                                         &steady_yd, steady_row_y, steady_col_y, periods,
                                         &block_structur, steady_state, block_decomposed, evaluate,
                                         block, &M_, &options_, print, &GlobalTemporaryTerms);
    }
  catch (GeneralException& feh)
    {
      mexErrMsgTxt(feh.message.c_str());
    }
#ifdef DEBUG
  mexPrintf("**************************************\n");
#endif

  BasicSymbolTable symbol_table {M_};

  if (!steady_state)
    {
      int field = mxGetFieldNumber(M_, "maximum_lag");
      if (field >= 0)
        y_kmin = static_cast<int>(floor(*(mxGetDoubles(mxGetFieldByNumber(M_, 0, field)))));
      else
        mexErrMsgTxt("maximum_lag is not a field of M_");
      field = mxGetFieldNumber(M_, "maximum_lead");
      if (field >= 0)
        y_kmax = static_cast<int>(floor(*(mxGetDoubles(mxGetFieldByNumber(M_, 0, field)))));
      else
        mexErrMsgTxt("maximum_lead is not a field of M_");
    }

  int field = mxGetFieldNumber(options_, "verbosity");
  if (field >= 0)
    verbosity = static_cast<int>(mxGetScalar(mxGetFieldByNumber(options_, 0, field)));
  else
    mexErrMsgTxt("verbosity is not a field of options_");

  if (!steady_state)
    {
      field = mxGetFieldNumber(options_, "simul");
      if (field < 0)
        mexErrMsgTxt("simul is not a field of options_");
    }
  else
    {
      field = mxGetFieldNumber(options_, "steady");
      if (field < 0)
        mexErrMsgTxt("steady is not a field of options_");
    }
  mxArray* temporarystruct {mxGetFieldByNumber(options_, 0, field)};

  field = mxGetFieldNumber(temporarystruct, "maxit");
  if (field < 0)
    {
      if (!steady_state)
        mexErrMsgTxt("maxit is not a field of options_.simul");
      else
        mexErrMsgTxt("maxit is not a field of options_.steady");
    }
  int maxit_
      = static_cast<int>(floor(*mxGetDoubles(mxGetFieldByNumber(temporarystruct, 0, field))));

  field = mxGetFieldNumber(options_, "markowitz");
  if (field < 0)
    mexErrMsgTxt("markowitz is not a field of options_");
  auto markowitz_c = static_cast<double>(*mxGetDoubles(mxGetFieldByNumber(options_, 0, field)));
  field = mxGetFieldNumber(options_, "minimal_solving_periods");
  if (field < 0)
    mexErrMsgTxt("minimal_solving_periods is not a field of options_");
  int minimal_solving_periods
      = static_cast<int>(*mxGetDoubles(mxGetFieldByNumber(options_, 0, field)));
  field = mxGetFieldNumber(options_, "stack_solve_algo");
  if (field < 0)
    mexErrMsgTxt("stack_solve_algo is not a field of options_");
  int stack_solve_algo = static_cast<int>(*mxGetDoubles(mxGetFieldByNumber(options_, 0, field)));

  field = mxGetFieldNumber(options_, "solve_algo");
  if (field < 0)
    mexErrMsgTxt("solve_algo is not a field of options_");
  int solve_algo = static_cast<int>(*mxGetDoubles(mxGetFieldByNumber(options_, 0, field)));

  /* Solver tolerance with respect to the residual. Equals options_.solve_tolf
     in the static case, or options_.dynatol.f in the dynamic case */
  double solve_tolf {[options_, steady_state] {
    if (steady_state)
      {
        int field {mxGetFieldNumber(options_, "solve_tolf")};
        if (field < 0)
          mexErrMsgTxt("solve_tolf is not a field of options_");
        return *mxGetDoubles(mxGetFieldByNumber(options_, 0, field));
      }
    else
      {
        int field {mxGetFieldNumber(options_, "dynatol")};
        if (field < 0)
          mexErrMsgTxt("dynatol is not a field of options_");
        mxArray* dynatol {mxGetFieldByNumber(options_, 0, field)};
        field = mxGetFieldNumber(dynatol, "f");
        if (field < 0)
          mexErrMsgTxt("f is not a field of options_.dynatol");
        return *mxGetDoubles(mxGetFieldByNumber(dynatol, 0, field));
      }
  }()};

  iter_solver_opts_t iter_solver_opts;
  if (steady_state)
    iter_solver_opts.set_static_values(options_);
  else
    {
      if (stack_solve_algo == 2 || stack_solve_algo == 3)
        {
          // Check that a preconditioner other than 'incomplete_lu' has not been requested
          int field {mxGetFieldNumber(options_, "simul")};
          if (field < 0)
            mexErrMsgTxt("simul is not a field of options_");
          mxArray* simul {mxGetFieldByNumber(options_, 0, field)};
          field = mxGetFieldNumber(simul, "preconditioner");
          if (field < 0)
            mexErrMsgTxt("preconditioner is not a field of options_.simul");
          mxArray* preconditioner {mxGetFieldByNumber(simul, 0, field)};
          if (!mxIsChar(preconditioner))
            mexErrMsgTxt("options_.simul.preconditioner should be a character array");
          char* preconditioner_str {mxArrayToString(preconditioner)};
          if (std::strcmp(preconditioner_str, "incomplete_lu"))
            mexErrMsgTxt("The 'bytecode' option is not compatible with a preconditioner other than "
                         "'incomplete_lu'");
          mxFree(preconditioner_str);
        }

      iter_solver_opts.set_dynamic_values(options_);
    }

  field = mxGetFieldNumber(M_, "fname");
  if (field < 0)
    mexErrMsgTxt("fname is not a field of M_");
  mxArray* mxa {mxGetFieldByNumber(M_, 0, field)};

  size_t buflen = mxGetM(mxa) * mxGetN(mxa) + 1;
  char* fname = static_cast<char*>(mxCalloc(buflen + 1, sizeof(char)));
  size_t status = mxGetString(mxa, fname, static_cast<int>(buflen));
  fname[buflen] = ' ';
  if (status != 0)
    mexWarnMsgTxt("Not enough space. Filename is truncated.");
  string file_name = fname;
  mxFree(fname);

  if (stack_solve_algo == 7 && !steady_state && !print)
    mexErrMsgTxt("Bytecode: Can't use option stack_solve_algo=7");

  if ((stack_solve_algo == 1 || stack_solve_algo == 6) && !steady_state && !print && !evaluate)
    mexErrMsgTxt("Bytecode: Can't use option stack_solve_algo=1 or 6");

  if (steady_state && !evaluate && !print && (solve_algo < 5 || solve_algo > 8))
    mexErrMsgTxt(
        "Bytecode: solve_algo must be between 5 and 8 when using the internal steady state solver");

  size_t size_of_direction = col_y * row_y * sizeof(double);
  auto* y = static_cast<double*>(mxMalloc(size_of_direction));
  test_mxMalloc(y, __LINE__, __FILE__, __func__, size_of_direction);
  auto* ya = static_cast<double*>(mxMalloc(size_of_direction));
  test_mxMalloc(ya, __LINE__, __FILE__, __func__, size_of_direction);
  direction = static_cast<double*>(mxMalloc(size_of_direction));
  test_mxMalloc(direction, __LINE__, __FILE__, __func__, size_of_direction);
  auto* x = static_cast<double*>(mxMalloc(col_x * row_x * sizeof(double)));
  test_mxMalloc(x, __LINE__, __FILE__, __func__, col_x * row_x * sizeof(double));

  ranges::fill_n(direction, row_y * col_y, 0);
  ranges::copy_n(xd, row_x * col_x, x);
  ranges::copy_n(yd, row_y * col_y, y);
  ranges::copy_n(yd, row_y * col_y, ya);

  const filesystem::path codfile {file_name + "/model/bytecode/"
                                  + (block_decomposed ? "block/" : "")
                                  + (steady_state ? "static" : "dynamic") + ".cod"};
  Evaluate evaluator {codfile, steady_state, symbol_table};

  Interpreter interprete {evaluator,
                          params,
                          y,
                          ya,
                          x,
                          steady_yd,
                          direction,
                          static_cast<int>(row_y),
                          static_cast<int>(row_x),
                          periods,
                          y_kmin,
                          y_kmax,
                          maxit_,
                          solve_tolf,
                          markowitz_c,
                          minimal_solving_periods,
                          stack_solve_algo,
                          solve_algo,
                          print,
                          GlobalTemporaryTerms,
                          steady_state,
                          block_decomposed,
                          static_cast<int>(col_x),
                          static_cast<int>(col_y),
                          symbol_table,
                          verbosity,
                          iter_solver_opts};
  bool r;
  vector<int> blocks;

  try
    {
      tie(r, blocks) = interprete.compute_blocks(file_name, evaluate, block);
    }
  catch (GeneralException& feh)
    {
      // Release the lock on dynamic.bin for MATLAB+Windows, see #1815
      interprete.Close_SaveCode();
      mexErrMsgTxt(feh.message.c_str());
    }

  if (nlhs > 0)
    {
      if (evaluate)
        {
          vector<double> residual = interprete.get_residual();
          plhs[0] = mxCreateDoubleMatrix(residual.size() / periods, periods, mxREAL);
          std::ranges::copy(residual, mxGetDoubles(plhs[0]));
        }
      else
        {
          plhs[0] = mxCreateDoubleMatrix(row_y, col_y, mxREAL);
          std::ranges::copy_n(y, row_y * col_y, mxGetDoubles(plhs[0]));
        }
      if (nlhs > 1)
        {
          if (evaluate)
            {
              int jacob_field_number = 0, jacob_exo_field_number = 0,
                  jacob_exo_det_field_number = 0;
              bool dont_store_a_structure {false};
              if (!block_structur)
                {
                  const char* field_names[] = {"g1", "g1_x", "g1_xd"};
                  jacob_field_number = 0;
                  jacob_exo_field_number = 1;
                  jacob_exo_det_field_number = 2;
                  mwSize dims[1] = {static_cast<mwSize>(blocks.size())};
                  plhs[1] = mxCreateStructArray(1, dims, std::extent_v<decltype(field_names)>,
                                                field_names);
                }
              else if (!mxIsStruct(block_structur))
                {
                  plhs[1] = interprete.get_jacob(blocks[0]);
                  dont_store_a_structure = true;
                }
              else
                {
                  plhs[1] = block_structur;
                  jacob_field_number = mxAddField(plhs[1], "g1");
                  if (jacob_field_number == -1)
                    mexErrMsgTxt("Fatal error in bytecode: in main, cannot add extra field jacob "
                                 "to the structArray");
                  jacob_exo_field_number = mxAddField(plhs[1], "g1_x");
                  if (jacob_exo_field_number == -1)
                    mexErrMsgTxt("Fatal error in bytecode: in main, cannot add extra field "
                                 "jacob_exo to the structArray");
                  jacob_exo_det_field_number = mxAddField(plhs[1], "g1_xd");
                  if (jacob_exo_det_field_number == -1)
                    mexErrMsgTxt("Fatal error in bytecode: in main, cannot add extra field "
                                 "jacob_exo_det to the structArray");
                }
              if (!dont_store_a_structure)
                for (size_t i {0}; i < blocks.size(); i++)
                  {
                    mxSetFieldByNumber(plhs[1], i, jacob_field_number,
                                       interprete.get_jacob(blocks[i]));
                    if (!steady_state)
                      {
                        mxSetFieldByNumber(plhs[1], i, jacob_exo_field_number,
                                           interprete.get_jacob_exo(blocks[i]));
                        mxSetFieldByNumber(plhs[1], i, jacob_exo_det_field_number,
                                           interprete.get_jacob_exo_det(blocks[i]));
                      }
                  }
            }
          else
            {
              plhs[1] = mxCreateDoubleMatrix(row_x, col_x, mxREAL);
              double* pind = mxGetDoubles(plhs[1]);
              for (i = 0; i < row_x * col_x; i++)
                pind[i] = x[i];
            }
          if (nlhs > 2)
            {
              plhs[2] = mxCreateDoubleMatrix(row_y, col_y, mxREAL);
              double* pind = mxGetDoubles(plhs[2]);
              for (i = 0; i < row_y * col_y; i++)
                pind[i] = y[i];
              if (nlhs > 3)
                plhs[3] = interprete.get_Temporary_Terms();
            }
        }
    }
  if (x)
    mxFree(x);
  if (y)
    mxFree(y);
  if (ya)
    mxFree(ya);
  if (direction)
    mxFree(direction);
}
