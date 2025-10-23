/*
 * Copyright © 2019-2025 Dynare Team
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

#ifndef DYNAMIC_MODEL_CALLER_HH
#define DYNAMIC_MODEL_CALLER_HH

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <type_traits>
#include <vector>

#include <dynmex.h>

#if defined(_WIN32) || defined(__CYGWIN32__)
# ifndef NOMINMAX
#  define NOMINMAX // Do not define "min" and "max" macros
# endif
# include <windows.h>
#else
# include <dlfcn.h> // unix/linux DLL (.so) handling routines
#endif

// Base class for all variants of dynamic model callers
class DynamicModelCaller
{
public:
  const bool linear;
  bool compute_jacobian; // Not constant, because will be changed from true to false for linear
                         // models after first evaluation

  // Used to store error messages (as exceptions cannot cross the OpenMP boundary)
  static std::string error_msg;
  static std::string error_id;
  static std::mutex error_mtx; // Guard for access in OpenMP context

  DynamicModelCaller(bool linear_arg, bool compute_jacobian_arg) :
      linear {linear_arg}, compute_jacobian {compute_jacobian_arg}
  {
  }
  virtual ~DynamicModelCaller() = default;
  [[nodiscard]] virtual double* y() = 0;
  [[nodiscard]] virtual double* x() = 0;
  [[nodiscard]] virtual double* T() = 0;
  /* Copy a column of the Jacobian to dest.
     Only copies non-zero elements, according to g1_sparse_{rowval,colval,colptr}. */
  virtual void copy_jacobian_column(mwIndex col, double* dest) const = 0;
  virtual void eval(double* resid) = 0;
  static void setErrMsg(std::string msg);
  static void setMException(const mxArray* exception);
};

// Base class for DLL-based dynamic model callers
class DynamicModelDllCaller : public DynamicModelCaller
{
protected:
  const double *params, *steady_state;
  std::vector<double> tt, y_p, x_p, jacobian_p;
  const int32_T* g1_sparse_colptr;

public:
  DynamicModelDllCaller(size_t ntt, mwIndex ny, mwIndex nx, const double* params_arg,
                        const double* steady_state_arg, const int32_T* g1_sparse_colptr_arg,
                        bool linear_arg, bool compute_jacobian_arg);
  [[nodiscard]] double*
  y() override
  {
    return y_p.data();
  }
  [[nodiscard]] double*
  x() override
  {
    return x_p.data();
  }
  [[nodiscard]] double*
  T() override
  {
    return tt.data();
  }
  void copy_jacobian_column(mwIndex col, double* dest) const override;
};

// Base class for MATLAB-based dynamic model callers
class DynamicModelMatlabCaller : public DynamicModelCaller
{
protected:
  std::string basename;
  mxArray *y_mx, *x_mx, *jacobian_mx, *params_mx, *steady_state_mx, *g1_sparse_rowval_mx,
      *g1_sparse_colval_mx, *g1_sparse_colptr_mx;

  /* Given a complex matrix (of double floats), returns a sparse matrix of same size.
     Real elements of the original matrix are copied as-is to the new one.
     Complex elements are replaced by NaNs.
     Destroys the original matrix.
     There are two versions, one for dense matrices, another for sparse
     matrices. */
  template<bool sparse>
  static mxArray* cmplxToReal(mxArray* cmplx_mx);

public:
  DynamicModelMatlabCaller(std::string basename_arg, mwIndex ny, mwIndex nx,
                           const mxArray* params_mx_arg, const mxArray* steady_state_mx_arg,
                           const mxArray* g1_sparse_rowval_mx_arg,
                           const mxArray* g1_sparse_colval_mx_arg,
                           const mxArray* g1_sparse_colptr_mx_arg, bool linear_arg,
                           bool compute_jacobian_arg);
  ~DynamicModelMatlabCaller() override;
  [[nodiscard]] double*
  y() override
  {
    return mxGetPr(y_mx);
  }
  [[nodiscard]] double*
  x() override
  {
    return mxGetPr(x_mx);
  }
  void copy_jacobian_column(mwIndex col, double* dest) const override;
};

template<bool sparse>
mxArray*
DynamicModelMatlabCaller::cmplxToReal(mxArray* cmplx_mx)
{
  mxArray* real_mx {
      sparse ? mxCreateSparse(mxGetM(cmplx_mx), mxGetN(cmplx_mx), mxGetNzmax(cmplx_mx), mxREAL)
             : mxCreateDoubleMatrix(mxGetM(cmplx_mx), mxGetN(cmplx_mx), mxREAL)};

  if constexpr (sparse)
    {
      std::ranges::copy_n(mxGetIr(cmplx_mx), mxGetNzmax(cmplx_mx), mxGetIr(real_mx));
      std::ranges::copy_n(mxGetJc(cmplx_mx), mxGetN(cmplx_mx) + 1, mxGetJc(real_mx));
    }

  mxComplexDouble* cmplx {mxGetComplexDoubles(cmplx_mx)};
  double* real {mxGetPr(real_mx)};
  for (std::conditional_t<sparse, mwSize, size_t> i {0};
       i <
       [&] {
         if constexpr (sparse)
           return mxGetNzmax(cmplx_mx);
         else
           return mxGetNumberOfElements(cmplx_mx);
       }(); // Use a lambda instead of the ternary operator to have the right type (there is no
            // constexpr ternary operator)
       i++)
    if (cmplx[i].imag == 0.0)
      real[i] = cmplx[i].real;
    else
      real[i] = std::numeric_limits<double>::quiet_NaN();

  mxDestroyArray(cmplx_mx);
  return real_mx;
}

class DynamicModelNoblockDllCaller : public DynamicModelDllCaller
{
private:
#if defined(_WIN32) || defined(__CYGWIN32__)
  static HINSTANCE resid_mex, g1_mex;
#else
  static void *resid_mex, *g1_mex;
#endif
  using dynamic_tt_fct = void (*)(const double* y, const double* x, const double* params,
                                  const double* steady_state, double* T);
  using dynamic_fct = void (*)(const double* y, const double* x, const double* params,
                               const double* steady_state, const double* T, double* value);
  static dynamic_tt_fct residual_tt_fct, g1_tt_fct;
  static dynamic_fct residual_fct, g1_fct;

public:
  DynamicModelNoblockDllCaller(size_t ntt, mwIndex ny, mwIndex nx, const double* params_arg,
                               const double* steady_state_arg, const int32_T* g1_sparse_colptr_arg,
                               bool linear_arg, bool compute_jacobian_arg);
  void eval(double* resid) override;
  static void load_dll(const std::string& basename);
  static void unload_dll();
};

class DynamicModelNoblockMatlabCaller : public DynamicModelMatlabCaller
{
public:
  DynamicModelNoblockMatlabCaller(std::string basename_arg, mwIndex ny, mwIndex nx,
                                  const mxArray* params_mx_arg, const mxArray* steady_state_mx_arg,
                                  const mxArray* g1_sparse_rowval_mx_arg,
                                  const mxArray* g1_sparse_colval_mx_arg,
                                  const mxArray* g1_sparse_colptr_mx_arg, bool linear_arg,
                                  bool compute_jacobian_arg);
  [[nodiscard]] double*
  T() override
  {
    std::cerr << "This should not happen" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  void eval(double* resid) override;
};

class DynamicModelBlockDllCaller : public DynamicModelDllCaller
{
private:
#if defined(_WIN32) || defined(__CYGWIN32__)
  static HINSTANCE mex;
#else
  static void* mex;
#endif
  using dynamic_resid_fct = void (*)(double* y, const double* x, const double* params,
                                     const double* steady_state, double* T, double* residual);
  using dynamic_g1_fct = void (*)(const double* y, const double* x, const double* params,
                                  const double* steady_state, double* T, double* g1_v);
  static dynamic_resid_fct residual_fct;
  static dynamic_g1_fct g1_fct;

public:
  DynamicModelBlockDllCaller(size_t ntt, size_t mfs, mwIndex ny, mwIndex nx,
                             const double* params_arg, const double* steady_state_arg,
                             const int32_T* g1_sparse_colptr_arg, bool linear_arg,
                             bool compute_jacobian_arg);
  void eval(double* resid) override;
  static void load_dll(const std::string& basename, int block_num);
  static void unload_dll();
};

class DynamicModelBlockMatlabCaller : public DynamicModelMatlabCaller
{
private:
  int block_num;
  mxArray* T_mx;

public:
  DynamicModelBlockMatlabCaller(std::string basename_arg, int block_num_arg, mwIndex ntt,
                                mwIndex ny, mwIndex nx, const mxArray* params_mx_arg,
                                const mxArray* steady_state_mx_arg,
                                const mxArray* g1_sparse_rowval_mx_arg,
                                const mxArray* g1_sparse_colval_mx_arg,
                                const mxArray* g1_sparse_colptr_mx_arg, bool linear_arg,
                                bool compute_jacobian_arg);
  ~DynamicModelBlockMatlabCaller() override;
  [[nodiscard]] double*
  T() override
  {
    return mxGetPr(T_mx);
  }
  void eval(double* resid) override;
};

#endif
