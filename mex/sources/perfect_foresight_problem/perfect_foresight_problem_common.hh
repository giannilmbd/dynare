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

#include <dynmex.h>
#include <utility>

/* In this file, njcol is the number of columns of a single-period Jacobian: ny in the no-block
   case, and mfs in the block case */

inline std::tuple<mxArray*, double*, mwIndex*, mwIndex*>
init_stacked_jacobian(mwIndex periods, mwIndex njcol, const int32_T* g1_sparse_rowval,
                      const int32_T* g1_sparse_colptr)
{
  /* Number of non-zero values in the stacked Jacobian.
     Contemporaneous derivatives appear at all periods, while lag or lead
     derivatives appear at all periods except one. */
  mwSize nzmax {static_cast<mwSize>((g1_sparse_colptr[3 * njcol] - 1) * (periods - 1)
                                    + (g1_sparse_colptr[2 * njcol] - g1_sparse_colptr[njcol]))};

  mxArray* stacked_jacobian_mx = mxCreateSparse(periods * njcol, periods * njcol, nzmax, mxREAL);
  double* stacked_jacobian = mxGetDoubles(stacked_jacobian_mx);
  mwIndex* ir = mxGetIr(stacked_jacobian_mx);
  mwIndex* jc = mxGetJc(stacked_jacobian_mx);

  /* Create the index vectors (IR, JC) of the sparse stacked jacobian. This
     makes parallelization across periods possible when evaluating the model,
     since all indices are known ex ante. */
  mwIndex k = 0;
  jc[0] = 0;
  for (mwIndex T = 0; T < periods; T++)
    for (mwIndex j {0}; j < njcol; j++) // Column within the period (i.e. variable)
      {
        if (T > 0)
          for (int32_T idx {g1_sparse_colptr[2 * njcol + j] - 1};
               idx < g1_sparse_colptr[2 * njcol + j + 1] - 1; idx++)
            ir[k++]
                = (T - 1) * njcol + g1_sparse_rowval[idx] - 1; // Derivatives w.r.t. y_{t+1} in T-1
        for (int32_T idx {g1_sparse_colptr[njcol + j] - 1};
             idx < g1_sparse_colptr[njcol + j + 1] - 1; idx++)
          ir[k++] = T * njcol + g1_sparse_rowval[idx] - 1; // Derivatives w.r.t. y_t in T
        if (T < periods - 1)
          for (int32_T idx {g1_sparse_colptr[j] - 1}; idx < g1_sparse_colptr[j + 1] - 1; idx++)
            ir[k++]
                = (T + 1) * njcol + g1_sparse_rowval[idx] - 1; // Derivatives w.r.t. y_{t-1} in T+1
        jc[T * njcol + j + 1] = k;
      }

  return {stacked_jacobian_mx, stacked_jacobian, ir, jc};
}

/* Fills vector of dynamic endogenous (y3n), given a period number (T), and paths for endogenous
   variables (y, y0, yT) */
inline void
fill_dynamic_endogenous(mwIndex T, mwIndex periods, mwIndex ny, const double* y, const double* y0,
                        const double* yT, double* y3n)
{
  if (T > 0 && T < periods - 1)
    std::ranges::copy_n(y + (T - 1) * ny, 3 * ny, y3n);
  else if (T > 0) // Last simulation period
    {
      std::ranges::copy_n(y + (T - 1) * ny, 2 * ny, y3n);
      std::ranges::copy_n(yT, ny, y3n + 2 * ny);
    }
  else if (T < periods - 1) // First simulation period
    {
      std::ranges::copy_n(y0, ny, y3n);
      std::ranges::copy_n(y + T * ny, 2 * ny, y3n + ny);
    }
  else // Special case: periods=1 (and so T=0)
    {
      std::ranges::copy_n(y0, ny, y3n);
      std::ranges::copy_n(y, ny, y3n + ny);
      std::ranges::copy_n(yT, ny, y3n + 2 * ny);
    }
}

inline void
fill_stacked_jacobian(mwIndex T, mwIndex periods, mwIndex njcol, double* stacked_jacobian,
                      mwIndex* ir, mwIndex* jc, auto& m)
{
  // Fill the stacked jacobian
  for (mwIndex col {T > 1 ? (T - 1) * njcol : 0}; // We can't use std::max() here, because
                                                  // mwIndex is unsigned under MATLAB
       col < std::min(periods * njcol, (T + 2) * njcol); col++)
    {
      mwIndex k = jc[col];
      while (k < jc[col + 1])
        {
          if (ir[k] >= (T + 1) * njcol)
            break; // Nothing to copy for this column
          if (ir[k] >= T * njcol)
            {
              /* Within the current column, this is the first line of
                 the stacked Jacobian that contains elements from the
                 (small) Jacobian just computed, so copy the whole
                 column of the latter to the former. */
              m->copy_jacobian_column(col - (T - 1) * njcol, stacked_jacobian + k);
              break;
            }
          k++;
        }
    }
}

/* The stacked jacobian so far constructed can still reference some zero
   elements, because some expressions that are symbolically different from
   zero may still evaluate to zero for some endogenous/parameter values. The
   following code further compresses the Jacobian by removing the references
   to those extra zeros. This makes a significant speed difference when
   inversing the Jacobian for some large models. */
inline void
compress_stacked_jacobian(mwIndex periods, mwIndex njcol, double* stacked_jacobian, mwIndex* ir,
                          mwIndex* jc)
{
  mwIndex k_orig = 0, k_new = 0;
  for (mwIndex col = 0; col < periods * njcol; col++)
    {
      while (k_orig < jc[col + 1])
        {
          if (stacked_jacobian[k_orig] != 0.0)
            {
              if (k_new != k_orig)
                {
                  stacked_jacobian[k_new] = stacked_jacobian[k_orig];
                  ir[k_new] = ir[k_orig];
                }
              k_new++;
            }
          k_orig++;
        }
      jc[col + 1] = k_new;
    }
}
