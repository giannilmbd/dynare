! Copyright © 2025 Dynare Team
!
! This file is part of Dynare.
!
! Dynare is free software: you can redistribute it and/or modify it under the terms of the
! GNU General Public License as published by the Free Software Foundation, either version 3 of
! the License, or (at your option) any later version.
!
! Dynare is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with Dynare. If not,
! see <https://www.gnu.org/licenses/>.
! This routine provides fast, vectorized bracket search and linear weight computation
! for 1D interpolation over a monotonic grid. It wraps the Fortran subroutine
! `bracket_linear_weight` and exposes it as a MATLAB MEX function.
!
! SYNTAX (in MATLAB)
!   [xqi, xqpi] = mexFunction(x, xq)
!
! INPUTS
!   x    [1×n or n×1 double] : Monotonically increasing vector of grid points.
!   xq   [1×nq or nq×1 double] : Vector of query points.
!
! OUTPUTS
!   xqi   [nq×1 int32]  : For each query point xq(i), the index ilow such that
!                         x(ilow) ≤ xq(i) ≤ x(ilow+1). Clamped if out-of-bounds.
!
!   xqpi  [nq×1 double] : Linear interpolation weights for the lower bracket point:
!                         xqpi(i) = (x(ilow+1) - xq(i)) / (x(ilow+1) - x(ilow))
!
! NOTES
! - Input vectors must be real, dense, and non-complex.
! - Extrapolation is clamped to the nearest interval: [x(1), x(2)] or [x(n-1), x(n)].
subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use iso_fortran_env, only: real64, int32
   use matlab_mex
   use interpolation
   implicit none (type, external)

   ! MATLAB MEX API
   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs
   type(c_ptr) :: x_mx, xq_mx ! Input arguments

   ! Fortran variables
   integer(int32) :: n, nq
   real(real64), pointer, contiguous :: x(:), xq(:), xqpi(:)
   integer(int32), pointer, contiguous :: xqi(:)

   ! Validate number of inputs/outputs
   if (nrhs /= 2_c_int .or. nlhs /= 2_c_int) then
      call mexErrMsgTxt("Rouwenhorst routine: incorrect number of inputs and/or outputs")
   end if

   ! Get input pointers
   x_mx = prhs(1)
   xq_mx = prhs(2)

   ! Validate input types
   if (.not. (mxIsDouble(x_mx) .and. (mxGetM(x_mx) == 1 .or. mxGetN(x_mx) == 1)) &
        .or. mxIsComplex(x_mx) .or. mxIsSparse(x_mx)) then
      call mexErrMsgTxt("1st argument (x) should be a real dense vector")
   end if
   if (.not. (mxIsDouble(xq_mx) .and. (mxGetM(xq_mx) == 1 .or. mxGetN(xq_mx) == 1)) &
        .or. mxIsComplex(xq_mx) .or. mxIsSparse(xq_mx)) then
      call mexErrMsgTxt("2nd argument (xq) should be a real dense vector")
   end if

   ! Convert MATLAB inputs to Fortran variables
   x => mxGetDoubles(x_mx)
   xq => mxGetDoubles(xq_mx)
   n = int(mxGetNumberOfElements(x_mx), int32)
   nq = int(mxGetNumberOfElements(xq_mx), int32)

   ! Allocate MATLAB output variables
   plhs(1) = mxCreateNumericMatrix(int(nq, mwSize), 1_mwSize, mxINT32_CLASS, mxREAL)
   plhs(2) = mxCreateDoubleMatrix(int(nq, mwSize), 1_mwSize, mxREAL)

   ! Copy results from Fortran arrays to MATLAB output
   xqi(1:nq) => mxGetInt32s(plhs(1))
   xqpi(1:nq) => mxGetDoubles(plhs(2))

   ! Call the Fortran function
   call bracket_linear_weight(x, n, xq, nq, xqi, xqpi)

end subroutine mexFunction
