! Initializes PARDISO for solving a real nonsymmetric sparse linear system with a direct solver.
!
! Synopsis:
!  [pt, iparm, dparm] = pardiso_init(number_of_threads);

! Copyright © 2025 Dynare Team
!
! This file is part of Dynare.
!
! Dynare is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Dynare is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
  use iso_fortran_env
  use matlab_mex
  use panua_pardiso
  implicit none (type, external)

  type(c_ptr), dimension(*), intent(in) :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  integer(int32) :: error
  integer(int32), pointer :: iparm(:)

  if (nrhs /= 1 .or. nlhs /= 3) &
       call mexErrMsgTxt("pardiso_init: takes one input argument and 3 output arguments.")

  plhs(1) = mxCreateNumericMatrix(64_mwSize, 1_mwSize, mxINT64_CLASS, mxREAL)
  plhs(2) = mxCreateNumericMatrix(64_mwSize, 1_mwSize, mxINT32_CLASS, mxREAL)
  plhs(3) = mxCreateDoubleMatrix(64_mwSize, 1_mwSize, mxREAL)

  iparm => mxGetInt32s(plhs(2))

  call pardisoinit(mxGetInt64s(plhs(1)), REAL_NONSYM, DIRECT_SOLVER, iparm, &
       mxGetDoubles(plhs(3)), error)
  if (error /= 0) call mexErrMsgTxt("pardisoinit failed: " // error_string(error))

  ! Set number of threads
  if (.not. (mxIsScalar(prhs(1)) .and. mxIsNumeric(prhs(1)))) &
       call mexErrMsgTxt("First input argument must be a numeric scalar")
  iparm(3) = int(mxGetScalar(prhs(1)), int32)
end subroutine mexFunction
