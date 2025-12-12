! Free the internal memory allocated by PARDISO
!
! Synopsis:
!  pardiso_free(pt, iparm, dparm);

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

  integer(int32), pointer, dimension(:) :: ia => null(), ja => null(), perm => null()
  real(real64), pointer :: a(:) => null(), b(:,:) => null(), x(:,:) => null()

  integer(int32) :: error

  if (.false.) print *, c_associated(plhs(1)) ! Unreachable code to silence the warning about unused plhs

  if (nrhs /= 3 .or. nlhs /= 0) &
       call mexErrMsgTxt("pardiso_free: takes three input arguments and no output argument.")

  if (.not. (mxIsInt64(prhs(1)) .and. mxGetNumberOfElements(prhs(1)) == 64)) &
       call mexErrMsgTxt("First argument must be an int64 array of size 64")

  if (.not. (mxIsInt32(prhs(2)) .and. mxGetNumberOfElements(prhs(2)) == 64)) &
       call mexErrMsgTxt("Second argument must be an int32 array of size 64")

  if (.not. (mxIsDouble(prhs(3)) .and. mxGetNumberOfElements(prhs(3)) == 64)) &
       call mexErrMsgTxt("Third argument must be an double precision array of size 64")

  call pardiso(mxGetInt64s(prhs(1)), 1, 1, REAL_NONSYM, RELEASE_ALL, 0_int32, a, ia, ja, &
       perm, 0_int32, mxGetInt32s(prhs(2)), MESSAGE_LEVEL_OFF, b, x, error, mxGetDoubles(prhs(3)))
  if (error /= 0) call mexErrMsgTxt("pardiso failed: " // error_string(error))
end subroutine mexFunction
