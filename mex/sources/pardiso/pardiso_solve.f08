! Solve a real nonsymmetric linear system A·x=b using PARDISO
!
! Synopsis:
!  [x, pt, iparm, dparm] = pardiso_solve(A, b, pt, iparm, dparm[, phase]);
!
! If omitted, phase defaults to 13 (analysis, numerical factorization, solve,
! iterative refinement).
!
! NB: the MEX takes care of flipping iparm(4) in order to handle the different
! sparse representations between MATLAB/Octave and PARDISO. iparm(4) should
! therefore be 0 on input (and thus output), unless one really wants to solve
! Aᵀ·x=b.
!
! NB: the present MEX and the pardiso_init MEX could easily be extended to
! solve other types of systems (symmetric, complex)

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

  type(c_ptr), dimension(*), intent(in), target :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  integer(int64), pointer :: pt(:)
  integer(int32), pointer :: iparm(:), perm(:) => null()
  real(real64), pointer :: dparm(:)

  integer(int32) :: phase, error, n, nunknowns, i
  integer(int32), allocatable, dimension(:) :: ia, ja
  real(real64), pointer :: a(:), b(:,:), x(:,:)

  ! Various sanity checks

  if ((nrhs /= 5 .and. nrhs /= 6) .or. nlhs /= 4) &
       call mexErrMsgTxt("pardiso_solve: takes 5 or 6 input arguments and exactly 4 output argument")

  if (.not. (mxIsDouble(prhs(1)) .and. mxIsSparse(prhs(1)))) &
       call mexErrMsgTxt("First argument must be a sparse real matrix")

  if (.not. (mxIsDouble(prhs(2)) .and. .not. mxIsSparse(prhs(2)))) &
       call mexErrMsgTxt("Second argument must be a dense real vector or matrix")

  if (.not. (mxIsInt64(prhs(3)) .and. mxGetNumberOfElements(prhs(3)) == 64)) &
       call mexErrMsgTxt("Third argument must be an int64 array of size 64")

  if (.not. (mxIsInt32(prhs(4)) .and. mxGetNumberOfElements(prhs(4)) == 64)) &
       call mexErrMsgTxt("Fourth argument must be an int32 array of size 64")

  if (.not. (mxIsDouble(prhs(5)) .and. mxGetNumberOfElements(prhs(5)) == 64)) &
       call mexErrMsgTxt("Fifth argument must be an double precision array of size 64")

  if (nrhs >= 6) then
     if (.not. (mxIsScalar(prhs(6)) .and. mxIsNumeric(prhs(6)))) &
          call mexErrMsgTxt("Sixth argument must be a numeric scalar")
  end if

  ! Map input and output arguments from MATLAB to PARDISO

  n = int(mxGetN(prhs(1)), int32)
  nunknowns = int(mxGetN(prhs(2)), int32)
  a => mxGetDoubles(prhs(1))
  associate (nnz => int(mxGetNzmax(prhs(1)), int32), ir => mxGetIr(prhs(1)), jc => mxGetJc(prhs(1)))
    allocate(ia(n+1), ja(nnz))
    do i = 1,n+1
       ia(i) = int(jc(i), int32) + 1
    end do
    do i = 1,nnz
       ja(i) = int(ir(i), int32) + 1
    end do
  end associate

  b(1:n,1:nunknowns) => mxGetDoubles(prhs(2))

  plhs(1) = mxCreateDoubleMatrix(int(n, mwSize), int(nunknowns, mwSize), mxREAL)
  x(1:n,1:nunknowns) => mxGetDoubles(plhs(1))

  plhs(2) = mxDuplicateArray(prhs(3))
  plhs(3) = mxDuplicateArray(prhs(4))
  plhs(4) = mxDuplicateArray(prhs(5))

  pt => mxGetInt64s(plhs(2))
  iparm => mxGetInt32s(plhs(3))
  dparm => mxGetDoubles(plhs(4))

  if (nrhs >= 6) then
     phase = int(mxGetScalar(prhs(6)), int32)
  else
     phase = ANALYSIS_NUM_FACT_SOLVE_REFINE
  end if

  if ((iparm(12) == 0 .and. mxGetN(prhs(1)) /= mxGetM(prhs(2))) .or. &
       (iparm(12) == 1 .and. mxGetM(prhs(1)) /= mxGetM(prhs(2)))) &
       call mexErrMsgTxt("First and second arguments have inconsistent dimensions")

  ! Flip the transpose flag, since PARDISO uses the CSR representation
  ! internally, and we will pass it a matrix in CSC representation (used by
  ! MATLAB/Octave)
  iparm(12) = 1 - iparm(12)

  ! Since perm is not initialized, check that PARDISO will not try to read it
  if (iparm(5) == 1) &
       call mexErrMsgTxt("Support for custom fill-in reducing ordering (iparm(5)=1) is not implemented")

  ! Compute the solution
  call pardiso(pt, 1, 1, REAL_NONSYM, phase, n, a, ia, ja, perm, nunknowns, iparm, &
       MESSAGE_LEVEL_OFF, b, x, error, dparm)
  if (error /= 0) call mexErrMsgTxt("pardiso failed: " // error_string(error))

  ! Restore the transpose flag
  iparm(12) = 1 - iparm(12)

end subroutine mexFunction
