! Copyright © 2019-2026 Dynare Team
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

module blas
  use iso_fortran_env
  implicit none (type, external)

#if defined(MATLAB_MEX_FILE) && __SIZEOF_POINTER__ == 8
  integer, parameter :: blint = int64
  integer, parameter :: bllog = 8 ! Logical kind, gfortran-specific
#else
  integer, parameter :: blint = int32
  integer, parameter :: bllog = 4 ! Logical kind, gfortran-specific
#endif

  interface
     subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
       import :: blint, real64
       implicit none
       character, intent(in) :: transa, transb
       integer(blint), intent(in) :: m, n, k, lda, ldb, ldc
       real(real64), dimension(*), intent(in) :: a, b
       real(real64), intent(in) :: alpha, beta
       real(real64), dimension(*), intent(inout) :: c
     end subroutine dgemm

     subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
       import :: blint, real64
       implicit none
       character, intent(in) :: trans
       integer(blint), intent(in) :: m, n, lda, incx, incy
       real(real64), dimension(*), intent(in) :: a, x
       real(real64), intent(in) :: alpha, beta
       real(real64), dimension(*), intent(inout) :: y
     end subroutine dgemv

     subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
       import :: blint, real64
       implicit none
       integer(blint), intent(in) :: m, n, lda, incx, incy
       real(real64), intent(in) :: alpha
       real(real64), dimension(*), intent(in) :: x, y
       real(real64), dimension(*), intent(inout) :: a
     end subroutine dger
  end interface

contains

  ! Updating a matrix using blas DGEMM
  ! C <- alpha*op(A)*op(B) + beta*C
  subroutine matmul_add(opA, opB, alpha, A, B, beta, C)
     character, intent(in) :: opA, opB
     ! The arrays used in BLAS/LAPACK calls are required to be contiguous, to
     ! avoid temporary copies before calling BLAS/LAPACK.
     real(real64), dimension(:,:), contiguous, intent(in) :: A, B
     real(real64), dimension(:,:), contiguous, intent(inout) :: C
     real(real64), intent(in) :: alpha, beta
     integer(blint) :: m, n, k
     if (opA == "N") then
        m = int(size(A,1), blint)
        k = int(size(A,2), blint)
     else
        m = int(size(A,2), blint)
        k = int(size(A,1), blint)
     end if 
     if (opB == "N") then
        n = int(size(B,2), blint)        
     else
        n = int(size(B,1), blint)
     end if
#ifdef DEBUG
     if ( (opA /= "N") .and. (opA /= "T") .and. (opA /= "C") ) then
        print *, "opA must be either N, T or C"
     end if 
     if ( (opB /= "N") .and. (opB /= "T") .and. (opB /= "C") ) then
        print *, "opB must be either N, T or C"
     end if 
     if (((opA == "N") .and. (opB == "N") .and. (size(A,2) /= size(B,1))) .or.&
        &((opA == "N") .and. (opB /= "N") .and. (size(A,2) /= size(B,2))) .or.&
        &((opA /= "N") .and. (opB == "N") .and. (size(A,1) /= size(B,1))) .or.&
        &((opA /= "N") .and. (opB /= "N") .and. (size(A,1) /= size(B,2))))    &
        then
        print *, "Inconsistent number of columns of op(A) and number of rows &
                 &of op(B)"
     end if
     if (m /= size(C,1)) then
        print *, "Inconsistent number of rows of op(A) and number of rows &
                 &of C"
     end if
#endif
     call dgemm(opA, opB, m, n, k, alpha, A, int(size(A, 1), blint), B, &
                int(size(B, 1), blint), beta, C, int(size(C, 1), blint))
  end subroutine

  ! Matrix/vector product update using blas DGEMV
  ! y <- alpha*op(A)*x + beta*y
  subroutine matvecmul_add(op, alpha, A, x, beta, y)
    character, intent(in) :: op
    real(real64), intent(in) :: alpha, beta
    ! The arrays used in BLAS/LAPACK calls are required to be contiguous, to
    ! avoid temporary copies before calling BLAS/LAPACK.
    real(real64), contiguous, intent(in) :: x(:), A(:,:)
    real(real64), contiguous, intent(inout) :: y(:)

    call dgemv(op, int(size(A, 1), blint), int(size(A, 2), blint), alpha, A, &
         int(size(A, 1), blint), x, 1_blint, beta, y, 1_blint)
  end subroutine matvecmul_add

  ! Rank-1 update using BLAS DGER
  ! A <- alpha * x * y^T + A
  subroutine rank_one_update(m, n, alpha, x, y, A)
    integer, intent(in) :: m, n
    real(real64), intent(in) :: alpha
    ! The arrays used in BLAS/LAPACK calls are required to be contiguous, to
    ! avoid temporary copies before calling BLAS/LAPACK.
    real(real64), contiguous, intent(in) :: x(:), y(:)
    real(real64), contiguous, intent(inout) :: A(:,:)

    call dger(int(m, blint), int(n, blint), alpha, x, 1_blint, y, 1_blint, &
              A, int(size(A, 1), blint))
  end subroutine rank_one_update

end module blas
