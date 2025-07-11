! Copyright © 2019-2023 Dynare Team
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

module lapack
  use blas
  implicit none (type, external)

  interface
     subroutine dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
       import :: blint, real64
       implicit none
       integer(blint), intent(in) :: m, n, nrhs, lda, ldb, lwork
       real(real64), dimension(*), intent(inout) :: a, b
       real(real64), dimension(*), intent(out) :: s, work
       real(real64), intent(in) :: rcond
       integer(blint), dimension(*), intent(out) :: iwork
       integer(blint), intent(out) :: rank, info
     end subroutine dgelsd

     subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
       import :: blint, real64
       implicit none
       integer(blint), intent(in) :: n, nrhs, lda, ldb
       real(real64), dimension(*), intent(inout) :: a, b
       integer(blint), dimension(*), intent(out) :: ipiv
       integer(blint), intent(out) :: info
     end subroutine dgesv

     subroutine dgges(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, &
          alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, &
          info)
       import :: blint, bllog, real64
       implicit none
       character, intent(in) :: jobvsl, jobvsr, sort
       interface
          logical(bllog) function selctg(alphar, alphai, beta)
            import :: bllog, real64
            real(real64), intent(in) :: alphar, alphai, beta
          end function selctg
       end interface
       integer(blint), intent(in) :: n, lda, ldb, ldvsl, ldvsr, lwork
       real(real64), dimension(*), intent(inout) :: a, b
       real(real64), dimension(*), intent(out) :: alphar, alphai, beta, vsl, vsr, work
       logical(bllog), dimension(*), intent(out) :: bwork
       integer(blint), intent(out) :: sdim, info
     end subroutine dgges

     subroutine dgecon(norm, n, a, lda, anorm, rcond, work, iwork, info)
       import :: blint, real64
       implicit none
       character, intent(in) :: norm
       integer(blint), intent(in) :: n, lda
       real(real64), dimension(*), intent(in) :: a
       real(real64), intent(in) :: anorm
       real(real64), intent(out) :: rcond
       real(real64), dimension(*), intent(out) :: work
       integer(blint), dimension(*), intent(out) :: iwork
       integer(blint), intent(out) :: info
     end subroutine dgecon

     subroutine dgetrf(m, n, a, lda, ipiv, info)
       import :: blint, real64
       implicit none
       integer(blint), intent(in) :: m, n, lda
       real(real64), dimension(*), intent(inout) :: a
       integer(blint), dimension(*), intent(out) :: ipiv
       integer(blint), intent(out) :: info
     end subroutine dgetrf

     subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
       import :: blint, real64
       implicit none
       integer(blint), intent(in) :: n, lda, lwork
       real(real64), dimension(*), intent(inout) :: a
       integer(blint), dimension(*), intent(in):: ipiv
       real(real64), dimension(*), intent(out) :: work
       integer(blint), intent(out) :: info
     end subroutine dgetri

     subroutine dpotrf(uplo, n, a, lda, info)
       import :: blint, real64
       implicit none
       character, intent(in) :: uplo
       integer(blint), intent(in) :: n, lda
       real(real64), dimension(*), intent(inout) :: a
       integer(blint), intent(out) :: info
     end subroutine dpotrf

     function dlange(norm, m, n, a, lda, work)
       import :: blint, real64
       implicit none
       character, intent(in) :: norm
       integer(blint), intent(in) :: m, n, lda
       real(real64), dimension(*), intent(in) :: a
       real(real64), dimension(*), intent(out) :: work
       real(real64) :: dlange
     end function dlange

     function ilaenv(ispec, name, opts, n1, n2, n3, n4)
       import :: blint  
       implicit none
       integer(blint) :: ispec, n1, n2, n3, n4
       character(len=*), intent(in) :: name, opts
       integer(blint) :: ilaenv
     end function ilaenv
  end interface

contains

  ! Computing inv(A)*B using lapack DGESV
  ! B <- inv(A)*B
  subroutine left_divide(A, B, ipiv, info)
     real(real64), dimension(:,:), intent(inout) :: A
     real(real64), dimension(:,:), intent(inout) :: B
     integer(blint), dimension(:), intent(inout) :: ipiv
     integer(blint), intent(inout) :: info
     integer(blint) :: m, n
#ifdef DEBUG
     if (size(ipiv) /= size(A,1)) then 
        print *, "Inconsistent number of rows of A and size of ipiv"
     end if
     if (size(A,2) /= size(B,1)) then
        print *, "Inconsistent number of columns of A and number of rows of B"
     end if
#endif
     m = int(size(A,1), blint)
     n = int(size(B,2), blint)
     call dgesv(m, n, A, m, ipiv, B, m, info)
  end subroutine left_divide

  ! Matlab-like norm function calling lapack DLANGE routine
  real(real64) function norm(A, c)
     real(real64), dimension(:,:), intent(in) :: A
     character, intent(in) :: c
     real(real64), dimension(:), allocatable :: work
     integer(blint) :: nrow
     nrow = int(size(A,1),blint)
     norm = dlange(c, nrow, int(size(A,2), blint), &
                  &A, nrow, work)
  end function norm

end module lapack
