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

module panua_pardiso
  use iso_fortran_env
  implicit none (type, external)

  ! Matrix types
  integer(int32), parameter :: REAL_SYM = 1, REAL_SYM_POSDEF = 2, REAL_SYM_INDEF = -2, &
    COMPLEX_STRUCT_SYM = 3, COMPLEX_HERM_POSDEF = 4, COMPLEX_HERM_INDEF = -4, COMPLEX_SYM = 6, &
    REAL_NONSYM = 11, COMPLEX_NONSYM = 13

  ! Verbosity levels
  integer(int32), parameter :: MESSAGE_LEVEL_OFF = 0, MESSAGE_LEVEL_ON = 1

  ! Available solvers
  integer(int32), parameter :: DIRECT_SOLVER = 0, ITERATIVE_SOLVER = 1

  ! Phases
  integer(int32), parameter :: ANALYSIS = 11, ANALYSIS_NUM_FACT = 12, &
       ANALYSIS_NUM_FACT_SOLVE_REFINE = 13, NUM_FACT = 22, SELECTED_INVERSION = -22, &
       NUM_FACT_SOLVE_REFINE = 23, SOLVE_ITERATIVE_REFINE = 33, &
       SOLVE_ITERATIVE_REFINE_ONLY_FORWARD = 331, SOLVE_ITERATIVE_REFINE_ONLY_DIAG = 332, &
       SOLVE_ITERATIVE_REFINE_ONLY_BACKWARD = 333, RELEASE_LU_MNUM = 0, RELEASE_ALL = -1

  interface
     subroutine pardisoinit(pt, mtype, solver, iparm, dparm, error)
       import :: int64, int32, real64
       implicit none
       integer(int64), intent(out) :: pt(64)
       integer(int32), intent(in) :: mtype, solver
       integer(int32), intent(out) :: iparm(64), error
       real(real64), intent(out) :: dparm(64)
     end subroutine pardisoinit

     subroutine pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, &
          x, error, dparm)
       import :: int64, int32, real64
       implicit none
       integer(int64), intent(inout) :: pt(64)
       integer(int32), intent(in) :: maxfct, mnum, mtype, phase, n, nrhs, ia(n+1), ja(*), perm(n), &
            msglvl
       integer(int32), intent(inout) :: iparm(64)
       integer(int32), intent(out) :: error
       real(real64), intent(in) :: a(*)
       real(real64), intent(inout) :: b(n,nrhs), dparm(64)
       real(real64), intent(out) :: x(n,nrhs)
     end subroutine pardiso
  end interface

contains

  function error_string(error)
    integer(int32), intent(in) :: error
    character(:), allocatable :: error_string
    select case (error)
    case (0)
       error_string = 'No error'
    case (-1)
       error_string = 'Input inconsistent'
    case (-2)
       error_string = 'Not enough memory'
    case (-3)
       error_string = 'Reordering problem'
    case (-4)
       error_string = 'Zero pivot, numerical factorization or iterative refinement problem'
    case (-5)
       error_string = 'Unclassified (internal) error'
    case (-6)
       error_string = 'Preordering failed'
    case (-7)
       error_string = 'Diagonal matrix problem'
    case (-8)
       error_string = '32-bit integer overflow problem'
    case (-10)
       error_string = 'No license file panua.lic found'
    case (-11)
       error_string = 'License is expired'
    case (-12)
       error_string = 'Wrong username or hostname'
    case (-100)
       error_string = 'Reached maximum number of Krylov-subspace iteration in iterative solver'
    case (-101)
       error_string = 'No sufficient convergence in Krylov-subspace iteration within 25 iterations'
    case (-102)
       error_string = 'Error in Krylov-subspace iteration'
    case (-103)
       error_string = 'Breakdown in Krylov-subspace iteration'
    case default
       error_string = 'unknown'
    end select
  end function error_string

end module panua_pardiso
