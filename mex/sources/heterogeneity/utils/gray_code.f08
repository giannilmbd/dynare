! Copyright © 2025-2026 Dynare Team
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
!
! Original author: Normann Rion <normann@dynare.org>
!
! Gray code utilities for efficient hypercube corner traversal

module gray_code
    use iso_fortran_env, only: real64, int32
    implicit none (type, external)

contains

    ! Generate Gray code sequence and flip indices
    ! This is used for efficient traversal of hypercube corners in multi-dimensional
    ! interpolation and expectation computations
    subroutine generate_flip_indices(n, flip_idx, Kcorn)
        integer(int32), intent(in) :: n
        integer(int32), allocatable, intent(out) :: flip_idx(:)
        integer(int32), intent(out), optional :: Kcorn

        ! Local variables
        integer(int32) :: Kcorn_local, t
        integer(int32), allocatable :: g(:)

        ! Number of corners in n-dimensional hypercube
        Kcorn_local = ishft(1, n)

        ! Generate Gray code sequence
        allocate(g(0:Kcorn_local-1))
        allocate(flip_idx(1:Kcorn_local-1))

        do t = 0, Kcorn_local-1
            g(t) = ieor(t, ishft(t,-1))
        end do

        ! Find which dimension flips between consecutive Gray codes
        do t = 2, Kcorn_local
            flip_idx(t-1) = trailz(ieor(g(t-1), g(t-2)))+1
        end do

        if (present(Kcorn)) Kcorn = Kcorn_local

    end subroutine generate_flip_indices

    ! Compute strides for state-only linear indices
    ! Converts multi-dimensional indices to linear indices
    subroutine compute_strides(dims, stride)
        integer(int32), intent(in) :: dims(:)
        integer(int32), allocatable, intent(out) :: stride(:)

        ! Local variables
        integer(int32) :: n, k

        n = size(dims, 1, int32)
        allocate(stride(n))

        stride(1) = 1_int32
        do k = 2, n
            stride(k) = stride(k-1) * dims(k-1)
        end do

    end subroutine compute_strides

    ! Compute low-corner linear indices from multi-dimensional indices
    subroutine compute_linear_indices(ind, stride, a)
        integer(int32), contiguous, intent(in) :: ind(:,:), stride(:)
        integer(int32), contiguous, intent(inout) :: a(:)

        ! Local variables
        integer(int32) :: n, k

        n = size(stride, 1, int32)
        a = 1_int32
        do concurrent (k=1:n)
            a = a + stride(k)*(ind(:, k)-1_int32)
        end do
    end subroutine compute_linear_indices

    ! Compute coefficient updates r_up, r_down and identify hard dimensions
    ! Hard dimensions are those where interpolation weight w = 1.0 (is_hard_one)
    ! or w = 0.0 (is_hard_zero), requiring special handling to avoid division by zero
    subroutine compute_coefficient_updates(w, r_up, r_down, is_hard_one, is_hard_zero)
        real(real64), contiguous, intent(in) :: w(:,:)
        real(real64), contiguous, intent(inout) :: r_up(:,:), r_down(:,:)
        logical, contiguous, intent(inout) :: is_hard_one(:,:), is_hard_zero(:,:)

        ! Local variables
        integer(int32) :: n, N_om, k, j

        N_om = size(w, 1, int32)
        n = size(w, 2, int32)

        ! Identify hard dimensions
        is_hard_one = (w == 1.0_real64)   ! Lower boundary: use x(ilow) only
        is_hard_zero = (w == 0.0_real64)  ! Upper boundary: use x(ilow+1) only

        ! Coefficient updates for moving between hypercube corners
        ! Safe computation avoiding division by zero
        do concurrent(k=1:n, j=1:N_om)
            if (is_hard_one(j,k)) then
                ! w=1: only lower corner contributes, never flip to upper
                r_up(j,k) = 0.0_real64    ! Not used, but safe default
                r_down(j,k) = 1.0_real64  ! Not used, but safe default
            else if (is_hard_zero(j,k)) then
                ! w=0: only upper corner contributes, start flipped
                r_up(j,k) = 1.0_real64    ! Not used, but safe default
                r_down(j,k) = 0.0_real64  ! Not used, but safe default
            else
                ! Normal interpolation: 0 < w < 1
                r_up(j,k) = (1.0_real64-w(j,k))/w(j,k)
                r_down(j,k) = w(j,k)/(1.0_real64-w(j,k))
            end if
        end do

    end subroutine compute_coefficient_updates

end module gray_code