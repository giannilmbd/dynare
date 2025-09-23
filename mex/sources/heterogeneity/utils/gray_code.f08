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
!
! Original author: Normann Rion <normann@dynare.org>

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
        integer(int32), intent(in) :: ind(:,:), stride(:)
        integer(int32), allocatable, intent(out) :: a(:)

        ! Local variables
        integer(int32) :: n, N_sp, k

        N_sp = size(ind, 1, int32)
        n = size(stride, 1, int32)
        allocate(a(N_sp))

        a = 1_int32
        do concurrent (k=1:n)
            a = a + stride(k)*(ind(:, k)-1_int32)
        end do

    end subroutine compute_linear_indices

    ! Compute coefficient updates r_up, r_down and identify hard dimensions
    ! Hard dimensions are those where interpolation weight w = 1.0
    subroutine compute_coefficient_updates(w, r_up, r_down, is_hard)
        real(real64), intent(in) :: w(:,:)
        real(real64), allocatable, intent(out) :: r_up(:,:), r_down(:,:)
        logical, allocatable, intent(out) :: is_hard(:,:)

        ! Local variables
        integer(int32) :: n, N_om, k, j

        N_om = size(w, 1, int32)
        n = size(w, 2, int32)

        allocate(r_up(N_om,n), r_down(N_om,n), is_hard(N_om,n))

        ! Hard dimensions (weight = 1, no interpolation needed)
        is_hard = (w == 1.0_real64)

        ! Coefficient updates for moving between hypercube corners
        r_up = (1-w)/w
        do concurrent(k=1:n, j=1:N_om)
            if (is_hard(j,k)) then
                r_down(j,k) = 1.0_real64
            else
                r_down(j,k) = w(j,k)/(1-w(j,k))
            end if
        end do

    end subroutine compute_coefficient_updates

end module gray_code