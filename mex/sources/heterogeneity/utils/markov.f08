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

module markov
    use iso_fortran_env, only: real64, int32
    implicit none (type, external)
contains

    ! Function to compute the Markov transition matrix using Rouwenhorst method
    subroutine markov_rouwenhorst(rho, sigma, N, tol, maxiter, y, p_vec, P_mat)
        integer(int32), intent(in) :: N, maxiter
        real(real64), intent(in) :: rho, sigma, tol
        real(real64), dimension(N), intent(inout) :: y, p_vec
        real(real64), dimension(N, N), intent(inout) :: P_mat
        
        integer(int32) :: i, k
        real(real64) :: p, scale_factor
        real(real64), allocatable :: P1(:, :), P2(:, :), P3(:, :), P4(:, :)
        real(real64), allocatable :: s(:)

        ! Step 1: Initialize transition matrix for N=2
        p = (1.0_real64 + rho) / 2.0_real64
        P_mat(1, 1) = p
        P_mat(1, 2) = 1.0_real64 - p
        P_mat(2, 1) = 1.0_real64 - p
        P_mat(2, 2) = p

        allocate(P1(N, N), P2(N, N), P3(N, N), P4(N, N))

        ! Step 2: Build transition matrix for N >= 3 using recursion
        do k = 3, N
            P1 = 0.0_real64
            P2 = 0.0_real64
            P3 = 0.0_real64
            P4 = 0.0_real64

            ! Construct new transition matrix
            P1(1:k-1, 1:k-1) = p * P_mat(1:k-1,1:k-1)
            P2(1:k-1, 2:k)   = (1.0_real64 - p) * P_mat(1:k-1,1:k-1)
            P3(2:k, 1:k-1)   = (1.0_real64 - p) * P_mat(1:k-1,1:k-1)
            P4(2:k, 2:k)     = p * P_mat(1:k-1,1:k-1)
            P_mat(1:k,1:k) = P1(1:k,1:k) + P2(1:k,1:k) + P3(1:k,1:k) + P4(1:k,1:k)

            ! Divide middle rows by 2
            do i = 2, k - 1
                P_mat(i, :) = P_mat(i, :) / 2.0_real64
            end do

        end do

        deallocate(P1, P2, P3, P4)

        ! Step 3: Compute stationary distribution p_vec
        call compute_stationary(P_mat, N, p_vec, tol, maxiter)

        ! Step 4: Construct the grid (Equivalent of `np.linspace(-1, 1, N)`)
        allocate(s(N))
        do i = 1, N
            s(i) = -1.0_real64 + (i - 1) * (2.0_real64 / (N - 1))
        end do

        ! Step 5: Scale `s` using sigma and variance correction
        scale_factor = sigma / sqrt(compute_variance(s, p_vec, N))
        s = s * scale_factor

        ! Step 6: Compute y values
        do i = 1, N
            y(i) = exp(s(i))
        end do

        ! Normalize y so that sum(p_vec * y) = 1
        y = y / sum(p_vec * y)

        deallocate(s)
    end subroutine markov_rouwenhorst

    ! Function to compute stationary distribution of a transition matrix
    subroutine compute_stationary(P_mat, N, p_vec, tol, maxiter)
        integer(int32), intent(in) :: N, maxiter
        real(real64), intent(in) :: P_mat(N, N), tol
        real(real64), intent(out) :: p_vec(N)
        integer(int32) :: iter
        real(real64) :: diff, temp(N)

        p_vec = 1.0_real64 / real(N, real64)  ! Start with a uniform guess

        do iter = 1, maxiter
            temp = matmul(p_vec, P_mat)
            diff = maxval(abs(temp - p_vec))
            p_vec = temp
            if (diff < tol) exit
        end do
    end subroutine compute_stationary

    ! Function to compute variance using a probability distribution
    function compute_variance(s, p_vec, N) result(var)
        integer(int32), intent(in) :: N
        real(real64), intent(in) :: s(N), p_vec(N)
        real(real64) :: mean_s, var

        mean_s = sum(s * p_vec)
        var = sum(p_vec * (s - mean_s)**2)
    end function compute_variance

end module markov