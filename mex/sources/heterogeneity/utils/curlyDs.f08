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

module curlyDs
    use iso_fortran_env, only: real64, int32
    use gray_code
    use blas
    implicit none (type, external)

contains

    ! Computes D(i, j', s) = sum_j (dB_i)_{j',i, j} * v_{i,j,k,s} * w_{i,j}
    ! and D(i',j',s) = \sum_{i} (Mu)ᵢ',ᵢ D(i, j', s)
    subroutine compute_curlyDs_1d(D, ind, inv_h, a_hat, om, Mu)
        ! Inputs/Outputs
        real(real64), contiguous, target, intent(out) :: D(:,:)
        real(real64), contiguous, intent(in) :: inv_h(:,:), a_hat(:,:,:), Mu(:,:), om(:)
        integer(int32), contiguous, intent(in) :: ind(:,:)

        ! Local variables
        real(real64) :: inc
        real(real64), allocatable, target :: acc(:,:)
        real(real64), pointer, contiguous :: D_mat(:,:), D_rhs(:,:)
        integer(int32) :: N_om, N_e, N_a, S, i, j, k, m, a, a_next

        ! Get dimensions
        N_e = size(Mu, 1, int32)
        N_om = size(D, 1, int32)
        S = size(D, 2, int32)
        N_a = N_om / N_e

        allocate(acc(N_om,S))
        acc = 0.0_real64
        do i=1,N_e
            do concurrent (k=1:S, j=1:N_a)
                m = (j-1)*N_e + i 
                a = (ind(i+(j-1)*N_e,1)-1)*N_e + i
                a_next = a+N_e
                inc = inv_h(m,1)*om(m)*a_hat(m,k,1)
                acc(a,k) = acc(a,k) - inc
                acc(a_next,k) = acc(a_next,k) + inc
            end do
        end do
        D_rhs(1:N_e,1:(N_a*S)) => acc
        D_mat(1:N_e,1:(N_a*S)) => D
        call matmul_add("T", "N", 1.0_real64, Mu, D_rhs, 0.0_real64, D_mat)

    end subroutine compute_curlyDs_1d

    ! Computes D(i, j', s) = sum_k sum_j (dB_i^k)_{j_k',j} * prod_{ℓ≠k}(B_i^ℓ)_{j_ℓ',j} * v_{i,j,k,s} * w_{i,j}
    subroutine compute_curlyDs_nd(D, ind, w, inv_h, a_hat, om, mu, dims)
        ! Inputs/Outpus
        real(real64), contiguous, target, intent(out) :: D(:,:)
        real(real64), contiguous, intent(in) :: a_hat(:,:,:), w(:,:), inv_h(:,:), mu(:,:), om(:)
        integer(int32), contiguous, intent(in) :: ind(:,:), dims(:)

        ! Local variables
        real(real64), allocatable, target :: acc(:,:) 
        real(real64), allocatable :: beta_k(:), beta_0(:), r_up(:,:), r_down(:,:), inc(:), gamma(:), rhs(:,:)
        real(real64), pointer, contiguous :: D_mat(:,:)
        integer(int32) :: N_om, N_e, N_a, S, n, Kcorn, t, i, j, k, l, m, kf, z, lin
        integer(int32), allocatable :: flip_idx(:), a(:)
        integer(int32), allocatable, target :: stride(:)
        logical, allocatable, target :: is_hard(:,:)
        logical, allocatable :: sk(:)

        N_om = size(D, 1, int32)
        S = size(D, 2, int32)
        n = size(dims, 1, int32)
        N_e = size(mu, 1, int32)
        N_a = N_om / N_e

        ! Generate Gray code flip indices
        call generate_flip_indices(n, flip_idx, Kcorn)

        ! ---- Strides for state-only linear indices (1..N_a) ----
        call compute_strides(dims, stride)

        ! Low-corner linear indices in 1, ..., N_a
        call compute_linear_indices(ind, stride, a)

        ! Fetch it into 1, ..., N_om
        do concurrent (j=1:N_a, i=1:N_e)
            a((j-1_int32)*N_e+i) = (a((j-1_int32)*N_e+i)-1)*N_e+i
        end do

        ! Low corner initial weight
        allocate(beta_0(N_om))
        beta_0 = 1.0_real64
        do k=1,n
            beta_0 = beta_0*w(:,k)
        end do

        ! Coefficient updates and hard dimensions
        call compute_coefficient_updates(w, r_up, r_down, is_hard)

        ! ---- Expectations ----
        allocate(sk(n), beta_k(N_om), acc(S,N_om), inc(S), gamma(N_om))
        acc = 0.0_real64
        stride = stride*N_e

        do k=1, n
            beta_k = inv_h(:,k)*beta_0*om / w(:,k)
            do i=1, N_e
                do j=1,N_a
                    sk = .false.
                    z = 0
                    lin = (j-1)*N_e + i
                    m = a(lin)
                    inc = -beta_k(lin)*a_hat(lin,:,k)
                    acc(:,m) = acc(:,m) + inc
                    do t=1, Kcorn-1
                        kf = flip_idx(t)
                        if (kf==k) then
                            inc = -inc
                            if (.not. sk(kf)) then
                                ! lower -> upper on dim kf
                                m = m+stride(kf)
                                sk(kf) = .true.
                            else
                                ! upper -> lower on dim kf
                                m = m-stride(kf)
                                sk(kf) = .false.
                            end if
                            if (z==0_int32) then
                                acc(:,m) = acc(:,m)+inc
                            end if
                        else
                            if (.not. sk(kf)) then
                                ! lower -> upper on dim kf
                                m = m + stride(kf)
                                sk(kf) = .true.
                                if (is_hard(lin,kf)) then
                                    ! No contribution. The number of
                                    ! hard dims flipped to upper increases
                                    z = z+1
                                else
                                    ! Both sides are higher than 0, but
                                    ! there is a contribution if and only
                                    ! if no hard dims is flipped to upper
                                    inc = inc*r_up(lin,kf)
                                    if (z==0_int32) then
                                        acc(:,m) = acc(:,m)+inc
                                    end if
                                end if
                            else
                                ! upper -> lower on dim kf
                                m = m - stride(kf)
                                sk(kf) = .false.
                                if (is_hard(lin,kf)) then
                                    z = z-1
                                    ! The number of hard dims flipped to upper decreases
                                    if (z==0) then
                                        acc(:,m) = acc(:,m)+inc
                                    end if
                                else
                                    ! Both sides are higher than 0, but
                                    ! there is a contribution if and only
                                    ! if no hard dims is flipped to upper
                                    inc = inc*r_down(lin,kf)
                                    if (z==0) then
                                        acc(:,m) = acc(:,m)+inc
                                   end if
                                end if
                            end if
                        end if
                    end do
                end do
            end do
        end do

        ! Weighting with the Markov kernel
        allocate(rhs(S*N_a,N_e))
        do concurrent(l=1:S,j=1:N_a,i=1:N_e)
            rhs((l-1)*N_a+j,i) = acc(l,(j-1)*N_e+i)
        end do
        D_mat(1:N_e,1:(N_a*S))=>D
        call matmul_add("T", "T", 1.0_real64, Mu, rhs, 0.0_real64, D_mat)

    end subroutine compute_curlyDs_nd

end module curlyDs