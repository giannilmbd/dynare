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
!
! Distribution derivative computation (curly D matrices)

module curlyDs
    use iso_fortran_env, only: real64, int32
    use gray_code
    use blas
    implicit none (type, external)

contains

    !---------------------------------------------------------------------------
    ! Compute distribution derivative matrix (curly D) for one-state case
    !
    ! Computes D(i,j',s) = Σⱼ (∂Bⁱ)ⱼ',ⱼ · vᵢ,ⱼ,ₖ,ₛ · ωᵢ,ⱼ
    ! then D(i',j',s) = Σᵢ μᵢ',ᵢ · D(i,j',s)
    !
    ! Arguments:
    !   D      [out]   : Distribution derivative [N_om × S]
    !   ind    [in]    : Bracket indices [N_om × 1]
    !   inv_h  [in]    : Inverse grid spacing [N_om × 1]
    !   a_hat  [in]    : Policy function derivatives [N_om × S × 1]
    !   om     [in]    : Distribution [N_om]
    !   Mu     [in]    : Transition matrix [N_e × N_e]
    !---------------------------------------------------------------------------
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

    !---------------------------------------------------------------------------
    ! Compute distribution derivative matrix (curly D) for n-state case
    !
    ! Uses Gray code traversal for efficient hypercube corner iteration.
    ! Handles hard boundaries (weight = 0 or 1) with special cases.
    !
    ! Arguments:
    !   D      [out]   : Distribution derivative [N_om × S]
    !   ind    [in]    : Multi-dimensional bracket indices [N_om × n]
    !   w      [in]    : Interpolation weights [N_om × n]
    !   inv_h  [in]    : Inverse grid spacing [N_om × n]
    !   a_hat  [in]    : Policy function derivatives [N_om × S × n]
    !   om     [in]    : Distribution [N_om]
    !   mu     [in]    : Transition matrix [N_e × N_e]
    !   dims   [in]    : Grid dimensions for each state [n]
    !---------------------------------------------------------------------------
    subroutine compute_curlyDs_nd(D, ind, w, inv_h, a_hat, om, mu, dims)
        ! Inputs/Outpus
        real(real64), contiguous, target, intent(out) :: D(:,:)
        real(real64), contiguous, intent(in) :: a_hat(:,:,:), w(:,:), inv_h(:,:), mu(:,:), om(:)
        integer(int32), contiguous, intent(in) :: ind(:,:), dims(:)

        ! Local variables
        real(real64), allocatable, target :: acc(:,:)
        real(real64), allocatable :: beta_k(:), beta_0(:), inc(:), gamma(:), rhs(:,:), r_up(:,:), r_down(:,:)
        real(real64), pointer, contiguous :: D_mat(:,:)
        integer(int32) :: N_om, N_e, N_a, S, n, Kcorn, t, i, j, k, l, m, kf, z, lin
        integer(int32), allocatable :: flip_idx(:), a(:)
        integer(int32), allocatable, target :: stride(:)
        logical, allocatable :: sk(:), is_hard_one(:,:), is_hard_zero(:,:)

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
        allocate(a(N_om))
        call compute_linear_indices(ind, stride, a)

        ! Fetch it into 1, ..., N_om
        do concurrent (j=1:N_a, i=1:N_e)
            a((j-1_int32)*N_e+i) = (a((j-1_int32)*N_e+i)-1)*N_e+i
        end do

        ! Coefficient updates and hard dimensions
        allocate(r_up(N_om, n), r_down(N_om, n), is_hard_one(N_om, n), is_hard_zero(N_om, n))
        call compute_coefficient_updates(w, r_up, r_down, is_hard_one, is_hard_zero)

        ! Low corner initial weight
        allocate(beta_0(N_om))
        beta_0 = product(w, dim=2, mask=(.not. is_hard_zero))

        ! ---- Expectations ----
        allocate(sk(n), beta_k(N_om), acc(S,N_om), inc(S), gamma(N_om))
        acc = 0.0_real64
        stride = stride*N_e

        do k=1, n
            beta_k = inv_h(:,k)*beta_0*om
            where (.not. is_hard_zero(:,k)) beta_k = beta_k / w(:,k)
            do i=1, N_e
                do j=1,N_a
                    sk = .false.
                    lin = (j-1)*N_e + i
                    z = count(is_hard_zero(lin,:))
                    m = a(lin)
                    inc = -beta_k(lin)*a_hat(lin,:,k)
                    if (z == 0_int32) acc(:,m) = acc(:,m) + inc
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
                        else
                            if (.not. sk(kf)) then
                                ! lower -> upper on dim kf
                                m = m + stride(kf)
                                sk(kf) = .true.
                                if (is_hard_one(lin,kf)) then
                                    ! No contribution. The number of
                                    ! hard-one dims flipped to upper increases
                                    z = z+1_int32
                                else if (is_hard_zero(lin,kf)) then
                                    ! No contribution. The number of
                                    ! hard-zero dims flipped to lower decreases
                                    z = z-1_int32
                                else
                                    ! Both sides are higher than 0, but
                                    ! there is a contribution if and only
                                    ! if no hard dims is flipped to upper
                                    inc = inc*r_up(lin,kf)
                                end if
                            else
                                ! upper -> lower on dim kf
                                m = m - stride(kf)
                                sk(kf) = .false.
                                if (is_hard_one(lin,kf)) then
                                    ! The number of hard-one dims flipped to upper decreases
                                    z = z-1_int32
                                else if (is_hard_zero(lin,kf)) then
                                    ! The number of hard-zero dims flipped to lower increases
                                    z = z+1_int32
                                else
                                    ! Both sides are higher than 0, but
                                    ! there is a contribution if and only
                                    ! if no hard dims is flipped to upper
                                    inc = inc*r_down(lin,kf)
                                end if
                            end if
                        end if
                        if (z==0_int32) then
                            acc(:,m) = acc(:,m)+inc
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