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
! Expected value and derivative computation for policy functions

module expectations
    use iso_fortran_env, only: real64, int32
    use blas
    use gray_code
    implicit none (type, external)

contains

    ! Compute 𝔼[x̄ₐ|θ,a], one-state case
    subroutine compute_expected_dx_1d(E, x, mu, ind, inv_h)
        ! Inputs/Outpus
        real(real64), contiguous, intent(out) :: E(:,:,:)
        real(real64), contiguous, intent(in), target :: x(:,:)
        real(real64), contiguous, intent(in) :: mu(:,:), inv_h(:,:)
        integer(int32), contiguous, intent(in) :: ind(:,:)

        ! Useful variables
        real(real64), allocatable :: Y_flat(:,:)
        real(real64), target, allocatable :: Dx_flat(:,:)
        real(real64), pointer, contiguous :: x_bar_dash(:,:,:)
        integer(int32) :: N_x, N_sp, N_a, N_e, i, a, j, j_a, j_e

        ! Sizes
        N_x = size(x,1,int32)
        N_sp = size(x,2,int32)
        N_e = size(mu,1,int32)
        N_a = N_sp / N_e

        ! Reshape of x
        x_bar_dash(1:N_x,1:N_e,1:N_a) => x

        ! First difference w.r.t the state
        allocate(Dx_flat(N_x*(N_a-1),N_e))
        do concurrent (j=1:N_a-1, i=1:N_x)
            Dx_flat((j-1)*N_x+i,:) = x_bar_dash(i,:,j+1)-x_bar_dash(i,:,j)
        end do

        ! Weighting with the transition matrix
        allocate(Y_flat(N_x*(N_a-1), N_e))
        call matmul_add("N", "T", 1.0_real64, Dx_flat, mu, 0.0_real64, Y_flat)

        ! Final calculation
        do concurrent (j_a=1:N_a, j_e=1:N_e)
            j = (j_a-1)*N_e+j_e
            a = ind(j,1)
            E(:,1,j) = inv_h(j,1) * Y_flat((a-1)*N_x+1:a*N_x, j_e)
        end do

    end subroutine compute_expected_dx_1d

    ! Compute 𝔼[x̄ₐ|θ,a], n-state case
    subroutine compute_expected_dx_nd(E, x, mu, ind, w, inv_h, dims)
        ! Inputs/Outpus
        real(real64), contiguous, intent(out) :: E(:,:,:)
        real(real64), contiguous, intent(in), target  :: x(:,:), w(:,:)
        real(real64), contiguous, intent(in) :: mu(:,:), inv_h(:,:)
        integer(int32), contiguous, intent(in), target :: ind(:,:)
        integer(int32), contiguous, intent(in) :: dims(:)

        ! Useful variables
        real(real64) :: beta
        real(real64), pointer, contiguous :: x_bar_dash(:,:,:)
        real(real64), allocatable :: Y(:,:), beta_0(:), beta_k(:), r_up(:,:), r_down(:,:), x_bar_dash_mat(:,:)
        integer(int32) :: N_x, N_sp, N_a, N_e, n, Kcorn, t, i, j_a, j_e, j, k, l_a, kf, z
        integer(int32), allocatable :: flip_idx(:), a(:), stride(:)
        logical, allocatable :: is_hard_one(:,:), is_hard_zero(:,:), sk(:)

        N_x = size(x,1,int32)
        N_sp = size(x,2,int32)
        N_e = size(mu,1,int32)
        N_a = N_sp / N_e
        n = size(E,2,int32)

        ! Useful reshape of x
        x_bar_dash(1:N_x,1:N_e,1:N_a) => x

        ! Generate flip indices using Gray code
        call generate_flip_indices(n, flip_idx)
        Kcorn = ishft(1, n)

        ! ---- Strides for state-only linear indices (1..N_a) ----
        call compute_strides(dims, stride)

        ! Low-corner linear indices in 1, ..., N_a
        allocate(a(N_sp))
        call compute_linear_indices(ind, stride, a)

        ! Coefficient updates and hard dimensions
        allocate(r_up(N_sp, n), r_down(N_sp, n), is_hard_one(N_sp, n), is_hard_zero(N_sp, n))
        call compute_coefficient_updates(w, r_up, r_down, is_hard_one, is_hard_zero)

        ! Low corner initial weight
        allocate(beta_0(N_sp))
        beta_0 = product(w, dim=2, mask=(.not. is_hard_zero))

        ! Weighting with the transition matrix
        allocate(Y(N_x*N_a,N_e), x_bar_dash_mat(N_x*N_a,N_e))
        do concurrent (j_a=1:N_a,j_e=1:N_e,i=1:N_x)
            x_bar_dash_mat((j_a-1)*N_x+i,j_e) = x_bar_dash(i,j_e,j_a)
        end do
        call matmul_add("N", "T", 1.0_real64, x_bar_dash_mat, mu, 0.0_real64, Y)

        ! ---- Expectations ----
        allocate(sk(n), beta_k(N_sp))
        do k=1, n
            beta_k = beta_0
            where (.not. is_hard_zero(:,k)) beta_k = beta_k / w(:,k)
            do j_a=1,N_a
                do j_e=1,N_e
                    ! do concurrent (j_a=1:N_a, j_e=1:N_e)
                    ! Linear index
                    j = (j_a-1)*N_e+j_e

                    ! Initialization
                    beta = -inv_h(j,k)*beta_k(j)
                    sk = .false.
                    z = count(is_hard_zero(j,:))
                    l_a = a(j)
                    if (z == 0_int32) then
                        E(:,k,j) = beta*Y((l_a-1)*N_x+1:l_a*N_x,j_e)
                    else
                        E(:,k,j) = 0.0_real64
                    end if

                    ! Remaining corners
                    do t=1, Kcorn-1
                        kf = flip_idx(t)
                        if (kf == k) then
                            beta = -beta
                            if (.not. sk(kf)) then
                                ! lower -> upper on dim kf
                                l_a = l_a + stride(kf)
                                sk(kf) = .true.
                            else
                                ! upper -> lower
                                l_a = l_a - stride(kf)
                                sk(kf) = .false.
                            end if
                        else
                            if (.not. sk(kf)) then
                                ! lower -> upper on dim kf
                                l_a = l_a + stride(kf)
                                sk(kf) = .true.
                                if (is_hard_one(j,kf)) then
                                    ! No contribution. The number of
                                    ! hard-one dims flipped to upper increases
                                    z = z+1_int32
                                else if (is_hard_zero(j,kf)) then
                                    ! No contribution. The number of
                                    ! hard-zero dims flipped to lower decreases
                                    z = z-1_int32
                                else
                                    ! Both sides are higher than 0, but
                                    ! there is a contribution if and only
                                    ! if no hard dims is flipped to upper
                                    beta = beta*r_up(j,kf)
                                end if
                            else
                                ! upper -> lower on dim kf
                                l_a = l_a - stride(kf)
                                sk(kf) = .false.
                                if (is_hard_one(j,kf)) then
                                    z = z-1_int32
                                else if (is_hard_zero(j,kf)) then
                                    z = z+1_int32
                                else
                                    ! Both sides are higher than 0, but
                                    ! there is a contribution if and only
                                    ! if no hard dims is flipped to upper
                                    beta = beta*r_down(j,kf)
                                end if
                            end if
                        end if
                        if (z==0_int32) E(:,k,j) = E(:,k,j)+beta*Y((l_a-1)*N_x+1:l_a*N_x,j_e)
                    end do
                end do
            end do
        end do
    end subroutine compute_expected_dx_nd

    subroutine compute_expected_y(E, y, mu, ind, w, dims)
        ! Inputs/Outpus
        real(real64), contiguous, intent(out) :: E(:,:)
        real(real64), contiguous, intent(in), target  :: y(:,:)
        real(real64), contiguous, intent(in) :: mu(:,:), w(:,:)
        integer(int32), contiguous, intent(in), target :: ind(:,:)
        integer(int32), contiguous, intent(in) :: dims(:)

        ! Local variables
        real(real64), pointer, contiguous :: y_ptr(:,:) => null()
        real(real64), allocatable, target :: y_mat(:,:)
        real(real64), allocatable :: beta(:), r_up(:,:), r_down(:,:)
        integer(int32) :: n, N_e, N_a, N_om, S, Kcorn, t, kf, i, j, m
        integer(int32), allocatable :: flip_idx(:), z(:), stride(:), a(:)
        logical, allocatable :: sk(:), is_hard_one(:,:), is_hard_zero(:,:)

        ! Size variables
        n = size(dims, 1, int32)
        N_e = size(mu, 1, int32)
        N_om = size(y, 1, int32)
        N_a = N_om/N_e
        S = size(y, 2, int32)

        ! Generate flip indices using Gray code
        call generate_flip_indices(n, flip_idx)
        Kcorn = ishft(1, n)

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
        allocate(beta(N_om))
        beta = product(w, dim=2, mask=(.not. is_hard_zero))

        ! Weighting with the Mu matrix
        y_ptr(1:N_e,1:(N_a*S))=>y
        allocate(y_mat(N_e,N_a*S))
        call matmul_add("N", "N", 1.0_real64, mu, y_ptr, 0.0_real64, y_mat)
        y_ptr(1:N_om,1:S)=>y_mat

        ! Computes ℰ(i, j, s) = sum_{i'} sum_{j'} [ ∏_ℓ (B^ℓ)_{j'_ℓ, i, j} ] *
        ! y_{i', j', s} * μ_{i,i'}. Note that the i and j dimensions are merged
        ! in the following calculation
        allocate(sk(n),z(N_om))

        ! Initialization at the low corner
        sk = .false.
        z = count(is_hard_zero, dim=2)
        do concurrent (j=1:N_om, m=1:S)
            if (z(j) == 0_int32) then
                E(j,m) = beta(j)*y_ptr(a(j),m)
            else
                E(j,m) = 0.0_real64
            end if
        end do
        stride = stride * N_e

        ! Going through the remaining corners
        do t=1,Kcorn-1
            kf = flip_idx(t)
            if (.not. sk(kf)) then
                ! lower -> upper on dim kf
                a = a + stride(kf)
                sk(kf) = .true.
                do concurrent (j=1:N_om)
                    if (is_hard_one(j,kf)) then
                        ! For nodes with a hard-one kf dim, the number of hard dims flipped
                        ! to upper increases
                        z(j) = z(j)+1_int32
                    else if (is_hard_zero(j,kf)) then
                        ! For nodes with a hard-zero kf dim, the number of hard dims flipped
                        ! to lower decreases
                        z(j) = z(j)-1_int32
                    else
                        ! For nodes with a soft kf dim, we update the beta
                        ! coefficient
                        beta(j) = beta(j)*r_up(j,kf)
                    end if
                end do
            else
                ! upper -> lower on dim kf
                a = a - stride(kf)
                sk(kf) = .false.
                do concurrent (j=1:N_om)
                    if (is_hard_one(j,kf)) then
                        ! For nodes with a hard-one kf dim, the number of hard dims flipped
                        ! to upper decreases
                        z(j) = z(j)-1_int32
                    else if (is_hard_zero(j,kf)) then
                        ! For nodes with a hard-zero kf dim, the number of hard dims flipped
                        ! to lower increases
                        z(j) = z(j)+1_int32
                    else
                        ! For nodes with a soft kf dim, we update the beta
                        ! coefficient
                        beta(j) = beta(j)*r_down(j,kf)
                    end if
                end do
            end if
            do concurrent (j=1:N_om)
                if (z(j) == 0_int32) then
                    do m=1,S
                        E(j,m) = E(j,m) + beta(j) * y_ptr(a(j),m)
                    end do
                end if
            end do
        end do

    end subroutine compute_expected_y

end module expectations
