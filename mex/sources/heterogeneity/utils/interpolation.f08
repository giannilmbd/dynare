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
! Linear interpolation with bracket search for policy functions

module interpolation
    use iso_fortran_env, only: real64, int32
    use calibration_types, only: gray_code_cache
    use gray_code, only: compute_linear_indices
    implicit none (type, external)
contains
    subroutine bracket_linear_weight(x, n, xq, nq, xqi, xqpi)
        !-----------------------------------------------------------------------
        ! Subroutine: bracket_linear_weight
        !
        ! Purpose:
        !   For a given monotonic increasing grid `x` and a vector of query
        !   points `xq`, computes for each query:
        !     - the index `xqi` such that xq ∈ [x(xqi), x(xqi+1)]
        !     - the linear interpolation weight `xqpi` for xqi
        !
        !   This allows for fast and robust piecewise linear interpolation.
        !
        ! Arguments:
        !   x     (in)  [real64, dimension(n)]
        !     Monotonically increasing 1D grid.
        !
        !   n     (in)  [int32]
        !     Number of elements in x (length of x).
        !
        !   xq    (in)  [real64, dimension(nq)]
        !     Query points to locate in the grid.
        !
        !   nq    (in)  [int32]
        !     Number of query points (length of xq).
        !
        !   xqi   (out) [int32, dimension(nq)]
        !     For each xq(i), the index `ilow` such that:
        !       x(ilow) ≤ xq(i) < x(ilow+1)
        !
        !   xqpi  (out) [real64, dimension(nq)]
        !     For each xq(i), the relative weight in the bracketed interval:
        !       xqpi(i) = (x(ilow+1) - xq(i)) / (x(ilow+1) - x(ilow))
        !
        ! Notes:
        !   - Extrapolation is clamped to the first or last interval.
        !   - Grid `x` must be strictly increasing for meaningful output.
        !
        !-----------------------------------------------------------------------

        integer(int32), intent(in) :: n, nq
        real(real64), intent(in)  :: x(n), xq(nq)
        integer(int32), intent(out) :: xqi(nq)
        real(real64), intent(out)   :: xqpi(nq)

        integer(int32) :: iq, ilow, ihigh, imid

        do iq = 1, nq
            if (xq(iq) < x(1)) then
                ilow = 1
                xqi(iq) = ilow
                xqpi(iq) = 1.0_real64
            else if (xq(iq) >= x(n)) then
                ! Clamp to last interval [x(n-1), x(n)] with weight 0 on lower point
                ! This ensures ilow+1 = n is valid for interpolation
                ilow = n - 1
                xqi(iq) = ilow
                xqpi(iq) = 0.0_real64
            else
                ilow = 1
                ihigh = n
                do while (ihigh - ilow > 1)
                    imid = (ihigh + ilow) / 2
                    if (xq(iq) >= x(imid)) then
                        ilow = imid
                    else
                        ihigh = imid
                    end if
                end do
                xqi(iq) = ilow
                xqpi(iq) = (x(ilow+1) - xq(iq)) / (x(ilow+1) - x(ilow))
            end if
        end do
    end subroutine bracket_linear_weight

    subroutine interpolate(out, xq)
        type(gray_code_cache), target, intent(inout) :: out
        real(real64), contiguous, target, intent(in) :: xq(:,:) ! N_x × (N_e ⋅ N_a_sp)

        ! Useful local variables
        integer(int32) :: n, M, N_a_om, N_a_sp, a, kf, t
        real(real64), pointer, contiguous :: mat(:,:) => null(), v(:,:) => null()
        N_a_om = size(out%is_hard_one, 1)
        n = size(out%is_hard_one, 2)
        M = size(out%acc, 1) / N_a_om ! M = N_x ⋅ N_e
        N_a_sp = size(xq, 2) / (M / size(xq,1))
        mat(1:M, 1:N_a_om) => out%acc
        v(1:M, 1:N_a_sp) => xq

        ! Initialization at the low corner
        out%linear_idx = out%corner_idx
        out%sk = .false.
        out%z = count(out%is_hard_zero, dim=2)
        out%beta = out%beta_0
        do concurrent (a=1:N_a_om)
            if (out%z(a) == 0_int32) then
                mat(:,a) = out%beta(a)*v(:,out%linear_idx(a))
            else
                mat(:,a) = 0.0_real64
            end if
        end do

        ! Gray code traversal over state dimensions for the other corners
        do t = 1, out%Kcorn - 1
            kf = out%flip_idx(t)
            if (.not. out%sk(kf)) then
                ! Lower -> upper on dimension kf
                out%linear_idx = out%linear_idx + out%stride_states(kf)
                out%sk(kf) = .true.
                do concurrent (a=1:N_a_om)
                    if (out%is_hard_one(a,kf)) then
                        ! For nodes with a hard-one kf dim, the number of hard-one dims flipped
                        ! to upper increases
                        out%z(a) = out%z(a)+1_int32
                    else if (out%is_hard_zero(a,kf)) then
                        ! For nodes with a hard-zero kf dim, the number of hard-zero dims flipped
                        ! to lower decreases
                        out%z(a) = out%z(a)-1_int32
                    else
                        ! For nodes with a soft kf dim, we update the beta
                        ! coefficient
                        out%beta(a) = out%beta(a)*out%r_up(a,kf)
                    end if
                end do
            else
                ! Upper -> lower on dimension kf
                out%linear_idx = out%linear_idx - out%stride_states(kf)
                out%sk(kf) = .false.
                do concurrent (a=1:N_a_om)
                    if (out%is_hard_one(a,kf)) then
                        ! For nodes with a hard-one kf dim, the number of hard-one dims flipped
                        ! to upper decreases
                        out%z(a) = out%z(a)-1_int32
                    else if (out%is_hard_zero(a,kf)) then
                        ! For nodes with a hard-zero kf dim, the number of hard-zero dims flipped
                        ! to lower increases
                        out%z(a) = out%z(a)+1_int32
                    else
                        ! For nodes with a soft kf dim, we update the beta
                        ! coefficient
                        out%beta(a) = out%beta(a)*out%r_down(a,kf)
                    end if
                end do
            end if
            do a = 1,N_a_om
                if (out%z(a) == 0_int32) mat(:,a) = mat(:,a) + out%beta(a)*v(:,out%linear_idx(a))
            end do
        end do

    end subroutine interpolate

end module interpolation