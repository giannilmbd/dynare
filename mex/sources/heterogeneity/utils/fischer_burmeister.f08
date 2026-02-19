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
! Fischer-Burmeister reformulation for complementarity problems
!
! For a complementarity problem F(x) = 0 ⟂ lb ≤ x ≤ ub, the FB function transforms
! the problem into an unconstrained system Φ(x) = 0 where:
!   - φ(a,b) = a + b - sqrt(a² + b² + ε)  satisfies φ(a,b) = 0 ⟺ a ≥ 0, b ≥ 0, ab = 0
!   - Lower bound only:  Φ = φ(x-lb, -F)
!   - Upper bound only:  Φ = φ(ub-x, F)
!   - Double bounds:     Φ = φ(x-lb, φ(ub-x, -F))  (nested formulation from LMMCP)

module fischer_burmeister
    use iso_fortran_env, only: real64, int32
    use calibration_types, only: BOUND_NONE, BOUND_LOWER, BOUND_UPPER, BOUND_BOTH
    implicit none (type, external)

    private
    public :: FB_EPS, fb_phi, fb_phi_derivative, fb_transform, fb_transform_jacobian

    ! Regularization parameter to avoid division by zero
    real(real64), parameter :: FB_EPS = 1.0e-14_real64

contains

    !---------------------------------------------------------------------------
    ! Core Fischer-Burmeister function
    ! φ(a,b) = a + b - sqrt(a² + b² + ε)
    ! Property: φ(a,b) = 0 ⟺ a ≥ 0, b ≥ 0, ab = 0
    !---------------------------------------------------------------------------
    pure elemental function fb_phi(a, b) result(phi)
        real(real64), intent(in) :: a, b
        real(real64) :: phi
        real(real64) :: r

        r = sqrt(a*a + b*b + FB_EPS)
        phi = a + b - r
    end function fb_phi

    !---------------------------------------------------------------------------
    ! Derivatives of the Fischer-Burmeister function
    ! ∂φ/∂a = 1 - a/r
    ! ∂φ/∂b = 1 - b/r
    ! where r = sqrt(a² + b² + ε)
    !---------------------------------------------------------------------------
    pure subroutine fb_phi_derivative(a, b, dphi_da, dphi_db)
        real(real64), intent(in) :: a, b
        real(real64), intent(out) :: dphi_da, dphi_db
        real(real64) :: r

        r = sqrt(a*a + b*b + FB_EPS)
        dphi_da = 1.0_real64 - a / r
        dphi_db = 1.0_real64 - b / r
    end subroutine fb_phi_derivative

    !---------------------------------------------------------------------------
    ! Transform residuals F(x) → Φ(x) using Fischer-Burmeister
    !
    ! For each component i:
    !   BOUND_NONE:  Φ_i = -F_i (just negate, no bounds)
    !   BOUND_LOWER: Φ_i = φ(x_i - lb_i, -F_i)
    !   BOUND_UPPER: Φ_i = φ(ub_i - x_i, F_i)
    !   BOUND_BOTH:  Φ_i = φ(x_i - lb_i, φ(ub_i - x_i, -F_i))  (nested)
    !---------------------------------------------------------------------------
    pure subroutine fb_transform(n, x, fvec, lb, ub, bound_type, phi_vec)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: x(:)
        real(real64), intent(in) :: fvec(:)
        real(real64), intent(in) :: lb(:)
        real(real64), intent(in) :: ub(:)
        integer(int32), intent(in) :: bound_type(:)
        real(real64), intent(out) :: phi_vec(:)

        integer(int32) :: i
        real(real64) :: phi_inner

        do i = 1, n
            select case (bound_type(i))
            case (BOUND_NONE)
                ! Unconstrained: just negate
                phi_vec(i) = -fvec(i)

            case (BOUND_LOWER)
                ! Lower bound only: φ(x-lb, -F)
                phi_vec(i) = fb_phi(x(i) - lb(i), -fvec(i))

            case (BOUND_UPPER)
                ! Upper bound only: φ(ub-x, F)
                phi_vec(i) = fb_phi(ub(i) - x(i), fvec(i))

            case (BOUND_BOTH)
                ! Double bounds (nested): φ(x-lb, φ(ub-x, -F))
                phi_inner = fb_phi(ub(i) - x(i), -fvec(i))
                phi_vec(i) = fb_phi(x(i) - lb(i), phi_inner)
            end select
        end do
    end subroutine fb_transform

    !---------------------------------------------------------------------------
    ! Transform Jacobian using FB chain rule
    !
    ! For Φ_i = φ(a_i, b_i), the chain rule gives:
    !   ∂Φ_i/∂x_j = (∂φ/∂a)·(∂a_i/∂x_j) + (∂φ/∂b)·(∂b_i/∂x_j)
    !
    ! For different bound types:
    !   BOUND_NONE:  ∂Φ_i/∂x_j = -∂F_i/∂x_j
    !   BOUND_LOWER: a=x-lb, b=-F → ∂Φ_i/∂x_j = (∂φ/∂a)·δ_ij + (∂φ/∂b)·(-∂F_i/∂x_j)
    !   BOUND_UPPER: a=ub-x, b=F  → ∂Φ_i/∂x_j = (∂φ/∂a)·(-δ_ij) + (∂φ/∂b)·(∂F_i/∂x_j)
    !   BOUND_BOTH:  nested chain rule (see implementation)
    !---------------------------------------------------------------------------
    pure subroutine fb_transform_jacobian(n, x, fvec, fjac, lb, ub, bound_type, phi_jac)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: x(:)
        real(real64), intent(in) :: fvec(:)
        real(real64), intent(in) :: fjac(:,:)  ! [n × n]
        real(real64), intent(in) :: lb(:)
        real(real64), intent(in) :: ub(:)
        integer(int32), intent(in) :: bound_type(:)
        real(real64), intent(out) :: phi_jac(:,:)  ! [n × n]

        integer(int32) :: i, j
        real(real64) :: a, b, dphi_da, dphi_db
        real(real64) :: a_inner, b_inner, dphi_inner_da, dphi_inner_db
        real(real64) :: phi_inner, dphi_outer_da, dphi_outer_db
        real(real64) :: dphiinner_dxj

        do i = 1, n
            select case (bound_type(i))
            case (BOUND_NONE)
                ! Unconstrained: ∂Φ_i/∂x_j = -∂F_i/∂x_j
                do j = 1, n
                    phi_jac(i, j) = -fjac(i, j)
                end do

            case (BOUND_LOWER)
                ! Lower bound only: Φ = φ(x-lb, -F)
                a = x(i) - lb(i)
                b = -fvec(i)
                call fb_phi_derivative(a, b, dphi_da, dphi_db)

                do j = 1, n
                    if (i == j) then
                        ! ∂a/∂x_i = 1, ∂b/∂x_j = -∂F_i/∂x_j
                        phi_jac(i, j) = dphi_da + dphi_db * (-fjac(i, j))
                    else
                        ! ∂a/∂x_j = 0
                        phi_jac(i, j) = dphi_db * (-fjac(i, j))
                    end if
                end do

            case (BOUND_UPPER)
                ! Upper bound only: Φ = φ(ub-x, F)
                a = ub(i) - x(i)
                b = fvec(i)
                call fb_phi_derivative(a, b, dphi_da, dphi_db)

                do j = 1, n
                    if (i == j) then
                        ! ∂a/∂x_i = -1, ∂b/∂x_j = ∂F_i/∂x_j
                        phi_jac(i, j) = dphi_da * (-1.0_real64) + dphi_db * fjac(i, j)
                    else
                        ! ∂a/∂x_j = 0
                        phi_jac(i, j) = dphi_db * fjac(i, j)
                    end if
                end do

            case (BOUND_BOTH)
                ! Double bounds (nested): Φ = φ_outer(x-lb, φ_inner(ub-x, -F))
                ! Inner: φ_inner = φ(ub-x, -F)
                a_inner = ub(i) - x(i)
                b_inner = -fvec(i)
                phi_inner = fb_phi(a_inner, b_inner)
                call fb_phi_derivative(a_inner, b_inner, dphi_inner_da, dphi_inner_db)

                ! Outer: Φ = φ(x-lb, φ_inner)
                a = x(i) - lb(i)
                b = phi_inner
                call fb_phi_derivative(a, b, dphi_outer_da, dphi_outer_db)

                do j = 1, n
                    ! ∂φ_inner/∂x_j = (∂φ_inner/∂a_inner)·(∂a_inner/∂x_j) + (∂φ_inner/∂b_inner)·(∂b_inner/∂x_j)
                    ! where ∂a_inner/∂x_j = -δ_ij, ∂b_inner/∂x_j = -∂F_i/∂x_j
                    if (i == j) then
                        dphiinner_dxj = dphi_inner_da * (-1.0_real64) + dphi_inner_db * (-fjac(i, j))
                    else
                        dphiinner_dxj = dphi_inner_db * (-fjac(i, j))
                    end if

                    ! ∂Φ/∂x_j = (∂φ_outer/∂a)·(∂a/∂x_j) + (∂φ_outer/∂b)·(∂φ_inner/∂x_j)
                    ! where ∂a/∂x_j = δ_ij
                    if (i == j) then
                        phi_jac(i, j) = dphi_outer_da + dphi_outer_db * dphiinner_dxj
                    else
                        phi_jac(i, j) = dphi_outer_db * dphiinner_dxj
                    end if
                end do
            end select
        end do
    end subroutine fb_transform_jacobian

end module fischer_burmeister
