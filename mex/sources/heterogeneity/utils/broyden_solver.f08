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
! Broyden quasi-Newton solver for heterogeneous-agent steady-state calibration
!
! Implements the Broyden method following the Python SSJ implementation
! (sequence_jacobian/utilities/solvers.py). Uses LAPACK dgesv for linear solves.

module broyden_solver
    use iso_fortran_env, only: real64, int32
    use ieee_arithmetic, only: ieee_is_finite
    use blas, only: blint, matvecmul_add, rank_one_update
    use lapack, only: dgesv
    use calibration_types, only: broyden_workspace
    use matlab_mex, only: mexPrintf

    implicit none (type, external)

    private
    public :: broyden_solve

    ! Default solver parameters
    real(real64), parameter :: DEFAULT_FD_STEP = 1.0e-5_real64      ! Finite difference step
    real(real64), parameter :: DEFAULT_BACKTRACK_C = 0.5_real64     ! Backtracking factor
    integer(int32), parameter :: DEFAULT_MAX_BACKTRACK = 30         ! Max backtracks per iteration

    ! Return codes
    integer(int32), parameter, public :: BROYDEN_SUCCESS = 0
    integer(int32), parameter, public :: BROYDEN_MAX_ITER = 1
    integer(int32), parameter, public :: BROYDEN_SINGULAR_JAC = 2
    integer(int32), parameter, public :: BROYDEN_MAX_BACKTRACK = 3
    integer(int32), parameter, public :: BROYDEN_FCN_ERROR = 4

    ! Abstract interface for residual function callback
    abstract interface
        subroutine residual_fcn(n, x, fvec, iflag)
            use iso_fortran_env, only: real64
            integer, intent(in) :: n
            real(real64), intent(in) :: x(n)
            real(real64), intent(out) :: fvec(n)
            integer, intent(inout) :: iflag
        end subroutine residual_fcn
    end interface

contains

    !---------------------------------------------------------------------------
    ! Broyden quasi-Newton solver
    !
    ! Solves f(x) = 0 using Broyden's method with finite-difference initial
    ! Jacobian and rank-1 updates.
    !
    ! Arguments:
    !   ws       : Pre-allocated Broyden workspace
    !   x        : On entry: initial guess; on exit: solution
    !   fvec     : On exit: residual at solution
    !   fcn      : Residual function f(x) to solve for f(x) = 0
    !   tol      : Convergence tolerance (max|fvec| < tol)
    !   max_iter : Maximum number of Broyden iterations
    !   nfev     : (out) Number of function evaluations
    !   fd_step  : (optional) Finite difference step size
    !   backtrack_c : (optional) Backtracking factor
    !   verbosity : (optional) Verbosity level (2 = verbose, 1 = silent)
    !   lower_bounds : (optional) Lower bounds for variables
    !   upper_bounds : (optional) Upper bounds for variables
    !
    ! Returns: status code (0 = success)
    !---------------------------------------------------------------------------
    function broyden_solve(ws, x, fvec, fcn, tol, max_iter, nfev, fd_step, backtrack_c, verbosity, &
                           lower_bounds, upper_bounds) result(info)
        type(broyden_workspace), intent(inout) :: ws
        real(real64), intent(inout) :: x(:)
        real(real64), intent(out) :: fvec(:)
        procedure(residual_fcn) :: fcn
        real(real64), intent(in) :: tol
        integer(int32), intent(in) :: max_iter
        integer(int32), intent(out) :: nfev
        real(real64), intent(in), optional :: fd_step, backtrack_c
        integer(int32), intent(in), optional :: verbosity
        real(real64), intent(in), optional :: lower_bounds(:), upper_bounds(:)
        integer(int32) :: info

        ! Local variables
        integer(int32) :: n, iter, bcount, iflag
        real(real64) :: h, bc, max_residual
        integer(blint) :: n_bl, nrhs, lda, ldb, lapack_info
        real(real64), allocatable :: x_new(:)
        logical :: has_bounds

        n = ws%n_vars
        h = DEFAULT_FD_STEP
        bc = DEFAULT_BACKTRACK_C
        if (present(fd_step)) h = fd_step
        if (present(backtrack_c)) bc = backtrack_c

        ! Check if bounds are provided
        has_bounds = present(lower_bounds) .and. present(upper_bounds)
        allocate(x_new(n))

        nfev = 0
        info = BROYDEN_SUCCESS

        ! Step 1: Evaluate initial residual f(x0)
        iflag = 0
        call fcn(n, x, fvec, iflag)
        nfev = nfev + 1

        if (iflag < 0) then
            info = BROYDEN_FCN_ERROR
            return
        end if

        ! Main iteration loop
        do iter = 1, max_iter
            ! Check convergence
            max_residual = maxval(abs(fvec))
            if (max_residual < tol) then
                info = BROYDEN_SUCCESS
                return
            end if

            ! Step 2: Compute Jacobian (first iteration only)
            if (iter == 1) then
                if (present(verbosity)) then
                    if (verbosity == 2) then
                        call mexPrintf('  Computing initial Jacobian via finite differences...'//NEW_LINE('A'))
                    end if
                end if
                call compute_fd_jacobian(ws, x, fvec, fcn, h, nfev, iflag, verbosity)
                if (iflag < 0) then
                    info = BROYDEN_FCN_ERROR
                    return
                end if
            end if

            ! Step 3: Solve J * dx = -fvec
            ! Copy J to J_copy (dgesv overwrites it)
            ws%J_copy = ws%J
            ! Set dx = -fvec (right-hand side)
            ws%dx = -fvec

            ! Call LAPACK dgesv
            n_bl = int(n, blint)
            nrhs = 1_blint
            lda = n_bl
            ldb = n_bl
            call dgesv(n_bl, nrhs, ws%J_copy, lda, ws%ipiv, ws%dx, ldb, lapack_info)

            if (lapack_info /= 0) then
                info = BROYDEN_SINGULAR_JAC
                return
            end if

            ! Step 4: Backtracking line search
            do bcount = 1, DEFAULT_MAX_BACKTRACK
                ! Compute new point: x + dx, then clamp to bounds if present
                x_new = x + ws%dx
                if (has_bounds) call clamp_to_bounds(x_new, lower_bounds, upper_bounds, n)

                ! Try new point
                iflag = 0
                call fcn(n, x_new, ws%fvec_new, iflag)
                nfev = nfev + 1

                if (iflag < 0) then
                    ! Function evaluation failed, backtrack
                    ws%dx = ws%dx * bc
                else
                    ! Success: update Jacobian with Broyden rank-1 update
                    ! df = fvec_new - fvec
                    ws%df = ws%fvec_new - fvec

                    ! Compute actual dx (may differ from ws%dx due to clamping)
                    ws%dx = x_new - x

                    ! Move to new point (must be done BEFORE broyden_update
                    ! because broyden_update overwrites ws%fvec_new as workspace)
                    fvec = ws%fvec_new
                    x = x_new

                    ! Broyden update: J = J + outer((df - J*dx) / ||dx||^2, dx)
                    call broyden_update(ws)

                    exit
                end if
            end do

            if (bcount > DEFAULT_MAX_BACKTRACK) then
                info = BROYDEN_MAX_BACKTRACK
                return
            end if
        end do

        ! Maximum iterations exceeded
        info = BROYDEN_MAX_ITER

    end function broyden_solve

    !---------------------------------------------------------------------------
    ! Clamp variables to bounds
    !
    ! For each variable, if the bound is finite (not -Inf or +Inf), clamp.
    !---------------------------------------------------------------------------
    subroutine clamp_to_bounds(x, lb, ub, n)
        real(real64), intent(inout) :: x(:)
        real(real64), intent(in) :: lb(:), ub(:)
        integer(int32), intent(in) :: n

        integer :: i

        do i = 1, n
            if (ieee_is_finite(lb(i))) x(i) = max(x(i), lb(i))
            if (ieee_is_finite(ub(i))) x(i) = min(x(i), ub(i))
        end do
    end subroutine clamp_to_bounds

    !---------------------------------------------------------------------------
    ! Compute initial Jacobian via forward finite differences
    !
    ! J(:,i) = (f(x + h*e_i) - f(x)) / h
    !---------------------------------------------------------------------------
    subroutine compute_fd_jacobian(ws, x, fvec, fcn, h, nfev, iflag, verbosity)
        type(broyden_workspace), intent(inout) :: ws
        real(real64), intent(in) :: x(:), fvec(:)
        procedure(residual_fcn) :: fcn
        real(real64), intent(in) :: h
        integer(int32), intent(inout) :: nfev
        integer, intent(out) :: iflag
        integer(int32), intent(in), optional :: verbosity

        integer(int32) :: i, n
        real(real64) :: x_orig
        character(len=256) :: msg

        n = ws%n_vars
        iflag = 0

        do i = 1, n
            ! Perturb x(i)
            x_orig = x(i)
            ws%fvec_new(1:n) = x  ! Use fvec_new as temporary x_pert
            ws%fvec_new(i) = x_orig + h

            ! Evaluate f(x + h*e_i)
            call fcn(n, ws%fvec_new, ws%df, iflag)
            nfev = nfev + 1

            if (iflag < 0) return

            ! Forward difference: J(:,i) = (f_pert - f) / h
            ws%J(:, i) = (ws%df - fvec) / h

            ! Debug output
            if (present(verbosity)) then
                if (verbosity == 2) then
                    write(msg, '(A,I2,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') &
                        '    FD[', i, ']: x=', x_orig, ', h=', h, &
                        ', f=', fvec(1), ', f_pert=', ws%df(1)
                    call mexPrintf(trim(msg)//NEW_LINE('A'))
                    write(msg, '(A,ES12.5)') '           J=', ws%J(1, i)
                    call mexPrintf(trim(msg)//NEW_LINE('A'))
                end if
            end if
        end do

    end subroutine compute_fd_jacobian

    !---------------------------------------------------------------------------
    ! Broyden rank-1 update
    !
    ! J = J + outer((df - J*dx) / ||dx||^2, dx)
    !
    ! where df = fvec_new - fvec (already computed in ws%df)
    !       dx = Newton step (in ws%dx)
    !---------------------------------------------------------------------------
    subroutine broyden_update(ws)
        type(broyden_workspace), intent(inout) :: ws

        integer(int32) :: n
        real(real64) :: dx_norm_sq, scale

        n = ws%n_vars

        ! Compute ||dx||^2
        dx_norm_sq = dot_product(ws%dx, ws%dx)

        if (dx_norm_sq < 1.0e-30_real64) then
            ! Skip update if dx is too small to avoid division by zero
            return
        end if

        ! Compute correction = J * dx using BLAS (store in fvec_new as temp)
        ws%fvec_new = 0.0_real64
        call matvecmul_add('N', 1.0_real64, ws%J, ws%dx, 0.0_real64, ws%fvec_new)

        ! fvec_new = (df - J*dx) / ||dx||^2
        scale = 1.0_real64 / dx_norm_sq
        ws%fvec_new = (ws%df - ws%fvec_new) * scale

        ! Rank-1 update: J = J + outer(fvec_new, dx) using BLAS dger
        ! dger: A = alpha * x * y^T + A
        call rank_one_update(n, n, 1.0_real64, ws%fvec_new, ws%dx, ws%J)

    end subroutine broyden_update

end module broyden_solver
