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
! Parameter calibration using Broyden quasi-Newton solver

module calibration
    use iso_fortran_env, only: real64, int32
    use iso_c_binding
    use matlab_mex               ! For mexPrintf
    use calibration_types
    use aggregate_residuals_tensor, only: agg_resid_tensor

    ! Broyden quasi-Newton solver for parameter calibration
    use broyden_solver, only: broyden_solve, &
        BROYDEN_SUCCESS, BROYDEN_MAX_ITER, BROYDEN_SINGULAR_JAC, &
        BROYDEN_MAX_BACKTRACK, BROYDEN_FCN_ERROR

    implicit none (type, external)

    private
    public :: calibrate

contains

    !---------------------------------------------------------------------------
    ! Calibrate parameters to match targets using Broyden quasi-Newton method
    !
    ! This implements the pattern from Python SSJ hank_1a.py lines 223-224:
    !   (beta, vphi), _ = utils.broyden_solver(res, initial_guess)
    !
    ! where res() calls household_trans.ss() which is like our solve_steady_state()
    !
    ! Returns: status (0=success, >0=error code)
    !---------------------------------------------------------------------------
    function calibrate(input, output) result(status)
        type(calibration_input), intent(inout) :: input
        type(calibration_output), intent(inout) :: output
        integer(int32) :: status

        ! Local variables for Broyden solver
        real(real64), allocatable :: x(:), fvec(:)
        integer(int32) :: info, nfev
        integer :: i
        character(len=256) :: msg

        ! Initialize
        status = 0

        ! Handle no-calibration case: just compute residuals
        if (input%dims%n_free_parameters == 0) then
            status = agg_resid_tensor(input, output)
            output%converged = (status == 0)
            output%iterations = 1
            ! residual_norm computed from ALL aggregate residuals
            output%residual_norm = maxval(abs(output%agg_output%residuals))
            if (input%cal_verbosity == 2) &
                 call mexPrintf('No calibration parameters - computing residuals only'//NEW_LINE('A'))
            return
        end if

        if (input%cal_verbosity == 2) &
             call mexPrintf('=== Parameter Calibration (Broyden) ==='//NEW_LINE('A'))

        ! Allocate local solution and residual vectors
        allocate(x(input%dims%n_free_parameters), fvec(input%dims%n_free_parameters))

        ! Set initial guess (work directly in θ-space, no transformation)
        do i = 1, input%dims%n_free_parameters
            x(i) = input%free_parameters_init(i)
        end do

        ! Print initial guesses and bounds
        if (input%cal_verbosity == 2) then
            call mexPrintf('  Initial parameter guesses:'//NEW_LINE('A'))
            do i = 1, input%dims%n_free_parameters
                write(msg, '(A,A,A,ES12.5,A,ES12.5,A,ES12.5,A)') &
                    '    ', trim(input%free_parameters_names(i)), ' = ', input%free_parameters_init(i), &
                    ' (lb=', input%free_parameters_bounds(1, i), &
                    ', ub=', input%free_parameters_bounds(2, i), ')'
                call mexPrintf(trim(msg)//NEW_LINE('A'))
            end do
        end if

        ! Solve using Broyden method with bounds (clamping)
        info = broyden_solve(input%broyden_ws, x, fvec, residual_wrapper, &
                             input%solver_opts%tol, input%solver_opts%max_iter, nfev, &
                             verbosity=input%cal_verbosity, &
                             lower_bounds=input%free_parameters_bounds(1,:), &
                             upper_bounds=input%free_parameters_bounds(2,:))

        ! Check convergence
        if (info == BROYDEN_SUCCESS) then
            output%converged = .true.
            if (input%cal_verbosity == 2) &
                 call mexPrintf('Calibration converged successfully!'//NEW_LINE('A'))
        else
            output%converged = .false.
            if (input%cal_verbosity == 2) then
                call mexPrintf('Calibration failed to converge'//NEW_LINE('A'))
                call print_broyden_error(info)
            end if
        end if

        ! Extract calibrated parameters (already in θ-space)
        if (.not. allocated(output%params)) allocate(output%params(input%dims%n_free_parameters))
        do i = 1, input%dims%n_free_parameters
            output%params(i) = x(i)
        end do

        ! Print calibrated parameters
        if (input%cal_verbosity == 2 .and. output%converged) then
            call mexPrintf('Calibrated parameters:'//NEW_LINE('A'))
            do i = 1, input%dims%n_free_parameters
                write(msg, '(A,A,A,ES12.5)') '  ', trim(input%free_parameters_names(i)), ' = ', output%params(i)
                call mexPrintf(trim(msg)//NEW_LINE('A'))
            end do
            write(msg, '(A,I4,A)') 'Converged in ', nfev, ' function evaluations'
            call mexPrintf(trim(msg)//NEW_LINE('A'))
        end if

        ! Extract calibration target equations residuals
        if (.not. allocated(output%residuals)) allocate(output%residuals(input%dims%n_free_parameters))
        output%residuals = fvec
        output%residual_norm = maxval(abs(output%residuals))

        ! Store iteration count
        output%iterations = nfev

        if (input%cal_verbosity == 2) call mexPrintf('Calibration complete'//NEW_LINE('A'))

        ! Set status based on convergence
        if (.not. output%converged) status = 10

    contains
        ! Residual wrapper (nested procedure to access input/output)
        subroutine residual_wrapper(nvar, x_params, resid, iflag)
            integer, intent(in) :: nvar
            real(real64), intent(in) :: x_params(nvar)
            real(real64), intent(out) :: resid(nvar)
            integer, intent(inout) :: iflag

            integer(int32) :: ss_status, k, n_target_eqs

            ! Update parameters with current guess (already in θ-space)
            do k = 1, nvar
                input%params(input%free_parameters_ind(k)) = x_params(k)
            end do

            ! Solve full steady state with these parameters (time iteration + distribution + aggregation)
            ss_status = agg_resid_tensor(input, output)

            if (ss_status /= 0) then
                iflag = -1  ! Signal error to Broyden solver
                return
            end if

            ! Populate resid with calibration target equations residuals
            n_target_eqs = size(input%target_equations, 1)
            do k = 1, n_target_eqs
                resid(k) = output%agg_output%residuals(input%target_equations(k))
            end do

            ! Print iteration info if verbosity == 2 (use parameter values from params)
            if (input%cal_verbosity == 2) &
                 call print_calibration_iteration(nvar, input%params, input%free_parameters_ind, &
                                                  input%free_parameters_names)

        end subroutine residual_wrapper

        subroutine print_broyden_error(code)
            integer(int32), intent(in) :: code
            select case (code)
            case (BROYDEN_MAX_ITER)
                call mexPrintf('    Broyden: Maximum iterations exceeded'//NEW_LINE('A'))
            case (BROYDEN_SINGULAR_JAC)
                call mexPrintf('    Broyden: Singular Jacobian encountered'//NEW_LINE('A'))
            case (BROYDEN_MAX_BACKTRACK)
                call mexPrintf('    Broyden: Too many backtracks'//NEW_LINE('A'))
            case (BROYDEN_FCN_ERROR)
                call mexPrintf('    Broyden: Function evaluation error'//NEW_LINE('A'))
            case default
                call mexPrintf('    Broyden: Unknown error code'//NEW_LINE('A'))
            end select
        end subroutine print_broyden_error

        subroutine print_calibration_iteration(n, params, free_parameters_ind, param_names)
            integer, intent(in) :: n
            real(real64), intent(in) :: params(:)  ! Full params array
            integer(int32), intent(in) :: free_parameters_ind(:)  ! Indices of free parameters in params
            character(len=*), intent(in) :: param_names(:)
            character(len=256) :: msg
            integer :: i

            do i = 1, n
                write(msg, '(A,A,A,ES12.5)') '  ', trim(param_names(i)), ' = ', &
                    params(free_parameters_ind(i))
                call mexPrintf(trim(msg)//NEW_LINE('A'))
            end do
        end subroutine print_calibration_iteration

    end function calibrate

end module calibration
