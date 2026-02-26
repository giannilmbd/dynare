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
! Main entry point for heterogeneous-agent steady-state computation

module aggregate_residuals_tensor
    use iso_fortran_env, only: real64, int32
    use iso_c_binding, only: C_DOUBLE, C_INT
    use calibration_types
    use matlab_mex
    use time_iteration_tensor
    use blas

    implicit none (type, external)

    private
    public :: agg_resid_tensor

contains

    !---------------------------------------------------------------------------
    ! Main solver function for tensor product grids (Cases 1 & 2)
    ! Simpler than solve_steady_state since grid is already constructed
    ! Returns: status (0=success, >0=error code)
    !---------------------------------------------------------------------------
    function agg_resid_tensor(input, output) result(status)
        type(calibration_input), intent(inout) :: input
        type(calibration_output), intent(inout) :: output
        integer(int32) :: status

        ! Initialize
        status = 0

        ! Phase 1: Time iteration on fixed tensor grid
        status = solve_time_iteration_tensor( &
            input%dims, input%ti_options, input%ti_verbosity, input%state_var, &
            input%tg_config, &
            input%old_pol, input%old_pol_array, input%new_pol, input%weighted_old_pol, &
            input%gc_cache_next, input%gc_cache_om, input%gc_cache_fi, &
            input%gp_ws, input%mex, &
            output%ti_output)

        ! Phase 2: Forward iteration for distribution
        status = solve_distribution_tensor(input, output%dist_output)

        ! Phase 3: Compute aggregated heterogeneous variables
        call compute_agg_het_variables(input, output)

        ! Phase 4: Compute aggregate residuals for calibration target equations
        call compute_aggregate_residuals(input, output)

        ! Only print "Calibration iteration complete" when there's actual calibration
        if (input%dims%n_free_parameters > 0 .and. input%cal_verbosity == 2) &
             call mexPrintf('  Calibration iteration complete'//NEW_LINE('A'))

    end function agg_resid_tensor

    function solve_distribution_tensor(input, dist_output) result(status)
        type(calibration_input), intent(inout) :: input
        type(distribution_output), intent(inout) :: dist_output
        integer(int32) :: status

        integer(int32) :: N_om, N_a, N_e, j, j_a, iter, m, kf, t, Kcorn, n_states
        real(real64) :: diff_norm
        real(real64), pointer, contiguous :: d_mat(:,:), d_mat_next(:,:)
        character(len=256) :: msg
        real(real64), pointer, contiguous :: inc(:), acc(:), r_up(:,:), r_down(:,:), beta_0(:)
        integer(int32), pointer, contiguous :: z(:), linear_idx(:), stride_shocks(:)
        logical, pointer, contiguous :: is_hard_one(:,:), is_hard_zero(:,:), sk(:)

        ! Initialize status
        status = 0
        dist_output%converged = .false.

        ! Useful size variables
        N_om = input%dims%N_om
        N_a = input%dims%N_a_om
        N_e = input%dims%N_e
        n_states = input%dims%n_states
        Kcorn = input%gc_cache_fi%Kcorn

        ! Create a few pointers to make the code more readable
        inc(1:N_om) => input%gc_cache_fi%inc
        acc(1:N_om) => input%gc_cache_fi%acc
        r_up(1:N_om, 1:n_states) => input%gc_cache_fi%r_up
        r_down(1:N_om, 1:n_states) => input%gc_cache_fi%r_down
        beta_0(1:N_om) => input%gc_cache_fi%beta_0
        z(1:N_om) => input%gc_cache_fi%z
        linear_idx(1:N_om) => input%gc_cache_fi%linear_idx
        stride_shocks(1:n_states) => input%gc_cache_fi%stride_shocks
        is_hard_one(1:N_om, 1:n_states) => input%gc_cache_fi%is_hard_one
        is_hard_zero(1:N_om, 1:n_states) => input%gc_cache_fi%is_hard_zero
        sk(1:n_states) => input%gc_cache_fi%sk

        ! Initialize the distribution to a uniform one
        d_mat(1:N_e, 1:N_a) => dist_output%distribution
        do j_a = 1, N_a
            d_mat(:,j_a) = input%tg_config%p / N_a
        end do
        input%d_next = dist_output%distribution
        d_mat_next(1:N_e, 1:N_a) => input%d_next

        if (input%forward_verbosity == 2) &
             call mexPrintf('  === Forward Iteration ==='//NEW_LINE('A'))

        ! Perform forward iterations
        do iter = 1, input%d_options%max_iter
            ! Perform the forward iteration algorithm
            ! Initialization
            linear_idx = input%gc_cache_fi%corner_idx  ! Reset to initial low-corner indices
            sk = .false.
            z = count(is_hard_zero, dim=2)

            ! Start at the low corner
            inc = beta_0*dist_output%distribution
            acc = 0.0_real64
            do j = 1, N_om
                if (z(j) == 0_int32) then
                    m = linear_idx(j)
                    acc(m) = acc(m) + inc(j)
                end if
            end do

            ! Go over the other corners
            do t = 1, input%gc_cache_fi%Kcorn-1
                kf = input%gc_cache_fi%flip_idx(t)
                if (.not. sk(kf)) then
                    ! lower -> upper on dim kf
                    linear_idx = linear_idx+stride_shocks(kf)
                    sk(kf) = .true.
                    do j = 1, N_om
                        if (is_hard_one(j,kf)) then
                            z(j) = z(j)+1_int32
                        else if (is_hard_zero(j,kf)) then
                            z(j) = z(j)-1_int32
                        else
                            inc(j) = inc(j)*r_up(j,kf)
                        end if
                    end do
                else
                    ! upper -> lower on dim kf
                    linear_idx = linear_idx-stride_shocks(kf)
                    sk(kf) = .false.
                    do j = 1, N_om
                        if (is_hard_one(j,kf)) then
                            z(j) = z(j)-1_int32
                        else if (is_hard_zero(j,kf)) then
                            z(j) = z(j)+1_int32
                        else
                            inc(j) = inc(j)*r_down(j,kf)
                        end if
                    end do
                end if
                do j=1, N_om
                    if (z(j) == 0_int32) then
                        m = linear_idx(j)
                        acc(m) = acc(m)+inc(j)
                    end if
                end do
            end do

            ! Weighting with the Markov kernel and update for the next distribution
            d_mat(1:N_e, 1:N_a) => acc
            call matmul_add("N", "N", 1.0_real64, input%tg_config%MuT, d_mat, &
                            0.0_real64, d_mat_next)

            ! Check convergence every check_every iterations
            if (mod(iter, input%d_options%check_every) == 0 .or. iter == input%d_options%max_iter) then
                diff_norm = maxval(abs(dist_output%distribution - input%d_next))

                ! Print per-iteration trace for verbosity == 2
                if (input%forward_verbosity == 2) then
                    write(msg, '(A,I6,A,ES10.3,A,ES10.3)') '  Forward iter ', iter, &
                        ': ||ΔD|| = ', diff_norm, ', tol = ', input%d_options%tol
                    call mexPrintf(trim(msg)//NEW_LINE('A'))
                end if

                if (diff_norm < input%d_options%tol) then
                    dist_output%converged = .true.
                    dist_output%iterations = iter
                    dist_output%residual_norm = diff_norm
                    if (input%forward_verbosity == 2) then
                        write(msg, '(A,I0,A,ES10.3,A)') '  Distribution converged in ', iter, &
                            ' iterations (norm: ', diff_norm, ')'
                        call mexPrintf(trim(msg)//NEW_LINE('A'))
                    end if
                    ! Final update
                    dist_output%distribution = input%d_next
                    exit
                end if
            end if
            ! Update for next iteration
            dist_output%distribution = input%d_next
        end do

        ! Check if converged
        if (.not. dist_output%converged) then
            status = 1
            dist_output%iterations = input%d_options%max_iter
            dist_output%residual_norm = diff_norm
            if (input%forward_verbosity >= 1) then
                write(msg, '(A,I0,A,ES10.3,A)') '    WARNING: Distribution did not converge after ', &
                    input%d_options%max_iter, ' iterations (norm: ', diff_norm, ')'
                call mexPrintf(trim(msg)//NEW_LINE('A'))
            end if
        end if

    end function solve_distribution_tensor

    subroutine compute_agg_het_variables(input, output)
        type(calibration_input), intent(inout) :: input
        type(calibration_output), intent(inout) :: output

        integer(int32) :: k, N_Ix

        if (input%cal_verbosity == 2) &
             call mexPrintf('  === Aggregate variables and residuals ==='//NEW_LINE('A'))

        ! Compute the aggregated heterogeneous policy functions
        ! Ix_j = sum_i D_i * policy_j(state_i)
        call matvecmul_add("N", 1.0_real64, output%ti_output%policies_om, &
                           output%dist_output%distribution, 0.0_real64, output%agg_output%Ix)

        ! Fill yagg - subset corresponding to SUM operators
        ! yagg = Ix(indices.Ix.in_x)
        input%yagg = output%agg_output%Ix(input%Ix_in_het)
        ! Fill the SUM_x variables in y
        N_Ix = size(input%Ix_in_agg,1)
        do k=1,N_Ix
            input%y(input%dims%n_agg_endo+input%Ix_in_agg(k)) = input%yagg(k)
        end do
    end subroutine compute_agg_het_variables

    subroutine compute_aggregate_residuals(input, output)
        type(calibration_input), intent(inout) :: input
        type(calibration_output), intent(inout) :: output
        character(len=256) :: msg
        integer :: i

        if (input%dims%n_agg_endo > 0_int32) &
             ! Use MATLAB MEX function
             call call_matlab_resid_tensor(output%agg_output%residuals, input%mex)

        ! Print all aggregate residuals if verbosity == 2
        if (input%cal_verbosity == 2) then
            call mexPrintf('  Aggregate residuals:'//NEW_LINE('A'))
            do i = 1, size(input%equation_names, 1)
                if (len_trim(input%equation_names(i)) > 0) then
                    ! Use equation name if available
                    write(msg, '(A,A,A,ES15.6)') '    ', trim(input%equation_names(i)), ': ', &
                        output%agg_output%residuals(i)
                else
                    ! Fall back to equation number
                    write(msg, '(A,I3,A,ES15.6)') '  ', i, ': ', &
                        output%agg_output%residuals(i)
                end if
                call mexPrintf(trim(msg)//NEW_LINE('A'))
            end do
            write(msg, '(A,ES15.6,A,ES15.6)') '  Max norm: ', &
                maxval(abs(output%agg_output%residuals(1:input%dims%n_agg_endo))), &
                ', tol = ', input%solver_opts%tol
            call mexPrintf(trim(msg)//NEW_LINE('A'))
        end if

    end subroutine compute_aggregate_residuals

    !---------------------------------------------------------------------------
    ! Call dynamic_het1_resid via mexCallMATLAB
    !---------------------------------------------------------------------------
    subroutine call_matlab_resid_tensor(residuals, input_mex)
        type(matlab_mex_handles), intent(in) :: input_mex
        real(real64), contiguous, intent(inout) :: residuals(:)

        ! MATLAB interface variables
        type(c_ptr), dimension(7) :: prhs
        type(c_ptr), dimension(1) :: plhs
        integer(C_INT) :: retval
        real(real64), pointer, contiguous :: resid(:)

        ! Build input array for mexCallMATLAB (7 inputs, no T management)
        prhs(1) = input_mex%y_mx
        prhs(2) = input_mex%x_mx
        prhs(3) = input_mex%params_mx
        prhs(4) = input_mex%y_mx
        prhs(5) = input_mex%yagg_mx

        ! Call MATLAB function: resid = model_name.dynamic_resid(y, x, params, ss, yagg)
        retval = mexCallMATLAB(1_C_INT, plhs, 5_C_INT, prhs, input_mex%agg_resid)

        if (retval /= 0) &
             call mexErrMsgTxt("MATLAB fallback: Failed to call " // input_mex%agg_resid)

        ! Extract residual from MATLAB output
        resid(1:size(residuals)) => mxGetDoubles(plhs(1))
        residuals = resid

    call mxDestroyArray(plhs(1))

    end subroutine call_matlab_resid_tensor

end module aggregate_residuals_tensor
