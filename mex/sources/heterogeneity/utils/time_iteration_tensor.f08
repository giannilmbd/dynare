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
! Time iteration for tensor product grids

module time_iteration_tensor
    use iso_fortran_env, only: real64, int32
    use iso_c_binding, only: C_DOUBLE, C_INT, C_INT64_T, c_ptr, c_f_procpointer, c_loc, &
                             c_f_pointer, c_associated, c_funloc
    use calibration_types
    use gray_code
    use matlab_mex
    use blas, only: matvecmul_add, matmul_add  ! BLAS matrix-vector and matrix-matrix multiplication
    use interpolation

    ! Trust region solver with Fischer-Burmeister bounds
    use trust_region, only: trust_region_solve
    use fischer_burmeister, only: fb_transform, fb_transform_jacobian, FB_EPS

    implicit none (type, external)

    private
    public :: solve_time_iteration_tensor, set_het_input, call_matlab_het_resid_tensor

contains

    !---------------------------------------------------------------------------
    ! Time iteration on fixed tensor product grid
    !---------------------------------------------------------------------------
function solve_time_iteration_tensor( &
        dims, ti_options, ti_verbosity, state_var, &
        tg_config, &
        old_pol, old_pol_array, new_pol, weighted_old_pol, &
        gc_cache_next, gc_cache_om, gc_cache_fi, &
        gp_ws, mex, &
        ti_output) result(status)

        ! Input parameters
        type(model_dimensions), intent(in) :: dims
        type(time_iteration_options), intent(in) :: ti_options
        integer(int32), intent(in) :: ti_verbosity
        integer(int32), intent(in), contiguous :: state_var(:)
        type(tensor_grid_config), pointer, intent(in) :: tg_config

        ! Policy arrays (inout - modified during iteration)
        real(real64), intent(inout), contiguous :: old_pol(:,:)
        real(real64), intent(inout), contiguous :: old_pol_array(:,:,:)
        real(real64), intent(inout), contiguous :: new_pol(:,:)
        real(real64), intent(inout), contiguous :: weighted_old_pol(:,:)

        ! Gray code caches (inout - modified for interpolation/forward iteration)
        type(gray_code_cache), intent(inout) :: gc_cache_next
        type(gray_code_cache), pointer, intent(inout) :: gc_cache_om
        type(gray_code_cache), pointer, intent(inout) :: gc_cache_fi

        ! Workspace and MATLAB handles
        type(grid_point_workspace), intent(inout) :: gp_ws
        type(matlab_mex_handles), intent(inout) :: mex

        ! Output
        type(time_iteration_output), target, intent(inout) :: ti_output

        ! Local variables
        integer(int32) :: status, iter, k, j, i
        integer(int32) :: n_consecutive_increases
        real(real64) :: diff_norm, prev_diff_norm
        real(real64), pointer :: policy(:)
        character(len=256) :: msg

        status = 0

        if (ti_verbosity == 2) call mexPrintf('  === Time Iteration ==='//NEW_LINE('A'))

        ! Initialize early stopping tracking
        prev_diff_norm = huge(1.0_real64)
        n_consecutive_increases = 0

        ! Time iteration loop
        do iter = 1, ti_options%max_iter

            ! Step 1: Solve household problem at all grid points using trust region solver
            call solve_household_all_points( &
                dims, &
                tg_config%mcp, &
                tg_config%sm, &
                tg_config%pol_grids, &
                tg_config%pol_dims, &
                tg_config%Mu, &
                ti_options, &
                state_var, &
                old_pol, &
                old_pol_array, &
                new_pol, &
                weighted_old_pol, &
                gc_cache_next, &
                gp_ws, &
                mex, &
                status)

            if (status /= 0) &
                 call mexErrMsgTxt('Time iteration failed: Household solver did not converge')

            ! Step 2: Check convergence
            diff_norm = maxval(abs(new_pol - old_pol))

            if (ti_verbosity == 2) then
                write(msg, '(A,I4,A,ES12.4,A,ES12.4)') '  Iter ', iter, &
                    ': ||Δpolicy|| = ', diff_norm, ', tol = ', ti_options%tol
                call mexPrintf(trim(msg)//NEW_LINE('A'))
            end if

            ! Early stopping: track consecutive increases in diff_norm
            if (diff_norm > prev_diff_norm) then
                n_consecutive_increases = n_consecutive_increases + 1
            else
                n_consecutive_increases = 0
            end if
            prev_diff_norm = diff_norm

            ! Check early stopping condition (if enabled)
            if (ti_options%early_stopping > 0 .and. &
                n_consecutive_increases >= ti_options%early_stopping) then
                if (ti_verbosity == 2) then
                    write(msg, '(A,I3,A)') '  Early stopping: ', n_consecutive_increases, &
                        ' consecutive increases in ||Δpolicy||'
                    call mexPrintf(trim(msg)//NEW_LINE('A'))
                end if
                ti_output%policies_sp = new_pol
                ti_output%converged = .false.
                ti_output%iterations = iter
                ti_output%residual_norm = diff_norm
                old_pol = new_pol
                exit
            end if

           if (diff_norm < ti_options%tol) then
                if (ti_verbosity == 2) then
                    write(msg, '(A,I4,A)') '  Converged in ', iter, ' iterations'
                    call mexPrintf(trim(msg)//NEW_LINE('A'))
                end if
                ti_output%policies_sp = new_pol
                ti_output%converged = .true.
                ti_output%iterations = iter
                ti_output%residual_norm = diff_norm
                old_pol = new_pol
                exit
            end if

            ! Update for next iteration (with learning rate dampening)
            old_pol = ti_options%learning_rate * new_pol &
                            + (1.0_real64 - ti_options%learning_rate) * old_pol
        end do

        if (.not. ti_output%converged) then
            if (ti_verbosity == 2) &
                 call mexPrintf('    Warning: Time iteration did not converge'//NEW_LINE('A'))
            ! If we exited the loop without convergence and without early stopping,
            ! we still need to store the last computed policy (not the dampened one)
            if (ti_output%iterations == 0) then
                ! Max iterations reached - store final policy
                ti_output%policies_sp = new_pol
                ti_output%iterations = ti_options%max_iter
                ti_output%residual_norm = diff_norm
                old_pol = new_pol
            end if
        end if

        ! Snap policy values to exact bounds where they are within FB tolerance.
        ! The Fischer-Burmeister regularization (FB_EPS > 0) prevents the solver
        ! from placing constrained variables exactly at their bounds.
        do i = 1, dims%n_orig
            do concurrent (j = 1:dims%N_sp)
                select case (tg_config%mcp%bound_type(i))
                case (BOUND_LOWER, BOUND_BOTH)
                    if (old_pol(i,j) - tg_config%mcp%lower_bounds(i) < sqrt(FB_EPS)) then
                        old_pol(i,j) = tg_config%mcp%lower_bounds(i)
                    end if
                end select
                select case (tg_config%mcp%bound_type(i))
                case (BOUND_UPPER, BOUND_BOTH)
                    if (tg_config%mcp%upper_bounds(i) - old_pol(i,j) < sqrt(FB_EPS)) then
                        old_pol(i,j) = tg_config%mcp%upper_bounds(i)
                    end if
                end select
            end do
        end do
        ti_output%policies_sp = old_pol

        ! Compute the interpolated values of the policy function on the dense grid
        call interpolate(gc_cache_om, old_pol)
        ti_output%policies_om(1:dims%n_het_endo,1:dims%N_om) => gc_cache_om%acc

        ! Compute bracketing indices and weights for each state dimension on the dense grid
        do k = 1, dims%n_states
            policy(1:dims%N_om) => ti_output%policies_om(state_var(k),:)
            call bracket_linear_weight(tg_config%d_grids(k)%nodes, & ! grid
                                       tg_config%d_dims(k), &        ! n_grid
                                       policy, &                           ! query points
                                       dims%N_om, &                  ! n_query
                                       tg_config%pol_ind_om(:, k), & ! indices (out)
                                       tg_config%pol_w_om(:, k))     ! weights (out)
        end do

        ! Update linear indices from multi-dimensional indices in 1,...,N_a
        call compute_linear_indices(tg_config%pol_ind_om, gc_cache_fi%stride_states, &
                                    gc_cache_fi%corner_idx)

        ! Cast the linear indices in 1,...,N_om
        do j=1, dims%N_a_om
            do i=1, dims%N_e
                gc_cache_fi%corner_idx((j-1)*dims%N_e+i) = (gc_cache_fi%corner_idx((j-1)*dims%N_e+i)-1)*dims%N_e+i
            end do
        end do

        ! Update coefficient updates (r_up, r_down, is_hard_one, is_hard_zero)
        call compute_coefficient_updates(tg_config%pol_w_om,     &
                                         gc_cache_fi%r_up,   &
                                         gc_cache_fi%r_down, &
                                         gc_cache_fi%is_hard_one, &
                                         gc_cache_fi%is_hard_zero)

        ! Update initial beta weights (product of all w's at each grid point)
        gc_cache_fi%beta_0 = product(tg_config%pol_w_om, dim=2, mask=(.not. gc_cache_fi%is_hard_zero))

    end function solve_time_iteration_tensor

    !---------------------------------------------------------------------------
    ! Set heterogeneous input: construct yh array and expected future policy
    !
    ! Prepares grid-point-specific inputs for household solver callback:
    ! - Builds yh array [yh(-1), yh(0), yh(+1)] from current policy and states
    ! - Computes E[yh(+1)] via multilinear interpolation over shock realizations
    ! - Computes dpol_next = ∂E[yh(+1)]/∂state for chain rule in Jacobian
    !---------------------------------------------------------------------------
    subroutine set_het_input(policy, j, j_e, dims, mcp, sm, state_var, &
                             pol_grids, pol_dims, weighted_old_pol, &
                             gc_cache_next, gp_ws, input_mex)
        use interpolation, only: bracket_linear_weight

        ! Grid point identification
        real(real64), intent(in) :: policy(:)               ! Current policy guess [n_orig]
        integer(int32), intent(in) :: j                     ! Linear grid point index
        integer(int32), intent(in) :: j_e                   ! Shock dimension index

        ! Model structure
        type(model_dimensions), intent(in) :: dims
        type(mcp_solver_config), intent(in) :: mcp
        real(real64), intent(in) :: sm(:,:)                 ! State matrix [(n_het_exo+n_states) × N_sp]
        integer(int32), intent(in) :: state_var(:)          ! Indices of state vars [n_states]

        ! Grid configuration
        type(grid), intent(in) :: pol_grids(:)              ! 1D state grids [n_states]
        integer(int32), intent(in) :: pol_dims(:)           ! Grid dimensions [n_states]

        ! Policy data
        real(real64), intent(in) :: weighted_old_pol(:,:)   ! Mu-weighted policy [N_sp × n_het_endo]

        ! Gray code cache for next-period interpolation
        type(gray_code_cache), intent(inout) :: gc_cache_next

        ! Workspace (outputs written here)
        type(grid_point_workspace), intent(inout) :: gp_ws

        ! MATLAB handles (inout because step_mx may be lazily created)
        type(matlab_mex_handles), intent(inout) :: input_mex

        ! Local variables
        integer(int32) :: l, k, t, kf, z, a(1), lv_idx
        integer(int32), allocatable :: dz(:)
        real(real64) :: beta
        real(real64), allocatable :: dbeta(:)
        real(real64), pointer :: pol_next(:)

        allocate(dbeta(dims%n_states))
        allocate(dz(dims%n_states))

        ! Build yh array: [yh(-1), yh(0), yh(+1)]
        gp_ws%yh = 0.0_real64
        ! Use state_var to place state values in correct positions
        gp_ws%yh(state_var) = sm(dims%n_het_exo+1:dims%n_het_exo+dims%n_states, j)
        ! Current policy guess
        gp_ws%yh(dims%n_het_endo+1:dims%n_het_endo+dims%n_orig) = policy
        ! Call dynamic_het1_set_auxiliary_variables to compute auxiliary variables
        ! This sets yh(0) auxiliaries from the originally declared variables
        ! Currently, the only level that's used is the level 0
        ! As we don't need the MCP multipliers value, we don't need the level 1
        if (mcp%set_auxiliary_variables) then
            lv_idx = 0
            call call_matlab_set_auxiliary_variables(gp_ws%yh, input_mex, l)
        end if
        if (dims%n_mult > 0) then
            ! Set MCP multipliers to zero for the FB function to work correctly
            gp_ws%yh(dims%n_het_endo+mcp%mult_in_het) = 0.0_real64
        end if

        ! Build xh array: current shocks (from grid point)
        gp_ws%xh = sm(1:dims%n_het_exo, j)

        ! Interpolate policy at t+1 for all shock realizations
        ! Extract next-period state values from policy
        do k = 1, dims%n_states
            gp_ws%states_next(k) = policy(state_var(k))
        end do

        ! Find bracketing indices and weights for next-period states
        do k = 1, dims%n_states
            call bracket_linear_weight(pol_grids(k)%nodes, &
                                       pol_dims(k), &
                                       gp_ws%states_next(k:k), 1_int32, &
                                       gp_ws%ind_next(1,k:k), gp_ws%w_next(1,k:k))
            ! Compute inverse grid spacing for derivative
            if (gp_ws%ind_next(1,k) == pol_dims(k)) then
                gp_ws%inv_h(k) = 1.0_real64 / (pol_grids(k)%nodes(pol_dims(k)) - pol_grids(k)%nodes(pol_dims(k)-1))
            else
                gp_ws%inv_h(k) = 1.0_real64 / (pol_grids(k)%nodes(gp_ws%ind_next(1,k)+1) - pol_grids(k)%nodes(gp_ws%ind_next(1,k)))
            end if
        end do

        ! Low-corner linear index in 1,...,N_a (state-space only)
        call compute_linear_indices(gp_ws%ind_next, gc_cache_next%stride_states, a)

        ! Fetch it into 1, ..., N_sp (including shock dimension)
        l = (a(1)-1)*dims%N_e+j_e

        ! Coefficient updates and hard dimensions
        call compute_coefficient_updates(gp_ws%w_next, &
                                         gc_cache_next%r_up, &
                                         gc_cache_next%r_down, &
                                         gc_cache_next%is_hard_one, &
                                         gc_cache_next%is_hard_zero)

        ! Low-corner initial weight
        beta = product(gp_ws%w_next(1,:), mask=(.not. gc_cache_next%is_hard_zero(1,:)))
        ! dbeta is a n_states vector where the index is the states w.r.t which
        ! we differentiate the policy function
        dbeta = -gp_ws%inv_h * beta
        where (.not. gc_cache_next%is_hard_zero(1,:)) dbeta = dbeta / gp_ws%w_next(1,:)

        ! Initialization at the low corner
        gc_cache_next%sk = .false.
        z = count(gc_cache_next%is_hard_zero(1,:))
        pol_next(1:dims%n_het_endo) => gp_ws%yh(2*dims%n_het_endo+1:3*dims%n_het_endo)
        pol_next = 0.0_real64
        if (z==0_int32) pol_next = pol_next + beta * weighted_old_pol(l,:)
        gp_ws%dpol_next = 0.0_real64
        do k=1, dims%n_states
            dz(k) = count(gc_cache_next%is_hard_zero(1,:))
            if (gc_cache_next%is_hard_zero(1,k)) then
                dz(k) = dz(k) - 1_int32
            end if
            if (dz(k) == 0_int32) then
                gp_ws%dpol_next(:,k) = gp_ws%dpol_next(:,k) + dbeta(k) * weighted_old_pol(l,:)
            end if
        end do

        ! Going through the remaining corners
        do t = 1, gc_cache_next%Kcorn - 1
            kf = gc_cache_next%flip_idx(t)
            ! For the expected derivative, the kf-state coefficient changes sign
            dbeta(kf) = -dbeta(kf)
            if (.not. gc_cache_next%sk(kf)) then
                ! lower -> upper on dimension kf
                l = l + gc_cache_next%stride_shocks(kf)
                gc_cache_next%sk(kf) = .true.
                ! If kf dim is hard one, the number of hard-one dims flipped to
                ! upper increases.
                if (gc_cache_next%is_hard_one(1,kf)) then
                    z = z + 1_int32
                    ! For the expected derivative, the number of hard-one dims
                    ! flipped to upper increases only for the dimensions k /= kf
                    do k=1, dims%n_states
                        if (kf /= k) dz(k) = dz(k)+1
                    end do
                else if (gc_cache_next%is_hard_zero(1,kf)) then
                    ! Similar case for hard-zero dims.
                    z = z-1_int32
                    do k=1, dims%n_states
                        if (kf /= k) dz(k) = dz(k)-1_int32
                    end do
                else
                    ! Otherwise the kf dim is soft and we update the beta
                    ! coefficient
                    beta = beta * gc_cache_next%r_up(1,kf)
                    ! For the expected derivative, we update the coefficient for
                    ! dimensions k /= kf
                    do k = 1, dims%n_states
                        if (kf /= k) dbeta(k) = dbeta(k) * gc_cache_next%r_up(1,kf)
                    end do
                end if
            else
                ! upper -> lower on dim kf
                l = l - gc_cache_next%stride_shocks(kf)
                gc_cache_next%sk(kf) = .false.
                ! If kf dim is hard-one, the number of hard dims flipped
                ! to upper decreases.
                if (gc_cache_next%is_hard_one(1,kf)) then
                    z = z - 1_int32
                    ! For the expected derivative, the number of hard-one dims
                    ! flipped to upper decreases only for the dimensions k /= kf
                    do k=1, dims%n_states
                        if (kf /= k) dz(k) = dz(k)-1_int32
                    end do
                else if (gc_cache_next%is_hard_zero(1,kf)) then
                    z = z + 1_int32
                    ! For the expected derivative, the number of hard-zero dims
                    ! flipped to lower increases only for the dimensions k /= kf
                    do k=1, dims%n_states
                        if (kf /= k) dz(k) = dz(k)+1_int32
                    end do
                else
                    ! Otherwise the kf dim is soft and we update the beta
                    ! coefficient
                    beta = beta * gc_cache_next%r_down(1,kf)
                    ! For the expected derivative, we update the coefficient for
                    ! dimensions k /= kf
                    do k = 1, dims%n_states
                        if (kf /= k) dbeta(k) = dbeta(k) * gc_cache_next%r_down(1,kf)
                    end do
                end if
            end if
            ! If the number of hard dims flipped to upper is zero, there is
            ! a contribution
            if (z == 0_int32) pol_next = pol_next + beta * weighted_old_pol(l,:)
            do k=1, dims%n_states
                if (dz(k) == 0_int32) gp_ws%dpol_next(:,k) = gp_ws%dpol_next(:,k) + dbeta(k)*weighted_old_pol(l,:)
            end do
        end do

    end subroutine set_het_input

    !---------------------------------------------------------------------------
    ! Solve household problem at all grid points using Trust Region + Fischer-Burmeister
    !
    ! This solver properly handles bound constraints by:
    ! 1. Eliminating multiplier variables (setting them to 0)
    ! 2. Skipping slackness equations
    ! 3. Applying Fischer-Burmeister transformation to handle bounds
    ! 4. Solving the reduced system with trust_region_solve
    !---------------------------------------------------------------------------
    subroutine solve_household_all_points( &
        dims,             &  ! Model dimensions (value type)
        mcp,              &  ! MCP solver config (bounds, reordering)
        sm,               &  ! State matrix [(n_het_exo+n_states) × N_sp]
        pol_grids,        &  ! 1D state grids [n_states]
        pol_dims,         &  ! Grid dimensions [n_states]
        Mu,               &  ! Transition matrix [N_e × N_e]
        ti_options,       &  ! Trust region solver options
        state_var,        &  ! State variable indices [n_states]
        old_pol,          &  ! Previous policy [n_het_endo × N_sp] (read)
        old_pol_array,    &  ! Previous policy [N_e × N_a × n_het_endo] (read)
        new_pol,          &  ! Updated policy [n_het_endo × N_sp] (write)
        weighted_old_pol, &  ! Mu-weighted policy [N_sp × n_het_endo] (write)
        gc_cache_next,    &  ! Gray code cache for interpolation
        gp_ws,            &  ! Grid point workspace
        mex,              &  ! MATLAB MEX handles
        status)

        ! Model dimensions
        type(model_dimensions), intent(in) :: dims

        ! MCP solver configuration (bounds, equation reordering)
        type(mcp_solver_config), intent(in) :: mcp

        ! State matrix [(n_het_exo+n_states) × N_sp]
        real(real64), intent(in), contiguous :: sm(:,:)

        ! Grid configuration
        type(grid), intent(in) :: pol_grids(:)              ! 1D state grids [n_states]
        integer(int32), intent(in), contiguous :: pol_dims(:)  ! Grid dimensions [n_states]

        ! Transition matrix [N_e × N_e]
        real(real64), intent(in), contiguous :: Mu(:,:)

        ! Trust region solver options
        type(time_iteration_options), intent(in) :: ti_options

        ! State variable indices [n_states]
        integer(int32), intent(in), contiguous :: state_var(:)

        ! Policy arrays
        real(real64), intent(in), contiguous, target :: old_pol(:,:) ! [n_het_endo × N_sp]
        real(real64), intent(inout), contiguous, target :: old_pol_array(:,:,:) ! [n_het_endo × N_e × N_a]
        real(real64), intent(inout), contiguous :: new_pol(:,:)      ! [n_het_endo × N_sp]
        real(real64), intent(inout), contiguous, target :: weighted_old_pol(:,:) ! [N_sp × n_het_endo]

        ! Gray code cache for next-period interpolation
        type(gray_code_cache), intent(inout) :: gc_cache_next

        ! Grid point workspace
        type(grid_point_workspace), intent(inout) :: gp_ws

        ! MATLAB MEX handles
        type(matlab_mex_handles), intent(inout) :: mex

        ! Output status
        integer(int32), intent(out) :: status

        ! Local variables
        integer(int32) :: j_a, j_e, j, i, info
        character(len=256) :: msg
        real(real64), pointer, contiguous :: old_pol_3d_ptr(:,:,:), weighted_old_pol_mat_ptr(:,:), old_pol_mat_ptr(:,:)

        ! Initialize status
        status = 0

        ! Weight the previous-iteration policy function with the Mu matrix
        ! (Useful to compute expectations)
        old_pol_3d_ptr(1:dims%n_het_endo,1:dims%N_e,1:dims%N_a_sp) => old_pol
        old_pol_array = reshape(old_pol_3d_ptr, [dims%N_e, dims%N_a_sp, dims%n_het_endo], order=[3, 1, 2])
        old_pol_mat_ptr(1:dims%N_e,1:(dims%N_a_sp*dims%n_het_endo)) => old_pol_array
        weighted_old_pol_mat_ptr(1:dims%N_e,1:(dims%N_a_sp*dims%n_het_endo))=>weighted_old_pol
        call matmul_add("N", "N", 1.0_real64, Mu, old_pol_mat_ptr, 0.0_real64, weighted_old_pol_mat_ptr)

        ! Loop over all grid points
        do j_a = 1, dims%N_a_sp
            do j_e = 1, dims%N_e
                j = dims%N_e*(j_a-1)+j_e

                ! Extract originally declared variables
                gp_ws%x_orig = old_pol(1:dims%n_orig, j)


                ! Solve using trust_region with FB transformation
                call trust_region_solve(gp_ws%x_orig, fb_callback, info, &
                                        tolx=ti_options%solver_tolx, &
                                        tolf=ti_options%solver_tolf, &
                                        maxiter=ti_options%solver_max_iter, &
                                        factor=ti_options%solver_factor)


                ! Check for errors (info=1 is success, info=-1 means initial
                ! guess is already solution). Otherwise, continue with current
                ! solution (x_orig was updated by trust_region_solve)
                if (info /= 1 .and. info /= -1) then
                    ! If stop_on_error is true, print diagnostics and fail
                    if (ti_options%solver_stop_on_error) then
                        write(msg, '(A,I0,A,I0,A)') &
                            '    Trust region failed at grid point (j_e=', j_e, ', j_a=', j_a, ')'
                        call mexPrintf(trim(msg)//NEW_LINE('A'))
                        select case (info)
                        case (0)
                            call mexPrintf('    Reason: Nonlinear system ill-behaved at the initial guess'//NEW_LINE('A'))
                        case (2)
                            call mexPrintf('    Reason: Maximum number of iterations reached'//NEW_LINE('A'))
                        case (3)
                            call mexPrintf('    Reason: Spurious convergence (trust region radius too small)'//NEW_LINE('A'))
                        case (4)
                            call mexPrintf('    Reason: Iteration not making good progress'//NEW_LINE('A'))
                        case (5)
                            call mexPrintf('    Reason: Tolerance too small, no further improvement possible'//NEW_LINE('A'))
                        end select
                        ! Print diagnostic information
                        call mexPrintf('    Initial guess (x_orig):'//NEW_LINE('A'))
                        do i = 1, dims%n_orig
                            write(msg, '(A,I0,A,ES15.6)') '      x(', i, ') = ', gp_ws%x_orig(i)
                            call mexPrintf(trim(msg)//NEW_LINE('A'))
                        end do
                        ! Also print grid point state values
                        call mexPrintf('    Grid point states:'//NEW_LINE('A'))
                        do i = 1, dims%n_states
                            write(msg, '(A,I0,A,ES15.6)') '      state(', i, ') = ', &
                                sm(dims%n_het_exo+i, j)
                            call mexPrintf(trim(msg)//NEW_LINE('A'))
                        end do
                        call mexPrintf('    Shocks:'//NEW_LINE('A'))
                        do i = 1, dims%n_het_exo
                            write(msg, '(A,I0,A,ES15.6)') '      shock(', i, ') = ', &
                                sm(i, j)
                            call mexPrintf(trim(msg)//NEW_LINE('A'))
                        end do
                        ! Print bounds info
                        call mexPrintf('    Bounds:'//NEW_LINE('A'))
                        do i = 1, dims%n_orig
                            write(msg, '(A,I0,A,ES12.4,A,ES12.4,A,I0)') '      var ', i, &
                                ': [', mcp%lower_bounds(i), &
                                ', ', mcp%upper_bounds(i), &
                                '], type=', mcp%bound_type(i)
                            call mexPrintf(trim(msg)//NEW_LINE('A'))
                        end do
                        ! Evaluate residuals at current point to show their values
                        call mexPrintf('    Evaluating residuals at failed point...'//NEW_LINE('A'))
                        call set_het_input(gp_ws%x_orig, j, j_e, dims, mcp, &
                                           sm, state_var, &
                                           pol_grids, pol_dims, &
                                           weighted_old_pol, gc_cache_next, &
                                           gp_ws, mex)
                        call call_matlab_het_resid_tensor(gp_ws%resid, mex)
                        call mexPrintf('    Full residuals (before extraction):'//NEW_LINE('A'))
                        do i = 1, dims%n_het_endo
                            write(msg, '(A,I0,A,ES15.6)') '      resid(', i, ') = ', gp_ws%resid(i)
                            call mexPrintf(trim(msg)//NEW_LINE('A'))
                        end do
                        status = 1
                        return
                    end if

                end if

                ! Store solution for originally declared variables
                ! We use yh as temporary storage
                gp_ws%yh(dims%n_het_endo+1:dims%n_het_endo+dims%n_orig) = gp_ws%x_orig
                ! Set the auxiliary variables (loop over topological levels)
                if (mcp%set_auxiliary_variables) then
                    do i = 0, dims%n_aux_levels - 1
                        call call_matlab_set_auxiliary_variables(gp_ws%yh, mex, i)
                    end do
                end if
                ! Copy back into new_pol
                new_pol(:, j) = gp_ws%yh(dims%n_het_endo+1:2*dims%n_het_endo)
                ! Set MCP multipliers to zero for the FB function to work correctly
                if (dims%n_mult > 0) new_pol(mcp%mult_in_het, j) = 0.0_real64
            end do
        end do

    contains
        ! Callback for trust_region_solve: computes FB-transformed residual and Jacobian
        subroutine fb_callback(x1, fvec, fjac)
            real(real64), dimension(:), intent(in) :: x1
            real(real64), dimension(size(x1)), intent(out) :: fvec
            real(real64), dimension(size(x1),size(x1)), intent(out), optional :: fjac

            integer(int32) :: ii, kk, eq_idx

            ! Set up heterogeneous input (yh, xh) and compute expected future policy
            ! yh(+1) = E[yh(+1)] and Ex = ∂E[yh(+1)]/∂state (for chain rule)
            call set_het_input(x1, j, j_e, dims, mcp, &
                               sm, state_var, &
                               pol_grids, pol_dims, &
                               weighted_old_pol, gc_cache_next, &
                               gp_ws, mex)


            ! Evaluate residual
            call call_matlab_het_resid_tensor(gp_ws%resid, mex)

            ! Extract equations for originally declared variables with MCP reordering
            ! Variable ii has its residual from equation mcp_eq_reordering(ii)
            ! This ensures bounds on variable i are paired with the correct equation's residual
            ! NOTE: Only extract n_orig equations, not n_het_endo (auxiliary vars are computed, not solved)
            do ii = 1, dims%n_orig
                eq_idx = mcp%eq_reordering(ii)
                gp_ws%resid_orig(ii) = gp_ws%resid(eq_idx)
            end do


            ! Apply Fischer-Burmeister transformation
            call fb_transform(dims%n_orig, x1, gp_ws%resid_orig, &
                              mcp%lower_bounds, &
                              mcp%upper_bounds, &
                              mcp%bound_type, &
                              fvec)

            ! Compute Jacobian if requested
            if (present(fjac)) then
                ! Evaluate Jacobian at the certainty equivalent point
                ! yh(+1) is already set to Eyh_next above
                ! jac and jac_next are extracted directly with MCP reordering
                call call_matlab_het_g1_tensor(gp_ws%jac, &
                                               gp_ws%jac_next, &
                                               mcp%eq_reordering_inv, &
                                               mex)


                ! Apply chain rule: jac(:, state_var(kk)) += jac_next * dpol_next(:, kk)
                ! Using BLAS DGEMV for each state dimension
                do kk = 1, dims%n_states
                    call matvecmul_add("N", 1.0_real64, gp_ws%jac_next, &
                                       gp_ws%dpol_next(:,kk), 1.0_real64, &
                                       gp_ws%jac(:, state_var(kk)))
                end do


                ! Get the jacobian over originally declared variables
                gp_ws%jac_orig = gp_ws%jac(1:dims%n_orig,1:dims%n_orig)


                ! Apply FB transformation to Jacobian
                call fb_transform_jacobian(dims%n_orig, x1, gp_ws%resid_orig, gp_ws%jac_orig, &
                                           mcp%lower_bounds, &
                                           mcp%upper_bounds, &
                                           mcp%bound_type, &
                                           fjac)
            end if

        end subroutine fb_callback

    end subroutine solve_household_all_points

    !---------------------------------------------------------------------------
    ! MATLAB MEX Fallback: Call dynamic_het1_resid via mexCallMATLAB
    ! Used when C DLL is not available
    ! Signature matches dynamic_het1_resid_mex.c with 7 inputs (no T management)
    !---------------------------------------------------------------------------
    subroutine call_matlab_het_resid_tensor(resid, input_mex)
        type(matlab_mex_handles), intent(inout) :: input_mex
        real(real64), contiguous, intent(inout) :: resid(:)

        ! MATLAB interface variables
        type(c_ptr), dimension(7) :: prhs
        type(c_ptr), dimension(1) :: plhs
        integer(C_INT) :: retval
        real(real64), pointer, contiguous :: resid_ptr(:)

        ! Build input array for mexCallMATLAB (7 inputs, no T management)
        prhs(1) = input_mex%y_mx
        prhs(2) = input_mex%x_mx
        prhs(3) = input_mex%params_mx
        prhs(4) = input_mex%y_mx         ! steady state (same as y for SS computation)
        prhs(5) = input_mex%yh_mx
        prhs(6) = input_mex%xh_mx
        prhs(7) = input_mex%params_mx    ! paramsh (same as params)

        ! Call MATLAB MEX function: resid = model_name.dynamic_het1_resid(y, x, params, ss, yh, xh, paramsh)
        retval = mexCallMATLAB(1_C_INT, plhs, 7_C_INT, prhs, input_mex%het_resid)
        if (retval /= 0) &
             call mexErrMsgTxt("MATLAB fallback: Failed to call " // input_mex%het_resid)

        ! Extract residual from MATLAB output
        resid_ptr(1:size(resid)) => mxGetDoubles(plhs(1))
        resid = resid_ptr

    end subroutine call_matlab_het_resid_tensor

    !---------------------------------------------------------------------------
    ! MATLAB MEX Fallback: Call dynamic_het1_g1 via mexCallMATLAB
    ! Used when C DLL is not available
    ! Signature matches dynamic_het1_g1_mex.c with 10 inputs (no T management)
    ! Returns sparse Jacobian values directly
    !---------------------------------------------------------------------------
    subroutine call_matlab_het_g1_tensor(jac, jac_next, eq_reordering_inv, input_mex)
        type(matlab_mex_handles), intent(inout) :: input_mex
        real(real64), contiguous, intent(inout) :: jac(:,:), jac_next(:,:)
        integer(int32), contiguous, intent(in) :: eq_reordering_inv(:)

        ! MATLAB interface variables
        type(c_ptr), dimension(10) :: prhs
        type(c_ptr), dimension(1) :: plhs
        type(c_ptr) :: g1_mx
        integer(C_INT) :: retval
        real(real64), pointer, contiguous :: pr(:)
        integer(mwSize), pointer, contiguous :: rowval(:), colptr(:)
        integer(mwSize) :: sparse_row, l, dynamic_g1_nnz
        integer(int32) :: ii, jj, n_het_endo

        n_het_endo = size(eq_reordering_inv)

        ! Build input array for mexCallMATLAB (10 inputs, no T management)
        ! Signature: (y, x, params, ss, yh, xh, paramsh, sparse_rowval, sparse_colval, sparse_colptr)
        prhs(1) = input_mex%y_mx
        prhs(2) = input_mex%x_mx
        prhs(3) = input_mex%params_mx
        prhs(4) = input_mex%y_mx
        prhs(5) = input_mex%yh_mx
        prhs(6) = input_mex%xh_mx
        prhs(7) = input_mex%params_mx
        prhs(8) = input_mex%rowval_mx
        prhs(9) = input_mex%colval_mx
        prhs(10) = input_mex%colptr_mx

        ! Call MATLAB MEX function: g1 = model_name.dynamic_het1_g1(...)
        retval = mexCallMATLAB(1_C_INT, plhs, 10_C_INT, prhs, input_mex%het_jac)
        if (retval /= 0) &
             call mexErrMsgTxt("MATLAB fallback: Failed to call " // input_mex%het_jac)
        g1_mx = plhs(1)

        ! Check that output is sparse
        if (.not. mxIsSparse(g1_mx)) &
             call mexErrMsgTxt("MATLAB fallback: " // input_mex%het_jac // " must return sparse matrix")

        ! Extract sparse matrix data
        colptr(1:input_mex%dynamic_g1_ncols) => mxGetJc(g1_mx)
        dynamic_g1_nnz = colptr(input_mex%dynamic_g1_ncols)
        rowval(1:dynamic_g1_nnz) => mxGetIr(g1_mx)
        pr(1:dynamic_g1_nnz) => mxGetDoubles(g1_mx)

        ! Extract yh(0) Jacobian directly into jac [n_orig x n_orig]
        ! Only extract rows eq_reordering(1:n_orig) and columns 1:n_orig
        ! Columns (n_het_endo+1) to (2*n_het_endo) of full Jacobian correspond to yh(0)
        jac = 0.0_real64
        do jj = 1, n_het_endo
            ! Column jj of yh(0) submatrix corresponds to column (n_het_endo + jj) of full Jacobian
            do l = colptr(n_het_endo + jj) + 1, colptr(n_het_endo + jj + 1)
                sparse_row = int(rowval(l), int32) + 1  ! Convert 0-based to 1-based
                ! Use inverse permutation for O(1) lookup
                ii = eq_reordering_inv(sparse_row)
                jac(ii, jj) = pr(l)
            end do
        end do

        ! Extract yh(+1) Jacobian directly into jac_next [n_orig x n_orig]
        ! Only extract rows eq_reordering(1:n_orig) and columns 1:n_orig
        ! Columns (2*n_het_endo+1) to (3*n_het_endo) of full Jacobian correspond to yh(+1)
        jac_next = 0.0_real64
        do jj = 1, n_het_endo
            ! Column jj of yh(+1) submatrix corresponds to column (2*n_het_endo + jj) of full Jacobian
            do l = colptr(2*n_het_endo + jj) + 1, colptr(2*n_het_endo + jj + 1)
                sparse_row = int(rowval(l), int32) + 1  ! Convert 0-based to 1-based
                ! Use inverse permutation for O(1) lookup
                ii = eq_reordering_inv(sparse_row)
                jac_next(ii, jj) = pr(l)
            end do
        end do

    end subroutine call_matlab_het_g1_tensor

    !---------------------------------------------------------------------------
    ! MATLAB MEX Fallback: Call dynamic_het1_set_auxiliary_variables
    ! Computes auxiliary variables from originally declared variables
    ! Signature: yh = set_aux_fn(y, x, params, ss, yh, xh, paramsh)
    !---------------------------------------------------------------------------
    subroutine call_matlab_set_auxiliary_variables(yh, input_mex, step)
        real(real64), dimension(:), intent(inout) :: yh
        type(matlab_mex_handles), intent(inout) :: input_mex
        integer(int32), intent(in) :: step

        ! MATLAB interface variables
        type(c_ptr), dimension(8) :: prhs
        type(c_ptr), dimension(1) :: plhs
        integer(C_INT) :: retval

        ! Pointer to the results
        real(real64), pointer, contiguous :: yh_ptr(:)

        ! Create mxArray for step parameter
        input_mex%step_mx = mxCreateDoubleScalar(real(step, c_double))

        ! Build input array for mexCallMATLAB (8 inputs)
        ! Signature: (y, x, params, ss, yh, xh, paramsh, step)
        prhs(1) = input_mex%y_mx
        prhs(2) = input_mex%x_mx
        prhs(3) = input_mex%params_mx
        prhs(4) = input_mex%y_mx         ! steady state (same as y for SS computation)
        prhs(5) = input_mex%yh_mx
        prhs(6) = input_mex%xh_mx
        prhs(7) = input_mex%params_mx    ! paramsh (same as params)
        prhs(8) = input_mex%step_mx

        ! Call MATLAB MEX function: yh = model_name.dynamic_het1_set_auxiliary_variables(...)
        retval = mexCallMATLAB(1_C_INT, plhs, 8_C_INT, prhs, input_mex%het_aux)
        if (retval /= 0) &
             call mexErrMsgTxt("MATLAB fallback: Failed to call " // input_mex%het_aux)

        ! Extract updated yh from MATLAB output (auxiliary variables have been set)
        yh_ptr(1:size(yh)) => mxGetDoubles(plhs(1))
        yh = yh_ptr

    end subroutine call_matlab_set_auxiliary_variables

end module time_iteration_tensor
