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
! Steady-state computation for heterogeneous-agent models (Cases 1 & 2: Tensor Grids)
!
! MATLAB SYNTAX:
!   output = compute_steady_state_tensor(M_fname, equation_names, ...
!       M_params, H_state_var, H_dynamic_g1_sparse_rowval, H_dynamic_g1_sparse_colval, ...
!       H_dynamic_g1_sparse_colptr, options_het, mat, indices)
!
! INPUTS:
!   M_fname                    [string]        : Model name (M_.fname)
!   equation_names             [cell]          : Equation names cell array for display
!   M_params                   [double vector] : Model parameters (M_.params)
!   H_state_var                [int32 vector]  : State variable indices (H_.state_var)
!   H_dynamic_g1_sparse_rowval [int32 vector]  : Jacobian sparse row indices
!   H_dynamic_g1_sparse_colval [int32 vector]  : Jacobian sparse column indices
!   H_dynamic_g1_sparse_colptr [int32 vector]  : Jacobian sparse column pointers
!   options_het                [struct]        : Heterogeneity options structure
!   mat                        [struct]        : Matrices structure containing:
!                                                - mat.Mu [N_e × N_e]: transition matrix
!                                                - mat.y, mat.x: aggregate variables
!                                                - mat.pol.*: policy function data
!                                                - mat.d.*: distribution data
!                                                - mat.unknowns.*: calibration parameters
!   indices                    [struct]        : Index structure containing:
!                                                - indices.Ix.*: aggregation indices
!                                                - indices.target_equations
!                                                - indices.unknowns.*
!                                                - indices.mult.* (optional)
!
! OUTPUTS:
!   output                     [struct]        : Results structure
!
subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
    use iso_c_binding
    use ieee_arithmetic, only: ieee_is_finite
    use matlab_mex
    use expectations
    use calibration_types
    use calibration
    use gray_code
    use time_iteration_tensor
    use pack_output
    use markov
    implicit none (type, external)
    integer(c_int), value :: nlhs, nrhs
    type(c_ptr) :: plhs(*), prhs(*)

    ! MATLAB input pointers
    type(c_ptr) :: M_fname_mx, equation_names_mx, M_params_mx, &
                   H_orig_endo_nbr, H_set_auxiliary_variables, H_dynamic_mcp_equations_ordering, &
                   H_state_var_mx, H_dynamic_g1_sparse_rowval_mx, H_dynamic_g1_sparse_colval_mx, &
                   H_dynamic_g1_sparse_colptr_mx, options_het_mx, mat_mx, &
                   indices_mx
    ! Useful intermediate MATLAB pointers
    type(c_ptr) :: pol_struct, d_struct, struct, cell_array, field, ti_mx, d_mx, cal_mx, Ix_mx

    ! Fortran data structures
    ! Input and output structures
    type(calibration_input) :: input
    type(calibration_output) :: output

    ! Sub-input instances (target attribute for pointer linkage)
    type(tensor_grid_config), target :: tg_config
    type(gray_code_cache), target :: gc_cache_fi, gc_cache_om, gc_cache_next
    type(time_iteration_options), target :: ti_options
    type(distribution_options), target :: d_options
    type(solver_options), target :: solver_opts

    ! Temporary pointers for extraction
    real(real64), pointer, contiguous :: temp_real(:), temp_real_2d(:,:)
    integer(int32), pointer, contiguous :: temp_int_2d(:,:)
    character(len=:), allocatable :: model_name
    ! Temporary pointers for full bounds (before reduction to original equations)
    real(real64), pointer, contiguous :: full_lower_bounds(:), full_upper_bounds(:)
    integer(int32), allocatable :: full_bound_type(:)

    ! Sub-output instances (target attribute for pointer linkage)
    type(time_iteration_output), target :: ti_output
    type(distribution_output), target :: dist_output
    type(aggregation_output), target :: agg_output

    ! Size variables
    integer(int32) :: i, status, n_states, N_sp, N_om, N_a_om, N_e, n_het_endo, n_agg_endo, &
                      total_cols, ntmp, n_params, n_yh, n_xh, n_unknowns, n_Ix, &
                      n_target_eqs, n_y, n_orig
    logical :: flag

    ! Check arguments
    if (nrhs /= 13) call mexErrMsgTxt("Need 13 inputs")
    if (nlhs < 1) call mexErrMsgTxt("Need 1 output")

    ! Assign input pointers
    M_fname_mx = prhs(1)
    equation_names_mx = prhs(2)
    M_params_mx = prhs(3)
    H_orig_endo_nbr = prhs(4)
    H_set_auxiliary_variables = prhs(5)
    H_dynamic_mcp_equations_ordering = prhs(6)
    H_state_var_mx = prhs(7)
    H_dynamic_g1_sparse_rowval_mx = prhs(8)
    H_dynamic_g1_sparse_colval_mx = prhs(9)
    H_dynamic_g1_sparse_colptr_mx = prhs(10)
    options_het_mx = prhs(11)
    mat_mx = prhs(12)
    indices_mx = prhs(13)

    ! ==================================================================
    ! PART 0: Put the direct MATLAB input under Fortran variables
    ! ==================================================================
    ! Extract model name from M_fname_mx for MEX function paths
    ! Heterogeneous functions: +model_name/dynamic_het1_resid, dynamic_het1_g1
    ! Aggregate functions: +model_name/dynamic_resid, dynamic_g1
    if (.not. c_associated(M_fname_mx)) call mexErrMsgTxt("M._fname not found")
    model_name = mxArrayToString(M_fname_mx)
    input%mex%het_aux = trim(model_name) // ".dynamic_het1_set_auxiliary_variables"
    input%mex%het_resid = trim(model_name) // ".dynamic_het1_resid"
    input%mex%het_jac = trim(model_name) // ".dynamic_het1_g1"
    input%mex%agg_resid = trim(model_name) // ".dynamic_resid"

    ! Parameters
    if (.not. c_associated(M_params_mx)) call mexErrMsgTxt("M_.params not found")
    n_params = int(mxGetNumberOfElements(M_params_mx), int32)
    temp_real(1:n_params) => mxGetDoubles(M_params_mx)
    ! Create params_mx and point input%params to it
    input%mex%params_mx = mxCreateDoubleMatrix(int(n_params, mwSize), 1_mwSize, mxREAL)
    input%params(1:n_params) => mxGetDoubles(input%mex%params_mx)
    ! Copy data from M_.params to input%params
    input%params = temp_real

    ! Number of originally declared heterogeneous endogenous variables
    if (.not. c_associated(H_orig_endo_nbr)) call mexErrMsgTxt("H_.orig_endo_nbr not found")
    n_orig = int(mxGetScalar(H_orig_endo_nbr), int32)
    input%dims%n_orig = n_orig

    ! Heterogeneous auxiliary variables flag
    if (.not. c_associated(H_set_auxiliary_variables)) call mexErrMsgTxt("H_.set_auxiliary_variables not found")
    tg_config%mcp%set_auxiliary_variables = logical(mxIsLogicalScalarTrue(H_set_auxiliary_variables))

    ! Heterogeneous equations MCP reordering
    if (.not. c_associated(H_dynamic_mcp_equations_ordering)) call mexErrMsgTxt("H_.dynamic_mcp_equations_ordering not found")
    n_het_endo =  int(mxGetNumberOfElements(H_dynamic_mcp_equations_ordering), int32)
    input%dims%n_het_endo = n_het_endo
    tg_config%mcp%eq_reordering(1:n_het_endo) => mxGetInt32s(H_dynamic_mcp_equations_ordering)

    ! Compute inverse permutation: eq_reordering_inv(eq_idx) = variable_idx
    allocate(tg_config%mcp%eq_reordering_inv(n_het_endo))
    do i = 1, n_het_endo
        tg_config%mcp%eq_reordering_inv(tg_config%mcp%eq_reordering(i)) = i
    end do

    ! Heterogeneous state variables
    if (.not. c_associated(H_state_var_mx)) call mexErrMsgTxt("H_.state_var not found")
    n_states = int(mxGetNumberOfElements(H_state_var_mx), int32)
    input%dims%n_states = n_states
    input%state_var(1:n_states) => mxGetInt32s(H_state_var_mx)

    ! Extract full sparse Jacobian row indices
    if (.not. c_associated(H_dynamic_g1_sparse_rowval_mx)) call mexErrMsgTxt("H_.dynamic_g1_sparse_rowval not found")
    input%mex%rowval_mx = H_dynamic_g1_sparse_rowval_mx

    ! Extract full sparse Jacobian column indices
    if (.not. c_associated(H_dynamic_g1_sparse_colval_mx)) call mexErrMsgTxt("H_.dynamic_g1_colval not found")
    input%mex%colval_mx = H_dynamic_g1_sparse_colval_mx

    ! Extract full sparse Jacobian column pointers
    ! Total columns = 3*n_het_endo + n_het_exo + 3*n_agg_endo + n_agg_exo
    if (.not. c_associated(H_dynamic_g1_sparse_colptr_mx)) call mexErrMsgTxt("H_.dynamic_g1_colptr not found")
    input%mex%colptr_mx = H_dynamic_g1_sparse_colptr_mx
    total_cols = int(mxGetNumberOfElements(input%mex%colptr_mx), int32)  ! colptr has ncols+1 elements
    input%mex%dynamic_g1_ncols = total_cols

    ! ==================================================================
    ! PART 1: Extract options from options_het
    ! ==================================================================

    ! Extract time_iteration options
    ti_mx = mxGetField(options_het_mx, 1_mwIndex, 'time_iteration')
    if (.not. c_associated(ti_mx)) call mexErrMsgTxt("options_het.time_iteration not found")

    field = mxGetField(ti_mx, 1_mwIndex, 'max_iter')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.max_iter not found")
    ti_options%max_iter = int(mxGetScalar(field), int32)

    field = mxGetField(ti_mx, 1_mwIndex, 'tol')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.tol not found")
    ti_options%tol = mxGetScalar(field)

    field = mxGetField(ti_mx, 1_mwIndex, 'solver_tolf')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.solver_tolf not found")
    ti_options%solver_tolf = mxGetScalar(field)

    field = mxGetField(ti_mx, 1_mwIndex, 'solver_tolx')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.solver_tolx not found")
    ti_options%solver_tolx = mxGetScalar(field)

    field = mxGetField(ti_mx, 1_mwIndex, 'solver_factor')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.solver_factor not found")
    ti_options%solver_factor = mxGetScalar(field)

    field = mxGetField(ti_mx, 1_mwIndex, 'solver_max_iter')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.solver_max_iter not found")
    ti_options%solver_max_iter = int(mxGetScalar(field), int32)

    field = mxGetField(ti_mx, 1_mwIndex, 'learning_rate')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.learning_rate not found")
    ti_options%learning_rate = mxGetScalar(field)

    field = mxGetField(ti_mx, 1_mwIndex, 'solver_stop_on_error')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.solver_stop_on_error not found")
    ti_options%solver_stop_on_error = mxGetScalar(field) > 0.5_real64

    field = mxGetField(ti_mx, 1_mwIndex, 'early_stopping')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.early_stopping not found")
    ti_options%early_stopping = int(mxGetScalar(field), int32)

    ! Extract forward distribution computation options
    d_mx = mxGetField(options_het_mx, 1_mwIndex, 'forward')
    if (.not. c_associated(d_mx)) call mexErrMsgTxt("options_het.forward not found")

    field = mxGetField(d_mx, 1_mwIndex, 'max_iter')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.forward.max_iter not found")
    d_options%max_iter = int(mxGetScalar(field), int32)

    field = mxGetField(d_mx, 1_mwIndex, 'tol')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.forward.tol not found")
    d_options%tol = mxGetScalar(field)

    field = mxGetField(d_mx, 1_mwIndex, 'check_every')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.forward.check_every not found")
    d_options%check_every = int(mxGetScalar(field), int32)

    ! Extract verbosity flags for each component
    field = mxGetField(d_mx, 1_mwIndex, 'verbosity')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.forward.verbosity not found")
    input%forward_verbosity = int(mxGetScalar(field), int32)

    field = mxGetField(ti_mx, 1_mwIndex, 'verbosity')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.time_iteration.verbosity not found")
    input%ti_verbosity = int(mxGetScalar(field), int32)

    ! Extract calibration options
    cal_mx = mxGetField(options_het_mx, 1_mwIndex, 'calibration')
    if (.not. c_associated(cal_mx)) call mexErrMsgTxt("options_het.calibration not found")

    field = mxGetField(cal_mx, 1_mwIndex, 'ftol')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.calibration.ftol not found")
    solver_opts%tol = mxGetScalar(field)

    field = mxGetField(cal_mx, 1_mwIndex, 'max_iter')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.calibration.max_iter not found")
    solver_opts%max_iter = int(mxGetScalar(field), int32)

    field = mxGetField(cal_mx, 1_mwIndex, 'verbosity')
    if (.not. c_associated(field)) call mexErrMsgTxt("options_het.calibration.verbosity not found")
    input%cal_verbosity = int(mxGetScalar(field), int32)

    ! Link options structures to input
    input%ti_options => ti_options
    input%d_options => d_options
    input%solver_opts => solver_opts
    input%tg_config => tg_config

    ! ==================================================================
    ! PART 2: Extract the relevant information from mat
    ! ==================================================================

    ! Extract transition matrix mat.Mu [N_e × N_e]
    field = mxGetField(mat_mx, 1_mwIndex, 'Mu')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.Mu not found")
    N_e = int(mxGetM(field), int32)
    input%dims%N_e = N_e
    tg_config%Mu(1:N_e, 1:N_e) => mxGetDoubles(field)

    ! Transpose Mu for efficient BLAS operations (column-major access)
    allocate(tg_config%MuT(N_e, N_e))
    tg_config%MuT = transpose(tg_config%Mu)

    ! Compute the stationnary distributions [N_e x 1]
    allocate(tg_config%p(N_e))
    call compute_stationary(tg_config%Mu, N_e, tg_config%p, d_options%tol, d_options%max_iter, flag)
    if (.not. flag) call mexErrMsgTxt("Stationary distribution of heterogeneous shocks: convergence failed!")

    ! MATLAB arrays y and x
    field = mxGetField(mat_mx, 1_mwIndex, 'y')
    n_y = int(mxGetNumberOfElements(field), int32)
    n_agg_endo = n_y / 3
    input%dims%n_agg_endo = n_agg_endo
    temp_real(1:n_y) => mxGetDoubles(field)
    ! Create y_mx and point input%y to it
    input%mex%y_mx = mxCreateDoubleMatrix(int(n_y, mwSize), 1_mwSize, mxREAL)
    input%y(1:n_y) => mxGetDoubles(input%mex%y_mx)
    ! Copy data from y_mx to input%y
    input%y = temp_real
    input%mex%x_mx = mxGetField(mat_mx, 1_mwIndex, 'x')

    ! Extract mat.pol and mat.d fields
    pol_struct = mxGetField(mat_mx, 1_mwIndex, 'pol')
    if (.not. c_associated(pol_struct)) call mexErrMsgTxt("mat.pol not found")
    d_struct = mxGetField(mat_mx, 1_mwIndex, 'd')
    if (.not. c_associated(d_struct)) call mexErrMsgTxt("mat.d not found")

    ! Dense-to-sparse interpolation indices [N_a_om × n_states]
    field = mxGetField(d_struct, 1_mwIndex, 'ind')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.d.ind not found")
    N_a_om = int(mxGetM(field), int32)
    input%dims%N_a_om = N_a_om
    temp_int_2d(1:N_a_om, 1:n_states) => mxGetInt32s(field)
    allocate(tg_config%d_ind(N_a_om, n_states))
    tg_config%d_ind = temp_int_2d

    ! Policy function values [n_het_endo × N_sp]
    field = mxGetField(pol_struct, 1_mwIndex, 'x_bar')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.pol.x_bar not found")
    N_sp = int(mxGetN(field), int32)
    input%dims%N_sp = N_sp
    input%dims%N_a_sp = N_sp / N_E
    temp_real_2d(1:n_het_endo, 1:N_sp) => mxGetDoubles(field)
    allocate(input%old_pol(n_het_endo, N_sp), input%new_pol(n_het_endo, N_sp))
    input%old_pol = temp_real_2d
    input%new_pol = temp_real_2d
    allocate(input%old_pol_array(1:N_e, 1:input%dims%N_a_sp, 1:n_het_endo), input%weighted_old_pol(1:N_sp,1:n_het_endo))

    ! Dense-to-sparse interpolation weights [N_a_om × n_states]
    field = mxGetField(d_struct, 1_mwIndex, 'w')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.d.w not found")
    temp_real_2d(1:N_a_om, 1:n_states) => mxGetDoubles(field)
    allocate(tg_config%d_w(N_a_om, n_states))
    tg_config%d_w = temp_real_2d

    ! Policy function interpolation indices [N_sp/N_om × n_states]
    N_om = N_a_om*N_e
    input%dims%N_om = N_om
    allocate(tg_config%pol_ind_sp(N_sp, n_states))
    allocate(tg_config%pol_ind_om(N_om, n_states))

    ! Policy function interpolation weights [N_sp/N_om × n_states]
    allocate(tg_config%pol_w_sp(N_sp, n_states))
    allocate(tg_config%pol_w_om(N_om, n_states))

    ! Grid dimensions [n_states]
    field = mxGetField(pol_struct, 1_mwIndex, 'dims')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.pol.dims not found")
    tg_config%pol_dims(1:n_states) => mxGetInt32s(field)
    field = mxGetField(d_struct, 1_mwIndex, 'dims')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.d.dims not found")
    tg_config%d_dims(1:n_states) => mxGetInt32s(field)

    ! Extract state matrix mat.pol.sm [(n_e + n_a) × N_sp]
    field = mxGetField(pol_struct, 1_mwIndex, 'sm')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.pol.sm not found")
    input%dims%n_het_exo = int(mxGetM(field), int32)-n_states
    tg_config%sm(1:input%dims%n_het_exo+n_states, 1:N_sp) => mxGetDoubles(field)

    ! Extract grids from cell array mat.pol.grids_array
    allocate(tg_config%pol_grids(n_states))
    cell_array = mxGetField(pol_struct, 1_mwIndex, 'grids_array')
    if (.not. c_associated(cell_array)) call mexErrMsgTxt("mat.pol.grids_array not found")
    if (.not. mxIsCell(cell_array)) call mexErrMsgTxt("mat.pol.grids_array is not a cell array")
    if (int(mxGetNumberOfElements(cell_array), int32) /= n_states) &
        call mexErrMsgTxt("mat.pol.grids_array has wrong number of elements")
    do i = 1, n_states
        field = mxGetCell(cell_array, int(i, mwIndex))
        if (.not. c_associated(field)) call mexErrMsgTxt("mat.pol.grids_array cell is null")
        tg_config%pol_grids(i)%nodes(1:tg_config%pol_dims(i)) => mxGetDoubles(field)
    end do
    allocate(tg_config%d_grids(n_states))
    cell_array = mxGetField(d_struct, 1_mwIndex, 'grids_array')
    if (.not. c_associated(cell_array)) call mexErrMsgTxt("mat.d.grids_array not found")
    if (.not. mxIsCell(cell_array)) call mexErrMsgTxt("mat.d.grids_array is not a cell array")
    if (int(mxGetNumberOfElements(cell_array), int32) /= n_states) &
        call mexErrMsgTxt("mat.d.grids_array has wrong number of elements")
    do i = 1, n_states
        field = mxGetCell(cell_array, int(i, mwIndex))
        if (.not. c_associated(field)) call mexErrMsgTxt("mat.d.grids_array cell is null")
        tg_config%d_grids(i)%nodes(1:tg_config%d_dims(i)) => mxGetDoubles(field)
    end do

    ! Extract mat.pol.bounds (complementarity conditions) into temporary full arrays
    struct = mxGetField(pol_struct, 1_mwIndex, 'bounds')
    if (.not. c_associated(struct)) call mexErrMsgTxt("mat.pol.bounds not found")

    ! Lower bounds [n_het_endo] - temporary full array
    field = mxGetField(struct, 1_mwIndex, 'lower_bounds')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.pol.bounds.lower_bounds not found")
    full_lower_bounds(1:n_het_endo) => mxGetDoubles(field)

    ! Upper bounds [n_het_endo] - temporary full array
    field = mxGetField(struct, 1_mwIndex, 'upper_bounds')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.pol.bounds.upper_bounds not found")
    full_upper_bounds(1:n_het_endo) => mxGetDoubles(field)

    ! Compute full bound types for Fischer-Burmeister transformation (temporary)
    allocate(full_bound_type(n_orig))
    do i = 1, n_orig
        if (ieee_is_finite(full_lower_bounds(i)) .and. &
            ieee_is_finite(full_upper_bounds(i))) then
            full_bound_type(i) = BOUND_BOTH
        else if (ieee_is_finite(full_lower_bounds(i))) then
            full_bound_type(i) = BOUND_LOWER
        else if (ieee_is_finite(full_upper_bounds(i))) then
            full_bound_type(i) = BOUND_UPPER
        else
            full_bound_type(i) = BOUND_NONE
        end if
    end do

    ! Extract mat.unknowns fields
    struct = mxGetField(mat_mx, 1_mwIndex, 'unknowns')
    if (.not. c_associated(struct)) call mexErrMsgTxt("mat.unknowns not found")

    ! Extract initial guesses for calibration parameters
    field = mxGetField(struct, 1_mwIndex, 'initial_values')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.unknowns.initial_values not found")
    n_unknowns = int(mxGetNumberOfElements(field), int32)
    input%dims%n_unknowns = n_unknowns
    input%unknowns_init(1:n_unknowns) => mxGetDoubles(field)

    ! Extract bounds for calibration parameters
    ! MATLAB stores as [2 × n_unknowns]: row 1 = lower bounds, row 2 = upper bounds
    field = mxGetField(struct, 1_mwIndex, 'bounds')
    if (.not. c_associated(field)) call mexErrMsgTxt("mat.unknowns.bounds not found")
    input%unknowns_bounds(1:2,1:n_unknowns) => mxGetDoubles(field)
    ! Note: bounds use -Inf/+Inf for unbounded dimensions (set in MATLAB)

    ! ==================================================================
    ! PART 3: Extract the relevant information from indices
    ! ==================================================================
    ! Extract indices.Ix.in_endo and indices.Ix.in_x
    Ix_mx = mxGetField(indices_mx, 1_mwIndex, 'Ix')
    if (.not. c_associated(Ix_mx)) call mexErrMsgTxt("indices.Ix not found")

    ! Get Ix.in_agg
    field = mxGetField(Ix_mx, 1_mwIndex, 'in_agg')
    if (.not. c_associated(field)) call mexErrMsgTxt("indices.Ix.in_agg not found")
    n_Ix = int(mxGetNumberOfElements(field), int32)
    input%Ix_in_agg(1:n_Ix) => mxGetInt32s(field)

    ! Get Ix.in_het
    field = mxGetField(Ix_mx, 1_mwIndex, 'in_het')
    if (.not. c_associated(field)) call mexErrMsgTxt("indices.Ix.in_het not found")
    input%Ix_in_het(1:n_Ix) => mxGetInt32s(field)

    ! Extract target_equations
    field = mxGetField(indices_mx, 1_mwIndex, 'target_equations')
    if (.not. c_associated(field)) call mexErrMsgTxt("indices.target_equations not found")
    n_target_eqs = int(mxGetNumberOfElements(field), int32)
    input%target_equations(1:n_target_eqs) => mxGetInt32s(field)

    ! Unknown parameter names from indices.unknowns.names
    struct = mxGetField(indices_mx, 1_mwIndex, 'unknowns')
    if (.not. c_associated(field)) call mexErrMsgTxt("indices.unknowns not found")
    cell_array = mxGetField(struct, 1_mwIndex, 'names')
    if (.not. c_associated(cell_array)) call mexErrMsgTxt("indices.unknowns.names not found")

    ! Allocate array for unknown parameter names
    allocate(character(len=256) :: input%unknowns_names(n_unknowns))

    ! Extract each parameter name from cell array
    do i = 1, n_unknowns
        field = mxGetCell(cell_array, int(i, mwIndex))
        if (.not. c_associated(field)) call mexErrMsgTxt("indices.unknowns.names cell is empty")
        input%unknowns_names(i) = mxArrayToString(field)
    end do

    ! Unknown parameters indices in M_.params
    field = mxGetField(struct, 1_mwIndex, 'ind')
    if (.not. c_associated(field)) call mexErrMsgTxt("indices.unknowns.ind not found")
    input%unknowns_ind(1:n_unknowns) => mxGetInt32s(field)

    ! Extract indices.mult.in_het (multiplier indices for FB solver)
    struct = mxGetField(indices_mx, 1_mwIndex, 'mult')
    if (c_associated(struct)) then
        field = mxGetField(struct, 1_mwIndex, 'in_het')
        if (c_associated(field)) then
            input%dims%n_mult = int(mxGetNumberOfElements(field), int32)
            if (input%dims%n_mult > 0) then
                tg_config%mcp%mult_in_het(1:input%dims%n_mult) => mxGetInt32s(field)
            end if
        end if
    end if

    ! Set the associated multiplier values to zero in old_pol
    if (input%dims%n_mult > 0) input%old_pol(tg_config%mcp%mult_in_het,:) = 0.0_real64

    ! Extract bounds for originally declared variables (1:orig_endo_nbr)
    ! Solver works over these variables; auxiliaries computed by set_auxiliary_variables
    allocate(tg_config%mcp%lower_bounds(n_orig))
    allocate(tg_config%mcp%upper_bounds(n_orig))
    allocate(tg_config%mcp%bound_type(n_orig))
    do i = 1, n_orig
        tg_config%mcp%lower_bounds(i) = full_lower_bounds(i)
        tg_config%mcp%upper_bounds(i) = full_upper_bounds(i)
        tg_config%mcp%bound_type(i) = full_bound_type(i)
    end do

    ! Deallocate temporary full bound type array
    deallocate(full_bound_type)

    ! ==================================================================
    ! PART 4: Precompute Gray code infrastructure for expectations
    ! ==================================================================
    ! Generate Gray code flip sequence
    call generate_flip_indices(n_states, gc_cache_fi%flip_idx, gc_cache_fi%Kcorn)
    call generate_flip_indices(n_states, gc_cache_om%flip_idx, gc_cache_om%Kcorn)
    call generate_flip_indices(n_states, gc_cache_next%flip_idx, gc_cache_next%Kcorn)

    ! Compute strides for linear indexing
    call compute_strides(tg_config%d_dims, gc_cache_fi%stride_states)
    call compute_strides(tg_config%pol_dims, gc_cache_om%stride_states)
    call compute_strides(tg_config%pol_dims, gc_cache_next%stride_states)
    allocate(gc_cache_fi%stride_shocks(n_states), gc_cache_om%stride_shocks(n_states), gc_cache_next%stride_shocks(n_states))
    gc_cache_om%stride_shocks = gc_cache_om%stride_states*N_e
    gc_cache_fi%stride_shocks = gc_cache_fi%stride_states*N_e
    gc_cache_next%stride_shocks = gc_cache_next%stride_states*N_e

    ! Work arrays for the interpolation of the policy functions from the sparse
    ! grid to the dense grid
    allocate(gc_cache_om%r_up(N_a_om, n_states))
    allocate(gc_cache_om%r_down(N_a_om, n_states))
    allocate(gc_cache_om%is_hard_one(N_a_om, n_states))
    allocate(gc_cache_om%is_hard_zero(N_a_om, n_states))
    call compute_coefficient_updates(tg_config%d_w,      &
                                     gc_cache_om%r_up,   &
                                     gc_cache_om%r_down, &
                                     gc_cache_om%is_hard_one, &
                                     gc_cache_om%is_hard_zero)
    allocate(gc_cache_om%beta_0(N_a_om))
    gc_cache_om%beta_0 = product(tg_config%d_w, dim=2, mask=(.not. gc_cache_om%is_hard_zero))
    allocate(gc_cache_om%beta(N_a_om))
    allocate(gc_cache_om%corner_idx(N_a_om))
    call compute_linear_indices(tg_config%d_ind, gc_cache_om%stride_states, gc_cache_om%corner_idx)
    allocate(gc_cache_om%linear_idx(N_a_om))
    allocate(gc_cache_om%z(N_a_om))
    allocate(gc_cache_om%acc(n_het_endo*N_om))

    ! Work arrays for forward iteration algorithm
    allocate(gc_cache_fi%linear_idx(N_om))
    allocate(gc_cache_fi%corner_idx(N_om))
    allocate(gc_cache_fi%r_up(N_om, n_states))
    allocate(gc_cache_fi%r_down(N_om, n_states))
    allocate(gc_cache_fi%is_hard_one(N_om, n_states))
    allocate(gc_cache_fi%is_hard_zero(N_om, n_states))
    allocate(gc_cache_fi%beta_0(N_om))
    allocate(gc_cache_fi%inc(N_om))
    allocate(gc_cache_fi%z(N_om))
    allocate(gc_cache_fi%acc(N_om))

    ! Work arrays for interpolated policy at t+1 (used in set_het_input)
    allocate(input%gp_ws%states_next(n_states))
    allocate(input%gp_ws%ind_next(1, n_states))
    allocate(input%gp_ws%w_next(1, n_states))
    allocate(gc_cache_next%r_up(1, n_states))
    allocate(gc_cache_next%r_down(1, n_states))
    allocate(gc_cache_next%is_hard_one(1, n_states))
    allocate(gc_cache_next%is_hard_zero(1, n_states))
    allocate(gc_cache_next%linear_idx(N_e))
    allocate(input%gp_ws%dpol_next(n_het_endo, n_states))
    allocate(input%gp_ws%inv_h(n_states))

    ! Common work arrays
    allocate(gc_cache_fi%sk(n_states))
    allocate(gc_cache_om%sk(n_states))
    allocate(gc_cache_next%sk(n_states))

    ! Link to input structure
    input%gc_cache_fi => gc_cache_fi
    input%gc_cache_om => gc_cache_om
    input%gc_cache_next => gc_cache_next

    ! ==================================================================
    ! PART 5: Create cached MATLAB arrays
    ! ==================================================================
    n_yh = 3*input%dims%n_het_endo
    n_xh = input%dims%n_het_exo

    ! Create yh_mx and point gp_ws%yh to it
    input%mex%yh_mx = mxCreateDoubleMatrix(int(n_yh, mwSize), 1_mwSize, mxREAL)
    input%gp_ws%yh(1:n_yh) => mxGetDoubles(input%mex%yh_mx)

    ! Create xh_mx and point gp_ws%xh to it
    input%mex%xh_mx = mxCreateDoubleMatrix(int(n_xh, mwSize), 1_mwSize, mxREAL)
    input%gp_ws%xh(1:n_yh) => mxGetDoubles(input%mex%xh_mx)

    ! Create yagg_mx and point input%yagg to it
    input%mex%yagg_mx = mxCreateDoubleMatrix(int(n_Ix, mwSize), 1_mwSize, mxREAL)
    input%yagg(1:n_Ix) => mxGetDoubles(input%mex%yagg_mx)

    ! Create empty paramsh array
    input%mex%paramsh_mx = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)

    ! ==================================================================
    ! PART 6: Initialize work arrays and solver workspace (once for entire computation)
    ! ==================================================================
    ! Allocate work arrays for the distribution
    allocate(input%d_next(N_om))

    ! Initialize grid point workspace for Fischer-Burmeister solver
    ! Solver works over originally declared variables (1:n_orig)
    allocate(input%gp_ws%x_orig(input%dims%n_orig))
    allocate(input%gp_ws%resid_orig(input%dims%n_orig))
    allocate(input%gp_ws%resid(n_het_endo))
    ! Working Jacobian arrays (extracted directly from sparse, chain rule applied to jac)
    allocate(input%gp_ws%jac_orig(input%dims%n_orig, input%dims%n_orig))
    allocate(input%gp_ws%jac(input%dims%n_het_endo, input%dims%n_het_endo))
    allocate(input%gp_ws%jac_next(input%dims%n_het_endo, input%dims%n_het_endo))
    input%gp_ws%initialized = .true.

    ! Initialize Broyden workspace for calibration
    if (n_unknowns > 0) then
        input%broyden_ws%n_vars = n_unknowns
        allocate(input%broyden_ws%J(n_unknowns, n_unknowns))
        allocate(input%broyden_ws%J_copy(n_unknowns, n_unknowns))
        allocate(input%broyden_ws%dx(n_unknowns))
        allocate(input%broyden_ws%fvec_new(n_unknowns))
        allocate(input%broyden_ws%df(n_unknowns))
        allocate(input%broyden_ws%ipiv(n_unknowns))
    end if

    ! Extract equation names from pre-built cell array (passed from MATLAB)
    ! Each cell contains a string: either the equation name or the index as fallback
    if (.not. c_associated(equation_names_mx)) call mexErrMsgTxt("equation_names not found")
    ntmp = int(mxGetNumberOfElements(equation_names_mx), int32)
    allocate(character(len=256) :: input%equation_names(ntmp))
    do i = 1, ntmp
        field = mxGetCell(equation_names_mx, int(i, mwIndex))
        if (c_associated(field)) then
            input%equation_names(i) = mxArrayToString(field)
        else
            input%equation_names(i) = ''
        end if
    end do

    ! ==================================================================
    ! PART 7: Initialize the output structures
    ! ==================================================================
    ! Allocate arrays in sub-output structures
    allocate(ti_output%policies_sp(n_het_endo, N_sp))
    allocate(dist_output%distribution(N_om))
    allocate(agg_output%Ix(n_het_endo))
    allocate(agg_output%residuals(n_agg_endo))

    ! Allocate calibration output arrays
    allocate(output%params(n_unknowns))
    allocate(output%residuals(n_unknowns))

    ! Link pointers from calibration_output to sub-outputs
    output%ti_output => ti_output
    output%dist_output => dist_output
    output%agg_output => agg_output

    ! ==================================================================
    ! PART 8: Call main solver
    ! ==================================================================
    status = calibrate(input, output)
    if (status /= 0) call mexErrMsgTxt("Steady-state computation failed")

    ! ==================================================================
    ! PART 9: Pack output into MATLAB struct
    ! ==================================================================
    call pack_calibration_output(output, input, plhs(1))

    ! ==================================================================
    ! PART 10: Cleanup cached MATLAB arrays
    ! ==================================================================
    if (c_associated(input%mex%params_mx)) call mxDestroyArray(input%mex%params_mx)
    if (c_associated(input%mex%y_mx)) call mxDestroyArray(input%mex%y_mx)
    if (c_associated(input%mex%yagg_mx)) call mxDestroyArray(input%mex%yagg_mx)
    if (c_associated(input%mex%yh_mx)) call mxDestroyArray(input%mex%yh_mx)
    if (c_associated(input%mex%xh_mx)) call mxDestroyArray(input%mex%xh_mx)
    if (c_associated(input%mex%paramsh_mx)) call mxDestroyArray(input%mex%paramsh_mx)

end subroutine mexFunction
