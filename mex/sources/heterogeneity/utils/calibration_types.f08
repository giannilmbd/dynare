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
! Derived types for heterogeneous-agent steady-state computation

module calibration_types
    use iso_fortran_env, only: real64, int32
    use iso_c_binding, only: c_ptr, c_null_ptr, C_FUNPTR, C_DOUBLE, C_INT64_T, c_int, c_associated
    use blas, only: blint

    implicit none (type, external)

    private
    public :: grid, tensor_grid_config, gray_code_cache, &
              time_iteration_options, distribution_options, solver_options, &
              mcp_solver_config, model_dimensions, broyden_workspace, &
              grid_point_workspace, matlab_mex_handles, &
              calibration_input, calibration_output, &
              time_iteration_output, distribution_output, aggregation_output, &
              BOUND_NONE, BOUND_LOWER, BOUND_UPPER, BOUND_BOTH

    ! Bound type enumeration for Fischer-Burmeister transformation
    integer(int32), parameter :: BOUND_NONE = 0   ! lb=-∞, ub=+∞ → Φ = -F (unconstrained)
    integer(int32), parameter :: BOUND_LOWER = 1  ! lb finite, ub=+∞ → Φ = φ(x-lb, -F)
    integer(int32), parameter :: BOUND_UPPER = 2  ! lb=-∞, ub finite → Φ = φ(ub-x, F)
    integer(int32), parameter :: BOUND_BOTH = 3   ! both finite → φ(x-lb, φ(ub-x, -F))

    ! 1D grid for cell array extraction
    type :: grid
        real(real64), pointer, contiguous :: nodes(:) => null()
    end type grid

    ! MCP solver configuration (Fischer-Burmeister transformation)
    ! Contains data arrays for complementarity constraints; dimensions are in model_dimensions
    ! NOTE: Defined before tensor_grid_config because it's used as a component
    type :: mcp_solver_config
        ! Complementarity bounds for original equations [n_orig]
        real(real64), allocatable :: lower_bounds(:)
        real(real64), allocatable :: upper_bounds(:)
        integer(int32), allocatable :: bound_type(:)

        ! MCP equation reordering [n_het_endo]
        integer(int32), pointer, contiguous :: eq_reordering(:) => null()
        ! Inverse of eq_reordering: maps equation index -> variable index [n_het_endo]
        integer(int32), allocatable :: eq_reordering_inv(:)

        ! Indices of original equations (rows to solve with FB) [n_orig]
        integer(int32), pointer, contiguous :: orig_eqs(:) => null()

        ! Indices of MULT_* variables in het endo [n_mult]
        integer(int32), pointer, contiguous :: mult_in_het(:) => null()

        ! Flag: presence of auxiliary het endo variables
        logical :: set_auxiliary_variables = .false.
    end type mcp_solver_config

    ! Tensor grid configuration
    ! Stores precomputed interpolation data from MATLAB
    type :: tensor_grid_config
        ! Interpolation matrices for the policy functions [N_sp × n_states]
        ! Sparse grid
        integer(int32), allocatable :: pol_ind_sp(:,:) ! Bracketing lower indices
        real(real64), allocatable :: pol_w_sp(:,:)     ! Interpolation weights [0,1]
        ! Dense grid
        integer(int32), allocatable :: pol_ind_om(:,:) ! Bracketing lower indices
        real(real64), allocatable :: pol_w_om(:,:)     ! Interpolation weights [0,1]

        ! Interpolation matrices for the dense-to-sparse grids [N_a × n_states]
        integer(int32), pointer, contiguous :: d_ind(:,:) ! Bracketing lower indices
        real(real64), pointer, contiguous :: d_w(:,:)     ! Interpolation weights [0,1]

        ! State grids (array of 1D grids)
        ! Corresponds to MATLAB cell array mat.pol.grids_array
        type(grid), allocatable :: pol_grids(:) ! [n_states]
        type(grid), allocatable :: d_grids(:)   ! [n_states]

        ! Grid dimensions
        ! Corresponds to MATLAB int32 array mat.pol.dims
        integer(int32), pointer, contiguous :: pol_dims(:) => null() ! [n_states] - length of each grid
        integer(int32), pointer, contiguous :: d_dims(:) => null() ! [n_states] - length of each grid

        ! State matrix [(n_e + n_a) × N_sp]
        ! Corresponds to MATLAB matrix mat.pol.sm
        ! Contains the state values (shocks and assets) at each grid point
        real(real64), pointer, contiguous :: sm(:,:) => null()

        ! Transition matrix for shocks [N_e × N_e]
        ! Corresponds to MATLAB matrix mat.Mu
        real(real64), pointer, contiguous :: Mu(:,:) => null()
        real(real64), allocatable :: p(:)

        ! Transpose of Mu for efficient BLAS operations [N_e × N_e]
        real(real64), allocatable :: MuT(:,:)

        ! MCP solver configuration (Fischer-Burmeister transformation)
        ! Used only in time iteration for complementarity constraints
        type(mcp_solver_config) :: mcp
    end type tensor_grid_config

    ! Gray code infrastructure for multi-dimensional interpolation
    ! Precomputed once for efficiency in expectation computations
    type :: gray_code_cache
        ! Gray code sequence for hypercube traversal
        integer(int32), allocatable :: flip_idx(:)      ! [Kcorn-1] - flip sequence
        integer(int32) :: Kcorn                         ! Number of corners: 2^n_states

        ! Grid geometry
        integer(int32), allocatable :: stride_states(:) ! [n_states] - strides for linear indexing
        integer(int32), allocatable :: stride_shocks(:) ! [n_states] - strides for linear indexing
        integer(int32), allocatable :: linear_idx(:)    ! [N_sp] - low-corner linear indices (working copy, modified during traversal)
        integer(int32), allocatable :: corner_idx(:)    ! [N_sp] - low-corner linear indices (initial values, never modified)

        ! Common work arrays
        logical, allocatable :: sk(:)

        ! Work arrays for forward iteration algorithm
        real(real64), allocatable :: r_up(:,:)   ! [N_om × n_states]
        real(real64), allocatable :: r_down(:,:) ! [N_om × n_states]
        real(real64), allocatable :: beta_0(:)   ! [N_om] initial value, never modified
        real(real64), allocatable :: beta(:)     ! [N_om]
        real(real64), allocatable :: inc(:)      ! [N_om]
        logical, allocatable :: is_hard_one(:,:)  ! [N_om × n_states] w=1 (lower boundary)
        logical, allocatable :: is_hard_zero(:,:) ! [N_om × n_states] w=0 (upper boundary)
        integer(int32), allocatable :: z(:)
        real(real64), allocatable :: acc(:)
    end type gray_code_cache

    ! Time iteration options
    type :: time_iteration_options
        integer(int32) :: max_iter                   ! Maximum iterations for time iteration
        real(real64) :: tol                          ! Convergence tolerance for time iteration
        real(real64) :: solver_tolf            ! Tolerance for pointwise nonlinear solver (tolf)
        real(real64) :: solver_factor          ! Initial trust region radius factor
        integer(int32) :: solver_max_iter      ! Max iterations for trust region solver
        real(real64) :: solver_tolx            ! Tolerance for spurious convergence check (tolx)
        logical :: solver_stop_on_error        ! If true, stop on solver failure; if false, continue with current solution
        real(real64) :: learning_rate                ! Policy update dampening factor (1 = no dampening)
        integer(int32) :: early_stopping             ! Stop after N consecutive diff increases (0 = disabled)
    end type time_iteration_options

    ! Forward iteration (distribution) options
    type :: distribution_options
        integer(int32) :: max_iter                   ! Maximum iterations
        real(real64) :: tol                         ! Convergence tolerance
        integer(int32) :: check_every               ! Check convergence every N iterations
    end type distribution_options

    ! Equilibrium solver options (generic for MINPACK or KINSOL)
    type :: solver_options
        real(real64) :: tol        ! Tolerance for equilibrium
        integer(int32) :: max_iter ! Max iterations
    end type solver_options

    ! Model dimensions (consolidated from scattered fields)
    ! Contains all dimension information used across the heterogeneous-agent solver
    type :: model_dimensions
        ! Grid sizes
        integer(int32) :: N_sp = 0       ! Total number of nodes in the sparse policy grid
        integer(int32) :: N_om = 0       ! Total number of nodes in the dense distribution grid
        integer(int32) :: N_e = 0        ! Total number of nodes of the shocks grids
        integer(int32) :: N_a_sp = 0     ! Nodes in sparse policy states grids (= N_sp / N_e)
        integer(int32) :: N_a_om = 0     ! Nodes in dense distribution states grids

        ! Variable counts
        integer(int32) :: n_het_endo = 0 ! Heterogeneous endogenous vars (full, including aux)
        integer(int32) :: n_het_exo = 0  ! Heterogeneous exogenous vars
        integer(int32) :: n_states = 0   ! Continuous state variables
        integer(int32) :: n_agg_endo = 0 ! Aggregate endogenous vars
        integer(int32) :: n_free_parameters = 0 ! Calibrated parameters

        ! Solver-related counts
        integer(int32) :: n_orig = 0     ! Originally declared het endo vars (solver optimizes over these)
        integer(int32) :: n_mult = 0     ! Number of multiplier variables
        integer(int32) :: n_aux_levels = 0 ! Number of topological levels for aux variable computation

        ! Level structure for aux variables (for time-shifting)
        integer(int32), allocatable :: het_aux_level_sizes(:)  ! Size of each level
        integer(int32), allocatable :: het_aux_level_vars(:)   ! Concatenated level vars (1-based indices)
    end type model_dimensions

    ! Broyden workspace for calibration solver
    ! Stores pre-allocated work arrays for the Broyden quasi-Newton method
    ! Following Python SSJ implementation (sequence_jacobian/utilities/solvers.py)
    type :: broyden_workspace
        integer(int32) :: n_vars = 0           ! Number of variables/equations

        ! Jacobian approximation and work copy for LU factorization
        real(real64), allocatable :: J(:,:)        ! [n x n] Approximate Jacobian
        real(real64), allocatable :: J_copy(:,:)   ! [n x n] Copy for dgesv (overwrites input)

        ! Newton step and residual work arrays
        real(real64), allocatable :: dx(:)         ! [n] Newton step
        real(real64), allocatable :: fvec_new(:)   ! [n] Residuals at trial point
        real(real64), allocatable :: df(:)         ! [n] Residual change (fvec_new - fvec)

        ! LAPACK pivot indices for dgesv
        integer(blint), allocatable :: ipiv(:)     ! [n] Pivot indices
    end type broyden_workspace

    ! Grid point workspace for time iteration solver
    ! Groups all mutable state needed for a single grid point solve
    ! Designed for future parallelization (one workspace per thread)
    type :: grid_point_workspace
        ! Trust region solver workspace (Fischer-Burmeister transformation)
        ! Dimensions n_orig and n_het_endo are in input%dims
        real(real64), allocatable :: x_orig(:)      ! [n_orig] - current solution
        real(real64), allocatable :: resid_orig(:)  ! [n_orig] - extracted residual (after MCP reordering)
        real(real64), allocatable :: resid(:)       ! [n_het_endo] - full residual from MATLAB

        ! Working Jacobian arrays (extracted directly from sparse, then chain rule applied)
        real(real64), allocatable :: jac_orig(:,:)       ! [n_orig x n_orig] - working Jacobian
        real(real64), allocatable :: jac(:,:)       ! [n_het_endo x n_het_endo] - working Jacobian
        real(real64), allocatable :: jac_next(:,:)  ! [n_het_endo x n_het_endo] - yh(+1) Jacobian for chain rule

        ! Heterogeneous variable buffers (for MATLAB callbacks)
        real(real64), pointer, contiguous :: yh(:)          ! [3*n_het_endo]
        real(real64), pointer, contiguous :: xh(:)          ! [n_het_exo]

        ! Next-period interpolation work arrays
        real(real64), allocatable :: states_next(:) ! [n_states]
        real(real64), allocatable :: w_next(:,:)    ! [1 x n_states]
        integer(int32), allocatable :: ind_next(:,:)! [1 x n_states]
        real(real64), allocatable :: dpol_next(:,:) ! [n_het_endo x n_states] - dE[yh(+1)]/d(state)
        real(real64), allocatable :: inv_h(:)       ! [n_states] - inverse grid spacing

        logical :: initialized = .false.
    end type grid_point_workspace

    ! Grouped MATLAB MEX array handles
    ! Separates MATLAB interface pointers from computational data
    type :: matlab_mex_handles
        ! Function names
        character(len=:), allocatable :: het_resid
        character(len=:), allocatable :: het_jac
        character(len=:), allocatable :: het_aux
        character(len=:), allocatable :: agg_resid

        ! Aggregate variables
        type(c_ptr) :: y_mx = c_null_ptr
        type(c_ptr) :: x_mx = c_null_ptr
        type(c_ptr) :: yagg_mx = c_null_ptr
        type(c_ptr) :: params_mx = c_null_ptr

        ! Heterogeneous variables
        type(c_ptr) :: yh_mx = c_null_ptr
        type(c_ptr) :: xh_mx = c_null_ptr
        type(c_ptr) :: paramsh_mx = c_null_ptr

        ! Auxiliary variable computation
        type(c_ptr) :: step_mx = c_null_ptr

        ! Sparse Jacobian structure
        type(c_ptr) :: rowval_mx = c_null_ptr
        type(c_ptr) :: colval_mx = c_null_ptr
        type(c_ptr) :: colptr_mx = c_null_ptr

        ! Total number of Jacobian columns
        integer(int32) :: dynamic_g1_ncols
    end type matlab_mex_handles

    ! Main input structure for steady-state computation
    ! NOTE: Most arrays use POINTER instead of ALLOCATABLE to point directly
    ! to MATLAB memory, avoiding unnecessary copies.
    type :: calibration_input
        type(tensor_grid_config), pointer :: tg_config  => null()

        ! Algorithm options (value types, assembled from MATLAB fields)
        type(time_iteration_options), pointer :: ti_options => null()
        type(distribution_options), pointer :: d_options => null()
        type(solver_options), pointer :: solver_opts => null()

        ! Calibration targets (unknown parameters)
        integer(int32), pointer, contiguous :: free_parameters_ind(:) => null()  ! Indices in params array
        real(real64), pointer, contiguous :: free_parameters_init(:) => null()  ! Initial guesses
        real(real64), pointer, contiguous :: free_parameters_bounds(:,:) => null()  ! Bounds [2 × n_free_parameters] (-Inf/+Inf for unbounded)

        ! Policy functions
        real(real64), allocatable :: old_pol(:,:)   ! Policy from previous time iteration [n_het_endo × N_sp]
        real(real64), allocatable :: old_pol_array(:,:,:) ! Transposed policy from previous time iteration [N_e × N_a × n_het_endo]
        real(real64), allocatable :: weighted_old_pol(:,:) ! Policy from previous time iteration [N_sp × n_het_endo]
        real(real64), allocatable :: new_pol(:,:)   ! Updated policy [n_het_endo × N_sp]

        ! Work arrays for the distribution forward iteration algorithm
        real(real64), pointer, contiguous :: d_next(:)

        ! Known aggregate steady-state values
        real(real64), pointer, contiguous :: y(:) => null()
        real(real64), pointer, contiguous :: x(:) => null()
        real(real64), pointer, contiguous :: yagg(:) => null()

        ! Grid point workspace (mutable state for time iteration)
        type(grid_point_workspace) :: gp_ws

        ! MATLAB MEX array handles (grouped for clarity)
        type(matlab_mex_handles) :: mex

        ! Aggregation indices from indices.Ix
        integer(int32), pointer, contiguous :: Ix_in_agg(:) => null() ! indices.Ix.in_agg
        integer(int32), pointer, contiguous :: Ix_in_het(:) => null() ! indices.Ix.in_het

        ! Indices of calibration target equations
        integer(int32), pointer, contiguous :: target_equations(:) => null() ! indices.target_equations

        ! Model parameters
        real(real64), pointer, contiguous :: params(:) => null()

        ! Names for calibration free_parameters (from indices.free_parameters.names)
        character(len=:), allocatable :: free_parameters_names(:)  ! [n_free_parameters] - parameter names

        ! Names for aggregate equations (from M_.equation_tags)
        character(len=:), allocatable :: equation_names(:)  ! [n_agg_endo] - equation names

        ! Verbosity: 2 = verbose (print all), 1 = silent (print nothing)
        integer(int32) :: forward_verbosity = 2
        integer(int32) :: ti_verbosity = 2
        integer(int32) :: cal_verbosity = 2

        ! Model dimensions (consolidated)
        type(model_dimensions) :: dims

        ! State variable names from M_.heterogeneity(1).state_var
        integer(int32), pointer, contiguous :: state_var(:) => null()  ! [n_states] - indices of state variables

        ! Sparse Jacobian structure from M_.heterogeneity(1)
        ! Used for efficient solver Jacobian evaluation

        ! Broyden solver workspace (for calibration)
        ! Initialized once in compute_steady_state_tensor, reused across all calls
        type(broyden_workspace) :: broyden_ws

        ! Gray code cache for efficient expectation computation
        ! - Forward-iteration algorithm
        type(gray_code_cache), pointer :: gc_cache_fi => null()
        ! - Interpolation of policy functions on the dense grid
        type(gray_code_cache), pointer :: gc_cache_om => null()
        ! - Next-period policy interpolation on the sparse grid
        type(gray_code_cache), pointer :: gc_cache_next => null()
    end type calibration_input

    !---------------------------------------------------------------------------
    ! Time iteration (backward) output
    !---------------------------------------------------------------------------
    type :: time_iteration_output
        ! Policy functions [n_het_endo × N_sp/N_om]
        real(real64), allocatable :: policies_sp(:,:)
        real(real64), pointer, contiguous :: policies_om(:,:)

        ! Convergence diagnostics
        logical :: converged
        integer(int32) :: iterations
        real(real64) :: residual_norm  ! Policy change norm (sup norm)
    end type time_iteration_output

    !---------------------------------------------------------------------------
    ! Distribution (forward iteration) output
    !---------------------------------------------------------------------------
    type :: distribution_output
        ! Stationary distribution [N_sp]
        real(real64), pointer, contiguous :: distribution(:) => null()

        ! Convergence diagnostics
        logical :: converged
        integer(int32) :: iterations
        real(real64) :: residual_norm  ! Distribution change norm
    end type distribution_output

    !---------------------------------------------------------------------------
    ! Aggregation output
    !---------------------------------------------------------------------------
    type :: aggregation_output
        ! Full aggregation of all heterogeneous endogenous variables [n_het_endo]
        ! Ix_j = sum_i D_i * policy_j(state_i)
        real(real64), allocatable :: Ix(:)

        ! Full aggregate residuals (all aggregate equations) [n_agg_endo]
        real(real64), allocatable :: residuals(:)
    end type aggregation_output

    !---------------------------------------------------------------------------
    ! Calibration output (master output structure)
    !---------------------------------------------------------------------------
    type :: calibration_output
        ! Sub-outputs (pointers to instances declared in calling scope)
        type(time_iteration_output), pointer :: ti_output => null()
        type(distribution_output), pointer :: dist_output => null()
        type(aggregation_output), pointer :: agg_output => null()

        ! Calibrated parameters [n_free_parameters]
        real(real64), allocatable :: params(:)

        ! Calibration target equations residuals [n_target_eqs]
        real(real64), allocatable :: residuals(:)

        ! Convergence info (for calibration loop)
        logical :: converged
        integer(int32) :: iterations
        real(real64) :: residual_norm
    end type calibration_output

end module calibration_types
