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
! Computes expected derivatives of policy functions for heterogeneous agent dynamics.
!
! MATLAB SYNTAX:
!   E = compute_expected_dx(x, ind, w, inv_h, mu, states, dims)
!
! INPUTS:
!   x       [N_x × N_sp double] : Policy function values on grid
!   ind     [struct]            : Structure with interpolation indices per state (N_sp × 1 per field)
!   w       [struct]            : Structure with interpolation weights per state (N_sp × 1 per field)
!   inv_h   [struct]            : Structure with inverse grid spacing per state (N_sp × 1 per field)
!   mu      [N_e × N_e double]  : Transition matrix for shocks
!   states  [cell array]        : State variable names (fastest dimension first)
!   dims    [struct]            : Structure with grid sizes per state
!
! OUTPUTS:
!   E       [N_x × n × N_sp double] : Expected derivatives (n = number of state variables)
!
! DESCRIPTION:
! Computes expected derivatives of policy functions with respect to state variables,
! incorporating shock transitions via the transition matrix mu. Optimized paths for
! 1-dimensional and n-dimensional state spaces.
!
subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
    use iso_c_binding
    use matlab_mex
    use expectations
    implicit none (type, external)
    integer(c_int), value :: nlhs, nrhs
    type(c_ptr) :: plhs(*), prhs(*)

    ! MATLAB inputs as raw pointers
    type(c_ptr) :: x_mx, ind_struct_mx, w_struct_mx, invh_struct_mx, mu_mx, states_mx, dims_struct_mx
    type(c_ptr) :: ind_field, w_field, invh_field, dim_field
    type(c_ptr) :: fieldname_cell
    character(len=:), allocatable :: state

    ! Fortran pointers (assigned via matlab_mex wrappers)
    real(real64), pointer, contiguous :: x(:,:), mu(:,:), w_col(:), inv_h_col(:)
    integer(int32), pointer, contiguous :: ind_col(:)
    real(real64), allocatable :: w(:,:), inv_h(:,:) 
    integer(int32), allocatable :: m(:), ind(:,:)

    ! Output
    integer(mwSize), dimension(3) :: dims
    real(real64), pointer, contiguous :: E(:,:,:)

    ! Useful size variables
    integer(int32) :: n, N_x, N_e, N_sp, s

    if (nrhs /= 7) call mexErrMsgTxt("Need 7 inputs: x, ind, w, inv_h, mu, states, dims")
    if (nlhs < 1) call mexErrMsgTxt("Need 1 output")

    ! Assign inputs
    x_mx           = prhs(1)
    ind_struct_mx  = prhs(2)
    w_struct_mx    = prhs(3)
    invh_struct_mx = prhs(4)
    mu_mx          = prhs(5)
    states_mx      = prhs(6)  ! cell array of fieldnames
    dims_struct_mx = prhs(7)

    ! Useful size variables
    N_x = int(mxGetM(x_mx), int32)
    N_sp = int(mxGetN(x_mx), int32)
    N_e = int(mxGetM(mu_mx), int32)
    n = int(mxGetNumberOfElements(states_mx), int32)

    ! Convert to Fortran arrays using matlab_mex wrappers
    x(1:N_x,1:N_sp) => mxGetDoubles(x_mx)
    mu(1:N_e,1:N_e) => mxGetDoubles(mu_mx)

    allocate(w(N_sp,n), inv_h(N_sp,n), ind(N_sp,n), m(n))

    ! Loop over fields of the structs (one per state variable)
    do s = 1, n
        ! Get fieldname from cell array
        fieldname_cell = mxGetCell(states_mx, int(s, mwIndex))
        state = mxArrayToString(fieldname_cell)

        ind_field  = mxGetField(ind_struct_mx,  1_mwIndex, state)
        w_field    = mxGetField(w_struct_mx,    1_mwIndex, state)
        invh_field = mxGetField(invh_struct_mx, 1_mwIndex, state)
        dim_field = mxGetField(dims_struct_mx, 1_mwIndex, state)

        ind_col(1:N_sp) => mxGetInt32s(ind_field)
        w_col(1:N_sp) => mxGetDoubles(w_field)
        inv_h_col(1:N_sp) => mxGetDoubles(invh_field)

        m(s) = int(mxGetScalar(dim_field), int32)
        ind(:,s) = ind_col
        w(:,s) = w_col
        inv_h(:,s) = inv_h_col
    end do

    ! Create output
    dims = [int(N_x, mwSize), int(n, mwSize), int(N_sp, mwSize)]
    plhs(1) = mxCreateNumericArray(3_mwSize, dims, mxDOUBLE_CLASS, mxREAL)
    E(1:N_x,1:n,1:N_sp) => mxGetDoubles(plhs(1))

    ! Call computational core
    if (n==1) then 
        call compute_expected_dx_1d(E, x, mu, ind, inv_h)
    else
        call compute_expected_dx_nd(E, x, mu, ind, w, inv_h, m)
    end if
end subroutine mexFunction