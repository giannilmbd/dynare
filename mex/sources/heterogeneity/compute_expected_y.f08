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
! Computes expectation operator for heterogeneous agent dynamics.
!
! MATLAB SYNTAX:
!   E = compute_expected_y(ind, w, y, mu, states, dims)
!
! INPUTS:
!   ind     [struct]      : Structure with interpolation indices per state (N_om × 1 per field)
!   w       [struct]      : Structure with interpolation weights per state (N_om × 1 per field)
!   y       [N_om × S double] : Data on tensor grid
!   mu      [N_e × N_e double] : Transition matrix for shocks
!   states  [cell array]  : State variable names (fastest dimension first)
!   dims    [struct]      : Structure with grid sizes per state
!
! OUTPUTS:
!   E       [N_om × S double] : Expectation values
!
! DESCRIPTION:
! Computes ℰ(i, j, s) = sum_{i'} sum_{j'} [ ∏_ℓ (B^ℓ)_{j'_ℓ, j} ] * y_{i', j', s} * μ_{i,i'}
! for expectation operations in heterogeneous agent models.
!
subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
    use iso_c_binding
    use matlab_mex
    use expectations
    implicit none (type, external)
    integer(c_int), value :: nlhs, nrhs
    type(c_ptr) :: plhs(*), prhs(*)

    ! MATLAB inputs as raw pointers
    type(c_ptr) :: y_mx, ind_struct_mx, w_struct_mx, mu_mx, states_mx, dims_struct_mx
    type(c_ptr) :: ind_field, w_field, dim_field
    type(c_ptr) :: fieldname_cell
    character(len=:), allocatable :: state

    ! Fortran pointers (assigned via matlab_mex wrappers)
    real(real64), pointer, contiguous :: y(:,:), mu(:,:), w_slice(:)
    integer(int32), pointer, contiguous :: ind_slice(:)
    real(real64), allocatable :: w(:,:)
    integer(int32), allocatable :: m(:), ind(:,:)

    ! Output
    real(real64), pointer, contiguous :: E(:,:)

    ! Useful size variables
    integer(int32) :: n, N_e, N_om, S, i

    if (nrhs /= 6) call mexErrMsgTxt("Need 6 inputs: ind, w, y, mu, states, dims")
    if (nlhs < 1) call mexErrMsgTxt("Need 1 output")

    ! Assign inputs
    ind_struct_mx  = prhs(1)
    w_struct_mx    = prhs(2)
    y_mx           = prhs(3)
    mu_mx          = prhs(4)
    states_mx      = prhs(5)
    dims_struct_mx = prhs(6)

    ! Useful size variables
    n = int(mxGetNumberOfElements(states_mx), int32)
    N_e = int(mxGetM(mu_mx), int32)
    N_om = int(mxGetM(y_mx), int32)
    S = int(mxGetN(y_mx), int32)

    ! Convert to Fortran arrays using matlab_mex wrappers
    y(1:N_om,1:S) => mxGetDoubles(y_mx)
    mu(1:N_e,1:N_e) => mxGetDoubles(mu_mx)

    allocate(w(N_om,n), ind(N_om,n), m(n))

    ! Loop over fields of the structs (one per state variable)
    do i = 1, n
        ! Get fieldname from cell array
        fieldname_cell = mxGetCell(states_mx, int(i, mwIndex))
        state = mxArrayToString(fieldname_cell)

        ind_field  = mxGetField(ind_struct_mx,  1_mwIndex, state)
        w_field    = mxGetField(w_struct_mx,    1_mwIndex, state)
        dim_field = mxGetField(dims_struct_mx, 1_mwIndex, state)

        ind_slice(1:N_om) => mxGetInt32s(ind_field)
        w_slice(1:N_om) => mxGetDoubles(w_field)

        m(i) = int(mxGetScalar(dim_field), int32)
        ind(:,i) = ind_slice
        w(:,i) = w_slice
    end do

    ! Create output
    plhs(1) = mxCreateDoubleMatrix(int(N_om, mwSize), int(S, mwSize), mxREAL)
    E(1:N_om,1:S) => mxGetDoubles(plhs(1))

    ! Call computational core
    call compute_expected_y(E, y, mu, ind, w, m)
end subroutine mexFunction