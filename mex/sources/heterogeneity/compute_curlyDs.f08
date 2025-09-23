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
! Computes distribution derivative coefficients for heterogeneous agents.
!
! MATLAB SYNTAX:
!   D = compute_curlyDs(Ind, W, inv_h, a_hat, om, mu, states, sizes)
!
! INPUTS:
!   Ind     [struct]             : Structure with interpolation indices per state (N_om × 1 per field)
!   W       [struct]             : Structure with interpolation weights per state (N_om × 1 per field)
!   inv_h   [struct]             : Structure with inverse grid spacing per state (N_om × 1 per field)
!   a_hat   [N_om × S × n double]: Policy function values for state variables
!   om      [N_om × 1 double]    : Distribution weights
!   mu      [N_e × N_e double]   : Transition matrix for shocks
!   states  [cell array]         : State variable names (fastest dimension first)
!   sizes   [struct]             : Structure with grid sizes per state
!
! OUTPUTS:
!   D       [N_om × S double]    : Distribution derivative coefficients
!
! DESCRIPTION:
! Computes D(i, s) = sum_k sum_j (dB_i^k)_{j_k',j} * prod_{ℓ≠k}(B_i^ℓ)_{j_ℓ',j} * a_hat_{i,s,k} * om_{i,j}
! for distribution derivatives in heterogeneous agent models. Optimized paths for
! 1-dimensional and n-dimensional state spaces.
!
subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
    use iso_c_binding
    use matlab_mex
    use curlyDs
    implicit none (type, external)
    integer(c_int), value :: nlhs, nrhs
    type(c_ptr) :: plhs(*), prhs(*)

    ! MATLAB inputs as raw pointers
    type(c_ptr) :: ind_struct_mx, w_struct_mx, inv_h_struct_mx, a_hat_mx, om_mx, mu_mx, states_mx, sizes_struct_mx
    type(c_ptr) :: ind_field, w_field, inv_h_field, size_field
    type(c_ptr) :: fieldname_cell
    character(len=:), allocatable :: state

    ! Fortran pointers (assigned via matlab_mex wrappers)
    real(real64), pointer, contiguous :: om(:), mu(:,:), w_col(:), inv_h_col(:)
    integer(int32), pointer, contiguous :: ind_col(:)
    real(real64), pointer, contiguous :: a_hat(:,:,:)
    real(real64), allocatable :: w(:,:), inv_h(:,:)
    integer(int32), allocatable :: ind(:,:), m(:)

    ! Output
    real(real64), pointer, contiguous :: D(:,:)

    ! Useful size variables
    integer(int32) :: n, N_e, N_om, S, k

    if (nrhs /= 8) call mexErrMsgTxt("Need 8 inputs: Ind, W, inv_h, a_hat, om, mu, states, sizes")
    if (nlhs < 1) call mexErrMsgTxt("Need 1 output")

    ! Assign inputs
    ind_struct_mx   = prhs(1)
    w_struct_mx     = prhs(2)
    inv_h_struct_mx = prhs(3)
    a_hat_mx        = prhs(4)
    om_mx           = prhs(5)
    mu_mx           = prhs(6)
    states_mx       = prhs(7)
    sizes_struct_mx = prhs(8)

    ! Get dimensions
    N_om = int(mxGetM(a_hat_mx), int32)
    N_e = int(mxGetM(mu_mx), int32)
    n = int(mxGetNumberOfElements(states_mx), int32)
    S = int(mxGetN(a_hat_mx) / n, int32)

    ! Convert distribution weights to Fortran array
    om(1:N_om) => mxGetDoubles(om_mx)
    a_hat(1:N_om,1:S,1:n) => mxGetDoubles(a_hat_mx)
    mu(1:N_e,1:N_e) => mxGetDoubles(mu_mx)

    ! Allocate arrays for struct data
    allocate(w(N_om,n), inv_h(N_om,n), ind(N_om,n), m(n))

    ! Loop over state fields to extract struct data
    do k = 1, n
        ! Get fieldname from cell array
        fieldname_cell = mxGetCell(states_mx, int(k, mwIndex))
        state = mxArrayToString(fieldname_cell)

        ! Extract fields from structures
        ind_field   = mxGetField(ind_struct_mx,   1_mwIndex, state)
        w_field     = mxGetField(w_struct_mx,     1_mwIndex, state)
        inv_h_field = mxGetField(inv_h_struct_mx, 1_mwIndex, state)
        size_field  = mxGetField(sizes_struct_mx, 1_mwIndex, state)

        ! Convert to Fortran pointers/arrays
        ind_col(1:N_om) => mxGetInt32s(ind_field)
        w_col(1:N_om) => mxGetDoubles(w_field)
        inv_h_col(1:N_om) => mxGetDoubles(inv_h_field)

        ! Copy to allocated arrays
        m(k) = int(mxGetScalar(size_field), int32)
        ind(:,k) = ind_col
        w(:,k) = w_col
        inv_h(:,k) = inv_h_col
    end do

    ! Create output array
    plhs(1) = mxCreateDoubleMatrix(int(N_om, mwSize), int(S, mwSize), mxREAL)
    D(1:N_om,1:S) => mxGetDoubles(plhs(1))

    ! Call computational core
    if (n == 1) then
        call compute_curlyDs_1d(D, ind, inv_h, a_hat, om, Mu)
    else
        call compute_curlyDs_nd(D, ind, w, inv_h, a_hat, om, Mu, m)
    end if

end subroutine mexFunction