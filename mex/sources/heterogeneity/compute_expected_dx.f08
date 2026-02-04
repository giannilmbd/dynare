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
!   E = compute_expected_dx(x, ind, w, inv_h, mu, dims)
!
! INPUTS:
!   x       [N_x × N_sp double]  : Policy function values on grid
!   ind     [N_sp × n int32]     : Interpolation indices per state (column per state)
!   w       [N_sp × n double]    : Interpolation weights per state (column per state)
!   inv_h   [N_sp × n double]    : Inverse grid spacing per state (column per state)
!   mu      [N_e × N_e double]   : Transition matrix for shocks
!   dims    [n int32]            : Grid sizes per state
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
    type(c_ptr) :: x_mx, ind_mx, w_mx, invh_mx, mu_mx, dims_mx

    ! Fortran pointers (assigned via matlab_mex wrappers)
    real(real64), pointer, contiguous :: x(:,:), mu(:,:), w(:,:), inv_h(:,:)
    integer(int32), pointer, contiguous :: ind(:,:), m(:)

    ! Output
    integer(mwSize), dimension(3) :: dims
    real(real64), pointer, contiguous :: E(:,:,:)

    ! Useful size variables
    integer(int32) :: n, N_x, N_e, N_sp

    if (nrhs /= 6) call mexErrMsgTxt("Need 6 inputs: x, ind, w, inv_h, mu, dims")
    if (nlhs < 1) call mexErrMsgTxt("Need 1 output")

    ! Assign inputs
    x_mx    = prhs(1)
    ind_mx  = prhs(2)
    w_mx    = prhs(3)
    invh_mx = prhs(4)
    mu_mx   = prhs(5)
    dims_mx = prhs(6)

    ! Useful size variables
    N_x = int(mxGetM(x_mx), int32)
    N_sp = int(mxGetN(x_mx), int32)
    N_e = int(mxGetM(mu_mx), int32)
    n = int(mxGetN(ind_mx), int32)

    ! Convert to Fortran arrays using matlab_mex wrappers
    x(1:N_x, 1:N_sp) => mxGetDoubles(x_mx)
    mu(1:N_e, 1:N_e) => mxGetDoubles(mu_mx)
    ind(1:N_sp, 1:n) => mxGetInt32s(ind_mx)
    w(1:N_sp, 1:n) => mxGetDoubles(w_mx)
    inv_h(1:N_sp, 1:n) => mxGetDoubles(invh_mx)
    m(1:n) => mxGetInt32s(dims_mx)

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