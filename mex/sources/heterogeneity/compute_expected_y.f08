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
!   E = compute_expected_y(ind, w, y, mu, dims)
!
! INPUTS:
!   ind     [N_om × n int32]   : Interpolation indices per state (column per state)
!   w       [N_om × n double]  : Interpolation weights per state (column per state)
!   y       [N_om × S double]  : Data on tensor grid
!   mu      [N_e × N_e double] : Transition matrix for shocks
!   dims    [n int32]          : Grid sizes per state
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
    type(c_ptr) :: ind_mx, w_mx, y_mx, mu_mx, dims_mx

    ! Fortran pointers (assigned via matlab_mex wrappers)
    real(real64), pointer, contiguous :: y(:,:), mu(:,:), w(:,:)
    integer(int32), pointer, contiguous :: ind(:,:), m(:)

    ! Output
    real(real64), pointer, contiguous :: E(:,:)

    ! Useful size variables
    integer(int32) :: n, N_e, N_om, S

    if (nrhs /= 5) call mexErrMsgTxt("Need 5 inputs: ind, w, y, mu, dims")
    if (nlhs < 1) call mexErrMsgTxt("Need 1 output")

    ! Assign inputs
    ind_mx  = prhs(1)
    w_mx    = prhs(2)
    y_mx    = prhs(3)
    mu_mx   = prhs(4)
    dims_mx = prhs(5)

    ! Useful size variables
    N_om = int(mxGetM(y_mx), int32)
    S = int(mxGetN(y_mx), int32)
    N_e = int(mxGetM(mu_mx), int32)
    n = int(mxGetN(ind_mx), int32)

    ! Convert to Fortran arrays using matlab_mex wrappers
    y(1:N_om, 1:S) => mxGetDoubles(y_mx)
    mu(1:N_e, 1:N_e) => mxGetDoubles(mu_mx)
    ind(1:N_om, 1:n) => mxGetInt32s(ind_mx)
    w(1:N_om, 1:n) => mxGetDoubles(w_mx)
    m(1:n) => mxGetInt32s(dims_mx)

    ! Create output
    plhs(1) = mxCreateDoubleMatrix(int(N_om, mwSize), int(S, mwSize), mxREAL)
    E(1:N_om,1:S) => mxGetDoubles(plhs(1))

    ! Call computational core
    call compute_expected_y(E, y, mu, ind, w, m)
end subroutine mexFunction