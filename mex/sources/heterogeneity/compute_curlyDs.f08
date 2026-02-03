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
!   D = compute_curlyDs(Ind, W, inv_h, a_hat, om, mu, dims)
!
! INPUTS:
!   Ind     [N_om × n int32]     : Interpolation indices per state (column per state)
!   W       [N_om × n double]    : Interpolation weights per state (column per state)
!   inv_h   [N_om × n double]    : Inverse grid spacing per state (column per state)
!   a_hat   [N_om × S × n double]: Policy function values for state variables
!   om      [N_om × 1 double]    : Distribution weights
!   mu      [N_e × N_e double]   : Transition matrix for shocks
!   dims    [n int32]            : Grid sizes per state
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
    type(c_ptr) :: ind_mx, w_mx, inv_h_mx, a_hat_mx, om_mx, mu_mx, dims_mx

    ! Fortran pointers (assigned via matlab_mex wrappers)
    real(real64), pointer, contiguous :: om(:), mu(:,:), w(:,:), inv_h(:,:)
    integer(int32), pointer, contiguous :: ind(:,:), m(:)
    real(real64), pointer, contiguous :: a_hat(:,:,:)

    ! Output
    real(real64), pointer, contiguous :: D(:,:)

    ! Useful size variables
    integer(int32) :: n, N_e, N_om, S

    if (nrhs /= 7) call mexErrMsgTxt("Need 7 inputs: Ind, W, inv_h, a_hat, om, mu, dims")
    if (nlhs < 1) call mexErrMsgTxt("Need 1 output")

    ! Assign inputs
    ind_mx    = prhs(1)
    w_mx      = prhs(2)
    inv_h_mx  = prhs(3)
    a_hat_mx  = prhs(4)
    om_mx     = prhs(5)
    mu_mx     = prhs(6)
    dims_mx   = prhs(7)

    ! Get dimensions
    N_om = int(mxGetM(a_hat_mx), int32)
    N_e = int(mxGetM(mu_mx), int32)
    n = int(mxGetN(ind_mx), int32)
    S = int(mxGetN(a_hat_mx) / n, int32)

    ! Convert to Fortran arrays using matlab_mex wrappers
    om(1:N_om) => mxGetDoubles(om_mx)
    a_hat(1:N_om, 1:S, 1:n) => mxGetDoubles(a_hat_mx)
    mu(1:N_e, 1:N_e) => mxGetDoubles(mu_mx)
    ind(1:N_om, 1:n) => mxGetInt32s(ind_mx)
    w(1:N_om, 1:n) => mxGetDoubles(w_mx)
    inv_h(1:N_om, 1:n) => mxGetDoubles(inv_h_mx)
    m(1:n) => mxGetInt32s(dims_mx)

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