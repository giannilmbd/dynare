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
! Computes the Rouwenhorst discretization of an AR(1) process.
!
! MATLAB SYNTAX:
!   [y, p_vec, P_mat] = rouwenhorst(rho, sigma, N, tol, maxiter)
!
! INPUTS:
!   rho     [double scalar]   : Persistence parameter of AR(1) process
!   sigma   [double scalar]   : Unconditional standard deviation
!   N       [int scalar]      : Number of grid points
!   tol     [double scalar]   : Tolerance for stationary distribution convergence
!   maxiter [int scalar]      : Maximum iterations for stationary distribution
!
! OUTPUTS:
!   y       [N × 1 double]    : Grid points for discretized process
!   p_vec   [N × 1 double]    : Stationary distribution probabilities
!   P_mat   [N × N double]    : Transition probability matrix
!
! DESCRIPTION:
! Constructs an N-point Markov chain approximation of an AR(1) process
! using the Rouwenhorst method. Returns the state grid, transition matrix,
! and stationary distribution.
!
subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
    use iso_fortran_env, only: real64, int32
    use matlab_mex
    use markov
    implicit none (type, external)

    ! MATLAB MEX API
    type(c_ptr), dimension(*), intent(in), target :: prhs
    type(c_ptr), dimension(*), intent(out) :: plhs
    integer(c_int), intent(in), value :: nlhs, nrhs
    type(c_ptr) :: rho_mx, sigma_mx, N_mx, tol_mx, maxiter_mx  ! Input arguments

    ! Fortran variables
    real(real64) :: rho, sigma, tol
    integer(int32) :: N, maxiter
    real(real64), pointer, contiguous :: y(:), p_vec(:), P_mat(:, :)

    ! Validate number of inputs/outputs
    if (nrhs /= 5_c_int .or. nlhs /= 3_c_int) &
         call mexErrMsgTxt("Rouwenhorst routine: incorrect number of inputs and/or outputs")

    ! Get input pointers
    rho_mx = prhs(1)
    sigma_mx = prhs(2)
    N_mx = prhs(3)
    tol_mx = prhs(4)
    maxiter_mx = prhs(5)

    ! Validate input types
    if (.not. (mxIsScalar(rho_mx) .and. mxIsDouble(rho_mx))) &
         call mexErrMsgTxt("1st argument (rho) must be a double scalar")
    if (.not. (mxIsScalar(sigma_mx) .and. mxIsDouble(sigma_mx))) &
         call mexErrMsgTxt("2nd argument (sigma) must be a double scalar")
    if (.not. (mxIsScalar(N_mx) .and. mxIsNumeric(N_mx))) &
         call mexErrMsgTxt("3rd argument (N) must be an integer scalar")
    if (.not. (mxIsScalar(tol_mx) .and. mxIsDouble(tol_mx))) &
         call mexErrMsgTxt("4th argument (tol) must be a double scalar")
    if (.not. (mxIsScalar(maxiter_mx) .and. mxIsNumeric(maxiter_mx))) &
         call mexErrMsgTxt("5th argument (maxiter) must be an integer scalar")

    ! Convert MATLAB inputs to Fortran variables
    rho = mxGetScalar(rho_mx)
    sigma = mxGetScalar(sigma_mx)
    N = int(mxGetScalar(N_mx))
    tol = mxGetScalar(tol_mx)
    maxiter = int(mxGetScalar(maxiter_mx))

    ! Allocate MATLAB output variables
    plhs(1) = mxCreateDoubleMatrix(int(N, mwSize), 1_mwSize, mxREAL)  ! y (Nx1)
    plhs(2) = mxCreateDoubleMatrix(int(N, mwSize), 1_mwSize, mxREAL)  ! p_vec (Nx1)
    plhs(3) = mxCreateDoubleMatrix(int(N, mwSize), int(N, mwSize), mxREAL)  ! P_mat (NxN)

    ! Copy results from Fortran arrays to MATLAB output
    y(1:N) => mxGetDoubles(plhs(1))
    p_vec(1:N) => mxGetDoubles(plhs(2))
    P_mat(1:N,1:N) => mxGetDoubles(plhs(3))

    ! Call the Fortran function
    call markov_rouwenhorst(rho, sigma, N, tol, maxiter, y, p_vec, P_mat)

end subroutine mexFunction
