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
! Computes transition matrices F for heterogeneous agent dynamics.
!
! MATLAB SYNTAX:
!   transition_matrices = compute_transition_matrices(curlyYs, expectations, curlyDs, T, N_Y, N_Ix)
!
! INPUTS:
!   curlyYs       [T × N_Y × N_Ix double]      : Aggregate impact coefficients
!   expectations  [N_om × N_Ix × (T-1) double] : Expectation matrices (N_om = N_e × N_a)
!   curlyDs       [N_om × T × N_Y double]      : Distribution derivative coefficients
!
! OUTPUTS:
!   transition_matrices [T × T × N_Y × N_Ix]   : Transition matrices F

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
    use iso_c_binding
    use matlab_mex
    use blas
    implicit none (type, external)

    ! MATLAB MEX API
    type(c_ptr), dimension(*), intent(in) :: prhs
    type(c_ptr), dimension(*), intent(out) :: plhs
    integer(c_int), intent(in), value :: nlhs, nrhs

    ! Pointers
    type(c_ptr) :: curlyYs_mx, expectations_mx, curlyDs_mx
    real(real64), pointer, contiguous :: curlyYs(:,:,:), expectations(:,:,:), &
                                         curlyDs(:,:,:), F(:,:,:,:),          &
                                         curlyDs_mat(:,:)

    ! Work arrays
    real(real64), allocatable :: F_mat(:,:)

    ! Dimensions
    integer(int32) :: T, N_Y, N_Ix, N_om, t0, t1, i, j
    integer(mwSize) :: dims(4)

    if (nrhs /= 3) call mexErrMsgTxt("Need 3 inputs: curlyYs, expectations, curlyDs")
    if (nlhs < 1) call mexErrMsgTxt("Need 1 output")

    ! Get inputs
    curlyYs_mx = prhs(1)
    expectations_mx = prhs(2)
    curlyDs_mx = prhs(3)

    ! Get dimensions from inputs
    T = int(mxGetM(curlyYs_mx), int32)
    N_om = int(mxGetM(expectations_mx), int32)
    N_Ix = int(mxGetN(expectations_mx) / (T-1), int32)
    N_Y = int(mxGetN(curlyYs_mx) / N_Ix, int32)

    ! Map input arrays
    curlyYs(1:T, 1:N_Y, 1:N_Ix) => mxGetDoubles(curlyYs_mx)
    expectations(1:N_om, 1:N_Ix, 1:T-1) => mxGetDoubles(expectations_mx)
    curlyDs(1:N_om, 1:T, 1:N_Y) => mxGetDoubles(curlyDs_mx)

    ! Create output array
    dims = [int(T, mwSize), int(T, mwSize), int(N_Y, mwSize), int(N_Ix, mwSize)]
    plhs(1) = mxCreateNumericArray(4_mwSize, dims, mxDOUBLE_CLASS, mxREAL)
    F(1:T, 1:T, 1:N_Y, 1:N_Ix) => mxGetDoubles(plhs(1))

    ! Allocate work arrays
    allocate(F_mat(T-1,T*N_Y))

    ! Put curlyDs in the matrix form
    curlyDs_mat(1:N_om,1:(T*N_Y)) => curlyDs

    ! First row of F: just curlyYs
    do concurrent (j=1:N_Ix, i=1:N_Y, t1=1:T)
        F(1, t1, i, j) = curlyYs(t1, i, j)
    end do

    ! Remaining rows
    do j=1, N_Ix
        call matmul_add("T", "N", 1.0_real64, expectations(:, j, :), curlyDs_mat, 0.0_real64, F_mat)
        do concurrent (i=1:N_Y, t0=2:T, t1=1:T)
            F(t0,t1,i,j) = F_mat(t0-1,t1+T*(i-1))
        end do
    end do

end subroutine mexFunction
