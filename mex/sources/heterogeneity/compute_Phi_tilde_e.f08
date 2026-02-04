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
! MEX function to build the interpolation matrix Φ̃ₑ for heterogeneous-agent models.
!
! This routine builds sparse triplets (I, J, V) for the interpolation matrix that
! combines multi-dimensional linear interpolation with shock transition mixing.
! The algorithm uses Gray code traversal for efficient hypercube corner enumeration.
!
! MATLAB Syntax:
!   [I, J, V] = compute_Phi_tilde_e(pol_ind, pol_w, pol_dims, Mu)
!
! Inputs:
!   pol_ind   [int32, N_sp × n]    : Lower bracket indices for each state dimension
!   pol_w     [double, N_sp × n]   : Interpolation weights [0,1] for lower corner
!   pol_dims  [int32, n × 1]       : Grid dimensions for each state variable
!   Mu        [double, N_e × N_e]  : Shock transition matrix
!
! Outputs:
!   I         [int32, TOT × 1]     : Row indices (1-based)
!   J         [int32, TOT × 1]     : Column indices (1-based)
!   V         [double, TOT × 1]    : Non-zero values
!
! where TOT = 2^n × N_e × N_sp (maximum number of triplets)

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
    use iso_fortran_env, only: real64, int32
    use iso_c_binding, only: c_ptr, c_int
    use matlab_mex
    use gray_code, only: generate_flip_indices, compute_strides, compute_linear_indices, &
                         compute_coefficient_updates
    implicit none (type, external)

    ! MEX interface
    type(c_ptr), dimension(*), intent(in) :: prhs
    type(c_ptr), dimension(*), intent(out) :: plhs
    integer(c_int), intent(in), value :: nlhs, nrhs

    ! Input MATLAB arrays
    type(c_ptr) :: pol_ind_mx, pol_w_mx, pol_dims_mx, Mu_mx

    ! Dimensions
    integer(int32) :: n              ! Number of state dimensions
    integer(int32) :: N_e            ! Number of shock grid points
    integer(int32) :: N_a            ! Product of pol_dims
    integer(int32) :: N_sp           ! N_e × N_a
    integer(int32) :: Kcorn          ! 2^n (number of hypercube corners)
    integer(int32) :: TOT            ! Total triplets: Kcorn × N_e × N_sp

    ! Input arrays (pointers to MATLAB data)
    integer(int32), pointer, contiguous :: pol_ind(:,:) => null()  ! [N_sp × n]
    real(real64), pointer, contiguous :: pol_w(:,:) => null()      ! [N_sp × n]
    integer(int32), pointer, contiguous :: pol_dims(:) => null()   ! [n]
    real(real64), pointer, contiguous :: Mu(:,:) => null()         ! [N_e × N_e]

    ! Gray code infrastructure
    integer(int32), allocatable :: flip_idx(:)           ! [Kcorn-1]
    integer(int32), allocatable :: stride(:)             ! [n]
    integer(int32), allocatable :: linear_idx(:)         ! [N_sp] linear indices in 1..N_a
    logical, allocatable :: sk(:)                        ! [n] corner state

    ! Coefficient updates
    real(real64), allocatable :: r_up(:,:)               ! [N_sp × n]
    real(real64), allocatable :: r_down(:,:)             ! [N_sp × n]
    logical, allocatable :: is_hard_one(:,:)              ! [N_sp × n] w=1 (lower boundary)
    logical, allocatable :: is_hard_zero(:,:)             ! [N_sp × n] w=0 (upper boundary)
    real(real64), allocatable :: beta(:)                 ! [N_sp] working coefficient
    integer(int32), allocatable :: z(:)                  ! [N_sp] count of hard dims at upper

    ! Work arrays for vectorized emission
    real(real64), allocatable, target :: V_block(:,:)    ! [N_sp × N_e]
    integer(int32), allocatable, target :: I_block(:,:)  ! [N_sp × N_e]
    integer(int32), allocatable, target :: J_block(:,:)  ! [N_sp × N_e]
    integer(int32), allocatable :: one_to_N_sp(:)        ! [N_sp]
    integer(int32), allocatable :: one_to_N_e(:)         ! [N_e]
    integer(int32), pointer, contiguous :: I(:) => null(), J(:) => null()

    ! Output arrays (pointers to MATLAB data)
    integer(int32), pointer, contiguous :: I_out(:) => null()
    integer(int32), pointer, contiguous :: J_out(:) => null()
    real(real64), pointer, contiguous :: V_out(:) => null()

    ! Loop variables
    integer(int32) :: ll, ii, jj, k, t, kf, ptr, block_size

    ! Input validation
    if (nrhs /= 4) call mexErrMsgTxt("Expected 4 inputs: pol_ind, pol_w, pol_dims, Mu")
    if (nlhs /= 3) call mexErrMsgTxt("Expected 3 outputs: I, J, V")

    ! Get input pointers
    pol_ind_mx = prhs(1)
    pol_w_mx = prhs(2)
    pol_dims_mx = prhs(3)
    Mu_mx = prhs(4)

    ! Validate input types
    if (.not. mxIsInt32(pol_ind_mx) .or. mxIsSparse(pol_ind_mx)) then
        call mexErrMsgTxt("1st argument (pol_ind) must be a dense int32 matrix")
    end if
    if (.not. mxIsDouble(pol_w_mx) .or. mxIsSparse(pol_w_mx)) then
        call mexErrMsgTxt("2nd argument (pol_w) must be a dense double matrix")
    end if
    if (.not. mxIsInt32(pol_dims_mx)) then
        call mexErrMsgTxt("3rd argument (pol_dims) must be an int32 vector")
    end if
    if (.not. mxIsDouble(Mu_mx) .or. mxIsSparse(Mu_mx)) then
        call mexErrMsgTxt("4th argument (Mu) must be a dense double matrix")
    end if

    ! Extract dimensions and create Fortran pointers
    N_sp = int(mxGetM(pol_ind_mx), int32)
    n = int(mxGetN(pol_ind_mx), int32)
    N_e = int(mxGetM(Mu_mx), int32)

    ! Verify consistency
    if (int(mxGetM(pol_w_mx), int32) /= N_sp .or. int(mxGetN(pol_w_mx), int32) /= n) then
        call mexErrMsgTxt("pol_w must have same dimensions as pol_ind")
    end if
    if (int(mxGetNumberOfElements(pol_dims_mx), int32) /= n) then
        call mexErrMsgTxt("pol_dims length must match number of columns in pol_ind")
    end if
    if (int(mxGetN(Mu_mx), int32) /= N_e) then
        call mexErrMsgTxt("Mu must be square")
    end if

    ! Create Fortran pointers to MATLAB data
    pol_ind(1:N_sp, 1:n) => mxGetInt32s(pol_ind_mx)
    pol_w(1:N_sp, 1:n) => mxGetDoubles(pol_w_mx)
    pol_dims(1:n) => mxGetInt32s(pol_dims_mx)
    Mu(1:N_e, 1:N_e) => mxGetDoubles(Mu_mx)

    ! Compute derived dimensions
    N_a = product(pol_dims)
    Kcorn = 2**n
    TOT = Kcorn * N_e * N_sp

    ! Verify N_sp = N_e × N_a
    if (N_sp /= N_e * N_a) then
        call mexErrMsgTxt("N_sp must equal N_e × product(pol_dims)")
    end if

    ! Allocate output arrays
    plhs(1) = mxCreateNumericMatrix(int(TOT, mwSize), 1_mwSize, mxINT32_CLASS, mxREAL)
    plhs(2) = mxCreateNumericMatrix(int(TOT, mwSize), 1_mwSize, mxINT32_CLASS, mxREAL)
    plhs(3) = mxCreateDoubleMatrix(int(TOT, mwSize), 1_mwSize, mxREAL)

    I_out(1:TOT) => mxGetInt32s(plhs(1))
    J_out(1:TOT) => mxGetInt32s(plhs(2))
    V_out(1:TOT) => mxGetDoubles(plhs(3))

    ! Allocate work arrays
    allocate(flip_idx(Kcorn - 1))
    allocate(stride(n))
    allocate(linear_idx(N_sp))
    allocate(sk(n))
    allocate(r_up(N_sp, n), r_down(N_sp, n), is_hard_one(N_sp, n), is_hard_zero(N_sp, n))
    allocate(beta(N_sp), z(N_sp))
    allocate(V_block(N_sp, N_e), I_block(N_sp, N_e), J_block(N_sp, N_e))
    allocate(one_to_N_sp(N_sp), one_to_N_e(N_e))
    one_to_N_sp = [(jj, jj=1,N_sp)]
    one_to_N_e = [(jj, jj=1,N_e)]

    ! Generate flip indices and strides
    call generate_flip_indices(n, flip_idx, Kcorn)
    call compute_strides(pol_dims, stride)

    ! Low-corner linear index in 1,...,N_a
    call compute_linear_indices(pol_ind, stride, linear_idx)

    ! Get coefficients
    call compute_coefficient_updates(pol_w, r_up, r_down, is_hard_one, is_hard_zero)
    beta = product(pol_w, dim=2, mask=(.not. is_hard_zero))

    ! Initialize Gray walk state
    sk = .false.
    z = count(is_hard_zero, dim=2)
    stride = stride*N_e

    ! Initialize blocks
    do concurrent (ii=1:N_e, jj=1:N_a, k=1:N_e)
        ll = (jj-1)*N_e+ii
        ! Initial weights for the low corner
        V_block(ll, k) = beta(ll) * Mu(ii, k)
    end do
    do concurrent (k=1:N_e)
        ! Span the destination index into 1,...,N_sp
        I_block(:, k) = (linear_idx - 1) * N_e + k
        ! We do this for all values on 1,...,N_sp
        J_block(:, k) = one_to_N_sp
    end do

    ! Save initial-corner block
    block_size = N_sp * N_e
    ptr = 0
    I(1:block_size) => I_block
    J(1:block_size) => J_block
    I_out(ptr+1:ptr+block_size) = I
    J_out(ptr+1:ptr+block_size) = J
    do concurrent (k=1:N_e)
        V_out(ptr+(k-1)*N_sp+1:ptr+k*N_sp) = merge(V_block(:,k), 0.0_real64, z == 0_int32)
    end do
    ptr = ptr + block_size

    ! Walk remaining 2^n - 1 corners
    do t = 1, Kcorn - 1
        kf = flip_idx(t)
        if (.not. sk(kf)) then
            ! Lower -> upper
            I_block = I_block + stride(kf)
            sk(kf) = .true.
            do concurrent (jj=1:N_sp)
                if (is_hard_one(jj, kf)) then
                    ! For nodes with a hard-one kf dim, the number of hard dims flipped
                    ! to upper increases
                    z(jj) = z(jj) + 1_int32
                else if (is_hard_zero(jj, kf)) then
                    ! For nodes with a hard-zero kf dim, the number of hard dims flipped
                    ! to lower decreases
                    z(jj) = z(jj) - 1_int32
                else
                    ! For nodes with a soft kf dim, there is a contribution iff no
                    ! hard dimension is flipped. We also update the
                    ! coefficient
                    do k = 1, N_e
                        V_block(jj,k) = V_block(jj,k) * r_up(jj, kf)
                    end do
                end if
            end do
        else
            ! Upper -> lower
            I_block = I_block - stride(kf)
            sk(kf) = .false.
            do concurrent (jj=1:N_sp)
                if (is_hard_one(jj, kf)) then
                    ! For nodes with a hard-one kf dim, the number of hard dims flipped
                    ! to upper decreases.
                    z(jj) = z(jj) - 1_int32
                else if (is_hard_zero(jj, kf)) then
                    ! For nodes with a hard-zero kf dim, the number of hard dims flipped
                    ! to lower increases.
                    z(jj) = z(jj) + 1_int32
                else
                    ! For nodes with a soft kf dim, we update the beta coefficient
                    do k = 1, N_e
                        V_block(jj,k) = V_block(jj,k) * r_down(jj, kf)
                    end do
                end if
            end do
        end if
        ! Store the contribution
        do concurrent (k=1:N_e)
            V_out(ptr+(k-1)*N_sp+1:ptr+k*N_sp) = merge(V_block(:,k), 0.0_real64, z == 0_int32)
        end do
        I_out(ptr+1:ptr+block_size) = I
        J_out(ptr+1:ptr+block_size) = J
        ptr = ptr+block_size
    end do

end subroutine mexFunction
