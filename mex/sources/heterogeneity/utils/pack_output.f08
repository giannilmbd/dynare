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
! Pack calibration output into MATLAB struct for return from MEX function

module pack_output
    use iso_c_binding
    use matlab_mex
    use calibration_types
    implicit none (type, external)

contains

    ! Pack calibration_output into MATLAB struct
    !
    ! Arguments:
    !   output     [in]  : Fortran calibration output structure
    !   input      [in]  : Fortran calibration input (for metadata)
    !   output_mx  [out] : MATLAB struct pointer
    !
    subroutine pack_calibration_output(output, input, output_mx)
        type(calibration_output), intent(in) :: output
        type(calibration_input), intent(in) :: input
        type(c_ptr), intent(out) :: output_mx

        ! MATLAB field names
        character(len=17), dimension(8) :: field_names = [ &
            character(len=17) :: 'converged', 'iterations', 'residual_norm', 'params', &
            'param_names', 'time_iteration', 'distribution', 'aggregates']

        ! Sub-struct field names
        character(len=13), dimension(4) :: ti_fields = [ &
            character(len=13) :: 'converged', 'iterations', 'residual_norm', 'policies']
        character(len=13), dimension(4) :: dist_fields = [ &
            character(len=13) :: 'converged', 'iterations', 'residual_norm', 'hist']
        character(len=9), dimension(2) :: agg_fields = [ &
            character(len=9) :: 'Ix', 'residuals']

        ! MATLAB pointers
        type(c_ptr) :: ti_struct, dist_struct, agg_struct
        type(c_ptr) :: field_ptr, names_cell
        type(c_ptr) :: str_ptr
        real(real64), pointer, contiguous :: temp_real(:), temp_real_2d(:,:)
        integer(int32) :: i, n_het_endo, N_sp, N_om, n_agg_endo, n_unknowns
        character(len=256), target :: str_buffer

        ! Extract dimensions
        n_unknowns = input%dims%n_unknowns
        n_het_endo = input%dims%n_het_endo
        N_sp = input%dims%N_sp
        N_om = input%dims%N_om
        n_agg_endo = input%dims%n_agg_endo

        ! Create main output struct using the wrapper (it handles the conversion)
        output_mx = mxCreateStructMatrix(1_mwSize, 1_mwSize, field_names)
        if (.not. c_associated(output_mx)) call mexErrMsgTxt("Failed to create output struct")

        ! ================================================================
        ! PART 1: Top-level convergence fields
        ! ================================================================

        ! converged (logical) - convert logical(4) to logical(mxLogical)
        field_ptr = mxCreateLogicalScalar(merge(.true._mxLogical, .false._mxLogical, output%converged))
        call mxSetField(output_mx, 1_mwIndex, 'converged', field_ptr)

        ! iterations (int32)
        field_ptr = mxCreateDoubleScalar(real(output%iterations, C_DOUBLE))
        call mxSetField(output_mx, 1_mwIndex, 'iterations', field_ptr)

        ! residual_norm (double)
        field_ptr = mxCreateDoubleScalar(output%residual_norm)
        call mxSetField(output_mx, 1_mwIndex, 'residual_norm', field_ptr)

        ! ================================================================
        ! PART 2: Calibrated parameters
        ! ================================================================

        ! params [n_unknowns × 1]
        if (n_unknowns > 0 .and. allocated(output%params)) then
            field_ptr = mxCreateDoubleMatrix(int(n_unknowns, mwSize), 1_mwSize, mxREAL)
            temp_real(1:n_unknowns) => mxGetDoubles(field_ptr)
            temp_real = output%params
        else
            field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
        end if
        call mxSetField(output_mx, 1_mwIndex, 'params', field_ptr)

        ! param_names {n_unknowns × 1} cell array of strings
        if (n_unknowns > 0) then
            names_cell = mxCreateCellMatrix(int(n_unknowns, mwSize), 1_mwSize)
            do i = 1, n_unknowns
                ! Convert Fortran string to C string (null-terminated)
                str_buffer = trim(input%unknowns_names(i)) // c_null_char
                str_ptr = mxCreateString(str_buffer)
                call mxSetCell(names_cell, int(i, mwIndex), str_ptr)
            end do
            call mxSetField(output_mx, 1_mwIndex, 'param_names', names_cell)
        else
            names_cell = mxCreateCellMatrix(0_mwSize, 0_mwSize)
            call mxSetField(output_mx, 1_mwIndex, 'param_names', names_cell)
        end if

        ! ================================================================
        ! PART 3: Time iteration sub-struct
        ! ================================================================

        ti_struct = mxCreateStructMatrix(1_mwSize, 1_mwSize, ti_fields)
        if (.not. c_associated(ti_struct)) call mexErrMsgTxt("Failed to create time_iteration struct")

        if (associated(output%ti_output)) then
            ! converged
            field_ptr = mxCreateLogicalScalar(merge(.true._mxLogical, .false._mxLogical, output%ti_output%converged))
            call mxSetField(ti_struct, 1_mwIndex, 'converged', field_ptr)

            ! iterations
            field_ptr = mxCreateDoubleScalar(real(output%ti_output%iterations, C_DOUBLE))
            call mxSetField(ti_struct, 1_mwIndex, 'iterations', field_ptr)

            ! residual_norm
            field_ptr = mxCreateDoubleScalar(output%ti_output%residual_norm)
            call mxSetField(ti_struct, 1_mwIndex, 'residual_norm', field_ptr)

            ! policies [n_het_endo × N_sp]
            if (allocated(output%ti_output%policies_sp)) then
                field_ptr = mxCreateDoubleMatrix(int(n_het_endo, mwSize), &
                                                 int(N_sp, mwSize), mxREAL)
                temp_real_2d(1:n_het_endo, 1:N_sp) => mxGetDoubles(field_ptr)
                temp_real_2d = output%ti_output%policies_sp
            else
                field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
            end if
            call mxSetField(ti_struct, 1_mwIndex, 'policies', field_ptr)
        else
            ! Fill with defaults if not associated
            field_ptr = mxCreateLogicalScalar(.false._mxLogical)
            call mxSetField(ti_struct, 1_mwIndex, 'converged', field_ptr)
            field_ptr = mxCreateDoubleScalar(0.0_C_DOUBLE)
            call mxSetField(ti_struct, 1_mwIndex, 'iterations', field_ptr)
            field_ptr = mxCreateDoubleScalar(0.0_C_DOUBLE)
            call mxSetField(ti_struct, 1_mwIndex, 'residual_norm', field_ptr)
            field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
            call mxSetField(ti_struct, 1_mwIndex, 'policies', field_ptr)
        end if

        call mxSetField(output_mx, 1_mwIndex, 'time_iteration', ti_struct)

        ! ================================================================
        ! PART 4: Distribution sub-struct
        ! ================================================================

        dist_struct = mxCreateStructMatrix(1_mwSize, 1_mwSize, dist_fields)
        if (.not. c_associated(dist_struct)) call mexErrMsgTxt("Failed to create distribution struct")

        if (associated(output%dist_output)) then
            ! converged
            field_ptr = mxCreateLogicalScalar(merge(.true._mxLogical, .false._mxLogical, output%dist_output%converged))
            call mxSetField(dist_struct, 1_mwIndex, 'converged', field_ptr)

            ! iterations
            field_ptr = mxCreateDoubleScalar(real(output%dist_output%iterations, C_DOUBLE))
            call mxSetField(dist_struct, 1_mwIndex, 'iterations', field_ptr)

            ! residual_norm
            field_ptr = mxCreateDoubleScalar(output%dist_output%residual_norm)
            call mxSetField(dist_struct, 1_mwIndex, 'residual_norm', field_ptr)

            ! hist [N_om × 1]
            if (associated(output%dist_output%distribution)) then
                field_ptr = mxCreateDoubleMatrix(int(N_om, mwSize), 1_mwSize, mxREAL)
                temp_real(1:N_om) => mxGetDoubles(field_ptr)
                temp_real = output%dist_output%distribution
            else
                field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
            end if
            call mxSetField(dist_struct, 1_mwIndex, 'hist', field_ptr)
        else
            ! Fill with defaults if not associated
            field_ptr = mxCreateLogicalScalar(.false._mxLogical)
            call mxSetField(dist_struct, 1_mwIndex, 'converged', field_ptr)
            field_ptr = mxCreateDoubleScalar(0.0_C_DOUBLE)
            call mxSetField(dist_struct, 1_mwIndex, 'iterations', field_ptr)
            field_ptr = mxCreateDoubleScalar(0.0_C_DOUBLE)
            call mxSetField(dist_struct, 1_mwIndex, 'residual_norm', field_ptr)
            field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
            call mxSetField(dist_struct, 1_mwIndex, 'hist', field_ptr)
        end if

        call mxSetField(output_mx, 1_mwIndex, 'distribution', dist_struct)

        ! ================================================================
        ! PART 5: Aggregates sub-struct
        ! ================================================================

        agg_struct = mxCreateStructMatrix(1_mwSize, 1_mwSize, agg_fields)
        if (.not. c_associated(agg_struct)) call mexErrMsgTxt("Failed to create aggregates struct")

        if (associated(output%agg_output)) then
            ! Ix [n_het_endo × 1]
            if (allocated(output%agg_output%Ix)) then
                field_ptr = mxCreateDoubleMatrix(int(n_het_endo, mwSize), 1_mwSize, mxREAL)
                temp_real(1:n_het_endo) => mxGetDoubles(field_ptr)
                temp_real = output%agg_output%Ix
            else
                field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
            end if
            call mxSetField(agg_struct, 1_mwIndex, 'Ix', field_ptr)

            ! residuals [n_agg_endo × 1]
            if (allocated(output%agg_output%residuals)) then
                field_ptr = mxCreateDoubleMatrix(int(n_agg_endo, mwSize), 1_mwSize, mxREAL)
                temp_real(1:n_agg_endo) => mxGetDoubles(field_ptr)
                temp_real = output%agg_output%residuals
            else
                field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
            end if
            call mxSetField(agg_struct, 1_mwIndex, 'residuals', field_ptr)
        else
            ! Fill with defaults if not associated
            field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
            call mxSetField(agg_struct, 1_mwIndex, 'Ix', field_ptr)
            field_ptr = mxCreateDoubleMatrix(0_mwSize, 0_mwSize, mxREAL)
            call mxSetField(agg_struct, 1_mwIndex, 'residuals', field_ptr)
        end if

        call mxSetField(output_mx, 1_mwIndex, 'aggregates', agg_struct)

    end subroutine pack_calibration_output

end module pack_output
