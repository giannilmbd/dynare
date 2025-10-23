! Fortran 2018 interface for a subset of MEX functions
!
! For some functions, exposing the C interface directly is not convenient (e.g.
! when they use C strings, since those need to be null-terminated, or when they
! return pointers to arrays of doubles/chars/integers…). A more Fortran-ish
! interface is provided for those, using a small glue code.
!
! Things to be aware of when adding new interfaces this file:
!
! – We only support the so-called interleaved API, available since MATLAB R2018a
!   and Octave 7 (but for the latter broken before 8.4.0).
! — The tricky part is to deal with API versioning.
!   For MATLAB ⩾ R2020a, all versions share the same API (embodied in API_VER
!   define), but in the past the API would differ across function, so there
!   was an API_VER2 define for those.
!   For each function, the information can be retrieved from either matrix.h or
!   mex.h.
!   Under Octave, some functions have a different symbol name with the
!   interleaved API, so we handle them through the API_VER_INTERLEAVED define;
!   see octave/mexproto.h for the list of affected functions.
! — C passes arguments by value, so the “value” keyword is often needed
! — Strings passed to C must be null terminated (hence a wrapper is needed to
!   append c_null_char)
! — We follow the Fortran convention that indices start at one. Hence, when
!   interfacing a function with a mwIndex argument, it is necessary to write
!   a glue code that substracts one to the index
! — When writing glue code, using the pure C interface as a starting point:
!   • remove the “use” declarations
!   • remove the “value” keywords
!   • convert input character arrays to character(kind=c_char, len=*)
!   • Fortran array pointers returned by the glue code must be marked
!     “contiguous” (which is always the case for C arrays). This will help the
!     Fortran compiler better optimize the code (in some cases, this will avoid
!     array copies)
! – If the function has no side-effect, mark it as “pure”, to avoid gfortran
!   warnings when the function may not be evaluated depending on the branch
!   (-Wfunction-elimination)

! Copyright © 2019-2025 Dynare Team
!
! This file is part of Dynare.
!
! Dynare is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Dynare is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

#ifdef MATLAB_MEX_FILE
# define API_VER "_800"
# define API_VER_INTERLEAVED API_VER
#else
! Octave
# define API_VER ""
# define API_VER_INTERLEAVED "_interleaved"
#endif

!!! C Matrix API
!!! Listed in same order as https://fr.mathworks.com/help/matlab/cc-mx-matrix-library.html
module matlab_mat
  use iso_fortran_env
  use iso_c_binding
  implicit none (type, external)

  !! C Data Types
  integer, parameter :: mwSize = c_size_t
  integer, parameter :: mwIndex = c_size_t
  integer, parameter :: mwSignedIndex = c_intptr_t
  integer, parameter :: mxLogical = c_bool
  integer, parameter :: mxComplexity = c_int
  integer, parameter :: mxClassID = c_int

  integer(mxComplexity), parameter :: mxREAL = 0
  integer(mxComplexity), parameter :: mxCOMPLEX = 1

  integer(mxClassID), parameter :: mxUNKNOWN_CLASS = 0
  integer(mxClassID), parameter :: mxCELL_CLASS = 1
  integer(mxClassID), parameter :: mxSTRUCT_CLASS = 2
  integer(mxClassID), parameter :: mxLOGICAL_CLASS = 3
  integer(mxClassID), parameter :: mxCHAR_CLASS = 4
  integer(mxClassID), parameter :: mxVOID_CLASS = 5
  integer(mxClassID), parameter :: mxDOUBLE_CLASS = 6
  integer(mxClassID), parameter :: mxSINGLE_CLASS = 7
  integer(mxClassID), parameter :: mxINT8_CLASS = 8
  integer(mxClassID), parameter :: mxUINT8_CLASS = 9
  integer(mxClassID), parameter :: mxINT16_CLASS = 10
  integer(mxClassID), parameter :: mxUINT16_CLASS = 11
  integer(mxClassID), parameter :: mxINT32_CLASS = 12
  integer(mxClassID), parameter :: mxUINT32_CLASS = 13
  integer(mxClassID), parameter :: mxINT64_CLASS = 14
  integer(mxClassID), parameter :: mxUINT64_CLASS = 15
  integer(mxClassID), parameter :: mxFUNCTION_CLASS = 16
  integer(mxClassID), parameter :: mxOPAQUE_CLASS = 17
  integer(mxClassID), parameter :: mxOBJECT_CLASS = 18
#if defined(_LP64) || defined(_WIN64)
  integer(mxClassID), parameter :: mxINDEX_CLASS = mxUINT64_CLASS
#else
  integer(mxClassID), parameter :: mxINDEX_CLASS = mxUINT32_CLASS
#endif

  interface
     !! mxArray attributes
     pure logical(c_bool) function mxIsNumeric(pm) bind(c, name="mxIsNumeric"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsNumeric

     pure logical(c_bool) function mxIsComplex(pm) bind(c, name="mxIsComplex"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsComplex

     pure integer(c_size_t) function mxGetNumberOfElements(pm) bind(c, name="mxGetNumberOfElements"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetNumberOfElements

     pure integer(c_size_t) function mxGetM(pm) bind(c, name="mxGetM"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetM

     pure integer(c_size_t) function mxGetN(pm) bind(c, name="mxGetN"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetN

     pure logical(c_bool) function mxIsEmpty(pm) bind(c, name="mxIsEmpty"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsEmpty

     !! Create, Query, and Access Data Types

     ! Numeric types
     type(c_ptr) function mxCreateDoubleMatrix(m, n, ComplexFlag) bind(c, name="mxCreateDoubleMatrix"//API_VER_INTERLEAVED)
       use iso_c_binding
       import :: mwSize, mxComplexity
       integer(mwSize), intent(in), value :: m, n
       integer(mxComplexity), intent(in), value :: ComplexFlag
     end function mxCreateDoubleMatrix

     type(c_ptr) function mxCreateDoubleScalar(value) bind(c, name="mxCreateDoubleScalar"//API_VER_INTERLEAVED)
       use iso_c_binding
       real(c_double), intent(in), value :: value
     end function mxCreateDoubleScalar

     type(c_ptr) function mxCreateNumericMatrix(m, n, ClassID, ComplexFlag) bind(c, name="mxCreateNumericMatrix"//API_VER_INTERLEAVED)
       use iso_c_binding
       import :: mwSize, mxClassID, mxComplexity
       integer(mwSize), intent(in), value :: m, n
       integer(mxClassID), intent(in), value :: ClassID
       integer(mxComplexity), intent(in), value :: ComplexFlag
     end

     ! Noncomplex Float
     pure logical(c_bool) function mxIsScalar(array_ptr) bind(c, name="mxIsScalar"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: array_ptr
     end function mxIsScalar

     pure real(c_double) function mxGetScalar(pm) bind(c, name="mxGetScalar"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetScalar

     pure logical(c_bool) function mxIsDouble(pm) bind(c, name="mxIsDouble"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsDouble

     pure type(c_ptr) function mxGetDoubles_internal(pm) bind(c, name="mxGetDoubles"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetDoubles_internal

     pure type(c_ptr) function mxGetPr_internal(pm) bind(c, name="mxGetPr"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetPr_internal

     ! Noncomplex integer

     pure logical(c_bool) function mxIsInt32(pm) bind(c, name="mxIsInt32"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsInt32

     pure type(c_ptr) function mxGetInt32s_internal(pa) bind(c, name="mxGetInt32s"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pa
     end function mxGetInt32s_internal

     pure logical(c_bool) function mxIsInt64(pm) bind(c, name="mxIsInt64"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsInt64

     pure type(c_ptr) function mxGetInt64s_internal(pa) bind(c, name="mxGetInt64s"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pa
     end function mxGetInt64s_internal

     ! Complex Float
     pure type(c_ptr) function mxGetComplexDoubles_internal(pa) bind(c, name="mxGetComplexDoubles"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pa
     end function mxGetComplexDoubles_internal

     ! Sparse
     type(c_ptr) function mxCreateSparse(m, n, nzmax, ComplexFlag) bind(c, name="mxCreateSparse"//API_VER_INTERLEAVED)
       use iso_c_binding
       import :: mwSize, mxComplexity
       integer(mwSize), intent(in), value :: m, n, nzmax
       integer(mxComplexity), intent(in), value :: ComplexFlag
     end function mxCreateSparse

     pure logical(c_bool) function mxIsSparse(pm) bind(c, name="mxIsSparse"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsSparse

     pure integer(c_size_t) function mxGetNzmax(pm) bind(c, name="mxGetNzmax"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetNzmax

     pure type(c_ptr) function mxGetIr_internal(pm) bind(c, name="mxGetIr"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetIr_internal

     pure type(c_ptr) function mxGetJc_internal(pm) bind(c, name="mxGetJc"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetJc_internal

     ! Nonnumeric types
     pure type(c_ptr) function mxGetData(pm) bind(c, name="mxGetData"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetData

     ! Character
     pure logical(c_bool) function mxIsChar(pm) bind(c, name="mxIsChar"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsChar

     ! Logical
     pure logical(c_bool) function mxIsLogical(pm) bind(c, name="mxIsLogical"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsLogical

     pure logical(c_bool) function mxIsLogicalScalar(array_ptr) bind(c, name="mxIsLogicalScalar"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: array_ptr
     end function mxIsLogicalScalar

     type(c_ptr) function mxCreateLogicalScalar(value) bind(c, name="mxCreateLogicalScalar"//API_VER_INTERLEAVED)
       use iso_c_binding
       import :: mxLogical
       logical(mxLogical), intent(in), value :: value
     end function mxCreateLogicalScalar

     pure type(c_ptr) function mxGetLogicals_internal(array_ptr) bind(c, name="mxGetLogicals"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: array_ptr
     end function mxGetLogicals_internal

     ! Object
     pure logical(c_bool) function mxIsClass_internal(pm, classname) bind(c, name="mxIsClass"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
       character(c_char), dimension(*), intent(in) :: classname
     end function mxIsClass_internal

     ! Structure
     type(c_ptr) function mxCreateStructMatrix_internal(m, n, nfields, fieldnames) bind(c, name="mxCreateStructMatrix" &
          //API_VER_INTERLEAVED)
       use iso_c_binding
       import :: mwSize
       integer(mwSize), intent(in), value :: m, n
       integer(c_int), intent(in), value :: nfields
       type(c_ptr), dimension(*), intent(in) :: fieldnames
     end function mxCreateStructMatrix_internal

     pure logical(c_bool) function mxIsStruct(pm) bind(c, name="mxIsStruct"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsStruct

     pure type(c_ptr) function mxGetField_internal(pm, index, fieldname) bind(c, name="mxGetField"//API_VER)
       use iso_c_binding
       import :: mwIndex
       type(c_ptr), intent(in), value :: pm
       integer(mwIndex), intent(in), value :: index
       character(c_char), dimension(*), intent(in) :: fieldname
     end function mxGetField_internal

     subroutine mxSetField_internal(pm, index, fieldname, pvalue) bind(c, name="mxSetField"//API_VER)
       use iso_c_binding
       import :: mwIndex
       type(c_ptr), intent(in), value :: pm
       type(c_ptr), intent(in), value :: pvalue
       integer(mwIndex), intent(in), value :: index
       character(c_char), dimension(*), intent(in) :: fieldname
     end subroutine mxSetField_internal

     pure integer(c_int) function mxGetNumberOfFields(pm) bind(c, name="mxGetNumberOfFields"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetNumberOfFields

     ! Cell
     pure logical(c_bool) function mxIsCell(pm) bind(c, name="mxIsCell"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsCell

     pure type(c_ptr) function mxGetCell_internal(pm, index) bind(c, name="mxGetCell"//API_VER)
       use iso_c_binding
       import :: mwIndex
       type(c_ptr), intent(in), value :: pm
       integer(mwIndex), intent(in), value :: index
     end function mxGetCell_internal

     !! Delete and Duplicate mxArray
     subroutine mxDestroyArray(pm) bind(c, name="mxDestroyArray"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end subroutine mxDestroyArray

     type(c_ptr) function mxDuplicateArray(in) bind(c, name="mxDuplicateArray"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: in
     end function mxDuplicateArray

     !! Convert mxArray

     ! Character
     type(c_ptr) function mxArrayToString_internal(array_ptr) bind(c, name="mxArrayToString"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: array_ptr
     end function mxArrayToString_internal


     !! Memory management

     subroutine mxFree(ptr) bind(c, name="mxFree"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in) :: ptr
     end subroutine mxFree
  end interface
contains
  ! Some helper functions to make the interface more Fortran-ish

  function mxGetDoubles(pm)
    type(c_ptr), intent(in) :: pm
    real(real64), dimension(:), pointer, contiguous :: mxGetDoubles
    call c_f_pointer(mxGetDoubles_internal(pm), mxGetDoubles, [ mxGetNumberOfElements(pm) ])
  end function mxGetDoubles

  function mxGetPr(pm)
    type(c_ptr), intent(in) :: pm
    real(real64), dimension(:), pointer, contiguous :: mxGetPr
    call c_f_pointer(mxGetPr_internal(pm), mxGetPr, [ mxGetNumberOfElements(pm) ])
  end function mxGetPr

  function mxGetInt32s(pa)
    type(c_ptr), intent(in) :: pa
    integer(int32), dimension(:), pointer, contiguous :: mxGetInt32s
    call c_f_pointer(mxGetInt32s_internal(pa), mxGetInt32s, [ mxGetNumberOfElements(pa) ])
  end function mxGetInt32s

  function mxGetInt64s(pa)
    type(c_ptr), intent(in) :: pa
    integer(int64), dimension(:), pointer, contiguous :: mxGetInt64s
    call c_f_pointer(mxGetInt64s_internal(pa), mxGetInt64s, [ mxGetNumberOfElements(pa) ])
  end function mxGetInt64s

  function mxGetComplexDoubles(pa)
    type(c_ptr), intent(in) :: pa
    complex(real64), dimension(:), pointer, contiguous :: mxGetComplexDoubles
    call c_f_pointer(mxGetComplexDoubles_internal(pa), mxGetComplexDoubles, [ mxGetNumberOfElements(pa) ])
  end function mxGetComplexDoubles

  function mxGetIr(pm)
    type(c_ptr), intent(in) :: pm
    integer(mwIndex), dimension(:), pointer, contiguous :: mxGetIr
    call c_f_pointer(mxGetIr_internal(pm), mxGetIr, [ mxGetNzmax(pm) ])
  end function mxGetIr

  function mxGetJc(pm)
    type(c_ptr), intent(in) :: pm
    integer(mwIndex), dimension(:), pointer, contiguous :: mxGetJc
    call c_f_pointer(mxGetJc_internal(pm), mxGetJc, [ mxGetN(pm)+1 ])
  end function mxGetJc

  function mxGetLogicals(array_ptr)
    type(c_ptr), intent(in) :: array_ptr
    logical(mxLogical), dimension(:), pointer, contiguous :: mxGetLogicals
    call c_f_pointer(mxGetLogicals_internal(array_ptr), mxGetLogicals, [ mxGetNumberOfElements(array_ptr) ])
  end function mxGetLogicals

  pure logical(c_bool) function mxIsClass(pm, classname)
    type(c_ptr), intent(in) :: pm
    character(kind=c_char, len=*), intent(in) :: classname
    mxIsclass = mxIsclass_internal(pm, classname // c_null_char)
  end function mxIsClass

  type(c_ptr) function mxCreateStructMatrix(m, n, fieldnames)
    integer(mwSize), intent(in) :: m, n
    character(kind=c_char, len=*), dimension(:), intent(in) :: fieldnames

    character(kind=c_char, len=len(fieldnames)+1), dimension(size(fieldnames)), target :: fieldnames_zero ! Stores zero-terminated strings
    type(c_ptr), dimension(size(fieldnames)) :: fieldnames_ptr ! C arrays of strings
    integer :: i
    do i = 1,size(fieldnames)
       fieldnames_zero(i) = trim(fieldnames(i)) // c_null_char
       fieldnames_ptr(i) = c_loc(fieldnames_zero(i))
    end do
    mxCreateStructMatrix = mxCreateStructMatrix_internal(m, n, int(size(fieldnames), c_int), fieldnames_ptr)
  end function mxCreateStructMatrix

  pure type(c_ptr) function mxGetField(pm, index, fieldname)
    type(c_ptr), intent(in) :: pm
    integer(mwIndex), intent(in) :: index
    character(kind=c_char, len=*), intent(in) :: fieldname
    mxGetField = mxGetField_internal(pm, index-1, fieldname // c_null_char)
  end function mxGetField

  subroutine mxSetField(pm, index, fieldname, pvalue)
       type(c_ptr), intent(in) :: pm
       type(c_ptr), intent(in) :: pvalue
       integer(mwIndex), intent(in) :: index
       character(kind=c_char, len=*), intent(in) :: fieldname
       call mxSetField_internal(pm, index-1, fieldname // c_null_char, pvalue)
  end subroutine mxSetField

  pure type(c_ptr) function mxGetCell(pm, index)
    type(c_ptr), intent(in) :: pm
    integer(mwIndex), intent(in) :: index
    mxGetCell = mxGetCell_internal(pm, index-1)
  end function mxGetCell

  function mxArrayToString(pm)
    type(c_ptr), intent(in) :: pm
    character(kind=c_char, len=:), allocatable :: mxArrayToString
    character(kind=c_char), dimension(:), pointer, contiguous :: chararray
    type(c_ptr) :: ptr
    integer :: i
    ptr = mxArrayToString_internal(pm)
    call c_f_pointer(ptr, chararray, [ mxGetNumberOfElements(pm) ])
    ! Convert the character array into a character scalar (of length > 1)
    allocate(character(kind=c_char, len=size(chararray)) :: mxArrayToString)
    do i=1,size(chararray)
       mxArrayToString(i:i) = chararray(i)
    end do
    call mxFree(ptr)
  end function mxArrayToString
end module matlab_mat


!!! C MEX API
!!! Listed in same order as https://fr.mathworks.com/help/matlab/call-mex-files-1.html
module matlab_mex
  use matlab_mat
  implicit none (type, external)

  interface
     integer(c_int) function mexCallMATLAB_internal(nlhs, plhs, nrhs, prhs, functionName) bind(c, name="mexCallMATLAB"//API_VER)
       use iso_c_binding
       integer(c_int), intent(in), value :: nlhs, nrhs
       type(c_ptr), dimension(*), intent(in) :: plhs, prhs
       character(c_char), dimension(*), intent(in) :: functionName
     end function mexCallMATLAB_internal

     type(c_ptr) function mexCallMATLABWithTrap_internal(nlhs, plhs, nrhs, prhs, functionName) &
          bind(c, name="mexCallMATLABWithTrap"//API_VER)
       use iso_c_binding
       integer(c_int), intent(in), value :: nlhs, nrhs
       type(c_ptr), dimension(*), intent(in) :: plhs, prhs
       character(c_char), dimension(*), intent(in) :: functionName
     end function mexCallMATLABWithTrap_internal

     subroutine mexErrMsgTxt_internal(msg) bind(c, name="mexErrMsgTxt"//API_VER)
       use iso_c_binding
       character(c_char), dimension(*), intent(in) :: msg
     end subroutine mexErrMsgTxt_internal

     subroutine mexErrMsgIdAndTxt_internal(id, msg) bind(c, name="mexErrMsgIdAndTxt"//API_VER)
       use iso_c_binding
       character(c_char), dimension(*), intent(in) :: id, msg
     end subroutine mexErrMsgIdAndTxt_internal

     subroutine mexPrintf_internal(message) bind(c, name="mexPrintf"//API_VER)
       use iso_c_binding
       character(c_char), dimension(*), intent(in) :: message
     end subroutine mexPrintf_internal
  end interface
contains
  ! Some helper functions to make the interface more Fortran-ish

  integer(c_int) function mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName)
    integer(c_int), intent(in) :: nlhs, nrhs
    type(c_ptr), dimension(*), intent(in) :: plhs, prhs
    character(kind=c_char, len=*), intent(in) :: functionName
    mexCallMATLAB = mexCallMATLAB_internal(nlhs, plhs, nrhs, prhs, functionName // c_null_char)
  end function mexCallMATLAB

  type(c_ptr) function mexCallMATLABWithTrap(nlhs, plhs, nrhs, prhs, functionName)
    integer(c_int), intent(in) :: nlhs, nrhs
    type(c_ptr), dimension(*), intent(in) :: plhs, prhs
    character(kind=c_char, len=*), intent(in) :: functionName
    mexCallMATLABWithTrap = mexCallMATLABWithTrap_internal(nlhs, plhs, nrhs, prhs, functionName // c_null_char)
  end function mexCallMATLABWithTrap

  subroutine mexErrMsgTxt(msg)
    character(kind=c_char, len=*), intent(in) :: msg
    call mexErrMsgTxt_internal(msg // c_null_char)
  end subroutine mexErrMsgTxt

  subroutine mexErrMsgIdAndTxt(id, msg)
    character(kind=c_char, len=*), intent(in) :: id, msg
    call mexErrMsgIdAndTxt_internal(id // c_null_char, msg // c_null_char)
  end subroutine mexErrMsgIdAndTxt

  subroutine mexPrintf(message)
    character(kind=c_char, len=*), intent(in) :: message
    call mexPrintf_internal(message // c_null_char)
  end subroutine mexPrintf

  ! Same as mexPrintf(), but trims trailing whitespace, and adds a new line
  subroutine mexPrintf_trim_newline(message)
    character(kind=c_char, len=*), intent(in) :: message
    call mexPrintf_internal(trim(message) // new_line(message) // c_null_char)
  end subroutine mexPrintf_trim_newline
end module matlab_mex
