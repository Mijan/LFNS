! This file was automatically generated by SWIG (http://www.swig.org).
! Version 4.0.0
!
! Do not make changes to this file unless you know what you are doing--modify
! the SWIG interface file instead.

! ---------------------------------------------------------------
! Programmer(s): Auto-generated by swig.
! ---------------------------------------------------------------
! SUNDIALS Copyright Start
! Copyright (c) 2002-2019, Lawrence Livermore National Security
! and Southern Methodist University.
! All rights reserved.
!
! See the top-level LICENSE and NOTICE files for details.
!
! SPDX-License-Identifier: BSD-3-Clause
! SUNDIALS Copyright End
! ---------------------------------------------------------------

module fsunlinsol_spfgmr_mod
 use, intrinsic :: ISO_C_BINDING
 use fsundials_linearsolver_mod
 use fsundials_types_mod
 use fsundials_nvector_mod
 use fsundials_types_mod
 use fsundials_matrix_mod
 use fsundials_nvector_mod
 use fsundials_types_mod
 implicit none
 private

 ! DECLARATION CONSTRUCTS
 integer(C_INT), parameter, public :: SUNSPFGMR_MAXL_DEFAULT = 5_C_INT
 integer(C_INT), parameter, public :: SUNSPFGMR_MAXRS_DEFAULT = 0_C_INT
 public :: FSUNLinSol_SPFGMR
 public :: FSUNLinSol_SPFGMRSetPrecType
 public :: FSUNLinSol_SPFGMRSetGSType
 public :: FSUNLinSol_SPFGMRSetMaxRestarts
 public :: FSUNSPFGMR
 public :: FSUNSPFGMRSetPrecType
 public :: FSUNSPFGMRSetGSType
 public :: FSUNSPFGMRSetMaxRestarts
 public :: FSUNLinSolGetType_SPFGMR
 public :: FSUNLinSolGetID_SPFGMR
 public :: FSUNLinSolInitialize_SPFGMR
 public :: FSUNLinSolSetATimes_SPFGMR
 public :: FSUNLinSolSetPreconditioner_SPFGMR
 public :: FSUNLinSolSetScalingVectors_SPFGMR
 public :: FSUNLinSolSetup_SPFGMR
 public :: FSUNLinSolSolve_SPFGMR
 public :: FSUNLinSolNumIters_SPFGMR
 public :: FSUNLinSolResNorm_SPFGMR
 public :: FSUNLinSolResid_SPFGMR
 public :: FSUNLinSolLastFlag_SPFGMR
 public :: FSUNLinSolSpace_SPFGMR
 public :: FSUNLinSolFree_SPFGMR

! WRAPPER DECLARATIONS
interface
function swigc_FSUNLinSol_SPFGMR(farg1, farg2, farg3) &
bind(C, name="_wrap_FSUNLinSol_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT), intent(in) :: farg3
type(C_PTR) :: fresult
end function

function swigc_FSUNLinSol_SPFGMRSetPrecType(farg1, farg2) &
bind(C, name="_wrap_FSUNLinSol_SPFGMRSetPrecType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSol_SPFGMRSetGSType(farg1, farg2) &
bind(C, name="_wrap_FSUNLinSol_SPFGMRSetGSType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSol_SPFGMRSetMaxRestarts(farg1, farg2) &
bind(C, name="_wrap_FSUNLinSol_SPFGMRSetMaxRestarts") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNSPFGMR(farg1, farg2, farg3) &
bind(C, name="_wrap_FSUNSPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT), intent(in) :: farg3
type(C_PTR) :: fresult
end function

function swigc_FSUNSPFGMRSetPrecType(farg1, farg2) &
bind(C, name="_wrap_FSUNSPFGMRSetPrecType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNSPFGMRSetGSType(farg1, farg2) &
bind(C, name="_wrap_FSUNSPFGMRSetGSType") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNSPFGMRSetMaxRestarts(farg1, farg2) &
bind(C, name="_wrap_FSUNSPFGMRSetMaxRestarts") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT), intent(in) :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolGetType_SPFGMR(farg1) &
bind(C, name="_wrap_FSUNLinSolGetType_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolGetID_SPFGMR(farg1) &
bind(C, name="_wrap_FSUNLinSolGetID_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolInitialize_SPFGMR(farg1) &
bind(C, name="_wrap_FSUNLinSolInitialize_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolSetATimes_SPFGMR(farg1, farg2, farg3) &
bind(C, name="_wrap_FSUNLinSolSetATimes_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
type(C_FUNPTR), value :: farg3
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolSetPreconditioner_SPFGMR(farg1, farg2, farg3, farg4) &
bind(C, name="_wrap_FSUNLinSolSetPreconditioner_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
type(C_FUNPTR), value :: farg3
type(C_FUNPTR), value :: farg4
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolSetScalingVectors_SPFGMR(farg1, farg2, farg3) &
bind(C, name="_wrap_FSUNLinSolSetScalingVectors_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
type(C_PTR), value :: farg3
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolSetup_SPFGMR(farg1, farg2) &
bind(C, name="_wrap_FSUNLinSolSetup_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolSolve_SPFGMR(farg1, farg2, farg3, farg4, farg5) &
bind(C, name="_wrap_FSUNLinSolSolve_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
type(C_PTR), value :: farg3
type(C_PTR), value :: farg4
real(C_DOUBLE), intent(in) :: farg5
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolNumIters_SPFGMR(farg1) &
bind(C, name="_wrap_FSUNLinSolNumIters_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolResNorm_SPFGMR(farg1) &
bind(C, name="_wrap_FSUNLinSolResNorm_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
real(C_DOUBLE) :: fresult
end function

function swigc_FSUNLinSolResid_SPFGMR(farg1) &
bind(C, name="_wrap_FSUNLinSolResid_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR) :: fresult
end function

function swigc_FSUNLinSolLastFlag_SPFGMR(farg1) &
bind(C, name="_wrap_FSUNLinSolLastFlag_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_LONG) :: fresult
end function

function swigc_FSUNLinSolSpace_SPFGMR(farg1, farg2, farg3) &
bind(C, name="_wrap_FSUNLinSolSpace_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
type(C_PTR), value :: farg2
type(C_PTR), value :: farg3
integer(C_INT) :: fresult
end function

function swigc_FSUNLinSolFree_SPFGMR(farg1) &
bind(C, name="_wrap_FSUNLinSolFree_SPFGMR") &
result(fresult)
use, intrinsic :: ISO_C_BINDING
type(C_PTR), value :: farg1
integer(C_INT) :: fresult
end function

end interface


contains
 ! MODULE SUBPROGRAMS
function FSUNLinSol_SPFGMR(y, pretype, maxl) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
type(SUNLinearSolver), pointer :: swig_result
type(N_Vector), target, intent(inout) :: y
integer(C_INT), intent(in) :: pretype
integer(C_INT), intent(in) :: maxl
type(C_PTR) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 
integer(C_INT) :: farg3 

farg1 = c_loc(y)
farg2 = pretype
farg3 = maxl
fresult = swigc_FSUNLinSol_SPFGMR(farg1, farg2, farg3)
call c_f_pointer(fresult, swig_result)
end function

function FSUNLinSol_SPFGMRSetPrecType(s, pretype) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT), intent(in) :: pretype
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = c_loc(s)
farg2 = pretype
fresult = swigc_FSUNLinSol_SPFGMRSetPrecType(farg1, farg2)
swig_result = fresult
end function

function FSUNLinSol_SPFGMRSetGSType(s, gstype) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT), intent(in) :: gstype
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = c_loc(s)
farg2 = gstype
fresult = swigc_FSUNLinSol_SPFGMRSetGSType(farg1, farg2)
swig_result = fresult
end function

function FSUNLinSol_SPFGMRSetMaxRestarts(s, maxrs) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT), intent(in) :: maxrs
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = c_loc(s)
farg2 = maxrs
fresult = swigc_FSUNLinSol_SPFGMRSetMaxRestarts(farg1, farg2)
swig_result = fresult
end function

function FSUNSPFGMR(y, pretype, maxl) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
type(SUNLinearSolver), pointer :: swig_result
type(N_Vector), target, intent(inout) :: y
integer(C_INT), intent(in) :: pretype
integer(C_INT), intent(in) :: maxl
type(C_PTR) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 
integer(C_INT) :: farg3 

farg1 = c_loc(y)
farg2 = pretype
farg3 = maxl
fresult = swigc_FSUNSPFGMR(farg1, farg2, farg3)
call c_f_pointer(fresult, swig_result)
end function

function FSUNSPFGMRSetPrecType(s, pretype) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT), intent(in) :: pretype
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = c_loc(s)
farg2 = pretype
fresult = swigc_FSUNSPFGMRSetPrecType(farg1, farg2)
swig_result = fresult
end function

function FSUNSPFGMRSetGSType(s, gstype) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT), intent(in) :: gstype
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = c_loc(s)
farg2 = gstype
fresult = swigc_FSUNSPFGMRSetGSType(farg1, farg2)
swig_result = fresult
end function

function FSUNSPFGMRSetMaxRestarts(s, maxrs) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT), intent(in) :: maxrs
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
integer(C_INT) :: farg2 

farg1 = c_loc(s)
farg2 = maxrs
fresult = swigc_FSUNSPFGMRSetMaxRestarts(farg1, farg2)
swig_result = fresult
end function

function FSUNLinSolGetType_SPFGMR(s) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(SUNLinearSolver_Type) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(s)
fresult = swigc_FSUNLinSolGetType_SPFGMR(farg1)
swig_result = fresult
end function

function FSUNLinSolGetID_SPFGMR(s) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(SUNLinearSolver_ID) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(s)
fresult = swigc_FSUNLinSolGetID_SPFGMR(farg1)
swig_result = fresult
end function

function FSUNLinSolInitialize_SPFGMR(s) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(s)
fresult = swigc_FSUNLinSolInitialize_SPFGMR(farg1)
swig_result = fresult
end function

function FSUNLinSolSetATimes_SPFGMR(s, a_data, atimes) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
type(C_PTR) :: a_data
type(C_FUNPTR), intent(in), value :: atimes
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 
type(C_FUNPTR) :: farg3 

farg1 = c_loc(s)
farg2 = a_data
farg3 = atimes
fresult = swigc_FSUNLinSolSetATimes_SPFGMR(farg1, farg2, farg3)
swig_result = fresult
end function

function FSUNLinSolSetPreconditioner_SPFGMR(s, p_data, pset, psol) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
type(C_PTR) :: p_data
type(C_FUNPTR), intent(in), value :: pset
type(C_FUNPTR), intent(in), value :: psol
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 
type(C_FUNPTR) :: farg3 
type(C_FUNPTR) :: farg4 

farg1 = c_loc(s)
farg2 = p_data
farg3 = pset
farg4 = psol
fresult = swigc_FSUNLinSolSetPreconditioner_SPFGMR(farg1, farg2, farg3, farg4)
swig_result = fresult
end function

function FSUNLinSolSetScalingVectors_SPFGMR(s, s1, s2) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
type(N_Vector), target, intent(inout) :: s1
type(N_Vector), target, intent(inout) :: s2
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 
type(C_PTR) :: farg3 

farg1 = c_loc(s)
farg2 = c_loc(s1)
farg3 = c_loc(s2)
fresult = swigc_FSUNLinSolSetScalingVectors_SPFGMR(farg1, farg2, farg3)
swig_result = fresult
end function

function FSUNLinSolSetup_SPFGMR(s, a) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
type(SUNMatrix), target, intent(inout) :: a
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 

farg1 = c_loc(s)
farg2 = c_loc(a)
fresult = swigc_FSUNLinSolSetup_SPFGMR(farg1, farg2)
swig_result = fresult
end function

function FSUNLinSolSolve_SPFGMR(s, a, x, b, tol) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
type(SUNMatrix), target, intent(inout) :: a
type(N_Vector), target, intent(inout) :: x
type(N_Vector), target, intent(inout) :: b
real(C_DOUBLE), intent(in) :: tol
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 
type(C_PTR) :: farg3 
type(C_PTR) :: farg4 
real(C_DOUBLE) :: farg5 

farg1 = c_loc(s)
farg2 = c_loc(a)
farg3 = c_loc(x)
farg4 = c_loc(b)
farg5 = tol
fresult = swigc_FSUNLinSolSolve_SPFGMR(farg1, farg2, farg3, farg4, farg5)
swig_result = fresult
end function

function FSUNLinSolNumIters_SPFGMR(s) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(s)
fresult = swigc_FSUNLinSolNumIters_SPFGMR(farg1)
swig_result = fresult
end function

function FSUNLinSolResNorm_SPFGMR(s) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
real(C_DOUBLE) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
real(C_DOUBLE) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(s)
fresult = swigc_FSUNLinSolResNorm_SPFGMR(farg1)
swig_result = fresult
end function

function FSUNLinSolResid_SPFGMR(s) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
type(N_Vector), pointer :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
type(C_PTR) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(s)
fresult = swigc_FSUNLinSolResid_SPFGMR(farg1)
call c_f_pointer(fresult, swig_result)
end function

function FSUNLinSolLastFlag_SPFGMR(s) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_LONG) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_LONG) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(s)
fresult = swigc_FSUNLinSolLastFlag_SPFGMR(farg1)
swig_result = fresult
end function

function FSUNLinSolSpace_SPFGMR(s, lenrwls, leniwls) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_LONG), dimension(*), target, intent(inout) :: lenrwls
integer(C_LONG), dimension(*), target, intent(inout) :: leniwls
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 
type(C_PTR) :: farg2 
type(C_PTR) :: farg3 

farg1 = c_loc(s)
farg2 = c_loc(lenrwls(1))
farg3 = c_loc(leniwls(1))
fresult = swigc_FSUNLinSolSpace_SPFGMR(farg1, farg2, farg3)
swig_result = fresult
end function

function FSUNLinSolFree_SPFGMR(s) &
result(swig_result)
use, intrinsic :: ISO_C_BINDING
integer(C_INT) :: swig_result
type(SUNLinearSolver), target, intent(inout) :: s
integer(C_INT) :: fresult 
type(C_PTR) :: farg1 

farg1 = c_loc(s)
fresult = swigc_FSUNLinSolFree_SPFGMR(farg1)
swig_result = fresult
end function


end module
