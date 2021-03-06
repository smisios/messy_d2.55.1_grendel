MODULE YOEPHY

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     -----------------------------------------------------------------
!*    ** *YOEPHY* - SWITCHES RELATED TO DIABATIC PROCESSES
!     -----------------------------------------------------------------

!        * E.C.M.W.F. PHYSICS PACKAGE *

LOGICAL LEPHYS
LOGICAL LECOND
LOGICAL LECUMF
LOGICAL LEDCLD
LOGICAL LEEVAP
LOGICAL LEGWDG
LOGICAL LEOZOC
LOGICAL LEQNGT
LOGICAL LERADI
LOGICAL LERADS
LOGICAL LESHCV
LOGICAL LESICE
LOGICAL LESURF
LOGICAL LEVDIF
LOGICAL LAGPHY
LOGICAL LEPCLD
LOGICAL LEO3CH
LOGICAL LBUD23
LOGICAL LEMETHOX
LOGICAL LERA40

!     J.-J. MORCRETTE       E.C.M.W.F.      91/07/14

!  NAME     TYPE     PURPOSE
!  ----  :  ----   : ---------------------------------------------------
! LEPHYS : LOGICAL : SWITCH THE FULL E.C.M.W.F. PHYSICS PACKAGE ON
! LAGPHY : LOGICAL : IF TRUE, PHYSICS PACKAGE CALLED IN LAGGED MODE
! LECOND : LOGICAL : TURN THE LARGE-SCALE CONDENSATION ON
! LECUMF : LOGICAL : TURN THE MASS-FLUX CUMULUS CONVECTION SCHEME ON
! LEDCLD : LOGICAL : TURN THE DIAGNOSTIC CLOUD SCHEME ON
! LEPCLD : LOGICAL : TURN THE PROGNOSTIC CLOUD SCHEME ON
! LEEVAP : LOGICAL : TURN THE EVAPORATION OF PRECIPITATION ON
! LEGWDG : LOGICAL : TURN THE GRAVITY WAVE DRAG ON
! LEOZOC : LOGICAL : TURN THE CLIMATOLOGICAL OZONE ON
! LEQNGT : LOGICAL : TURN THE NEGATIVE HUMIDITY FIXER ON
! LERADI : LOGICAL : TURN THE RADIATION SCHEME ON
! LERADS : LOGICAL : TURN THE INTERACTIVE SURFACE RADIATIVE PROPERTIESON
! LESHCV : LOGICAL : TURN THE SHALLOW CONV. IN THE MASS-FLUX SCHEME ON
! LESICE : LOGICAL : TURN THE INTERACTIVE SEA ICE PROCESSES ON
! LESURF : LOGICAL : TURN THE INTERACTIVE SURFACE PROCESSES ON
! LEVDIF : LOGICAL : TURN THE VERTICAL DIFFUSION ON
! LEO3CH : LOGICAL : TURN THE O3 CHEMISTRY ON (for EC prog. ozone)
! LBUD23 : LOGICAL : SWITCH FOR 3 AND 2 DIMENSIONAL BUDGETS 
! LEMETHOX: LOGICAL : TURN THE METHANE OXIDATION ON
! LERA40 : LOGICAL  : EXTRA PHYSICS DIAGNOSTICS FOR ERA40
! LTE... : LOGICAL : CHECK OF I/O ARGUMENTS FOR EVERY PHYSICS ROUTINE
!     -----------------------------------------------------------------
END MODULE YOEPHY
