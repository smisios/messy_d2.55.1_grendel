! ====================================================================
!
! TITLE    : H2OEMIS
! 
! PROJECT  : STRATOFLY
!
! Module   : smcl/messy_h2oemis
!
! DATE     : 01.10.2019
!
! NOTE     : Comment format '!>' is used for generation of code 
!              documentation from source file via 'doxygen'.
! 
!> \brief
!> This is  documentation for the  Submodel Core Layer (SMCL) routine. 
!> Here, the submodel is integrated in the SMCL. \n
!
!> \authors Johannes Emmerig, DLR Oberpfaffenhofen, 2019,
!>          (Johannes.Emmerig@dlr.de)
!>          - Development of 'H2OEMIS'
!
!> \authors Patrick Jöckel, DLR Oberpfaffenhofen, 2019,
!>          (Patrick.Joeckel@dlr.de)
!>          - Mentoring the development
!
!> \version 1.0
!>     - Integration of 'H2OEMIS' in SMCL
!
!> \todo
!>     - ...
!
! ====================================================================


! **********************************************************************
MODULE messy_h2oemis

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! GLOBAL PARAMETERS
  !> Unique submodel string for 'H2OEMIS'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'h2oemis'
  !> Submodel version 
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'
 
END MODULE messy_h2oemis
! **********************************************************************
