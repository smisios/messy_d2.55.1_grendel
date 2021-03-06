MODULE messy_convect_tiedtke_param

  USE messy_main_constants_mem,   ONLY: dp
  USE messy_main_tools,           ONLY: t_reset_par ! fb_mk_20120116

  IMPLICIT NONE
  INTRINSIC :: NULL ! op_pj_20100122

  SAVE

  ! ----------------------------------------------------------------
  !
  ! module *mo_cumulus_flux* - parameters for cumulus massflux scheme
  !
  ! ----------------------------------------------------------------

  REAL(dp) :: entrpen      !    entrainment rate for penetrative convection
  REAL(dp) :: entrscv      !    entrainment rate for shallow convection
  REAL(dp) :: entrmid      !    entrainment rate for midlevel convection
  REAL(dp) :: entrdd       !    entrainment rate for cumulus downdrafts
  REAL(dp) :: centrmax     !
  REAL(dp) :: cmfctop      !    relat. cloud massflux at level above nonbuoyanc
  REAL(dp) :: cmfcmax      !    maximum massflux value allowed for
  REAL(dp) :: cmfcmin      !    minimum massflux value (for safety)
  REAL(dp) :: cmfdeps      !    fractional massflux for downdrafts at lfs
  REAL(dp) :: rhcdd        !    relative saturation in downdrafts
! op_pj_20100122+
!!$  REAL(dp) :: cprcon       !    coefficients for determining conversion
!!$                       !    from cloud water to rain
  REAL(dp), POINTER :: cprcon => NULL()
                       !    coefficients for determining conversion
                       !    from cloud water to rain
! op_pj_20100122-
  INTEGER :: nmctop    !    max. level for cloud base of mid level conv.
  LOGICAL :: lmfpen    !    true if penetrative convection is switched on
  LOGICAL :: lmfscv    !    true if shallow     convection is switched on
  LOGICAL :: lmfmid    !    true if midlevel    convection is switched on
  LOGICAL :: lmfdd     !    true if cumulus downdraft      is switched on
  LOGICAL :: lmfdudv   !    true if cumulus friction       is switched on
! fb_mk_20120116+
  ! Switches to reset the corresponding parameters in Tietke scheme
  TYPE(t_reset_par) :: rset_cmfctop = t_reset_par(.FALSE.,0.0_dp)
  TYPE(t_reset_par) :: rset_cprcon  = t_reset_par(.FALSE.,0.0_dp)
  TYPE(t_reset_par) :: rset_entrscv = t_reset_par(.FALSE.,0.0_dp)
! fb_mk_20120116-
! fb_mk_20140210+
  TYPE(t_reset_par) :: rset_entrpen = t_reset_par(.FALSE.,0.0_dp)
  TYPE(t_reset_par) :: rset_entrmid = t_reset_par(.FALSE.,0.0_dp)
! fb_mk_20140210-

CONTAINS

SUBROUTINE cuparam_init(status, nlev, nlevp1, nvclev, vct,  &
                        nn, lmidatm, p_parallel_io)

  ! Description:
  !
  ! Defines disposable parameters for massflux scheme
  !
  ! Method:
  !
  ! This routine is called from init_tracer from messy_main_control_e5
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, February 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, Jan 1999, subroutine cuparam -> module mo_cumulus_flux
  ! M. Esch, MPI, July 1999, modifications for ECHAM5
  ! U. Schlese, MPI, August 2000, mid level cloud base *nmctop*
  ! H. Tost, MPICH, April 2004, adapted for convection module
  ! 
  ! for more details see file AUTHORS
  ! 

  IMPLICIT NONE

  INTRINSIC :: REAL

  INTEGER, INTENT(IN) :: nlev, nlevp1, nvclev, nn
  REAL(dp),INTENT(IN) :: vct(nvclev*2)
  LOGICAL, INTENT(IN) :: lmidatm, p_parallel_io

  INTEGER, INTENT(INOUT) :: status

! local variables
  REAL(dp)    :: za, zb
  REAL(dp)    :: zph(nlevp1), zp(nlev)

  INTEGER :: jk

  !  Executable Statements 

  status = 0

!-- 1. Specify parameters for massflux-scheme

  entrpen = 1.0E-4_dp !

  entrscv = 1.0E-3_dp ! Average entrainment rate for shallow convection

  entrmid = 1.0E-4_dp ! Average entrainment rate for midlevel convection

  entrdd = 2.0E-4_dp ! Average entrainment rate for downdrafts

  centrmax= 3.E-4_dp !

  cmfcmax = 1.0_dp ! Maximum massflux value allowed for updrafts etc

  cmfcmin = 1.0E-10_dp ! Minimum massflux value (for safety)

  cmfdeps = 0.3_dp ! Fractional massflux for downdrafts at lfs

! mz_jb_20031216+
! The following section was rewritten in order to make the use of new vertical
! model resolutions possible. All parameter values are based on the model
! version echam5.2.02. Changes compared to echam5.1.07_mz10i are
! - new values for cmfctop
! - new values for cprcon

IF (lmidatm) THEN
  SELECT CASE (nn)
  CASE (21)
    IF (20<=nlev .AND. nlev<=38) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (31)
    IF (20<=nlev .AND. nlev<=38) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSEIF (39<=nlev .AND. nlev<=90) THEN
      cmfctop = 0.25_dp+0.05_dp*REAL(nlev-39,dp)/51.0_dp
      cprcon  = 3.0E-4_dp-1.5E-4_dp*REAL(nlev-39,dp)/51.0_dp
    ELSEIF (90<nlev) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE
       status = 1
    END IF
  CASE (42)
    IF (20<=nlev .AND. nlev<=38) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSEIF (39<=nlev .AND. nlev<=90) THEN
      cmfctop = 0.25_dp+0.05_dp*REAL(nlev-39,dp)/51.0_dp
      cprcon  = 3.0E-4_dp-1.5E-4_dp*REAL(nlev-39,dp)/51.0_dp
    ELSEIF (90<nlev) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE
       status = 1
    END IF
  CASE (63)
    IF (20<=nlev .AND. nlev<=38) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSEIF (39<=nlev .AND. nlev<=90) THEN
      cmfctop = 0.25_dp+0.05_dp*REAL(nlev-39,dp)/51.0_dp
      cprcon  = 3.0E-4_dp-2.0E-4_dp*REAL(nlev-39,dp)/51.0_dp
    ELSEIF (90<nlev) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (85)
    IF (20<=nlev .AND. nlev<=38) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSEIF (39<=nlev .AND. nlev<=90) THEN
      cmfctop = 0.25_dp+0.10_dp*REAL(nlev-39,dp)/51.0_dp
      cprcon  = 2.0E-4_dp-1.0E-4_dp*REAL(nlev-39,dp)/51.0_dp
    ELSEIF (90<nlev) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (106)
    IF (20<=nlev .AND. nlev<=38) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSEIF (39<=nlev .AND. nlev<=90) THEN
      cmfctop = 0.25_dp+0.10_dp*REAL(nlev-39,dp)/51.0_dp
      cprcon  = 2.0E-4_dp-1.0E-4_dp*REAL(nlev-39,dp)/51.0_dp
    ELSEIF (90<nlev) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (159)
    IF (20<=nlev .AND. nlev<=38) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSEIF (39<=nlev .AND. nlev<=90) THEN
      cmfctop = 0.25_dp+0.10_dp*REAL(nlev-39,dp)/51.0_dp
      cprcon  = 2.0E-4_dp-1.0E-4_dp*REAL(nlev-39,dp)/51.0_dp
    ELSEIF (90<nlev) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (255)
    IF (39<=nlev .AND. nlev<=90) THEN
      cmfctop = 0.25_dp+0.10_dp*REAL(nlev-39,dp)/51.0_dp
      cprcon  = 2.0E-4_dp-1.0E-4_dp*REAL(nlev-39,dp)/51.0_dp
    ELSEIF (90<nlev) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE DEFAULT
     status = 1
  END SELECT
ELSE
  SELECT CASE (nn)
  CASE (10, 21)
    IF (10<=nlev .AND. nlev<=19) THEN
      cmfctop = 0.10_dp
      cprcon  = 8.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (31)
    IF (10<=nlev .AND. nlev<=19) THEN
      cmfctop = 0.10_dp
      cprcon  = 4.0E-4_dp
    ELSEIF (20<=nlev .AND. nlev<=31) THEN
      cmfctop = 0.10_dp+0.20_dp*REAL(nlev-19,dp)/12.0_dp
      cprcon  = 4.0E-4_dp-2.5E-4_dp*REAL(nlev-19,dp)/12.0_dp
    ELSEIF (31<nlev) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE
       status = 1
    END IF
  CASE (42)
    IF (10<=nlev .AND. nlev<=19) THEN
      cmfctop = 0.12_dp
      cprcon  = 4.0E-4_dp
    ELSEIF (20<=nlev .AND. nlev<=31) THEN
      cmfctop = 0.12_dp+0.18_dp*REAL(nlev-19,dp)/12.0_dp
      cprcon  = 4.0E-4_dp-2.5E-4_dp*REAL(nlev-19,dp)/12.0_dp
    ELSEIF (31<nlev) THEN
      ! op_ck_20070904+
      IF (nlev == 41) then
         cmfctop = 0.22_dp
      ELSE
         cmfctop = 0.30_dp
      ENDIF
      ! op_ck_20070904-
      cprcon  = 1.5E-4_dp
    ELSE
       status = 1
    END IF
  CASE (63)
    IF (10<=nlev .AND. nlev<=19) THEN
      cmfctop = 0.15_dp
      cprcon  = 3.0E-4_dp
    ELSEIF (20<=nlev .AND. nlev<=31) THEN
      cmfctop = 0.15_dp+0.15_dp*REAL(nlev-19,dp)/12.0_dp
      cprcon  = 3.0E-4_dp-2.0E-4_dp*REAL(nlev-19,dp)/12.0_dp
    ELSEIF (31<nlev) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (85)
    IF (10<=nlev .AND. nlev<=19) THEN
      cmfctop = 0.20_dp
      cprcon  = 2.5E-4_dp
    ELSEIF (20<=nlev .AND. nlev<=31) THEN
      cmfctop = 0.20_dp+0.15_dp*REAL(nlev-19,dp)/12.0_dp
      cprcon  = 2.5E-4_dp-1.5E-4_dp*REAL(nlev-19,dp)/12.0_dp
    ELSEIF (31<nlev) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (106)
    IF (10<=nlev .AND. nlev<=19) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.5E-4_dp
    ELSEIF (20<=nlev .AND. nlev<=31) THEN
      cmfctop = 0.25_dp+0.10_dp*REAL(nlev-19,dp)/12.0_dp
      cprcon  = 2.5E-4_dp-1.5E-4_dp*REAL(nlev-19,dp)/12.0_dp
    ELSEIF (31<nlev) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (159)
    IF (10<=nlev .AND. nlev<=19) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.5E-4_dp
    ELSEIF (20<=nlev .AND. nlev<=31) THEN
      cmfctop = 0.25_dp+0.10_dp*REAL(nlev-19,dp)/12.0_dp
      cprcon  = 2.5E-4_dp-1.5E-4_dp*REAL(nlev-19,dp)/12.0_dp
    ELSEIF (31<nlev) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE (255)
    IF (19<=nlev) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    END IF
  CASE DEFAULT
     status = 1
  END SELECT
END IF
! mz_jb_20031216-

! fb_mk_20120116+
  !
  ! Overwrite the initialization of cmfctop, cprcon, ans entrscv with the
  ! the value provided by the namelist.
  IF (rset_cmfctop%l) cmfctop = rset_cmfctop%v
  IF (rset_cprcon%l)  cprcon  = rset_cprcon%v
  IF (rset_entrscv%l) entrscv = rset_entrscv%v
  IF (rset_entrpen%l) entrpen = rset_entrpen%v
  IF (rset_entrmid%l) entrmid = rset_entrmid%v
  IF (p_parallel_io) THEN
     WRITE (*,*) 'SUBROUTINE cuparam_init:'
     WRITE (*,*) ''
     WRITE (*,*) ' lmidatm = ',lmidatm,', nlev = ',nlev,' nn = ',nn
     WRITE (*,*) ''
     WRITE (*,*) ' rset_cmfctop:',rset_cmfctop
     WRITE (*,*) ' relat. cloud massflux at level above nonbuoyanc: cmfctop = ',cmfctop
     WRITE (*,*) ''
     WRITE (*,*) ' rset_cprcon:',rset_cprcon
     WRITE (*,*) ' coefficients for determining conversion from'
     WRITE (*,*) ' cloud water to rain:                             cprcon = ',cprcon
     WRITE (*,*) ''
     WRITE (*,*) ' rset_entrscv:',rset_entrscv
     WRITE (*,*) ' entrainment rate for shallow convection:         entrscv = ',entrscv
     WRITE (*,*) ''
     WRITE (*,*) ' rset_entrpen:',rset_entrpen
     WRITE (*,*) ' entrainment rate for penetrative convection:     entrpen = ',entrpen
     WRITE (*,*) ''
     WRITE (*,*) ' rset_entrmid:',rset_entrmid
     WRITE (*,*) ' entrainment rate for midlevel convection:        entrmid = ',entrmid
  END IF
! fb_mk_20120116-

  ! Next value is relative saturation in downdrafts
  ! but is no longer used ( formulation implies saturation)

  rhcdd = 1.0_dp


! Determine highest level *nmctop* for cloud base of midlevel convection
! assuming nmctop=9 (300 hPa) for the standard 19 level model

!-- half level pressure values, assuming 101320. Pa surface pressure 

  DO jk=1,nlevp1
    za=vct(jk)
    zb=vct(jk+nvclev)
    zph(jk)=za+zb*101320.0_dp
  END DO
!
! -- full level pressure
!
  DO jk = 1, nlev
    zp(jk)=(zph(jk)+zph(jk+1))*0.5_dp
  END DO
!
! -- search for 300 hPa level
!
  DO jk = 1, nlev
    nmctop=jk
    IF(zp(jk).GE.30000.0_dp) EXIT
  END DO
!
  IF (p_parallel_io) THEN
    WRITE (*,*) &
         'max. level for cloud base of mid level convection: nmctop= ',nmctop
  END IF

  RETURN
END SUBROUTINE cuparam_init

!===========================================================================================================

SUBROUTINE cuparam_init_1(status, nlev, nlevp1, nvclev, vct,  &
                        nn, lmidatm, p_parallel_io, lcouple, lipcc)

  ! Description:
  !
  ! Defines disposable parameters for massflux scheme
  !
  ! Method:
  !
  ! This routine is called from *iniphy*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, February 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, Jan 1999, subroutine cuparam -> module mo_cumulus_flux
  ! M. Esch, MPI, July 1999, modifications for ECHAM5
  ! U. Schlese, MPI, August 2000, mid level cloud base *nmctop*
  ! 
  ! for more details see file AUTHORS
  ! 

  IMPLICIT NONE

  INTRINSIC :: REAL

  INTEGER, INTENT(IN) :: nlev, nlevp1, nvclev, nn
  REAL(dp),INTENT(IN) :: vct(nvclev*2)
  LOGICAL, INTENT(IN) :: lmidatm, p_parallel_io, lcouple, lipcc

  INTEGER, INTENT(INOUT) :: status

! local variables
  REAL(dp)    :: za, zb, zph(nlevp1), zp(nlev)
  INTEGER :: jk

  !  Executable Statements 

  status = 0

!-- 1. Specify parameters for massflux-scheme

  IF ( nlev == 11 ) THEN
     entrpen = 1.0E-4_dp !
     entrscv = 3.0E-4_dp !
  ELSE
     entrpen = 1.0E-4_dp !
     entrscv = 1.0E-3_dp !
  ENDIF

!
! coupled model:
  IF(lcouple .OR. lipcc) entrscv = 3.0E-4_dp  ! mz_ap_20090519 

  entrmid = 1.0E-4_dp ! Average entrainment rate for midlevel convection

  entrdd = 2.0E-4_dp ! Average entrainment rate for downdrafts

  centrmax= 3.E-4_dp !

  cmfcmax = 1.0_dp ! Maximum massflux value allowed for updrafts etc

  cmfcmin = 1.E-10_dp ! Minimum massflux value (for safety)

  cmfdeps = 0.3_dp ! Fractional massflux for downdrafts at lfs

!
!                19 Level, no middle atmosphere
!
  IF (nlev == 11 .AND. .NOT. lmidatm) THEN
#ifndef MESSY
     IF (nn == 21) THEN
#else
     IF ((nn == 21) .OR. (nn == 10)) THEN
#endif
      cmfctop = 0.3_dp
      cprcon  = 1.0E-3_dp
    ELSE IF (nn == 31) THEN
      cmfctop = 0.3_dp
      cprcon  = 1.0E-3_dp
    ELSE
       status = 1
    ENDIF
  ELSE IF (nlev == 19  .AND. .NOT. lmidatm) THEN
#ifndef MESSY
    IF (nn == 21) THEN
#else
    IF ((nn == 21) .OR. (nn == 10)) THEN
#endif
      cmfctop = 0.1_dp 
      cprcon  = 8.0E-4_dp
    ELSE IF (nn == 31) THEN
      cmfctop = 0.1_dp 
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.20_dp  ! mz_ap_20090519 
      cprcon  = 4.0E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.12_dp
      cprcon  = 4.0E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.15_dp 
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.20_dp
      cprcon  = 2.5E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.5E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.5E-4_dp
    ELSE
       status = 1
    ENDIF
!
!                31 Level, no middle atmosphere
!
  ELSE IF (nlev == 31 .AND. .NOT. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
! op_ck_20070907+
!
!                41 Level (DLR)
!
 ELSE IF (nlev == 41) THEN
    IF (nn == 42) THEN
       cmfctop = 0.22_dp
       cprcon  = 1.5E-4_dp
    ELSE
       status = 1
    ENDIF
! op_ck_20070907-
!
!                60 Level
!
  ELSE IF (nlev == 60) THEN
    IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
!
!                39 Level, middle atmosphere
!
  ELSE IF (nlev == 39 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.25_dp  
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.25_dp  
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.25_dp 
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.25_dp 
      cprcon  = 2.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.25_dp
      cprcon  = 2.0E-4_dp
    ELSE
       status = 1
    ENDIF
!
!                47 Level, middle atmosphere
!
  ELSE IF (nlev == 47 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
      ! mz_bk_20110429+
      ! coupled model:
      ! mz_bk_20120820+
      !!$ IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      IF(lcouple .OR. lipcc) THEN
         cmfctop = 0.21_dp
         cprcon  = 3.0E-4_dp
      END IF
      ! mz_bk_20120820-
      ! mz_bk_20110429-
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
!
!                49 Level, middle atmosphere
!
  ELSE IF (nlev == 49 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
!
!                87 Level, middle atmosphere
!
  ELSE IF (nlev == 87 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
!
!                90 Level, middle atmosphere
!
  ELSE IF (nlev == 90 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
!
!                95 Level, middle atmosphere
!
  ELSE IF (nlev == 95 .AND. lmidatm) THEN
    IF (nn == 31) THEN
      cmfctop = 0.30_dp  
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 42) THEN
      cmfctop = 0.30_dp
      cprcon  = 1.5E-4_dp
    ELSE IF (nn == 63) THEN
      cmfctop = 0.30_dp
    ! coupled model:
      IF(lcouple .OR. lipcc) cmfctop = 0.22_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 85) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
!
!                191 Level, middle atmosphere
!
  ELSE IF (nlev == 191 .AND. lmidatm) THEN
    IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
! ka_sv_20160406+
!                74 Level, middle atmosphere and extended above
!
  ELSE IF (nlev == 74 .AND. lmidatm) THEN
    IF (nn==42) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
!                84 Level, middle atmosphere and extended above
!
   ELSE IF (nlev == 84 .AND. lmidatm) THEN
    IF (nn==42) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
!                94 Level, middle atmosphere and extended above
!
   ELSE IF (nlev == 94 .AND. lmidatm) THEN
    IF (nn==42) THEN
      cmfctop = 0.25_dp
      cprcon  = 3.0E-4_dp
    ELSE IF (nn == 106) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 159) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 255) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE IF (nn == 319) THEN
      cmfctop = 0.35_dp
      cprcon  = 1.0E-4_dp
    ELSE
       status = 1
    ENDIF
! ka_sv_20180112-
  ELSE
       status = 1
  ENDIF

! op_pj_20140204+
  !
  ! Overwrite the initialization of cmfctop, cprcon, ans entrscv with the
  ! the value provided by the namelist.
  IF (rset_cmfctop%l) cmfctop = rset_cmfctop%v
  IF (rset_cprcon%l)  cprcon  = rset_cprcon%v
  IF (rset_entrscv%l) entrscv = rset_entrscv%v
  IF (rset_entrpen%l) entrpen = rset_entrpen%v
  IF (rset_entrmid%l) entrmid = rset_entrmid%v
  IF (p_parallel_io) THEN
     WRITE (*,*) 'SUBROUTINE cuparam_init_1:'
     WRITE (*,*) ''
     WRITE (*,*) ' lmidatm = ',lmidatm,', nlev = ',nlev,' nn = ',nn
     WRITE (*,*) ''
     WRITE (*,*) ' rset_cmfctop:',rset_cmfctop
     WRITE (*,*) ' relat. cloud massflux at level above nonbuoyanc: cmfctop = ',cmfctop
     WRITE (*,*) ''
     WRITE (*,*) ' rset_cprcon:',rset_cprcon
     WRITE (*,*) ' coefficients for determining conversion from'
     WRITE (*,*) ' cloud water to rain:                             cprcon = ',cprcon
     WRITE (*,*) ''
     WRITE (*,*) ' rset_entrscv:',rset_entrscv
     WRITE (*,*) ' entrainment rate for shallow convection:         entrscv = ',entrscv
     WRITE (*,*) ''
     WRITE (*,*) ' rset_entrpen:',rset_entrpen
     WRITE (*,*) ' entrainment rate for penetrative convection:     entrpen = ',entrpen
     WRITE (*,*) ''
     WRITE (*,*) ' rset_entrmid:',rset_entrmid
     WRITE (*,*) ' entrainment rate for midlevel convection:        entrmid = ',entrmid
  END IF
! op_pj_20140204-

  ! Next value is relative saturation in downdrafts
  ! but is no longer used ( formulation implies saturation)

  rhcdd = 1.0_dp


! Determine highest level *nmctop* for cloud base of midlevel convection
! assuming nmctop=9 (300 hPa) for the standard 19 level model

!-- half level pressure values, assuming 101320. Pa surface pressure 

  DO jk=1,nlevp1
    za=vct(jk)
    zb=vct(jk+nvclev)
    zph(jk)=za+zb*101320.0_dp
  END DO
!
! -- full level pressure
!
  DO jk = 1, nlev
    zp(jk)=(zph(jk)+zph(jk+1))*0.5_dp
  END DO
!
! -- search for 300 hPa level
!
  DO jk = 1, nlev
    nmctop=jk
    IF(zp(jk).GE.30000.0_dp) EXIT
  END DO
!
  IF (p_parallel_io) THEN
    WRITE (*,*) &
         'max. level for cloud base of mid level convection: nmctop= ',nmctop
  END IF

  RETURN
END SUBROUTINE cuparam_init_1

!******************************************************************************
! op_mm_20140327 +
! DEFINING PARAMETER FOR COSMO
! USING SAME PARAMETERS AS IN src_conv_tiedtke.f90 

SUBROUTINE cuparam_init_2(status, nlev,hyai,hybi,p_parallel_io) 


  ! Description:
  !
  ! Defines disposable parameters for massflux scheme
  !
  ! Method:
  !
  ! This routine is called from *iniphy*
  !
  ! Authors:
  !
  ! M. Tiedtke, ECMWF, February 1989, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A. Rhodin, MPI, Jan 1999, subroutine cuparam -> module mo_cumulus_flux
  ! M. Esch, MPI, July 1999, modifications for ECHAM5
  ! U. Schlese, MPI, August 2000, mid level cloud base *nmctop*
  ! 
  ! for more details see file AUTHORS
  ! 

  IMPLICIT NONE

  INTRINSIC :: REAL

  INTEGER, INTENT(IN) :: nlev
  REAL(dp),INTENT(IN) :: hyai(nlev+1), hybi(nlev+1)
  LOGICAL, INTENT(IN) ::  p_parallel_io

  INTEGER, INTENT(INOUT) :: status

! local variables
  REAL(dp)    :: za, zb, zph(nlev+1), zp(nlev)
  INTEGER :: jk

  !  Executable Statements 

  status = 0

  entrpen       = 0.00010
  entrscv  = 0.00030  !definied in cosmo/data_convection.f90
  entrmid       = 0.00010  
  entrdd        = 0.00020
  centrmax= 3.E-4_dp   ! is not used in cosmo. Defining same as for echam
  cmfcmax = 1.0_dp 
  cmfcmin       = 1.E-10 
  cmfdeps       = 0.3 
  cmfctop       = 0.33
  cprcon        = 0.0002


  ! Next value is relative saturation in downdrafts
  ! but is no longer used ( formulation implies saturation)

  rhcdd = 1.0_dp


! Determine highest level *nmctop* for cloud base of midlevel convection
! assuming nmctop=9 (300 hPa) for the standard 19 level model

!-- half level pressure values, assuming 101320. Pa surface pressure 

  DO jk=1,nlev+1
    za=hyai(jk)
    zb=hybi(jk)
    zph(jk)=za+zb*101320.0_dp
  END DO
!
! -- full level pressure
!
  DO jk = 1, nlev
    zp(jk)=(zph(jk)+zph(jk+1))*0.5_dp
  END DO
!
! -- search for 300 hPa level
!
  DO jk = 1, nlev
    nmctop=jk
    IF(zp(jk).GE.30000.0_dp) EXIT
  END DO
!
  IF (p_parallel_io) THEN
    WRITE (*,*) &
         'max. level for cloud base of mid level convection: nmctop= ',nmctop
  END IF

  RETURN
END SUBROUTINE cuparam_init_2
!===========================================================================================================
END MODULE messy_convect_tiedtke_param
