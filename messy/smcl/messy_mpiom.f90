
MODULE messy_mpiom

  ! MESSy
  USE messy_main_constants_mem,  ONLY: DP, SP, WP, STRLEN_MEDIUM

  IMPLICIT NONE 
  PUBLIC

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'mpiom'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.0'

  PUBLIC :: dp, sp, STRLEN_MEDIUM

  ! CTRL-NAMELIST PARAMETERS

  INTEGER*8, PUBLIC :: IDATE,III
  INTEGER  , PUBLIC :: NACTYEAR
  REAL(dp) , PUBLIC :: ABACK,CAULAPTS,CAULAPUV,CWT,DBACK,GFDL_DIFF
#if defined MPIOM_13B
  REAL(dp) , PUBLIC :: CRELSAL,CRELTEM
#elif defined MPIOM_2000
   REAL(dp) , PUBLIC :: cwa
#endif
  REAL(dp) , PUBLIC :: PIBOG,RRELSAL,RRELTEM,DIST
  REAL(dp) , PUBLIC :: CDZW(500)    
#ifdef MPIOM_13B
  INTEGER  , PUBLIC :: NDTDAY
#endif
  INTEGER  , PUBLIC :: I,ICONVA,IWIRB,J,K,M,N,NANF,ENDDAY
  INTEGER  , PUBLIC :: ICOU,IJJU,JMANF,JMEND,JMM,JMMM,LDTRUN
  INTEGER  , PUBLIC :: IERR, LDTDAY, LEN_FEB, LMON1, LMON2, LREAD
  INTEGER  , PUBLIC :: LREADMAX, MONMON, NACYEAR, NDTYEAR, IMAL
  INTEGER  , PUBLIC :: L,IREADC
  REAL*8   , PUBLIC :: ttts, tttr, ttt
  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC :: GRID = ''
  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC :: VLEVELS = ''
  INTEGER, PUBLIC                      :: delta_time_mpiom
  INTEGER, PUBLIC                      :: NPROCA
  INTEGER, PUBLIC                      :: NPROCB
  INTEGER, PUBLIC                      :: IE_G_nml
  INTEGER, PUBLIC                      :: JE_G_nml
  LOGICAL, PUBLIC                      :: L_TIDES 
  LOGICAL, PUBLIC                      :: L_COUPLING 
  LOGICAL, PUBLIC                      :: L_HAMOCC_COUPLING 
  LOGICAL, PUBLIC                      :: L_CHECK 

#ifdef MPIOM_2000
  ! time step variables
  INTEGER LDTYEAR, LDTMONTH, ldt
  INTEGER :: NYEARS, NMONTS, NDAYS

  INTEGER :: ll, lmont, ntracerloop
  INTEGER :: ihalo_sor, imod
  ! write grid information back to file
  LOGICAL :: lgridinfo = .FALSE.
  ! check floating-point operations for exceptions
  LOGICAL :: fp_tracing_enabled = .FALSE.
  ! apply sea level correction?
  LOGICAL :: lzo_correct = .TRUE.
  LOGICAL ::lundelayed_momentum_advection = .FALSE.
#endif

  PUBLIC :: mpiom_read_nml_ctrl

CONTAINS

   ! ---------------------------------------------------------------------------
 
   SUBROUTINE mpiom_read_nml_ctrl(status, iou)
 
     !  MPIOM MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
     !
     ! read namelist for 'coupling' to ECHAM5
     !
     ! Author: Pozzer Andrea, MPICH, Oct 2004
 
 
     USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
 
     IMPLICIT NONE
 
     ! I/O
     INTEGER, INTENT(OUT) :: status     ! error status
     INTEGER, INTENT(IN)  :: iou        ! I/O unit
 
   NAMELIST /CTRL/ GRID, VLEVELS, delta_time_mpiom, NPROCA,NPROCB, IE_G_nml, JE_G_nml, &
                    L_TIDES, L_COUPLING, L_HAMOCC_COUPLING, L_CHECK   !declared in mo_param1.f90
 
     ! LOCAL
     CHARACTER(LEN=*), PARAMETER :: substr='mpiom_read_nml_ctrl'
     LOGICAL              :: lex      ! file exists ?
     INTEGER              :: fstat    ! file status
 
     status = 1
 
!----------------------------------------------------------------------
! DEFAULT PARAMETER SETTINGS - CAN BE OVERWRITEN BY THE NAMELIST
     delta_time_mpiom = 0
     IE_G_nml = 0
     JE_G_nml = 0
     L_TIDES          =.TRUE.
     L_COUPLING       =.TRUE.
     L_HAMOCC_COUPLING=.TRUE.
     L_CHECK          =.FALSE.
!----------------------------------------------------------------------
     
     CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
     IF (.not.lex) RETURN    ! <modstr>.nml does not exist
 
     READ(iou, NML=CTRL, IOSTAT=fstat)

     CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
     IF (fstat /= 0) RETURN  ! error while reading namelist
 
     CALL read_nml_close(substr, iou, modstr)
 
     status = 0 ! NO ERROR
 
   END SUBROUTINE mpiom_read_nml_ctrl
 
   ! ---------------------------------------------------------------------------

END MODULE
