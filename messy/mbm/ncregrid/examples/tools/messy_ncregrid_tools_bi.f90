! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_TOOLS_BI
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, October 2002
! ******************************************************************

#ifndef ICON

  USE messy_main_constants_mem, ONLY: dp
  USE messy_main_blather_bi,    ONLY: error_bi
  ! MESSy
  USE messy_main_timer_event,   ONLY: time_event, io_time_event  &
                                    , TIME_INC_MONTHS, TRIG_EXACT
  USE messy_ncregrid_tools,     ONLY: NCRGCNT, LMODE2DH

  IMPLICIT NONE
  PRIVATE

  ! MAX. NUMBER OF RGTEVENTS PER SUB-MODEL
  INTEGER, PARAMETER         :: NMAXRGTE     = 500
  INTEGER, PARAMETER, PUBLIC :: RGTMAXACTSTR = 200

  ! READ IN OPTION
  INTEGER, PARAMETER, PUBLIC :: RGREAD_NCVAR =1 ! USE RGTOOL_BI_READ_NCVAR
  INTEGER, PARAMETER, PUBLIC :: RGREAD_NCFILE=2 ! USE RGTOOL_BI_READ_NCFILE

  ! USED FOR NAMELIST-INPUT (ONLY INTERNALLY)
  TYPE RGTEVENT_IO
     TYPE(io_time_event)         :: evt = &
          io_time_event(1, TIME_INC_MONTHS,TRIG_EXACT,0)
     TYPE(NCRGCNT)               :: cnt = &
          NCRGCNT('', -1, -1, -1, -1)
     CHARACTER(LEN=RGTMAXACTSTR) :: act = ''
  END TYPE RGTEVENT_IO

  ! FOR USER INTERFACE
  TYPE RGTEVENT
     TYPE(RGTEVENT_IO) :: io
     TYPE(time_event)  :: event
  END TYPE RGTEVENT

  ! PUBLIC, BUT NOT DIRECT USER INTERFACE
  ! NEEDED TO ACCESS SUB-STRUCTURES OF RGTEVENT
  PUBLIC :: RGTEVENT_IO ! (only public, because component of RGTEVENT !)
  PUBLIC :: NCRGCNT     ! (only public, because component of RGTEVENT !)

  ! USER INTERFACE
  ! - DATA IMPORT
  ! op_pj_20120510+
  INTERFACE RGTOOL_BI_READ_NCVAR
     MODULE PROCEDURE RGTOOL_BI_READ_NCVAR_4D
     MODULE PROCEDURE RGTOOL_BI_READ_NCVAR_5D
  END INTERFACE
  ! op_pj_20120510-
  PUBLIC  :: RGTOOL_BI_READ_NCVAR  ! READ ONE FIELD WITH NCREGRID
  ! op_pj_20120510+
  INTERFACE RGTOOL_BI_READ_NCFILE
     MODULE PROCEDURE RGTOOL_BI_READ_NCFILE_4D
     MODULE PROCEDURE RGTOOL_BI_READ_NCFILE_5D
  END INTERFACE
  ! op_pj_20120510-
  PUBLIC  :: RGTOOL_BI_READ_NCFILE ! READ ONE FILE WITH NCREGRID
  PUBLIC  :: RGTEVENT_READ         ! READ DATA FIELD(S) ON EVENT
  ! - EVENT/COUNTER HANDLING
  PUBLIC  :: RGTEVENT              ! RG-TRIGGER EVENT STRUCTURE
  PUBLIC  :: RGTEVENT_INIT_NML     ! READ RGT-EVENT INFORMATION FROM NAMELIST
  PUBLIC  :: RGTEVENT_INIT_ONL     ! INIT RGT-EVENT INFORMATION ONLINE
  PUBLIC  :: RGTEVENT_STATUS       ! GET STATUS OF NAMED/INDEXED EVENT AND
                                   ! UPDATE RESTART FILE
  ! - COUNTER RESTART HANDLING
  PUBLIC  :: messy_ncregrid_write_restart
  PUBLIC  :: messy_ncregrid_read_restart
  PUBLIC  :: messy_ncregrid_free_memory
  ! - SWITCH MODI (BASEMODEL DEPENDENT)
  PUBLIC  :: messy_ncregrid_initialize

  INTERFACE RGTEVENT_READ
     MODULE PROCEDURE RGTEVENT_READ_2D
     MODULE PROCEDURE RGTEVENT_READ_3D
     MODULE PROCEDURE RGTEVENT_READ_4D
  END INTERFACE

  ! PRIVATE HELPER ROUTINES
!  PRIVATE :: RGTEVENT_INIT        ! INITIALIZE ONE RGT-EVENT FROM I/O
!  PRIVATE :: RGTEVENT_STAT        ! CHECK STATUS OF ONE RGT-EVENT
!  PRIVATE :: RGTEVENT_INDEX       ! CALCULATES EVENT-INDEX FROM NAME

CONTAINS

! ####### PUBLIC ROUTINES ##########################################

! ------------------------------------------------------------------
SUBROUTINE messy_ncregrid_initialize

#ifndef HOMMESE
  LMODE2DH = .TRUE.
#else
  LMODE2DH = .FALSE.
#endif

END SUBROUTINE messy_ncregrid_initialize
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTOOL_BI_READ_NCVAR_4D(modstr, vname, t, dat       &
                               , lrg, lrgx, lrgy, lrgz, lok &  ! OPTIONAL
                               , hyam, hybm, p0, ps         &  ! OPTIONAL
                               , hyai, hybi                 &  ! OPTIONAL
                               , latm, lonm, lati, loni     &  ! OPTIONAL
                               )

  ! PERFORMS ONE NCREGRID STEP FOR ONE FIELD
  ! AND RETURNS DATA AS 4D ARRAY (x,z,n,y),
  ! AND OPTIONALLY GRID-STRUCTURE AS ARRAYS
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast &
                                 , dcg, scatter_gp
  USE messy_main_data_bi,    ONLY: ngpblks, nproma
  ! MESSY
  USE messy_main_tools,      ONLY: find_next_free_unit
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf, ONLY: ncvar, init_ncvar
  USE messy_ncregrid_geohyb, ONLY: geohybgrid, init_geohybgrid
  USE messy_ncregrid_tools,  ONLY: rgtool_read_ncvar, rgtool_convert &
                                 , rgtool_g2c

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, PRESENT, SIZE, TRIM

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)             :: modstr    ! calling module
  CHARACTER(LEN=*), INTENT(IN)             :: vname     ! name of variable
  INTEGER,          INTENT(IN)             :: t         ! netCDF time step
  REAL(dp), DIMENSION(:,:,:,:), POINTER    :: dat       ! (local) data field
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrg       ! regrid really ?
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgx      ! regrid in x
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgy      ! regrid in y
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgz      ! regrid in z
  LOGICAL,          INTENT(OUT),  OPTIONAL :: lok       ! OK?
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyam  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybm  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: p0    ! reference pressure
  REAL(dp), DIMENSION(:,:), POINTER,  OPTIONAL :: ps    ! (local) surf. press.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyai  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybi  ! hybrid-B-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: latm  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lonm  ! longitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lati  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: loni  ! longitude

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_BI_READ_NCVAR'
  INTEGER                            :: iou       ! logical I/O unit
!!$  INTEGER                            :: status    ! memory status
  TYPE(ncvar)                        :: var       ! nc-variable
  TYPE(geohybgrid)                   :: grid      ! grid-structure
  REAL(dp), DIMENSION(:,:,:,:), POINTER  :: gdat  ! (global) data field
  REAL(dp), DIMENSION(:,:), POINTER      :: gps   ! surface pressure
  !
  LOGICAL                            :: lsc       ! scatter results ?
  LOGICAL                            :: llrg, llrgx, llrgy, llrgz, llok
  INTEGER                            :: klat
  INTEGER                            :: klon
  INTEGER                            :: klev
  INTEGER                            :: kparam
  INTEGER                            :: s1d        ! size of 1D field

  ! INITIALIZE
  ! I/O
  IF (ASSOCIATED(dat)) THEN
!!$     DEALLOCATE(dat, STAT=status)
!!$     CALL ERRMSG(substr, status, 1)
     DEALLOCATE(dat)
  END IF
  NULLIFY(dat)
  !
  IF (PRESENT(ps)) THEN
     IF (ASSOCIATED(ps)) THEN
!!$        DEALLOCATE(ps, STAT=status)
!!$        CALL ERRMSG(substr, status, 2)
        DEALLOCATE(ps)
     END IF
     NULLIFY(ps)
  END IF
  !
  IF (PRESENT(hyam)) THEN
     IF (ASSOCIATED(hyam)) THEN
!!$        DEALLOCATE(hyam, STAT=status)
!!$        CALL ERRMSG(substr, status, 3)
        DEALLOCATE(hyam)
     END IF
     NULLIFY(hyam)
  END IF
  !
  IF (PRESENT(hybm)) THEN
     IF (ASSOCIATED(hybm)) THEN
!!$        DEALLOCATE(hybm, STAT=status)
!!$        CALL ERRMSG(substr, status, 4)
        DEALLOCATE(hybm)
     END IF
     NULLIFY(hybm)
  END IF
  !
  IF (PRESENT(hyai)) THEN
     IF (ASSOCIATED(hyai)) THEN
!!$        DEALLOCATE(hyai, STAT=status)
!!$        CALL ERRMSG(substr, status, 5)
        DEALLOCATE(hyai)
     END IF
     NULLIFY(hyai)
  END IF
  !
  IF (PRESENT(hybi)) THEN
     IF (ASSOCIATED(hybi)) THEN
!!$        DEALLOCATE(hybi, STAT=status)
!!$        CALL ERRMSG(substr, status, 6)
        DEALLOCATE(hybi)
     END IF
     NULLIFY(hybi)
  END IF
  !
  IF (PRESENT(p0)) THEN
     IF (ASSOCIATED(p0)) THEN
!!$        DEALLOCATE(p0, STAT=status)
!!$        CALL ERRMSG(substr, status, 7)
        DEALLOCATE(p0)
     END IF
     NULLIFY(p0)
  END IF
  !
  IF (PRESENT(latm)) THEN
     IF (ASSOCIATED(latm)) THEN
!!$        DEALLOCATE(latm, STAT=status)
!!$        CALL ERRMSG(substr, status, 8)
        DEALLOCATE(latm)
     END IF
     NULLIFY(latm)
  END IF
  !
  IF (PRESENT(lonm)) THEN
     IF (ASSOCIATED(lonm)) THEN
!!$        DEALLOCATE(lonm, STAT=status)
!!$        CALL ERRMSG(substr, status, 9)
        DEALLOCATE(lonm)
     END IF
     NULLIFY(lonm)
  END IF
  !
  IF (PRESENT(lati)) THEN
     IF (ASSOCIATED(lati)) THEN
!!$        DEALLOCATE(lati, STAT=status)
!!$        CALL ERRMSG(substr, status, 10)
        DEALLOCATE(lati)
     END IF
     NULLIFY(lati)
  END IF
  !
  IF (PRESENT(loni)) THEN
     IF (ASSOCIATED(loni)) THEN
!!$        DEALLOCATE(loni, STAT=status)
!!$        CALL ERRMSG(substr, status, 11)
        DEALLOCATE(loni)
     END IF
     NULLIFY(loni)
  END IF
  !
  ! LOCAL
  IF (PRESENT(lrg)) THEN
     llrg = lrg
  ELSE
     llrg = .true.
  END IF
  IF (PRESENT(lrgx)) THEN
     llrgx = lrgx
  ELSE
     llrgx = .true.
  END IF
  IF (PRESENT(lrgy)) THEN
     llrgy = lrgy
  ELSE
     llrgy = .true.
  END IF
  IF (PRESENT(lrgz)) THEN
     llrgz = lrgz
  ELSE
     llrgz = .true.
  END IF
  !
  NULLIFY(gdat)       ! (global) data field
  NULLIFY(gps)        ! (global) surface pressure field

#ifdef HOMMESE
  llrg = .false.
  llrgx = .false.
  llrgy = .false.
  llrgz = .false.
#endif

  IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL RGTOOL_READ_NCVAR(iou, TRIM(modstr)//'.nml', vname, t, var &
                           ,grid, llrg, llrgx, llrgy, llrgz, llok)
     IF (llok) THEN
        CALL RGTOOL_CONVERT(var, gdat, grid, order='xzny')
        klon = SIZE(gdat, 1)
        klev = SIZE(gdat, 2)
        kparam = SIZE(gdat, 3)
        klat = SIZE(gdat, 4)
        CALL RGTOOL_G2C(grid, hyam, hybm, p0, gps      &
                        ,hyai, hybi                    &
                        ,latm, lonm                    &
                        ,lati, loni, klat, klon, klev)
     END IF
     CALL INIT_NCVAR(var)
     CALL INIT_GEOHYBGRID(grid)
  END IF ! p_parallel_io

  ! BROADCAST SUCCESS
  CALL p_bcast(llok, p_io)
  IF (PRESENT(lok)) THEN
     lok = llok
  END IF

  ! BROADCST DIMENSIONS
  IF (llok) THEN
     CALL p_bcast(klon, p_io)
     CALL p_bcast(klev, p_io)
     CALL p_bcast(klat, p_io)
     CALL p_bcast(kparam, p_io)
  END IF

  ! NOTE: SCATTERING THE RESULTS (dat) TO THE PROCESSORS IS ONLY POSSIBLE
  !       IF REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy.AND.lrgz)
#ifndef HOMMESE
  lsc = llrg.AND.llok.AND.llrgx.AND.llrgy.AND.llrgz
#else
  lsc = llok
#endif
  !
  IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
     ! ALLOCATE SPACE FOR LOCAL DATA FIELD
!!$     ALLOCATE(dat(nproma,klev,kparam,ngpblks), STAT=status)
!!$     CALL ERRMSG(substr, status, 12)
     ALLOCATE(dat(nproma,klev,kparam,ngpblks))
     ! SCATTER GLOBAL DATA TO LOCAL FIELDS
     CALL scatter_gp(gdat, dat, dcg)
  ELSE            ! COPY AND BROADCAST ENTIRE FIELD
     IF (llok) THEN
        ! ALLOCATE SPACE FOR DATA FIELD
!!$        ALLOCATE(dat(klon,klev,kparam,klat),STAT=status)
!!$        CALL ERRMSG(substr, status, 13)
        ALLOCATE(dat(klon,klev,kparam,klat))
        ! COPY FIRST ON I/O-PROCESSOR
        IF (p_parallel_io) THEN
           dat(:,:,:,:) = gdat(:,:,:,:)
        END IF
        ! BROADCAST RESULT
        CALL p_bcast(dat, p_io)
     END IF ! llok
  END IF ! lsc

  ! NOTE: SCATTERING THE RESULTS (ps) TO THE PROCESSORS IS ONLY POSSIBLE
  !       IF HORIZONTAL REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy)
  IF (PRESENT(ps)) THEN
#ifndef HOMMESE
     lsc = llrg.AND.llok.AND.llrgx.AND.llrgy
#else
     lsc = llok
#endif
     !
     IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
        ! ALLOCATE SPACE FOR LOCAL DATA FIELD
!!$        ALLOCATE(ps(nproma,ngpblks), STAT=status)
!!$        CALL ERRMSG(substr, status, 14)
        ALLOCATE(ps(nproma,ngpblks))
        ! SCATTER GLOBAL DATA TO LOCAL FIELDS
        CALL scatter_gp(gps, ps, dcg)
     ELSE
        IF (llok) THEN
           ! ALLOCATE SPACE FOR DATA FIELD
!!$           ALLOCATE(ps(klon,klat),STAT=status)
!!$           CALL ERRMSG(substr, status, 15)
           ALLOCATE(ps(klon,klat))
           IF (p_parallel_io) THEN
              ps(:,:) = gps(:,:)
           END IF
           ! BROADCAST RESULT
           CALL p_bcast(ps , p_io)
        END IF ! llok
     END IF ! lsc
  END IF ! ps present

  ! 1D ARRAYS NEVER SCATTERED
  IF (llok) THEN
     IF (PRESENT(hyam)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hyam)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hyam(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 16)
           ALLOCATE(hyam(s1d))
        END IF
        CALL p_bcast(hyam ,p_io)
     END IF
     !
     IF (PRESENT(hybm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hybm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hybm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 17)
           ALLOCATE(hybm(s1d))
        END IF
        CALL p_bcast(hybm ,p_io)
     END IF
     !
     IF (PRESENT(hyai)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hyai)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hyai(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 18)
           ALLOCATE(hyai(s1d))
        END IF
        CALL p_bcast(hyai ,p_io)
     END IF
     !
     IF (PRESENT(hybi)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hybi)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hybi(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 19)
           ALLOCATE(hybi(s1d))
        END IF
        CALL p_bcast(hybi ,p_io)
     END IF
     !
     IF (PRESENT(p0)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(p0)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(p0(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 20)
           ALLOCATE(p0(s1d))
        END IF
        CALL p_bcast(p0 ,p_io)
     END IF
     !
     IF (PRESENT(latm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(latm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(latm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 21)
           ALLOCATE(latm(s1d))
        END IF
        CALL p_bcast(latm ,p_io)
     END IF
     !
     IF (PRESENT(lonm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(lonm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(lonm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 22)
           ALLOCATE(lonm(s1d))
        END IF
        CALL p_bcast(lonm ,p_io)
     END IF
     !
     IF (PRESENT(lati)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(lati)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(lati(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 23)
           ALLOCATE(lati(s1d))
        END IF
        CALL p_bcast(lati ,p_io)
     END IF
     !
     IF (PRESENT(loni)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(loni)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(loni(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 24)
           ALLOCATE(loni(s1d))
        END IF
        CALL p_bcast(loni ,p_io)
     END IF
  END IF

  ! CLEAN UP: FREE MEMORY OF GLOBAL FIELDS
  IF (llok) THEN
     IF (ASSOCIATED(gdat)) THEN
!!$        DEALLOCATE(gdat, STAT=status)
!!$        CALL ERRMSG(substr, status, 25)
        DEALLOCATE(gdat); NULLIFY(gdat)
     END IF
!!$     IF (PRESENT(ps)) THEN ! op_pj_20130806 memory leak
        IF (ASSOCIATED(gps)) THEN
!!$           DEALLOCATE(gps, STAT=status)
!!$           CALL ERRMSG(substr, status, 26)
           DEALLOCATE(gps) ; NULLIFY(gps)
        END IF
!!$     END IF ! op_pj_20130806 memory leak
  END IF

END SUBROUTINE RGTOOL_BI_READ_NCVAR_4D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTOOL_BI_READ_NCVAR_5D(modstr, vname, t, dat       &
                               , lrg, lrgx, lrgy, lrgz, lok &  ! OPTIONAL
                               , hyam, hybm, p0, ps         &  ! OPTIONAL
                               , hyai, hybi                 &  ! OPTIONAL
                               , latm, lonm, lati, loni     &  ! OPTIONAL
                               )

  ! PERFORMS ONE NCREGRID STEP FOR ONE FIELD
  ! AND RETURNS DATA AS 4D ARRAY (x,z,n,y),
  ! AND OPTIONALLY GRID-STRUCTURE AS ARRAYS
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast &
                                 , dcg, scatter_gp
  USE messy_main_data_bi,    ONLY: ngpblks, nproma
  ! MESSY
  USE messy_main_tools,      ONLY: find_next_free_unit
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf, ONLY: ncvar, init_ncvar
  USE messy_ncregrid_geohyb, ONLY: geohybgrid, init_geohybgrid
  USE messy_ncregrid_tools,  ONLY: rgtool_read_ncvar, rgtool_convert &
                                 , rgtool_g2c

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, PRESENT, SIZE, TRIM

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)             :: modstr    ! calling module
  CHARACTER(LEN=*), INTENT(IN)             :: vname     ! name of variable
  INTEGER,          INTENT(IN)             :: t         ! netCDF time step
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER  :: dat       ! (local) data field
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrg       ! regrid really ?
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgx      ! regrid in x
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgy      ! regrid in y
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgz      ! regrid in z
  LOGICAL,          INTENT(OUT),  OPTIONAL :: lok       ! OK?
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyam  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybm  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: p0    ! reference pressure
  REAL(dp), DIMENSION(:,:), POINTER,  OPTIONAL :: ps    ! (local) surf. press.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyai  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybi  ! hybrid-B-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: latm  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lonm  ! longitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lati  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: loni  ! longitude

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_BI_READ_NCVAR'
  INTEGER                            :: iou       ! logical I/O unit
!!$  INTEGER                            :: status    ! memory status
  TYPE(ncvar)                        :: var       ! nc-variable
  TYPE(geohybgrid)                   :: grid      ! grid-structure
  REAL(dp), DIMENSION(:,:,:,:), POINTER  :: gdat  ! (global) data field
  REAL(dp), DIMENSION(:,:), POINTER      :: gps   ! surface pressure
  REAL(dp), DIMENSION(:,:,:,:), POINTER  :: zlptr => NULL() ! (local) data pointer
  !
  LOGICAL                            :: lsc       ! scatter results ?
  LOGICAL                            :: llrg, llrgx, llrgy, llrgz, llok
  INTEGER                            :: klat
  INTEGER                            :: klon
  INTEGER                            :: klev
  INTEGER                            :: kparam
  INTEGER                            :: s1d        ! size of 1D field

  ! INITIALIZE
  ! I/O
  IF (ASSOCIATED(dat)) THEN
!!$     DEALLOCATE(dat, STAT=status)
!!$     CALL ERRMSG(substr, status, 1)
     DEALLOCATE(dat)
  END IF
  NULLIFY(dat)
  !
  IF (PRESENT(ps)) THEN
     IF (ASSOCIATED(ps)) THEN
!!$        DEALLOCATE(ps, STAT=status)
!!$        CALL ERRMSG(substr, status, 2)
        DEALLOCATE(ps)
     END IF
     NULLIFY(ps)
  END IF
  !
  IF (PRESENT(hyam)) THEN
     IF (ASSOCIATED(hyam)) THEN
!!$        DEALLOCATE(hyam, STAT=status)
!!$        CALL ERRMSG(substr, status, 3)
        DEALLOCATE(hyam)
     END IF
     NULLIFY(hyam)
  END IF
  !
  IF (PRESENT(hybm)) THEN
     IF (ASSOCIATED(hybm)) THEN
!!$        DEALLOCATE(hybm, STAT=status)
!!$        CALL ERRMSG(substr, status, 4)
        DEALLOCATE(hybm)
     END IF
     NULLIFY(hybm)
  END IF
  !
  IF (PRESENT(hyai)) THEN
     IF (ASSOCIATED(hyai)) THEN
!!$        DEALLOCATE(hyai, STAT=status)
!!$        CALL ERRMSG(substr, status, 5)
        DEALLOCATE(hyai)
     END IF
     NULLIFY(hyai)
  END IF
  !
  IF (PRESENT(hybi)) THEN
     IF (ASSOCIATED(hybi)) THEN
!!$        DEALLOCATE(hybi, STAT=status)
!!$        CALL ERRMSG(substr, status, 6)
        DEALLOCATE(hybi)
     END IF
     NULLIFY(hybi)
  END IF
  !
  IF (PRESENT(p0)) THEN
     IF (ASSOCIATED(p0)) THEN
!!$        DEALLOCATE(p0, STAT=status)
!!$        CALL ERRMSG(substr, status, 7)
        DEALLOCATE(p0)
     END IF
     NULLIFY(p0)
  END IF
  !
  IF (PRESENT(latm)) THEN
     IF (ASSOCIATED(latm)) THEN
!!$        DEALLOCATE(latm, STAT=status)
!!$        CALL ERRMSG(substr, status, 8)
        DEALLOCATE(latm)
     END IF
     NULLIFY(latm)
  END IF
  !
  IF (PRESENT(lonm)) THEN
     IF (ASSOCIATED(lonm)) THEN
!!$        DEALLOCATE(lonm, STAT=status)
!!$        CALL ERRMSG(substr, status, 9)
        DEALLOCATE(lonm)
     END IF
     NULLIFY(lonm)
  END IF
  !
  IF (PRESENT(lati)) THEN
     IF (ASSOCIATED(lati)) THEN
!!$        DEALLOCATE(lati, STAT=status)
!!$        CALL ERRMSG(substr, status, 10)
        DEALLOCATE(lati)
     END IF
     NULLIFY(lati)
  END IF
  !
  IF (PRESENT(loni)) THEN
     IF (ASSOCIATED(loni)) THEN
!!$        DEALLOCATE(loni, STAT=status)
!!$        CALL ERRMSG(substr, status, 11)
        DEALLOCATE(loni)
     END IF
     NULLIFY(loni)
  END IF
  !
  ! LOCAL
  IF (PRESENT(lrg)) THEN
     llrg = lrg
  ELSE
     llrg = .true.
  END IF
  IF (PRESENT(lrgx)) THEN
     llrgx = lrgx
  ELSE
     llrgx = .true.
  END IF
  IF (PRESENT(lrgy)) THEN
     llrgy = lrgy
  ELSE
     llrgy = .true.
  END IF
  IF (PRESENT(lrgz)) THEN
     llrgz = lrgz
  ELSE
     llrgz = .true.
  END IF
  !
  NULLIFY(gdat)       ! (global) data field
  NULLIFY(gps)        ! (global) surface pressure field

#ifdef HOMMESE
  llrg = .false.
  llrgx = .false.
  llrgy = .false.
  llrgz = .false.
#endif

  IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL RGTOOL_READ_NCVAR(iou, TRIM(modstr)//'.nml', vname, t, var &
                           ,grid, llrg, llrgx, llrgy, llrgz, llok)
     IF (llok) THEN
        CALL RGTOOL_CONVERT(var, gdat, grid, order='xzny')
        klon = SIZE(gdat, 1)
        klev = SIZE(gdat, 2)
        kparam = SIZE(gdat, 3)
        klat = SIZE(gdat, 4)
        CALL RGTOOL_G2C(grid, hyam, hybm, p0, gps      &
                        ,hyai, hybi                    &
                        ,latm, lonm                    &
                        ,lati, loni, klat, klon, klev)
     END IF
     CALL INIT_NCVAR(var)
     CALL INIT_GEOHYBGRID(grid)
  END IF ! p_parallel_io

  ! BROADCAST SUCCESS
  CALL p_bcast(llok, p_io)
  IF (PRESENT(lok)) THEN
     lok = llok
  END IF

  ! BROADCST DIMENSIONS
  IF (llok) THEN
     CALL p_bcast(klon, p_io)
     CALL p_bcast(klev, p_io)
     CALL p_bcast(klat, p_io)
     CALL p_bcast(kparam, p_io)
  END IF

  ! NOTE: SCATTERING THE RESULTS (dat) TO THE PROCESSORS IS ONLY POSSIBLE
  !       IF REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy.AND.lrgz)
#ifndef HOMMESE
  lsc = llrg.AND.llok.AND.llrgx.AND.llrgy.AND.llrgz
#else
  lsc = llok
#endif
  !
  IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
     ! ALLOCATE SPACE FOR LOCAL DATA FIELD
!!$     ALLOCATE(dat(nproma,klev,kparam,ngpblks), STAT=status)
!!$     CALL ERRMSG(substr, status, 12)
     ALLOCATE(dat(nproma,klev,kparam,ngpblks,1))
     ! SCATTER GLOBAL DATA TO LOCAL FIELDS
     zlptr => dat(:,:,:,:,1)
     CALL scatter_gp(gdat, zlptr, dcg)
  ELSE            ! COPY AND BROADCAST ENTIRE FIELD
     IF (llok) THEN
        ! ALLOCATE SPACE FOR DATA FIELD
!!$        ALLOCATE(dat(klon,klev,kparam,klat),STAT=status)
!!$        CALL ERRMSG(substr, status, 13)
        ALLOCATE(dat(klon,klev,kparam,klat,1))
        ! COPY FIRST ON I/O-PROCESSOR
        IF (p_parallel_io) THEN
           dat(:,:,:,:,1) = gdat(:,:,:,:)
        END IF
        ! BROADCAST RESULT
        CALL p_bcast(dat(:,:,:,:,1), p_io)
     END IF ! llok
  END IF ! lsc

  ! NOTE: SCATTERING THE RESULTS (ps) TO THE PROCESSORS IS ONLY POSSIBLE
  !       IF HORIZONTAL REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy)
  IF (PRESENT(ps)) THEN
#ifndef HOMMESE
     lsc = llrg.AND.llok.AND.llrgx.AND.llrgy
#else
     lsc = llok
#endif
     !
     IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
        ! ALLOCATE SPACE FOR LOCAL DATA FIELD
!!$        ALLOCATE(ps(nproma,ngpblks), STAT=status)
!!$        CALL ERRMSG(substr, status, 14)
        ALLOCATE(ps(nproma,ngpblks))
        ! SCATTER GLOBAL DATA TO LOCAL FIELDS
        CALL scatter_gp(gps, ps, dcg)
     ELSE
        IF (llok) THEN
           ! ALLOCATE SPACE FOR DATA FIELD
!!$           ALLOCATE(ps(klon,klat),STAT=status)
!!$           CALL ERRMSG(substr, status, 15)
           ALLOCATE(ps(klon,klat))
           IF (p_parallel_io) THEN
              ps(:,:) = gps(:,:)
           END IF
           ! BROADCAST RESULT
           CALL p_bcast(ps , p_io)
        END IF ! llok
     END IF ! lsc
  END IF ! ps present

  ! 1D ARRAYS NEVER SCATTERED
  IF (llok) THEN
     IF (PRESENT(hyam)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hyam)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hyam(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 16)
           ALLOCATE(hyam(s1d))
        END IF
        CALL p_bcast(hyam ,p_io)
     END IF
     !
     IF (PRESENT(hybm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hybm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hybm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 17)
           ALLOCATE(hybm(s1d))
        END IF
        CALL p_bcast(hybm ,p_io)
     END IF
     !
     IF (PRESENT(hyai)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hyai)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hyai(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 18)
           ALLOCATE(hyai(s1d))
        END IF
        CALL p_bcast(hyai ,p_io)
     END IF
     !
     IF (PRESENT(hybi)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hybi)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hybi(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 19)
           ALLOCATE(hybi(s1d))
        END IF
        CALL p_bcast(hybi ,p_io)
     END IF
     !
     IF (PRESENT(p0)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(p0)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(p0(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 20)
           ALLOCATE(p0(s1d))
        END IF
        CALL p_bcast(p0 ,p_io)
     END IF
     !
     IF (PRESENT(latm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(latm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(latm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 21)
           ALLOCATE(latm(s1d))
        END IF
        CALL p_bcast(latm ,p_io)
     END IF
     !
     IF (PRESENT(lonm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(lonm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(lonm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 22)
           ALLOCATE(lonm(s1d))
        END IF
        CALL p_bcast(lonm ,p_io)
     END IF
     !
     IF (PRESENT(lati)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(lati)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(lati(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 23)
           ALLOCATE(lati(s1d))
        END IF
        CALL p_bcast(lati ,p_io)
     END IF
     !
     IF (PRESENT(loni)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(loni)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(loni(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 24)
           ALLOCATE(loni(s1d))
        END IF
        CALL p_bcast(loni ,p_io)
     END IF
  END IF

  ! CLEAN UP: FREE MEMORY OF GLOBAL FIELDS
  IF (llok) THEN
     IF (ASSOCIATED(gdat)) THEN
!!$        DEALLOCATE(gdat, STAT=status)
!!$        CALL ERRMSG(substr, status, 25)
        DEALLOCATE(gdat) ; NULLIFY(gdat)
     END IF
!!$     IF (PRESENT(ps)) THEN ! op_pj_20130806 memory leak
        IF (ASSOCIATED(gps)) THEN
!!$           DEALLOCATE(gps, STAT=status)
!!$           CALL ERRMSG(substr, status, 26)
           DEALLOCATE(gps) ; NULLIFY(gps)
        END IF
!!$     END IF ! op_pj_20130806 memory leak
  END IF

END SUBROUTINE RGTOOL_BI_READ_NCVAR_5D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTOOL_BI_READ_NCFILE_4D(modstr, fname, t, dat, vars   &
                               , lrg, lrgx, lrgy, lrgz, lok    &  ! OPTIONAL
                               , hyam, hybm, p0, ps            &  ! OPTIONAL
                               , hyai, hybi                    &  ! OPTIONAL
                               , latm, lonm, lati, loni        &  ! OPTIONAL
                               )

  ! PERFORMS ONE NCREGRID STEP FOR ONE FILE
  ! AND RETURNS m DATA FIELDS AS 4D ARRAY (x,z,m,y),
  ! AND OPTIONALLY GRID-STRUCTURE AS ARRAYS
  ! NOTE: THIS DOES ONLY WORK FOR (x,z,n=1,y) VARIABLES
  !       (i.e., 3-D VAIABLES)
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast &
                                , dcg, scatter_gp
  USE messy_main_data_bi,   ONLY: ngpblks, nproma
  ! MESSY
  USE messy_main_tools,     ONLY: find_next_free_unit
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf, ONLY: ncvar, init_ncvar,GRD_MAXSTRLEN
  USE messy_ncregrid_geohyb, ONLY: geohybgrid, init_geohybgrid
  USE messy_ncregrid_tools,  ONLY: rgtool_read_ncfile, rgtool_convert &
                                 , rgtool_g2c


  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, PRESENT, SIZE, TRIM

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)             :: modstr    ! calling module
  CHARACTER(LEN=*), INTENT(IN)             :: fname     ! file name
  INTEGER,          INTENT(IN)             :: t         ! netCDF time step
  REAL(dp), DIMENSION(:,:,:,:), POINTER    :: dat       ! (local) data field
  CHARACTER(LEN=GRD_MAXSTRLEN), DIMENSION(:), &
                                   POINTER :: vars      ! variable names
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrg       ! regrid really ?
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgx      ! regrid in x
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgy      ! regrid in y
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgz      ! regrid in z
  LOGICAL,          INTENT(OUT),  OPTIONAL :: lok       ! OK?
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyam  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybm  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: p0    ! reference pressure
  REAL(dp), DIMENSION(:,:), POINTER,  OPTIONAL :: ps    ! (local) surf. press.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyai  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybi  ! hybrid-B-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: latm  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lonm  ! longitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lati  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: loni  ! longitude

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_BI_READ_NCFILE_4D'
  INTEGER                             :: iou        ! logical I/O unit
!!$  INTEGER                             :: status     ! memory status
  TYPE(ncvar), DIMENSION(:), POINTER  :: var        ! nc-variable list
  TYPE(geohybgrid)                    :: grid       ! grid-structure
  REAL(dp), DIMENSION(:,:,:,:), POINTER   :: hdat   ! temporary data field
  REAL(dp), DIMENSION(:,:,:,:), POINTER   :: gdat   ! (global) data field
  REAL(dp), DIMENSION(:,:), POINTER       :: gps    ! (global) surface pressure
  !
  LOGICAL                             :: lsc        ! scatter results ?
  LOGICAL                             :: llrg, llrgx, llrgy, llrgz, llok
  INTEGER                             :: i          ! counter
  INTEGER                             :: nvar       ! number of variables
  INTEGER                             :: klat       ! number of latitudes
  INTEGER                             :: klon       ! number of longitudes
  INTEGER                             :: klev       ! number of levels
  INTEGER                             :: kparam     ! number of parameters
  INTEGER                             :: s1d        ! size of 1D field

  ! INITIALIZE
  ! I/O
  IF (ASSOCIATED(dat)) THEN
!!$     DEALLOCATE(dat, STAT=status)
!!$     CALL ERRMSG(substr, status, 1)
     DEALLOCATE(dat)
  END IF
  NULLIFY(dat)
  !
  IF (ASSOCIATED(vars)) THEN
!!$     DEALLOCATE(vars, STAT=status)
!!$     CALL ERRMSG(substr, status, 2)
     DEALLOCATE(vars)
  END IF
  NULLIFY(vars)
  !
  IF (PRESENT(ps)) THEN
     IF (ASSOCIATED(ps)) THEN
!!$        DEALLOCATE(ps, STAT=status)
!!$        CALL ERRMSG(substr, status, 3)
        DEALLOCATE(ps)
     END IF
     NULLIFY(ps)
  END IF
  !
  IF (PRESENT(hyam)) THEN
     IF (ASSOCIATED(hyam)) THEN
!!$        DEALLOCATE(hyam, STAT=status)
!!$        CALL ERRMSG(substr, status, 4)
        DEALLOCATE(hyam)
     END IF
     NULLIFY(hyam)
  END IF
  !
  IF (PRESENT(hybm)) THEN
     IF (ASSOCIATED(hybm)) THEN
!!$        DEALLOCATE(hybm, STAT=status)
!!$        CALL ERRMSG(substr, status, 5)
        DEALLOCATE(hybm)
     END IF
     NULLIFY(hybm)
  END IF
  !
  IF (PRESENT(hyai)) THEN
     IF (ASSOCIATED(hyai)) THEN
!!$        DEALLOCATE(hyai, STAT=status)
!!$        CALL ERRMSG(substr, status, 6)
        DEALLOCATE(hyai)
     END IF
     NULLIFY(hyai)
  END IF
  !
  IF (PRESENT(hybi)) THEN
     IF (ASSOCIATED(hybi)) THEN
!!$        DEALLOCATE(hybi, STAT=status)
!!$        CALL ERRMSG(substr, status, 7)
        DEALLOCATE(hybi)
     END IF
     NULLIFY(hybi)
  END IF
  !
  IF (PRESENT(p0)) THEN
     IF (ASSOCIATED(p0)) THEN
!!$        DEALLOCATE(p0, STAT=status)
!!$        CALL ERRMSG(substr, status, 8)
        DEALLOCATE(p0)
     END IF
     NULLIFY(p0)
  END IF
  !
  IF (PRESENT(latm)) THEN
     IF (ASSOCIATED(latm)) THEN
!!$        DEALLOCATE(latm, STAT=status)
!!$        CALL ERRMSG(substr, status, 9)
        DEALLOCATE(latm)
     END IF
     NULLIFY(latm)
  END IF
  !
  IF (PRESENT(lonm)) THEN
     IF (ASSOCIATED(lonm)) THEN
!!$        DEALLOCATE(lonm, STAT=status)
!!$        CALL ERRMSG(substr, status, 10)
        DEALLOCATE(lonm)
     END IF
     NULLIFY(lonm)
  END IF
  !
  IF (PRESENT(lati)) THEN
     IF (ASSOCIATED(lati)) THEN
!!$        DEALLOCATE(lati, STAT=status)
!!$        CALL ERRMSG(substr, status, 11)
        DEALLOCATE(lati)
     END IF
     NULLIFY(lati)
  END IF
  !
  IF (PRESENT(loni)) THEN
     IF (ASSOCIATED(loni)) THEN
!!$        DEALLOCATE(loni, STAT=status)
!!$        CALL ERRMSG(substr, status, 12)
        DEALLOCATE(loni)
     END IF
     NULLIFY(loni)
  END IF
  !
  ! LOCAL
  IF (PRESENT(lrg)) THEN
     llrg = lrg
  ELSE
     llrg = .true.
  END IF
  IF (PRESENT(lrgx)) THEN
     llrgx = lrgx
  ELSE
     llrgx = .true.
  END IF
  IF (PRESENT(lrgy)) THEN
     llrgy = lrgy
  ELSE
     llrgy = .true.
  END IF
  IF (PRESENT(lrgz)) THEN
     llrgz = lrgz
  ELSE
     llrgz = .true.
  END IF
  !
  NULLIFY(hdat)       ! temporary data field
  NULLIFY(gdat)       ! (global) data field
  NULLIFY(gps)        ! (global) surface pressure field
  NULLIFY(var)        ! variable names

#ifdef HOMMESE
  llrg = .false.
  llrgx = .false.
  llrgy = .false.
  llrgz = .false.
#endif

  IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL RGTOOL_READ_NCFILE(iou, TRIM(modstr)//'.nml', fname, t, var &
                            ,grid, llrg, llrgx, llrgy, llrgz, llok)
     IF (llok) THEN
        nvar = SIZE(var)
!!$        ALLOCATE(vars(nvar),STAT=status)
!!$        CALL ERRMSG(substr, status, 13)
        ALLOCATE(vars(nvar))
        ! 1st STEP
        vars(1) = TRIM(var(1)%name)
        CALL RGTOOL_CONVERT(var(1), hdat, grid, order='xzny')
        CALL INIT_NCVAR(var(1))
        klon = SIZE(hdat, 1)
        klev = SIZE(hdat, 2)
        klat = SIZE(hdat, 4)
        kparam = SIZE(hdat, 3)
        IF (kparam > 1) THEN
           CALL RGMSG(substr, RGMLE, &
                'FILE REGRIDDING NOT POSSIBLE FOR PARAMETER VARIABLES !')
        END IF
        ! ALLOCATE SPACE FOR GLOBAL DATA FIELD
!!$        ALLOCATE(gdat(klon, klev, nvar, klat), STAT=status)
!!$        CALL ERRMSG(substr, status, 14)
        ALLOCATE(gdat(klon, klev, nvar, klat))
        ! COPY VARIABLE-DATA TO GLOBAL DATA FIELD
        gdat(:,:,1,:) = hdat(:,:,1,:)
        ! FREE MEMORY OF TEMPORARY-DATA
!!$        DEALLOCATE(hdat, STAT=status)
!!$        CALL ERRMSG(substr, status, 15)
        DEALLOCATE(hdat)
        NULLIFY(hdat)
        ! 2nd ... nvar
        DO i=2, nvar
           vars(i) = TRIM(var(i)%name)
           CALL RGTOOL_CONVERT(var(i), hdat, grid, order='xzny')
           CALL INIT_NCVAR(var(i))
           kparam = SIZE(hdat, 3)
           IF (kparam > 1) THEN
              CALL RGMSG(substr, RGMLE, &
                   'FILE REGRIDDING NOT POSSIBLE FOR PARAMETER VARIABLES !')
           END IF
           ! COPY VARIABLE-DATA TO GLOBAL DATA FIELD
           gdat(:,:,i,:) = hdat(:,:,1,:)
           ! FREE MEMORY OF VARIABLE-DATA
!!$           DEALLOCATE(hdat, STAT=status)
!!$           CALL ERRMSG(substr, status, 16)
           DEALLOCATE(hdat)
           NULLIFY(hdat)
        END DO
!!$        DEALLOCATE(var, STAT=status)
!!$        CALL ERRMSG(substr, status, 17)
        DEALLOCATE(var)
        NULLIFY(var)
        ! GRID
        CALL RGTOOL_G2C(grid, hyam, hybm, p0, gps     &
                       ,hyai, hybi                    &
                       ,latm, lonm                    &
                       ,lati, loni, klat, klon, klev)
        CALL INIT_GEOHYBGRID(grid)
     END IF ! llok
  END IF ! p_parallel_io
  !
  ! BROADCAST SUCCESS
  CALL p_bcast(llok, p_io)
  IF (PRESENT(lok)) THEN
     lok = llok
  END IF

  ! BROADCST DIMENSIONS
  IF (llok) THEN
     CALL p_bcast(nvar, p_io)
     CALL p_bcast(klon, p_io)
     CALL p_bcast(klev, p_io)
     CALL p_bcast(klat, p_io)
     CALL p_bcast(kparam, p_io)
  END IF

  ! BROADCAST VARIABLE NAMES
  IF (llok) THEN
     IF (.NOT.p_parallel_io) THEN
!!$        ALLOCATE(vars(nvar), STAT=status)
!!$        CALL ERRMSG(substr, status, 18)
        ALLOCATE(vars(nvar))
     END IF
     DO i=1, nvar
        CALL p_bcast(vars(i), p_io)
     END DO
  END IF

  ! NOTE: SCATTERING THE RESULTS (dat, ps) TO THE PROCESSORS IS ONLY POSSIBLE
  !       IF REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy.AND.lrgz)
#ifndef HOMMESE
  lsc = llrg.AND.llok.AND.llrgx.AND.llrgy.AND.llrgz
#else
  lsc = llok
#endif
  !
  IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
     ! ALLOCATE SPACE FOR LOCAL DATA FIELD
!!$    ALLOCATE(dat(nproma,klev,nvar,ngpblks), STAT=status)
!!$     CALL ERRMSG(substr, status, 19)
    ALLOCATE(dat(nproma,klev,nvar,ngpblks))
     ! SCATTER GLOBAL DATA TO LOCAL FIELDS
     CALL scatter_gp(gdat, dat, dcg)
  ELSE            ! COPY AND BROADCAST ENTIRE FIELD
     IF (llok) THEN
        ! ALLOCATE SPACE FOR DATA FIELD
!!$        ALLOCATE(dat(klon,klev,nvar,klat),STAT=status)
!!$        CALL ERRMSG(substr, status, 20)
        ALLOCATE(dat(klon,klev,nvar,klat))
        ! COPY FIRST ON I/O-PROCESSOR
        IF (p_parallel_io) THEN
           dat(:,:,:,:) = gdat(:,:,:,:)
        END IF
        ! BROADCAST RESULT
        CALL p_bcast(dat, p_io)
     END IF ! llok
  END IF ! lsc

  ! NOTE: SCATTERING THE RESULTS (ps) TO THE PROCESSORS IS ONLY POSSIBLE
  !       IF HORIZONTAL REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy)
  IF (PRESENT(ps)) THEN
     lsc = llrg.AND.llok.AND.llrgx.AND.llrgy
     !
     IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
        ! ALLOCATE SPACE FOR LOCAL DATA FIELD
!!$        ALLOCATE(ps(nproma,ngpblks), STAT=status)
!!$        CALL ERRMSG(substr, status, 21)
        ALLOCATE(ps(nproma,ngpblks))
        ! SCATTER GLOBAL DATA TO LOCAL FIELDS
        CALL scatter_gp(gps, ps, dcg)
     ELSE
        IF (llok) THEN
           ! ALLOCATE SPACE FOR DATA FIELD
!!$           ALLOCATE(ps(klon,klat),STAT=status)
!!$           CALL ERRMSG(substr, status, 22)
           ALLOCATE(ps(klon,klat))
           IF (p_parallel_io) THEN
              ps(:,:) = gps(:,:)
           END IF
           ! BROADCAST RESULT
           CALL p_bcast(ps , p_io)
        END IF ! llok
     END IF ! lsc
  END IF ! ps present

  ! 1D ARRAYS NEVER SCATTERED
  IF (llok) THEN
     IF (PRESENT(hyam)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hyam)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hyam(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 23)
           ALLOCATE(hyam(s1d))
        END IF
        CALL p_bcast(hyam ,p_io)
     END IF
     !
     IF (PRESENT(hybm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hybm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hybm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 24)
           ALLOCATE(hybm(s1d))
        END IF
        CALL p_bcast(hybm ,p_io)
     END IF
     !
     IF (PRESENT(hyai)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hyai)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hyai(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 25)
           ALLOCATE(hyai(s1d))
        END IF
        CALL p_bcast(hyai ,p_io)
     END IF
     !
     IF (PRESENT(hybi)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hybi)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hybi(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 26)
           ALLOCATE(hybi(s1d))
        END IF
        CALL p_bcast(hybi ,p_io)
     END IF
     !
     IF (PRESENT(p0)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(p0)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(p0(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 27)
           ALLOCATE(p0(s1d))
        END IF
        CALL p_bcast(p0 ,p_io)
     END IF
     !
     IF (PRESENT(latm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(latm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(latm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 28)
           ALLOCATE(latm(s1d))
        END IF
        CALL p_bcast(latm ,p_io)
     END IF
     !
     IF (PRESENT(lonm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(lonm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(lonm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 29)
           ALLOCATE(lonm(s1d))
        END IF
        CALL p_bcast(lonm ,p_io)
     END IF
     !
     IF (PRESENT(lati)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(lati)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(lati(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 30)
           ALLOCATE(lati(s1d))
        END IF
        CALL p_bcast(lati ,p_io)
     END IF
     !
     IF (PRESENT(loni)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(loni)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(loni(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 31)
           ALLOCATE(loni(s1d))
        END IF
        CALL p_bcast(loni ,p_io)
     END IF
  END IF

  ! CLEAN UP: FREE MEMORY OF GLOBAL FIELDS ON I/O PROCESSOR
  IF (llok) THEN
     IF (ASSOCIATED(gdat)) THEN
!!$        DEALLOCATE(gdat, STAT=status)
!!$        CALL ERRMSG(substr, status, 32)
        DEALLOCATE(gdat) ; NULLIFY(gdat)
     END IF
!!$     IF (PRESENT(ps)) THEN ! op_pj_20130806 memory leak
        IF (ASSOCIATED(gps)) THEN
!!$           DEALLOCATE(gps, STAT=status)
!!$           CALL ERRMSG(substr, status, 33)
           DEALLOCATE(gps) ; NULLIFY(gps)
        END IF
!!$     END IF ! op_pj_20130806 memory leak
  END IF

END SUBROUTINE RGTOOL_BI_READ_NCFILE_4D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTOOL_BI_READ_NCFILE_5D(modstr, fname, t, dat, vars   &
                               , lrg, lrgx, lrgy, lrgz, lok    &  ! OPTIONAL
                               , hyam, hybm, p0, ps            &  ! OPTIONAL
                               , hyai, hybi                    &  ! OPTIONAL
                               , latm, lonm, lati, loni        &  ! OPTIONAL
                               )

  ! PERFORMS ONE NCREGRID STEP FOR ONE FILE
  ! AND RETURNS m DATA FIELDS AS 5D ARRAY (x,z,m,y,n),
  ! AND OPTIONALLY GRID-STRUCTURE AS ARRAYS
  !
  ! Author: Patrick Joeckel, DLR, May 2012

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast &
                                , dcg, scatter_gp
  USE messy_main_data_bi,   ONLY: ngpblks, nproma
  ! MESSY
  USE messy_main_tools,     ONLY: find_next_free_unit
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf, ONLY: ncvar, init_ncvar,GRD_MAXSTRLEN
  USE messy_ncregrid_geohyb, ONLY: geohybgrid, init_geohybgrid
  USE messy_ncregrid_tools,  ONLY: rgtool_read_ncfile, rgtool_convert &
                                 , rgtool_g2c


  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, PRESENT, SIZE, TRIM

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)             :: modstr    ! calling module
  CHARACTER(LEN=*), INTENT(IN)             :: fname     ! file name
  INTEGER,          INTENT(IN)             :: t         ! netCDF time step
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER  :: dat       ! (local) data field
  CHARACTER(LEN=GRD_MAXSTRLEN), DIMENSION(:), &
                                   POINTER :: vars      ! variable names
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrg       ! regrid really ?
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgx      ! regrid in x
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgy      ! regrid in y
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgz      ! regrid in z
  LOGICAL,          INTENT(OUT),  OPTIONAL :: lok       ! OK?
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyam  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybm  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: p0    ! reference pressure
  REAL(dp), DIMENSION(:,:), POINTER,  OPTIONAL :: ps    ! (local) surf. press.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyai  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybi  ! hybrid-B-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: latm  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lonm  ! longitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lati  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: loni  ! longitude

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_BI_READ_NCFILE_5D'
  INTEGER                             :: iou        ! logical I/O unit
!!$  INTEGER                             :: status     ! memory status
  TYPE(ncvar), DIMENSION(:), POINTER  :: var        ! nc-variable list
  TYPE(geohybgrid)                    :: grid       ! grid-structure
  REAL(dp), DIMENSION(:,:,:,:), POINTER   :: hdat   ! temporary data field
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: gdat   ! (global) data field
  REAL(dp), DIMENSION(:,:), POINTER       :: gps    ! (global) surface pressure
  REAL(dp), DIMENSION(:,:,:,:), POINTER   :: zgptr => NULL() ! (global) data pointer
  REAL(dp), DIMENSION(:,:,:,:), POINTER   :: zlptr => NULL() ! (local) data pointer
  !
  LOGICAL                             :: lsc        ! scatter results ?
  LOGICAL                             :: llrg, llrgx, llrgy, llrgz, llok
  INTEGER                             :: i          ! counter
  INTEGER                             :: nvar       ! number of variables
  INTEGER                             :: klat       ! number of latitudes
  INTEGER                             :: klon       ! number of longitudes
  INTEGER                             :: klev       ! number of levels
  INTEGER                             :: kparam     ! number of parameters
  INTEGER                             :: kparam2    ! number of parameters
  INTEGER                             :: s1d        ! size of 1D field

  ! INITIALIZE
  ! I/O
  IF (ASSOCIATED(dat)) THEN
!!$     DEALLOCATE(dat, STAT=status)
!!$     CALL ERRMSG(substr, status, 1)
     DEALLOCATE(dat)
  END IF
  NULLIFY(dat)
  !
  IF (ASSOCIATED(vars)) THEN
!!$     DEALLOCATE(vars, STAT=status)
!!$     CALL ERRMSG(substr, status, 2)
     DEALLOCATE(vars)
  END IF
  NULLIFY(vars)
  !
  IF (PRESENT(ps)) THEN
     IF (ASSOCIATED(ps)) THEN
!!$        DEALLOCATE(ps, STAT=status)
!!$        CALL ERRMSG(substr, status, 3)
        DEALLOCATE(ps)
     END IF
     NULLIFY(ps)
  END IF
  !
  IF (PRESENT(hyam)) THEN
     IF (ASSOCIATED(hyam)) THEN
!!$        DEALLOCATE(hyam, STAT=status)
!!$        CALL ERRMSG(substr, status, 4)
        DEALLOCATE(hyam)
     END IF
     NULLIFY(hyam)
  END IF
  !
  IF (PRESENT(hybm)) THEN
     IF (ASSOCIATED(hybm)) THEN
!!$        DEALLOCATE(hybm, STAT=status)
!!$        CALL ERRMSG(substr, status, 5)
        DEALLOCATE(hybm)
     END IF
     NULLIFY(hybm)
  END IF
  !
  IF (PRESENT(hyai)) THEN
     IF (ASSOCIATED(hyai)) THEN
!!$        DEALLOCATE(hyai, STAT=status)
!!$        CALL ERRMSG(substr, status, 6)
        DEALLOCATE(hyai)
     END IF
     NULLIFY(hyai)
  END IF
  !
  IF (PRESENT(hybi)) THEN
     IF (ASSOCIATED(hybi)) THEN
!!$        DEALLOCATE(hybi, STAT=status)
!!$        CALL ERRMSG(substr, status, 7)
        DEALLOCATE(hybi)
     END IF
     NULLIFY(hybi)
  END IF
  !
  IF (PRESENT(p0)) THEN
     IF (ASSOCIATED(p0)) THEN
!!$        DEALLOCATE(p0, STAT=status)
!!$        CALL ERRMSG(substr, status, 8)
        DEALLOCATE(p0)
     END IF
     NULLIFY(p0)
  END IF
  !
  IF (PRESENT(latm)) THEN
     IF (ASSOCIATED(latm)) THEN
!!$        DEALLOCATE(latm, STAT=status)
!!$        CALL ERRMSG(substr, status, 9)
        DEALLOCATE(latm)
     END IF
     NULLIFY(latm)
  END IF
  !
  IF (PRESENT(lonm)) THEN
     IF (ASSOCIATED(lonm)) THEN
!!$        DEALLOCATE(lonm, STAT=status)
!!$        CALL ERRMSG(substr, status, 10)
        DEALLOCATE(lonm)
     END IF
     NULLIFY(lonm)
  END IF
  !
  IF (PRESENT(lati)) THEN
     IF (ASSOCIATED(lati)) THEN
!!$        DEALLOCATE(lati, STAT=status)
!!$        CALL ERRMSG(substr, status, 11)
        DEALLOCATE(lati)
     END IF
     NULLIFY(lati)
  END IF
  !
  IF (PRESENT(loni)) THEN
     IF (ASSOCIATED(loni)) THEN
!!$        DEALLOCATE(loni, STAT=status)
!!$        CALL ERRMSG(substr, status, 12)
        DEALLOCATE(loni)
     END IF
     NULLIFY(loni)
  END IF
  !
  ! LOCAL
  IF (PRESENT(lrg)) THEN
     llrg = lrg
  ELSE
     llrg = .true.
  END IF
  IF (PRESENT(lrgx)) THEN
     llrgx = lrgx
  ELSE
     llrgx = .true.
  END IF
  IF (PRESENT(lrgy)) THEN
     llrgy = lrgy
  ELSE
     llrgy = .true.
  END IF
  IF (PRESENT(lrgz)) THEN
     llrgz = lrgz
  ELSE
     llrgz = .true.
  END IF
  !
  NULLIFY(hdat)       ! temporary data field
  NULLIFY(gdat)       ! (global) data field
  NULLIFY(gps)        ! (global) surface pressure field
  NULLIFY(var)        ! variable names

#ifdef HOMMESE
  llrg = .false.
  llrgx = .false.
  llrgy = .false.
  llrgz = .false.
#endif

  IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL RGTOOL_READ_NCFILE(iou, TRIM(modstr)//'.nml', fname, t, var &
                            ,grid, llrg, llrgx, llrgy, llrgz, llok)
     IF (llok) THEN
        nvar = SIZE(var)
!!$        ALLOCATE(vars(nvar),STAT=status)
!!$        CALL ERRMSG(substr, status, 13)
        ALLOCATE(vars(nvar))
        ! 1st STEP
        vars(1) = TRIM(var(1)%name)
        CALL RGTOOL_CONVERT(var(1), hdat, grid, order='xzny')
        CALL INIT_NCVAR(var(1))
        klon = SIZE(hdat, 1)
        klev = SIZE(hdat, 2)
        klat = SIZE(hdat, 4)
        kparam = SIZE(hdat, 3)
        ! ALLOCATE SPACE FOR GLOBAL DATA FIELD
!!$        ALLOCATE(gdat(klon, klev, kparam, klat, nvar), STAT=status)
!!$        CALL ERRMSG(substr, status, 14)
        ALLOCATE(gdat(klon, klev, kparam, klat, nvar))
        ! COPY VARIABLE-DATA TO GLOBAL DATA FIELD
        gdat(:,:,:,:,1) = hdat(:,:,:,:)
        ! FREE MEMORY OF TEMPORARY-DATA
!!$        DEALLOCATE(hdat, STAT=status)
!!$        CALL ERRMSG(substr, status, 15)
        DEALLOCATE(hdat)
        NULLIFY(hdat)
        ! 2nd ... nvar
        DO i=2, nvar
           vars(i) = TRIM(var(i)%name)
           CALL RGTOOL_CONVERT(var(i), hdat, grid, order='xzny')
           CALL INIT_NCVAR(var(i))
           kparam2 = SIZE(hdat, 3)
           IF (kparam2 /= kparam) THEN
              CALL RGMSG(substr, RGMLE, &
                   'PARAMETER DIMENSION MISMATCH OF VARIABLES !')
           END IF
           ! COPY VARIABLE-DATA TO GLOBAL DATA FIELD
           gdat(:,:,:,:,i) = hdat(:,:,:,:)
           ! FREE MEMORY OF VARIABLE-DATA
!!$           DEALLOCATE(hdat, STAT=status)
!!$           CALL ERRMSG(substr, status, 16)
           DEALLOCATE(hdat)
           NULLIFY(hdat)
        END DO
!!$        DEALLOCATE(var, STAT=status)
!!$        CALL ERRMSG(substr, status, 17)
        DEALLOCATE(var)
        NULLIFY(var)
        ! GRID
        CALL RGTOOL_G2C(grid, hyam, hybm, p0, gps     &
                       ,hyai, hybi                    &
                       ,latm, lonm                    &
                       ,lati, loni, klat, klon, klev)
        CALL INIT_GEOHYBGRID(grid)
     END IF ! llok
  END IF ! p_parallel_io
  !
  ! BROADCAST SUCCESS
  CALL p_bcast(llok, p_io)
  IF (PRESENT(lok)) THEN
     lok = llok
  END IF

  ! BROADCST DIMENSIONS
  IF (llok) THEN
     CALL p_bcast(nvar, p_io)
     CALL p_bcast(klon, p_io)
     CALL p_bcast(klev, p_io)
     CALL p_bcast(klat, p_io)
     CALL p_bcast(kparam, p_io)
  END IF

  ! BROADCAST VARIABLE NAMES
  IF (llok) THEN
     IF (.NOT.p_parallel_io) THEN
!!$        ALLOCATE(vars(nvar), STAT=status)
!!$        CALL ERRMSG(substr, status, 18)
        ALLOCATE(vars(nvar))
     END IF
     DO i=1, nvar
        CALL p_bcast(vars(i), p_io)
     END DO
  END IF

  ! NOTE: SCATTERING THE RESULTS (dat, ps) TO THE PROCESSORS IS ONLY POSSIBLE
  !       IF REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy.AND.lrgz)
#ifndef HOMMESE
  lsc = llrg.AND.llok.AND.llrgx.AND.llrgy.AND.llrgz
#else
  lsc = llok
#endif
  !
  IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
     ! ALLOCATE SPACE FOR LOCAL DATA FIELD
!!$    ALLOCATE(dat(nproma,klev,kparam,ngpblks,nvar), STAT=status)
!!$     CALL ERRMSG(substr, status, 19)
    ALLOCATE(dat(nproma,klev,kparam,ngpblks,nvar))
     ! SCATTER GLOBAL DATA TO LOCAL FIELDS
    DO i=1, nvar
       IF (p_parallel_io) zgptr => gdat(:,:,:,:,i)
       zlptr => dat(:,:,:,:,i)
       CALL scatter_gp(zgptr, zlptr, dcg)
    END DO
  ELSE            ! COPY AND BROADCAST ENTIRE FIELD
     IF (llok) THEN
        ! ALLOCATE SPACE FOR DATA FIELD
!!$        ALLOCATE(dat(klon,klev,kparam,klat,nvar),STAT=status)
!!$        CALL ERRMSG(substr, status, 20)
        ALLOCATE(dat(klon,klev,kparam,klat,nvar))
        ! COPY FIRST ON I/O-PROCESSOR
        IF (p_parallel_io) THEN
           dat(:,:,:,:,:) = gdat(:,:,:,:,:)
        END IF
        ! BROADCAST RESULT
        DO i=1, nvar
           CALL p_bcast(dat(:,:,:,:,i), p_io)
        END DO
     END IF ! llok
  END IF ! lsc

  ! NOTE: SCATTERING THE RESULTS (ps) TO THE PROCESSORS IS ONLY POSSIBLE
  !       IF HORIZONTAL REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy)
  IF (PRESENT(ps)) THEN
     lsc = llrg.AND.llok.AND.llrgx.AND.llrgy
     !
     IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
        ! ALLOCATE SPACE FOR LOCAL DATA FIELD
!!$        ALLOCATE(ps(nproma,ngpblks), STAT=status)
!!$        CALL ERRMSG(substr, status, 21)
        ALLOCATE(ps(nproma,ngpblks))
        ! SCATTER GLOBAL DATA TO LOCAL FIELDS
        CALL scatter_gp(gps, ps, dcg)
     ELSE
        IF (llok) THEN
           ! ALLOCATE SPACE FOR DATA FIELD
!!$           ALLOCATE(ps(klon,klat),STAT=status)
!!$           CALL ERRMSG(substr, status, 22)
           ALLOCATE(ps(klon,klat))
           IF (p_parallel_io) THEN
              ps(:,:) = gps(:,:)
           END IF
           ! BROADCAST RESULT
           CALL p_bcast(ps , p_io)
        END IF ! llok
     END IF ! lsc
  END IF ! ps present

  ! 1D ARRAYS NEVER SCATTERED
  IF (llok) THEN
     IF (PRESENT(hyam)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hyam)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hyam(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 23)
           ALLOCATE(hyam(s1d))
        END IF
        CALL p_bcast(hyam ,p_io)
     END IF
     !
     IF (PRESENT(hybm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hybm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hybm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 24)
           ALLOCATE(hybm(s1d))
        END IF
        CALL p_bcast(hybm ,p_io)
     END IF
     !
     IF (PRESENT(hyai)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hyai)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hyai(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 25)
           ALLOCATE(hyai(s1d))
        END IF
        CALL p_bcast(hyai ,p_io)
     END IF
     !
     IF (PRESENT(hybi)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(hybi)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(hybi(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 26)
           ALLOCATE(hybi(s1d))
        END IF
        CALL p_bcast(hybi ,p_io)
     END IF
     !
     IF (PRESENT(p0)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(p0)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(p0(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 27)
           ALLOCATE(p0(s1d))
        END IF
        CALL p_bcast(p0 ,p_io)
     END IF
     !
     IF (PRESENT(latm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(latm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(latm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 28)
           ALLOCATE(latm(s1d))
        END IF
        CALL p_bcast(latm ,p_io)
     END IF
     !
     IF (PRESENT(lonm)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(lonm)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(lonm(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 29)
           ALLOCATE(lonm(s1d))
        END IF
        CALL p_bcast(lonm ,p_io)
     END IF
     !
     IF (PRESENT(lati)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(lati)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(lati(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 30)
           ALLOCATE(lati(s1d))
        END IF
        CALL p_bcast(lati ,p_io)
     END IF
     !
     IF (PRESENT(loni)) THEN
        IF (p_parallel_io) THEN
           s1d = SIZE(loni)
        END IF
        CALL p_bcast(s1d ,p_io)
        IF (.NOT.p_parallel_io) THEN
!!$           ALLOCATE(loni(s1d), STAT=status)
!!$           CALL ERRMSG(substr, status, 31)
           ALLOCATE(loni(s1d))
        END IF
        CALL p_bcast(loni ,p_io)
     END IF
  END IF

  ! CLEAN UP: FREE MEMORY OF GLOBAL FIELDS ON I/O PROCESSOR
  IF (llok) THEN
     IF (ASSOCIATED(gdat)) THEN
!!$        DEALLOCATE(gdat, STAT=status)
!!$        CALL ERRMSG(substr, status, 32)
        DEALLOCATE(gdat) ; NULLIFY(gdat)
     END IF
!!$     IF (PRESENT(ps)) THEN ! op_pj_20130806 memory leak
        IF (ASSOCIATED(gps)) THEN
!!$           DEALLOCATE(gps, STAT=status)
!!$           CALL ERRMSG(substr, status, 33)
           DEALLOCATE(gps) ; NULLIFY(gps)
        END IF
!!$     END IF ! op_pj_20130806 memory leak
  END IF

END SUBROUTINE RGTOOL_BI_READ_NCFILE_5D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_READ_2D(rgt, modstr, name, RGREAD_TYPE   &
                           ,efield, lrg, lstop, vars, levent)

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,        ONLY: p_parallel_io !!$, finish
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf,    ONLY: GRD_MAXSTRLEN

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, NULL, PRESENT, SIZE, TRIM

  ! I/O
  TYPE (RGTEVENT), DIMENSION(:), POINTER  :: rgt
  CHARACTER(LEN=*), INTENT(IN)            :: modstr       ! submodel-string
  CHARACTER(LEN=*), INTENT(IN)            :: name         ! rgtevent name
  INTEGER,          INTENT(IN)            :: RGREAD_TYPE  ! NCVAR or NCFILE
  REAL(dp), DIMENSION(:,:), INTENT(OUT)   :: efield       ! 2D field
  CHARACTER(LEN=GRD_MAXSTRLEN),  OPTIONAL &
                   ,DIMENSION(:), POINTER :: vars    ! variable names
  LOGICAL,          INTENT(IN),  OPTIONAL :: lrg     ! regrid or raw data
  LOGICAL,          INTENT(IN),  OPTIONAL :: lstop   ! stop on error
  LOGICAL,          INTENT(OUT), OPTIONAL :: levent  ! return event status

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER       :: substr='RGTEVENT_READ_2D'
  REAL(dp), POINTER, DIMENSION(:,:,:,:) :: field => NULL()
  LOGICAL                           :: evtstat  ! rgtevent status
!!$  INTEGER                           :: status   ! memory status
  LOGICAL                           :: lok      ! regrid ok?
  CHARACTER(LEN=RGTMAXACTSTR)       :: actstr   ! rgtevent action
  INTEGER                           :: tnc      ! time step in netCDF file
  CHARACTER(LEN=GRD_MAXSTRLEN), &
              DIMENSION(:), POINTER :: hvars => NULL() ! variable names

  ! get/update event status
  CALL rgtevent_status(evtstat, tnc, rgt, modstr, &
    name=name, action=actstr,lstop=lstop)

  IF (PRESENT(levent)) levent = evtstat

  IF (.NOT.evtstat) RETURN

  ! read file
  SELECT CASE(RGREAD_TYPE)

  CASE(RGREAD_NCVAR)
     CALL RGTOOL_BI_READ_NCVAR(modstr                   &
          ,TRIM(actstr), tnc, field, lrg=lrg, lok=lok)

  CASE(RGREAD_NCFILE)
     IF (p_parallel_io) THEN
        CALL RGMSG(substr, RGMLW, &
             '2D-FIELD READ VIA RGTOOL_BI_READ_NCFILE!')
     END IF
     !
     CALL RGTOOL_BI_READ_NCFILE(modstr                         &
          ,TRIM(actstr), tnc, field, hvars, lrg=lrg, lok=lok)
     IF (lok) THEN
        IF (PRESENT(vars)) THEN
           IF (ASSOCIATED(vars)) THEN
!!$              DEALLOCATE(vars, STAT=status)
!!$              CALL ERRMSG(substr, status, 1)
              DEALLOCATE(vars)
           END IF
           NULLIFY(vars)
           IF (ASSOCIATED(hvars)) THEN
!!$              ALLOCATE(vars(SIZE(hvars)), STAT=status)
!!$              CALL ERRMSG(substr, status, 2)
              ALLOCATE(vars(SIZE(hvars)))
              vars(:) = hvars(:)
           END IF
        END IF
     END IF
     IF (ASSOCIATED(hvars)) THEN
!!$        DEALLOCATE(hvars, STAT=status)
!!$        CALL ERRMSG(substr, status, 3)
        DEALLOCATE(hvars)
     END IF
     NULLIFY(hvars)

  CASE DEFAULT
     IF (p_parallel_io) THEN
        CALL RGMSG(substr, RGMLE, &
             'UNKNOWN READ-TYPE OF RGTEVENT !',.false.)
!!$        CALL finish(substr)
        CALL error_bi('',substr)
     END IF
!     CALL finish(substr)

  END SELECT

  IF (lok) THEN
     IF ( (SIZE(efield,1) /= SIZE(field,1)) .OR. &
          (SIZE(efield,2) /= SIZE(field,4))      &
          ) THEN
        IF (p_parallel_io) THEN
           CALL RGMSG(substr, RGMLE, &
                'ARRAY SIZE MISMATCH !',.false.)
!!$           CALL finish(substr)
           CALL error_bi('',substr)
        END IF
!        CALL finish(substr)
     END IF
     efield(:,:) = field(:,1,1,:)
!!$     DEALLOCATE(field, STAT=status)
!!$     CALL ERRMSG(substr, status, 4)
     DEALLOCATE(field)
     NULLIFY(field)
  ELSE
     IF (p_parallel_io) THEN
        CALL RGMSG(substr, RGMLE, &
             'UPDATE OF/FROM '''//TRIM(actstr)//''' FAILED !',.false.)
!!$        CALL finish(substr)
        CALL error_bi('',substr)
     END IF
!     CALL finish(substr)
  END IF

END SUBROUTINE RGTEVENT_READ_2D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_READ_3D(rgt, modstr, name, RGREAD_TYPE   &
                           ,efield, lrg, lstop, vars, levent)

  ! ECHAM5
  USE messy_main_mpi_bi,        ONLY: p_parallel_io !!$, finish
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf,    ONLY: GRD_MAXSTRLEN

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, NULL, PRESENT, SIZE, TRIM

  ! I/O
  TYPE (RGTEVENT), DIMENSION(:), POINTER  :: rgt
  CHARACTER(LEN=*), INTENT(IN)            :: modstr        ! submodel-string
  CHARACTER(LEN=*), INTENT(IN)            :: name          ! rgtevent name
  INTEGER,          INTENT(IN)            :: RGREAD_TYPE   ! NCVAR or NCFILE
  REAL(dp), DIMENSION(:,:,:), INTENT(OUT) :: efield        ! 3D field
  CHARACTER(LEN=GRD_MAXSTRLEN),  OPTIONAL &
                   ,DIMENSION(:), POINTER :: vars    ! variable names
  LOGICAL,          INTENT(IN),  OPTIONAL :: lrg     ! regrid or raw data
  LOGICAL,          INTENT(IN),  OPTIONAL :: lstop   ! stop on error
  LOGICAL,          INTENT(OUT), OPTIONAL :: levent  ! return event status

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER       :: substr='RGTEVENT_READ_3D'
  REAL(dp), POINTER, DIMENSION(:,:,:,:) :: field => NULL()
  LOGICAL                           :: evtstat  ! rgtevent status
!!$  INTEGER                           :: status   ! memory status
  LOGICAL                           :: lok      ! regrid ok?
  CHARACTER(LEN=RGTMAXACTSTR)       :: actstr   ! rgtevent action
  INTEGER                           :: tnc      ! time step in netCDF file
  CHARACTER(LEN=GRD_MAXSTRLEN), &
              DIMENSION(:), POINTER :: hvars => NULL() ! variable names

  ! get/update event status
  CALL rgtevent_status(evtstat, tnc, rgt, modstr, &
    name=name, action=actstr,lstop=lstop)

  IF (PRESENT(levent)) levent = evtstat

  IF (.NOT.evtstat) RETURN

  ! read file
  SELECT CASE(RGREAD_TYPE)

  CASE(RGREAD_NCVAR)
     CALL RGTOOL_BI_READ_NCVAR(modstr                   &
          ,TRIM(actstr), tnc, field, lrg=lrg, lok=lok)

  CASE(RGREAD_NCFILE)
     CALL RGTOOL_BI_READ_NCFILE(modstr                         &
          ,TRIM(actstr), tnc, field, hvars, lrg=lrg, lok=lok)
     IF (lok) THEN
        IF (PRESENT(vars)) THEN
           IF (ASSOCIATED(vars)) THEN
!!$              DEALLOCATE(vars, STAT=status)
!!$              CALL ERRMSG(substr, status, 1)
              DEALLOCATE(vars)
           END IF
           NULLIFY(vars)
           IF (ASSOCIATED(hvars)) THEN
!!$              ALLOCATE(vars(SIZE(hvars)), STAT=status)
!!$              CALL ERRMSG(substr, status, 2)
              ALLOCATE(vars(SIZE(hvars)))
              vars(:) = hvars(:)
           END IF
        END IF
     END IF
     IF (ASSOCIATED(hvars)) THEN
!!$        DEALLOCATE(hvars, STAT=status)
!!$        CALL ERRMSG(substr, status, 3)
        DEALLOCATE(hvars)
     END IF
     NULLIFY(hvars)

  CASE DEFAULT
     IF (p_parallel_io) THEN
        CALL RGMSG(substr, RGMLE, &
             'UNKNOWN READ-TYPE OF RGTEVENT !',.false.)
!!$        CALL finish(substr)
        CALL error_bi('',substr)
     END IF
!     CALL finish(substr)

  END SELECT

  IF (lok) THEN
     IF (SIZE(efield,2) > 1) THEN
        IF ( (SIZE(efield,1) /= SIZE(field,1)) .OR. &
             (SIZE(efield,2) /= SIZE(field,2)) .OR. &
             (SIZE(efield,3) /= SIZE(field,4))      &
             ) THEN
           IF (p_parallel_io) THEN
              CALL RGMSG(substr, RGMLE, &
                   'ARRAY SIZE MISMATCH !',.false.)
!!$              CALL finish(substr)
              CALL error_bi('',substr)
           END IF
!           CALL finish(substr)
        END IF
        efield(:,:,:) = field(:,:,1,:)
     ELSE
        IF ( (SIZE(efield,1) /= SIZE(field,1)) .OR. &
             (SIZE(efield,2) /= SIZE(field,3)) .OR. &
             (SIZE(efield,3) /= SIZE(field,4))      &
             ) THEN
           IF (p_parallel_io) THEN
              CALL RGMSG(substr, RGMLE, &
                   'ARRAY SIZE MISMATCH !',.false.)
!!$              CALL finish(substr)
              CALL error_bi('',substr)
           END IF
!           CALL finish(substr)
        END IF
        efield(:,:,:) = field(:,1,:,:)
     END IF
!!$     DEALLOCATE(field, STAT=status)
!!$     CALL ERRMSG(substr, status, 4)
     DEALLOCATE(field)
     NULLIFY(field)
  ELSE
     IF (p_parallel_io) THEN
        CALL RGMSG(substr, RGMLE, &
             'UPDATE OF/FROM '''//TRIM(actstr)//''' FAILED !',.false.)
!!$        CALL finish(substr)
        CALL error_bi('',substr)
     END IF
!     CALL finish(substr)
  END IF

END SUBROUTINE RGTEVENT_READ_3D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_READ_4D(rgt, modstr, name, RGREAD_TYPE   &
                           ,efield, lrg, lstop, vars, levent)

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,        ONLY: p_parallel_io !!$, finish
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf,    ONLY: GRD_MAXSTRLEN

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, NULL, PRESENT, SIZE, TRIM

  ! I/O
  TYPE (RGTEVENT), DIMENSION(:), POINTER  :: rgt
  CHARACTER(LEN=*), INTENT(IN)            :: modstr        ! submodel-string
  CHARACTER(LEN=*), INTENT(IN)            :: name          ! rgtevent name
  INTEGER,          INTENT(IN)            :: RGREAD_TYPE   ! NCVAR or NCFILE
  REAL(dp), DIMENSION(:,:,:,:), INTENT(OUT) :: efield        ! 4D field
  CHARACTER(LEN=GRD_MAXSTRLEN),  OPTIONAL &
                   ,DIMENSION(:), POINTER :: vars    ! variable names
  LOGICAL,          INTENT(IN),  OPTIONAL :: lrg     ! regrid or raw data
  LOGICAL,          INTENT(IN),  OPTIONAL :: lstop   ! stop on error
  LOGICAL,          INTENT(OUT), OPTIONAL :: levent  ! return event status

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER       :: substr='RGTEVENT_READ_4D'
  REAL(dp), POINTER, DIMENSION(:,:,:,:) :: field => NULL()
  LOGICAL                           :: evtstat  ! rgtevent status
!!$  INTEGER                           :: status   ! memory status
  LOGICAL                           :: lok      ! regrid ok?
  CHARACTER(LEN=RGTMAXACTSTR)       :: actstr   ! rgtevent action
  INTEGER                           :: tnc      ! time step in netCDF file
  CHARACTER(LEN=GRD_MAXSTRLEN), &
              DIMENSION(:), POINTER :: hvars => NULL() ! variable names

  ! get/update event status
  CALL rgtevent_status(evtstat, tnc, rgt, modstr, &
    name=name, action=actstr,lstop=lstop)

  IF (PRESENT(levent)) levent = evtstat

  IF (.NOT.evtstat) RETURN

  ! read file
  SELECT CASE(RGREAD_TYPE)

  CASE(RGREAD_NCVAR)
     CALL RGTOOL_BI_READ_NCVAR(modstr                   &
          ,TRIM(actstr), tnc, field, lrg=lrg, lok=lok)

  CASE(RGREAD_NCFILE)
     CALL RGTOOL_BI_READ_NCFILE(modstr                         &
          ,TRIM(actstr), tnc, field, hvars, lrg=lrg, lok=lok)
     IF (lok) THEN
        IF (PRESENT(vars)) THEN
           IF (ASSOCIATED(vars)) THEN
!!$              DEALLOCATE(vars, STAT=status)
!!$              CALL ERRMSG(substr, status, 1)
              DEALLOCATE(vars)
           END IF
           NULLIFY(vars)
           IF (ASSOCIATED(hvars)) THEN
!!$              ALLOCATE(vars(SIZE(hvars)), STAT=status)
!!$              CALL ERRMSG(substr, status, 2)
              ALLOCATE(vars(SIZE(hvars)))
              vars(:) = hvars(:)
           END IF
        END IF
     END IF
     IF (ASSOCIATED(hvars)) THEN
!!$        DEALLOCATE(hvars, STAT=status)
!!$        CALL ERRMSG(substr, status, 3)
        DEALLOCATE(hvars)
     END IF
     NULLIFY(hvars)

  CASE DEFAULT
     IF (p_parallel_io) THEN
        CALL RGMSG(substr, RGMLE, &
             'UNKNOWN READ-TYPE OF RGTEVENT !',.false.)
!!$        CALL finish(substr)
        CALL error_bi('',substr)
     END IF
!     CALL finish(substr)

  END SELECT

  IF (lok) THEN
     IF ( (SIZE(efield,1) /= SIZE(field,1)) .OR. &
          (SIZE(efield,2) /= SIZE(field,2)) .OR. &
          (SIZE(efield,3) /= SIZE(field,3)) .OR. &
          (SIZE(efield,4) /= SIZE(field,4))      &
          ) THEN
        IF (p_parallel_io) THEN
           CALL RGMSG(substr, RGMLE, &
                'ARRAY SIZE MISMATCH !',.false.)
!!$           CALL finish(substr)
        CALL error_bi('',substr)
        END IF
!        CALL finish(substr)
     END IF
     efield(:,:,:,:) = field(:,:,:,:)
!!$     DEALLOCATE(field, STAT=status)
!!$     CALL ERRMSG(substr, status, 4)
     DEALLOCATE(field)
     NULLIFY(field)
  ELSE
     IF (p_parallel_io) THEN
        CALL RGMSG(substr, RGMLE, &
             'UPDATE OF/FROM '''//TRIM(actstr)//''' FAILED !',.false.)
!!$        CALL finish(substr)
        CALL error_bi('',substr)
     END IF
!     CALL finish(substr)
  END IF

END SUBROUTINE RGTEVENT_READ_4D
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_INIT_NML (rgt, modstr)

  ! INITIALIZES RGT-EVENTS FROM NAMELIST modstr.nml
  ! NAMELIST:
  !   RG_TRIG(1) = 1,'months','first',0,  'name_1',1, 1, 12, 1,'action1',
  !   RG_TRIG(2) = 1,'months','first',0,  'name_2',1, 1, 12, 1,'action2',
  !   RG_TRIG(3) = 1,'months','first',0,  'name_3',1, 1, 12, 1,'action3',
  !   ...
  !                |____________________| |___________________||________|
  !                          V                     V               V
  !                  I/O EVENT SYNTAX       COUNTER SYNTAX       ACTION
  !
  ! ACTION= variable-name or filename for regridding
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast !!$, finish
  USE messy_main_timer_bi,   ONLY: p_bcast_event
  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
  ! MESSY
  USE messy_main_tools,      ONLY: find_next_free_unit
  ! NCREGRID
  USE messy_ncregrid_base

  IMPLICIT NONE

  INTRINSIC :: TRIM

  ! I/O
  TYPE (RGTEVENT), DIMENSION(:), POINTER :: rgt
  CHARACTER(LEN=*), INTENT(IN)           :: modstr  ! submodel-string

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTEVENT_INIT_NML'
  INTEGER            :: iou                  ! logical I/O unit
  TYPE (RGTEVENT_IO) :: RG_TRIG(NMAXRGTE)    ! I/O RGT-events
  INTEGER            :: status               ! status
  INTEGER            :: n                    ! number of EVENT-COUNTERS
  INTEGER            :: i                    ! counter
  LOGICAL            :: lex                  ! file exists ?
  INTEGER            :: fstat                ! file status
  LOGICAL            :: lstat                ! event status
  INTEGER            :: cpos                 ! current counter position

  IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL rgtevent_read_nml_cpl(status, iou)
!!$     IF (status /= 0) CALL finish(substr)
        IF (status /= 0) CALL error_bi('',substr)
  END IF

  CALL start_message_bi(modstr,'RGTEVENT INITIALISATION',substr)
  
  ! GET NUMBER OF EVENTS
  IF (p_parallel_io) THEN
     n = 0
     DO i=1, NMAXRGTE
        IF (RG_TRIG(i)%cnt%start == -1) CYCLE
        n = n+1
        ! 
        IF (TRIM(RG_TRIG(i)%cnt%name) == '') THEN
           WRITE(RG_TRIG(i)%cnt%name(1:7),'(a3,i4.4)') 'RGT',i
        END IF
     END DO
     CALL RGMSG(substr, RGMLIC, &
          ' ',N,' RGTEVENTS RECOGNIZED IN '''//TRIM(modstr)//'.nml''')
  END IF
  ! BROADCAST RESULT
  CALL p_bcast(n, p_io)

  ! ALLOCATE SPACE
!!$  ALLOCATE(RGT(n), STAT=status)
!!$  CALL ERRMSG(substr,status,1)
  ALLOCATE(RGT(n))
  ! TRANSFER INFORMATION
  IF (p_parallel_io) THEN
     n = 0
     DO i=1, NMAXRGTE
        IF (RG_TRIG(i)%cnt%start == -1) CYCLE
        n = n+1
        RGT(n)%io = RG_TRIG(i)
     END DO
  END IF
  ! BROADCAST RESULTS
  DO i=1, n
     CALL p_bcast(rgt(i)%io%cnt%start,   p_io)
     CALL p_bcast(rgt(i)%io%cnt%step,    p_io)
     CALL p_bcast(rgt(i)%io%cnt%reset,   p_io)
     CALL p_bcast(rgt(i)%io%cnt%current, p_io)
     CALL p_bcast(rgt(i)%io%cnt%name,    p_io)
     CALL p_bcast(rgt(i)%io%act,         p_io)
     CALL p_bcast_event(rgt(i)%io%evt,   p_io)
     ! INITIALIZE EVENTS
     CALL RGTEVENT_INIT(rgt(i))
     ! UPDATE TO/FROM EVENT COUNTER LIST
     CALL RGTEVENT_STATUS(lstat, cpos, rgt, modstr, index=i, linit=.TRUE.)
  END DO

  CALL end_message_bi(modstr,'RGTEVENT INITIALISATION',substr)

  CONTAINS
    ! ------------------------------------------------------------------
    SUBROUTINE RGTEVENT_read_nml_cpl(status, iou)

      ! read namelist with RGT-events
      !
      ! Author: Patrick Joeckel, MPICH, Mar 2004
      
      ! MESSy
      USE messy_main_tools, ONLY: read_nml_open, read_nml_check &
                                , read_nml_close
      
      IMPLICIT NONE
      
      ! I/O
      INTEGER, INTENT(OUT) :: status     ! error status
      INTEGER, INTENT(IN)  :: iou        ! I/O unit
      
      ! (LOCAL) NAMELIST VARIABLES
      CHARACTER(LEN=*), PARAMETER :: substr = 'rgtevent_read_nml_cpl'
      
      NAMELIST /RGTEVENTS/ RG_TRIG
      
      ! LOCAL
      LOGICAL :: lex      ! file exists ?
      INTEGER :: fstat    ! file status
      
      status = 1
      
      ! INITIALIZE NAMELIST VARIABLES
!      DO i=1, NMAXRGTE
!         RG_TRIG(i)%cnt = NCRGCNT('', -1, -1, -1, -1)
!      END DO
      
      CALL read_nml_open(lex, substr, iou, 'RGTEVENTS', modstr)
      IF (.not.lex) RETURN    ! <modstr>.nml does not exist
      
      READ(iou, NML=RGTEVENTS, IOSTAT=fstat)
      CALL read_nml_check(fstat, substr, iou, 'RGTEVENTS', modstr)
      IF (fstat /= 0) RETURN  ! error while reading namelist
      
      CALL read_nml_close(substr, iou, modstr)
      
      status = 0  ! no ERROR
      
    END SUBROUTINE RGTEVENT_read_nml_cpl
    ! ------------------------------------------------------------------

END SUBROUTINE RGTEVENT_INIT_NML
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_INIT_ONL (rgt, modstr, index &
     , counter, unit, adjustment, offset   &
     , name, start, step, reset, current, action)

  ! INITIALIZES RGT-EVENT ONLINE
  !   RG_TRIG(1) = 1,'months','first',0,  'name_1',1, 1, 12, 1,'action1',
  !   RG_TRIG(2) = 1,'months','first',0,  'name_2',1, 1, 12, 1,'action2',
  !   RG_TRIG(3) = 1,'months','first',0,  'name_3',1, 1, 12, 1,'action3',
  !   ...
  !                |____________________| |___________________||________|
  !                          V                     V               V
  !                  I/O EVENT SYNTAX       COUNTER SYNTAX       ACTION
  !
  ! ACTION= variable-name or filename for regridding
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, June 2006

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast !!$, finish
  USE messy_main_timer_bi,   ONLY: p_bcast_event
  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
  ! MESSY
  USE messy_main_tools,      ONLY: find_next_free_unit
  ! NCREGRID
  USE messy_ncregrid_base

  IMPLICIT NONE

  INTRINSIC :: TRIM

  ! I/O
  TYPE (RGTEVENT),  DIMENSION(:), POINTER :: rgt
  CHARACTER(LEN=*), INTENT(IN) :: modstr  ! submodel-string
  INTEGER,          INTENT(IN) :: index
  !
  INTEGER,          INTENT(IN) :: counter    ! No. of steps in given unit
  CHARACTER(len=*), INTENT(IN) :: unit       ! counter unit type
  CHARACTER(len=*), INTENT(IN) :: adjustment ! adjustment in side the unit
  INTEGER,          INTENT(IN) :: offset
  !
  CHARACTER(LEN=*), INTENT(IN) :: name
  INTEGER,          INTENT(IN) :: start, step, reset, current
  CHARACTER(LEN=*), INTENT(IN) :: action

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTEVENT_INIT_ONL'
  INTEGER            :: status               ! status
  LOGICAL            :: lstat                ! event status
  INTEGER            :: cpos                 ! current counter position

  CALL start_message_bi(modstr,'RGTEVENT INITIALISATION',substr)

  rgt(index)%io%evt = io_time_event(counter, unit, adjustment, offset)
  rgt(index)%io%cnt%name    = TRIM(name)   ! name for identification of event
  rgt(index)%io%cnt%start   = start
  rgt(index)%io%cnt%step    = step
  rgt(index)%io%cnt%reset   = reset
  rgt(index)%io%cnt%current = current
  rgt(index)%io%act         = TRIM(action)

  ! INITIALIZE EVENTS
  CALL RGTEVENT_INIT(rgt(index))
  ! UPDATE TO/FROM EVENT COUNTER LIST
  CALL RGTEVENT_STATUS(lstat, cpos, rgt, modstr &
       , index=index, linit=.TRUE., lstop=.TRUE.)

  CALL end_message_bi(modstr,'RGTEVENT INITIALISATION',substr)

END SUBROUTINE RGTEVENT_INIT_ONL
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_STATUS(status, cpos, rgt, modstr, name, index, action &
                          , lstop, linit)

  ! RETURNS THE EVENT STATUS (flag) OF A NAMED (name)
  ! OR INDEXED (index) RGT-EVENT IN A LIST (rgt)
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  ! ECHAM5/MESSy
  USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_bcast !!$, finish
  USE messy_main_timer,     ONLY: lresume, lstart
  ! NCREGRID
  USE messy_ncregrid_tools, ONLY: rgtool_ncrgcnt_rst
  USE messy_ncregrid_base

  IMPLICIT NONE

  INTRINSIC :: PRESENT, SIZE, TRIM

  ! I/O
  LOGICAL,                      INTENT(OUT)           :: status ! event status
  INTEGER,                      INTENT(OUT)           :: cpos   ! counter pos.
  TYPE(RGTEVENT), DIMENSION(:), POINTER               :: rgt    ! RGT-event list
  CHARACTER(LEN=*),             INTENT(IN)            :: modstr ! calling module
  CHARACTER(LEN=*),             INTENT(IN), OPTIONAL  :: name   ! name of event
  INTEGER,                      INTENT(IN), OPTIONAL  :: index  ! index of event
  CHARACTER(LEN=RGTMAXACTSTR),  INTENT(OUT), OPTIONAL :: action ! action string
  LOGICAL,                      INTENT(IN), OPTIONAL  :: lstop  ! stop on error
  LOGICAL,                      INTENT(IN), OPTIONAL  :: linit  ! initialize?

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTEVENT_STATUS'
  INTEGER :: ix     ! INDEX
  LOGICAL :: llstop ! stop on error
  INTEGER :: RGMLX  ! RGML(W,E)
  LOGICAL :: zlinit

  ! INIT
  status = .false.
  cpos = 0
  IF (PRESENT(action)) action = ''

  ! STOP ON ERROR ?
  IF (PRESENT(lstop)) THEN
    llstop = lstop
  ELSE
    llstop = .true.  ! DEFAULT
  END IF
  !
  IF (llstop) THEN
    RGMLX = RGMLE
  ELSE
    RGMLX = RGMLW
  END IF
  !
  IF (PRESENT(linit)) THEN
     zlinit = linit
  ELSE
     zlinit = .FALSE. ! DEFAULT
  END IF
  
  ! CHECK SUBROUTINE CALL
  IF (((.NOT.PRESENT(name)).AND.(.NOT.PRESENT(index))) .OR. &
       (PRESENT(name).AND.PRESENT(index))) THEN
     IF (p_parallel_io) THEN
        CALL RGMSG(substr, RGMLE, &
             'EITHER NAME OR INDEX MUST BE GIVEN !', .false.)
!        CALL finish(substr)
     END IF
!!$     CALL finish(substr)
     CALL error_bi('',substr)
  END IF

  ! CALLED BY NAME: LOOK FOR INDEX
  IF (PRESENT(name)) THEN
     CALL RGTEVENT_INDEX(rgt, name, ix)
     IF (ix < 0) THEN
        IF (p_parallel_io) THEN
           CALL RGMSG(substr, RGMLX, &
                'RGTEVENT WITH NAME '''//TRIM(name)//&
                &''' NOT FOUND IN LIST !', .false.)
!           IF (llstop) CALL finish(substr)
        END IF
!!$        IF (llstop) CALL finish(substr)
        IF (llstop) CALL error_bi('',substr)
        RETURN
     END IF
  END IF

  ! CALLED BY INDEX
  IF (PRESENT(index)) THEN
     IF (index > SIZE(rgt)) THEN
        IF (p_parallel_io) THEN
           CALL RGMSG(substr, RGMLX, 'INDEX OUT OF RANGE !', .false.)
!           IF (llstop) CALL finish(substr)
        END IF
!!$        IF (llstop) CALL finish(substr)
        IF (llstop) CALL error_bi('',substr)
        RETURN
     ELSE
        ix = index
     END IF
  END IF

  ! GET EVENT STATUS
  IF (.NOT. zlinit) CALL RGTEVENT_STAT(rgt(ix), status)

  ! UPDATE TO/FROM COUNTER LIST
  CALL RGTOOL_NCRGCNT_RST(modstr, lstart, lresume, &
       status, rgt(ix)%io%cnt, p_parallel_io, linit)

  ! RETURN VALUES
  cpos = rgt(ix)%io%cnt%current
  status = status.OR.lstart.OR.lresume
  IF (PRESENT(action)) THEN
     action = TRIM(rgt(ix)%io%act)
  END IF

END SUBROUTINE RGTEVENT_STATUS
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE messy_ncregrid_write_restart

  ! ECHAM5/MESSy
  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
  USE messy_main_mpi_bi,     ONLY: p_parallel_io !!$, finish
  ! MESSy
  USE messy_main_tools,      ONLY: find_next_free_unit
  USE messy_ncregrid_tools,  ONLY: write_ncrgcnt_list

  IMPLICIT NONE

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'messy_ncregrid_write_restart'
  INTEGER          :: iou
  INTEGER          :: status
  INTEGER, SAVE    :: nrstcount     = 0
  CHARACTER(LEN=4) :: nrstcount_str = ''

  CALL start_message_bi('NCREGRID_TOOLS','WRITING RESTART FILE',substr)

  nrstcount = nrstcount + 1
  write(nrstcount_str,'(i4.4)') nrstcount

  IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL write_ncrgcnt_list(status, &
          'restart_'//nrstcount_str//'_ncregrid.rst', iou)
!!$     IF (status /= 0) CALL finish(substr, &
!!$          'WRITE_NCRGCNT_LIST REPORTED AN ERROR!')
     IF (status /= 0) CALL error_bi(&
          'WRITE_NCRGCNT_LIST REPORTED AN ERROR!', substr)
  END IF

  CALL end_message_bi('NCREGRID_TOOLS','WRITING RESTART FILE',substr)

END SUBROUTINE messy_ncregrid_write_restart
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE messy_ncregrid_read_restart

  ! ECHAM5/MESSy
  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
  USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast !!$, finish
  ! MESSy
  USE messy_main_tools,      ONLY: find_next_free_unit
  USE messy_ncregrid_tools,  ONLY: read_ncrgcnt_list, get_next_ncrgcnt, ncrgcnt

  IMPLICIT NONE

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'messy_ncregrid_read_restart'
  INTEGER                :: iou
  TYPE(ncrgcnt), POINTER :: cntptr
  LOGICAL                :: last
  INTEGER                :: status

  CALL start_message_bi('NCREGRID_TOOLS','READING RESTART FILE',substr)

  IF (p_parallel_io) THEN
     iou = find_next_free_unit(100,200)
     CALL read_ncrgcnt_list(status, 'restart_ncregrid.rst', iou)
!!$     IF (status /= 0) CALL finish(substr, &
!!$          'READ_NCRGCNT_LIST REPORTED AN ERROR!')
     IF (status /= 0) CALL error_bi(&
          'READ_NCRGCNT_LIST REPORTED AN ERROR!', substr)
  END IF

  DO
     CALL get_next_ncrgcnt(last, cntptr)
     IF (last) EXIT
     CALL p_bcast(cntptr%start,   p_io)
     CALL p_bcast(cntptr%step,    p_io)
     CALL p_bcast(cntptr%reset,   p_io)
     CALL p_bcast(cntptr%current, p_io)
  END DO

  CALL end_message_bi('NCREGRID_TOOLS','READING RESTART FILE',substr )

END SUBROUTINE messy_ncregrid_read_restart
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE messy_ncregrid_free_memory

  ! ECHAM5/MESSy
!!$  USE messy_main_mpi_bi,    ONLY: finish
  ! MESSy
  USE messy_ncregrid_tools, ONLY: clean_ncrgcnt_list
  
  IMPLICIT NONE

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'messy_ncregrid_free_memory'
  INTEGER :: status

  CALL clean_ncrgcnt_list(status)
!!$  IF (status /= 0) CALL finish(substr, &
!!$       'clean_ncrgcnt_list reported an error')
  IF (status /= 0) CALL error_bi(&
       'clean_ncrgcnt_list reported an error', substr)

END SUBROUTINE messy_ncregrid_free_memory
! ------------------------------------------------------------------

! ####### PRIVATE ROUTINES #########################################

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_INIT(rgt)

  ! INITIALIZES ONE INTERNAL RGT-EVENT (rgt) FROM
  ! I/O NAMELIST INFORMATION
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  USE messy_main_timer_bi,          ONLY: timer_event_init

  IMPLICIT NONE

  INTRINSIC :: TRIM

  ! I/O
  TYPE (RGTEVENT), INTENT(INOUT) :: rgt

  CALL timer_event_init(rgt%event, rgt%io%evt, &
       TRIM(rgt%io%cnt%name), 'next')

END SUBROUTINE RGTEVENT_INIT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_STAT(rgt, flag)

  ! RETURNS EVENT STATUS (flag) OF ONE RGT-EVENT
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  ! BM/MESSy
  USE messy_main_timer_bi,   ONLY: event_state
  ! MESSy
  USE messy_main_timer,      ONLY: next_date

  IMPLICIT NONE

  ! I/O
  TYPE (RGTEVENT), INTENT(INOUT) :: rgt
  LOGICAL,         INTENT(OUT)   :: flag

  flag = event_state(rgt%event, next_date)

END SUBROUTINE RGTEVENT_STAT
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTEVENT_INDEX(rgt, name, ix)

  ! RETURNS INDEX (ix) OF NAMED (name) RGT-EVENT IN LIST (rgt)
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  IMPLICIT NONE

  INTRINSIC :: SIZE, TRIM

  ! I/O
  TYPE (RGTEVENT),DIMENSION(:), INTENT(IN)   :: rgt   ! counter-list
  CHARACTER(LEN=*),             INTENT(IN)   :: name  ! counter-name
  INTEGER,                      INTENT(OUT)  :: ix    ! index

  ! LOCAL
  INTEGER :: i

  ix = -1
  DO i=1, SIZE(rgt)
     IF (TRIM(name) == TRIM(rgt(i)%io%cnt%name)) THEN
        ix = i
        EXIT
     END IF
  END DO

END SUBROUTINE RGTEVENT_INDEX
! ------------------------------------------------------------------

#endif
! ... of #ifndef ICON

! ******************************************************************
! ------------------------------------------------------------------
END MODULE MESSY_NCREGRID_TOOLS_BI
! ------------------------------------------------------------------
! ******************************************************************
