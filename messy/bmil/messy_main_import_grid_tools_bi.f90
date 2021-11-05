#include "messy_main_ppd_bi.inc"

! ******************************************************************
MODULE messy_main_import_grid_tools_bi
! ******************************************************************

  ! ******************************************************************
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002
  !         extended for curvilinear ocean grids by
  !         Astrid  Kerkweg, UniMz, Mainz, 2012-2013
  !         -> re-structured/expanded for more general GRID application
  ! ******************************************************************

  USE messy_main_timer_event,   ONLY: time_event, io_time_event   &
                                    , TIME_INC_MONTHS, TRIG_EXACT
  USE messy_main_blather_bi,    ONLY: ERROR_BI
  USE messy_main_import_grid,   ONLY: t_ncrgcnt
  USE messy_main_channel,       ONLY: STRLEN_OBJECT
  USE messy_main_grid_trafo,    ONLY: GTRF_NONE &
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
                                    , GTRF_NRGD &
#endif
                                    , GTRF_SCRP
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  ! MAX. NUMBER OF RGTEVENTS PER SUB-MODEL
  INTEGER, PARAMETER         :: NMAXRGTE     = 500
  INTEGER, PARAMETER, PUBLIC :: RGTMAXACTSTR = 800

  ! READ IN OPTION
  INTEGER, PARAMETER, PUBLIC :: RGREAD_NCVAR =1 ! USE RGTOOL_BI_READ_NCVAR
  INTEGER, PARAMETER, PUBLIC :: RGREAD_NCFILE=2 ! USE RGTOOL_BI_READ_NCFILE

  ! DEFINE DEFAULT order
  CHARACTER(LEN=4), PARAMETER :: def_order = 'xyzn'

  ! USED FOR NAMELIST-INPUT RGTINIT (ONLY INTERNALLY)
  TYPE t_rgtinit_io
     CHARACTER(LEN=STRLEN_OBJECT) :: name = ''
     INTEGER                      :: tstep
     CHARACTER(LEN=RGTMAXACTSTR)  :: act = ''
  END TYPE t_rgtinit_io

  TYPE t_rgtinit
     CHARACTER(LEN=STRLEN_OBJECT) :: name = ''
     CHARACTER(LEN=RGTMAXACTSTR)  :: act = ''
     INTEGER                      :: domain_idx = 0
     INTEGER                      :: tstep
  END TYPE t_rgtinit

  ! USED FOR NAMELIST-INPUT RGTEVENT (ONLY INTERNALLY)
  TYPE t_rgtevent_io
     TYPE(io_time_event)         :: evt = &
          io_time_event(1, TIME_INC_MONTHS,TRIG_EXACT,0)
     TYPE(t_ncrgcnt)             :: cnt = &
          T_NCRGCNT('', -1, -1, -1, -1)
     CHARACTER(LEN=RGTMAXACTSTR) :: act = ''
  END TYPE t_rgtevent_io

  ! INTERNAL WORKSPACE
  TYPE t_rgtevent
     TYPE(t_rgtevent_io) :: io
     TYPE(time_event)    :: event
     INTEGER             :: domain_idx = 0
     REAL(DP), POINTER   :: n => NULL()  ! pointer to SCALAR channel object
  END TYPE t_rgtevent

  ! PUBLIC, BUT NOT DIRECT USER INTERFACE
  ! NEEDED TO ACCESS SUB-STRUCTURES OF RGTEVENT
  PUBLIC :: t_rgtevent_io ! (only public, because component of RGTEVENT !)
  PUBLIC :: t_ncrgcnt     ! (only public, because component of RGTEVENT !)

  ! USER INTERFACE
  ! - DATA IMPORT
  INTERFACE RGTOOL_BI_READ_NCVAR
     MODULE PROCEDURE RGTOOL_BI_READ_NCVAR_5D
  END INTERFACE
  PUBLIC  :: RGTOOL_BI_READ_NCVAR  ! READ ONE FIELD WITH NCREGRID

  INTERFACE RGTOOL_BI_READ_NCFILE
     MODULE PROCEDURE RGTOOL_BI_READ_NCFILE_5D
  END INTERFACE
  PUBLIC  :: RGTOOL_BI_READ_NCFILE ! READ ONE FILE WITH NCREGRID

  ! - INIT HANDLING
  PUBLIC  :: T_RGTINIT             ! RG-TRIGGER EVENT STRUCTURE
  PUBLIC  :: RGTINIT_INIT_NML      ! READ NAMELIST FOR STATIC DATA (read in
                                   ! INITAL-PHASE ONLY
  ! - EVENT/COUNTER HANDLING
  PUBLIC  :: T_RGTEVENT            ! RG-TRIGGER EVENT STRUCTURE
  PUBLIC  :: RGTEVENT_INIT_NML     ! READ RGT-EVENT INFORMATION FROM NAMELIST
  PUBLIC  :: RGTEVENT_STATUS       ! GET STATUS OF NAMED/INDEXED EVENT AND
  ! UPDATE RESTART FILE
  ! - COUNTER RESTART HANDLING
  PUBLIC  :: RGTEVENT_TO_CHANNEL   ! SAVE RGT-EVENT INFORMATION TO CHANNEL
  PUBLIC  :: RGTEVENT_FROM_CHANNEL ! INIT COUNTER INFORMATION FROM CHANNEL

  ! PRIVATE HELPER ROUTINES
  !  PRIVATE :: RGTEVENT_INIT        ! INITIALIZE ONE RGT-EVENT FROM I/O
  !  PRIVATE :: RGTEVENT_STAT        ! CHECK STATUS OF ONE RGT-EVENT
  !  PRIVATE :: ncvar2choatt         ! convert netCDF variable attributes
  !                                  ! into channel object attribute list

CONTAINS

  ! ####### PUBLIC ROUTINES ##########################################

  ! ------------------------------------------------------------------
  SUBROUTINE RGTOOL_BI_READ_NCVAR_5D(modstr, vname, t, dat    &
       , iipol, rg_range, ogrid_ID, igrid_id  &
       , SDID                       &  ! OPTIONAL
       , lrg, lrgx, lrgy, lrgz, lok &  ! OPTIONAL
       , hyam, hybm, p0, ps         &  ! OPTIONAL
       , hyai, hybi                 &  ! OPTIONAL
       , latm, lonm, lati, loni     &  ! OPTIONAL
       , oarea                      &  ! OPTIONAL
       , att                        &  ! OPTIONAL
       , ldatreq                    &  ! OPTIONAL
       , nlon, nlat, nlev, npar     &  ! OPTIONAL
       , mnc)                          ! OPTIONAL

    ! PERFORMS ONE INTERPOLATION STEP FOR ONE FIELD
    ! AND RETURNS m DATA AS 5D ARRAY (x,z,m,y,n),
    ! AND OPTIONALLY GRID-STRUCTURE AS ARRAYS
    !
    ! Author: Patrick Joeckel, MPICH, Mainz, October 2002    (NCREGRID)
    !         Astrid Kerkweg, Uni Mainz, 2012/2013 (expanded for SCRIP)

    USE messy_main_mpi_bi,          ONLY: p_pe, p_nprocs, p_all_comm
    USE messy_main_tools,           ONLY: FIND_NEXT_FREE_UNIT
    USE messy_main_grid_netcdf,     ONLY: t_ncvar, INIT_NCVAR, t_multinc
    USE messy_main_grid,            ONLY: t_geohybgrid &
                                        , GRID_ERROR   &
                                        , INIT_GEOHYBGRID
    USE messy_main_grid_tools,      ONLY: RGTOOL_CONVERT, RGTOOL_G2C
    USE messy_main_import_grid,     ONLY: RGTOOL_READ_NCVAR
    USE messy_main_import_grid_par, ONLY: INIT_PARALLEL
    USE messy_main_channel_attributes, ONLY: t_att_ptr
    USE messy_main_channel_repr,       ONLY: repr_def_axes
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma
    USE messy_main_mpi_bi,          ONLY: dcg, P_BCAST, p_io, p_parallel_io &
                                        , scatter_gp
    USE messy_main_grid_bi,         ONLY: P_BCAST_GRID
    USE messy_main_channel_bi,      ONLY: P_BCAST_ATTRIBUTE
    USE messy_main_grid,            ONLY: NEW_GEOHYBGRID
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT, SIZE, TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)             :: modstr    ! calling module
    CHARACTER(LEN=*), INTENT(IN)             :: vname     ! name of variable
    INTEGER,          INTENT(IN)             :: t         ! netCDF time step
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER  :: dat       ! (local) data field
    INTEGER,          INTENT(IN)             :: iipol     ! interpolation method
    INTEGER                       , POINTER  :: rg_range(:,:)
    INTEGER,          INTENT(INOUT)          :: ogrid_ID  ! output grid ID
    INTEGER,          INTENT(INOUT)          :: igrid_ID  !  input grid ID
    INTEGER,          INTENT(INOUT),OPTIONAL :: SDID      ! SCRIP grid ID
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
    REAL(dp), DIMENSION(:,:), POINTER, OPTIONAL :: oarea  ! box area output grid
    TYPE(t_att_ptr), DIMENSION(:), POINTER, OPTIONAL :: att ! var attributes
    LOGICAL, INTENT(IN),  OPTIONAL :: ldatreq
    INTEGER, INTENT(OUT), OPTIONAL :: nlon
    INTEGER, INTENT(OUT), OPTIONAL :: nlat
    INTEGER, INTENT(OUT), OPTIONAL :: nlev
    INTEGER, INTENT(OUT), OPTIONAL :: npar
    TYPE(t_multinc), INTENT(INOUT), OPTIONAL :: mnc

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER            :: substr = 'RGTOOL_BI_READ_NCVAR_5D'
    INTEGER                                :: iou       ! logical I/O unit
    INTEGER                                :: status    ! memory status
    TYPE(t_ncvar)                          :: var       ! nc-variable
    REAL(dp), DIMENSION(:,:,:,:), POINTER  :: gdat      ! (global) data field
    REAL(dp), DIMENSION(:,:),     POINTER  :: gps       ! surface pressure
    ! (local) data pointer
    REAL(dp), DIMENSION(:,:,:,:), POINTER  :: zlptr => NULL()
    TYPE(t_geohybgrid)                     :: convgrid  ! grid-structure
    LOGICAL                                :: llrg, llrgx, llrgy, llrgz, llok
    INTEGER                                :: klat
    INTEGER                                :: klon
    INTEGER                                :: klev
    INTEGER                                :: kparam
    CHARACTER(LEN=4)                       :: dat_order = def_order
    INTEGER                                :: idx_lon, idx_lat, idx_lev
    INTEGER                                :: idx_param
    INTEGER                                :: ix
    INTEGER                                :: SCRIP_ID
    INTEGER , SAVE                         :: callnum = 0
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    INTEGER                                :: s1d  ! size of 1D field
    INTEGER                                :: natt ! number of attributes
    LOGICAL                                :: lsc  ! scatter results ?
    INTEGER                                :: ia   ! attribute counter
#endif
    INTEGER                                :: i
    LOGICAL                                :: zldatreq

    ! INITIALIZE
    ! I/O
    callnum = callnum + 1

    IF (ASSOCIATED(dat)) THEN
       DEALLOCATE(dat)
    END IF
    NULLIFY(dat)
    !
    IF (ASSOCIATED(rg_range)) THEN
       DEALLOCATE(rg_range)
    END IF
    NULLIFY(rg_range)
    ALLOCATE(rg_range(1,2))
    !
    IF (PRESENT(ps)) THEN
       IF (ASSOCIATED(ps)) THEN
          DEALLOCATE(ps)
       END IF
       NULLIFY(ps)
    END IF
    !
    IF (PRESENT(hyam)) THEN
       IF (ASSOCIATED(hyam)) THEN
          DEALLOCATE(hyam)
       END IF
       NULLIFY(hyam)
    END IF
    !
    IF (PRESENT(hybm)) THEN
       IF (ASSOCIATED(hybm)) THEN
          DEALLOCATE(hybm)
       END IF
       NULLIFY(hybm)
    END IF
    !
    IF (PRESENT(hyai)) THEN
       IF (ASSOCIATED(hyai)) THEN
          DEALLOCATE(hyai)
       END IF
       NULLIFY(hyai)
    END IF
    !
    IF (PRESENT(hybi)) THEN
       IF (ASSOCIATED(hybi)) THEN
          DEALLOCATE(hybi)
       END IF
       NULLIFY(hybi)
    END IF
    !
    IF (PRESENT(p0)) THEN
       IF (ASSOCIATED(p0)) THEN
          DEALLOCATE(p0)
       END IF
       NULLIFY(p0)
    END IF
    !
    IF (PRESENT(latm)) THEN
       IF (ASSOCIATED(latm)) THEN
          DEALLOCATE(latm)
       END IF
       NULLIFY(latm)
    END IF
    !
    IF (PRESENT(lonm)) THEN
       IF (ASSOCIATED(lonm)) THEN
          DEALLOCATE(lonm)
       END IF
       NULLIFY(lonm)
    END IF
    !
    IF (PRESENT(lati)) THEN
       IF (ASSOCIATED(lati)) THEN
          DEALLOCATE(lati)
       END IF
       NULLIFY(lati)
    END IF
    !
    IF (PRESENT(loni)) THEN
       IF (ASSOCIATED(loni)) THEN
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
    IF (PRESENT(SDID)) THEN
       SCRIP_ID = SDID
    ELSE
       SCRIP_ID = -99
    ENDIF
    !
    IF (PRESENT(ldatreq)) THEN
       zldatreq = ldatreq
    ELSE
       zldatreq = .TRUE.
    END IF
    !
    NULLIFY(gdat)       ! (global) data field
    NULLIFY(gps)        ! (global) surface pressure field

    IF (PRESENT(att)) THEN
       IF (ASSOCIATED(att)) THEN
          DO i=1, SIZE(att)
             IF (ASSOCIATED(att(i)%ptr)) DEALLOCATE(att(i)%ptr)
             NULLIFY(att(i)%ptr)
          END DO
          DEALLOCATE(att)
          NULLIFY(att)
       END IF
    END IF

    CALL INIT_PARALLEL(p_pe, p_nprocs, p_all_comm)

#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    IF ( p_parallel_io .OR. (iipol == GTRF_NONE)) THEN
#endif
       iou = find_next_free_unit(100,200)
       CALL RGTOOL_READ_NCVAR(status, iou &
            , TRIM(modstr)//'.nml', vname, t, var, iipol &
            , ogrid_ID, igrid_ID, SCRIP_ID               &
            , llrg, llrgx, llrgy, llrgz, llok, oarea,convgrid &
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
            , ldompar = .FALSE., lvarpar = .FALSE. &
#else
            , ldompar = .TRUE. , lvarpar = .FALSE. &
#endif
            , mnc=mnc)
       IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)

       IF (llok) THEN
          rg_range(1,:) = 0
          DO ix=1,var%natts
             IF (TRIM(var%att(ix)%name) == 'mmig_ixf_range') THEN
                rg_range(1,1:2) = INT(var%att(ix)%dat%vi(1:2))
             END IF
          END DO
          IF (.NOT. llrg) THEN
             IF (iipol == GTRF_NONE) THEN
                ogrid_ID = igrid_ID
             ENDIF
          ENDIF

          IF (iipol == GTRF_NONE) THEN
             dat_order = 'xyzn'
             idx_lon   = 1
             idx_lat   = 2
             idx_lev   = 3
             idx_param = 4
          ELSE
             dat_order=repr_def_axes(_RI_XYZN_('x','y','z','n'))
             idx_lon   = _IX_XYZN_
             idx_lat   = _IY_XYZN_
             idx_lev   = _IZ_XYZN_
             idx_param = _IN_XYZN_
          END IF

          !CALL PRINT_GEOHYBGRID(convgrid, 'CONVGRID')

          CALL RGTOOL_CONVERT(var, gdat, convgrid, order=dat_order)

          !write (0,*) 'IG: NCVAR: GDAT', UBOUND(gdat)

          klon   = SIZE(gdat, idx_lon)
          klat   = SIZE(gdat, idx_lat)
          klev   = SIZE(gdat, idx_lev)
          kparam = SIZE(gdat, idx_param)
          IF (PRESENT(nlon)) nlon = klon
          IF (PRESENT(nlat)) nlat = klat
          IF (PRESENT(nlev)) nlev = klev
          IF (PRESENT(npar)) npar = kparam

          CALL RGTOOL_G2C(convgrid, hyam, hybm, p0, gps      &
               ,hyai, hybi                    &
               ,latm, lonm                    &
               ,lati, loni, klat, klon, klev, order=dat_order)

          IF (PRESENT(att)) THEN
             ALLOCATE(att(1))
             CALL NCVAR2CHOATT(var, att(1)%ptr)
          ENDIF

          CALL INIT_NCVAR(var)
       END IF ! llok

#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    END IF ! p_parallel_io

    ! BROADCAST SUCCESS
    CALL P_BCAST(llok, p_io)
    IF (.NOT. llrg .AND. iipol /= GTRF_NONE) THEN
       CALL P_BCAST_GRID(convgrid, p_io)
       IF (.NOT. p_parallel_io) THEN
          CALL NEW_GEOHYBGRID(status, igrid_ID, convgrid)
       END IF
    END IF

    IF (llok) THEN
       IF (PRESENT(att)) THEN
          IF (.NOT. p_parallel_io) THEN
             ALLOCATE(att(1))
          ELSE
             IF (ASSOCIATED(att(1)%ptr)) THEN
                natt = SIZE(att(1)%ptr)
             ELSE
                natt = 0
             END IF
          END IF
          CALL P_BCAST(natt, p_io)
          IF (.NOT. p_parallel_io) THEN
             IF (natt > 0) THEN
                ALLOCATE(att(1)%ptr(natt))
             END IF
          END IF
          DO ia=1, natt
             CALL P_BCAST_ATTRIBUTE(att(1)%ptr(ia), p_io)
          END DO
       ENDIF
    ENDIF
#endif

    CALL INIT_GEOHYBGRID(convgrid) ! memory leak fixed

    IF (PRESENT(lok)) THEN
       lok = llok
    END IF

    SELECT CASE (iipol)
    CASE (GTRF_NONE)
       ! NOTHING TO DO
       IF (llok) THEN
          IF (zldatreq) THEN
             ALLOCATE(dat(SIZE(gdat,1),SIZE(gdat,2),SIZE(gdat,3) &
                  , SIZE(gdat,4),1))
             ! COPY FIRST ON I/O-PROCESSOR
             dat(:,:,:,:,1) = gdat(:,:,:,:)
          END IF
       END IF

#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
     CASE (GTRF_NRGD, GTRF_SCRP)
        ! BROADCST DIMENSIONS
       IF (llok) THEN
          CALL P_BCAST(klon, p_io)
          CALL P_BCAST(klev, p_io)
          CALL P_BCAST(klat, p_io)
          CALL P_BCAST(kparam, p_io)
          IF (PRESENT(nlon)) nlon = klon
          IF (PRESENT(nlat)) nlat = klat
          IF (PRESENT(nlev)) nlev = klev
          IF (PRESENT(npar)) npar = kparam
       END IF

       IF (zldatreq) THEN

          ! NOTE: SCATTERING THE RESULTS (dat) TO THE PROCESSORS IS ONLY
          !       POSSIBLE
          !       IF REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy.AND.lrgz)
          lsc = (llrg.AND.llok.AND.llrgx.AND.llrgy)
          !
          IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
             ALLOCATE(dat(nproma,klev,kparam,ngpblks,1))
             ! SCATTER GLOBAL DATA TO LOCAL FIELDS
             zlptr => dat(:,:,:,:,1)
             CALL SCATTER_GP(gdat, zlptr, dcg)
          ELSE            ! COPY AND BROADCAST ENTIRE FIELD
             IF (llok) THEN
                ! ALLOCATE SPACE FOR DATA FIELD
                ALLOCATE(dat(klon,klev,kparam,klat,1))
                ! COPY FIRST ON I/O-PROCESSOR
                IF (p_parallel_io) THEN
                   dat(:,:,:,:,1) = gdat(:,:,:,:)
                END IF
                ! BROADCAST RESULT
                CALL P_BCAST(dat(:,:,:,:,1), p_io)
             END IF ! llok
          END IF ! lsc

       END IF

       ! NOTE: SCATTERING THE RESULTS (ps) TO THE PROCESSORS IS ONLY POSSIBLE
       !       IF HORIZONTAL REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy)
       IF (PRESENT(ps)) THEN
          lsc = llrg.AND.llok.AND.llrgx.AND.llrgy
          !
          IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
             ! ALLOCATE SPACE FOR LOCAL DATA FIELD
             ALLOCATE(ps(nproma,ngpblks))
             ! SCATTER GLOBAL DATA TO LOCAL FIELDS
             CALL SCATTER_GP(gps, ps, dcg)
          ELSE
             IF (llok) THEN
                ! ALLOCATE SPACE FOR DATA FIELD
                ALLOCATE(ps(klon,klat))
                IF (p_parallel_io) THEN
                   ps(:,:) = gps(:,:)
                END IF
                ! BROADCAST RESULT
                DO ix = 1, SIZE(ps,2)
                   CALL P_BCAST(ps(:,ix) , p_io)
                END DO
             END IF ! llok
          END IF ! lsc
       END IF ! ps present

       ! 1D ARRAYS NEVER SCATTERED
       IF (llok) THEN
          IF (.NOT. p_parallel_io) THEN
             ALLOCATE(rg_range(1,2))
          END IF
          CALL P_BCAST(rg_range, p_io)
          IF (PRESENT(hyam)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(hyam)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(hyam(s1d))
             END IF
             CALL P_BCAST(hyam ,p_io)
          END IF
          !
          IF (PRESENT(hybm)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(hybm)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(hybm(s1d))
             END IF
             CALL P_BCAST(hybm ,p_io)
          END IF
          !
          IF (PRESENT(hyai)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(hyai)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(hyai(s1d))
             END IF
             CALL P_BCAST(hyai ,p_io)
          END IF
          !
          IF (PRESENT(hybi)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(hybi)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(hybi(s1d))
             END IF
             CALL P_BCAST(hybi ,p_io)
          END IF
          !
          IF (PRESENT(p0)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(p0)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(p0(s1d))
             END IF
             CALL P_BCAST(p0 ,p_io)
          END IF
          !
          IF (PRESENT(latm)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(latm)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(latm(s1d))
             END IF
             CALL P_BCAST(latm ,p_io)
          END IF
          !
          IF (PRESENT(lonm)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(lonm)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(lonm(s1d))
             END IF
             CALL P_BCAST(lonm ,p_io)
          END IF
          !
          IF (PRESENT(lati)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(lati)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(lati(s1d))
             END IF
             CALL P_BCAST(lati ,p_io)
          END IF
          !
          IF (PRESENT(loni)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(loni)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN

                ALLOCATE(loni(s1d))
             END IF
             CALL P_BCAST(loni ,p_io)
          END IF
       END IF
#else
    CASE (GTRF_SCRP)
       ! ALLOCATE SPACE FOR DATA FIELD
       IF (zldatreq) THEN
          ALLOCATE(dat(SIZE(gdat,1),SIZE(gdat,2),SIZE(gdat,3),SIZE(gdat,4),1))
          dat(:,:,:,:,1) = gdat(:,:,:,:)
       END IF

       IF (PRESENT(ps)) THEN
          ALLOCATE(ps(klon,klat))
          ps(:,:) = gps(:,:)
       END IF
#endif
    END SELECT

    IF (PRESENT(SDID)) THEN
       SDID = SCRIP_ID
    ENDIF

    ! CLEAN UP: FREE MEMORY OF GLOBAL FIELDS
    IF (llok) THEN
       IF (ASSOCIATED(gdat)) THEN
          DEALLOCATE(gdat)
       END IF
       IF (ASSOCIATED(gps)) THEN
          DEALLOCATE(gps)
       END IF
    END IF

    NULLIFY(zlptr)

  END SUBROUTINE RGTOOL_BI_READ_NCVAR_5D
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE RGTOOL_BI_READ_NCFILE_5D(modstr, fname, t, dat      &
       , iipol, rg_range, ogrid_ID, igrid_id     &
       , vars                          &
       , SDID                          &  ! OPTIONAL
       , lrg, lrgx, lrgy, lrgz, lok    &  ! OPTIONAL
       , hyam, hybm, p0, ps            &  ! OPTIONAL
       , hyai, hybi                    &  ! OPTIONAL
       , latm, lonm, lati, loni        &  ! OPTIONAL
       , oarea                         &  ! OPTIONAL
       , att                           &  ! OPTIONAL
       , ldatreq                       &  ! OPTIONAL
       , nlon, nlat, nlev, npar        &  ! OPTIONAL
       , mnc)                             ! OPTIONAL

    ! PERFORMS ONE INTERPOLATION STEP FOR ONE FILE
    ! AND RETURNS m DATA FIELDS AS 5D ARRAY (x,z,m,y,n),
    ! AND OPTIONALLY GRID-STRUCTURE AS ARRAYS
    !
    ! Author: Patrick Joeckel, DLR, May 2012                 (NCREGRID)
    !         Astrid Kerkweg, Uni Mainz, 2012/2013 (expanded for SCRIP)

    USE messy_main_mpi_bi,          ONLY: p_pe, p_nprocs, p_all_comm
    USE messy_main_tools,           ONLY: FIND_NEXT_FREE_UNIT
    USE messy_main_grid_netcdf,     ONLY: t_ncvar, INIT_NCVAR,GRD_MAXSTRLEN &
                                        , RGMLE, RGMSG, t_multinc
    USE messy_main_grid,            ONLY: t_geohybgrid &
                                        , GRID_ERROR   &
                                        , INIT_GEOHYBGRID
    USE messy_main_grid_tools,      ONLY: RGTOOL_CONVERT, RGTOOL_G2C
    USE messy_main_import_grid,     ONLY: RGTOOL_READ_NCFILE
    USE messy_main_import_grid_par, ONLY: INIT_PARALLEL
    USE messy_main_channel_attributes, ONLY: t_att_ptr
    USE messy_main_channel_repr,       ONLY: repr_def_axes
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma
    USE messy_main_mpi_bi,          ONLY: p_parallel_io, p_io, P_BCAST   &
                                        , SCATTER_GP, dcg
    USE messy_main_grid_bi,         ONLY: P_BCAST_GRID
    USE messy_main_channel_bi,      ONLY: P_BCAST_ATTRIBUTE
    USE messy_main_grid,            ONLY: NEW_GEOHYBGRID
#endif

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, PRESENT, SIZE, TRIM

    ! I/O
    CHARACTER(LEN=*),     INTENT(IN)             :: modstr ! calling module
    CHARACTER(LEN=*),     INTENT(IN)             :: fname ! file name
    INTEGER,              INTENT(IN)             :: t     ! netCDF time step
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER      :: dat   ! (local) data field
    INTEGER,              INTENT(IN)             :: iipol ! interpolation method
    INTEGER                       , POINTER      :: rg_range(:,:)
    INTEGER,              INTENT(INOUT)          :: ogrid_ID ! output grid ID
    INTEGER,              INTENT(INOUT)          :: igrid_ID !  input grid ID
    CHARACTER(LEN=GRD_MAXSTRLEN), DIMENSION(:), &
                                         POINTER :: vars  ! variable names
    INTEGER,              INTENT(INOUT),OPTIONAL :: SDID  ! SCRIP grid ID
    LOGICAL,              INTENT(IN),   OPTIONAL :: lrg   ! regrid really ?
    LOGICAL,              INTENT(IN),   OPTIONAL :: lrgx  ! regrid in x
    LOGICAL,              INTENT(IN),   OPTIONAL :: lrgy  ! regrid in y
    LOGICAL,              INTENT(IN),   OPTIONAL :: lrgz  ! regrid in z
    LOGICAL,              INTENT(OUT),  OPTIONAL :: lok   ! OK?
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
    REAL(dp), DIMENSION(:,:), POINTER,  OPTIONAL :: oarea ! box area output grid
    TYPE(t_att_ptr), DIMENSION(:), POINTER, OPTIONAL :: att ! var attributes
    LOGICAL, INTENT(IN),  OPTIONAL :: ldatreq
    INTEGER, INTENT(OUT), OPTIONAL :: nlon
    INTEGER, INTENT(OUT), OPTIONAL :: nlat
    INTEGER, INTENT(OUT), OPTIONAL :: nlev
    INTEGER, INTENT(OUT), OPTIONAL :: npar
    TYPE(t_multinc), INTENT(INOUT), OPTIONAL :: mnc

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_BI_READ_NCFILE_5D'
    INTEGER                                 :: iou        ! logical I/O unit
    INTEGER                                 :: status     ! memory status
    TYPE(t_ncvar), DIMENSION(:),    POINTER :: var        ! nc-variable list
    REAL(dp), DIMENSION(:,:,:,:),   POINTER :: hdat       ! temporary data field
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: gdat       ! (global) data field
    ! (global) surface pressure
    REAL(dp), DIMENSION(:,:),       POINTER :: gps
    ! (global) data pointer
    REAL(dp), DIMENSION(:,:,:,:),   POINTER :: zgptr => NULL()
    ! (local) data pointer
    REAL(dp), DIMENSION(:,:,:,:),   POINTER :: zlptr => NULL()
    TYPE(t_geohybgrid)                      :: convgrid   ! grid-structure
    LOGICAL                                 :: llrg, llrgx, llrgy, llrgz, llok
    INTEGER                                 :: i          ! counter
    INTEGER                                 :: nvar       ! number of variables
    INTEGER                                 :: klat       ! number of latitudes
    INTEGER                                 :: klon       ! number of longitudes
    INTEGER                                 :: klev       ! number of levels
    INTEGER                                 :: kparam     ! number of parameters
    INTEGER                                 :: kparam2    ! number of parameters
    CHARACTER(LEN=4)                        :: dat_order = def_order
    INTEGER                                 :: idx_lon, idx_lat, idx_lev
    INTEGER                                 :: idx_param
    INTEGER                                 :: ix
    INTEGER                                 :: SCRIP_ID
    INTEGER , SAVE                          :: callnum = 0
    LOGICAL                                 :: zldatreq
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    INTEGER                                 :: s1d        ! size of 1D field
    INTEGER                                 :: natt ! number of attributes
    INTEGER                                 :: ia   ! attribute counter
    LOGICAL                                 :: lsc        ! scatter results ?
#endif

    ! INITIALIZE
    ! I/O
    callnum = callnum + 1

    IF (ASSOCIATED(dat)) THEN
       DEALLOCATE(dat)
    END IF
    NULLIFY(dat)
    !
    IF (ASSOCIATED(rg_range)) THEN
       DEALLOCATE(rg_range)
    END IF
    NULLIFY(rg_range)
    !
    IF (ASSOCIATED(vars)) THEN
       DEALLOCATE(vars)
    END IF
    NULLIFY(vars)
    !
    IF (PRESENT(ps)) THEN
       IF (ASSOCIATED(ps)) THEN
          DEALLOCATE(ps)
       END IF
       NULLIFY(ps)
    END IF
    !
    IF (PRESENT(hyam)) THEN
       IF (ASSOCIATED(hyam)) THEN
          DEALLOCATE(hyam)
       END IF
       NULLIFY(hyam)
    END IF
    !
    IF (PRESENT(hybm)) THEN
       IF (ASSOCIATED(hybm)) THEN
          DEALLOCATE(hybm)
       END IF
       NULLIFY(hybm)
    END IF
    !
    IF (PRESENT(hyai)) THEN
       IF (ASSOCIATED(hyai)) THEN
          DEALLOCATE(hyai)
       END IF
       NULLIFY(hyai)
    END IF
    !
    IF (PRESENT(hybi)) THEN
       IF (ASSOCIATED(hybi)) THEN
          DEALLOCATE(hybi)
       END IF
       NULLIFY(hybi)
    END IF
    !
    IF (PRESENT(p0)) THEN
       IF (ASSOCIATED(p0)) THEN
          DEALLOCATE(p0)
       END IF
       NULLIFY(p0)
    END IF
    !
    IF (PRESENT(latm)) THEN
       IF (ASSOCIATED(latm)) THEN
          DEALLOCATE(latm)
       END IF
       NULLIFY(latm)
    END IF
    !
    IF (PRESENT(lonm)) THEN
       IF (ASSOCIATED(lonm)) THEN
          DEALLOCATE(lonm)
       END IF
       NULLIFY(lonm)
    END IF
    !
    IF (PRESENT(lati)) THEN
       IF (ASSOCIATED(lati)) THEN
          DEALLOCATE(lati)
       END IF
       NULLIFY(lati)
    END IF
    !
    IF (PRESENT(loni)) THEN
       IF (ASSOCIATED(loni)) THEN
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
    IF (PRESENT(SDID)) THEN
       SCRIP_ID = SDID
    ELSE
       SCRIP_ID = -99
    ENDIF
    !
    IF (PRESENT(ldatreq)) THEN
       zldatreq = ldatreq
    ELSE
       zldatreq = .TRUE.
    END IF
    !
    NULLIFY(hdat)       ! temporary data field
    NULLIFY(gdat)       ! (global) data field
    NULLIFY(gps)        ! (global) surface pressure field
    NULLIFY(var)        ! variable names

    IF (PRESENT(att)) THEN
       IF (ASSOCIATED(att)) THEN
          DO i=1, SIZE(att)
             IF (ASSOCIATED(att(i)%ptr)) DEALLOCATE(att(i)%ptr)
             NULLIFY(att(i)%ptr)
          END DO
          DEALLOCATE(att)
          NULLIFY(att)
       END IF
    END IF

    CALL INIT_PARALLEL(p_pe, p_nprocs, p_all_comm)

#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
     IF ( p_parallel_io .OR. (iipol == GTRF_NONE)) THEN
#endif
        iou = find_next_free_unit(100,200)
        CALL RGTOOL_READ_NCFILE( status,iou, TRIM(modstr)//'.nml' &
             , fname, t, var &
             , iipol, ogrid_ID, igrid_ID, SCRIP_ID &
             , llrg, llrgx, llrgy, llrgz, llok, oarea &
             , convgrid &
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
             , ldompar = .FALSE., lvarpar = .FALSE. &
#else
             , ldompar = .TRUE. , lvarpar = .FALSE. &
#endif
             , mnc=mnc)

        IF (status /= 0) CALL ERROR_BI(grid_error(status), substr)

        IF (llok) THEN
           IF (.NOT. llrg) THEN
              IF (iipol == GTRF_NONE) THEN
                 ogrid_ID = igrid_ID
              ENDIF
           ENDIF
           nvar = SIZE(var)
           ALLOCATE(vars(nvar))

           ALLOCATE(rg_range(nvar,2))
           rg_range(:,:) = 0
           DO i=1, nvar
              DO ix=1,var(i)%natts
                 IF (TRIM(var(i)%att(ix)%name) == 'mmig_ixf_range') THEN
                    rg_range(i,1:2) = INT(var(i)%att(ix)%dat%vi(1:2))
                 END IF
              END DO
           END DO

           ! 1st STEP
           vars(1) = TRIM(var(1)%name)

           IF (iipol == GTRF_NONE) THEN
              dat_order = 'xyzn'
              idx_lon   = 1
              idx_lat   = 2
              idx_lev   = 3
              idx_param = 4
           ELSE
              dat_order = repr_def_axes(_RI_XYZN_('x','y','z','n'))
              idx_lon   = _IX_XYZN_
              idx_lat   = _IY_XYZN_
              idx_lev   = _IZ_XYZN_
              idx_param = _IN_XYZN_
           END IF

           CALL RGTOOL_CONVERT(var(1), hdat, convgrid, order=dat_order)

           klon   = SIZE(hdat, idx_lon)
           klat   = SIZE(hdat, idx_lat)
           klev   = SIZE(hdat, idx_lev)
           kparam = SIZE(hdat, idx_param)
           IF (PRESENT(nlon)) nlon = klon
           IF (PRESENT(nlat)) nlat = klat
           IF (PRESENT(nlev)) nlev = klev
           IF (PRESENT(npar)) npar = kparam

           ! ALLOCATE SPACE FOR GLOBAL DATA FIELD
           IF (zldatreq) THEN
              ALLOCATE(gdat(SIZE(hdat,1),SIZE(hdat,2),SIZE(hdat,3)&
                   ,SIZE(hdat,4),nvar))
              gdat(:,:,:,:,1) = hdat(:,:,:,:)
           END IF
           ! COPY VARIABLE-DATA TO GLOBAL DATA FIELD
           ! FREE MEMORY OF TEMPORARY-DATA
           DEALLOCATE(hdat)
           NULLIFY(hdat)

           IF (PRESENT(att)) THEN
              ALLOCATE(att(nvar))
              DO i=1, nvar
                 CALL NCVAR2CHOATT(var(i), att(i)%ptr)
              END DO
           END IF
           CALL INIT_NCVAR(var(1))

           ! 2nd ... nvar
           DO i=2, nvar
              vars(i) = TRIM(var(i)%name)

              CALL RGTOOL_CONVERT(var(i), hdat, convgrid, order=dat_order)

              CALL INIT_NCVAR(var(i))
              kparam2 = SIZE(hdat, idx_param)
              IF (kparam2 /= kparam) THEN
                 CALL RGMSG(substr, RGMLE, &
                      'PARAMETER DIMENSION MISMATCH OF VARIABLES !')
              END IF
              IF (zldatreq) THEN
                 ! COPY VARIABLE-DATA TO GLOBAL DATA FIELD
                 gdat(:,:,:,:,i) = hdat(:,:,:,:)
              END IF
              ! FREE MEMORY OF VARIABLE-DATA
              DEALLOCATE(hdat)
              NULLIFY(hdat)
           END DO

           DEALLOCATE(var)
           NULLIFY(var)

           ! GRID
           CALL RGTOOL_G2C(convgrid, hyam, hybm, p0, gps     &
                ,hyai, hybi                    &
                ,latm, lonm                    &
                ,lati, loni, klat, klon, klev, order=dat_order)
        ENDIF ! llok
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
     END IF ! p_parallel_io

     ! BROADCAST SUCCESS
     CALL P_BCAST(llok, p_io)
     IF (.NOT. llrg .AND. iipol /= GTRF_NONE) THEN
        CALL P_BCAST_grid(convgrid, p_io)
        IF (.NOT. p_parallel_io) THEN
           CALL NEW_GEOHYBGRID(status, igrid_ID, convgrid)
        END IF
     END IF

     IF (llok) THEN
        IF (PRESENT(att)) THEN
           CALL P_BCAST(nvar, p_io)
           IF (.NOT. p_parallel_io) THEN
              ALLOCATE(att(nvar))
           ENDIF
           DO i=1, nvar
              IF (p_parallel_io) THEN
                 IF (ASSOCIATED(att(i)%ptr)) THEN
                    natt = SIZE(att(i)%ptr)
                 ELSE
                    natt = 0
                 END IF
              END IF
              CALL P_BCAST(natt, p_io)
              IF (.NOT. p_parallel_io) THEN
                 IF (natt > 0) THEN
                    ALLOCATE(att(i)%ptr(natt))
                 END IF
              END IF
              DO ia=1, natt
                 CALL P_BCAST_ATTRIBUTE(att(i)%ptr(ia), p_io)
              END DO
           END DO
        ENDIF
     END IF

#endif

    CALL INIT_GEOHYBGRID(convgrid) ! memory leak fixed

    IF (PRESENT(lok)) THEN
       lok = llok
    END IF

    SELECT CASE (iipol)
    CASE (GTRF_NONE)
       ! NOTHING TO DO
       IF (llok) THEN
          IF (zldatreq) THEN
             ALLOCATE(dat(SIZE(gdat,1),SIZE(gdat,2),SIZE(gdat,3) &
                  , SIZE(gdat,4),SIZE(gdat,5)))
             ! COPY FIRST ON I/O-PROCESSOR
             dat(:,:,:,:,:) = gdat(:,:,:,:,:)
          END IF
       END IF ! llok
#if defined(ECHAM5) || (defined(CESM1) && !defined(HOMMESE))
    CASE (GTRF_NRGD, GTRF_SCRP)
       ! BROADCST DIMENSIONS
       IF (llok) THEN
          CALL P_BCAST(nvar, p_io)
          CALL P_BCAST(klon, p_io)
          CALL P_BCAST(klev, p_io)
          CALL P_BCAST(klat, p_io)
          CALL P_BCAST(kparam, p_io)
          IF (PRESENT(nlon)) nlon = klon
          IF (PRESENT(nlat)) nlat = klat
          IF (PRESENT(nlev)) nlev = klev
          IF (PRESENT(npar)) npar = kparam
       END IF

       ! BROADCAST VARIABLE NAMES
       IF (llok) THEN
          IF (.NOT.p_parallel_io) THEN
             ALLOCATE(vars(nvar))
          END IF
          DO i=1, nvar
             CALL P_BCAST(vars(i), p_io)
          END DO
       END IF

       IF (zldatreq) THEN

          ! NOTE: SCATTERING THE RESULTS (dat, ps) TO THE PROCESSORS IS ONLY
          !       POSSIBLE
          !       IF REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy.AND.lrgz)
          lsc = (llrg.AND.llok.AND.llrgx.AND.llrgy)
          !
          IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
             ! ALLOCATE SPACE FOR LOCAL DATA FIELD
             ALLOCATE(dat(nproma,klev,kparam,ngpblks,nvar))
             ! SCATTER GLOBAL DATA TO LOCAL FIELDS
             DO i=1, nvar
                IF (p_parallel_io) zgptr => gdat(:,:,:,:,i)
                zlptr => dat(:,:,:,:,i)
                CALL SCATTER_GP(zgptr, zlptr, dcg)
             END DO
          ELSE            ! COPY AND BROADCAST ENTIRE FIELD
             IF (llok) THEN
                ! ALLOCATE SPACE FOR DATA FIELD
                ALLOCATE(dat(klon,klev,kparam,klat,nvar))
                ! COPY FIRST ON I/O-PROCESSOR
                IF (p_parallel_io) THEN
                   dat(:,:,:,:,:) = gdat(:,:,:,:,:)
                END IF
                ! BROADCAST RESULT
                DO i=1, nvar
                   CALL P_BCAST(dat(:,:,:,:,i), p_io)
                END DO
             END IF ! llok
          END IF ! lsc

       END IF

       ! NOTE: SCATTERING THE RESULTS (ps) TO THE PROCESSORS IS ONLY POSSIBLE
       !       IF HORIZONTAL REGRIDDING WAS PERFORMED (lrg.AND.lrgx.AND.lrgy)
       IF (PRESENT(ps)) THEN
          lsc = llrg.AND.llok.AND.llrgx.AND.llrgy
          !
          IF (lsc) THEN   ! SCATTER DATA TO PROCESSORS
             ! ALLOCATE SPACE FOR LOCAL DATA FIELD
             ALLOCATE(ps(nproma,ngpblks))
             ! SCATTER GLOBAL DATA TO LOCAL FIELDS
             CALL SCATTER_GP(gps, ps, dcg)
          ELSE
             IF (llok) THEN
                ! ALLOCATE SPACE FOR DATA FIELD
                ALLOCATE(ps(klon,klat))
                IF (p_parallel_io) THEN
                   ps(:,:) = gps(:,:)
                END IF
                ! BROADCAST RESULT
                DO ix = 1,SIZE(ps,2)
                   CALL P_BCAST(ps(:,ix) , p_io)
                END DO
             END IF ! llok
          END IF ! lsc
       END IF ! ps present

       ! 1D ARRAYS NEVER SCATTERED
       IF (llok) THEN
          IF (.NOT. p_parallel_io) THEN
             ALLOCATE(rg_range(nvar,2))
          END IF
          CALL P_BCAST(rg_range, p_io)
          IF (PRESENT(hyam)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(hyam)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(hyam(s1d))
             END IF
             CALL P_BCAST(hyam ,p_io)
          END IF
          !
          IF (PRESENT(hybm)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(hybm)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(hybm(s1d))
             END IF
             CALL P_BCAST(hybm ,p_io)
          END IF
          !
          IF (PRESENT(hyai)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(hyai)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(hyai(s1d))
             END IF
             CALL P_BCAST(hyai ,p_io)
          END IF
          !
          IF (PRESENT(hybi)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(hybi)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(hybi(s1d))
             END IF
             CALL P_BCAST(hybi ,p_io)
          END IF
          !
          IF (PRESENT(p0)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(p0)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(p0(s1d))
             END IF
             CALL P_BCAST(p0 ,p_io)
          END IF
          !
          IF (PRESENT(latm)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(latm)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(latm(s1d))
             END IF
             CALL P_BCAST(latm ,p_io)
          END IF
          !
          IF (PRESENT(lonm)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(lonm)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(lonm(s1d))
             END IF
             CALL P_BCAST(lonm ,p_io)
          END IF
          !
          IF (PRESENT(lati)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(lati)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(lati(s1d))
             END IF
             CALL P_BCAST(lati ,p_io)
          END IF
          !
          IF (PRESENT(loni)) THEN
             IF (p_parallel_io) THEN
                s1d = SIZE(loni)
             END IF
             CALL P_BCAST(s1d ,p_io)
             IF (.NOT.p_parallel_io) THEN
                ALLOCATE(loni(s1d))
             END IF
             CALL P_BCAST(loni ,p_io)
          END IF
       END IF
#else
    CASE (GTRF_SCRP)
       ! ALLOCATE SPACE FOR DATA FIELD
       IF (zldatreq) THEN
          ALLOCATE(dat(SIZE(gdat,1),SIZE(gdat,2) &
               ,SIZE(gdat,3),SIZE(gdat,4),SIZE(gdat,5)))

          dat(:,:,:,:,:) = gdat(:,:,:,:,:)
       END IF

       IF (PRESENT(ps)) THEN
          ALLOCATE(ps(klon,klat))
          ps(:,:) = gps(:,:)
       END IF
#endif

    END SELECT

    IF (PRESENT(SDID)) THEN
       SDID = SCRIP_ID
    ENDIF

    ! CLEAN UP: FREE MEMORY OF GLOBAL FIELDS ON I/O PROCESSOR
    IF (llok) THEN
       IF (ASSOCIATED(gdat)) THEN
          DEALLOCATE(gdat)
       END IF
       IF (ASSOCIATED(gps)) THEN
          DEALLOCATE(gps)
       END IF
    END IF

    NULLIFY(zlptr)
    NULLIFY(zgptr)

  END SUBROUTINE RGTOOL_BI_READ_NCFILE_5D
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE RGTINIT_INIT_NML(rgtinit, modstr)

    ! INITIALIZES RGTINITS FROM NAMELIST modstr.nml
    ! NAMELIST:
    !   RG_INIT(1) = 'action1',
    !   RG_INIT(2) = 'action2',
    !   RG_INIT(3) = 'action3',
    !   ...
    !                |________|
    !                    V
    !                  ACTION
    !
    ! ACTION= variable-name or filename for regridding
    !
    ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

    USE messy_main_channel_bi,           ONLY: n_dom
    USE messy_main_mpi_bi,               ONLY: p_parallel_io, p_io, P_BCAST
    USE messy_main_blather_bi,           ONLY: START_MESSAGE_BI, END_MESSAGE_BI
    USE messy_main_tools,                ONLY: find_next_free_unit &
                                             , domains_from_string &
                                             , MISSING
    USE messy_main_grid_netcdf,          ONLY: RGMSG, RGMLIC

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    TYPE (t_rgtinit), DIMENSION(:), POINTER :: rgtinit
    CHARACTER(LEN=*), INTENT(IN)            :: modstr  ! submodel-string

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTINIT_INIT_NML'
    INTEGER             :: iou                  ! logical I/O unit
    TYPE(t_rgtinit_io)  :: RG_INIT(NMAXRGTE)    ! I/O RGT-events
    INTEGER             :: n                    ! number of EVENT-COUNTERS
    INTEGER             :: i                    ! counter
    CHARACTER(LEN=32)   :: varname = '  '
    INTEGER             :: nd, num, status
    INTEGER, DIMENSION(:), POINTER :: domnum  => NULL()

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL RGTINIT_READ_NML_CPL(status, iou)
       IF (status /= 0 .AND. status /= MISSING) &
          CALL ERROR_BI('rgtinit_read_nml_cpl reported an error',substr)
    END IF

    CALL START_MESSAGE_BI(modstr,'RGTINIT INITIALISATION',substr)

    ! GET NUMBER OF EVENTS
    IF (p_parallel_io) THEN
       n = 0
       DO i=1, NMAXRGTE
          IF (TRIM(RG_INIT(i)%act) == '') CYCLE
           !
          CALL domains_from_string(status,TRIM(RG_INIT(i)%name),n_dom, num)
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          n = n+num
       END DO
       CALL RGMSG(substr, RGMLIC, &
            ' ',N,' RGTINIT RECOGNIZED IN '''//TRIM(modstr)//'.nml''')
    END IF
    ! BROADCAST RESULT
    CALL P_BCAST(n, p_io)

    ! ALLOCATE SPACE
    ALLOCATE(RGTINIT(n))
    ! TRANSFER INFORMATION
    IF (p_parallel_io) THEN
       n = 0
       rgtrig_loop: DO i=1, NMAXRGTE
          IF (TRIM(RG_INIT(i)%act) == '') CYCLE
          CALL domains_from_string(status,TRIM(RG_INIT(i)%name),n_dom, num &
                ,varname,dnums=domnum )
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          domain_loop: DO nd = 1, SIZE(domnum)
             n = n+1
             RGTINIT(n)%act = RG_INIT(i)%act
             IF (TRIM(varname) == '') THEN
                WRITE(varname(1:7),'(a3,i4.4)') 'RGT',i
             END IF
             RGTINIT(n)%name        = TRIM(ADJUSTL(varname))
             RGTINIT(n)%domain_idx  = domnum(nd)
             RGTINIT(n)%tstep        = RG_INIT(i)%tstep
          END DO domain_loop
          DEALLOCATE(domnum);  NULLIFY(domnum)
       END DO rgtrig_loop
    END IF
    ! BROADCAST RESULTS
    DO i=1, n
       CALL P_BCAST(rgtinit(i)%name,       p_io)
       CALL P_BCAST(rgtinit(i)%act,        p_io)
       CALL P_BCAST(rgtinit(i)%domain_idx, p_io)
       CALL P_BCAST(rgtinit(i)%tstep,      p_io)
    END DO

    status = 0

    CALL END_MESSAGE_BI(modstr,'RGTINIT INITIALISATION',substr)

  CONTAINS

    ! ------------------------------------------------------------------
    SUBROUTINE RGTINIT_READ_NML_CPL(status, iou)

      ! read namelist with RGT-events
      !
      ! Author: Patrick Joeckel, MPICH, Mar 2004

      USE messy_main_tools, ONLY: READ_NML_OPEN, READ_NML_CHECK &
                                , READ_NML_CLOSE, MISSING

      IMPLICIT NONE

      ! I/O
      INTEGER, INTENT(OUT) :: status     ! error status
      INTEGER, INTENT(IN)  :: iou        ! I/O unit

      ! (LOCAL) NAMELIST VARIABLES
      CHARACTER(LEN=*), PARAMETER :: substr = 'rgtinit_read_nml_cpl'
      LOGICAL                     :: lex      ! file exists ?
      INTEGER                     :: fstat    ! file status

      NAMELIST /RGTINIT/ RG_INIT

      status = 10

      CALL READ_NML_OPEN(lex, substr, iou, 'RGTINIT', modstr, status=status)
      IF (.not.lex) RETURN    ! <modstr>.nml does not exist

      IF (status /= MISSING) THEN
         status = 10
         READ(iou, NML=RGTINIT, IOSTAT=fstat)
         CALL READ_NML_CHECK(fstat, substr, iou, 'RGTINIT', modstr)
         IF (fstat /= 0) RETURN  ! error while reading namelist

         CALL READ_NML_CLOSE(substr, iou, modstr)
         status = 0  ! no ERROR
      ELSE
         write (*,*) ' WARNING: no &RGTINIT namelist in import.nml'
         CALL READ_NML_CLOSE(substr, iou, modstr)
      END IF

    END SUBROUTINE RGTINIT_READ_NML_CPL
    ! ------------------------------------------------------------------

  END SUBROUTINE RGTINIT_INIT_NML
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE RGTEVENT_INIT_NML(rgt, modstr)

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

    USE messy_main_channel_bi,  ONLY: n_dom
    USE messy_main_mpi_bi,      ONLY: p_parallel_io, p_io, P_BCAST
    USE messy_main_timer_bi,    ONLY: P_BCAST_EVENT
    USE messy_main_blather_bi,  ONLY: START_MESSAGE_BI, END_MESSAGE_BI
    USE messy_main_tools,       ONLY: find_next_free_unit &
                                    , domains_from_string
    USE messy_main_grid_netcdf, ONLY: RGMSG,  RGMLIC

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    TYPE (t_rgtevent), DIMENSION(:), POINTER :: rgt
    CHARACTER(LEN=*),  INTENT(IN)            :: modstr  ! submodel-string

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTEVENT_INIT_NML'
    INTEGER             :: iou                  ! logical I/O unit
    TYPE(t_rgtevent_io) :: RG_TRIG(NMAXRGTE)    ! I/O RGT-events
    INTEGER             :: n                    ! number of EVENT-COUNTERS
    INTEGER             :: i                    ! counter
    LOGICAL             :: lstat                ! event status
    INTEGER             :: cpos                 ! current counter position
    CHARACTER(LEN=32)   :: varname = '  '
    INTEGER             :: nd, num, status
    INTEGER, DIMENSION(:), POINTER :: domnum  => NULL()

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL RGTEVENT_READ_NML_CPL(status, iou)
       IF (status /= 0) &
          CALL ERROR_BI('rgtevent_read_nml_cpl reported an error',substr)
    END IF

    CALL START_MESSAGE_BI(modstr,'RGTEVENT INITIALISATION',substr)

    ! GET NUMBER OF EVENTS
    IF (p_parallel_io) THEN
       n = 0
       DO i=1, NMAXRGTE
          IF (RG_TRIG(i)%cnt%start == -1) CYCLE
           !
          CALL domains_from_string(status,TRIM(RG_TRIG(i)%cnt%name),n_dom, num)
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          n = n+num
       END DO
       CALL RGMSG(substr, RGMLIC, &
            ' ',N,' RGTEVENTS RECOGNIZED IN '''//TRIM(modstr)//'.nml''')
    END IF
    ! BROADCAST RESULT
    CALL P_BCAST(n, p_io)

    ! ALLOCATE SPACE
    ALLOCATE(RGT(n))
    ! TRANSFER INFORMATION
    IF (p_parallel_io) THEN
       n = 0
       rgtrig_loop: DO i=1, NMAXRGTE
          IF (RG_TRIG(i)%cnt%start == -1) CYCLE
          CALL domains_from_string(status,TRIM(RG_TRIG(i)%cnt%name),n_dom, num &
                ,varname,dnums=domnum )
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)
          domain_loop: DO nd = 1, SIZE(domnum)
             n = n+1
             RGT(n)%io = RG_TRIG(i)
             IF (TRIM(varname) == '') THEN
                WRITE(varname(1:7),'(a3,i4.4)') 'RGT',i
             END IF
             rgt(n)%io%cnt%name = TRIM(ADJUSTL(varname))
             rgt(n)%domain_idx  = domnum(nd)
          END DO domain_loop
          DEALLOCATE(domnum);  NULLIFY(domnum)
       END DO rgtrig_loop
    END IF
    ! BROADCAST RESULTS
    DO i=1, n
       CALL P_BCAST(rgt(i)%io%cnt%start,   p_io)
       CALL P_BCAST(rgt(i)%io%cnt%step,    p_io)
       CALL P_BCAST(rgt(i)%io%cnt%reset,   p_io)
       CALL P_BCAST(rgt(i)%io%cnt%current, p_io)
       CALL P_BCAST(rgt(i)%io%cnt%name,    p_io)
       CALL P_BCAST(rgt(i)%io%act,         p_io)
       CALL P_BCAST(rgt(i)%domain_idx,     p_io)
       CALL P_BCAST_event(rgt(i)%io%evt,   p_io)

       ! INITIALIZE EVENTS
       CALL RGTEVENT_INIT(rgt(i))

       ! UPDATE TO/FROM EVENT COUNTER LIST
       CALL RGTEVENT_STATUS(lstat, cpos, rgt, i, linit=.TRUE.)
    END DO

    status = 0

    CALL END_MESSAGE_BI(modstr,'RGTEVENT INITIALISATION',substr)

  CONTAINS

    ! ------------------------------------------------------------------
    SUBROUTINE RGTEVENT_READ_NML_CPL(status, iou)

      ! read namelist with RGT-events
      !
      ! Author: Patrick Joeckel, MPICH, Mar 2004

      USE messy_main_tools, ONLY: READ_NML_OPEN, READ_NML_CHECK &
                                , READ_NML_CLOSE

      IMPLICIT NONE

      ! I/O
      INTEGER, INTENT(OUT) :: status     ! error status
      INTEGER, INTENT(IN)  :: iou        ! I/O unit

      ! (LOCAL) NAMELIST VARIABLES
      CHARACTER(LEN=*), PARAMETER :: substr = 'rgtevent_read_nml_cpl'
      LOGICAL                     :: lex      ! file exists ?
      INTEGER                     :: fstat    ! file status

      NAMELIST /RGTEVENTS/ RG_TRIG

      status = 1

      CALL READ_NML_OPEN(lex, substr, iou, 'RGTEVENTS', modstr)
      IF (.not.lex) RETURN    ! <modstr>.nml does not exist

      READ(iou, NML=RGTEVENTS, IOSTAT=fstat)
      CALL READ_NML_CHECK(fstat, substr, iou, 'RGTEVENTS', modstr)
      IF (fstat /= 0) RETURN  ! error while reading namelist

      CALL READ_NML_CLOSE(substr, iou, modstr)

      status = 0  ! no ERROR

    END SUBROUTINE RGTEVENT_READ_NML_CPL
    ! ------------------------------------------------------------------

  END SUBROUTINE RGTEVENT_INIT_NML
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE RGTEVENT_STATUS(status, cpos, rgt, ix, action &
       , lstop, linit, lftrig)

    ! RETURNS THE EVENT STATUS (flag) OF A NAMED (name)
    ! OR INDEXED (index) RGT-EVENT IN A LIST (rgt)
    !
    ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

    USE messy_main_mpi_bi,      ONLY: p_parallel_io
    USE messy_main_blather_bi,  ONLY: ERROR_BI
#if !defined(ICON) && !defined(COSMO)
    USE messy_main_timer,       ONLY: lresume
#endif
    USE messy_main_timer,       ONLY: lstart
    USE messy_main_import_grid, ONLY: RGTOOL_NCRGCNT_RST
    USE messy_main_grid_netcdf, ONLY: RGMSG, RGMLE, RGMLW

    IMPLICIT NONE

    INTRINSIC :: PRESENT, SIZE, TRIM

    ! I/O
    LOGICAL,                     INTENT(OUT)        :: status ! event status
    INTEGER,                     INTENT(OUT)        :: cpos   ! counter pos.
    TYPE(t_rgtevent), DIMENSION(:), POINTER         :: rgt    ! RGT-event list
    INTEGER,                     INTENT(IN)         :: ix     ! index of event
    CHARACTER(LEN=RGTMAXACTSTR), INTENT(OUT), OPTIONAL :: action ! action string
    LOGICAL,                     INTENT(IN),  OPTIONAL :: lstop  ! stop on error
    LOGICAL,                     INTENT(IN),  OPTIONAL :: linit  ! initialize?
    LOGICAL,                     INTENT(IN),  OPTIONAL :: lftrig ! force event

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'RGTEVENT_STATUS'
    LOGICAL                     :: llstop ! stop on error
    INTEGER                     :: RGMLX  ! RGML(W,E)
    LOGICAL                     :: zlinit
#if defined(ICON) || defined(COSMO)
    LOGICAL                     :: zlftrig
#endif

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

#if defined(ICON) || defined(COSMO)
    IF (PRESENT(lftrig)) THEN
       zlftrig = lftrig
    ELSE
       zlftrig = .FALSE. ! DEFAULT
    END IF
#endif

    IF (ix > SIZE(rgt)) THEN
       IF (p_parallel_io) THEN
          CALL RGMSG(substr, RGMLX, 'INDEX OUT OF RANGE !', .false.)
       END IF
       IF (llstop) CALL ERROR_BI(' ',substr)
       RETURN
    END IF

    ! GET EVENT STATUS
    IF (.NOT. zlinit) CALL RGTEVENT_STAT(rgt(ix), status)

    ! UPDATE TO/FROM COUNTER LIST
#if defined(ICON) || defined(COSMO)
    CALL RGTOOL_NCRGCNT_RST(lstart, zlftrig, status, &
         rgt(ix)%io%cnt, rgt(ix)%n, p_parallel_io, linit)
    status = status .OR. lstart .OR. zlftrig
#else
    CALL RGTOOL_NCRGCNT_RST(lstart, lresume, status, &
         rgt(ix)%io%cnt, rgt(ix)%n, p_parallel_io, linit)
    status = status .OR. lstart .OR. lresume
#endif

    ! RETURN VALUES
    cpos = rgt(ix)%io%cnt%current

    IF (PRESENT(action)) THEN
       action = TRIM(rgt(ix)%io%act)
    END IF

  END SUBROUTINE RGTEVENT_STATUS
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE RGTEVENT_TO_CHANNEL(cha, rgt)

    USE messy_main_channel,            ONLY: new_channel_object, new_attribute &
                                           , AF_RST_CMP
    USE messy_main_channel_error_bi,   ONLY: CHANNEL_HALT
    USE messy_main_channel_bi,         ONLY: SCALAR
    USE messy_main_channel_mem,        ONLY: dom_curid

    IMPLICIT NONE
    INTRINSIC :: SIZE, REAL

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: cha
    TYPE(t_rgtevent), DIMENSION(:), POINTER :: rgt ! (IN) ! RGT-event list

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'rgtevent_to_channel'
    CHARACTER(LEN=*), PARAMETER :: pre = 'RGT_'
    INTEGER                     :: status
    INTEGER :: i

    DO i=1, SIZE(rgt)

       IF (rgt(i)%domain_idx /= dom_curid) CYCLE

       ! prepend RGT_ to enable output suppression with wildcards
       CALL new_channel_object(status, TRIM(cha)           &
            , pre//TRIM(rgt(i)%io%cnt%name), p0 = rgt(i)%n &
            , reprid=SCALAR, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       !
       rgt(i)%n = REAL(rgt(i)%io%cnt%current, DP)
       !
       CALL new_attribute(status, TRIM(cha)            &
            , pre//TRIM(rgt(i)%io%cnt%name)            &
            , 'event_counter', i=rgt(i)%io%evt%counter &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, TRIM(cha)      &
            , pre//TRIM(rgt(i)%io%cnt%name)      &
            , 'event_unit', c=rgt(i)%io%evt%unit &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, TRIM(cha)                  &
            , pre//TRIM(rgt(i)%io%cnt%name)                  &
            , 'event_adjustment', c=rgt(i)%io%evt%adjustment &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, TRIM(cha)          &
            , pre//TRIM(rgt(i)%io%cnt%name)          &
            , 'event_offset', r=rgt(i)%io%evt%offset &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, TRIM(cha)         &
            , pre//TRIM(rgt(i)%io%cnt%name)         &
            , 'stepper_min', i=rgt(i)%io%cnt%start  &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, TRIM(cha)         &
            , pre//TRIM(rgt(i)%io%cnt%name)         &
            , 'stepper_step', i=rgt(i)%io%cnt%step  &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, TRIM(cha)         &
            , pre//TRIM(rgt(i)%io%cnt%name)         &
            , 'stepper_max', i=rgt(i)%io%cnt%reset  &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, TRIM(cha)             &
            , pre//TRIM(rgt(i)%io%cnt%name)             &
            , 'stepper_start', i=rgt(i)%io%cnt%current  &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, TRIM(cha)     &
            , pre//TRIM(rgt(i)%io%cnt%name)     &
            , 'stepper_action', c=rgt(i)%io%act &
            , iflag = AF_RST_CMP)
       CALL channel_halt(substr, status)
       !
    END DO

  END SUBROUTINE RGTEVENT_TO_CHANNEL
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE RGTEVENT_FROM_CHANNEL(rgt)

    USE messy_main_channel_mem, ONLY: dom_curid

    IMPLICIT NONE
    INTRINSIC :: SIZE, NINT

    ! I/O
    TYPE(t_rgtevent), DIMENSION(:), POINTER :: rgt ! (IN) ! RGT-event list

    ! LOCAL
    INTEGER :: i

    DO i=1, SIZE(rgt)
       IF (rgt(i)%domain_idx /= dom_curid) CYCLE
       rgt(i)%io%cnt%current = NINT(rgt(i)%n)
    END DO

  END SUBROUTINE RGTEVENT_FROM_CHANNEL
  ! ------------------------------------------------------------------

  ! ####### PRIVATE ROUTINES #########################################

  ! ------------------------------------------------------------------
  SUBROUTINE RGTEVENT_INIT(rgt)

    ! INITIALIZES ONE INTERNAL RGT-EVENT (rgt) FROM
    ! I/O NAMELIST INFORMATION
    !
    ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

    USE MESSY_MAIN_TIMER_BI,          ONLY: TIMER_EVENT_INIT

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    TYPE (t_rgtevent), INTENT(INOUT)          :: rgt

    CALL TIMER_EVENT_INIT(rgt%event, rgt%io%evt, &
         TRIM(rgt%io%cnt%name), 'next')

  END SUBROUTINE RGTEVENT_INIT
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE RGTEVENT_STAT(rgt, flag)

    ! RETURNS EVENT STATUS (flag) OF ONE RGT-EVENT
    !
    ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

    USE messy_main_timer_bi,   ONLY: EVENT_STATE
    USE messy_main_timer,      ONLY: next_date

    IMPLICIT NONE

    ! I/O
    TYPE (t_rgtevent), INTENT(INOUT) :: rgt
    LOGICAL,           INTENT(OUT)   :: flag

    flag = event_state(rgt%event, next_date)

  END SUBROUTINE RGTEVENT_STAT
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE NCVAR2CHOATT(ncvar, att)

    USE messy_main_grid_netcdf,        ONLY: t_ncvar, EXTRACT_NCATT
    USE messy_main_channel_attributes, ONLY: t_attribute, AF_NONE &
                                           , TYPE_INTEGER, TYPE_STRING &
                                           , TYPE_REAL_DP

    IMPLICIT NONE

    ! I/O
    TYPE(t_ncvar),                   INTENT(IN)  :: ncvar   ! netCDF variable
    ! INTENT(OUT)
    TYPE(t_attribute), DIMENSION(:), POINTER     :: att ! attribute list

    ! LOCAL
    INTEGER :: n, i
    INTEGER :: type

    IF (ASSOCIATED(att)) THEN
       DEALLOCATE(att); NULLIFY(att)
    ENDIF

    n = ncvar%natts

    ! NO ATTRIBUTES -> NOTHING ELSE TO DO
    IF (n == 0) RETURN

    ALLOCATE(att(n))

    DO i=1, n

       CALL EXTRACT_NCATT(ncvar%att(i) &
            , type        &
            , att(i)%name &
            , att(i)%i    &
            , att(i)%c    &
            , att(i)%r)

       att(i)%iflag = AF_NONE

       SELECT CASE(type)
       CASE(1)
          att(i)%type = TYPE_INTEGER
       CASE(2)
          att(i)%type = TYPE_STRING
       CASE(3)
          att(i)%type = TYPE_REAL_DP
       END SELECT

    END DO

  END SUBROUTINE NCVAR2CHOATT
  ! ------------------------------------------------------------------

! ******************************************************************
END MODULE messy_main_import_grid_tools_bi
! ******************************************************************
