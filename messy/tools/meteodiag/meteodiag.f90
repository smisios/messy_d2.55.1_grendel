! ******************************************************************
! ------------------------------------------------------------------
PROGRAM METEODIAG
  ! ------------------------------------------------------------------
  ! Author: Andreas Baumgaertner, MPICH, Mainz, July 2010
  !         Heini Wernli, ETH Zuerich, July 2010
  ! based on NCREGRID by Patrick Joeckel, MPICH, Mainz, June 2002
  ! ------------------------------------------------------------------

    !- öffnen der relevanten files für einen Monat
    !- loop über alle Zeitschritte auf einem file
    !    - einlesen SLP als 2d Feld: slp(i,j)
    !    - einlesen U und V als 3d Felder
    !    - interpolation von U und V auf Druckflächen 100, 150, ..., 500 hPa
    !    - Berechnung der vertikal gemittelten 2d Windgeschwindigkeit auf
    !      diesen Schichten: vel(i,j)
    !
    !    Hier werden dann unsere Programmteile kommen
    !    slp(i,j) --> cyclone identification --> cyc(i,j)   0/1 field of cyclones
    !    vel(i,j) --> jet identification --> jet(i,j)               0/1 field
    !            of jets
    !
    !    Dann wieder von dir
    !    - output der 2d-Felder cyc(i,j) und jet(i,j) im von dir gewünschten
    !      fileformat
    !- endloop


  USE mo_f2kcli                    ! command line interface
  USE messy_ncregrid_tools_meteodiag
  USE messy_ncregrid_geohyb,    ONLY: geohybgrid, COPY_GEOHYBGRID
  USE messy_ncregrid_netcdf,    ONLY: INIT_NCVAR
  USE messy_main_constants_mem, ONLY: dp, pi, STRLEN_LONG

  IMPLICIT NONE


  INTRINSIC :: TRIM, COS

  CHARACTER(LEN=*), PARAMETER  :: modstr = 'meteodiag'
 
  ! FOR COMMAND LINE
  CHARACTER(LEN=256) :: EXE          ! program name
  CHARACTER(LEN=80)  :: CMD          ! argument
  INTEGER            :: NARG         ! number of arguments

  LOGICAL                               :: lok
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: dat ! data field
  REAL(dp), DIMENSION(:,:,:),   POINTER :: wind_u, wind_v  ! zonal and meridional wind   
  TYPE(geohybgrid)                      :: grid_in, grid_out ! grid-structure
  REAL(dp), DIMENSION(:),   POINTER     :: hyam, hybm, p0    &
                                         , hyai, hybi, latm, lonm &
                                         , latfactor

  REAL(dp), DIMENSION(:,:),   POINTER   :: ps

  ! NAMELIST
  CHARACTER(STRLEN_LONG) :: namelist_U, varname_U, namelist_V, varname_V, outfile
  LOGICAL :: coslat = .TRUE.   

  REAL(dp), DIMENSION(:,:),   POINTER   :: cyc,jet
  INTEGER :: tstep
  INTEGER :: klon, klat, klev
  INTEGER :: lon, lat, lev

  INTEGER :: status

  NARG = COMMAND_ARGUMENT_COUNT()    ! number of arguments
  CALL GET_COMMAND_ARGUMENT(0,EXE)   ! program name

  IF (NARG > 1) THEN 
     WRITE(*,*) 'Too many arguments !'
     CALL USAGE(TRIM(EXE)) 
     STOP
  END IF

  IF (NARG == 0) THEN 
     CALL USAGE(TRIM(EXE)) 
     STOP
  END IF

  CALL GET_COMMAND_ARGUMENT(1,CMD)  

  CALL meteodiag_read_nml_ctrl(status, CMD)
  IF (status /= 0) THEN
     WRITE(*,*) 'NAMELIST READ ERROR'
     STOP
  END IF
 

  DO tstep=1,2 ! for testing
 ! DO ! endless DO loop, must be terminated with exit

     !***********************************************************
     ! LOAD
     !***********************************************************
     CALL RGTOOL_METEODIAG_READ_NCVAR(TRIM(namelist_U),TRIM(varname_U) &
          , tstep, dat, lrg = .true.   &
          , lok = lok,  grid=grid_in &
          , hyam=hyam, hybm=hybm, p0=p0, ps=ps        &
          , hyai=hyai, hybi=hybi, latm=latm, lonm=lonm)
     IF (.NOT.lok) THEN
        EXIT
     END IF
     klon = SIZE(lonm)
     klat = SIZE(latm)
     klev = SIZE(hyam)
     ALLOCATE(wind_u(klon,klev,klat) &
             ,wind_v(klon,klev,klat) & 
             ,cyc(klon, klat)        &
             ,jet(klon, klat))

     ALLOCATE(latfactor(klat))
     IF (coslat) THEN 
        latfactor(:) = COS(latm(:)/180.*pi)
     ELSE
        latfactor(:) = 1.
     END IF

     DO lev=1,klev
        DO lon=1,klon
           wind_u(lon,lev,:) = dat(lon,lev,1,:)/latfactor(:)
        END DO
     END DO

     CALL RGTOOL_METEODIAG_READ_NCVAR(TRIM(namelist_V),TRIM(varname_V), &
          tstep, dat, lrg = .true. , lok = lok)
      IF (.NOT.lok) THEN
        EXIT
     END IF

     DO lev=1,klev
        DO lon=1,klon
           wind_v(lon,lev,:) = dat(lon,lev,1,:)/latfactor(:)
        END DO
     END DO


     !***********************************************************
     ! PROCESSING
     !***********************************************************

     ! example:
     cyc(:,:)=0.
     cyc(1,1)=4. ! 88deg N, 1deg E
     cyc(:,:) = wind_u(:,1,:)+wind_v(:,2,:) ! (xzy)
     jet(:,:) = ps(:,:) ! (xy)
     
     !***********************************************************
     ! SAVE
     !***********************************************************
     CALL COPY_GEOHYBGRID(grid_out,grid_in)
     ! "Delete" these variables
     CALL INIT_NCVAR(grid_out%hyam)
     CALL INIT_NCVAR(grid_out%hybm)
     CALL INIT_NCVAR(grid_out%hyai)
     CALL INIT_NCVAR(grid_out%hybi)
     CALL RGTOOL_METEODIAG_SAVE_NCVAR(outfile,'cyc',tstep,cyc,grid_out)
     CALL RGTOOL_METEODIAG_SAVE_NCVAR(outfile,'jet',tstep,jet,grid_out)
     DEALLOCATE(wind_u,wind_v,cyc,jet,latfactor) 
 
  END DO



CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE USAGE(EXE)
    CHARACTER (LEN=*) :: EXE
    WRITE(*,*) '------------------------------------------------------'
    WRITE(*,*) 'METEODIAG Version 0.1'
    WRITE(*,*) 'Author: Andreas Baumgaertner, MPICH, July 2010'
    WRITE(*,*) '        Heini Wernli, ETH Zuerich, July 2010' 
    WRITE(*,*) 'based on NCREGRID by Patrick Joeckel, MPICH, June 2002'
    WRITE(*,*) '-----------------------------------------'
    WRITE(*,*) 'Usage: '//TRIM(EXE)//' <namelist-file>'
    WRITE(*,*) '------------------------------------------------------'
  END SUBROUTINE USAGE
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE meteodiag_read_nml_ctrl(status, nmlfile)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    CHARACTER(LEN=*) , INTENT(IN)  :: nmlfile

    ! mz_pj_20080118: OUT_PREC added
    NAMELIST /CTRL/ namelist_U, varname_U  &
                  , namelist_V, varname_V, coslat, outfile

    INTEGER                     :: iou = 17

    status = 1 ! ERROR

    OPEN(iou,file=TRIM(nmlfile))
    READ(iou, NML=CTRL, IOSTAT=status)
    IF (status /= 0) THEN
       RETURN
    END IF

    CLOSE(iou)

    status = 0  ! no ERROR

  END SUBROUTINE meteodiag_read_nml_ctrl
  ! -------------------------------------------------------------------

END PROGRAM METEODIAG
