! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL MLOCEAN
!
! Author : Markus Kunze, FUB-IFM, January-February  2010
!          
! References:
!
! * P. Joeckel, R. Sander, A. Kerkweg, H. Tost, and J. Lelieveld,
!   Technical Note: The Modular Earth Submodel System (MESSy) - a new
!   approach towards Earth System Modeling,
!   Atmos. Chem. Phys., 5, 433-444, 2005.
!   http://www.atmos-chem-phys.net/5/433 
!
! **********************************************************************

! **********************************************************************
MODULE messy_mlocean
  ! **********************************************************************

  ! ----------- >
  
  USE messy_main_constants_mem, ONLY : dp  &
       , rho_sea  &   ! density of sea water in kg m^-3
       , csw          ! specific heat of sea waterJ/K/kg

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: dp

  ! ----------- <
  !
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'mlocean'  ! submodel name
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.4'      ! submodel version
  !
  ! CTRL-NAMELIST PARAMETERS  public for broadcast in SMIL
  INTEGER,  PUBLIC :: mloswitch ! 1 -> mlo    : use mixed layer ocean from ECHAM
  !                               2 -> plasim : use mixed layer ocean from Planet Sim.
  INTEGER,  PUBLIC :: nflxcorr  ! 0 -> no flux correction
  !                               1 -> use prescribed flux correction from input file
  !                               2 -> input of the net surface energy budget
  REAL(dp), PUBLIC :: flxscale  ! scaling factor for flux correction
  !
  REAL(dp), PUBLIC :: mldmix       ! depth of the mixed layer
  REAL(dp), PUBLIC :: fbase_north  ! basic flux correction for sea ice (north)
  REAL(dp), PUBLIC :: fbase_south  ! basic flux correction for sea ice (south)
  !
  ! PUBLIC SUBROUTINES (to be called from messy_mlocean_e5.f90)
  !
  PUBLIC :: mlocean_read_nml_ctrl
  PUBLIC :: mlocean_mlflux
  PUBLIC :: mlocean_mlfluxp
  PUBLIC :: mlocean_readflux

CONTAINS

  ! =========================================================================
  ! ### add your own public subroutines here
  ! =========================================================================

  ! =========================================================================
  SUBROUTINE mlocean_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools,     ONLY : read_nml_open, read_nml_check, read_nml_close

    USE messy_mlocean_plasim, ONLY : ndiag, nout, ntspd, nocean, nprint, nprhor, &
                                     nperpetual_ocean, taunc, dlayer,            &
                                     vdiffk, newsurf
    
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/         mloswitch
    NAMELIST /CTRL_MLO/     mldmix,  nflxcorr, flxscale, &
                            fbase_north, fbase_south
    NAMELIST /CTRL_PLASIM/  ndiag, nout, ntspd, nocean, nprint, nprhor, &
                            nperpetual_ocean, taunc, mldmix,            &
                            vdiffk, newsurf
    ! LOCAL

    CHARACTER(LEN=*), PARAMETER :: substr='mlocean_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE with default values; can be changed via namelist
    mloswitch   = 1       ! use mixed layer ocean from ECHAM
    nflxcorr    = 2       ! input of the net surface energy budget
    flxscale    = 1._dp   ! scaling factor for the flux correction
    
    mldmix      = 50      ! default depth of mixed layer
    fbase_north = 0.0_dp  ! basic flux correction for sea ice (north)
    fbase_south = 20.0_dp ! basic flux correction for sea ice (south)

    ndiag            = 480 ! diagnostics each ndiag timesteps
    nout             = 32  ! afterburner output each nout timesteps
    nocean           = 1   ! compute ocean yes/no
    newsurf          = 0   ! update surface arrays at restart
    ntspd            = 32  ! ocean timesteps per day
    nperpetual_ocean = 0   ! perpetual climate conditions
    nprint           = 0   ! print debug information
    nprhor           = 0   ! gp to print debug information
    !
    taunc            =  0.   ! newtonian cooling timescale (d)
    vdiffk           = 1.E-4 ! vertikal diffusion coeff. [m**2/s]
    status           = 1       ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    WRITE (*,NML=CTRL)
    WRITE (*,*) TRIM(substr)//':mloswitch:',mloswitch
    !
    SELECT CASE(mloswitch)
    CASE(1)
       READ(iou, NML=CTRL_MLO, IOSTAT=fstat)
       CALL read_nml_check(fstat, substr, iou, 'CTRL_MLO', modstr)
       IF (fstat /= 0) RETURN  ! error while reading namelist
       !
       WRITE (*,NML=CTRL_MLO) 
       WRITE (*,*) TRIM(substr)//':mldmix:  ',mldmix
       WRITE (*,*) TRIM(substr)//':nflxcorr:',nflxcorr
       WRITE (*,*) TRIM(substr)//':flxscale:',flxscale
       WRITE (*,*) TRIM(substr)//':fbase_north:',fbase_north
       WRITE (*,*) TRIM(substr)//':fbase_south:',fbase_south
       !
    CASE(2)
       READ(iou, NML=CTRL_PLASIM, IOSTAT=fstat)
       CALL read_nml_check(fstat, substr, iou, 'CTRL_PLASIM', modstr)
       IF (fstat /= 0) RETURN  ! error while reading namelist
       !
       WRITE (*,'(/," *******************************************")')
       WRITE (*,'(" * '//TRIM(modstr)//' ",a30," *")') TRIM(modver)
       WRITE (*,'(" *******************************************")')
       WRITE (*,'(" * Namelist CTRL_PLASIM from <mlocean.nl>  *")')
       WRITE (*,'(" *******************************************")')
       WRITE (*,NML=CTRL_PLASIM)
       !
    CASE DEFAULT
       WRITE (*,'(/," *******************************************")')
       WRITE (*,'(" * '//TRIM(modstr)//' ",a30," *")') TRIM(modver)
       WRITE (*,'(" *******************************************")')
       WRITE (*,'(" * Namelist CTRL from <mlocean.nl>         *")')
       WRITE (*,'(" * mloswitch ",i4,"is not supported.     *")') mloswitch 
       WRITE (*,'(" *******************************************")')
       RETURN  ! ERROR
    END SELECT
    !
    CALL read_nml_close(substr, iou, modstr)
    !
    dlayer(:) = mldmix   ! layer depth (m)
    status    = 0        ! NO ERROR
    !
    RETURN
  END SUBROUTINE mlocean_read_nml_ctrl
  ! =========================================================================
  !  -----------------------------------------------------------------------
  SUBROUTINE mlocean_mlflux (kproma, nmw1, nmw2, wgt1, wgt2 &
                    , pmldmix, pslf, pamlcorr, palake   &
                    , psst1, psst2, paice1, paice2, paflux1, paflux2)
    !
    ! Description:
    !
    ! Passes flux corrections to mixed layer ocean
    !
    ! Method:
    !
    ! This subroutine uses the climatology of the net ocean surface energy budget,
    ! provided by an external input file.  The Q-flux correction of
    ! the mixed layer ocean is calculated using the climatological sea surface 
    ! temperatures.
    ! The Q-flux correction is then interpolated appropriate.
    !
    ! mlocean_mlflux() is called from messy_mlocean_e5.f90:mlocean_physc().
    !
    ! Authors (ml_flux.f90): 
    !
    ! U. Schlese, DKRZ, January 1993, original source (clsst)
    ! M. Esch, MPI, September 2002, rearranged for mixed layer ocean's
    !                               flux correction
    ! for more details see file AUTHORS
    ! 
    !
    ! Arguments INTENT(in)
    ! --------------------
    INTEGER,                INTENT(in) :: kproma ! number of longitudes on PE
    INTEGER,                INTENT(in) :: nmw1, nmw2
    REAL(dp),               INTENT(in) :: wgt1, wgt2
    REAL(dp),               INTENT(in) :: pmldmix           ! mixed layer depth
    REAL(dp), DIMENSION(:), INTENT(in) :: pslf              ! sea-land mask
    REAL(dp), DIMENSION(:), INTENT(in) :: palake            ! lake mask
    REAL(dp), DIMENSION(:), INTENT(in) :: psst1, psst2      ! sea surface temperature
    REAL(dp), DIMENSION(:), INTENT(in) :: paice1, paice2    ! sea ice
    REAL(dp), DIMENSION(:), INTENT(in) :: paflux1, paflux2  ! net ocean surface energy budget
    !
    ! Arguments INTENT(out)
    ! --------------------
    REAL(dp), DIMENSION(:), INTENT(out) :: pamlcorr ! updated flux correction
    !
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr='mlocean_mlflux'
    INTEGER  :: jl      ! longitude loop index
    REAL(dp) :: zmixcap ! heat capacity of the ocean layer
    REAL(dp) :: zmonlen ! length of month in seconds 
    REAL(dp) :: zic     ! sea ice fraction
    
    !  Executable statements
    !
    ! 1. Set up constants
    !
    zmixcap = rho_sea * csw * pmldmix
    zmonlen = 2592000._dp
    !
    !-- 2. Update mixed layer correction
    !
    DO jl = 1, kproma
       
       IF (palake(jl) == 0._dp .AND. pslf(jl) < 1._dp) THEN
          !
          IF (nmw2 < nmw1) THEN       ! first part of month
             pamlcorr(jl) = wgt1 * paflux1(jl) + wgt2 * paflux2(jl) &
                           -zmixcap * (psst1(jl) - psst2(jl))       &
                           /zmonlen
          ELSE                        !second half of month
             pamlcorr(jl) = wgt1 * paflux1(jl) + wgt2 * paflux2(jl) &
                           -zmixcap * (psst2(jl) - psst1(jl))        &
                           /zmonlen
          END IF
          !
          ! no flux correction on climatological ice points
          !
          zic = (wgt1 * paice1(jl) + wgt2 * paice2(jl)) * 0.01_dp
          zic = MAX(0._dp, MIN(1._dp,zic))
          IF (zic > 0.9_dp) pamlcorr(jl) = 0._dp
          !
       ELSE                          ! land or lake
          pamlcorr(jl) = 0._dp
       END IF
    END DO
    !
    ! apply the scaling factor for the flux correction
    !
    pamlcorr(:) = pamlcorr(:) * flxscale
    
    RETURN
  END SUBROUTINE mlocean_mlflux
  !  -----------------------------------------------------------------------
  !  -----------------------------------------------------------------------
  SUBROUTINE mlocean_mlfluxp (kproma, wgt1, wgt2 &
                    , pslf, palake     &
                    , paflux1, paflux2 &
                    , pamlcorr)
    !
    ! Description:
    !
    ! Passes flux corrections to mixed layer ocean
    !
    ! Method:
    !
    ! This subroutine interpolates the flux corrections .
    !
    ! mlocean_mlflux() is called from messy_mlocean_e5.f90:mlocean_mlflx().
    !
    ! Authors (ml_flux.f90): 
    ! for more details see file AUTHORS
    ! 
    !
    ! Arguments INTENT(in)
    ! --------------------
    INTEGER,                INTENT(in) :: kproma ! number of longitudes on PE
    REAL(dp),               INTENT(in) :: wgt1, wgt2
    REAL(dp), DIMENSION(:), INTENT(in) :: pslf
    REAL(dp), DIMENSION(:), INTENT(in) :: palake
    REAL(dp), DIMENSION(:), INTENT(in) :: paflux1, paflux2 ! flux corr. for 2 months
    !
    ! Arguments INTENT(out)
    ! --------------------
    REAL(dp), DIMENSION(:), INTENT(out) :: pamlcorr ! updated flux correction
    !
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr='mlocean_mlfluxp'
    INTEGER  :: jl      ! longitude loop index
    !
    !  Executable statements
    !
    !-- Update mixed layer correction
    !
    DO jl = 1, kproma
       
       IF (palake(jl) == 0._dp .AND. pslf(jl) < 1._dp) THEN
          !
          ! interpolate in time and apply the scaling factor 
          ! for the flux correction
          !
          pamlcorr(jl) = (wgt1 * paflux1(jl) + wgt2 * paflux2(jl)) &
                         * flxscale
          !
       ELSE                          ! land or lake
          pamlcorr(jl) = 0._dp
       END IF
    END DO
    RETURN
  END SUBROUTINE mlocean_mlfluxp
  !  -----------------------------------------------------------------------
  !  -----------------------------------------------------------------------
  SUBROUTINE mlocean_readflux(status, fname, fluxcorr)

    ! U. Schlese, DKRZ,  May 1993, original version (readsst)
    ! M. Esch,    MPI,   Sep 2002, modified for flux correction
    ! M. Kunze,   FUB,   Feb 2010, adapted for MESSy SMCL
    
!!$    USE messy_main_tools, ONLY: start_message, end_message
    USE messy_main_blather, ONLY: start_message, end_message
    !
    USE netcdf
    !
    ! I/O
    INTEGER,                     INTENT(out) :: status     ! error status
    REAL(dp), DIMENSION(:,:,0:), INTENT(out) :: fluxcorr
    CHARACTER(LEN=*),            INTENT(in)  :: fname
    !
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='mlocean_readflux'
    INTEGER, PARAMETER    :: nt = 12
    INTEGER               :: nlon, nlat
    INTEGER               :: ncid, varid, latid, i
    INTEGER, DIMENSION(3) :: start       ! start index for NetCDF-read
    INTEGER, DIMENSION(3) :: count       ! number of iterations for NetCDF-read
    REAL(dp), DIMENSION(:),     ALLOCATABLE :: lat
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: tmp 
    LOGICAL               :: lex
    
    CALL start_message(modstr, 'DATA IMPORT', substr)
    !
    nlon = SIZE(fluxcorr,dim=1)
    nlat = SIZE(fluxcorr,dim=2)
    !
    status = 1
    INQUIRE (file=fname, exist=lex)
    IF (lex) THEN
       status = NF90_OPEN (TRIM(fname), NF90_NOWRITE, ncid)
       IF (status /= NF90_NOERR) THEN
          WRITE(*,*) 'Could not open flux file '//TRIM(fname)
          WRITE(*,*) NF90_STRERROR(status)
          RETURN
       ELSE
          WRITE(*,*) 'Open flux file '//TRIM(fname)
       END IF
    ELSE
       WRITE(*,*) 'flux file '//TRIM(fname)//' does not exist.'
       RETURN
    ENDIF
    
    status = NF90_INQ_VARID (ncid, 'aflux', varid)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netcdf error.'
       WRITE(*,*) NF90_STRERROR(status)
       RETURN
    END IF
    
    start = (/    1,    1,  1 /)
    count = (/ nlon, nlat, nt /)
    status = NF90_GET_VAR (ncid, varid, fluxcorr(:,:,1:nt), start, count)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netcdf error.'
       WRITE(*,*) NF90_STRERROR(status)
       RETURN
    ELSE
       WRITE(*,*) 'read flux from file '//TRIM(fname)
    END IF
    !
    ! test the order of the latitudes
    !
    ALLOCATE (lat(nlat))
    status = NF90_INQ_VARID (ncid, 'lat', latid)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netcdf error.'
       WRITE(*,*) NF90_STRERROR(status)
       RETURN
    END IF
    status = NF90_GET_VAR   (ncid, latid, lat(:), (/ 1 /), (/ nlat /))
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netcdf error.'
       WRITE(*,*) NF90_STRERROR(status)
       RETURN
    END IF
    IF (lat(1) < lat(nlat)) THEN
       !
       ! S -> N: COARDS (reorder required)
       !
       ALLOCATE (tmp(nlon,nlat,nt))
       tmp(:,:,:) = fluxcorr(:,:,1:nt)
       DO i=1, nlat
          fluxcorr(:,i,1:nt) = tmp(:,nlat+1-i,:)
       END DO
       DEALLOCATE(tmp)
    END IF
    DEALLOCATE(lat)
    !
    fluxcorr(:,:,0)  = fluxcorr(:,:,12)
    fluxcorr(:,:,13) = fluxcorr(:,:,1)
    !
    status = NF90_CLOSE(ncid)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'SUBROUTINE  readflux netcdf error.'
       WRITE(*,*) NF90_STRERROR(status)
       RETURN
    END IF
    CALL end_message(modstr, 'DATA IMPORT', substr)
    status = 0
    RETURN
  END SUBROUTINE mlocean_readflux
  ! **********************************************************************
END MODULE messy_mlocean
! **********************************************************************
