
MODULE tools 

  !-------------------------------------------------------------------------
  ! files with subroutine used by biogen.f90
  ! 
  !     Laurens Ganzeveld, MPIC-Mainz, 2002
  !     Final code develpment:
  !     Andrea Pozzer, MPIC-Mainz, 2006-2009
  !     Patrick Joeckel, DLR, 2010: minor modifications
  !-------------------------------------------------------------------------

  IMPLICIT NONE

  ! VERSION
  CHARACTER(LEN=*), PARAMETER :: VERSION = '0.6'
  ! ... INTERNAL
  INTEGER, PARAMETER :: iou  = 21       ! I/O unit
  INTEGER, PARAMETER :: str_short = 30  ! length of short strings
  INTEGER, PARAMETER :: str_long  = 100 ! length of long strings
  INTEGER, PARAMETER :: str_vlong = 500 ! length of very long strings
  INTEGER,  PARAMETER :: dp       = SELECTED_REAL_KIND(12,307)
  ! PARAMETER
  CHARACTER(LEN=*), PARAMETER :: hline='==========================================================='
  INTEGER, PARAMETER, PUBLIC :: nmonth = 12
  CHARACTER(LEN=2), DIMENSION(nmonth)  :: mon = (/'01','02','03','04','05','06',  &
                                              '07','08','09','10','11','12'/)
  INTEGER, PARAMETER  :: nlon = 360   ! number of longitude intervals (1 deg)
  INTEGER, PARAMETER  :: nlat = 180   ! number of latitude  intervals (1 deg)
  INTEGER, PARAMETER, PUBLIC :: ILON=720, ILAT=360 ! special for Olson vegetation 
  INTEGER, PARAMETER, PUBLIC :: IILON=2160, IILAT=1080
  ! FROM NAMELIST 
  CHARACTER(LEN=str_vlong), PUBLIC,SAVE :: inputpath_olson = ''
  CHARACTER(LEN=str_vlong), PUBLIC,SAVE :: inputpath_geia = ''
  CHARACTER(LEN=str_vlong), PUBLIC,SAVE :: output_dir = ''
  LOGICAL, PUBLIC,SAVE :: l_output_olson 
  LOGICAL, PUBLIC,SAVE :: l_flux
  LOGICAL, PUBLIC,SAVE :: l_ocean_alt ! alternative input file for oceanic emission
  REAL, PUBLIC, SAVE :: olson(nlon,nlat,nmonth)          
  REAL, PUBLIC, SAVE :: land(nlon,nlat,nmonth)            
  INTEGER, PUBLIC, SAVE :: mask
  ! VARIABLE TO PRODUCE
  REAL,PUBLIC,SAVE :: landmask(nlon,nlat,nmonth)
  INTEGER, PUBLIC, PARAMETER :: NMAXTRAC =100        !maximum number of tracer possible 
  INTEGER, PUBLIC, SAVE ::  TRAC_TOT                 !total number of tracer defined 
  CHARACTER(LEN=str_short), PUBLIC, SAVE     :: TRAC_NAME(NMAXTRAC)     !name of tracer
  REAL,DIMENSION(NMAXTRAC), PUBLIC,SAVE :: TRAC_MOLAR_MASS              ! 
  REAL,DIMENSION(NMAXTRAC), PUBLIC,SAVE :: TRAC_LAND                    ! 
  REAL,DIMENSION(NMAXTRAC), PUBLIC,SAVE :: TRAC_OCEAN                   ! 
  INTEGER, DIMENSION(NMAXTRAC), PUBLIC,SAVE ::  TRAC_YEAR               ! year emissions 
  LOGICAL,DIMENSION(NMAXTRAC), PUBLIC, SAVE     :: L_TERP            ! basic distribution used 
  LOGICAL,DIMENSION(NMAXTRAC), PUBLIC, SAVE     :: L_NVOC            
  LOGICAL,DIMENSION(NMAXTRAC), PUBLIC, SAVE     :: L_ISOP            
  ! FOR NETCDF
  !CHARACTER(LEN=str_long)         :: s_freq    = ''
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: output_field
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: ocean_field
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: land_field
  REAL, DIMENSION(:),SAVE, ALLOCATABLE :: xlon    ! longitudes [deg]
  REAL, DIMENSION(:),SAVE, ALLOCATABLE :: xlat    ! latitudes  [deg]
  REAL, DIMENSION(:),SAVE, ALLOCATABLE :: xlev    ! levels [1]
  REAL, DIMENSION(:),SAVE, ALLOCATABLE :: xtime   ! time [1]
  INTEGER, PARAMETER  :: nlev  = 1 ! levels
  INTEGER, SAVE       :: ntime = 12 ! time steps
  CHARACTER(LEN=str_vlong), PUBLIC :: OUTPUT      = ''     ! netCDF-file
  CHARACTER(LEN=str_vlong) , PUBLIC :: SPECIES     = ''
  CHARACTER(LEN=str_long) , PUBLIC :: UNIT_OUTPUT = ''
  CHARACTER(LEN=str_vlong) , PUBLIC ::LONG_NAME_OUTPUT = ''
  INTEGER, PUBLIC, SAVE     :: YEAR      = 0        
  REAL,PUBLIC, SAVE     :: TOTAL = 0.0         
  REAL,PUBLIC, SAVE     :: MOLARMASS = 0.0 ! [g/mol]        
  REAL, PUBLIC, SAVE :: sland = 0.0
  REAL, PUBLIC, SAVE :: socean = 0.0
  LOGICAL, PUBLIC, SAVE :: lterp = .FALSE.
  LOGICAL, PUBLIC, SAVE :: lnvoc = .FALSE.
  LOGICAL, PUBLIC, SAVE :: lisop = .FALSE.
  ! CONSTANTS
  REAL, PARAMETER :: r_earth  = 6371000.0! radius of the earth in m
  REAL, PARAMETER :: pi       = 3.14159265358979323846
  REAL, PARAMETER :: N_A      = 6.022045E23 ! Avogadro constant [1/mol]
  INTEGER  :: spy    = 365 * 24 * 3600  ! seconds per year
  ! BASIC DISTRIBUTION
  REAL, PUBLIC :: surfarea(nlon,nlat)       ! surface area of grid boxes (m^2)  ! rvk now in surfco.com
  REAL,PUBLIC,SAVE :: base_dist(nlon,nlat,nmonth)   ! biogenic  surface emission rate, kg/m2/s 
  PUBLIC :: read_olson
  PUBLIC :: calcland
  PUBLIC :: nc_creation 
  PUBLIC :: NFERR
  PUBLIC :: areacal
  PUBLIC :: land_distr
  PUBLIC :: ocean_distr
  PUBLIC :: add_geia
  PUBLIC :: add_distr
  PUBLIC :: mask_distr

CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE read_olson(status, iou, path, olson, land)

     IMPLICIT NONE

     REAL :: tmp(ilon,ilat,nmonth)
     INTEGER :: olsveg(ilon,ilat,nmonth)                 !0.5x0.5
     REAL :: lai(ilon,ilat,nmonth)
     REAL :: laimax(ilon,ilat,nmonth)
     REAL :: dm(ilon,ilat,nmonth)
     REAL, INTENT(inout) :: olson(:,:,:)
     REAL, INTENT(inout) :: land(:,:,:)
     CHARACTER(LEN=*), INTENT(IN)  :: path
     INTEGER,          INTENT(IN)  :: iou
     INTEGER,          INTENT(OUT) :: status
     !local
     CHARACTER(LEN=str_vlong)  :: fname 
     INTEGER :: month, IP, I, J
  
     status = 0
     write (*,*) hline
     write (*,*) " READING OLSON VEGETATION INDEX"
     write (*,*) " resolution : LON",ILON,"LAT",ILAT
     write (*,*) hline
     DO month=1,nmonth
     WRITE (*,'(1a,i3)')' Month: ',month
     write (*,*) MON(month)
     FNAME=path//'/veg_Ols_'
     IP=INDEX(FNAME,' ')
     FNAME=FNAME(1:IP-1)//MON(month)
#if defined(__ibm__) || defined(__G95__) || defined(__GFORTRAN__)
     WRITE(*,*) 'OPEN(...,FORM=BINARY,...) NOT ALLOWED.'
     STOP
#else
!!$ NAG compliance +
     ! see
     ! https://stackoverflow.com/questions/36717653/standard-equivalent-for-intel-fortran-form-binary
     ! or
     ! http://fortranwiki.org/fortran/show/Stream+Input+Output
!!$     OPEN(UNIT=iou,FILE=FNAME,FORM='BINARY',                        &
!!$              STATUS='UNKNOWN')
     OPEN(UNIT=iou,FILE=FNAME,ACCESS='STREAM',FORM='UNFORMATTED',   &
          STATUS='UNKNOWN')
!!$ NAG compliance -
#endif
     !     reading of file with vegetation data derived from the
     !     Olson ecosystem database and the satellite data, if this
     !     file is not available for the the specific month
     !     this file is created (ERR=999)
     REWIND(iou)
     READ(iou) ((OLSVEG(I,J,month),I=1,ilon),J=1,ilat)
     READ(iou) ((LAI(I,J,month),I=1,ilon),J=1,ilat)
     READ(iou) ((LAIMAX(I,J,month),I=1,ilon),J=1,ilat)
     READ(iou) ((DM(I,J,month),I=1,ilon),J=1,ilat)
     CLOSE(iou)
     END DO

     !------------------------OLSON VEGETATION-------------------
     olson(:,:,:)=0
     tmp(:,:,:)=0
     do i=1,ilon
     do j=1,ilat
     do month=1,nmonth
       tmp(i,j,month)=real(olsveg(i,ilat+1-j,month))
     end do
     end do
     end do
     do i=1,nlon
     do j=1,nlat
       land(i,j,:)=(tmp(i*2,j*2,:)+tmp(i*2-1,j*2,:)+tmp(i*2,j*2-1,:)+tmp(i*2-1,j*2-1,:))/4 
       olson(i,j,:)=tmp(i*2-1,j*2,:) ! approx
     end do
     end do

  END SUBROUTINE read_olson

!############################################################################

  SUBROUTINE calcland(field, landmask, mask)

     IMPLICIT NONE
     REAL, INTENT(IN) :: field(:,:,:)
     INTEGER, INTENT(IN) :: mask
     REAL, INTENT(OUT) :: landmask(:,:,:)
     INTEGER :: i,j,m

     do i=1,nlon
       do j=1,nlat
         do m=1,nmonth
           IF (field(i,j,m).GT.0) THEN
              IF (mask.eq.1) THEN
                landmask(i,j,:) = 1. 
              ELSE
                landmask(i,j,:) = 0. 
              ENDIF
           ELSE 
              IF (mask.eq.1) THEN
                 landmask(i,j,:) = 0.
              ELSE
                landmask(i,j,:) = 1. 
              ENDIF
           ENDIF
         end do
       end do
     end do

  END SUBROUTINE calcland

  !##############################################################################

  SUBROUTINE land_distr(iou,inclterp,inclnvoc,inclisop,coemo )
  

    !----------------------------------------------------------------------
    !     Get the distribution of NMHC effective surface emissions from GEIA
    !     based on RvKuhlmann code!
    !----------------------------------------------------------------------
    !-----------------------------------------------------------------------


    IMPLICIT NONE

    !     Input Parameters

    integer, INTENT(IN) :: iou
    logical, INTENT(IN) :: inclterp          ! should the terpene distribution be included ?
    logical, INTENT(IN) :: inclisop          ! should the isoprene distribution be included ?
    logical, INTENT(IN) :: inclnvoc          ! should the voc distribution be included ?
    REAL, INTENT(OUT) :: coemo(nlon,nlat,nmonth)
    !     Local Workspace
    integer ::i,j, month, m  ! longitude, latitude, month indices
    ! GAIA GRID SPECIFICATION
    INTEGER, PARAMETER ::  glon = 360
    INTEGER, PARAMETER ::  glat = 180
    REAL :: conv ! conversion factor
    REAL :: coemi(360,180,nmonth)

    integer :: mday(12)        ! days per month
    REAL :: secpermon  


    !     initialize

    coemi(:,:,:) = 0.
    coemo(:,:,:) = 0.

    ! only include this if you want a mixture of both distributions (for CO,Acetone)
    IF ( inclterp ) THEN
       CALL add_geia( 'terp.1a', iou, glon, glat, coemi )
    ENDIF

    ! mainly for isoprene (see onlem)
    IF ( inclisop ) THEN
       CALL add_geia( 'isop90.1a', iou, glon, glat, coemi )
    ENDIF

    IF ( inclnvoc ) THEN
       CALL add_geia( 'nvoc90.1a', iou, glon, glat, coemi )
    ENDIF

    !***************************************************************

    !     Next convert input data from mgC/m2/mo to kg(CO)/m2/s

    conv = 1.e-6*(28./12.)
    DATA mday /31,28,31,30,31,30,31,31,30,31,30,31/
    do month = 1,nmonth
       secpermon = 86400.*mday(month)
       do j = 1,glat
          do i = 1,glon
             coemi(i,j,month) = ( conv/secpermon ) * coemi(i,j,month) 
             ! LG- since the interpolation is not required the value 
             !     coemi is assigned to coemo
             coemo(i,j,month) = coemi(i,j,month)
          end do
       end do
    end do

  END SUBROUTINE land_distr



    !***************************************************************


  SUBROUTINE ocean_distr(coemo)

    !----------------------------------------------------------------------
    !     Set the distribution of oceanic CO surface emissions
    !----------------------------------------------------------------------

    IMPLICIT NONE

    !     Input Parameters

    real, intent(out) ::  coemo(nlon,nlat,nmonth)  ! output surface emission rate, kg/m2/s

    !     Local Workspace

    integer :: j, month,m,i,n          ! indices

    integer,parameter ::nregs=10

    real :: bounds(nregs+1)      ! boundaries for regions
    real :: area(nregs)          ! ocean area of regions (10^6 km2)
    real :: annem(nregs)         ! annual emissions (Gmol/yr) in region
    real :: regflux(nregs)       ! flux (kg/m2/s) in region

    real :: latval

    !----------------------------------------------------------------------
    !     program to set Ocean surface CO surface emissions rate
    !     the data are from Bates et al., JGR, 1995, p. 23093, table 4
    !     note that the initial orography field is used, which could perhaps
    !     miss some ocean regions due to sea ice, but this is doubtfully 
    !     important because the bounds only go to +/- 75 deg lat

    !----------------------------------------------------------------------

    !     data from Bates et al.

    DATA bounds / -75,-60,-45,-30,-15,0,15,30,45,60,75 /

    DATA area / 19.3,39.0,46.8,46.1,49.8,49.5,39.0,27.9,16.6,8.5 /

    DATA annem / 22,63,90,44,95,64,34,32,15,5 /

    !     convert to kg/m2/s

    do n = 1,nregs
       regflux(n) = ( annem(n)/(area(n)*1.e6*1.e6) )       &! Gmol/m2/yr
            &        * (1./spy)                         &! Gmol/m2/s                  
            &        * 1.e9*(28.e-3)                          ! kg/m2/s
       !debug         write(*,*) n,regflux(n)
    end do

    !     fill into model fluxes

    coemo(:,:,:)=0
    do j = 1,nlat
       latval = j-90. 
       do n = 1,nregs
          if (latval.gt.bounds(n) .and. latval.le.bounds(n+1)) then
             do m = 1,12
                do i = 1,nlon
                   coemo(i,j,m) = regflux(n)
                end do
             end do
          end if
       end do
    end do

  END SUBROUTINE ocean_distr

  !##############################################################################

  SUBROUTINE mask_distr (distr, mask)

    IMPLICIT NONE

    REAL, INTENT(INOUT) :: distr(nlon,nlat,nmonth)
    REAL, INTENT(IN) :: mask(nlon,nlat,nmonth)

    distr=distr*mask 

  END SUBROUTINE mask_distr

  !##############################################################################

  !##############################################################################

  SUBROUTINE add_geia( filename, iou, glon, glat, coemi )

    !----------------------------------------------------------------------
    !     Add in the distribution of NMHC surface emissions from GEIA
    !----------------------------------------------------------------------

    IMPLICIT NONE

    character(LEN=*), intent(in) :: filename
    integer, intent(in) ::iou 
    integer, intent(in) ::glon 
    integer, intent(in) ::glat 
    real, intent(inout)::coemi(glon,glat,nmonth)  ! input surface emission rate, mg(C)/m2/mo
    integer i,j,jj,n,month
    integer igrid
    real x(nmonth)

    !----------------------------------------------------------------------

    !write(*,*) TRIM(inputpath_geia)//'/'//TRIM(filename)
    open( iou, status='old',file= TRIM(inputpath_geia)//'/'//TRIM(filename))

    !     read header
    do n = 1,10
       read(iou,*)
    end do

    !     read in data
    do n = 1,nlon*nlat+1
       read(iou,'(i6,12(1x,e10.5))',end=100) igrid,x
       !         write(*,*) i,x
       !     set grid
       j = int(float(igrid)/1000.)

       i = igrid - j*1000
       do month = 1,nmonth
          coemi(i,j,month) = coemi(i,j,month) + x(month)
       end do
       !         write(*,*) i,j,x
       !         stop
    end do

    stop 'something wrong reading geia data'

100 continue
    close( iou )
    !     write(55,*) 'added in emissions from ',filename,', nlocs=',n-1
    return

  END SUBROUTINE add_geia

  !##############################################################################

  SUBROUTINE add_distr ( destot, distr, addon )

    !----------------------------------------------------------------------
    !  add any distribution of a model grid to the array "addon", with a
    !  specified total of destot ( arbitry units, conversion has to be done elsewhere )
    !----------------------------------------------------------------------

    implicit none

    !----------------------------------------------------------------------

    !     Input Parameters
    real, intent(in) :: destot             ! desired total of the field
    real, intent(in) :: distr(nlon,nlat,nmonth)! the distribution to be used, unit: arbitrary
    !     Input / output
    real, intent(inout) :: addon(nlon,nlat,nmonth) ! the field where to add on the scaled field

    !     locals
    integer i,j,m            ! indices
    real tot                 ! total of raw data
    integer :: mday(nmonth)        ! days per month
    REAL :: secpermon  
    !-----------------------------------------------------------------------
    !     written by Rolf von Kuhlmann (April 1999, last changes April 2000) for 
    !     MATCH-MPIC uersion 2.0, Author: mgl, March 1998
    !     corrected : Andrea Pozzer MPIC 2006
    !-----------------------------------------------------------------------
    !     calculate the raw total of the sample distribution (units arbitrary)
    tot = 0.
    DATA mday /31,28,31,30,31,30,31,31,30,31,30,31/
    do m=1,nmonth
       do j=1,nlat
          do i=1,nlon           ! convert from
             secpermon = 86400.*mday(m)
             tot = tot +  distr(i,j,m)*surfarea(i,j)*secpermon 
          enddo
       enddo
    enddo

    ! convert to Tg/y ! 1E12 is to convert to Tg
    if (tot.gt.0.) addon = addon + (destot/tot) * distr*1E12

  END SUBROUTINE add_distr
  ! ------------------------------------------------------------------------
  SUBROUTINE nc_creation

    USE netcdf

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, CHAR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'nc_creation'
    INTEGER :: ncid      ! netCDF-ID
    INTEGER :: dimid_lat, dimid_lon, dimid_lev, dimid_time
    INTEGER :: varid_lat, varid_lon, varid_lev, varid_time
    INTEGER :: varid_mass, varid_flux, varid_height
    INTEGER                  :: jc
    CHARACTER(LEN=str_long)  :: timestr = ''
    CHARACTER(LEN=4)         :: yrstr = '    '
    CHARACTER(LEN=1000 + NMAXTRAC*50) :: nmlstr = ''
    CHARACTER(LEN=4)         :: jcstr = '    '
    CHARACTER(LEN=4)         :: levstr = '    '
    CHARACTER(LEN=str_long)  :: scalestr = ''
    CHARACTER(LEN=str_long)  :: mmassstr = ''
    !
    CHARACTER(LEN=str_long)  :: landstr = ''
    CHARACTER(LEN=str_long)  :: oceanstr = ''
    CHARACTER(LEN=str_long)  :: lterpstr = ''
    CHARACTER(LEN=str_long)  :: lnvocstr = ''
    CHARACTER(LEN=str_long)  :: lisopstr = ''
    !
    CHARACTER(LEN=8)         :: date
    CHARACTER(LEN=10)        :: time
    CHARACTER(LEN=5)         :: zone

    ! INIT
    ! CONVERT TO STRING
    WRITE(yrstr,'(i4)') YEAR
    CALL DATE_AND_TIME(date, time, zone)
    WRITE(mmassstr,*)  MOLARMASS
    WRITE(lterpstr,*) lterp
    WRITE(lisopstr,*) lisop
    WRITE(lnvocstr,*) lnvoc
    WRITE(landstr,*) sland
    WRITE(oceanstr,*) socean

    ! OUTPUT NAME
    !OUTPUT=TRIM(OUTPUT)//TRIM(SPECIES)//'_biogen_'//yrstr//'.nc'
    OUTPUT=TRIM(OUTPUT)//'_'//yrstr//'01_'//yrstr//'12'//'.nc'
    !write (*,*) 'output file : ',OUTPUT
    ! CREATE NEW FILE
    CALL NFERR( &
         nf90_create(TRIM(OUTPUT), NF90_CLOBBER, ncid) &
         ,1)

    ! ADD GLOBALE ATTRIBUTES
    ! - VERSION
    CALL NFERR( &
         nf90_put_att(ncid, NF90_GLOBAL, 'created_with',     &
         'biogen2nc Version '//TRIM(VERSION)) &
         ,2)
    ! - DATE AND TIME
    CALL NFERR( &
         nf90_put_att(ncid, NF90_GLOBAL, 'date', date) &
         ,3)
    CALL NFERR( &
         nf90_put_att(ncid, NF90_GLOBAL, 'time', TRIM(time)//TRIM(zone)) &
         ,4)

    ! - NAMELIST
    nmlstr = 'FILE= '''//TRIM(OUTPUT)//'''; '//CHAR(10)//'NML: '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'TRAC_NAME = '''//TRIM(SPECIES)//''', '//CHAR(10)
 nmlstr = TRIM(nmlstr)//'TRAC_MOLAR_MASS = '''//TRIM(mmassstr)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'TRAC_YEAR = '//yrstr//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'TRAC_LAND ='//TRIM(landstr)//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'TRAC_OCEAN ='//TRIM(oceanstr)//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'L_TERP = '//TRIM(lterpstr)//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'L_NVOC = '//TRIM(lnvocstr)//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'L_ISOP = '//TRIM(lisopstr)//', '//CHAR(10)

    CALL NFERR( &
         nf90_put_att(ncid, NF90_GLOBAL, 'namelist', nmlstr) &
         ,5)

    ! DEFINE DIMENSIONS
    CALL NFERR( &
         nf90_def_dim(ncid, 'lon', nlon, dimid_lon) &
         ,6)
    CALL NFERR( &
         nf90_def_dim(ncid, 'lat', nlat, dimid_lat) &
         ,7)
    CALL NFERR( &
         nf90_def_dim(ncid, 'lev', nlev, dimid_lev) &
         ,5)
    CALL NFERR( &
         nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid_time) &
         ,8)

    ! DEFINE COORDINATE VARIABLES WITH ATTRIBUTES
    CALL NFERR( &
         nf90_def_var(ncid, 'lon', NF90_FLOAT, (/ dimid_lon /), varid_lon) &
         ,9)
    CALL NFERR( &
         nf90_put_att(ncid, varid_lon, 'long_name', 'longitude') &
         ,10)
    CALL NFERR( &
         nf90_put_att(ncid, varid_lon, 'units', 'degrees_east') &
         ,11)

    CALL NFERR( &
         nf90_def_var(ncid, 'lat', NF90_FLOAT, (/ dimid_lat /), varid_lat) &
         ,12)
    CALL NFERR( &
         nf90_put_att(ncid, varid_lat, 'long_name', 'latitude') &
         ,13)
    CALL NFERR( &
         nf90_put_att(ncid, varid_lat, 'units', 'degrees_north') &
         ,14)

    CALL NFERR( &
         nf90_def_var(ncid, 'lev', NF90_FLOAT, (/ dimid_lev /), varid_lev) &
         ,15)
    CALL NFERR( &
         nf90_put_att(ncid, varid_lev, 'long_name', 'level index') &
         ,16)
    CALL NFERR( &
         nf90_put_att(ncid, varid_lev, 'units', 'level') &
         ,17)

    CALL NFERR( &
         nf90_def_var(ncid, 'time', NF90_FLOAT, (/ dimid_time /), varid_time) &
         ,18)
    CALL NFERR( &
         nf90_put_att(ncid, varid_time, 'long_name', 'time') &
         ,19)
    !
    !SELECT CASE(TRIM(ADJUSTL(s_freq)))
    !CASE('annual')
    !   timestr = 'year since '//yrstr//'-01-01 00:00:00'
    !CASE('seasonal')
    !   timestr = 'season since '//yrstr//'-01-01 00:00:00'
    !CASE('monthly')
       timestr = 'month since '//yrstr//'-01-01 00:00:00'
    !CASE DEFAULT
    !   WRITE(*,*) substr,': UNKNOWN DATA FREQUENCY ',TRIM(ADJUSTL(s_freq))
    !   STOP
    !END SELECT
    CALL NFERR( &
         nf90_put_att(ncid, varid_time, 'units', timestr) &
         ,20)

!    ! DEFINE VARIABLES
!    ! - emission height
!    CALL NFERR( &
!         nf90_def_var(ncid, 'height', NF90_FLOAT  &
!         , (/ dimid_lev /), varid_height) &
!         ,21)
!    CALL NFERR( &
!         nf90_put_att(ncid, varid_height, 'long_name' &
!         , 'emission height') &
!         ,22)
!    CALL NFERR( &
!         nf90_put_att(ncid, varid_height, 'units', 'm') &
!         ,23)
!
!    ! - mass flux
!    IF (L_MASSFLUX) THEN
!       CALL NFERR( &
!            nf90_def_var(ncid, TRIM(SPECIES)//'_mflux', NF90_FLOAT  &
!            , (/ dimid_lon, dimid_lat, dimid_lev, dimid_time /), varid_mass) &
!            ,24)
!       CALL NFERR( &
!            nf90_put_att(ncid, varid_mass, 'long_name' &
!            , 'mass flux of '//TRIM(s_species)) &
!            ,25)
!       CALL NFERR( &
!            nf90_put_att(ncid, varid_mass, 'units', TRIM(s_unit)) &
!            ,26)
!       CALL NFERR( &
!            nf90_put_att(ncid, varid_mass, 'molar_mass', MOLARMASS) &
!            ,27)
!    END IF
!
!    ! - flux
    IF (l_flux) THEN
       CALL NFERR( &
            nf90_def_var(ncid, TRIM(SPECIES)//'_flux', NF90_FLOAT  &
            , (/ dimid_lon, dimid_lat, dimid_lev, dimid_time /), varid_flux) &
            ,28)
    ELSE
       CALL NFERR( &
            nf90_def_var(ncid, TRIM(SPECIES), NF90_FLOAT  &
            , (/ dimid_lon, dimid_lat, dimid_lev, dimid_time /), varid_flux) &
            ,28)
    ENDIF
    CALL NFERR( &
         nf90_put_att(ncid, varid_flux, 'long_name' &
         , TRIM(LONG_NAME_OUTPUT)) &
         ,29)
    CALL NFERR( &
         nf90_put_att(ncid, varid_flux, 'units', TRIM(UNIT_OUTPUT)) &
         ,30)
    CALL NFERR( &
         nf90_put_att(ncid, varid_flux, 'molar_mass', MOLARMASS) &
         ,31)

    ! SWITCH MODUS
    CALL NFERR( &
         nf90_enddef(ncid) &
         ,32)

    ! SAVE COORDINATE VARIBLES
    CALL NFERR( &
         nf90_put_var(ncid, varid_lon, xlon) &
         ,33)
    CALL NFERR( &
         nf90_put_var(ncid, varid_lat, xlat) &
         ,34)
    CALL NFERR( &
         nf90_put_var(ncid, varid_lev, xlev) &
         ,35)
    CALL NFERR( &
         nf90_put_var(ncid, varid_time, xtime) &
         ,36)
!    CALL NFERR( &
!         nf90_put_var(ncid, varid_height, HEIGHT(1:nlev)) &
!         ,37)
!
!    ! SAVE VARIABLES
!    IF (L_MASSFLUX) THEN
!       CALL NFERR( &
!            nf90_put_var(ncid, varid_mass, emismass) &
!            ,38)
!    END IF
!
    CALL NFERR( &
         nf90_put_var(ncid, varid_flux, output_field) &
         ,39)

    ! CLOSE FILE
    CALL NFERR( &
         nf90_close(ncid) &
         ,40)


!--------------------------------------------------------------------------------
! DIAGNOSTIC OUTPUT
!--------------------------------------------------------------------------------
  WRITE(*,*) '==========================================================='
  WRITE(*,*) 'OUTPUT         : ', TRIM(OUTPUT)
  WRITE(*,*) 'SPECIES        : ', TRIM(SPECIES)
  WRITE(*,*) 'MOLAR MASS     : ', MOLARMASS
!  WRITE(*,*) 'GLOBAL SCALING : ', GLOBALSCALE
!  WRITE(*,*) 'MASS FLUX      : ', L_MASSFLUX
  WRITE(*,*) 'YEAR           : ', YEAR
!  WRITE(*,*) 'HEIGHT         : ', HEIGHT(1:nlev)
!  WRITE(*,*) 'SPECIES NAME   : ',TRIM(s_species)
!  WRITE(*,*) 'UNIT           : ',TRIM(s_unit)
!  WRITE(*,*) 'FREQUENCY      : ',TRIM(s_freq)
 ! WRITE(*,*)
 ! WRITE(*,*) 'MASS FLUX      : '
 ! WRITE(*,*) ' LBOUND        : ',LBOUND(emismass)
 ! WRITE(*,*) ' UBOUND        : ',UBOUND(emismass)
 ! WRITE(*,*) ' MIN           : ',MINVAL(emismass)
 ! WRITE(*,*) ' MAX           : ',MAXVAL(emismass)
 ! WRITE(*,*) ' SUM           : ',SUM(emismass)
 ! WRITE(*,*)
 ! WRITE(*,*) 'FLUX           : '
 ! WRITE(*,*) ' AREA CHECKSUM : ',af_sum
  WRITE(*,*) ' LBOUND        : ',LBOUND(output_field)
  WRITE(*,*) ' UBOUND        : ',UBOUND(output_field)
  WRITE(*,*) ' MIN           : ',MINVAL(output_field)
  WRITE(*,*) ' MAX           : ',MAXVAL(output_field)
  WRITE(*,*) '==========================================================='

  END SUBROUTINE nc_creation
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE NFERR(status, pos)

    USE netcdf, ONLY: NF90_NOERR, nf90_strerror
    
    IMPLICIT NONE
    
    ! I/O
    INTEGER,          INTENT(IN) :: status
    INTEGER,          INTENT(IN) :: pos
    
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netCDF ERROR at position: ', pos
       WRITE(*,*) 'netCDF ERROR status     : ',status
       WRITE(*,*) 'netCDF ERROR            : ',nf90_strerror(status)
    END IF
  
  END SUBROUTINE NFERR
  ! ------------------------------------------------------------------
  SUBROUTINE areacal(im,jm,dxyp)

    integer,INTENT(IN) :: im
    integer,INTENT(IN) :: jm
    real, INTENT(OUT) :: dxyp(im,jm)
    REAL :: twopi, dlat, fjeq, dd
    INTEGER :: j,jj

    twopi=2.*pi
    dlat=.5*twopi/(jm-1)
    fjeq=.5*(1+jm)
    dd=2.*(r_earth**2)*(twopi/im)*sin(0.5*dlat)
    !areag=0.

    do j=1,jm
       jj=int(j-0.5*jm)
       dxyp(:,j)=dd*cos(dlat*(j-fjeq))
    enddo

  END SUBROUTINE areacal


 END MODULE tools
