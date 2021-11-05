! ******************************************************************
! ------------------------------------------------------------------
PROGRAM EDGAR2NC
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, October 2004
! ******************************************************************
!
! CONVERT EDGAR-DATABASE EMISSION FILES FROM GEIA-FORMAT TO netCDF
!
! -----------------------------------------------------------------

  USE mo_f2kcli                    ! command line interface

  IMPLICIT NONE

  INTRINSIC :: ADJUSTL, ALLOCATED, ASSOCIATED, MAXVAL, MINVAL, SIZE &
             , SUM, TRIM, UBOUND, LBOUND, REAL, SIN, SELECTED_REAL_KIND

  ! VERSION
  CHARACTER(LEN=*), PARAMETER :: VERSION = '1.4b'

  ! FOR COMMAND LINE
  CHARACTER(LEN=256) :: EXE          ! program name
  CHARACTER(LEN=80)  :: CMD          ! argument
  INTEGER            :: NARG         ! number of arguments

  ! CONSTANTS
  INTEGER,  PARAMETER :: dp       = SELECTED_REAL_KIND(12,307)
  REAL(DP), PARAMETER :: r_earth  = 6371000.0_dp! radius of the earth in m
  REAL(DP), PARAMETER :: pi       = 3.14159265358979323846_dp
  REAL(DP), PARAMETER :: N_A      = 6.022045E23_dp ! Avogadro constant [1/mol]

  ! PARAMETER
  INTEGER, PARAMETER  :: MAX_NCLASS  = 10000 ! max. number of emission-classes
  INTEGER, PARAMETER  :: MAX_HEIGHTS =   100 ! max. number of emission heights
  INTEGER, PARAMETER  :: nlon = 360   ! number of longitude intervals (1 deg)
  INTEGER, PARAMETER  :: nlat = 180   ! number of latitude  intervals (1 deg)
  REAL(DP), PARAMETER :: dhlat = 0.5_DP ! half latitude width  [deg]
  ! ... INTERNAL
  INTEGER, PARAMETER :: str_short = 30  ! length of short strings
  INTEGER, PARAMETER :: str_long  = 100 ! length of long strings
  INTEGER, PARAMETER :: str_vlong = 500 ! length of very long strings
  INTEGER, PARAMETER :: iou  = 21       ! I/O unit

  ! NAMELIST
  TYPE CLASS_IO
     CHARACTER(LEN=str_long) :: fname = ''
     INTEGER                 :: level = 1
     REAL(DP)                :: scale = 1.0_DP
  END TYPE CLASS_IO
  !
  CHARACTER(LEN=str_long)               :: OUTPUT    = ''     ! netCDF-file
  CHARACTER(LEN=str_long)               :: SPECIES   = ''
  REAL(DP)                              :: MOLARMASS = 0.0_DP ! [g/mol]
  ! global scaling factor
  REAL(DP)                              :: GLOBALSCALE = 0.0_DP
  LOGICAL                               :: L_MASSFLUX = .FALSE.
  INTEGER                               :: YEAR      = 0        
  REAL(DP), DIMENSION(MAX_HEIGHTS)      :: HEIGHT             ! [m]
  CHARACTER(LEN=str_vlong)              :: INPUTPATH = ''
  TYPE(CLASS_IO), DIMENSION(MAX_NCLASS) :: CLASS              ! emission class
  LOGICAL,        DIMENSION(MAX_NCLASS) :: lfile_ok

  ! VARIABLES
  INTEGER                         :: status       ! status flag
  INTEGER                         :: jc           ! class counter
  INTEGER                         :: ji, jj, jk, jl
  INTEGER                         :: nc           ! act. number of classes
  ! ... PER DATA-FILE
  CHARACTER(LEN=str_long)             :: x_species = ''
  CHARACTER(LEN=str_long)             :: x_freq = ''
  CHARACTER(LEN=str_long)             :: x_unit = ''
  INTEGER                             :: x_year
  REAL(DP), DIMENSION(:,:,:), POINTER :: x_data
  ! ... FOR CHECKING
  CHARACTER(LEN=str_long)         :: s_species = ''
  CHARACTER(LEN=str_long)         :: s_freq    = ''
  CHARACTER(LEN=str_long)         :: s_unit    = ''
  ! ... DATA AND GRID
  INTEGER                                   :: nlev  = 1 ! levels
  INTEGER                                   :: ntime = 1 ! time steps
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: emismass
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: emisflux
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xlon    ! longitudes [deg]
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xlat    ! latitudes  [deg]
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xlev    ! levels [1]
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xtime   ! time [1]
  ! ... CONVERSION
  REAL(DP) :: af_lat = 0.0_DP
  REAL(DP) :: af_lon = 0.0_DP
  REAL(DP) :: af     = 0.0_DP
  REAL(DP) :: af_sum = 0.0_DP
  REAL(DP) :: area   = 0.0_DP
  INTEGER  :: spy    = 365 * 24 * 3600  ! seconds per year
  REAL(DP) :: conv                      ! conversion

  ! (1) READ COMMAND LINE
  NARG = COMMAND_ARGUMENT_COUNT()    ! number of arguments
  CALL GET_COMMAND_ARGUMENT(0,EXE)   ! program name
  !
  IF (NARG > 1) THEN
     WRITE(*,*) 'COMMAND-LINE ERROR: TOO MANY ARGUMENTS !'
     CALL USAGE(TRIM(EXE))
     STOP
  END IF
  !
  IF (NARG == 0) THEN
     CALL USAGE(TRIM(EXE))
     STOP
  END IF
  !
  CALL GET_COMMAND_ARGUMENT(1,CMD)

  ! (2) INIT
  nc = 0                 ! act. number of classes
  HEIGHT(:) = 0.0        ! emission height [m]
  ALLOCATE(xlon(nlon))
  DO ji=1, nlon
     xlon(ji) = -180.5_DP + REAL(ji, DP)
  END DO
  ALLOCATE(xlat(nlat))
  DO jj=1, nlat
     xlat(jj) = -90.5_DP + REAL(jj, DP)
  END DO
  lfile_ok(:) = .FALSE.  ! NO INPUT FILE, OR NOT PRESENT

  ! (3) READ NAMELIST-FILE
  CALL read_nml(status, iou, TRIM(CMD))
  IF (status /= 0) STOP

  ! (4) GET NUMBER OF EMISSION LEVELS
  DO jc = 1, MAX_NCLASS
     IF (class(jc)%level > nlev) nlev = class(jc)%level
  END DO
  ALLOCATE(xlev(nlev))
  DO jk=1, nlev
     xlev(jk) = REAL(jk, DP)
  END DO

  ! (5) LOOP OVER EMISSION-CLASSES
  class_loop: DO jc = 1, MAX_NCLASS
     
     ! SKIP IF NAME IS EMPTY
     IF (TRIM(class(jc)%fname) == '') CYCLE

     CALL read_geia(status, iou                                  &
          , TRIM(INPUTPATH)//'/'//TRIM(ADJUSTL(class(jc)%fname)) &
          , x_species, x_year, x_freq, x_unit, x_data)
     IF (status < 0) CYCLE  ! FILE DOES NOT EXIST
     IF (status > 0) STOP   ! ERROR
     ! STATUS == 0 -> CONTINUE
     lfile_ok(jc) = .TRUE.  ! INPUT FILE FOR CLASS jc IS PRESENT

     nc = nc + 1

     IF (nc == 1) THEN
        ! FIRST FILE
        s_species = TRIM(x_species)
        s_freq    = TRIM(x_freq)
        s_unit    = TRIM(x_unit)
        ntime     = SIZE(x_data,3)
        ALLOCATE(emismass(nlon, nlat, nlev, ntime))
        emismass(:,:,:,:) = 0.0
        ALLOCATE(emisflux(nlon, nlat, nlev, ntime))
        emisflux(:,:,:,:) = 0.0
        ALLOCATE(xtime(ntime))
        DO jl=1, ntime
           xtime(jl) = REAL(jl-1, DP)
        END DO
     END IF

     IF (TRIM(s_species) /= TRIM(x_species))        &
          WRITE(*,*) ' WARNING : SPECIES MISMATCH ' &
          ,TRIM(s_species),' <-> ',TRIM(x_species)

     IF (TRIM(s_unit) /= TRIM(x_unit))           &
          WRITE(*,*) ' WARNING : UNIT MISMATCH ' &
          ,TRIM(s_unit),' <-> ',TRIM(x_unit)

     IF (year /= x_year) &
          WRITE(*,*) ' WARNING : YEAR MISMATCH ' &
          ,year,' <-> ',x_year

     IF (TRIM(s_freq) /= TRIM(x_freq)) THEN
          WRITE(*,*) ' ERROR   : FREQUENCY MISMATCH ' &
             ,TRIM(s_freq),' <-> ',TRIM(x_freq)
        STOP
     END IF

     IF (SIZE(x_data,3) /= ntime) THEN
        WRITE(*,*) ' ERROR   : TIMESTEP MISMATCH ' &
             ,ntime,' <-> ',SIZE(x_data,3)
        STOP
     END IF

     emismass(:,:,class(jc)%level,:) = emismass(:,:,class(jc)%level,:) &
          + x_data(:,:,:) * class(jc)%scale

  END DO class_loop

  ! (6) CONVERT MASS FLUX TO FLUX
  af_lon = 1.0_DP / REAL(nlon, DP)  ! area fraction (longitude)

  lat_loop: DO jj = 1, nlat
     ! area fraction (latitude)
     af_lat = (0.5_dp * ( sin(((xlat(jj)+dhlat)/180._DP)*pi)     &
                        - sin(((xlat(jj)-dhlat)/180._DP)*pi) ) )

     lon_loop: DO ji=1, nlon

        af     = af_lat * af_lon                ! area fraction
        af_sum = af_sum + af                    ! check sum
        area   = af * 4.0_DP * pi * r_earth**2  ! area in [m^2]

!qqq CHECK (1/ntime) WITH NEW EDGAR4 monthly, seasonal files
        conv   = ( N_A * 1000.0_DP / ( MOLARMASS * area  &
             * (REAL(spy, DP)/REAL(ntime, DP)) ) )

        emismass(ji,jj,:,:) = emismass(ji,jj,:,:) * GLOBALSCALE
        emisflux(ji,jj,:,:) = emismass(ji,jj,:,:) * conv
        
     END DO lon_loop
  END DO lat_loop

  ! (7) DIAGNOSTIC OUTPUT
  WRITE(*,*) '==========================================================='
  WRITE(*,*) 'OUTPUT         : ', TRIM(OUTPUT)
  WRITE(*,*) 'SPECIES        : ', TRIM(SPECIES)
  WRITE(*,*) 'MOLAR MASS     : ', MOLARMASS
  WRITE(*,*) 'GLOBAL SCALING : ', GLOBALSCALE
  WRITE(*,*) 'MASS FLUX      : ', L_MASSFLUX
  WRITE(*,*) 'YEAR           : ', YEAR
  WRITE(*,*) 'HEIGHT         : ', HEIGHT(1:nlev)
  WRITE(*,*) 'INPUTPATH      : ', TRIM(INPUTPATH)
  WRITE(*,*) 'SPECIES NAME   : ',TRIM(s_species)
  WRITE(*,*) 'UNIT           : ',TRIM(s_unit)
  WRITE(*,*) 'FREQUENCY      : ',TRIM(s_freq)
  WRITE(*,*)
  WRITE(*,*) 'MASS FLUX      : '
  WRITE(*,*) ' LBOUND        : ',LBOUND(emismass)
  WRITE(*,*) ' UBOUND        : ',UBOUND(emismass)
  WRITE(*,*) ' MIN           : ',MINVAL(emismass)
  WRITE(*,*) ' MAX           : ',MAXVAL(emismass)
  WRITE(*,*) ' SUM           : ',SUM(emismass)
  WRITE(*,*)
  WRITE(*,*) 'FLUX           : '
  WRITE(*,*) ' AREA CHECKSUM : ',af_sum
  WRITE(*,*) ' LBOUND        : ',LBOUND(emisflux)
  WRITE(*,*) ' UBOUND        : ',UBOUND(emisflux)
  WRITE(*,*) ' MIN           : ',MINVAL(emisflux)
  WRITE(*,*) ' MAX           : ',MAXVAL(emisflux)
  WRITE(*,*) '==========================================================='

  ! (8) OUTPUT TO NETCDF-FILE
  CALL nc_dump

  ! (9) CLEAN
  IF (ALLOCATED(emismass)) DEALLOCATE(emismass)
  IF (ALLOCATED(emisflux)) DEALLOCATE(emisflux)
  IF (ALLOCATED(xlon)) DEALLOCATE(xlon)
  IF (ALLOCATED(xlat)) DEALLOCATE(xlat)
  IF (ALLOCATED(xlev)) DEALLOCATE(xlev)
  IF (ALLOCATED(xtime)) DEALLOCATE(xtime)
  IF (ASSOCIATED(x_data))  DEALLOCATE(x_data)

! --------------------------------------------------------------------------
! ##########################################################################
! --------------------------------------------------------------------------
CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE USAGE(EXE)
    CHARACTER (LEN=*) :: EXE
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'EDGAR2NC Version ',VERSION
    WRITE(*,*) 'Author: Patrick Joeckel, MPICH, October 2004'
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'Usage: '//TRIM(EXE)//' <namelist-file>'
    WRITE(*,*) '--------------------------------------------'
  END SUBROUTINE USAGE
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE read_nml(status, iou, fname)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: iou
    CHARACTER(LEN=*), INTENT(IN)  :: fname

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_nml'
    LOGICAL :: lex       ! file exists
    INTEGER :: fstat = 0 ! file status

    NAMELIST /CTRL/ OUTPUT, SPECIES, MOLARMASS, GLOBALSCALE, L_MASSFLUX  & 
                  , YEAR, HEIGHT, INPUTPATH, CLASS
  
    status = 1 ! ERROR

    WRITE(*,*) '==========================================================='

    ! CHECK IF FILE EXISTS
    INQUIRE(file=TRIM(fname), exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) substr,': FILE DOES NOT EXIST (',TRIM(fname),')'
       status = 1
       RETURN
    END IF

    ! OPEN FILE
    OPEN(iou,file=TRIM(fname))

    ! READ NEMELIST
    WRITE(*,*) 'READING NAMELIST ''CTRL'''//&
         &' FROM '''//TRIM(fname),''' (unit ',iou,') ...'
    !
    READ(iou, NML=CTRL, IOSTAT=fstat)
    !
    IF (fstat /= 0) THEN
       WRITE(*,*) substr,': READ ERROR IN NAMELIST ''CTRL'' (',TRIM(fname),')'
       status = 3  ! READ ERROR IN NAMELIST
       RETURN
    END IF

    WRITE(*,*) ' OUTPUT   : ', TRIM(OUTPUT)
    WRITE(*,*) ' SPECIES  : ', TRIM(SPECIES)
    WRITE(*,*) ' YEAR     : ', YEAR
    WRITE(*,*) ' INPUTPATH: ', TRIM(INPUTPATH)

    ! CLOSE FILE
    CLOSE(iou)

    WRITE(*,*) '==========================================================='

    status = 0

  END SUBROUTINE read_nml
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE read_geia(status, iou, fname  &
       , species, year, freq, unit, data)

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL, NINT

    ! I/O
    INTEGER, INTENT(OUT)                   :: status
    INTEGER, INTENT(IN)                    :: iou     ! logical I/O unit
    CHARACTER(LEN=*),          INTENT(IN)  :: fname   ! filename
    CHARACTER(LEN=str_short),  INTENT(OUT) :: species
    INTEGER,                   INTENT(OUT) :: year
    CHARACTER(LEN=str_short),  INTENT(OUT) :: unit
    CHARACTER(LEN=str_short),  INTENT(OUT) :: freq
    REAL(DP), DIMENSION(:,:,:), POINTER    :: data    ! INTENT(OUT)
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_geia'
    LOGICAL                     :: lex    ! file exists
    LOGICAL                     :: lopn   ! file is already opened
    INTEGER                     :: iout   ! unit, if already opened
    INTEGER                     :: fstat  ! file status
    CHARACTER(LEN=str_vlong)    :: line   ! line
    INTEGER                     :: nl     ! line counter
    INTEGER                     :: nhl, ndl, nel  ! header, data, empty lines
    CHARACTER(LEN=str_vlong), DIMENSION(:), POINTER :: word => NULL()
    INTEGER                     :: nw     ! word counter
    INTEGER                     :: jw,jw2 ! loop counter for words
    REAL                        :: lon, lat
    REAL(DP), DIMENSION(:), POINTER :: value => NULL()
    INTEGER                     :: i
    INTEGER                     :: ilon, ilat
    LOGICAL                     :: ldl    ! line is data line; not header line
    INTEGER                     :: ntime  ! number of time-steps
    
    ! INITIALIZE
    species = ''
    year = 0
    unit = ''
    freq = ''
    ntime = 1   ! DEFAULT: annual
    nl  = 0     ! number of lines
    nhl = 0     ! number of header lines
    ndl = 0     ! number of data lines
    nel = 0     ! number of empty lines
    IF (ASSOCIATED(data)) DEALLOCATE(data)

    WRITE(*,*) '-----------------------------------------------------------'

    ! CHECK IF FILE EXISTS
    INQUIRE(file=TRIM(fname), exist=lex, opened=lopn, number=iout)
    !
    IF (.NOT.lex) THEN
       WRITE(*,*) substr,': FILE DOES NOT EXIST (',TRIM(fname),')'
       status = -1   ! FILE DOES NOT EXIST
       RETURN ! ERROR
    END IF
    IF (lopn) THEN
       WRITE(*,*) substr,': FILE ALREADY OPEN ON UNIT ',iout
       status = 2   ! FILE IS ALREADY OPEN
       RETURN ! ERROR
    END IF
    
    ! OPEN FILE
    OPEN(unit=iou,file=TRIM(fname))
    WRITE(*,*) 'FILE           : ',TRIM(fname)

    ! READ FILE
    line_loop: DO
       READ(iou,'(a)',IOSTAT=fstat) line
       IF (fstat /= 0) EXIT
       nl = nl + 1
       IF (TRIM(line) == '') THEN
          nel = nel + 1
          CYCLE
       END IF
       ! READ LINE
       ldl = is_numeric(line)
       data_line: IF (ldl) THEN
          ! DATA LINE
          ndl = ndl + 1
          CALL strcrack(line, ',' ,word, nw)
          IF (nw /= ntime + 2) THEN
             WRITE(*,*) substr,': READ ERROR IN DATA; LINE=',nl
             WRITE(*,*) word
             status = 4 ! READ ERROR IN DATA LINE
             RETURN             
          END IF
          IF (.NOT. ASSOCIATED(value)) THEN
             ALLOCATE(value(ntime))
          END IF
          value(:) = 0.0_DP
          IF (.NOT. ASSOCIATED(data)) THEN
             ALLOCATE(data(nlon,nlat,ntime))
             data(:,:,:) = 0.0_DP
          END IF
          !
          ! CRACK DATA LINE
          READ(word(1),*) lon
          READ(word(2),*) lat
          DO i=1, ntime
             READ(word(i+2),*) value(i)
          END DO
          ilon = NINT(lon + 180.0) + 1
          ilat = NINT(lat + 90.0)  + 1
          data(ilon, ilat, :) = data(ilon, ilat, :) + value(:)
       ELSE
          ! HEADER LINE
          nhl = nhl + 1
          CALL strcrack(line, ' ',word, nw)
          word_loop: DO jw=1, nw
             ! UNIT
             IF (word(jw) == 'Units') THEN
                DO jw2 = jw+2, nw
                   IF (word(jw2) == 'For') EXIT
                   unit = TRIM(unit)//' '//TRIM(word(jw2))
                END DO
                unit = TRIM(ADJUSTL(unit))
             END IF
             ! SPECIES, YEAR, AND FREQ
             IF (word(jw) == '#') THEN
                READ(word(jw-2),*) year
                !
                DO jw2=1, jw-3
                   species = TRIM(species)//' '//TRIM(word(jw2))
                END DO
                species = TRIM(ADJUSTL(species))
                !
                freq = TRIM(word(jw-1))
                SELECT CASE(TRIM(ADJUSTL(freq)))
                CASE('annual')
                   ntime = 1
                CASE('seasonal')
                   ntime = 4
                CASE('monthly')
                   ntime = 12
                CASE DEFAULT
                   WRITE(*,*) substr, &
                        ': UNKNOWN DATA FREQUENCY',TRIM(ADJUSTL(freq))
                   status = 5  ! UNKNOWN FREQUENCY
                   RETURN
                END SELECT
             END IF
          END DO word_loop
       END IF data_line
    END DO line_loop
    
    IF (fstat > 0) THEN
       WRITE(*,*) substr,': READ ERROR (IOSTAT=',fstat,') ON UNIT '&
            ,iou,' (FILE=',TRIM(fname),')'
       status = 3  ! READ ERROR
    ELSE
       ! END OF FILE REACHED
    END IF
    
    ! CLOSE FILE
    CLOSE(iou)

    WRITE(*,*) '# LINES        : ',nl
    WRITE(*,*) '# EMPTY LINES  : ',nel
    WRITE(*,*) '# HEADER-LINES : ',nhl
    WRITE(*,*) '# DATA-LINES   : ',ndl

    WRITE(*,*) 'SPECIES        : ',TRIM(species)
    WRITE(*,*) 'UNIT           : ',TRIM(unit)
    WRITE(*,*) 'FREQUENCY      : ',TRIM(freq)
    WRITE(*,*) 'YEAR           : ',year

    WRITE(*,*) 'MIN            : ',MINVAL(data)
    WRITE(*,*) 'MAX            : ',MAXVAL(data)
    WRITE(*,*) 'SUM            : ',SUM(data)
    WRITE(*,*) '-----------------------------------------------------------'

    ! CLEAN
    IF (ASSOCIATED(word)) DEALLOCATE(word)
    IF (ASSOCIATED(value)) DEALLOCATE(value)

    ! RETURN
    status = 0
    
  END SUBROUTINE read_geia
  ! ------------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE strcrack(str, ch, el, n)

    IMPLICIT NONE
    
    INTRINSIC :: INDEX, LEN_TRIM

    ! I/O
    CHARACTER(LEN=*),               INTENT(IN)  :: str
    CHARACTER,                      INTENT(IN)  :: ch
    CHARACTER(LEN=*), DIMENSION(:), POINTER     :: el       
    INTEGER,                        INTENT(OUT) :: n
    
    ! LOCAL
    INTEGER :: idx1, idx2
    
    ! INIT
    IF (ASSOCIATED(el)) DEALLOCATE(el)
    NULLIFY(el)
    n = 0
    
    ! EMPTY STRING
    IF ( (TRIM(str) == '') .OR. (TRIM(str) == ch) ) RETURN
    
    idx1 = 0
    idx2 = 0
    DO 
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE
       
       n = n + 1
       
    END DO
    
    ! ALLOCATE SPACE
    ALLOCATE(el(n))
    
    n = 0
    idx1 = 0
    idx2 = 0
    DO     
       idx1 = idx2 + 1
       IF (idx1 > LEN_TRIM(str(:))) EXIT
       IF (INDEX(TRIM(str(idx1:)), ch) == 0) THEN
          idx2 = LEN_TRIM(str(:)) + 1
       ELSE
          idx2 = idx2 + INDEX(TRIM(str(idx1:)), ch)
       END IF
       IF (idx1 == idx2) CYCLE
       
       n = n + 1
       
       el(n) = str(idx1:idx2-1)
       
    END DO

  END SUBROUTINE strcrack
  ! ---------------------------------------------------------------------  

  ! ------------------------------------------------------------------------
  FUNCTION is_numeric(string)

    IMPLICIT NONE

    INTRINSIC :: INDEX, LEN_TRIM

    ! I/O
    LOGICAL                      :: is_numeric
    CHARACTER(LEN=*), INTENT(IN) :: string

    ! LOCAL
    INTEGER                     :: n, i

    is_numeric = .true.

    n  = LEN_TRIM(string)

    IF (INDEX(string,'--') /= 0) THEN
       is_numeric = .false.
       RETURN
    END IF

    DO i=1, n
       SELECT CASE(string(i:i))
       CASE ('0','1','2','3','4','5','6','7','8','9')
       CASE ('E','e','.','-','+',',',' ')
       CASE DEFAULT
          is_numeric = .false.
          EXIT
       END SELECT
    END DO
    
  END FUNCTION is_numeric
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE nc_dump

    USE netcdf

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, CHAR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'nc_dump'
    INTEGER :: ncid      ! netCDF-ID
    INTEGER :: dimid_lat, dimid_lon, dimid_lev, dimid_time
    INTEGER :: varid_lat, varid_lon, varid_lev, varid_time
    INTEGER :: varid_mass, varid_flux, varid_height
    INTEGER                  :: jc
    CHARACTER(LEN=str_long)  :: timestr = ''
    CHARACTER(LEN=4)         :: yrstr = '    '
    CHARACTER(LEN=1000 + MAX_NCLASS*50) :: nmlstr = ''
    CHARACTER(LEN=4)         :: jcstr = '    '
    CHARACTER(LEN=4)         :: levstr = '    '
    CHARACTER(LEN=str_long)  :: scalestr = ''
    CHARACTER(LEN=str_long)  :: mmassstr = ''
    CHARACTER(LEN=str_long)  :: glscstr = ''
    CHARACTER(LEN=str_vlong) :: heightstr = ''
    !
    CHARACTER(LEN=8)         :: date
    CHARACTER(LEN=10)        :: time
    CHARACTER(LEN=5)         :: zone

    ! INIT
    ! CONVERT TO STRING
    WRITE(yrstr,'(i4)') YEAR
    CALL DATE_AND_TIME(date, time, zone)
    WRITE(mmassstr,*)  MOLARMASS
    WRITE(glscstr,*)   GLOBALSCALE
    WRITE(heightstr,*) HEIGHT(1:nlev)

    ! CREATE NEW FILE
    CALL NFERR( &
         nf90_create(TRIM(OUTPUT), NF90_CLOBBER, ncid) &
         ,1)

    ! ADD GLOBALE ATTRIBUTES
    ! - VERSION
    CALL NFERR( &
         nf90_put_att(ncid, NF90_GLOBAL, 'created_with',     &
         'edgar2nc Version '//TRIM(VERSION)) &
         ,2)
    ! - DATE AND TIME
    CALL NFERR( &
         nf90_put_att(ncid, NF90_GLOBAL, 'date', date) &
         ,3)
    CALL NFERR( &
         nf90_put_att(ncid, NF90_GLOBAL, 'time', TRIM(time)//TRIM(zone)) &
         ,4)

    ! - NAMELIST
    nmlstr = CHAR(10)//'OUTPUT= '''//TRIM(OUTPUT)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'SPECIES= '''//TRIM(SPECIES)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'MOLARMASS= '''//TRIM(mmassstr)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'GLOBALSCALE= '''//TRIM(glscstr)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'YEAR= '//yrstr//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'HEIGHT= '//TRIM(heightstr)//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'INPUTPATH= '''//TRIM(INPUTPATH)//''', '//CHAR(10)
    class_loop: DO jc = 1, MAX_NCLASS
     
       ! SKIP IF NAME IS EMPTY
       IF (TRIM(class(jc)%fname) == '') CYCLE

       ! SKIP IF FILE WAS NOT PRESENT
       IF (.NOT. lfile_ok(jc)) CYCLE

       WRITE(jcstr,'(i4)') jc
       WRITE(levstr,'(i4)') class(jc)%level
       WRITE(scalestr,*) class(jc)%scale

       nmlstr=TRIM(nmlstr)//'CLASS('//jcstr//') = '''&
            &//TRIM(class(jc)%fname)//''', '//levstr//', '//TRIM(scalestr)&
            &//', '//CHAR(10)

    END DO class_loop
    !
    WRITE(*,*) TRIM(nmlstr)
    !
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
    SELECT CASE(TRIM(ADJUSTL(s_freq)))
    CASE('annual')
       timestr = 'year since '//yrstr//'-01-01 00:00:00'
    CASE('seasonal')
       timestr = 'season since '//yrstr//'-01-01 00:00:00'
    CASE('monthly')
       timestr = 'month since '//yrstr//'-01-01 00:00:00'
    CASE DEFAULT
       WRITE(*,*) substr,': UNKNOWN DATA FREQUENCY ',TRIM(ADJUSTL(s_freq))
       STOP
    END SELECT
    CALL NFERR( &
         nf90_put_att(ncid, varid_time, 'units', timestr) &
         ,20)

    ! DEFINE VARIABLES
    ! - emission height
    CALL NFERR( &
         nf90_def_var(ncid, 'height', NF90_FLOAT  &
         , (/ dimid_lev /), varid_height) &
         ,21)
    CALL NFERR( &
         nf90_put_att(ncid, varid_height, 'long_name' &
         , 'emission height') &
         ,22)
    CALL NFERR( &
         nf90_put_att(ncid, varid_height, 'units', 'm') &
         ,23)

    ! - mass flux
    IF (L_MASSFLUX) THEN
       CALL NFERR( &
            nf90_def_var(ncid, TRIM(SPECIES)//'_mflux', NF90_FLOAT  &
            , (/ dimid_lon, dimid_lat, dimid_lev, dimid_time /), varid_mass) &
            ,24)
       CALL NFERR( &
            nf90_put_att(ncid, varid_mass, 'long_name' &
            , 'mass flux of '//TRIM(s_species)) &
            ,25)
       CALL NFERR( &
            nf90_put_att(ncid, varid_mass, 'units', TRIM(s_unit)) &
            ,26)
       CALL NFERR( &
            nf90_put_att(ncid, varid_mass, 'molar_mass', MOLARMASS) &
            ,27)
    END IF

    ! - flux
    CALL NFERR( &
         nf90_def_var(ncid, TRIM(SPECIES)//'_flux', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_lev, dimid_time /), varid_flux) &
         ,28)
    CALL NFERR( &
         nf90_put_att(ncid, varid_flux, 'long_name' &
         , 'flux of '//TRIM(s_species)) &
         ,29)
    CALL NFERR( &
         nf90_put_att(ncid, varid_flux, 'units', 'molecules m-2 s-1') &
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
    CALL NFERR( &
         nf90_put_var(ncid, varid_height, HEIGHT(1:nlev)) &
         ,37)

    ! SAVE VARIABLES
    IF (L_MASSFLUX) THEN
       CALL NFERR( &
            nf90_put_var(ncid, varid_mass, emismass) &
            ,38)
    END IF

    CALL NFERR( &
         nf90_put_var(ncid, varid_flux, emisflux) &
         ,39)

    ! CLOSE FILE
    CALL NFERR( &
         nf90_close(ncid) &
         ,40)

  END SUBROUTINE nc_dump
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

! --------------------------------------------------------------------------
! ##########################################################################
END PROGRAM EDGAR2NC
! ##########################################################################
! --------------------------------------------------------------------------
