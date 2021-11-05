  PROGRAM biogen

    ! ------------------------------------------------------------------------
    !     This program is mostly consistent with the routines being developed and
    !     used by Rolf Von Kuhlmann, MPIC-Mainz for the offline chemistry and
    !     tracer and transport model MATCH. The code produces some default set
    !     of emission files in the netCDF format at the highest resolution 
    !     available which can then be used in a model like echam5, to define
    !     the emission fluxes at the grid-resolution of the model using the
    !     regridding routines (2002-05).
    !
    !     The code has been extended such that other users can easily 
    !     manipulate the code in order to produce alternative emission fields
    !     for use in specific case studies.
    !
    !     Laurens Ganzeveld, MPIC-Mainz, 2002
    !   
    !     Final code develpment:
    !     Andrea Pozzer, MPIC-Mainz, 2006
    !     Patrick Joeckel, DLR, 2010: minor modifications
    ! ----------------------------------------------------------------------

    USE tools 
      
  IMPLICIT NONE


  ! variables are defined in tools.f90
  ! internal
  CHARACTER(LEN=80)  :: file_nml="biogen.nml"          ! argument
  INTEGER  :: status       ! status flag
  INTEGER :: jt,i    !counter

!--------------------------------------------------------------------------------
! 1) INIT
!--------------------------------------------------------------------------------

  DO jt=1, NMAXTRAC
    TRAC_NAME(jt) = ''
    TRAC_MOLAR_MASS(jt) = 0.0
    TRAC_LAND(jt)        = 0.0 
    TRAC_OCEAN(jt)      = 0.0 
    TRAC_YEAR(jt)       = 9999 
    L_NVOC(jt)          = .FALSE.
    L_TERP(jt)          = .FALSE.
    L_ISOP(jt)          = .FALSE.
  END DO
  
!--------------------------------------------------------------------------------
! 2) READ NAMELIST-FILE
!--------------------------------------------------------------------------------

  CALL read_nml(status, iou, TRIM(file_nml))
  IF (status /= 0) STOP
  
  write(*,*) hline
  write(*,*) "WARNINGS:"
  TRAC_TOT = 0 
  DO jt=1, NMAXTRAC
     IF (TRIM(TRAC_NAME(jt)) == '') CYCLE
     TRAC_TOT = TRAC_TOT + 1
     !
     TRAC_NAME(TRAC_TOT)       = TRIM(TRAC_NAME(jt))
     TRAC_MOLAR_MASS(TRAC_TOT) = TRAC_MOLAR_MASS(jt)
     TRAC_LAND(TRAC_TOT)        = TRAC_LAND(jt)
     TRAC_OCEAN(TRAC_TOT)      = TRAC_OCEAN(jt)
     TRAC_YEAR(TRAC_TOT)       = TRAC_YEAR(jt)
     L_NVOC(TRAC_TOT)          = L_NVOC(jt)
     L_TERP(TRAC_TOT)          = L_TERP(jt)
     L_ISOP(TRAC_TOT)          = L_ISOP(jt)
     ! WARNINGS
     IF (TRAC_LAND(jt).eq.0.) write(*,*)       " 0.0 LAND  (NATURAL) EMISSION FOR : ", TRAC_NAME(jt) 
     IF (TRAC_OCEAN(jt).eq.0.) write(*,*)      " 0.0 OCEAN (NATURAL) EMISSION FOR : ", TRAC_NAME(jt) 
     IF (TRAC_MOLAR_MASS(jt).eq.0.) write(*,*) " 0.0 MOLARMASS         FOR :        ", TRAC_NAME(jt) 
     IF (TRAC_YEAR(jt).eq.9999) write(*,*)     " DEFAULT TRAC_YEAR     FOR :        ", TRAC_NAME(jt) 
  END DO
  write(*,*) hline
  write (*,*) " Total number of tracer  = " , trac_tot 
  write(*,*) hline
  write(*,*) '    OUTPUT  :  '
  write(*,*) '    LON       = ',nlon
  write(*,*) '    LAT       = ',nlat 
  write(*,*) '    LEV       = ',nlev
  write(*,*) '    DIRECTORY = ',output_dir     
  
!--------------------------------------------------------------------------------
! 3) INPUT FILES AND SOME CALCULATION 
!--------------------------------------------------------------------------------

 ! Olson vegetation input needed for landmask (land) 
 CALL read_olson(status, iou, TRIM(INPUTPATH_OLSON), olson, land)  ! reading Olson data
 ! grid calculation
 call areacal(nlon,nlat,surfarea)

!--------------------------------------------------------------------------------
! 4) PREPARATION OUTPUT
!--------------------------------------------------------------------------------

    ALLOCATE (xlon(nlon), xlat(nlat), xlev(nlev))! coordinates
    ALLOCATE(xtime(ntime))
    ALLOCATE(output_field(nlon, nlat, nlev, ntime))
    ALLOCATE(ocean_field(nlon, nlat, ntime))
    ALLOCATE(land_field(nlon, nlat, ntime))
    output_field(:,:,:,:) = 0.0
    ! define the grid
    xlon = (/ (-179.5 + (360. * REAL(i-1) / REAL(nlon)),i=1,nlon) /) !
    xlat = (/ (- 89.5 + (180. * REAL(i-1) / REAL(nlat)),i=1,nlat) /)   ! 
    xlev = (/(i,i=1,nlev)/) 
    xtime = (/ (REAL(i-1),i=1,ntime) /)

!--------------------------------------------------------------------------------
! 5) OUTPUT
!--------------------------------------------------------------------------------

  IF (l_output_olson) THEN
     MOLARMASS=1.
     TOTAL=1.
     YEAR=2000
     l_flux= .FALSE.
     output_field(:,:,1,:)=olson(:,:,:)
     OUTPUT=TRIM(output_dir)//'/Olson_biome_'
     SPECIES=TRIM('biome')
     UNIT_OUTPUT='biome index'
     LONG_NAME_OUTPUT = 'biome'
     CALL nc_creation
  END IF

     !--------------------------------------------------------
     !                  EMISSIONS FIELDS 
     !---------------------------------------------------------
  DO jt=1,TRAC_TOT
     !--------------------------------------------------------
     !                       GENERAL
     !---------------------------------------------------------
     l_flux= .TRUE.
     output_field(:,:,:,:)=0.
     ocean_field(:,:,:)=0.
     land_field(:,:,:)=0.
     ! specific for tracer
     MOLARMASS=TRAC_MOLAR_MASS(jt)
     YEAR = TRAC_YEAR(jt)
     OUTPUT=TRIM(output_dir)//'/'//'biogen_'//TRIM(TRAC_NAME(jt))
     SPECIES=TRIM(TRAC_NAME(jt))
     LONG_NAME_OUTPUT='flux of '//TRIM(SPECIES)
     UNIT_OUTPUT='molecules m-2 s-1'
     lterp = L_TERP(jt)
     lisop = L_ISOP(jt)
     lnvoc = L_NVOC(jt)
     sland = TRAC_LAND(jt)
     socean = TRAC_OCEAN(jt)
     !--------------------------------------------------------
     !                       LAND
     !---------------------------------------------------------
     ! READ BASIC DISTRIBUTION land and scale it!
     call land_distr(iou,l_terp(jt),l_nvoc(jt),l_isop(jt),base_dist)
     ! land sea mask 
     CALL CALCLAND(land,landmask,1) ! ocean=0 land=1 
     ! mask the distribution
      call mask_distr(base_dist,landmask) ! not done in the original code.... but better
     ! scale the distribution
     call add_distr (TRAC_LAND(jt), base_dist, land_field(:,:,:))
     !--------------------------------------------------------
     !                       OCEAN
     !---------------------------------------------------------
     ! READ BASIC DISTRIBUTION ocean and scale it!
     IF (l_ocean_alt) THEN
       call land_distr(iou,.FALSE.,.TRUE.,.FALSE.,base_dist)
     ELSE
       call ocean_distr(base_dist)
     ENDIF
     ! calculate mask for ocean
     ! land sea mask 
     CALL CALCLAND(land,landmask,0) ! ocean=1 land=0 
     ! mask the distribution
     call mask_distr(base_dist,landmask)
     ! scale the distribution
     call add_distr (TRAC_OCEAN(jt), base_dist, ocean_field(:,:,:))
     !--------------------------------------------------------
     !                       SUM
     !---------------------------------------------------------
     ! Add the two distribution
     ! in g s-1 m-2 !! we multiply with (N_A/MOLARMASS) --> transform in mcl s-1 m-2
     output_field(:,:,1,:)=(ocean_field(:,:,:)+land_field(:,:,:))*(N_A/MOLARMASS)
     !--------------------------------------------------------
     !             netcdf files creation
     !---------------------------------------------------------
     CALL nc_creation
  END DO

!--------------------------------------------------------------------------------
! 6) CLEAN MEMORY
!--------------------------------------------------------------------------------

  DEALLOCATE (xlon, xlat, xlev)
  DEALLOCATE (xtime)
  DEALLOCATE (output_field)
  DEALLOCATE (land_field)
  DEALLOCATE (ocean_field)

! --------------------------------------------------------------------------
! ##########################################################################
! --------------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE read_nml(status, iou, fname)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: iou
    CHARACTER(LEN=*), INTENT(IN)  :: fname

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_nml'
    LOGICAL :: lex   ! file exists
    INTEGER :: fstat ! file status
    INTEGER :: jt 

    NAMELIST /CTRL/ INPUTPATH_OLSON,INPUTPATH_GEIA, output_dir, l_output_olson, &
                    l_ocean_alt,                                      &
                    TRAC_NAME,TRAC_MOLAR_MASS,TRAC_LAND,TRAC_OCEAN, TRAC_YEAR,   &
                    L_TERP, L_NVOC, L_ISOP
  
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

    WRITE(*,*) ' INPUTPATH_OLSON: ', TRIM(INPUTPATH_OLSON)
    WRITE(*,*) ' INPUTPATH_GEIA: ', TRIM(INPUTPATH_GEIA)

    IF (l_ocean_alt) THEN
       WRITE(*,*) ' ALTERNATIVE OCEANIC EMISSION DISTRIBUTION : ON'
    ELSE
       WRITE(*,*) ' ALTERNATIVE OCEANIC EMISSION DISTRIBUTION : OFF'
    END IF
    
    DO jt=1, NMAXTRAC
       IF (TRIM(TRAC_NAME(jt)) == '') CYCLE
       write (*,*) hline 
       write (*,*) "         TRACER DEFINED:       "
       write (*,*) "                               "
       write (*,*) "         TRAC_NAME       = ", TRIM(TRAC_NAME(jt)) 
       write (*,*) "         TRAC_MOLAR_MASS = ", TRAC_MOLAR_MASS(jt)
       write (*,*) "         TRAC_LAND        = ", TRAC_LAND(jt) 
       write (*,*) "         TRAC_OCEAN      = ", TRAC_OCEAN(jt) 
       write (*,*) "         TRAC_YEAR       = ", TRAC_YEAR(jt) 
       write (*,*) "          L_TERP         = ", L_TERP(jt) 
       write (*,*) "          L_ISOP         = ", L_ISOP(jt) 
       write (*,*) "          L_NVOC         = ", L_NVOC(jt) 
       write (*,*) hline 
    END DO


    ! CLOSE FILE
    CLOSE(iou)
    status = 0

  END SUBROUTINE read_nml
  ! ------------------------------------------------------------------------


  END PROGRAM biogen 
