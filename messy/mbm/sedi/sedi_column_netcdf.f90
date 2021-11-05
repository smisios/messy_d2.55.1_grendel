! -*- f90 -*-
!*******************************************************************************
!                Time-stamp: <2006-08-11 15:31:27 akerkweg>
!*******************************************************************************

! This file produces netcdf output. It can only be used if the netcdf
! library is available.

MODULE sedi_column_netcdf

  USE messy_main_constants_mem, ONLY: dp
  USE netcdf ! for nf90_* functions

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: open_nc_file, write_nc_file, close_nc_file

  ! choose a precision:
  INTEGER, PARAMETER :: PREC = nf90_float
  ! INTEGER, PARAMETER :: PREC = nf90_double

  INTEGER :: dimids2d(3), dimids3d(4)
  INTEGER :: londimid, latdimid, levdimid, tdimid
  INTEGER :: lonid, latid, levid, tid
  INTEGER :: start3d(4), cnt3d(4)
  INTEGER :: nlon, ngl, nlev
  INTEGER :: specid(1500) ! same as MAX_EQN in gdata.h

CONTAINS

  !*****************************************************************************

  SUBROUTINE nf(status) ! turns nf90_* function into subroutine + checks status
    INTEGER :: status
    IF (status /= nf90_noerr) THEN
      WRITE (*,*) 'netcdf error: ', nf90_strerror(status)
      STOP
    ENDIF
  END SUBROUTINE nf

  !*****************************************************************************

  SUBROUTINE open_nc_file (ncid, filename, &
       species, units, long_name, klon, klat, klev)

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: filename, species(:), units(:)
    CHARACTER(*), OPTIONAL, INTENT(IN) :: long_name(:)
    INTEGER, OPTIONAL, INTENT(IN) :: klev
    INTEGER, OPTIONAL, INTENT(IN) :: klon
    INTEGER, OPTIONAL, INTENT(IN) :: klat
    INTEGER, INTENT(OUT) :: ncid
    REAL, DIMENSION(:), ALLOCATABLE :: nlondata, ngldata, nlevdata

    INTEGER :: i

    ! define the grid size
    IF (PRESENT(klon)) THEN
       nlon=klon
    ELSE
       nlon = 1
    ENDIF
    IF (PRESENT(klat)) THEN
       ngl=klat
    ELSE
       ngl = 1
    ENDIF
    IF (PRESENT(klev)) THEN
       nlev=klev
    ELSE
       nlev = 1
    ENDIF

    ALLOCATE(nlondata(nlon))
    ALLOCATE(ngldata(ngl))
    ALLOCATE(nlevdata(nlev))

    ! define the grid
    nlondata(:) = 0. !(/ 0. /)
    ngldata(:)  = 0. !(/ 0. /)
    DO i=1, nlev
       nlevdata(i) = REAL(i)
    END DO
 
    ! open the netcdf file (nf90_clobber = overwrite existing file)
    CALL nf(nf90_create(filename//'.nc', nf90_clobber, ncid))

    ! global attributes
    CALL nf(nf90_put_att(ncid, nf90_global, 'title', 'sedi_column'))

    ! definition of the dimensions
    ! syntax: nf90_def_dim(IN:ncid, IN:name, IN:len, OUT:dimid)
    CALL nf(nf90_def_dim(ncid, 'lon',  nlon,           londimid))
    CALL nf(nf90_def_dim(ncid, 'lat',  ngl,            latdimid))
    CALL nf(nf90_def_dim(ncid, 'lev',  nlev,           levdimid))
    CALL nf(nf90_def_dim(ncid, 'time', nf90_unlimited, tdimid))

    dimids2d(:) = (/ londimid, latdimid,           tdimid /) ! 2d (3rd dim=t)
    dimids3d(:) = (/ londimid, latdimid, levdimid, tdimid /) ! 3d (4th dim=t)
    cnt3d(:)    = (/ nlon,     ngl,      nlev,     1 /)

    ! definition of variables
    ! syntax: nf90_def_var(IN:ncid, IN:name, IN:xtype, IN:dimids, OUT:varid)
    ! coordinate variables
    CALL nf(nf90_def_var(ncid, 'lon',   PREC,        londimid,  lonid))
    CALL nf(nf90_def_var(ncid, 'lat',   PREC,        latdimid,  latid))
    CALL nf(nf90_def_var(ncid, 'lev',   PREC,        levdimid,  levid))
    CALL nf(nf90_def_var(ncid, 'time',  nf90_double, tdimid,    tid))

    DO i = 1, SIZE(species)
      CALL nf(nf90_def_var(ncid, TRIM(species(i)), PREC, dimids3d, specid(i)))
      IF (PRESENT(long_name)) THEN
        CALL nf(nf90_put_att(ncid, specid(i), 'long_name', TRIM(long_name(i))))
      ELSE
        CALL nf(nf90_put_att(ncid, specid(i), 'long_name', TRIM(species(i))))
      ENDIF
      CALL nf(nf90_put_att(ncid, specid(i), 'units',     TRIM(units(i))))
    END DO

    ! assign attributes
    ! syntax: nf90_put_att(IN:ncid, IN:vid, IN:name, IN:values)
    ! longitude
    CALL nf(nf90_put_att(ncid, lonid,  'long_name', 'longitude'))
    CALL nf(nf90_put_att(ncid, lonid,  'units',     'degrees_east'))
    ! latitude
    CALL nf(nf90_put_att(ncid, latid,  'long_name', 'latitude'))
    CALL nf(nf90_put_att(ncid, latid,  'units',     'degrees_north'))
    ! levels
    CALL nf(nf90_put_att(ncid, levid,  'long_name', 'level index'))
    CALL nf(nf90_put_att(ncid, levid,  'units',     'level'))
    CALL nf(nf90_put_att(ncid, levid,  'positive',  'down'))
    ! time
    CALL nf(nf90_put_att(ncid, tid,    'long_name', 'time'))
    CALL nf(nf90_put_att(ncid, tid,    'units',     'seconds since 2000-01-01 00:00:00'))

    ! end of the definitions, switch to data mode
    CALL nf(nf90_enddef(ncid))

    ! syntax: nf90_put_var(IN:ncid, IN:varid, IN:values)
    ! write the data of the grid
    CALL nf(nf90_put_var(ncid, lonid,   nlondata))
    CALL nf(nf90_put_var(ncid, latid,   ngldata))
    CALL nf(nf90_put_var(ncid, levid,   nlevdata))

  END SUBROUTINE open_nc_file

  !*****************************************************************************

  SUBROUTINE write_nc_file (ncid, time, x, klev)

    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: ncid
    REAL(dp), INTENT(in) :: time
    REAL(dp), INTENT(in) :: x(:)
    INTEGER,  INTENT(IN) :: klev

    INTEGER :: i, timestep, limit, indu,indl

    ! write timestep
    CALL nf(nf90_inquire_dimension(ncid, tid, len=timestep))
    timestep = timestep + 1

    ! syntax: nf90_put_var(ncid, varid, values, start, cnt)
    ! start:  start in netcdf variable
    ! cnt:    number of netcdf variable points
    ! values: starting point of the fortran variable
    start3d = (/ 1,    1,   1,    timestep /)

    CALL nf(nf90_put_var(ncid, tid, time, (/timestep/) ))
    limit =  size(x)/klev 
    species_loop: DO I = 1, limit
       indu=i*klev
       indl=(i-1)*klev+1
!       write (*,*) 'write_nc_file', klev, limit, indl,indu,x(indl:indu)
       CALL nf(nf90_put_var(ncid, specid(i), (/x(indl:indu)/), start3d, cnt3d))
    END DO species_loop

    CALL nf(nf90_sync(ncid)) ! write buffer to file

  END SUBROUTINE write_nc_file

  !*****************************************************************************

  SUBROUTINE close_nc_file (ncid)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ncid

    CALL nf(nf90_close(ncid))

  END SUBROUTINE close_nc_file

  !*****************************************************************************

END MODULE sedi_column_netcdf

!*******************************************************************************
