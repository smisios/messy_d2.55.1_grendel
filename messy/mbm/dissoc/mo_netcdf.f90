!*******************************************************************************
!                Time-stamp: <2018-05-07 20:29:39 sander>
!*******************************************************************************

MODULE mo_netcdf

  USE netcdf ! for nf90_* functions
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: open_dissoc_nc_file, write_dissoc_nc_file, close_dissoc_nc_file

  ! choose a precision:
  INTEGER, PARAMETER :: PREC = nf90_float
  ! INTEGER, PARAMETER :: PREC = nf90_double

  INTEGER :: dimids2d(2)
  INTEGER :: szadimid, pressdimid
  INTEGER :: szaid, levid, pressid
  INTEGER :: specid

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

  SUBROUTINE open_dissoc_nc_file (ncid, nlev, press, nsza, sza)

    IMPLICIT NONE

    INTEGER, INTENT(OUT)              :: ncid
    INTEGER, INTENT(IN)               :: nlev, nsza
    REAL,    INTENT(IN), DIMENSION(:) :: sza, press

    INTEGER :: i
    REAL :: nlevdata(nlev)

    DO i = 1, nlev
      nlevdata(i) = REAL(i)
    END DO

    ! open the netcdf file (nf90_clobber = overwrite existing file)
    CALL nf(nf90_create('dissoc.nc', nf90_clobber, ncid))

    ! global attributes
    CALL nf(nf90_put_att(ncid, nf90_global, 'title', 'dissoc-box'))

    ! definition of the dimensions
    ! syntax: nf90_def_dim(IN:ncid, IN:name, IN:len, OUT:dimid)
    CALL nf(nf90_def_dim(ncid, 'sza',  nsza, szadimid))
    CALL nf(nf90_def_dim(ncid, 'press',  nlev, pressdimid))

    dimids2d(:) = (/ szadimid, pressdimid /)

    ! definition of variables
    ! syntax: nf90_def_var(IN:ncid, IN:name, IN:xtype, IN:dimids, OUT:varid)
    ! coordinate variables
    CALL nf(nf90_def_var(ncid, 'sza',   PREC, szadimid,  szaid))
    !CALL nf(nf90_def_var(ncid, 'lev',   PREC, levdimid,  levid))
    CALL nf(nf90_def_var(ncid, 'press',   PREC, pressdimid,  pressid))

    ! assign attributes
    ! syntax: nf90_put_att(IN:ncid, IN:vid, IN:name, IN:values)
    ! sza
    CALL nf(nf90_put_att(ncid, szaid,  'long_name', 'solar zenith angle'))
    CALL nf(nf90_put_att(ncid, szaid,  'units',     'degrees'))
    ! levels
    !CALL nf(nf90_put_att(ncid, levid,  'long_name', 'level index'))
    !CALL nf(nf90_put_att(ncid, levid,  'units',     'level'))
    !CALL nf(nf90_put_att(ncid, levid,  'positive',  'down'))

    CALL nf(nf90_put_att(ncid, pressid,  'long_name', 'Pressure'))
    CALL nf(nf90_put_att(ncid, pressid,  'units',     'hPa'))
    CALL nf(nf90_put_att(ncid, pressid,  'positive',  'down'))

    ! end of the definitions, switch to data mode
    CALL nf(nf90_enddef(ncid))

    ! syntax: nf90_put_var(IN:ncid, IN:varid, IN:values)
    ! write the data of the grid
    CALL nf(nf90_put_var(ncid, szaid, sza))
    !CALL nf(nf90_put_var(ncid, levid, nlevdata))
    CALL nf(nf90_put_var(ncid, pressid, press))

  END SUBROUTINE open_dissoc_nc_file

  !*****************************************************************************

  SUBROUTINE write_dissoc_nc_file (ncid, species, x)

    USE messy_main_constants_mem, ONLY: DP
    IMPLICIT NONE

    INTEGER,                      INTENT(IN) :: ncid
    CHARACTER(*),                 INTENT(IN) :: species
    REAL(DP),     DIMENSION(:,:), INTENT(in) :: x

    ! back to definition mode
    CALL nf(nf90_redef(ncid))

    CALL nf(nf90_def_var(ncid, TRIM(species), PREC, dimids2d, specid))
    CALL nf(nf90_put_att(ncid, specid, 'long_name', TRIM(species)))
    CALL nf(nf90_put_att(ncid, specid, 'units',     '1/s'))

    ! end of the definitions, switch to data mode
    CALL nf(nf90_enddef(ncid))

    ! syntax: nf90_put_var(IN:ncid, IN:varid, IN:values)
    CALL nf(nf90_put_var(ncid, specid, x))

    CALL nf(nf90_sync(ncid)) ! write buffer to file

  END SUBROUTINE write_dissoc_nc_file

  !*****************************************************************************

  SUBROUTINE close_dissoc_nc_file (ncid)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ncid
    CALL nf(nf90_close(ncid))

  END SUBROUTINE close_dissoc_nc_file

  !*****************************************************************************

END MODULE mo_netcdf

!*******************************************************************************
