!> module mo_wncdf
!>
!> contains methods to write simple netCDF files for a single array
MODULE mo_wncdf
  USE mo_kind, ONLY: dp!, sp
  USE mo_commo1, ONLY: alat_g, alon_g, alatpsi_g, alonpsi_g
!  USE mo_units, ONLY: io_stderr
!  use mo_parallel, only: p_abort
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  PRIVATE
  PUBLIC ::  write_ncdf_single_array
CONTAINS
  SUBROUTINE write_ncdf_single_array(dump_fname, a, a_shape, a_name, &
       a_valid_range, a_fill_value, lat_lon_geometry)
    CHARACTER(len=*), intent(in) :: dump_fname
    REAL(dp), INTENT(in) :: a(*)
    INTEGER, INTENT(in) :: a_shape(:)
    CHARACTER(len=*), INTENT(in) :: a_name
    REAL(dp), OPTIONAL, INTENT(in) :: a_valid_range(2), a_fill_value
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: lat_lon_geometry
    ! FIXME: optional lat, lon args


    CHARACTER(len=3), PARAMETER :: dim_name(7) = (/ 'x  ', 'y  ', 'z  ', &
         'xx ', 'yy ', 'zz ', 'xxx' /)
    INTEGER :: ndims
    INTEGER :: ncid, varid_a, varid_lat, varid_lon
    INTEGER :: i
    INTEGER, ALLOCATABLE :: a_dimids(:)
    INTEGER, PARAMETER :: ntrials = 2
    CHARACTER(*), PARAMETER :: title = 'MPIOM graph decomposition diagnostics'
    CHARACTER(*), PARAMETER :: institution = &
         'Deutsches Klimarechenzentrum GmbH'
    CHARACTER(*), PARAMETER :: source = &
         'decomposition computations'
    CHARACTER(*), PARAMETER :: longitude = 'longitude', latitude = 'latitude', &
         degrees_east = 'degrees_east', degrees_north = 'degrees_north', &
         lat_lon_coord = 'lat lon'

    ndims = SIZE(a_shape)
    ALLOCATE(a_dimids(ndims))

    CALL handle_ncdf_err(nf_create(dump_fname, NF_CLOBBER, ncid), &
         __LINE__)
    DO i = 1, ndims
      CALL handle_ncdf_err(nf_def_dim(ncid, TRIM(dim_name(i)), a_shape(i), &
           a_dimids(i)), __LINE__)
    END DO
    CALL handle_ncdf_err(nf_def_var(ncid, a_name, NF_DOUBLE, ndims, &
         a_dimids, varid_a), __LINE__)
!      CALL handle_ncdf_err(nf_put_att_text(ncid, varid_doms, &
!           'long_name', LEN(doms_lname), doms_lname), __LINE__)



    ! write global cf information
    CALL handle_ncdf_err(nf_put_att_text(ncid, NF_GLOBAL, &
           'title', LEN(title), title), __LINE__)
    CALL handle_ncdf_err(nf_put_att_text(ncid, NF_GLOBAL, &
         'institution', LEN(institution), institution), __LINE__)
    CALL handle_ncdf_err(nf_put_att_text(ncid, NF_GLOBAL, &
         'source', LEN(source), source), __LINE__)
    IF (PRESENT(a_valid_range)) THEN
      CALL handle_ncdf_err(nf_put_att_double(ncid, varid_a, &
         'valid_range', NF_DOUBLE, 2, a_valid_range(1)), __LINE__)
    END IF
    IF (PRESENT(a_fill_value)) THEN
      CALL handle_ncdf_err(nf_put_att_double(ncid, varid_a, &
           '_FillValue', NF_DOUBLE, 1, a_fill_value), __LINE__)
    END IF
    IF (PRESENT(lat_lon_geometry)) THEN
      SELECT CASE (lat_lon_geometry)
      CASE ('s')
        ! FIXME: check for the following
        !  - ndims == 2
        !  - a_shape(1) == ie_g
        !  - a_shape(2) == je_g
        CALL handle_ncdf_err(nf_put_att_text(ncid, varid_a, &
             'coordinates', LEN(lat_lon_coord), lat_lon_coord), __LINE__)
        CALL handle_ncdf_err(nf_def_var(ncid, 'lat', NF_DOUBLE, ndims, &
             a_dimids, varid_lat), __LINE__)
        CALL handle_ncdf_err(nf_def_var(ncid, 'lon', NF_DOUBLE, ndims, &
             a_dimids, varid_lon), __LINE__)
      CASE ('psi')
        ! FIXME: check for the following
        !  - ndims == 2
        !  - a_shape(1) == ie_g
        !  - a_shape(2) == je_g + 1
        CALL handle_ncdf_err(nf_put_att_text(ncid, varid_a, &
             'coordinates', LEN(lat_lon_coord), lat_lon_coord), __LINE__)
        CALL handle_ncdf_err(nf_def_var(ncid, 'lat', NF_DOUBLE, ndims, &
             a_dimids, varid_lat), __LINE__)
        CALL handle_ncdf_err(nf_def_var(ncid, 'lon', NF_DOUBLE, ndims, &
             a_dimids, varid_lon), __LINE__)
      END SELECT
    END IF
    CALL handle_ncdf_err(nf_enddef(ncid), __LINE__)

    CALL handle_ncdf_err(nf_put_var_double(ncid, varid_a, a), &
         __LINE__)
    IF (PRESENT(lat_lon_geometry)) THEN
      SELECT CASE (lat_lon_geometry)
      CASE ('s')
        CALL handle_ncdf_err(nf_put_var_double(ncid, varid_lat, alat_g), &
             __LINE__)
        CALL handle_ncdf_err(nf_put_var_double(ncid, varid_lon, alon_g), &
             __LINE__)
      CASE ('psi')
        CALL handle_ncdf_err(nf_put_var_double(ncid, varid_lat, alatpsi_g), &
             __LINE__)
        CALL handle_ncdf_err(nf_put_var_double(ncid, varid_lon, alonpsi_g), &
             __LINE__)
      END SELECT
    END IF

    CALL handle_ncdf_err(nf_close(ncid), __LINE__)

  END SUBROUTINE write_ncdf_single_array

  SUBROUTINE handle_ncdf_err(errcode, line)
    INTEGER, INTENT(in) :: errcode, line
    IF (errcode .NE. nf_noerr) THEN
      WRITE (0, *) 'Error at line ', line, ': ', nf_strerror(errcode)
!      CALL p_abort
       STOP 1
    END IF
  END SUBROUTINE handle_ncdf_err

END MODULE mo_wncdf
