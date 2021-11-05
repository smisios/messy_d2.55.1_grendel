PROGRAM anta2nc

  ! syntax:
  !   anta2nc < inputlist
  ! compile on LINUX:
  !   lf95 -g --chk asue -I/client/include -o anta2nc anta2nc.f90 -L/client/lib -lnetcdf
  !    f95 -C=all -g -gline -I/client/include -o anta2nc anta2nc.f90 -L/client/lib -lnetcdf
  ! compile on cross:
  ! efc anta2nc.f90 -I/pool/ia64/netcdf/netcdf-3.6.0-p1/include -L/pool/ia64/netcdf/netcdf-3.6.0-p1/lib -lnetcdf
  !
  ! syntax for bigendian extra file compile lf95 on LINUX
  !   anta2nc -Wl,-T < inputlist

  ! Convert MPIOM1 grid information from EXTRA file anta.ext
  ! to SCRIP compatible netCDF format.

  ! Most of the grid conversion is from mo_couple:grids_writing
  USE iso_varying_string
  USE mo_kind, ONLY: i8, sp
  IMPLICIT NONE

  INTEGER :: ie, je

  INTEGER, PARAMETER :: io_in_anta=14
  TYPE(varying_string) :: ifilename, ofilename

  LOGICAL :: lex
  INTEGER(i8) :: ext_hdr(4)
  INTEGER :: ipoints, flen, ioerrstat, verbose

  INTEGER       ::   i,j, ii, jj         ! looping indicees
  INTEGER       ::   ip1,im1             ! i+1, i-1
  REAL          ::   pi                  ! PI
  REAL          ::   rad2deg             ! 180/PI
  INTEGER, ALLOCATABLE ::   mask(:,:)    ! (ie-2,je) inverse land sea mask
  REAL, ALLOCATABLE    ::   lon(:,:)     ! (2*ie,2*je) longitudes of doubled array
  REAL, ALLOCATABLE    ::   lat(:,:)     ! (2*ie,2*je) latitudes of doubled array
  REAL(sp), ALLOCATABLE    ::   lons(:,:)    ! (ie-2,je) longitudes of scalars
  REAL(sp), ALLOCATABLE    ::   lats(:,:)    ! (ie-2,je) latitudes of scalars
  REAL(sp), ALLOCATABLE    ::   lonu(:,:)    ! (ie-2,je) longitudes of vector u-components
  REAL(sp), ALLOCATABLE    ::   latu(:,:)    ! (ie-2,je) latitudes of vector u-components
  REAL(sp), ALLOCATABLE    ::   lonv(:,:)    ! (ie-2,je) longitudes of vector v-components
  REAL(sp), ALLOCATABLE    ::   latv(:,:)    ! (ie-2,je) latitudes of vector v-components
  REAL(sp), ALLOCATABLE    ::   lonc(:,:)    ! (ie-2,je) corner grid points of scalar grid
  REAL(sp), ALLOCATABLE    ::   latc(:,:)    ! (ie-2,je) corner grid points of scalar grid
  REAL(sp), ALLOCATABLE    ::   clons(:,:,:) ! (4,ie-2,je) corner longitudes of scalars
  REAL(sp), ALLOCATABLE    ::   clats(:,:,:) ! (4,ie-2,je) corner latitudes of scalars
  REAL(sp), ALLOCATABLE    ::   clonu(:,:,:) ! (4,ie-2,je) corner longitudes of vector u-components
  REAL(sp), ALLOCATABLE    ::   clatu(:,:,:) ! (4,ie-2,je) corner latitudes of vector u-components
  REAL(sp), ALLOCATABLE    ::   clonv(:,:,:) ! (4,ie-2,je) corner longitudes of vector v-components
  REAL(sp), ALLOCATABLE    ::   clatv(:,:,:) ! (4,ie-2,je) corner latitudes of vector v-components
  NAMELIST /gridparams/ ie, je, verbose
  verbose = 0
  ie = -1
  je = -1
  READ (*, nml=gridparams)

  WRITE(*,*) 'Enter name of the input file (ANTA file in EXTRA format):'
  CALL get(ifilename, iostat=ioerrstat)
  IF (ioerrstat .gt. 0) CALL err_exit


  INQUIRE (FILE=CHAR(ifilename), EXIST=lex)
  IF (.NOT. lex) THEN
    WRITE(0,*) 'Could not open file <', CHAR(ifilename), '>'
    STOP 'Open failed'
  END IF

  WRITE(*,*) 'Enter prefix of the output file (e.g. grob3):'
  CALL get(ofilename, iostat=ioerrstat)
  IF (ioerrstat .gt. 0) CALL err_exit

  OPEN(io_in_anta, FILE=CHAR(ifilename), FORM='UNFORMATTED')

  READ(io_in_anta) ext_hdr

  ipoints = INT(ext_hdr(4)) / 4

  WRITE(*,*) 'The input file has ', ipoints, ' gridpoints.'

  IF ( ie*je .NE. ipoints ) THEN
    STOP 'Number of longitudes and latitudes does not match input!'
  END IF

  WRITE(*,*) 'o.k.'

  ALLOCATE(mask(ie-2,je))
  ALLOCATE(lon(2*ie,2*je))
  ALLOCATE(lat(2*ie,2*je))
  ALLOCATE(lons(ie-2,je))
  ALLOCATE(lats(ie-2,je))
  ALLOCATE(lonu(ie-2,je))
  ALLOCATE(latu(ie-2,je))
  ALLOCATE(lonv(ie-2,je))
  ALLOCATE(latv(ie-2,je))
  ALLOCATE(lonc(ie-2,je))
  ALLOCATE(latc(ie-2,je))
  ALLOCATE(clons(4,ie-2,je))
  ALLOCATE(clats(4,ie-2,je))
  ALLOCATE(clonu(4,ie-2,je))
  ALLOCATE(clatu(4,ie-2,je))
  ALLOCATE(clonv(4,ie-2,je))
  ALLOCATE(clatv(4,ie-2,je))

  READ(io_in_anta) lon
  READ(io_in_anta) ext_hdr
  READ(io_in_anta) lat

  CLOSE (io_in_anta)

  !
  !-- create arrays of longitudes and latitudes
  !
  !--  follow the OASIS conventions: S -> N
  !
!!  DO j = 1, 2*je
!!    lon(:,j)=gila(:,2*je+1-j)
!!    lat(:,j)=giph(:,2*je+1-j)
!!  ENDDO
  !
  !--  convert from radiant to degree
  !
  pi = ATAN(1.)*4.0
  rad2deg = 180./pi
  lon(:,:) = lon(:,:) * rad2deg
  lat(:,:) = lat(:,:) * rad2deg

  WHERE (lon(:,:) < 0.)
    lon(:,:) = lon(:,:) + 360.
  END WHERE

  IF (verbose > 0) THEN
    DO J = 1, 6
      WRITE(*,*) 'lon: ', j, lon(1:6, j)
    ENDDO
    DO J = 2*je-6, 2*je
      WRITE(*,*) 'lon: ', j, lon(1:6, j)
    ENDDO
    DO J = 1, 6
      WRITE(*,*) 'lat: ', j, lat(1:6, j)
    ENDDO
    DO J = 2*je-6, 2*je
      WRITE(*,*) 'lat: ', j, lat(1:6, j)
    ENDDO
  END IF
  !
  !--  extract scalar/vector grid points
  !
  !          1  2  3  4  5  6  ... 2*ie
  !
  !      1   c  v  c  v  c  v
  !      2   u  s  u  s  u  s
  !      3   c  v  c  v  c  v
  !      4   u  s  u  s  u  s
  !      5   c  v  c  v  c  v           c: grid cell corners of scalars
  !      6   u  s  u  s  u  s           v: vector-v
  !      :                              u: vector-u
  !      :                              s: scalar
  !     2*ij
  !
  DO i = 1, ie-2
    ii = (i*2)+2
    DO j = 1, je
      jj = j*2
      !--     scalar
      lats(i,j) = real(lat(ii,jj), sp)
      lons(i,j) = real(lon(ii,jj), sp)
      !--     vector - u
      latu(i,j) = real(lat(ii+1,jj), sp)
      lonu(i,j) = real(lon(ii+1,jj), sp)
    ENDDO
    !--     vector - v
    DO j = 1, je-1
      jj = j*2
      latv(i,j) = real(lat(ii,jj+1), sp)
      lonv(i,j) = real(lon(ii,jj+1), sp)
    ENDDO
    latv(i,je) = real(lat(ii,je*2), sp)
    lonv(i,je) = real(lon(ii,je*2), sp)
    !--     corners of scalar grid cells
    latc(i,1) = real(lat(ii+1,1), sp)
    lonc(i,1) = real(lon(ii+1,1), sp)
    DO j = 2, je
      jj = j*2
      latc(i,j) = real(lat(ii+1,jj-1), sp)
      lonc(i,j) = real(lon(ii+1,jj-1), sp)
    ENDDO
  ENDDO
  !
  !--  create corner arrays for SCRIP interpolation
  !
  DO i = 1, ie-2
    ii = (i*2)+2
    im1 = i-1
    ip1 = i+1
    IF (im1 == 0) im1 = ie-2
    IF (ip1 == ie-1) ip1 = 1
    !
    !--     scalar
    !
    DO j = 1, je-1
      clons(1,i,j) = lonc(im1,j  )
      clons(2,i,j) = lonc(im1,j+1)
      clons(3,i,j) = lonc(i  ,j+1)
      clons(4,i,j) = lonc(i  ,j  )
      clats(1,i,j) = latc(im1,j  )
      clats(2,i,j) = latc(im1,j+1)
      clats(3,i,j) = latc(i  ,j+1)
      clats(4,i,j) = latc(i  ,j  )
    ENDDO
    clons(1,i,je) = lonc(im1,je)
    clons(2,i,je) = lonu(im1,je)
    clons(3,i,je) = lonu(i  ,je)
    clons(4,i,je) = lonc(i  ,je)
    clats(1,i,je) = latc(im1,je)
    clats(2,i,je) = latu(im1,je)
    clats(3,i,je) = latu(i  ,je)
    clats(4,i,je) = latc(i  ,je)
    !
    !--     vector - u
    !
    DO j = 2, je
      clonu(1,i,j) = lonv(i  ,j-1)
      clonu(2,i,j) = lonv(i  ,j  )
      clonu(3,i,j) = lonv(ip1,j  )
      clonu(4,i,j) = lonv(ip1,j-1)
      clatu(1,i,j) = latv(i  ,j-1)
      clatu(2,i,j) = latv(i  ,j  )
      clatu(3,i,j) = latv(ip1,j  )
      clatu(4,i,j) = latv(ip1,j-1)
    ENDDO
    clonu(1,i,1) = real(lon(ii  ,1), sp)
    clonu(2,i,1) = lonv(i  ,1)
    clonu(3,i,1) = lonv(ip1,1)
    clonu(4,i,1) = real(lon(ii+2,1), sp)
    clatu(1,i,1) = real(lat(ii  ,1), sp)
    clatu(2,i,1) = latv(i  ,1)
    clatu(3,i,1) = latv(ip1,1)
    clatu(4,i,1) = real(lat(ii+2,1), sp)
    !
    !--     vector - v
    !
    DO j = 1, je-1
      clonv(1,i,j) = lonu(im1,j  )
      clonv(2,i,j) = lonu(im1,j+1)
      clonv(3,i,j) = lonu(i  ,j+1)
      clonv(4,i,j) = lonu(i  ,j  )
      clatv(1,i,j) = latu(im1,j  )
      clatv(2,i,j) = latu(im1,j+1)
      clatv(3,i,j) = latu(i  ,j+1)
      clatv(4,i,j) = latu(i  ,j  )
    ENDDO
    clonv(1,i,je) = lonu(im1,je)
    clonv(2,i,je) = lonu(im1,je)
    clonv(3,i,je) = lonu(i  ,je)
    clonv(4,i,je) = lonu(i  ,je)
    clatv(1,i,je) = latu(im1,je)
    clatv(2,i,je) = latu(im1,je)
    clatv(3,i,je) = latu(i  ,je)
    clatv(4,i,je) = latu(i  ,je)
  ENDDO

  mask(:,:) = 1

  ! write scalar grid to ofilenames.nc

  flen = LEN_TRIM(ofilename)

  ofilename = replace(ofilename, flen+1, 's.nc')
  CALL write_grid(CHAR(ofilename), je, ie-2, lats, lons, clats, clons)

  ! write vector-u grid to ofilenameu.nc

  ofilename = replace(ofilename, flen+1, 'u.nc')
  CALL write_grid(CHAR(ofilename), je, ie-2, latu, lonu, clatu, clonu)

  ! write vector-v grid to ofilenamev.nc

  ofilename = replace(ofilename, flen+1, 'v.nc')
  CALL write_grid(CHAR(ofilename), je, ie-2, latv, lonv, clatv, clonv)

END PROGRAM anta2nc

SUBROUTINE nfce(status)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: status

  INCLUDE 'netcdf.inc'

  IF ( status .NE. NF_NOERR ) THEN
    WRITE(0,*) NF_STRERROR(status)
    STOP 'netCDF error'
  END IF

END SUBROUTINE nfce

SUBROUTINE err_exit
  IMPLICIT NONE
  stop 1
END SUBROUTINE err_exit

SUBROUTINE write_grid(ofilename, nlat, nlon, grid_center_lat, grid_center_lon, &
                      grid_corner_lat, grid_corner_lon)

  USE mo_kind, ONLY: sp
  IMPLICIT NONE

  CHARACTER (*) :: ofilename

  INTEGER :: nlon, nlat

  REAL(sp) :: grid_center_lat(nlon, nlat), grid_center_lon(nlon, nlat), grid_corner_lat(4,nlon, nlat), grid_corner_lon(4,nlon, nlat)


  INTEGER :: grid_size

  INTEGER :: grid_corners, grid_rank, grid_dims(2)
  INTEGER :: fileid, nc_dims_id(3)
  INTEGER :: nc_gridsize_id, nc_gridcorn_id, nc_gridrank_id, nc_griddims_id
  INTEGER :: nc_gridxsize_id, nc_gridysize_id
  INTEGER :: nc_grdcntrlat_id, nc_grdcntrlon_id, nc_grdimask_id
  INTEGER :: nc_grdcrnrlat_id, nc_grdcrnrlon_id
  INTEGER :: grid_imask(nlon*nlat)

  INCLUDE 'netcdf.inc'

  grid_size    = nlon*nlat
  grid_dims(1) = nlon
  grid_dims(2) = nlat
  grid_corners = 4
  grid_rank    = 2
  grid_imask(:) = 1

  WRITE(*,*) 'write grid to ', ofilename
  !***
  !*** create netCDF dataset for this grid
  !***
  CALL nfce(nf_create (ofilename, NF_CLOBBER, fileid))

  CALL nfce(nf_put_att_text (fileid, NF_GLOBAL, 'title', len_trim(ofilename)-3, ofilename))

  !***
  !*** define grid size dimension
  !***
  CALL nfce(nf_def_dim (fileid, 'grid_size', grid_size, nc_gridsize_id))
  CALL nfce(nf_def_dim (fileid, 'grid_xsize', nlon, nc_gridxsize_id))
  CALL nfce(nf_def_dim (fileid, 'grid_ysize', nlat, nc_gridysize_id))

  !***
  !*** define grid corner dimension
  !***
  CALL nfce(nf_def_dim (fileid, 'grid_corners', grid_corners, nc_gridcorn_id))

  !***
  !*** define grid rank dimension
  !***
  CALL nfce(nf_def_dim (fileid, 'grid_rank', grid_rank, nc_gridrank_id))

  !***
  !*** define grid dimension size array
  !***
  nc_dims_id(1) = nc_gridrank_id
  CALL nfce(nf_def_var (fileid, 'grid_dims', NF_INT, 1, nc_dims_id, nc_griddims_id))

  !***
  !*** define grid center latitude array
  !***
  nc_dims_id(1) = nc_gridxsize_id
  nc_dims_id(2) = nc_gridysize_id

  CALL nfce(nf_def_var (fileid, 'grid_center_lat', NF_REAL, 2, nc_dims_id, nc_grdcntrlat_id))

  CALL nfce(nf_put_att_text (fileid, nc_grdcntrlat_id, 'units', 7, 'degrees'))
  CALL nfce(nf_put_att_text (fileid, nc_grdcntrlat_id, 'bounds', 15, 'grid_corner_lat'))

  !***
  !*** define grid center longitude array
  !***
  CALL nfce(nf_def_var (fileid, 'grid_center_lon', NF_REAL, 2, nc_dims_id, nc_grdcntrlon_id))

  CALL nfce(nf_put_att_text (fileid, nc_grdcntrlon_id, 'units', 7, 'degrees'))
  CALL nfce(nf_put_att_text (fileid, nc_grdcntrlon_id, 'bounds', 15, 'grid_corner_lon'))

  !***
  !*** define grid mask
  !***
  CALL nfce(nf_def_var (fileid, 'grid_imask', NF_INT, 2, nc_dims_id, nc_grdimask_id))

  CALL nfce(nf_put_att_text (fileid, nc_grdimask_id, 'units', 8, 'unitless'))
  CALL nfce(nf_put_att_text (fileid, nc_grdimask_id, 'coordinates', 31, 'grid_center_lon grid_center_lat'))

  !***
  !*** define grid corner latitude array
  !***
  nc_dims_id(1) = nc_gridcorn_id
  nc_dims_id(2) = nc_gridxsize_id
  nc_dims_id(3) = nc_gridysize_id

  CALL nfce(nf_def_var (fileid, 'grid_corner_lat', NF_REAL, 3, nc_dims_id, nc_grdcrnrlat_id))

  CALL nfce(nf_put_att_text (fileid, nc_grdcrnrlat_id, 'units', 7, 'degrees'))

  !***
  !*** define grid corner longitude array
  !***
  CALL nfce(nf_def_var (fileid, 'grid_corner_lon', NF_REAL, 3, nc_dims_id, nc_grdcrnrlon_id))

  CALL nfce(nf_put_att_text (fileid, nc_grdcrnrlon_id, 'units', 7, 'degrees'))

  !***
  !*** end definition stage
  !***
  CALL nfce(nf_enddef(fileid))

  !-----------------------------------------------------------------------
  !
  !     write grid data
  !
  !-----------------------------------------------------------------------

  CALL nfce(nf_put_var_int(fileid, nc_griddims_id, grid_dims))

  CALL nfce(nf_put_var_int(fileid, nc_grdimask_id, grid_imask))

  CALL nfce(nf_put_var_real(fileid, nc_grdcntrlat_id, grid_center_lat))

  CALL nfce(nf_put_var_real(fileid, nc_grdcntrlon_id, grid_center_lon))

  CALL nfce(nf_put_var_real(fileid, nc_grdcrnrlat_id, grid_corner_lat))

  CALL nfce(nf_put_var_real(fileid, nc_grdcrnrlon_id, grid_corner_lon))

  CALL nfce(nf_close(fileid))

END SUBROUTINE write_grid
