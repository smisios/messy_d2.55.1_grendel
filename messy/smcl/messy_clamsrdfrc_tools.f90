MODULE messy_clamsrdfrc_tools
!***********************************************************************
! contains tools for CLaMS radiative forcing in EMAC
! - SUBROUTINE nc_read_all_ap
! - SUBROUTINE get_apg
! - SUBROUTINE set_coor
! - SUBROUTINE get_aps_and_datas
! - SUBROUTINE nc_write_parameters 
! - SUBROUTINE nc_read_ap_s_info
! - SUBROUTINE ncdf_error
!***********************************************************************
  USE messy_clams_global, ONLY: PREC

  IMPLICIT NONE

  PUBLIC :: nc_read_all_ap
  PUBLIC :: get_apg
  PUBLIC :: set_coor
  PUBLIC :: get_aps_and_datas
  PUBLIC :: nc_write_parameters 
  PUBLIC :: nc_read_ap_s_info
  PRIVATE :: ncdf_error

CONTAINS
  
!*************************************************************************
! read points from netcdf-file (lat, lon, lev, time)
!*************************************************************************
SUBROUTINE nc_read_all_ap (LAT, LON, LEV, ap_all)

  USE messy_clams_global, ONLY: nparts_max
  USE messy_clamsrdfrc_global, ONLY: ap

  IMPLICIT NONE

  REAL(PREC), DIMENSION(:), POINTER :: LAT, LON, LEV
  TYPE(ap), DIMENSION(:), POINTER :: ap_all

  INTEGER :: nparts_local

  nparts_local = nparts_max

  ALLOCATE (ap_all(nparts_local))
  ap_all(:)%lat = LAT
  ap_all(:)%lon = LON 
  ap_all(:)%lev = LEV

END SUBROUTINE nc_read_all_ap

!*******************************************************************************
!  transform a layer between lev_min and lev_max to the structure ap_g
!*******************************************************************************
SUBROUTINE get_apg (ap_g, indexarr, lev_min, lev_max, data_present)

  USE messy_clamsrdfrc_global, ONLY: ap, npart_g, ap_g_all
  
  IMPLICIT NONE
  
  TYPE(ap), DIMENSION(:), POINTER  :: ap_g
  INTEGER,  DIMENSION(:), POINTER  :: indexarr
  REAL(PREC), INTENT(in)           :: lev_min, lev_max
  LOGICAL, INTENT(out)             :: data_present
  
  LOGICAL, DIMENSION(:), POINTER   :: mask
  INTEGER                          :: nparts_level, i, k

  ALLOCATE (mask(npart_g))

  mask = .FALSE.
  WHERE (lev_min<=ap_g_all(:)%lev .AND. ap_g_all(:)%lev<lev_max)
     mask = .TRUE.
  END WHERE
  
  nparts_level = COUNT(mask)
  
  IF (nparts_level == 0) THEN
     data_present = .FALSE.
  ELSE
     data_present = .TRUE.
     ALLOCATE (ap_g(nparts_level))
     ALLOCATE (indexarr(nparts_level))

     i = 0
     DO k = 1, npart_g
        IF (mask(k)) THEN
           i = i+1
           indexarr(i) = k
        ENDIF
     ENDDO
    
     ap_g = PACK(ap_g_all, mask)
  ENDIF
 
  DEALLOCATE (mask)
  
END SUBROUTINE get_apg

!****************************************************************************
! Transform ap-coordinates (lon, lat) on kart. coord. on a unit sphere
!****************************************************************************
SUBROUTINE set_coor (ap_array, vert_array)

  USE messy_clamsrdfrc_global, ONLY: ap

  IMPLICIT NONE
    
  TYPE(ap), DIMENSION(:),   POINTER :: ap_array
  REAL(PREC),     DIMENSION(:,:), POINTER :: vert_array
 
  REAL(PREC)  :: pi
  REAL(PREC), DIMENSION(:), ALLOCATABLE :: lons_rad, lats_rad

  pi=4.*ATAN(1.)

  ALLOCATE (lons_rad(SIZE(ap_array)))
  ALLOCATE (lats_rad(SIZE(ap_array)))

  lons_rad = ap_array%lon * pi / 180.          ! transform from deg to rad
  lats_rad = (90. - ap_array%lat) * pi / 180.  ! use spherical coordinates
                                               !   0. < lons < 2*PI
                                               !   0. < lats < PI
  vert_array(1,:) = COS(lons_rad)*SIN(lats_rad) ! transf from spheric
  vert_array(2,:) = SIN(lons_rad)*SIN(lats_rad) ! coord. to  kart.
  vert_array(3,:) = COS(lats_rad)  ! coord on a unit sphere (r=1)

  DEALLOCATE (lons_rad, lats_rad)

END SUBROUTINE set_coor

!*******************************************************************************
!  transform a layer between lev_min and lev_max to the structure ap_s
!  and write appropriate data on structure data_s
!*******************************************************************************
SUBROUTINE get_aps_and_datas (ap_s, data_s, lev_min, lev_max, data_present)
   
  USE messy_clamsrdfrc_global, ONLY: ap, npart_s, ap_s_all, ntags, data_s_all
  USE messy_clams_global, ONLY: nparts_max

  IMPLICIT NONE
  
  TYPE(ap), DIMENSION(:),   POINTER :: ap_s
  REAL(PREC),     DIMENSION(:,:), POINTER :: data_s
  REAL(PREC)                              :: lev_min, lev_max
  LOGICAL, INTENT(out)              :: data_present
  
  LOGICAL, DIMENSION(:), POINTER  :: mask
  INTEGER                         :: nparts_level, i

  npart_s = nparts_max
  ALLOCATE (mask(npart_s))

  mask = .FALSE.
  WHERE (lev_min<=ap_s_all(:)%lev .AND. ap_s_all(:)%lev<lev_max)
     mask = .TRUE.
  END WHERE
  
  nparts_level = COUNT(mask)
  write(*,*) 'nparts_level', nparts_level
  IF (nparts_level == 0) THEN
     data_present = .FALSE.
  ELSE
     data_present = .TRUE.
     
     ALLOCATE (ap_s(nparts_level))
     ap_s = PACK(ap_s_all, mask)
     
     ALLOCATE (data_s(ntags, nparts_level))
     DO i = 1, ntags
        data_s(i,:) = PACK(data_s_all(i,:), mask)
     ENDDO
     
  ENDIF
  
  DEALLOCATE (mask)
  
END SUBROUTINE get_aps_and_datas
!*************************************************************************
! 
!*************************************************************************
SUBROUTINE nc_write_parameters (nlev,nlat,nlon, ECHAM_DATA)

  USE messy_clamsrdfrc_global, ONLY: data_g_tot

  IMPLICIT NONE

  INTEGER      :: nlev, nlat, nlon
  REAL(PREC), DIMENSION(:,:,:) :: ECHAM_DATA  

  REAL(PREC), DIMENSION(:,:,:), POINTER :: arr3d
  
  INTEGER :: ilon, ilat, ilev, ipart
  
  ! write new variables
  ALLOCATE (arr3d(nlev,nlat,nlon))

     ! data_g -> arr3d
     ipart = 0
     DO ilon = 1, nlon
        DO ilat = 1, nlat
           DO ilev = 1, nlev
              ipart = ipart + 1
              arr3d(ilev,ilat,ilon) = data_g_tot(1,ipart)
           ENDDO
        ENDDO
       ENDDO
       
       DO ilon = 1, nlon
          DO ilat = 1, nlat
             DO ilev = 1, nlev
                ECHAM_DATA(ilon,ilev,ilat) = arr3d(ilev,ilat,ilon)
             ENDDO
          ENDDO
       ENDDO       
       
!!$write(*,*) 'SIZE(ECHAM_DATA,1)',SIZE(ECHAM_DATA,1)
!!$write(*,*) 'SIZE(ECHAM_DATA,2)',SIZE(ECHAM_DATA,2)
!!$write(*,*) 'SIZE(ECHAM_DATA,3)',SIZE(ECHAM_DATA,3)

    DEALLOCATE (arr3d)
    
  END SUBROUTINE nc_write_parameters
  
!*************************************************************************
! read information from netcdf-file:
! number of points, chemical species, lev_grid, some global variables
!*************************************************************************
SUBROUTINE nc_read_ap_s_info (file_init, lev_window, dir)

  USE messy_clams_tools_dateconv
  USE messy_clams_tools_utils,  ONLY: lowercase, uppercase
  USE messy_clams_tools_ncutils, ONLY: nc_get_vertcoorname
  USE netcdf
  USE messy_clamsrdfrc_global, ONLY: datetype, lat_down, lat_up, r_coarse, r_high

  IMPLICIT NONE 

  CHARACTER*(*), INTENT(in)                        :: file_init
  REAL(PREC), DIMENSION(:),POINTER                       :: lev_grid 
  REAL(PREC), DIMENSION(:,:),POINTER                     :: lev_window
  TYPE(datetype)                                   :: date
  CHARACTER*(*), INTENT(in), OPTIONAL              :: dir   

  INTEGER                                          :: nparts_local
  INTEGER                                          :: nlevs

  CHARACTER(100)               :: file_init_tot
  CHARACTER(40)                :: vertcoorname
  REAL                         :: lev_min, lev_max
  REAL, DIMENSION(:),POINTER   :: lev_delta
  REAL                         :: dummy
  REAL                         :: lev_act

  INTEGER,DIMENSION(2)         :: tstart,tcount,dim_array
  INTEGER, DIMENSION(1)        :: data_int
  INTEGER                      :: rcode, ncid, dummy_id
  INTEGER                      :: i,dim_len,text_len,varid
  INTEGER                      :: nlevs_value
  REAL, DIMENSION(1)           :: data_real
  REAL(KIND(0d0)),DIMENSION(1) :: jultime


  file_init_tot=file_init
  IF (PRESENT(dir)) file_init_tot=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(file_init))    

  PRINT *, 'File: ', TRIM(file_init_tot)

  rcode= nf90_open(file_init_tot,nf90_nowrite,ncid)
  IF (rcode /= nf90_noerr) THEN
     WRITE(*,*) 'Error on open file ',file_init_tot
     IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)
  ENDIF

  CALL nc_get_vertcoorname(ncid,vertcoorname)

  ! get time
  rcode = nf90_inq_varid (ncid,'time',varid)
  rcode = nf90_get_var (ncid,varid,jultime)
  CALL js2ymds (jultime(1),date%year,date%month,date%day,date%sec)
  WRITE (*,'(A,4I6)') ' time (year, month, day, hour) : ', &
       date%year, date%month, date%day, date%sec/3600 

  ! Read nparts (# of APs)
  rcode= nf90_inq_dimid(ncid,'NPARTS',dummy_id)     
  IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)
  rcode= nf90_inquire_dimension(ncid,dummy_id,len=nparts_local)
  IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)

  ! Read global parameters describing initial distribution
  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_'//TRIM(lowercase(vertcoorname))//'_min', data_real)
  lev_min=data_real(1)
  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_'//TRIM(lowercase(vertcoorname))//'_max', data_real)
  lev_max=data_real(1)
  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_n'//TRIM(lowercase(vertcoorname))//'s', data_int)
  nlevs=data_int(1)
  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_lat_down', data_real)
  lat_down=data_real(1)
  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_lat_up', data_real)
  lat_up=data_real(1)

  ! lat_min und lat_max werden aus i3d.inp gelesen (NICHT ueberschreiben!)
  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_lat_min', data_real)
  dummy=data_real(1)
  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_lat_max', data_real)
  dummy=data_real(1)

  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_r_coarse', data_real)
  r_coarse=data_real(1)
  rcode=nf90_get_att(ncid,NF90_GLOBAL,'exp_POS_r_high', data_real)
  r_high=data_real(1)

  ! Read nlevs (# of lev niveaus)
  rcode= nf90_inq_dimid(ncid,'N'//TRIM(uppercase(vertcoorname))//'S',dummy_id)     
  IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)

  ! read lev information 
  rcode = nf90_inquire_dimension(ncid,dummy_id,len=nlevs_value)
  IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)
  ALLOCATE(lev_grid(nlevs_value),lev_delta(nlevs_value))
  
  rcode=nf90_inq_varid(ncid,TRIM(uppercase(vertcoorname))//'_GRID',varid)     
  IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)
  rcode=nf90_get_var(ncid,varid,lev_grid)
  IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)
  rcode=nf90_inq_varid(ncid,TRIM(uppercase(vertcoorname))//'_DELTA',varid)
  IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)
  rcode=nf90_get_var(ncid,varid,lev_delta)
  IF (rcode /= nf90_noerr) CALL ncdf_error(rcode)

  ALLOCATE (lev_window(nlevs,2))
  lev_act = lev_min
  IF (nlevs==1) THEN
     WRITE (*,*) '2-D-Interpolation!!!'
     lev_window(1,1) = 0
     lev_window(1,2) = 10000.
  ELSE
     DO i=1, nlevs  
        lev_window(i,1) = lev_act
        lev_window(i,2) = lev_act+lev_delta(i)
        lev_act = lev_act+lev_delta(i)
     END DO
  ENDIF

  WRITE (*,*) 'n',TRIM(lowercase(vertcoorname)),'s=',nlevs
  WRITE (*,*)

  rcode=nf90_close(ncid)

END SUBROUTINE nc_read_ap_s_info

!*************************************************************************
! handels NetCDF errors 
!*************************************************************************
SUBROUTINE ncdf_error (error)
  
  USE netcdf

  IMPLICIT NONE

  INTEGER :: error

  WRITE (*,*)
  WRITE (*,*) '***********************************************************'

  WRITE (*,*) 'NetCDF Error:'
  
  WRITE (*,*) nf90_strerror(error)

  WRITE (*,*) '***********************************************************'
  WRITE (*,*)
  WRITE (*,*) 'THE PROGRAM WILL STOP !!!'   
  WRITE (*,*)
  WRITE (*,*)
  STOP

END SUBROUTINE ncdf_error

END MODULE messy_clamsrdfrc_tools
