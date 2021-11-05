!*******************************************************************************
!
! Mask routine used to generate masks taking into account the land-sea 
! distribution. It uses land-sea mask and orography at 0,125 degree from
! ECMWF. 
!
! Written originally by Ann'Sophie Tissier (matlab)
! 04/2014 Modified by B. Legras 
! 02/2016 New refined regions in Asia 
! 04/2016 Converted to FORTRAN90, Nicole Thomas
!
! edge rules: long >= and <
!             lat  >= and <
! except for N boundaries where <= can be used
!
!*******************************************************************************
!
! Module messy_clamstracer_regional_masks contains the following subroutines:
!
!    Subroutine read_era_data (filename, data, varname)
!    Function set_lsm_mask (latitude, longitude, lsm)
!    Function set_oro_mask (latitude, longitude, gph)
!         contains: function interpol_lin (x1, y1, x2, y2, x)
!    Subroutine define_regional_mask (longitude, latitude, level, levelbound, &
!                                     maskname, maskval, lsm, gph)
!
!*******************************************************************************

module messy_clamstracer_regional_masks

  USE messy_clams_global, ONLY: PREC, rank

  implicit none

  type eradata_type
     character(30)                           :: name
     integer                                 :: nlon, nlat
     real(prec)                              :: dlon, dlat, lon0,lat0
     real(prec), dimension(:),   allocatable :: lon, lat
     real(prec), dimension(:,:), allocatable :: values
  end type eradata_type


contains 


   !*******************************************************************
   ! Read 2D-variable from ECMWF file
   !*******************************************************************
   subroutine read_era_data (filename, data, varname)

     use netcdf
     use messy_clams_tools_ncutils, only: nc_check_error

     implicit none

     character(*)        :: filename, varname
     type (eradata_type) :: data

     integer :: status, ncid, dimid, varid
     integer :: ind(1)

     ! open netcdf file
     if (rank==0) write (*,*) 'read file ',trim(filename)
     status = nf90_open (filename, nf90_nowrite, ncid)
     call nc_check_error (status, "Cannont open file "//trim(filename))

     ! set name
     data%name = varname
     
     ! read dimensions lon, lat
     status = nf90_inq_dimid (ncid, "longitude", dimid)
     call nc_check_error (status, "Cannont find dimension longitude")
     status = nf90_inquire_dimension (ncid, dimid, len=data%nlon)
     call nc_check_error (status, "Cannont read dimension longitude")
     status = nf90_inq_dimid (ncid, "latitude", dimid)
     call nc_check_error (status, "Cannont find dimension latitude")
     status = nf90_inquire_dimension (ncid, dimid, len=data%nlat)
     call nc_check_error (status, "Cannont read dimension latitude")

     ! read coordinate variables lon, lat
     allocate (data%lon(data%nlon))
     allocate (data%lat(data%nlat))
     status = nf90_inq_varid (ncid, "longitude", varid)
     call nc_check_error (status, "Cannont find variable longitude")
     status = nf90_get_var (ncid, varid, data%lon)
     call nc_check_error (status, "Cannont read variable longitude")
     status = nf90_inq_varid (ncid, "latitude", varid)
     call nc_check_error (status, "Cannont find variable latitude")
     status = nf90_get_var (ncid, varid, data%lat)
     call nc_check_error (status, "Cannont read variable latitude")
     data%lon0 = data%lon(1)
     data%lat0 = data%lat(1)
     data%dlon = data%lon(2)-data%lon(1)
     data%dlat = data%lat(2)-data%lat(1)

     ! read parameter
     allocate (data%values(data%nlon,data%nlat))
     status = nf90_inq_varid (ncid, TRIM(data%name), varid)
     call nc_check_error (status, "Cannont find variable "//trim(data%name))
     status = nf90_get_var (ncid, varid, data%values)
     call nc_check_error (status, "Cannont read variable "//trim(data%name))
     
     ! longitudes [-180,180] => [0,360]
     if (data%lon0 < 0.) then
        
        ind = minloc(abs(data%lon))
        if (data%lon(ind(1)) < 0.) ind(1) = ind(1)+1
        if (ind(1) > 1) then
           
           data%lon0 = data%lon(ind(1))
           
           ! shift longitudes
           if (rank==0) write (*,*) 'shift longitudes'
           data%lon = cshift(data%lon,ind(1)-1)
           where (data%lon < 0.) 
              data%lon = data%lon + 360.
           end where
           
           ! shift parameter
           if (rank==0) write (*,*) 'shift ',trim(data%name)             
           data%values = cshift(data%values,ind(1)-1,1)
           
        endif

     endif

     if (rank==0) write (*,*) 'nlon, nlat, lon0, lat0, dlon, dlat=', &
          data%nlon, data%nlat, data%lon0, data%lat0, data%dlon, data%dlat

     ! close netcdf file
     status = nf90_close (ncid)

  end subroutine read_era_data

  !*******************************************************************
  ! Get land-sea mask for current latitudes and longitudes 
  !*******************************************************************
  Function set_lsm_mask (latitude, longitude, lsm)

    use messy_clams_global,  ONLY: mdi, eps

    implicit none
    
    real(prec), dimension(:)  :: latitude, longitude
    type (eradata_type)       :: lsm ! land-sea mask from ECMWF
    real(prec), dimension(size(latitude)) :: set_lsm_mask
    
    integer :: npoints, i
    integer :: ilon, ilat
    
    !set_lsm_mask = -1.
    set_lsm_mask = mdi
    
    npoints = size(latitude)

    do i = 1, npoints

       
       if (ABS((latitude(i)-mdi)/mdi)>eps) then
    
          !------------------------------------------------
          ! get value of nearest grid point
          !------------------------------------------------
          
          ilon = (longitude(i)-lsm%lon0)/lsm%dlon+1
          if (ilon==lsm%nlon) then
             if (abs(longitude(i)-lsm%lon(1)+360.) < abs(longitude(i)-lsm%lon(lsm%nlon))) then
                ilon = 1
             endif
          else
             if (abs(longitude(i)-lsm%lon(ilon+1)) < abs(longitude(i)-lsm%lon(ilon))) then
                ilon = ilon+1
             endif
          endif
          
          ilat = (latitude(i)-lsm%lat0)/lsm%dlat+1
          if (abs(latitude(i)-lsm%lat(ilat+1)) < abs(latitude(i)-lsm%lat(ilat))) then
             ilat = ilat+1
          endif
          
          set_lsm_mask(i) = lsm%values(ilon,ilat)
          
       endif

    enddo

  End Function set_lsm_mask

  !*********************************************************************
  ! Get orographie (gph on surface) for current latitudes and longitudes
  !*********************************************************************
  Function set_oro_mask (latitude, longitude, gph)

    use messy_clams_global,  ONLY: mdi, eps

    implicit none
    
    real(prec), dimension(:)              :: latitude, longitude
    type (eradata_type)                   :: gph ! gph (on surface) from ECMWF
    real(prec), dimension(size(latitude)) :: set_oro_mask

    real(prec)    :: val1, val2
    integer       :: npoints
    integer       :: i, ilon1, ilon2, ilat1, ilat2


    !set_oro_mask = 0.
    set_oro_mask = mdi

    npoints = size(latitude)

    do i = 1, npoints

       if (ABS((latitude(i)-mdi)/mdi)>eps) then
    
          !------------------------------------------------
          ! 2d linear interpolation
          !------------------------------------------------
          
          ilon1 = (longitude(i)-gph%lon0)/gph%dlon+1
          ilon2 = mod(ilon1,gph%nlon)+1

          ilat1 = (latitude(i)-gph%lat0)/gph%dlat+1
          ilat2 = ilat1 +1
          
          val1 = interpol_lin (gph%lon(ilon1),gph%values(ilon1,ilat1), &
                               gph%lon(ilon2),gph%values(ilon2,ilat1),longitude(i))
          val2 = interpol_lin (gph%lon(ilon1),gph%values(ilon1,ilat2), &
                               gph%lon(ilon2),gph%values(ilon2,ilat2),longitude(i))
          set_oro_mask(i) = interpol_lin (gph%lat(ilat1),val1, &
                                          gph%lat(ilat2),val2,latitude(i))

       endif

    enddo

    contains

      function interpol_lin (x1, y1, x2, y2, x)
        
        implicit none
        
        real(prec) :: interpol_lin, x1, y1, x2, y2, x
        
        if (x1==x2) then
           interpol_lin = y1
        else     
           interpol_lin = y1 + (y2-y1)/(x2-x1)*(x-x1)
        endif
        
      end function interpol_lin

  End Function set_oro_mask

  !*******************************************************************
  ! Define regional masks
  !*******************************************************************
  Subroutine define_regional_mask (longitude, latitude, level, levelbound, &
                                   maskname, maskval, lsm, gph)

    use messy_clams_global,  ONLY: mdi, eps

    implicit none

    real(prec), dimension(:) :: longitude, latitude, level
    logical,    dimension(:) :: maskval
    real(prec)               :: levelbound
    character(*)             :: maskname
    type (eradata_type)      :: lsm ! land-sea mask from ECMWF
    type (eradata_type)      :: gph ! gph (on surface) from ECMWF

    real(prec), dimension(:), pointer :: lsm_mask, oro_mask

    allocate (lsm_mask(size(longitude)))
    allocate (oro_mask(size(longitude)))

    maskval = .false.

    select case (maskname)

    case ('l20')
       where (latitude<=20. .and. latitude>=-20) 
          maskval = .true.
       end where

    case ('NAf') ! North Africa

       where ((longitude<51 .or. longitude>=340) .and. &
              (latitude<=20 .and. latitude>=0))
          maskval = .true.
       endwhere
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       where (lsm_mask<1) 
          maskval = .false.
       end where
   
       ! Big lakes
       where (longitude<=37.5 .and. longitude>=30 .and. &
              latitude<=3.5 .and. latitude>=0) 
          maskval = .true.
       end where

       ! Other lakes
       where (longitude>=15 .and. longitude<=38 .and. &
              latitude>=0 .and. latitude<=15)
          maskval = .true.
       end where
       where (longitude>=14 .and. longitude<=15 .and. &
              latitude>=12.5 .and. latitude<=13.5)
          maskval = .true.
       end where
       where (longitude>=4 .and. longitude<=5 .and. &
              latitude>=10 .and. latitude<=11)
          maskval = .true.
       end where
       where (longitude>=0 .and. longitude<=1 .and. &
              latitude>=6.5 .and. latitude<=8) 
          maskval = .true.
       end where
       where (longitude>=349 .and. longitude<=350 .and. &
              latitude>=12.5 .and. latitude<=14)
          maskval = .true.
       end where
       where (longitude>=349 .and. longitude<=360 .and. &
              latitude>=6 .and. latitude<=7.5)
          maskval = .true.
       end where

    case ('SAf') ! South Africa 

       where ((longitude<51 .or. longitude>=340) .and. &
              (latitude<0 .and. latitude>=-20))
          maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       where (lsm_mask<1) 
          maskval = .false.
       end where
       
       ! Big lakes
       where (longitude<=37.5 .and. longitude>=30 .and. &
              latitude<0 .and. latitude>=-15)
          maskval = .true.
       end where
       ! Islands
       where (longitude<=44 .and. longitude>=43 .and. &
              latitude>=-12 .and. latitude<=-11)
          maskval = .false.
       end where
       ! Other lakes
       where (longitude>=15 .and. longitude<=38 .and. &
              latitude>=-15 .and. latitude<0)
          maskval = .true.
       end where
       where (longitude>=27 .and. longitude<=30 .and. &
              latitude>=-18 .and. latitude<=-16)
          maskval = .true.
       end where
       where (longitude>=32 .and. longitude<=36 .and. &
              latitude>=-16 .and. latitude<=-14)
          maskval = .true.
       end where

    case ('ITA') ! Inter Tropical Atlantic, previously Africa Ocean

       where (longitude<325 .and. longitude>=310 .and. &
              latitude<=20 .and. latitude>=0)
          maskval = .true.
       end where
       where ((longitude<15 .or. longitude>=325) .and. &
              (latitude<=20  .and. latitude>=-20))
          maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       !where (lsm_mask==1)    
       where (lsm_mask/=0)  ! = 1 or mdi    
          maskval = .false.
       end where
       
       ! Islands
       where (longitude>=335 .and. longitude<=337 .and. &
              latitude>=14 .and. latitude<=16)
          maskval = .true.
       end where
       ! Lakes
       where (longitude>=0 .and. longitude<=16 .and. &
              latitude>=9 .and. latitude<=46)
          maskval = .false.
       end where
       where (longitude>=0 .and. longitude<=1 .and. &
              latitude>=6.4 .and. latitude<=8) 
           maskval = .false.
       end where

    case ('NAP') ! North Asia Pacific

       where (longitude>=66  .and. longitude<=150 .and. &
              latitude<=35 .and. latitude>=5)
           maskval = .true.
       end where

    case ('WPac') !  West Pacific

       where (longitude>150 .and. longitude<210 .and. &
              latitude<5 .and. latitude>=-20)
           maskval = .true.
       end where

    case ('IndMal') ! Indonesie + Malaysia

       where (longitude>=95 .and. longitude<150 .and. &
              latitude<5 .and. latitude>=-11)
           maskval = .true.
       end where
       ! include north Sumatra
       where (longitude>=95 .and. longitude<=99 .and. &
              latitude<=7 .and. latitude>=5)
           maskval = .true.
       end where
       ! include north Borneo
       where (longitude>114 .and. longitude<=119.5 .and. &
              latitude>=5 .and. latitude<=7.5)
           maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       !where (lsm_mask==0)    
       where (lsm_mask/=1)  ! 0 or mdi    
          maskval = .false.
       end where

    case ('WPool') ! WarmPool

       where (longitude>=95 .and. longitude<150 .and. &
              latitude<5 .and. latitude>=-11)
           maskval = .true.
       end where
       where (longitude>=113 .and. longitude<150 .and. &
              latitude<-11 .and. latitude>=-20)
           maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       !where (lsm_mask==1)    
       where (lsm_mask/=0)    
          maskval = .false.
       end where

    case ('NAus') ! North Australia

       where (longitude>=113 .and. longitude<150 .and. &
               latitude<-11 .and. latitude>=-20)
           maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       !where (lsm_mask==0)    
       where (lsm_mask/=1)    
          maskval = .false.
       end where
       where (longitude<=130 .and. longitude>=128 .and. &
              latitude>=-15.5 .and. latitude<=-17.5)
           maskval = .true.
       end where
    
    case ('SAus') ! South Australia

       where (longitude>=113 .and. longitude<154 .and. &
              latitude<=-20 .and. latitude>=-39)
           maskval = .true.
       end where

    case ('IO') ! Indian Ocean

       where (longitude<95 .and. longitude>=35 .and. &
              latitude<5 .and. latitude>=-20)
           maskval = .true.
       end where
       where (longitude<78 .and. longitude>=43 .and. &
              latitude<26 .and. latitude>=5) 
           maskval = .true.
       end where
       where (longitude<113 .and. longitude>=95 .and. &
              latitude<-11 .and. latitude>=-20)
           maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       !where (lsm_mask==1)    
       where (lsm_mask/=0)    
          maskval = .false.
       end where
       
       ! Islands
       where (longitude<=44 .and. longitude>=43 .and. &
              latitude>=-12 .and. latitude<=-11)
           maskval = .true.
       end where
       where (longitude<=56 .and. longitude>=52 .and. &
              latitude>=11.5 .and. latitude<=13.5)
           maskval = .true.
       end where
       
       ! Lakes
       where (longitude<=38 .and. longitude>=34 .and. &
              latitude>=1.5 .and. latitude<=5.5)
           maskval = .false.
       end where
       where (longitude<=36 .and. longitude>=34 .and. &
              latitude>=-15 .and. latitude<=-13)
           maskval = .false.
       end where
       
    case ('India') ! Indian peninsula + western regions 

       where (longitude>=66 .and. longitude<92 .and. &
               latitude<=35 .and. latitude>=5)
           maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       where (lsm_mask<1)    
          maskval = .false.
       end where

       oro_mask = set_oro_mask(latitude, longitude, gph)
       !where (oro_mask>=3500.*9.81) 
       where (oro_mask>=3500.*9.81 .or. ABS((oro_mask-mdi)/mdi)<=eps) 
          maskval = .false.
       end where

       ! Lakes
!!$       where (longitude>=80 .and. longitude<=92 .and. &
!!$              latitude>=30 .and. latitude<=35)
!!$           maskval = .true.
!!$       end where
       where (longitude>=75.5 .and. longitude<=76.5 .and. &
              latitude>=31.5 .and. latitude<=32.5)
           maskval = .true.
       end where

    case ('China') ! Central China

       where (longitude>=92 .and. longitude<125 .and. &
              latitude<=35 .and. latitude>=23)
           maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       where (lsm_mask<1)    
          maskval = .false.
       end where

       oro_mask = set_oro_mask(latitude, longitude, gph)
       where (oro_mask>=3500.*9.81 .or. ABS((oro_mask-mdi)/mdi)<=eps) 
          maskval = .false.
       end where
       
       ! Lakes
       where (longitude>=115 .and. longitude<=117 .and. &
              latitude>=28 .and. latitude<=32)
           maskval = .true.
       end where
       where (longitude>=118 .and. longitude<=120 .and. &
              latitude>=33 .and. latitude<=34)
           maskval = .true.
       end where
       where (longitude>=119.10245 .and. longitude<=120.5 .and. &
              latitude>=30.5 .and. latitude<=32)
           maskval = .true.
       end where

    case ('Pen') ! 

       where (longitude>=92 .and. longitude<=115 .and. &
              latitude<=23 .and. latitude>=5)
           maskval = .true.
       end where
       ! mask north Sumatra
       where (longitude>=95 .and. longitude<=99 .and. &
              latitude<=7 .and. latitude>=5)
           maskval = .false.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       where (lsm_mask<1)    
          maskval = .false.
       end where
    
       ! Erase Andaman Islands
       where (longitude>=92 .and. longitude<=94 .and. &
              latitude>=11 .and. latitude<=14)
           maskval = .false.
       end where
       ! Include Lakes
       where (longitude>=103.5 .and. longitude<=105 .and. &
              latitude>=12 .and. latitude<=14)
           maskval = .true.
       end where
       where (longitude>=103.5 .and. longitude<=104 .and. &
              latitude>=10.8 .and. latitude<=11.2)
           maskval = .true.
       end where
       where (longitude>=102 .and. longitude<=103 .and. &
              latitude>=16 .and. latitude<=20)
           maskval = .true.
       end where
       where (longitude>=102.3 .and. longitude<=103 .and. &
              latitude>=5 .and. latitude<=5.3)
           maskval = .true.
       end where


!!$*****************************************************+    
!!$    case 'AML' % Asia Main Land previously NAPL3
!!$     ilon = (longitude>=66) & (longitude<=150);
!!$     ilat = (latitude<=35) & (latitude>=5);   
!!$     maskbox(ilat,ilon)=true;    
!!$     interpmask_sea=masksea(latitude,longitude);
!!$     interpmask_geopot=maskorog(latitude,longitude);
!!$     maskbox(interpmask_s1024ea<1)=false;
!!$     
!!$     % Andaman Islands
!!$     ilon1=(longitude>=92) & (longitude<=94);
!!$     ilat1=(latitude>=11) & (latitude<=14);
!!$     maskbox(ilat1,ilon1)=false;   
!!$     % Philippines, Borneo, Malaysia
!!$     ilon1=(longitude>=115) & (longitude<=130);
!!$     ilat1=(latitude>=5) & (latitude<=20);
!!$     maskbox(ilat1,ilon1)=false;
!!$     ilon1=(longitude>=94) & (longitude<=103.5);
!!$     ilat1=(latitude>=5) & (latitude<=9.5);
!!$     maskbox(ilat1,ilon1)=false;
!!$     % Lakes
!!$     ilon1=(longitude>=80) & (longitude<=92);
!!$     ilat1=(latitude>=30) & (latitude<=35);
!!$     maskbox(ilat1,ilon1)=true;
!!$     ilon1=(longitude>=103.5) & (longitude<=105);
!!$     ilat1=(latitude>=12) & (latitude<=14);
!!$     maskbox(ilat1,ilon1)=true;
!!$     ilon1=(longitude>=102) & (longitude<=103);
!!$     ilat1=(latitude>=16) & (latitude<=20);
!!$     maskbox(ilat1,ilon1)=true;
!!$     ilon1=(longitude>=115) & (longitude<=117);
!!$     ilat1=(latitude>=28) & (latitude<=32);
!!$     maskbox(ilat1,ilon1)=true;
!!$     ilon1=(longitude>=118) & (longitude<=120);
!!$     ilat1=(latitude>=33) & (latitude<=34);
!!$     maskbox(ilat1,ilon1)=true;
!!$     ilon1=(longitude>=119.5) & (longitude<=120.5);
!!$     ilat1=(latitude>=30.5) & (latitude<=32);
!!$     maskbox(ilat1,ilon1)=true;
!!$     ilon1=(longitude>=75.5) & (longitude<=76.5);
!!$     ilat1=(latitude>=31.5) & (latitude<=32.5);
!!$     maskbox(ilat1,ilon1)=true;
!!$     
!!$     maskbox(interpmask_geopot>=3500.*9.81)=false;
!!$*****************************************************+    

    case ('BoB') ! Bay of Bengal

       where (longitude>=78 .and. longitude<99 .and. &
              latitude<25 .and. latitude>=5)
           maskval = .true.
       end where
       where (longitude>=98 .and. longitude < 100.5 .and. &
              latitude<7.5 .and. latitude >= 5)
           maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       !where (lsm_mask==1)    
       where (lsm_mask/=0)    
          maskval = .false.
       end where

       ! Andaman Islands
       where (longitude>=92 .and. longitude<=94 .and. &
              latitude>=11 .and. latitude<=14)
           maskval = .true.
       end where

    case ('SCSPhi') ! South China Sea & Philippines Sea

       where (longitude>=99 .and. longitude<150 .and. &
              latitude<=23 .and. latitude>=5)
           maskval = .true.
       end where
       where (longitude>=98 .and. longitude<100.5 .and. &
              latitude<7.5 .and. latitude>=5)
           maskval = .false.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       !where (lsm_mask==1)    
       where (lsm_mask/=0)    
          maskval = .false.
       end where
    
       ! Erase Lakes on Pen
       where (longitude>=103.5 .and. longitude<=105 .and. &
              latitude>=12 .and. latitude<=14)
           maskval = .false.
       end where
       where (longitude>=103.5 .and. longitude<=104 .and. &
              latitude>=10.8 .and. latitude<=11.2)
           maskval = .false.
       end where
       where (longitude>=102 .and. longitude<=103 .and. &
              latitude>=16 .and. latitude<=20)
           maskval = .false.
       end where
       where (longitude>=102.4 .and. longitude<=103 .and. &
              latitude>=4.7 .and. latitude<=5.3)
           maskval = .false.
       end where

    case ('Phi') ! Philippines

       where (longitude>=116 .and. longitude<=130 .and. &
              latitude>7.5 .and. latitude<=20)
           maskval = .true.
       end where
       where (longitude>=119.5 .and. longitude<=130 .and. &
              latitude>5 .and. latitude<=20)
           maskval = .true.
       end where
       lsm_mask = set_lsm_mask(latitude, longitude, lsm)
       !where (lsm_mask==0)    
       where (lsm_mask/=1)    
          maskval = .false.
       end where

!!$   case 'NAPO' %North Asia Pacific Ocean
!!$     clear ilon ilat
!!$     ilon = (longitude>=78) & (longitude<=150);
!!$     ilat = (latitude<=35) & (latitude>5);         
!!$     maskbox(ilat,ilon)=true;
!!$     interpmask_sea=masksea(latitude,longitude);
!!$     maskbox(interpmask_sea==1)=false;
!!$     % Andaman Islands
!!$     ilon1=(longitude>=92) & (longitude<=94);
!!$     ilat1=(latitude>=11) & (latitude<=14);
!!$     maskbox(ilat1,ilon1)=true;   
!!$     % Philippines, Borneo, Malaysia
!!$     ilon1=(longitude>=115) & (longitude<=130);
!!$     ilat1=(latitude>5) & (latitude<=20);
!!$     maskbox(ilat1,ilon1)=true;
!!$     ilon1=(longitude>=94) & (longitude<=103.5);
!!$     ilat1=(latitude>5) & (latitude<=9.5);
!!$     maskbox(ilat1,ilon1)=true;

    case ('CAm') ! Central America

       where (longitude<310 .and. longitude>=235 .and. &
              latitude<30 .and. latitude>=0)
           maskval = .true.
       end where
 
    case ('SAm') ! South America

       where (longitude<325 .and. longitude>=279 .and. &
              latitude<0 .and. latitude>=-35)
           maskval = .true.
       end where

    case ('NCP') ! North Central Pacific

       where (longitude>=150 .and. longitude<210 .and. &
              latitude>=5 .and. latitude<=20)
           maskval = .true.
       end where
       where (longitude>=210 .and. longitude<235 .and. &
              latitude>=0 .and. latitude<=20)
           maskval = .true.
       end where

    case ('SEP') ! South East Pacific

       where (longitude>=210 .and. longitude<235 .and. &
              latitude>=-20 .and. latitude<0)
           maskval = .true.
       end where
       where (longitude>=235 .and. longitude<279 .and. &
              latitude>=-20 .and. latitude<0)
           maskval = .true.
       end where

    case ('Tibet') ! Tibetan plateau
       
       where (longitude<=120 .and. longitude>=50 .and. &
              latitude<=50 .and. latitude>=25)
           maskval = .true.
       end where
       oro_mask = set_oro_mask(latitude, longitude, gph)
       where (oro_mask<3500.*9.81) 
          maskval = .false.
       end where

    case ('HimSlope') ! Himalaya slope

       where (longitude<=95 .and. longitude>=72 .and. &
              latitude<=35 .and. latitude>=22)
           maskval = .true.
       end where

       oro_mask = set_oro_mask(latitude, longitude, gph)
       where (oro_mask>=9.81*3500. .or. oro_mask<=9.81*600.) 
          maskval = .false.
       end where

    end select

    ! level below levelbound
    where (level > levelbound) 
       maskval = .false.
    endwhere
    

    deallocate (lsm_mask, oro_mask)

  End Subroutine define_regional_mask


end module messy_clamstracer_regional_masks
