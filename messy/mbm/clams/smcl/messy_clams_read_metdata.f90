Module messy_clams_read_metdata

  USE netcdf

  PRIVATE

  PUBLIC :: clams_read_ecmwf_files
  PUBLIC :: nc_read_corrfile

  INTERFACE calcthetap
     MODULE PROCEDURE calcthetap_theta1d
     MODULE PROCEDURE calcthetap_theta3d
  END INTERFACE calcthetap
    
Contains

!--------------------------------------------------------------------
! Reads winddata from ECMWF data files
!--------------------------------------------------------------------

  SUBROUTINE clams_read_ecmwf_files (status, USE_CLAMSSEDI) 
    ! THIS subroutine belongs in SMIL!!

    USE messy_clams_global,     ONLY: YEAR_START, MONTH_START, DAY_START,     &
                                      HOUR_START, MINUTE_START, SECOND_START, &
                                      YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, &
                                      lfirst_cycle, lstart, lresume

    USE messy_clams_global,     ONLY: rank, prec, dp, &
                                      met_dir, met_prefix, met_freq, &
                                      theta_dir, theta_prefix, &
                                      predata, futdata, irdday, irdsec, &
                                      nx, ny, nz, &
                                      pre_year_co, pre_month_co, &
                                      pre_day_co, pre_sec_co, &
                                      pre_metfile, fut_metfile, &
                                      pre_thetafile, fut_thetafile, &
                                      fut_year, fut_month, fut_day, fut_sec, &
                                      pre_year, pre_month, pre_day, pre_sec, &
                                      udt, vdt, wdt, ufut, vfut, wfut, &
                                      leveldt, levelfut, dlevdzdt, dlevdzfut, &
                                      nparams, levdotname, use_3d, buffersize, &
                                      filenamelen
    USE messy_clamssedi_global, ONLY: UDT_sedi, VDT_sedi, WDT_sedi, &
                                      leveldt_sedi, dlevdzdt_sedi
    
    USE messy_clams_tools_dateconv,  ONLY: incdat, datsec
    USE messy_clams_tools_utils,     ONLY: get_metfilename
 
    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL

    INTEGER                 :: status
    LOGICAL, INTENT(IN)     :: USE_CLAMSSEDI

    CHARACTER(LEN=*), PARAMETER :: substr='clams_read_ecmwf_files'

    character(40)               :: data_orig=' ', data_desc=' '

    INTEGER                 :: ncuv
    INTEGER                 :: rcode
    INTEGER                 :: sod
    CHARACTER(filenamelen)  :: srcfile
    INTEGER                 :: pos, i
 
    status = 0 ! no error 

    !**************************************
    ! Read WINDFILES
    !**************************************
    IF (met_freq<1 .or. met_freq>24) then
       status = 6
       return
    endif
    SELECT CASE(met_freq)
    CASE (24)
       irdday = 1
       irdsec = 0
    CASE (1,2,3,4,6,8,12)
       irdday = 0
       irdsec = met_freq*3600
    CASE DEFAULT
       WRITE (*,*) 'Frequency of windfiles must be 1,2,3,4,6,8,12 or 24  !!!' 
       status = 1
       return
    END SELECT

    IF (lfirst_cycle) THEN
       !---------------------------------
       ! First windfile (at start time or restart time)
       !---------------------------------
       if(lstart) then
          if (rank==0) write(*,*) 'lstart in ecmwf reading'
          pre_year  = YEAR_START
          pre_month = MONTH_START
          pre_day   = DAY_START
          pre_sec   = HOUR_START*3600 + MINUTE_START*60 + SECOND_START
       else !restart
          if (rank==0) then
             write(*,*) 'pre_year, pre_month, pre_day, pre_sec'
             write(*,*) pre_year, pre_month, pre_day, pre_sec
          endif
       end if
          
       !---------------------------------------------------
       ! Get data_description and data_origin
       !---------------------------------------------------
       rcode = nf90_open(pre_metfile, NF90_NOWRITE, ncuv, buffersize)
       if (rcode /= 0) then
          WRITE (*,*) 'Wind file could not be opened ! '
          WRITE (*,*) 'The program will stop.'
          status = -1
          return
       ENDIF
       rcode = nf90_get_att (ncuv, NF90_GLOBAL, 'data_description', data_desc)
       if (rcode /= 0) then
          rcode = nf90_get_att (ncuv, NF90_GLOBAL, 'source_file_name', srcfile)
          pos = index(srcfile,'ppass')
          if (pos/=0) then
             data_desc = 'UKMO'
          else
             data_desc = 'noname'
          endif
       endif
       rcode = nf90_get_att (ncuv, NF90_GLOBAL, 'data_origin', data_orig)
       if (rcode /= 0)   data_orig = 'noname'
       rcode = nf90_close(ncuv)
       IF(rank==0)THEN
          write (*,*) 'Data origin: ', TRIM(data_orig)
          write (*,*) 'Data description: ', TRIM(data_desc)
          write (*,*) 
       ENDIF

        ! 3d or 2d
       IF (rank==0) THEN
          WRITE(*,*) '3D run: ',use_3d
          WRITE(*,*) 
       ENDIF

      ! short name of variable containing THETA-DOT information
       IF(rank==0) THEN
          WRITE (*,*) 'variable containing lev-dot information: ', levdotname
          write (*,*) 
       ENDIF

       IF (lstart) THEN

          CALL nc_readday (status, PREDATA, udt, vdt, wdt, leveldt, &
               pre_metfile, pre_thetafile, &
               pre_year, pre_month,  pre_day, pre_sec/3600)
          IF (status /= 0) RETURN
          
          dlevdzdt = -1

          IF (USE_CLAMSSEDI) THEN
             udt_sedi = udt
             vdt_sedi = vdt
             wdt_sedi = wdt
             leveldt_sedi = leveldt
             call get_dlevdz (status, dlevdzdt_sedi, leveldt_sedi, pre_metfile)
             
             IF (status /= 0) RETURN
          ENDIF
       
       END IF

       !---------------------------------
       ! Next windfile
       !---------------------------------
       ! Calculate date for second windfile:
       
       fut_year = pre_year
       fut_month= pre_month
       fut_day = pre_day
       fut_sec = pre_sec

       CALL incdat (fut_sec, fut_day, fut_month, fut_year,  &
            irdsec,irdday,0,0)              

       ! Name of second windfile:
       fut_metfile = get_metfilename (met_prefix, met_dir, &
             fut_year, fut_month, fut_day, fut_sec/3600)
       if (theta_prefix=='') then
          fut_thetafile = ''
       else
          fut_thetafile = get_metfilename (theta_prefix, theta_dir, &
             fut_year, fut_month, fut_day, fut_sec/3600)
       endif

       IF (lstart) THEN

          CALL nc_readday (status, FUTDATA, ufut, vfut, wfut, levelfut, &
               fut_metfile, fut_thetafile, &
               fut_year, fut_month,  fut_day, fut_sec/3600)
         
          IF (USE_CLAMSSEDI)  THEN
             call get_dlevdz (status, dlevdzfut, levelfut, fut_metfile)
             IF (status /= 0) RETURN
          ELSE
             dlevdzfut = -1.
          ENDIF

       END IF

    ELSE ! not first cycle

       sod = HOUR*3600 +MINUTE*60 + SECOND

       IF (datsec(sod, DAY, MONTH, YEAR) ==  &
            datsec(fut_sec, fut_day, fut_month, fut_year)) THEN  

          if (rank==0) write(*,*) 
          if (rank==0) write(*,*) 'Next windfile reached'
!!!!!
!          write(140+rank,'(I,I,I,I)') fut_year,fut_month,fut_day,fut_sec
          pre_year  = fut_year
          pre_month = fut_month 
          pre_day   = fut_day
          pre_sec   = fut_sec
          pre_year_co  = fut_year
          pre_month_co = fut_month 
          pre_day_co   = fut_day
          pre_sec_co   = fut_sec

          CALL incdat (fut_sec, fut_day, fut_month, fut_year,  &
               irdsec,irdday,0,0)              

          pre_metfile = fut_metfile
          pre_thetafile = fut_thetafile

          DO i = 1, nparams
             PREDATA(i)%values = FUTDATA(i)%values
          ENDDO
 
          fut_metfile = get_metfilename (met_prefix, met_dir, &
               fut_year, fut_month, fut_day, fut_sec/3600)
          if (theta_prefix == '') then
             fut_thetafile = ''
          else
             fut_thetafile = get_metfilename (theta_prefix, theta_dir, &
                  fut_year, fut_month, fut_day, fut_sec/3600)
          endif

          if (rank==0) then
             write(*,*) 'pre_date:',pre_year,pre_month,pre_day,pre_sec
             write(*,*) 'fut_date:',fut_year,fut_month,fut_day,fut_sec
             write(*,*) 'pre_metfile: ', trim(pre_metfile)
             write(*,*) 'fut_metfile: ', trim(fut_metfile)
             if (pre_thetafile/='') then
                write(*,*) 'pre_thetafile: ', trim(pre_thetafile)
                write(*,*) 'fut_thetafile: ', trim(fut_thetafile)
             endif
          endif

          CALL nc_readday (status, FUTDATA, ufut, vfut, wfut, levelfut, &
               fut_metfile, fut_thetafile, &
               fut_year, fut_month,  fut_day, fut_sec/3600)
         

          IF (USE_CLAMSSEDI) THEN
             call get_dlevdz (status, dlevdzfut, levelfut, fut_metfile)
             IF (status /= 0) RETURN
          ELSE
             dlevdzfut = -1.
          ENDIF

       ENDIF ! next windfile?
    ENDIF ! lfirst_cycle?

    IF (lresume) THEN
       sod = HOUR*3600 +MINUTE*60 + SECOND

       IF (datsec(sod, DAY, MONTH, YEAR) ==  &
            datsec(fut_sec, fut_day, fut_month, fut_year)) THEN  

          if (rank==0) write(*,*) 
          if (rank==0) write(*,*) 'Next windfile reached'
!!!!!
!          write(140+rank,'(I,I,I,I)') fut_year,fut_month,fut_day,fut_sec
          pre_year  = fut_year
          pre_month = fut_month 
          pre_day   = fut_day
          pre_sec   = fut_sec
          pre_year_co  = fut_year
          pre_month_co = fut_month 
          pre_day_co   = fut_day
          pre_sec_co   = fut_sec

          CALL incdat (fut_sec, fut_day, fut_month, fut_year,  &
               irdsec,irdday,0,0)              

          pre_metfile = fut_metfile
          pre_thetafile = fut_thetafile

          DO i = 1, nparams
             PREDATA(i)%values = FUTDATA(i)%values
          ENDDO
 
          fut_metfile = get_metfilename (met_prefix, met_dir, &
               fut_year, fut_month, fut_day, fut_sec/3600)
          if (theta_prefix == '') then
             fut_thetafile = ''
          else
             fut_thetafile = get_metfilename (theta_prefix, theta_dir, &
                  fut_year, fut_month, fut_day, fut_sec/3600)
          endif

          if (rank==0) then
             write(*,*) 'pre_metfile: ', trim(pre_metfile)
             write(*,*) 'fut_metfile: ', trim(fut_metfile)
             if (pre_thetafile/='') then
                write(*,*) 'pre_thetafile: ', trim(pre_thetafile)
                write(*,*) 'fut_thetafile: ', trim(fut_thetafile)
             endif
          endif

          CALL nc_readday (status, FUTDATA, ufut, vfut, wfut, levelfut, &
               fut_metfile, fut_thetafile, &
               fut_year, fut_month,  fut_day, fut_sec/3600)
          
          IF (USE_CLAMSSEDI)  THEN
             call get_dlevdz (status, dlevdzfut, levelfut, fut_metfile)
             if (status /= 0) return
          ELSE
             dlevdzfut = -1.
          ENDIF

       END IF
    END IF

  END SUBROUTINE clams_read_ecmwf_files

  !*********************************************************************                  
  ! SUBROUTINE nc_readday(status, data, uwind, vwind, wwind, levdata, file)
  !
  ! subroutine to get data of all parameters for one time and all
  ! levels from the netCDF- file
  !
  !**********************************************************************
  SUBROUTINE nc_readday (status, data, u, v, w, lev, file, thetafile, &
                         year, month,  day, hour)
    
    USE messy_clams_global,         ONLY: prec, dp, rank, datatype, &
                                          nx, ny, nz, nparams, &
                                          use_3d, level_is_vertcoor, levelgrid, levdotname, &
                                          mdi, eps, buffersize, &
                                          init_vertcoorname, met_vertcoorname
    USE messy_clams_tools_utils,    ONLY: uppercase, lowercase
    USE messy_clams_tools_ncutils,  ONLY: nc_check_error, nc_get_param_2dint
    USE messy_clams_tools_dateconv, ONLY: ymds2js

    IMPLICIT NONE

    TYPE(datatype), dimension(:), pointer :: data
    INTEGER     :: status
    character(*):: file, thetafile

    REAL(PREC) :: u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz),lev(nx,ny,nz)
    INTEGER    :: year,month,day,hour


    !Local Variables
    real(prec),  dimension(:,:,:), allocatable :: param3d
    integer,     dimension(:,:),   allocatable :: indarr
    integer                                    :: nin, nin2, varid
    integer                                    :: pos
    integer                                    :: iparam, ix, iy, iz
    real(dp)                                   :: jultime

    
    status = 0

    !open actual netCDF- file
    status = nf90_open (file,NF90_NOWRITE,nin, buffersize)
    call nc_check_error (status,'Can not open file '//TRIM(file)//' !',abort=.false.)
    if (status /= nf90_noerr) return

    if (thetafile /= '') then
       status = nf90_open (thetafile,nf90_nowrite,nin2, buffersize)
       call nc_check_error (status,'Can not open file '//TRIM(thetafile)//' !',abort=.false.)
       if (status /= nf90_noerr) return
    endif

    ! Read U
    call get_species (status,nin,'U',u)

    ! Read V
    call get_species (status,nin,'V',v)

    ! model level data: read vertical coordinate (used in init-file)
    if (.not. level_is_vertcoor) then
       call get_species (status,nin,trim(uppercase(init_vertcoorname)),lev)
    else
       lev = 0.0
    endif

    ! Read W (THETA_DOT/ZETA_DOT/OMEGA)
    if (use_3d) then
       call get_species (status,nin,trim(levdotname),w)
       ! change units 
       if (trim(init_vertcoorname)=='press') then
          ! changes units of OMEGA from Pa/s to hPa/s 
          where (ABS((w-MDI)/MDI)>eps)
             w = w/100.
          endwhere
       elseif (.not. level_is_vertcoor) then
          jultime = ymds2js (year,month,day,hour*3600)
          DO iz=1,nz
             call calcthetap(w(1:nx,1:ny,iz),nx,ny,lev(1:nx,1:ny,iz),jultime)    
          ENDDO
       else
          jultime = ymds2js (year,month,day,hour*3600)
          DO iz=1,nz
             call calcthetap(w(1:nx,1:ny,iz),nx,ny,levelgrid(iz),jultime)    
          ENDDO
       endif
    else  ! set d theta/dt=0
       w = 0.0
    endif

    ! Read other output parameters
    do iparam = 1, nparams

       ! if parameter is vertical coordinate in windfile:
       if (data(iparam)%name == uppercase(met_vertcoorname)) then
          call get_vertcoor3d (status,nin,trim(met_vertcoorname),data(iparam)%values)
          
       else
!!$       if (trim(data(iparam)%name)=='SH_REPLACE') then
!!$          call get_species (status,nin,'SH',data(iparam)%values)
!!$       else
          
          ! read EQLAT and PV from isentropic file
          if (data(iparam)%name=='EQLAT' .or. data(iparam)%name=='PV') then
             call get_species (status, nin2, trim(data(iparam)%name), data(iparam)%values)
             if (status /= nf90_noerr) return
             
          else
             
             pos = INDEX(data(iparam)%name,'_',back=.true.)

             ! no underscore in parameter name: read parameter
             if (pos == 0) then
                call get_species (status,nin,trim(data(iparam)%name),data(iparam)%values)
                if (status /= nf90_noerr) return

             ! underscore in parameter name, but suffix not "TROP1" or "TROP2": read parameter
             elseif (data(iparam)%name(pos+1:)/='TROP1' .and. data(iparam)%name(pos+1:)/='TROP2') then
                call get_species (status,nin,trim(data(iparam)%name),data(iparam)%values)
                if (status /= nf90_noerr) return

             ! suffix of parameter name "TROP1" or "TROP2"
             else

                allocate (param3d(nx,ny,nz))
                allocate (indarr(nx,ny))

                ! read parameter
                call get_species (status,nin,trim(data(iparam)%name(1:pos-1)),data(iparam)%values)
                if (status /= nf90_noerr) return
                param3d = data(iparam)%values

                ! read TROP1/TROP2
                call nc_get_param_2dint (status,nin,trim(data(iparam)%name(pos+1:)),indarr)
                if (status /= nf90_noerr) return

                ! get parameter values on TROP1 or TROP2 and copy to every level
                ! (for following interpolation routines)
                do iy = 1, ny
                   do ix = 1, nx
                      if (indarr(ix,iy)<=0) then
                         data(iparam)%values(ix,iy,:) = mdi
                      else
                         data(iparam)%values(ix,iy,:) = param3d(ix,iy,indarr(ix,iy))
                      endif
                   enddo
                enddo

                deallocate (param3d, indarr)
                
             endif
             
          endif

       endif

    enddo

    status = nf90_close(nin)

    if (thetafile /= '') status = nf90_close (nin2)
  
 end subroutine nc_readday


  !*************************************************************
  !
  !*************************************************************
  SUBROUTINE get_dlevdz (status, dlevdz, level, metfile)

    USE messy_clams_global,        ONLY: prec, mdi, eps, nx, ny, nz, &
                                         level_is_vertcoor, levelgrid, rank, &
                                         buffersize
    USE messy_clams_tools_ncutils, ONLY: nc_check_error
    
    IMPLICIT NONE
    
    real(prec), dimension(:,:,:), pointer :: dlevdz, level
    integer      :: status
    character(*) :: metfile
    
    
    real(prec), dimension(:,:,:), pointer :: gph
    real(prec) :: dlev, dz
    integer    :: ncid, rcode, varid
    integer    :: ix, iy, iz
    
    
    status = 0

    ! open meteorological file
    status = nf90_open (metfile,NF90_NOWRITE,ncid, buffersize)
    call nc_check_error (status,'Can not open file '//TRIM(metfile)//' !',abort=.false.)
    if (status /= nf90_noerr) return
    
    ! read GPH
    allocate (gph(nx,ny,nz))
    call get_species(status, ncid, 'GPH', gph)
    if (status /= 0) then
       deallocate (gph)
       return
    endif
     
    ! calculate dlev / dz 
    do iz = 2, nz-1
       do iy = 1, ny
          do ix = 1, nx
             if (ABS((gph(ix,iy,iz+1)-mdi)/mdi)<=eps .or. &
                 ABS((gph(ix,iy,iz-1)-mdi)/mdi)<=eps) then
                dlevdz(ix,iy,iz) = mdi
             else
                if (.not. level_is_vertcoor) then
                   dlev = level(ix,iy,iz+1)-level(ix,iy,iz-1)
                else
                   dlev = levelgrid(iz+1) - levelgrid(iz-1)
                endif
                dz = gph(ix,iy,iz+1)/9.80665 - gph(ix,iy,iz-1)/9.80665
                dlevdz(ix,iy,iz) = dlev / dz
             endif
          enddo
       enddo
    enddo
    dlevdz(:,:,1)  = mdi
    dlevdz(:,:,nz) = mdi

    ! close met. file
    rcode = nf90_close (ncid)

    deallocate (gph)

  END SUBROUTINE get_dlevdz

  !*********************************************************************                  
  !
  !*********************************************************************                  
  subroutine get_species(status, nin, varname, species)

    USE messy_clams_global,          ONLY: prec
    USE messy_clams_tools_ncutils,   ONLY: nc_check_error, nc_get_var_cf

    IMPLICIT NONE
    
    integer                       :: status,nin
    character(*)                  :: varname         
    real(prec), dimension(:,:,:)  :: species

    integer                       :: varid

    status = 0
    
    status = nf90_inq_varid (nin, trim(varname), varid)
    call nc_check_error (status, &
         'ERROR - '//trim(varname)//' not found in windfile !!!', abort=.false.)
    if (status /= nf90_noerr) return

    call nc_get_var_cf (status, nin, varid, varname, species)


  end subroutine get_species

  !*************************************************************
  !
  !*************************************************************                  
  subroutine get_vertcoor3d (status, nin,varname,species)

    use netcdf

    use messy_clams_global,         only: prec, nz
    use messy_clams_tools_ncutils,  only: nc_check_error
    
    implicit none
    
    real(prec), dimension(:,:,:) :: species
    integer                      :: status, nin
    character(*)                 :: varname         

    real, dimension(nz) :: values1d
    integer             :: varid
    integer             :: iz
    
    status = 0

    status = nf90_inq_varid(nin,trim(varname),varid)
    call nc_check_error (status, &
         'ERROR - '//trim(varname)//' not found in windfile !!!', abort=.false.)
    if (status /= nf90_noerr) return

    status = nf90_get_var(nin,varid,values1d)
    call nc_check_error (status, &
         'Coordinate variable '//trim(varname)//' could not be read !!!', abort=.false.)
    if (status /= nf90_noerr) return

    do iz = 1, nz
       species(:,:,iz) = values1d(iz)
    enddo
        
  end subroutine get_vertcoor3d

!**************************************************************
! - changes units of theta_dot from K/day to K/s 
! - add correction (if correction file is specified)
!**************************************************************
SUBROUTINE calcthetap_theta1d (theta_dot,nlon,nlat,theta,jultime)

  use messy_clams_global,     only: prec, dp, mdi, eps, latgrid, &
                                    corr_thetadot,corr_3d

  implicit none

  integer     :: nlon, nlat
  real(prec)  :: theta_dot(nlon,nlat),theta
  real(dp)    :: jultime

  ! local variables
  integer    :: ilat
  real(prec) :: pi, corr
  
  pi=4.*atan(1.)                                                        

  ! if no correction file is specified:
  if (.not. corr_thetadot) then
     where (ABS((theta_dot-mdi)/mdi)>eps) 
        theta_dot = theta_dot/86400. 
     end where
     ! if a correction file is specified: add correction
  else
     ! if CORR is two-dimensional (level,time)
     if (.not. corr_3d) then
        corr = get_thetadot_corr (theta,jultime)    
        where (ABS((theta_dot-mdi)/mdi)>eps) 
           theta_dot = theta_dot/86400. + corr/86400.
        end where
        ! if CORR is three-dimensional (level,lat,time)
     else
        do ilat = 1, nlat
           corr = get_thetadot_corr3d (theta,jultime,latgrid(ilat))
           where (ABS((theta_dot(:,ilat)-mdi)/mdi)>eps) 
              theta_dot(:,ilat) = theta_dot(:,ilat)/86400. + corr/86400.
           end where
        end do
     end if
  end if
  
END SUBROUTINE calcthetap_theta1d

SUBROUTINE calcthetap_theta3d (theta_dot,nlon,nlat,theta,jultime)

  use messy_clams_global,     only: prec, dp, mdi, eps, latgrid, &
                                    corr_thetadot,corr_3d

  implicit none

  integer      :: nlon, nlat
  real(prec)   :: theta_dot(nlon,nlat),theta(nlon,nlat)
  real(dp)     :: jultime

  ! local variables
  integer     :: ilat, ilon
  real(prec)  :: pi, corr
  
  pi=4.*atan(1.)                                                        

  ! if no correction file is specified:
  if (.not. corr_thetadot) then
     where (ABS((theta_dot-mdi)/mdi)>eps) 
        theta_dot = theta_dot/86400. 
     end where
     ! if a correction file is specified: add correction
  else
     ! if CORR is two-dimensional (level,time)
     if (.not. corr_3d) then
        do ilat = 1, nlat
           do ilon = 1, nlon
              corr = get_thetadot_corr (theta(ilon,ilat),jultime)    
              if (ABS((theta_dot(ilon,ilat)-mdi)/mdi)>eps) then
                 theta_dot(ilon,ilat) = theta_dot(ilon,ilat)/86400. + corr/86400.
              endif
           enddo
        enddo

     ! if CORR is three-dimensional (level,lat,time)
     else
        do ilat = 1, nlat
           do ilon = 1, nlon
              corr = get_thetadot_corr3d (theta(ilon,ilat),jultime,latgrid(ilat))
              if (ABS((theta_dot(ilon,ilat)-mdi)/mdi)>eps) then
                 theta_dot(ilon,ilat) = theta_dot(ilon,ilat)/86400. + corr/86400.
              endif
           enddo
        enddo
     end if
  end if
  
END SUBROUTINE calcthetap_theta3d


!*************************************************************
!
!*************************************************************
function get_thetadot_corr (theta,jultime)
 
  use messy_clams_global,         only: dp, prec, &
                                        theta_corr, time_corr, thetadot_corr
  use messy_clams_tools_interpol, only: interpol_lin, interpol_lin_double

  implicit none

  real(prec)  :: get_thetadot_corr, theta
  real(dp)    :: jultime
  
  ! local variables:
  integer    :: ntheta, ntime, i, itheta1, itheta2, itime1, itime2
  real(prec) :: val1, val2, corr
  
  ntheta = size(theta_corr)
  ntime = size(time_corr)

  ! get theta indices for interpolation
  if (theta <= theta_corr(1)) then
     itheta1 = 1
     itheta2 = 1
  elseif (theta >= theta_corr(ntheta)) then 
     itheta1 = ntheta
     itheta2 = ntheta
  else
     do i=1,ntheta-1
        if (theta_corr(i)<=theta .and. theta<=theta_corr(i+1)) itheta1=i 
     end do
     itheta2 = itheta1+1
  endif
  
  ! get time indices for interpolation
  if (jultime <= time_corr(1)) then
     itime1 = 1
     itime2 = 1
  elseif (jultime >= time_corr(ntime)) then
     itime1 = ntime
     itime2 = ntime
  else
     do i=1,ntime-1
        if (time_corr(i)<=jultime .and. jultime<=time_corr(i+1)) itime1=i
     end do
     itime2 = itime1+1
  end if
  
  ! interpolate between theta levels for first time
  val1 = interpol_lin (theta_corr(itheta1),thetadot_corr(itime1,  itheta1), &
                       theta_corr(itheta2),thetadot_corr(itime1,  itheta2),theta)

  ! interpolate between theta levels for second time
  val2 = interpol_lin (theta_corr(itheta1),thetadot_corr(itime2,itheta1), &
                       theta_corr(itheta2),thetadot_corr(itime2,itheta2),theta)
  
  ! interpolate between two times
  corr = interpol_lin_double (time_corr(itime1),val1,time_corr(itime2),val2,jultime)

  get_thetadot_corr = corr

end function get_thetadot_corr


function get_thetadot_corr3d (theta,jultime,lat)

  use messy_clams_global,         only: dp, prec, &
                                        theta_corr, time_corr, lat_corr, thetadot_corr3d
  use messy_clams_tools_interpol, only: interpol_lin, interpol_lin_double

  implicit none

  real(prec)  :: get_thetadot_corr3d, theta, lat
  real(dp)    :: jultime
  
  ! local variables:
  integer    :: ntheta, ntime, nlat, i, itheta1, itheta2, ilat1, ilat2, itime1, itime2
  real(prec) :: val_lat1, val_lat2, val_time1, val_time2, corr
  
  ntheta = size(theta_corr)
  ntime  = size(time_corr)
  nlat   = size(lat_corr)
  
  ! get theta indices for interpolation
  if (theta <= theta_corr(1)) then
     itheta1 = 1
     itheta2 = 1
  elseif (theta >= theta_corr(ntheta)) then
     itheta1 = ntheta
     itheta2 = ntheta
  else
     do i=1,ntheta-1
        if (theta_corr(i)<=theta .and. theta<=theta_corr(i+1)) itheta1=i 
     end do
     itheta2 = itheta1+1
  endif
     
  ! get lat indices for interpolation
  if (lat <= lat_corr(1)) then
     ilat1 = 1
     ilat2 = 1
  elseif (lat >= lat_corr(nlat)) then
     ilat1 = nlat
     ilat2 = nlat
  else
     do i=1,nlat-1
        if (lat_corr(i)<=lat .and. lat<=lat_corr(i+1)) ilat1=i
     end do
     ilat2 = ilat1+1
  endif
  
  ! get time indices for interpolation
  itime1=0
  if (jultime <= time_corr(1)) then
     itime1 = 1
     itime2 = 1
  elseif (jultime >= time_corr(ntime)) then
     itime1 = ntime
     itime2 = ntime
  else
     do i=1,ntime-1
        if (time_corr(i)<=jultime .and. jultime<=time_corr(i+1)) itime1=i
     end do
     itime2 = itime1+1
  end if

  ! interpolate vertical (between itheta1 and itheta2) for first latitude and first time
  val_lat1 = interpol_lin (theta_corr(itheta1),thetadot_corr3d(itime1,ilat1,itheta1), &
                           theta_corr(itheta2),thetadot_corr3d(itime1,ilat1,itheta2),theta)
  
  ! interpolate vertical (between itheta1 and itheta2) for second latitude and first time
  val_lat2 = interpol_lin (theta_corr(itheta1),thetadot_corr3d(itime1,ilat2,itheta1), &
                           theta_corr(itheta2),thetadot_corr3d(itime1,ilat2,itheta2),theta)

  ! interpolate horizontal (between lat1 and lat2) for first time
  val_time1 = interpol_lin (lat_corr(ilat1),val_lat1,lat_corr(ilat2),val_lat2,lat)


  ! interpolate vertical (between itheta1 and itheta2) for first latitude and second time
  val_lat1 = interpol_lin (theta_corr(itheta1),thetadot_corr3d(itime2,ilat1,itheta1), &
                           theta_corr(itheta2),thetadot_corr3d(itime2,ilat1,itheta2),theta)

  ! interpolate vertical (between itheta1 and itheta2) for second latitude and second time
  val_lat2 = interpol_lin (theta_corr(itheta1),thetadot_corr3d(itime2,ilat2,itheta1), &
                           theta_corr(itheta2),thetadot_corr3d(itime2,ilat2,itheta2),theta)

  ! interpolate horizontal (between lat1 and lat2) for second time
  val_time2 = interpol_lin (lat_corr(ilat1),val_lat1,lat_corr(ilat2),val_lat2,lat)

  ! interpolate time 
  corr = interpol_lin_double (time_corr(itime1),val_time1,time_corr(itime2),val_time2,jultime)


  get_thetadot_corr3d = corr

end function get_thetadot_corr3d

!**************************************************************
! - changes units of OMEGA from Pa/s to hPa/s 
!**************************************************************
!!!!! wird nicht mehr genutzt !?!
! SUBROUTINE calc_omega (P,LEN,press)
! use messy_clams_global, only: mdi,eps,prec
! IMPLICIT NONE

! INTEGER    :: LEN,I
! REAL(PREC) :: P(LEN),press

! DO i=1,LEN
!    if (ABS((P(i)-MDI)/MDI)>eps) then
!       P(i)=P(i)/100. 
!    endif
! ENDDO

! END SUBROUTINE calc_omega


  !*********************************************************************
  !
  !*********************************************************************                  
  subroutine nc_read_corrfile (status,filename,vertcoorname)

    USE netcdf
    use messy_clams_global,        only: time_corr,theta_corr,lat_corr, &
                                         thetadot_corr,thetadot_corr3d,corr_3d, &
                                         buffersize
    use messy_clams_tools_ncutils, only: nc_check_error
    use messy_clams_tools_utils,   only: lowercase

    implicit none

    integer      :: status
    character(*) :: filename, vertcoorname

    integer :: ncid, dimid, varid
    integer :: ntime, nlevel, nlat

    status = 0 ! no error

    ! open netcdf file
    status = nf90_open (filename,nf90_nowrite,ncid, buffersize)
    call nc_check_error (status,'Can not open file '//TRIM(filename)//' !',abort=.false.)
    if (status/=0) return
      
    !-------------------------------------------------------
    ! read dimensions
    !-------------------------------------------------------

    ! read time
    status = nf90_inq_dimid (ncid,'time',dimid)
    call nc_check_error (status,'Cannot find dimension time',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=ntime)
    call nc_check_error (status,'Cannot read dimension time',abort=.false.)
    if (status/=0) return

    ! read level
    status = nf90_inq_dimid (ncid,trim(lowercase(vertcoorname)),dimid)
    call nc_check_error (status,'Cannot find dimension '//trim(lowercase(vertcoorname)),abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=nlevel)
    call nc_check_error (status,'Cannot read dimension '//trim(lowercase(vertcoorname)),abort=.false.)
    if (status/=0) return

    ! read lat (if existing)
    status = nf90_inq_dimid (ncid,'lat',dimid)
    if (status /= nf90_noerr) then
       corr_3d = .false.
    else
       status = nf90_inquire_dimension (ncid,dimid,len=nlat)
       call nc_check_error (status,'Cannot read dimension lat',abort=.false.)
       if (status/=0) return
       corr_3d = .true.
    end if

    !-------------------------------------------------------
    ! read variables 
    !-------------------------------------------------------

    ! allocate arrays
    allocate (time_corr(ntime))
    allocate (theta_corr(nlevel))
    if (corr_3d) then
       allocate (lat_corr(nlat))
       allocate (thetadot_corr3d(ntime,nlat,nlevel))
    else
       allocate (thetadot_corr(ntime,nlevel))
    endif

    ! read coordinate time
    status = nf90_inq_varid (ncid,'time',varid)
    call nc_check_error (status,'Cannot find variable time',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid,varid,time_corr)
    call nc_check_error (status,'Cannot read variable time',abort=.false.)
    if (status/=0) return

    ! read coordinate level
    status = nf90_inq_varid (ncid,trim(lowercase(vertcoorname)),varid)
    call nc_check_error (status,'Cannot find variable '//trim(lowercase(vertcoorname)),abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid,varid,theta_corr)
    call nc_check_error (status,'Cannot read variable '//trim(lowercase(vertcoorname)),abort=.false.)
    if (status/=0) return

    ! read coordinate lat 
    if (corr_3d) then
       status = nf90_inq_varid (ncid,'lat',varid)
       call nc_check_error (status,'Cannot find variable lat',abort=.false.)
       if (status/=0) return
       status = nf90_get_var (ncid,varid,lat_corr)
       call nc_check_error (status,'Cannot read variable lat',abort=.false.)
       if (status/=0) return
    endif

    ! read variable CORR
    status = nf90_inq_varid (ncid,'CORR',varid)
    call nc_check_error (status,'Cannot find variable CORR',abort=.false.)
    if (status/=0) return
    if (corr_3d) then
       status = nf90_get_var (ncid,varid,thetadot_corr3d)
     else
       status = nf90_get_var (ncid,varid,thetadot_corr)
    endif
    call nc_check_error (status,'Cannot read variable CORR',abort=.false.)
    if (status/=0) return

    ! close netcdf file
    status = nf90_close (ncid)

  end subroutine nc_read_corrfile

End Module messy_clams_read_metdata
