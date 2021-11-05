      SUBROUTINE WRITE_BGCMEAN_2D(kpie,kpje,kpke,pddpo)

!****************************************************************
!
!**** *WRITE_BGCMEAN_2D* - calculate and write 2-dimensional bgc mean data.
!
!     Patrick Wetzel,    *MPI-Met, HH*    24.05.04
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Write bgc mean data.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *WRITE_BGCMEAN_2D(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
!                  ptiestu) *
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL :: of model grid.
!     *INTEGER* *kpje*    - 2nd REAL :: of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL :: of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd REAL ::) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st REAL ::) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd REAL ::) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!     *REAL*    *ptiestu* - depth of layers [m].
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

      use mo_commo1, only: lyears,lmonts,ldays,ldtdayc,dt

      USE mo_bgcmean
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      use mo_param1_bgc 

      USE mo_control_bgc

      use mo_parallel

      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,nk,ndtday
      REAL :: pddpo(kpie,kpje,kpke)

      real :: time2d,t_hour


#ifndef __cpl_co2
#ifdef PANTHROPOCO2
#ifdef DIFFAT
      REAL :: emission_p_s(nmeantime_2d)
#endif
#endif
#endif


      INCLUDE 'netcdf.inc'
      INTEGER :: ncvarid,ncstat,ncoldmod


      INTEGER :: START_1D,COUNT_1D       ! size of internal array 
      INTEGER :: START_2D(3),COUNT_2D(3) ! size of internal array 
      INTEGER :: START_3D(4),COUNT_3D(4) ! size of internal array 

      
!-----------------------------------------------------------------------

      START_1D =meantime_2d
      COUNT_1D =1

      START_2D(1) =1  
      START_2D(2) =1
      START_2D(3) =meantime_2d
      COUNT_2D(1) =kpie  
      COUNT_2D(2) =kpje
      COUNT_2D(3) =1

      START_3D(1) =1  
      START_3D(2) =1
      START_3D(3) =1
      START_3D(4) =meantime_2d
      COUNT_3D(1) =kpie  
      COUNT_3D(2) =kpje
      COUNT_3D(3) =kwrbioz
      COUNT_3D(4) =1      
      
!-----------------------------------------------------------------------


!
!  Masking bgcmean data.
!
!-----------------------------------------------------------------------
      
      ! atm(l=12:14) is not masked !!
      DO j=1,kpje
      DO i=1,kpie

        IF(pddpo(i,j,1) .LT. 0.5) THEN

          bgcm2d(i,j,jkwco2  )=rmasko      
          bgcm2d(i,j,jpco2   )=rmasko

!!$#ifdef PCOMPONENT_ANALYSIS          
!!$        bgcm2d(i,j,jdpco2_dalk)=rmasko
!!$        bgcm2d(i,j,jdpco2_ddic)=rmasko
!!$        bgcm2d(i,j,jdpco2_dsst)=rmasko
!!$        bgcm2d(i,j,jdpco2_dsss)=rmasko
!!$#endif /* PCOMPONENT_ANALYSIS */      

          bgcm2d(i,j,jdms)    =rmasko       
          bgcm2d(i,j,jdmsflux)=rmasko      
          bgcm2d(i,j,jdmsprod)=rmasko      
          bgcm2d(i,j,jdms_bac)=rmasko      
          bgcm2d(i,j,jdms_uv )=rmasko      
          bgcm2d(i,j,jco2fxd )=rmasko      
          bgcm2d(i,j,jco2fxu )=rmasko      
          bgcm2d(i,j,joxflux )=rmasko      
          bgcm2d(i,j,jniflux )=rmasko      

!js    add jexport etc

          bgcm2d(i,j,jexport  )=rmasko
          bgcm2d(i,j,jexpoca  )=rmasko
          bgcm2d(i,j,jexposi  )=rmasko

#ifdef AGG         
          bgcm2d(i,j,jnos     )=rmasko       
#endif        

#ifdef ANTC14
          bgcm2d(i,j,jac14fx  )=rmasko
#endif       

#ifdef PCFC
          bgcm2d(i,j,jcfc11fx )=rmasko
        bgcm2d(i,j,jcfc12fx )=rmasko
          bgcm2d(i,j,jpcfc11  )=rmasko
        bgcm2d(i,j,jpcfc12  )=rmasko           
#endif              
        ENDIF
      ENDDO
      ENDDO

!js masking of coex90 etc
! attention: kbo is zero everywhere if called from end_bgc.f90 !!!!
!      if(p_pe==p_io) THEN
!       write(io_stdo_bgc,*) 'in write_bgcmean at month ' 
!       write(io_stdo_bgc,*) 'kbo= ', kbo(45,45),kbo(100,100)
!       write(io_stdo_bgc,*) 'n90depth= ',n90depth
!      ENDIF

      do j=1,kpje
       do i=1,kpie

       if (kbo (i,j).lt.n90depth) then
          bgcm2d(i,j,jcoex90  )=rmasko       
          bgcm2d(i,j,jopex90  )=rmasko       
          bgcm2d(i,j,jcaex90  )=rmasko       
       endif
       if (kbo (i,j).lt.n1000depth) then
          bgcm2d(i,j,jcoex1000)=rmasko       
          bgcm2d(i,j,jopex1000)=rmasko       
          bgcm2d(i,j,jcaex1000)=rmasko       
       endif
       if (kbo (i,j).lt.n2000depth) then
          bgcm2d(i,j,jcoex2000)=rmasko       
          bgcm2d(i,j,jopex2000)=rmasko       
          bgcm2d(i,j,jcaex2000)=rmasko      
       endif

       enddo
      enddo
      
!
! Set fill mode


!
!-----------------------------------------------------------------------
      IF (p_pe==p_io) THEN

      ncstat = NF_SET_FILL(nc_2d_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_2D: Problem with netCDF97')

!
! Write bgcmean data : 1-D bgc mean data
!
!-----------------------------------------------------------------------
!      ncstat = NF_INQ_VARID(nc_2d_id,'steps_p_m',ncvarid )
!      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_2D: Problem with netCDF102d')
!      ncstat = NF_PUT_VARA_DOUBLE (nc_2d_id,ncvarid,start_1d,count_1d,stepspm(meantime_2d) )
!      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_2D: Problem with netCDF103d')


      if(mean_2D_freq.ne.4) then
         t_hour=0.0
      else
         NDTDAY=NINT(86400./DT) 
         t_hour=real(ldtdayc)/real(ndtday)
      endif
      time2d=LDAYS+LMONTS*100+LYEARS*10000+t_hour

      ncstat = NF_INQ_VARID(nc_2d_id,'time',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_2D: Problem with netCDF102d')

!      write(0,*) 'in bgcmean 2d',nc_2d_id,ncvarid,start_1d,count_1d,time2d

      ncstat = NF_PUT_VARA_DOUBLE (nc_2d_id,ncvarid,start_1d,count_1d,time2d)


!      write(0,*) 'in bgcmean 2d',ncstat,NF_STRERROR(ncstat)


      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_2D: Problem with netCDF103d')

   END IF


!
! Write bgcmean data : 2-D bgc mean data
!
!-----------------------------------------------------------------------
                
!      ncstat = NF_INQ_VARID(nc_2d_id,'co2fluxdown_mean',ncvarid )
!      IF ( ncstat .NE. NF_NOERR ) STOP 'WRITE_BGCMEAN_2D: Problem with netCDF204'
!      ncstat = NF_PUT_VARA_DOUBLE                                         &
!     &         (nc_2d_id,ncvarid,start_2d,count_2d,bgcm2d(1,1,jco2fxd) )
!      IF ( ncstat .NE. NF_NOERR ) STOP 'WRITE_BGCMEAN_2D: Problem with netCDF205'

      nk=1
      CALL write_netcdf_var(nc_2d_id,'co2fluxdown_mean',bgcm2d(1,1,jco2fxd),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'co2fluxup_mean',bgcm2d(1,1,jco2fxu),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'pco2',bgcm2d(1,1,jpco2),nk,meantime_2d)

      CALL write_netcdf_var(nc_2d_id,'pco2',bgcm2d(1,1,jpco2),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'dms',bgcm2d(1,1,jdms),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'dmsflux',bgcm2d(1,1,jdmsflux),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'dmsprod',bgcm2d(1,1,jdmsprod),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'dms_bac',bgcm2d(1,1,jdms_bac),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'dms_uv',bgcm2d(1,1,jdms_uv),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'oxflux',bgcm2d(1,1,joxflux),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'kwco2',bgcm2d(1,1,jkwco2),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'niflux',bgcm2d(1,1,jniflux),nk,meantime_2d)

#ifdef DIFFAT      
      CALL write_netcdf_var(nc_2d_id,'atmco2',bgcm2d(1,1,jatmco2),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'atmo2',bgcm2d(1,1,jatmo2),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'atmn2',bgcm2d(1,1,jatmn2),nk,meantime_2d)
#endif

      CALL write_netcdf_var(nc_2d_id,'opex90',bgcm2d(1,1,jopex90),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'opex1000',bgcm2d(1,1,jopex1000),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'opex2000',bgcm2d(1,1,jopex2000),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'caex90',bgcm2d(1,1,jcaex90),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'caex1000',bgcm2d(1,1,jcaex1000),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'caex2000',bgcm2d(1,1,jcaex2000),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'coex90',bgcm2d(1,1,jcoex90),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'coex1000',bgcm2d(1,1,jcoex1000),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'coex2000',bgcm2d(1,1,jcoex2000),nk,meantime_2d)

      CALL write_netcdf_var(nc_2d_id,'export',bgcm2d(1,1,jexport),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'exposi',bgcm2d(1,1,jexposi),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'expoca',bgcm2d(1,1,jexpoca),nk,meantime_2d)

#ifdef ANTC14
      CALL write_netcdf_var(nc_2d_id,'ac14fx',bgcm2d(1,1,jac14fx),nk,meantime_2d)
#endif       
#ifdef PCFC
      CALL write_netcdf_var(nc_2d_id,'cfc11fx',bgcm2d(1,1,jcfc11fx),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'cfc12fx',bgcm2d(1,1,jcfc12fx),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'pcfc11',bgcm2d(1,1,jpcfc11),nk,meantime_2d)
      CALL write_netcdf_var(nc_2d_id,'pcfc12',bgcm2d(1,1,jpcfc12),nk,meantime_2d)
#endif      

!
! Reset mean fields
!
!-----------------------------------------------------------------------
      DO k=1,nbgcm2d
      DO j=1,kpje
      DO i=1,kpie
          bgcm2d(i,j,k) = 0.      
      ENDDO
      ENDDO
      ENDDO
    

      RETURN
      END








