      SUBROUTINE WRITE_BGCMEAN_BIOZ(kpie,kpje,kpke,pddpo)

!****************************************************************
!
!**** *WRITE_BGCMEAN_BIOZ* - calculate and write 2-dimensional bgc mean data.
!
!     Patrick Wetzel,    *MPI-Met, HH*    24.05.04
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Write and mask bgc mean data.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *WRITE_BGCMEAN_BIOZ(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
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
!     NOTE:
!
!     output is written with frequency of 2D-fields, even though fields are 3D
!
!**************************************************************************

      USE mo_bgcmean
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      use mo_param1_bgc 

      USE mo_control_bgc

      USE mo_parallel
 
      use mo_commo1, only: lyears,lmonts,ldays,ldtdayc,dt

      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,l,nk,ndtday
      REAL :: pddpo(kpie,kpje,kpke)

      real :: timebioz,t_hour

      INCLUDE 'netcdf.inc'
      INTEGER :: ncvarid,ncstat,ncoldmod


      
      INTEGER :: START_1D,COUNT_1D       ! size of internal array 
      INTEGER :: START_3D(4),COUNT_3D(4) ! size of internal array 

      
!-----------------------------------------------------------------------

      START_1D =meantime_2d
      COUNT_1D =1

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

      DO l=1,nbgcm3d
      DO k=1,kwrbioz
      DO j=1,kpje
      DO i=1,kpie
        IF(pddpo(i,j,k) .LT. 0.5) THEN
          bgcm3d(i,j,k,l)=rmasko
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
!
! Set fill mode
!
!-----------------------------------------------------------------------
!

      if (p_pe==p_io) then
      ncstat = NF_SET_FILL(nc_bioz_id,NF_NOFILL, ncoldmod)
      IF ( ncstat .NE. NF_NOERR ) call stop_all('WRITE_BGCMEAN_BIOZ: Problem with netCDF97')
     
!
! Write bgcmean data : 1-D bgc mean data
!
!-----------------------------------------------------------------------
!      ncstat = NF_INQ_VARID(nc_bioz_id,'steps_p_m',ncvarid )
!      IF ( ncstat .NE. NF_NOERR ) call stop_all('WRITE_BGCMEAN_BIOZ: Problem with netCDF102d')
!      ncstat = NF_PUT_VARA_DOUBLE (nc_bioz_id,ncvarid,start_1d,count_1d,stepspm(meantime_2d) )
!      IF ( ncstat .NE. NF_NOERR ) call stop_all('WRITE_BGCMEAN_BIOZ: Problem with netCDF103d')


      if(mean_2D_freq.ne.4) then
         t_hour=0.0
      else
         NDTDAY=NINT(86400./DT) 
         t_hour=real(ldtdayc)/real(ndtday)
      endif

      timebioz=LDAYS+LMONTS*100+LYEARS*10000+t_hour


      ncstat = NF_INQ_VARID(nc_bioz_id,'time',ncvarid )
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_BIOZ: Problem with netCDF102d')

!       write(0,*) 'in bgcmean bioz',nc_bioz_id,ncvarid,start_1d,count_1d,timebioz

      ncstat = NF_PUT_VARA_DOUBLE (nc_bioz_id,ncvarid,start_1d,count_1d,timebioz)
      IF ( ncstat .NE. NF_NOERR ) call STOP_ALL('WRITE_BGCMEAN_BIOZ: Problem with netCDF103d')

      end if


!
! Write bgcmean data : 3-D euphotic layer (here k=1:8) 
!
!-----------------------------------------------------------------------
      nk = kwrbioz
      CALL write_netcdf_var(nc_bioz_id,'phosph',bgcm3d(1,1,1,jphosph),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'phyto',bgcm3d(1,1,1,jphyto),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'grazer',bgcm3d(1,1,1,jgrazer),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'phosy',bgcm3d(1,1,1,jphosy),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'oxygen',bgcm3d(1,1,1,joxygen),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'iron',bgcm3d(1,1,1,jiron),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'nitrate',bgcm3d(1,1,1,jano3),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'alkali',bgcm3d(1,1,1,jalkali),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'silica',bgcm3d(1,1,1,jsilica),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'dic',bgcm3d(1,1,1,jdic),nk,meantime_2d)
      CALL write_netcdf_var(nc_bioz_id,'doc',bgcm3d(1,1,1,jdoc),nk,meantime_2d)

!     CALL write_netcdf_var(nc_bioz_id,'dic13',bgcm3d(1,1,1,isco213),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'dic14',bgcm3d(1,1,1,isco214),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'poc',bgcm3d(1,1,1,jpoc),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'atten',bgcm3d(1,1,1,jatten),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'dms',bgcm3d(1,1,1,jdms),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'calc',bgcm3d(1,1,1,jcalc),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'opal',bgcm3d(1,1,1,jopal),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'export',bgcm3d(1,1,1,jexport),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'expoca',bgcm3d(1,1,1,jexpoca),nk,meantime_2d)
!     CALL write_netcdf_var(nc_bioz_id,'exposi',bgcm3d(1,1,1,jexposi),nk,meantime_2d)

 
#ifdef AGG
      CALL write_netcdf_var(nc_bioz_id,'numbers',bgcm3d(1,1,1,jnos),nk,meantime_2d)
#endif


!
! Reset mean fields
!
!-----------------------------------------------------------------------
!

      DO l=1,nbgcm3d
      DO k=1,kwrbioz
      DO j=1,kpje
      DO i=1,kpie
        bgcm3d(i,j,k,l) = 0.
      ENDDO
      ENDDO
      ENDDO
      ENDDO


      RETURN
      END
