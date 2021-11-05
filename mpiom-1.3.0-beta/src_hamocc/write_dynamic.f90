      SUBROUTINE WRITE_DYNAMIC(kpie,kpje,kpke,pddpo)

!****************************************************************
!
!**** *WRITE_DYNAMIC* - calculate and write 2-dimensional bgc mean data.
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
!     *CALL*       *WRITE_DYNAMIC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,
!                  ptiestu) *
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL :: of model grid.
!     *INTEGER* *kpje*    - 2nd REAL :: of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL :: of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd REAL ::) [m].
!     *REAL*    *pgila*   - geographical longitude of grid points [degree E].
!     *REAL*    *pgiph*   - geographical latitude  of grid points [degree N].
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

      USE mo_bgcmean
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      use mo_param1_bgc 
      USE mo_dynamic

      USE mo_control_bgc

      USE mo_parallel
 
      implicit none
      INTEGER :: kpie,kpje,kpke,i,j,k,l,nk
      REAL :: pddpo(kpie,kpje,kpke)


      INCLUDE 'netcdf.inc'
      INTEGER :: ncvarid,ncstat


      INTEGER :: START_2D(3),COUNT_2D(3) ! size of internal array 
      INTEGER :: START_3D(4),COUNT_3D(4) ! size of internal array 

      
!-----------------------------------------------------------------------

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
      COUNT_3D(3) =nbgcdyn
      COUNT_3D(4) =1      
      
!-----------------------------------------------------------------------
!
!  Masking bgcmean data.
!
!-----------------------------------------------------------------------
      
      DO l=1,kdtot
      DO k=1,nbgcdyn
      DO j=1,kpje
      DO i=1,kpie
        IF(pddpo(i,j,1) .LT. 0.5) THEN
        bgcdyn(i,j,k,l)=rmasko
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      DO j=1,kpje
      DO i=1,kpie
        IF(pddpo(i,j,1) .LT. 0.5) THEN
        bgc_nmld(i,j)=rmasko
        bgc_zmld(i,j)=rmasko
        ENDIF
      ENDDO
      ENDDO    

!
! Write WRITE_DYNAMIC data : 2-D mean data
!
!-----------------------------------------------------------------------
      nk=1
      CALL write_netcdf_var(nc_dyn_id,'nmld',bgc_nmld(1,1),nk,meantime_2d)
      CALL write_netcdf_var(nc_dyn_id,'zmld',bgc_zmld(1,1),nk,meantime_2d)

! Write WRITE_DYNAMIC data : 2-D mean data
!
!-----------------------------------------------------------------------

      nk=nbgcdyn
      CALL write_netcdf_var(nc_dyn_id,'adv',bgcdyn(1,1,1,kdadv),nk,meantime_2d)
      CALL write_netcdf_var(nc_dyn_id,'dif',bgcdyn(1,1,1,kddif),nk,meantime_2d)
      CALL write_netcdf_var(nc_dyn_id,'pre',bgcdyn(1,1,1,kdpre),nk,meantime_2d)
      CALL write_netcdf_var(nc_dyn_id,'gmp',bgcdyn(1,1,1,kdgmp),nk,meantime_2d)
      CALL write_netcdf_var(nc_dyn_id,'bio',bgcdyn(1,1,1,kdbio),nk,meantime_2d)
!
! Reset mean fields
!
!-----------------------------------------------------------------------
      DO l=1,kdtot
      DO k=1,nbgcdyn
      DO j=1,kpje
      DO i=1,kpie
        bgcdyn(i,j,k,l)=0.
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      
      DO j=1,kpje
      DO i=1,kpie
        bgc_nmld(i,j)=0.
        bgc_zmld(i,j)=0.
      ENDDO
      ENDDO    
    

      RETURN
      END
