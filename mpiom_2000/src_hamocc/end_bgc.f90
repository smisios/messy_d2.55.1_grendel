      SUBROUTINE END_BGC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp,      &
     &        pgila,pgiph,ptiestu,kplyear,kplmon,kplday,kpldtoce)
!****************************************************************
!
!**** *END_BGC* - finish with marine bio-geo-chemistry module.
!
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     - call inventory
!     - save time series
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!     called by mpiom.f90
!
!     *CALL*       *END_BGC(list....)*
!
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
!
!     Externals
!     ---------
!     INVENTORY_BGC SAVE_TIMESER_BGC
!
!**********************************************************************
      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_bgcmean
      USE mo_control_bgc
      use mo_param1_bgc

      use mo_parallel
implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k
      INTEGER :: kplyear,kplmon,kplday,kpldtoce
      REAL(wp) :: pddpo(kpie,kpje,kpke)
      REAL(wp) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL(wp) :: pgila(kpie*2,kpje*2)
      REAL(wp) :: pgiph(kpie*2,kpje*2)
      REAL(wp) :: ptiestu(kpke)

!
! BGC mean for last timestep
!
!!$#ifdef PDYNAMIC_BGC
!!$      call avrg_dynamic(kpie,kpje,kpke) ! has to be called before AVRG_BGCMEAN !
!!$      CALL WRITE_DYNAMIC(kpie,kpje,kpke,pddpo)
!!$      CALL close_dynamic
!!$#endif
!!$      if(p_pe==p_io) THEN
!!$       write(io_stdo_bgc,*) 'before avrg_bgcmean in end_bgc'
!!$!       write(io_stdo_bgc,*) 'kbo= ', kbo(45,45),kbo(100,100)
!!$       write(io_stdo_bgc,*) 'n90depth= ',n90depth
!!$      ENDIF
!!$
!!$      call avrg_bgcmean_2d(kpie,kpje,kpke)
!!$      call avrg_bgcmean_3d(kpie,kpje,kpke)
      k=kpke
      DO 1321 j=1,kpje
      DO 1321 i=1,kpie
         kbo(i,j)=1
         bolay(i, j) = 0._wp
         IF (pddpo(i, j, k) .GT. 0.5_wp) THEN
            bolay(i,j)=pddpo(i,j,k)
            kbo(i,j)=k
         ENDIF
1321  CONTINUE
! evaluate min depth of last layer
      bolaymin = 8000._wp
      DO 1322 k=kpke-1,1,-1
      DO 1322 j=1,kpje
      DO 1322 i=1,kpie
         IF (pddpo(i, j, k) .GT. 0.5_wp .AND. pddpo(i, j, k+1) .LE. 0.5_wp) THEN
            bolay(i,j)=pddpo(i,j,k)
            kbo(i,j)=k
            bolaymin = min(bolaymin,bolay(i,j))
         ENDIF
1322  CONTINUE

!!$         if(p_pe==p_io) THEN
!!$!       write(io_stdo_bgc,*) 'kbo= ', kbo(45,45),kbo(100,100)
!!$       write(io_stdo_bgc,*) 'n90depth= ',n90depth
!!$      ENDIF
!!$
!!$      CALL WRITE_BGCMEAN_2D(kpie,kpje,kpke,pddpo)
!!$      CALL WRITE_BGCMEAN_3D(kpie,kpje,kpke,pddpo)
!!$      CALL WRITE_BGCMEAN_BIOZ(kpie,kpje,kpke,pddpo)
!!$      CALL WRITE_BGCMEAN_SED(kpie,kpje,kpke,pddpo)


#ifdef OLD_IO
!       write(0,*) 'calling close_bgcmean_bioz from end_bgc'
      CALL CLOSE_BGCMEAN_BIOZ
!       write(0,*) 'calling close_bgcmean_2d from end_bgc'
      CALL CLOSE_BGCMEAN_2D
!       write(0,*) 'calling close_bgcmean_3d from end_bgc'
      CALL CLOSE_BGCMEAN_3D
!       write(0,*) 'calling close_bgcmean_sed from end_bgc'
      CALL CLOSE_BGCMEAN_SED
#endif/*def OLD_IO */


!
! Global inventory of all tracers
!
      CALL INVENTORY_BGC(kpie,kpje,kpke)

!
! save the time series of bgc
!
      CALL SAVE_TIMESER_BGC

#ifdef DMSASSIM
      CALL DMS_ASSIM(kpie,kpje)
#endif

      RETURN
      END
