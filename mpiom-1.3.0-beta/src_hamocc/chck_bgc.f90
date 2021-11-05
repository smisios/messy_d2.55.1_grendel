      SUBROUTINE CHCK_BGC(kunit,kpicycli,ystring,kpie,kpje,kpke,pddpo)

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/chck_bgc.f90,v $\\
!$Revision: 1.2.2.1.4.1.2.2.4.1.2.2.2.3.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!$Name: mpiom_1_2_0 $\\
!
!*********************************************************************
!
!**** *CHCK_BGC* - check marine bio-geo-chemistry fields
!
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Modified
!     --------
!     js: this routine is very time consuming, needs to
!         be thinned out (e.g. do k=1,kpke,3) or so
!     
!     Purpose
!     -------
!     - check max/min values on wet/dry cells.
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *CHCK_BGC(kunit,kpicycli,ystring,kpie,kpje,kpke,pddpo)*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!
!     called by:
!
!               aufw_bgc
!               carchm
!               ini_bgc
!               ocprod
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kunit*    - stdout logical unit.
!     *INTEGER* *kpicycli* - flag for cyclicity.
!     *CHARACTER* *ystring*- .
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!
!     Externals
!     ---------
!     EXTR
!
!***********************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      use mo_param1_bgc 

      USE mo_control_bgc

      USE mo_param1, ONLY: ie_g

implicit none

!mz_ap_20070626+
      EXTERNAL extr
!mz_ap_20070626-

      CHARACTER*(*) ystring
      INTEGER :: kunit,kpicycli,kpie,kpje,kpke,i,j,k,l
      REAL ::  pddpo(kpie,kpje,kpke)

      WRITE(kunit,*) ystring

!                        
! check for invalid values in dry and wet cells:
!
      WRITE(kunit,*)'Check values of ocean tracer hi :'
      DO k=1,kpke,5
         WRITE(kunit,*)'      layer ',k,' :'
         CALL EXTR(kpie,kpje,hi(1,1,k),pddpo(1,1,k),rmasko,kunit)
      ENDDO
      WRITE(kunit,*)'Check values of ocean tracer CO3 :'
      DO k=1,kpke,5
         WRITE(kunit,*)'      layer ',k,' :'
         CALL EXTR(kpie,kpje,co3(1,1,k),pddpo(1,1,k),rmasko,kunit)
      ENDDO

                         
! check for invalid values of oceanic tracers in dry and wet cells
!
      DO l=1,nocetra
         WRITE(kunit,*)                                                &
     &   'Check values of ocean tracer no. ',l,' :'
      DO k=1,kpke,5
         WRITE(kunit,*)'      layer ',k,' :'
         CALL EXTR(kpie,kpje,ocetra(1,1,k,l),pddpo(1,1,k),rmasko,kunit)
      ENDDO
      ENDDO

!                        
! check for invalid values in dry and wet cells of sediment tracer.
!
      DO l=1,npowtra
         WRITE(kunit,*)                                                &
     &   'Check values of pore water sediment tracer no. ',l,' :'
      DO k=1,ks,3
         WRITE(kunit,*)'      layer ',k,' :'
         CALL EXTR(kpie,kpje,powtra(1,1,k,l),bolay,rmasks,kunit)
      ENDDO
      ENDDO
!     write(0,291)issster
! 291   format('issster=',i5)
      DO l=1,nsedtra
         WRITE(kunit,*)                                                &
     &   'Check values of solid sediment tracer no. ',l,' :'
      DO k=1,ks,3
         WRITE(kunit,*)'      layer ',k,' :'
         CALL EXTR(kpie,kpje,sedlay(1,1,k,l),bolay,rmasks,kunit)
      ENDDO
      ENDDO

      WRITE(kunit,*)                                                   &
     &'Check values of accumulated sediment tracer sedhpl (hi):'
      DO k=1,ks,3
         WRITE(kunit,*)'      layer ',k,' :'
         CALL EXTR(kpie,kpje,sedhpl(1,1,k),bolay,rmasks,kunit)
      ENDDO

!                        
! check for cyclicity: sediment tracer.
! RJ: To do: Implement the following check also for the distributed case
!
      IF ( kpicycli .EQ. 1 .AND. kpie .EQ. ie_g) THEN
      WRITE(kunit,*)'doing cyclicity check...'
      DO l=1,npowtra
      DO k=1,ks,3
      DO j=1,kpje
      IF (ABS(powtra(1,j,k,l)-powtra(kpie-1,j,k,l)).GT.1.e-15 .OR.     &
     &    ABS(powtra(2,j,k,l)-powtra(kpie  ,j,k,l)).GT.1.e-15     ) THEN
         WRITE(kunit,*)                                                &
     &   'Sediment tracer no. ',l,' is not cyclic at j,k=',j,k,' : ',  &
     &    powtra(1,j,k,l),powtra(kpie-1,j,k,l),                        &
     &    powtra(2,j,k,l),powtra(kpie  ,j,k,l)
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!                        
! check for cyclicity: ocean tracer.
!
      DO l=1,nocetra
      DO k=1,kpke
      DO j=1,kpje
      IF (ABS(ocetra(1,j,k,l)-ocetra(kpie-1,j,k,l)).GT.1.e-15 .OR.     &
     &    ABS(ocetra(2,j,k,l)-ocetra(kpie  ,j,k,l)).GT.1.e-15     ) THEN
         WRITE(kunit,*)                                                &
     &   'Ocean tracer no. ',l,' is not cyclic at j,k=',j,k,' : ',     &
     &    ocetra(1,j,k,l),ocetra(kpie-1,j,k,l),                        &
     &    ocetra(2,j,k,l),ocetra(kpie  ,j,k,l)
      ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDIF


      RETURN
      END
