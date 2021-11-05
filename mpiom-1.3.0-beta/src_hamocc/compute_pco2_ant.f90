      SUBROUTINE COMPUTE_PCO2_ANT                                               &
     &           (kpie,kpje,kpke,pddpo,psao,                          &
     &            pco2,kplmon,kplday,kmonlen)
!**********************************************************************
!
!**** *COMPUTE_PCO2_ANT* - .
!
!     Patrick Wetzel,    *MPI-Met, HH*    20.07.04
!
!     Purpose
!     -------
!     Inorganic carbon cycle.
!
!     Method
!     -------
!
!     *CALL* *COMPUTE_PCO2_ANT(kpie,kpje,kpke,pddpo,
!                          kplmon,kplday,kmonlen)*
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL of model grid.
!     *INTEGER* *kpje*    - 2nd REAL of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd REAL) [m].
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *pco2*    - sea level pressure [Pascal].
!
!     Externals
!     ---------
!     .
!
!**********************************************************************
#ifdef PANTHROPOCO2

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_timeser_bgc
      USE mo_control_bgc
      USE mo_bgcmean
      use mo_param1_bgc 

      implicit none
      INTEGER :: i,j,k,kpie,kpje,kpke
      INTEGER :: kplmon,kplday,kmonlen,laumo1

      REAL psao(kpie,kpje,kpke)
      REAL pddpo(kpie,kpje,kpke)
      REAL pco2(kpie,kpje)
      
      REAL :: dddhhh,dadh,a,h,c,alk,t1,t2
      REAL :: akbi,ak2i,ak1i
      REAL :: contppm
      REAL :: AHI, ABE,RMONLEN,RPLDAY
      REAL :: AK0,AK1,AK2,AKB,AKW,BT,oxysa,anisa


      laumo1=kplmon+1
      IF(laumo1.GT.12) laumo1=1
      
      RMONLEN=kmonlen
      RPLDAY=kplday
      AHI=RPLDAY/RMONLEN
      ABE=1.-AHI

      k=1

      DO 1 j=1,kpje
      DO 1 i=1,kpie

      IF(pddpo(i,j,1).GT.0.5) THEN

!*       21.11 SET CHEMICAL CHEMICAL CONSTANTS
!      AK0  =AHI*CHEMCM(i,j,5,LAUMO1)+ABE*CHEMCM(i,j,5,kplmon)
!      AK1  =AHI*CHEMCM(i,j,4,LAUMO1)+ABE*CHEMCM(i,j,4,kplmon)
!      AK2  =AHI*CHEMCM(i,j,3,LAUMO1)+ABE*CHEMCM(i,j,3,kplmon)
!      AKB  =AHI*CHEMCM(i,j,1,LAUMO1)+ABE*CHEMCM(i,j,1,kplmon)
!      AKW  =AHI*CHEMCM(i,j,2,LAUMO1)+ABE*CHEMCM(i,j,2,kplmon)
!      BT   =AHI*CHEMCM(i,j,6,LAUMO1)+ABE*CHEMCM(i,j,6,kplmon)
!      oxysa=AHI*CHEMCM(i,j,7,LAUMO1)+ABE*CHEMCM(i,j,7,kplmon)
!      anisa=AHI*CHEMCM(i,j,8,LAUMO1)+ABE*CHEMCM(i,j,8,kplmon)
! CHEMCM is now computed daily, so no time interpolation !!!!!
      AK0  =CHEMCM(i,j,5,kplmon)
      AK1  =CHEMCM(i,j,4,kplmon)
      AK2  =CHEMCM(i,j,3,kplmon)
      AKB  =CHEMCM(i,j,1,kplmon)
      AKW  =CHEMCM(i,j,2,kplmon)
      BT   =CHEMCM(i,j,6,kplmon)
      oxysa=CHEMCM(i,j,7,kplmon)
      anisa=CHEMCM(i,j,8,kplmon)

      ak1i=1./ak1
      ak2i=1./ak2
      akbi=1./akb
      
      contppm=1./0.35e-3

! Calculate new hi concentration

         h=hi_ant(i,j,k)
         c=ocetra(i,j,k,isco2_ant)
         akw=akw3(i,j,k)
         bt=rrrcl*psao(i,j,k)
         akb=akb3(i,j,k)
         alk=ocetra(i,j,k,ialk_ant)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2

         dddhhh=a/dadh
         hi_ant(i,j,k)=hi_ant(i,j,k)-dddhhh
         h=hi_ant(i,j,k)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2 
         dddhhh=a/dadh
         hi_ant(i,j,k)=hi_ant(i,j,k)-dddhhh
         h=hi_ant(i,j,k)

         t1=h*ak1i
         t2=h*ak2i
         a=c*(2.+t2)/(1.+t2+t2*t1)+akw/h-h+bt/(1.+h*akbi)-alk
         dadh=c*(1./(ak2*(1.+t2+t2*t1))-(2.+t2)*((1.+2.*t1)*ak2i)/  &
     &       (1.+t2+t2*t1)**2)                                      &
     &       -akw/h**2-1.-(bt*akbi)/(1.+h*akbi)**2
         dddhhh=a/dadh
         hi_ant(i,j,k)=hi_ant(i,j,k)-dddhhh
         h=hi_ant(i,j,k)

!
! Calculate pCO2 [ppmv] from total dissolved inorganic carbon (DIC: SCO212)
! the calculation also includes solubility
!
         pco2(i,j)=c/((1.+ak1*(1.+ak2/h)/h)*ak0)


      ENDIF
1     CONTINUE

#endif /*PANTHROPOCO2*/
      RETURN
      END
