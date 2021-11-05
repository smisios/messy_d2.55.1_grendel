      SUBROUTINE DILUTE_BGC(kpie,kpje,l1_bgc,l1_new)
!*******************************************************************
!
!**** *DILUTE_BGC* - dilute tracer concentration.
!
!     S.Legutke,        *MPI-MaD, HH*    08.08.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    04.03.02
!     - only advected tracer may be diluted.
!     P.Wetzel ,        *MPI-Met, HH*    18.10.02
!     - dilution by change of ZO  ??
!
!     Purpose
!     -------
!     - new tracer concentration due to freshwater flux.
!     js: in the call to dilute distinction between l1_bgc and l1_new
!         is only from vertical velocity component * time step ????:
!   excerpt from mpiom.f90:
!           layer1_new(i,j)=DDPO(I,J,1)+ZO(I,J)                    &
!    &       -WO(I,J,1)*DT-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA
!           layer1_bgc(i,j)=DDPO(I,J,1)+ZO(I,J)                    &
!                         -SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA
!
!
!     Method
!     -------
!     Total freshwater flux is used to calculate new
!     tracer concentration. The update accounts for all fluxes of
!     water relevant for tracer concentration changes (e.g.
!     P-E, runoff, snowmelt (pfresh), and icemelt/growth (pbrine)).
!     Tracer concentration in all water fluxes is assumed to be zero.
!
!
!**   Interface.
!     ----------
!
!     *CALL*       *DILUTE_BGC(kpie,kpje,l1_bgc,l1_new)*
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *REAL*    *pfresh*  - 1st layer height change due to atm. fluxes [m/s].
!     *REAL*    *pbrine*  - 1st layer height change due to ice growth [m/s].
!     *REAL*    *pzold*   - old height of upper layer [m].
!     *REAL*    *pddpo*   - land/sea mask.
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      use mo_param1_bgc

implicit none

      INTEGER :: kpie,kpje,i,j,l
      REAL(wp) :: l1_bgc(kpie,kpje)
      REAL(wp) :: l1_new(kpie,kpje)

      REAL(wp) zfac(kpie,kpje)

      DO j=1,kpje
      DO i=1,kpie
         zfac(i, j) = 1.0_wp
         if (l1_new(i, j) .GT. 0.0_wp) then
            zfac(i,j)=l1_bgc(i,j)/l1_new(i,j)
         endif
      ENDDO
      ENDDO

      DO l=1,nocetra

      DO j=1,kpje
      DO i=1,kpie
        ocetra(i,j,1,l) = ocetra(i,j,1,l)*zfac(i,j)
      ENDDO
      ENDDO

      ENDDO

      RETURN
      END
