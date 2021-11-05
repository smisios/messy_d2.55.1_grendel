      SUBROUTINE ATMOTR(kpie,kpje,kpke,kplmon,pdlxp,pdlyp)
!*******************************************************************
!
!**** *ATMOTR* - calculate the atmospheric diffusion
!
!     Maier-Reimer,          *MPI-Met, HH*
!
!     P.Wetzel,              *MPI-Met, HH*    27.03.02
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     - calculate the atmospheric diffusion for surface CO2 exchange.
!
!     Method
!     -------
!     -
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *CALL*       *ATMOTR(kpie,kpje,kpke)*
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_control_bgc
      USE mo_bgcmean
      USE mo_param1_bgc

      USE mo_bgc_diagnostic                                  !CDI tinka

      USE mo_boundsexch, ONLY : bounds_exch
      IMPLICIT NONE

      INTEGER :: kpie,kpje,kpke,kplmon
      INTEGER :: i,j,l,is
      REAL(wp), DIMENSION(kpie,kpje) :: ucor, vcor
      REAL(wp), DIMENSION(kpie,kpje) :: pdlxp,pdlyp
#if !defined(__cpl_co2) && defined(PANTHROPOCO2)
      REAL(wp) :: ppm2con
#endif

      do j=1,kpje
      do i=1,kpie
        ucor(i, j) = 0._wp
        vcor(i, j) = 0._wp
      enddo
      enddo

!      do l=1,natm ! outer loop over all tracers

!     call perio2(123,atm(1,1,l))
      call bounds_exch(1,'p',atm)

      do l=1,3     ! 3 hard coded !!!!!!!!!!!!!!!!!!!!!!!!!
      do is=1,4 ! speed up

      do j=2,kpje-1
      do i=2,kpie-1
!        ucor(i,j)=0.01*(atm(i+1,j)-atm(i,j))*dlxu(i,j)*dlyu(i,j)
!        vcor(i,j)=0.01*(atm(i,j)-atm(i,j+1))*dlxv(i,j)*dlyv(i,j)*atdifv(i,j)
        ucor(i, j) = 0.2_wp * (atm(i+1, j, l) - atm(i, j, l)) &
             * pdlxp(i, j) * pdlyp(i, j)
        vcor(i, j) = 0.2_wp * (atm(i, j, l) - atm(i, j+1, l)) &
             * pdlxp(i, j) * pdlyp(i, j) * atdifv(i, j)
      enddo
      enddo
!     call perio2(124,ucor)
!     call perio2(125,vcor)
      call bounds_exch(1,'u+',ucor)
      call bounds_exch(1,'v+',vcor)

      do j=2,kpje
      do i=2,kpie
        atm(i,j,l)=atm(i,j,l)+(ucor(i,j)-ucor(i-1,j)                &
     &           +vcor(i,j-1)-vcor(i,j))/(pdlxp(i,j)*pdlyp(i,j))
      enddo
      enddo
!     call perio2(126,atm(1,1,l))

      enddo  ! speed up
      enddo  ! outer loop over all tracers
      call bounds_exch(1,'p',atm)

      do j=1,kpje
      do i=1,kpie

        bgcm2d(i,j,jatmco2)=bgcm2d(i,j,jatmco2)+atm(i,j,iatmco2)
        bgcm2d(i,j,jatmo2) =bgcm2d(i,j,jatmo2) +atm(i,j,iatmo2)
        bgcm2d(i,j,jatmn2) =bgcm2d(i,j,jatmn2) +atm(i,j,iatmn2)
      enddo
      enddo

      RETURN
      END
