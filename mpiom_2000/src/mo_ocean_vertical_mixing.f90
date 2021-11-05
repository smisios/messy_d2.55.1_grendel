!> Calculation of vertical eddy viscosity and vertical eddy diffusivity.
!> @author Ernst Maier-Reimer, Uwe Mikolajewicz, ... 1999
!> @date Last modified on 09.03.2009 by Helmuth Haak
!> @todo Documentation , add interface to GOTM
!> This module holds subroutines for calculation of vertical eddy viscosity and vertical eddy diffusivity .
!! Included subroutines are:
!! calc_rinum    : Calculation of vertical eddy viscosity avo and vertical eddy diffusivity dvo
!! following the Philander and Pacanowski scheme

MODULE mo_ocean_vertical_mixing
  USE mo_kind, ONLY: wp
  USE mo_param1, ONLY: kep
  USE mo_constants, ONLY: api
  USE mo_units, ONLY: io_stdout
  USE mo_commo1, ONLY: tiestw
  IMPLICIT NONE
  PRIVATE
  REAL(wp) :: av0, dv0, wt, wa
  REAL(wp), ALLOCATABLE :: dbackv(:), abackv(:)
  PUBLIC :: av0, dv0
  PUBLIC :: setup_ocean_vertical_mixing, calc_rinum
CONTAINS

  SUBROUTINE setup_ocean_vertical_mixing(cwt, cwa, aback, dback)
    REAL(wp), INTENT(in) :: cwt, cwa, aback, dback

#ifdef DBACKPROFIL
    INTEGER :: idback
#endif
    INTEGER :: k
    REAL(wp) :: gfdl_diff

    wt = cwt / 6._wp**3
    wa = cwa / 6._wp**3
    ALLOCATE(dbackv(kep), abackv(kep))
    !HH   background diffusion
    abackv(:) = aback
    dbackv(:) = dback
#ifdef DBACKPROFIL
    IDBACK = 0
    DO K=1,KEP
      ABACKV(K) = ABACK
      dbackv(k) = dback                                                &
           + (1 - idback) * 1.e-4_wp                                   &
           & * (0.8_wp                                                 &
           &    + 1.05_wp/api                                          &
           &      * ATAN(4.5_wp * 1.e-3_wp * (tiestw(k) - 2500._wp)))  &
           & * (0.5_wp + SIGN(0.5_wp, tiestw(k) - 500._wp))            &
           & * SQRT(ABS((tiestw(k) - 500._wp) / (3500._wp - 500._wp)))
    END DO
#endif

    DO K=1,KEP
      gfdl_diff = 1.e-4_wp                                                     &
           * (0.8_wp + 1.05_wp                                                 &
           &           / api * ATAN(4.5_wp * 1.e-3_wp * (tiestw(k) - 2500._wp)))
#ifdef DBACKGFDL
      dbackv(k) = gfdl_diff
#endif
    END DO

#ifdef DBACKGFDL2
    ! FIXME: shouldn't this be inside the above loop?
    dbackv(k) = 0.5_wp * (dback + gfdl_diff)
#endif

    DO K=1,KEP
      WRITE(IO_STDOUT,6002)'BACKGROUND DIFFUSIVITY AT '              &
           ,INT(TIESTW(K))                           &
           ,'M : HOPE : ',DBACKV(K),' GFDL : ',GFDL_DIFF
      WRITE(IO_STDOUT,6002)'BACKGROUND VISCOSITY AT ',INT(TIESTW(K)) &
           ,'M : HOPE : ',ABACKV(K)
    END DO

6002 FORMAT(1X,A27,I5,A10,E12.3,A8,E12.3)

  END SUBROUTINE setup_ocean_vertical_mixing

   !> calculation of vertical eddy viscosity avo and vertical eddy diffusivity dvo
   !! following the Philander and Pacanowski scheme
   !!
   !! this requires the calculation of the Richardson number, hence the name of the routine
   !!
   !! in case of the convection parametrizations ioconv = 1 and ioconv = 4, the
   !! vertical eddy viscosity avo and vertical eddy diffusivity dvo are
   !! increased to the namelist parameters cavocon and cdvocon if they are
   !! smaller than those
   !!
   !! rewritten by Aiko Voigt, Max Planck Institute for Meteorology (Jan 13, 2009)

  SUBROUTINE calc_rinum

    USE mo_mpi, only : p_abort
    USE mo_mean
    USE mo_commo1
    USE mo_planetary_constants, ONLY: g,rhoref_water
    USE mo_boundsexch, ONLY : bounds_exch

    !! local scalars

    INTEGER i         !< loop variable in W-E direction
    INTEGER j         !< loop variable in N-S direction
    INTEGER k         !< loop variable for vertical levels (k=1 top, k=ke bottom)
    INTEGER ko
    INTEGER ku

    REAL(wp)    rinumo    !< Richardson number
    REAL(wp)    stabeps
    REAL(wp)    avohelp   !< vert. eddy viscosity
    REAL(wp)    dvohelp   !< veri. eddy diffusivity
    REAL(wp)    dudo
    REAL(wp)    gor
    REAL(wp)    hho
    REAL(wp)    topnah
    REAL(wp)    wpendep
    REAL(wp)    wtdecay
    REAL(wp)    relne
    REAL(wp)    relax
    REAL(wp)    dudz
    REAL(wp)    drdz0
    REAL(wp)    cra
    REAL(wp)    crd

    !! local array

    REAL(wp)    zsurf(ke) !< mask for top layer, i.e. zsurf=1 for k=1 and 0 otherwise

!$omp parallel private(i,j,k, rinumo,ku,ko,stabeps,wtdecay,dudo,hho    &
!$omp                 ,avohelp,dvohelp)

    !======================================================================
    !
    !     d)
    !
    !                   calculation of richardson number dependent
    !                      vertical eddy viscosity   av     and
    !                      vertical eddy diffusivity dv
    !
    !
    !         rinum : ( (g/rho)*d(rho)/dz ) / ( (d(velocity)/dz)**2 )
    !                  richardson number (even,odd ==> rinume,rinumo)
    !         gor   : g/rho
    !
    !         av0   : numerical value of vertical eddy viscosity in case
    !                 of neutral stability, i.e. free turbulence
    !
    !---------------------------------------------------------------------
    zsurf(1) = 1._wp
    DO k=2,ke
       zsurf(k) = 0._wp
    ENDDO

    !
    relne = 0.4_wp
    relax = 1._wp - relne
    !
    dudz = 1.e4_wp
    drdz0 = 1.e-3_wp
    !
    gor=g/rhoref_water
    !
    !--------------------------------------------------------------------
    !
    !     c.1)
    !
    !     vertical eddy viscosity  (for momentum equation)
    !
    !--------------------------------------------------------------------
    !
    !     d.1)    mixed-layer turbulence
    !
    !  amplitude of turbulence decays by factor wtdecay every model level.
    !  turbulence stops once density difference reaches the equivalent of
    !   wtdt
    !  temperature difference.
    !  turbulence under ice is / is not  enhanced.
    !                             ==
    cra = 5._wp
    crd = 5._wp

    wpendep = 40._wp
    !FIXME: duplicate from lines 97
    relne = 0.4_wp
    relax = 1._wp - relne
    wtdecay=EXP(-dzw(1)/wpendep)

!$omp do
    DO j=1,je
       DO i=1,ie
          !hh       t1o(i,j,1)=wt*2.*(1.-sicomo(i,j))*fu10(i,j)**3*wtdecay
          !hh       s1o(i,j,1)=wt*(1.-sicomo(i,j))*fu10(i,j)**3*wtdecay
#ifdef REDWMICE
          t1o(i, j, 1) = wa * (1._wp - sicomo(i, j))**2 * fu10(i, j)**3
          s1o(i, j, 1) = wt * (1._wp - sicomo(i, j))**2 * fu10(i, j)**3
#else
          t1o(i, j, 1) = wa * (1._wp - sicomo(i, j)) * fu10(i, j)**3
          s1o(i, j, 1) = wt * (1._wp - sicomo(i, j)) * fu10(i, j)**3
#endif /*REDWMICE*/

       ENDDO
    ENDDO
!$omp end do

    if (iMEAN.ne.0)then
       if (LDIFFDIAG) then
          wtmix(:,:,1)=s1o(:,:,1)
          rinu(:,:, 1) = 0._wp
       endif
    endif


!$omp do
    DO j=2,je1

       DO k=2,ke
          ku=MIN(k+1,ke)
          ko=MAX(k-1,1)
          stabeps=cstabeps/dzw(ko)
          wtdecay=EXP(-dzw(ko)/wpendep)

          topnah = 0.0_wp
          DO i=2,ie1

             t1o(i,j,1)=t1o(i,j,1)*wtdecay*stabeps                            &
                  / (stabeps + 0.5_wp * (stabio(i, j, k) + (1._wp - zsurf(ko)) &
                  &                      * stabio(i, j, ko)))

             s1o(i,j,1)=s1o(i,j,1)*wtdecay*stabeps                            &
                  / (stabeps + 0.5_wp * (stabio(i, j, k) + (1._wp - zsurf(ko)) &
                  &                      * stabio(i, j, ko)))

             if (iMEAN.ne.0)then
                if (LDIFFDIAG) then
                   wtmix(i,j,k)=s1o(i,j,1)
                endif
             endif


             dudo=almzer + weto(i,j,k) * di(k)**2                                  &
                  * (   ( uko(i-1,j,k) - uko(i-1,j,ko) )**2                        &
                  + ( uko(i,j,k)   - uko(i,j,ko)   )**2                        &
                  + ( vke(i,j-1,k) - vke(i,j-1,ko) )**2                        &
                  + ( vke(i,j,k)   - vke(i,j,ko)   )**2   ) * 0.5_wp

             hho=weto(i,j,k)*(amsuo(i,j,k)+amsuo(i-1,j,k)+amsue(i,j-1,k)       &
                  + amsue(i, j, k)) * 0.25_wp

             rinumo = hho * MAX(gor * stabio(i, j, k) / dudo, 0._wp)

             if (iMEAN.ne.0)then
                if (LDIFFDIAG) then
                   rinu(i,j,k)=rinumo
                endif
             endif

             avohelp = avo(i,j,k)
             avohelp = (relax*MIN(avohelp,av0+abackv(k))+relne*(t1o(i,j,1) &
                  + av0 / ((1._wp + cra*rinumo)**2) + abackv(k) + topnah)) &
                  & * weto(i,j,k)

             dvohelp = dvo(i,j,k)

             dvohelp = (relax*MIN(dvohelp,dv0+s1o(i,j,1))+relne*(s1o(i,j,1) &
                  + dv0 / ((1._wp + crd*rinumo)**3) + dbackv(k) + topnah))  &
                  & * weto(i,j,k)

             !! enhance vertical eddy diffusivity and viscosity for ioconv = 1,4
             !! do not enhance vert. eddy diffusivity and viscosity  for ioconv = 2,3
             !! cdvocon and cavocon can be set in the namelist


             SELECT CASE (ioconv)

             CASE(1)     !! convection by vertical diffusion (1)

               ! FIXME: why not use almzer here?
                dvo(i,j,k) = weto(i, j, k) &
                     * MAX(cdvocon * (1.e-11_wp - stabio(i,j,k)) &
                     &     / (1.e-11_wp + ABS(stabio(i, j, k))), dvohelp)
                avo(i,j,k) = weto(i, j, k) &
                     * MAX(cavocon * (1.e-11_wp - stabio(i,j,k)) &
                     &     / (1.e-11_wp + ABS(stabio(i, j, k))), avohelp)

             CASE(2,3,4)     !! complete mixing (2), interchange of tho and sao (3), plume convection (4)

                dvo(i,j,k) = dvohelp
                avo(i,j,k) = avohelp

             CASE default !! unsupported value for ioconv

                WRITE(0,*) 'ioconv value unsupported:', ioconv
                CALL P_ABORT

             END SELECT

             stabio(i, j, k) = MAX(stabio(i, j, k), 0._wp)

          ENDDO
       ENDDO
    ENDDO ! j-loop
!$omp end do


!$omp single
    CALL bounds_exch(1,'p',s1o(:,:,1),'mo_octher 10')
    CALL bounds_exch(1,'p',t1o(:,:,1),'mo_octher 11')
!#ifdef bounds_exch_save
    CALL bounds_exch(1,'p',dvo,'mo_octher 12')
    CALL bounds_exch(1,'p',avo,'mo_octher 13')
    CALL bounds_exch(1,'p',stabio,'mo_octher 14')
!#endif
!$omp end single

    avo(:,:,kep) = 0._wp
    dvo(:,:,kep) = 0._wp

!$omp end parallel

  END SUBROUTINE calc_rinum





end module mo_ocean_vertical_mixing
