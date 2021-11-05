!>
 !! Routines for ocean thermodynamics. This includes the following subroutines (in order they are called):
 !! dilcor_gtrf2  : dilutes tracers, only if coupled to biogeochemistry model HAMOCC
 !! dilcor_ptrf2  : dilutes tracers, only if coupled to biogeochemistry model HAMOCC
 !! calc_dens     : calculates density
 !!
 !! @author: Helmuth Haak, Max Planck Institute for Meteorology
 !!

 MODULE mo_octher

  USE mo_param1
  USE mo_param3      ! This is not actually used here but only in subroutine adisitj which is called from octher
                     ! and is supposed to be expanded inline. The purpose of this dependency is to force compilation
                     ! of mo_param3.f90 before compilation of mo_octher.f90, because mo_param3.mod is needed for
                     ! the inline expansion of adisitj.

  USE mo_parallel
  USE mo_commo1
  USE mo_commoau1
  USE mo_commoau2
  USE mo_levitus, only: relax_surf
  USE mo_units, only: io_stdout

  USE mo_runoff, ONLY : river_runoff_stations,river_runoff_omip,glac_calv, &
                        luse_river_runoff_stations,luse_glac_calv

#ifdef PBGC
  USE mo_tro ,only: dilcor_gtrf2 , dilcor_ptrf2
  USE mo_param1_bgc, only: nocetra
  USE mo_carbch, only: ocetra
#endif


  IMPLICIT NONE

 CONTAINS

  SUBROUTINE octher

  USE mo_diagnosis, only: calc_potential_energy_release
  USE mo_convection, only : convection
  USE mo_ocean_vertical_mixing, only : calc_rinum


#ifdef PBGC
  INTEGER :: l
#endif
  !-----------------------------------------------------------------------
  !
  !     sbr octher computes
  !
  !          river_runoff -  only in the uncoupled model
  !          convection   -  baroclinic pressure in each layer +  convective adjustment
  !          calc_rinum   -  richardson-number depending coefficients for
  !                          vertical diffusion of momentum (avo) and
  !                          temperature and salinity (dvo)
  !          calc_den     -  update of density field
  !
  !
  !
  !-----------------------------------------------------------------------

#ifndef __coupled

#ifdef PBGC
    call dilcor_gtrf2
#endif

    CALL relax_surf      ! only in the uncoupled model

    IF ( luse_river_runoff_stations ) THEN
      CALL river_runoff_stations
    ELSE
      CALL river_runoff_omip
    ENDIF

    IF ( luse_glac_calv ) THEN
      CALL glac_calv
    ENDIF

#ifdef PBGC
    do l=1,nocetra
       call dilcor_ptrf2(ocetra(:,:,1,l))
    enddo
#endif

#endif

    IF (LCONVDIAG .AND. ioconv .NE. 1 ) CALL calc_potential_energy_release(1)

    CALL convection

    IF (LCONVDIAG .AND. ioconv .NE. 1 ) CALL calc_potential_energy_release(2)

    CALL calc_rinum

    IF  (ioconv .NE. 1 ) THEN
      CALL calc_dens
    ENDIF

  END SUBROUTINE octher




  SUBROUTINE calc_dens

    REAL(wp) :: shelp(ie,je),thelp(ie,je),rhelp(ie,je)
    INTEGER :: i,j,k

!$omp parallel private(i,j,k)
!$omp do
    DO j=1,je

!       DO i=1,ie
!          shelp(i,j)=0.0
!          thelp(i,j)=0.0
!       ENDDO

       !uwe  include downward propagation of tho and sao in land

       DO k=2,ke
          DO i=1,ie
             IF(.NOT. lweto(i,j,k)) THEN
                tho(i,j,k)=tho(i,j,k-1)
                sao(i,j,k)=sao(i,j,k-1)
             ENDIF
          ENDDO
       ENDDO

       !uwe  compute pressure and density new

       DO k=1,ke
          DO i=1,ie
             thelp(i,j)=tho(i,j,k)
             shelp(i,j)=sao(i,j,k)
          ENDDO
          CALL adisitj(thelp,shelp,preff(k),j)
          CALL rho1j(thelp,shelp,preff(k),rhelp,j)

          !uwe   po noch sauber bestimmen, aus rhoo noch machen!!!

          DO i=1,ie
             rhoo(i,j,k)=rhelp(i,j)
          ENDDO
       ENDDO
    ENDDO
!$omp end do
!$omp end parallel

  END SUBROUTINE calc_dens

 END MODULE mo_octher
