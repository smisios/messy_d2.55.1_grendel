!>
!! SWR-ABSORPTION SCHEME
!! Shortwave radiation is done in three steps
!! As first step the relative swr absorption factor (swrab) in each level is calculated
!! either for
!! a) a static jerlov type of water (default Type II)
!!    sbr jerlov_swr_absorption
!!
!! or
!! b) derived from a dynamic 3d cholophyll map taken from hamocc
!!    sbr swr_frac_from_chl
!!    sbr dynamic_swr_absorption
!!
!! The actual absorption and the update of the temperature is afterwards done
!! again in two steps.
!! a) for the uppermost level in sbr growth together with the rest of the
!!    thermodynamic forcing.
!! b) for the subsurface levels in subsurface_swr_absorption
!!    Heat is distributed over the water column in proportion to the absorption.
!!    No Heat is transpered into sea floor, but is added to the bottom level.
!!
!! @author H.Haak, MPI-M
!! @author U.Mikolajewicz, MPI-M
!!
!!
!! @par Revision History
!! Introduced by H.Haak, MPI-M (2010-03-12)
!!
!!
MODULE mo_swr_absorption

  USE mo_kind, ONLY: wp, dp
  USE mo_param1, ONLY : ie,je,ke
  USE mo_commo1, ONLY : tho, dzw, tiestw, weto, dt, dpio, zo, ddpo
  USE mo_commoau1, ONLY : cc
  USE mo_grid, ONLY : lwith_one_layer_shelfs



  IMPLICIT NONE

  PRIVATE


! Jerlov Water Type values from Kara et. al. (2005)
!   - namelist parameter
!   - default Type II
!
! type   jerlov_atten  jerlov_bluefrac
! ------------------------------------
!  I       0.05            0.45
!  IA      0.06            0.41
!  IB      0.08            0.36
!  II      0.12            0.28
!  III     0.17            0.27


  REAL(wp) ::  jerlov_atten       !<  attenuation of photosyntheic active radiation (kpar) [m-1]
  REAL(wp) ::  jerlov_bluefrac     !< blue fraction of par

  REAL(wp), ALLOCATABLE ,TARGET :: swr_frac(:,:,:) !<  available light fraction passed to hamocc
  REAL(wp), ALLOCATABLE ,TARGET :: swsum(:,:)      !<  swr fraction absorbed in the surface layer
  REAL(wp), ALLOCATABLE ,TARGET :: swrab(:,:,:)    !<  relative swr absorption factor
  REAL(wp), ALLOCATABLE ,TARGET :: heatabs(:,:)    !<  heating due to absorption [J m-2]

  LOGICAL :: lfb_bgc_oce=.FALSE.


  PUBLIC :: swsum, swrab, heatabs, swr_frac         &
       , jerlov_swr_absorption                      &
       , old_swr_absorption                         &
       , dynamic_swr_absorption                     &
       , subsurface_swr_absorption                  &
       , alloc_mem_swr_absorption                   &
       , jerlov_atten,jerlov_bluefrac,lfb_bgc_oce   &
       , swr_frac_from_chl


CONTAINS

  SUBROUTINE alloc_mem_swr_absorption

    !>
    !! This is the initialisation part of the sw-absorption scheme
    !! ( should be called before the timeloop )

    ALLOCATE(swsum(ie,je), swrab(ie,je,ke), &
         heatabs(ie,je))

    swsum = 0._wp
    swrab = 0._wp
    heatabs = 0._wp

#ifdef PBGC
    ALLOCATE(swr_frac(ie,je,ke+1))
    swr_frac(:,:,:)=0.0_wp
#endif

  END SUBROUTINE alloc_mem_swr_absorption


  !>
  !! This is the old sw-absorption scheme and is only here for backward compatibility
  !! It should not be used anaymore for prduction

  SUBROUTINE  old_swr_absorption

    INTEGER :: k

    REAL(wp), PARAMETER :: opendep=11._wp
    REAL(wp) :: swsumi(ie,je)

    swsum(:,:)=EXP(-tiestw(2)/opendep)
    swsumi(:,:) = 1._wp / swsum(:,:)

    DO k=1,ke
      swrab(:,:,k)=swsumi(:,:)*(EXP(-tiestw(k)/opendep)    &
           -EXP(-tiestw(k+1)/opendep))
    ENDDO

  END SUBROUTINE old_swr_absorption




  !>
  !! This is the first part of the sw-absorption scheme
  !! ( can to be called before the time loop )
  !! This sbr calculates the relative swr absorption factor (swrab) in each level.
  !!
  !! Numbers are derived for jerlov  type of water
  !!    N.G. Jerlov: Optical Studies of Ocean Waters.
  !!             In: Reports of the Swedish Deep-Sea Expedition, 3:1-59; 1951.
 

  SUBROUTINE  jerlov_swr_absorption

    REAL(wp) :: swsumi(ie,je)

    INTEGER :: k

    ! length scale for penetration of sw-radiation [m]
    !! new sw-penetration
    !! default parameters for JERLOV II

    swsum(:,:)=jerlov_bluefrac*EXP(-tiestw(2)*jerlov_atten)
    swsumi(:,:) = 1._wp / swsum(:,:)

    DO k=1,ke
      swrab(:,:,k)=swsumi(:,:)*jerlov_bluefrac*(EXP(-tiestw(k)*jerlov_atten)-EXP(-tiestw(k+1)*jerlov_atten))
    ENDDO

  END SUBROUTINE jerlov_swr_absorption


  !>
  !! This is the first part of the sw-absorption scheme
  !! ( needs to be called directly before upper layer thermodynamics -> growth )
  !! This sbr calculates the relative swr absorption factor (swrab) in each level.
  !!
  !! Availble light fraction (swr_frac) is derived from a 3d chlorophyll map taken 
  !! from hamocc ( sbr swr_frac_from_chl ; LFB_BGC_OCE=.true.)


  SUBROUTINE  dynamic_swr_absorption

    REAL(wp) :: swsumi(ie,je)

    INTEGER :: k

    swsum(:,:)=swr_frac(:,:,2)
    swsumi(:,:)=1._wp /swsum(:,:)

    DO k=1,ke
      swrab(:,:,k)=swsumi(:,:) * (swr_frac(:,:,k) - swr_frac(:,:,k+1))
    ENDDO


  END SUBROUTINE dynamic_swr_absorption

  SUBROUTINE swr_frac_from_chl

#ifdef PBGC
    USE mo_carbch, ONLY : ocetra
    USE mo_param1_bgc, ONLY : iphy


    !> parameters for sw-radiation fraction
    !! Analogue to Zielinski et al., Deep-Sea Research II 49 (2002), 3529-3542

    REAL(wp), PARAMETER :: redfrac=0.4_wp !< red fraction of the spectral domain (> 580nm)

    REAL(wp), PARAMETER :: c_to_chl=12.0_wp/60.0_wp   !< ration Carbon to Chlorophyll
    REAL(wp), PARAMETER :: r_car=122.0_wp   !< Redfield ratio
    REAL(wp), PARAMETER :: pho_to_chl=r_car*c_to_chl*1.e6_wp !< 1 kmolP = (122*12/60)*10^6 mg[Chlorophyll]

    REAL(wp), PARAMETER :: atten_r=0.35_wp !< attenuation of red light [m-1]
    REAL(wp), PARAMETER :: atten_w=0.03_wp !< attenuation of blue/green light
                                           !! in clear water between 400nm and 580nm [m-1]
    REAL(wp), PARAMETER :: atten_c=0.04_wp !< attenuation of blue/green light
                                           !! by chlorophyll [m-1]


    INTEGER :: k

    swr_frac(:,:,1) = 1.0_wp

    DO k=2,ke
      swr_frac(:,:,k) = swr_frac(:,:,k-1) * (                         &
          redfrac         * EXP(-dzw(k-1) *  atten_r)                 &
      +  (1.0_wp-redfrac) * EXP(-dzw(k-1) * (atten_w +                &
           atten_c*pho_to_chl*MAX(0.0_wp,ocetra(:,:,k-1,iphy)))))
    END DO

#endif

  END SUBROUTINE swr_frac_from_chl


  SUBROUTINE subsurface_swr_absorption

    !>
    !! This is part two of the sw-absorption scheme and should be called after
    !! the surface thermodynamics -> growth )
    !! This sbr calculates the subsurface absorption from the botton to level 2.
    !! Absorption in the surface level is done in sbr growth.
    !! SWR does not penetrate into Land. Heat from SWR that would penetrate theoreticly
    !! below the bottom is added to the bottom layer.


    INTEGER :: k

    REAL(wp) :: heatabb(ie,je) !> absorbed heat, accumulates
                               !! the heat that would be absorbed in
                               !! grid points below the bottom
    REAL(wp) :: heatabs_t(ie,je)

    heatabs_t(:,:)=weto(:,:,1)*heatabs(:,:)*dt/cc
    heatabb(:,:) = 0._wp

    DO k=ke,2,-1
      heatabb(:,:)=heatabb(:,:)+swrab(:,:,k)*heatabs_t(:,:)
      tho(:,:,k)=tho(:,:,k)+heatabb(:,:)*dpio(:,:,k)
      heatabb(:,:) = heatabb(:,:) * (1._wp - weto(:,:,k))
    END DO

    IF ( lwith_one_layer_shelfs ) THEN
      WHERE (weto(:,:,1) > 0.5_wp) tho(:,:,1) = tho(:,:,1) &
           +heatabb(:,:)/(ddpo(:,:,1)+zo(:,:))
    ENDIF



  END SUBROUTINE subsurface_swr_absorption


END MODULE mo_swr_absorption

