!***********************************************************************
!
!**** *MODULE mo_bgc_diagnostic* -
!          allocates variables for bgcflux field (2d, for surface fluxes, as in bgcmean_2d)
!          allocates variables for bgcprod field (3d, for fluxes within euphotic zone only)
!
!
!     Katharina Six    *MPI-Met, HH*    9.12.09
!
!     Purpose
!     -------
!     - declaration and memory allocation for output fields (CDI Output)
!
!**********************************************************************
MODULE mo_bgc_diagnostic

  USE mo_carbch
  USE mo_commo1
  USE mo_param1
  USE mo_param1_bgc


  IMPLICIT NONE

  ! 3d fields -fluxes within euph. zone-

  INTEGER, PARAMETER :: &
       ibgcprod = 8, &
       kphosy   = 1, &
       kgraz    = 2, &
       kexport  = 3, &
       kdelcar  = 4, &
       kdelsil  = 5, &
       kdmsprod = 6, &
       kdms_bac = 7, &
       kdms_uv  = 8

  INTEGER, PARAMETER :: &
       ibgc_rate = 6, &
       ksco212  =1,  &
       kalkali  =2,  &
       kphosph  =3,  &
       kano3    =4,  &
       ksilica  =5,  &
       kiron   = 6


  ! 2d fields -fluxes-

  INTEGER, PARAMETER :: &
       ibgcflux  = 38, &
       kco2flux  =  1, &
       kco214f   =  2, &
       ko2flux   =  3, &
       kn2flux   =  4, &
       kn2oflux  =  5, &
       kdmsflux  =  6, &
       kprorca   =  7, &
       kprcaca   =  8, &
       ksilpro   =  9, &
       kprodus   = 10, &
       kkwco2    = 11, &
       kpco2     = 12, &
       kco2fxd   = 13, &
       kco2fxu   = 14, &
       kcoex90   = 17, &
       kopex90   = 18, &
       kcaex90   = 19, &
       kcoex1000 = 20, &
       kopex1000 = 21, &
       kcaex1000 = 22, &
       kcoex2000 = 23, &
       kopex2000 = 24, &
       kcaex2000 = 25, &
  !#ifdef ANTC14
       kac14fx   = 26, &
  !#ifdef __c_isotopes
       kc13flux  = 27, &
       kc14flux  = 28, &
  !#ifdef PCFC
       kcfc11fx  = 29, &
       kcfc12fx  = 30, &
       kpcfc11   = 31, &
       kpcfc12   = 32, &
  !#ifdef AMMO
       knh3flux  = 33, &
       kn2fix    = 34, &
       klysokl   = 35, &
       kfeatm    = 36, &
       kdpco2    = 37, &
       kdpo2     = 38

  INTEGER, PARAMETER :: &
       ibgcomz  = 2, &
       komz      = 1,  &
       komz_depth= 2

  ! 3d fields -fluxes whole ocean-

  INTEGER, PARAMETER :: &
       ibgc_o_pro = 2, &
       kdissol   = 1, &
       kdenit  = 2

  REAL(wp), DIMENSION (:,:,:,:), ALLOCATABLE, TARGET :: bgcprod    ! 3d concentration EU
  REAL(wp), DIMENSION (:,:,:), ALLOCATABLE, TARGET :: bgcflux      ! 2d flux field
  REAL(wp), DIMENSION (:,:,:,:), ALLOCATABLE, TARGET :: bgc_o_pro  ! 3d fluxes whole ocean
  REAL(wp), DIMENSION (:,:,:), ALLOCATABLE, TARGET :: bgcomz     ! 2d O2-concentration and depth of OMZ
  REAL(wp), DIMENSION (:,:,:,:), ALLOCATABLE, TARGET :: bgcrate     ! 3d rate change within EU

CONTAINS

  SUBROUTINE bgc_diagnostic_init(io_stdo_bgc, kpie, kpje, kpke, kwrbioz)
    INTEGER, INTENT(in) :: io_stdo_bgc, kpie, kpje, kpke, kwrbioz

    WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcprod ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',kwrbioz
    WRITE(io_stdo_bgc,*)'Fourth dimension    : ',ibgcprod

    ALLOCATE(bgcprod(kpie, kpje, kwrbioz, ibgcprod))
    bgcprod(:,:,:,:) = 0.0_wp


    WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcflux ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',ibgcflux

    ALLOCATE(bgcflux(kpie,kpje,ibgcflux))
    bgcflux(:,:,:) = 0.0_wp


    WRITE(io_stdo_bgc,*)'Memory allocation for variable bgc_o_pro ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',kpke
    WRITE(io_stdo_bgc,*)'Fourth dimension    : ',ibgc_o_pro

    ALLOCATE(bgc_o_pro(kpie,kpje,kpke,ibgc_o_pro))
    bgc_o_pro(:,:,:,:) = 0.0_wp


    WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcomz ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',ibgcomz

    ALLOCATE(bgcomz(kpie,kpje,ibgcomz))
    bgcomz(:,:,:) = 0.0_wp

    WRITE(io_stdo_bgc,*)'Memory allocation for variable bgcrate ...'
    WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
    WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
    WRITE(io_stdo_bgc,*)'Third dimension    : ',kwrbioz
    WRITE(io_stdo_bgc,*)'Fourth dimension    : ',ibgc_rate

    ALLOCATE(bgcrate(kpie,kpje,kwrbioz,ibgc_rate))
    bgcrate(:,:,:,:) = 0.0_wp

  END SUBROUTINE bgc_diagnostic_init

  !>
  !! Compute minimum O2 concentration in water column and corresponding depth.
  !!
  SUBROUTINE calc_omz_depth(pddpo)

    REAL(wp), INTENT(IN) :: pddpo(ie,je,ke)
    REAL(wp) :: sumdep

    INTEGER :: i, j, k

    bgcomz = 100._wp

    DO j = 1,je
       DO i = 1,ie
          sumdep = 0._wp
          DO k = 1,ke
             IF (pddpo(i, j, k) > 0.5_wp) THEN
                IF(ocetra(i,j,k,ioxygen) < bgcomz(i,j,komz))  THEN
                   bgcomz(i,j,komz) = ocetra(i,j,k,ioxygen)
                   sumdep = sumdep + pddpo(i,j,k)
                END IF
             END IF
          END DO
          IF(pddpo(i,j,1) > 0.5_wp) bgcomz(i,j,komz_depth) = sumdep
       END DO
    END DO

  END SUBROUTINE calc_omz_depth

  !>
  !! Initializer tracer concentration for rate computation.
  !!
  SUBROUTINE store_tracer(pddpo)

    REAL(wp), INTENT(IN) :: pddpo(ie,je,ke)
    INTEGER :: i,j,k

    DO k = 1,kwrbioz
       DO j = 1,je
          DO i = 1,ie
             IF(pddpo(i,j,k) > 0.5_wp) THEN
                bgcrate(i,j,k,ksco212) = ocetra(i,j,k,isco212)
                bgcrate(i,j,k,kalkali) = ocetra(i,j,k,ialkali)
                bgcrate(i,j,k,kphosph) = ocetra(i,j,k,iphosph)
                bgcrate(i,j,k,kano3) = ocetra(i,j,k,iano3)
                bgcrate(i,j,k,ksilica) = ocetra(i,j,k,isilica)
                bgcrate(i,j,k,kiron) = ocetra(i,j,k,iiron)
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE store_tracer

  !>
  !! Get rate of tracer concentration as simple difference.
  !!
  SUBROUTINE rate_tracer(pddpo)

    REAL(wp), INTENT(IN) :: pddpo(ie,je,ke)
    INTEGER :: i,j,k

    DO k=1,kwrbioz
       DO j=1,je
          DO i=1,ie
             IF (pddpo(i, j, k) .GT. 0.5_wp) THEN
                bgcrate(i,j,k,ksco212)= ocetra(i,j,k,isco212) - bgcrate(i,j,k,ksco212)
                bgcrate(i,j,k,kalkali)= ocetra(i,j,k,ialkali) - bgcrate(i,j,k,kalkali)
                bgcrate(i,j,k,kphosph)= ocetra(i,j,k,iphosph) - bgcrate(i,j,k,kphosph)
                bgcrate(i,j,k,kano3)  = ocetra(i,j,k,iano3) - bgcrate(i,j,k,kano3)
                bgcrate(i,j,k,ksilica)= ocetra(i,j,k,isilica) - bgcrate(i,j,k,ksilica)
                bgcrate(i,j,k,kiron)  = ocetra(i,j,k,iiron) - bgcrate(i,j,k,kiron)
             END IF
          END DO
       END DO
    END DO

  END SUBROUTINE rate_tracer

END MODULE mo_bgc_diagnostic
