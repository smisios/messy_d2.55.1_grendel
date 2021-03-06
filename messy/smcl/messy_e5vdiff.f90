MODULE MESSY_E5VDIFF
  ! AUTHOR:
  !  H.G. Ouwersloot, MPIC, April 2016

  USE messy_main_constants_mem, ONLY: dp
!!$#ifdef ECHAM5
!!$  USE messy_main_data_bi,       ONLY: cvdifts
  USE messy_main_constants_mem, ONLY: cvdifts
!!$#endif
  USE messy_main_constants_mem, ONLY: rd, rv

  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL(dp) :: clam    = 150._dp        !  *asymptotic mixing length for momentum.
  REAL(dp) :: ckap    = 0.4_dp         !  *karman constant.
  REAL(dp) :: cb      = 5._dp          !  *stability parameter near neutrality.
  REAL(dp) :: cc      = 5._dp          !  *stability parameter for unstable cases.
  REAL(dp) :: cchar   = 0.018_dp       !  *charnock constant.
!!$#ifndef ECHAM5
!!$  REAL(dp) :: cvdifts = 1.5_dp         !  *factor for timestep weighting
!!$                                       !   in *rhs* of *vdiff* and *scv*.
!!$#endif
  REAL(dp) :: cfreec  = 0.001_dp       !  *free convection parameter
  REAL(dp) :: cgam    = 1.25_dp        !  *free convection parameter
  REAL(dp) :: cz0ice  = 0.001_dp       !  *roughness over sea-ice
  REAL(dp) :: csncri  = 5.85036E-3_dp  !  *critical snow depth for soil computations

  ! Constants used for computation of saturation mixing ratio
  ! over liquid water (*c_les*) or ice(*c_ies*)
  ! Based on Murray (Journal of Applied Meteorology, 1967)
  REAL(dp), PARAMETER :: c1es  = 610.78_dp
  REAL(dp), PARAMETER :: c2es  = c1es*rd/rv
  REAL(dp), PARAMETER :: c3les = 17.269_dp
  REAL(dp), PARAMETER :: c3ies = 21.875_dp
  REAL(dp), PARAMETER :: c4les = 35.86_dp
  REAL(dp), PARAMETER :: c4ies = 7.66_dp

  ! Variables for computation of evapotranspiration
  REAL(dp), PARAMETER :: cva   = 5000._dp
  REAL(dp), PARAMETER :: cvb   = 10._dp
  REAL(dp), PARAMETER :: cvc   = 100._dp
  REAL(dp), PARAMETER :: cvbc  = cvb*cvc
  REAL(dp), PARAMETER :: cvabc = (cva+cvbc)/cvc
  REAL(dp), PARAMETER :: cvk   = .9_dp
  REAL(dp), PARAMETER :: cvkc  = cvk*cvc
  REAL(dp), PARAMETER :: cvrad = 0.55_dp

  PUBLIC   :: dp

  PUBLIC   :: clam
  PUBLIC   :: ckap
  PUBLIC   :: cb
  PUBLIC   :: cc
  PUBLIC   :: cchar
  PUBLIC   :: cvdifts
  PUBLIC   :: cfreec
  PUBLIC   :: cgam
  PUBLIC   :: cz0ice
  PUBLIC   :: csncri
  PUBLIC   :: c2es
  PUBLIC   :: c3les
  PUBLIC   :: c3ies
  PUBLIC   :: c4les
  PUBLIC   :: c4ies
  PUBLIC   :: cva
  PUBLIC   :: cvb
  PUBLIC   :: cvc
  PUBLIC   :: cvbc
  PUBLIC   :: cvabc
  PUBLIC   :: cvk
  PUBLIC   :: cvkc
  PUBLIC   :: cvrad

  PUBLIC   :: e5vdiff_surftemp

  CHARACTER(len=*), PUBLIC, PARAMETER :: MODSTR='e5vdiff'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver='1.1'

  CONTAINS

!=========================================================================
  SUBROUTINE e5vdiff_surftemp(klon, pdt,                                    &
! CONSTANTS
          pemi, pboltz, pcp, pc16, platev, platsu,                        &
! COEFFICIENTS FROM THE ELIMINATION
          pfscoe, pescoe, pfqcoe, peqcoe,                                 &
! OLD VALUES AT THE SURFACE
          psold, pqsold, pdqsold,                                         &
! OTHER FLUXES
          pnetrad, pgrdfl,                                                &
! DIFFUSION COEFFICIENTS, CAIR AND CSAT FOR EVAP AND SOIL HEAT CAPACITY
          pcfh, pcair, pcsat, pfracsu, pgrdcap,                           &
! Logical land mask
          lpland,                                                         &
! OUTPUT
          psnew, pqsnew)
!
! COMPUTES THE ENERGY BALANCE AT THE SURFACE WITH AN IMPLICIT SCHEME
! THAT IS CONNECTED TO THE RICHTMYER AND MORTON ALGORITHM OF THE PBL.
!
! INPUT
! -----
! KLON     : LENGTH OF ARRAYS TO BE USED
! PDT      : LEAP-FROG TIMESTEP IN SECONDS TIMES ALPHA (1.5)
!            ALPHA USED TO COMPENSATE FOR EVALUATING X* INSTEAD OF
!            X^(t+1); X* = ALPHA X^(t+1) + ( 1 - ALPHA ) X^(t-1)
!
! PEMI     : SURFACE EMISSIVITY
! PBOLTZ   : STEFAN-BOLTZMANN CONSTANT
! PCP      : SPECIFIC HEAT OF AIR
! PC16     : CPD*VTMPC2=CPD*(DELTA-1) (CF. VDIFF),FOR SENS.HEAT FL.
! PLATEV   : LATENT HEAT OF EVAPORATION
! PLATSU   : LATENT HEAT OF SUBLIMATION
!
! PFSCOE, PESCOE : COEFFICIENTS OF THE RICHTMYER AND MORTON SCHEME
!                  FOR DRY STATIC ENERGY
! PFQCOE, PEQCOE : AS ABOVE BUT FOR SPECIFIC HUMIDITY
!
! PSOLD   : OLD SURFACE DRY STATIC ENERGY (TS * CP)
! PQSOLD  : SATURATED  SPECIFIC HUMIDITY FOR OLD TEMPERATURE
! PDQSOLD : DERIVATIVE OF SATURATED  SPECIFIC HUMIDITY AT THE
!           OLD TEMPERATURE
!
! PNETRAD : NET RADIATION AT THE SURFACE (UPWARD LONGWAVE IS
!           INCLUDED BUT FOR THE OLD SURFACE TEMPERATURE)
! PGRDFL  : GROUND HEAT FLUX
!
! PCFH    : DIFFUSION COEFFICIENT FOR STATIC ENERGY AND MOISTURE
! PCAIR   : COEFFICIENT IN LATENT HEAT FLUX FORMULA (SEE VDIFF)
! PCSAT   : COEFFICIENT IN LATENT HEAT FLUX FORMULA (SEE VDIFF)
! PFRACSU : FRACTION OF SURFACE FOR SUBLIMATION
! PGRDCAP : SURFACE HEAT CAPACITY
!
! OUTPUT
! ------
! PSNEW   : NEW SURFACE STATIC ENERGY
! PQSNEW  : NEW SATURATED SURFACE AIR MOISTURE
!
!
! AUTHOR.
! -------
!
! J. POLCHER  *LMD*  AND  J.-P. SCHULZ  *MPI*,  MAY 1995
! H.G. OUWERSLOOT *MPIC*, 2016
!
!
! MODIFICATIONS.
! --------------
!
! J.-P. SCHULZ  *MPI*,  OCTOBER 1997:
!    MODIFY ACCORDING TO LATENT HEAT FLUX FORMULATION IN VDIFF
!    USING ZCAIR AND ZCSAT COEFFICIENTS.
!
! J.-P. SCHULZ  *MPI*,  AUGUST 1998:
!    MODIFY ACCORDING TO SENSIBLE HEAT FLUX FORMULATION IN VDIFF.
!
! PF.COE IS ADAPTED TO EVALUATE X INSTEAD OF X/ALPHA IN EQUATIONS
!
! OUWERSLOOT *MPIC*, 2016:
! TRANSFER FROM ECHAM TO MESSY AND CLARIFICATIONS
!
    IMPLICIT NONE

!   ARGUMENTS
    INTEGER :: klon
    REAL(dp):: pdt, pemi, pboltz, pcp(klon), pc16, platev, platsu
    REAL(dp):: pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
    REAL(dp):: psold(klon), pqsold(klon), pdqsold(klon)
    REAL(dp):: pnetrad(klon), pgrdfl(klon)
    REAL(dp):: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
    REAL(dp):: pgrdcap(klon)
    REAL(dp):: psnew(klon), pqsnew(klon)
    LOGICAL :: lpland(klon)

    INTEGER :: jl
    REAL(dp):: zcolin, zcohfl, zcoind, zicp, zca, zcs

!   MAIN ROUTINE
    DO jl = 1,klon
      IF (lpland(jl)) THEN
      ! Defining zcolin, zcohfl and zcoind such that at the surface
      ! s^{new} = ( zcolin s^{old} + zcoind) / ( zcolin + zcohfl )
      ! (Cs/cp)*(s^{new}-s^{old})/(alpha*2*Delta{T} =
      !   R_n^{impl.} + LE^{impl.} + H^{impl.} + G + X, (fluxes pointing to the surface)
      ! See Eq. (24) of Schulz et al. (2001); extra term X
      ! X = - ( LE^{impl.} / L ) * (cpv - cpd) T_S^{old}
      ! Compensates the use of dry static energy flux, H, (based on s_S and s_air) instead of
      ! directly using SH itself (based on T_S and T_air):
      !   w's' = w'(cp T)' => rho cp w'T' = rho w's' - rho T w'cp'
      !   w'cp' = w'q' (cpv - cpd)
      !   => SH = H - LE / L * T * (cpv - cpd)

        zicp = 1._dp/pcp(jl)                                                ! 1 / cp

        zca    = platsu*pfracsu(jl) +  platev*(pcair(jl) - pfracsu(jl))     ! L * beta
        zcs    = platsu*pfracsu(jl) +  platev*(pcsat(jl) - pfracsu(jl))     ! L * beta * h

        zcolin = pgrdcap(jl)*zicp +                                       & ! Cs / cp + 2 alpha Delta{t} * ( 1/cp * 4*epsilon*sigma*(s_old/c_p)^3
                          pdt*(zicp*4._dp*pemi*pboltz*                    & !        - rho C_H |U| ( L beta E_q - L beta h
                          ((zicp*psold(jl))**3._dp) -                     & !             - 1/cp * (cpv - cpd) s^{old} ( beta E_q - beta h) )
                          pcfh(jl)*(zca*peqcoe(jl) - zcs -                & !        * 1/cp d{q_s}/d{T} )
                                    zicp*pc16*psold(jl)*                  &
                                    (pcair(jl)*peqcoe(jl) - pcsat(jl)))*  &
                               zicp*pdqsold(jl))

        zcohfl = -pdt*pcfh(jl)*(pescoe(jl)-1._dp)                           ! - 2 alpha Delta{t} rho C_H |U| ( E_s - 1)

        zcoind = pdt * (pnetrad(jl) + pcfh(jl)*pfscoe(jl) +  pcfh(jl)*    & ! 2 alpha Delta{t} * (R_n^{expl.} + rho C_H |U| F_s +
                          ((zca*peqcoe(jl)-zcs)*pqsold(jl) +              & !     rho C_H |U| ( ( L beta E_q - L beta h ) q_s^{old}
                          zca*pfqcoe(jl) - zicp*pc16*psold(jl)*           & !         + L beta F_q - 1/cp * (cpv - cpd) s^{old} *
                         ((pcair(jl)*peqcoe(jl) - pcsat(jl))*pqsold(jl) + & !             ( ( beta E_q - beta h) q_s^{old} + beta F_q ) )
                             pcair(jl)*pfqcoe(jl))) + pgrdfl(jl))           !     + G )

        psnew(jl) = (zcolin * psold(jl) + zcoind) / (zcolin + zcohfl)
        pqsnew(jl) = pqsold(jl) + zicp * pdqsold(jl) * (psnew(jl) -       & ! q_s,new = q_s,old + (1 / cp) * d{q_s}/d{T} * (s_new - s_old)
                                                              psold(jl))
      ELSE
        psnew(jl) = psold(jl)
        pqsnew(jl) = pqsold(jl)
      END IF
    END DO
    RETURN
  END SUBROUTINE e5vdiff_surftemp
!=========================================================================
END MODULE MESSY_E5VDIFF
