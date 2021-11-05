! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL H2OISO
!
! THIS SUBMODEL IS TO SIMULATE WATER ISOTOPES
! 
! Authors: * Roland Eichinger, DLR-IPA, 2013-2015
!          * Franziska Frank, DLR-IPA, 20??
!            - chemistry moved to submodel CH4
!            - submodel TRSYNC
!            - corrections for MSBM coupling
!          * Patrick Joeckel, DLR-IPA, 2020:
!            - messy_main_data_wiso_bi.g90
!            - messy_main_tools_wiso.f90
!            - cloud processes moved to CLOUD submodel
!
! References:
!
! * Eichinger, R., Jöckel, P., & Lossow, S.: Simulation of the isotopic
!   composition of stratospheric water vapour – Part 2: Investigation of
!   HDO / H2O variations, Atmospheric Chemistry and Physics, 15, 7003–7015,
!   doi: 10.5194/acp-15-7003-2015,
!   URL http://www.atmos-chem-phys.net/15/7003/2015/ (2015b)
!
! * Eichinger, R., Jöckel, P., Brinkop, S., Werner, M., & Lossow, S.:
!   Simulation of the isotopic composition of stratospheric water vapour
!   – Part 1: Description and evaluation of the EMAC model, Atmospheric
!   Chemistry and Physics, 15, 5537–5555, doi: 10.5194/acp-15-5537-2015,
!   URL http://www.atmos-chem-phys.net/15/5537/2015/ (2015a)
!
! **********************************************************************

! **********************************************************************
MODULE messy_h2oiso
! **********************************************************************

   USE messy_main_constants_mem,    ONLY: DP, tmelt, alf, cpd=>cp_air, &
                                          vtmpc2, vtmpc1, rd,          &
                                          rhoh2o=>rho_H2O
   USE messy_main_tools_wiso

   IMPLICIT NONE
   !PRIVATE
   SAVE
 
   PUBLIC :: DP

   CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'h2oiso'
   CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.2'

  interface fractcal
     module procedure fractcal
  end interface

  interface wiso_calc_delta
     module procedure wiso_calc_delta
  end interface

  !  SUBROUTINES
  ! VDIFF SUBROUTINES CALLED FROM BMIL
  PUBLIC :: h2oiso_vdiff_surfhum
  PUBLIC :: h2oiso_vdiff_vapcoeff
  PUBLIC :: h2oiso_vdiff_kinfac
  PUBLIC :: h2oiso_vdiff_RMscheme
  PUBLIC :: h2oiso_vdiff_zqdif
  PUBLIC :: h2oiso_vdiff_mflux
  ! SURF SUBROUTINES CALLED FROM BMIL
  PUBLIC :: h2oiso_surf_convflux
  PUBLIC :: h2oiso_surf_snowc
  PUBLIC :: h2oiso_surf_snsulw
  PUBLIC :: h2oiso_surf_snglme
  PUBLIC :: h2oiso_surf_snbmew
  PUBLIC :: h2oiso_surf_skinres
  PUBLIC :: h2oiso_surf_soilres
  PUBLIC :: h2oiso_surf_corrmp
  ! WISO SUBROUTINES, DIRECTLY COPIED
  PUBLIC :: fractcal
  PUBLIC :: wiso_calc_delta ! that is a function! (wiso_)

 CONTAINS

   ! =========================================================================
   ! ### own public subroutines
   ! =========================================================================

   ! =========================================================================

SUBROUTINE h2oiso_vdiff_surfhum( kproma                       &
                               , pwsmx   , ptsw    , ptsi     &
                               , zqsl    , zqsw    , zqsi     &
                               , zwisoqsw, zwisoqsi, zwisoqsl &
                               , pwisows , pwisosw_d   )

  IMPLICIT NONE


  INTEGER, INTENT(IN) :: kproma


  ! Local wiso variables
  REAL(dp):: zwisofracl ! Eq. fractionation 
  REAL(dp):: zwisofracw ! coefficients,
  REAL(dp):: zwisofraci ! eff. frac for ice
  REAL(dp):: zocean
  REAL(dp):: pwisosw_d(kproma,mwiso) ! from file
  REAL(dp):: zice
  REAL(dp):: zsatval

  ! Final results of block
  REAL(dp), INTENT(OUT) :: zwisoqsw(kproma,mwiso) ! iso. hum. sat. over water
  REAL(dp), INTENT(OUT) :: zwisoqsi(kproma,mwiso) ! iso. hum. sat. over ice
  REAL(dp), INTENT(OUT) :: zwisoqsl(kproma,mwiso) ! iso. hum. sat. over land - 
                                    ! not needed yet because fractionation
                                    ! processes over land not considered

  ! Wiso variables needed/coming from other processes 
  REAL(dp), INTENT(INOUT) :: pwisows(kproma,mwiso) ! from/to "surf"

  ! Variables from vdiff
    ! calculated in vdiff
  REAL(dp), INTENT(IN):: zqsl(kproma)  ! surf. sat. spec. humidity land
  REAL(dp), INTENT(IN):: zqsw(kproma)  ! hum. sat. over water
  REAL(dp), INTENT(IN):: zqsi(kproma)  ! hum. sat. over ice
    ! assigned to vdiff
  REAL(dp), INTENT(IN) :: pwsmx(kproma) ! Water holding capacity [m] of soil
  REAL(dp), INTENT(IN) :: ptsw(kproma)  ! surface temperature over water
  REAL(dp), INTENT(IN) :: ptsi(kproma)  ! surface temperature over ice

  INTEGER :: jt, jl
  
!    surface humidity for land, water and ice - water isotopes
     DO jt=1,mwiso
        DO jl=1,kproma !kidia,kfdia!
!    land ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!          no fractionation for land processes
           zwisofracl=1.0_dp
           zwisoqsl(jl,jt)=zwisofracl*zqsl(jl)
           pwisows(jl,jt)=MIN(pwisows(jl,jt),zwisofracl*pwsmx(jl))

!    water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!          fractionation over water
           zwisofracw=exp(talph1(jt)/(ptsw(jl)**2)+talph2(jt)/ptsw(jl)+talph3(jt))
!          "corrections" of ocean surface isotope values by prescribed seasurface delta values 
           zocean=toce(jt)* (1.0_dp+(pwisosw_d(jl,jt)/1000.0_dp))
           zwisofracw=zocean/zwisofracw
!          calculation of tracer humidity saturation
           zwisoqsw(jl,jt)=zwisofracw*zqsw(jl)
!
!    ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
           zwisofraci=1.0_dp
           IF (jt == i_HH18O) THEN ! vapour-liquid fractionation for H2-18O
              zwisofraci=exp(talps1(jt)/ptsi(jl)+talps2(jt))
           ENDIF
           IF (jt == i_HDO) THEN ! vapour-liquid fractionation for HD-16O
              zwisofraci=exp(talps1(jt)/(ptsi(jl)**2)+talps2(jt))
           ENDIF
!          effective fractionation over ice if necessary (isotopes only)
           IF ( (jt /= i_HHO) .and. (ptsi(jl) < tmelt)) THEN
              zsatval=tsatbase - tsatfac*(ptsi(jl)-tmelt)
              zwisofraci=zwisofraci*(zsatval/(1.0_dp+zwisofraci*(zsatval-1.0_dp)*tdifrel(jt)))
           ENDIF
!          "corrections" of seaice surface isotope values by prescribed seasurface delta values 
           zice=toce(jt)* (1.0_dp+(pwisosw_d(jl,jt)/1000.0_dp))
           zwisofraci=zice/zwisofraci
!          calculation of tracer humidity saturation
           zwisoqsi(jl,jt)=zwisofraci*zqsi(jl)

!if (isnan(zwisoqsi(jl,jt)))zwisoqsi(jl,jt)=zqsi(jl)*tnat(jt)!re

        END DO
     END DO


END SUBROUTINE h2oiso_vdiff_surfhum

   ! =================================

   ! =================================


   ! =================================


SUBROUTINE h2oiso_vdiff_vapcoeff( kproma, klev         &
                                , pum1, pvm1  , pqm1   &
                                , zqsl, zhum  , psn    &
                                , pws , pwsmx , pcvs   &
                                , pcvw, pvgrat, pwl    &
                                , zchl, zwet_tmp       &
                                , pwisows              &
                                , pwisosn              &
                                , pwisowl              &
                                , zwisocsat, zwisocair &
                                , lpland   )

  IMPLICIT NONE

  ! usable
  INTEGER, INTENT(IN)  :: kproma, klev
  LOGICAL, INTENT(IN)  :: lpland(kproma)       ! : land-sea flag
  REAL(dp), INTENT(IN) :: pum1(kproma,klev)
  REAL(dp), INTENT(IN) :: pvm1(kproma,klev)
  REAL(dp), INTENT(IN) :: pqm1(kproma,klev)

!--------
! preliminary result !!!!!
  REAL(dp), INTENT(IN)::  zwet_tmp(kproma)  ! dummy field to store a preliminary value of zwet
!--------

  ! Local variables
  REAL(dp) :: zwslim
  REAL(dp) :: zepdu2
  REAL(dp) :: zdu2(kproma)
  LOGICAL  :: lo,lo2
  REAL(dp) :: zepevap

  ! Local wiso variables
  REAL(dp) :: zdeltasn(kproma,mwiso) !isotopic ratio (delta-value)in snow reservoir
  REAL(dp) :: zdeltawl(kproma,mwiso) !isotopic ratio (delta-value)in soil water reservoir
  REAL(dp) :: zdeltaws(kproma,mwiso) !isotopic ratio (delta-value)in skin water reservoir
  REAL(dp) :: zwisowet(kproma,mwiso)
  REAL(dp) :: zwisosec

  ! Wiso constants
  REAL(dp) :: zwisosurfmin, zwisoroundoff

  ! Variables from vdiff
    ! calculated in vdiff
  REAL(dp), INTENT(IN) :: zchl(kproma) 
  REAL(dp), INTENT(IN) :: zhum(kproma)
  REAL(dp), INTENT(IN) :: zqsl(kproma) ! surf. sat. spec. humidity
    !assigned to vdiff
  REAL(dp), INTENT(IN) :: psn(kproma)  ! snow depth [m water equivalent] at the ground
  REAL(dp), INTENT(IN) :: pws(kproma)  ! Soil water content [m]
  REAL(dp), INTENT(IN) :: pwsmx(kproma)! Water holding capacity [m] of soil
  REAL(dp), INTENT(IN) :: pcvs(kproma) ! fractional snow cover (defined in *physc*)
  REAL(dp), INTENT(IN) :: pcvw(kproma) ! wet skin fraction
  REAL(dp), INTENT(IN) :: pvgrat(kproma)! vegetation ratio

  REAL(dp), INTENT(IN) :: pwl(kproma)  ! Water content [m] in skin reservoir
                          ! (vegetation and bare soil)
                          ! this one isnt really assigned, just in wiso

  ! Wiso variables coming from other processes (here:surf)
  REAL(dp), INTENT(IN) :: pwisosn(kproma,mwiso) ! Iso. snow depth [m water equivalent] at the ground
  REAL(dp), INTENT(IN) :: pwisowl(kproma,mwiso) ! Water content [m] in skin reservoir
                                   ! (vegetation and bare soil)
  REAL(dp), INTENT(IN) :: pwisows(kproma,mwiso) ! Soil water content [m]

  ! Final result of block
  REAL(dp),INTENT(OUT) :: zwisocsat(kproma,mwiso) ! iso. vapour
  REAL(dp), INTENT(OUT) :: zwisocair(kproma,mwiso) ! coefficients 

  INTEGER :: jt, jl
  
  ! Security Paramters for calculation of delta values
  zwisosurfmin=1.0e-12_dp
  zwisoroundoff=1.0e-8_dp  

  ! Computational constants from vdiff
  zwslim=0.35_dp
  zepevap=1.0e-10_dp

  ! precalculating these values
  zepdu2=1.0_dp
  DO jl=1,kproma
     zdu2(jl)=MAX(zepdu2,pum1(jl,klev)**2+pvm1(jl,klev)**2)
!!$     lo_dew(jl)=pqm1(jl,klev).GT.zqsl(jl) !---wiso-code: remember if dew occurs
  ENDDO

!   equivalent evapotranspiration efficiency coeff. over land - water isotopes
    DO jt=1,mwiso
       DO jl=1,kproma

!         calculate isotope ratio of snow reservoir
          zdeltasn(jl,jt)=1.0_dp
          IF ((pwisosn(jl,jt).gt.zwisosurfmin).and.(psn(jl).gt.zwisosurfmin)) zdeltasn(jl,jt)=(pwisosn(jl,jt)/psn(jl))/tnat(jt)
          IF (abs(1.0_dp-zdeltasn(jl,jt)).lt.zwisoroundoff) zdeltasn(jl,jt)=1.0_dp    ! cut off rounding errors
IF(jt.eq.1)zdeltasn(jl,jt)=1.0_dp!rehc
!         calculate isotope ratio of soil reservoir !re no, thats skin
          zdeltawl(jl,jt)=1.0_dp
          IF ((pwisowl(jl,jt).gt.zwisosurfmin).and.(pwl(jl).gt.zwisosurfmin)) zdeltawl(jl,jt)=(pwisowl(jl,jt)/pwl(jl))/tnat(jt)
          IF (abs(1.0_dp-zdeltawl(jl,jt)).lt.zwisoroundoff) zdeltawl(jl,jt)=1.0_dp    ! cut off rounding errors
IF(jt.eq.1)zdeltawl(jl,jt)=1.0_dp!rehc
!         calculate isotope ratio of skin reservoir !re no, thats soil
          zdeltaws(jl,jt)=1.0_dp
          IF ((pwisows(jl,jt).gt.zwisosurfmin).and.(pws(jl).gt.zwisosurfmin)) zdeltaws(jl,jt)=(pwisows(jl,jt)/pws(jl))/tnat(jt)
          IF (abs(1.0_dp-zdeltaws(jl,jt)).lt.zwisoroundoff) zdeltaws(jl,jt)=1.0_dp    ! cut off rounding errors
IF(jt.eq.1)zdeltaws(jl,jt)=1.0_dp!rehc

! calculation of tracer coefficients *zwisocair* and *zwisocsat* analog to
! normal vapour coefficients *zcair* and *zcsat*
! *zwisocsat* and *zwisocair* are used to calculate surface flux j:
! j = constant * ( zwisocair * tracer(nlev) - zwisocsat * tracer_saturation )
! (see also: g. hoffmann, dissertation, p.29ff)

          IF (pws(jl) .GT. zwslim*pwsmx(jl)) THEN
              zwisowet(jl,jt)= pcvs(jl)*zdeltasn(jl,jt)                        &
                      + (1.0_dp-pcvs(jl))*(pcvw(jl)*zdeltawl(jl,jt)      &
                      + (1.0_dp-pcvw(jl))*zdeltaws(jl,jt)                & 
                      / (1.0_dp+zchl(jl)*SQRT(zdu2(jl))*zwet_tmp(jl)))
          ELSE
             zwisowet(jl,jt) = pcvs(jl)*zdeltasn(jl,jt)                       &
                      + (1.0_dp-pcvs(jl))* pcvw(jl)*zdeltawl(jl,jt) 
          ENDIF

          lo=zhum(jl).LE.pqm1(jl,klev)/zqsl(jl) .OR. zhum(jl).LT.zepevap

          zwisocsat(jl,jt)=pcvs(jl)*zdeltasn(jl,jt)+(1.0_dp-pcvs(jl))                     &
                         *(pcvw(jl)*zdeltawl(jl,jt)+(1.0_dp-pcvw(jl))*zdeltaws(jl,jt)     &
                         *MERGE(0.0_dp,zhum(jl),lo))
          zwisocair(jl,jt)=pcvs(jl)*zdeltasn(jl,jt)+(1.0_dp-pcvs(jl))                     &
                         *(pcvw(jl)*zdeltawl(jl,jt)+(1.0_dp-pcvw(jl))*zdeltaws(jl,jt)     &
                         *MERGE(0.0_dp,1.0_dp,lo))

          zwisosec=MERGE(pcvs(jl)*zdeltasn(jl,jt)+((1.0_dp-pcvs(jl))*(pcvw(jl)*zdeltawl(jl,jt)  &
                  + (1.0_dp-pcvw(jl))*zdeltaws(jl,jt))),1.0_dp,lpland(jl))

          lo2=pqm1(jl,klev).GT.zqsl(jl)

          zwisocsat(jl,jt)=MERGE(zwisosec,zwisocsat(jl,jt),lo2)
          zwisocair(jl,jt)=MERGE(zwisosec,zwisocair(jl,jt),lo2)

          zwisocsat(jl,jt)=pvgrat(jl)*zwisowet(jl,jt)+(1.0_dp-pvgrat(jl))*zwisocsat(jl,jt)
          zwisocair(jl,jt)=pvgrat(jl)*zwisowet(jl,jt)+(1.0_dp-pvgrat(jl))*zwisocair(jl,jt)

          zwisocsat(jl,jt)=zwisocsat(jl,jt)*tnat(jt)
          zwisocair(jl,jt)=zwisocair(jl,jt)*tnat(jt)

        END DO
      END DO


END SUBROUTINE h2oiso_vdiff_vapcoeff

   ! =================================

   ! =================================

   ! =================================
SUBROUTINE h2oiso_vdiff_kinfac( kproma, klev     &
                              , pum1    , pvm1   &
                              , zwisokin         &
                              , lpland           )

  IMPLICIT NONE

  ! local
  REAL(dp):: zwspeed     ! wind speed for kin. fraction
  REAL(dp):: zwspeedmin  ! minimum wind speed for kin. fraction

  ! usable
  INTEGER,  INTENT(IN) :: kproma, klev
  REAL(dp), INTENT(IN) :: pum1(kproma,klev)
  REAL(dp), INTENT(IN) :: pvm1(kproma,klev)
  LOGICAL, INTENT(IN)  :: lpland(kproma)

  ! Result of block
  REAL(dp), INTENT(OUT) :: zwisokin(kproma,mwiso)  ! kinetic fractionation factor at air/sea interface

  INTEGER :: jt, jl
  
! Minimum Windspeed for calculation of kinetic fractionation effects
  zwspeedmin=7.0_dp

!     calculation of the kinetic fract. factor *zwisokin*
      DO jt=1,mwiso
         DO jl=1,kproma
!           absolute value of windspeed
            zwspeed=sqrt(pum1(jl,klev)**2+pvm1(jl,klev)**2)
!           kin. fractionation according to formula of brutsaert, 1975
!           (see also: g. hoffmann, dissertation, p.21)
            IF (zwspeed.le.zwspeedmin) THEN
               zwisokin(jl,jt)=tkinsl(jt)
            ELSE
               zwisokin(jl,jt)=tkinfa1(jt)*zwspeed+tkinfa2(jt)
            ENDIF
!           no kin. fractionation over land surface
            zwisokin(jl,jt)=MERGE(1.0_dp,1.0_dp-zwisokin(jl,jt),lpland(jl))
         END DO
      END DO

END SUBROUTINE h2oiso_vdiff_kinfac


   ! =================================

   ! =================================

SUBROUTINE h2oiso_vdiff_RMscheme( kproma, klev, klevp1, klevm1 &
                                , paphm1, cvdifts              &
                                , zcfhl, zcfhw, zcfhi          &
                                , zcfh, zebsh, zwisokin        &
                                , zwisoeqnl, zwisofqnl         &
                                , zwisoeqnw, zwisofqnw         &
                                , zwisoeqni, zwisofqni         &
                                , pwisoqm1                     &
                                , zwisoqdif )


  IMPLICIT NONE

  ! usable
  INTEGER, INTENT(IN) :: kproma, klev, klevp1, klevm1
  REAL(dp), INTENT(IN) :: cvdifts !->cvdifts->mo_physc2

  ! GETS PRESET HERE
  !REAL(dp) :: zwisoqdif(kproma,klev,mwiso)

  REAL(dp), INTENT(OUT) :: zwisoqdif(kproma,klev,mwiso)

  ! variable
  REAL(dp),INTENT(IN):: pwisoqm1(kproma,klev,mwiso), zwisokin(kproma,mwiso)

  ! local variable
  REAL(dp):: zfac, zdisc
  REAL(dp):: zwisodisql(kproma,mwiso)
  REAL(dp):: zwisodisqw
  REAL(dp):: zwisodisqi
  REAL(dp):: ztpfac1, ztpfac2!->cvdifts->mo_physc2
  INTEGER:: itop, itopp1, jk

  ! variables from vdiff
    ! calculated in vdiff
  REAL(dp):: zqdp
  REAL(dp), INTENT(IN):: zcfh(kproma,klev) ! heat transfer coeff.
  REAL(dp), INTENT(IN):: zebsh(kproma,klev)
  REAL(dp), INTENT(IN):: zcfhl(kproma)
  REAL(dp), INTENT(IN):: zcfhw(kproma)
  REAL(dp), INTENT(IN):: zcfhi(kproma)
    ! assigned to vdiff
  REAL(dp), INTENT(IN):: paphm1(kproma,klevp1) ! half level pressure (t-dt)

  ! Result of block
  REAL(dp), INTENT(OUT) :: zwisoeqnl(kproma,mwiso),zwisofqnl(kproma,mwiso) ! RM-coeffs
  REAL(dp), INTENT(OUT) :: zwisoeqnw(kproma,mwiso),zwisofqnw(kproma,mwiso) ! for moisture
  REAL(dp), INTENT(OUT) :: zwisoeqni(kproma,mwiso),zwisofqni(kproma,mwiso) ! land,water,ice

  INTEGER :: jt, jl
  
  itop = 1
  itopp1 = itop+1
  ztpfac1=cvdifts
  ztpfac2=1.0_dp/ztpfac1
!----------------------------------------------------------------
  DO jt=1,mwiso
     DO jk=1,klev
        DO jl=1,kproma
           zwisoqdif(jl,jk,jt)=0.0_dp
        END DO
     END DO
  END DO

  DO jt=1,mwiso
     DO jk=itop,klev
        DO jl=1,kproma
           zwisoqdif(jl,jk,jt)=ztpfac2*pwisoqm1(jl,jk,jt)
        END DO
     ENDDO
  END DO

  DO jt=1,mwiso
     DO jl=1,kproma
        zqdp=1.0_dp/(paphm1(jl,itopp1)-paphm1(jl,itop))
        zdisc=1.0_dp/(1.0_dp+zcfh(jl,itop)*zqdp)
        zwisoqdif(jl,itop,jt)=zdisc*zwisoqdif(jl,itop,jt)
     ENDDO
  ENDDO

  DO jt=1,mwiso
     DO jk=itopp1,klevm1
        DO jl=1,kproma
           zqdp=1.0_dp/(paphm1(jl,jk+1)-paphm1(jl,jk))
           zfac=zcfh(jl,jk-1)*zqdp
           zdisc=1.0_dp/(1.0_dp+zfac*(1.0_dp-zebsh(jl,jk-1))              &
                +zcfh(jl,jk)*zqdp)
           zwisoqdif(jl,jk,jt) = zdisc*(zwisoqdif(jl,jk,jt)+zfac*zwisoqdif(jl,jk-1,jt))
        END DO
     END DO
  END DO
!----------------------------------------------------------------

     DO jt=1,mwiso
       DO jl=1,kproma
        zqdp=1.0_dp/(paphm1(jl,klevp1)-paphm1(jl,klev))
        zfac=zcfh(jl,klevm1)*zqdp
!
! calculation of the EN and FN coefficients of the Richtmyer-Morton-Scheme - water isotopes
!
!       land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zwisodisql=1.0_dp/(1.0_dp+zfac*(1.0_dp-zebsh(jl,klevm1)))
        zwisoeqnl(jl,jt)=zwisodisql(jl,jt)*zcfhl(jl)*zqdp
        zwisofqnl(jl,jt)=zwisodisql(jl,jt)*(zwisoqdif(jl,klev,jt)+zfac*zwisoqdif(jl,klevm1,jt))*ztpfac1
!
!       water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zwisodisqw=1.0_dp/(1.0_dp+zfac*(1.0_dp-zebsh(jl,klevm1))              &
                  +zwisokin(jl,jt)*zcfhw(jl)*zqdp)
        zwisoeqnw(jl,jt)=zwisodisqw*zwisokin(jl,jt)*zcfhw(jl)*zqdp
        zwisofqnw(jl,jt)=zwisodisqw*(zwisoqdif(jl,klev,jt)+zfac*zwisoqdif(jl,klevm1,jt))*ztpfac1
!

!       ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zwisodisqi=1.0_dp/(1.0_dp+zfac*(1.0_dp-zebsh(jl,klevm1))              &
                  +zcfhi(jl)*zqdp)
        zwisoeqni(jl,jt)=zwisodisqi*zcfhi(jl)*zqdp
        zwisofqni(jl,jt)=zwisodisqi*(zwisoqdif(jl,klev,jt)+zfac*zwisoqdif(jl,klevm1,jt))*ztpfac1

        END DO
     END DO

END SUBROUTINE h2oiso_vdiff_RMscheme

   ! =================================

   ! =================================

SUBROUTINE h2oiso_vdiff_zqdif(kproma, klev, klevm1, ztmst      &
                             , zcfhl, zcfhw, zcfhi             &
                             , pfrl , pfrw  , pfri             &
                             , cvdifts, zqdif_jl               &
                             , zqslnew, zwisoqdif              &
                             , zwisoeqnl, zwisofqnl, zwisoeqnw &
                             , zwisofqnw, zwisoeqni, zwisofqni &
                             , zwisoqsw , zwisoqsi, zwisocsat  &
                             , zwisocair, zebsh                &
                             , zwisodqdt2, pwisoqm1 )

  IMPLICIT NONE


  ! usable
  INTEGER, INTENT(IN) :: kproma, klev, klevm1
  REAL(dp), INTENT(IN):: cvdifts, ztmst

  REAL(dp), INTENT(IN):: pwisoqm1(kproma,klev,mwiso)

  ! RESULT
  REAL(dp), INTENT(INOUT) :: zwisoqdif(kproma,klev,mwiso)
  REAL(dp), INTENT(OUT):: zwisodqdt2(kproma,klev,mwiso)

  ! from vdiff
  ! assigned to vdiff
  REAL(dp), INTENT(IN):: pfrl(kproma) ! land cover -> slm
  REAL(dp), INTENT(IN):: pfrw(kproma) ! sea cover  -> seacov
  REAL(dp), INTENT(IN):: pfri(kproma) ! ice cover  -> icecov
  ! calculated in vdiff
  REAL(dp), INTENT(IN):: zcfhl(kproma)
  REAL(dp), INTENT(IN):: zcfhw(kproma)
  REAL(dp), INTENT(IN):: zcfhi(kproma)
  REAL(dp), INTENT(IN):: zqdif_jl(kproma,klev)
  REAL(dp), INTENT(IN):: zqslnew(kproma)
  REAL(dp), INTENT(IN):: zebsh(kproma,klev)

  ! local
  INTEGER:: jk, itop
  REAL(dp):: zwisoqklevi(kproma,mwiso), zwisoqklevw(kproma,mwiso), zwisoqklevl(kproma,mwiso)
  REAL(dp):: zgqsum,zgtl,zgtw,zgti,zgqw,zgqi,zgql
  REAL(dp):: ztpfac1, ztpfac2, ztpfac3
  REAL(dp):: zwisoqslnew(kproma,mwiso)
!  REAL(dp), INTENT(IN):: zwisoqsl(kproma,mwiso)
!  REAL(dp):: zwisodqdt
  REAL(dp):: zcons13

  ! from RMscheme block
  REAL(dp), INTENT(IN):: zwisoeqnl(kproma,mwiso)
  REAL(dp), INTENT(IN):: zwisofqnl(kproma,mwiso)
  REAL(dp), INTENT(IN):: zwisoeqnw(kproma,mwiso)
  REAL(dp), INTENT(IN):: zwisofqnw(kproma,mwiso)
  REAL(dp), INTENT(IN):: zwisoeqni(kproma,mwiso)
  REAL(dp), INTENT(IN):: zwisofqni(kproma,mwiso)
  ! from surfhum block
  REAL(dp), INTENT(IN):: zwisoqsw(kproma,mwiso)
  REAL(dp), INTENT(IN):: zwisoqsi(kproma,mwiso)
  REAL(dp), INTENT(IN):: zwisocsat(kproma,mwiso)
  REAL(dp), INTENT(IN):: zwisocair(kproma,mwiso)

  INTEGER :: jt, jl
  
  ztpfac1=cvdifts
  ztpfac2=1.0_dp/ztpfac1
  ztpfac3=1.0_dp-ztpfac2
  itop = 1
  zcons13=1.0_dp/ztmst

!   DO jl = 1, kproma
!      DO jk = 1, klev
!         DO jt = 1, mwiso
!            zwisoqdif(jl,jk,jt)=0.0_dp
!         ENDDO
!      ENDDO
!   ENDDO

  DO jt=1,mwiso
     DO jl=1,kproma
        zwisoqslnew(jl,jt)=zqslnew(jl)  ! no fractionation for land surface processes
     END DO
  END DO

!* Calculation of qklev using the new surface value ZQSLNEW - water isotopes
!* (assuming no fractionation)
!
     DO jt=1,mwiso
      DO jl=1,kproma
!      land   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zwisoqklevl(jl,jt)=zwisoeqnl(jl,jt)                                      &
                           *(zwisocsat(jl,jt)*zwisoqslnew(jl,jt)                 &
                            -zwisocair(jl,jt)*zqdif_jl(jl,klev)*ztpfac1)            &
                           +zwisofqnl(jl,jt)
!
!
!      water   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zwisoqklevw(jl,jt)=zwisoeqnw(jl,jt)*zwisoqsw(jl,jt)+zwisofqnw(jl,jt)
!
!      ice   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
        zwisoqklevi(jl,jt)=zwisoeqni(jl,jt)*zwisoqsi(jl,jt)+zwisofqni(jl,jt)
!
!    Grid-mean specific humidity at the
!    'blending height' (here: lowest model level 'klev') - water isotopes
!
        zgtl=pfrl(jl)*zcfhl(jl)
        zgtw=pfrw(jl)*zcfhw(jl)
        zgti=pfri(jl)*zcfhi(jl)
        IF (pfrl(jl).LT.1.0_dp) THEN
           zgql=zgtl*zwisocair(jl,jt)
        ELSE
           zgql=zgtl
        ENDIF
        zgqw=zgtw
        zgqi=zgti
        zgqsum=(zgql+zgqw+zgqi)/ztpfac2
        zwisoqdif(jl,klev,jt)=(zgql*zwisoqklevl(jl,jt)+zgqw*zwisoqklevw(jl,jt)+zgqi*zwisoqklevi(jl,jt))/zgqsum

       END DO
     END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! BACK-SUBSTITUTION.
     DO jt=1,mwiso
        DO jk=klevm1,itop,-1
           DO jl=1,kproma
              zwisoqdif(jl,jk,jt)=zwisoqdif(jl,jk,jt)+zebsh(jl,jk)*zwisoqdif(jl,jk+1,jt)
           ENDDO
        ENDDO
     ENDDO

     ! INCREMENTATION OF Q TENDENCIES.
     DO jt=1,mwiso
        DO jk=itop,klev
           DO jl=1,kproma
              zwisoqdif(jl,jk,jt)=zwisoqdif(jl,jk,jt)+ztpfac3*pwisoqm1(jl,jk,jt)

              ! Dimenson change in result because of outsourcing 
              ! of the tendency computation
              !zwisodqdt=(zwisoqdif(jl,jk,jt)-pwisoqm1(jl,jk,jt))*zcons13
              zwisodqdt2(jl,jk,jt)=(zwisoqdif(jl,jk,jt)-pwisoqm1(jl,jk,jt))*zcons13
              ! The actual tendency calculation takes place in the SMIL!!
              !pwisoqte(jl,jk,jt)=pwisoqte(jl,jk,jt)+zwisodqdt

           ENDDO
        ENDDO
      ENDDO

!       
END SUBROUTINE h2oiso_vdiff_zqdif
   ! =================================

   ! =================================
SUBROUTINE h2oiso_vdiff_mflux( kproma, klev, g                 &
                             , zcfhl, zcfhw, zcfhi             &
                             , cvdifts, pqm1                   &
                             , ztmst, zqdif                    &
                             , zqslnew, lpland                 &
                             , pwisoevapl, pwisoevapot         &
                             , zwisoqdif, zwisocsat, zwisocair &
                             , pwisoqm1, delta_time            &
                             , pfrl, pfrw, pfri                &
                             , zwisoqsw, zwisoqsi              &
                             , pwisoevap, pwisoevapw           &
                             , pwisoevapi, pwisoevaplac        &
                             , pwisoevapwac, pwisoevapiac )
    
  IMPLICIT NONE


  ! usable
  INTEGER, INTENT(IN) :: klev, kproma
  REAL(dp), INTENT(IN):: cvdifts, delta_time
  REAL(dp), INTENT(IN):: pqm1(kproma,klev)
  !variable
  REAL(DP), INTENT(IN):: pwisoqm1(kproma,klev,mwiso)
  LOGICAL, INTENT(IN) :: lpland(kproma)
  REAL(DP), INTENT(IN):: g
  REAL(DP), INTENT(IN):: ztmst

  REAL(dp), INTENT(IN):: pfrl(kproma), pfri(kproma), pfrw(kproma)

  ! results of block
   REAL(dp), INTENT(INOUT):: pwisoevap(kproma,mwiso)
   REAL(dp), INTENT(INOUT):: pwisoevapl(kproma,mwiso)
   REAL(dp), INTENT(INOUT):: pwisoevapw(kproma,mwiso)
   REAL(dp), INTENT(INOUT):: pwisoevapi(kproma,mwiso)
   REAL(dp), INTENT(INOUT):: pwisoevaplac(kproma,mwiso)
   REAL(dp), INTENT(INOUT):: pwisoevapwac(kproma,mwiso)
   REAL(dp), INTENT(INOUT):: pwisoevapiac(kproma,mwiso)

  REAL(dp), INTENT(OUT):: pwisoevapot(kproma,mwiso)

  ! local
  REAL(dp):: zcons15, zdtime
  REAL(dp):: ztpfac1, ztpfac2, ztpfac3
  REAL(dp):: zwisozqs
  REAL(dp):: zwisoqnlev
  REAL(dp):: zqnlev
  REAL(dp):: zcoefl, zcoefw, zcoefi
  REAL(dp):: zwisoqhfll(kproma,mwiso), zwisoqhflw(kproma,mwiso), zwisoqhfli(kproma,mwiso)
  REAL(dp):: zwisoqslnew(kproma,mwiso)

  ! variable from vdiff
    ! calculated in vdiff
  REAL(dp), INTENT(IN):: zqdif(kproma,klev)
  REAL(dp), INTENT(IN):: zcfhl(kproma), zcfhw(kproma), zcfhi(kproma)
  REAL(dp), INTENT(IN):: zqslnew(kproma)
  REAL(dp), INTENT(IN):: zwisoqsw(kproma,mwiso), zwisoqsi(kproma,mwiso)

  ! from zqdif block
  REAL(dp), INTENT(IN) :: zwisoqdif(kproma,klev,mwiso)
  ! from surfhum block
  REAL(dp), INTENT(IN) :: zwisocair(kproma,mwiso)
  REAL(dp), INTENT(IN) :: zwisocsat(kproma,mwiso)

  INTEGER :: jt, jl
  
  DO jt=1,mwiso
     DO jl=1,kproma
        zwisoqslnew(jl,jt)=zqslnew(jl)  ! no fractionation for land surface processes
     END DO
  END DO

  zdtime = delta_time

  ztpfac1=cvdifts
  ztpfac2=1.0_dp/ztpfac1
  ztpfac3=1.0_dp-ztpfac2
  zcons15=1.0_dp/(g*ztmst)

  !    Surface fluxes of moisture - water isotopes

  DO jt=1,mwiso
     DO jl=1,kproma
       zcoefl=zcons15*zcfhl(jl)
       zcoefw=zcons15*zcfhw(jl)
       zcoefi=zcons15*zcfhi(jl)

!    Moisture fluxes - water isotopes
        zwisoqnlev=zwisoqdif(jl,klev,jt)-ztpfac3*pwisoqm1(jl,klev,jt)
        zqnlev=zqdif(jl,klev)-ztpfac3*pqm1(jl,klev)
        zwisoqnlev=MERGE(zqnlev,zwisoqnlev,lpland(jl))
!
!       land   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        zwisozqs=ztpfac2*zwisoqslnew(jl,jt)
        zwisoqhfll(jl,jt)=zcoefl*(zwisocair(jl,jt)*zwisoqnlev-zwisocsat(jl,jt)*zwisozqs)
!
!       water   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        zwisoqhflw(jl,jt)=zcoefw*(zwisoqnlev-ztpfac2*zwisoqsw(jl,jt))
!       ice   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        zwisoqhfli(jl,jt)=zcoefi*(zwisoqnlev-ztpfac2*zwisoqsi(jl,jt))
!
!    Accumulated evaporation - water isotopes
        pwisoevapl(jl,jt)=zwisoqhfll(jl,jt)
        pwisoevapw(jl,jt)=zwisoqhflw(jl,jt)
        pwisoevapi(jl,jt)=zwisoqhfli(jl,jt)

        pwisoevap(jl,jt)=pwisoevap(jl,jt)+(pfrl(jl)*zwisoqhfll(jl,jt)              &
                                          +pfrw(jl)*zwisoqhflw(jl,jt)              &
                                          +pfri(jl)*zwisoqhfli(jl,jt))*zdtime
!
!    Potential evaporation over land (for snow and skin reservoir) - water isotopes
          pwisoevapot(jl,jt)=zcoefl*(zwisoqnlev-zwisozqs)

!     Accumulated area weighted evaporations for land/water/ice - water isotopes
       pwisoevaplac(jl,jt)=pwisoevaplac(jl,jt)+pfrl(jl)*pwisoevapl(jl,jt)*zdtime
       pwisoevapwac(jl,jt)=pwisoevapwac(jl,jt)+pfrw(jl)*pwisoevapw(jl,jt)*zdtime
       pwisoevapiac(jl,jt)=pwisoevapiac(jl,jt)+pfri(jl)*pwisoevapi(jl,jt)*zdtime
     END DO
  END DO

END SUBROUTINE h2oiso_vdiff_mflux
   ! =================================



   ! =================================
   ! =================================
   ! =================================

! All the ISO CONVECT computation is in external files
! (messy_h2oiso_convect.f90 and messy_h2oiso_tiedtke_convect.f90)

   ! =================================
   ! =================================
   ! =================================

! All the ISO CLOUD computation is in an external file
! due to clarity of the file (messy_h2oiso_cloud.f90)

   ! =================================
   ! =================================
   ! =================================



   ! =================================
SUBROUTINE h2oiso_surf_convflux( kproma, delta_time, psn             &
                               , pcvs, pcvw, pwisoevapot             &
                               , pwl, pwisoevapl                     &
                               , zwisoevsnd, zwisoraind              &
                               , zwisosnowd, zwisoevwld, zwisoevwsd  &
                               , pwisorsfl, pwisorsfc                &
                               , pwisossfl, pwisossfc                &
                               , pwisosn, pwisowl                    &
                               , zwisoevttd            )

  IMPLICIT NONE


  ! local
  INTEGER, INTENT(IN) :: kproma
  REAL(dp), INTENT(IN) :: delta_time
  REAL(dp) :: zdtime
  REAL(dp) :: zdelta
  REAL(dp) :: zwisomin
  REAL(dp) :: zwisosec

  ! assigned to surf
  REAL(dp), INTENT(IN) :: pcvs(kproma)
  REAL(dp), INTENT(IN) :: pcvw(kproma)

  REAL(dp), INTENT(OUT) :: zwisoevttd(kproma,mwiso)

  ! Results of block
  REAL(dp), INTENT(OUT):: zwisoevsnd(kproma,mwiso)
  REAL(dp), INTENT(OUT):: zwisoevwld(kproma,mwiso)
  REAL(dp), INTENT(OUT):: zwisoraind(kproma,mwiso)
  REAL(dp), INTENT(OUT):: zwisosnowd(kproma,mwiso)
  REAL(dp), INTENT(OUT):: zwisoevwsd(kproma,mwiso)

  ! Coming from vdiff
  REAL(dp), INTENT(IN):: pwisoevapot(kproma,mwiso)
  REAL(dp), INTENT(IN):: pwisoevapl(kproma,mwiso)

  ! Coming from cloud
  REAL(dp), INTENT(IN):: pwisorsfl(kproma,mwiso)
  REAL(dp), INTENT(IN):: pwisorsfc(kproma,mwiso)
  REAL(dp), INTENT(IN):: pwisossfl(kproma,mwiso)
  REAL(dp), INTENT(IN):: pwisossfc(kproma,mwiso)

  ! Coming from last TS
  REAL(dp), INTENT(IN):: pwisosn(kproma,mwiso)
  REAL(dp), INTENT(IN):: pwisowl(kproma,mwiso)

  ! assigned to surf
  REAL(dp), INTENT(IN) :: psn(kproma)
  REAL(dp), INTENT(IN) :: pwl(kproma)

  INTEGER :: jt, jl
  
  zwisosec= 1.0e-10_dp
  zwisomin = 1.0e-12_dp
  zdtime = delta_time

! Convert water fluxes to [m water equivalent * timestep] - water isotopes
!
 DO jt=1,mwiso
  DO jl=1,kproma
     zwisoraind(jl,jt)=(pwisorsfl(jl,jt)+pwisorsfc(jl,jt))*zdtime/twisorhoh2o(jt)
     zwisosnowd(jl,jt)=(pwisossfl(jl,jt)+pwisossfc(jl,jt))*zdtime/twisorhoh2o(jt)

!    calculate isotope ratio of snow reservoir
     zdelta=tnat(jt)
     IF (psn(jl).gt.zwisomin .and. pwisosn(jl,jt).gt.zwisomin) zdelta = MIN(pwisosn(jl,jt)/psn(jl),1.0_dp)!remin
     IF (abs(1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp  ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc

     zwisoevsnd(jl,jt)=pcvs(jl)*zdelta*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)
     
!    calculate isotope ratio of skin reservoir
     zdelta=tnat(jt)
     IF (pwl(jl).gt.zwisomin .and. pwisowl(jl,jt).gt.zwisomin) zdelta = MIN(pwisowl(jl,jt)/pwl(jl),1.0_dp)!remin
     IF (abs(1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp  ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
     zwisoevwld(jl,jt)=(1.0_dp-pcvs(jl))*pcvw(jl)*zdelta*pwisoevapot(jl,jt)*zdtime/twisorhoh2o(jt)

     zwisoevttd(jl,jt)=pwisoevapl(jl,jt)*zdtime/twisorhoh2o(jt)
     zwisoevwsd(jl,jt)=zwisoevttd(jl,jt)-zwisoevsnd(jl,jt)
  END DO
 END DO

END SUBROUTINE h2oiso_surf_convflux
   ! =================================

   ! =================================
SUBROUTINE h2oiso_surf_snowc( kproma, klev, rhoh2o, lpglac        &
                            , cwlmax, cvinter, delta_time         &
                            , zwisosnowd, zwisoevsnd, zevsnd      &
                            , zwisosncmelt, psnc, zsnowd          &
                            , ptm1, tte, pssfl, pssfc             &
                            , lpland, pevapot                     &
                            , pcvs, pu10, pv10                    &
                            , vlt, psn, zsn_tmp1, zsn_tmp2        &
                            , pwisosnc                            &
                            , zsncmelt, pwlmx, time_step_len      &
                            , zwisosn        )


  IMPLICIT NONE
  
  ! local
  INTEGER, INTENT(IN) :: kproma, klev
  LOGICAL, INTENT(IN) :: lpland(kproma)
  REAL(dp), INTENT(IN) :: ptm1(kproma,klev)
  REAL(dp), INTENT(IN) :: delta_time, cvinter, time_step_len
  REAL(dp), INTENT(INOUT) :: psnc(kproma)
  REAL(dp), INTENT(IN) :: pssfl(kproma) 
  REAL(dp), INTENT(IN) :: pssfc(kproma)
  REAL(dp), INTENT(IN) :: rhoh2o, cwlmax
  REAL(dp):: zmprcp1_tmp(kproma)
  REAL(dp), INTENT(IN) :: pevapot(kproma)
  REAL(dp), INTENT(IN) :: pcvs(kproma)
  REAL(dp), INTENT(IN) :: pu10(kproma), pv10(kproma)
  REAL(dp), INTENT(IN) :: vlt(kproma), tte(kproma,klev)
  REAL(dp), INTENT(INOUT) :: psn(kproma)
  LOGICAL, INTENT(IN)  :: lpglac(kproma)
  ! results
  REAL(dp), INTENT(OUT):: zwisosncmelt(kproma,mwiso)
  ! Coming from convflux block
  REAL(dp), INTENT(INOUT):: zwisosnowd(kproma,mwiso)
  REAL(dp), INTENT(INOUT):: zwisoevsnd(kproma,mwiso)

  ! local
  REAL(dp), INTENT(INOUT):: zwisosn(kproma,mwiso)
  REAL(dp):: zdelta
  REAL(dp):: zwisomin, zwisosec
  REAL(dp):: zwisomprcp
  REAL(dp):: zwisosncp
  REAL(dp), INTENT(INOUT):: pwisosnc(kproma,mwiso)
  REAL(dp):: zwisosubl
  REAL(dp):: zc2, zexpt, zdtime, zc3, zexpw
  REAL(dp):: zwisosncwind, zsncmax
  REAL(dp):: zsnowd_tmp(kproma), zsncp_tmp(kproma)
  REAL(dp),INTENT(OUT):: zsnowd(kproma)
  REAL(dp):: zsubl_tmp(kproma), zsnc_tmp(kproma)
  REAL(dp):: zsnc_tmp2(kproma), zsncwind_tmp(kproma)
  LOGICAL :: lo_sub(kproma)
  REAL(dp):: zevsnd_tmp(kproma), lo_tm1(kproma,klev)
  REAL(DP),INTENT(OUT)::zevsnd(kproma)

  REAL(dp), INTENT(OUT):: zsncmelt(kproma)
  REAL(dp), INTENT(OUT):: zsn_tmp1(kproma), zsn_tmp2(kproma)
  REAL(dp), INTENT(OUT):: pwlmx(kproma)
  INTEGER :: jt, jl
  
  zwisomin= 1.0e-12_dp
  zwisosec= 1.0e-10_dp
  zc2=1.87E5_dp
  zc3=1.56E5_dp
  zdtime = delta_time

  lo_tm1(1:kproma,:) = ptm1(1:kproma,:) +  tte(1:kproma,:) * time_step_len

  DO jt=1,mwiso
     DO jl=1,kproma
        zwisosn(jl,jt)=0._dp
     ENDDO
  ENDDO

  ! -----------------
  ! This is to create some of the temporary values right here
  DO jl=1,kproma

     zsn_tmp1(jl)= 0.0_dp
     zsn_tmp2(jl)= 0.0_dp

     zsncmelt(jl)=0.0_dp
     IF (.NOT.lpglac(jl)) THEN
        pwlmx(jl)=cwlmax*(1.0_dp+vlt(jl))
     ELSE
        pwlmx(jl)=0.0_dp
     END IF

        zevsnd(jl) = pcvs(jl) * pevapot(jl) * zdtime/rhoh2o
        zsnowd(jl) = (pssfl(jl) + pssfc(jl)) * zdtime / rhoh2o
        zsnowd_tmp(jl) = zsnowd(jl)        !---wiso-code: remember total snow

     IF (lpland(jl).AND..NOT.lpglac(jl)) THEN

        zsncmax=MAX(0.0_dp,pwlmx(jl)-cwlmax)
        zmprcp1_tmp(jl)=MIN(zsnowd(jl)*cvinter,zsncmax-psnc(jl))
        zsncp_tmp(jl) = psnc(jl) + zmprcp1_tmp(jl)
        zsnowd(jl)=zsnowd(jl)-zmprcp1_tmp(jl)
        zevsnd_tmp(jl) = zevsnd(jl)
        lo_sub(jl) = (zevsnd_tmp(jl).GT.0.0_dp)  
        psnc(jl)=MIN(MAX(0.0_dp,zsncp_tmp(jl)+zevsnd_tmp(jl)),zsncmax)
        zsubl_tmp(jl) = psnc(jl) - zsncp_tmp(jl)
        zsnc_tmp(jl) = psnc(jl)
        !zexpt=MAX(0.0_dp,ptm1(jl,klev)+3.0_dp-tmelt)*zdtime/zc2
        zexpt=MAX(0.0_dp,lo_tm1(jl,klev)+3.0_dp-tmelt)*zdtime/zc2
        zsncmelt(jl)=psnc(jl)*(1.0_dp-EXP(-zexpt))
        zevsnd(jl)=zevsnd_tmp(jl)-(psnc(jl)-zsncp_tmp(jl))
        psnc(jl)=psnc(jl)-zsncmelt(jl)
        zsnc_tmp2(jl)=psnc(jl)
        zexpw = SQRT(pu10(jl)**2+pv10(jl)**2)*zdtime/zc3
        zsncwind_tmp(jl) = psnc(jl)*(1.0_dp-EXP(-zexpw))
        psnc(jl)=psnc(jl)-zsncwind_tmp(jl)

        !zsn_tmp1 is for snsulw, zsn_tmp2 for snglme -------

        zsnowd(jl)=zsnowd(jl)+zsncwind_tmp(jl)

        !zsn_tmp1(jl)=psn(jl)+zsnowd_tmp2(jl)+zevsnd_tmp(jl)
!        zsn_tmp2(jl)=psn(jl)+zsnowd_tmp2(jl)+zevsnd_tmp(jl)
!        IF (psn(jl).LT.0.0_dp) THEN
!           zsn_tmp2(jl)=0.0_dp
!        END IF

     ELSE
        psnc(jl)=0.0_dp

!        zsn_tmp2(jl)=0.0_dp
        !---------------------------------------------------------

     ENDIF

  ENDDO
  ! -----------------


      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            psn(jl)=psn(jl)+zsnowd(jl)+zevsnd(jl)
            zsn_tmp1(jl)=psn(jl) ! wiso-code: remember snow fall
            IF (psn(jl).LT.0.0_dp) THEN
               psn(jl)=0.0_dp
            END IF
         ELSE
            psn(jl)=0.0_dp
         END IF
         zsn_tmp2(jl) = psn(jl)
      END DO



   ! Snow changes in the canopy (interception of snowfall,
   ! sublimation, melting, unloading due to wind) - water isotopes (no additional fractionation)
   DO jt=1,mwiso
      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            zwisosn(jl,jt)=zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)

            zdelta=tnat(jt)
            IF (zsnowd_tmp(jl).gt.zwisomin .and. zwisosnowd(jl,jt).gt.zwisomin) zdelta=MIN(zwisosnowd(jl,jt)/zsnowd_tmp(jl),1.0_dp)!remin
            IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp   ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
            zwisomprcp=zmprcp1_tmp(jl) * zdelta
              zwisosncp=pwisosnc(jl,jt)+zwisomprcp
            zwisosnowd(jl,jt)=zwisosnowd(jl,jt)-zwisomprcp
            zdelta=tnat(jt)
            IF (lo_sub(jl)) THEN                           ! deposition, assume identical delta value as the atmospheric flux
              IF (zevsnd_tmp(jl).gt.zwisomin .and. zwisoevsnd(jl,jt).gt.zwisomin) &
                   zdelta=MIN(zwisoevsnd(jl,jt)/zevsnd_tmp(jl),1.0_dp)!remin
            ELSE                                           ! sublimation, assume identical delta value as the canopy layer
              IF (zsncp_tmp(jl).gt.zwisomin .and. zwisosncp.gt.zwisomin) &
                   zdelta=MIN(zwisosncp/zsncp_tmp(jl),1.0_dp)!remin
            ENDIF  
            IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp   ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc

            zwisosubl=zsubl_tmp(jl) * zdelta
            pwisosnc(jl,jt)=zwisosncp+zwisosubl
            zwisoevsnd(jl,jt)=zwisoevsnd(jl,jt)-zwisosubl
            
            zdelta=tnat(jt)
            IF (zsnc_tmp(jl).gt.zwisomin .and. pwisosnc(jl,jt).gt.zwisomin) zdelta=MIN(pwisosnc(jl,jt)/zsnc_tmp(jl),1.0_dp)!remin
            IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp   ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
            zwisosncmelt(jl,jt)=zsncmelt(jl) * zdelta
            pwisosnc(jl,jt)=pwisosnc(jl,jt)-zwisosncmelt(jl,jt)

            zdelta=tnat(jt)
            IF (zsnc_tmp2(jl).gt.zwisomin .and. pwisosnc(jl,jt).gt.zwisomin) zdelta=MIN(pwisosnc(jl,jt)/zsnc_tmp2(jl),1.0_dp)!remin
            IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp   ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
            zwisosncwind=zsncwind_tmp(jl) * zdelta
            pwisosnc(jl,jt)=pwisosnc(jl,jt)-zwisosncwind
            zwisosnowd(jl,jt)=zwisosnowd(jl,jt)+zwisosncwind
         ELSE
            pwisosnc(jl,jt)=0.0_dp
         END IF
      END DO
   END DO




END SUBROUTINE h2oiso_surf_snowc
   ! =================================


   ! =================================

SUBROUTINE h2oiso_surf_snsulw( kproma                   &
                             , pwisosn, zwisosnowd      &
                             , zwisoevsnd, zwisoevwsd   &
                             , pwisogld                 &
                             , lpland, lpglac, zsn_tmp1 &
                             , pwisoalac, zwisoevttd    &
                             , zwisorogl, zwisoraind    )

  IMPLICIT NONE

  INTEGER, INTENT(IN):: kproma
  LOGICAL, INTENT(IN):: lpland(kproma), lpglac(kproma)

  REAL(DP), INTENT(INOUT):: pwisoalac(kproma,mwiso)

  ! from previous and on
  REAL(DP), INTENT(INOUT):: pwisosn(kproma,mwiso)
  REAL(DP), INTENT(INOUT):: zwisoevwsd(kproma,mwiso)

  ! from snowc block
  REAL(DP), INTENT(IN):: zwisosnowd(kproma,mwiso)
  REAL(DP), INTENT(IN):: zwisoevsnd(kproma,mwiso)
  REAL(DP), INTENT(IN):: zwisoevttd(kproma,mwiso)

  REAL(DP), INTENT(IN):: zwisoraind(kproma,mwiso)

  ! for next block
  REAL(DP), INTENT(INOUT):: pwisogld(kproma,mwiso)

  ! calculated in snowc
  REAL(DP), INTENT(IN):: zsn_tmp1(kproma)
  REAL(DP), INTENT(INOUT):: zwisorogl(kproma,mwiso)
  INTEGER :: jt, jl

   ! Snowfall and sublimation on land (excluding glaciers) - water isotopes
   DO jt=1,mwiso
      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            pwisosn(jl,jt)=pwisosn(jl,jt)+zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)
            IF (zsn_tmp1(jl).LT.0.0_dp) THEN
               zwisoevwsd(jl,jt)=zwisoevwsd(jl,jt)+pwisosn(jl,jt)
               pwisosn(jl,jt)=0.0_dp
            END IF
         ELSE
            pwisosn(jl,jt)=0.0_dp
         END IF
      END DO
   END DO

   ! Snowfall and sublimation on glaciers and diagnostics - water isotopes
   DO jt=1,mwiso
      DO jl=1,kproma
         IF (lpglac(jl)) THEN
            pwisogld(jl,jt)=pwisogld(jl,jt)+zwisosnowd(jl,jt)+zwisoevsnd(jl,jt)
            pwisoalac(jl,jt)=zwisoraind(jl,jt)+zwisosnowd(jl,jt)+zwisoevttd(jl,jt)
            zwisorogl(jl,jt)=zwisoraind(jl,jt)
         END IF
      END DO
   END DO



END SUBROUTINE h2oiso_surf_snsulw

   ! =================================


   ! =================================
SUBROUTINE h2oiso_surf_snglme( kproma, lstart, rhoh2o          &
                             , pwisosn                         &
                             , pwisogld                        &
                             , lpland, lpglac                  &
                             , pgld, ptsl, psn                 &
                             , zwisosnmel, zsnmel, zsn_tmp2    &
                             , pgrndcapc, zsnowd, zevsnd       &
                             , zwisorogl)

  IMPLICIT NONE

  INTEGER, INTENT(IN):: kproma
  LOGICAL, INTENT(IN):: lpland(kproma), lpglac(kproma)
  LOGICAL, INTENT(IN):: lstart  
  REAL(dp), INTENT(INOUT) :: psn(kproma)
  REAL(dp), INTENT(IN) :: pgrndcapc(kproma), zsnowd(kproma)

  REAL(DP), INTENT(IN):: zsn_tmp2(kproma)
  REAL(DP), INTENT(IN):: ptsl(kproma)
  REAL(dp), INTENT(IN):: rhoh2o, zevsnd(kproma)
  REAL(dp), INTENT(INOUT) :: pgld(kproma)
  REAL(dp), INTENT(INOUT) :: zwisorogl(kproma,mwiso)

  ! local
  REAL(DP):: zdelta, zsmelt
  REAL(DP):: zwisomin, zwisosec
  REAL(DP):: ztsl_tmp(kproma),zgld_tmp(kproma)

  !
  REAL(DP),INTENT(OUT):: zsnmel(kproma)!re this is now a temp value!

  ! results
  REAL(DP), INTENT(OUT):: zwisosnmel(kproma,mwiso)
  ! from previous and on
  REAL(DP), INTENT(INOUT):: pwisosn(kproma,mwiso)
  REAL(DP), INTENT(INOUT):: pwisogld(kproma,mwiso)
  INTEGER :: jl, jt
  
  zwisomin= 1.0e-12_dp
  zwisosec= 1.0e-10_dp

  DO jl = 1, kproma
     zsnmel(jl) = 0.0_dp
     DO jt = 1, mwiso
        zwisosnmel(jl,jt) = 0.0_dp
     ENDDO
     IF (lpglac(jl)) THEN
        pgld(jl)=pgld(jl)+zsnowd(jl)+zevsnd(jl)
     ENDIF
  END DO

  IF (.NOT. lstart) THEN
     DO jl=1,kproma
      zgld_tmp(jl) = pgld(jl)
      ztsl_tmp(jl) = ptsl(jl)
        IF (lpland(jl).AND.ptsl(jl).GT.tmelt) THEN
           IF (lpglac(jl)) THEN
              zsnmel(jl)=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
              pgld(jl)=pgld(jl)-zsnmel(jl)
           ELSE IF (psn(jl).GT.0.0_dp) THEN
              zsmelt=pgrndcapc(jl)*(ptsl(jl)-tmelt)/(alf*rhoh2o)
              zsnmel(jl)=MIN(psn(jl),zsmelt)
              psn(jl)=MAX(psn(jl)-zsnmel(jl),0.0_dp)
           END IF
        END IF
     END DO
   END IF

  ! Snow and glacier melt - water isotopes (no additional fractionation)
   IF (.NOT. lstart) THEN
     DO jt=1,mwiso
       DO jl=1,kproma
         IF (lpland(jl).AND.ztsl_tmp(jl).GT.tmelt) THEN
            IF (lpglac(jl)) THEN
               zdelta=tnat(jt)
               IF (zgld_tmp(jl).gt.zwisomin .and. pwisogld(jl,jt).gt.zwisomin) zdelta=MIN(pwisogld(jl,jt)/zgld_tmp(jl),1.0_dp)!remin
               IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
                 zwisosnmel(jl,jt)=zsnmel(jl) * zdelta
               pwisogld(jl,jt)=pwisogld(jl,jt)-zwisosnmel(jl,jt)
               zwisorogl(jl,jt)=zwisorogl(jl,jt)+zwisosnmel(jl,jt)
            ELSE IF (zsn_tmp2(jl).GT.0.0_dp) THEN
               zdelta=tnat(jt)
               IF (zsn_tmp2(jl).gt.zwisomin .and. pwisosn(jl,jt).gt.zwisomin) zdelta=MIN(pwisosn(jl,jt)/zsn_tmp2(jl),1.0_dp)!remin
               IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
                 zwisosnmel(jl,jt)=zsnmel(jl) * zdelta
               pwisosn(jl,jt)=MAX(pwisosn(jl,jt)-zwisosnmel(jl,jt),0.0_dp)
            END IF
         END IF
       END DO
     END DO
   END IF





END SUBROUTINE h2oiso_surf_snglme
   ! =================================


   ! =================================
SUBROUTINE h2oiso_surf_snbmew( kproma                   &
                             , pwisowl, zwisomlres      &
                             , zwisosnmel, zwisosncmelt &
                             , lpland, lpglac, zsncmelt &
                             , pwl, pwlmx               &
                             , zsnmlt, zsnmel           &
                             , zwisosn )

  IMPLICIT NONE

  INTEGER, INTENT(IN):: kproma
  LOGICAL, INTENT(IN):: lpland(kproma), lpglac(kproma)
  REAL(DP),INTENT(IN):: zsncmelt(kproma), pwlmx(kproma)
  REAL(dp),INTENT(INOUT):: pwl(kproma)

  ! local
  REAL(DP):: zdelta
  REAL(DP):: zwisomin, zwisosec
  REAL(DP):: zwisosnmlt

  ! results of block
  REAL(DP),INTENT(INOUT):: pwisowl(kproma,mwiso)
  REAL(DP),INTENT(OUT):: zwisomlres(kproma,mwiso)
  REAL(DP),INTENT(INOUT):: zwisosnmel(kproma,mwiso)
  REAL(DP),INTENT(INOUT):: zsnmlt(kproma), zsnmel(kproma)

  ! from snowc
  REAL(DP),INTENT(IN):: zwisosncmelt(kproma,mwiso)
  REAL(DP), INTENT(INOUT):: zwisosn(kproma,mwiso)

  ! local temps
  REAL(DP):: zwl_tmp(kproma)
  REAL(DP):: zsnmlt_tmp(kproma)!beware theres another one in cloud
  INTEGER :: jl, jt

  zwisomin= 1.0e-12_dp
  zwisosec= 1.0e-10_dp

  DO jl=1,kproma
     IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
        pwl(jl)=pwl(jl)+zsncmelt(jl)
        zwl_tmp(jl)=pwl(jl)                                                !---wiso-code: remember skin reservoir
        zsnmlt_tmp(jl)=MAX(0.0_dp,pwl(jl)-pwlmx(jl))                        !---wiso-code: remember additional meltwater
        zsnmlt(jl)=zsnmel(jl)+MAX(0.0_dp,pwl(jl)-pwlmx(jl))
        pwl(jl)=MIN(pwlmx(jl),pwl(jl))
        zsnmel(jl)=zsnmel(jl)+zsncmelt(jl)
     END IF
  END DO


   ! Snow budget and meltwater (glacier-free land only) - water isotopes (no additional fractionation)
   DO jt=1,mwiso
      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
            pwisowl(jl,jt)=pwisowl(jl,jt)+zwisosncmelt(jl,jt)
            zdelta=tnat(jt)
            IF (zwl_tmp(jl).gt.zwisomin .and. pwisowl(jl,jt).gt.zwisomin) zdelta=MIN(pwisowl(jl,jt)/zwl_tmp(jl),1.0_dp)!remin
            IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
            zwisosnmlt=zwisosnmel(jl,jt) + (zsnmlt_tmp(jl) * zdelta)
            pwisowl(jl,jt)=pwisowl(jl,jt) - (zsnmlt_tmp(jl) * zdelta)
            zwisomlres(jl,jt)=zwisosnmlt
            zwisosnmel(jl,jt)=zwisosnmel(jl,jt)+zwisosncmelt(jl,jt)
            zwisosn(jl,jt)=zwisosn(jl,jt)-zwisosnmel(jl,jt)
         END IF
      END DO
   END DO


END SUBROUTINE h2oiso_surf_snbmew
   ! =================================


   ! ==================================
SUBROUTINE h2oiso_surf_skinres( kproma, delta_time, rhoh2o &
                              , zwisoraind, pwisowl                &
                              , zwisoevwld, zwisoevwsd             &
                              , lpland, lpglac                     &
                              , pcvs, pcvw, pevapot                &
                              , pwlmx_tmp                          &
                              , pwl, cvinter, zraind_tmp           &
                              , zmprcp2_tmp                        )


  IMPLICIT NONE

  INTEGER, INTENT(IN):: kproma
  LOGICAL, INTENT(IN):: lpland(kproma), lpglac(kproma)
  REAL(dp), INTENT(IN):: pcvs(kproma), pcvw(kproma)
  REAL(dp), INTENT(IN) :: rhoh2o, delta_time, cvinter
  REAL(dp), INTENT(IN) :: pevapot(kproma)
  REAL(dp), INTENT(IN) :: pwlmx_tmp(kproma)
  REAL(dp), INTENT(INOUT) :: pwl(kproma)

  ! results of block
  REAL(DP), INTENT(INOUT):: zwisoraind(kproma,mwiso)
  REAL(DP), INTENT(INOUT):: pwisowl(kproma,mwiso)
  REAL(DP), INTENT(INOUT):: zwisoevwsd(kproma,mwiso)

  ! from previous block
  REAL(DP), INTENT(IN):: zwisoevwld(kproma,mwiso)

  ! local
  REAL(DP):: zdelta, zwisoevap
  REAL(DP):: zwisomin, zwisosec
  REAL(DP):: zwisomprcp, zwisowlp
  REAL(DP):: zdtime, zevwld_tmp(kproma)
  REAL(DP), INTENT(IN):: zraind_tmp(kproma)
  REAL(DP), INTENT(OUT):: zmprcp2_tmp(kproma)
  REAL(DP):: zwlp_tmp(kproma), zevap_tmp(kproma)
  LOGICAL:: lo_evap(kproma)
  INTEGER :: jl, jt
  
  !------------------
  zdtime = delta_time
  DO jl=1,kproma
     zevwld_tmp(jl) = (1.0_dp-pcvs(jl))*pcvw(jl)*pevapot(jl)*zdtime/rhoh2o

     IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
        zmprcp2_tmp(jl) = MIN(zraind_tmp(jl)*cvinter,pwlmx_tmp(jl)-pwl(jl))
        zwlp_tmp(jl) = pwl(jl)+zmprcp2_tmp(jl)
        lo_evap(jl) = (zevwld_tmp(jl).GT.0.0_dp)
        pwl(jl)=MIN(MAX(0.0_dp,zwlp_tmp(jl)+zevwld_tmp(jl)),pwlmx_tmp(jl))
        zevap_tmp(jl)=pwl(jl)-zwlp_tmp(jl)
     ELSE
        pwl(jl)=0.0_dp
     ENDIF

  ENDDO

  zwisomin= 1.0e-12_dp
  zwisosec= 1.0e-10_dp
  !------------------

   ! Water budget - water isotopes
   ! Skin reservoir (vegetation and bare soil) - water isotopes
   DO jt=1,mwiso
      DO jl=1,kproma
         IF (lpland(jl).AND..NOT.lpglac(jl)) THEN

            ! Interception of rain - water isotopes

            zdelta=tnat(jt)
            IF (zraind_tmp(jl).gt.zwisomin .and. zwisoraind(jl,jt).gt.zwisomin) zdelta=MIN(zwisoraind(jl,jt)/zraind_tmp(jl),1.0_dp)!remin
            IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp   ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
            zwisomprcp=zmprcp2_tmp(jl) * zdelta
            zwisowlp=pwisowl(jl,jt)+zwisomprcp
            zwisoraind(jl,jt)=zwisoraind(jl,jt)-zwisomprcp

            ! Evaporation or dew collection - water isotopes

            zdelta=tnat(jt)
            IF (lo_evap(jl)) THEN                          ! dew collection, assume identical delta value as the atmospheric flux
              IF (zevwld_tmp(jl).gt.zwisomin .and. zwisoevwld(jl,jt).gt.zwisomin) &
                   zdelta=MIN(zwisoevwld(jl,jt)/zevwld_tmp(jl),1.0_dp)!remin
            ELSE                                           ! evaporation, assume identical delta value as the skin reservoir
              IF (zwlp_tmp(jl).gt.zwisomin .and. zwisowlp.gt.zwisomin) &
                   zdelta=MIN(zwisowlp/zwlp_tmp(jl),1.0_dp)!remin
            ENDIF  
            IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp   ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc
            zwisoevap=zevap_tmp(jl) * zdelta
            pwisowl(jl,jt)=zwisowlp+zwisoevap
              zwisoevwsd(jl,jt)=zwisoevwsd(jl,jt)-zwisoevap
          ELSE
            pwisowl(jl,jt)=0.0_dp
          END IF
      END DO
   END DO


END SUBROUTINE h2oiso_surf_skinres
   ! =================================


   ! ==================================
SUBROUTINE h2oiso_surf_soilres( kproma, jpgrnd                 &
                              , zwisomlres, zwisoraind         &
                              , zwisoevwsd, pwisows            &
                              , zevwsd, zsnmlt, pwsmx          &
                              , lpland, lpglac, zraind_tmp     &
                              , delta_time, ptsoil             &
                              , zmprcp2_tmp, ngl, porostd, pws &
                              , zwisoros, zwisodrain )

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: kproma, jpgrnd
  LOGICAL, INTENT(IN) :: lpland(kproma), lpglac(kproma)
  REAL(DP), INTENT(IN):: pwsmx(kproma)
  REAL(DP), INTENT(IN):: delta_time
  REAL(DP), INTENT(IN):: ptsoil(kproma,jpgrnd)
  REAL(DP):: zdrain(kproma)
  REAL(DP), INTENT(IN):: zevwsd(kproma)

  ! local
  REAL(DP):: zwslim, zwisoinfil, zdelta
  REAL(DP):: zwisoprfl, zwisomin, zwisosec
  REAL(DP), INTENT(OUT):: zwisoros(kproma,mwiso) 
  REAL(DP), INTENT(OUT):: zwisodrain(kproma,mwiso)
  REAL(DP):: zwisowsup

  ! result of block
  REAL(DP), INTENT(INOUT):: pwisows(kproma,mwiso)

  ! from previous block
  REAL(DP), INTENT(IN):: zwisomlres(kproma,mwiso)
  REAL(DP), INTENT(IN):: zwisoraind(kproma,mwiso)
  REAL(DP), INTENT(IN):: zwisoevwsd(kproma,mwiso)
  REAL(DP), INTENT(IN):: zsnmlt(kproma)

  REAL(DP),INTENT(IN):: zraind_tmp(kproma), zmprcp2_tmp(kproma)
  REAL(DP):: zraind2_tmp(kproma)
  REAL(DP):: zmlres(kproma)
  REAL(DP):: zprfl_tmp(kproma)
  REAL(DP):: zros_tmp(kproma)
  LOGICAL:: lo_zros(kproma)
  REAL(DP):: zws_tmp(kproma), zwsup_tmp(kproma)
  REAL(DP):: zdtime
  !new
  REAL(DP):: zvol, zlysic, zros(kproma), zconw1, zbm, zlyeps
  REAL(DP), INTENT(INOUT):: pws(kproma)
  REAL(DP):: zb1, zbws, zroeff, zorvari, zorvars
  INTEGER, INTENT(IN):: ngl
  REAL(DP), INTENT(IN):: porostd(kproma)
  REAL(DP):: zinfil, zwsup, zdrmin, zwdtr, zconw2, zconw3, zdrexp, zdrmax
  INTEGER :: jt, jl
  
  zdtime = delta_time
  zdrmin=0.001_dp/(3600.0_dp*1000.0_dp)
  zwisomin= 1.0e-12_dp
  zwisosec= 1.0e-10_dp
  zorvari=100.0_dp
  zorvars=1000.0_dp*64.0_dp/ngl
  zdrmax=0.1_dp/(3600.0_dp*1000.0_dp)
  zdrexp=1.5_dp

  DO jt=1,mwiso
     DO jl=1,kproma
        zwisoros(jl,jt)   = 0.0_dp
        zwisodrain(jl,jt) = 0.0_dp
     ENDDO
  ENDDO

  DO jl=1,kproma
     IF (lpland(jl).AND..NOT.lpglac(jl)) THEN
        
        zwdtr=0.90_dp*pwsmx(jl)
        zmlres(jl)=zsnmlt(jl)
        zconw2=pwsmx(jl)-zwdtr
        zconw3=zdrmax-zdrmin



        zraind2_tmp(jl) = zraind_tmp(jl)-zmprcp2_tmp(jl)

        IF (zevwsd(jl) >= 0.0_dp) THEN
           zprfl_tmp(jl) = zmlres(jl)+zraind2_tmp(jl)+zevwsd(jl)
        ELSE
           pws(jl)=pws(jl)+zevwsd(jl)
           zprfl_tmp(jl) = zmlres(jl)+zraind2_tmp(jl)
        END IF

        zroeff=MAX(0.0_dp, porostd(jl)-zorvari)/(porostd(jl)+zorvars)
        zbws=MAX(MIN(zroeff,0.5_dp),0.01_dp)
        zb1=1.0_dp+zbws
        zbm=1.0_dp/zb1
        zconw1=pwsmx(jl)*zb1
        zwslim=0.05_dp*pwsmx(jl)
        zvol = 0.0_dp
        zros(jl) = 0.0_dp
        zinfil=0.0_dp

        IF (ptsoil(jl,1).LT.tmelt) THEN
           zros(jl)=zprfl_tmp(jl)
        ELSE
           lo_zros(jl) = (zprfl_tmp(jl).GT.0.0_dp.AND.pws(jl).GT.zwslim) !---wiso-code: remember if runoff occurs
           IF (zprfl_tmp(jl).GT.0.0_dp.AND.pws(jl).GT.zwslim) THEN
              IF (pws(jl).GT.pwsmx(jl)) THEN
                 zlyeps=pws(jl)-pwsmx(jl)
              ELSE
                 zlyeps=0.0_dp
              END IF
              zlysic=(pws(jl)-zlyeps)/pwsmx(jl)
              zlysic=MIN(zlysic,1.0_dp)
              zvol=(1.0_dp-zlysic)**zbm-zprfl_tmp(jl)/zconw1
              zros(jl)=zprfl_tmp(jl)-(pwsmx(jl)-pws(jl))
              IF (zvol.GT.0.0_dp) THEN
                 zros(jl)=zros(jl)+pwsmx(jl)*zvol**zb1
              END IF
              zros(jl)=MAX(zros(jl),0.0_dp)
              zros_tmp(jl) = zros(jl)           !---wiso-code: remember runoff value
              zinfil=zprfl_tmp(jl)-zros(jl)
           ELSE
              zros(jl)=0.0_dp
              zinfil=zprfl_tmp(jl)
           END IF
           pws(jl)=pws(jl)+zinfil
        END IF
        
!*    4.2.2  Drainage and total runoff
!
             zws_tmp(jl)=pws(jl)                                               !---wiso-code: remember value of pws before drainage
             IF (pws(jl).LE.zwslim) THEN
                zdrain(jl)=0.0_dp
             ELSE
                IF (ptsoil(jl,1).GT.tmelt) THEN
                   zdrain(jl)=zdrmin*pws(jl)/pwsmx(jl)
                   IF (pws(jl).GT.zwdtr) THEN
                      zdrain(jl)=zdrain(jl)+zconw3*                  &
                                ((pws(jl)-zwdtr)/zconw2)**zdrexp
                   END IF
                   zdrain(jl)=zdrain(jl)*zdtime
                   zdrain(jl)=MIN(zdrain(jl),pws(jl)-zwslim)
                   pws(jl)=pws(jl)-zdrain(jl)
                ELSE
                   zdrain(jl)=0.0_dp
                END IF
             END IF
             zwsup=MAX(pws(jl)-pwsmx(jl),0.0_dp)
             zwsup_tmp(jl) = zwsup                                             !---wiso-code: remember soil water surplus
             pws(jl)=pws(jl)-zwsup
             zros(jl)=zros(jl)+zdrain(jl)+zwsup
          ELSE
             pws(jl)=0.0_dp

     END IF
  END DO

!------------------------------------
!---until here is variables fetching
!------------------------------------

   ! Soil reservoir - water isotopes
   DO jt=1,mwiso
      DO jl=1,kproma
          IF (lpland(jl).AND..NOT.lpglac(jl)) THEN

             zwslim=0.05_dp*pwsmx(jl)
             zwisoinfil=0.0_dp

             ! Surface runoff, infiltration and evaporation from soil - water isotopes
             IF (zevwsd(jl) >= 0.0_dp) THEN
               zwisoprfl=zwisomlres(jl,jt)+zwisoraind(jl,jt)+zwisoevwsd(jl,jt)
             ELSE
               pwisows(jl,jt)=pwisows(jl,jt)+zwisoevwsd(jl,jt)
               zwisoprfl=zwisomlres(jl,jt)+zwisoraind(jl,jt)
             END IF

             ! assume that surface runoff has same delta value as precipitation (no mixing with soil water pool ws)
             zdelta=tnat(jt)
             IF (zprfl_tmp(jl).gt.zwisomin .and. zwisoprfl.gt.zwisomin) zdelta=MIN(zwisoprfl/zprfl_tmp(jl),1.0_dp)!remin
             IF ((1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc

             IF (ptsoil(jl,1).LT.tmelt) THEN
                zwisoros(jl,jt)=zwisoprfl
             ELSE
                IF (lo_zros(jl)) THEN
                   zwisoros(jl,jt)=zros_tmp(jl)*zdelta
                   zwisoinfil=zwisoprfl-zwisoros(jl,jt)
                ELSE
                   zwisoros(jl,jt)=0.0_dp
                   zwisoinfil=zwisoprfl
                END IF
                pwisows(jl,jt)=pwisows(jl,jt)+zwisoinfil
             END IF

             ! Drainage and total runoff - water isotopes

             ! drainage water and soil water surplus have same delta value as soil water ws
             zdelta=tnat(jt)
             IF (zws_tmp(jl).gt.zwisomin .and. pwisows(jl,jt).gt.zwisomin) zdelta=MIN(pwisows(jl,jt)/zws_tmp(jl),1.0_dp)!remin
             IF (ABS(1.0_dp-zdelta).lt.zwisosec) zdelta=1.0_dp ! cut off rounding errors
IF(jt.eq.1)zdelta=1.0_dp!rehc

             IF (zws_tmp(jl).LE.zwslim) THEN
                zwisodrain(jl,jt)=0.0_dp
             ELSE
                IF (ptsoil(jl,1).GT.tmelt) THEN
                   zwisodrain(jl,jt)=zdrain(jl)*zdelta
                   pwisows(jl,jt)=pwisows(jl,jt)-zwisodrain(jl,jt)
                ELSE
                   zwisodrain(jl,jt)=0.0_dp
                END IF
             END IF
             zwisowsup=zwsup_tmp(jl)*zdelta
             pwisows(jl,jt)=pwisows(jl,jt)-zwisowsup
             zwisoros(jl,jt)=zwisoros(jl,jt)+zwisodrain(jl,jt)+zwisowsup
          ELSE
             pwisows(jl,jt)=0.0_dp
          END IF

      END DO
   END DO



END SUBROUTINE h2oiso_surf_soilres

   ! =======================================

   ! ========================================
SUBROUTINE h2oiso_surf_corrmp( kproma                 &
                             , pwl, pws, psn          &
                             , psnc, pgld             &
                             , pwisowl                &
                             , pwisows                &
                             , pwisosn                &
                             , pwisosnc               &
                             , pwisogld               &
                             , pslm, zwisoros         &
                             , zwisosnmel, pwisoalac  &
                             , zwisodrain, zwisosn    &
                             , zwisorogl, pwisorunoff &
                             , pwisosnmel, pwisoapmegl&
                             , pwisodrain, pwisosnacl &
                             , pwisorogl)

  IMPLICIT NONE

  ! local 
  INTEGER,INTENT(IN):: kproma
  REAL(DP):: zwisomin

  REAL(DP),INTENT(IN):: pslm(kproma)
  REAL(DP),INTENT(IN):: zwisoros(kproma,mwiso), zwisosnmel(kproma,mwiso)

  ! results
  REAL(DP),INTENT(INOUT)::pwisowl(kproma,mwiso), pwisogld(kproma,mwiso)
  REAL(DP),INTENT(INOUT)::pwisows(kproma,mwiso), pwisosnc(kproma,mwiso)
  REAL(DP),INTENT(INOUT)::pwisosn(kproma,mwiso)

  REAL(dp), INTENT(IN) :: pwisoalac(kproma,mwiso), zwisodrain(kproma,mwiso)
  REAL(dp), INTENT(IN) :: zwisosn(kproma,mwiso), zwisorogl(kproma,mwiso)

  REAL(dp),INTENT(OUT):: pwisorunoff(kproma,mwiso), pwisorogl(kproma,mwiso)
  REAL(dp),INTENT(OUT):: pwisosnmel(kproma,mwiso), pwisoapmegl(kproma,mwiso)
  REAL(dp),INTENT(OUT):: pwisodrain(kproma,mwiso), pwisosnacl (kproma,mwiso)

  ! assigned to surf
  REAL(dp), INTENT(IN) :: pwl(kproma), pgld(kproma)
  REAL(dp), INTENT(IN) :: pws(kproma)
  REAL(dp), INTENT(IN) :: psn(kproma)
  REAL(dp), INTENT(IN) :: psnc(kproma)

  INTEGER :: jt, jl
  
  zwisomin= 1.0e-12_dp

   ! Water fluxes - water isotopes
   DO jt=1,mwiso
      DO jl=1,kproma
         pwisorunoff(jl,jt)= pwisorunoff(jl,jt) +zwisoros(jl,jt)   *twisorhoh2o(jt)*pslm(jl)
         pwisosnmel(jl,jt) = pwisosnmel(jl,jt)  +zwisosnmel(jl,jt) *twisorhoh2o(jt)*pslm(jl)
         pwisoapmegl(jl,jt)= pwisoapmegl(jl,jt) +pwisoalac(jl,jt)  *twisorhoh2o(jt)*pslm(jl)
         pwisodrain(jl,jt) = pwisodrain(jl,jt)  +zwisodrain(jl,jt) *twisorhoh2o(jt)*pslm(jl)
         pwisosnacl(jl,jt) = pwisosnacl(jl,jt)  +zwisosn(jl,jt)    *twisorhoh2o(jt)*pslm(jl)
         pwisorogl(jl,jt)  = pwisorogl(jl,jt)   +zwisorogl(jl,jt)  *twisorhoh2o(jt)*pslm(jl)
      END DO
   END DO


  !Corrections of minor water pools
   DO jt = 1,mwiso
     DO jl = 1,kproma

       IF (pwl(jl).LT.zwisomin) pwisowl(jl,jt)=pwl(jl)*tnat(jt)
       IF (pws(jl).LT.zwisomin) pwisows(jl,jt)=pws(jl)*tnat(jt)
       IF (psn(jl).LT.zwisomin) pwisosn(jl,jt)=psn(jl)*tnat(jt)
       IF (psnc(jl).LT.zwisomin) pwisosnc(jl,jt)=psnc(jl)*tnat(jt)
       IF (pgld(jl).LT.zwisomin) pwisogld(jl,jt)=pgld(jl)*tnat(jt)

       IF (pwisowl(jl,jt).LT.0.0_dp)  pwisowl(jl,jt)=pwl(jl)*tnat(jt)!retest
       IF (pwisows(jl,jt).LT.0.0_dp)  pwisows(jl,jt)=pws(jl)*tnat(jt)!retest
       IF (pwisosn(jl,jt).LT.0.0_dp)  pwisosn(jl,jt)=psn(jl)*tnat(jt)!retest
       IF (pwisosnc(jl,jt).LT.0.0_dp) pwisosnc(jl,jt)=psnc(jl)*tnat(jt)!retest
       IF (pwisogld(jl,jt).LT.0.0_dp) pwisogld(jl,jt)=pgld(jl)*tnat(jt)!retest


     END DO
   END DO


END SUBROUTINE h2oiso_surf_corrmp

   ! =======================================

  ! =========================================================================

  SUBROUTINE fractcal(kproma,kbdim,kwiso,pwisofrac,pwisofracice,pt,pmelt)
  ! ---------------------------------------------------
  !
  ! fractcal calculates fractionation coefficients
  ! for a band of longitudinal grid points, simultaneously 
  !
  ! ---------------------------------------------------
  
    IMPLICIT NONE

  ! input arguments
    INTEGER, INTENT(IN)     :: kproma, kbdim, kwiso
    REAL(dp), INTENT(IN)    :: pt(kbdim), pmelt
  
  ! input/output arguments
    REAL(dp), INTENT(INOUT) :: pwisofrac(kbdim,kwiso), pwisofracice(kbdim,kwiso)
  
  ! local variables
    INTEGER     :: jl,jt
    REAL(dp)    :: zsatval
    
  ! fractionation over water
    DO jt=1,kwiso
    DO jl=1,kproma
      pwisofrac(jl,jt)=exp(talph1(jt)/(pt(jl)**2)+talph2(jt)/pt(jl)+talph3(jt))
    END DO
    END DO
  
  ! fractionation over ice
    DO jt=1,kwiso
  
    IF (jt == i_HHO) THEN
      DO jl=1,kproma
      pwisofracice(jl,jt)=1.0_dp
      END DO
    ELSEIF (jt == i_HH18O) THEN
      DO jl=1,kproma
      pwisofracice(jl,jt)=exp(talps1(jt)/pt(jl)+talps2(jt))
      END DO
    ELSEIF (jt == i_HDO) THEN
      DO jl=1,kproma
      pwisofracice(jl,jt)=exp(talps1(jt)/(pt(jl)**2)+talps2(jt))
      END DO
    ENDIF
  
  ! effective fractionation over ice if necessary
    IF (jt /= 1) THEN
      DO jl=1,kproma
      IF (pt(jl).lt.pmelt) THEN
        zsatval=tsatbase-tsatfac*(pt(jl)-tmelt)
        pwisofracice(jl,jt)=pwisofracice(jl,jt)*(zsatval/(1.0_dp+pwisofracice(jl,jt)*(zsatval-1.0_dp)*tdifrel(jt)))
      ENDIF
      END DO
    ENDIF
  
    END DO
   
    RETURN
  END SUBROUTINE fractcal

  ! =========================================================================
  ! ### own private subroutines
  ! =========================================================================


  FUNCTION wiso_calc_delta(jt,pwiso,pdefault)
  
  IMPLICIT NONE
  
  REAL(dp), INTENT(IN) :: pwiso, pdefault
  INTEGER,  INTENT(IN) :: jt
  
  REAL(dp) :: wiso_calc_delta
  
  IF (ABS(pdefault).gt.cwisomin) THEN
     wiso_calc_delta = ((pwiso/pdefault)/tnat(jt)-1.0_dp)*1000.0_dp
  ELSE
     wiso_calc_delta = -9999.9_dp
  ENDIF
    
  END FUNCTION wiso_calc_delta
  
  ! =========================================================================

!--------------------------------------------------------------
!--------------------------------------------------------------
!--------------------------------------------------------------



SUBROUTINE physcadjwiso (kproma,kbdim,klev,kwiso,time_step_len, &
                         pqm1,     pxlm1,     pxim1,        &
                         pqte,     pxlte,     pxite,        &
                         pwisoqm1, pwisoxlm1, pwisoxim1,    &
                         pwisoqte, pwisoxlte, pwisoxite)
!
!      M. WERNER           AWI, BREMERHAVEN      2009
!
!      PURPOSE
!      -------
!      ADJUSTMENT OF WATER ISOTOPE TRACERS TO DEFAULT MODEL WATER VARIABLES Q, XL AND XI
!      (ATTENTION: IT IS ASSUMEND, THAT ISOTOPE TRACER #1 IS SET TO H2-16O)
!
!      INTERFACE
!      ---------
!      THIS SUBROUTINE IS CALLED FROM
!        *PHYSC*
!
!      INPUT  NORMAL AND ISOTOPE VALUES AND TENDENCIES OF
!             VAPOUR, CLOUD WATER AND CLOUD ICE, 
!      OUTPUT NEW ISOTOPE VALUES AND TENDENCIES 
!
!      EXTERNALS
!      ---------
!      NONE


  IMPLICIT NONE

! input arguments  
  INTEGER, INTENT(IN) :: kproma,kbdim,klev,kwiso

  REAL(dp), INTENT(IN) :: time_step_len
  REAL(dp), INTENT(IN) :: pqm1(kbdim,klev),                       &
                          pxlm1(kbdim,klev),                      &
                          pxim1(kbdim,klev),                      &
                          pqte(kbdim,klev),                       &
                          pxlte(kbdim,klev),                      &
                          pxite(kbdim,klev)

! input/output arguments  
  REAL(dp), INTENT(INOUT) :: pwisoqm1(kbdim,klev,kwiso),                       &
                             pwisoxlm1(kbdim,klev,kwiso),                      &
                             pwisoxim1(kbdim,klev,kwiso),                      &
                             pwisoqte(kbdim,klev,kwiso),                       &
                             pwisoxlte(kbdim,klev,kwiso),                      &
                             pwisoxite(kbdim,klev,kwiso)
  
! local variables
  REAL(dp) :: zcorr(kbdim)
  
  LOGICAL  :: lcorr(kbdim)

  REAL(dp) :: zqpone,     zxlpone,     zxipone,           &
              zwisoqpone, zwisoxlpone, zwisoxipone
              
  REAL(dp) :: ztwodt

  INTEGER :: jl,jk,jt
    

  ztwodt=time_step_len

!
! time step m1
  DO jk=1,klev
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    jt=1
    DO jl = 1,kproma
      IF (pqm1(jl,jk) .GT.cwisomin .AND. pwisoqm1(jl,jk,jt) .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =pqm1(jl,jk) /pwisoqm1(jl,jk,1)
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        IF (lcorr(jl)) THEN 
          pwisoqm1 (jl,jk,jt)=zcorr(jl) *pwisoqm1 (jl,jk,jt)
        ELSE
          pwisoqm1(jl,jk,jt) = pqm1(jl,jk) * tnat(jt)  
        ENDIF
      ENDDO
    ENDDO
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    jt=1
    DO jl = 1,kproma
      IF (pxlm1(jl,jk) .GT.cwisomin .AND. pwisoxlm1(jl,jk,jt) .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =pxlm1(jl,jk) /pwisoxlm1(jl,jk,1)
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        IF (lcorr(jl)) THEN 
          pwisoxlm1 (jl,jk,jt)=zcorr(jl) *pwisoxlm1 (jl,jk,jt)
        ELSE
          pwisoxlm1(jl,jk,jt) = pxlm1(jl,jk) * tnat(jt)  
        ENDIF
      ENDDO
    ENDDO
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    jt=1
    DO jl = 1,kproma
      IF (pxim1(jl,jk) .GT.cwisomin .AND. pwisoxim1(jl,jk,jt) .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =pxim1(jl,jk) /pwisoxim1(jl,jk,1)
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        IF (lcorr(jl)) THEN 
          pwisoxim1 (jl,jk,jt)=zcorr(jl) *pwisoxim1 (jl,jk,jt)
        ELSE
          pwisoxim1(jl,jk,jt) = pxim1(jl,jk) * tnat(jt)  
        ENDIF
      ENDDO
    ENDDO    
  ENDDO

! time step p1
  DO jk=1,klev
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    DO jl = 1,kproma
      jt=1
      zqpone=pqm1(jl,jk)+pqte(jl,jk)*ztwodt
      zwisoqpone=pwisoqm1(jl,jk,jt)+pwisoqte(jl,jk,jt)*ztwodt
      IF (zqpone .GT.cwisomin .AND. zwisoqpone .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =zqpone /zwisoqpone
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        zqpone=pqm1(jl,jk)+pqte(jl,jk)*ztwodt
        zwisoqpone=pwisoqm1(jl,jk,jt)+pwisoqte(jl,jk,jt)*ztwodt
        IF (lcorr(jl)) THEN 
          zwisoqpone=zcorr(jl) *zwisoqpone
        ELSE
          zwisoqpone = zqpone * tnat(jt)  
        ENDIF
        pwisoqte(jl,jk,jt)=(zwisoqpone-pwisoqm1(jl,jk,jt))/ztwodt
      ENDDO
    ENDDO
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    DO jl = 1,kproma
      jt=1
      zxlpone=pxlm1(jl,jk)+pxlte(jl,jk)*ztwodt
      zwisoxlpone=pwisoxlm1(jl,jk,jt)+pwisoxlte(jl,jk,jt)*ztwodt
      IF (zxlpone .GT.cwisomin .AND. zwisoxlpone .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =zxlpone /zwisoxlpone
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        zxlpone=pxlm1(jl,jk)+pxlte(jl,jk)*ztwodt
        zwisoxlpone=pwisoxlm1(jl,jk,jt)+pwisoxlte(jl,jk,jt)*ztwodt
        IF (lcorr(jl)) THEN 
          zwisoxlpone=zcorr(jl) *zwisoxlpone
        ELSE
          zwisoxlpone = zxlpone * tnat(jt)  
        ENDIF
        pwisoxlte(jl,jk,jt)=(zwisoxlpone-pwisoxlm1(jl,jk,jt))/ztwodt
      ENDDO
    ENDDO
    zcorr(:) = 1._dp
    lcorr(:) = .FALSE.
    DO jl = 1,kproma
      jt=1
      zxipone=pxim1(jl,jk)+pxite(jl,jk)*ztwodt
      zwisoxipone=pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztwodt
      IF (zxipone .GT.cwisomin .AND. zwisoxipone .GT.cwisomin) THEN
        lcorr(jl)=.TRUE.
        zcorr(jl) =zxipone /zwisoxipone
      ENDIF
    ENDDO
    DO jt = 1,kwiso
      DO jl = 1,kproma
        zxipone=pxim1(jl,jk)+pxite(jl,jk)*ztwodt
        zwisoxipone=pwisoxim1(jl,jk,jt)+pwisoxite(jl,jk,jt)*ztwodt
        IF (lcorr(jl)) THEN 
          zwisoxipone=zcorr(jl) *zwisoxipone
        ELSE
          zwisoxipone = zxipone * tnat(jt)  
        ENDIF
        pwisoxite(jl,jk,jt)=(zwisoxipone-pwisoxim1(jl,jk,jt))/ztwodt
      ENDDO
    ENDDO
  ENDDO
  
  
!      Adjustment of water isotope tracers for negative default model water variables q, xl and xi
!      - for negative default water values the tracers are set to negative values, too, with a delta value equal to SMOW
!
! time step m1
  DO jt = 1, kwiso
    DO jk=1,klev
      DO jl = 1,kproma
        IF (pqm1(jl,jk).LT.0.0_dp) THEN
          pwisoqm1(jl,jk,jt)=tnat(jt)*pqm1(jl,jk)
        ENDIF
        IF (pxlm1(jl,jk).LT.0.0_dp) THEN
          pwisoxlm1(jl,jk,jt)=tnat(jt)*pxlm1(jl,jk)
        ENDIF
        IF (pxim1(jl,jk).LT.0.0_dp) THEN
          pwisoxim1(jl,jk,jt)=tnat(jt)*pxim1(jl,jk)
        ENDIF
      END DO
    END DO
  END DO
! time step p1
  DO jt = 1, kwiso
    DO jk=1,klev
      DO jl = 1,kproma
        zqpone=pqm1(jl,jk)+pqte(jl,jk)*ztwodt
        IF (zqpone.LT.0.0_dp) THEN
          zwisoqpone=tnat(jt)*zqpone
          pwisoqte(jl,jk,jt)=(zwisoqpone-pwisoqm1(jl,jk,jt))/ztwodt
        ENDIF
        zxlpone=pxlm1(jl,jk)+pxlte(jl,jk)*ztwodt
        IF (zxlpone.LT.0.0_dp) THEN
          zwisoxlpone=tnat(jt)*zxlpone
          pwisoxlte(jl,jk,jt)=(zwisoxlpone-pwisoxlm1(jl,jk,jt))/ztwodt
        ENDIF
        zxipone=pxim1(jl,jk)+pxite(jl,jk)*ztwodt
        IF (zxipone.LT.0.0_dp) THEN
          zwisoxipone=tnat(jt)*zxipone
          pwisoxite(jl,jk,jt)=(zwisoxipone-pwisoxim1(jl,jk,jt))/ztwodt
        ENDIF
      END DO
    END DO
  END DO
  

  RETURN
END SUBROUTINE physcadjwiso

!--------------------------------------------------------------
!--------------------------------------------------------------
!--------------------------------------------------------------

SUBROUTINE h2oiso_constdelta_tracer( X_rare, &
                                     R_s, &
                                     X_total, &
                                     M_total, &
                                     M_rare)
!
!      F. FRANK            DLR, OBERPFAFFENHOFEN 2017
!
!      PURPOSE
!      -------
!      CALCULATES THE RARE ISOTOPOLOGUE SO THAT THE SAMPLE ISOTOPIC RATIO IS
!      PRESERVED.
!
!      INTERFACE
!      ---------
!      THIS SUBROUTINE IS CALLED FROM
!        *GLOBAL_START*
!
!      INPUT  ISOTOPIC RATIO OF DESIRED SAMPLE (R_s), TOTAL MASTER ISOTOPOLOGUE (X_total)
!             MOLARMASS OF MASTER (M_total) AND MOLARMASS OF RARE ISOTOPOLOGUE (M_rare)
!      OUTPUT NEW ISOTOPOLOGUE VALUE (X_rare)
!
!      EXTERNALS
!      ---------
!      NONE
! 
! Note: This formula is only applicable in the h2oiso submodel where
!         1. H2OISOHHOvap is the total mixing ratio, hence HHO + HDO
!         2. the tracer H2OISOHDOvap = 0.5 HDO at all times, and the other 
!            isotopologues are with respect to Oxygen, which is included in H2O only once

  IMPLICIT NONE

! input arguments
  REAL(dp), INTENT(IN) :: R_s
  REAL(dp), INTENT(IN) :: M_total
  REAL(dp), INTENT(IN) :: M_rare
! array
  REAL(dp), INTENT(IN) :: X_total
! output arguments
  REAL(dp), INTENT(OUT) :: X_rare
  
  X_rare = ( R_s / (1.0_dp + R_s) ) * ( M_rare / M_total ) * X_total

END SUBROUTINE h2oiso_constdelta_tracer

!--------------------------------------------------------------

SUBROUTINE h2oiso_constdelta_tendency( X_rare_te, &
                                     R_s, &
                                     X_total, &
                                     X_rare, &
                                     M_total, &
                                     M_rare, &
                                     timestep)
!
!      F. FRANK            DLR, OBERPFAFFENHOFEN 2017
!
!      PURPOSE
!      -------
!      CALCULATES THE TENDENCY OF RARE ISOTOPOLOGUE NEEDED SO THAT THE SAMPLE RATIO IS
!      PRESERVED.
!
!      INTERFACE
!      ---------
!      THIS SUBROUTINE IS CALLED FROM
!        *PHYSC*
!
!      INPUT  RATIO OF DESIRED SAMPLE, TOTAL MASTER ISOTOPOLOGUE
!             MOLARMASS OF MASTER AND MOLARMASS OF RARE ISOTOPOLOGUE
!      OUTPUT NEW ISOTOPOLOGUE VALUE
!
!      EXTERNALS
!      ---------
!      NONE
!
! Note: This formula is only applicable in the h2oiso submodel where
!         1. H2OISOHHOvap is the total mixing ratio, hence HHO + HDO
!         2. the tracer H2OISOHDOvap = 0.5 HDO at all times, and the other 
!            isotopologues are with respect to Oxygen, which is included in H2O only once


  IMPLICIT NONE

! input arguments
  REAL(dp), INTENT(IN) :: R_s
  REAL(dp), INTENT(IN) :: M_total
  REAL(dp), INTENT(IN) :: M_rare
  REAL(dp), INTENT(IN) :: timestep
  
! array
  REAL(dp), INTENT(IN) :: X_total
  REAL(dp), INTENT(IN) :: X_rare
! output arguments
  REAL(dp), INTENT(OUT) :: X_rare_te
  

   X_rare_te = ( ( R_s / (1.0_dp + R_s) ) * ( M_rare / M_total ) * X_total - X_rare ) / timestep


END SUBROUTINE h2oiso_constdelta_tendency

! **********************************************************************
END MODULE messy_h2oiso
! **********************************************************************

