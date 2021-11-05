! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL MLOCEAN
!
! Authors:  Markus Kunze, FUB-IFM, January-February  2010 (MESSy Interface)
!           Erich Roeckner, MPI-M, 1995, (mixed layer ocean: ml_ocean.f90)
!
! References:
!
! * Roeckner, E., T. Siebert und J. Feichter, 
!   Climate response to anthropogenic sulfate forcing simulated with a general 
!   circulation model, 
!   In: Aerosol forcing of climate, Verlag JohnWiley & Sons, Chichester, 
!   New England, USA, 349-362, 1995.
!
! * P. Joeckel, R. Sander, A. Kerkweg, H. Tost, and J. Lelieveld,
!   Technical Note: The Modular Earth Submodel System (MESSy) - a new
!   approach towards Earth System Modeling,
!   Atmos. Chem. Phys., 5, 433-444, 2005.
!   http://www.atmos-chem-phys.net/5/433 
!
! **********************************************************************

! **********************************************************************
MODULE messy_mlocean_mlo
  ! **********************************************************************

  ! ----------- >
  
  USE messy_main_constants_mem, ONLY : dp  &
       , stbo       & ! Stephan-Boltzmann constant [W/m2/K4]
       , tmelt      & ! melting temperature of ice/snow
       , alf        & ! latent heat for fusion in J/kg (L_f)
       , rho_H2O    & ! density of H2O [kg/m3]
       , rho_sea    & ! density of sea water in kg m^-3
       , csw        & ! specific heat of sea waterJ/K/kg
       , ctfreez      ! temperature at which sea starts freezing/melting

  USE messy_mlocean, ONLY : mldmix &
                          , fbase_north &
                          , fbase_south
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: dp

  ! ----------- <
  !
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'mlocean_mlo' ! submodel name
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '0.2'         ! submodel version
  !
  REAL(dp), PARAMETER :: zalpha    = 2.1656_dp ! thermal conductivity of ice in J m^-1 K^-1
  REAL(dp), PARAMETER :: zalphas   = 0.31_dp   ! thermal conductivity of snow in J m^-1 K^-1
  REAL(dp), PARAMETER :: zrho_sn   = 330._dp   ! density of snow in kg m^-3
  REAL(dp), PARAMETER :: zrhoice   = 910._dp   ! density of ice in kg m^-3
  REAL(dp), PARAMETER :: zcpice    = 2106._dp  ! specific heat of ice Ws kg^-1 K^-1
  REAL(dp), PARAMETER :: zdice     = 0.10_dp   ! minimum thickness of ice slab in meter
  !
  ! PUBLIC SUBROUTINES (to be called from messy_mlocean_e5.f90)
  !
  PUBLIC :: mlocean_mlo
  PUBLIC :: mlocean_mlicetemp
  
CONTAINS
  !
  ! =========================================================================
  SUBROUTINE mlocean_mlo (kproma, lonorth, delta_time     &
                        , pslm,       palake,    pfsnet   &
                        , pevapi,     pcvsi,     pfri     & 
                        , pamlcorr,   pseaice,   psiced   &
                        , ptsi,       ptsw                &
                        , pfluxres,   pahfres,   psni     &
!!$                     , pamlcorac,  pamlheatac )          ! op_pj_20160617
                        , pamlheat )                        ! op_pj_20160617
    !
    !  ---------------------------------------------------------------------
    !
    ! Calculate parameters of a mixed layer ocean.
    ! (original file: ml_ocean.f90)
    !
    ! Author: Erich Roeckner, MPI-M, 1995
    !
    ! References:
    !
    ! * Roeckner, E., T. Siebert und J. Feichter, 
    !   Climate response to anthropogenic sulfate forcing simulated with a general 
    !   circulation model, 
    !   In: Aerosol forcing of climate, Verlag JohnWiley & Sons, Chichester, 
    !   New England, USA, 349-362, 1995.
    !
    ! Arguments INTENT(in)
    ! --------------------
    INTEGER,                INTENT(in) :: kproma
    LOGICAL,  DIMENSION(:), INTENT(in) :: lonorth    ! .TRUE. for northern latitude
    REAL(dp),               INTENT(in) :: delta_time 
    REAL(dp), DIMENSION(:), INTENT(in) :: pslm       ! sea-land mask
    REAL(dp), DIMENSION(:), INTENT(in) :: palake     ! lake fraction     (lake == 1)
    REAL(dp), DIMENSION(:), INTENT(in) :: pfsnet     ! net ocean surface energy budget W m^-2
    REAL(dp), DIMENSION(:), INTENT(in) :: pevapi     ! evaporation over ice
    REAL(dp), DIMENSION(:), INTENT(in) :: pcvsi      ! snow cover over ice (fraction of grid box)
    REAL(dp), DIMENSION(:), INTENT(in) :: pfri       ! ice cover (fraction of grid box)
    !
    ! Arguments INTENT(in out)
    ! --------------------
    REAL(dp), DIMENSION(:), INTENT(in out) :: pamlcorr   ! flux correction for mixed layer ocean
    REAL(dp), DIMENSION(:), INTENT(in out) :: pseaice    ! seaice fraction rel to ocean
    REAL(dp), DIMENSION(:), INTENT(in out) :: psiced     ! seaice thickness
    REAL(dp), DIMENSION(:), INTENT(in out) :: ptsi       ! surface temperature of ice
    REAL(dp), DIMENSION(:), INTENT(in out) :: ptsw       ! surface temperature of water
    REAL(dp), DIMENSION(:), INTENT(in out) :: pfluxres   ! residual heat flux in W m-2
    REAL(dp), DIMENSION(:), INTENT(in out) :: pahfres    ! residual heat flux over ice in W m-2
    REAL(dp), DIMENSION(:), INTENT(in out) :: psni       ! water equivalent of snow on ice in m
! op_pj_20160617+
!!$    REAL(dp), DIMENSION(:), INTENT(out)    :: pamlcorac  ! accumulated mixed layer fluc correction
!!$    REAL(dp), DIMENSION(:), INTENT(out)    :: pamlheatac ! accumulated mixed layer heat ??
    REAL(dp), DIMENSION(:), INTENT(out)    :: pamlheat ! mixed layer heat ??
! op_pj_20160617-
    !
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr='mlocean_mlo'
    INTEGER :: jl
    REAL(dp) :: ziscond
    REAL(dp) :: zfw
    REAL(dp) :: zheat
    REAL(dp) :: zfbase
    REAL(dp) :: zrhoilf, zicecap
    REAL(dp) :: zdtrilf
    REAL(dp) :: zfreez              ! threshold for ice formation
    REAL(dp) :: zmixcap, zmcapdt
    REAL(dp) :: zmcaprilf
    REAL(dp) :: zfluxw              ! net surface heat flux
    REAL(dp) :: zts                 ! surface temperature/ preliminary
    REAL(dp) :: zfres               ! 
    REAL(dp) :: zconhflx            ! conductive heat flux (H_c)
    REAL(dp) :: zhi                 ! preliminary ice thickness
    REAL(dp) :: zsubice             ! sublimation/ deposition (E_i)
    REAL(dp) :: zcpdt
    REAL(dp) :: zsniced             ! effective ice thickness (h_eff)
    !
    ! Executable statements
    !
    ! 1. Set up constants
    !
    ziscond   = zalpha / zalphas * rho_sea/ zrho_sn
    zrhoilf   = zrhoice * alf
    zicecap   = zrhoice * zcpice * zdice  ! heat capacity of ice slab 
    zdtrilf   = delta_time / zrhoilf
    zfreez    = -zdice / zdtrilf          ! threshold for ice formation
    zmixcap   = rho_sea * csw * mldmix    ! heat capacity of the mixed layer C_w
    zmcapdt   = delta_time / zmixcap      ! delta_t / C_w
    zmcaprilf = zmixcap / zrhoilf
    zcpdt     = zicecap / delta_time
    !
    ! 2. Mixed layer ocean temperature and ice thickness
    !
    DO jl = 1, kproma
       !
       IF (palake(jl) == 0._dp .AND. pslm(jl) < 0.5_dp) THEN  ! ocean points; no lake points
          !
          IF (pseaice(jl) < 0.5_dp) THEN                      ! open water
             !
             !  Ocean temperature (ptsw) calculated from 
             !    - the net surface heat flux (zfluxw);
             !    - the residual flux from last time step (pfluxres -> unrealised freezing);
             !    - the flux correction (pamlcorr -> provided from input file).
             !
             zfluxw = pfsnet(jl) - pamlcorr(jl)
             zts    = ptsw(jl)   + zmcapdt * (zfluxw + pfluxres(jl))
             IF (zts < ctfreez) THEN
                pamlcorr(jl) = MIN(0._dp, pamlcorr(jl))
                zfluxw = pfsnet(jl) - pamlcorr(jl)
                zts    = ptsw(jl)   + zmcapdt * (zfluxw + pfluxres(jl))
             END IF
             ptsi(jl)     = ctfreez
             pfluxres(jl) = 0._dp
             psiced(jl)   = 0._dp
             IF (zts >= ctfreez) THEN                     ! open water (unchanged)
                ptsw(jl) = zts                            ! case 1: no ice possible
             ELSE                                         ! check ice formation
                ptsw(jl) = ctfreez
                zfres    = (zts - ctfreez)/ zmcapdt       ! < 0.
                IF (zfres <= zfreez) THEN                 ! case 3: ice formation
                   psiced(jl)   = zmcaprilf*(ctfreez-zts) ! > zdice
                   pseaice(jl)  = 1._dp
                ELSE                                      ! case 2: no ice formation yet
                   pfluxres(jl) = zfres
                END IF
             END IF
             
          ELSE IF (psiced(jl) >= zdice) THEN              ! sea ice is present
          
             IF (lonorth(jl)) THEN
                zfbase = fbase_north
             ELSE
                zfbase = fbase_south
             END IF
             !
             !       Ice thickness (psiced)
             !
             zsniced  = psiced(jl) + ziscond * psni(jl) ! h_eff - effective ice thickness
             zconhflx = zalpha * (ptsi(jl) - ctfreez)/ zsniced
             zsubice  = (1._dp - pcvsi(jl)) * pevapi(jl) * delta_time/ zrhoice
             pamlcorr(jl) = MIN(0._dp, pamlcorr(jl)) - zfbase
             zhi      = psiced(jl) - zdtrilf * (zconhflx + pfluxres(jl) - pamlcorr(jl)) + zsubice
             ptsw(jl) = ctfreez
             IF (zhi >= zdice) THEN
                psiced(jl)   = zhi
                pseaice(jl)  = 1._dp
                pfluxres(jl) = 0._dp
             ELSE IF (zhi <= 0._dp) THEN                 ! complete melting
                ptsw(jl)     = ctfreez - zhi/ zmcaprilf  ! ptsw > ctfreez
                psiced(jl)   = 0._dp
                pseaice(jl)  = 0._dp
                pfluxres(jl) = -rho_H2O * alf * psni(jl)/ delta_time
                psni(jl)     = 0._dp
             ELSE                                        ! incomplete melting
                psiced(jl)   = zdice
                pseaice(jl)  = 1._dp
                pfluxres(jl) = (zdice - zhi)/ zdtrilf
                pahfres(jl)  = pahfres(jl) - delta_time * pfri(jl) * pfluxres(jl)
             END IF
          END IF
          zfw = 1._dp - pfri(jl) - pslm(jl)
          ! accumulate mixed layer variables
          zheat = zmixcap * zfw * ptsw(jl) - zrhoilf * pfri(jl) * psiced(jl)

! op_pj_20160617+
!!$       pamlheatac(jl) = pamlheatac(jl) + zheat * delta_time
!!$       pamlcorac(jl)  = pamlcorac (jl) + pamlcorr(jl) * delta_time
          pamlheat(jl) = zheat
! op_pj_20160617-
       END IF
    END DO
    !  ---------------------------------------------------------------------
    RETURN
  END SUBROUTINE mlocean_mlo
  !  -----------------------------------------------------------------------
  !  -----------------------------------------------------------------------
  SUBROUTINE mlocean_mlicetemp (kproma, delta_time    &
                      , psiced,     palake,  pslf     &
                      , ptrfli,     psofli            &
                      , pahfsi,     pahfli,  pfri     &
                      , ptsi,       psni,    pfluxres &
                      , pahfres,    pqres,   pahfice, pahfcon ) 
    ! Description:
    !
    ! Prognostic calculation of sea-ice temperature
    !
    ! Method:
    ! (taken from sicetemp.f90 
    !        IF (lmlo) THEN
    ! )
    !
    ! Authors:
    !
    ! U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
    ! U. Schlese, MPI December 2002, ice thickness over ocean
    !
    ! for more details see file AUTHORS
    !
    ! Arguments INTENT(in)
    ! --------------------
    INTEGER,                INTENT(in) :: kproma
    REAL(dp),               INTENT(in) :: delta_time 
    REAL(dp), DIMENSION(:), INTENT(in) :: psiced   ! ice depth in m
    REAL(dp), DIMENSION(:), INTENT(in) :: palake   ! lake fraction of grid box
    REAL(dp), DIMENSION(:), INTENT(in) :: pslf     ! sea land fraction (land == 1)
    REAL(dp), DIMENSION(:), INTENT(in) :: ptrfli   ! LW flux over ice in W m-2
    REAL(dp), DIMENSION(:), INTENT(in) :: psofli   ! SW flux over ice in W m-2
    REAL(dp), DIMENSION(:), INTENT(in) :: pahfsi   ! sensible heat flux over ice in W m-2
    REAL(dp), DIMENSION(:), INTENT(in) :: pahfli   ! latent heat flux over ice in W m-2
    REAL(dp), DIMENSION(:), INTENT(in) :: pfri     ! ice cover (fraction of grid box)
    !
    ! Arguments INTENT(in out)
    ! --------------------
    REAL(dp), DIMENSION(:), INTENT(in out) :: ptsi     ! surface temperature of ice in K
    REAL(dp), DIMENSION(:), INTENT(in out) :: psni     ! water equivalent of snow on ice in m
    REAL(dp), DIMENSION(:), INTENT(in out) :: pfluxres ! residual heat flux in W m-2
    REAL(dp), DIMENSION(:), INTENT(in out) :: pahfres  ! residual heat flux over ice in W m^-2
    REAL(dp), DIMENSION(:), INTENT(in out) :: pqres    ! residual heat flux for melting sea ice in W m^-2
    REAL(dp), DIMENSION(:), INTENT(in out) :: pahfice  ! conductive heat flux in W m^-2
    REAL(dp), DIMENSION(:), INTENT(in out) :: pahfcon  ! conductive heat flux through ice in W m^-2
    !
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr='mlocean_mlicetemp'
    INTEGER  :: jl
    REAL(dp) :: ziscond
    REAL(dp) :: zcpcon
    REAL(dp) :: zcpdt
    REAL(dp) :: zsniced  ! effective ice thickness (h_eff)
    REAL(dp) :: zicefl
    REAL(dp) :: zsflx    ! net ocean surface energy budget + residual term (H + R_m)
    !
    !  Executable statements
    !
    !-- 1. Set up constants
    !
    ziscond = zalpha / zalphas * rho_sea / zrho_sn
    zcpcon  = zrhoice * zcpice * zdice
    zcpdt   = zcpcon / delta_time
    !
    !-- 2. Compute new skin-temperature
    !    
    DO jl = 1, kproma
    
       IF (palake(jl) == 0._dp) THEN          ! ocean points
          IF (psiced(jl) >= zdice) THEN       ! ice covered ocean points
 
             zsniced  = psiced(jl) + ziscond * psni(jl)
             zicefl   = zalpha * ctfreez / zsniced 
             zsflx    = ptrfli(jl) + psofli(jl) + pahfsi(jl) + pahfli(jl) &
                      + pfluxres(jl)
             pfluxres(jl) = 0._dp
             ptsi(jl) = (zcpdt * ptsi(jl) + zsflx + zicefl)/  &
                        (zcpdt + zalpha/zsniced)

             IF (ptsi(jl) > tmelt) THEN
                pfluxres(jl) = (zcpdt + zalpha/zsniced) * (ptsi(jl) - tmelt)
                ptsi(jl)     = tmelt
             END IF
             pahfres(jl) = pahfres(jl) + delta_time * pfri(jl) * pfluxres(jl)
             pahfice(jl) = zalpha * (ptsi(jl) - ctfreez) / zsniced
          ELSE                                ! water
            pahfice(jl) = 0._dp
            ptsi(jl)    = tmelt
            psni(jl)    = 0._dp
         END IF
         pahfcon(jl)    = pahfcon(jl) + delta_time * pfri(jl) * pahfice(jl)
       END IF
    END DO
    !
    RETURN
  END SUBROUTINE mlocean_mlicetemp
  !  -----------------------------------------------------------------------
END MODULE messy_mlocean_mlo
! **********************************************************************
