MODULE MESSY_SCAV_INTER

  USE messy_scav_mem
  USE messy_scav_l_kpp,             ONLY: str_field_kpp_l => spc_names
  USE messy_scav_i_kpp,             ONLY: str_field_kpp_i => spc_names
  USE messy_main_constants_mem,     ONLY: dp
  USE messy_main_tools,             ONLY: PTR_2D_ARRAY

  IMPLICIT NONE
  REAL(dp), PUBLIC, DIMENSION(:,:), POINTER   :: gboxarea_2d => NULL()
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:)  :: idx_evap_num
  REAL(dp), PUBLIC                            :: thres_val = 1.e-6_dp

  TYPE proc
    TYPE(PTR_2D_ARRAY), DIMENSION(:,:), POINTER :: frac_eva => NULL()
  END TYPE proc
  TYPE(proc), DIMENSION(:), POINTER, PUBLIC :: process => NULL()
  INTEGER, PUBLIC, PARAMETER :: proc_idx_rain = 1
  INTEGER, PUBLIC, PARAMETER :: proc_idx_snow = 2
  INTEGER, PUBLIC, PARAMETER :: proc_idx_lwc  = 3
  INTEGER, PUBLIC, PARAMETER :: proc_idx_iwc  = 4

  TYPE, PUBLIC :: t_aermod_made3
    INTEGER          :: ks        = 0
    INTEGER          :: km        = 0
    INTEGER          :: ki        = 0
    INTEGER          :: as        = 0
    INTEGER          :: am        = 0
    INTEGER          :: ai        = 0
    INTEGER          :: cs        = 0
    INTEGER          :: cm        = 0
    INTEGER          :: ci        = 0
    INTEGER          :: nmod      = 0
    INTEGER          :: i_mode    = 0
    CHARACTER(LEN=2) :: csubname  = ''
    INTEGER          :: evap_mode = 0
    INTEGER          :: n_resmod  = 0
  END TYPE t_aermod_made3
  TYPE(t_aermod_made3), PUBLIC, ALLOCATABLE, DIMENSION(:) :: made3

  PRIVATE
  SAVE
  PUBLIC :: in_flux, out_flux, in_cloud, out_cloud
  PUBLIC :: melting, sedi_cloud
  PUBLIC :: evapo, evapo_aer
  PUBLIC :: scav_cvdep
!  PUBLIC :: pack_aerosol, unpack_aerosol
  PUBLIC :: aer_scav_drv
  PUBLIC :: scav_main
  PUBLIC :: cloudtrac2aer, cloud2trac, aer2cloudtrac, calc_evapo_frac
  PUBLIC :: trac2l, trac2i
  PUBLIC :: proc
  PUBLIC :: get_core_mode_made3

CONTAINS
!===============================================================================

  SUBROUTINE in_flux(input_field, pxtp1, aerflux, vol, oldcov, &
                     ntrac, lproma, js, phase, nspec)

    USE messy_scav_l_kpp,   ONLY: lspec => nspec
    USE messy_scav_i_kpp,   ONLY: ispec => nspec

    IMPLICIT NONE
    
    INTRINSIC :: MAX
    INTEGER,  INTENT(IN)    :: ntrac, lproma, js, phase, nspec
    REAL(dp), INTENT(IN)    :: input_field(nspec,lproma), pxtp1(ntrac,lproma), &
                               vol(lproma), oldcov(lproma)
    REAL(dp), INTENT(INOUT) :: aerflux(aer_count(js)%numbers(c_all), lproma)

    REAL(dp)                :: ppmx(0:ntrac, lproma), vol2cm(lproma)
    INTEGER                 :: jt

    vol2cm(1:lproma) = vol(1:lproma) * 1.e6_dp
    ppmx(0,1:lproma) = 0.0_dp 
    ppmx(1:ntrac,1:lproma) = pxtp1(1:ntrac,1:lproma)
!   eliminate negative values before KPP

    SELECT CASE (PHASE)

    CASE(1)
      L_SPEC(:,:) = 0.0_dp
      DO jt=1,lspec
        L_SPEC(jt,1:lproma) = MAX(input_field(jt,1:lproma)  / &
                            ( vol2cm(1:lproma) * oldcov(1:lproma) ), 0.0_dp)
      ENDDO
      CALL in_kpp_l(ppmx, L_SPEC(0:lspec,1:lproma), ntrac, js, lproma, lspec)
    CASE(2)
      I_SPEC(:,:) = 0.0_dp
      DO jt=1,ispec
        I_SPEC(jt,1:lproma) = MAX(input_field(jt,1:lproma)  / &
                            ( vol2cm(1:lproma) * oldcov(1:lproma) ), 0.0_dp)
      ENDDO
      CALL in_kpp_i(ppmx, I_SPEC(0:ispec,1:lproma), ntrac, js, lproma, ispec)
    END SELECT

    IF (cpl_aerosol > 0) THEN
      DO jt=1,aer_count(js)%numbers(c_all)
        aerflux(jt,1:lproma) = aerflux(jt,1:lproma) / &
                             ( vol2cm(1:lproma) * oldcov(1:lproma) )
      ENDDO

      SELECT CASE (PHASE)
      CASE(1)
        CALL aer2l(aerflux, L_SPEC(0:lspec,1:lproma), js, lproma, lspec)
      CASE(2)
        CALL aer2i(aerflux, I_SPEC(0:ispec,1:lproma), js, lproma, ispec)
      END SELECT

      DO jt=1,aer_count(js)%numbers(c_all)
        aerflux(jt,1:lproma) = aerflux(jt,1:lproma) * &
                               vol2cm(1:lproma) * oldcov(1:lproma) 
      ENDDO

    ENDIF
    RETURN
  END SUBROUTINE in_flux
!-------------------------------------------------------------------------------

  SUBROUTINE out_flux(output_field, pxtp1, pHrain, lwc, vol, oldcov, &
                      ntrac, lproma, js, phase, kspec, cprod_arr)

    USE messy_scav_l_kpp,     ONLY: lspec => nspec, IND_Hp_l
    USE messy_scav_i_kpp,     ONLY: ispec => nspec
    USE messy_main_constants_mem,   ONLY: avo => N_A

    IMPLICIT NONE

    INTRINSIC :: MAX, LOG10
    INTEGER,  INTENT(IN)    :: ntrac, lproma, js, kspec
    INTEGER,  INTENT(IN)    :: phase   ! 1 = liquid, 2 = ice
    REAL(dp), INTENT(INOUT) :: output_field(kspec, lproma), pxtp1(ntrac, lproma)
    REAL(dp), INTENT(INOUT) :: pHrain(lproma)
    REAL(dp), INTENT(IN)    :: vol(lproma), oldcov(lproma), lwc(lproma)
    ! compiler workaround for g95 (0.92, 0.93) allocate from 0
    REAL(dp), INTENT(INOUT) :: cprod_arr(lproma,0:nprod)

    REAL(dp) :: ppmx(0:ntrac,lproma), zfw2a, vol2cm(lproma)
    INTEGER  :: jt, jl, idx2

    L_SPEC(0,1:lproma) = 0.0_dp
    vol2cm(1:lproma) = vol(1:lproma) * 1.e6_dp
    ppmx(:,1:lproma) = 0.0_dp 

    ppmx(1:ntrac,1:lproma) = pxtp1(1:ntrac,1:lproma)

    SELECT CASE (PHASE)

    CASE(1)
      CALL out_kpp_l(ppmx, L_SPEC(0:lspec,1:lproma), ntrac, js, lproma, lspec)
      DO jt=1,lspec
        output_field(jt,1:lproma) = L_SPEC(jt,1:lproma) * &
                                    oldcov(1:lproma) * vol2cm(1:lproma)
      ENDDO
    CASE(2)
      CALL out_kpp_i(ppmx, I_SPEC(0:ispec,1:lproma), ntrac, js, lproma, ispec)
      DO jt=1,ispec
        output_field(jt,1:lproma) = I_SPEC(jt,1:lproma) * &
                                    oldcov(1:lproma) * vol2cm(1:lproma)
      ENDDO
    END SELECT

    DO jt=1,ntrac
      pxtp1(jt,1:lproma) = (1._dp - oldcov(1:lproma)) * pxtp1(jt,1:lproma) + &
                            oldcov(1:lproma) * ppmx(jt,1:lproma) 
    ENDDO

    if (phase == 1) then
!   calculation of rainpH
      DO jl=1,lproma
        zfw2a = lwc(jl) * avo * 1.e-6
        pHrain(jl) = -LOG10(MAX(L_SPEC(IND_hp_l,jl)/zfw2a, 1.e-20_dp))
      ENDDO
! calculation of production rates
      DO jt=1,nprod
        idx2 = PR(js)%PROD_IDX(jt)
        DO jl=1,lproma
          cprod_arr(jl,jt) =  cprod_arr(jl,jt) + &
            MAX(L_SPEC(idx2,jl),0._dp)
        END DO
      END DO
    END if

    RETURN
  END SUBROUTINE out_flux

!-------------------------------------------------------------------------------
! Used to release chemically active species (input_field) from the cloud phase
! to the gas and aerosol phases (pxtp1).
! NOTE: In case MADE3 is used as aerosol microphysics submodel, only gas species
! are released here! Aerosol species are released via evapo_aer.
  SUBROUTINE evapo(input_field, pxtp1, vol, frac, &
                   ntrac, lproma, lidx, js, jk, jrow, phase, kspec)

    USE messy_scav_l_kpp,  ONLY: lspec => nspec,           &
                                       IND_CO2_l, IND_H2O_l
    USE messy_scav_i_kpp,  ONLY: ispec => nspec, &
                                      ! IND_CO2_i, 
                                 IND_H2O_i !, kpp_hno3_i
    USE messy_scav,        ONLY: NCB

    IMPLICIT NONE

    INTEGER,  INTENT(IN)    :: ntrac, lproma, js, jrow, jk, phase, kspec
    INTEGER,  INTENT(IN)    :: lidx(lproma)
    REAL(dp), INTENT(IN)    :: vol(lproma)
    REAL(dp), INTENT(IN)    :: frac(lproma)
    REAL(dp), INTENT(INOUT) :: pxtp1(ntrac, lproma), input_field(kspec, lproma)

    REAL(dp)                :: ppmx(0:ntrac, lproma), vol2cm(lproma)
    REAL(dp)                :: inp_field(0:kspec, lproma)
    INTEGER                 :: jt, jl, idx2

    vol2cm(1:lproma) = vol(1:lproma) * 1.e6_dp
    ppmx(0,1:lproma) = 0.0_dp
    inp_field(:,:) = 0.0_dp
    ppmx(1:ntrac,1:lproma) = pxtp1(1:ntrac,1:lproma)

    SELECT CASE(PHASE)
    CASE(1)
      DO jt=1,lspec
        inp_field(jt,1:lproma)   = input_field(jt,1:lproma) / vol2cm(1:lproma)
      END DO
     
!  as CO2 is a fixed species its liquid concentration should not influence the gas phase and is set to 0.
      inp_field(IND_CO2_l,1:lproma) = 0._dp
      inp_field(IND_H2O_l,1:lproma) = 0._dp

!   redistributing soluted trace species into gas and aerosol species

!     fraction of material evaporating into evap_max_mode, 
!     rest in the next smaller mode

      IF (cpl_aerosol > 0) THEN
        CALL l2evap(inp_field, ppmx, ntrac, js, lproma, kspec, 1, frac)

        CALL l2evap(inp_field, ppmx, ntrac, js, lproma, kspec, 2, frac)
      ENDIF
    CASE(2)
      DO jt=1,ispec
        inp_field(jt,1:lproma)   = input_field(jt,1:lproma) / vol2cm(1:lproma)
      END DO
!  as CO2 is a fixed species its liquid concentration should not influence the gas phase and is set to 0.
!      inp_field(IND_CO2_i,1:lproma) = 0.
      inp_field(IND_H2O_i,1:lproma) = 0.
!   redistributing soluted trace species into gas and aerosol species
      IF (cpl_aerosol > 0) &
        CALL i2evap(inp_field, ppmx, ntrac, js, lproma, kspec)

    END SELECT

    IF (l_lg) THEN
      DO jl=1,lproma
        idx2 = lidx(jl)
        IF ( NINT(NCB(idx2,jk,jrow)) >= 1 ) THEN
          input_field(1:kspec,jl) = 0._dp
          pxtp1(1:ntrac,jl) = ppmx(1:ntrac,jl)
        ENDIF
      ENDDO
    ELSE
      DO jl=1,lproma
        input_field(1:kspec,jl) = 0._dp
        pxtp1(1:ntrac,jl) = ppmx(1:ntrac,jl)
      ENDDO
    ENDIF

    RETURN
  END SUBROUTINE evapo

!-------------------------------------------------------------------------------

  SUBROUTINE cloud2trac(input_field, pxtp1, &
                        ntrac, lproma, js, phase, kspec)

    USE messy_scav_l_kpp,  ONLY: lspec => nspec,           &
                                 IND_CO2_l, IND_H2O_l
    USE messy_scav_i_kpp,  ONLY: ispec => nspec, &
                                 IND_CO2_i, IND_H2O_i !, kpp_hno3_i

    IMPLICIT NONE

    INTEGER,  INTENT(IN)    :: ntrac, lproma, js,  phase, kspec
    REAL(dp), INTENT(INOUT) :: pxtp1(ntrac, lproma), input_field(kspec, lproma)

    REAL(dp)                :: ppmx(0:ntrac, lproma)
    REAL(dp)                :: inp_field(0:kspec, lproma)
    INTEGER                 :: jt

    ppmx(0,1:lproma) = 0.0_dp
    inp_field(:,:) = 0.0_dp
    ppmx(1:ntrac,1:lproma) = pxtp1(1:ntrac,1:lproma)

    SELECT CASE(PHASE)
    CASE(1)
      DO jt=1,lspec
        inp_field(jt,1:lproma)   = input_field(jt,1:lproma) 
      END DO
     
!  as CO2 and H2O are fixed species, their liquid concentration should not influence
!  the gas phase and is set to 0.
      inp_field(IND_CO2_l,1:lproma) = 0._dp
      inp_field(IND_H2O_l,1:lproma) = 0._dp

      CALL l2trac(inp_field, ppmx, ntrac, js, lproma, kspec)

    CASE(2)
      DO jt=1,ispec
        inp_field(jt,1:lproma)   = input_field(jt,1:lproma) 
      END DO
!  as CO2 and H2O are fixed species, their ice concentration should not influence 
!  the gas phase and is set to 0.

      inp_field(IND_CO2_i,1:lproma) = 0.
      inp_field(IND_H2O_i,1:lproma) = 0.
!   redistributing soluted trace species into gas and aerosol species

      CALL i2trac(inp_field, ppmx, ntrac, js, lproma, kspec)

    END SELECT
    
    pxtp1(1:ntrac,1:lproma) = ppmx(1:ntrac,1:lproma)
    input_field(1:kspec,1:lproma) = inp_field(1:kspec,1:lproma)
    
    RETURN
  END SUBROUTINE cloud2trac
!-------------------------------------------------------------------------------

! Used to redistribute chemically inactive species (aerflux) to the aerosol
! phase (pxtp1).
! NOTE: In case MADE3 is used as aerosol microphysics submodel, chemically
! active aerosol species are also released here.
  SUBROUTINE evapo_aer(pxtp1, aerflux, vol, frac, &
                      ntrac, lproma, lidx, js, jk, jrow)

    USE messy_scav,      ONLY: NCB
    IMPLICIT NONE

    INTEGER,  INTENT(IN)    :: ntrac, js, lproma, jrow, lidx(lproma), jk
    REAL(dp), INTENT(IN)    :: vol(lproma)
    REAL(dp), INTENT(IN)    :: frac(lproma,ntrac,n_resmod)
    REAL(dp), INTENT(INOUT) :: aerflux(aer_count(js)%numbers(c_all), lproma), &
                               pxtp1(ntrac, lproma)

    REAL(dp)                :: ppmx(0:ntrac, lproma), vol2cm(lproma), &
                               aer_field(0:ntrac, lproma)

    INTEGER                 :: jt, jl, idx, idx2
    INTEGER                 :: jmod
    REAL(dp)                :: hlp(lproma)

    vol2cm(1:lproma)             = vol(1:lproma) * 1.e6_dp
    ppmx(0,1:lproma)             = 0.0_dp
    ppmx(1:ntrac, 1:lproma)      = pxtp1(1:ntrac,1:lproma)
    aer_field(0:ntrac, 1:lproma) = 0.0_dp
    CALL aerspec2aertrac(ntrac, js, lproma, aerflux, aer_field)

    IF (i_evap == 2) THEN
      DO jt=1,ntrac
!    aerosol numbers should be decreased by cloud processing due to 
!    coagulation of cloud droplets -> after evaporation less aerosol particles
!    should remain
        IF (attr(js)%evap_trac(n_resmod,jt).EQ.idx_evap_num(js)) THEN
          hlp(:) = frac_resnum * aer_field(jt,1:lproma) / vol2cm(1:lproma)
        ELSE
          hlp(:) = aer_field(jt,1:lproma) / vol2cm(1:lproma)
        ENDIF

!    ! each tracer is evaporated to the same species of the
!    ! mode(s) assigned in scav_init_coupling

        DO jmod = 1, n_resmod
          idx = attr(js)%evap_trac(jmod,jt)

!    aerosol numbers should be decreased by cloud processing due to 
!    coagulation of cloud droplets -> after evaporation less aerosol particles
!    should remain
!    reduction to 25% is too strong and yields low AOD
!    reduction to 50% is not too bad with OC SOA hydrophlic particle emissions

          DO jl=1,lproma
             ppmx(idx,jl) = ppmx(idx,jl) + hlp(jl) * frac(jl,jt,jmod)
          ENDDO

        ENDDO


      ENDDO
    ELSE 
      DO jt=1,ntrac
        !    ! each tracer is evaporated to the one where it comes from
        DO jl=1,lproma
          ppmx(jt,jl) = ppmx(jt,jl) + aer_field(jt,jl) / vol2cm(jl)
        ENDDO
      ENDDO
    ENDIF
    
    IF (l_lg) THEN
      DO jl=1,lproma
        ppmx(jt,jl) = ppmx(jt,jl) + aer_field(jt,jl) / vol2cm(jl)
      ENDDO
    ENDIF

    IF (l_lg) THEN
      DO jl=1,lproma
        idx2 = lidx(jl)
        IF ( NINT(NCB(idx2,jk,jrow)) >= 1 ) THEN
          aerflux(1:aer_count(js)%numbers(c_all),jl) = 0.0_dp
          pxtp1(1:ntrac,jl) = ppmx(1:ntrac,jl)
        ENDIF
      ENDDO
    ELSE
      aerflux(1:aer_count(js)%numbers(c_all),1:lproma) = 0.0_dp
      pxtp1(1:ntrac,1:lproma) = ppmx(1:ntrac,1:lproma)
    ENDIF
    
    RETURN
    
  END SUBROUTINE evapo_aer
!-------------------------------------------------------------------------------

  SUBROUTINE in_cloud(inp_field, pxtp1, aerflux, vol, cover, &
                      ntrac, lproma, js, phase, kspec)

    USE messy_scav_l_kpp,         ONLY: lspec => nspec
    USE messy_scav_i_kpp,         ONLY: ispec => nspec

    IMPLICIT NONE

    INTRINSIC :: MAX
    INTEGER,  INTENT(IN)    :: ntrac, js, lproma, phase, kspec
    REAL(dp), INTENT(IN)    :: inp_field(kspec, lproma), pxtp1(ntrac, lproma), &
                               vol(lproma), cover(lproma)
    REAL(dp), INTENT(INOUT) :: aerflux(aer_count(js)%numbers(c_all), lproma)

    REAL(dp)                :: ppmx(0:ntrac, lproma), vol2cm(lproma)
    INTEGER                 :: jt
    
!!!
!    if (cover.lt.1e-7_dp) RETURN
!!!
    vol2cm(1:lproma) = vol(1:lproma) * 1.e6_dp
    ppmx(0,1:lproma) = 0.0_dp 
    ppmx(1:ntrac,1:lproma) = pxtp1(1:ntrac,1:lproma)
!   eliminate negative values before KPP  
    SELECT CASE (PHASE)
    CASE(1)
      L_SPEC(:,:) = 0.0_dp
      DO jt=1,lspec
        L_SPEC(jt,1:lproma) = MAX(inp_field(jt,1:lproma)  / &
                            ( vol2cm(1:lproma) * cover(1:lproma) ), 0.0_dp)
      ENDDO
      CALL in_kpp_l(ppmx, L_SPEC(0:lspec,1:lproma), ntrac, js, lproma, lspec)
    CASE(2)
      I_SPEC(:,:) = 0.0_dp
      DO jt=1,ispec
        I_SPEC(jt,1:lproma) = MAX(inp_field(jt,1:lproma)  / &
                            ( vol2cm(1:lproma) * cover(1:lproma) ), 0.0_dp)
      ENDDO
      CALL in_kpp_i(ppmx, I_SPEC(0:ispec,1:lproma), ntrac, js, lproma, ispec)
    END SELECT

    IF (cpl_aerosol.GT.0) THEN
      DO jt=1,aer_count(js)%numbers(c_all)
        aerflux(jt,1:lproma) = aerflux(jt,1:lproma) / &
                             ( vol2cm(1:lproma) * cover(1:lproma))
      ENDDO
      
      SELECT CASE (PHASE)
      CASE(1)
        CALL aer2l(aerflux, L_SPEC(0:lspec,1:lproma), js, lproma, lspec)
      CASE(2)
        CALL aer2i(aerflux, I_SPEC(0:ispec,1:lproma), js, lproma, ispec)
      END SELECT
      

      DO jt=1,aer_count(js)%numbers(c_all)
        aerflux(jt,1:lproma) = aerflux(jt,1:lproma) * &
                               vol2cm(1:lproma) * cover(1:lproma)
      ENDDO
    ENDIF
    RETURN

  END SUBROUTINE in_cloud
!-------------------------------------------------------------------------------
  SUBROUTINE out_cloud(cloud_field, rain_field, pxtp1, clpH, Hp_l,&
                       cover, lwc, pmrat, vol, ntrac, lproma, js, &
                       phase, kspec, cprod_arr)

    USE messy_scav_l_kpp,           ONLY: lspec => nspec, IND_Hp_l
    USE messy_scav_i_kpp,           ONLY: ispec => nspec
    USE messy_main_constants_mem,   ONLY: avo => N_A

    INTRINSIC :: MIN, MAX, LOG10
    INTEGER,  INTENT(IN)    :: ntrac, js, lproma, phase, kspec
    REAL(dp), INTENT(IN)    :: cover(lproma), lwc(lproma), &
                               pmrat(lproma), vol(lproma)
    REAL(dp), INTENT(INOUT) :: cloud_field(kspec, lproma), &
                               rain_field(kspec, lproma),  &
                               pxtp1(ntrac, lproma),       &
                               clpH(lproma)
    ! compiler workaround for g95 (0.92, 0.93) allocate from 0:
    REAL(dp), INTENT(INOUT)   :: Hp_l(lproma), cprod_arr(lproma,0:nprod)

    REAL(dp)                :: ppmx(0:ntrac, lproma)
    REAL(dp)                :: zfac(lproma), vol2cm(lproma), zfw2a
    INTEGER :: jt, jl, idx2

!!!
!    if (cover.lt.1e-7_dp) RETURN
!!!
    L_SPEC(0,1:lproma) = 0.0_dp
    zfac(1:lproma)     = MAX(MIN(pmrat(1:lproma) / lwc(1:lproma), 1.0_dp),0._dp)
    vol2cm(1:lproma)   = vol(1:lproma) * 1.e6_dp
    ppmx(0,1:lproma)   = 0.0_dp 
    ppmx(1:ntrac,1:lproma) = pxtp1(1:ntrac,1:lproma)

    SELECT CASE (PHASE)
    CASE(1)
      CALL out_kpp_l(ppmx, L_SPEC(0:lspec,1:lproma), ntrac, js, lproma, lspec)
      DO jt=1,lspec
        cloud_field(jt,1:lproma) = L_SPEC(jt,1:lproma) * cover(1:lproma)
        rain_field(jt,1:lproma)  = rain_field(jt,1:lproma) + &
                                   cloud_field(jt,1:lproma) * zfac(1:lproma) * &
                                   vol2cm(1:lproma) 
        cloud_field(jt,1:lproma) = cloud_field(jt,1:lproma) * &
                                   (1._dp - zfac(1:lproma)) * vol2cm(1:lproma)
      ENDDO
    CASE(2)
      CALL out_kpp_i(ppmx, I_SPEC(0:ispec,1:lproma), ntrac, js, lproma, ispec)
      DO jt=1,ispec
        cloud_field(jt,1:lproma) = I_SPEC(jt,1:lproma) * cover(1:lproma)
        rain_field(jt,1:lproma)  = rain_field(jt,1:lproma) + &
                                   cloud_field(jt,1:lproma) * zfac(1:lproma) * &
                                   vol2cm(1:lproma) 
        cloud_field(jt,1:lproma) = cloud_field(jt,1:lproma) * &
                                   (1._dp - zfac(1:lproma)) * vol2cm(1:lproma)
      ENDDO
    END SELECT

    DO jt=1,ntrac
      pxtp1(jt,1:lproma) = (1.0_dp - cover(1:lproma)) * pxtp1(jt,1:lproma) + &
                            cover(1:lproma) * ppmx(jt,1:lproma)
    ENDDO
  
    IF (phase == 1) THEN
    !    calculation of the cloud pH
      DO jl=1,lproma
        zfw2a   = lwc(jl) * avo * 1.e-6_dp
        clpH(jl) = -LOG10(MAX(L_SPEC(IND_Hp_l,jl) / zfw2a, 1.e-20_dp))
        Hp_l(jl) =  L_SPEC(IND_Hp_l,jl)
      ENDDO
      ! calculation of production rates
      DO jt=1,nprod
        idx2 = PR(js)%PROD_IDX(jt)
        DO jl=1,lproma
          cprod_arr(jl,jt) =  cprod_arr(jl,jt) + &
            MAX(L_SPEC(idx2,jl),0._dp) * cover(jl)
        END DO
      END DO
    ENDIF
    RETURN
  END SUBROUTINE out_cloud

!-------------------------------------------------------------------

 SUBROUTINE sedi_cloud(lkppfld, lnewflx, lkppsediflx, laersediflx, &
                       pmisedi, iwc, lproma, js, kspec)

    USE messy_scav_i_kpp,   ONLY: ispec => nspec


    INTRINSIC :: MIN
    INTEGER,  INTENT(IN)    :: js, lproma, kspec
    REAL(dp), INTENT(IN)    :: iwc(lproma), pmisedi(lproma)
    REAL(dp), INTENT(INOUT) ::    &
      lkppsediflx(kspec, lproma), &
      lkppfld(kspec, lproma)
    REAL(dp), INTENT(INOUT) ::    &
      laersediflx(aer_count(js)%numbers(c_all), lproma), &
      lnewflx(aer_count(js)%numbers(c_all), lproma)

    REAL(dp)                :: zfac(lproma)
    INTEGER :: jt

    zfac(1:lproma)     = MAX(0._dp, &
                         MIN(pmisedi(1:lproma) / iwc(1:lproma), 1.0_dp))
    DO jt=1,ispec
      lkppsediflx(jt,1:lproma) = lkppfld(jt,1:lproma) * zfac(1:lproma) 
      lkppfld(jt,1:lproma)     = lkppfld(jt,1:lproma) * (1._dp - zfac(1:lproma))
    enddo

    DO jt = 1,aer_count(js)%numbers(c_all)
      laersediflx(jt,1:lproma) = lnewflx(jt,1:lproma) * zfac(1:lproma) 
      lnewflx(jt,1:lproma)     = lnewflx(jt,1:lproma) * (1._dp - zfac(1:lproma))
    ENDDO

    RETURN
  END SUBROUTINE sedi_cloud

!-------------------------------------------------------------------------------

  SUBROUTINE MELTING(lproma,   js,      kspec,    lmelt,   lsnow, &
                     lkppiflx, lkppflx, laeriflx, laerflx )

    IMPLICIT NONE

    INTEGER,  INTENT(IN)    :: lproma, js, kspec
    REAL(dp), INTENT(IN)    :: lmelt(lproma), lsnow(lproma)
    REAL(dp), INTENT(INOUT) :: lkppiflx(kspec,lproma), lkppflx(kspec,lproma)
    REAL(dp), INTENT(INOUT) :: laeriflx(aer_count(js)%numbers(c_all),lproma)
    REAL(dp), INTENT(INOUT) :: laerflx(aer_count(js)%numbers(c_all),lproma)

    REAL(dp) :: frac(lproma)
    INTEGER  :: jl, jt, idx1, idx2

    frac(:) = 0._dp
    DO jl=1,lproma
      IF (lsnow(jl) > 1e-15_dp) &
        frac(jl) = MAX(0._dp,MIN(1._dp, lmelt(jl)/lsnow(jl)))
    enddo

    DO jt=1,ispec_ice
      idx2 = KPP_I_IDX(js)%ice_spec(jt,ice2liq)
      idx1 = KPP_I_IDX(js)%ice_spec(jt,ice_idx)
      do jl=1,lproma
        lkppflx(idx2,jl)  = lkppflx(idx2,jl) + frac(jl) * lkppiflx(idx1,jl)
        lkppiflx(idx1,jl) = (1._dp - frac(jl)) * lkppiflx(idx1,jl) 
      enddo
    enddo
    DO jt = 1,aer_count(js)%numbers(c_all)
      do jl=1,lproma
        laerflx(jt,jl)  = laerflx(jt,jl) + frac(jl) * laeriflx(jt,jl)
        laeriflx(jt,jl) = (1._dp - frac(jl)) * laeriflx(jt,jl) 
      enddo
    enddo
    

  END SUBROUTINE MELTING
!===============================================================================

  SUBROUTINE scav_cvdep (kproma, kbdim,  klev,   jk,  jrow, ztmst, &
       ktrac,  ktopm2, ldcum,  kctop,  &
       kcbot,  kdtop,  lddraf, pmfu,   pmfd,   pdmfup, pdmfdp,                 &
       pmfuxt, pmfdxt, zrainh, zcover_conv, pten, zmfu)

    ! ----------------------------------------------------------------------
    ! mz_ht_20030303+, This subroutines does some of the wet deposition 
    !     calculations normally being included directly in the cumulus 
    !     subroutine cuflx but has moved to this subroutine to facilitate
    !     the exchange of this routine. The subroutine requires some input
    !     parameters from cuflx to calculate the in-cloud convective 
    !     scavenging. The tendencies are then subsequently being used in the 
    !     subroutine scav_physc to update the tracer concentrations. The 
    !     main reason why these calculations are not performed in scav_physc
    !     but in cuflx.f90 is that not only the input parameters from this 
    !     routine are required for the calculations but also that the 
    !     convective tracer mass fluxes must be updated within the vertical
    !     loop for the gas-aqueous phase transfer and subsequent aqueous-
    !     phase chemical transformations. 
    !     ==================================================================
    !     NOTE that the calculations in this routine are generally based on
    !     tracer concentrations in molecules cm-3 and consequently when the
    !     concentrations are given in volume mixing ratio, like in MECCA
    !     some recalculations are introduced in this routine, which should
    !     be removed/changed when other concentrations are applied 
    !     ==================================================================
    ! 
    !     First efforts by Geert-Jan Roelofs
    !     Developed and implemented by Holger Tost, 2003 to May 2004
    ! ----------------------------------------------------------------------         

    USE messy_main_constants_mem,       ONLY: avo => N_a, amd => M_air !, g
    USE messy_scav                  
    USE messy_scav_aer,                 ONLY: aer_scav, aer_scav_cloud, &
                                              aermodel, max_mode
    USE messy_scav_l_kpp,               ONLY: nspec
   
    IMPLICIT NONE 

    SAVE

    INTRINSIC :: ABS, ASSOCIATED, MAX, MIN, TINY, REAL

    INTEGER  :: kproma, kbdim, klev, ktopm2, &
                jrow, jk, jl, i, jt, jh, ktrac, js
    REAL(dp) :: pmfu(kbdim,klev),   zmfu(kbdim), pmfd(kbdim,klev),   &
                pdmfup(kbdim,klev), pdmfdp(kbdim,klev)

    INTEGER :: kctop(kbdim),         kcbot(kbdim),    kdtop(kbdim)
    LOGICAL :: ldcum(kbdim),         lddraf(kbdim)

    REAL(dp) :: pmfuxt(kbdim,klev,ktrac),pmfdxt(kbdim,klev,ktrac)
    REAL(dp) :: zrainh(kbdim),  zcover_conv(kbdim)
    
    REAL(dp) ::   pxtu(kbdim,ktrac),  pxtd(kbdim,ktrac)
  

    REAL(dp) :: zt(kproma,klev),    zrhoa(kproma,klev)
    REAL(dp) :: zrlwc(kproma),                        &
                mc(kbdim,klev),       cm(kbdim,klev)

    REAL(dp) :: zrdrad(kproma)
    REAL(dp) :: pten(kbdim,klev)
  
    REAL(dp) :: ztmst,  zevapco

    REAL(dp) :: zwashfrac(kbdim)

    REAL(dp) :: snowflux(kbdim), zriwc(kbdim), zlwc(kbdim)

    REAL(dp), PARAMETER :: epsilon = TINY(0.0_dp)

    ! =========================================================================
    ! variables for aerosol scav

    INTEGER  :: l
    REAL(dp) :: aer_rad(max_mode)
    REAL(dp) :: dummy_field_cloud(nspec), kpp_field_cloud2(nspec,kproma), &
                new_aerflx(kbdim, aer_count(1)%numbers(c_all)) 
!----------------------------------------------------------------------
! new packing for scavenging
    INTEGER  :: lproma, mproma, naermo, nmode, num
    INTEGER  :: lidx(kproma), midx(kproma)
    REAL(dp) :: llwc(kproma), llwc2(kproma), liwc(kproma)
    REAL(dp) :: lrdrad(kproma), lpress(kproma), ltemp(kproma)
    REAL(dp) :: lprec(kproma), lsnow(kproma), lcov(kproma), lrho(kproma)
    REAL(dp) :: lgrvol(kproma), lprod(kproma)
    REAL(dp) :: loc_aer_rad(max_mode, kproma)
    REAL(dp) :: lkppflx(nspec,kproma), lkppfld(nspec,kproma)
    REAL(dp) :: lph(kproma), lhp(kproma)
    REAL(dp) , ALLOCATABLE :: xt_vec(:,:), laerflx(:,:), mxt_vec(:,:), l_xt(:,:)
    REAL(dp) :: l_aer_rad(max_mode,kproma)
    REAL(dp) :: l_philfrac(max_mode,kproma)
    ! compiler workaround for g95 (0.92, 0.93) allocate from 0
    REAL(dp) :: cprod_array(kproma,0:nprod)
    REAL(dp) :: lsteps(kproma), lrsteps(kproma)
    !+ for output of partitioning parameters

    REAL(dp) :: mkppflx(nspec,kproma), kpp_rain_cv_help(kproma)
    REAL(dp) :: mgrvol(kproma), mcov(kproma)
    REAL(dp) :: evfrac(kproma)
    INTEGER  :: lwork(kproma), mwork(kproma)
    REAL(dp) :: lfrac(max_mode,kproma)

!----------------------------------------------------------------------
    IF (.NOT.lscav_cv) RETURN
    js = 1
    IF (ALLOCATED(l_xt)) DEALLOCATE(l_xt)
    ALLOCATE(l_xt(aer_count(js)%numbers(c_all),kproma))
    IF (ALLOCATED(xt_vec)) DEALLOCATE(XT_vec)
    ALLOCATE(XT_VEC(ktrac,kproma))
    IF (ALLOCATED(mxt_vec)) DEALLOCATE(MXT_vec)
    ALLOCATE(MXT_VEC(ktrac,kproma))
    IF (ALLOCATED(laerflx)) DEALLOCATE(laerflx)
    ALLOCATE(laerflx(aer_count(js)%numbers(c_all), kproma))

    dummy_field_cloud(:) = 0.0_dp
    new_aerflx(:,:)      = 0.0_dp
    zrhoa(:,:)           = 0.0_dp
    zriwc(:)             = 0.0_dp
    zlwc(:)              = 0.0_dp
    snowflux(:)          = 0.0_dp

    pxtu(:,:)            = 0.0_dp
    pxtd(:,:)            = 0.0_dp
    aer_rad(:)           = 2.0e-6_dp       ! pseudo aerosol radius in case of aerosol tracers,
                                           ! but no aerosol module running
    l_philfrac(:,:)      = 1.0_dp

    zrlwc(:)             = 0.0_dp
    zwashfrac(:)         = 0.0_dp
    zrdrad(:)            = 0.0_dp
    lfrac(:,:)           = 0.0_dp

    IF (jk.EQ.ktopm2) zrainh(:)=0.

    ! HT- calculating the concentration changes due to convective transport
    !     with cumulus clouds if any of the nconv tracer switches are set to
    !     1

    IF (l_anynconvect) THEN

      DO jl=1,kproma
        zt(jl,jk)   = pten(jl,jk) 
        kconbot(jl,jrow)= 0.
      ENDDO
      DO jh=1,klev
        DO jl=1,kproma
          IF (grvol(jl,jh,jrow) > 0._dp) &
            zrhoa(jl,jh)= grmass(jl,jh,jrow)/grvol(jl,jh,jrow) ! kg m-3
        ENDDO
      ENDDO

      DO jl=1,kproma
      ! convert mixing ratio [mol/mol] to concentration [mcl/cc]
      ! mc = mixing ratio to concentration
!       this should work if calculation of grmass and grvol use actual pressure
!       and not as it is at the moment the values from timestep before
        mc(jl,jk) = avo * zrhoa(jl,jk) * 1.e-3_dp /amd

      ! convert concentration [mcl/cc] to mixing ratio [mol/mol]
      ! cm = concentration to mixing ratio
        cm(jl,jk)= 1._dp/mc(jl,jk)
     
        IF (ldcum(jl))&
          kconbot(jl,jrow) = REAL(kcbot(jl),dp)

      ENDDO
 
!       print*,'vor der schleife der jl',p_pe
      DO jl=1,kproma

          ! mz_HT_20020131 in-cloud scavenging calculated here, the budget
          !     calculations are performed in cudtdq
         
!          print*, 'vor IF 1 bei jl= ',jl, p_pe
        IF (ldcum(jl).AND.jk.GE.kctop(jl)-1.AND.jk.LE.kcbot(jl)) THEN

             ! calculating the re-evaporation of dissolved 
             !     gases in the downdrafts. 
!             print*, "vor IF 2 bei jl= ",jl, p_pe
          IF (lddraf(jl).AND.jk.GE.kdtop(jl).AND.pdmfdp(jl,jk-1).LT.0.) THEN
            zevapco=0.
!                print*, "vor IF 3 bei jl= ",jl, p_pe
            IF (zrainh(jl).GT.epsilon) zevapco=pdmfdp(jl,jk-1)/zrainh(jl)

            zevapco=MIN(zevapco,0._dp)
            zevapco=MAX(zevapco,-1._dp)
!                print*, 'zevapco ',zevapco, p_pe
            zrainh(jl)=zrainh(jl)+pdmfdp(jl,jk-1)
!                print*, zrainh(jl),pdmfdp(jl,jk-1), 'zrainh, pdmfdp', p_pe


! ===============================================================================
          ENDIF
        ENDIF
      ENDDO


!   calculate cloud cover
       DO jl=1,kproma
          CALL calc_covercv (conv_cover(jl,jk,jrow), zmfu(jl), &
                             zrhoa(jl,jk), zcover_conv(jl))
       ENDDO
     
!  calc of actual tracer concentrations
       do jt=1, ktrac
         do jl=1,kproma
           if (abs(zmfu(jl)).gt.0_dp) &
             pxtu(jl,jt) = pmfuxt(jl,jk,jt) * mc(jl,jk) / zmfu(jl) 
           if (abs(pmfd(jl,jk)).gt.0_dp) &
             pxtd(jl,jt) = pmfdxt(jl,jk,jt) * mc(jl,jk) / pmfd(jl,jk) 
         ENDDO
       ENDDO

!------------------------------------------------

!   impaction scavenging :
       IF (lscav_imp) THEN

!    calculation of the rain liquid water content of the box and the 
!    input parameters for the chemistry interactions      

       DO jl=1,kproma
          zrainh(jl)=zrainh(jl)+pdmfup(jl,jk)
          CALL  CALC_LWC_BC_CV (zrainh(jl), zcover_conv(jl), &
               zrlwc(jl), zwashfrac(jl), zrdrad(jl), jk, kconbot(jl,jrow))
          IF (zrainh(jl) > 1.e-10_dp) &
          zlwc(jl) = zrainh(jl)/zcover_conv(jl) * ztmst *&
                     gboxarea_2d(jl,jrow)/grvol(jl,jk,jrow)
       ENDDO
!===============================================================================

!   aerosol scavenging included here
        IF (lscav_aer) THEN

        !------------------------------------
! pack for rain and snow
! into a vector that contains all the boxes with liquid impaction scavenging
          lproma = 0
          llwc(:)          = 0._dp             
          lrdrad(:)        = 0._dp           
          lpress(:)        = 0._dp       
          ltemp(:)         = 0._dp        
          lprec(:)         = 0._dp        
          lcov(:)          = 0._dp         
          lrho(:)          = 0._dp         
          lgrvol(:)        = 0._dp       
          lsnow(:)         = 0._dp   
          xt_vec(:,:)      = 0._dp     
          loc_aer_rad(:,:) = 0._dp
          laerflx(:,:)     = 0._dp
          lidx(:)          = 0

          lwork(:) = 0
          DO jl=1,kproma
              IF ( (zlwc(jl) > thres_val) .AND. &
                 (zwashfrac(jl) > 1.e-7_dp) ) THEN
              lproma = lproma + 1
              lwork(lproma)=jl
            ENDIF
          ENDDO
          DO jl=1,lproma
            i=lwork(jl)
            llwc(jl)    = zrlwc(i)
            lrdrad(jl)  = zrdrad(i)
            lpress(jl)  = press_3d(i,jk,jrow)
            ltemp(jl)   = zt(i,jk)
            lprec(jl)   = zrainh(i)
            lcov(jl)    = zwashfrac(i)
            lrho(jl)    = zrhoa(i,jk)
            lgrvol(jl)  = grvol(i,jk,jrow)
            lsnow(jl)   = snowflux(i)
            lidx(jl)    = i
          END DO
          DO jt=1,ktrac
            DO jl=1,lproma
              i=lwork(jl)
              xt_vec(jt,jl)   = pxtu(i,jt) 
            ENDDO
          ENDDO
          DO jt=1,aer_count(js)%numbers(c_all)
            DO jl=1,lproma
              i=lwork(jl)
              laerflx(jt,jl) = aero_flx(js)%aer_flx_cv(i,1,jt,jrow)
            ENDDO
          ENDDO

          IF (lproma > 0) THEN

            DO naermo = 1, aermodel(js)%aeromodnum
              nmode = aermodel(js)%aer_input(naermo)%lmode
              num   = aermodel(js)%aer_mass_attr(naermo)%number
              !    scavenging of aerosol masses
              l_xt(:,1:lproma) = 0._dp
              CALL pack_aerosol(num, lproma, &
                       l_xt(1:num,1:lproma), & ! zx
                       xt_vec(:,1:lproma),                             & ! pmx
                       ktrac, 1, js, naermo)
              IF (ASSOCIATED(aermodel(js)%aer_input(naermo)%wetradius)) THEN 
                do l=1,nmode
                  DO jl=1,lproma
                    loc_aer_rad(l,jl)  = &
                      aermodel(js)%aer_input(naermo)%wetradius(lidx(jl),jk,&
                      l,jrow) 
                  ENDDO
                ENDDO
              ELSE
                DO l=1,nmode
                  loc_aer_rad(l,1:lproma) = aer_rad(l)
                ENDDO
              ENDIF
              DO l=1,nmode
                l_aer_rad(l,1:lproma) = loc_aer_rad(l,1:lproma) * &
                  aermodel(js)%aer_input(naermo)%cmr2mmr(l)
              ENDDO
              
              CALL aer_scav(num, lproma, js,      &
                           nmode, naermo, 1,                                &
                           l_xt(1:num,1:lproma),  &
                           laerflx(1:aer_count(js)%numbers(c_all),1:lproma),&
                           ztmst, lprec(1:lproma), lcov(1:lproma),          &
                           lrdrad(1:lproma), l_aer_rad(1:nmode,1:lproma),   &
                           ltemp(1:lproma), lpress(1:lproma),               &
                           lrho(1:lproma), lsnow(1:lproma),                 &
                           lgrvol(1:lproma), 1)

              CALL unpack_aerosol(num, lproma,                  &
                     l_xt(1:num,1:lproma), xt_vec(:,1:lproma),  &
                     ktrac, 1, js, naermo)

!    scavenging of aerosol numbers
              l_xt(:,1:lproma) = 0._dp
              num   = aermodel(js)%aer_num_attr(naermo)%number
              CALL pack_aerosol(num, lproma, &
                       l_xt(1:num,1:lproma), & ! zx
                       xt_vec(:,1:lproma),                            & ! pmx
                       ktrac, 2, js, naermo)
              DO jt=1,aer_count(js)%numbers(c_all)
                DO jl=1,lproma
                  l_xt(jt,jl) = l_xt(jt,jl) / avo
                ENDDO
              ENDDO
              DO l=1,nmode
                l_aer_rad(l,1:lproma) = loc_aer_rad(l,1:lproma) 
              ENDDO

              CALL aer_scav(num, lproma, js,         &
                            nmode, naermo, 2,                                 &
                            l_xt(1:num,1:lproma),    &
                            laerflx(1:aer_count(js)%numbers(c_all),1:lproma), &
                            ztmst, lprec, lcov, lrdrad,                       &
                            l_aer_rad(1:nmode,1:lproma), ltemp,               &
                            lpress, lrho, lsnow, lgrvol, 1)

              DO jt=1,aer_count(js)%numbers(c_all)
                DO jl=1,lproma
                  l_xt(jt,jl) = l_xt(jt,jl) * avo
                ENDDO
              ENDDO

              CALL unpack_aerosol(num, lproma,&
                     l_xt(1:num,1:lproma),    &! zx
                     xt_vec(:,1:lproma),                               &! pmx
                     ktrac, 2, js, naermo)
            ENDDO
          ENDIF

          DO jt=1,ktrac
            DO jl=1,lproma
              i=lwork(jl)
              pxtu(i,jt) = xt_vec(jt,jl) 
            ENDDO
          ENDDO
          DO jt=1,aer_count(js)%numbers(c_all)
            DO jl=1,lproma
              i=lwork(jl)
              aero_flx(js)%aer_flx_cv(i,1,jt,jrow) =          &
                laerflx(jt,jl) 
            ENDDO
          ENDDO

          
        ENDIF        ! scav_aer
!===============================================================================
        IF (lscav_gas) THEN
           
            lproma = 0
            llwc(:)          = 0._dp             
            lrdrad(:)        = 0._dp           
            lpress(:)        = 0._dp       
            ltemp(:)         = 0._dp         
            lprec(:)         = 0._dp        
            lcov(:)          = 0._dp         
            lrho(:)          = 0._dp         
            lgrvol(:)        = 0._dp       
            lsnow(:)         = 0._dp   
            xt_vec(:,:)      = 0._dp 
            laerflx(:,:)     = 0._dp
            lkppflx(:,:)     = 0._dp
            lidx(:)          = 0
            cprod_array(:,:) = 0._dp
            lsteps(:)        = 0._dp
            lrsteps(:)       = 0._dp
            mproma = 0
            mkppflx(:,:)     = 0._dp
            mxt_vec(:,:)     = 0._dp 
            mgrvol(:)        = 0._dp
            mcov(:)          = 0._dp
            midx(:)          = 0

            kpp_rain_cv_help(:)=0.0_dp
            DO jt=1,nspec
              DO jl=1,kproma
                kpp_rain_cv_help(jl)=MAX(kpp_rain_cv_help(jl), &
                  kpp_rain_cv(jl,js,jt,jrow))
              ENDDO
            ENDDO

            lwork(:) = 0
            mwork(:) = 0
            DO jl=1,kproma
              IF ( (zlwc(jl) > thres_val) .AND. &
                 (zwashfrac(jl) > 1.e-7_dp) ) THEN
                lproma = lproma + 1               
                lwork(lproma)=jl
              ELSE
                IF (kpp_rain_cv_help(jl) > 0.0_dp) THEN
                  mproma = mproma + 1
                  mwork(mproma)=jl
                ENDIF
              ENDIF
            ENDDO
            DO jt=1,ktrac
              DO jl=1,lproma
                i=lwork(jl)
                xt_vec(jt,jl) = pxtu(i,jt)
              ENDDO
            ENDDO
            DO jt=1,nspec
              DO jl=1,lproma
                i=lwork(jl)
                lkppflx(jt,jl)  = kpp_rain_cv(i,js,jt,jrow)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                laerflx(jt,jl) = aero_flx(js)%aer_flx_cv(i,1,jt,jrow)
              ENDDO
            ENDDO
            DO jl=1,lproma
              i=lwork(jl)
              lgrvol(jl)  = grvol(i,jk,jrow)
              lcov(jl)    = zwashfrac(i)
              llwc(jl)    = zlwc(i)
              llwc2(jl)   = zrlwc(i)
              ltemp(jl)   = zt(i,jk)
              lpress(jl)  = press_3d(i,jk,jrow)
              lrdrad(jl)  = zrdrad(i)
              lidx(jl)    = i
              lph(jl)     = ph(js)%rainpH_cv(i,jk,jrow)
            ENDDO
            DO jt=1,nspec
              DO jl=1,mproma
                i=mwork(jl)
                mkppflx(jt,jl) = kpp_rain_cv(i,js,jt,jrow)
              ENDDO
            ENDDO
            DO jt=1,ktrac
              DO jl=1,mproma
                i=mwork(jl)
                mxt_vec(jt,jl) = pxtd(i,jt)
              ENDDO
            ENDDO
            DO jl=1,mproma
             i=mwork(jl)
             mgrvol(jl) = grvol(i,jk,jrow)
             mcov(jl)   = zwashfrac(i)
             midx(jl)   = i
           ENDDO

           IF (lproma > 0) THEN
             CALL in_flux(lkppflx, xt_vec, laerflx, lgrvol, lcov, &
                          ktrac, lproma, js, 1, nspec)

             IF (lscav_easy) THEN
               CALL scav_easy_liq(llwc2, ltemp, lproma, js)
             ELSE
               CALL scav_aqchem_kpp(ltemp, lpress, llwc2, llwc, lrdrad, &
                                    ztmst, 2, lproma, lidx, jk, jrow, js, &
                                    lsteps, lrsteps)
             ENDIF
              
             CALL out_flux(lkppflx, xt_vec, lph, llwc, lgrvol, lcov, &
                           ktrac, lproma, js, 1, nspec,              &
                           cprod_array(1:lproma,:))
           ENDIF

           IF (mproma > 0) CALL evapo(mkppflx, mxt_vec, mgrvol, evfrac, &
                                      ktrac, mproma, midx, js, jk, jrow,&
                                      1, nspec)
           DO jt=1,ktrac
             DO jl=1,lproma
               i=lwork(jl)
               pxtu(i,jt) = xt_vec(jt,jl) 
             ENDDO
           ENDDO
           DO jt=1,nspec
             DO jl=1,lproma
               i=lwork(jl)
               kpp_rain_cv(i,js,jt,jrow) = lkppflx(jt,jl)
             ENDDO
           ENDDO
           DO jt=1,aer_count(js)%numbers(c_all)
             DO jl=1,lproma
               i=lwork(jl)
               aero_flx(js)%aer_flx_cv(i,1,jt,jrow) = laerflx(jt,jl)
             ENDDO
           ENDDO
           DO jl=1,lproma
             i=lwork(jl)
             ph(js)%rainpH_cv(i,jk,jrow) = lph(jl)
           ENDDO
           DO jt=1,ktrac
             DO jl=1,mproma
               i=mwork(jl)
               pxtd(i,jt) = mxt_vec(jt,jl) 
             ENDDO
           ENDDO
           DO jt=1,nspec
             DO jl=1,mproma
               i=mwork(jl)
               kpp_rain_cv(i,js,jt,jrow) = mkppflx(jt,jl)
             ENDDO
           ENDDO

!==============================
    
          ENDIF    ! scav_gas
!-------------------------------------------------------------------------------
!     convective in cloud scavenging is considered in this section
!  nucleation scavenging from freshly formed convective precipitation
        ENDIF    ! scav_imp
        
        IF (lscav_nuc) THEN
!   nucleation scavenging :
!===============================================================================
        DO jl=1,kproma
          IF (pdmfup(jl,jk).GT.0._dp) &
            CALL calc_lwc_ic_cv(pdmfup(jl,jk), zcover_conv(jl),          &
                                zrlwc(jl), zwashfrac(jl), zrdrad(jl), jk,&
                                kconbot(jl,jrow))
          zrdrad(jl) = 1.75e-2_dp
          zlwc(jl)   = zrlwc(jl)
        ENDDO
!===============================================================================
       !   aerosol scavenging included here
        IF (lscav_aer) THEN
! pack for rain and snow
! into a vector that contains all the boxes with liquid nucleation scavenging
            lproma = 0
            llwc(:)          = 0._dp
            liwc(:)          = 0._dp
            lrdrad(:)        = 0._dp
            lpress(:)        = 0._dp
            ltemp(:)         = 0._dp
            lcov(:)          = 0._dp
            lgrvol(:)        = 0._dp
            xt_vec(:,:)      = 0._dp
            loc_aer_rad(:,:) = 0._dp
            laerflx(:,:)     = 0._dp
            lidx(:)          = 0
            
            lwork(:) = 0
            DO jl=1,kproma
              IF ( (zrlwc(jl) > thres_val) .AND. &
                 (zwashfrac(jl) > 1.e-7_dp) ) THEN
                lproma = lproma + 1
                lwork(lproma)=jl
              ENDIF
            ENDDO
            Do jl=1,lproma
              i=lwork(jl)
              llwc(jl)    = zrlwc(i)
              liwc(jl)    = zriwc(i)
              lrdrad(jl)  = zrdrad(i)
              lpress(jl)  = press_3d(i,jk,jrow)
              ltemp(jl)   = zt(i,jk)
              lcov(jl)    = conv_cover(i,jk,jrow)
              lgrvol(jl)  = grvol(i,jk,jrow)
              lidx(jl)    = i
            enddo
            DO jt=1,ktrac
              DO jl=1,lproma
                i=lwork(jl)
                xt_vec(jt,jl)   = pxtu(i,jt)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                laerflx(jt,jl) = new_aerflx(i,jt)
              ENDDO
            ENDDO
            
            IF (lproma > 0) THEN
             
              DO naermo = 1, aermodel(js)%aeromodnum
                nmode = aermodel(js)%aer_input(naermo)%lmode
                !    scavenging of aerosol masses
                num   = aermodel(js)%aer_mass_attr(naermo)%number
                l_xt(:,1:lproma) = 0._dp
                CALL pack_aerosol(num, lproma, &
                         l_xt(1:num,1:lproma), & ! zx
                         xt_vec(:,1:lproma),                             & ! pmx
                         ktrac, 1, js, naermo)
                IF (ASSOCIATED(aermodel(js)%aer_input(naermo)%wetradius)) THEN 
                  DO l=1,nmode
                    DO jl=1,lproma
                      loc_aer_rad(l,jl)  = &
                        aermodel(js)%aer_input(naermo)%wetradius(lidx(jl),jk,&
                        l,jrow) 
                    enddo
                  ENDDO
                ELSE
                  DO l=1,nmode
                    loc_aer_rad(l,1:lproma) = aer_rad(l)
                  ENDDO
                ENDIF
                DO l=1,nmode
                  l_aer_rad(l,1:lproma) = loc_aer_rad(l,1:lproma) * &
                    aermodel(js)%aer_input(naermo)%cmr2mmr(l)
                ENDDO
                IF (ASSOCIATED(aermodel(js)%aer_input(naermo)%mfracnuc)) THEN 
                  DO l=1,SIZE(aermodel(js)%aer_input(naermo)%mfracnuc,3)
                    DO jl=1,lproma
                      lfrac(l,jl) = aermodel(js)%aer_input(naermo)%mfracnuc(lidx(jl),jk,l,jrow)
                    END DO
                  END DO
                END IF

                CALL aer_scav_cloud(num, lproma,     &
                   js, nmode, naermo, 1, 1,                                    &
                   l_xt(1:num,1:lproma),             &
                   laerflx(1:aer_count(js)%numbers(c_all),1:lproma),           &
                   ztmst, llwc(1:lproma), liwc(1:lproma),                      &
                   lcov(1:lproma), lrdrad(1:lproma),                           &
                   l_aer_rad(1:nmode,1:lproma), ltemp(1:lproma),               &
                   lpress(1:lproma), lgrvol(1:lproma), lfrac(:,1:lproma),      &
                   l_philfrac(1:nmode,1:lproma) )
               
                CALL unpack_aerosol(num, lproma, &
                          l_xt(1:num,1:lproma),  & 
                          xt_vec(:,1:lproma),                              & 
                          ktrac, 1, js, naermo)

                !    scavenging of aerosol numbers
                l_xt(:,1:lproma) = 0._dp
                num   = aermodel(js)%aer_num_attr(naermo)%number
                CALL pack_aerosol(num, lproma, &
                         l_xt(1:num,1:lproma), & ! zx
                         xt_vec(:,1:lproma),                            & ! pmx
                         ktrac, 2, js, naermo)
                
                DO jt=1,aer_count(js)%numbers(c_all)
                  do jl=1,lproma
                    l_xt(jt,jl) = l_xt(jt,jl) / avo
                  enddo
                ENDDO
                DO l=1,nmode
                  l_aer_rad(l,1:lproma) = loc_aer_rad(l,1:lproma)
                ENDDO
                IF (ASSOCIATED(aermodel(js)%aer_input(naermo)%nfracnuc)) THEN 
                  DO l=1,SIZE(aermodel(js)%aer_input(naermo)%nfracnuc,3)
                    DO jl=1,lproma
                      lfrac(l,jl) = aermodel(js)%aer_input(naermo)%nfracnuc(lidx(jl),jk,l,jrow)
                    END DO
                  END DO
                END IF
                CALL aer_scav_cloud(num, lproma,     &
                  js, nmode, naermo, 2, 1,                                    &
                  l_xt(1:num,1:lproma),              &
                  laerflx(1:aer_count(js)%numbers(c_all),1:lproma),           &
                  ztmst, llwc(1:lproma), liwc(1:lproma),                      &
                  lcov(1:lproma), lrdrad(1:lproma),                           &
                  l_aer_rad(1:nmode,1:lproma), ltemp(1:lproma),               &
                  lpress(1:lproma), lgrvol(1:lproma), lfrac(1:nmode,1:lproma),&
                  l_philfrac(1:nmode,1:lproma) )
           
                DO jt=1,aer_count(js)%numbers(c_all)
                  do jl=1,lproma
                    l_xt(jt,jl) = l_xt(jt,jl) * avo
                  enddo
                ENDDO
                                
                CALL unpack_aerosol(num, lproma, &
                         l_xt(1:num,1:lproma),   & 
                         xt_vec(:,1:lproma),                              & 
                         ktrac, 2, js, naermo)
              ENDDO
            ENDIF

            DO jt=1,ktrac
              DO jl=1,lproma
                i=lwork(jl)
                pxtu(i,jt) = xt_vec(jt,jl)   
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                new_aerflx(i,jt) = laerflx(jt,jl) 
              ENDDO
            ENDDO

          ENDIF          ! on lscav_aer
!===============================================================================
        IF (lscav_gas) THEN
             
           kpp_field_cloud2(1:nspec,:) = 0.0_dp
      
           lproma = 0
           llwc(:)          = 0._dp             
           lrdrad(:)        = 0._dp           
           lpress(:)        = 0._dp       
           ltemp(:)         = 0._dp        
           lprec(:)         = 0._dp        
           lprod(:)         = 0._dp        
           lcov(:)          = 0._dp         
           lrho(:)          = 0._dp         
           lgrvol(:)        = 0._dp       
           xt_vec(:,:)      = 0._dp 
           laerflx(:,:)     = 0._dp
           lkppflx(:,:)     = 0._dp
           lkppfld(:,:)     = 0._dp
           lidx(:)          = 0
           cprod_array(:,:) = 0._dp

           lwork(:) = 0    
           DO jl=1,kproma
             IF ( (zlwc(jl) > thres_val) .AND. &
                  (conv_cover(jl,jk,jrow) > 1.e-7_dp) ) THEN
               lproma = lproma + 1
               lwork(lproma)=jl
             ENDIF
           ENDDO
           DO jt=1,ktrac
             DO jl=1,lproma
               i=lwork(jl)
               xt_vec(jt,jl)   = pxtu(i,jt)
             ENDDO
           ENDDO
           DO jt=1,nspec
             DO jl=1,lproma
               i=lwork(jl)
               lkppflx(jt,jl)  = kpp_rain_cv(i,js,jt,jrow)
               lkppfld(jt,jl)  = kpp_field_cloud2(jt,i)
             ENDDO
           ENDDO
           DO jt=1,aer_count(js)%numbers(c_all)
             DO jl=1,lproma
               i=lwork(jl)
               laerflx(jt,jl) = new_aerflx(i,jt)
             ENDDO
           ENDDO
           DO jl=1,lproma
             i=lwork(jl)
             lgrvol(jl)  = grvol(i,jk,jrow)
             lcov(jl)    = conv_cover(i,jk,jrow)
             llwc(jl)    = zlwc(i)
             lprod(jl)   = zlwc(i)
             ltemp(jl)   = zt(i,jk)
             lpress(jl)  = press_3d(i,jk,jrow)
             lrdrad(jl)  = zrdrad(i)
             lrho(jl)    = zrhoa(i,jk)
             lidx(jl)    = i
             lph(jl)     = ph(js)%cloudpH_cv(i,jk,jrow)
           ENDDO
         
           IF (lproma > 0) THEN
             CALL in_cloud(lkppfld, xt_vec, laerflx, lgrvol, lcov, &
                            ktrac, lproma, js, 1, nspec)
             IF (lscav_easy) THEN
               CALL scav_easy_liq(llwc, ltemp, lproma, js)
             ELSE
               CALL scav_aqchem_kpp(ltemp, lpress, llwc, llwc, lrdrad, &
                                     ztmst, 1, lproma, lidx, jk, jrow, js, &
                                     lsteps, lrsteps)
             ENDIF
             CALL out_cloud(lkppfld, lkppflx, xt_vec, lph, lhp, lcov, &
                             llwc, lprod, lgrvol, ktrac, lproma, js, 1, &
                             nspec, cprod_array(1:lproma,:) )
             lhp(1:lproma) = lhp(1:lproma) / (1.0e-6_dp*lrho(1:lproma))

             CALL evapo(lkppfld, xt_vec, lgrvol, evfrac, ktrac, lproma, &
                        lidx, js, jk, jrow, 1, nspec)
           ENDIF

           DO jt=1,ktrac
             DO jl=1,lproma
               i=lwork(jl)
               pxtu(i,jt) = xt_vec(jt,jl)
             ENDDO
           ENDDO
           DO jt=1,nspec
             DO jl=1,lproma
               i=lwork(jl)
               kpp_rain_cv(i,js,jt,jrow) = lkppflx(jt,jl)
             ENDDO
           ENDDO
           DO jt=1,aer_count(js)%numbers(c_all)
             DO jl=1,lproma
               i=lwork(jl)
               new_aerflx(i,jt) = laerflx(jt,jl) 
             ENDDO
           ENDDO
           DO jl=1,lproma
             i=lwork(jl)
             ph(js)%cloudpH_cv(i,jk,jrow)   = lph(jl)
             ph(js)%Hp_cloud_cv(i,jk,jrow)  = lhp(jl)
           ENDDO
       
         ENDIF   ! on scav_gas
      
       
        ! update of aerosol wet deposition flux by new_aer_flx

         DO jl=1,kproma
           IF (zlwc(jl) > 0._dp) THEN
             aero_flx(js)%aer_flx_cv(jl,1,                            &
               1:aer_count(js)%numbers(c_all),jrow) =                 &
               aero_flx(js)%aer_flx_cv(jl,1,                          &
               1:aer_count(js)%numbers(c_all),jrow) +                 &
               new_aerflx(jl,1:aer_count(js)%numbers(c_all)) 
           ENDIF
         ENDDO

!-------------------------------------------------------------------------------
        ENDIF   ! scav_nuc
          
        DO jt=1,ktrac
          DO jl=1,kproma
            if (abs(zmfu(jl)).gt.0._dp) &
              pmfuxt(jl,jk,jt) = zmfu(jl) * pxtu(jl,jt) * cm(jl,jk)
            if (abs(pmfd(jl,jk)).gt.0._dp) &
              pmfdxt(jl,jk,jt) = pmfd(jl,jk) * pxtd(jl,jt) * cm(jl,jk)
          ENDDO
        ENDDO

      ENDIF  ! on l_anyconvect

   RETURN

 END SUBROUTINE scav_cvdep

! ==============================================================================
! This subroutine calculates a convective cloud cover estimate 

  SUBROUTINE calc_covercv(conv_cover, updraft, &
                          zrhoa, zcover_conv)


    IMPLICIT NONE

    INTRINSIC :: MAX, MIN

    REAL(dp), INTENT(INOUT) :: conv_cover, zcover_conv
    REAL(dp), INTENT(IN)    :: updraft, zrhoa

!  local variables
    REAL(dp) :: vup_min

!      estimated convective cloud cover adapted from MATCH from Mark Lawrence

!     first determine the fraction swept out by precip from the core cloud;
!        the precipitating part is the cloud width needed to support the 
!        updraft mass flux with a nominal updraft velocity (of 1 m/s);
!        this might miss out on some diffs betw land and sea, but should be 
!        good for a first attempt; 

!     ATTENTION: Do not sum up pclcover and conv_cover for a total cover!
!                They might get higher values than 1 !!!!!!

              
          vup_min = 1._dp      ! [m/s] necessary updraft velocity to support convective rainfall
                               ! still to be tested or modified

          conv_cover = updraft / (vup_min * zrhoa) 
          zcover_conv = MAX(conv_cover,zcover_conv)
       

!     next compute a "coverage adjustment"; this is to help reduce the
!     time-step dependence of the process, since it is expected that for longer
!     time-steps the convective precip will tend to strip out the entire
!     precipitating column (so that without this adjustment, in twice as 
!     much time, the same total amount will still be scavenged); at the 
!     other end (short time steps), the tendency will not be availability-
!     limited, so that the deltat term in the tracer update will take care
!     of appropriate time-dependence.  (Could collect statistics in a future
!     study to see if the threshold chosen here is really appropriate...)
!     Assume:
!     1) the characteristic horiz dimension of any individual tower is 10 km;
!     2) the cells propogate horizontally with a characteristic speed of 10 m/s
   
      zcover_conv = MIN(1._dp, zcover_conv)

!     zcover_conv = 0.05   !old parametrization assuming a 5% convective cloud cover

  END SUBROUTINE calc_covercv



!==============================================================================
  SUBROUTINE AERSPEC2aertrac(ntrac, js, lproma, infield, outfield)

   USE messy_main_constants_mem,       ONLY: avo => N_a
   IMPLICIT NONE

   INTEGER,  INTENT(IN) :: ntrac, js, lproma
   REAL(dp), INTENT(IN) :: infield(aer_count(js)%numbers(c_all),lproma)
   REAL(dp), INTENT(INOUT) :: outfield(0:ntrac, lproma)

   INTEGER   :: jt, idx, jl
   INTRINSIC :: MAX
   
   DO jt = 1,ntrac
     IF ( attr(js)%log_att(jt,laerosol).AND. &
          attr(js)%log_att(jt,lwetdep) ) THEN
       idx = attr(js)%int_att(jt,aerosol_index)
       IF (attr(js)%int_att(jt,lquantity) == 2) THEN
         DO jl=1,lproma
           outfield(jt,jl) = MAX(infield(idx,jl) * avo,0.0_dp)
         ENDDO
       ELSE
         DO jl=1,lproma
           outfield(jt,jl) = MAX(infield(idx,jl),0.0_dp)
         ENDDO
       ENDIF
     ENDIF
   ENDDO

 END SUBROUTINE AERSPEC2aertrac
!-------------------------------------------------------------------------------

 SUBROUTINE in_kpp_l(pmx, d, ntrac, js, lproma, kspec)
  
   IMPLICIT NONE

   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(IN)    :: pmx(0:ntrac,lproma)
   REAL(dp), INTENT(INOUT) :: d(0:kspec,lproma)
   INTEGER :: idx1, idx2, jt, jl

   DO jt=1,LSPEC_GAS
     idx1 = KPP_L_IDX(js)%gas_spec(jt,gas_idx)
     idx2 = KPP_L_IDX(js)%gas_spec(jt,gas2trac) 
     DO jl=1,lproma
       D(idx1,jl) = pmx(idx2,jl)
     ENDDO
   ENDDO

 END SUBROUTINE in_kpp_l
!------------------------------------------------
 SUBROUTINE in_kpp_i(pmx, d, ntrac, js, lproma, kspec)
   
   IMPLICIT NONE

   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(IN)    :: pmx(0:ntrac,lproma)
   REAL(dp), INTENT(INOUT) :: d(0:kspec,lproma)
   INTEGER :: idx1, idx2, jt, jl

   DO jt=1,ISPEC_GAS
     idx1 = KPP_I_IDX(js)%gas_spec(jt,gas_idx)
     idx2 = KPP_I_IDX(js)%gas_spec(jt,gas2trac) 
     DO jl=1,lproma
       D(idx1,jl) = pmx(idx2,jl)
     ENDDO
   ENDDO

 END SUBROUTINE in_kpp_i
!-------------------------------------------------------------------------------

 SUBROUTINE out_kpp_l(pmx, d, ntrac, js, lproma, kspec)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: pmx(0:ntrac, lproma)
   REAL(dp), INTENT(INOUT) :: d(0:kspec, lproma)
   INTEGER :: idx1, idx2, jt, jl

   DO jt = 1,LSPEC_GAS
     idx1 = KPP_L_IDX(js)%gas_spec(jt,gas_idx)
     idx2 = KPP_L_IDX(js)%gas_spec(jt,gas2trac) 
     DO jl=1,lproma
       pmx(idx2,jl) = D(idx1,jl)
       D(idx1,jl)   = 0._dp
    ENDDO
   ENDDO

 END SUBROUTINE out_kpp_l
!-------------------------------------------------
 SUBROUTINE out_kpp_i(pmx, d, ntrac, js, lproma, kspec)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: pmx(0:ntrac, lproma)
   REAL(dp), INTENT(INOUT) :: d(0:kspec, lproma)
   INTEGER :: idx1, idx2, jt, jl

   DO jt = 1,ISPEC_GAS
     idx1 = KPP_I_IDX(js)%gas_spec(jt,gas_idx)
     idx2 = KPP_I_IDX(js)%gas_spec(jt,gas2trac) 
     DO jl=1,lproma
       pmx(idx2,jl) = D(idx1,jl)
       D(idx1,jl)   = 0._dp
     ENDDO
   ENDDO

 END SUBROUTINE out_kpp_i
!-------------------------------------------------------------------------------

 SUBROUTINE aer2l(aerflux, D, js, lproma, kspec) 

   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: aerflux(aer_count(js)%numbers(c_all),lproma)
   REAL(dp), INTENT(INOUT) :: D(0:kspec,lproma)
   INTEGER :: jt, idx1, jl

   DO jt = 1,aer_count(js)%numbers(c_all)
     idx1 = aer_attr(js)%aer_spec(jt,aer2l_idx)
     IF (IDX1 == 0) CYCLE
     DO jl=1,lproma
       D(idx1,jl) = D(idx1,jl) + MAX(aerflux(jt,jl),0._dp)
       aerflux(jt,jl) = 0._dp !MAX(aerflux(jt,jl),0._dp) * REAL(1-MIN(1,idx1),dp)
     ENDDO
   ENDDO
      
 END SUBROUTINE aer2l
!-------------------------------------------------------------------------------
 SUBROUTINE aer2i(aerflux, D, js, lproma, kspec)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: aerflux(aer_count(js)%numbers(c_all),lproma)
   REAL(dp), INTENT(INOUT) :: D(0:kspec,lproma)
   INTEGER :: jt, idx1, jl

   DO jt = 1,aer_count(js)%numbers(c_all)
     idx1 = aer_attr(js)%aer_spec(jt,aer2i_idx)
     IF (IDX1 == 0) CYCLE
     DO jl=1,lproma
       D(idx1,jl) = D(idx1,jl) + MAX(aerflux(jt,jl),0._dp)
       aerflux(jt,jl) = 0._dp !MAX(aerflux(jt,jl),0._dp) * REAL(1-MIN(1,idx1),dp)
     ENDDO
   ENDDO
      
 END SUBROUTINE aer2i
!-------------------------------------------------------------------------------
 SUBROUTINE l2evap(D, pmx, ntrac, js, lproma, kspec, idx0, frac)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec, idx0
   REAL(dp), INTENT(INOUT) :: pmx(0:ntrac,lproma)
   REAL(dp), INTENT(IN)    :: d(0:kspec, lproma)
   REAL(dp), INTENT(IN)    :: frac(lproma)
   INTEGER  :: idx1, idx2, jt, jl, idx
   REAL(dp) :: fract(lproma)

   IF (idx0 == 1) THEN
     idx = liq2evap
     fract(1:lproma) = frac(1:lproma)
   ELSE
     idx = liq2eva2
     fract(1:lproma) = 1._dp - frac(1:lproma)
   ENDIF

   DO jt=1,lspec_liq
     idx1 = KPP_L_IDX(js)%liq_spec(jt,liq_idx)
     idx2 = KPP_L_IDX(js)%liq_spec(jt,idx)
     DO jl=1,lproma
       pmx(idx2,jl) = pmx(idx2,jl) + D(idx1,jl) * fract(jl)
     ENDDO
   ENDDO

 END SUBROUTINE l2evap
!-------------------------------------------------------------------------------
 SUBROUTINE i2evap(D, pmx, ntrac, js, lproma, kspec)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: pmx(0:ntrac,lproma)
   REAL(dp), INTENT(IN)    :: d(0:kspec, lproma)
   INTEGER :: idx1, idx2, jt, jl

   DO jt=1,ispec_ice
     idx1 = KPP_I_IDX(js)%ice_spec(jt,ice_idx)
     idx2 = KPP_I_IDX(js)%ice_spec(jt,ice2evap)
     DO jl=1,lproma
       pmx(idx2,jl) = pmx(idx2,jl) + D(idx1,jl)
     ENDDO
   ENDDO

 END SUBROUTINE i2evap
!-------------------------------------------------------------------------------
 SUBROUTINE l2trac(D, pmx, ntrac, js, lproma, kspec)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: pmx(0:ntrac,lproma)
   REAL(dp), INTENT(INOUT) :: d(0:kspec, lproma)
   INTEGER :: idx1, idx2, jt, jl

   DO jt=1,lspec_liq
     idx1 = KPP_L_IDX(js)%liq_spec(jt,liq_idx)
     idx2 = KPP_L_IDX(js)%liq_spec(jt,liq2trac)
     IF (IDX2 == 0) CYCLE
     DO jl=1,lproma
       pmx(idx2,jl) = pmx(idx2,jl) + D(idx1,jl)
       D(idx1,jl) = 0._dp
     ENDDO
   ENDDO

 END SUBROUTINE l2trac
!-------------------------------------------------------------------------------
SUBROUTINE trac2l(D, pmx, ntrac, js, lproma, kspec)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: pmx(lproma,ntrac)
   REAL(dp), INTENT(INOUT) :: d(lproma,1:kspec)
   INTEGER :: idx1, idx2, jt, jl

   DO jt=1,lspec_liq
     idx1 = KPP_L_IDX(js)%liq_spec(jt,liq_idx)
     idx2 = KPP_L_IDX(js)%liq_spec(jt,liq2trac)
     IF (IDX2 == 0) CYCLE
     DO jl=1,lproma
       D(jl,idx1) = D(jl,idx1) + pmx(jl,idx2)
       pmx(jl,idx2) = 0._dp
     ENDDO
   ENDDO

 END SUBROUTINE trac2l
!-------------------------------------------------------------------------------
 SUBROUTINE i2trac(D, pmx, ntrac, js, lproma, kspec)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: pmx(0:ntrac,lproma)
   REAL(dp), INTENT(INOUT) :: d(0:kspec, lproma)
   INTEGER :: idx1, idx2, jt, jl

   DO jt=1,ispec_ice
     idx1 = KPP_I_IDX(js)%ice_spec(jt,ice_idx)
     idx2 = KPP_I_IDX(js)%ice_spec(jt,ice2trac)
     IF (IDX2 == 0) CYCLE
     DO jl=1,lproma
       pmx(idx2,jl) = pmx(idx2,jl) + D(idx1,jl)
       D(idx1,jl) = 0._dp
     ENDDO
   ENDDO

 END SUBROUTINE i2trac
!-------------------------------------------------------------------------------
SUBROUTINE trac2i(D, pmx, ntrac, js, lproma, kspec)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: ntrac, js, lproma, kspec
   REAL(dp), INTENT(INOUT) :: pmx(lproma,ntrac)
   REAL(dp), INTENT(INOUT) :: d(lproma,1:kspec)
   INTEGER :: idx1, idx2, jt, jl

   DO jt=1,ispec_ice
     idx1 = KPP_I_IDX(js)%ice_spec(jt,ice_idx)
     idx2 = KPP_I_IDX(js)%ice_spec(jt,ice2trac)
     IF (IDX2 == 0) CYCLE
     DO jl=1,lproma
       D(jl,idx1) = D(jl,idx1) + pmx(jl,idx2)
       pmx(jl,idx2) = 0._dp
     ENDDO
   ENDDO

 END SUBROUTINE trac2i
!-------------------------------------------------------------------------------
SUBROUTINE cloudtrac2aer(aerflux, pmx, js, lproma, ntrac, phase)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: js, lproma, ntrac
   REAL(dp), INTENT(INOUT) :: aerflux(lproma, aer_count(js)%numbers(c_all))
   REAL(dp), INTENT(INOUT) :: pmx(lproma,ntrac)
   INTEGER :: jt, idx1, jl, phase, idx

   SELECT CASE (PHASE)
   CASE(1)
     idx = aer2tracl
   CASE(2)
     idx = aer2traci
   END SELECT
   DO jt = 1,aer_count(js)%numbers(c_all)
     idx1 = aer_attr(js)%aer_spec(jt,idx)
     IF (IDX1 == 0) CYCLE
     DO jl=1,lproma
       aerflux(jl,jt) = aerflux(jl,jt) + pmx(jl,idx1)
       pmx(jl,idx1) = 0._dp
     ENDDO
   ENDDO
      
 END SUBROUTINE cloudtrac2aer
!-------------------------------------------------------------------------------
SUBROUTINE aer2cloudtrac(aerflux, pmx, js, lproma, ntrac, phase)
   
   IMPLICIT NONE
   INTEGER,  INTENT(IN)    :: js, lproma, ntrac, phase
   REAL(dp), INTENT(INOUT) :: aerflux(aer_count(js)%numbers(c_all),lproma)
   REAL(dp), INTENT(INOUT) :: pmx(ntrac,lproma)
   INTEGER :: jt, idx1, jl, idx

   REAL(dp) :: lpmx(0:ntrac,lproma)

   lpmx(0:ntrac,1:lproma) = 0._dp

   SELECT CASE (PHASE)
   CASE(1)
     idx = aer2tracl
   CASE(2)
     idx = aer2traci
   END SELECT
   DO jt = 1,aer_count(js)%numbers(c_all)
     idx1 = aer_attr(js)%aer_spec(jt,idx)
     IF (idx1 == 0 ) CYCLE
     DO jl=1,lproma
       lpmx(idx1,jl) = lpmx(idx1,jl) + aerflux(jt,jl)
       aerflux(jt,jl) = 0._dp
     ENDDO
   ENDDO
   DO jt=1,ntrac
     pmx(jt,1:lproma) = pmx(jt,1:lproma) + lpmx(jt,1:lproma)
   ENDDO


 END SUBROUTINE aer2cloudtrac
!-------------------------------------------------------------------------------
 SUBROUTINE pack_aerosol(kspec, lproma, zx, pmx, ntrac, TYPE, set, model) 

   USE messy_scav_aer,   ONLY: aermodel, spec_idx
   
   IMPLICIT NONE
   INTEGER, INTENT(IN)     :: ntrac, TYPE, set, kspec, lproma, model
   REAL(dp), INTENT(IN)    :: pmx(ntrac, lproma)
   REAL(dp), INTENT(INOUT) :: zx(kspec, lproma)
   INTEGER :: jt, idx, jl
   INTEGER, DIMENSION(:), POINTER :: index_type
   
   zx(:,:) = 0._dp
   SELECT CASE (TYPE)
     CASE (1)
       index_type => aermodel(set)%aer_mass_attr(model)%mass_int_att(:,spec_idx)
     CASE (2)
       index_type => aermodel(set)%aer_num_attr(model)%num_int_att(:,spec_idx)
   END SELECT

   DO jt = 1, kspec
     idx = index_type(jt)
     DO jl=1,lproma
       zx(jt,jl) = pmx(idx,jl)
     ENDDO
   ENDDO
   
 END SUBROUTINE pack_aerosol

!-------------------------------------------------------------------------------

 SUBROUTINE unpack_aerosol(kspec, lproma, zx, pmx, ntrac, TYPE, set, model) 

   USE messy_scav_aer,   ONLY: aermodel, spec_idx

   IMPLICIT NONE
   INTEGER, INTENT(IN)     :: ntrac, TYPE, set, kspec, lproma, model
   REAL(dp), INTENT(INOUT) :: pmx(ntrac, lproma)
   REAL(dp), INTENT(IN)    :: zx(kspec, lproma)
   INTEGER :: jt, idx, jl
   INTEGER, DIMENSION(:), POINTER :: index_type
   
   SELECT CASE (TYPE)
     CASE (1)
       index_type => aermodel(set)%aer_mass_attr(model)%mass_int_att(:,spec_idx)
     CASE (2)
       index_type => aermodel(set)%aer_num_attr(model)%num_int_att(:,spec_idx)
   END SELECT
  
   DO jt = 1, kspec
     idx = index_type(jt)
     DO jl=1,lproma
       pmx(idx,jl) = zx(jt,jl) 
     ENDDO
   ENDDO
   
 END SUBROUTINE unpack_aerosol
 
!===============================================================================

 SUBROUTINE AER_SCAV_DRV(lproma, ntrac, jk, js, itype, &
                         phase,  lwork,  time_step_len,      &
                         xt_vec, lnewflx,                    &
                         llwc,   liwc,    lcov,    lrdrad,   &
                         ltemp,  lpress,  lgrvol,  lrho,     &
                         lprec,  lsnow )

   USE messy_scav_aer,           ONLY: aermodel, max_mode, &
                                       aer_scav2, aer_scav_cloud
   USE messy_main_constants_mem, ONLY: avo => N_A

   INTEGER,  INTENT(IN)    :: lproma, ntrac, jk, js, itype, phase
   INTEGER,  INTENT(IN)    :: lwork(lproma)
   REAL(dp), INTENT(IN)    :: time_step_len
   REAL(dp), INTENT(IN)    :: llwc(lproma),  liwc(lproma)
   REAL(dp), INTENT(IN)    :: lcov(lproma),  lrdrad(lproma)
   REAL(dp), INTENT(IN)    :: lprec(lproma), lsnow(lproma)
   REAL(dp), INTENT(IN)    :: ltemp(lproma), lpress(lproma)
   REAL(dp), INTENT(IN)    :: lgrvol(lproma), lrho(lproma)
   REAL(dp), INTENT(INOUT) :: xt_vec(ntrac,lproma)
   REAL(dp), INTENT(INOUT) :: lnewflx(aer_count(js)%numbers(c_all),lproma)
   
   INTEGER  :: naermo, nmode, num
   INTEGER  :: jl, l, jt
   REAL(dp) :: l_xt(aer_count(js)%numbers(c_all),lproma)
   REAL(dp) :: loc_aer_rad(max_mode, lproma)
   REAL(dp) :: l_aer_rad(max_mode, lproma)
   REAL(dp) :: l_philfrac(max_mode, lproma)
   REAL(dp) :: aer_rad(max_mode)
   REAL(dp) :: lfrac(max_mode,lproma)

   aer_rad(:) = 2.0e-6_dp
   lfrac(:,:) = 0._dp
   l_philfrac(:,:) = 1.0_dp

   do naermo = 1,aermodel(js)%aeromodnum
     nmode = aermodel(js)%aer_input(naermo)%lmode
     num   = aermodel(js)%aer_mass_attr(naermo)%number
!    scavenging of aerosol masses
     l_xt(:,1:lproma) = 0._dp
     CALL pack_aerosol(num, lproma, &
       l_xt(1:num,1:lproma),        & ! zx
       xt_vec(:,1:lproma),          & ! pmx
       ntrac, 1, js, naermo)
     if (associated(aermodel(js)%aer_input(naermo)%wetrad)) then
!CDIR NOLOOPCHG
       do l=1,nmode
         do jl=1,lproma
           loc_aer_rad(l,jl)  = &
             aermodel(js)%aer_input(naermo)%wetrad(lwork(jl),l,jk)
         enddo
       ENDDO
     else
!CDIR NOLOOPCHG
       do l=1,nmode
         do jl=1,lproma
           loc_aer_rad(l,jl) = aer_rad(l)
         enddo
       enddo
     endif
!CDIR NOLOOPCHG
     do l=1,nmode
       do jl=1,lproma
         l_aer_rad(l,jl) = loc_aer_rad(l,jl) * &
           aermodel(js)%aer_input(naermo)%cmr2mmr(l)
       enddo
     enddo

     lfrac(:,:) = 0._dp
     aermodel(js)%aer_input(naermo)%lcalc_nucl_aer = .TRUE.
     IF (ASSOCIATED(aermodel(js)%aer_input(naermo)%mfrac)) THEN 
       aermodel(js)%aer_input(naermo)%lcalc_nucl_aer = .FALSE.
       DO l=1,SIZE(aermodel(js)%aer_input(naermo)%mfrac,2)
         DO jl=1,lproma
           lfrac(l,jl) = aermodel(js)%aer_input(naermo)%mfrac(lwork(jl),l,jk)
         END DO
       END DO
     END IF
!     print*, "mfrac ", minval(lfrac(:,:)), MAXVAL(lfrac(:,:)), &
!      SIZE(aermodel(js)%aer_input(naermo)%mfrac,2)

     SELECT CASE (itype)
     CASE(1) ! impaction scavenging
       call aer_scav2(num, lproma, js,                      &
         nmode, naermo, 1,                                 &
         l_xt(1:num,1:lproma),                             &
         lnewflx(1:aer_count(js)%numbers(c_all),1:lproma), &
!           laeriflx(1:aer_count(js)%numbers(c_all),1:lproma),&
         time_step_len, lprec(1:lproma), lcov(1:lproma),   &
         lrdrad(1:lproma), l_aer_rad(1:nmode,1:lproma),    &
         ltemp(1:lproma), lpress(1:lproma),                &
         lrho(1:lproma), lsnow(1:lproma),                  &
         lgrvol(1:lproma), phase)
     CASE(2) ! nucleation scavenging
       CALL aer_scav_cloud(num, lproma,                      &
         js, nmode, naermo, 1, phase,                        &
         l_xt(1:num,1:lproma),                               &
         lnewflx(1:aer_count(js)%numbers(c_all),1:lproma),   &
         time_step_len, llwc(1:lproma), liwc(1:lproma),      &
         lcov(1:lproma), lrdrad(1:lproma),                   &
         l_aer_rad(1:nmode,1:lproma), ltemp(1:lproma),       &
         lpress(1:lproma), lgrvol(1:lproma),                 &
         lfrac(1:nmode,1:lproma), l_philfrac(1:nmode,1:lproma) ) 
     END SELECT

     CALL unpack_aerosol(num, lproma, &
       l_xt(1:num,1:lproma),          &
       xt_vec(:,1:lproma),            &
       ntrac, 1, js, naermo)
     
!    scavenging of aerosol numbers
     l_xt(:,1:lproma) = 0._dp
     num   = aermodel(js)%aer_num_attr(naermo)%number
     CALL pack_aerosol(num, lproma, &
       l_xt(1:num,1:lproma),        & ! zx
       xt_vec(:,1:lproma),          & ! pmx
       ntrac, 2, js, naermo)
!CDIR NOLOOPCHG
     DO jt=1,num
       DO jl=1,lproma
         l_xt(jt,jl) = l_xt(jt,jl) / avo
       enddo
     enddo
!CDIR NOLOOPCHG
     do l=1,nmode
       do jl=1,lproma
         l_aer_rad(l,jl) = loc_aer_rad(l,jl)
       end do
     enddo

     aermodel(js)%aer_input(naermo)%lcalc_nucl_aer = .TRUE.
     IF (ASSOCIATED(aermodel(js)%aer_input(naermo)%nfrac)) THEN 
       aermodel(js)%aer_input(naermo)%lcalc_nucl_aer = .FALSE.
       lfrac(:,:) = 0._dp
       DO l=1,SIZE(aermodel(js)%aer_input(naermo)%nfrac,2)
         DO jl=1,lproma
           lfrac(l,jl) = aermodel(js)%aer_input(naermo)%nfrac(lwork(jl),l,jk)
         END DO
       END DO
     END IF

!     print*, "nfrac ", minval(lfrac(:,:)), MAXVAL(lfrac(:,:))
     SELECT CASE (itype)
     CASE(1) ! impaction scavenging
       call aer_scav2(num, lproma, js,                      &
         nmode, naermo, 2,                                 &
         l_xt(1:num,1:lproma),                             &
!           laeriflx(1:aer_count(js)%numbers(c_all),1:lproma),&
         lnewflx(1:aer_count(js)%numbers(c_all),1:lproma), &
         time_step_len, lprec(1:lproma), lcov(1:lproma),   &
         lrdrad(1:lproma), l_aer_rad(1:nmode,1:lproma),    &
         ltemp(1:lproma), lpress(1:lproma),                &
         lrho(1:lproma), lsnow(1:lproma),                  &
         lgrvol(1:lproma), phase)
     CASE(2) ! nucleation scavenging
       if (associated(aermodel(js)%aer_input(naermo)%philfrac)) then
!CDIR NOLOOPCHG
         do l=1,nmode
           do jl=1,lproma
             l_philfrac(l,jl)  = &
                  aermodel(js)%aer_input(naermo)%philfrac(lwork(jl),l,jk)
           enddo
         ENDDO
       endif
       CALL aer_scav_cloud(num, lproma,                      &
         js, nmode, naermo, 2, phase,                        &
         l_xt(1:num,1:lproma),                               &
         lnewflx(1:aer_count(js)%numbers(c_all),1:lproma),   &
         time_step_len, llwc(1:lproma), liwc(1:lproma),      &
         lcov(1:lproma), lrdrad(1:lproma),                   &
         l_aer_rad(1:nmode,1:lproma), ltemp(1:lproma),       &
         lpress(1:lproma), lgrvol(1:lproma),                 &
         lfrac(1:nmode,1:lproma),l_philfrac(1:nmode,1:lproma) )
     END SELECT

!CDIR NOLOOPCHG
     DO jt=1,num
       DO jl=1,lproma
         l_xt(jt,jl) = l_xt(jt,jl) * avo
       enddo
     enddo

     CALL unpack_aerosol(num, lproma, &
       l_xt(1:num,1:lproma),          &
       xt_vec(:,1:lproma),            &
       ntrac, 2, js, naermo)
   enddo
   
 END SUBROUTINE AER_SCAV_DRV

!===============================================================================


!===============================================================================

! Calculates the factor (frac_evapo) that has to be applied (evapo_aer) to the
! tracer concentrations for each grid cell (1..lproma), each tracer (1..ntrac),
! and each residual mode (1..n_resmod) to assign the cloud or precipitation
! residuals to the appropriate tracers.
  SUBROUTINE CALC_EVAPO_FRAC(status, aer_field, kpp_field, lproma, kspec, &
                             phase, js, ntrac, frac_evapo, resimp)

    USE messy_scav_aer
    USE messy_main_constants_mem,   ONLY: avo => N_A, pi

    INTRINSIC :: ALLOCATED, TRIM, SUM, MAX, MIN, REAL

    INTEGER,  INTENT(IN)    :: lproma, js, kspec, phase &
                               , ntrac
    INTEGER,  INTENT(OUT)   :: status
    REAL(dp), INTENT(IN)    :: aer_field(aer_count(js)%numbers(c_all),lproma)
    REAL(dp), INTENT(IN)    :: kpp_field(kspec,lproma)
    REAL(dp), INTENT(INOUT) :: frac_evapo(lproma,ntrac,n_resmod)
    REAL(dp), INTENT(IN), OPTIONAL :: &
         resimp(aer_count(js)%numbers(c_all),lproma)

    INTEGER                 :: idx, jl, jt, ji
    INTEGER                 :: jg, idx_mode, idx_core, idx_trac, idx_aero &
                               , nspec, ntrac_made3, j, i_made3
    INTEGER                 :: idx_made3(ntrac)
    INTEGER,  ALLOCATABLE   :: idx_aer2cld(:), idx_incl(:)
    REAL(dp), ALLOCATABLE   :: molm(:)
    REAL(dp)                :: nucnum(n_resmod,lproma)
    REAL(dp)                :: resnuc(aer_count(js)%numbers(c_all),lproma)
    REAL(dp)                :: f(n_resmod)
    REAL(dp)                :: gam, gamma, nuc
    REAL(dp)                :: Nf_as_am, Nf_as_cs, Nf_as_cm, Nf_asp2_cm &
                               , Nf_am_cm, Nf_cs_cm
    REAL(dp)                :: sum_spec
    REAL(dp)                :: zhelp2
    REAL(dp), ALLOCATABLE   :: impnum(:,:), chem(:,:)
    REAL(dp), ALLOCATABLE   :: g(:)
    CHARACTER(LEN=5)        :: str_h2o, str_so4mm, str_hso4m
    REAL(dp)                :: scavenged_num(lproma), scavenged_mass(lproma)
    REAL(dp), PARAMETER     :: num_max=1._dp
    REAL(dp)                :: zhelp

    frac_evapo(:,:,:)   = 0._dp
    IF (ALLOCATED(made3)) THEN
      ALLOCATE(impnum(made3(1)%nmod,lproma))
      ALLOCATE(g(made3(1)%nmod))
      nucnum(:,:)       = 0._dp
      impnum(:,:)       = 0._dp
      ntrac_made3       = 0
      idx_made3(:)      = 0
      resnuc(:,:)       = aer_field(:,:)
    ELSE
      scavenged_num(:)  = 0._dp
      scavenged_mass(:) = 0._dp
    END IF
    i_made3             = 0
    status              = 19 ! unknown error

    ! Store number concentrations of scavenged particles
    DO ji = 1,aermodel(js)%aeromodnum
      IF ((TRIM(aermodel(js)%aermodname(ji)) == 'made3') .OR. &
             (TRIM(aermodel(js)%aermodname(ji)) == 'MADE3')) THEN
        i_made3 = ji
      END IF
      DO jt = nucl_mode + 1,aermodel(js)%aer_num_attr(ji)%number
        idx = aermodel(js)%aer_num_attr(ji)%num_int_att(jt,spec_idx_comp)
        IF (ji .EQ. i_made3) THEN
          ! Determine mode and core mode (i.e. mode to which this tracer is
          ! assigned upon nucleation) for current tracer
          idx_mode = aermodel(js)%aer_num_attr(ji)%num_int_att(jt,spec_mode)
          idx_core = get_core_mode_made3(idx_mode)
          IF (idx_core < 0) THEN
            status = 11
            RETURN
          END IF
          ! Store number of nuclei in core mode and number of impacted
          ! particles from original mode
          DO jl = 1, lproma
            ! Exclude N_ks from the sum
            IF (idx_mode .NE. made3(1)%ks) &
                 nucnum(idx_core,jl) = nucnum(idx_core,jl) + resnuc(idx,jl)
            IF (PRESENT(resimp)) &
                 impnum(idx_mode,jl) = impnum(idx_mode,jl) + resimp(idx,jl)
          END DO
        ELSE
          IF (.NOT.aermodel(js)%aer_num_attr(ji)%num_log_att(jt,spec_sol)) CYCLE
        DO jl = 1,lproma
          scavenged_num(jl) = scavenged_num(jl) + aer_field(idx,jl) 
        END DO
        END IF
      END DO
   ENDDO

   ! from molecules per box to g / cm??
   ! NOTE: mw is 0 for number tracers
   DO jt= 1,aer_count(js)%numbers(c_all)
     IF (ALLOCATED(made3)) THEN
       idx_trac = aer_attr(js)%aer_spec(jt,aer_idx)
       ! Store number and indices of MADE3 tracers
       IF (attr(js)%int_att(idx_trac,aer_mod_num) .EQ. i_made3) THEN
         ntrac_made3 = ntrac_made3 + 1
         idx_made3(ntrac_made3) = idx_trac
       END IF
       ! Determine core mode for current tracer (i.e. mode to which this
       ! tracer is assigned upon nucleation)
       idx_core = get_core_mode_made3(attr(js)%int_att(idx_trac,tmode))
       IF (idx_core < 0) THEN
         status = 12
         RETURN
       END IF
     ELSE
     DO jl = 1,lproma
       scavenged_mass(jl) = scavenged_mass(jl) + aer_attr(js)%mw(jt) * &
            aer_field(jt,jl) / avo  
     ENDDO
     END IF
   END DO

   IF (PHASE == 2) THEN
      DO jt = 1,ispec_ice
         idx = kpp_i_idx(js)%ice_spec(jt,ice_idx)
         if (TRIM(str_field_kpp_i(idx)) == "H2O_i")  CYCLE
         DO jl=1,lproma
            scavenged_mass(jl) = scavenged_mass(jl) + &
                 kpp_i_idx(js)%ice_attr(jt,ice_mw)  * &
                 kpp_field(idx,jl) / avo
         ENDDO   
      END DO
   ELSEIF (PHASE == 1) THEN
      DO jt = 1,lspec_liq
         idx = kpp_l_idx(js)%liq_spec(jt,liq_idx)
         if (TRIM(str_field_kpp_l(idx)) == "H2O_l")  CYCLE
         DO jl=1,lproma
            scavenged_mass(jl) = scavenged_mass(jl) + &
                 kpp_l_idx(js)%liq_attr(jt,liq_mw)  * &
                 kpp_field(idx,jl) / avo
         ENDDO
      END DO
   ENDIF

   IF (n_resmod .EQ. 2) THEN

   DO jl=1,lproma

!      frac_evapo(jl) = 1._dp - MAX(0._dp, MIN(1._dp, &
!                       scavenged_num(jl) / num_max) )
!      frac_evapo(jl) = scavenged_num(jl) 
     zhelp = 0._dp

     IF (scavenged_num(jl) > 0._dp ) THEN
! g/cm?? / #/cm?? -> g
       zhelp = scavenged_mass(jl) / scavenged_num(jl) 
! divided by an approximated aerosol density of 1500 kg/m?? = 1.5e6 g/m??
! -> volume m??
       If (zhelp > 0._dp) THEN
         zhelp = zhelp / 1.5e6_dp
! -> divide by 4/3 pi and take the 3rd root
! -> radius in m
         zhelp = (zhelp * 0.75_dp / pi)**(1._dp/3._dp)

! -> from radius of average mass to Count Median Radius
!    assumed sigma of 1.59 (minimum of all existing submodels,
!    conservative approach)
!         zhelp = zhelp / EXP(1.5_dp*(LOG(1.59))**2)
         zhelp = zhelp / EXP(1.5_dp*(LOG(2.0_dp))**2._dp)
       ENDIF
         IF (zhelp < 3.e-7_dp) THEN
           frac_evapo(jl,:,2) = 1.e-9_dp
         ELSEIF ( (zhelp >= 3.e-7_dp) .AND. &
                  (zhelp <  5.e-7_dp) ) THEN
           frac_evapo(jl,:,2) = 1._dp - ( (5.e-7_dp - zhelp) / 2.e-7_dp)**4._dp
         ELSE
           frac_evapo(jl,:,2) = 1._dp
         ENDIF

       ELSE 
         frac_evapo(jl,:,2) = 1._dp
     ENDIF

   END DO

   frac_evapo(:,:,1)   = 1._dp - frac_evapo(:,:,2)

   END IF

   use_made3: IF (ALLOCATED(made3)) THEN

     ! The terminology in this block closely follows Kaiser et al., to be
     ! submitted to GMD, 2016, which is also the reference for the idea behind
     ! the algorithm. Note, however, that ice and liquid scavenging are treated
     ! by separate calls of the present subroutine.

     ! For in-cloud/in-precip generated mass calculations
     ALLOCATE(idx_aer2cld(aer_count(js)%numbers(c_all)))
     IF (PHASE .EQ. 2) THEN
       nspec             = ispec_ice
       ALLOCATE(idx_incl(nspec))
       ALLOCATE(molm(nspec))
       idx_aer2cld(1:aer_count(js)%numbers(c_all)) = &
            aer_attr(js)%aer_spec(1:aer_count(js)%numbers(c_all),aer2i_idx)
       idx_incl(1:nspec) = kpp_i_idx(js)%ice_spec(1:nspec,ice_idx)
       str_h2o           = "H2O_i"
       str_so4mm         = "SO4mm_i"
       str_hso4m         = "HSO4m_i"
       molm(1:nspec)     = kpp_i_idx(js)%ice_attr(1:nspec,ice_mw)
     ELSE IF (PHASE .EQ. 1) THEN
       nspec             = lspec_liq
       ALLOCATE(idx_incl(nspec))
       ALLOCATE(molm(nspec))
       idx_aer2cld(1:aer_count(js)%numbers(c_all)) = &
            aer_attr(js)%aer_spec(1:aer_count(js)%numbers(c_all),aer2l_idx)
       idx_incl(1:nspec) = kpp_l_idx(js)%liq_spec(1:nspec,liq_idx)
       str_h2o           = "H2O_l"
       str_so4mm         = "SO4mm_l"
       str_hso4m         = "HSO4m_l"
       molm(1:nspec)     = kpp_l_idx(js)%liq_attr(1:nspec,liq_mw)
     END IF
! FIXME: It would be more elegant to do the following via
! something like the liq2evap indices that are used in the default case (to add
! the KPP species to the residuals via `evapo'), but this should do the trick in
! our case as only one aerosol species corresponds to two KPP species.
     ALLOCATE(chem(kspec,lproma))
     chem(:,:) = kpp_field(:,:)
     DO jg = 1, nspec
       IF (TRIM(str_field_kpp_l(idx_incl(jg))) .EQ. str_so4mm) THEN
         DO j = 1, nspec
           IF (TRIM(str_field_kpp_l(idx_incl(j))) .EQ. str_hso4m) THEN
             chem(idx_incl(jg),:) = chem(idx_incl(jg),:) &
                  + chem(idx_incl(j),:)
             EXIT
           END IF
         END DO
         EXIT
       END IF
     END DO

     loop_cells: DO jl = 1, lproma

       ! Set coefficients
       zhelp = SUM(nucnum(:,jl))
       IF (zhelp .GT. 1.0e-30_dp) THEN
         zhelp = 1.0_dp / zhelp
       ELSE
         zhelp = 0.0_dp
       END IF
       f     = nucnum(:,jl) * zhelp
       g     = impnum(:,jl) * zhelp
       gam   = MIN(g(made3(1)%km)+g(made3(1)%ki)+g(made3(1)%am)+g(made3(1)%ai) &
                   , 1._dp)
       gamma = MIN(g(made3(1)%cm)+g(made3(1)%ci), 1._dp)

       ! Compute transferred number fractions
       Nf_as_am   = MAX(MIN(gam, 1._dp-gamma)-g(made3(1)%cs), 0._dp)
       Nf_as_cs   = MAX(MIN(g(made3(1)%cs), 1._dp-gamma)-gam, 0._dp)
       Nf_as_cm   = gamma
       Nf_asp2_cm = MIN(MIN(gam, g(made3(1)%cs)), 1._dp-gamma)
       IF (nucnum(2,jl) .GT. 1.e-30_dp) THEN
         Nf_am_cm = MIN(f(2)*(impnum(&
              made3(1)%cs,jl)+impnum(made3(1)%cm,jl)+impnum(made3(1)%ci,jl))&
              /nucnum(2,jl), 1._dp)
       ELSE
         Nf_am_cm = 0._dp
       END IF
       IF (nucnum(3,jl) .GT. 1.e-30_dp) THEN
         Nf_cs_cm = MIN(f(3)*(impnum(&
              made3(1)%km,jl)+impnum(made3(1)%ki,jl)+impnum(made3(1)%am,jl)+&
              impnum(&
              made3(1)%ai,jl)+impnum(made3(1)%cm,jl)+impnum(made3(1)%ci,jl))&
              /nucnum(3,jl), 1._dp)
       ELSE
         Nf_cs_cm = 0._dp
       END IF

       ! Loop over wet deposited MADE3 aerosol tracers
       loop_tracers: DO jt = 1, ntrac_made3

         idx_trac = idx_made3(jt)
         idx_aero = attr(js)%int_att(idx_trac,aerosol_index)

         ! Determine core mode for current tracer (i.e. mode to which this
         ! tracer is assigned upon nucleation)
         idx_core = get_core_mode_made3(attr(js)%int_att(idx_trac,tmode))
         IF (idx_core < 0) THEN
           status = 12
           RETURN
         END IF

         ! Set variables `nuc' and `zhelp' depending on whether current tracer
         ! is a mass or a number tracer
         nuc = resnuc(idx_aero,jl)
         IF (attr(js)%int_att(idx_trac,lquantity) == 1) THEN ! mass tracer
           nuc   = nuc * aer_attr(js)%mw(idx_aero) / avo
           zhelp = nuc
           IF (idx_aer2cld(idx_aero) .NE. 0) THEN ! KPP
             ! Add in-cloud/in-precip generated (or subtract lost) mass
             DO jg = 1, nspec
               IF ((idx_incl(jg) .EQ. idx_aer2cld(idx_aero)) &
                    .AND. (TRIM(str_field_kpp_l(idx_incl(jg))) /= str_h2o)) &
                    THEN
                 IF (chem(idx_incl(jg),jl) .GE. 0._dp) THEN
                   ! Generated mass is distributed equally across residual modes
                   IF ((attr(js)%int_att(idx_trac,tmode) == idx_core) &
                        ) THEN ! residual mode
                     nuc = nuc + molm(jg) * chem(idx_incl(jg),jl) &
                          / (avo * REAL(n_resmod, kind=dp))
                   END IF
                 ELSE
                   ! Calculate sum of concentrations of MADE3 tracers associated
                   ! with currently processed KPP species
                   sum_spec = 0._dp
                   DO j = 1, ntrac_made3
                     idx = attr(js)%int_att(idx_made3(j),aerosol_index)
                     IF (idx_aer2cld(idx) .EQ. idx_incl(jg)) THEN
                       sum_spec = sum_spec &
                            + resnuc(idx,jl) * aer_attr(js)%mw(idx) / avo
                       IF (PRESENT(resimp)) THEN
                         sum_spec = sum_spec &
                              + resimp(idx,jl) * aer_attr(js)%mw(idx) / avo
                       END IF
                     END IF
                   END DO
                   ! Lost mass is taken proportionally from each tracer
                   ! Note: `nuc' can become negative here (but in that case
                   ! should become positive again below)!
                   IF (sum_spec .GT. 1.e-30_dp) THEN
                     IF (PRESENT(resimp)) THEN
                        nuc = nuc + (nuc + resimp(idx_aero,jl) &
                            * aer_attr(js)%mw(idx_aero) / avo) &
                            * molm(jg) * chem(idx_incl(jg),jl) &
                            / (avo * sum_spec)
                     ELSE
                       nuc = nuc &
                            * (1._dp + molm(jg) * chem(idx_incl(jg),jl)& 
                            / (avo * sum_spec))
                     END IF
                   ELSE
                     ! No aerosol corresponding to currently processed KPP
                     ! species available
                     nuc = 0._dp
                   END IF
                 END IF
               END IF
             END DO
           END IF
         ELSE IF (attr(js)%int_att(idx_trac,lquantity) == 2) THEN
           zhelp = nuc
         END IF
         IF (PRESENT(resimp)) THEN
           IF (attr(js)%int_att(idx_trac,lquantity) == 1) THEN  ! mass tracer
             zhelp = zhelp &
                  + resimp(idx_aero,jl) * aer_attr(js)%mw(idx_aero) / avo
             ! Add impacted mass that stays within the mode of the core
             nuc = nuc + f(idx_core) * resimp(idx_aero,jl) &
                         * aer_attr(js)%mw(idx_aero) / avo
           ELSE IF (attr(js)%int_att(idx_trac,lquantity) == 2) THEN
             zhelp = zhelp + resimp(idx_aero,jl)
           END IF
         END IF
         IF (zhelp .GT. 1.e-30_dp) THEN
           zhelp = 1._dp / zhelp
         ELSE
           zhelp = 0._dp
         END IF

         ! Perform actual mode assignment, i.e. set frac_evapo.

         ! The idea behind the implementation was to first consider one residual
         ! mode after the other, except `idx_core', and perform the following
         ! three steps for each of them:
         ! 1. Calculate transfers from `idx_core' to the mode under
         !    consideration.
         ! 2. (mass only) Calculate the fraction of the impacted mass that
         !    remains in the mode under consideration.
         ! 3. (mass only) Calculate the transfer from the other modes to the
         !    mode under consideration.
         ! Subsequently, the remaining fraction is assigned to `idx_core'.

         ! ... code that is used for both number and mass tracers
         frac_evapo(jl,idx_trac,:) = nuc * zhelp
         SELECT CASE (idx_core)
         CASE (1) ! as
           frac_evapo(jl,idx_trac,2) = Nf_as_am * frac_evapo(jl,idx_trac,2)
           frac_evapo(jl,idx_trac,3) = Nf_as_cs * frac_evapo(jl,idx_trac,3)
           frac_evapo(jl,idx_trac,4) = (Nf_as_cm + Nf_asp2_cm) &
                * frac_evapo(jl,idx_trac,4)
           frac_evapo(jl,idx_trac,1) = frac_evapo(jl,idx_trac,1) &
                - SUM(frac_evapo(jl,idx_trac,2:4))
         CASE (2) ! am
           frac_evapo(jl,idx_trac,1) = 0._dp
           frac_evapo(jl,idx_trac,3) = 0._dp
           frac_evapo(jl,idx_trac,4) = Nf_am_cm * frac_evapo(jl,idx_trac,4)
           frac_evapo(jl,idx_trac,2) = frac_evapo(jl,idx_trac,2) &
                - frac_evapo(jl,idx_trac,4)
         CASE (3) ! cs
           frac_evapo(jl,idx_trac,1) = 0._dp
           frac_evapo(jl,idx_trac,2) = 0._dp
           frac_evapo(jl,idx_trac,4) = Nf_cs_cm * frac_evapo(jl,idx_trac,4)
           frac_evapo(jl,idx_trac,3) = frac_evapo(jl,idx_trac,3) &
                - frac_evapo(jl,idx_trac,4)
         CASE (4) ! cm
           frac_evapo(jl,idx_trac,1:3) = 0._dp
           ! frac_evapo(jl,idx_trac,4): all residuals from cloud constituents
           ! whose core belongs to mode `cm' are assigned to residual mode `cm'
         CASE DEFAULT
           status = 12
           RETURN
         END SELECT

         ! ... code that is only used for number tracers
         IF (attr(js)%int_att(idx_trac,lquantity) .EQ. 2) THEN
           ! If specified by the user (i.e. namelist parameter frac_resnum < 1),
           ! reduce number conc. of residual aerosol
           frac_evapo(jl,idx_trac,:) = frac_evapo(jl,idx_trac,:) * frac_resnum

           ! Set N_ks fraction to zero
           IF (attr(js)%int_att(idx_trac,tmode) .EQ. made3(1)%ks) &
                frac_evapo(jl,idx_trac,:) = 0._dp
         END IF

         ! ... code that is only used for mass tracers
         mass_tracer: IF (PRESENT(resimp) &
              .AND. (attr(js)%int_att(idx_trac,lquantity) .EQ. 1)) THEN
           zhelp2 = resimp(idx_aero,jl) * aer_attr(js)%mw(idx_aero) * zhelp &
                 / avo
           SELECT CASE (idx_core)
           CASE (1)
             frac_evapo(jl,idx_trac,2) = frac_evapo(jl,idx_trac,2) &
                  + f(2) * (1._dp - Nf_am_cm) * zhelp2
             frac_evapo(jl,idx_trac,3) = frac_evapo(jl,idx_trac,3) &
                  + f(3) * (1._dp - Nf_cs_cm) * zhelp2
             frac_evapo(jl,idx_trac,4) = frac_evapo(jl,idx_trac,4) &
                  + (f(4) + f(2) * Nf_am_cm + f(3) * Nf_cs_cm) * zhelp2
             frac_evapo(jl,idx_trac,1) = nuc * zhelp &
                  + (1._dp - f(1)) * zhelp2 - SUM(frac_evapo(jl,idx_trac,2:4))
           CASE (2)
             frac_evapo(jl,idx_trac,4) = frac_evapo(jl,idx_trac,4) &
                   + (f(4) + f(3) + f(1)) * zhelp2
             IF (gam .GT. 1.e-30_dp) frac_evapo(jl,idx_trac,4) = &
                  frac_evapo(jl,idx_trac,4) - f(1) * Nf_as_am * zhelp2 / gam
             frac_evapo(jl,idx_trac,2) = nuc * zhelp &
                  + (1._dp - f(2)) * zhelp2 - frac_evapo(jl,idx_trac,4)
           CASE (3)
             frac_evapo(jl,idx_trac,4) = frac_evapo(jl,idx_trac,4) &
                   + (f(4) + f(2) + f(1)) * zhelp2
             IF (g(made3(1)%cs) .GT. 1.e-30_dp) &
                  frac_evapo(jl,idx_trac,4) = frac_evapo(jl,idx_trac,4) &
                  - f(1) * Nf_as_cs * zhelp2 / g(made3(1)%cs)
             frac_evapo(jl,idx_trac,3) = nuc * zhelp &
                  + (1._dp - f(3)) * zhelp2 - frac_evapo(jl,idx_trac,4)
           CASE (4)
             frac_evapo(jl,idx_trac,4) = nuc * zhelp + (1._dp - f(4)) * zhelp2
           CASE DEFAULT
             status = 12
             RETURN
           END SELECT
         END IF mass_tracer

       END DO loop_tracers

     END DO loop_cells

     DEALLOCATE(impnum)
     DEALLOCATE(g)
     DEALLOCATE(idx_aer2cld)
     DEALLOCATE(idx_incl)
     DEALLOCATE(molm)
     DEALLOCATE(chem)

   END IF use_made3

   status = 0 ! no error

  END SUBROUTINE CALC_EVAPO_FRAC
!===============================================================================

  SUBROUTINE SCAV_MAIN(status, kproma, nlev, ntrac, jrow, js, kspec, &
                      max_lev_scav, time_step_len,          &
                      rain,   snow, rainprod,    snowprod,  &
                      lwc,    iwc,  temp,        press,     &
                      cover,  rcover, vol,       rho,       &
                      imelt,        isedi,                  &
                      bc_lwc,       gboxarea,               &
                      xt,                                   &
                      cloudl_field, cloudi_field,           &
                      cloudl_aer,   cloudi_aer,             &
                      rain_field,   snow_field,             &
                      rain_aer,     snow_aer,               &
                      rainpH,       cloudpH,                &
                      Hp_cloud,                             &
                      L_LS,   lnucl_i,  lnucl_l,            &
                      xsteps_rain,  xsteps_cloud,           &
                      xrsteps_rain, xrsteps_cloud,          &
                      prod_array,                           &
                      phi_i, mju_i,                         &
                      iwc_i, iwc_T_i )

   USE messy_scav_l_kpp,     ONLY: lspec => nspec
   USE messy_scav_i_kpp,     ONLY: ispec => nspec

   USE messy_scav,                 ONLY: scav_aqchem_kpp,    &
                                         scav_easy_liq,      &
                                         scav_easy_ice,      &
                                         calc_lwc_bc_n

   IMPLICIT NONE

   INTRINSIC :: ALLOCATED, REAL

   INTEGER,  INTENT(IN) :: js, jrow, kproma, nlev, ntrac
   INTEGER,  INTENT(IN) :: kspec, max_lev_scav
   INTEGER,  INTENT(OUT):: status
   REAL(dp), INTENT(IN) :: time_step_len

   ! liquid precipitation flux
   REAL(dp), INTENT(IN) :: rain(kproma,nlev)      
   ! frozen precipitation flux
   REAL(dp), INTENT(IN) :: snow(kproma,nlev)     
   ! liquid precipitation production 
   REAL(dp), INTENT(IN) :: rainprod(kproma,nlev)
   ! frozen precipitation production 
   REAL(dp), INTENT(IN) :: snowprod(kproma,nlev)   
   ! liquid water content
   REAL(dp), INTENT(IN) :: lwc(kproma,nlev)
   ! ice water content
   REAL(dp), INTENT(IN) :: iwc(kproma,nlev)  
   ! air temperature
   REAL(dp), INTENT(IN) :: temp(kproma,nlev)  
   ! air pressure
   REAL(dp), INTENT(IN) :: press(kproma,nlev)  
   ! cloud cover
   REAL(dp), INTENT(IN) :: cover(kproma,nlev)  
   ! precipitating cover
   REAL(dp), INTENT(IN) :: rcover(kproma,nlev)  
   ! grid volume
   REAL(dp), INTENT(IN) :: vol(kproma,nlev)  
   ! air density
   REAL(dp), INTENT(IN) :: rho(kproma,nlev)
   ! melting of ice
   REAL(dp), INTENT(IN) :: imelt(kproma,nlev)  
   ! sedimentation of ice
   REAL(dp), INTENT(IN) :: isedi(kproma,nlev)  
   ! grid box area
   REAL(dp), INTENT(IN) :: gboxarea(kproma)
   ! large_scale or convective switch
   LOGICAL,  INTENT(IN) :: L_LS

   ! "tracers"
   REAL(dp), INTENT(INOUT) :: xt(kproma,nlev,ntrac)
!  WARNING: The MADE3 scheme assumes that the following `.*_field' and `.*_aer'
!  arrays are initialized and contain only zeros upon entry to this subroutine
   ! "in cloud liquid chemistry species"
   REAL(dp), INTENT(INOUT) :: cloudl_field(kproma,lspec,nlev)
   ! "in cloud ice chemistry species"
   REAL(dp), INTENT(INOUT) :: cloudi_field(kproma,ispec,nlev)
   ! "in rain liquid chemistry species"
   REAL(dp), INTENT(INOUT) :: rain_field(kproma,lspec)
   ! "in rain ice chemistry species"
   REAL(dp), INTENT(INOUT) :: snow_field(kproma,ispec)
   ! "in liquid cloud aerosols
   REAL(dp), INTENT(INOUT) :: &
     cloudl_aer(kproma,aer_count(js)%numbers(c_all),nlev)
   ! "in ice cloud aerosols
   REAL(dp), INTENT(INOUT) :: &
     cloudi_aer(kproma,aer_count(js)%numbers(c_all),nlev)
   ! "in rain flux aerosols"
   REAL(dp), INTENT(INOUT) :: rain_aer(kproma,aer_count(js)%numbers(c_all))
   ! "in snow flux aerosols"
   REAL(dp), INTENT(INOUT) :: snow_aer(kproma,aer_count(js)%numbers(c_all))
   ! rain pH
   REAL(dp), INTENT(INOUT) :: rainpH(kproma,nlev)
   ! cloud pH
   REAL(dp), INTENT(INOUT) :: cloudpH(kproma,nlev)
   ! cloud H+ concentration
   REAL(dp), INTENT(INOUT) :: Hp_cloud(kproma,nlev)
! for output of partitioning parameters
   REAL(dp), INTENT(INOUT) :: phi_i(kproma,nlev)
   REAL(dp), INTENT(INOUT) :: mju_i(kproma,nlev)
   REAL(dp), INTENT(INOUT) :: iwc_i(kproma,nlev)
   REAL(dp), INTENT(INOUT) :: iwc_T_i(kproma,nlev)

   ! logicals, if nucleation scavenging has been active in a grid box
   ! to determine, whether the cloud tracer has to be resettet
   LOGICAL, INTENT(INOUT)  :: lnucl_i(kproma,nlev), lnucl_l(kproma,nlev)
   ! production rates
   REAL(dp), INTENT(INOUT) :: prod_array(kproma,nlev,nprod,2)
   ! number of kpp_steps
   REAL(dp), INTENT(INOUT) :: xsteps_rain(kproma,nlev), xsteps_cloud(kproma,nlev)
   REAL(dp), INTENT(INOUT) :: xrsteps_rain(kproma,nlev), &
                              xrsteps_cloud(kproma,nlev)
   REAL(dp), INTENT(INOUT) :: bc_lwc(kproma,nlev)
   

   INTEGER  :: phase ! 1=water, 2=ice
   INTEGER  :: jl, jk, jt, i
   ! mean droplet radius
   REAL(dp) :: rdrad(kproma), col_cov(kproma)  
   ! precipitation liquid water content
   REAL(dp) :: lwc_bc(kproma), zlwc(kproma), ziwc(kproma)
   REAL(dp) :: frac
   REAL(dp) :: zfac(kproma)
   INTEGER  :: idx1, idx2

! new packing for scavenging
   INTEGER  :: lproma, mproma
   REAL(dp) :: llwc(kproma), llwc2(kproma), liwc(kproma)
   REAL(dp) :: lrdrad(kproma), lpress(kproma), ltemp(kproma)
   REAL(dp) :: lprec(kproma), lsnow(kproma), lcov(kproma), lrho(kproma)
   REAL(dp) :: lgrvol(kproma), lprod(kproma), lsedi(kproma)
   REAL(dp) :: lkppflx(kspec,kproma), lkppfld(kspec,kproma)
   REAL(dp) :: lkppiflx(kspec,kproma), mkppiflx(kspec,kproma)
   REAL(dp) :: lph(kproma), lhp(kproma), lmelt(kproma), lsnow_up(kproma)
   REAL(dp) :: lsteps(kproma), lrsteps(kproma)
! for output of partitioning parameters
   REAL(dp) :: lphi_i(kproma), lmju_i(kproma)
   REAL(dp) :: liwc_i(kproma), liwc_T_i(kproma)

   REAL(dp) :: xt_vec(ntrac,kproma), mxt_vec(ntrac,kproma)
   REAL(dp) :: laerflx(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: laeriflx(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: maerflx(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: maeriflx(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: lnewflx(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: laersediflx(aer_count(js)%numbers(c_all), kproma)
   
   REAL(dp) :: mkppflx(kspec,kproma), mkppfld(kspec,kproma)
   REAL(dp) :: lkppsediflx(kspec,kproma)

   REAL(dp) :: aerfield_evapo(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: impfield_evapo(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: kppfield_evapo(kspec, kproma)
   
   REAL(dp) :: mgrvol(kproma), mmelt(kproma), msnow_up(kproma)
   REAL(dp) :: kpp_rain_help(kproma), kpp_snow_help(kproma)
   REAL(dp) :: aer_rain_help(kproma), aer_snow_help(kproma)

   INTEGER  :: lwork(kproma),mwork(kproma)
   REAL(dp) :: kpp_field_cloud2(kspec,kproma)
   REAL(dp) :: new_aerflx(kproma, aer_count(js)%numbers(c_all))

   ! sedimented ice residuals from in-cloud scavenging, aerosol
   REAL(dp) :: rescldinuc(kproma, aer_count(js)%numbers(c_all))
   ! snow residuals from in-cloud scavenging, aerosol
   REAL(dp) :: ressnownuc(kproma, aer_count(js)%numbers(c_all))
   ! rain residuals from in-cloud scavenging, aerosol
   REAL(dp) :: resrainnuc(kproma, aer_count(js)%numbers(c_all))
   ! snow residuals from impaction scavenging, aerosol
   REAL(dp) :: ressnowimp(kproma, aer_count(js)%numbers(c_all))
   ! rain residuals from impaction scavenging, aerosol
   REAL(dp) :: resrainimp(kproma, aer_count(js)%numbers(c_all))
   ! sedimented ice residuals from chemistry, KPP species
   REAL(dp) :: rescldichm(kproma, kspec)
   ! snow residuals from chemistry, KPP species
   REAL(dp) :: ressnowchm(kproma, kspec)
   ! rain residuals from chemistry, KPP species
   REAL(dp) :: resrainchm(kproma, kspec)
   ! arrays that are used for the actual computations, correspond to the above
   ! (l.* for cells with clouds, m.* for those without)
   REAL(dp) :: lrescldinuc(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: lressedinuc(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: lressnownuc(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: lrescldlnuc(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: lresrainnuc(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: lressnowimp(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: lresrainimp(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: lrescldichm(kspec, kproma)
   REAL(dp) :: lressedichm(kspec, kproma)
   REAL(dp) :: lressnowchm(kspec, kproma)
   REAL(dp) :: lrescldlchm(kspec, kproma)
   REAL(dp) :: lresrainchm(kspec, kproma)
   REAL(dp) :: mrescldinuc(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: mressnownuc(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: mresrainnuc(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: mressnowimp(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: mresrainimp(aer_count(js)%numbers(c_all), kproma)
   REAL(dp) :: mrescldichm(kspec, kproma)
   REAL(dp) :: mressnowchm(kspec, kproma)
   REAL(dp) :: mresrainchm(kspec, kproma)

   REAL(dp) :: tmpbuffer(kproma,ntrac,n_resmod)
   INTEGER  :: idt, j

! compiler workaround for g95 (0.92, 0.93) allocate from 0
   REAL(dp) :: cprod_array(kproma,0:nprod)
   REAL(dp) :: evfrac(kproma)
   REAL(dp) :: evfrac_aer(kproma,ntrac,n_resmod)

   kpp_field_cloud2(:,:) = 0._dp
   new_aerflx(:,:)       = 0._dp
   rescldinuc(:,:)       = 0._dp
   ressnownuc(:,:)       = 0._dp
   resrainnuc(:,:)       = 0._dp
   ressnowimp(:,:)       = 0._dp
   resrainimp(:,:)       = 0._dp
   rescldichm(:,:)       = 0._dp
   ressnowchm(:,:)       = 0._dp
   resrainchm(:,:)       = 0._dp
   status                = 9 ! unknown error

   phi_i(:,:)            = 0._dp
   mju_i(:,:)            = 0._dp
   iwc_i(:,:)            = 0._dp
   iwc_T_i(:,:)          = 0._dp
   col_cov(:)            = 0._dp

   level_loop: do jk=max_lev_scav,nlev

!-------------------------------------------------------------------------
!--- Nucleation scavenging --------
!---  Tansfer from cloud  ---------
!---     to rain water    ---------
      nucleation_scavenging: if (lscav_nuc) then
!       fixed cloud droplet mean radius
        rdrad(:) = 1.75e-2_dp
!-----------------------------------------
        ice_nucl_scav: if (lscav_i) THEN
! pack for snow and ice
! into a vector that contains all the boxes with 
! solid/frozen nucleation scavenging
          phase  = 2
            
          lproma = 0
          liwc(:)          = 0._dp
          llwc(:)          = 0._dp
          lrdrad(:)        = 0._dp
          lsedi(:)         = 0._dp
          lpress(:)        = 0._dp
          ltemp(:)         = 0._dp
          lcov(:)          = 0._dp
          lgrvol(:)        = 0._dp
          xt_vec(:,:)      = 0._dp
          lprod(:)         = 0._dp
          laerflx(:,:)     = 0._dp
          lkppflx(:,:)     = 0._dp
          lkppfld(:,:)     = 0._dp
          lkppsediflx(:,:) = 0._dp
          laersediflx(:,:) = 0._dp

          lrescldinuc(:,:) = 0._dp
          lressedinuc(:,:) = 0._dp
          lressnownuc(:,:) = 0._dp
          lrescldichm(:,:) = 0._dp
          lressedichm(:,:) = 0._dp
          lressnowchm(:,:) = 0._dp

! for output of partitioning parameters
          lphi_i(:)        = 0._dp
          lmju_i(:)        = 0._dp
          liwc_i(:)        = 0._dp
          liwc_T_i(:)      = 0._dp
!
          lhp(:)           = 0._dp
          lph(:)           = 0._dp
          lrho(:)          = 0._dp
          
          mproma           = 0
          mxt_vec(:,:)     = 0._dp
          maerflx(:,:)     = 0._dp

          mrescldinuc(:,:) = 0._dp
          mrescldichm(:,:) = 0._dp

          mkppflx(:,:)     = 0._dp
          mkppfld(:,:)     = 0._dp
          mgrvol(:)        = 0._dp
          cprod_array(:,:) = 0._dp
          lwork(:) = 0
          mwork(:) = 0
          kpp_snow_help(:) = 0._dp
          aer_snow_help(:) = 0._dp
          evfrac(:)        = 0._dp

          evfrac_aer(:,:,:)= 0._dp

          ! replaced maxval statement
          DO jt=1,ispec
            DO jl=1,kproma
              IF (L_LS) THEN
                kpp_field_cloud2(jt,jl) = kpp_field_cloud2(jt,jl) &
                                        + cloudi_field(jl,jt,jk)
                cloudi_field(jl,jt,jk) = 0._dp
              ENDIF
              kpp_snow_help(jl) = MAX(kpp_snow_help(jl),kpp_field_cloud2(jt,jl))
            ENDDO
          ENDDO
          DO jt=1,aer_count(js)%numbers(c_all)
            DO jl=1,kproma
              IF (L_LS) THEN
                new_aerflx(jl,jt) = new_aerflx(jl,jt) &
                                  + cloudi_aer(jl,jt,jk)
                cloudi_aer(jl,jt,jk) = 0._dp
              ENDIF
              aer_snow_help(jl) = MAX(aer_snow_help(jl), new_aerflx(jl,jt))
            ENDDO
          ENDDO
          IF (L_LS) lnucl_i(1:kproma,jk) = .TRUE.

          do jl=1,kproma
            if ( (iwc(jl,jk)   > thres_val) .AND. &
                 (cover(jl,jk) > 0._dp) ) THEN
              lproma = lproma + 1
              lwork(lproma) = jl
              lnucl_i(jl,jk)   = .TRUE.
            ELSE
              IF ( (kpp_snow_help(jl) > 0.0_dp) .or. &
                   (aer_snow_help(jl) > 0.0_dp) ) THEN
                mproma = mproma + 1
                mwork(mproma)=jl
              ENDIF
            ENDIF
          enddo
          
          IF (LPROMA > 0) THEN
            DO jl=1,lproma
              i = lwork(jl)
              liwc(jl)         = iwc(i,jk)
              lrdrad(jl)       = rdrad(i)
              lpress(jl)       = press(i,jk)
              ltemp(jl)        = temp(i,jk)
              lcov(jl)         = cover(i,jk)
              lgrvol(jl)       = vol(i,jk)
              lprod(jl)        = snowprod(i,jk)
              lrho(jl)         = rho(i,jk)
              lsedi(jl)        = isedi(i,jk)
! for output of partitioning parameters
              lphi_i(jl)       = phi_i(i,jk)
              lmju_i(jl)       = mju_i(i,jk)
              liwc_i(jl)       = iwc_i(i,jk)
              liwc_T_i(jl)     = iwc_T_i(i,jk) 
!
            ENDDO
            DO jt=1,ntrac
              DO jl=1,lproma
                i=lwork(jl)
                xt_vec(jt,jl)   = xt(i,jk,jt)
              ENDDO
            ENDDO
            DO jt=1,ispec
              DO jl=1,lproma
                i=lwork(jl)
                lkppflx(jt,jl)  = snow_field(i,jt)
                lkppfld(jt,jl)  = kpp_field_cloud2(jt,i) &
                                + cloudi_field(i,jt,jk)
                cloudi_field(i,jt,jk) = 0._dp
                lrescldichm(jt,jl) = rescldichm(i,jt)
                lressnowchm(jt,jl) = ressnowchm(i,jt)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                laerflx(jt,jl) = snow_aer(i,jt)
                lnewflx(jt,jl) = new_aerflx(i,jt) &
                               + cloudi_aer(i,jt,jk)
                cloudi_aer(i,jt,jk) = 0._dp
                lrescldinuc(jt,jl) = rescldinuc(i,jt)
                lressnownuc(jt,jl) = ressnownuc(i,jt)
              ENDDO
            ENDDO

          ENDIF  ! on lproma
          IF (mproma > 0) THEN
            DO jt=1,ntrac
              DO jl=1,mproma
                i=mwork(jl)
                mxt_vec(jt,jl)   = xt(i,jk,jt)
              ENDDO
            ENDDO
            DO jl=1,mproma
              i=mwork(jl)
              mgrvol(jl) = vol(i,jk)
            ENDDO
            DO jt=1,ispec
              DO jl=1,mproma
                i=mwork(jl)
                mkppfld(jt,jl)  = kpp_field_cloud2(jt,i)
                mrescldichm(jt,jl) = rescldichm(i,jt)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,mproma
                i=mwork(jl)
                maerflx(jt,jl) = new_aerflx(i,jt)
                mrescldinuc(jt,jl) = rescldinuc(i,jt)
              ENDDO
            ENDDO
          ENDIF ! on mproma

!-----------------------------------------------------------------------
! aerosol scavenging
          if (lscav_aer) then
            if (lproma > 0) THEN

              CALL AER_SCAV_DRV(lproma, ntrac, jk, js, 2,             &
                phase, lwork(1:lproma), time_step_len,                &
                xt_vec(:,1:lproma),                                   &
                lnewflx(1:aer_count(js)%numbers(c_all),1:lproma),     &
                llwc(1:lproma),   liwc(1:lproma),   lcov(1:lproma),   &
                lrdrad(1:lproma), ltemp(1:lproma),  lpress(1:lproma), &
                lgrvol(1:lproma), lrho(1:lproma),   lprec(1:lproma),  &
                lsnow(1:lproma)  )

              ! store new ice nuclei
              DO jt=1,aer_count(js)%numbers(c_all)
                DO jl=1,lproma
                  i=lwork(jl)
                  lrescldinuc(jt,jl) = lrescldinuc(jt,jl) + &
                    MAX(lnewflx(jt,jl) - new_aerflx(i,jt), 0._dp)
                ENDDO
              ENDDO

            END IF  ! on lproma > 0
          ENDIF  ! on lscav_aer

!-----------------------------------------------------------------------
! gas scavenging and chemistry
          IF (lscav_gas) THEN

            IF (lproma > 0) then

              CALL in_cloud(lkppfld(:,1:lproma), xt_vec(:,1:lproma), &
                            lnewflx(:,1:lproma), lgrvol(1:lproma),   &
                            lcov(1:lproma),                          &
                            ntrac, lproma, js, 2, kspec)

              ! Subtract current concentrations of KPP species in order to
              ! calculate increment (below)
              DO jt = 1, ispec
                lrescldichm(jt,1:lproma) = lrescldichm(jt,1:lproma) &
                     - I_SPEC(jt,1:lproma) * lgrvol(1:lproma) * 1.e6_dp &
                     * lcov(1:lproma)
              END DO

              if (iscav_easy) then
                call scav_easy_ice(liwc(1:lproma), ltemp(1:lproma), &
                                lpress(1:lproma),                   &
! for output of partitioning parameters
                                lphi_i(1:lproma),                   &
                                lmju_i(1:lproma), liwc_i(1:lproma), &
                                liwc_T_i(1:lproma),                 &
!
                                lproma, js)
              else
 !                 call scav_icechem_kpp(ltemp, lpress, llwc, llwc, lrdrad, &
 !                                       time_step_len, 1, lproma, lwork,   &
 !                                       jk, jrow)
                print*, "icechem"
              endif
              
              call out_cloud(lkppfld(:,1:lproma), lkppflx(:,1:lproma),      &
                             xt_vec(:,1:lproma), lph(1:lproma),             &
                             lhp(1:lproma), lcov(1:lproma), liwc(1:lproma), &
                             lprod(1:lproma), lgrvol(1:lproma),             &
                             ntrac, lproma, js, 2, kspec,                   &
                             cprod_array(1:lproma,:) )

              ! Calculate increment of KPP species concentrations due to
              ! chemistry
              DO jt = 1, ispec
                lrescldichm(jt,1:lproma) = lrescldichm(jt,1:lproma) &
                     + I_SPEC(jt,1:lproma) * lgrvol(1:lproma) * 1.e6_dp &
                     * lcov(1:lproma)
              END DO
            ENDIF
          ENDIF   ! lscav_gas

!------------------------------------------------------
!sedimentation calculation of both dissolved gases and aerosols

          IF (lproma > 0) then
    
              CALL sedi_cloud(lkppfld(:,1:lproma), lnewflx(:,1:lproma), &
                              lkppsediflx(:,1:lproma),                  &
                              laersediflx(:,1:lproma),                  &
                              lsedi(1:lproma), liwc(1:lproma),          &
                              lproma, js, kspec)

              ! store sedimenting ice residuals
              zfac(1:lproma) = &
                   MAX(MIN(lsedi(1:lproma) / liwc(1:lproma), 1.0_dp), 0._dp)
              DO jt = 1, aer_count(js)%numbers(c_all)
                lressedinuc(jt,1:lproma) = &
                     lrescldinuc(jt,1:lproma) * zfac(1:lproma)
              END DO
              DO jt = 1, ispec
                lressedichm(jt,1:lproma) = &
                     lrescldichm(jt,1:lproma) * zfac(1:lproma)
              END DO

              ! As `sedi_cloud' operates on the aerosol array `lnewflx', which
              ! still contains the residuals that are removed by precipitation,
              ! and since it uses `lsedi'/`liwc' (instead of
              ! `lsedi'/(`liwc'+`lsedi')) as the fraction of sedimenting in-ice
              ! aerosol, I assume that `zfac' (below) should also correspond to
              ! those arrays here (instead of to those with sedimenting ice
              ! already removed).

              ! split residual between in-cloud ice crystals and snow
              zfac(1:lproma) = &
                   MAX(MIN(lprod(1:lproma) / liwc(1:lproma), 1.0_dp), 0._dp)
              DO jt = 1, aer_count(js)%numbers(c_all)
                lressnownuc(jt,1:lproma) = lressnownuc(jt,1:lproma) &
                     + lrescldinuc(jt,1:lproma) * zfac(1:lproma)
                lrescldinuc(jt,1:lproma) = &
                     lrescldinuc(jt,1:lproma) * (1._dp - zfac(1:lproma))
              END DO
              DO jt = 1, ispec
                lressnowchm(jt,1:lproma) = lressnowchm(jt,1:lproma) &
                     + lrescldichm(jt,1:lproma) * zfac(1:lproma)
                lrescldichm(jt,1:lproma) = &
                     lrescldichm(jt,1:lproma) * (1._dp - zfac(1:lproma))
              END DO

              ! adjust sedimenting residual if necessary (this gives precedence
              ! to snow production over sedimentation)
              DO jl = 1, lproma
                DO jt = 1, aer_count(js)%numbers(c_all)
                  lressedinuc(jt,jl) = &
                       MIN(lressedinuc(jt,jl), lrescldinuc(jt,jl))
                END DO
                DO jt = 1, ispec
                  lressedichm(jt,jl) = &
                       MIN(lressedichm(jt,jl), lrescldichm(jt,jl))
                END DO
              END DO

              ! split residual between in-cloud and sedimenting ice
              lrescldinuc(1:aer_count(js)%numbers(c_all),1:lproma) = &
                   lrescldinuc(1:aer_count(js)%numbers(c_all),1:lproma) &
                   - lressedinuc(1:aer_count(js)%numbers(c_all),1:lproma)
              lrescldichm(1:ispec,1:lproma) = lrescldichm(1:ispec,1:lproma) &
                   - lressedichm(1:ispec,1:lproma)
          ENDIF
!-------------------------------------
           aerfield_evapo(:,:) = 0._dp
           kppfield_evapo(:,:) = 0._dp

           IF (ALLOCATED(made3)) THEN
             DO jl = 1, lproma
               i = lwork(jl)
               aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
                    lrescldinuc(1:aer_count(js)%numbers(c_all),jl)
               kppfield_evapo(1:ispec,i) = lrescldichm(1:ispec,jl)
             END DO
             DO jl = 1, mproma
               i = mwork(jl)
               aerfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                    mrescldinuc(1:aer_count(js)%numbers(c_all),jl)
               kppfield_evapo(1:ispec,i) = mrescldichm(1:ispec,jl)
             END DO
           ELSE

           DO jl=1,lproma
             i =lwork(jl)
             aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
               lnewflx(1:aer_count(js)%numbers(c_all),jl)
             kppfield_evapo(1:ispec,i) = lkppfld(1:ispec,jl)
           ENDDO
           DO jl=1,mproma
             i =mwork(jl)
             aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
               maerflx(1:aer_count(js)%numbers(c_all),jl)
             kppfield_evapo(1:ispec,i) = mkppfld(1:ispec,jl)
           ENDDO

           END IF

           IF (n_resmod .GT. 1) THEN

              DO idt = 1, SIZE(process(4)%frac_eva, 1)
                 DO j = 1, SIZE(process(4)%frac_eva, 2)
                    tmpbuffer(1:kproma,idt,j) = &
                         process(4)%frac_eva(idt,j)%ptr(1:kproma,jk)  
                 END DO
              END DO

             CALL calc_evapo_frac(status, aerfield_evapo(:,1:kproma), &
                  kppfield_evapo(1:ispec,1:kproma), kproma, ispec, 2, &
                  js, ntrac, tmpbuffer(1:kproma,:,:))
             IF (status .NE. 0) RETURN
             DO idt = 1, SIZE(process(4)%frac_eva, 1)
                DO j = 1, SIZE(process(4)%frac_eva, 2)
                   process(4)%frac_eva(idt,j)%ptr(1:kproma,jk) = &
                        tmpbuffer(1:kproma,idt,j)
                END DO
             END DO
           ELSE
              DO idt = 1, ntrac
                 process(4)%frac_eva(idt,1)%ptr(1:kproma,jk) = 1._dp
              END DO
           END IF

! evaporation for gaseous compounds sedimented in cloud species
            if (mproma > 0) THEN
              DO jl=1,mproma
                i = mwork(jl)
                evfrac(jl) = 0.0_dp
                DO idt = 1, ntrac
                   ! All (ntrac) elements of
                   ! procces(4)%frac_eva(i,jk,:,n_resmod) should have the same
                   ! value here. To avoid arbitrarily selecting one of these
                   ! elements, they are summed and divided by ntrac.
                   evfrac(jl) = evfrac(jl) + &
                   process(4)%frac_eva(idt,n_resmod)%ptr(i,jk)
                END DO
                evfrac(jl) = evfrac(jl) / REAL(ntrac, kind=dp)
                DO idt = 1, ntrac
                   DO j = 1, n_resmod
                      evfrac_aer(jl,idt,j) = &
                           process(4)%frac_eva(idt,j)%ptr(i,jk)
                   END DO
                END DO
             ENDDO

              IF (ALLOCATED(made3)) evfrac(1:mproma) = 1._dp

              CALL evapo(mkppfld(:,1:mproma), mxt_vec(:,1:mproma), &
                         mgrvol(1:mproma), evfrac(1:mproma),       &
                         ntrac, mproma, mwork(1:mproma), js, jk,   &
                         jrow, 2, kspec)
            END IF
                               
! evaporation of remaining in-cloud species are done in scav_physc(2)
!   at the end of the timestep
!--------------------------------------

!--------------------------------------
!  update aerosol deposition flux
          if (lproma > 0) THEN
            do jl=1,lproma
              if (lprod(jl) > 0) THEN
                frac = MIN(lprod(jl) / liwc(jl), 1._dp)
                laerflx(1:aer_count(js)%numbers(c_all),jl) =   &
                  laerflx(1:aer_count(js)%numbers(c_all),jl) + &
                  lnewflx(1:aer_count(js)%numbers(c_all),jl) * frac
                lnewflx(1:aer_count(js)%numbers(c_all),jl) =   &
                  lnewflx(1:aer_count(js)%numbers(c_all),jl) * (1._dp - frac)
              ENDIF
            enddo
          endif


! evaporation of sedimenting aerosol species
          IF (MPROMA > 0) THEN

            IF (.NOT. ALLOCATED(made3)) THEN
!
               CALL evapo_aer(mxt_vec(:,1:mproma), maerflx(:,1:mproma), &
                    mgrvol(1:mproma), evfrac_aer(1:mproma,:,:),         &
                    ntrac, mproma, mwork(1:mproma), js, jk, jrow)

            ELSE
               CALL evapo_aer(mxt_vec(:,1:mproma), mrescldinuc(:,1:mproma), &
                    mgrvol(1:mproma), evfrac_aer(1:mproma,:,:),             &
                    ntrac, mproma, mwork(1:mproma), js, jk, jrow)
               mrescldichm(:,1:mproma) = 0._dp
            END IF

          END IF

! evaporation of in-cloud species are done in scav_physc(2)

!-------------------------------------------------------------------------
! unpack for snow and ice
          IF (lproma > 0) THEN
            DO jt=1,ntrac
              DO jl=1,lproma
                i=lwork(jl)
                xt(i,jk,jt) = xt_vec(jt,jl)
              ENDDO
            ENDDO
! for output of partitioning parameters
            DO jl=1,lproma
              i=lwork(jl)
              phi_i(i,jk) = lphi_i(jl)
              mju_i(i,jk) = lmju_i(jl)
              iwc_i(i,jk) = liwc_i(jl)
              iwc_T_i(i,jk) = liwc_T_i(jl)
            ENDDO
!
            DO jt=1,ispec
              DO jl=1,lproma
                i=lwork(jl)
                snow_field(i,jt)             = lkppflx(jt,jl)
                kpp_field_cloud2(jt,i)       = lkppsediflx(jt,jl)
                cloudi_field(i,jt,jk)        = lkppfld(jt,jl)
!
                ! store residuals for other processes and for the layer below
                ressnowchm(i,jt)    = lressnowchm(jt,jl)
                rescldichm(i,jt)    = lressedichm(jt,jl)
!
              ENDDO
            ENDDO

            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                new_aerflx(i,jt)    = laersediflx(jt,jl)
                snow_aer(i,jt)      = laerflx(jt,jl)

                IF (.NOT. ALLOCATED(made3)) THEN
                   cloudi_aer(i,jt,jk) = lnewflx(jt,jl)
                ELSE
                   cloudi_aer(i,jt,jk) = lrescldinuc(jt,jl)
                END IF
                ! store residuals for other processes and for the layer below
                ressnownuc(i,jt)    = lressnownuc(jt,jl)
                rescldinuc(i,jt)    = lressedinuc(jt,jl)
              ENDDO
            ENDDO


            IF (jk == nlev) THEN
              DO jl=1,lproma
                i=lwork(jl)
                snow_aer(i,1:aer_count(js)%numbers(c_all)) =   &
                  snow_aer(i,1:aer_count(js)%numbers(c_all)) + & 
                  new_aerflx(i,1:aer_count(js)%numbers(c_all))
                snow_field(i,1:ispec) = snow_field(i,1:ispec) + &
                  kpp_field_cloud2(1:ispec,i)
              END DO
              new_aerflx(:,:) = 0._dp
              kpp_field_cloud2(:,:) = 0._dp
            ENDIF
                
          ENDIF
          IF (mproma > 0) THEN
            DO jt=1,ntrac
              DO jl=1,mproma
                i=mwork(jl)
                xt(i,jk,jt) = mxt_vec(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,ispec
              DO jl=1,mproma
                i=mwork(jl)
                kpp_field_cloud2(jt,i)    = mkppfld(jt,jl)
                rescldichm(i,jt) = mrescldichm(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,mproma
                i=mwork(jl)
                new_aerflx(i,jt) = maerflx(jt,jl)
                rescldinuc(i,jt) = mrescldinuc(jt,jl)
              ENDDO
            ENDDO
          ENDIF

        ENDIF ice_nucl_scav
!---------------------------------------------------------------------------
! liquid phase
        liq_nucl_scav: if (LSCAV_L) THEN

! pack for liquid water clouds
! into a vector that contains all the boxes 
! with liquid nucleation scavenging
          phase  = 1
          lproma = 0
          llwc(:)          = 0._dp
          liwc(:)          = 0._dp
          lrdrad(:)        = 0._dp
          lpress(:)        = 0._dp
          ltemp(:)         = 0._dp
          lcov(:)          = 0._dp
          lgrvol(:)        = 0._dp
          xt_vec(:,:)      = 0._dp
          laerflx(:,:)     = 0._dp
          lnewflx(:,:)     = 0._dp
          lkppflx(:,:)     = 0._dp
          lkppfld(:,:)     = 0._dp

          lrescldlnuc(:,:) = 0._dp
          lresrainnuc(:,:) = 0._dp
          lrescldlchm(:,:) = 0._dp
          lresrainchm(:,:) = 0._dp

          lrho(:)          = 0._dp
          lprod(:)         = 0._dp
          lhp(:)           = 0._dp
          lph(:)           = 0._dp
          lsteps(:)        = 0._dp
          lrsteps(:)       = 0._dp
          lwork(:) = 0
          cprod_array(:,:) = 0._dp

          kpp_rain_help(:) = 0._dp
          aer_rain_help(:) = 0._dp 

          mwork(:) = 0
          mproma   = 0
          
          mkppfld(:,:)     = 0._dp
          maerflx(:,:)     = 0._dp
          mxt_vec(:,:)     = 0._dp
          evfrac(:)        = 0._dp

          evfrac_aer(:,:,:)= 0._dp

          IF (L_LS) THEN
            DO jt=1,lspec
              DO jl=1,kproma
                kpp_rain_help(jl) = MAX(kpp_rain_help(jl), &
                                        cloudl_field(jl,jt,jk))
              END DO
            END DO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,kproma
                aer_rain_help(jl) = MAX(aer_rain_help(jl), &
                                        cloudl_aer(jl,jt,jk) )
              END DO
            END DO
            lnucl_l(1:kproma,jk) = .TRUE.
          END IF

          DO jl=1,kproma
            IF ( (lwc(jl,jk)   > thres_val) .AND. &
                 (cover(jl,jk) > 1.e-7_dp) ) then
              lproma = lproma + 1
              lwork(lproma)=jl
              lnucl_l(jl,jk) = .TRUE.
            ELSE
              IF ( (kpp_rain_help(jl) > 0.0_dp) .or. &
                   (aer_rain_help(jl) > 0.0_dp) ) THEN
                mproma = mproma + 1
                mwork(mproma)=jl
              ENDIF
            ENDIF
          ENDDO
          
          if (lproma > 0) THEN
            DO jl=1,lproma
              i=lwork(jl)
              llwc(jl)       = lwc(i,jk)
              lrdrad(jl)     = rdrad(i)
              lpress(jl)     = press(i,jk)
              ltemp(jl)      = temp(i,jk)
              lrho(jl)       = rho(i,jk)
              lcov(jl)       = cover(i,jk)
              lgrvol(jl)     = vol(i,jk)
              lprod(jl)      = rainprod(i,jk)
              lph(jl)        = cloudpH(i,jk)
            ENDDO
            DO jt=1,ntrac
              DO jl=1,lproma
                i=lwork(jl)
                xt_vec(jt,jl)   = xt(i,jk,jt)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                laerflx(jt,jl) = rain_aer(i,jt)
                lnewflx(jt,jl) = cloudl_aer(i,jt,jk)
                cloudl_aer(i,jt,jk) = 0._dp

                ! rain residual aerosol
                lresrainnuc(jt,jl) = resrainnuc(i,jt)
              ENDDO
            ENDDO
            DO jt=1,lspec
              DO jl=1,lproma
                i=lwork(jl)
                lkppflx(jt,jl)  = rain_field(i,jt)
                lkppfld(jt,jl)  = cloudl_field(i,jt,jk)
                cloudl_field(i,jt,jk) = 0._dp
                lresrainchm(jt,jl) = resrainchm(i,jt)
              ENDDO
            ENDDO
            DO jt=1,nprod
              DO jl=1,lproma
                i=lwork(jl)
                cprod_array(jl,jt) = prod_array(i,jk,jt,1) * time_step_len
              ENDDO
            ENDDO
          ENDIF ! lproma
          
          IF (mproma > 0) THEN
            DO jt=1,ntrac
              DO jl=1,mproma
                i=mwork(jl)
                mxt_vec(jt,jl)   = xt(i,jk,jt)
              ENDDO
            ENDDO
            DO jl=1,mproma
              i=mwork(jl)
              mgrvol(jl) = vol(i,jk)
            ENDDO
            DO jt=1,lspec
              DO jl=1,mproma
                i=mwork(jl)
                mkppfld(jt,jl)  = cloudl_field(i,jt,jk)
                cloudl_field(i,jt,jk) = 0._dp
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,mproma
                i=mwork(jl)
                maerflx(jt,jl) = cloudl_aer(i,jt,jk)
              ENDDO
            ENDDO
          ENDIF ! on mproma

!---------------------------------------------------------------
! aerosol nucleation scavenging
          if (lscav_aer) THEN
           
            if (lproma > 0) THEN

              CALL AER_SCAV_DRV(lproma, ntrac, jk, js, 2,             &
                phase, lwork(1:lproma), time_step_len,                &
                xt_vec(:,1:lproma),                                   &
                lnewflx(1:aer_count(js)%numbers(c_all),1:lproma),     &
                llwc(1:lproma),   liwc(1:lproma),   lcov(1:lproma),   &
                lrdrad(1:lproma), ltemp(1:lproma),  lpress(1:lproma), &
                lgrvol(1:lproma), lrho(1:lproma),   lprec(1:lproma),  &
                lsnow(1:lproma)  )

              DO jt=1,aer_count(js)%numbers(c_all)
                DO jl=1,lproma
                  i=lwork(jl)
                  lrescldlnuc(jt,jl) = &
                    MAX(lnewflx(jt,jl) - cloudl_aer(i,jt,jk), 0._dp)
                ENDDO
              ENDDO

            END IF ! on lproma > 0

          endif ! on lscav_aer
!-------------------------------------------------------------------------------
! gas scavenging and liquid chemistry
          if (lscav_gas) then

            IF (lproma > 0) then
              CALL in_cloud(lkppfld(:,1:lproma), xt_vec(:,1:lproma), &
                            lnewflx(:,1:lproma), lgrvol(1:lproma),   &
                            lcov(1:lproma),                          &
                            ntrac, lproma, js, 1, kspec)

              ! Subtract current concentrations of KPP species in order to
              ! calculate increment (below)
              DO jt = 1, lspec
                lrescldlchm(jt,1:lproma) = lrescldlchm(jt,1:lproma) &
                   - L_SPEC(jt,1:lproma) * lgrvol(1:lproma) * 1.e6_dp &
                   * lcov(1:lproma)
              END DO

              if (lscav_easy) then
                call scav_easy_liq(llwc(1:lproma), ltemp(1:lproma), lproma,js)
              else
                call scav_aqchem_kpp(ltemp(1:lproma),  lpress(1:lproma),   &
                                     llwc(1:lproma),   llwc(1:lproma),     &
                                     lrdrad(1:lproma), time_step_len, 1,   &
                                     lproma, lwork, jk, jrow, js,          &
                                     lsteps(1:lproma), lrsteps(1:lproma) )
              endif
              call out_cloud(lkppfld(:,1:lproma), lkppflx(:,1:lproma),      &
                             xt_vec(:,1:lproma), lph(1:lproma),             &
                             lhp(1:lproma), lcov(1:lproma), llwc(1:lproma), &
                             lprod(1:lproma), lgrvol(1:lproma),             &
                             ntrac, lproma, js, 1, kspec,                   &
                             cprod_array(1:lproma,:) )

              ! Calculate increment of KPP species concentrations due to
              ! chemistry
              DO jt = 1, lspec
                lrescldlchm(jt,1:lproma) = lrescldlchm(jt,1:lproma) &
                     + L_SPEC(jt,1:lproma) * lgrvol(1:lproma) * 1.e6_dp &
                     * lcov(1:lproma)
              END DO

              lhp(1:lproma) = lhp(1:lproma) / (1.0e-6_dp*lrho(1:lproma))
                
            ENDIF
              
          ENDIF   ! on lscav_gas          

          ! split residual between in-cloud droplets and rain
          IF (lproma .GT. 0) THEN
            zfac(1:lproma) = &
                 MAX(MIN(lprod(1:lproma) / llwc(1:lproma), 1.0_dp), 0._dp)
            DO jt = 1, aer_count(js)%numbers(c_all)
              lresrainnuc(jt,1:lproma) = lresrainnuc(jt,1:lproma) &
                   + lrescldlnuc(jt,1:lproma) * zfac(1:lproma)
              lrescldlnuc(jt,1:lproma) = &
                   lrescldlnuc(jt,1:lproma) * (1._dp - zfac(1:lproma))
            END DO
            DO jt = 1, lspec
              lresrainchm(jt,1:lproma) = lresrainchm(jt,1:lproma) &
                   + lrescldlchm(jt,1:lproma) * zfac(1:lproma)
              lrescldlchm(jt,1:lproma) = &
                   lrescldlchm(jt,1:lproma) * (1._dp - zfac(1:lproma))
            END DO
          END IF

          aerfield_evapo(:,:) = 0._dp
          kppfield_evapo(:,:) = 0._dp

          IF (ALLOCATED(made3)) THEN
            DO jl = 1, lproma
              i = lwork(jl)
              aerfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                   lrescldlnuc(1:aer_count(js)%numbers(c_all),jl)
              kppfield_evapo(1:lspec,i) = lrescldlchm(1:lspec,jl)
            END DO
          ELSE

             DO jl=1,lproma
                i =lwork(jl)
                aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
                     lnewflx(1:aer_count(js)%numbers(c_all),jl)
                kppfield_evapo(1:lspec,i) = lkppfld(1:lspec,jl)
             ENDDO
             DO jl=1,mproma
                i =mwork(jl)
                aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
                     maerflx(1:aer_count(js)%numbers(c_all),jl)
                kppfield_evapo(1:lspec,i) = mkppfld(1:lspec,jl)
             ENDDO

          END IF

           IF (n_resmod .GT. 1) THEN
              DO idt = 1, SIZE(process(3)%frac_eva, 1)
                 DO j = 1, SIZE(process(3)%frac_eva, 2)
                    tmpbuffer(1:kproma,idt,j) = &
                         process(3)%frac_eva(idt,j)%ptr(1:kproma,jk)  
                 END DO
              END DO

             CALL calc_evapo_frac(status, aerfield_evapo(:,1:kproma), &
                  kppfield_evapo(1:lspec,1:kproma), kproma, lspec, 1, &
                  js, ntrac, tmpbuffer(1:kproma,:,:))
             IF (status .NE. 0) RETURN

             DO idt = 1, SIZE(process(3)%frac_eva, 1)
                DO j = 1, SIZE(process(3)%frac_eva, 2)
                   process(3)%frac_eva(idt,j)%ptr(1:kproma,jk) = &
                        tmpbuffer(1:kproma,idt,j)
                END DO
             END DO
           ELSE
              DO idt = 1, ntrac
                 process(3)%frac_eva(idt,1)%ptr(1:kproma,jk) = 1._dp
              END DO
           END IF

!--------------------------------------------------
! update aerosol deposition flux
          if (lproma > 0) THEN
            do jl=1,lproma
              if (lprod(jl) > 0) THEN
                frac = MIN(lprod(jl) / llwc(jl), 1._dp)
                laerflx(1:aer_count(js)%numbers(c_all),jl) =   &
                  laerflx(1:aer_count(js)%numbers(c_all),jl) + &
                  lnewflx(1:aer_count(js)%numbers(c_all),jl) * frac
                lnewflx(1:aer_count(js)%numbers(c_all),jl) =   &
                  lnewflx(1:aer_count(js)%numbers(c_all),jl) * (1._dp - frac)
              ENDIF
            enddo
          endif
!------------------------------------------------------------------
! evaporation of species in clouds (from previous time steps)
          IF (MPROMA > 0) THEN
            DO jl=1,mproma
              i = mwork(jl)
                evfrac(jl) = 0.0_dp
                DO idt = 1, ntrac
                   ! For explanation, see above (process(4))
                   evfrac(jl) = evfrac(jl) + &
                   process(3)%frac_eva(idt,n_resmod)%ptr(i,jk)
                END DO
                evfrac(jl) = evfrac(jl) / REAL(ntrac, kind=dp)
                DO idt = 1, ntrac
                   DO j = 1, n_resmod
                      evfrac_aer(jl,idt,j) = &
                           process(3)%frac_eva(idt,j)%ptr(i,jk)
                   END DO
                END DO
            ENDDO
            IF (ALLOCATED(made3)) evfrac(1:mproma) = 1._dp

            CALL evapo(mkppfld(:,1:mproma), mxt_vec(:,1:mproma), &
                       mgrvol(1:mproma), evfrac(1:mproma),       &
                       ntrac, mproma, mwork(1:mproma), js, jk,   &
                       jrow, 1, kspec)
!!            print*, "MXTVEC_before: ", jk, MAXVAL(mxt_vec(153,1:mproma)),&
!!              MAXVAL(mxt_vec(250,1:mproma))
!!            print*, "incloud field before: ",jk, MAXVAL(maerflx(:,:))

            ! Not (currently) necessary for MADE3, as no cloud phase aerosol
            ! tracers defined. Hence, maerflx should not contain any aerosol (in
            ! contrast to ice_nucl_scav, where it contains sedimented ice).
            IF (.NOT. ALLOCATED(made3)) THEN
               CALL evapo_aer(mxt_vec(:,1:mproma), maerflx(:,1:mproma), &
                    mgrvol(1:mproma), evfrac_aer(1:mproma,:,:),         &
                    ntrac, mproma, mwork(1:mproma), js, jk, jrow)
            END IF
          ENDIF
!------------------------------------------------------------------
! unpacking for liquid water clouds
          if (lproma > 0) THEN
            DO jt=1,ntrac
              DO jl=1,lproma
                i=lwork(jl)
                xt(i,jk,jt) = xt_vec(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,lspec
              DO jl=1,lproma
                i=lwork(jl)
                rain_field(i,jt)      = lkppflx(jt,jl)
                cloudl_field(i,jt,jk) = lkppfld(jt,jl)
                resrainchm(i,jt)      = lresrainchm(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                rain_aer(i,jt)      = laerflx(jt,jl) 

                IF (.NOT. ALLOCATED(made3)) THEN
                   cloudl_aer(i,jt,jk) = lnewflx(jt,jl)
                ELSE
                   cloudl_aer(i,jt,jk) = lrescldlnuc(jt,jl)
                END IF
                ! store residuals for other processes and for the layer below
                resrainnuc(i,jt)    = lresrainnuc(jt,jl)
              ENDDO
            ENDDO
            DO jl=1,lproma
              i=lwork(jl)
              cloudpH(i,jk)  = lph(jl)
              Hp_cloud(i,jk) = lhp(jl)
              xsteps_cloud(i,jk)   = lsteps(jl)
              xrsteps_cloud(i,jk)  = lrsteps(jl)
            ENDDO
            DO jt=1,nprod
              DO jl=1,lproma
                i=lwork(jl)
                prod_array(i,jk,jt,1) = cprod_array(jl,jt) / time_step_len
              ENDDO
            ENDDO
          ENDIF  ! on lproma

          IF (MPROMA > 0) THEN
            DO jt=1,ntrac
              DO jl=1,mproma
                i=mwork(jl)
                xt(i,jk,jt) = mxt_vec(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,ispec
              DO jl=1,mproma
                i=mwork(jl)
                cloudl_field(i,jt,jk) = mkppfld(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,mproma
                i=mwork(jl)
                cloudl_aer(i,jt,jk) = maerflx(jt,jl)
              ENDDO
            ENDDO
          ENDIF
!-------------------------------------------
        ENDIF liq_nucl_scav
!-------------------------------------------
      END if nucleation_scavenging
!-------------------------------------------
!-------------------------------------------------
     ziwc(:)   = 0._dp
     zlwc(:)   = 0._dp
     rdrad(:)  = 0._dp
     lwc_bc(:) = 0._dp

     impaction_scavenging: if (lscav_imp) then

       do jl=1,kproma
         ! determine mean droplet radius and precipitation lwc 
         ! in terms of falling precip 
         ! (used for calculation of the transfer velocity)
         call calc_lwc_bc_n(rain(jl,jk),   snow(jl,jk),  rcover(jl,jk), &
                            lwc_bc(jl),    rdrad(jl),    col_cov(jl)    )
         bc_lwc(jl,jk) = lwc_bc(jl)

         ! determine total precipitation lwc and iwc
         if ( (rain(jl,jk) > 1.e-12_dp) .and. (col_cov(jl) >  1.e-10_dp) )  &
           zlwc(jl) = rain(jl,jk) / col_cov(jl) * time_step_len *    &
                      gboxarea(jl) / vol(jl,jk)

         if ( (snow(jl,jk) > 1.e-12_dp) .and. (col_cov(jl) >  1.e-10_dp) )  &
           ziwc(jl) = snow(jl,jk) / col_cov(jl) * time_step_len *    &
                      gboxarea(jl) / vol(jl,jk)

       enddo

       ice_imp_scav: if (lscav_i) then

! pack for snow and ice
! into a vector that contains all the boxes with frozen impaction scavenging
         phase = 2 ! ice

         lproma = 0
         llwc(:)          = 0._dp
         liwc(:)          = 0._dp             
         lrdrad(:)        = 0._dp           
         lpress(:)        = 0._dp       
         ltemp(:)         = 0._dp        
         lcov(:)          = 0._dp         
         lrho(:)          = 0._dp         
         lgrvol(:)        = 0._dp       
         lprec(:)         = 0._dp
         lsnow(:)         = 0._dp   
         xt_vec(:,:)      = 0._dp     
         laerflx(:,:)     = 0._dp

         lressnownuc(:,:) = 0._dp
         lresrainnuc(:,:) = 0._dp
         lressnowimp(:,:) = 0._dp
         lresrainimp(:,:) = 0._dp
         lressnowchm(:,:) = 0._dp
         lresrainchm(:,:) = 0._dp

         lprod(:)         = 0._dp
         lkppflx(:,:)     = 0._dp
         lkppiflx(:,:)    = 0._dp
! for output of partitioning parameters
         lphi_i(:)        = 0._dp
         lmju_i(:)        = 0._dp
         liwc_i(:)        = 0._dp
         liwc_T_i(:)      = 0._dp
!
         lmelt(:)         = 0._dp
         lhp(:)           = 0._dp
         lph(:)           = 0._dp
         lwork(:) = 0
         cprod_array(:,:) = 0._dp
         
         mproma = 0
         mkppiflx(:,:)    = 0._dp
         maeriflx(:,:)    = 0._dp
         maerflx(:,:)     = 0._dp
         mressnownuc(:,:) = 0._dp
         mresrainnuc(:,:) = 0._dp
         mressnowimp(:,:) = 0._dp
         mresrainimp(:,:) = 0._dp
         mressnowchm(:,:) = 0._dp
         mresrainchm(:,:) = 0._dp

         mxt_vec(:,:)     = 0._dp
         mgrvol(:)        = 0._dp
         mmelt(:)         = 0._dp
         mwork(:) = 0
         evfrac(:)        = 0._dp
         evfrac_aer(:,:,:)= 0._dp

         ! replaced maxval statement
         kpp_snow_help(:) = 0._dp
         aer_snow_help(:) = 0._dp
         DO jt=1,ispec
           DO JL=1,kproma
             kpp_snow_help(jl) = MAX(kpp_snow_help(jl), snow_field(jl,jt))
           ENDDO
         ENDDO
         DO jt=1,aer_count(js)%numbers(c_all)
           DO jl=1,kproma
             aer_snow_help(jl) = MAX(aer_snow_help(jl), snow_aer(jl,jt))
           ENDDO
         ENDDO

         ! determine mask
         do jl=1,kproma
           if (ziwc(jl) > thres_val) THEN
             lproma = lproma + 1
             lwork(lproma) = jl
           ELSE
             IF ( (kpp_snow_help(jl) > 0.0_dp) .or. &
                  (aer_snow_help(jl) > 0.0_dp) ) THEN
               mproma = mproma + 1
               mwork(mproma)=jl
             ENDIF
           ENDIF
         enddo
         ! packing for lproma
         if (lproma > 0) THEN
           do jl=1,lproma
             i=lwork(jl)
             liwc(jl)               = ziwc(i)
             lrdrad(jl)             = rdrad(i)
             lpress(jl)             = press(i,jk)
             ltemp(jl)              = temp(i,jk)
             lcov(jl)               = col_cov(i)
             lrho(jl)               = rho(i,jk)
             lgrvol(jl)             = vol(i,jk)
             lsnow(jl)              = snow(i,jk)
             lmelt(jl)              = imelt(i,jk)
             lsnow_up(jl)           = snow(i,jk-1)
! for output of partitioning parameters
             lphi_i(jl)             = phi_i(i,jk)
             lmju_i(jl)             = mju_i(i,jk)
             liwc_i(jl)             = iwc_i(i,jk)
             liwc_T_i(jl)           = iwc_T_i(i,jk) 
!
           ENDDO
           DO jt=1,ntrac
             DO jl=1,lproma
               i=lwork(jl)
               xt_vec(jt,jl)   = xt(i,jk,jt)
             ENDDO
           ENDDO
           DO jt=1,aer_count(js)%numbers(c_all)
             DO jl=1,lproma
               i=lwork(jl)
               laerflx(jt,jl)  = rain_aer(i,jt)
               laeriflx(jt,jl) = snow_aer(i,jt)

               lressnownuc(jt,jl) = ressnownuc(i,jt)
               lresrainnuc(jt,jl) = resrainnuc(i,jt)
               lressnowimp(jt,jl) = ressnowimp(i,jt)
               lresrainimp(jt,jl) = resrainimp(i,jt)
             ENDDO
           ENDDO
           DO jt=1,lspec
             DO jl=1,lproma
               i=lwork(jl)
               lkppflx(jt,jl)   = rain_field(i,jt)
               lresrainchm(jt,jl) = resrainchm(i,jt)
             ENDDO
           ENDDO
           DO jt=1,ispec
             DO jl=1,lproma
               i=lwork(jl)
               lkppiflx(jt,jl)  = snow_field(i,jt)
               lressnowchm(jt,jl) = ressnowchm(i,jt)
             ENDDO
           ENDDO
         endif  ! packing for lproma

         if (mproma > 0) THEN   ! packing for mproma
           DO jt=1,lspec
             DO jl=1,mproma
               i=mwork(jl)
               mkppflx(jt,jl)  = rain_field(i,jt)
               mresrainchm(jt,jl) = resrainchm(i,jt)
             ENDDO
           ENDDO
           DO jt=1,ispec
             DO jl=1,mproma
               i=mwork(jl)
               mkppiflx(jt,jl) = snow_field(i,jt)
               mressnowchm(jt,jl) = ressnowchm(i,jt)
             ENDDO
           ENDDO
           DO jt=1,aer_count(js)%numbers(c_all)
             DO jl=1,mproma
               i=mwork(jl)
               maerflx(jt,jl)  = rain_aer(i,jt)
               maeriflx(jt,jl) = snow_aer(i,jt)
               mresrainnuc(jt,jl) = resrainnuc(i,jt)
               mresrainimp(jt,jl) = resrainimp(i,jt)
               mressnownuc(jt,jl) = ressnownuc(i,jt)
               mressnowimp(jt,jl) = ressnowimp(i,jt)
             ENDDO
           ENDDO
           DO jt=1,ntrac
             DO jl=1,mproma
               i=mwork(jl)
               mxt_vec(jt,jl) = xt(i,jk,jt)
             ENDDO
           ENDDO
           DO jl=1,mproma
             i=mwork(jl)
             mgrvol(jl)          = vol(i,jk)
             mmelt(jl)           = imelt(i,jk)
             msnow_up(jl)        = snow(i,jk-1)
           ENDDO
         ENDIF ! packing for mproma

!-----------------------------------------------------------
!   aerosol scavenging
         if (lscav_aer) then
           if (lproma > 0) THEN

             CALL AER_SCAV_DRV(lproma, ntrac, jk, js, 1,           &
             phase, lwork(1:lproma), time_step_len,                &
             xt_vec(:,1:lproma),                                   &
             laeriflx(1:aer_count(js)%numbers(c_all),1:lproma),    &
             llwc(1:lproma),   liwc(1:lproma),   lcov(1:lproma),   &
             lrdrad(1:lproma), ltemp(1:lproma),  lpress(1:lproma), &
             lgrvol(1:lproma), lrho(1:lproma),   lprec(1:lproma),  &
             lsnow(1:lproma)  )

             DO jt = 1, aer_count(js)%numbers(c_all)
               DO jl = 1, lproma
                 i = lwork(jl)
                 lressnowimp(jt,jl) = lressnowimp(jt,jl) + &
                   MAX(laeriflx(jt,jl) - snow_aer(i,jt), 0._dp)
               ENDDO
             ENDDO

           END IF     ! lproma > 0
         endif     ! lscav_aer
!-----------------------------------------------------------
! gas scavenging and aqueous phase chemistry

         if (lscav_gas) THEN

           IF (lproma > 0) then
              CALL in_flux(lkppiflx(:,1:lproma), xt_vec(:,1:lproma),  &
                           laeriflx(:,1:lproma), lgrvol(1:lproma),    &
                           lcov(1:lproma),                            &
                           ntrac, lproma, js, 2, kspec)
              ! Subtract current concentrations of KPP species in order to
              ! calculate increment (below)
              DO jt = 1, ispec
                lressnowchm(jt,1:lproma) = lressnowchm(jt,1:lproma) &
                     - I_SPEC(jt,1:lproma) * lgrvol(1:lproma) * 1.e6_dp &
                     * lcov(1:lproma)
              END DO

              if (iscav_easy) then
                call scav_easy_ice(liwc(1:lproma), ltemp(1:lproma),    &
                                   lpress(1:lproma),                   &
! for output of partitioning parameters 
                                   lphi_i(1:lproma),                   &
                                   lmju_i(1:lproma), liwc_i(1:lproma), &
                                   liwc_T_i(1:lproma),                 &
!
                                   lproma, js )

              else
!                call scav_icechem_kpp(ltemp, lpress, llwc2, llwc, lrdrad, &
!                                      time_step_len, 2, lproma, lwork,    &
!                                      jk, jrow)
                print*, "icechem"
              endif

              call out_flux(lkppiflx(:,1:lproma), xt_vec(:,1:lproma), &
                            lph(1:lproma), liwc(1:lproma),            &
                            lgrvol(1:lproma), lcov(1:lproma),         &
                            ntrac, lproma, js, 2, kspec,              &
                            cprod_array(1:lproma,:) )
              ! Calculate increment of KPP species concentrations due to
              ! chemistry
              DO jt = 1, ispec
                lressnowchm(jt,1:lproma) = lressnowchm(jt,1:lproma) &
                     + I_SPEC(jt,1:lproma) * lgrvol(1:lproma) * 1.e6_dp &
                     * lcov(1:lproma)
              END DO

            ENDIF

          ENDIF   ! scav_gas


!-----------------------------------------------------------
! check for melting transitions
          if (lproma > 0) THEN

            CALL MELTING(lproma,  js, kspec,                        &
                         lmelt(1:lproma), lsnow_up(1:lproma),       &
                         lkppiflx(:,1:lproma), lkppflx(:,1:lproma), &
                         laeriflx(:,1:lproma), laerflx(:,1:lproma))

            DO jl = 1, lproma
              IF (lsnow_up(jl) > 1e-15_dp) THEN
                zfac(jl) = MAX(MIN(lmelt(jl)/lsnow_up(jl), 1._dp), 0._dp)
              ELSE
                zfac(jl) = 0._dp
              END IF
            END DO
            DO jt = 1, aer_count(js)%numbers(c_all)
              lresrainnuc(jt,1:lproma) = lresrainnuc(jt,1:lproma) &
                   + lressnownuc(jt,1:lproma) * zfac(1:lproma)
              lressnownuc(jt,1:lproma) = &
                   lressnownuc(jt,1:lproma) * (1._dp - zfac(1:lproma))
              lresrainimp(jt,1:lproma) = lresrainimp(jt,1:lproma) &
                   + lressnowimp(jt,1:lproma) * zfac(1:lproma)
              lressnowimp(jt,1:lproma) = &
                   lressnowimp(jt,1:lproma) * (1._dp - zfac(1:lproma))
            END DO
            DO jt = 1, ispec_ice
              idx2 = KPP_I_IDX(js)%ice_spec(jt,ice2liq)
              idx1 = KPP_I_IDX(js)%ice_spec(jt,ice_idx)
              lresrainchm(idx2,1:lproma) = lresrainchm(idx2,1:lproma) &
                   + lressnowchm(idx1,1:lproma) * zfac(1:lproma)
              lressnowchm(idx1,1:lproma) = &
                   lressnowchm(idx1,1:lproma) * (1._dp - zfac(1:lproma))
            END DO
          END IF

          if (mproma > 0) THEN

            CALL MELTING(mproma,  js, kspec,                        &
                         mmelt(1:mproma), msnow_up(1:mproma),       &
                         mkppiflx(:,1:mproma), mkppflx(:,1:mproma), &
                         maeriflx(:,1:mproma), maerflx(:,1:mproma))

            DO jl = 1, mproma
              IF (msnow_up(jl) > 1e-15_dp) THEN
                zfac(jl) = MAX(MIN(mmelt(jl)/msnow_up(jl), 1._dp), 0._dp)
              ELSE
                zfac(jl) = 0._dp
              END IF
            END DO
            DO jt = 1, aer_count(js)%numbers(c_all)
              mresrainnuc(jt,1:mproma) = mresrainnuc(jt,1:mproma) &
                   + mressnownuc(jt,1:mproma) * zfac(1:mproma)
              mressnownuc(jt,1:mproma) = &
                   mressnownuc(jt,1:mproma) * (1._dp - zfac(1:mproma))
              mresrainimp(jt,1:mproma) = mresrainimp(jt,1:mproma) &
                   + mressnowimp(jt,1:mproma) * zfac(1:mproma)
              mressnowimp(jt,1:mproma) = &
                   mressnowimp(jt,1:mproma) * (1._dp - zfac(1:mproma))
            END DO
            DO jt = 1, ispec_ice
              idx2 = KPP_I_IDX(js)%ice_spec(jt,ice2liq)
              idx1 = KPP_I_IDX(js)%ice_spec(jt,ice_idx)
              mresrainchm(idx2,1:mproma) = mresrainchm(idx2,1:mproma) &
                   + mressnowchm(idx1,1:mproma) * zfac(1:mproma)
              mressnowchm(idx1,1:mproma) = &
                   mressnowchm(idx1,1:mproma) * (1._dp - zfac(1:mproma))
            END DO
          END IF

! evaporation fraction per mode
           aerfield_evapo(:,:) = 0._dp
           kppfield_evapo(:,:) = 0._dp

           impfield_evapo(:,:) = 0._dp
           IF (ALLOCATED(made3)) THEN
             DO jl = 1, lproma
               i = lwork(jl)
               aerfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                    lressnownuc(1:aer_count(js)%numbers(c_all),jl)
               kppfield_evapo(1:ispec,i) = lressnowchm(1:ispec,jl)
               impfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                    lressnowimp(1:aer_count(js)%numbers(c_all),jl)
             END DO
             DO jl = 1, mproma
               i = mwork(jl)
               aerfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                    mressnownuc(1:aer_count(js)%numbers(c_all),jl)
               kppfield_evapo(1:ispec,i) = mressnowchm(1:ispec,jl)
               impfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                    mressnowimp(1:aer_count(js)%numbers(c_all),jl)
             END DO

              DO idt = 1, SIZE(process(2)%frac_eva, 1)
                 DO j = 1, SIZE(process(2)%frac_eva, 2)
                    tmpbuffer(1:kproma,idt,j) = &
                         process(2)%frac_eva(idt,j)%ptr(1:kproma,jk)  
                 END DO
              END DO

             CALL calc_evapo_frac(status, aerfield_evapo(:,1:kproma), &
                  kppfield_evapo(1:ispec,1:kproma), kproma, ispec, 2, &
                  js, ntrac, tmpbuffer(1:kproma,:,:),    & 
                  resimp=&
                  impfield_evapo(1:aer_count(js)%numbers(c_all),1:kproma))
             IF (status .NE. 0) RETURN

             DO idt = 1, SIZE(process(2)%frac_eva, 1)
                DO j = 1, SIZE(process(2)%frac_eva, 2)
                   process(2)%frac_eva(idt,j)%ptr(1:kproma,jk) = &
                        tmpbuffer(1:kproma,idt,j)
                END DO
             END DO

           ELSE

           DO jl=1,lproma
             i =lwork(jl)
             aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
               laeriflx(1:aer_count(js)%numbers(c_all),jl)
             kppfield_evapo(1:ispec,i) = lkppiflx(1:ispec,jl)
           ENDDO
           DO jl=1,mproma
             i =mwork(jl)
             aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
               maeriflx(1:aer_count(js)%numbers(c_all),jl)
             kppfield_evapo(1:ispec,i) = mkppiflx(1:ispec,jl)
           ENDDO
           IF (n_resmod .GT. 1) THEN

              DO idt = 1, SIZE(process(2)%frac_eva, 1)
                 DO j = 1, SIZE(process(2)%frac_eva, 2)
                    tmpbuffer(1:kproma,idt,j) = &
                         process(2)%frac_eva(idt,j)%ptr(1:kproma,jk)  
                 END DO
              END DO

              CALL calc_evapo_frac(status, aerfield_evapo(:,1:kproma), &
                   kppfield_evapo(1:ispec,1:kproma), kproma, ispec, 2, &
                   js, ntrac, tmpbuffer(1:kproma,:,:))
              IF (status .NE. 0) RETURN

              DO idt = 1, SIZE(process(2)%frac_eva, 1)
                 DO j = 1, SIZE(process(2)%frac_eva, 2)
                    process(2)%frac_eva(idt,j)%ptr(1:kproma,jk) = &
                         tmpbuffer(1:kproma,idt,j)
                 END DO
              END DO
           ELSE
              DO idt = 1, ntrac
                 process(2)%frac_eva(idt,1)%ptr(1:kproma,jk) = 1._dp
              END DO
           END IF
        END IF

! evaporation
          if (mproma > 0) THEN
            DO jl=1,mproma
              i = mwork(jl)
                evfrac(jl) = 0.0_dp
                DO idt = 1, ntrac
                   ! For explanation, see above (process(4))
                   evfrac(jl) = evfrac(jl) + &
                   process(2)%frac_eva(idt,n_resmod)%ptr(i,jk)
                END DO
                evfrac(jl) = evfrac(jl) / REAL(ntrac, kind=dp)
                DO idt = 1, ntrac
                   DO j = 1, n_resmod
                      evfrac_aer(jl,idt,j) = &
                           process(2)%frac_eva(idt,j)%ptr(i,jk)
                   END DO
                END DO
            ENDDO

            IF (ALLOCATED(made3)) evfrac(1:mproma) = 1._dp

            call evapo(mkppiflx(:,1:mproma), mxt_vec(:,1:mproma), &
                       mgrvol(1:mproma), evfrac(1:mproma),        &
                       ntrac, mproma, mwork(1:mproma), js, jk,    &
                       jrow, phase, kspec)

            IF (.NOT. ALLOCATED(made3)) THEN

               call evapo_aer(mxt_vec(:,1:mproma), maeriflx(:,1:mproma), &
                    mgrvol(1:mproma), evfrac_aer(1:mproma,:,:),          &
                    ntrac, mproma, mwork(1:mproma), js, jk, jrow)
            ELSE
               CALL evapo_aer(mxt_vec(:,1:mproma), mressnownuc(:,1:mproma), &
                    mgrvol(1:mproma), evfrac_aer(1:mproma,:,:),             &
                    ntrac, mproma, mwork(1:mproma), js, jk, jrow)
               mressnowimp(:,1:mproma) = 0._dp
               mressnowchm(:,1:mproma) = 0._dp
            END IF

         endif

!-----------------------------------------------------------
! unpack for snow and ice
          IF (lproma > 0 ) THEN
            DO jt=1,ntrac
              DO jl=1,lproma
                i=lwork(jl)
                xt(i,jk,jt) = xt_vec(jt,jl)
              ENDDO
            ENDDO
! for output of partitioning parameters
            DO jl=1,lproma
              i=lwork(jl)
              phi_i(i,jk) = lphi_i(jl)
              mju_i(i,jk) = lmju_i(jl)
              iwc_i(i,jk) = liwc_i(jl)
              iwc_T_i(i,jk) = liwc_T_i(jl)
            ENDDO

            DO jt=1,lspec
              DO jl=1,lproma
                i=lwork(jl)
                rain_field(i,jt) = lkppflx(jt,jl)

                ! store residuals for liq. imp. scav. and for the layers below
                resrainchm(i,jt) = lresrainchm(jt,jl)

              ENDDO
            ENDDO
            DO jt=1,ispec
              DO jl=1,lproma
                i=lwork(jl)
                snow_field(i,jt) = lkppiflx(jt,jl)

                ! store residuals for the layers below
                ressnowchm(i,jt) = lressnowchm(jt,jl)

              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                rain_aer(i,jt) = laerflx(jt,jl)
                snow_aer(i,jt) = laeriflx(jt,jl)

                ! store residuals for liq. imp. scav. and for the layers below
                resrainnuc(i,jt) = lresrainnuc(jt,jl)
                resrainimp(i,jt) = lresrainimp(jt,jl)
                ressnownuc(i,jt) = lressnownuc(jt,jl)
                ressnowimp(i,jt) = lressnowimp(jt,jl)

              ENDDO
            ENDDO
          ENDIF ! for lproma

          IF (mproma > 0) THEN
            DO jt=1,ntrac
              DO jl=1,mproma
                i=mwork(jl)
                xt(i,jk,jt) = mxt_vec(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,lspec
              DO jl=1,mproma
                i=mwork(jl)
                rain_field(i,jt) = mkppflx(jt,jl)
                resrainchm(i,jt) = mresrainchm(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,ispec
              DO jl=1,mproma
                i=mwork(jl)
                snow_field(i,jt) = mkppiflx(jt,jl)
                ressnowchm(i,jt) = mressnowchm(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,mproma
                i=mwork(jl)
                rain_aer(i,jt) = maerflx(jt,jl)
                snow_aer(i,jt) = maeriflx(jt,jl)

                resrainnuc(i,jt) = mresrainnuc(jt,jl)
                resrainimp(i,jt) = mresrainimp(jt,jl)
                ressnownuc(i,jt) = mressnownuc(jt,jl)
                ressnowimp(i,jt) = mressnowimp(jt,jl)
              ENDDO
            ENDDO
          ENDIF ! for mproma

        END if ice_imp_scav
!----------------------------------------------------
!----------------------------------------------------
        liq_imp_scav:  IF (LSCAV_L) THEN

! pack for rain liquid water
! into a vector that contains all the boxes with liquid impaction scavenging
          phase  = 1  ! water 

          lproma = 0
          liwc(:)          = 0._dp             
          llwc(:)          = 0._dp             
          llwc2(:)         = 0._dp             
          lrdrad(:)        = 0._dp           
          lpress(:)        = 0._dp       
          ltemp(:)         = 0._dp        
          lprec(:)         = 0._dp        
          lsnow(:)         = 0._dp        
          lcov(:)          = 0._dp         
          lrho(:)          = 0._dp         
          lgrvol(:)        = 0._dp       
          xt_vec(:,:)      = 0._dp     
          laerflx(:,:)     = 0._dp
          lkppflx(:,:)     = 0._dp

          lresrainnuc(:,:) = 0._dp
          lresrainimp(:,:) = 0._dp
          lresrainchm(:,:) = 0._dp

          lph(:)           = 0._dp
          lhp(:)           = 0._dp
          lsteps(:)        = 0._dp
          lrsteps(:)       = 0._dp
          cprod_array(:,:) = 0._dp

          mproma = 0
          mkppflx(:,:)     = 0._dp

          maerflx(:,:)     = 0._dp
          mresrainnuc(:,:) = 0._dp
          mresrainimp(:,:) = 0._dp
          mresrainchm(:,:) = 0._dp

          mxt_vec(:,:)     = 0._dp 
          mgrvol(:)        = 0._dp
          
          kpp_rain_help(:) = 0.0_dp
          aer_rain_help(:) = 0.0_dp
          evfrac(:)        = 0.0_dp

          evfrac_aer(:,:,:)= 0.0_dp

          ! replaced maxval statement
          DO jt=1,lspec
            DO JL=1,kproma
              kpp_rain_help(jl) = MAX(kpp_rain_help(jl), rain_field(jl,jt)) 
            ENDDO
          ENDDO
          DO jt=1,aer_count(js)%numbers(c_all)
            DO jl=1,kproma
              aer_rain_help(jl) = MAX(aer_rain_help(jl), rain_aer(jl,jt))
            ENDDO
          ENDDO
          
          lwork(:) = 0
          mwork(:) = 0
          ! determine mask
          DO jl=1,kproma
            IF (zlwc(jl) > thres_val) THEN
              IF (lwc_bc(jl) > 0._dp .AND. rdrad(jl) > 0._dp) THEN
                lproma = lproma + 1
                lwork(lproma)=jl
              ENDIF
            ELSE
              IF ( (kpp_rain_help(jl) > 0.0_dp) .or.      &
                   (aer_rain_help(jl) > 0.0_dp) ) THEN
                mproma = mproma + 1
                mwork(mproma)=jl
              ENDIF
            ENDIF
          ENDDO

          !packing for lproma
          IF (LPROMA > 0) THEN
            DO jl=1,lproma
              i=lwork(jl)
              llwc(jl)               = zlwc(i)
              llwc2(jl)              = lwc_bc(i)
              lrdrad(jl)             = rdrad(i)
              lpress(jl)             = press(i,jk)
              ltemp(jl)              = temp(i,jk)
              lprec(jl)              = rain(i,jk)
              lcov(jl)               = col_cov(i)
              lrho(jl)               = rho(i,jk)
              lgrvol(jl)             = vol(i,jk)
              lph(jl)                = rainph(i,jk)
            ENDDO
            DO jt=1,ntrac
              DO jl=1,lproma
                i=lwork(jl)
                xt_vec(jt,jl)   = xt(i,jk,jt)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                laerflx(jt,jl) = rain_aer(i,jt)

                lresrainnuc(jt,jl) = resrainnuc(i,jt)
                lresrainimp(jt,jl) = resrainimp(i,jt)
              ENDDO
            ENDDO
            DO jt=1,lspec
              do jl=1,lproma
                i=lwork(jl)
                lkppflx(jt,jl) = rain_field(i,jt)
                lresrainchm(jt,jl) = resrainchm(i,jt)
              enddo
            enddo
            DO jt=1,nprod
              DO jl=1,lproma
                i=lwork(jl)
                cprod_array(jl,jt) = prod_array(i,jk,jt,2) * time_step_len
              ENDDO
            ENDDO
            
          ENDIF !for lproma

          !packing for mproma
          IF (mproma > 0) THEN
            DO jt=1,lspec
              DO jl=1,mproma
                i=mwork(jl)
                mkppflx(jt,jl) = rain_field(i,jt)
                mresrainchm(jt,jl) = resrainchm(i,jt)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,mproma
                i=mwork(jl)
                maerflx(jt,jl)  = rain_aer(i,jt)
                mresrainnuc(jt,jl) = resrainnuc(i,jt)
                mresrainimp(jt,jl) = resrainimp(i,jt)
              ENDDO
            ENDDO
            
            DO jt=1,ntrac
              DO jl=1,mproma
                i=mwork(jl)
                mxt_vec(jt,jl) = xt(i,jk,jt)
              ENDDO
            ENDDO
            DO jl=1,mproma
              i=mwork(jl)
              mgrvol(jl)          = vol(i,jk)
            ENDDO
          ENDIF ! for mproma
!----------------------------------------------------------------
!   aerosol scavenging
          if (lscav_aer) then
           
            if (lproma > 0) THEN

              CALL AER_SCAV_DRV(lproma, ntrac, jk, js, 1,             &
                phase, lwork(1:lproma), time_step_len,                &
                xt_vec(:,1:lproma),                                   &
                laerflx(1:aer_count(js)%numbers(c_all),1:lproma),     &
                llwc(1:lproma),   liwc(1:lproma),   lcov(1:lproma),   &
                lrdrad(1:lproma), ltemp(1:lproma),  lpress(1:lproma), &
                lgrvol(1:lproma), lrho(1:lproma),   lprec(1:lproma),  &
                lsnow(1:lproma)  )

              DO jt=1,aer_count(js)%numbers(c_all)
                DO jl=1,lproma
                  i=lwork(jl)
                  lresrainimp(jt,jl) = lresrainimp(jt,jl) + &
                       MAX(laerflx(jt,jl) - rain_aer(i,jt), 0._dp)
                ENDDO
              ENDDO

            END IF    ! lproma > 0

          ENDIF    ! lscav_aer

!--------------------------------------------------------------------

          if (lscav_gas) then
 
            IF (lproma > 0) then
              CALL in_flux(lkppflx(:,1:lproma), xt_vec(:,1:lproma),  &
                           laerflx(:,1:lproma), lgrvol(1:lproma),    &
                           lcov(1:lproma),                           &
                           ntrac, lproma, js, 1, kspec)

              ! Subtract current concentrations of KPP species in order to
              ! calculate increment (below)
              DO jt = 1, lspec
                lresrainchm(jt,1:lproma) = lresrainchm(jt,1:lproma) &
                     - L_SPEC(jt,1:lproma) * lgrvol(1:lproma) * 1.e6_dp &
                     * lcov(1:lproma)
              END DO

              if (lscav_easy) then
                call scav_easy_liq(llwc2(1:lproma), ltemp(1:lproma), lproma, js)
              else
                call scav_aqchem_kpp(ltemp(1:lproma), lpress(1:lproma),    &
                                     llwc2(1:lproma), llwc(1:lproma),      &
                                     lrdrad(1:lproma), time_step_len, 2,   &
                                     lproma, lwork(1:lproma), jk, jrow, js,&
                                     lsteps(1:lproma), lrsteps(1:lproma) )
              endif
              
              call out_flux(lkppflx(:,1:lproma), xt_vec(:,1:lproma),  &
                            lph(1:lproma), llwc(1:lproma),            &
                            lgrvol(1:lproma), lcov(1:lproma),         &
                            ntrac, lproma, js, 1, kspec,              &
                            cprod_array(1:lproma,:))

              ! Calculate increment of KPP species concentrations due to
              ! chemistry
              DO jt = 1, lspec
                lresrainchm(jt,1:lproma) = lresrainchm(jt,1:lproma) &
                     + L_SPEC(jt,1:lproma) * lgrvol(1:lproma) * 1.e6_dp &
                     * lcov(1:lproma)
              END DO

            ENDIF
            
          ENDIF   ! on lscav_gas
!----------------------------------------------------------
          aerfield_evapo(:,:) = 0._dp
          kppfield_evapo(:,:) = 0._dp

          impfield_evapo(:,:) = 0._dp
          IF (ALLOCATED(made3)) THEN
             DO jl = 1, lproma
                i = lwork(jl)
                aerfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                     lresrainnuc(1:aer_count(js)%numbers(c_all),jl)
                kppfield_evapo(1:lspec,i) = lresrainchm(1:lspec,jl)
                impfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                     lresrainimp(1:aer_count(js)%numbers(c_all),jl)
             END DO
             DO jl = 1, mproma
                i = mwork(jl)
                aerfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                     mresrainnuc(1:aer_count(js)%numbers(c_all),jl)
               kppfield_evapo(1:lspec,i) = mresrainchm(1:lspec,jl)
               impfield_evapo(1:aer_count(js)%numbers(c_all),i) = &
                    mresrainimp(1:aer_count(js)%numbers(c_all),jl)
            END DO

            DO idt = 1, SIZE(process(1)%frac_eva, 1)
               DO j = 1, SIZE(process(1)%frac_eva, 2)
                  tmpbuffer(1:kproma,idt,j) = &
                       process(1)%frac_eva(idt,j)%ptr(1:kproma,jk)
               END DO
            END DO

            CALL calc_evapo_frac(status, aerfield_evapo(:,1:kproma), &
                 kppfield_evapo(1:lspec,1:kproma), kproma, lspec, 1, &
                 js, ntrac, tmpbuffer(1:kproma,:,:),    &
                 resimp=&
                 impfield_evapo(1:aer_count(js)%numbers(c_all),1:kproma))
            IF (status .NE. 0) RETURN

            DO idt = 1, SIZE(process(1)%frac_eva, 1)
               DO j = 1, SIZE(process(1)%frac_eva, 2)
                  process(1)%frac_eva(idt,j)%ptr(1:kproma,jk) = &
                       tmpbuffer(1:kproma,idt,j)
               END DO
            END DO

         ELSE
            DO jl=1,lproma
               i =lwork(jl)
               aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
                    laerflx(1:aer_count(js)%numbers(c_all),jl)
               kppfield_evapo(1:lspec,i) = lkppflx(1:lspec,jl)
            ENDDO
            DO jl=1,mproma
               i =mwork(jl)
               aerfield_evapo(1:aer_count(js)%numbers(c_all),i) =  &
                    maerflx(1:aer_count(js)%numbers(c_all),jl)
               kppfield_evapo(1:lspec,i) = mkppflx(1:lspec,jl)
            ENDDO

            IF (n_resmod .GT. 1) THEN

               DO idt = 1, SIZE(process(1)%frac_eva, 1)
                  DO j = 1, SIZE(process(1)%frac_eva, 2)
                     tmpbuffer(1:kproma,idt,j) = &
                          process(1)%frac_eva(idt,j)%ptr(1:kproma,jk)
                  END DO
               END DO

               CALL calc_evapo_frac(status, aerfield_evapo(:,1:kproma), &
                    kppfield_evapo(1:lspec,1:kproma), kproma, lspec, 1, &
                    js, ntrac, tmpbuffer(1:kproma,:,:))

               IF (status .NE. 0) RETURN

               DO idt = 1, SIZE(process(1)%frac_eva, 1)
                  DO j = 1, SIZE(process(1)%frac_eva, 2)
                     process(1)%frac_eva(idt,j)%ptr(1:kproma,jk) = &
                          tmpbuffer(1:kproma,idt,j)
                  END DO
               END DO
            ELSE
               DO idt = 1, ntrac
                  process(1)%frac_eva(idt,1)%ptr(1:kproma,jk) = 1._dp
               END DO
            END IF
         END IF
         ! evaporation
              
         if (mproma > 0) THEN
            DO jl=1,mproma
               i = mwork(jl)
                evfrac(jl) = 0.0_dp
                DO idt = 1, ntrac
                   ! For explanation, see above (process(4))
                   evfrac(jl) = evfrac(jl) + &
                   process(1)%frac_eva(idt,n_resmod)%ptr(i,jk)
                END DO
                evfrac(jl) = evfrac(jl) / REAL(ntrac, kind=dp)
                DO idt = 1, ntrac
                   DO j = 1, n_resmod
                      evfrac_aer(jl,idt,j) = &
                           process(1)%frac_eva(idt,j)%ptr(i,jk)
                   END DO
                END DO
            ENDDO
            IF (ALLOCATED(made3)) evfrac(1:mproma) = 1._dp
            call evapo(mkppflx(:,1:mproma), mxt_vec(:,1:mproma), &
                       mgrvol(1:mproma), evfrac(1:mproma),       &
                       ntrac, mproma, mwork(1:mproma), js, jk,   &
                       jrow, 1, kspec)

            IF (.NOT. ALLOCATED(made3)) THEN

               call evapo_aer(mxt_vec(:,1:mproma), maerflx(:,1:mproma), &
                    mgrvol(1:mproma), evfrac_aer(1:mproma,:,:),         &
                    ntrac, mproma, mwork(1:mproma), js, jk, jrow)
            ELSE
               CALL evapo_aer(mxt_vec(:,1:mproma), mresrainnuc(:,1:mproma), &
                    mgrvol(1:mproma), evfrac_aer(1:mproma,:,:),       &
                    ntrac, mproma, mwork(1:mproma), js, jk, jrow)
               mresrainimp(:,1:mproma) = 0._dp
               mresrainchm(:,1:mproma) = 0._dp
            END IF
          ENDIF

!----------------------------------------------------------
! unpack for liquid water scavenging
          IF( lproma > 0) THEN
            DO jt=1,ntrac
              DO jl=1,lproma
                i=lwork(jl)
                xt(i,jk,jt) = xt_vec(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,lspec
              DO jl=1,lproma
                i=lwork(jl)
                rain_field(i,jt) = lkppflx(jt,jl)
                resrainchm(i,jt) = lresrainchm(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,lproma
                i=lwork(jl)
                rain_aer(i,jt) = laerflx(jt,jl)
                resrainnuc(i,jt) = lresrainnuc(jt,jl)
                resrainimp(i,jt) = lresrainimp(jt,jl)
              ENDDO
            ENDDO
            DO jl=1,lproma
              i=lwork(jl)
              rainpH(i,jk) = lph(jl)
              xsteps_rain(i,jk)  = lsteps(jl)
              xrsteps_rain(i,jk) = lrsteps(jl)
            ENDDO
            DO jt=1,nprod
              DO jl=1,lproma
                i=lwork(jl)
                prod_array(i,jk,jt,2) = cprod_array(jl,jt) / time_step_len
              ENDDO
            ENDDO
          ENDIF
          IF (mproma > 0) THEN
            DO jt=1,ntrac
              DO jl=1,mproma
                i=mwork(jl)
                xt(i,jk,jt) = mxt_vec(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,lspec
              DO jl=1,mproma
                i=mwork(jl)
                rain_field(i,jt) = mkppflx(jt,jl)
                resrainchm(i,jt) = mresrainchm(jt,jl)
              ENDDO
            ENDDO
            DO jt=1,aer_count(js)%numbers(c_all)
              DO jl=1,mproma
                i=mwork(jl)
                rain_aer(i,jt) = maerflx(jt,jl)
                resrainnuc(i,jt) = mresrainnuc(jt,jl)
                resrainimp(i,jt) = mresrainimp(jt,jl)
              ENDDO
            ENDDO
          ENDIF
!----------------------------------
        END IF liq_imp_scav

!----------------------------------
      END if impaction_scavenging

      ! For MADE3, overwrite diagnostic output using appropriate fields
      IF (jk == nlev .AND. ALLOCATED(made3)) THEN

         ! Aerosol species and chemically active species
         DO jt=1,aer_count(js)%numbers(c_all)
            DO jl=1,kproma
               rain_aer(jl,jt) = resrainnuc(jl,jt) + resrainimp(jl,jt)
               snow_aer(jl,jt) = ressnownuc(jl,jt) + ressnowimp(jl,jt) + &
                    rescldinuc(jl,jt)
            END DO
         ENDDO

         ! Chemically inactive species - liquid
         DO jt=1,lspec
            DO jl=1,kproma
               rain_field(jl,jt) = resrainchm(jl,jt)
            END DO
         END DO

         ! Chemically inactive species - ice
         DO jt=1,ispec
            DO jl=1,kproma
               snow_field(jl,jt) = ressnowchm(jl,jt) + rescldichm(jl,jt)
            END DO
         END DO

      ENDIF

    END do level_loop

    status = 0

  END SUBROUTINE SCAV_MAIN
 
!===============================================================================

  INTEGER FUNCTION get_core_mode_made3(i_mode)

    INTRINSIC :: ALLOCATED

    ! I/O
    INTEGER, INTENT(in) :: i_mode

    ! Local
    INTEGER :: idx_core

    IF (.NOT. ALLOCATED(made3)) THEN
      get_core_mode_made3 = -1
      RETURN
    END IF

    IF ((i_mode .EQ. made3(1)%ks) .OR. (i_mode .EQ. made3(1)%as)) THEN
      ! ks, as -> as
      idx_core = 1
    ELSE IF ((i_mode .EQ. made3(1)%km) .OR. (i_mode .EQ. made3(1)%ki) &
         .OR. (i_mode .EQ. made3(1)%am) .OR. (i_mode .EQ. made3(1)%ai)) THEN
      ! km, ki, am, ai -> am
      idx_core = 2
    ELSE IF (i_mode .EQ. made3(1)%cs) THEN
      ! cs -> cs
      idx_core = 3
    ELSE IF ((i_mode .EQ. made3(1)%cm) .OR. (i_mode .EQ. made3(1)%ci)) THEN
      ! cm, ci -> cm
      idx_core = 4
    ELSE
      idx_core = -2
    END IF

    get_core_mode_made3 = idx_core

  END FUNCTION get_core_mode_made3

!=================================================================================================
END MODULE MESSY_SCAV_INTER
