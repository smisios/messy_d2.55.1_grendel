MODULE MESSY_GMXE_AERCHEM

  USE MESSY_GMXE_MEM

  IMPLICIT NONE

CONTAINS


!=============================== KPP LIQUID chemistry ========================

  SUBROUTINE AERCHEM_AQCHEM_KPP (PTEMP, PPRESS, AQRLWC,         &
                                 AQRDRAD, TIME_STEP_LEN, LPROMA, &
                                 LIDX, JK)

    ! NEW SUBROUTINE IN WHICH THE AQUEOUS CHEMISTRY IS BEING
    ! CALCULATED PRODUCED BY KPP 
    USE MESSY_MAIN_CONSTANTS_MEM, ONLY: R_GAS, N_A
    USE MESSY_GMXE_AERCHEM_LIQ
! KPP VARIABLES
    USE MESSY_GMXE_AERCHEM_KPP

    IMPLICIT NONE
    SAVE
    
    INTEGER,  INTENT(IN) :: LPROMA, JK, LIDX(LPROMA)
    REAL(DP), INTENT(IN) :: TIME_STEP_LEN
    REAL(DP), INTENT(IN) :: PTEMP(LPROMA), AQRLWC(LPROMA), &
                            AQRDRAD(LPROMA), PPRESS(LPROMA)

    REAL(dp) :: NSTEPS(lproma), RSTEPS(lproma)

    INTEGER  :: I, JL, IDX2
    LOGICAL  :: LOC(LPROMA)
    CHARACTER(LEN=*), PARAMETER :: MODSTR ='GMXe_aerchem'

    REAL(dp) :: cair(lproma)
    REAL(dp) :: conc(lproma,1:nspec)

    INTEGER  :: kpp_stat(lproma), kpp_steps(lproma), kpp_rsteps(lproma)
    INTEGER  :: status

!!$    REAL(DP) :: N_BEFORE, N_AFTER
!!$    REAL(DP) :: HNO3_BEFORE, HNO3_AFTER
!!$    REAL(DP) :: NO3M_BEFORE, NO3M_AFTER
!!$    REAL(DP) :: HNO3L_BEFORE, HNO3L_AFTER
!!$    REAL(DP) :: N2O5_BEFORE, N2O5_AFTER

    INTRINSIC :: REAL, SUM
    
    CALL ALLOC_CHEM_FIELDS(LPROMA)
! SETTING UP PHYSICAL PARAMETERS FOR EACH BOX
    CALL LIQ_PHYSC_INIT(AQRLWC, AQRDRAD, PTEMP, PPRESS)

    CALL SET_PHOTO(LIDX, JK)
! CALCULATES THE EQUILIBRIUM COEFFICIENTS / TESTFACS
    CALL EQUIL_DISSOCIATION_COEFF
! CALCULATES THE TRANSFER COEFFICIENTS
    CALL TRANSFER_COEFFICIENT_LIQ

    conc(:,:) = 0._dp
    cair(:) = (N_A/1.E6_dp) * ppress(:) / (R_gas*ptemp(:))
    IF (.NOT.ivec) THEN       ! solver is not created with kp4 and uses a non-vector format
       DO jl=1,lproma

          CALL fill_jx (status, jx_vec(jl:jl,:))
          IF (status /= 0) STOP "fill_jx array size"
          
          CALL fill_temp (status,ptemp(jl:jl))
          IF (status /= 0) STOP "fill_temp array size"
          
          CALL fill_press (status, ppress(jl:jl))
          IF (status /= 0) STOP "fill_press array size"
          
          CALL fill_cair (status, cair(jl:jl))
          IF (status /= 0) STOP "fill_cair array size"
          
          CALL fill_lwc (status, vliq_lwc(jl:jl))
          IF (status /= 0) STOP "fill_lwc array size"
          
          CALL fill_cv_l (status, vcv_l(jl:jl))
          IF (status /= 0) STOP "fill_cv_l"
          
          CALL fill_k_exf(status, vk_exf(jl:jl,1:nspec))
          IF (status /= 0) STOP "fill_k_exf array size"
          CALL fill_k_exb(status, vk_exb(jl:jl,1:nspec))
          IF (status /= 0) STOP "fill_k_exb array size"
          CALL fill_k_exf_N2O5(status, vk_exf_N2O5(jl:jl))
          IF (status /= 0) STOP "fill_k_exf_N2O5 array size"
          CALL fill_k_exf_ClNO3(status, vk_exf_ClNO3(jl:jl))
          IF (status /= 0) STOP "fill_k_exf_ClNO3 array size"
          CALL fill_k_exf_BrNO3(status, vk_exf_BrNO3(jl:jl))
          IF (status /= 0) STOP "fill_k_exf_BrNO3 array size"

          conc(1,1:NSPEC) = L_SPEC(1:NSPEC,JL)
    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
          CALL kpp_integrate(TIME_STEP_LEN, conc, ierrf=kpp_stat, xNacc=kpp_steps, xNrej=kpp_rsteps)
    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
          conc(:,:) = MAX(conc(:,:), 0.0_DP) ! force positive results
          L_SPEC(1:NSPEC,JL) = conc(1,1:NSPEC)

          nsteps(jl) = kpp_steps(1)
          rsteps(jl) = kpp_rsteps(1)
       END DO

    ELSE  ! solver uses a vectorised format
       CALL fill_jx (status, jx_vec)
       IF (status /= 0) STOP "fill_jx array size"

       CALL fill_temp (status,ptemp(1:LPROMA))
       IF (status /= 0) STOP "fill_temp array size"
       
       CALL fill_press (status, ppress)
       IF (status /= 0) STOP "fill_press array size"
       
       CALL fill_cair (status, cair)
       IF (status /= 0) STOP "fill_cair array size"
       
       CALL fill_lwc (status, vliq_LWC)
       IF (status /= 0) STOP "fill_lwc array size"
       
       CALL fill_cv_l (status, VCV_L)
       IF (status /= 0) STOP "fill_cv_l"
       
       CALL fill_k_exf(status, vk_exf(:,1:NSPEC))
       IF (status /= 0) STOP "fill_k_exf array size"
       CALL fill_k_exb(status, vk_exb(:,1:NSPEC))
       IF (status /= 0) STOP "fill_k_exb array size"
       CALL fill_k_exf_N2O5(status, vk_exf_N2O5)
       IF (status /= 0) STOP "fill_k_exf_N2O5 array size"
       CALL fill_k_exf_ClNO3(status, vk_exf_ClNO3)
       IF (status /= 0) STOP "fill_k_exf_ClNO3 array size"
       CALL fill_k_exf_BrNO3(status, vk_exf_BrNO3)
       IF (status /= 0) STOP "fill_k_exf_BrNO3 array size"

       DO jl=1,lproma
          conc(jl,1:NSPEC) = L_SPEC(1:NSPEC,JL)
       ENDDO

    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
       !        CALL kpp_integrate(TIME_STEP_LEN, conc)
       CALL kpp_integrate(TIME_STEP_LEN, conc, ierrf=kpp_stat, xNacc=kpp_steps, xNrej=kpp_rsteps)
    ! KPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPPKPP
       conc(:,:) = MAX(conc(:,:), 0.0_DP) ! force positive results
       DO jl=1,lproma
          L_SPEC(1:NSPEC,JL) = conc(jl,1:NSPEC)
          IF (kpp_stat(jl) /=1 ) print*, "Error in kpp in aerchem! ", &
               kpp_stat(jl), "Box: ", jl
       ENDDO
       nsteps(:) = kpp_steps
       rsteps(:) = kpp_rsteps

    END IF


!!$    ! mz_ht_20080918+
!!$    ! special for Hg
!!$    DO jl=1,lproma
!!$      L_SPEC(IND_RGM_l,jl)  =  L_SPEC(IND_RGM_l,jl) + L_SPEC(IND_Hg_l,jl)
!!$      L_SPEC(IND_Hg_l,jl)   = 0._dp
!!$    ENDDO
!!$    ! mz_ht_20080918-
    CALL DEALLOC_CHEM_FIELDS

    RETURN

  END SUBROUTINE AERCHEM_AQCHEM_KPP
!==============================================================================
  SUBROUTINE AERCHEM_DRIVER(kproma, nlev, jm, iwat, ncomp_aer,      &
                            xaerml, prwet, temp, press,             &
                            time_step_len)!,                        &
                            !pH, klwc, &
                            !prod_1, prod_2, prod_3,                 &
                            !prod_4, prod_5, prod_6, prod_7, prod_8, &
                            !prod_9, prod_10, prod_11, prod_12)
                            

    USE MESSY_MAIN_CONSTANTS_MEM,       ONLY: dp, AVO => N_A
    USE MESSY_GMXE_AERCHEM_LIQ,         ONLY: diag_aerchem, lspec, l_spec,  &
                                              nprod
    USE MESSY_GMXE_AERCHEM_KPP,         ONLY: IND_HP_L
    USE MESSY_GMXE_AERCHEM_INP_KPP,     ONLY: thres_lwc

!!$    USE MESSY_GMXE_AERCHEM_KPP_G_MEM,   ONLY: KPP_Prod_01_l, &
!!$                                              KPP_Prod_02_l, KPP_Prod_03_l, &
!!$                                              KPP_Prod_04_l, KPP_Prod_05_l, &
!!$                                              KPP_Prod_06_l, KPP_Prod_07_l, &
!!$                                              KPP_Prod_08_l, KPP_Prod_09_l, &
!!$                                              KPP_Prod_10_l, KPP_Prod_11_l, &
!!$                                              KPP_Prod_12_l
    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: iwat         ! aerosol water index in this mode
    INTEGER,  INTENT(IN) :: ncomp_aer    ! number of compounds of paerml 
    INTEGER,  INTENT(IN) :: kproma, nlev, jm
    REAL(dp), INTENT(IN) :: prwet(kproma,nlev) !aerosol radius of this mode [m]
    REAL(dp), INTENT(IN) :: temp(kproma,nlev)  !air temperature [K]
    REAL(dp), INTENT(IN) :: press(kproma,nlev) !air pressure [Pa]

    REAL(dp), INTENT(IN) :: time_step_len

!    !liquid fraction of compounds
!    REAL(dp), INTENT(IN)   :: liq_frac(kproma,nlev,ncomp_aer) 
    !tracer array 
    REAL(dp), INTENT(INOUT):: xaerml(kproma,nlev,0:ncomp_aer)   
!!$    REAL(dp), INTENT(INOUT):: pH(kproma,nlev)   
!!$    REAL(dp), INTENT(INOUT):: klwc(kproma,nlev)
!!$    REAL(dp), INTENT(INOUT):: prod_1(kproma,nlev), prod_2(kproma,nlev),&
!!$                              prod_3(kproma,nlev), prod_4(kproma,nlev),&
!!$                              prod_5(kproma,nlev), prod_6(kproma,nlev),&
!!$                              prod_7(kproma,nlev), prod_8(kproma,nlev),&
!!$                              prod_9(kproma,nlev), prod_10(kproma,nlev),&
!!$                              prod_11(kproma,nlev), prod_12(kproma,nlev)
    ! local
    INTEGER :: lproma
    INTEGER :: jl, jk, jt, idx1, idx2, idx

    REAL(dp) :: lwc(kproma)
    REAL(dp) :: lrdrad(kproma), ltemp(kproma), lpress(kproma), llwc(kproma)
    INTEGER  :: lidx(kproma)
    REAL(dp) :: zfw2a
    
    ALLOCATE (L_SPEC(0:lspec,kproma))

!!$    pH(:,:) = -1.0e34_dp
!!$    prod_1(:,:) = 0._dp
!!$    prod_2(:,:) = 0._dp
!!$    prod_3(:,:) = 0._dp
!!$    prod_4(:,:) = 0._dp
!!$    prod_5(:,:) = 0._dp
!!$    prod_6(:,:) = 0._dp
!!$    prod_7(:,:) = 0._dp
!!$    prod_8(:,:) = 0._dp
!!$    prod_9(:,:) = 0._dp
!!$    prod_10(:,:) = 0._dp
!!$    prod_11(:,:) = 0._dp
!!$    prod_12(:,:) = 0._dp
    diag_aerchem(jm)%ph_2d(:,:) = -1.0e34_dp

    level_loop: DO jk=1,nlev

      lproma = 0
      L_SPEC(:,:) = 0._dp
      packing_loop_0: DO jl=1,kproma
        ! lwc in kg/m^3
        ! xaerml is in molecules/cm^3 -> *1/avo   => mol/cm^3
        !                             -> *1/1e-6  => mol/m^3
        !                             -> *1/55.5 = *1/(1000/M_h2o) => kg/m^3
        lwc(jl)   = xaerml(jl,jk,iwat) / (55.5_dp * avo * 1e-6_dp)
!!$        klwc(jl,jk) = lwc(jl)
        diag_aerchem(jm)%lwc_2d(jl,jk) = lwc(jl)
        ! calculate aerosol chemistry only if sufficient liquid water available
        IF ( (prwet(jl,jk) > 1.e-10_dp) .AND. &
             (lwc(jl) > thres_lwc) ) THEN

          lproma = lproma + 1
          lidx(lproma)  = jl
        ENDIF
      END DO packing_loop_0
      packing_loop_1: DO jl=1,lproma
        idx = lidx(jl)
        lrdrad(jl) = prwet(idx,jk)
        ltemp(jl)  = temp(idx,jk)
        lpress(jl) = press(idx,jk)
        llwc(jl)   = lwc(idx)
      ENDDO packing_loop_1
      species_packing_loop_liquid:  DO jt=1, spec_number
        idx1 = species(jt)%aermlidx(jm)
        IF (idx1 == 0) CYCLE
        idx2 = species(jt)%kppidx(jm)
        DO jl=1,lproma
          idx = lidx(jl)
          ! use only the liquid fraction for the partitioning
          !            l_xt(lproma,jt) = xaerml(jl,jk,jt) * liq_frac(jl,jk,jt)
!          L_SPEC(idx2,jl) = MAX(0._dp, xaerml(idx,jk,idx1))
          L_SPEC(idx2,jl) = xaerml(idx,jk,idx1)
        ENDDO
      END DO species_packing_loop_liquid
      species_packing_loop_gas:  DO jt=1, spec_number
        idx1 = species(jt)%aermlidx(0)
        IF (idx1 == 0) CYCLE
        idx2 = species(jt)%kppidx(0)
        DO jl=1,lproma
          idx = lidx(jl)
          ! use only the liquid fraction for the partitioning
          !            l_xt(lproma,jt) = xaerml(jl,jk,jt) * liq_frac(jl,jk,jt)
!          L_SPEC(idx2,jl) = MAX(0._dp,xaerml(idx,jk,idx1))
          L_SPEC(idx2,jl) = xaerml(idx,jk,idx1)
        ENDDO
      END DO species_packing_loop_gas

     
      CALL AERCHEM_AQCHEM_KPP(LTEMP(1:lproma),     lPRESS(1:lproma),  &
                              llwc(1:lproma),      lrdrad(1:lproma),  &
                              TIME_STEP_LEN, LPROMA, LIDX(1:lproma), JK)
      
      unpacking_loop_liquid: DO jt=1, spec_number
        idx1 = species(jt)%aermlidx(jm)
        IF (idx1 == 0) CYCLE
        idx2 = species(jt)%kppidx(jm)
        IF (idx2 == 0) CYCLE
        DO jl=1,lproma
          ! add the new liquid fraction to the solid fraction
!!$!          xaerml(lidx(jl),jk,jt) = (1._dp - liq_frac(jl,jk,jt)) *   &
!!$!                                    xaerml(lidx(jl),jk,jt) + l_xt(jl,jt)
          xaerml(lidx(jl),jk,idx1) = L_SPEC(idx2,jl)
        END DO
      END DO unpacking_loop_liquid
      unpacking_loop_gas: DO jt=1, spec_number
        idx1 = species(jt)%aermlidx(0)
        IF (idx1 == 0) CYCLE
        idx2 = species(jt)%kppidx(0)
        IF (idx2 == 0) CYCLE
        DO jl=1,lproma
          ! add the new liquid fraction to the solid fraction
!!$!          xaerml(lidx(jl),jk,jt) = (1._dp - liq_frac(jl,jk,jt)) *   &
!!$!                                    xaerml(lidx(jl),jk,jt) + l_xt(jl,jt)
          xaerml(lidx(jl),jk,idx1) = L_SPEC(idx2,jl)
        END DO
      END DO unpacking_loop_gas
      DO jl=1,lproma
        idx = lidx(jl)
        zfw2a   = llwc(jl) * avo * 1.e-6_dp

        diag_aerchem(jm)%ph_2d(idx,jk) = -LOG10 ( L_SPEC(IND_Hp_l, jl) / zfw2a ) 
      END DO
!!$      DO jl=1,lproma
!!$        idx = lidx(jl)
!!$        prod_1(idx,jk)  = L_SPEC(KPP_Prod_01_l,jl) / time_step_len
!!$        prod_2(idx,jk)  = L_SPEC(KPP_Prod_02_l,jl) / time_step_len
!!$        prod_3(idx,jk)  = L_SPEC(KPP_Prod_03_l,jl) / time_step_len
!!$        prod_4(idx,jk)  = L_SPEC(KPP_Prod_04_l,jl) / time_step_len
!!$        prod_5(idx,jk)  = L_SPEC(KPP_Prod_05_l,jl) / time_step_len
!!$        prod_6(idx,jk)  = L_SPEC(KPP_Prod_06_l,jl) / time_step_len
!!$        prod_7(idx,jk)  = L_SPEC(KPP_Prod_07_l,jl) / time_step_len
!!$        prod_8(idx,jk)  = L_SPEC(KPP_Prod_08_l,jl) / time_step_len
!!$        prod_9(idx,jk)  = L_SPEC(KPP_Prod_09_l,jl) / time_step_len
!!$        prod_10(idx,jk) = L_SPEC(KPP_Prod_10_l,jl) / time_step_len
!!$        prod_11(idx,jk) = L_SPEC(KPP_Prod_11_l,jl) / time_step_len
!!$        prod_12(idx,jk) = L_SPEC(KPP_Prod_12_l,jl) / time_step_len
!!$      END DO

      DO jt=1,nprod
        idx2 = diag_aerchem(jm)%PROD_IDX(jt)
        DO jl=1,lproma
          idx = lidx(jl)
          diag_aerchem(jm)%Prod(jt)%ptr_2d(idx,jk) = &
            L_SPEC(idx2,jl) / time_step_len
        END DO
      END DO
      

    END DO level_loop
    IF (ALLOCATED(L_SPEC)) DEALLOCATE(L_SPEC)


  END SUBROUTINE AERCHEM_DRIVER
!==============================================================================

  SUBROUTINE INIT_AERCHEM

    USE MESSY_MAIN_TOOLS,             ONLY: strcrack, match_wild
    USE MESSY_GMXE_AERCHEM_LIQ
!    USE MESSY_GMXE_AERCHEM_AUTO_KPP
    USE messy_gmxe_AERCHEM_kpp,       ONLY: INITVAL => INITIALIZE, ICNTRL
    USE MESSY_GMXE_AERCHEM_KPP,       ONLY: rtol, atol, str_field_kpp_l => spc_names

    INTEGER :: jt, count_g, count_l, dummy, status, IDX
    CHARACTER(LEN=26), POINTER     :: strname(:) => NULL()
    REAL(dp)  :: rtols

    LSPEC_GAS = 0
    LSPEC_LIQ = 0
    
    DO jt=1,lspec
      if (associated(strname)) DEALLOCATE (strname)
      NULLIFY(strname)     
      call strcrack(STR_FIELD_KPP_L(jt),'_', strname, dummy)
      IF (dummy==1) THEN 
        LSPEC_GAS = LSPEC_GAS + 1
      ELSE
        LSPEC_LIQ = LSPEC_LIQ + 1
      ENDIF
    ENDDO

    ALLOCATE(kpp_l_idx%gas_spec(lspec_gas,gas_max))
    ALLOCATE(kpp_l_idx%gas_attr(lspec_gas,GATT_MAX))
    ALLOCATE(kpp_l_idx%liq_spec(lspec_liq,liq_max))
    ALLOCATE(kpp_l_idx%liq_attr(lspec_liq,LATT_MAX))
    count_g = 0
    count_l = 0

    DO jt=1,lspec
      dummy=0
      if (associated(strname)) DEALLOCATE (strname)
      NULLIFY(strname)     
      call strcrack(STR_FIELD_KPP_L(jt),'_', strname, dummy)
      IF (dummy==1) THEN 
        count_g = count_g + 1
        kpp_l_idx%gas_spec(count_g,gas_idx) = jt
      ELSE
        count_l = count_l + 1
        kpp_l_idx%liq_spec(count_l,liq_idx) = jt

      ENDIF
    ENDDO
    
    CALL INITVAL
    ! set solver
    ICNTRL(3) = 2               ! =ROS3
!    ICNTRL(3) = 3               ! =ROS4
 !   ICNTRL(3) = 4               ! =RODAS3
    ! solver accuracy for liquid kpp
    rtols = 5.e-2_dp
    DO jt=1,LSPEC
       rtol(jt) = rtols      ! must be defined for ros3
       atol(jt) = 1.e2_dp    ! must be defined for ros3
    END DO
    ICNTRL(4) = 1e7    ! maximum number of steps

    nprod = 0
    DO jt=1,lspec
      IF (MATCH_WILD( '*Prod*', str_field_kpp_l(jt) )) THEN
        nprod = nprod + 1
      END IF
    END DO

  END SUBROUTINE INIT_AERCHEM

!-------------------------------------------------------------------------------

END MODULE MESSY_GMXE_AERCHEM
