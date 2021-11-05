PROGRAM MESSY_SCAV_BOX
!     This routine is intended for small testing purposes of the 
!     main SCAV routines.
!     Author Holger Tost, May 2005 , MPI-Chemie, Mainz, Germany

  USE messy_main_constants_mem, ONLY: g, M_air, R_gas, avo => N_A
  USE messy_scav_l_kpp,         ONLY: dp, nspec, IND_HNO3, &
                                      IND_no3m_l, IND_HNO3_l,&
                                      str_field_kpp_l => spc_names
  USE messy_scav,               ONLY: calc_lwc_bc, scav_aqchem_kpp
  USE messy_scav_mem
  USE messy_scav_inter,         ONLY: in_flux, out_flux, in_cloud, &
                                      out_cloud, evapo


  IMPLICIT NONE

  INTEGER, PARAMETER  :: nlev  = 10
  INTEGER, PARAMETER  :: ntrac = 21
  INTEGER, PARAMETER  :: js = 1
  
  INTEGER, PARAMETER  :: idt_HNO3    = 1
  INTEGER, PARAMETER  :: idt_H2O2    = 2
  INTEGER, PARAMETER  :: idt_O3      = 3
  INTEGER, PARAMETER  :: idt_SO2     = 4
  INTEGER, PARAMETER  :: idt_H2SO4   = 5
  INTEGER, PARAMETER  :: idt_HCHO    = 6
  INTEGER, PARAMETER  :: idt_HCOOH   = 7
  INTEGER, PARAMETER  :: idt_NH3     = 8
  INTEGER, PARAMETER  :: idt_N2O5    = 9
  INTEGER, PARAMETER  :: idt_CH3OOH  = 10
  INTEGER, PARAMETER  :: idt_CH3CO2H = 11
  INTEGER, PARAMETER  :: idt_OH      = 12
  INTEGER, PARAMETER  :: idt_HO2     = 13
  INTEGER, PARAMETER  :: idt_PAN     = 14
  INTEGER, PARAMETER  :: idt_HNO4    = 15
  INTEGER, PARAMETER  :: idt_HONO    = 16
  INTEGER, PARAMETER  :: idt_CO2     = 17
  INTEGER, PARAMETER  :: idt_CH3O2   = 18
  INTEGER, PARAMETER  :: idt_CH3OH   = 19
  INTEGER, PARAMETER  :: idt_CH3COCH3= 20
  INTEGER, PARAMETER  :: idt_O2      = 21
  INTEGER, PARAMETER  :: aero_count = 1
  INTEGER, PARAMETER  :: kproma = 1

  CHARACTER(LEN=30)   :: str_field(21)

  INTEGER             :: nsteps
  INTEGER             :: t, jk, jt, jl, lproma
  INTEGER             :: kcltop

  REAL(dp) :: TEMPI(NLEV+1,kproma), PRESSI(NLEV+1,kproma), zi(NLEV+1,kproma)
  REAL(dp) :: press(NLEV,kproma), TEMP(NLEV,kproma), RHO(NLEV,kproma)
  REAL(dp) :: Z(NLEV,kproma), DELTA_Z, GRVOL(NLEV,kproma), DELTA_P(NLEV,kproma)
  REAL(dp) :: LWC(NLEV,kproma), IWC(NLEV,kproma), COVER(NLEV,kproma)
  REAL(dp) :: RAINPROD(NLEV,kproma), SNOWPROD(NLEV,kproma)
  REAL(dp) :: RAINFLUX(NLEV,kproma), SNOWFLUX(NLEV,kproma)
  REAL(dp) :: lwc_bc(NLEV,kproma), lwc_bc_eff(NLEV,kproma), washfrac(NLEV,kproma)
  REAL(dp) :: mc(NLEV,kproma), CM(NLEV,kproma)
  REAL(dp) :: ZRDRAD(kproma), Hp_l(kproma)
  REAL(dp) :: time_step_len
  REAL(dp) :: area
  REAL(dp) :: pxtp1(nlev,ntrac,kproma)
  REAL(dp) :: wetflx(nspec,kproma)
  REAL(dp) :: dummy_field(nspec,kproma)
  REAL(dp) :: ph_rain(nlev,kproma), ph_cloud(nlev,kproma), Hp_cloud(NLEV,kproma)
  REAL(dp), ALLOCATABLE :: aero_flux(:,:), l_aerflx(:,:)

  REAL(dp) :: l_wetflx(nspec,kproma), l_xt(ntrac,kproma)
  REAL(dp) :: l_temp(kproma), l_press(kproma), l_lwc(kproma), l_rainprod(kproma)
  REAL(dp) :: l_ph(kproma), l_hp(kproma), l_rad(kproma), l_rho(kproma)
  REAL(dp) :: l_dummy(1:nspec,kproma), l_cov(kproma), l_vol(kproma)
  REAL(dp) :: l_lwc1(kproma), l_lwc2(kproma), l_hp_cloud(kproma), lsteps(1:kproma)
  INTEGER  :: lidx(kproma)
  REAL(dp) :: l_steps(kproma), r_steps(kproma)
! op_pj_20130424+
!!$  REAL(dp) :: lprod_array(kproma,nprod)
  REAL(dp), ALLOCATABLE :: lprod_array(:,:)
! op_pj_20130424-  

  ALLOCATE(aero_flux(aero_count,kproma))
  ALLOCATE(l_aerflx(aero_count,kproma))
  nsteps = 10
  area = 2.5e9_dp  ! in m²
!  time_step_len= 1._dp    ! in s
!  time_step_len= 60._dp    ! in s
!  time_step_len= 600._dp    ! in s
  time_step_len= 3600._dp    ! in s
!  time_step_len= 43200._dp    ! in s
!  time_step_len= 10000._dp    ! in s

  dummy_field(:,:) = 0._dp
  lidx(:) = 1
  CALL INITIALIZE

!    built a column grid of 10 levels
!    level 1 is top and 10 bottom layer
  DELTA_Z  = 500._dp

  do jk=1,nlev+1
    do jl=1,kproma
      zi(jk,jl)        = (nlev + 1 - jk) * DELTA_Z
      tempi(jk,jl)     = 288._dp - (nlev-jk) * 6.8_dp                ! K
      pressi(jk,jl)    = 101300._dp * exp(-zi(jk,jl) / &
                        (R_gas * tempi(jk,jl)/(g*1.e-3_dp*M_air)))
    enddo
  enddo
  do jk=1,nlev
    do jl=1,kproma
      temp(jk,jl)     = (tempi(jk+1,jl) - tempi(jk,jl))/2._dp + tempi(jk,jl)
      delta_p(jk,jl)  = pressi(jk+1,jl) - pressi(jk,jl)
      press(jk,jl)    = delta_p(jk,jl)/2._dp + pressi(jk,jl)
      rho(jk,jl)      = press(jk,jl) * M_air * &
                        1.e-3_dp/ (R_gas*temp(jk,jl))         ! kg/m³
      grvol(jk,jl)    = area * DELTA_Z                        ! m³
      print*, jk, temp(jk,jl), press(jk,jl), delta_p(jk,jl), rho(jk,jl)
          ! convert mixing ratio [mol/mol] to concentration [mcl/cc]
          ! mc = mixing ratio to concentration
      mc(jk,jl) = (avo/1.e6_dp) * press(jk,jl) / (R_gas * temp(jk,jl))
          ! convert concentration [mcl/cc] to mixing ratio [mol/mol]
          ! cm = concentration to mixing ratio
      cm(jk,jl)= 1._dp/mc(jk,jl)
    enddo
  enddo
! set rain and cloud parameters
  
  zrdrad(:)       = 0._dp
  lwc(:,:)        = 0._dp
  iwc(:,:)        = 0._dp
  lwc_bc(:,:)     = 0._dp 
  lwc_bc_eff(:,:) = 0._dp
  rainprod(:,:)   = 0._dp
  rainflux(:,:)   = 0._dp
  snowflux(:,:)   = 0._dp

  lwc(6,1)       = 1.e-3_dp
  rainprod(6,1)  = 0.5e-3_dp
  lwc(7,1)       = 2.e-3_dp
  rainprod(7,1)  = 1.e-3_dp
  lwc(8,1)       = 2.e-3_dp
  rainprod(8,1)  = 1.e-3_dp
  cover(:,1)      = 1._dp
  do jl=2,kproma
     lwc(6,jl)       = 1.e-3_dp + 1.e-3_dp*(jl-1)/(kproma-1)        ! kg/m³
     rainprod(6,jl)  = 0.5e-3_dp + 0.5e-3_dp*(jl-1)/(kproma-1)       ! kg/m³
     lwc(7,jl)       = 2.e-3_dp + 1.e-3_dp*(jl-1)/(kproma-1)        ! kg/m³
     rainprod(7,jl)  = 1.e-3_dp + 0.5e-3_dp*(jl-1)/(kproma-1)       ! kg/m³
     lwc(8,jl)       = 2.e-3_dp + 1.e-3_dp*(jl-1)/(kproma-1)        ! kg/m³
     rainprod(8,jl)  = 1.e-3_dp + 0.5e-3_dp*(jl-1)/(kproma-1)       ! kg/m³
     cover(:,jl)      = 1._dp
  end do
 do jl=1,kproma

    do jk=2,nlev
      rainflux(jk,jl) = rainflux(jk-1,jl) + &
                        rainprod(jk-1,jl) * delta_p(jk-1,jl) / &
                        (g * time_step_len)
    enddo

!!$    rainflux(8,1)  = 4._dp/3600._dp
!!$    rainflux(9,1)  = 4.5_dp/3600._dp
!!$    rainflux(10,1) = 4.5_dp/3600._dp
!!$    rainflux(8,2)  = 4.1_dp/3600._dp
!!$    rainflux(9,2)  = 4.25_dp/3600._dp
!!$    rainflux(10,2) = 4.35_dp/3600._dp

  enddo
  ph_rain(:,:)    = 19.99_dp
  ph_cloud(:,:)   = 19.99_dp
  aero_flux(:,:) = 0._dp
  do jl=1,kproma
    CALL INIT_CONC(pxtp1(1:nlev,1:ntrac,jl), ntrac, nlev)

    do jk=1,nlev
      write(50+jl,'(21e25.15)' ) pxtp1(jk,1:ntrac,jl)
    enddo
  enddo
!    loop over number of timesteps
  do t=1,nsteps

    do jk=1,nlev
      do jt=1,ntrac
        do jl=1,kproma
          pxtp1(jk,jt,jl) = pxtp1(jk,jt,jl) * mc(jk,jl)
        enddo
      enddo
    enddo

    wetflx(1:nspec,:) = 0._dp

!    loop over levels
    do jk=1,nlev
      do jl=1,kproma
!    if rainflux above threshold calc imp_scav      
        call calc_lwc_bc(rainflux(jk,jl), snowflux(jk,jl), cover(jk,jl),     &
                       lwc_bc_eff(jk,jl), washfrac(jk,jl), zrdrad(jl),       &
                       nlev, jk, kcltop)

        lwc_bc(jk,jl) = rainflux(jk,jl) / cover(jk,jl) * &
                        time_step_len * area / grvol(jk,jl)
!!$
!!$!      aerosol_impaction_scavenging
!!$      call aer_scav()
!!$

        print*, "column:", jl, "level: ", jk, "LWC: ",lwc(jk,jl),"RAINFLUX: ", &
              rainflux(jk,jl), "RAIN_RATE: ", rainflux(jk,jl) * 3600._dp,      &
              "LWC_BC: ", lwc_bc(jk,jl), "RADIUS: ",zrdrad(jl)
      enddo

!      pack for local
      l_wetflx(:,:) = 0._dp
      l_xt(:,:)     = 0._dp
      l_aerflx(:,:) = 0._dp
      l_vol(:)      = 0._dp
      l_cov(:)      = 0._dp
      l_temp(:)     = 0._dp
      l_press(:)    = 0._dp
      l_lwc1(:)     = 0._dp
      l_lwc2(:)     = 0._dp
      l_rad(:)      = 0._dp
      l_ph(:)       = 0._dp
      l_steps(:)    = 0._dp
      r_steps(:)    = 0._dp
      lprod_array(:,:) = 0._dp
      lproma = 0
      do jl=1,kproma
        if (lwc_bc(jk,jl) >= 1.e-7_dp) then
          lproma = lproma + 1
          l_wetflx(1:nspec,lproma) = wetflx(1:nspec,jl)
          l_xt(1:ntrac,lproma)     = pxtp1(jk,1:ntrac,jl)
          l_aerflx(:,lproma)       = aero_flux(:,jl)
          l_vol(lproma)            = grvol(jk,jl)
          l_cov(lproma)            = cover(jk,jl)
          l_temp(lproma)           = temp(jk,jl)
          l_press(lproma)          = press(jk,jl)
          l_lwc1(lproma)           = lwc_bc_eff(jk,jl)
          l_lwc2(lproma)           = lwc_bc(jk,jl)
          l_rad(lproma)            = zrdrad(jl)
          l_ph(lproma)             = ph_rain(jk,jl)
        endif
      enddo
      if (lproma > 0) then

        call in_flux(l_wetflx(1:nspec,1:lproma), l_xt(1:ntrac,1:lproma),&
                     l_aerflx(:,1:lproma), l_vol(1:lproma),             &
                     l_cov(1:lproma), ntrac, lproma, js, 1, nspec)
        call scav_aqchem_kpp(l_temp(1:lproma), l_press(1:lproma), &
                             l_lwc1(1:lproma), l_lwc2(1:lproma), &
                             l_rad(1:lproma), time_step_len, 2, &
                             lproma, lidx, jk, 1, js,           &
                             l_steps(1:lproma), r_steps(1:lproma) )
        call out_flux(l_wetflx(1:nspec,1:lproma), l_xt(1:ntrac,1:lproma), &
                      l_pH(1:lproma), l_lwc2(1:lproma), &
                      l_vol(1:lproma), l_cov(1:lproma), ntrac, lproma, js,&
                      1, nspec, lprod_array(1:lproma,:) )

        print*,"steps: ", l_steps, r_steps
      endif
      lproma = 0
      do jl=1,kproma
        if (lwc_bc(jk,jl) >= 1.e-7_dp) then
          lproma = lproma + 1
          wetflx(1:nspec,jl)   = l_wetflx(1:nspec,lproma) 
          pxtp1(jk,1:ntrac,jl) = l_xt(1:ntrac,lproma)  
          ph_rain(jk,jl)       = l_ph(lproma)     
        endif
      enddo
!!$
!!$!    if lwc above threshold calc nuc_scav
!!$      call aer_scav_cloud()
!!$
      l_wetflx(:,:) = 0._dp
      l_dummy(:,:)  = 0._dp
      l_xt(:,:)     = 0._dp
      l_aerflx(:,:) = 0._dp
      l_vol(:)      = 0._dp
      l_cov(:)      = 0._dp
      l_temp(:)     = 0._dp
      l_press(:)    = 0._dp
      l_lwc(:)      = 0._dp
      l_rainprod(:) = 0._dp
      l_rad(:)      = 0._dp
      l_ph(:)       = 0._dp
      l_rho(:)      = 0._dp
      l_steps(:)    = 0._dp
      r_steps(:)    = 0._dp
      lprod_array(:,:) = 0._dp
      lproma = 0
      do jl=1,kproma
        if (lwc(jk,jl) >= 1.e-7_dp) then
          lproma = lproma + 1
          l_wetflx(1:nspec,lproma) = wetflx(1:nspec,jl)
          l_dummy(:,lproma)        = 0._dp
          l_xt(1:ntrac,lproma)     = pxtp1(jk,1:ntrac,jl)
          l_aerflx(:,lproma)       = aero_flux(:,jl)
          l_vol(lproma)            = grvol(jk,jl)
          l_cov(lproma)            = cover(jk,jl)
          l_temp(lproma)           = temp(jk,jl)
          l_press(lproma)          = press(jk,jl)
          l_lwc(lproma)            = lwc(jk,jl)
          l_rainprod(lproma)       = rainprod(jk,jl)
          l_rad(lproma)            = 1.75e-2_dp
          l_ph(lproma)             = ph_cloud(jk,jl)
          l_rho(lproma)            = rho(jk,jl)
        endif
      enddo
      
      if (lproma > 0) then

        call in_cloud(l_dummy(1:nspec,1:lproma), l_xt(1:ntrac,1:lproma),&
                      l_aerflx(:,1:lproma), l_vol(1:lproma), l_cov(1:lproma), &
                      ntrac, lproma, js, 1, nspec)
        

        call scav_aqchem_kpp(l_temp(1:lproma), l_press(1:lproma), &
                             l_lwc(1:lproma), l_lwc(1:lproma), l_rad(1:lproma),&
                             time_step_len, 1, lproma, lidx, jk, 1, js, &
                             l_steps(1:lproma), r_steps(1:lproma) )

        call out_cloud(l_dummy(1:nspec,1:lproma), l_wetflx(1:nspec,1:lproma),&
                       l_xt(1:ntrac,1:lproma), l_pH(1:lproma), l_hp(1:lproma),&
                       l_cov(1:lproma), l_lwc(1:lproma), l_rainprod(1:lproma),&
                       l_vol(1:lproma), ntrac, lproma, js, 1, nspec, &
                       lprod_array(1:lproma,:) )
        l_Hp_cloud(1:lproma) = l_hp(1:lproma) /(1.0e-6_dp * l_rho(1:lproma))

      !  call evapo(l_dummy(1:nspec,1:lproma), l_xt(1:ntrac,1:lproma), &
      !             l_vol(1:lproma), l_cov(1:lproma), ntrac, lproma, js)
      endif
      
      lproma = 0
      do jl=1,kproma
        if (lwc(jk,jl) >= 1.e-7_dp) then
          lproma = lproma + 1
          wetflx(1:nspec,jl)   = l_wetflx(1:nspec,lproma) 
          pxtp1(jk,1:ntrac,jl) = l_xt(1:ntrac,lproma) 
          ph_cloud(jk,jl)      = l_ph(lproma)  
          Hp_cloud(jk,jl)      = l_Hp_cloud(lproma)
        endif
      enddo

   enddo

!!$!     output (ASCII only)
    do jk=1,nlev
      do jt=1,ntrac
        do jl=1,kproma
          pxtp1(jk,jt,jl) = pxtp1(jk,jt,jl) * cm(jk,jl)
        enddo
      enddo
    enddo
    do jk=1,nlev
      do jl=1,kproma
        write(50+jl,'(21e25.15)') pxtp1(jk,1:ntrac,jl)
      enddo
    enddo
    
    do jl=1,kproma
      write(20+jl,'(21e25.15)') wetflx(1:nspec,jl)/(area*time_step_len)
    enddo
    print*, "rain pH:" ,ph_rain
    
  enddo

!----------------------------------
CONTAINS
!-------------------------------------------------------------------------------------
  SUBROUTINE INITIALIZE

    USE messy_scav_liq,          ONLY: use_schwartz, ALLOC_SCAV_VALUES_L
    USE messy_scav_l_kpp,        only: INITVAL => initialize_l, rtol, atol, icntrl
    USE messy_scav_mem,          ONLY: cpl_aerosol
    USE messy_main_tools,        ONLY: match_wild, strcrack  ! op_pj_20130424
    
    IMPLICIT NONE

    INTEGER           :: i, j, ji, count_g, count_l, kprod, dummy
    INTEGER           :: scav2trac_idx(nspec)
    CHARACTER(LEN=30) :: string
    REAL(dp)          :: mw_trac(21)
    CHARACTER(LEN=26), POINTER     :: strname(:) => NULL()
    REAL(dp)          :: rtols

    USE_schwartz=.true.
!    USE_schwartz=.false.
    cpl_aerosol = 0
    CALL INITVAL  
! set solver
    ICNTRL(3) = 2               ! =ROS3
!   ICNTRL(3) = 4               ! =RODAS3

!    CALL STR_KPP(str_field_kpp_L)
    CALL INDEX_FIELD
    CALL ALLOC_SCAV_VALUES_L
    CALL DEFINE_MW(mw_trac)
    ! accuracy
    rtols = 1.e-2_dp
    DO i=1,nspec
       rtol(i) = rtols  ! must be defined for ros3
       atol(i) = 1.e-3_dp   ! must be defined for ros3
    END DO

    scav2trac_idx(:) = 0
    count_g = 0
    count_l = 0
    
    print*,"Gas - Liquid - Interaction"

    do jt=1,ntrac
      do ji=1,nspec
        if ( TRIM(str_field_kpp_l(ji)) == TRIM(str_field(jt)) ) then
          scav2trac_idx(ji)=jt
          print*, scav2trac_idx(ji),str_field(jt),str_field_kpp_l(ji)
          endif
      enddo
    enddo
    do jt=1,nspec
      if (scav2trac_idx(jt) /= 0) then
        count_g = count_g + 1
      else
        count_l = count_l + 1
      end if
    enddo
    lspec_gas = count_g
    lspec_liq = count_l
      
    ALLOCATE(kpp_l_idx(1))
    ALLOCATE(kpp_l_idx(js)%gas_spec(lspec_gas, gas_max))
    ALLOCATE(kpp_l_idx(js)%liq_spec(lspec_liq, liq_max))
    ALLOCATE(kpp_l_idx(js)%gas_attr(lspec_gas, gatt_max))
    i = 0 
    j = 0

    do jt=1,nspec
      if (scav2trac_idx(jt)  /= 0 ) then
        i = i + 1
        kpp_l_idx(js)%gas_spec(i,gas_idx)  = jt 
        kpp_l_idx(js)%gas_spec(i,gas2trac) = scav2trac_idx(jt)
        kpp_l_idx(js)%gas_attr(i,gas_mw)   = mw_trac(scav2trac_idx(jt))
      else
        j = j + 1
        kpp_l_idx(js)%liq_spec(j,liq_idx) = jt
      end if
    end do
    print*, "nspecs", nspec, lspec_gas, lspec_liq 
    do jt=1,lspec_gas
      print*, "Species ", kpp_l_idx(js)%gas_spec(jt,gas_idx), &
        str_field_kpp_l(kpp_l_idx(js)%gas_spec(jt,gas_idx)), &
              " interacts with tracer ", kpp_l_idx(js)%gas_spec(jt,gas2trac),&
              str_field(kpp_l_idx(js)%gas_spec(jt,gas2trac))
    enddo
 
    kpp_l_idx(js)%liq_spec(:,liq2evap) = 0
    do j=1,lspec_liq
      do jt=1,ntrac
        string=TRIM(str_field(jt))//'_l'
        if ( TRIM(str_field_kpp_l(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
          == TRIM(string) )    kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
      enddo
    enddo
    ALLOCATE(aer_count(1))
    ALLOCATE(aer_count(js)%numbers(aercount_max))
    aer_count(js)%numbers(:) = 0
    ALLOCATE(L_spec(0:nspec,kproma))

    ! op_pj_20130424+
    nprod = 0
    DO ji=1,nspec
      IF (MATCH_WILD( '*Prod*', str_field_kpp_l(ji) )) THEN
        nprod = nprod + 1
      END IF
    END DO
    
    ALLOCATE(PR(js))

    ALLOCATE(PR(js)%PROD(nprod))
    ALLOCATE(PR(js)%PROD_IDX(nprod))
    ALLOCATE(PR(js)%PROD_CHAR(nprod))

    PR(js)%PROD_CHAR(:) = ""
    PR(js)%PROD_IDX(:)  = 0

    kprod = 0
    DO jt=1,nspec
      IF (ASSOCIATED(strname)) DEALLOCATE (strname)
      NULLIFY(strname)
      !        print*, "in specloop ", jt, str_field_kpp_l(jt)
      IF (MATCH_WILD( '*Prod*', str_field_kpp_l(jt) )) THEN
        kprod = kprod + 1
        CALL strcrack(str_field_kpp_l(jt), '_', strname, dummy)
!          print*, "correctly found a prod rate ", kprod, strname, strname(3), dummy
        PR(js)%PROD_CHAR(kprod) = strname(2)
        PR(js)%PROD_IDX(kprod) = jt
      END IF
    END DO
    ! op_mm_20131106+
!!$    ALLOCATE(lprod_array(kproma,nprod))
    ALLOCATE(lprod_array(kproma,0:nprod))
    ! op_mm_20131106-
    ! op_pj_20130424-


  END SUBROUTINE INITIALIZE
!--------------------------------------------------------------------------------------
  SUBROUTINE index_field

    str_field(1)='HNO3'
    str_field(2)='H2O2'
    str_field(3)='O3'
    str_field(4)='SO2'
    str_field(5)='H2SO4'
    str_field(6)='HCHO'
    str_field(7)='HCOOH'
    str_field(8)='NH3'
    str_field(9)='N2O5'
    str_field(10)='CH3OOH'
    str_field(11)='CH3CO2H'
    str_field(12)='OH'
    str_field(13)='HO2'
    str_field(14)='PAN'
    str_field(15)='HNO4'
    str_field(16)='HONO'
    str_field(17)='CO2'
    str_field(18)='CH3O2'
    str_field(19)='CH3OH'
    str_field(20)='CH3COCH3'
    str_field(21)='O2'
    
  END SUBROUTINE index_field
!--------------------------------------------------------------------------------------
 
!--------------------------------------------------------------------------------------
  SUBROUTINE INIT_CONC(pxtp1, ntrac, nlev)
    
    IMPLICIT NONE

    INTEGER,INTENT(IN)      :: ntrac, nlev
    REAL(dp), INTENT(INOUT) :: pxtp1(NLEV, NTRAC)

    INTEGER :: jk, jt
! in mol/mol
    pxtp1(:,:) = 0._dp

    pxtp1(1:nlev,idt_HNO3)     = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_HNO3)     = 1.e-1_dp *  pxtp1(1:nlev,idt_HNO3) 
!    pxtp1(1:nlev,idt_HNO3)     = (/1.e-9_dp,2.e-9_dp,3.e-9_dp,4.e-9_dp,&
!                                       5.e-9_dp,6.e-9_dp,7.e-9_dp,8.e-9_dp,9.e-9_dp,1.e-8_dp/)
    pxtp1(1:nlev,idt_H2O2)     = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_O3)       = (/1.e-8_dp,1.e-8_dp,1.e-8_dp,1.e-8_dp,1.e-8_dp,1.e-8_dp,1.e-8_dp,1.e-8_dp,1.e-8_dp,2.e-8_dp/)
    pxtp1(1:nlev,idt_SO2)      = &
      (/1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-11_dp,1.e-8_dp,1.e-8_dp,1.e-8_dp,1.e-8_dp/)
    pxtp1(1:nlev,idt_H2SO4)    = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-8_dp,1.e-8_dp/)
    pxtp1(1:nlev,idt_HCHO)     = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_HCOOH)    = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_NH3)      = &
      (/1.e-10_dp,1.e-10_dp,1.e-10_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-10_dp,1.e-10_dp/)
    pxtp1(1:nlev,idt_N2O5)     = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_CH3OOH)   = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_CH3CO2H)  = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_OH)       = &
      (/1.e-14_dp,1.e-14_dp,1.e-14_dp,1.e-14_dp,1.e-14_dp,1.e-14_dp,1.e-14_dp,1.e-14_dp,1.e-14_dp,1.e-14_dp/)
    pxtp1(1:nlev,idt_HO2)      = &
     (/1.e-13_dp,1.e-13_dp,1.e-13_dp,1.e-13_dp,1.e-13_dp,1.e-13_dp,1.e-13_dp,1.e-13_dp,1.e-13_dp,1.e-13_dp/)
    pxtp1(1:nlev,idt_O2)       = (/ 0.78_dp, 0.78_dp, 0.78_dp, 0.78_dp, 0.78_dp, 0.78_dp, 0.78_dp, 0.78_dp, 0.78_dp, 0.78_dp/)
    pxtp1(1:nlev,idt_HNO4)     = &
      (/1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp,1.e-12_dp/)
    pxtp1(1:nlev,idt_HONO)     = &
      (/1.e-10_dp,1.e-10_dp,1.e-10_dp,1.e-10_dp,1.e-10_dp,1.e-10_dp,1.e-10_dp,1.e-10_dp,1.e-10_dp,1.e-10_dp/)
    pxtp1(1:nlev,idt_CO2)      = &
      (/3.5e-4_dp,3.5e-4_dp,3.5e-4_dp,3.5e-4_dp,3.5e-4_dp,3.5e-4_dp,3.5e-4_dp,3.5e-4_dp,3.5e-4_dp,3.5e-4_dp/)
    pxtp1(1:nlev,idt_CH3O2)    = &
      (/1.e-11_dp,1.e-11_dp,1.e-11_dp,1.e-11_dp,1.e-11_dp,1.e-11_dp,1.e-11_dp,1.e-11_dp,1.e-11_dp,1.e-11_dp/)
    pxtp1(1:nlev,idt_CH3OH)    = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_CH3COCH3) = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)
    pxtp1(1:nlev,idt_PAN)      = (/1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp,1.e-9_dp/)

    
  END SUBROUTINE INIT_CONC
  
!--------------------------------------------------------------------------------------
  SUBROUTINE DEFINE_MW(mw)

    IMPLICIT NONE
    INTEGER, PARAMETER  :: idt_HNO3    = 1
    INTEGER, PARAMETER  :: idt_H2O2    = 2
    INTEGER, PARAMETER  :: idt_O3      = 3
    INTEGER, PARAMETER  :: idt_SO2     = 4
    INTEGER, PARAMETER  :: idt_H2SO4   = 5
    INTEGER, PARAMETER  :: idt_HCHO    = 6
    INTEGER, PARAMETER  :: idt_HCOOH   = 7
    INTEGER, PARAMETER  :: idt_NH3     = 8
    INTEGER, PARAMETER  :: idt_N2O5    = 9
    INTEGER, PARAMETER  :: idt_CH3OOH  = 10
    INTEGER, PARAMETER  :: idt_CH3CO2H = 11
    INTEGER, PARAMETER  :: idt_OH      = 12
    INTEGER, PARAMETER  :: idt_HO2     = 13
    INTEGER, PARAMETER  :: idt_PAN     = 14
    INTEGER, PARAMETER  :: idt_HNO4    = 15
    INTEGER, PARAMETER  :: idt_HONO    = 16
    INTEGER, PARAMETER  :: idt_CO2     = 17
    INTEGER, PARAMETER  :: idt_CH3O2   = 18
    INTEGER, PARAMETER  :: idt_CH3OH   = 19
    INTEGER, PARAMETER  :: idt_CH3COCH3= 20
    INTEGER, PARAMETER  :: idt_O2      = 21
    REAL(dp) :: mw(21)
    mw(idt_O2)      = 3.2E-2_DP
    mw(idt_O3)      = 4.8E-2_DP
    
    mw(idt_OH)      = 1.7E-2_DP
    mw(idt_HO2)     = 3.3E-2_DP
    mw(idt_H2O2)    = 3.4E-2_DP
    
    mw(idt_NH3)     = 1.7E-2_DP
    mw(idt_HONO)    = 4.7E-2_DP
    mw(idt_HNO3)    = 6.3E-2_DP
    mw(idt_N2O5)    = 1.08E-1_DP
    mw(idt_HNO4)    = 7.9E-2_DP
    
    mw(idt_HCHO)    = 3.E-2_DP
    mw(idt_HCOOH)   = 4.6E-2_DP
    mw(idt_CH3OH)   = 3.2E-2_DP
    mw(idt_CH3OOH)  = 4.8E-2_DP
    mw(idt_CO2)     = 4.4E-2_DP
    mw(idt_CH3O2)   = 4.7E-2_DP
    
    mw(idt_CH3CO2H) = 6.E-2_DP
    mw(idt_PAN)     = 1.21E-1_DP
    
    mw(idt_CH3COCH3)= 5.8E-2_DP      
    
    
    mw(idt_SO2)     = 6.4E-2_DP
    mw(idt_H2SO4)   = 9.8E-2_DP

  END SUBROUTINE DEFINE_MW

!-----------------------------------------------------------------------------

END PROGRAM MESSY_SCAV_BOX










