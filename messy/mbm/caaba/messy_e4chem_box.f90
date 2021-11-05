!*****************************************************************************

MODULE messy_e4chem_box

  USE caaba_mem,                  ONLY: model_start_day, model_time, &
                                        timesteplen, init_spec,         &
                                        c, cair, temp, press, relhum,     & 
                                        photrat_channel
  USE messy_cmn_photol_mem        ! IP_MAX, ip_*, jname
  USE messy_e4chem ! ind_*, update_rconst, 
                      ! initialize_indexarrays, SPC_NAMES, EQN_NAMES,
                      ! NSPEC,
  USE messy_main_constants_mem,   ONLY: dp,N_A, rho_H2O, M_H2O, R_gas,     &
                                        OneDay, STRLEN_VLONG,              &
                                        pi, HLINE2
  USE messy_e4chem,               ONLY: modstr

  USE caaba_io,                   ONLY: nf, open_input_file, nf90_inquire, &
                                        nf90_inquire_variable,             &
                                        nf90_get_var, open_output_file,    &
                                        write_output_file, close_file

  IMPLICIT NONE

  REAL(DP) :: jx(IP_MAX) = 0.
  REAL(DP), DIMENSION(:), ALLOCATABLE :: lwc
  ! (for C, cair, and temp, see caaba_mem.f90)

  INTEGER :: ncid_tracer, ncid_spec, ncid_diag ! op_pj_20110217
  ! op_pj_20110217+
  INTEGER, PARAMETER :: NDIAG = 14
  CHARACTER(LEN=32), PARAMETER :: DIAG_NAMES(NDIAG) = &
       (/ 'ZDELTAO3_BRV                    ' , &
          'ZPRODO2                         ' , &
          'ZPRODCO                         ' , &
          'ZPRODCH4                        ' , &
          'ZDESTH12                        ' , &
          'ZDESTH14                        ' , &
          'ZDESTN13                        ' , &
          'ZDESTC1                         ' , &
          'ZDESTCL2O2                      ' , &
          'ZDESTCLOH                       ' , &
          'ZDESTH8                         ' , &
          'ZDESTH4                         ' , &
          'ZDESTH11                        ' , &
          'ZDESTO1                         ' /)
  CHARACTER(LEN=20) :: unit_diag(NDIAG)
  REAL(dp) :: zdiag(NDIAG)
  ! op_pj_20110217+

  REAL(dp), DIMENSION(:), ALLOCATABLE :: DANI, DANIM

  PRIVATE
  PUBLIC :: e4chem_init   ! initialize chemistry
  PUBLIC :: e4chem_physc  ! calculate chemistry
  PUBLIC :: e4chem_result ! print results
  PUBLIC :: e4chem_finish ! close files

CONTAINS

  !***************************************************************************

  SUBROUTINE e4chem_init

    USE messy_main_tools,         ONLY: str, spechum2mr, rel2spechum
    USE messy_e4chem,             ONLY: e4chem_read_nml_ctrl
    USE caaba_mem,                ONLY: degree_lat, &
                                        l_hum_emac, l_psat_liquid, l_relhum_wmo

    IMPLICIT NONE

    INTRINSIC :: TRIM

    CHARACTER(LEN=20)   :: unit(NSPEC)
    INTEGER, PARAMETER  :: iou = 999   ! I/O unit
    INTEGER             :: status ! error status
    INTEGER             :: i,jr
    REAL(DP)            :: latitude(1,1)
    INTEGER             :: ilat(1,1)

    ! read e4chem ctrl namelist:
    CALL e4chem_read_nml_ctrl(status, iou)
    IF (status /= 0) STOP 1

    C(:) = 0. ! default value unless explicitly initialized

    ! mz_rs_20160610+
    !c(ind_H2O) = rh2mr(status,relhum,temp,press,l_hum_emac,l_relhum_wmo) * cair
    c(ind_H2O) = cair * spechum2mr(status, &
      rel2spechum(status, relhum, temp, press, &
      l_hum_emac, l_psat_liquid, l_relhum_wmo))
    ! mz_rs_20160610-
    IF (status /= 0) STOP 1

    CALL x0 ! initial gas phase mixing ratios

    ! ------------------------------------------------------------------------

    ! open output files and write headers:

    ! open output file e4chem_tracer_gp.nc:
    unit(:) = 'mol/mol'    ! define species' units
    CALL open_output_file(ncid_tracer, 'caaba_e4chem', SPC_NAMES, unit)
    ! op_pj_20110217+
    unit_diag(:) = 'molec/s'
    CALL open_output_file(ncid_diag,   'caaba_e4chem_diag' &
         , DIAG_NAMES, unit_diag)
    zdiag(:) = 0.0_dp
    ! op_pj_20110217-

    ALLOCATE(RCGAS(NUMTEM, NUMRAT))
    ALLOCATE(sulook(ksul, 1,1))
    ALLOCATE(DANI(1),DANIM(1))
    DANI(1)=-1.1_dp
    DANIM(1)=-1.1_dp

    ! PRECALCULATE REACTION RATES (LOOKUP TABLE)
    CALL INRCGAS

    latitude(1,1)=degree_lat
! op_pj_20101216+
!!$    CALL inisulnew(latitude,1,1)
    CALL inisulnew(latitude,1,1,1)
! op_pj_20101216-
 
    ! ------------------------------------------------------------------------

  END SUBROUTINE e4chem_init

  !***************************************************************************


  !***************************************************************************

  SUBROUTINE x0

    USE caaba_mem,         ONLY: degree_lat
    USE messy_main_tools,  ONLY: ucase   ! conversion to uppercase

    INTRINSIC :: TRIM

    INTEGER                     :: i, n_var, n_dim, varid_x, ct_spc
    REAL(DP)                    :: mr_x ! mixing ration of species x
    CHARACTER(LEN=STRLEN_VLONG) :: name_x, name_spc

    ! ------------------------------------------------------------------------

    CALL x0_strato
    !CALL x0_mtchem

    ! ------------------------------------------------------------------------

    ! external chemical species' initialization
    IF (TRIM(init_spec)/="") THEN
      CALL open_input_file(ncid_spec, init_spec) ! get ID for input file
      CALL nf(nf90_inquire(ncid_spec, n_dim, n_var)) ! no. dims, vars (=specs)
      !print *, 'mm_box: ncid_spec = ', ncid_spec
      !print *, 'mm_box: no. vars = ', n_var

      varid_x = 1
      DO WHILE (varid_x <= n_var) ! loop over ext init species
        ! get names of variables:
        CALL nf(nf90_inquire_variable(ncid_spec, varid_x, name_x))
        ! convert to uppercase for comparison:
        CALL ucase(name_x)
        !print *, 'mm_box: ext spec no. = ', varid_x,'*'
        !print *, 'mm_box: name = ', name_x,'*'
        ct_spc = 1
        DO WHILE (ct_spc <= NSPEC) ! loop over chemical species
          name_spc = SPC_NAMES(ct_spc)
          ! convert to uppercase for comparison:
          CALL ucase(name_spc)
          !print *, 'mm_box: ct_spc = ',ct_spc
          !print *, 'mm_box: name_spc = ', name_spc,'*'
          !print *, 'mm_box: name_x   = ', name_x,'*'
          IF (TRIM(name_spc) == TRIM(name_x)) THEN
            CALL nf(nf90_get_var(ncid_spec, varid_x, mr_x))
            c(ct_spc) = mr_x * cair
            print *, 'Chemical species initialized:      ', TRIM(name_x)
            !print *, 'mm_box: mr(', ct_spc, ') = ', mr_x
            !print *, 'mm_box: c(', ct_spc, ')  = ', c(ct_spc)
            !ct_spc = ct_spc + 1
            EXIT
          ELSE
            ct_spc = ct_spc + 1
            IF (ct_spc .GT. NSPEC) THEN
              print *, 'Chemical species NOT initialized:  ', TRIM(name_x)
            ENDIF
          ENDIF
        ENDDO ! species loop
        varid_x = varid_x + 1
      ENDDO ! extloop
      CALL close_file(ncid_spec)
    ENDIF

    ! ------------------------------------------------------------------------

    PRINT *, HLINE2
    PRINT *, 'Initial gas-phase mixing ratios and concentrations:'
    PRINT *, HLINE2
    DO i = 1,NSPEC
      IF (c(i)>0) THEN
        WRITE(*,'(2A,ES10.2,A,ES10.2,A)') ' ', SPC_NAMES(i), &
          c(i)/cair, ' mol/mol   = ', c(i), ' mcl/cm3'
      ENDIF
    ENDDO

    ! ------------------------------------------------------------------------

  CONTAINS

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_strato

!!$      ! stratosphere, 20 hPa
!!$      ! from scout02
!!$      IF (ind_H        /= 0) c(ind_H)      =   1.E-12 * cair
!!$      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
!!$      IF (ind_HO2      /= 0) c(ind_HO2)    =   1.E-15 * cair
!!$      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
!!$      IF (ind_NO3      /= 0) c(ind_NO3)    =   1.E-12 * cair
!!$      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
!!$      IF (ind_HNO4     /= 0) c(ind_HNO4)   =   1.E-10 * cair
!!$      IF (ind_CL       /= 0) c(ind_CL)     =   1.E-21 * cair
!!$      IF (ind_CLO      /= 0) c(ind_CLO)    =   1.E-15 * cair
!!$      IF (ind_HOCl     /= 0) c(ind_HOCl)   =  40.E-12 * cair !1.E-15 * cair
!!$      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  =   1.E-13 * cair
!!$      IF (ind_CL2      /= 0) c(ind_CL2)    =   1.E-13 * cair
!!$      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  =   1.E-12 * cair
!!$      IF (ind_N2O      /= 0) c(ind_N2O)    =  1.3E-07 * cair
!!$      IF (ind_CO       /= 0) c(ind_CO)     =  1.4E-08 * cair
!!$      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) =  12.E-12 * cair
!!$      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  =   6.E-10 * cair
!!$      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  =  1.4E-11 * cair
!!$      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) =   1.E-12 * cair
!!$      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  =   1.E-12 * cair
!!$      IF (ind_CCL4     /= 0) c(ind_CCL4)   =   1.E-12 * cair
!!$      IF (ind_CH3CCL3  /= 0) c(ind_CH3CCL3)=   1.E-12 * cair
!!$      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
!!$      IF (ind_H2O      /= 0) c(ind_H2O)    =   1.E-12 * cair
!!$      IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
!!$      IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
!!$      IF (ind_H2       /= 0) c(ind_H2)     =   5.E-07 * cair
!!$      IF (ind_O3       /= 0) c(ind_O3)     =   4.E-06 * cair
!!$      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
!!$      IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
!!$      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.8E-06 * cair
!!$      IF (ind_HCHO     /= 0) c(ind_HCHO)   =   7.E-11 * cair
!!$      IF (ind_CO       /= 0) c(ind_CO)     =  70.E-09 * cair
!!$      IF (ind_CO2      /= 0) c(ind_CO2)    = 350.E-06 * cair
!!$      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 450.E-12 * cair
!!$      IF (ind_HCl      /= 0) c(ind_HCl)    = 400.E-12 * cair
!!$
!!$      IF (ind_OHAB       /= 0) c(ind_OHAB)     =   1.E-16 * cair
!!$      IF (ind_HO2AB      /= 0) c(ind_HO2AB)    =   1.E-15 * cair
!!$   
      ! stratosphere, 10 hPa
      IF (ind_H        /= 0) c(ind_H)      =   1.E-16 * cair
      IF (ind_OH       /= 0) c(ind_OH)     =   1.E-16 * cair
      IF (ind_HO2      /= 0) c(ind_HO2)    =   1.E-15 * cair
      IF (ind_N        /= 0) c(ind_N)      =   1.E-12 * cair
      IF (ind_NO3      /= 0) c(ind_NO3)    =   1.E-12 * cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   =   1.E-10 * cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   =   1.5E-10 * cair
      IF (ind_CL       /= 0) c(ind_CL)     =   1.E-21 * cair
      IF (ind_CLO      /= 0) c(ind_CLO)    =   1.E-15 * cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   =  40.E-12 * cair !1.E-15 * cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  =   1.E-13 * cair
      IF (ind_CL2      /= 0) c(ind_CL2)    =   1.E-13 * cair
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  =   1.E-12 * cair
      IF (ind_N2O      /= 0) c(ind_N2O)    =  1.3E-07 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  1.4E-08 * cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) =  12.E-12 * cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  =   6.E-10 * cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  =  1.4E-11 * cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) =   1.E-12 * cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  =   1.E-12 * cair
      IF (ind_CCL4     /= 0) c(ind_CCL4)   =   1.E-12 * cair
      IF (ind_CH3CCL3  /= 0) c(ind_CH3CCL3)=   1.E-12 * cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   =   5.E-09 * cair
      IF (ind_H2O      /= 0) c(ind_H2O)    = 4.251E-06 * cair
      IF (ind_O3P      /= 0) c(ind_O3P)    =   9.E-34 * cair
      IF (ind_O1D      /= 0) c(ind_O1D)    =   1.E-16 * cair
      IF (ind_H2       /= 0) c(ind_H2)     =   5.E-07 * cair
      IF (ind_O3       /= 0) c(ind_O3)     =   8.E-06 * cair
      IF (ind_NO       /= 0) c(ind_NO)     =   1.E-24 * cair
      IF (ind_NO2      /= 0) c(ind_NO2)    =   1.E-09 * cair
      IF (ind_CH4      /= 0) c(ind_CH4)    =  1.8E-06 * cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   =   7.E-11 * cair
      IF (ind_CO       /= 0) c(ind_CO)     =  70.E-09 * cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 350.E-06 * cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 180.E-12 * cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 400.E-12 * cair

      IF (ind_OHAB       /= 0) c(ind_OHAB)     =   1.E-16 * cair
      IF (ind_HO2AB      /= 0) c(ind_HO2AB)    =   1.E-15 * cair
   
    END SUBROUTINE x0_strato

    ! ------------------------------------------------------------------------

    SUBROUTINE x0_mtchem

      ! mesosphere, 0.01 hPa
      IF (ind_H        /= 0) c(ind_H)      = 1.876E-07* cair
      IF (ind_OH       /= 0) c(ind_OH)     = 1.240E-08* cair
      IF (ind_HO2      /= 0) c(ind_HO2)    = 5.012E-09* cair
      IF (ind_N        /= 0) c(ind_N)      = 2.114E-10* cair
      IF (ind_NO3      /= 0) c(ind_NO3)    = 3.083E-21* cair
      IF (ind_N2O5     /= 0) c(ind_N2O5)   = 9.072E-27* cair
      IF (ind_HNO4     /= 0) c(ind_HNO4)   = 4.754E-18* cair
      IF (ind_CL       /= 0) c(ind_CL)     = 8.001E-11* cair
      IF (ind_CLO      /= 0) c(ind_CLO)    = 1.564E-13* cair
      IF (ind_HOCl     /= 0) c(ind_HOCl)   = 7.015E-15* cair
      IF (ind_CL2O2    /= 0) c(ind_CL2O2)  = 1.558E-25* cair
      IF (ind_CL2      /= 0) c(ind_CL2)    = 1.E-13* cair ! ????
      IF (ind_CH3O2    /= 0) c(ind_CH3O2)  = 1.E-12* cair! ????
      IF (ind_N2O      /= 0) c(ind_N2O)    = 7.077E-11* cair
      IF (ind_CO       /= 0) c(ind_CO)     = 9.862E-07* cair
      IF (ind_CH3OOH   /= 0) c(ind_CH3OOH) = 3.558E-13* cair
      IF (ind_ClNO3    /= 0) c(ind_ClNO3)  = 1.460E-22* cair
      IF (ind_CFCl3    /= 0) c(ind_CFCl3)  = 2.854E-38* cair
      IF (ind_CF2Cl2   /= 0) c(ind_CF2Cl2) = 2.776E-17* cair
      IF (ind_CH3CL    /= 0) c(ind_CH3CL)  = 2.212E-14* cair
!!$   IF (ind_CCL4     /= 0) c(ind_CCL4)   = 5.605E-45* cair
      IF (ind_CH3CCL3  /= 0) c(ind_CH3CCL3)= 1.401E-44* cair
      IF (ind_HNO3     /= 0) c(ind_HNO3)   = 2.084E-13* cair
      IF (ind_H2O      /= 0) c(ind_H2O)    = 4.567E-06* cair
      IF (ind_O3P      /= 0) c(ind_O3P)    = 5.350E-06* cair
      IF (ind_O1D      /= 0) c(ind_O1D)    = 3.143E-14* cair
      IF (ind_H2       /= 0) c(ind_H2)     = 1.310E-06* cair
      IF (ind_O3       /= 0) c(ind_O3)     = 7.823E-08* cair
      IF (ind_NO       /= 0) c(ind_NO)     = 2.454E-09* cair
      IF (ind_NO2      /= 0) c(ind_NO2)    = 1.685E-12* cair
      IF (ind_CH4      /= 0) c(ind_CH4)    = 1.113E-07* cair
      IF (ind_HCHO     /= 0) c(ind_HCHO)   = 1.417E-12* cair
      IF (ind_CO       /= 0) c(ind_CO)     = 9.862E-07* cair
      IF (ind_CO2      /= 0) c(ind_CO2)    = 3.641E-04* cair
      IF (ind_H2O2     /= 0) c(ind_H2O2)   = 1.169E-10* cair
      IF (ind_HCl      /= 0) c(ind_HCl)    = 3.342E-09* cair
    END SUBROUTINE x0_mtchem

  END SUBROUTINE x0

  !***************************************************************************

  SUBROUTINE e4chem_physc

    USE messy_e4chem,             ONLY: CHEMICS
    USE messy_main_tools,         ONLY: spechum2mr, rel2spechum
    ! um_ak_20110711+    
    !USE messy_jval,               ONLY: jval_gp
    USE messy_jval,               ONLY: jval_2d
    ! um_ak_20110711-    
    USE caaba_mem,                ONLY: l_hum_emac, l_psat_liquid, l_relhum_wmo, &
                                        photol_clev, degree_lat, &
         month => lmonth, &
         cossza



    REAL(DP) :: henry(0:NSPEC)=0.     ! Inverse Henry constant, dimensionless
    INTEGER :: ip, status

    ! tracer mixing ratio before chemistry
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zmr
    ! zmr after tracer scaling
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zmrbc
    ! tracer mixing ratio after chemistry
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zmrac

    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: temp_2D
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: pressi_2D
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: pressm_2D

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Conc
    ! specific humidity = m(H2O)/m(air) [kg/kg]
    REAL(DP), DIMENSION(:,:), ALLOCATABLE   :: sphum_2D

! op_pj_20101216+
!!$    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: tp_i
    REAL(DP), DIMENSION(:), ALLOCATABLE :: tp_i
! op_pj_20101216-
    ! op_pj_20110217+
    ! ozone molecules destroyed by BrO+ClO cycle
    REAL(DP), ALLOCATABLE :: ZDELTAO3_BRV(:,:)
    ! ozone budget
    REAL(DP), ALLOCATABLE :: ZPRODO2(:,:)
    REAL(DP), ALLOCATABLE :: ZPRODCO(:,:)
    REAL(DP), ALLOCATABLE :: ZPRODCH4(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTH12(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTH14(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTN13(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTC1(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTCL2O2(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTCLOH(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTH8(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTH4(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTH11(:,:)
    REAL(DP), ALLOCATABLE :: ZDESTO1(:,:)
    ! op_pj_20110217-

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE  :: JX

    REAL(DP) :: zlat(1)

    REAL(dp), PARAMETER :: U0LIM = pi/2.+0.052335956 ! 3deg past sunset
    REAL(dp) :: RAMU02 

    INTEGER :: KLEV, NLON, KLON
    KLEV=1
    NLON=1
    KLON=1

    ALLOCATE(pressi_2D(NLON,KLEV+1))
    ALLOCATE(pressm_2D(NLON,KLEV))
    ALLOCATE(temp_2D(NLON,KLEV))
    ALLOCATE(sphum_2D(NLON,KLEV))
    ALLOCATE(Conc(NLON,KLEV,NSPEC))
    ALLOCATE(JX(NLON,KLEV,IP_MAX))
! op_pj_20101216+
!!$ ALLOCATE(tp_i(NLON,KLEV))
    ALLOCATE(tp_i(NLON))
! op_pj_20101216-
    ! op_pj_20110217+
    ALLOCATE(ZDELTAO3_BRV(NLON,KLEV))
    ! ozone budget
    ALLOCATE(ZPRODO2(NLON,KLEV))
    ALLOCATE(ZPRODCO(NLON,KLEV))
    ALLOCATE(ZPRODCH4(NLON,KLEV))
    ALLOCATE(ZDESTH12(NLON,KLEV))
    ALLOCATE(ZDESTH14(NLON,KLEV))
    ALLOCATE(ZDESTN13(NLON,KLEV))
    ALLOCATE(ZDESTC1(NLON,KLEV))
    ALLOCATE(ZDESTCL2O2(NLON,KLEV))
    ALLOCATE(ZDESTCLOH(NLON,KLEV))
    ALLOCATE(ZDESTH8(NLON,KLEV))
    ALLOCATE(ZDESTH4(NLON,KLEV))
    ALLOCATE(ZDESTH11(NLON,KLEV))
    ALLOCATE(ZDESTO1(NLON,KLEV))
    ! op_pj_20110217-

    ! mz_rs_20160610+
    !c(ind_H2O) = rh2mr(status,relhum,temp,press,l_hum_emac,l_relhum_wmo) * cair
    c(ind_H2O) = cair * spechum2mr(status, &
      rel2spechum(status, relhum, temp, press, &
      l_hum_emac, l_psat_liquid, l_relhum_wmo))
    ! mz_rs_20160610-
    IF (status /= 0) STOP 1

    DO ip=1, IP_MAX
       ! um_ak_20110711+
       !jx(1,1,ip) = jval_gp(ip)%ptr(1,photol_clev,1)
       jx(1,1,ip) = jval_2d(ip)%ptr(1,photol_clev)
       ! um_ak_20110711-
    ENDDO

    RAMU02=ACOS(cossza)
    IF (RAMU02.le.U0LIM) THEN
       DANI(1) = 1.1_dp
    ELSE
       DANI(1) = -1.1_dp
       ! sunset
       IF (DANIM(1).GT.1.) DANI(1) = 0.1_dp
    ENDIF
  
    CALL check_range('before chem:',c(:))
    c(:) = MAX(c(:),0._DP) ! set negative values to zero
    conc(1,1,:) = c(:)
    temp_2D(1,1) = temp
    pressi_2D(1,1) = press-1./100.*press
    pressi_2D(1,2) = press+1./100.*press
    pressm_2D(1,1) = press
    sphum_2D(1,1) = rel2spechum(status,relhum,temp,press)
! op_pj_20101216+
!!$ tp_i(:,:)=2._dp ! --> 2: always above tropopause
    tp_i(:)=2._dp ! --> 2: always above tropopause
! op_pj_20101216-
    zlat(:)=degree_lat/180._dp*pi
    ! JX: CFC11 => CFCl3, CFC12 => CF2Cl2

! op_pj_20101216+
!!$ CALL CHEMICS(1,    1,        1,            1,   
    CALL CHEMICS(1,    1,        1,                                                &
! op_pj_20101216-
         temp_2D(:,:),      sphum_2D(:,:),    Conc(:,:,:),                           & 
         pressm_2D(:,:),    pressi_2D(:,:),   DANI(:), DANIM(:),                     &
         JX(:,:,ip_O3P),    JX(:,:,ip_O1D),   JX(:,:,ip_NO2),   JX(:,:,ip_HNO3),     &
         JX(:,:,ip_COH2),   JX(:,:,ip_CHOH),  JX(:,:,ip_N2O5),  JX(:,:,ip_HNO4),     &
         JX(:,:,ip_NO2O),   JX(:,:,ip_NOO2),  JX(:,:,ip_H2O2),  JX(:,:,ip_CH3OOH),   &
         JX(:,:,ip_O2),     JX(:,:,ip_CFCl3), JX(:,:,ip_CF2Cl2), JX(:,:,ip_N2O),     &
         JX(:,:,ip_CLONO2),                                                          &
         JX(:,:,ip_CL2O2),  JX(:,:,ip_HOCL),                                         &
         JX(:,:,ip_CCL4),   JX(:,:,ip_CH3CL), JX(:,:,ip_CH3CCL3), JX(:,:,ip_HCL),    &
         JX(:,:,ip_H2O),    JX(:,:,ip_NO),    JX(:,:,ip_CO2),                        &
         zlat,   MONTH,     timesteplen,      tp_i    &  
         ! op_pj_20110217+ ozone diagnostics/budget
         , ZDELTAO3_BRV(:,:)                          &
         , ZPRODO2(:,:),   ZPRODCO(:,:)               &
         , ZPRODCH4(:,:),  ZDESTH12(:,:)              &
         , ZDESTH14(:,:),  ZDESTN13(:,:)              &
         , ZDESTC1(:,:),   ZDESTCL2O2(:,:)            &
         , ZDESTCLOH(:,:), ZDESTH8(:,:)               &
         , ZDESTH4(:,:),   ZDESTH11(:,:)              &
         , ZDESTO1(:,:) )

    zdiag( 1) = ZDELTAO3_BRV(1,1)
    zdiag( 2) = ZPRODO2(1,1)
    zdiag( 3) = ZPRODCO(1,1)
    zdiag( 4) = ZPRODCH4(1,1)
    zdiag( 5) = ZDESTH12(1,1)
    zdiag( 6) = ZDESTH14(1,1)
    zdiag( 7) = ZDESTN13(1,1)
    zdiag( 8) = ZDESTC1(1,1)
    zdiag( 9) = ZDESTCL2O2(1,1)
    zdiag(10) = ZDESTCLOH(1,1)
    zdiag(11) = ZDESTH8(1,1)
    zdiag(12) = ZDESTH4(1,1)
    zdiag(13) = ZDESTH11(1,1)
    zdiag(14) = ZDESTO1(1,1)
    ! op_pj_20110217-

    c(:) = Conc(1,1,:)          ! remove two dimensions

    DANIM(1) = DANI(1)

    ! DEALLOCATE MEMORY
    DEALLOCATE(Conc)
    DEALLOCATE(temp_2D,pressi_2D,pressm_2D)
    DEALLOCATE(JX)
    DEALLOCATE(sphum_2D)
    DEALLOCATE(tp_i)
    ! op_pj_20110217+
    DEALLOCATE(ZDELTAO3_BRV)
    ! ozone budget
    DEALLOCATE(ZPRODO2)
    DEALLOCATE(ZPRODCO)
    DEALLOCATE(ZPRODCH4)
    DEALLOCATE(ZDESTH12)
    DEALLOCATE(ZDESTH14)
    DEALLOCATE(ZDESTN13)
    DEALLOCATE(ZDESTC1)
    DEALLOCATE(ZDESTCL2O2)
    DEALLOCATE(ZDESTCLOH)
    DEALLOCATE(ZDESTH8)
    DEALLOCATE(ZDESTH4)
    DEALLOCATE(ZDESTH11)
    DEALLOCATE(ZDESTO1)
    ! op_pj_20110217-


  CONTAINS

    !-------------------------------------------------------------------------

    SUBROUTINE check_range(infostring,conc)

      ! print a warning if a concentration is not in the correct range

      CHARACTER(LEN=*), INTENT(IN) :: infostring
      REAL(DP),         INTENT(IN) :: conc(:) ! tracer concentration
      INTEGER :: jt

      INTRINSIC :: SIZE

      tracer_loop: DO jt=1,SIZE(conc)
         wrong_conc: IF ((conc(jt)<0._DP).OR.(conc(jt)>cair)) THEN
            WRITE(*,'(2A,F10.0,A,1PG12.3E3,2A)') infostring, &
                 ' time =', model_time, &
                 ', c =', conc(jt), ' mcl/cm3 for ', TRIM(SPC_NAMES(jt))
         ENDIF wrong_conc
      ENDDO tracer_loop

    END SUBROUTINE check_range

    !-------------------------------------------------------------------------

  END SUBROUTINE e4chem_physc

  !***************************************************************************

  SUBROUTINE e4chem_result

    CALL write_output_file(ncid_tracer, model_time, c/cair)
    CALL write_output_file(ncid_diag, model_time, zdiag) ! op_pj_20110217

  END SUBROUTINE e4chem_result

  !***************************************************************************

  SUBROUTINE e4chem_finish

    CALL close_file(ncid_tracer)
    CALL close_file(ncid_diag)   ! op_pj_20110217
    DEALLOCATE(RCGAS, sulook)
    DEALLOCATE(DANI,DANIM)

  END SUBROUTINE e4chem_finish

  !***************************************************************************

END MODULE messy_e4chem_box

!*****************************************************************************
