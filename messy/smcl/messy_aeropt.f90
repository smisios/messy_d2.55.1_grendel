MODULE MESSY_AEROPT

  USE MESSY_AEROPT_MEM

  REAL(dp), PARAMETER :: MIN_VAL       = 1.e-10_dp
  REAL(dp), PARAMETER :: MIN_VAL_OMEGA = 0.9999_dp

CONTAINS

!===============================================================================

  SUBROUTINE aeropt_read_nml_ctrl(status, iou)

    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    SAVE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    !--- Local variables:

    CHARACTER(LEN=*), PARAMETER :: substr = 'aeropt_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE

    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    !--- 1) Read CTRL namelist:

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR


  END SUBROUTINE aeropt_read_nml_ctrl
!===============================================================================
  SUBROUTINE aeropt_initialise(L_IO, IDX)

    USE MESSY_AEROPT_INPUT, ONLY: lookup_initialize_rad0, &
                                  lookup_initialize_rad1, &
                                  lookup_initialize_rad2, &
                                  lookup_initialize_rad3
    USE MESSY_AEROPT_SETS,  ONLY: allocate_lut_struct, &
                                  map_lut_struct

    LOGICAL, INTENT(IN) :: L_IO
    INTEGER, INTENT(IN) :: IDX

    CALL ALLOCATE_LUT_STRUCT(IDX)

    CALL MAP_LUT_STRUCT(IDX)

    ! get dimension of fields in NetCDF files

    CALL lookup_initialize_rad0(l_io, IDX)
    CALL MAP_LUT_STRUCT(IDX)

    ! allocate memory for lookup tables
    CALL lookup_initialize_rad1(idx)
    CALL MAP_LUT_STRUCT(IDX)

    ! read NetCDF files
    CALL lookup_initialize_rad2(l_io)
    CALL MAP_LUT_STRUCT(IDX)

    ! finish rad initialization
    CALL lookup_initialize_rad3(l_io)
    CALL MAP_LUT_STRUCT(IDX)


  END SUBROUTINE aeropt_initialise

!===============================================================================

  SUBROUTINE calc_aero_properties(kproma, klev, paernl, paernl2, prwet, &
                                  ppressi, aot_opt,                     &
                                  extcoeff_opt,                         & 
                                  aot_sw, aot_lw,                       &
                                  omega_sw, gamma_sw,                   &
                                  jv_sca, jv_abs, jv_ga                 &
                                  )
  !
  !  Author:
  !  --------
  !  H. Tost, MPIC (adopted MADE aerosol radiation for use within EMAC),  2010
  !  Original MADE implementation in ECHAM5/MESSy by Axel Lauer, DLR, 2001-2003
  !
  !  Purpose:
  !  ---------
  !  Calculate aerosol optical properties (extinction cross section
  !  (sigma), single scattering albedo (omega), and assymetry factor
  !  (gamma)) for the current aerosol composition and size distribution
  !  for each mode from look-up tables.
  !
  !  Interface:
  !  ----------
  !  calc_aero_properties is called from aeropt_e5
  !
  !  Externals:
  !  ----------
  !  none
  !
  !--- Parameter list:
  !
  ! paerml(kproma,klev,0:naertot) = aerosol mass for each compound per mode [umol m-3]
  ! paernl(kproma,klev,nmod)      = aerosol number for each mode [N cm-3]
  ! paernl2(kproma,klev,nmod)     = aerosol number for each mode [N mol-1]
  ! prwet (kproma,klev,nmod)      = aerosol wet radius  [cm]
  ! znr   (kproma,klev,nmod)      = aerosol refractive index (Re)
  ! zni   (kproma,klev,nmod)      = aerosol refractive index (Im)
  ! zmsp  (kproma,klev,nmod)      = aerosol mie-size-parameter
  ! ppressi(kproma,klev+1)        = air pressure (interface) [Pa]
  ! aot_opt(kproma,klev,num_diag_elements,max_diag_wavelens) = aerosol optical thicknes, AOT (optional output)
  ! extcoeff_opt(kproma,klev,num_diag_elements,max_diag_wavelens) = aerosol extinction coefficient (optional output)
  ! aot_sw(kproma,klev,nsw)       = shortwave (sw) AOT (for ECHAM5 radiation scheme)
  ! aot_lw(kproma,klev,jpband)    = longwave  (lw) AOT (for ECHAM5 radiation scheme)
  ! omega_sw(kproma,klev,jpband)  = shortwave (sw) aerosol single scattering albedo (for ECHAM5 radiation scheme)
  ! gamma_sw(kproma,klev,jpband)  = shortwave (sw) aerosol assymetry factor         (for ECHAM5 radiation scheme)

    USE MESSY_MAIN_CONSTANTS_MEM,    ONLY: pi, g, M_air
    


    IMPLICIT NONE

    INTEGER :: kproma, klev

    REAL(dp), DIMENSION(kproma,klev,nmod)        :: paernl, paernl2, prwet
    REAL(dp), DIMENSION(kproma,klev+1)           :: ppressi

    ! aot optional wavelength for diagnostics
    REAL(dp), &
      DIMENSION(kproma,klev,num_diag_elements,max_diag_wavelens,nmod+1) :: &
      aot_opt
    ! extinction coeff. optional wavelength
    REAL(dp), &
      DIMENSION(kproma,klev,num_diag_elements,max_diag_wavelens,nmod+1) :: &
      extcoeff_opt
    ! shortwave
    REAL(dp), DIMENSION(kproma,klev,nsw)         :: aot_sw, omega_sw, gamma_sw 
    ! longwave
    REAL(dp), DIMENSION(kproma,klev,jpband)      :: aot_lw               
    ! jval
    REAL(dp), DIMENSION(kproma,klev,n_jv)        :: jv_sca, jv_abs, jv_ga
    REAL(dp), DIMENSION(kproma,klev,n_jv_calc)   :: zjv_sca, zjv_abs, zjv_ga
    
    !--- Local variables:

    INTEGER :: jm, jk, jl, jmod, jint, jt            ! loop indices 
                                                     !(jint = wavelength index)
    INTEGER :: n, idx1, n2                           ! help indices

    INTEGER :: nnr, nni, nmsp                        ! look-up table indices
    INTEGER :: nswlw                                 ! longwave/shortwave index
    ! index of ECHAM5 sw/lw band for current sub-interval
    INTEGER :: band                                  

    LOGICAL :: lshortwave

    REAL(dp), DIMENSION(kproma,klev,nmod)        :: zwetvol

    ! aerosol optical properties of each mode for ECHAM5 bands
    ! extinction cross section divided by wavelength squared [part-1]

    REAL(dp), DIMENSION(kproma,klev,nmod,nsw)    :: e5_sw_sigma   ! shortwave
    REAL(dp), DIMENSION(kproma,klev,nmod,jpband) :: e5_lw_sigma   ! longwave
    ! single scattering albedo [1]
    REAL(dp), DIMENSION(kproma,klev,nmod,nsw)    :: e5_sw_omega   ! shortwave
    ! assymetry factor [1]
    REAL(dp), DIMENSION(kproma,klev,nmod,nsw)    :: e5_sw_gamma   ! shortwave

    REAL(dp)                                     :: aottotal      ! total AOT
    REAL(dp)                                     :: extcoefftotal ! total extinction coeff.
    REAL(dp)                                     :: nr, ni, msp, zdpg, zfac


    REAL(dp), DIMENSION(kproma,klev,nmod)        :: znr   ! Re(refractive index)
    REAL(dp), DIMENSION(kproma,klev,nmod)        :: zni   ! Im(refractive index)
    REAL(dp), DIMENSION(kproma,klev,nmod)        :: zmsp  ! mie-size-parameter
    ! total volume [m3/kg(air)]
    REAL(dp), DIMENSION(kproma,klev,nmod)        :: zvoltot  
    ! component vol. [m3/kg(air)]
    REAL(dp), DIMENSION(kproma,klev,nmod,num_rad_types) :: zvol
    ! fraction of component vol. to total volume [0-1]
    REAL(dp), DIMENSION(kproma,klev,nmod,num_rad_types) :: zvolfrac 
    ! vert. integrated particle number conc. in _each_ level [#/m2]
    REAL(dp), DIMENSION(kproma,klev,nmod)        :: zsaer    
    REAL(dp), DIMENSION(kproma,klev,nmod)        :: zsigma, zomega, zgamma
    REAL(dp) :: zsig, zomg, zgam, dummy

    
    REAL(dp), PARAMETER :: ZEPS          = 1.e-20_dp
    
    REAL(dp) :: test_var
    ! 1. Initialization


    ZSIGMA(:,:,:) = 0._dp
    Zomega(:,:,:) = 0._dp
    Zgamma(:,:,:) = 0._dp

    ! 1.1 local arrays

    e5_sw_sigma(:,:,:,:) = 0._dp
    e5_lw_sigma(:,:,:,:) = 0._dp
    e5_sw_omega(:,:,:,:) = 0._dp
    e5_sw_gamma(:,:,:,:) = 0._dp

    zvolfrac(:,:,:,:)    = 0._dp
    ! 1.2 arrays of output stream
    
    aot_sw(:,:,:)   = 0._dp
    omega_sw(:,:,:) = 0._dp
    gamma_sw(:,:,:) = 0._dp
    aot_lw(:,:,:)   = 0._dp
    
    zvoltot(:,:,:)  = 0._dp
    zvol(:,:,:,:)   = 0._dp
    zsaer(:,:,:)    = 0._dp
    
    zjv_sca=0._dp
    zjv_abs=0._dp
    zjv_ga =0._dp

    ! factor for calculating the number of moles (air) per m2 from delta_p
    zfac = 1.0e3_dp / M_air / g
    DO jm=ns+1,nmod
      DO jk=1,klev
        DO jl=1,kproma
          
          ! 2.1 calculate vertically integrated number of particles 
          !     in grid cell [#/m2]
          zdpg = (ppressi(jl,jk+1) - ppressi(jl,jk)) * zfac
!          zsaer(jl,jk,jm) = paernl(jl,jk,jm) * 1.e6_dp * zdpg
          zsaer(jl,jk,jm) = paernl2(jl,jk,jm) * zdpg

          ! 2.2 Calculate volume of aerosol components and total volume.

          ! wet volume [cm3 m-3(air)]
          zwetvol(jl,jk,jm) = 4._dp * pi/3._dp * &
            (prwet(jl,jk,jm)/aerspec(jm)%ram2cmr)**3._dp
          zwetvol(jl,jk,jm) = zwetvol(jl,jk,jm) * paernl(jl,jk,jm) * 1.e6_dp
 
          ! *** calculate species bulk volume and fractions ***
          zvoltot(jl,jk,jm) = MAX(zwetvol(jl,jk,jm), 1.e-30_dp)


! new version start


        ENDDO
      ENDDO
    ENDDO

    zvoltot(:,:,:) = 0._dp
    DO jm = ns+1,nmod
      DO jt = 1,aerspec(jm)%count
        IDX1 = aerspec(jm)%rad_spec(jt)
        DO jk = 1,klev
          DO jl = 1,kproma
            zvol(jl,jk,jm,IDX1) = zvol(jl,jk,jm,IDX1) +                     &
                                  aerspec(jm)%paerml(jl,jk,jt) * 1.e-3_dp * &
                                  aerspec(jm)%molarmass(jt) /               &
                                  aerspec(jm)%density(jt)
          END DO
        END DO
      END DO

      DO jt=1,num_rad_types
         DO jk=1,klev
            DO jl=1,kproma
               zvoltot(jl,jk,jm)   = zvoltot(jl,jk,jm) + zvol(jl,jk,jm,jt)
            END DO
         ENDDO
      ENDDO

      DO jk = 1,klev
        DO jl = 1,kproma
          IF (zvoltot(jl,jk,jm) < 9.9e-29_dp) CYCLE
          DO jt=1,num_rad_types
            zvolfrac(jl,jk,jm,jt) = MIN( 1._dp, &
                                    zvol(jl,jk,jm,jt) / zvoltot(jl,jk,jm) )
          ENDDO
        ENDDO
      ENDDO
      test_var = 0._dp
      do jk=1,klev
        DO jl=1,kproma
          test_var = MAX(test_var, sum(zvolfrac(jl,jk,jm,:)))
          IF (paernl(jl,jk,jm) > 1.e-15_dp) THEN
            dummy = zvoltot(jl,jk,jm)/(paernl(jl,jk,jm)*1e6_dp)
            dummy = ( dummy * 0.75_dp / pi )**(1._dp/3._dp) * aerspec(jm)%ram2cmr
            prwet(jl,jk,jm) = dummy
          ENDIF
        enddo
      enddo
    ENDDO


!!! new version end


   ! 3. Get aerosol optical properties (extinction cross section, single
   !    scattering albedo, assymetry factor) for each mode and each ECHAM5
   !    shortwave band and the extinction cross section for each ECHAM5
   !    longwave band from look-up tables.
   !    Input for look-up tables: Re(refractive index), Im(refractive index),
   !    mie-size-parameter (=2*pi*r/lambda)

    lshortwave = .true.

    DO jint = 1, num_sw_intervals+num_opt_wavelens+n_jv_bands+num_lw_intervals

      ! 3.1 Set shortwave/longwave switch according to current sub-interval.
      
      IF (lambda(jint) > lambda(num_sw_intervals) &
           .OR. (jint > num_sw_intervals+num_opt_wavelens+n_jv_bands) ) THEN
         lshortwave = .FALSE.
      ELSE
         lshortwave = .TRUE.
      ENDIF

      ! 3.2 Calculate volume averaged real (nr) and imaginary (ni) part of
      !     the refractive index for an internal mixture (Ouimetter and Flagan,
      !     Atmos. Environm., 16, 1982). This method follows Stier et al., 
      !     Atmos. Chem. Phys., 2005.

      znr(:,:,:)     = 0._dp
      zni(:,:,:)     = 0._dp

      DO jm = ns+1, nmod
      ! 3.2 Calculate mie-size-parameter (msp).
        DO jk = 1, klev
          DO jl = 1, kproma
            zmsp(jl,jk,jm) = 2.0_dp * pi * prwet(jl,jk,jm)*1.e-2_dp / lambda(jint)
          END DO ! jl-loop
        END DO ! jk-loop (levels)

      ! 3.3a Calculate volume averaged real (nr) and imaginary (ni) part of
      !      the refractive index for an internal mixture (Ouimetter and Flagan,
      !      Atmos. Environm., 16, 1982). This method follows Stier et al., 
      !      Atmos. Chem. Phys., 2005.

      ! special treatmen if hydrophobic modes are externally mixed    
        IF (L_EXTMIXT) THEN
          IF (jm > CS) CYCLE
        ENDIF

        DO jk = 1, klev
          DO jl = 1, kproma

            ! "WASO" ---> SO4 + NH4 + NO3

            znr(jl,jk,jm)     = znr(jl,jk,jm) + ref_re(ri_waso,jint) &
                              * zvolfrac(jl,jk,jm,ri_waso)
            zni(jl,jk,jm)     = zni(jl,jk,jm) + ref_im(ri_waso,jint) &
                              * zvolfrac(jl,jk,jm,ri_waso)

            ! BC

            znr(jl,jk,jm)     = znr(jl,jk,jm) + ref_re(ri_bc,jint)   &
                              * zvolfrac(jl,jk,jm,ri_bc)
            zni(jl,jk,jm)     = zni(jl,jk,jm) + ref_im(ri_bc,jint)   &
                              * zvolfrac(jl,jk,jm,ri_bc)

            ! POM

            znr(jl,jk,jm)     = znr(jl,jk,jm) + ref_re(ri_oc,jint)   &
                              * zvolfrac(jl,jk,jm,ri_oc)
            zni(jl,jk,jm)     = zni(jl,jk,jm) + ref_im(ri_oc,jint)   &
                              * zvolfrac(jl,jk,jm,ri_oc)

            ! mineral dust (DU)

            znr(jl,jk,jm)     = znr(jl,jk,jm) + ref_re(ri_du,jint)   &
                              * zvolfrac(jl,jk,jm,ri_du)
            zni(jl,jk,jm)     = zni(jl,jk,jm) + ref_im(ri_du,jint)   &
                              * zvolfrac(jl,jk,jm,ri_du)

            ! sea salt (SS)

            znr(jl,jk,jm)     = znr(jl,jk,jm) + ref_re(ri_ss,jint)   &
                              * zvolfrac(jl,jk,jm,ri_ss)
            zni(jl,jk,jm)     = zni(jl,jk,jm) + ref_im(ri_ss,jint)   &
                              * zvolfrac(jl,jk,jm,ri_ss)

            ! aerosol liquid water (H2O)

            znr(jl,jk,jm)     = znr(jl,jk,jm) + ref_re(ri_h2o,jint)  &
                              * zvolfrac(jl,jk,jm,ri_h2o)
            zni(jl,jk,jm)     = zni(jl,jk,jm) + ref_im(ri_h2o,jint)  &
                              * zvolfrac(jl,jk,jm,ri_h2o)

          END DO ! jl-loop
        END DO ! jk-loop (levels)

      END DO ! jm-loop (aerosol modes)

  ! 3.4a Get aerosol optical properties for each mode from look-up tables
  !      for ECHAM5 shortwave and longwave bands.
  !      Multiply sigma by wavelength squared to get exinction cross
  !      section [m2/particle].

      IF (lshortwave) THEN
        nswlw  = 1   ! shortwave
      ELSE
        nswlw  = 2   ! longwave
      END IF

      DO jm = ns + 1, nmod

        jmod = radmod(jm)

        ! special treatmen if hydrophobic modes are externally mixed    
        IF (L_EXTMIXT) THEN
          IF (jm > CS) CYCLE
        ENDIF

        DO jk = 1, klev
          DO jl = 1, kproma
              
            ! 1. Limit input values to valid range.

            nr  = MAX(nr_min(nswlw,jmod), &
                  MIN(znr(jl,jk,jm), nr_max(nswlw,jmod)))
            ni  = MAX(ni_min(nswlw,jmod), &
                  MIN(zni(jl,jk,jm), ni_max(nswlw,jmod)))
            msp = MAX(msp_min(nswlw,jmod),&
                  MIN(zmsp(jl,jk,jm), msp_max(nswlw,jmod)))
            ! 2. Get look-up table indices.

            !    1st dimension "nr"  = linear
            nnr  = 1 + NINT((nr - nr_min(nswlw,jmod)) / nr_step(nswlw,jmod))
            !    2nd dimension "ni"  = logarithmic
            nni  = 1 + NINT((LOG(ni)  - log_ni_min(nswlw,jmod)) &
                 / ni_step(nswlw,jmod))
            !    3rd dimension "msp" = logarithmic
            nmsp = 1 + NINT((LOG(msp) - log_msp_min(nswlw,jmod)) &
                 / msp_step(nswlw,jmod))

            ! 3. Do table look-up.
            IF (nswlw == 1) THEN ! shortwave
!do only for aitken and larger particles
              if(prwet(jl,jk,jm) > 6.e-7_dp) then
                zsigma(jl,jk,jm) = lut_sw_sigma(nnr,nni,nmsp,jmod) &
                                 * lambda_squared(jint)
                zomega(jl,jk,jm) = lut_sw_omega(nnr,nni,nmsp,jmod)
                zgamma(jl,jk,jm) = lut_sw_gamma(nnr,nni,nmsp,jmod)
              else
                zsigma(jl,jk,jm) = zeps
                zomega(jl,jk,jm) = MIN_VAL_OMEGA
                zgamma(jl,jk,jm) = zeps
              endif
            ELSE                ! longwave
              if(prwet(jl,jk,jm) > 6.e-7_dp) then
                zsigma(jl,jk,jm) = lut_lw_sigma(nnr,nni,nmsp,jmod) &
                                 * lambda_squared(jint)
              ELSE
                zsigma(jl,jk,jm) = zeps
              END if
            END IF
            
          END  DO ! jl-loop
        END DO ! jk-loop (levels)

      END DO ! jmod-loop (aerosol modes)

      ! special treatment for externally mixed modes
      ! weight the resulting extinction and not the refractive indices
      ! with the volume of the individual compounds
      IF (L_EXTMIXT) THEN
        IF (lshortwave) THEN
          nswlw  = 1   ! shortwave
        ELSE
          nswlw  = 2   ! longwave
        END IF

        DO jm = CS+1, nmod
          zsigma(:,:,jm) = 0._dp
          zomega(:,:,jm) = 0._dp
          zgamma(:,:,jm) = 0._dp
          ! map KI -> KS, AI -> AS, CI -> CS
          jmod = jm - CS

          DO jt=1,num_rad_types
            ! skip hydrophilic material since it will not occur 
            ! in the hydrophobic modes
            IF ( (jt == ri_ss)  .OR. (jt == ri_h2o)  .OR. &
                 (jt == ri_so4) .OR. (jt == ri_waso) ) CYCLE

            DO jl = 1,kproma
              DO jk=1,klev
                zsig = 0._dp
                zomg = 0._dp
                zgam = 0._dp
      ! 3.3b Determine the refractive indices per species
                znr(jl,jk,jm) = ref_re(jt,jint)
                zni(jl,jk,jm) = ref_im(jt,jint)

      ! 3.4b Get aerosol optical properties for each mode from look-up tables
      !      for ECHAM5 shortwave and longwave bands.
      !      Multiply sigma by wavelength squared to get exinction cross
      !      section [m2/particle].
                ! 1. Limit input values to valid range.
                nr  = MAX(nr_min(nswlw,jmod), &
                      MIN(znr(jl,jk,jm), nr_max(nswlw,jmod)))
                ni  = MAX(ni_min(nswlw,jmod), &
                      MIN(zni(jl,jk,jm), ni_max(nswlw,jmod)))
                msp = MAX(msp_min(nswlw,jmod),&
                      MIN(zmsp(jl,jk,jm), msp_max(nswlw,jmod)))
                ! 2. Get look-up table indices.

                !    1st dimension "nr"  = linear
                nnr  = 1 + NINT((nr - nr_min(nswlw,jmod)) / nr_step(nswlw,jmod))
                !    2nd dimension "ni"  = logarithmic
                nni  = 1 + NINT((LOG(ni)  - log_ni_min(nswlw,jmod)) &
                     / ni_step(nswlw,jmod))
                !    3rd dimension "msp" = logarithmic
                nmsp = 1 + NINT((LOG(msp) - log_msp_min(nswlw,jmod)) &
                     / msp_step(nswlw,jmod))

                ! 3. Do table look-up.
                IF (nswlw == 1) THEN ! shortwave
                  !do only for aitken and larger particles
                  if(prwet(jl,jk,jm) > 6.e-7_dp) then
                    zsig = lut_sw_sigma(nnr,nni,nmsp,jmod) &
                         * lambda_squared(jint)
                    zomg = lut_sw_omega(nnr,nni,nmsp,jmod)
                    zgam = lut_sw_gamma(nnr,nni,nmsp,jmod)
                  endif
                ELSE                ! longwave
                  if(prwet(jl,jk,jm) > 6.e-7_dp) then
                    zsig = lut_lw_sigma(nnr,nni,nmsp,jmod) &
                         * lambda_squared(jint)
                  END if
                END IF
                zsigma(jl,jk,jm) = zsigma(jl,jk,jm) + &
                                   zsig * zvolfrac(jl,jk,jm,jt)
                zomega(jl,jk,jm) = zomega(jl,jk,jm) + &
                                   zomg * zvolfrac(jl,jk,jm,jt)
                zgamma(jl,jk,jm) = zgamma(jl,jk,jm) + &
                                   zgam * zvolfrac(jl,jk,jm,jt)
              END DO
            END DO
          END DO
        END DO
      END IF
      

      ! 4. Map sub-intervals to ECHAM5 bands (I),
      !    calculate optional diagnostics.


      band = int2band(jint) ! index of ECHAM5 sw/lw band for current
                            ! sub-interval

      IF (band > 0) THEN  ! wavelength relevant to ECHAM5 radiation scheme
        IF (lshortwave) THEN  ! sub-interval of a ECHAM5 shortwave band
       
          DO jmod = ns+1, nmod
            DO jk = 1, klev
              DO jl = 1, kproma
                e5_sw_sigma(jl,jk,jmod,band) = e5_sw_sigma(jl,jk,jmod,band) &
                                             + zsigma(jl,jk,jmod)           &
                                             * weight(jint)
                e5_sw_omega(jl,jk,jmod,band) = e5_sw_omega(jl,jk,jmod,band) &
                                             + zsigma(jl,jk,jmod)           &
                                             * zomega(jl,jk,jmod)           &
                                             * weight(jint)
                e5_sw_gamma(jl,jk,jmod,band) = e5_sw_gamma(jl,jk,jmod,band) &
                                             + zsigma(jl,jk,jmod)           &
                                             * zomega(jl,jk,jmod)           &
                                             * zgamma(jl,jk,jmod)           &
                                             * weight(jint)
              END DO ! jl-loop
            END DO ! jk-loop (levels)
          END DO ! jmod-loop (aerosol modes)

        ELSE                  ! sub-interval of a ECHAM5 longwave band
          DO jmod = ns+1, nmod
            DO jk = 1, klev
              DO jl = 1, kproma
                e5_lw_sigma(jl,jk,jmod,band) = e5_lw_sigma(jl,jk,jmod,band) &
                                             + zsigma(jl,jk,jmod)           &
                                             * weight(jint)
              END DO ! jl-loop
            END DO ! jk-loop (levels)
          END DO ! jmod-loop (aerosol modes)

        END IF

      ELSE ! 4.2 optional wavelength, for diagnostics only

        ! calculate vertically integrated total/component
        ! aerosol optical thickness (AOT)

        n = jint - num_sw_intervals

        IF(n < 1) CYCLE

        IF (n <= num_opt_wavelens) THEN

          aot_opt(:,:,:,n,:) = 0._dp
          extcoeff_opt(:,:,:,n,:) = 0._dp

          DO jmod = ns+1, nmod
            DO jk = 1, klev
              DO jl = 1, kproma
                aottotal = zsigma(jl,jk,jmod) * zsaer(jl,jk,jmod)

                aot_opt(jl,jk,diag_tot,n,jmod) = aot_opt(jl,jk,diag_tot,n,jmod) &
                  + aottotal
                aot_opt(jl,jk,diag_tot,n,nmod+1) = &
                  aot_opt(jl,jk,diag_tot,n,nmod+1) + aottotal
                
                aot_opt(jl,jk,diag_bc,n,jmod)   = aot_opt(jl,jk,diag_bc,n,jmod)  &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_bc)
                aot_opt(jl,jk,diag_bc,n,nmod+1)  = &
                  aot_opt(jl,jk,diag_bc,n,nmod+1)  &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_bc)
                
                aot_opt(jl,jk,diag_oc,n,jmod)   = aot_opt(jl,jk,diag_oc,n,jmod)  &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_oc)
                aot_opt(jl,jk,diag_oc,n,nmod+1)  = &
                  aot_opt(jl,jk,diag_oc,n,nmod+1)  &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_oc)
                
                aot_opt(jl,jk,diag_du,n,jmod)   = aot_opt(jl,jk,diag_du,n,jmod)  &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_du)
                aot_opt(jl,jk,diag_du,n,nmod+1)  = &
                  aot_opt(jl,jk,diag_du,n,nmod+1)  &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_du)
                
                aot_opt(jl,jk,diag_ss,n,jmod)   = aot_opt(jl,jk,diag_ss,n,jmod)  &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_ss)
                aot_opt(jl,jk,diag_ss,n,nmod+1)  = &
                  aot_opt(jl,jk,diag_ss,n,nmod+1)  &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_ss)
                
                aot_opt(jl,jk,diag_h2o,n,jmod)  = aot_opt(jl,jk,diag_h2o,n,jmod) &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_h2o)
                aot_opt(jl,jk,diag_h2o,n,nmod+1) = &
                  aot_opt(jl,jk,diag_h2o,n,nmod+1) &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_h2o)
                
                aot_opt(jl,jk,diag_waso,n,jmod) = aot_opt(jl,jk,diag_waso,n,jmod)&
                  + aottotal * zvolfrac(jl,jk,jmod,ri_waso)
                aot_opt(jl,jk,diag_waso,n,nmod+1) = &
                  aot_opt(jl,jk,diag_waso,n,nmod+1) &
                  + aottotal * zvolfrac(jl,jk,jmod,ri_waso)

               extcoefftotal = zsigma(jl,jk,jmod) * paernl(jl,jk,jmod) * 1.e6_dp

                extcoeff_opt(jl,jk,diag_tot,n,jmod) = extcoeff_opt(jl,jk,diag_tot,n,jmod) &
                  + extcoefftotal
                extcoeff_opt(jl,jk,diag_tot,n,nmod+1) = &
                  extcoeff_opt(jl,jk,diag_tot,n,nmod+1) + extcoefftotal

                extcoeff_opt(jl,jk,diag_bc,n,jmod) = extcoeff_opt(jl,jk,diag_bc,n,jmod) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_bc)
                extcoeff_opt(jl,jk,diag_bc,n,nmod+1) = &
                 extcoeff_opt(jl,jk,diag_bc,n,nmod+1) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_bc)

               extcoeff_opt(jl,jk,diag_oc,n,jmod) = extcoeff_opt(jl,jk,diag_oc,n,jmod) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_oc)
                extcoeff_opt(jl,jk,diag_oc,n,nmod+1) = &
                  extcoeff_opt(jl,jk,diag_oc,n,nmod+1) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_oc)

                extcoeff_opt(jl,jk,diag_du,n,jmod) = extcoeff_opt(jl,jk,diag_du,n,jmod) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_du)
                extcoeff_opt(jl,jk,diag_du,n,nmod+1) = &
                  extcoeff_opt(jl,jk,diag_du,n,nmod+1) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_du)

                extcoeff_opt(jl,jk,diag_ss,n,jmod) = extcoeff_opt(jl,jk,diag_ss,n,jmod) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_ss)
                extcoeff_opt(jl,jk,diag_ss,n,nmod+1) = &
                  extcoeff_opt(jl,jk,diag_ss,n,nmod+1) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_ss)

                extcoeff_opt(jl,jk,diag_h2o,n,jmod) = extcoeff_opt(jl,jk,diag_h2o,n,jmod) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_h2o)
                extcoeff_opt(jl,jk,diag_h2o,n,nmod+1) = &
                  extcoeff_opt(jl,jk,diag_h2o,n,nmod+1) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_h2o)

                extcoeff_opt(jl,jk,diag_waso,n,jmod) = extcoeff_opt(jl,jk,diag_waso,n,jmod) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_waso)
                extcoeff_opt(jl,jk,diag_waso,n,nmod+1) = &
                  extcoeff_opt(jl,jk,diag_waso,n,nmod+1) &
                  + extcoefftotal * zvolfrac(jl,jk,jmod,ri_waso)
              END DO ! jl-loop
            END DO ! jk-loop (levels)
          END DO ! jmod-loop (aerosol modes)

        END IF
        ! for JVAL
        n2 = jint - num_sw_intervals - num_opt_wavelens
        IF ((n2 > 0) .AND. (n2 <= n_jv_bands)) THEN
          DO jmod = ns+1, nmod
            DO jk = 1, klev
              DO jl = 1, kproma
                zjv_sca(jl,jk,n2) = zjv_sca(jl,jk,n2) + zsaer(jl,jk,jmod) * &
                  zsigma(jl,jk,jmod)
                zjv_abs(jl,jk,n2) = zjv_abs(jl,jk,n2) + zsaer(jl,jk,jmod) * &
                  zsigma(jl,jk,jmod) * zomega(jl,jk,jmod)
                zjv_ga(jl,jk,n2)  = zjv_ga(jl,jk,n2)  + zsaer(jl,jk,jmod) * &
                  zsigma(jl,jk,jmod) * zomega(jl,jk,jmod) * zgamma(jl,jk,jmod)
              END DO
            END DO
          END DO
          DO jk=1,klev
            DO jl=1,kproma
              IF (zjv_abs(jl,jk,n2) > ZEPS) THEN
                zjv_ga(jl,jk,n2)  = zjv_ga(jl,jk,n2)  / zjv_abs(jl,jk,n2)
                zjv_abs(jl,jk,n2) = zjv_abs(jl,jk,n2) / zjv_sca(jl,jk,n2)
              ELSE
                zjv_sca(jl,jk,n2) = MIN_VAL
                zjv_abs(jl,jk,n2) = MIN_VAL_OMEGA
                zjv_ga(jl,jk,n2)  = MIN_VAL
              END IF
            END DO
          END DO
          
          DO jk=1,klev
            DO jl=1,kproma
              jv_sca(jl,jk,n2) = zjv_sca(jl,jk,n2) * zjv_abs(jl,jk,n2)
              jv_abs(jl,jk,n2) = zjv_sca(jl,jk,n2) * (1._dp -zjv_abs(jl,jk,n2))
              jv_ga(jl,jk,n2)  = zjv_ga(jl,jk,n2)
            END DO
          END DO
        END IF
      END IF
    END DO ! jint-loop (all sub-intervals of ECHAM5 sw+lw bands)

    
    ! 5. Map sub-intervals to ECHAM5 bands (II)

    ! 5.1 shortwave bands

    DO band = 1, nsw
      DO jmod = ns+1, nmod
        DO jk = 1, klev
          DO jl = 1, kproma
            IF (e5_sw_omega(jl,jk,jmod,band) > ZEPS) THEN
              e5_sw_gamma(jl,jk,jmod,band) = e5_sw_gamma(jl,jk,jmod,band) &
                                           / e5_sw_omega(jl,jk,jmod,band)
              e5_sw_omega(jl,jk,jmod,band) = e5_sw_omega(jl,jk,jmod,band) &
                                           / e5_sw_sigma(jl,jk,jmod,band)
              e5_sw_sigma(jl,jk,jmod,band) = e5_sw_sigma(jl,jk,jmod,band) &
                                           /sumweight_sw(band)
            END IF
          END DO ! jl-loop
        END DO ! jk-loop (levels)
      END DO ! jmod-loop (aerosol modes)
    END DO ! band-loop (ECHAM sw bands)

   ! 5.2 longwave bands
    DO band = 1, jpband
      DO jmod = ns+1, nmod
        DO jk = 1, klev
          DO jl = 1, kproma
            e5_lw_sigma(jl,jk,jmod,band) = e5_lw_sigma(jl,jk,jmod,band) &
                                         / sumweight_lw(band)
          END DO ! jl-loop
        END DO ! jk-loop (levels)
      END DO ! jmod-loop (aerosol modes)
   END DO ! band-loop (ECHAM lw bands)
   
   ! 6. Calculate aerosol optical properties for ECHAM5 radiation scheme

   !    arrays for output stream (aot_sw, omega_sw, gamma_sw, aot_lw)
   !    have already been set to 0._dp.

   DO band = 1, nsw
     DO jk = 1, klev
       DO jmod = ns+1, nmod
         DO jl = 1, kproma
           aot_sw(jl,jk,band)   = aot_sw(jl,jk,band)         &
                                + e5_sw_sigma(jl,jk,jmod,band) &
                                * zsaer(jl,jk,jmod)

           omega_sw(jl,jk,band) = omega_sw(jl,jk,band)       &
                                + e5_sw_sigma(jl,jk,jmod,band) &
                                * e5_sw_omega(jl,jk,jmod,band) &
                                * zsaer(jl,jk,jmod)

           gamma_sw(jl,jk,band) = gamma_sw(jl,jk,band)       &
                                + e5_sw_sigma(jl,jk,jmod,band) &
                                * e5_sw_omega(jl,jk,jmod,band) &
                                * e5_sw_gamma(jl,jk,jmod,band) &
                                * zsaer(jl,jk,jmod)

         END DO ! jl-loop
       END DO ! jmod-loop (aerosol modes)

       DO jl = 1, kproma
         IF(omega_sw(jl,jk,band) > ZEPS) THEN
           gamma_sw(jl,jk,band) = gamma_sw(jl,jk,band) / omega_sw(jl,jk,band)
           omega_sw(jl,jk,band) = omega_sw(jl,jk,band) / aot_sw  (jl,jk,band)
         ELSE
           aot_sw  (jl,jk,band) = MIN_VAL
           omega_sw(jl,jk,band) = MIN_VAL_OMEGA
           gamma_sw(jl,jk,band) = MIN_VAL
         END IF
       END DO ! jl-loop
     END DO ! jk-loop (levels)
   END DO ! band-loop (ECHAM sw bands)



   ! 6.2 ECHAM5 longwave bands
   DO band = 1, jpband
     DO jmod = ns+1, nmod
       DO jk = 1, klev
         DO jl = 1, kproma
           aot_lw(jl,jk,band)   = aot_lw(jl,jk,band)        &
                                + e5_lw_sigma(jl,jk,jmod,band) &
                                    * zsaer(jl,jk,jmod)
         END DO ! jl-loop
       END DO ! jk-loop (levels)
     END DO ! jmod-loop (aerosol modes)
   END DO ! band-loop (ECHAM lw bands)

   !--------------------------------------------------------------------

 END SUBROUTINE calc_aero_properties

!===============================================================================

 ! =======================================================================
 FUNCTION vrev2(p)

   IMPLICIT NONE
   INTRINSIC :: SIZE

   ! inverts array p(kproma,klev) in the second index
   
   REAL(dp), INTENT(in), DIMENSION(:,:)     :: p
   REAL(dp), DIMENSION(SIZE(p,1),SIZE(p,2)) :: vrev2
   
   INTEGER :: klev, jk

   klev=SIZE(p,2)
   
   DO jk=1,klev
      vrev2(:,jk)=p(:,klev+1-jk)
   END DO
   
 END FUNCTION vrev2
 ! =======================================================================

END MODULE MESSY_AEROPT
