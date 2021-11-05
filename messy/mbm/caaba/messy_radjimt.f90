!*****************************************************************************
MODULE messy_radjimt

  USE messy_main_constants_mem, ONLY: DP, k_B, HLINE2, STRLEN_ULONG, R_gas, pi
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_3D_ARRAY &
                                    , Spline1d, Splint1d

  USE messy_cmn_photol_mem      ! IP_MAX, ip_*, jname

  IMPLICIT NONE
  PRIVATE
  ! NAME OF SUBMODEL
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'radjimt'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'

  PUBLIC :: radjimt_read_nml_ctrl           ! read CTRL namelist and initialize
  PUBLIC :: MesoSol
  PUBLIC :: ThermoSol
  PUBLIC :: CO2ChapManHeating 
  PUBLIC :: radjimt_readfluxes
  PUBLIC :: Branchingratios

  ! Number of heating and cooling rates
  INTEGER, PARAMETER,PUBLIC :: heatcool_rates  = 9

  INTEGER, PARAMETER :: total_heat_l    = 1
  INTEGER, PARAMETER :: o3_chap_heat_l  = 2
  INTEGER, PARAMETER :: o3_hart_heat_l  = 3
  INTEGER, PARAMETER :: o3_hugg_heat_l  = 4 !
  INTEGER, PARAMETER :: o3_herz_heat_l  = 5
  INTEGER, PARAMETER :: o2_meso_heat_l  = 6 !
  INTEGER, PARAMETER :: euv_heat_l      = 7 !
  INTEGER, PARAMETER :: uv_heat_l       = 8 !

  INTEGER, PARAMETER :: total_cool_l    = 1
  INTEGER, PARAMETER :: co2_cool_l      = 2
  INTEGER, PARAMETER :: o3_cool_l       = 3
  INTEGER, PARAMETER :: o_cool_l        = 4 
  INTEGER, PARAMETER :: no_cool_l       = 5

  INTEGER, PARAMETER, PUBLIC :: ht_branches = 17                 
  REAL(dp), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: branch

  INTEGER :: b_aurqo2_b1   = 1    ! O2 + e* -> O2+ + 2*e- 
  INTEGER :: b_aurqo2_b2   = 2    ! O2 + e* -> O+ + O + e- 
  INTEGER :: b_aurqn2_b1   = 3    ! N2 + e* -> N2+ + 2*e- 
  INTEGER :: b_aurqn2_b2   = 4    ! N2 + e* -> N+ + N{4s} + 2*e- 
  INTEGER :: b_aurqn2_b3   = 5    ! N2 + e* -> N+ + N{2d} + 2*e- 
  INTEGER :: b_aurqn2_b4   = 6    ! N2 + e* -> 0.47*N{4s} + 0.53*N{2d} + e- 
  INTEGER :: b_jio2_b1     = 7    ! O2 + hv(<105nm) -> O2+ + e* 
  INTEGER :: b_jio2_b2     = 8    ! O2 + hv(<105nm) -> O+ + O + e* 
  INTEGER :: b_jin2_b1     = 9    ! N2 + hv(<105nm) -> N2+ + e* 
  INTEGER :: b_jin2_b2     = 10   ! N2 + hv(<105nm) -> N+ + N{4s} + e* 
  INTEGER :: b_jin2_b3     = 11   ! N2 + hv(<105nm) -> N+ + N{2d} + e* 
  INTEGER :: b_sen2_b1     = 12   ! N2 + e* -> N2+ + 2*e- 
  INTEGER :: b_sen2_b2     = 13   ! N2 + e* -> N+ + N{4s} + 2e-* 
  INTEGER :: b_sen2_b3     = 14   ! N2 + e* -> N+ + N{2d} + 2e-* 
  INTEGER :: b_sen2_b4     = 15   ! N2 + e* -> 0.47*N{4s} + 0.53*N{2d} + e-* 
  INTEGER :: b_seo2_b1     = 16   ! O2 + e* -> O2+ + 2*e-       
  INTEGER :: b_seo2_b2     = 17   ! O2 + e* -> O+ + O + e- 

  ! Flag to use exothermic chemical heating and associated
  ! efficiencies, or not (see MesoThermoheat for details)
  INTEGER, PUBLIC :: chem_heat_on = 0      &
                         , chem_heat_on_last = 0

  ! Bounds variables for efficiency interpolation
  INTEGER   :: maxn=-1, minn=1000

  REAL(dp), POINTER, PUBLIC :: jx(:,:,:,:)

  ! Mesospheric heating efficiencies  
  REAL(dp), ALLOCATABLE, DIMENSION(:), PUBLIC :: O2eff_bulk, O3eff_bulk &
       , htlyeff, euveff, src1eff, PEOP, PEN2, PEO2

  ! NWAVES_h2o2 = no. of wavelengths in range 260 - 350nm
  INTEGER, PARAMETER :: NWAVES_FUV = 92

  ! s2k flux and wavelengths
  REAL(dp) :: S2K_h2o2_ALLFLUXES(NWAVES_fuv,2)
  REAL(dp) :: S2K_FLUX_FUV(NWAVES_FUV), WAVELS_FUV(NWAVES_FUV)

  ! Cross sections of O3(main absorbing species)
  REAL(dp) :: CSAO3(NWAVES_FUV)

  REAL(dp), ALLOCATABLE, DIMENSION(:,:), PUBLIC :: CSH2O2

  REAL(dp), PARAMETER :: GSCON = R_gas*1.e3_dp

  ! GLOBAL CTRL-NAMELIST
  CHARACTER(LEN=STRLEN_ULONG), PUBLIC :: file_sk2_h2o2_fluxes = ''  &
                                       , file_O3_Xsections  = ''    &
                                       , file_s2k_fluxes    = ''    &
                                       , file_all_Xsections = ''
CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE radjimt_read_nml_ctrl(status, iou)

    ! READ RADJIMT NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    ! LOCAL
    LOGICAL :: lex   ! file exists?
    INTEGER :: fstat ! file status
    CHARACTER(LEN=*), PARAMETER :: substr = 'radjimt_read_nml_ctrl'

    NAMELIST /CTRL/ file_sk2_h2o2_fluxes, file_O3_Xsections &
                  , file_s2k_fluxes, file_all_Xsections


    ! INITIALIZE
    status = 1 ! error

    ! INPUT NAMELIST
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! no error

  END SUBROUTINE radjimt_read_nml_ctrl

  ! --------------------------------------------------------------------------

  SUBROUTINE radjimt_readfluxes(s2k_h2o2_fluxes_f,O3_Xsections_f)

    IMPLICIT NONE

    ! file handles
    INTEGER, INTENT(IN) :: s2k_h2o2_fluxes_f, O3_Xsections_f

    ! LOCAL
    ! Cross sections of O3(main absorbing species)
    REAL(dp) :: CSAO3_all(NWAVES_FUV,2)
    INTEGER    :: countw, countc

    ! Initialise solar fluxes and cross section
    ! Read the solar flux file to get wavelengths and fluxes,
    ! and the O3 absorption cross sections file
    ! at wavelength bins of 1 Angstrom.

    OPEN (s2k_h2o2_fluxes_f,                                            &
         FILE=TRIM(file_sk2_h2o2_fluxes),STATUS='old')

    OPEN (O3_Xsections_f,                                               &
         FILE=TRIM(file_O3_Xsections),STATUS='old')

    WAVELOOP: DO countw=1,NWAVES_fuv

          READ(s2k_h2o2_fluxes_f,*)                                     &
               (S2K_h2o2_ALLFLUXES(countw,countc),countc=1,2)

          READ(O3_Xsections_f,*)                                        &
               (CSAO3_all(countw,countc),countc=1,2)

          ! convert wavelengths from m to nm
          WAVELS_fuv(countw) = S2K_h2o2_ALLFLUXES(countw,1)             &
                               *1.e9
          S2K_FLUX_fuv(countw) = S2K_h2o2_ALLFLUXES(countw,2)

          CSAO3(countw) = CSAO3_all(countw,2)

    ENDDO WAVELOOP

    CLOSE(s2k_h2o2_fluxes_f)
    CLOSE(O3_Xsections_f)


  END SUBROUTINE radjimt_readfluxes


  !========================================================
  !=                        MesoSol                       =
  !========================================================
  
  SUBROUTINE MesoSol(m,l,ht_dim, heatrates,                             &
       gas1_1d, gas2_1d, gas3_1d, gas_O3_1d, gas_H2O_1d, gas_H2O2_1d, &
       nden1d, csza,                               &
       r0_eff, scht1d, temp1d, ht1d, cp1d, den1d, eccentric,                  &
       f107, pres,  mmw_gas1, mmw_gas2, GRAV0,l_reduceO3heating)

    INTEGER,  INTENT(IN) :: m, l, ht_dim
    REAL(dp), INTENT(IN) :: csza, f107, eccentric, mmw_gas1, mmw_gas2
    REAL(dp), INTENT(IN) :: ht1d(ht_dim), temp1d(ht_dim)
    REAL(dp), INTENT(IN) :: gas1_1d(ht_dim), gas2_1d(ht_dim), gas3_1d(ht_dim)
    REAL(dp), INTENT(IN) :: gas_O3_1d(ht_dim), gas_H2O_1d(ht_dim), gas_H2O2_1d(ht_dim)
    REAL(dp), INTENT(IN) :: nden1d(ht_dim), den1d(ht_dim), pres(ht_dim)
    REAL(dp), INTENT(IN) :: cp1d(ht_dim), R0_eff(ht_dim), scht1d(ht_dim)
    REAL(dp), POINTER :: heatrates(:,:,:,:)
    REAL(dp), INTENT(IN) :: GRAV0
    LOGICAL,  INTENT(IN) :: l_reduceO3heating
    REAL(dp) :: colo2(ht_dim), o2cm3, o3cm3, colo3(ht_dim)
    REAL(dp) :: colh2o(ht_dim), h2ocm3
    REAL(dp) :: Z, csza_ch 

    REAL(dp) :: scht1d_o3(ht_dim), csza_o3
    REAL(dp) :: scht1d_o2(ht_dim), csza_o2
    REAL(dp) :: scht1d_h2o(ht_dim), csza_h2o
    
    REAL(dp) :: srband, hartley, chappuis, huggins, meso_heat
    REAL(dp) :: ozheat, o2heat
    REAL(dp) :: srband_fac, hartley_fac, huggins_fac, chappuis_fac
    REAL(dp) :: herzberg
    
    REAL(dp) :: O3eff63(63), O2eff63(63), pres63(63), splft1(63)
    REAL(dp) :: splft2(63)
    REAL(dp) :: eff1, eff2, pval, ex

    REAL(dp) :: fmxfmn

    INTEGER    :: initialise = 1

    INTEGER    :: n

    INTEGER :: status = 1


    ! Hartley band absorption cross sections (*1e-19) and 200-300nm fluxes

    REAL(dp) :: htly(11)=(/1.5,5.0,18.,46.,82.,110., 110.,75.,40.,13.,2.66/)
    REAL(dp) :: fluxh(11)=(/6.7E11,2.7E12,5.1E12,6.1E12,5.3E12,7.6E12,       &
                  1.3E13,4.1E13,1.3E13,6.2E13,5.9E13/)

    REAL(dp) :: section, dfluxh(11), tmp
    INTEGER    :: nn

    ! Level for crude stratospheric heatin
    INTEGER    :: fade = 0

    ! Mylynczak and Solomon bulk O3 heating efficienies
    DATA O3EFF63/1.000,1.000,1.000,1.000,1.000,                               &
         1.000,1.000,.999,.992,.980,.960,.946,                                &
         .914,.870,.846,.814,.783,.753,.726,                                  &
         .704,.690,.685,.687,.690,.691,.691,                                  &
         .690,.690,.690,.690,.690,.690,.690,                                  &
         .690,.691,.691,.691,.691,.691,.691,                                  &
         .691,.691,.691,.691,.691,.691,.691,                                  &
         .691,.691,.691,.691,.691,.691,.691,                                  &
         .691,.691,.691,.691,.691,.691,.691,                                  &
         .691,.691/
    
    ! Mylynczak and Solomon bulk O2 heating efficienies
    DATA O2EFF63/1.000,1.000,1.000,1.000,1.000,                               &
         1.000,1.000,1.000,1.000,1.000,1.000,                                 &
         1.000,1.000,1.000,.980,.983,.982,                                    &
         .970,.945,.913,.880,.852,.832,.820,                                  &
         .817,.816,.810,.795,.777,.765,.764,                                  &
         .759,.730,.664,.579,.500,.446,.416,                                  &
         .400,.393,.390,.390,.391,.391,.390,                                  &
         .388,.384,.380,.375,.366,.350,.324,                                  &
         .291,.260,.234,.214,.200,.190,.184,                                  &
         .180,.176,.173,.170/
    
    ! Log(pres) of efficiencies grid
    DATA pres63/6.90775528,  6.57442194,                                      &
         6.24108859,  5.90775525,  5.57442191,                                &
         5.24108856,  4.90775522,  4.57442188,                                &
         4.24108853,  3.90775519,  3.57442185,                                &
         3.2410885,   2.90775516,  2.57442182,                                &
         2.24108847,  1.90775513,  1.57442179,                                &
         1.24108844,  0.9077551,   0.574421757,                               &
         0.241088414, -0.0922449296,-0.425578273,                             &
         -0.758911616,-1.09224496,  -1.4255783,                               &
         -1.75891165, -2.09224499,  -2.42557833,                              &
         -2.75891168, -3.09224502,  -3.42557836,                              &
         -3.75891171, -4.09224505,  -4.42557839,                              &
         -4.75891174, -5.09224508,  -5.42557842,                              &
         -5.75891177, -6.09224511,  -6.42557845,                              &
         -6.75891179, -7.09224514,  -7.42557848,                              &
         -7.75891182, -8.09224517,  -8.42557851,                              &
         -8.75891185, -9.0922452,   -9.42557854,                              &
         -9.75891188, -10.0922452,  -10.4255786,                              &
         -10.7589119, -11.0922453,  -11.4255786,                              &
         -11.7589119, -12.0922453,  -12.4255786,                              &
         -12.758912,  -13.0922453,  -13.4255787,                              &
         -13.758912/
    
    ! Mlynczak and Solomon 1993 heating efficiency coefficients
    ! equations given only calculate up to 10-4 mb (~130km), therefore
    ! O3 efficiencies kept constant above this height.
    ! O2 efficiencies above this come from Roble.
    
    ! If neutral-neutral chemical heating is off, bulk efficiencies
    ! should be used, else individual o2o3 values applied to solar
    ! energy minus bond energy. See M+S.
    
    ! Table 7 in M+S paper for individual efficiencies
    
    ! Hartley band    - see polynomial
    ! Huggins         - unit
    ! Chappuis        - unit
    ! SRC             - see polynomial
    ! Ly alpha        - 0.95
    ! SRB             - unit
    ! HERZBERG        - unit
    ! H + O3          - 0.6

    ! Calculate heating efficiencies first time through

    ! If you are introducing chemical heating mid-run, you must
    ! re- initialse the efficiencies here.

    IF(chem_heat_on == 1 .AND. chem_heat_on_last == 0)             &
       initialise = 1

    IF(initialise == 1 ) THEN

!!$#ifdef MESSY
!!$       IF (p_parallel_io) &
!!$#endif
!!$       WRITE(6,*) "Initialising Mesospheric heating ..."
       
       ! If self-consistent composition not on, use bulk efficiencies
       !IF(SW(COMPO) /2) THEN
       IF(chem_heat_on == 0) THEN
          
          ! Points have to be ascending for spline to work
          DO n=1,63
             pres63(n)=-pres63(n)
          ENDDO
          
          ! Fit spline for efficiencies on old grid
          CALL SPLINE1D(pres63,o2eff63,63,0._dp,0._dp,splft1,.TRUE.)
          CALL SPLINE1D(pres63,o3eff63,63,0._dp,0._dp,splft2,.TRUE.)
          
          minn=-1
          maxn=-1
          
          DO n = ht_dim, 1, -1
             
             eff1=1.
             eff2=1.
             
             ! Points have to be ascending for spline to work
             pval=-log(pres(n))
             
             IF((pres63(63)-pval) > 0. .AND.                                  &
                  pval > pres63(1)) THEN
                
                CALL SPLINT1D(pres63,o2eff63,splft1,63,pval,eff1,status)
#ifndef MESSY
                IF (status > 0) RETURN
#else
                IF (status > 0) THEN
                   write(*,*) 'SPLINT1D',m,l,status
                   RETURN
                ENDIF
#endif
                CALL SPLINT1D(pres63,o3eff63,splft2,63,pval,eff2,status)
#ifndef MESSY
                IF (status > 0) RETURN
#else
                IF (status > 0) THEN
                   write(*,*) 'SPLINT1D',m,l,status
                   RETURN
                ENDIF
#endif
               
             ENDIF
             
             O2eff_bulk(n)=eff1
             O3eff_bulk(n)=eff2
             
             ! Make individual efficiency 1.
             htlyeff(n)=1.
             
             ! Find the points where we went out of range
             IF(pval > pres63(63)) maxn=n
             IF(pval < pres63(1) .AND. minn == -1) minn=n
             
             
          ENDDO

          ! Lets tidy up where we went out of bounds.
          IF(maxn > 1) THEN
             DO n=maxn, ht_dim
                O2eff_bulk(n)=O2eff_bulk(n-1)
                O3eff_bulk(n)=O3eff_bulk(n-1)
             ENDDO
          ENDIF
          IF(minn > 1) THEN
             DO n=minn, 1, -1
                O2eff_bulk(n)=O2eff_bulk(n+1)
                O3eff_bulk(n)=O3eff_bulk(n+1)
             ENDDO
          ENDIF

          ! Use individual efficiencies and don't include bond energy if 
          ! composition on ...

       ELSE
          
          ! Evaluate hartley efficiency polynomial fit
          DO n=1,ht_dim
             pval=(pres(n))*(1.e-2) ! in mb
             IF(pval <= 10. .AND. pval > 1.) THEN
                htlyeff(n)=1.0
             ELSE
                IF(pval <= 1 .AND. pval > 1.e-2) THEN ! M+S eff
                   ex=LOG10(pval)+1.
                   htlyeff(n)=0.926210+0.133960*ex-0.076863*(ex**2)+          &
                        0.006897*(ex**3)
                ELSE
                   ex=LOG10(pval)+3
                   htlyeff(n)=0.926210+0.133960*ex-                           &
                        0.076863*(ex**2)+0.006897*(ex**3)
                   htlyeff(n)=0.669650-0.009682*ex+                           &
                        0.033093*(ex**2)+0.017938*(ex**3)
                   ! srceff(n)=0.75349+0.0036*(ex)+0.059468*(ex**2)-
                   ! &                    0.022795*(ex**3)
                ENDIF
                
                IF(pval <= 1.e-4) THEN
                   htlyeff(n)=htlyeff(n-1)
                ENDIF
                
             ENDIF
             
          ENDDO
          
          O2eff_bulk(:)=1.
          O3eff_bulk(:)=1.
          
       ENDIF

! ALD test: reduce O3 heating efficiency abve 70km if SW(compo) = 2
! Self consistent O3 concentrations in the mesosphere compare well with
! lit but result in unrealisticly large heating. Needs investigating.

       IF (l_reduceO3heating) THEN   
         DO n = 1,ht_dim
            if(n.ge.22) then
              if(n.eq.22)write(6,*)"ALD TEST ** Reducing O3 heating "
              O3eff_bulk(n) = O3eff_bulk(n)*0.2
            endif
         ENDDO
       ENDIF
 
       initialise = 0
             
    ENDIF ! initialise

    ! Solar variability factor
    fmxfmn=(F107-71.)/(220-71.)
    
    ! Strobez parameterisations don't include bond energy. Therefore for 
    ! each range we need to add this back in if we're not calculating 
    ! composition to calculate chemical heating. To estimate the factor 
    ! in each, we compare the bond energy with that of an average photon in 
    ! that range.
    ! ALD apr 08: I'm not convinced about this. The Strobel 78 paper and 
    ! B&S p138 suggest that bond energy IS included if using the
    ! calculation of total heating rate (rather than net heating rate)
    ! which we are. I'm keeping the factors in though as the profile do
    ! look better when they're included.
    ! NB also that cmat1 has a factor of 1.8 which accounts for 'not
    ! considering bond energy' - i think this is .....wrong.

   
    ! We also need to apply Mlynczak and Solomon bulk heating efficiencies 
    ! rather than efficiencies for individual processes. If composition is 
    ! on, need to use individual efficiencies, not subtract bond energy, 
    ! and include exothermic chemical heating.
    
    ! fac=bond_energy/(h*c/wavelength)
    
    ! O3BOND=1.10*1.602e-19
    ! O2BOND=5.12*1.602e-19
    
    ! srband 175-205nm   average=190nm      
    !                    fac= 5.12*1.602e-19/(6.6256e-34*2.9979e8/190.e-9)
    ! hartley 200-300nm  average=250nm      
    !                    fac= 1.10*1.602e-19/(6.6256e-34*2.9979e8/250.e-9)
    ! huggins 300-335nm  average=317.5nm    
    !                    fac= 1.10*1.602e-19/(6.6256e-34*2.9979e8/317.5e-9)
    ! chappuis 450-850nm average=650nm      
    !                    fac= 1.10*1.602e-19/(6.6256e-34*2.9979e8/650.e-9)
    ! herzberg (205 - 242nm) includes bond energy, so no factor required
    
    ! We need to add back in bond energy if no composition
    ! IF(SW(COMPO) /= 2) THEN
    IF(chem_heat_on == 0)THEN     
       
       srband_fac=0.784591675
       hartley_fac=0.221795559
       huggins_fac=0.281680375
       chappuis_fac=0.576668501
       
       ! We calculate without bond energy when composition on, since bond 
       ! energy included in chemical heating
    ELSE
      
       srband_fac=0.
       hartley_fac=0.
       huggins_fac=0.
       chappuis_fac=0.
       
    ENDIF

    ! Initialise flux
    do nn = 1,11
       dfluxh(nn)= fluxh(nn)*eccentric*(1.+0.03*fmxfmn)
    enddo

    
    DO n = ht_dim, 1, -1
       
       heatrates(n, m,O3_chap_heat_l,l)=0.
       heatrates(n, m,O3_hart_heat_l,l)=0.
       heatrates(n, m,O3_hugg_heat_l,l)=0.
       heatrates(n, m,O3_herz_heat_l,l)=0.
       heatrates(n, m,O2_meso_heat_l,l)=0.
       jx(n,m,ip_O1D,l) = 0.   

       srband = 0.
       chappuis = 0.
       hartley = 0.
       huggins = 0.
       herzberg = 0.
       meso_heat = 0.
       ozheat = 0.
       o2heat = 0.
       
       ! Find csza, with Chapmann grazing incidence
       csza_ch=CHAPMANN(csza, scht1d(n), r0_eff(n))

       ! ALD : Find csza, with Chapmann grazing incidence for each
       ! constituent.
       ! first scale heights m:

       scht1d_o3(n)=((GSCON*temp1d(n))/(mmw_gas1*GRAV0))/3
       csza_o3=CHAPMANN(csza,scht1d_o3(n), r0_eff(n))
       
       scht1d_o2(n)=((GSCON*temp1d(n))/(mmw_gas2*GRAV0))
       csza_o2=CHAPMANN(csza,scht1d_o2(n), r0_eff(n))

       scht1d_h2o(n)=((GSCON*temp1d(n))/(18.*GRAV0))
       csza_h2o=CHAPMANN(csza,scht1d_h2o(n), r0_eff(n))

       ! Is the sun up? Note, this is inside height loop due to Chapmann 
       ! function
       IF(csza_ch > 0.) THEN
          
          ! O2 concentration per cm3
          o2cm3=gas2_1d(n)*1.e-6
          
          ! O3 concentration per cm3
          o3cm3=gas_O3_1d(n)*nden1d(n)*1.e-6
          
          ! H2O concentration per cm3
          h2ocm3=gas_H2O_1d(n)*nden1d(n)*1.e-6

          IF(n == ht_dim) THEN
             
             ! Column approximation above top of atmosphere
             ! (scale height in cm)
             colo2(n)=o2cm3*scht1d_o2(n)*1.e2*(1./csza_o2) 
             
             ! There's no O3 up here
             colo3(n)=1.
             
             ! There's no H2O up here
             colh2o(n)=1.

              Z=(ht1d(n)-ht1d(n-1))*1.e2
         ELSE
             
             ! Level thickness in cm
             Z=(ht1d(n+1)-ht1d(n))*1.e2
             
             ! Integrate Column O2
             colo2(n)=colo2(n+1)+o2cm3*Z*(1./csza_o2)
             
             ! Integrate Column O3
             colo3(n)=colo3(n+1)+o3cm3*Z*(1./csza_o3)
             
             ! Integrate Column H2O
             colh2o(n)=colh2o(n+1)+h2ocm3*Z*(1./csza_h2o)

          ENDIF
          
          ! Calculate O2 Schumunn Runge band
          IF (colo2(n) < 1E18) THEN
             srband=o2cm3*2.43E-19
          ELSE
             srband=(o2cm3)/(0.66*colo2(n)+3.44E9*SQRT(colo2(n)))
          ENDIF
          
          ! Calculate O3 heating. Note these are with bond energy subtracted.
          IF (colo3(n) > 0.01) THEN
             
             chappuis=O3cm3*1.0545E-15*                                       &
                  EXP(-2.85E-21*colo3(n))
             
             hartley=O3cm3*4.8048E-14*                                        &
                  EXP(-8.8E-18*colo3(n))
             
             huggins=O3cm3/colo3(n)*(4.66E3-7.8E2*                            &
                  EXP(-1.77E-19*colo3(n))                                     &
                  - 3.88E3*EXP(-4.22E-18*colo3(n)))

             herzberg=1.5e3*(6.6e-24*o2cm3 + 4.9e-18*o3cm3)*                  &
                  EXP(-6.6e-24*colo2(n) - 4.9e-18*colo3(n))
            
          ENDIF
          
          ! Add back bond energy
          ! - Convert to J/s/m3 rather than 0.1J/s/m-3 (/10. as per B+S) 
          ! - Solar variability
          ! - Orbit eccentricity

          srband =(1.+srband_fac)*                                            &
               (srband/10.)*(1.+0.11*fmxfmn)*eccentric
          
          chappuis =(1.+chappuis_fac)*                                        &
               (chappuis/10.)*(1+0.03*fmxfmn)*eccentric
          
          hartley = htlyeff(n)*(1.+hartley_fac)*                              &
               (hartley/10.)*(1+0.03*fmxfmn)*eccentric
          
          huggins = (1.+huggins_fac)*                                         &
               (huggins/10.)*(1+0.03*fmxfmn)*eccentric
          
          herzberg = (herzberg/10.)*(1+0.04*fmxfmn)*eccentric

          if(chappuis.lt.0.) chappuis = 0. 
          if(hartley.lt.0.)  hartley  = 0.
          if(huggins.lt.0.)  huggins  = 0.
          if(herzberg.lt.0.) herzberg = 0.
          if(srband.lt.0.)   srband   = 0.

! ALD TEST: try setting herzberg as it is done in cmat1
!          herzberg = herzberg/(1.e-4*86400)
! this is done differently in cmat1 but herberg doesn't contribute 
! much to heating anyway so I'm not going to worry about it now

          ! ozheat=o3eff_bulk(n)*(chappuis+hartley+huggins) + herzberg
          chappuis = o3eff_bulk(n)*chappuis
          hartley  = o3eff_bulk(n)*hartley
          huggins  = o3eff_bulk(n)*huggins
 
          o2heat   = o2eff_bulk(n)*srband

          ! ozheat=o3eff_bulk(n)*(chappuis+hartley+huggins+herzberg)

          !  meso_heat=chappuis + hartley + huggins + herzberg + o2heat

          !  if(m.eq.46.and.l.eq.10.or.l.eq.1.and.(n.ge.15.or.n.le.55))&
          !  write(6,*)' o2 o3 heat ',n,l, o2heat,meso_heat

          ! Very Crude parameterized heating for 15-30km to account for 
          ! scattering and H2O and IR. Basically rate is held constant
          ! at 30km. This means the heating rate in Kday-1 falls off with
          ! density. This is offset a little with an exponent, being
          ! a function of the diostance from the fade off height. Very
          ! crude, but it gives a resonably smooth transition, and global
          ! heating/cooling profile in line with the literature.
          !

          IF(ht1d(n) < 30e3 .AND. n < ht_dim) THEN

             ! The factor after csza_ch is peak K/day
             !ozheat = ozheat +                                                &
             !         csza_ch*3.*EXP(-((ht1d(n) - 18e3)**2.)/(2.*(5.e3)**2.))
             !ozheat = ozheat*den1d(n)/82.

             IF(fade == 0) fade = n+1

             chappuis = heatrates(fade, m,O3_chap_heat_l,l)*den1d(n+1)
             hartley  = heatrates(fade, m,O3_hart_heat_l,l)*den1d(n+1)
             huggins  = heatrates(fade, m,O3_hugg_heat_l,l)*den1d(n+1)
             herzberg = heatrates(fade, m,O3_herz_heat_l,l)*den1d(n+1)
             chappuis = chappuis*EXP(-(ht1d(fade-1) - ht1d(n))/10.e3)
             hartley  = hartley*EXP(-(ht1d(fade-1) - ht1d(n))/10.e3)
             huggins  = huggins*EXP(-(ht1d(fade-1) - ht1d(n))/10.e3)
             herzberg = herzberg*EXP(-(ht1d(fade-1) - ht1d(n))/10.e3)
          ENDIF
          ! 
          ! Dissociation rate coefficients. See rates.txt.README for 
          ! details

          IF (colo2(n).LT.1.E19) THEN
             jx(n,m,ip_O2,l)=(1.+0.11*fmxfmn)*1.1E-7*              &
                  EXP(-1.97E-10*(colo2(n)**0.522))
          ELSE
             jx(n,m,ip_O2,l)=(1.+0.11*fmxfmn)*1.45E8*              &
                  (colo2(n)**(-0.83))
          ENDIF

          ! ALD I've multiplied jherz_pr photodissociation rate by cos(solar 
          ! zenith angle) to give some transition between day and night.
          ! The value of 8.e-10 is from B&S and is for and overhead sun
          ! i.e. cos(sza) of 1. Also changed eqn for tail below 60km to match
          ! profile in B&S

          IF (ht1d(n) > 60.e3) THEN
             jx(n,m,ip_O2,l)= jx(n,m,ip_O2,l) + 8.e-10*(csza_ch)
          ELSE
             jx(n,m,ip_O2,l)= jx(n,m,ip_O2,l) + (8.e-10 - ((60.e3-ht1d(n))  &
                                        /15.e3)*(3.9e-10))*(csza_ch)  
          ENDIF


          DO  nn=1,11
             section=htly(nn)*1.E-19*o3cm3*(1./csza_o3)*Z
             tmp = dfluxh(nn)*EXP(-section) 
             jx(n,m,ip_O1D,l) =jx(n,m,ip_O1D,l) +          &
                  tmp*htly(nn)*1.e-19*10.
             dfluxh(nn)=tmp
          ENDDO


          !   Thermosphere/mesosphere optically thin for Chappuis and Huggins
          !   bands, therefore assumed as zero optical depth values given in B+S
          !   Scattering left unaccounted for now (enhances bottom 15kms)
          !   JHUG=1.24e-4     JChap=4.4e-4
          
          jx(n,m,ip_O3P,l) = 1.24e-4+4.4e-4     

          !   First part Ly alpha, second SRB (see B+S,nicolet 84)
          !
          !   JH2Oa  -> H + OH     LYA
          !   JH2Ob  -> H2 + (O1d) SRB
          !
          ! Paramaterization different (factor -3.4e-6 instead of 
          ! 1.e-7 in SRB, and -1.e-14 instead of -4.4e-19), to 
          ! get match with B+S profile.
          !
          ! ALD: I think the column density should be that of O2 here. 
          ! as in B&S. Previously colh2o(n) was used here

          ! ALD NB: jh2oa looks too small.
          ! may be due to low O2 col restulting from too much O again?

          jx(n,m,ip_H2O,l) =1.e-6*(4.+ 2.5*fmxfmn)*              &
                                      EXP(-4.4e-19*colo2(n)**0.917)

          jx(n,m,ip_H2O1D,l) =(1.2e-6*EXP(-1.e-7*colo2(n)**0.35))

          ! calculate H2O2 rate: H2O2 + hv = 2OH for wavelength range
          ! 260nm to 350nm, temps 200 to 400K.

          CALL H2O2_photorate(n, m, l, ht_dim, temp1d, eccentric,                  &
                              gas2_1d, nden1d, R0_eff, csza,               &
                              ht1d, gas_O3_1d, gas_H2O2_1d, GRAV0)


          ! ALD temporarily setting these jx here so  can plot them.
          ! BUT they are 1d so can return them to chemistry_gen when finished
          jx(n,m,ip_NO2,l) =  8.50d-03 ! NO2 + hv (190-410nm) -> NO + O
          jx(n,m,ip_N2O,l) =  9.00d-07 ! N2O + hv (>175nm) -> N2 + O{1d}
          jx(n,m,ip_NO2O,l) = 1.56d-01 ! NO3 + hv (<580nm) -> NO2 + O
          jx(n,m,ip_NOO2,l) = 2.01d-02 ! NO3 + hv (>580nm) -> NO + O2 

       ELSE

          jx(n,m,ip_H2O,l) = 0.
          jx(n,m,ip_H2O1D,l) = 0.
          jx(n,m,ip_O3P,l)  = 0.
          jx(n,m,ip_O1D,l)  = 0.
          jx(n,m,ip_O2,l)  = 0.
          jx(n,m,ip_H2O2,l) = 0.

          jx(n,m,ip_NO2,l) =  0.
          jx(n,m,ip_N2O,l) =  0. 
          jx(n,m,ip_NO2O,l) = 0. 
          jx(n,m,ip_NOO2,l) = 0.  

       ENDIF


       ! check for any negative values

          IF(jx(n,m,ip_H2O,l).lt.1.e-20)                          &
             jx(n,m,ip_H2O,l) = 0.
          IF(jx(n,m,ip_H2O1D,l).lt.0.)                              &
             jx(n,m,ip_H2O1D,l) = 0.
          IF(jx(n,m,ip_O3P,l).lt.0.)                               &
             jx(n,m,ip_O3P,l) = 0.
          IF(jx(n,m,ip_O1D,l).lt.0.)                               &
             jx(n,m,ip_O1D,l) = 0.
          IF(jx(n,m,ip_O2,l).lt.0.)                               &
             jx(n,m,ip_O2,l) = 0.
!!$          IF(jx(n,m,jherz_pr,l).lt.0.)                              &
!!$             jx(n,m,jherz_pr,l) = jx(n+1,m,jherz_pr,l)*0.4

       heatrates(n,m,O3_chap_heat_l,l) = chappuis/(den1d(n))
       heatrates(n,m,O3_hart_heat_l,l) = hartley/(den1d(n))
       heatrates(n,m,O3_hugg_heat_l,l) = huggins/(den1d(n))
       heatrates(n,m,O3_herz_heat_l,l) = herzberg/(den1d(n))
       heatrates(n,m,O2_meso_heat_l,l) = o2heat/(den1d(n))       

     !  if(m.eq.46.and.l.eq.10) &
     !         write(6,*)'o3heat',n,(heatrates(n,m,l,O3_hart_heat_l)+ &
     !                    heatrates(n,m,l,O3_hugg_heat_l) + &
     !                    heatrates(n,m,l,O3_herz_heat_l) + &
     !                    heatrates(n,m,l,O3_chap_heat_l)) !*86400/cp1d(n) 

    ENDDO   

    status = 0 ! no error 

  END SUBROUTINE MesoSol

  !========================================================
  !=                   H2O2_photorate                     =
  !========================================================

  !r
  !r H2O2_photorate  Routine to calculate photodissociation
  !r                  rate for H2O2 + hv = 2OH
  !r                  For wavelength range 260nm to 350nm
  !r                  and temperatures from 200 to 400K
  !r
  !r                  Solar flux data is from V24a of SOLAR2000
  !r                  and does not vary with F107 at these
  !r                  wavelengths
  !r

  SUBROUTINE H2O2_photorate(n, m, l, ht_dim, temp1d, eccentric,               &
                           gas2_1d, nden1d, R0_eff, csza,   &
                           ht1d, gas_O3_1d, gas_H2O2_1d, GRAV0)

    INTEGER,  INTENT(IN)  :: n, m, l, ht_dim
    REAL(dp), INTENT(IN)  :: eccentric , csza
    REAL(dp), INTENT(IN)  :: temp1d(ht_dim), ht1d(ht_dim)
    REAL(dp), INTENT(IN)  :: gas2_1d(ht_dim), nden1d(ht_dim)
    REAL(dp), INTENT(IN)  :: R0_eff(ht_dim)
    REAL(dp), INTENT(IN)  :: gas_O3_1d(ht_dim), gas_H2O2_1d(ht_dim)
    REAL(dp), INTENT(IN)  :: GRAV0

    !Local variables
    INTEGER    :: countw, countc

    ! Generic flux array used for calculation
    REAL(dp) :: local_flux_fuv

    ! Cross sections of O2
    REAL(dp) :: CSAO2_fuv(NWAVES_FUV)

    ! Ratio of nighttime to daytime ionization
    REAL(dp) :: nightfac = 1.e-6
    REAL(dp) :: rnight_h2o2

    ! column densities of level
    REAL(dp) :: WO2, WO3, WH2O2 

    ! Integrated column densities
    REAL(dp) :: TAUO2, TAUO3, TAUH2O2, TAU, attenuation

    ! Consituent scale heights
    REAL(dp) :: HO(Ht_dim), HO2(Ht_dim), HO3(Ht_dim),          &
                  HH2O2(Ht_dim)

    ! Photodissociation rate coeff for H2O2 + hv = 2OH
    REAL(dp) :: JH2O2_loc
    REAL(dp) :: jh2o2(ht_dim)

    ! Total height R0+ht1d for Chapmann function
    REAL(dp) :: Z(ht_dim)

    ! Constiuent densities cm-3
    REAL(dp) :: O2ALL(ht_dim),O3ALL(ht_dim), H2O2ALL(Ht_dim)

    ! acceleration due to gravity, seczenang for o3, o2 and h2o2
    REAL(dp) :: gravy(ht_dim), SECO2, SECO3, SECH2O2

    ! Height between levels
    REAL(dp) :: zht


    ! Calculate H2O2 absorption cross section
    CALL H2O2_CROSS_SECTION(n, NWAVES_fuv, WAVELS_fuv, temp1d, CSH2O2)

    ! Reset columns
    WO2=0.
    WO3=0.
    WH2O2=0.

    ! number densities cm-3
    O2all(n)=gas2_1d(n)*1.e-6
    O3all(n)=gas_O3_1d(n)*nden1d(n)*1.e-6
    H2O2all(n)=gas_H2O2_1d(n)*nden1d(n)*1.e-6

    ! Set gravity
    gravy(n) = GRAV0

    ! Scale heights, in cm
    HO(n)=(GSCON*temp1d(n)/16./gravy(n))*1.E2
    HO2(n)=HO(n)/2.
    HO3(n)=HO(n)/3.
    HH2O2(n)=HO(n)/7.*8.

    ! Set total height in cm
    z(n)=(R0_eff(n))*1.E2

    jh2o2(n)  = 0.

    ! calculate sec zenith angle, incorporating Chapmann grazing incidence
    ! function
    SECO2=1./CHAPMANN(csza, HO2(n), Z(n))
    SECO3=1./CHAPMANN(csza, HO3(n), Z(n))
    SECH2O2=1./CHAPMANN(csza, HH2O2(n), Z(n))

    rnight_h2o2=1.

    ! The stars are out ...
    IF(SECO2 < 0) THEN
       rnight_h2o2=nightfac
    ENDIF


    ! Integrate Column densities
    IF(n == ht_dim) THEN
       WO2=O2all(n)*HO2(n)*SECO2
       WO3=O3all(n)*HO3(n)*SECO3
       WH2O2=H2O2all(n)*HH2O2(n)*SECH2O2
    ELSE
       ! distance between 2 levels in cm
       ZHT=(ht1d(n+1)-ht1d(n))*1.E2
       WO2=WO2+O2ALL(n)*ZHT*SECO2
       WO3=WO3+O3ALL(n)*ZHT*SECO3
       WH2O2=WH2O2+H2O2ALL(n)*ZHT*SECH2O2
    ENDIF

    ! Only bother if the sun's up
    IF(SECO2 >0 .OR. SECO3 >0 .OR. SECH2O2 >0) THEN

        JH2O2_loc   = 0.

          !$OMP PARALLEL PRIVATE(countw,TAUO3,TAUO2,TAUH2O2,TAU,attenuation,    &
          !$OMP local_flux_fuv)                                                 &
          !$OMP SHARED(CSAO3, CSAO2_fuv, CSH2O2,WO3,WO2,WH2O2,                  &
          !$OMP rnight_h2o2,WAVELS_fuv)                                         &
          !$OMP REDUCTION(+:JH2O2_loc)
          !$OMP DO

        ! Loop over wavelength
        WAVELENGTH1: DO countw=1, NWAVES_fuv

            CSAO2_fuv(countw) = 1e-24
            TAUO2=(CSAO2_fuv(countw))*WO2
            TAUO3=(CSAO3(countw))*WO3
            TAUH2O2=(CSH2O2(n,countw))*WH2O2
            TAU=(TAUO2+TAUO3+TAUH2O2)

            attenuation=0.
            IF (tau <= 100. .AND. tau /= 0) attenuation = exp(-tau)

            local_flux_fuv=attenuation*S2K_FLUX_fuv(countw)*eccentric
            IF(local_flux_fuv < 1 .OR. tau == 0) local_flux_fuv=0.

            JH2O2_loc = JH2O2_loc+(local_flux_fuv)*CSH2O2(n,countw)

        ENDDO WAVELENGTH1
        !$OMP END DO
        !$OMP END PARALLEL

    ELSE

        ! Nightime rates
        JH2O2_loc   = 0.

    ENDIF

    ! Populate photorate coeffcient arrays
    JH2O2(n) = JH2O2_loc

   jx(n,m,ip_H2O2,l) =  JH2O2_loc*rnight_h2o2

    RETURN
  END SUBROUTINE H2O2_photorate


  !========================================================
  !                   H2O2_cross_section
  !========================================================

  !r ALD 26 Nov 2007
  !r Routine to calculate absorption cross section for h2o2
  !r Valid for (260nm <= hv <= 350nm) and
  !r 200K <= T <+ 400K
  !r Numerical expression for Absorption cross section of h2o2
  !r taken from JPL 2003, Sanders et al. 2003.


  SUBROUTINE H2O2_cross_section(n, NWAVES_fuv, WAVELS_fuv,         &
                                                temp1d, CSH2O2)

    ! passed in: number of wavelengths, 
    !            wavelengths(m), temperatures(K)
    INTEGER,  INTENT(IN)  :: n
    INTEGER,  INTENT(IN)  :: NWAVES_fuv
    REAL(dp), INTENT(IN)  :: WAVELS_fuv(:)
    REAL(dp), INTENT(IN)  :: temp1d(:)

    ! output h2o2 cross section (cm-2)
    REAL(dp), INTENT(INOUT) :: CSH2O2(:,:)

    ! local variables
    INTEGER                 :: i, j
    REAL(dp)              :: Chi, sumA, sumB
    REAL(dp)              :: constant_A(8)
    REAL(dp)              :: constant_B(5)

    DATA constant_A/ 6.4761e4, -9.2170972e2, 4.535649,             &
                     -4.4589016e-3, -4.035101e-5, 1.6878206e-7,    &
                     -2.652014e-10, 1.5534675e-13/

    DATA constant_B/ 6.8123e3, -5.1351e1, 1.1522e-1,               &
                     -3.0493e-5, -1.0924e-7/


    Chi = 1/(1 + exp(-1265./temp1d(n)))

    DO i = 1, NWAVES_fuv

            CSH2O2(n,i) = 0.
            sumA = 0.

            DO j = 1,8
                  sumA = sumA + constant_A(j)*(WAVELS_fuv(i))**(j-1)
            ENDDO

            sumB = 0.

            DO j = 1,5
                 sumB = sumB + constant_B(j)*(WAVELS_fuv(i))**(j-1)
            ENDDO

            CSH2O2(n,i) = (chi*sumA + (1-chi)*sumB)/1.e21

    ENDDO

  END SUBROUTINE H2O2_cross_section

  !========================================================
  !=                   ThermoHeat                         =
  !========================================================

  !r
  !r thermo_sol    Routine to calculate thermospheric heating and jx. 
  !r               Can be run using solar2000 data, or data used in CMAT1. See
  !r               comments in S2K_fluxes for data description, and CMAT thesis
  !r               for CMAT data. S2k_fluxes formulated and coded by Alison 
  !r               Dobbin.
  !r
  !r This routine could be tuned. It integrates over spectrum even at night.
  !r Probably overkill, one point would do. mjh
  !r
  !r

  SUBROUTINE ThermoSol(m,l,ht_dim, gas1_1d,gas2_1d, gas3_1d, gas_O3_1d,    &
       F107,csza, ht1d, heatrates,temp1d,                                  &
       eccentric, cp1d, scht1d, pres, den1d, R0_eff,                       &
       nden1d, GRAV0, iou1, iou2)
    
    
    INTEGER,  INTENT(IN)  :: m, l, ht_dim       
    REAL(dp), INTENT(IN)  :: csza, f107, eccentric
    REAL(dp), INTENT(IN)  :: ht1d(ht_dim), temp1d(ht_dim)
    REAL(dp), INTENT(IN)  :: gas1_1d(ht_dim), gas2_1d(ht_dim)
    REAL(dp), INTENT(IN)  :: gas3_1d(ht_dim), gas_O3_1d(ht_dim)
    REAL(dp), INTENT(IN)  :: den1d(ht_dim), pres(ht_dim), cp1d(ht_dim)
    REAL(dp), INTENT(IN)  :: scht1d(ht_dim), R0_eff(ht_dim)
    REAL(dp), INTENT(IN)  :: nden1d(ht_dim)
    REAL(dp), INTENT(IN)  :: GRAV0
    INTEGER,  INTENT(IN)  :: iou1, iou2

    REAL(dp), POINTER  :: heatrates(:,:,:,:)
 
    INTEGER :: initialise = 1
    INTEGER :: status = 1
    
    ! See data statement below
    INTEGER :: FLUX_SOURCE
    
    ! NWAVES is largey than the number of wave bins returrned from
    ! cmat1_fluxes, or s2k_fluxes.
    ! the arrays are populated as
    ! Xrays in 1 -> NWAVES_XRAYS,
    ! EUV in NWAVES_XRAYS+1->NWAVES_EUV+NWAVES_XRAYS,
    ! SRC in NWAVES_EUV+NWAVES_XRAYS + 1->
    ! NWAVES_EUV+NWAVES_XRAYS+NWAVES_SRC
    
    ! NWAVES_EUV, NWAVES_SRC, NWAVES_XRAY are returned from
    ! flux routines
    
    INTEGER :: NWAVES, NWAVES_EUV, NWAVES_SRC, NWAVES_XRAY
    PARAMETER(NWAVES=200)
    
    ! s2k flux and wavelengths
    REAL(dp) :: S2K_FLUX(NWAVES), WAVELS(NWAVES)
    
    ! s2k cross sections
    REAL(dp) :: CSAO(NWAVES),CSIO(NWAVES), CSAN2(NWAVES),                   &
         CSIN2(NWAVES), CSAO2(NWAVES), CSIO2(NWAVES)
    
    ! Generic flux array used for calculation
    REAL(dp) :: FLUX(NWAVES)
    REAL(dp) :: local_flux
    
    ! Lyman-a flux
    REAL(dp) :: lyman_a_flux
    INTEGER :: lyman_a_num
    
    ! Bounds variables for efficiency interpolation
    INTEGER :: maxn, minn
    
    ! Save these variables between calls. Need this if not compiled with
    ! static flag.
    SAVE CSAO,CSIO, CSAN2,CSIN2, CSAO2,                                       &
         CSIO2, FLUX, WAVELS, lyman_a_flux,                                   &
         lyman_a_num, minn, maxn, NWAVES_EUV, NWAVES_SRC, NWAVES_XRAY
    
    ! Raw s2k fluxes and cross-sections in input file
    REAL(dp) :: S2K_ALLFLUXES(NWAVES,12), XSECTIONS(NWAVES,7)
    
    ! Units conversion for s2k
    REAL(dp) :: unit_conv
    
    ! Old CMAT pressure levels
    INTEGER :: TOPLVL
    PARAMETER (TOPLVL=63)
    
    ! Plank*c. See below.
    REAL(dp) :: PCC
    
    ! column densities of level
    REAL(dp) :: WO, WO2, WN2, WO3
    
    ! Integrated column densities
    REAL(dp) :: TAUO, TAUO2, TAUN2, TAU, attenuation
    
    ! Some countv variables
    INTEGER :: wav, countv
    
    ! Consituent scale heights
    REAL(dp) :: HO(ht_dim),HO2(ht_dim),HN2(ht_dim),HO3(ht_dim)
    
    ! Total height R0+ht1d for Chapmann function
    REAL(dp) :: Z(ht_dim)
    
    ! EUV total heating efficiency
    REAL(dp) :: euveff63(TOPLVL)
    
    ! Pressure coords of efficiencies grid
    REAL(dp) :: pres63(TOPLVL), pval
    
    ! Efficiencies fit variables
    REAL(dp) :: splft1(63), splft2(63)
    
    ! ionisation rate coefficients
    REAL(dp) :: PEUVO(ht_dim), PEUVO2(ht_dim), PEUVN2(ht_dim)
    
    ! N2 pre-dissociation rate coefficient
    REAL(dp) :: pre(ht_dim)
    
    ! Lyman-a and SRC O2 dissociation rate coefficients
    REAL(dp) :: JLYA(ht_dim), jsrc(ht_dim)

    ! EUV and SRC heating
    REAL(dp) :: paeuv(ht_dim), src1(ht_dim)
    
    ! Constiuent densities cm-3
    REAL(dp) ::  OALL(ht_dim),O2ALL(ht_dim),N2ALL(ht_dim), O3ALL(ht_dim)
    
    ! SRC bulk heating efficiency
    REAL(dp) :: src1eff63(TOPLVL)
    
    ! Temporary efficiency variables
    REAL(dp) :: eff, eff1, eff2
    
    ! acceleration due to gravity, seczenang for o3,o2,o, n2 
    REAL(dp) :: gravy(ht_dim),seco3,seco2,seco, secn2
    
    ! Ratio of nighttime to daytime ionization
    REAL(dp) :: rnight_o,rnight_o2, rnight_n2, nightfac
    
    ! Height between levels
    REAL(dp) :: zht

    INTEGER    :: n

    REAL(dp) PEUVO2_loc, PEUVO_loc, PEUVN2_loc, pre_loc, PAEUV_loc, SRC1_loc 
    REAL(dp) JSRC_loc, JLYA_loc

    REAL(dp) :: last_f107 = 0.
 
    ! Log(pres) of efficiencies grid
    DATA pres63/6.90775528,  6.57442194,                                      &
         6.24108859,  5.90775525,  5.57442191,                                &
         5.24108856,  4.90775522,  4.57442188,                                &
         4.24108853,  3.90775519,  3.57442185,                                &
         3.2410885,   2.90775516,  2.57442182,                                &
         2.24108847,  1.90775513,  1.57442179,                                &
         1.24108844,  0.9077551,   0.574421757,                               &
         0.241088414, -0.0922449296,-0.425578273,                             &
         -0.758911616,-1.09224496,  -1.4255783,                               &
         -1.75891165, -2.09224499,  -2.42557833,                              &
         -2.75891168, -3.09224502,  -3.42557836,                              &
         -3.75891171, -4.09224505,  -4.42557839,                              &
         -4.75891174, -5.09224508,  -5.42557842,                              &
         -5.75891177, -6.09224511,  -6.42557845,                              &
         -6.75891179, -7.09224514,  -7.42557848,                              &
         -7.75891182, -8.09224517,  -8.42557851,                              &
         -8.75891185, -9.0922452,   -9.42557854,                              &
         -9.75891188, -10.0922452,  -10.4255786,                              &
         -10.7589119, -11.0922453,  -11.4255786,                              &
         -11.7589119, -12.0922453,  -12.4255786,                              &
         -12.758912,  -13.0922453,  -13.4255787,                              &
         -13.758912/
    
    ! Heating efficiencies, after Roble
    DATA euveff63/1.00000,   1.00000,   1.00000,   1.00000,   1.00000,        &
         1.00000,   1.00000,   1.00000,   1.00000,   1.00000,                 &
         1.00000,   1.00000,   1.00000,   1.00000,   1.00000,                 &
         1.00000,   1.00000,   1.00000,   1.00000,   1.00000,                 &
         1.00000,   1.00000,   1.00000,   1.00000,   1.00000,                 &
         1.00000,   1.00000,   1.00000,   1.00000,   1.00000,                 &
         1.00000,   1.00000,   1.00000,   1.00000,   1.00000,                 &
         1.00000,   1.00000,   1.00000,   1.03040,   0.93341,                 &
         0.77392,   0.62373,   0.55243,   0.56369,   0.60571,                 &
         0.62715,   0.61404,   0.58380,   0.55422,   0.53229,                 &
         0.51595,   0.50270,   0.48777,   0.46444,   0.42657,                 &
         0.37873,   0.33454,   0.30718,   0.29809,   0.29883,                 &
         0.30095,   0.30112,   0.30034/
    
    
    
    ! SRC efficiency only need be used if composition off. The 0.3-0.4
    ! comes from loss due to diss (and subsequent transport down) and airglow.
    ! whilst the airglow isn't represented brilliantly, bond subtraction
    ! is covered if composition on
    
    DATA src1eff63/  0.28249,   0.28249,   0.28249,   0.28249,                &
         0.28249,   0.28249,   0.28249,   0.28249,                            &
         0.28249,   0.28249,   0.28249,   0.28249,                            &
         0.28249,   0.28249,   0.28249,   0.28249,                            &
         0.28249,   0.28249,   0.28249,   0.28249,                            &
         0.28249,   0.28249,   0.28249,   0.28249,                            &
         0.28249,   0.28249,   0.28249,   0.28249,                            &
         0.28249,   0.28249,   0.28249,   0.28776,                            &
         0.29549,   0.30716,   0.32284,   0.34137,                            &
         0.36150,   0.38078,   0.39573,   0.40311,                            &
         0.40356,   0.40102,   0.39934,   0.39951,                            &
         0.40004,   0.39944,   0.39697,   0.39255,                            &
         0.38608,   0.37734,   0.36590,   0.35136,                            &
         0.33326,   0.31104,   0.28436,   0.25557,                            &
         0.22937,   0.21026,   0.19879,   0.19213,                            &
         0.18745,   0.18358,   0.18080/
    
    ! Which fluxes/cross sections to use 1=CMAT1, 2=SOLAR2000
    DATA FLUX_SOURCE/2/
    
    ! Plank constant*c, with a units conversion factor
    PCC=1.985E-25*1.e6   !*1.e9
    
    ! Ratio of nightime ionisation to sec=1. This is to represent ionisation by
    ! Galactic cosmic rays
    ! ald: changed from 1.e-3 to 1.e-6 in line with cmat1
    nightfac=1.e-6
    
    ! --------------------- Initialise begin -----------------------------
    
    IF(f107 /= last_f107) THEN

!!$#ifdef MESSY
!!$       IF (p_parallel_io) &
!!$#endif
!!$       write(6,*) "F107 changed, recalculating thermospheric fluxes ..."

       ! Get appropriate flux/cross-section data
       IF(FLUX_SOURCE == 2) THEN

          CALL S2K_FLUXES(f107, eccentric, CSAO, CSIO,                        &
               CSAO2, CSIO2, CSAN2, CSIN2, FLUX, WAVELS,                      &
               lyman_a_flux,                                                  &
               NWAVES_XRAY, NWAVES_SRC, NWAVES_EUV, lyman_a_num,              &
               iou1, iou2)
       ELSE

          CALL CMAT1_FLUXES(f107, eccentric, CSAO, CSIO,                      &
               CSAO2, CSIO2, CSAN2, CSIN2, FLUX, WAVELS,                      &
               lyman_a_flux,                                                  &
               NWAVES_XRAY, NWAVES_SRC, NWAVES_EUV, lyman_a_num)

       ENDIF

    ENDIF


    ! Initialise efficiencies
    IF(initialise == 1) THEN
       
!!$#ifdef MESSY
!!$       IF (p_parallel_io) &
!!$#endif
!!$       write(6,*) "Initialising Thermospheric heating efficiencies ..."

       
       ! Interpolate efficiencies onto current grid
       ! Points have to be ascending for spline to work
       DO n=1,TOPLVL
          pres63(n)=-pres63(n)
       ENDDO
       
       ! Fit spline for efficiencies on old grid
       CALL SPLINE1D(pres63,euveff63,TOPLVL,0._dp,0._dp,splft1,.TRUE.)
       CALL SPLINE1D(pres63,src1eff63,TOPLVL,0._dp,0._dp,splft2,.TRUE.)
       
       minn=-1
       maxn=-1
       
       DO n= ht_dim,1, -1
          
          eff1=1.
          eff2=1.
          
          ! Points have to be ascending for spline to work
          pval=-log(pres(n))
          
          IF((pres63(TOPLVL)-pval) > 0. .AND.                                 &
               pval > pres63(1)) THEN
             
             CALL SPLINT1D(pres63,euveff63,splft1,TOPLVL,pval,eff1,status)
#ifndef MESSY
                IF (status > 0) RETURN
#else
                IF (status > 0) THEN
                   write(*,*) 'SPLINT1D',m,l,status
                   RETURN
                ENDIF
#endif
             CALL SPLINT1D(pres63,src1eff63,splft2,TOPLVL,pval,eff2,status)
#ifndef MESSY
                IF (status > 0) RETURN
#else
                IF (status > 0) THEN
                   write(*,*) 'SPLINT1D',m,l,status
                   RETURN
                ENDIF
#endif
             
          ENDIF
          
          ! Find the points where we went out of range
          if(pval > pres63(TOPLVL)) maxn=n
          if(pval < pres63(1) .AND. minn == -1) minn=n
          
          euveff(n)=eff1
          src1eff(n)=eff2
          
       ENDDO
       
       ! Lets tidy up where we went out of bounds
       IF(maxn > 1) THEN
          DO n=maxn, ht_dim
             euveff(n)=euveff(n-1)
             src1eff(n)=src1eff(n-1)
          ENDDO
       ENDIF
       IF(minn > 1) THEN
          DO n=minn, 1, -1
             euveff(n)=euveff(n+1)
             src1eff(n)=src1eff(n+1)
          ENDDO
       ENDIF
       
       initialise = 0
       
       !RETURN  ! NO DO NOT RETURN, no init call exists mz_ab_20100910
       
    ENDIF
    ! --------------------- Initialise end -----------------------------


    ! Reset columns
    WO=0.
    WO2=0.
    WN2=0.
    WO3=0.

  
    ! Main height loop from top of atmosphere
    HEIGHT: DO n=ht_dim,2,-1
       
       ! number densities cm-3
       Oall(n)=gas1_1d(n)*1.e-6
       O2all(n)=gas2_1d(n)*1.e-6
       N2all(n)=gas3_1d(n)*1.e-6
       O3all(n)=gas_O3_1d(n)*nden1d(n)*1.e-6

       ! Set gravity
       gravy(n) = GRAV0
       
       ! Scale heights, in cm
       HO(n)=(GSCON*temp1d(n)/16./gravy(n))*1.E2
       HO2(n)=HO(n)/2.
       HO3(n)=HO(n)/3.
       HN2(n)=HO(n)/7.*4.
       
       ! Set total height in cm
       z(n)=(R0_eff(n))*1.E2
       
       PEUVO(n)  = 0.
       PEUVO2(n) = 0.
       PEUVN2(n) = 0.
       pre(n)    = 0.0
       paeuv(n)  = 0.0
       src1(n)   = 0.
       jlya(n)   = 0.
       jsrc(n)   = 0.
       
       ! calculate sec zenith angle, incorporating Chapmann grazing incidence
       ! function
       SECO=1./CHAPMANN(csza, HO(n), Z(n))
       SECO2=1./CHAPMANN(csza, HO2(n), Z(n))
       SECN2=1./CHAPMANN(csza, HN2(n), Z(n))
       SECO3=1./CHAPMANN(csza, HO3(n), Z(n))
        
       rnight_o=1.
       rnight_o2=1.
       rnight_n2=1.
       

       ! The stars are out ...
       IF(SECO < 0) THEN
          rnight_o=nightfac
          SECO=1.
       ENDIF
       IF(SECO2 < 0) THEN
          rnight_o2=nightfac
          SECO2=1.
       ENDIF
       IF(SECN2 < 0) THEN
          rnight_n2=nightfac
          SECN2=1.
       ENDIF
       IF(SECO3 < 0) THEN
          SECO3=1.
       ENDIF

       
       ! Integrate Column densities
       IF(n == ht_dim) THEN
          WO=Oall(n)*HO(n)*SECO
          WO2=O2all(n)*HO2(n)*SECO2
          WN2=N2all(n)*HN2(n)*SECN2
          WO3=O3all(n)*HO3(n)*SECO3
       ELSE
          ! distance between 2 levels in cm
          ZHT=(ht1d(n+1)-ht1d(n))*1.E2
          WO=WO+OALL(n)*ZHT*SECO
          WO2=WO2+O2ALL(n)*ZHT*SECO2
          WN2=WN2+N2ALL(n)*ZHT*SECN2
          WO3=WO3+O3ALL(n)*ZHT*SECO3
       ENDIF
       
       ! If not using composition, and subbing bond calculating airglow etc,
       ! use bulk efficiencies (eff = src1eff(n))
       eff=1.
     
       ! IF(SW(COMPO) /= 2) eff = src1eff(n)
       IF(chem_heat_on == 0) eff = src1eff(n)  
       
       ! Only bother if the sun's up
       IF(SECO >0 .OR. SECO2 >0 .OR. SECN2 > 0) THEN

          ! Loop over wavelength
      
          PEUVO2_loc  = 0.
          PEUVO_loc   = 0.
          PEUVN2_loc  = 0.
          pre_loc     = 0.
          PAEUV_loc   = 0.
          SRC1_loc    = 0.
          JSRC_loc    = 0.
          JLYA_loc    = 0.

          !$OMP PARALLEL PRIVATE(wav,TAUO,TAUO2,TAUN2,TAU,attenuation,        &
          !$OMP local_flux)                                                   &
          !$OMP SHARED(CSAO,CSAO2,CSAN2,WO,WO2,WN2,CSIO,CSIO2,CSIN2,N2all,    &
          !$OMP rnight_n2,rnight_o2,rnight_o,euveff,WAVELS,pcc,eff,           &
          !$OMP lyman_a_num)                                                  &
          !$OMP REDUCTION(+:PEUVO_loc,PEUVO2_loc,PEUVN2_loc,pre_loc,paeuv_loc,&
          !$OMP SRC1_loc,JSRC_loc)
          !$OMP DO
          WAVELENGTH1: DO wav=1, NWAVES_EUV+NWAVES_XRAY+NWAVES_SRC

             TAUO=(CSAO(wav))*WO
             TAUO2=(CSAO2(wav))*WO2
             TAUN2=(CSAN2(wav))*WN2
             TAU=(TAUO+TAUO2+TAUN2)
             
             attenuation=0.
             IF (tau <= 100. .AND. tau /= 0) attenuation = exp(-tau)
             
             local_flux=attenuation*FLUX(wav)
             IF(local_flux < 1 .OR. tau == 0.) local_flux=0.
             
             ! Calculate EUV ionization rate coefficients
             PEUVO_loc  = PEUVO_loc+local_flux*CSIO(wav)
             PEUVO2_loc = PEUVO2_loc+local_flux*CSIO2(wav)
             PEUVN2_loc = PEUVN2_loc+local_flux*CSIN2(wav)
             
             ! Calculate pre-dissociation of N2
             !pre_loc=pre_loc+(csan2(wav)-csin2(wav))*                         &
             !     local_flux*N2all(n)*rnight_n2

             ! ald; something funny here with pren2?
             pre_loc=pre_loc+(csan2(wav)-csin2(wav))*                         &
                  local_flux
             
             ! EUV heating calculation
             paeuv_loc=paeuv_loc+euveff(n)*(local_flux*(                      &
                  csio(wav)*oall(n)*rnight_o+                                 &
                  csio2(wav)*o2all(n)*rnight_o2+                              &
                  csin2(wav)*N2all(n)*rnight_n2)                              &
                  *(pcc))/WAVELS(wav)
             
             ! Schumann Runge Continuum, if the sun's up.
             IF(rnight_o2 == 1) THEN
                
                IF(wav > NWAVES_XRAY+NWAVES_EUV+1) THEN
                   
                   ! Calculate SRC heating
                   SRC1_loc = SRC1_loc + (eff*local_flux*                     &
                        (CSAO2(wav)*O2ALL(n))*pcc)/WAVELS(wav)
                   
                   ! Calculate JSRC, excluding JLYA
                   IF(wav /= lyman_a_num) JSRC_loc=JSRC_loc+                  &
                        (local_flux)*CSAO2(wav)

                   ! Calculate JLYA, which is at position 18
                   IF(wav == lyman_a_num)                                    &
                        JLYA_loc=local_flux*CSAO2(wav)
                   !IF(wav == lyman_a_num) JLYA_loc=FLUX(lyman_a_num)

                ENDIF

             ENDIF
             
          ENDDO WAVELENGTH1
          !$OMP END DO
          !$OMP END PARALLEL
          

       ELSE

       ! Nightime rates are either 0. or a fraction of the daytime rate

          PEUVO2_loc = 0.
          PEUVO_loc  = 0.
          PEUVN2_loc = 0.
          pre_loc    = 0.
          PAEUV_loc  = 0.
          SRC1_loc   = 0.
          JSRC_loc   = 0.
          JLYA_loc   = 0. 
          
       ENDIF

       PEUVO2(n) = PEUVO2_loc
       PEUVO(n)  = PEUVO_loc
       PEUVN2(n) = PEUVN2_loc
       pre(n)    = pre_loc
       PAEUV(n)  = PAEUV_loc
       SRC1(n)   = SRC1_loc
       JSRC(n)   = JSRC_loc
       JLYA(n)   = JLYA_loc

       ! Calculate EUV ionisation rates
       ! PEUVO(n)=PEUVO(n)*Oall(n)*rnight_o
       ! PEUVO2(n)=PEUVO2(n)*O2all(n)*rnight_o2
       ! PEUVN2(n)=PEUVN2(n)*N2all(n)*rnight_n2
       
       ! Populate photorate coeffcient arrays
       jx(n,m,ip_O3Pp,l)     = PEUVO_loc*rnight_o
       jx(n,m,ip_O2_b1b2,l)  = PEUVO2_loc*rnight_o2 
       jx(n,m,ip_N2,l)       = PEUVN2_loc*rnight_n2
       jx(n,m,ip_O3PO1D,l)   = (JSRC_loc+JLYA_loc)*rnight_o2

       ! NO + hv (Lya) = NO+ + e
       ! where 2.e-18 and 1.e-20 are the cross sections of NO and 
       ! O2 in cm-2. B&S.
       jx(n,m,ip_NOp,l) = LYMAN_A_FLUX*2.0e-18*        &
                                   exp(-1.0e-20*WO2)*rnight_n2

       ! Dissociation of NO : NO + hv(>150nm) = N(4S) + O
       ! Murrey et al 1994 say this should be 7.e-6
       jx(n,m,ip_NO,l) = 4.5e-6*exp(-1.e-8*(WO2)**.38) &
                                  *exp(-5.e-19*WO3)*rnight_n2

       ! Pre-dissociation N2
       ! jx(n,m,l,ip_NN2D)=pre(n)*1.e-30
       jx(n,m,ip_NN2D,l)= pre(n)*rnight_n2

       ! Update heating 
       ! ALD: If composition on, need to subtract bond energy.
       ! 5.12 is O2 bond energy, and 1.96 O1d energy.

       ! IF(SW(COMPO) == 2) THEN
       IF(chem_heat_on == 1) THEN 
          src1(n)=src1(n)-jx(n,m,ip_O3PO1D,l)*      &
                               ((5.12+1.96)*1.602e-19)*      &
                               gas2_1d(n)
       ENDIF

       heatrates(n,m, EUV_heat_l,l)=paeuv(n)/den1d(n)
       heatrates(n,m, UV_heat_l,l)=src1(n)/den1d(n)

       if(heatrates(n,m, UV_heat_l,l) .lt. 0.)                            &
          heatrates(n,m, UV_heat_l,l) = 0.

       
       !if(m.eq.45.and.l.eq.1.or.l.eq.10)write(6,*)                        &
       !  'SRC heating after ', n, heatrates(n,m,l, UV_heat_l)

    ENDDO HEIGHT   ! End height loop

    !* ALD 18/01/06: Call routine 'Secondaries'
    ! to get the ionisation rates due to secondary electrons. 
    ! seo_pr, seo2_pr, sen2_pr.

    CALL Secondaries(m,l,ht_dim,ht1d,pres)

    last_f107 = f107

    status = 0

    RETURN
  END SUBROUTINE ThermoSol
  
  !r ---------------------------------------------------------------------
  !r
  !r S2K_FLUXES    Routine to return Solar2000 fluxes and cross sections to
  !r               thermo_sol.
  !r
  !r  Inputs :
  !r
  !r  f107, eccentric
  !r
  !r  Outputs :
  !r
  !r  CSAO, CSIO,
  !r  CSAO2, CSIO2, CSAN2, CSIN2, FLUX, WAVELS, lyman_a_flux
  !r  lyman_a_num,
  !r  NWAVES_XRAY_OUT, NWAVES_SRC_OUT, NWAVES_EUV_OUT
  !r
  !r Alison's comments in all_Xsections.txt datafile -
  !r
  !r ALD 15/9/03:
  !r For wavelnegth bins 2-4,4-8,8-18,  the cross sections and fluxes
  !r are from Stan Solomon, based on Henke. They are not listed in this 
  !r table and are in the routine thermoheat.f as a data table.
  !r I only have flux values at 2 F107s and a simple extrapolation is done.
  !r
  !r For wavelengths 18-23, 23-32, 32-44, 33-60, 60-70 etc up to 1040-1050,
  !r cross sections are from Stan solomon based on Henke below 50A and 
  !r Fennelly and Torr above.For wavelengths 1050-1060 up to 1390-1400A, 
  !r values of O2 Absorption come from Rod Viereck thermsopheric heating 
  !r routine (in ../ROD/solflux_RAV.f) For wavelengths from 1400-1410 up 
  !r to 1790-1800, values are simple averages from BRUNOUT_O2.dat
  !r (file from Rod with references inside file), stored on apl machine in 
  !r /alison/Cross_sections/ROD/
  !r
  !r I have added lyman alpha as a seperate line at 1215.67A
  !r
  !r ALD: AUG 2003: DAT file of solar fluxes from Solar 2000 model (V24a) - 
  !r see below for explanation
  !r ALD: 20/8/03: This file contains lists of solar flux 
  !r (photons cm-2s-1nm-1) at specific wavelengths (Angstroms).
  !r Column 1 is wavelength bin - the first 5 bins are 18-23, 23-32, 32-44, 
  !r 44-60 Angstrom bins, the others are 10 Angstrom bins (i.e.60-70, 
  !r 70-80 etc). (These bins have been chosen as I have corresponding 
  !r photoabsorption and ionisation cross sections.)
  !r Columns 2 to 12 are photon fluxes at 1AU for the f107's listed.
  !r These numbers have been obtained by running the solar2000 model 
  !r avaiable at the following address:
  !r http://www.sel.noaa.gov/spacewx/index.html
  !r (Stored at UCL under /usr/users/alison/SOLAR2000/SOLAR2000). Results 
  !r from runs with full date info etc are in 
  !r /usr/users/alison/SOLAR2000/SOLAR2000/Results (listed by date), 
  !r or listed by F107,
  !r /usr/users/alison/imre-home/IDLprograms/CMAT_view/
  !r Files_to_view/Solar_Flux/s2k_out_*.
  !r
  !r
  !r
  
  SUBROUTINE S2K_FLUXES(f107, eccentric, CSAO, CSIO,                          &
       CSAO2, CSIO2, CSAN2, CSIN2, FLUX, WAVELS, lyman_a_flux,                &
       NWAVES_XRAY_OUT, NWAVES_SRC_OUT, NWAVES_EUV_OUT,                       &
       lyman_a_num, s2k_fluxes_f, all_Xsections_f)
    
    ! XRAYS are 1->3, EUV 4->103, SRC 104->183
    INTEGER :: NWAVES, NWAVES_EUV, NWAVES_SRC, NWAVES_XRAY
    PARAMETER(NWAVES=182, NWAVES_XRAY=3,                                      &
         NWAVES_SRC=76, NWAVES_EUV = 103)

    ! file handles
    INTEGER, INTENT(IN) :: s2k_fluxes_f, all_Xsections_f

    ! F107
    REAL(dp) :: f107
    
    ! Lymann-a flux
    REAL(dp) :: lyman_a_flux
    ! Lymann-a record number
    INTEGER :: lyman_a_num
    
    ! Orbital eccentricity factor
    REAL(dp) :: eccentric
    
    ! s2k flux and wavelengths
    REAL(dp) :: S2K_FLUX(NWAVES), WAVELS(NWAVES)
    
    ! FLUX, which is returned
    REAL(dp) :: FLUX(NWAVES)
    
    ! s2k cross sections
    REAL(dp) :: CSAO(NWAVES),CSIO(NWAVES), CSAN2(NWAVES),                   &
         CSIN2(NWAVES), CSAO2(NWAVES),                                        &
         CSIO2(NWAVES)
    
    ! X ray fluxes, cross-sections, and wavelength
    REAL(dp) :: XRFL(NWAVES_XRAY),XRFH(NWAVES_XRAY), XRF(NWAVES_XRAY),      &
         XCSAO(NWAVES_XRAY),XCSAO2(NWAVES_XRAY),XCSAN2(NWAVES_XRAY),          &
         XCSIO(NWAVES_XRAY),XCSIO2(NWAVES_XRAY),XCSIN2(NWAVES_XRAY),          &
         XWAVE(NWAVES_XRAY)
    
    ! Raw s2k fluxes and cross-sections in input file
    REAL(dp) :: S2K_ALLFLUXES(NWAVES-NWAVES_XRAY,12),                       &
         XSECTIONS(NWAVES-NWAVES_XRAY,7)
    
    ! Units conversion factor
    REAL(dp) :: unit_conv
    
    ! Some count variables
    INTEGER :: wav, countv, J
    
    ! Output versions, to be returned to heating code
    INTEGER :: NWAVES_EUV_OUT, NWAVES_SRC_OUT, NWAVES_XRAY_OUT
    
    ! ALD: Stan Solomons glow model Xray wavelength bands and fluxes :
    ! Wavelengths bands 2-4A, 4-8A, 8-18A.
    ! Photon flux in 1E9 cm-2 s-1
    DATA XWAVE/ 3., 6., 13.5/
    DATA XRFL/ 5.E-8, 1.E-5, 2.E-3/
    DATA XRFH/ 3.E-5, 8.E-4, 5.E-2/
    
    ! ALD: new xray cross sections from Stan Solomon (file sigma_1nm.dat)
    DATA XCSAO/ 0.001, 0.005, 0.038/
    DATA XCSAO2/ 0.001, 0.006, 0.046/
    DATA XCSAN2/ 0.001, 0.011, 0.076/
    DATA XCSIO/ 0.001, 0.005, 0.038/
    DATA XCSIO2/ 0.001, 0.006, 0.046/
    DATA XCSIN2/ 0.001, 0.011, 0.076/
    
!!$#ifdef MESSY
!!$    IF (p_parallel_io) &
!!$#endif
!!$    WRITE(6,*) "Initialising SOLAR2000 fluxes and cross-sections ..."
    
    unit_conv=1.E-18
    
    ! ALD SEPT2003:  Read the solar flux and cross section files
    ! at wavelength bins of 1 Angstrom. These come from Stan Solomon
    ! (Xsections), and the solar200 model.
    
    OPEN (s2k_fluxes_f,                                                       &
         FILE=TRIM(file_s2k_fluxes),STATUS='old')
    OPEN (all_Xsections_f,                                                    &
         FILE=TRIM(file_all_Xsections),STATUS='old')
    
    DO 2020 wav=1,NWAVES
       
       J=wav-NWAVES_XRAY
       
       IF(wav > NWAVES_XRAY) THEN
          
          READ(s2k_fluxes_f,*)                                                &
               (S2K_ALLFLUXES(J,countv),countv=1,12)
          READ(all_Xsections_f,*)                                             &
               (XSECTIONS(J,countv),countv=1,7)
          
          WAVELS(wav) = S2K_ALLFLUXES(J,1)
          CSAO(wav)   = XSECTIONS(J,2)
          CSIO(wav)   = XSECTIONS(J,3)
          CSAN2(wav)  = XSECTIONS(J,4)
          CSIN2(wav)  = XSECTIONS(J,5)
          CSAO2(wav)  = XSECTIONS(J,6)
          CSIO2(wav)  = XSECTIONS(J,7)
          
       ELSE
          
          ! Put xrays in
          CSAO(wav)   = XCSAO(wav)
          CSIO(wav)   = XCSIO(wav)
          CSAN2(wav)  = XCSAN2(wav)
          CSIN2(wav)  = XCSIN2(wav)
          CSAO2(wav)  = XCSAO2(wav)
          CSIO2(wav)  = XCSIO2(wav)
          WAVELS(wav) = XWAVE(wav)
          
          ! Get flux for current f107
          FLUX(wav)=(XRFH(wav) - XRFL(wav))*((F107-67.)/(243.-67.))           &
               + XRFL(wav) ! bn_ab_20110723 bug fix from CMAT v1.5
          IF(FLUX(wav) < 0.0) FLUX(wav)=0.
          FLUX(wav)=FLUX(wav)
          
       ENDIF
       
       
2020 ENDDO
    
    CLOSE(s2k_fluxes_f)
    CLOSE(all_Xsections_f)
    
    ! Now extrapolate the S2K fluxes for the specified F107.
    ! NB to use E10.7 instead of F107 you need to WRITE a new
    ! interpolation routine and check what the E107's are for
    ! the chosen dates.
    
    CALL S2K_INTERP(S2K_ALLFLUXES, F107, S2K_FLUX)
    
    ! Now put these in flux array
    DO 30 wav=NWAVES_XRAY+1,NWAVES
       J=wav-NWAVES_XRAY
       FLUX(wav) = S2K_FLUX(J)
30  ENDDO
    
    
    ! Convert to appropriate units
    DO 35 wav=1,NWAVES
       
       CSAO(wav)=CSAO(wav)*unit_conv
       CSAO2(wav)=CSAO2(wav)*unit_conv
       CSAN2(wav)=CSAN2(wav)*unit_conv
       CSIO(wav)=CSIO(wav)*unit_conv
       CSIO2(wav)=CSIO2(wav)*unit_conv
       CSIN2(wav)=CSIN2(wav)*unit_conv
       WAVELS(wav)=WAVELS(wav)*1.e-10
       FLUX(wav)=FLUX(wav)*eccentric
       
       
       ! WRITE(6,"(4E10.2, A)") WAVELS(wav), CSIO(wav),
       ! &        CSIO2(wav),CSIN2(wav), "s2k"
       ! WRITE(6,"(4E10.2,A)") WAVELS(wav), CSAO(wav),
       ! &        CSAO2(wav),CSAN2(wav) , "s2k"
       ! WRITE(6,"(2E10.2,A)") WAVELS(wav),
       ! &        FLUX(wav), "s2k"
       
35  ENDDO
    
    ! set Lya flux used in chemistry. We should really
    ! indentify this by wavelength.
    lyman_a_num = NWAVES_XRAY + NWAVES_EUV + 18
    lyman_a_flux = FLUX(lyman_a_num)
    
    NWAVES_XRAY_OUT=NWAVES_XRAY
    NWAVES_SRC_OUT=NWAVES_SRC
    NWAVES_EUV_OUT=NWAVES_EUV
    
    RETURN
  END SUBROUTINE S2K_FLUXES
  
  
  !r  -------------------------------------------------------------------
  !r
  !r  CMAT1_FLUXES   routine to return CMAT1 fluxes. It's a bit rough and
  !r                 ready this routine, only for regression testing. New
  !r                 Solar2000 (s2k) data is preferable.
  !r
  !r  Inputs :
  !r
  !r  f107, eccentric
  !r
  !r  Outputs :
  !r
  !r  CSAO, CSIO,
  !r  CSAO2, CSIO2, CSAN2, CSIN2, FLUX, WAVELS, lyman_a_flux
  !r  lyman_a_num
  !r  NWAVES_XRAY, NWAVES_SRC, NWAVES_EUV)
  !r
  
  SUBROUTINE CMAT1_FLUXES(f107, eccentric, CSAO, CSIO,                        &
       CSAO2, CSIO2, CSAN2, CSIN2, FLUX, WAVELS, lyman_a_flux,                &
       NWAVES_XRAY_OUT, NWAVES_SRC_OUT, NWAVES_EUV_OUT, lyman_a_num)

    ! XRAYS are 1->7, EUV 8->65, SRC 66->91
    INTEGER :: NWAVES, NWAVES_EUV, NWAVES_SRC, NWAVES_XRAY
    PARAMETER(NWAVES=91, NWAVES_XRAY=7,                                       &
         NWAVES_SRC=27, NWAVES_EUV = 57)
    
    ! Output versions, to be returned to heating code
    INTEGER :: NWAVES_EUV_OUT, NWAVES_SRC_OUT, NWAVES_XRAY_OUT
    
    ! FLUX, which is returned
    REAL(dp) :: FLUX(NWAVES)
    
    ! cross sections
    REAL(dp) :: CSAO(NWAVES),CSIO(NWAVES), CSAN2(NWAVES),                   &
         CSIN2(NWAVES), CSAO2(NWAVES),                                        &
         CSIO2(NWAVES), WAVELS(NWAVES)
    
    ! Lyman-a flux
    REAL(dp) :: lyman_a_flux
    ! Lyman-a record number
    INTEGER :: lyman_a_num
    
    REAL(dp) :: f107, FMXFMN, eccentric
    INTEGER :: wav, j
    
    REAL(dp) :: SFL(NWAVES_EUV),CSAO_C1(NWAVES_EUV),CSAO2_C1(NWAVES_EUV),   &
         CSAN2_C1(NWAVES_EUV),CSIO_C1(NWAVES_EUV),                            &
         CSIO2_C1(NWAVES_EUV),CSIN2_C1(NWAVES_EUV),                           &
         SFH(NWAVES_EUV),SF(NWAVES_EUV),RLAM(65)
    
    REAL(dp) :: xrfl(NWAVES_XRAY),xrfh(NWAVES_XRAY),xrf(NWAVES_XRAY),       &
         xcsao(NWAVES_XRAY),xcsao2(NWAVES_XRAY),                              &
         xcsan2(NWAVES_XRAY),xcsio(NWAVES_XRAY),xcsio2(NWAVES_XRAY),          &
         xcsin2(NWAVES_XRAY)
    
    REAL(dp) :: O2CROSS(NWAVES_SRC),FLUX_C1(NWAVES_SRC),                    &
         LAMBDA(NWAVES_SRC),TMP1(NWAVES_SRC)
    REAL(dp) :: SOLVAR(NWAVES_SRC)
    
    
    ! EUV and X-ray cross-sections and fluxes -- START -------
    
    ! Bottom 20 are soft-xrays,up until last 8 EUV, last 8 are O2 SRC
    DATA RLAM/18.6,19.0,21.6,21.8,22.1,28.5,28.8,29.5,30.0,                   &
         30.4,33.7,41.0,43.8,44.0,44.2,45.7,46.4,46.7,47.9,49.2,              &
         75.,125.,175.,225.,256.3,284.15,275.,303.31,303.78,   & ! <-Torr
         325.,368.07,375.,425.,465.22,475.,525.,554.37,584.33,                &
         575.,609.76,629.73,625.,675.,730.36,725.,765.15,770.41,              &
         789.36,775.,825.,875.,925.,977.62,975.,1025.72,1031.91,              &
         1025.,1387.5,1425.,1475.,1525.,1575.,1625.,1675.,1725./
    
    DATA SFL/                                                                 &
         .0001,.0001,.0003,.0001,.0003,.0005,.0025,.0022,.0012,               &
         .0006,.0011,.0006,.0021,.0008,.0009,.0005,.0027,.0052,               &
         .0059,.0043,                                                         &
         .38,.13,1.84,.92,.27,.1,.84,.24,6.,.87,.74,.21,.39,.18,              &
         .31,.51,.80,1.58,.48,.45,1.5,.17,.22,.39,.17,.2,.24,                 &
         .79,.87,1.93,4.43,4.22,5.96,1.79,4.38,3.18,3.64/
    
    DATA SFH/                                                                 &
         .0016,.0053,.0048,.0016,.0048,.0072,.0211,.0186,.0024,               &
         .0104,.0158,.0073,.0130,.0097,.0109,.0061,.0168,.0107,               &
         .0121,.0267,                                                         &
         1.37,.468,5.7,7.14,1.08,5.72,12.16,4.69,14.39,6.83,1.53,             &
         2.54,1.53,.736,1.82,1.64,1.52,4.3,1.048,2.48,3.87,1.37,              &
         .539,.746,.429,.439,1.19,1.514,2.454,4.85,12.219,9.85,               &
         10.217,4.078,11.85,6.1,6.09/


    DATA CSAO_C1/                                                             &
         .34,.36,.5,.51,.52,.05,.05,.06,.06,.06,.08,.13,.15,                  &
         .15,.16,.17,.18,.18,.19,.21,                                         &
         .06,3.53,5.96,7.55,8.43,9.26,8.78,9.7,9.72,10.03,10.84,              &
         10.7,11.21,11.25,11.64,11.91,12.13,12.17,11.9,12.23,12.22,           &
         12.21,10.04,11.35,8.0,4.18,4.18,4.28,4.23,4.38,4.18,2.12,            &
         0.,0.,0.,0.,0./
    
    DATA CSAO2_C1/                                                            &
         .69,.72,.99,1.01,1.05,.10,.11,.11,.12,.12,.16,.26,.3,.31,            &
         .31,.34,.35,.36,.38,.41,                                             &
         1.18,3.61,7.27,10.5,12.8,14.8,13.65,15.98,16.,17.19,18.40,           &
         18.17,19.39,20.4,21.59,24.06,25.59,22.0,25.04,26.1,25.8,             &
         26.02,26.27,25.,29.05,21.98,25.18,26.66,27.09,20.87,9.85,            &
         15.54,4.0,16.53,1.6,1.0,1.1/
    
    DATA CSAN2_C1/                                                            &
         .44,.47,.65,.67,.69,1.13,1.13,1.12,1.11,1.10,.1,.16,.19,             &
         .19,.19,.21,.22,.22,.24,.25,                                         &
         .6,2.32,5.4,8.15,9.65,10.6,10.8,11.58,11.6,                          &
         14.6,18.0,17.51,21.07,21.8,                                          &
         21.85,24.53,24.69,23.2,22.38,23.1,23.2,23.22,29.75,26.3,             &
         30.94,35.46,26.88,19.26,30.71,15.05,46.63,16.99,.7,                  &
         36.16,0.,0.,0./
    
    DATA CSIO_C1/                                                             &
         .34,.36,.5,.51,.52,.05,.05,.06,.06,.06,.08,.13,.15,                  &
         .15,.16,.17,.18,.18,.19,.21,                                         &
         1.06,3.53,5.96,7.55,8.43,9.26,8.78,9.7,9.72,10.03,10.84,             &
         10.7,11.21,11.25,11.64,11.91,12.13,12.17,11.9,12.23,12.22,           &
         12.21,10.04,11.35,8.0,4.18,4.18,4.28,4.23,4.38,4.18,2.12,            &
         0.,0.,0.,0.,0./
    
    DATA CSIO2_C1/                                                            &
         .69,.72,.99,1.01,1.05,.10,.11,.11,.12,.12,.16,.26,.3,.31,            &
         .31,.34,.35,.36,.38,.41,                                             &
         1.18,3.61,7.27,10.5,12.8,14.8,13.65,15.98,16.,17.19,18.40,           &
         18.17,19.39,20.4,21.59,24.06,25.59,22.0,25.04,26.1,25.8,             &
         25.94,22.05,23.0,23.81,8.59,9.69,11.05,9.39,6.12,4.69,9.34,          &
         2.5,12.22,1.,0.,.27/

    DATA CSIN2_C1/                                                            &
         .44,.47,.65,.67,.69,1.13,1.13,1.12,1.11,1.10,.1,.16,.19,& ! Soft X-ray
         .19,.19,.21,.22,.22,.24,.25,                            & ! Soft X-ray
         .6,2.32,5.4,8.15,9.65,10.6,10.08,11.58,11.6,            & ! EUV
         14.6,18.0,17.51,21.07,21.8,                             & ! EUV
         21.85,24.53,24.69,23.2,22.38,23.1,23.2,23.22,25.06,23., & ! EUV
         23.2,23.77,18.39,10.18,16.75,0.,0.,0.,0.,0.,0.,0.,0./    ! EUV+SRC pad
    

    !c    xray fluxes based on the SOLRAD data analysis performed by
    !c    Frank Eparvier, May 94.
    data xrfl/.000000026,.00000027,.0000015,.00004,.00021,.00037,            &
         .00063/
    data xrfh/.000028,.000068,.00038,.0086,.025,.074,.20/
    
    data xcsao/.002,.009,.03,.09,.14,.2,.25/
    data xcsao2/.004,.018,.05,.18,.27,.35,.45/
    data xcsan2/.0025,.01,.03,.12,.17,.25,.35/
    data xcsio/.002,.009,.03,.09,.14,.2,.25/
    data xcsio2/.004,.018,.05,.18,.27,.35,.45/
    data xcsin2/.0025,.01,.03,.12,.17,.25,.35/
    
    ! EUV and X-ray cross-sections and fluxes -- END   -------
    
    ! SRC cross-sections and fluxes ------------ START -------
    
    DATA LAMBDA/ 1075E-8,1125E-8,1175E-8,1215.67E-8,1225E-8,1275E-8,          &
         1302.17E-8,1304.86E-8,1306.03E-8,1334.53E-8,                         &
         1335.70E-8,1325E-8,1393.76E-8,                                       &
         1375E-8,1402.77E-8,1425E-8,                                          &
         1475E-8,1548.19E-8,1525E-8,1550.77E-8,1561.0E-8,                     &
         1575E-8,1625E-8,1657.2E-8,1675E-8,1725E-8,                           &
         1775E-8/
    
    DATA FLUX_C1/2.9E9,0.091E9,4.4E9,251.0E9,8.0E9,4.1E9,1.1E9,1.13E9,        &
         1.23E9,1.8E9,2.5E9,4.65E9,1.3E9,6.1E9,0.91E9,9.5E9,                  &
         16.2E9,3.8E9,25.2E9,1.9E9,2.5E9,35.6E9,56.0E9,8.5E9,                 &
         121.5E9,225.0E9,357.0E9/
    
    ! Solar variance of above fluxes
    DATA SOLVAR/1.,1.,1.,1.8,1.1,1.3,0.57,0.57,0.57,1.2,1.2,                  &
         1.2,1.2,1.2,1.2,1.1,1.1,0.8,0.9,0.8,0.8,                             &
         0.7,0.7,0.6,0.7,0.6,0.3/
    
    ! Mass ab cross section 8e-20
    DATA O2CROSS/ 7.387E-19,4.755E-19,1.673E-18,1E-20,6.304E-18,              &
         9.778E-19,4.5E-19,4.2E-19,4.0E-19,2.082E-18,                         &
         1.8E-18,1.8E-18,1.419E-17,1.36E-17,1.4E-17,1.50E-17,                 &
         1.239E-17,9.833E-18,1.06E-17,9.2E-18,8.1E-18,                        &
         4.551E-18,3.468E-18,2.0E-18,1.600E-18,                               &
         3.149E-19,1.051E-19/
    
    ! SRC cross-sections and fluxes ------------ END    -------
    
!!$#ifdef MESSY
!!$    IF (p_parallel_io) &
!!$#endif
!!$    WRITE(6,*) "Initialising CMAT1 fluxes and cross-sections ...."
    
    
    FMXFMN = (F107-71.)/(220.-71.)     ! Used for SRC
    
    DO 964 j=1,NWAVES_XRAY    ! Hard X-rays
       xrFL(J)=xrFL(J)*1.E9
       xrFH(J)=xrFH(J)*1.E9
       xCSAO(J)=xCSAO(J)*1.E-18
       xCSAO2(J)=xCSAO2(J)*1.E-18
       xCSAN2(J)=xCSAN2(J)*1.E-18
       xCSIO(J)=xCSIO(J)*1.E-18
       xCSIO2(J)=xCSIO2(J)*1.E-18
       xCSIN2(J)=xCSIN2(J)*1.E-18
964 ENDDO
    
    DO 1 J=1,NWAVES_EUV    ! Soft X-rays(1-20), and EUV(20-37)
       SFL(J)=SFL(J)*1.E9
       SFH(J)=SFH(J)*1.E9
       CSAO_C1(J)=CSAO_C1(J)*1.E-18
       CSAO2_C1(J)=CSAO2_C1(J)*1.E-18
       CSAN2_C1(J)=CSAN2_C1(J)*1.E-18
       CSIO_C1(J)=CSIO_C1(J)*1.E-18
       CSIO2_C1(J)=CSIO2_C1(J)*1.E-18
       CSIN2_C1(J)=CSIN2_C1(J)*1.E-18
       
       
       IF(j <= NWAVES_SRC) THEN
          TMP1(j) = FLUX_C1(j)*ECCENTRIC*(1.+SOLVAR(j)*                       &
               FMXFMN)
       ENDIF
       
1   ENDDO
    
    ! set Lya flux used in chemistry. We should really
    ! indentify this by wavelength.
    lyman_a_num = NWAVES_XRAY + NWAVES_EUV + 4
    lyman_a_flux = FLUX(lyman_a_num)
    
    DO 1922 j=1,20     ! Soft X-ray Solar flux
       
       sf(j)=sfl(j)*(0.82*f107-49.86)/1.5
       !c
       !c       normalize soft X-ray flux to measured values at F10.7=127
       !c       observed on day 79 1998.  SNOE observed 7.83x10(7),
       !c       reference spectrum is 0.32x10(8) which is scaled by 35.09
       !c       correction factor = 0.78/(35.09x0.32) = 0.07
       sf(j)=sfl(j)*(0.82*f107-49.86)*0.07/1.5
       
1922 ENDDO
    
    
    
    DO 7 J=21,NWAVES_EUV              ! EUV Solar flux Flux
       ! SF(J)=(SFH(J)-SFL(J))*F107/172.-0.413*SFH(J)+1.413*SFL(J)
       ! SF(J)=(SFH(J)-SFL(J))*F107/176. -0.38*SFH(J)+1.38*SFL(J)
       SF(J)=(SFH(J)-SFL(J))*(F107-67)/176. + SFL(J)
       if(sf(j) <= 0.0)sf(j)=sfl(j)
7   ENDDO
    
    DO 966 J=1,NWAVES_XRAY              ! Hard x-ray flux
       xrF(J)=(xrFH(J)-xrFL(J))*F107/176.-0.38*xrFH(J)+1.38*xrFL(J)
       if(xrf(j) <= 0.0)xrf(j)=xrfl(j)
966 ENDDO
    
    
    ! Populate data for return
    DO 35 wav=1,NWAVES
       
       ! X rays
       IF(wav <= NWAVES_XRAY) THEN
          
          CSAO(wav)=xCSAO(wav)
          CSAO2(wav)=xCSAO2(wav)
          CSAN2(wav)=xCSAN2(wav)
          CSIO(wav)=xCSIO(wav)
          CSIO2(wav)=xCSIO2(wav)
          CSIN2(wav)=xCSIN2(wav)
          WAVELS(wav)=RLAM(wav)*1.e-10
          FLUX(wav)=xrF(wav)
          
       ENDIF
       
       ! EUV
       IF(wav > NWAVES_XRAY .AND. wav <= (NWAVES_EUV+NWAVES_XRAY)) THEN
          
          CSAO(wav)=CSAO_C1(wav)
          CSAO2(wav)=CSAO2_C1(wav)
          CSAN2(wav)=CSAN2_C1(wav)
          CSIO(wav)=CSIO_C1(wav)
          CSIO2(wav)=CSIO2_C1(wav)
          CSIN2(wav)=CSIN2_C1(wav)
          WAVELS(wav)=RLAM(wav)*1.e-10
          FLUX(wav)=SF(wav)
          
       ENDIF
       
       ! SRC
       IF(wav > (NWAVES_EUV+NWAVES_XRAY)) THEN
          
          CSAO(wav)=0.
          CSAN2(wav)=0.
          CSIO(wav)=0.
          CSIO2(wav)=0.
          CSIN2(wav)=0.
          
          CSAO2(wav)=O2CROSS(wav-(NWAVES_EUV+NWAVES_XRAY))
          WAVELS(wav)=LAMBDA(wav-(NWAVES_EUV+NWAVES_XRAY))*1.e-2
          FLUX(wav)=TMP1(wav-(NWAVES_EUV+NWAVES_XRAY))
          
          
       ENDIF
       
       ! WRITE(6,"(4E10.2,A)") WAVELS(wav), CSIO(wav),
       ! &        CSIO2(wav),CSIN2(wav), "cmat"
       ! WRITE(6,"(4E10.2,A)") WAVELS(wav), CSAO(wav),
       ! &        CSAO2(wav),CSAN2(wav), "cmat"
       ! WRITE(6,"(2E10.2,A)") WAVELS(wav),
       ! &        FLUX(wav), "cmat"
       
35  ENDDO
    
    NWAVES_XRAY_OUT=NWAVES_XRAY
    NWAVES_SRC_OUT=NWAVES_SRC
    NWAVES_EUV_OUT=NWAVES_EUV
    
    RETURN
    
  END SUBROUTINE CMAT1_FLUXES
  
  
  SUBROUTINE S2K_INTERP(S2K_ALLFLUXES, F107_in, S2K_FLUX)
    
    IMPLICIT NONE
    
    ! --------------------------------------------------------------------
    ! ALD: SEPT 2003:
    ! Routine to interpolate solar flux values onto chosen F107.
    ! Input:  array of fluxes from solar 2000 model (s2k_allfluxes(179,12))
    ! Units: photons (photons cm-2s-1). See /input_params/s2k_fluxes
    ! for explanation.
    ! Input:  F107
    ! Output: s2k_flux(179) - an array of photon fluxes for 179
    ! wavelength bins, linearly interpolated to the correct F107
    ! by looking at the photon fluxes at F107 values above and below
    ! Other:  S2K_F107s - data table of the F107s at which S2K model was run.
    ! SFH, SFL: Interpolation values of high and low flux
    ! F107H, F107L: as above
    ! NWAVES: No of wavelength bins
    ! NF107s: No of F107s in S2K_ALLFLUXES data file
    ! ----------------------------------------------------------------------
    
    INTEGER :: J,K,L, NWAVES, NF107s
    PARAMETER (NWAVES = 179, NF107s = 11)

    REAL(dp) :: f107_in 
    REAL(dp) :: S2K_ALLFLUXES(NWAVES,NF107s+1), F107, S2K_FLUX(NWAVES)
    REAL(dp) :: S2K_F107s(NF107s)
    REAL(dp) :: SFH(NWAVES), SFL(NWAVES), F107H, F107L

    DATA  S2K_F107s/  67., 85., 102., 120., 137., 155.,                       &
         172., 190., 208., 225., 243/
    

    f107 = f107_in

    IF (F107 < 67.) F107 = 67.
    IF (F107 > 243.) F107 = 243.
    
    DO 10 J=1,(NF107s-1)
       IF (F107 >= S2K_F107s(J)) THEN
          IF (F107 == 243.) THEN
             DO K=1,NWAVES
                S2K_FLUX(K)=S2K_ALLFLUXES(K,NF107s+1)
                ! WRITE(6,*)'flux',J,K,S2K_FLUX(K)
             ENDDO
             GOTO 20
          ELSE
             IF  (F107 < S2K_F107s(J+1)) THEN
                DO 30 L=1,NWAVES
                   
                   ! Don't forget S2K_ALFLUXES has a colmn for wavelengths so 
                   ! you need to add an extra 1 to the column number
                   
                   SFH(L) = S2K_ALLFLUXES(L,J+2)
                   SFL(L) = S2K_ALLFLUXES(L,J+1)
                   ! WRITE(6,*)'SFH, SFL', SFH(L), SFL(L)
                   F107H = S2K_F107s(J+1)
                   F107L = S2K_F107s(J)
                   
                   S2K_FLUX(L)=(((F107-F107L)/(F107H-F107L))                  &
                        *(SFH(L)-SFL(L)))+SFL(L)
                   ! WRITE(6,*)'Flux', S2K_FLUX(L)
                   
30              ENDDO
                
             ENDIF
          ENDIF
       ENDIF
10  ENDDO
20  CONTINUE
    
  END SUBROUTINE S2K_INTERP

  ! Not implemented yet.         
  SUBROUTINE H2O_IR(T,HT,Z)


    !**********************************************************************
    !*
    !*                 SUBROUTINE MODH2O
    !*
    !**********************************************************************
    !
    ! Routine to derive simple H2O thermal infrared heating rates
    ! - uses coefficients calculated for Met.Office 32-level grid
    ! - hence must interpolate on and off this old grid
    !
    ! Called by MODIR
    ! Calls     nothing
    !

    INTEGER, PARAMETER   :: NLONG=72+2,NLAT=36,NHGHT=30,                      &
         NP1=NHGHT+1, NP2=NP1+1, NP3=NP1+2, NP4=NP1+3

    !----------------------------------------------------------------------
    
    ! input temperatures and vertical grid
    
    REAL T(NLONG,NLAT,NHGHT), Z(NHGHT)
    
    ! Output heating rates
    
    REAL HT(NLONG,NLAT,NHGHT)
    
    ! internal arrays, new grid
    
    REAL DZD2(NHGHT), TT
    
    INTEGER ILEV2(2,NHGHT)
    
    ! internal array, old grid
    
    REAL ZOLD(32), HTOL(32), DZD(32)
    
    INTEGER ILEV(2,32)
    
    ! Constants in quadratic fit of heating rate to local temperature
    
    REAL F1(32), F2(32), F3(32)
    
    ! Loop variables
    
    INTEGER LNG, LAT, LHT
    
    ! Other integer variables
    
    INTEGER IFL, LOLD, LNEW, JJ
    
    DATA IFL/0/
    
    ! Constants in quadratic fit at 32 levels
    
    DATA F1/-0.00636,-0.007687,-0.0198,-0.01476,-0.0117,-0.01080, &
         &-0.01004,-0.01030,-0.009887,-0.01038,-0.01005,-0.009818,-0.01004, &
         &-0.01035,-0.01077,-0.01114,-0.01163,-0.01147,-0.01167, &
         &-0.01121,-0.01147,-0.01106,-0.01125,-0.01067,-0.01048,-0.01065, &
         &-0.01004,-0.009944,-0.00933,-0.008312,-0.007283,-0.007146/
    
    DATA F2/4.094E-5,4.73E-5,1.197E-4,8.544E-5,6.5E-5,6.002E-5, &
         &5.357E-5,5.592E-5,5.196E-5,5.535E-5,5.168E-5,4.867E-5,4.881E-5, &
         &4.923E-5,4.965E-5,5.052E-5,5.354E-5,5.372E-5,5.772E-5,5.627E-5, &
         &6.02E-5,5.871E-5,6.246E-5,6.007E-5, &
         &6.003E-5,6.196E-5,5.761E-5,5.83E-5,5.449E-5,4.767E-5,4.072E-5, &
         &4.101E-5/
    
    DATA F3/-7.933E-8,-8.841E-8,-2.218E-7,-1.538E-7,-1.134E-7, &
         &-1.048E-7,-9.055E-8,-9.581E-8,-8.699E-8,-9.336E-8,-8.508E-8, &
         &-7.831E-8,-7.753E-8,-7.722E-8,-7.63E-8,-7.706E-8,-8.236E-8, &
         &-8.352E-8,-9.257E-8,-9.116E-8,-1.002E-7, &
         &-9.865E-8,-1.08E-7,-1.051E-7,-1.064E-7,-1.109E-7,-1.023E-7, &
         &-1.052E-7,-9.803E-8,-8.503E-8,-7.152E-8,-7.356E-8/
    


    ! First interpolate from CMAT2 grid to Met office new grid

    ! This first section is called only on the first entry: It
    ! sets up the old 32-level grid and calculates the interpolation
    !
    IF(IFL.EQ.0)THEN

       WRITE(6,1000)
1000   FORMAT(' SUBROUTINE MODH2O VERSION 3.2 : 26 APR 1988 ')
       IFL=1

       DO 10 LHT=1,32
          ZOLD(LHT)=2.608+(LHT-1)*0.2878
10     CONTINUE


       ! first calculate interpolation constants to get from temperatures
       ! on new grid to 32 level grid

        LOLD=1
       DO 20 LHT=1,NHGHT-1
25     CONTINUE
       IF(ZOLD(LOLD).LT.Z(1))THEN
          DZD(LOLD)=0.0
          ILEV(1,LOLD)=1
          ILEV(2,LOLD)=1
          LOLD=LOLD+1
          GOTO 25
       ENDIF

       IF(ZOLD(LOLD).GE.Z(LHT).AND.ZOLD(LOLD).LE.Z(LHT+1))THEN
          ILEV(1,LOLD)=LHT
          ILEV(2,LOLD)=LHT+1
          DZD(LOLD)=(ZOLD(LOLD)-Z(LHT))/(Z(LHT+1)-Z(LHT))
          LOLD=LOLD+1
          IF(LOLD.GT.32)GOTO 30
          GOTO 25
       ENDIF
20     CONTINUE
30     CONTINUE

       IF(LOLD.LE.32)THEN
          DO 40 LHT=LOLD,32
             ILEV(1,LHT)=NHGHT
             ILEV(2,LHT)=NHGHT
             DZD(LHT)=0.0
40        CONTINUE
       ENDIF

       ! and then calculate interpolation constants to go from old grid
       ! heating rates to new grid heating rates

       LNEW=1
       DO 50 LHT=1,31
55        CONTINUE
          IF(Z(LNEW).LT.ZOLD(1))THEN
             DZD2(LNEW)=0.0
             ILEV2(1,LNEW)=1
             ILEV2(2,LNEW)=1
             LNEW=LNEW+1
             GOTO 55
          ENDIF
          IF(Z(LNEW).GE.ZOLD(LHT).AND.Z(LNEW).LE.ZOLD(LHT+1))THEN
             ILEV2(1,LNEW)=LHT
             ILEV2(2,LNEW)=LHT+1
             DZD2(LNEW)=(Z(LNEW)-ZOLD(LHT))/0.2878
             LNEW=LNEW+1
             IF(LNEW.GT.NHGHT)GOTO 60
             GOTO 55
          ENDIF
50     CONTINUE
60     CONTINUE

       IF(LNEW.LE.NHGHT)THEN
          DO 70 LHT=LNEW,NHGHT
             ILEV2(1,LHT)=32
             ILEV2(2,LHT)=32
             DZD2(LHT)=0.0
70        CONTINUE
       ENDIF

    ENDIF

    ! Main calculation of heating rates - interpolate temperatures
    ! onto old model grid to calculate heating rates

      DO 801 LAT=1,NLAT
        DO 802 LNG=1,NLONG
          DO 90 LHT=1,32
            JJ=33-LHT
            TT=(1.0-DZD(LHT))*T(LNG,LAT,ILEV(1,LHT)) &
     &        +DZD(LHT)*T(LNG,LAT,ILEV(2,LHT))

!  To convet to K/s
!   6.56226E-13 = 5.6698E-8/86400.0

            HTOL(LHT)= &
     &        ((F3(JJ)*TT+F2(JJ))*TT+F1(JJ))*6.56226E-13*TT**4
   90     CONTINUE

! and then interpolate heating rates onto new grid

          DO 100 LHT=1,NHGHT
             HT(LNG,LAT,LHT)=(1.0-DZD2(LHT))*HTOL(ILEV2(1,LHT)) &
     &         +DZD2(LHT)*HTOL(ILEV2(2,LHT))
  100     CONTINUE
   802 CONTINUE
   801 CONTINUE


  END SUBROUTINE H2O_IR


  !========================================================
  !=                      Secondaries                     =
  !========================================================
  !r Routine to get the ionisation rate due to secondary 
  !r electrons 
  !r
  !r Inputs:
  !r m and l   : latitude and longitude grid counters
  !r ht1d      : 1d array of heights at m and l
  !r jx        : array of photoionisation rates
  !r
  !r Outputs:
  !r jx        : array of photoionisation rates including
  !r             secondary electron ionisation
  
  SUBROUTINE Secondaries(m,l,ht_dim,ht1d,pres)
    
    INTEGER,    INTENT(IN)  :: m, l, ht_dim
    INTEGER                 :: i, n
    
    REAL(dp), INTENT(IN)  :: ht1d(ht_dim), pres(ht_dim)
    
    ! PEOP, PEO2, PEN2 : Ratio's of Photoelectron ionisation to photon
    !  ionisation for each major
    
    ! PEOP,PEN2,PEO2 defined globally for saving

    REAL(dp) :: PEO2POL(11), PEN2POL(11), PEOPPOL(11)
    REAL(dp) :: EXEX
    
    ! variables for smoothing
    ! REAL(dp) :: sm_dim1=7., sm_dim2=3.
    
    DATA PEOPPOL/5.6585493, -0.21288884, 0.26014609,                        &
         0.12820983, -0.20868770, 0.072935121,                      &
         -0.012377555, 0.0011888792, -6.6231050e-05,                &
         2.0024786e-06, -2.5484502e-08/
    
    DATA PEO2POL/0.0094207625, -0.036520917, 0.039762641,                   &
         0.015361183, -0.033185664, 0.016094355,                    &
         -0.0037619671, 0.00048778766, -3.5949878e-05,              &
         1.4130961e-06, -2.3047003e-08/
    
    DATA PEN2POL/51.365471, -16.681270, 1.3220374,                          &
         -4.1530864, 3.6937091, -1.3627621,                         &
         0.27327087, -0.032298698, 0.0022525608,                    &
         -8.5878012e-05,1.3813762e-06/
    
    INTEGER :: initialise = 1
    
    ! Set up Secondary electron profiles
    ! See Fuller-Rowell, modelling changes in NO in the Thermosphere
    ! and upper Mesosphere JGR, 1992
    
    ! Only need to set up photo / photoelectron ionisation ratios
    ! once
    
    ! IF(m == 45 .and. l==1) then
    !  DO N=1,HT_DIM
    !    write(6,"(a30,i3,5e12.4)") "jx 1, N2, O2, O2 ht, P ", n,    &
    !                       ht1d(n), pres(n),                                &
    !                       jx(n,m,l,ip_N2P),                       &
    !                       jx(n,m,l,ip_O2P),                       &
    !                       jx(n,m,l,ip_O3Pp)
    !  ENDDO
    ! ENDIF
    
    IF(initialise == 1) THEN
       
       DO N=1,HT_DIM
          
          peop(n)=0.
          peo2(n)=0.
          pen2(n)=0.
          if(pres(n) <= 1.0376) then
             DO i=1,11
                EXEX=+1.*log(1.0376)-1.*log(pres(n))
                PEOP(n)=PEOP(n)+(PEOPPOL(i)*(exex**(i-1.)))
                PEO2(n)=PEO2(n)+(PEO2POL(i)*(exex**(i-1.)))
                PEN2(n)=PEN2(n)+(PEN2POL(i)*(exex**(i-1.)))
             ENDDO
             if(exex > 12.96) then
                PEOP(n)=PEOP(n-1)
                PEO2(n)=PEO2(n-1)
             endif
          endif
          IF(pres(n) < 0.4e-4) PEN2(n)=PEN2(n-1)
          !IF(pres(n) < 0.4e-4) write(6,*) 'here', n, pres(n)
          
       ENDDO
       
       !  call Smooft(peo2,ht_dim,sm_dim1)
       !  call Smooft(peop,ht_dim,sm_dim2)
       !  call Smooft(pen2,ht_dim,sm_dim2)
       
       DO N=HT_DIM,1,-1   ! Tidy up bottom.
          if(ht1d(n) < 100.e3) then
             
             peop(n)=2.*peop(n+1)-peop(n+2)
             if(peop(n).le.0.) peop(n)=0.

             peo2(n)=2.*peo2(n+1)-peo2(n+2)
             if(peo2(n).le.0.) peo2(n)=0.
             
             pen2(n)=2.*pen2(n+1)-pen2(n+2)
             if(pen2(n) .le.0.) pen2(n)=0.
             
          endif
       ENDDO
    ENDIF
   
    ! Add to the jx array
 
    DO N=1, HT_DIM
       
       ! photolysis
       jx(n,m,ip_O2_b1,l) = jx(n,m,ip_O2_b1b2,l)*branch(b_jio2_b1,n)
       jx(n,m,ip_O2_b2,l) = jx(n,m,ip_O2_b1b2,l)*branch(b_jio2_b2,n)

       jx(n,m,ip_N2_b1,l) = jx(n,m,ip_N2,l)*branch(b_jin2_b1,n)
       jx(n,m,ip_N2_b2,l) = jx(n,m,ip_N2,l)*branch(b_jin2_b2,n)
       jx(n,m,ip_N2_b3,l) = jx(n,m,ip_N2,l)*branch(b_jin2_b3,n)

       ! secondary electrons (se)
       jx(n,m,ip_se_Op_em,l) = jx(n,m,ip_O3Pp,l)*PEOP(n)

       jx(n,m,ip_se_O2_b1,l) = jx(n,m,ip_O2_b1b2,l)*PEO2(n)*branch(b_seo2_b1,n)
       jx(n,m,ip_se_O2_b2,l) = jx(n,m,ip_O2_b1b2,l)*PEO2(n)*branch(b_seo2_b2,n)

       jx(n,m,ip_se_N2_b1,l) = jx(n,m,ip_N2,l)*PEN2(n)*branch(b_sen2_b1,n)
       jx(n,m,ip_se_N2_b2,l) = jx(n,m,ip_N2,l)*PEN2(n)*branch(b_sen2_b2,n)
       jx(n,m,ip_se_N2_b3,l) = jx(n,m,ip_N2,l)*PEN2(n)*branch(b_sen2_b3,n)
       jx(n,m,ip_se_N2_b4,l) = jx(n,m,ip_N2,l)*PEN2(n)*branch(b_sen2_b4,n)

       ! auroral electrons e*
       ! aurqo: O + e* -> O+ + 2*e-
       ! not necessary here because can already be done in ionparams
       ! jx(n,m,ip_Op_em,l) = jx(n,m,ip_Op_em,l) 

       ! rt(aurqo2_b1) O2 + e* -> O2+ + 2*e- ht_dep_branch*                         
       ! rt(aurqo2_b2) O2 + e* -> O+ + O + e- ht_dep_branch*                        
       jx(n,m,ip_O2p_em,l)  = jx(n,m,ip_O2_aurq,l)*branch(b_aurqo2_b1,n) ! 
       jx(n,m,ip_Op_O_em,l) = jx(n,m,ip_O2_aurq,l)*branch(b_aurqo2_b2,n) ! 

       ! rt(aurqn2_b1) N2 + e* -> N2+ + 2*e- ht_dep_branch*                         
       ! rt(aurqn2_b2) N2 + e* -> N+ + N{4s} + 2*e- ht_dep_branch*                  
       ! rt(aurqn2_b3) N2 + e* -> N+ + N{2d} + 2*e- ht_dep_branch*                  
       ! rt(aurqn2_b4) N2 + e* -> 0.47*N{4s} + 0.53*N{2d} + e- ht_dep_branch*       
       jx(n,m,ip_N2p_em,l)    = jx(n,m,ip_N2_aurq,l)*branch(b_aurqn2_b1,n) ! 
       jx(n,m,ip_Np_N_em,l)   = jx(n,m,ip_N2_aurq,l)*branch(b_aurqn2_b2,n) ! 
       jx(n,m,ip_Np_N2D_em,l) = jx(n,m,ip_N2_aurq,l)*branch(b_aurqn2_b3,n) ! 
       jx(n,m,ip_N_N2D_em,l)  = jx(n,m,ip_N2_aurq,l)*branch(b_aurqn2_b4,n) ! 


       ! mz_ab_20110612 deleted
!!$       if (jx(n,m,ip_se_Op_em,l) < 1.e-18) jx(n,m,ip_se_Op_em,l)=0.
!!$       if (jx(n,m,ip_se_O2,l) < 1.e-18) jx(n,m,ip_se_O2,l)=0.
!!$       if (jx(n,m,ip_se_N2,l) < 1.e-18) jx(n,m,ip_se_N2,l)=0.
       
       
       ! write out jx
       !  IF(m == 45 .and. l==1) then
       !   write(6,"(a26,i3,3e12.4)") "jx 2, N2, O2, O2  ", n,        &
       !                      jx(n,m,l,ip_N2P),                      &
       !                      jx(n,m,l,ip_O2P),                      &
       !                      jx(n,m,l,ip_O3Pp)
       !  ENDIF
       
    ENDDO
    
    initialise = 0
    
  END SUBROUTINE Secondaries


  !=======================================================
  !=                   CO2ChapManHeating                 =
  !=======================================================

  ! INPUT:

  ! den1d(ht_dim)               1-dimensional density profile
  ! cp1d(ht_dim)                1-dimensional heat capacity profile
  ! ht1d(ht_dim)                array with heights

  ! OUPUT:

  ! co2_chap_heat               heating rate due absorption of CO2
  !                             modeled as Chapman layer
  !                             actually now used as O3 heating

  SUBROUTINE CO2ChapManHeating(m,l,den1d, cp1d, ht1d, ht_dim, csza, co2_chap_heat)

      INTEGER,  INTENT(IN)  :: m, l, ht_dim
      REAL(dp), INTENT(IN)  :: cp1d(ht_dim)
      REAL(dp), INTENT(IN)  :: den1d(ht_dim)
      REAL(dp), INTENT(IN)  :: ht1d(ht_dim)
      REAL(dp), INTENT(IN)  :: csza

      REAL                       qm0,hscale         ! maximum production rate, scale height
      REAL                       hmax1,hmax2        ! heights of max. production
      REAL                       en1,en2            ! energies available for heating
      INTEGER                    nlay               ! number of Chapman layers
      REAL                    :: hmax(1:2)          ! array with heights of max. production
      REAL                    :: energy_in(1:2)     ! array with the input energy for each layer
      REAL                    :: z(1:2,1:ht_dim)    ! altitude (dimensionless)
      real                    :: q(1:2,1:ht_dim)    ! heat production (sza,z), dimensionless
      real                    :: qnew(1:2,1:ht_dim) ! halfway scaled heat production for
      real                    :: qtot(1:ht_dim)     ! qnew summed over all layers
      INTEGER                    n,i
      REAL                       csza_loc, MAXVAL, MINVAL
      REAL                       func_n,func_np1       ! chapman production function at n and n+1
      REAL                       zdum,dz               ! dimensionless height + interval
      REAL                       local_integral        ! integral of chapman function from z(n)
      REAL                       local_integral2       ! to z(n+1)
      REAL                       integral,integral2    ! total integral of chapman function from
                                                       ! top of atmosphere to lower boundary

      REAL(dp), INTENT(OUT) :: co2_chap_heat(ht_dim) ! total scaled heat production converted to
                                                       ! units of K day-1

      !WRITE(6,*) "DEBUG 0 - subroutine CO2ChapManHeating starting..."

      ! define heights of maximum production rate (hmax), maximum production rate
      ! (qm0), scale height (hscale), available energies for each layer, and number
      ! of layers
      nlay      = 2             ! number of Chapman layers
      qm0       = 1.0
      hscale    =  7500.        ! scale height (m)
      hmax1     = 50000.        ! height of maximum production of Chapman layer 1 (m)
      hmax2     = 90000.        ! height of maximum production of Chapman layer 2 (m)
      ! mjh changed
      !hmax      = [hmax1,hmax2]
      hmax(1)   = hmax1
      hmax(2)   = hmax2
      en1       = 5.0           ! input energy of Chapman layer 1 (Wm-2)
      en2       = 5.0           ! input energy of Chapman layer 2 (Wm-2) ----> BOTH CHANGED HERE
                                ! ----> FOR OZONE TESTING!!!

      ! mjh changed.
      ! energy_in = [en1,en2]
      energy_in(1) = en1
      energy_in(1) = en2

      ! Make sure model doesn't crash for small csza

      csza_loc = csza
      IF(csza .lt. 0.01) THEN 
         csza_loc = 0.01
      ENDIF
 
      ! calculate production rate as function of height
      do i=1,nlay
         do n=1,ht_dim          !loop over altitudes (m)
            z(i,n)=(ht1d(n)-hmax(i))/hscale
            q(i,n)=qm0*exp(1-z(i,n)-(1/csza_loc)*exp(-z(i,n)))
         enddo
      enddo

!!$      IF(m == INT(lat_dim/2.) .AND. l == 10) THEN
!!$
!!$#ifdef MESSY
!!$         IF (p_parallel_io) &
!!$#endif
!!$         WRITE(6,*) 'm = ',m,' l = ',l,'ht_dim = ',ht_dim,' hmax1 = ',hmax(1), &
!!$                   ' z(layer1) = ',z(1,20:25),' q(layer1) = ',q(1,20:25), "DEBUG1"
!!$      ENDIF

      ! define scaling factor for production rate to heating rate
      !      according to 5 Wm-2 used for heating (function 1 - 200-300 nm)
      !      and to 9.0e-2 Wm-2 (function 2 - 1-20 microns)
      ! CHANGED HERE FOR OZONE TESTING TO BOTH 5 W/M2

      ! calculate integral of standard production rate curves (chapman functions)
      ! integrate from dimensionless zdum = -25 to +25 (max. prod. at z = 0)

      integral = 0.

      do n=1,101
         
         zdum = (n - 51.)/2
         dz = 0.5
         ! func=exp(1-zdum-(1/cos(sza*DTR))*exp(-zdum))

         func_n    = exp(1.-zdum-(1./csza_loc)*exp(-zdum))
         func_np1  = exp(1.-(zdum+dz)-(1./csza_loc)*exp(-zdum+dz))
         
         local_integral = func_n*dz + (func_np1-func_n)*dz/2.

         integral = integral + local_integral

      enddo

!!$      IF(m == INT(lat_dim/2.) .AND. l == 10) THEN
!!$
!!$#ifdef MESSY
!!$         IF (p_parallel_io) &
!!$#endif
!!$         WRITE(6,*) 'csza = ',csza_loc,' integral = ',integral, "DEBUG3"
!!$
!!$      ENDIF

      ! For a CO2 density of 1.6 kg/m3 at the surface and exponentially falling of,
      ! there is roughly 5 Wm-2 availabe for heating in the band from 200 to 300 nm
      ! due to absorption by CO2.
      ! For the same CO2 density profile as above, there is roughly 9.0e-2 Wm-2 available
      ! for heating in the band from 1 to 20 microns due to absorption by CO2
      ! BUT WE'RE DOING O3 TESTING NOW WITH 2 CHAPMAN LAYERS OF 5 W/M2

      ! Apply scaling factors of energy_in(i)/(hscale*integral value)
      ! The factor 1/hscale comes from scaling from dimensionless quantity to
      ! a quantity with a dimension: dz = (1/hscale)*dh 
      ! Units of qnew are W/m2
      do i=1,nlay
         do n=1,ht_dim
            qnew(i,n)=q(i,n)*energy_in(i)/(hscale*integral)
            !print,'q=',q(i,n),' qnew=',qnew(i,n)
         enddo
      enddo

       
      ! Test if integral of qnew gives back the energy input
      
      integral2 = 0.

      do n=1,101

         zdum = (n - 51.)/2
         dz = 0.5
         ! func=exp(1-zdum-(1/cos(sza*DTR))*exp(-zdum))

         func_n    = energy_in(1)*exp(1.-zdum-(1./csza_loc)*exp(-zdum))/integral
         func_np1  = energy_in(2)*exp(1.-(zdum+dz)-(1./csza_loc)*exp(-zdum+dz))/integral

         local_integral2 = func_n*dz + (func_np1-func_n)*dz/2.

         integral2 = integral2 + local_integral2

      enddo

!!$      IF(m == INT(lat_dim/2.) .AND. l == 10) THEN
!!$
!!$#ifdef MESSY
!!$         IF (p_parallel_io) &
!!$#endif
!!$         WRITE(6,*) m,l,'energy(layer1) = ',energy_in(1),' qnew(layer1) = ',qnew(1,20:25),&
!!$                    ' qnew(layer2) = ',qnew(2,20:25),' integral2 = ',integral2, "DEBUG4"
!!$
!!$      ENDIF

      ! Add Chapman production functions for 200-300 nm and 1-20 microns
      ! IN THIS CASE NOW TWO PRODUCTION FUNCTIONS OF 5 W/M2 AT 50 AND 90 KM
      qtot(1:ht_dim) = 0.
      do n=1,ht_dim
         do i=1,nlay
            qtot(n)=qtot(n)+qnew(i,n)
         enddo
      enddo

!!$      IF(m == INT(lat_dim/2.) .AND. l == 10) THEN
!!$
!!$         WRITE(6,*) m,l,'qtot = ',qtot(20:25), "DEBUG5"
!!$
!!$      ENDIF

      ! Convert from W/m2 to units of W/kg = J/(kg*s):
      ! Divide by density (kg/m3) and column height(m)
      ! den1d is density profile of mixture
      ! test for strangely high or low values with max and min val

      MAXVAL=-99999.
      MINVAL=99999.

      do n=1,ht_dim-1
         co2_chap_heat(n)=qtot(n)/(den1d(n)*(ht1d(n+1)-ht1d(n)))
         !print,'co2_chap_heat=',co2_chap_heat(n)
         
         IF (co2_chap_heat(n) .gt. MAXVAL) THEN
            MAXVAL=co2_chap_heat(n)
         ENDIF

         IF (co2_chap_heat(n) .lt. MINVAL) THEN
            MINVAL=co2_chap_heat(n)
         ENDIF

      enddo

      ! co2_chap_heat is the total heating rate in J/(kg*s) due to absorption by CO2 in the
      ! 200-300 nm and 1-20 microns wavelength ranges
      ! HERE CHANGED FOR O3 HEATING DUE TO 2 CHAPMAN LAYERS AT 50 AND 90 KM

!!$      IF(m == INT(lat_dim/2.) .AND. l == 10) THEN
!!$
!!$#ifdef MESSY
!!$         IF (p_parallel_io) &
!!$#endif
!!$         WRITE(6,*) m,l,'final heating = ',co2_chap_heat,' maxval = ',MAXVAL,' minval = ',MINVAL, &
!!$          ' density = ',den1d(20:25),' heat capacity = ',cp1d(20:25),'height = ',ht1d(20:25), "DEBUG6"
!!$
!!$      ENDIF

    END SUBROUTINE CO2ChapManHeating

  ! **************************************************************************

    !r========================================================
    !r=                   Chapmann                           =
    !r========================================================
    !r
    !r    Chapmann     Chapmann grazing incidence function from
    !r                 Risbeth and Garriot
    !r
    !r    Input :
    !r
    !r    csza         cos zenith angle
    !r    scht         scht of constituent in question
    !r    ht           height
    !r
    !r    Output :
    !r
    !r    Chapmann cos zenith angle
    !r
    
    REAL(dp) FUNCTION Chapmann(coszen, scht, ht)
      
      REAL(dp), INTENT(IN) :: coszen, ht, scht

      REAL(dp) :: zenangdeg, chap, ht_fac, fc,                              &
                    maxang, tee, seczenang, errfc, zenang

      !IF ( ABS(coszen).LE.0.01 ) THEN
      !   coszen = coszen + 0.015
      !ENDIF
      !seczenang=1./coszen
      !CHAPMANN=1./seczenang
      !RETURN
      
      seczenang=1./coszen
      zenang=ACOS(coszen)
      zenangdeg=180.*(zenang/PI)
      chap=-99999.
      maxang=92.
      
      ! At 500km altitude, zenang is about 95 degrees for sunrise/set
      IF(((zenangdeg) > 85.) .AND. (zenangdeg < maxang)) THEN

         ht_fac=(ht)/scht
         chap=SQRT(0.5*PI*ht_fac*SIN(zenang))
         chap=chap*EXP(0.5*ht_fac*((COS(zenang))**2))
         
         ! error function  of sqrt(0.5*ht_fac*(COS(zenang)**2))
         fc=SQRT(0.5*ht_fac*((COS(zenang))**2))
         tee=1./(1.+0.5*fc)
         errfc=tee*EXP(-fc*fc-1.26551223+tee*(1.00002638+tee*                 &
              (0.37409196+tee*(0.09678418+tee*(-.18628806+tee*                &
              (.27886807+tee*                                                 &
              (-1.13520398+tee*(1.48851587+tee*(-.82215223+tee*               &
              .17087277)))))))))
         
         IF(zenangdeg <= 90.)  THEN
            errfc=(1.)-errfc
            chap=chap*(1.-errfc)
         ELSE
            errfc=1.-errfc
            chap=chap*(1.+errfc)
         ENDIF
         IF(chap > 0.)  THEN
            seczenang=chap
         ENDIF

      ENDIF

      CHAPMANN=1./seczenang

    END FUNCTION Chapmann
    
 ! ****************************************************************************

  !r========================================================
  !r=                    BranchingRatios                   =
  !r========================================================
  !r
  !r This routine sets up any height dependent branching
  !r ratios, as defined in rates sheet by ht_dep_branch in
  !r in the first branching rate column. As well as the rates
  !r sheet, you can check the generated Chemistry_gen.h (see
  !r declaration for array "branch".
  !r
  !r Inputs
  !r
  !r pres(ht_dim)   : pressure height array
  !r htav(ht_dim)   : Global mean height
  !r
  !r Output
  !r
  !r branch(ht_branches,ht_dim) : branching ratios (passed via module scope.
  !r 
  !r
  !r


  SUBROUTINE BranchingRatios(ht_dim, pres, htav)


    ! Polynomial coefficients for O2 and N2 ionization partitioning.
    ! Taken from CMAT1, Dobbin version, where
    !
    ! EUVNP  ! Fraction of N+ created from photon ionisation of N2
    ! EUVOP  !      "      O+                  "                O2
    !
    ! The original source was Tim Fuller-Rowell "Modelling NO
    ! in the thermosphere" JGR. 1992.
    !
    
    INTEGER, INTENT(IN) :: ht_dim
    REAL(dp) :: pres(ht_dim), htav(ht_dim)

    REAL(dp):: EUVOPPOL(11) = (/0.00032620530,-0.013361594,0.0075514941,&
         0.0091426697,-0.0080997479,0.0021838333,&
         -0.00019717393,-9.7139409e-06,3.0529147e-06&
         ,-2.0378266e-07,4.5798682e-09/)    

    REAL(dp):: EUVNPPOL(6) = (/ 0.38940824,0.0023363129,0.011755605 &
         ,-0.0044112183,0.00040532805,-1.1299813e-05/)


    ! Some local variables

    REAL(dp) :: euvop(ht_dim), exex
    REAL(dp) :: euvnp(ht_dim)
    INTEGER    :: i,n 

    ! The ratio N2d to N4s when ion product of N2 ionization is N+

    REAL(dp) :: ratio_n2d_to_n4s

    euvop = 0.
    euvnp = 0.

   ! WRITE(6,*) "Setting up height dependent ionization branching ratios ..."

    DO n=1,ht_dim

       IF(pres(n) <= 1.0376) THEN
          DO i=1,11
             exex=+1.*LOG(1.0376)-1.*LOG(pres(n))
             EUVOP(n)=EUVOP(n)+(EUVOPPOL(i)*(exex**(i-1.)))
             IF(i < 7) EUVNP(n)=EUVNP(n)+(EUVNPPOL(i)*(exex**(i-1.)))
          ENDDO
          IF(exex.gt.12.96) THEN
             EUVOP(n)=EUVOP(n-1)
             EUVNP(n)=EUVNP(n-1)
          ENDIF
       ENDIF

    ENDDO

    DO n=ht_dim,1,-1   ! Tidy up bottom.
       IF(htav(n) < 100.e3) THEN
          EUVNP(n)=EUVNP(n+1)
          EUVOP(n)=EUVOP(n+1)
       ENDIF
    ENDDO

    ! ------------------------------------

    ! Ratio of O+ to total ionization of O2
    ! (by auroral particles and photons)
    ! b_aurqo2_b2 O2 + e* -> O+ + O + 2*e-
    ! b_jio2_b2   O2 + hv(<105nm) -> O+ + O + e*

     branch(b_aurqo2_b2,:)  = 0.33
     branch(b_jio2_b2,:)    = euvop(:)

    ! -------------------------------------

    ! Ratio O2+ to total ionization of O2
    ! (by auroral particles and photons)
    ! b_aurqo2_b1 O2 + e* -> O2+ + 2*e-
    ! b_jio2_b1   O2 + hv(<105nm) -> O2+ + e*

     branch(b_aurqo2_b1,:)  = 1. - branch(b_aurqo2_b2,:)
     branch(b_jio2_b1,:)  = 1. - euvop(:)

    ! -------------------------------------

    ! Ratio N2+ to ionization of N2
    ! b_aurqn2_b1 N2 + e* -> N2+ + 2*e-
    ! b_jin2_b1   N2 + hv(<105nm) -> N2+ + e*

     branch(b_aurqn2_b1,:) = 0.76
     branch(b_jin2_b1,:) = 1. - euvnp(:)

    ! -------------------------------------

    ! Ratio of N+ to ionization of N2. Note this is
    ! further partitioned to N4s and N2d
    ! (by auroral particles and photons)
    ! b_aurqn2_b2 N2 + e* -> N+ + N{4s} + 2*e-
    ! b_aurqn2_b3 N2 + e* -> N+ + N{2d} + 2*e-
    ! b_jin2_b2   N2 + hv(<105nm) -> N+ + N{4s} + e*
    ! b_jin2_b3   N2 + hv(<105nm) -> N+ + N{2d} + e*

     ratio_n2d_to_n4s = 0.5

     branch(b_aurqn2_b2,:) = (1 - ratio_n2d_to_n4s)*(1-branch(b_aurqn2_b1,:))
     branch(b_aurqn2_b3,:) = ratio_n2d_to_n4s* (1-branch(b_aurqn2_b1,:))
     branch(b_jin2_b2,:)   = (1 - ratio_n2d_to_n4s)*euvnp(:)
     branch(b_jin2_b3,:)   = ratio_n2d_to_n4s*euvnp(:)

    ! -------------------------------------

    ! Ratio of O2+ to total ionisation of O2 by secondary electrons
    ! b_seo2_b2  O2 + e* -> O2+ + 2*e-

      branch(b_seo2_b1,:) = 0.67
     
    ! Ratio of O+ to total ionisation of O2 by secondary electrons
    ! b_seo2_b1  O2 + e* -> O+ + O + 2*e-

      branch(b_seo2_b2,:) = 1. - branch(b_seo2_b1,:)

    ! Ratio N2+ to ionization of N2 by secondary electrons
    ! b_sen2_b1 N2 + e* -> N2+ + 2*e-

      branch(b_sen2_b1,:) = 0.76

    ! -------------------------------------

    ! Ratio of N+ to ionization of N2 by secondary electrons. Note this is
    ! further partitioned to N4s and N2d
    ! b_sen2_b2 N2 + e* -> N+ + N{4s} + 2*e-
    ! b_sen2_b3 N2 + e* -> N+ + N{2d} + 2*e-

       branch(b_sen2_b2,:) = (1 - ratio_n2d_to_n4s)*(1-branch(b_sen2_b1,:))
       branch(b_sen2_b3,:) = ratio_n2d_to_n4s* (1-branch(b_sen2_b1,:))

    ! -------------------------------------

    ! Ratio of N(2D:4S) to dissociation of N2 by auroral particle precipitation
    ! and secondary electrons. 
    ! b_aurqn2_b4 N2 + e* -> N{4s} + N{2d} + *e-
    ! b_sen2_b4 N2 + e* -> N{4s} + N{2d} + *e-
    ! NB the photon equivalent of this reaction is predissociation of N2

      branch(b_aurqn2_b4,:) = 1.34
      branch(b_sen2_b4,:) = 1.34


  END SUBROUTINE BranchingRatios

END MODULE messy_radjimt

! ****************************************************************************
