! ***********************************************************************
MODULE messy_spe
  ! ***********************************************************************
  ! Solar proton event parameterization
  !
  ! Authors: Andreas Baumgaertner,
  ! MPI-CH Mainz (abaumg@mpch-mainz.mpg.de), Jan 2007
  ! Authors: Stefan Versick,
  ! KIT-SCC & KIT-IMK-ASF (stefan.versick@kit.edu), Jul 2016
  ! Sabine Barthlott
  ! KIT-IMK-ASF (sabine.barthlott@kit.edu), Feb 2016
  ! Literature:
  ! - Vitt and Jackman (1996), "A comparison of sources of odd nitrogen production from
  !   1974 through 1993 in the Earth's middle atmosphere as calculated using a
  !   two-dimensional model", JGR 101(D3), 6729-6739
  ! - Armstrong et al. (1989), "Interplanetary Energetic Ions and Polar Radio Wave 
  !   Absorption", JGR 94(A4), 3543-3554
  ! - Reid et al. (1991), "Reponse of the middle atmosphere to the solar proton events
  !   of August-December, 1989), GRL 18(6), 1019-1022
  ! - Jackman, C. H., et al. (2005), Neutral atmospheric influences of the solar proton events
  !   in October-November 2003, Journal of Geophysical Research-Space Physics, 110(A9), 10.
  ! - Bethe, H.A.; Ashkin, J.; in: Segre, E., Experimental Nuclear Physics, Volume I, 
  !   Chapt II Experimental Nuclear Physics, 166-, John Wiley and Sons, New York
  ! - Nieder, H.: Modellstudien zur Untersuchung des Einflusses solarer Prozesse auf die mittlere Atmosphaere, 
  !   phd-Thesis, Karlsruhe Institute of Technology, 2015 .

  ! TODO:
  !    see TODO list in messy_spe_e5.f90

  USE messy_main_constants_mem,  ONLY: DP
 

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'spe'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.1'

  INTEGER,   PUBLIC :: spe_method = 0       ! internal (0) or external (1) ionization rates
  ! ka_sv_20171213+
  LOGICAL,   PUBLIC :: calc_pos_ions = .FALSE.       ! calculate ions directly? default: no; if true be sure that you have selected ion reactions in your meccanism
  ! ka_sv_20171213+
  ! ka_sv_20180327+
  LOGICAL,   PUBLIC :: calc_three_o_ions = .FALSE.       ! calculate only one O+ ion or all three seperate?
  REAL(dp), PUBLIC :: br_n4s, br_n2d
  ! ka_sv_20180327-
  ! FOR METHOD 0
  REAL(dp),   DIMENSION(:), POINTER,  PUBLIC :: chrig  => NULL() ! rigidity of each channel [MeV]
  INTEGER,   PUBLIC :: npch = 7       ! number of proton flux channels of the instrument
  REAL(DP),  PUBLIC :: rigres ! rigidity spectrum resolution

  ! energies for interpolated spectrum
  REAL(dp),  DIMENSION(:),  ALLOCATABLE,  PUBLIC :: energies   
  ! range for energies at STP 
  REAL(dp),  DIMENSION(:),  ALLOCATABLE,  PUBLIC :: rangestp
  ! proton flux (instrument channels)
  REAL(dp),  DIMENSION(:), POINTER,       PUBLIC :: pflx => NULL()
  ! interpolated proton flux spectrum
  REAL(dp),  DIMENSION(:), ALLOCATABLE,   PUBLIC :: pflxspc 

  INTEGER, PUBLIC :: k_max = 25 ! number of pitch angles

  ! FOR METHOD 1
  ! ionization rates (calculated externally)
  REAL(dp),  DIMENSION(:), POINTER,       PUBLIC :: ions_ext => NULL()
  INTEGER,   PUBLIC :: naltitudes=0

  ! FOR METHOD 2
  REAL(dp),  DIMENSION(:,:,:), POINTER,     PUBLIC :: AIMOS_p => NULL() ! protons
  REAL(dp),  DIMENSION(:,:,:), POINTER,     PUBLIC :: AIMOS_e => NULL() ! electrons
  REAL(dp),  DIMENSION(:,:,:), POINTER,     PUBLIC :: AIMOS_a => NULL() ! alpha particles
  REAL(dp),  DIMENSION(:,:,:), POINTER,     PUBLIC :: AIMOS_ionrate   => NULL() ! total

  ! FOR METHOD 3

!ka_sb_20160722+
  CHARACTER(LEN=200), PUBLIC :: input_spe3_file = ''
!ka_sb_20160722-
  
!ka_sv_20160615+  
  ! FOR METHOD 4
  REAL(dp),  DIMENSION(:,:,:), POINTER,     PUBLIC :: CMIP6_p => NULL() ! protons
  REAL(dp),  DIMENSION(:,:,:), POINTER,     PUBLIC :: CMIP6_e => NULL() ! electrons
  REAL(dp),  DIMENSION(:,:,:), POINTER,     PUBLIC :: CMIP6_g => NULL() ! galactic cosmic rays particles
  REAL(dp),  DIMENSION(:,:,:), POINTER,     PUBLIC :: CMIP6_ionrate   => NULL() ! total
!ka_sv_20160615-

!ka_sb_20160222+
  INTEGER,   PUBLIC :: npe_method = 0       ! constant (0), internally set (1) or holger nieders (2) NOX production rates
  CHARACTER(LEN=200), PUBLIC :: input_sap_file = ''
  INTEGER,   PUBLIC :: switchphi = 0        ! photoionization on (1) or off (0)
!ka_sb_20160222-
!ka_sv_20170426+
  INTEGER, PUBLIC :: n2oprod = 0
!ka_sv_20170426-

!ka_sb_20160222 FOR METHOD 0 and METHOD 1
!ka_sb_20160222 internally set prodution rates (actually: Jackman)
  ! MAX LENGTH OF CTRL NML INPUT ion_km, Nperion_km, NOperion_km  
  INTEGER, PARAMETER, PUBLIC :: NMAXIONKM = 50

  ! GLOBAL CTRL NAMELIST VARIABLES
  ! ionization rate data file
!csv  CHARACTER(LEN=200), PUBLIC :: ion_data_file = ''   ! kit_sv_20140320: never used
  ! spe data file
!csv  CHARACTER(LEN=200), PUBLIC :: spe_data_file = ''   ! kit_sv_20140320: never used
  ! full ionization for |geomagnlatitude| >= r_lat1 [deg]
  REAL(dp), PUBLIC :: r_lat1 = 60.0_dp
  ! no   ionization for |geomagnlatitude| <  r_lat2 [deg]
  REAL(dp), PUBLIC :: r_lat2 = 55.0_dp
  INTEGER,  PUBLIC :: nbin = 50 ! number of bins
  INTEGER,  PUBLIC :: minoffs=5  ! determines smallest energy of interpolated spectrum
  INTEGER,  PUBLIC :: spect_interp_method=0
  ! parameters of the range-energy relationship by Bethe (1953)
  REAL(dp), PUBLIC :: rangerela=9.3_dp, rangerelb=1.8_dp
  REAL(dp), PUBLIC, DIMENSION(NMAXIONKM), SAVE :: ion_km=-999., Nperion_km=-999., NOperion_km= -999.
  
!ka_sb_20160222 FOR METHOD 2
!ka_sb_20160222 linearly interpolated SAP NOx production rates by Holger Nieder

  ! op_pj_20170220+
  REAL(dp), DIMENSION(:), POINTER :: SAP_Np  => NULL()  !DIONP-CONCM-TEMP-H2O-O-N-NO-NOpf-Nf-NMf-OMf
  REAL(dp), DIMENSION(:), POINTER :: SAP_NOp => NULL() !DIONP-CONCM-TEMP-H2O-O-N-NO-NOpf-Nf-NMf-OMf
  ! op_pj_20170220-
  
!ka_sv_20180514+
  CHARACTER(LEN=200), PUBLIC :: spe_aimos_dir, spe_aimos_prefix
  INTEGER, PUBLIC :: spe_nlat, spe_nlon, spe_nlev, aimos_time
  LOGICAL, PUBLIC :: aimpro, aimele, aimalp
  LOGICAL, PUBLIC :: spe_hox
!ka_sv_20180514-
  
  ! END GLOBAL CTRL NAMELIST VARIABLES

  PUBLIC :: SPE_read_nml_ctrl
  PUBLIC :: spe_clean
  PUBLIC :: spe_interp
  PUBLIC :: SPE_NPE
  PUBLIC :: SPE_MultiLinPhion
  PUBLIC :: SPE_ReadNPE
!ka_sb_20160222+
  PUBLIC :: SPE_IONS
  PUBLIC :: SPE_PROD_XNOX
!ka_sb_20160222- 
!ka_sb_20160706+
  PUBLIC :: SPE_READ_IONRATE
  PUBLIC :: HN_ION_PROVIDE_DATA
  PUBLIC :: SPE_INTERP_EXT_IONRATE
!ka_sb_20160706-
!ka_sb_20161114+
  PUBLIC :: SPE_PHIONIZ
!ka_sb_20161114-
!ka_sv_20180514+
  PUBLIC :: spe_read_aimos_orig
!ka_sv_20180514-

CONTAINS


! =========================================================================
 SUBROUTINE spe_interp

    ! ****************************************************************
    ! Interpolate measured integral proton fluxes to a spectrum
    ! ****************************************************************

    IMPLICIT NONE

    ! local variables
    INTEGER :: pch     ! index of measurement channel
    INTEGER :: bin
    REAL(DP) :: maxrig=850._dp ! maximum rigidity in spectrum
    REAL(DP) :: rig1,rig2,p,p0,p1,p2,j0,j1,j2,p0m,j0m
    REAL(dp),  DIMENSION(:),  ALLOCATABLE :: gamma,A ! for method 2  

    INTRINSIC EXP, LOG, REAL, SQRT

    ! INITIALIZE

    rigres=LOG(maxrig)/REAL(nbin,dp)

    !***********************************************************
    ! Calculate energy levels 
    !***********************************************************
    DO bin=1, nbin  ! loop over all energies
       rig1=EXP(REAL(bin-1-minoffs,dp)*rigres)
       rig2=EXP(REAL(bin-minoffs,dp)*rigres)
       energies(bin)=SQRT(rig1*rig2) 
    END DO

    ! Calculate range for energies at STP, after Bethe (1953)
    rangestp=(energies/rangerela)**rangerelb ! [m]

    !**************************************************
    ! Calculate spectra pflxspc
    !**************************************************
    ! Initialize final array
    DO bin=1, nbin
        pflxspc(bin)=0._dp
    END DO
    ! loop over all measurements
    SELECT CASE(spect_interp_method)
    CASE(0)
       ! extrapolate and interpolate between measurements, smooth spectrum
       p0m=0._dp
       j0m=0._dp
       DO pch=1, npch-1
          p1=REAL(chrig(pch),dp)
          p2=REAL(chrig(pch+1),dp)
          j1=pflx(pch)
          j2=pflx(pch+1)
          IF ((j2/j1).eq.1._dp) j2=j1*1.01_dp ! modify to avoid div by 0
          p0=(p1-p2)/LOG(j2/j1)
          p0m=p0m+p0
          j0=j1*EXP(p1/p0)
          j0m=j0m+j0
       END DO
       p0m=p0m/REAL(npch-1,dp) ! arithmetic mean of p0
       j0m=j0m/REAL(npch-1,dp)
       DO bin=1, nbin
          p=energies(bin)
          pflxspc(bin)=j0m*EXP(-p/p0m)
       END DO
       ! convert from integral flux spectrum to flux spectrum
       DO bin=1, nbin-1
          pflxspc(bin)=(pflxspc(bin)-pflxspc(bin+1))
          IF (pflxspc(bin)<0._dp) pflxspc(bin)=0._dp
       END DO
       pflxspc(nbin)=pflxspc(nbin-1)
    CASE(1)
       ! extrapolate and interpolate between measurements, rough spectrum, still experimental
       pch=1
       DO bin=1, nbin
          p=energies(bin)
          IF ((p>REAL(chrig(pch+1))).AND.(pch<npch-1))  pch=pch+1
          p1=REAL(chrig(pch),dp)
          p2=REAL(chrig(pch+1),dp)
          j1=pflx(pch)
          j2=pflx(pch+1)
          IF ((j2/j1).eq.1._dp) j2=j1*1.01_dp ! modify to avoid div by 0
          p0=(p1-p2)/(LOG(j2/j1))
          j0=j1*EXP(p1/p0)
          pflxspc(bin)=j0*EXP(-p/p0)
       END DO
       ! convert from integral flux spectrum to flux spectrum
       DO bin=1, nbin-1
          pflxspc(bin)=(pflxspc(bin)-pflxspc(bin+1))
          IF (pflxspc(bin)<0._dp) pflxspc(bin)=0._dp
       END DO
       pflxspc(nbin)=pflxspc(nbin-1)
    CASE(2)
       ! extrapolate and interpolate between measurements, Vitt and Jackman (1996) method
       ALLOCATE(gamma(npch))
       ALLOCATE(A(npch))
       DO pch=1, npch-1
          p1=REAL(chrig(pch),dp)
          p2=REAL(chrig(pch+1),dp)
          j1=pflx(pch)
          j2=pflx(pch+1)
!          IF ((j1==0._dp).OR.(j2==0._dp)) THEN
          IF (j1==0._dp) THEN
             gamma(pch)=0
             A(pch)=0
             continue
          END IF
          IF ((j1/j2).eq.1._dp) j1=j2*1.01_dp ! modify to avoid div by 0
          gamma(pch)=1-LOG(j1/j2)/log(p1/p2) 
          IF ((gamma(pch).le.1.001).AND.(gamma(pch).ge.0.999)) gamma(pch)=0
          A(pch)=j1/(p1**(1-gamma(pch)))*(gamma(pch)-1)
       END DO
       DO bin=1, nbin
          rig1=EXP(REAL(bin-1-minoffs,dp)*rigres)
          rig2=EXP(REAL(bin-minoffs,dp)*rigres)
          DO pch=npch-1, 2, -1
             IF (energies(bin)>=chrig(pch)) EXIT
          END DO
          pflxspc(bin)=A(pch)/(1-gamma(pch))* (rig2**(1-gamma(pch))-rig1**(1-gamma(pch)))
       END DO
       DEALLOCATE(gamma)
       DEALLOCATE(A)
    END SELECT


  END SUBROUTINE spe_interp
  ! =========================================================================


  ! =========================================================================
  SUBROUTINE SPE_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/  r_lat1, r_lat2 &
         , spe_method &
!ka_sv_20171213+
         , calc_pos_ions &
!ka_sv_20171213-
!ka_sv_20180327+
         , calc_three_o_ions &
         , br_n4s, br_n2d &
!ka_sv_20180327-
         , nbin, minoffs, spect_interp_method, rangerela, rangerelb &
!ka_sb_20160722+ 
         , npe_method, input_sap_file & 
         , switchphi  &
         , input_spe3_file &
!ka_sv_20170426+
         , n2oprod &
!ka_sv_20170426-
!ka_sb_20160722-         
         , ion_km, Nperion_km, NOperion_km &
!ka_sv_20180514+
         , spe_aimos_dir, spe_aimos_prefix, spe_nlat, spe_nlon, spe_nlev, aimos_time, aimpro, aimele, aimalp, spe_hox
!ka_sv_20180514-
         

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='spe_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    WRITE(*,*) 'SPE_METHOD            : ',spe_method
    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY      
       SELECT CASE(spect_interp_method)
       CASE(0)   
          WRITE(*,*) 'METHOD: extrapolate and interpolate between measurements, smooth spectrum'
       CASE(1) 
          WRITE(*,*) 'METHOD: extrapolate and interpolate between measurements, rough spectrum, still experimental'
       CASE(2) 
          WRITE(*,*) 'METHOD: extrapolate and interpolate between measurements, Vitt and Jackman'
       CASE DEFAULT
          WRITE(*,*) ' ERROR: unknown spect_interp_method!'
          RETURN
       END SELECT
    CASE(1) ! EXTERNAL IONIZATION RATES
       WRITE(*,*) ' EXTERNAL IONIZATION RATES'       
    CASE(2) ! EXTERNAL IONIZATION RATES FROM AIMOS / JAN MAIK WISSING
       WRITE(*,*) ' USING AIMOS IONIZATION RATES'
! ka_sb_20160705+
    CASE(3) ! EXTERNAL IONIZATION RATES FROM AIMOS / HOLGER NIEDER
        WRITE(*,*) ' USING H. NIEDER AIMOS IONIZATION RATES'
        WRITE(*,*) 'file used: ',input_spe3_file
!ka_sv_20171213+
        IF (calc_pos_ions) THEN
          WRITE(*,*) 'N2D, N4S, N2+, N+, O2+, O+ and NO+ are calculated directly'
          IF (calc_three_o_ions) THEN
            WRITE(*,*) 'Instead of O+ O(4S)+, O(2D)+ and O(2P)+ are calculated'
          ELSE
            WRITE(*,*) 'Only O+ is calculated'
          END IF
        ELSE
          WRITE(*,*) 'N and NO are parametrized. No Ions are calculated'
        END IF
!ka_sv_20171213-
! ka_sb_20160705-
! ka_sv_20160615+
    CASE(4) ! EXTERNAL FORM CMIP6
       WRITE(*,*) ' USING CMIP6 IONIZATION RATES'
! ka_sv_20160615-
! ka_sv_20180514+
    CASE(5) ! AIMOS original files
       WRITE(*,*) ' USING ORIGINAL AIMOS FILES'
! ka_sv_20180514-
    CASE DEFAULT
       WRITE(*,*) ' ERROR: unknown spe_method!'
       RETURN
    END SELECT
    
!ka_sb_20160222+
    WRITE(*,*) 'NPE_METHOD            : ',npe_method
    SELECT CASE (npe_method) 
    CASE(0) ! CONSTANT RATES
      WRITE(*,*) 'USING CONSTANT VALUES'
    CASE(1) ! USE values set in namelist
      WRITE(*,*) 'USING values set in namelist'
    CASE(2) ! USE NOx PRODUCTION RATES FROM HOLGER NIEDER
      WRITE(*,*) 'USING NOx production rates by HOLGER NIEDER'   
      WRITE(*,*) 'file used: ',input_sap_file
    CASE DEFAULT
      WRITE(*,*) 'ERROR: unknown npe_method!'
      RETURN
    END SELECT
!ka_sb_20160222-  
!!ka_sb_20161213+
    WRITE(*,*) 'switchphi          :',switchphi
    SELECT CASE (switchphi)
    CASE(0) !no photoionization
      WRITE(*,*) 'NO PHOTOIONIZATION'
    CASE(1) !photoionization
      WRITE(*,*) 'WITH PHOTOIONIZATION'
    CASE DEFAULT
      WRITE(*,*) 'ERROR: unknown switchphi value!'
    RETURN
    END SELECT
!!ka_sb_20161213-    
!ka_sv_20170426+
    WRITE(*,*) 'N2O PRODUCTION DUE TO EEP:',n2oprod
!ka_sv_20170426-

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE SPE_read_nml_ctrl
  ! =========================================================================

  !--------------------------------------------------------------------------
  SUBROUTINE spe_clean

    IMPLICIT NONE

    ! DATA
    IF (ALLOCATED(pflxspc))    DEALLOCATE(pflxspc)
    IF (ALLOCATED(energies))   DEALLOCATE(energies)
    IF (ALLOCATED(rangestp))   DEALLOCATE(rangestp)

    IF (ASSOCIATED(pflx))      DEALLOCATE(pflx)
    NULLIFY(pflx)
    IF (ASSOCIATED(ions_ext))  DEALLOCATE(ions_ext)
    NULLIFY(ions_ext)

    ! op_pj_20170220+
    IF (ASSOCIATED(SAP_Np))  DEALLOCATE(SAP_Np)
    NULLIFY(SAP_Np)
    IF (ASSOCIATED(SAP_NOp)) DEALLOCATE(SAP_NOp)
    NULLIFY(SAP_NOp)
    ! op_pj_20170220-

  END SUBROUTINE spe_clean
  !--------------------------------------------------------------------------

!ka_sb_20160216+
  SUBROUTINE SPE_NPE(TEMP       & !input: temperature [K]
                    ,CONCM        & !input: air density [1/cm3]
                    ,O3P,O1D,NATOM,NO,H2O & !input: vmr
                    ,DIONP        & !input: ionization rate [1/cm3 s] = ions(jp,jk,jrow)
                    ,RNO_POS,RN_4S_POS & !OUTPUT [#/ionpair]
                    ,N2,O2        & !input: vmr
                    ,RO2p_phi,ROp4S_phi,ROp2D_phi,ROp2P_phi,RN2p_phi,RNp_phi,RN4S_phi,RNOp_phi & !input: (values for photoionization) calculated in subroutine spe_phioniz for spe_method = 3 
                    , RO2p, ROp4S, ROp2D, ROp2P, RN2p, RNp, RN4S_out, RNOp, RN2D_out &  ! output
                    ,PRate_phi,PRate_part)     !output

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Subroutine to read precomputed rates
!     and return interpolated values

!     Written in Karlsruhe, April 2012, Holger Nieder
!     Adapted to K3dCTM, March 2013, Holger Nieder
!     Updated for Photoionization, February 2014, Holger Nieder
!     Adapted to form ions directly, February 2016, Holger Nieder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE


  REAL(dp)  :: TEMP,F,CONCM,O3P,O1D,NATOM,NO, &
         H2O,DIONP,JNOion,RNO_POS, RN_4S_POS,DIONP_new

! op_pj_20170220+
!!$  REAL(dp), DIMENSION(1548288) :: SAP_Np !DIONP-CONCM-TEMP-H2O-O-N-NO-NOpf-Nf-NMf-OMf 4*6*7*1*8*6*8*3*2*2*2
!!$  REAL(dp),  DIMENSION(1548288) :: SAP_NOp !DIONP-CONCM-TEMP-H2O-O-N-NO-NOpf-Nf-NMf-OMf
! op_pj_20170220-

  REAL(dp) :: VAL_RNO
  REAL(dp) :: VAL_RN4S

  REAL(dp) :: IND_D,IND_C,IND_N,IND_O,IND_NO,IND_T,IND_H2O
  INTEGER  :: J !Actual Height Level
  INTEGER  :: K !Dummy index
  REAL(dp) :: inGLOBALTIME,inZLEV,inLAT,inLON, &
            inRN2D,inRN4S,inRNO,XX,inDIONP,inJNOi,inNO,inNATOM
  REAL(dp) :: NOpf,Nf,NMf,OMf !Factors of primary ion distribution, important for photoionization and low-EPE
  REAL(dp) :: IND_NOpf,IND_Nf,IND_NMf,IND_OMf

  REAL(dp) :: NOpfac,Nfac,NMolfac,OMolfac !different names for ion distribution inside routine
  REAL(dp) :: N2,O2,EPENRM
  REAL(dp) :: RN2p,RNp,RO2p,ROp4S,ROp2D,ROp2P,RNOp,RN4S !RN4S=RN
  REAL(dp) :: RN2p_phi,RNp_phi,RO2p_phi,ROp4S_phi,ROp2D_phi,ROp2P_phi,RNOp_phi,RN4S_phi !RN4S=RN !rates due to photoionization
  REAL(dp) :: RN2p_part,RNp_part,RO2p_part,ROp4S_part,ROp2D_part,ROp2P_part,RNOp_part,RN4S_part  !rates due to particles
  REAL(dp) :: RN2D_phi,RN2D_part,RN2D,PRate_phi,PRate_part,PRate
  REAL(dp) :: RN4S_out, RN2D_out


      IF(DIONP .LE. 0.0) THEN
       IND_D = 1.0
      ELSE
       IND_D = (LOG10(DIONP/5.0)/LOG10(20.0) + 1.0)

       IND_D = MAX(IND_D,1.0)
       IND_D = MIN(IND_D,4.0)
      ENDIF


      IF (CONCM .LT. 2.0e9) THEN
       IND_C = 1.0
      ELSE
       IND_C = (LOG10(CONCM/(2.0e18-2.0e9)*2.0e18/2.0e9)*9.0 &
               /LOG10(2.0e18/2.0e9)+1.0)
      ENDIF


      IND_C = MAX(IND_C,1.0)
      IND_C = MIN(IND_C,6.0) !Easy way: crop the CONCM-Index here, CONCM <= 6 with phioniz!

      IND_T = ((TEMP-20.0)/60.0)
      IF(TEMP .GT. 320.0) THEN
       IND_T = (TEMP - 320.0)/120.0 + 5.0
      ENDIF
      IF(TEMP .GT. 440.0) THEN
       IND_T = (TEMP - 440.0)/300.0 + 6.0
      ENDIF
      
      IND_T = MAX(IND_T,1.0)
      IND_T = MIN(IND_T,7.0)

      IF((O3P+O1D) .LT. 1.0e-12) THEN
       IND_O = 1.0

      ELSE IF((O3P+O1D) .LT. 4.0e-8) THEN
       IND_O = LOG10((O3P+O1D)/1.0e-12) &
               *1.0/LOG10(4.0e-8/1.0e-12) + 1.0

      ELSE IF ((O3P+O1D) .LE. 0.086) THEN
       IND_O = (LOG10((O3P+O1D)/(0.086-4.0e-8)*0.086 &
               /4.0e-8)*8.0/LOG10(0.086/4.0e-8)+2.0)

       IND_O = MAX(IND_O,1.0)
       IND_O = MIN(IND_O,10.0)
      ELSE IF ((O3P+O1D) .LE. 0.16) THEN
       IND_O = LOG10((O3P+O1D)/0.086) &
               *1.0/LOG10(0.16/0.086) + 10.0

      ELSE IF ((O3P+O1D) .LE. 0.3) THEN
       IND_O = LOG10((O3P+O1D)/0.16) &
               *1.0/LOG10(0.3/0.16) + 11.0

      ELSE
       IND_O = 12.0

      ENDIF

      IND_O = IND_O - 4.0 !1-8 instead of 5-12 for photoionization
      IND_O = MAX(IND_O,1.0)



      IF(NO .LT. 4.0e-15) THEN
       IND_NO = 1.0

      ELSE IF(NO .LT. 4.0e-11) THEN
       IND_NO = LOG10(NO/4.0e-15) &
               *1.0/LOG10(4.0e-11/4.0e-15) + 1.0

      ELSE
       IND_NO = (LOG10(NO/(0.02-4.0e-11)*0.02/4.0e-11) &
                *6.0/LOG10(0.02/4.0e-11)+2.0)

      ENDIF
      IND_NO = MAX(IND_NO,1.0)
      IND_NO = MIN(IND_NO,8.0)


      IF(NATOM .LT. 5.0e-13) THEN
       IND_N = 1.0

      ELSE IF(NATOM .LT. 8.0e-10) THEN
       IND_N = LOG10(NATOM/5.0e-13) &
               *1.0/LOG10(8.0e-10/5.0e-13) + 1.0

      ELSE
       IND_N = (LOG10(NATOM/(8.0e-4-8.0e-10)*8.0e-4 &
               /8.0e-10)*4.0/LOG10(8.0e-4/8.0e-10)+2.0)

      ENDIF
      IND_N = MAX(IND_N,1.0)
      IND_N = MIN(IND_N,6.0)

!ka_sb_20161109+
          RN2D_phi = 0.0 !RNO

          RN2p_part = 0.0_dp
          RNp_part = 0.0
          RO2p_part = 0.0
          ROp4S_part= 0.0
          ROp2D_part= 0.0
          ROp2P_part= 0.0
          RNOp_part = 0.0

          RN4S_part = 0.0 !RN
          RN2D_part = 0.0 !RNO

          RO2p = 0.0
          ROp4S = 0.0
          ROp2D = 0.0
          ROp2P = 0.0

! summing up all rates due to photoionization
          PRate_phi = RN2p_phi + RNp_phi + RO2p_phi + ROp4S_phi + ROp2D_phi + ROp2P_phi + RNOp_phi 

! ka_sv_20170714+
!          EPENRM=0.8978*N2+1.0*O2+0.56*(O3P+O1D) !considers ionization potential of trace gases
          EPENRM=max(0.8978*N2+1.0*O2+0.56*(O3P+O1D),1.e-30) !considers ionization potential of trace gases
! ka_sv_20170714-

 ! particle ionization
          RN2p_part = (0.76*0.8978*N2)/EPENRM*DIONP
          RNp_part  = (0.24*0.8978*N2)/EPENRM*DIONP
          RO2p_part = (0.67*O2)/EPENRM*DIONP
          ROp4S_part= (0.33*O2)/EPENRM*DIONP+(0.56*(O3P+O1D))/EPENRM*DIONP

!! ka_sv_20170915+
!! quick hack
!! ONLY FOR TESTS
!          RN2p_part = 0.0
!          RNp_part = 0.0
!          RO2p_part = 0.0
!          ROp4S_part= 0.0
!          ROp2D_part = 0.0
!          RNOp_part=0.0
!!ka_sv_20170915-

! sum up individual rates of photoionization and particle ionization
          RN2p = RN2p_phi  +RN2p_part
          RNp  = RNp_phi   +RNp_part
          RO2p = RO2p_phi  +RNp_part
          ROp4S= ROp4S_phi +ROp4S_part
          ROp2D= ROp2D_phi +ROp2D_part
! ka_sv_20170502+
          RNOp = RNOp_phi + RNOp_part
! ka_sv_20170502-

 ! primary production rate due to particle ionization 
!csv         RN4S_part = DIONP* (0.538*1.27/1.2) ! Porter 1MeV ohne N+ auf asymptotic hochskaliert, FIXME: fuer "air"
!csv         RN2D_part = DIONP* (0.506*1.27/1.2) ! Porter 1MeV ohne N+ auf asymptotic hochskaliert, FIXME: fuer "air"
    RN4S_part = DIONP*br_n4s
    RN2D_part = DIONP*br_n2d
!! ka_sv_20170915+
!! quick hack
!! ONLY FOR TESTS
!         RN4S_part = 0.0 !RN
!         RN2D_part = 0.0 !RNO
!!ka_sv_20170915-

 ! add primary production rate due to photoionisation
       RN4S = RN4S_part + 1.0*RN4S_phi ! photoionisation produziert N4S
       RN2D = RN2D_part + 0.0*RN4S_phi! photoionisation produziert N2D  ! eventuell durch RN2D_phi ersetzen und Vorfaktor aendern (Stellschraube für NO-Produktion)
       
       ! ka_sv_20180326+
       ! keep correct numbers for direct production; RN4S and RN2D are changed later on in the code
       RN4S_out = RN4S
       RN2D_out = RN2D
       ! ka_sv_20180326-

       PRate_part= RN2p_part + RNp_part + RO2p_part + ROp4S_part + ROp2D_part + ROp2P_part + RNOp_part


       PRate  = PRate_phi + PRate_part


       RN4S = RN4S / MAX(PRate,0.01)
       RN2D = RN2D / MAX(PRate,0.01)

                  NOpfac  = RNOp/MAX(PRate,0.01)
                  Nfac = (RN2p+RNp)/MAX(PRate,0.01)

                  IF((RN2p+RNp) .GT. 0.0) THEN
                  NMolfac = RN2p/(RN2p+RNp)
                  ELSE
                  NMolfac = 0.0
                  ENDIF

                  IF(RO2p+ROp4S+ROp2D+ROp2P .GT. 0.0) THEN
                  OMolfac = RO2p/( RO2p+ROp4S+ROp2D+ROp2P )
                  ELSE
                  OMolfac = 0.0
                  ENDIF

! compute Indices for ion distribution factors
      IND_NOpf = 2.0*NOpfac + 1.0
      IND_Nf   = Nfac + 1.0
      IND_NMf  = NMolfac+ 1.0
      IND_OMf  = OMolfac+ 1.0

! Multilinear interpolation
        CALL SPE_MultiLinPhion(RN_4S_POS, &
             IND_D,IND_C,IND_N,IND_O,IND_NO,IND_T,IND_H2O, &
             IND_NOpf,IND_Nf,IND_NMf,IND_OMf,SAP_Np)

        CALL SPE_MultiLinPhion(RNO_POS, &
             IND_D,IND_C,IND_N,IND_O,IND_NO,IND_T,IND_H2O, &
             IND_NOpf,IND_Nf,IND_NMf,IND_OMf,SAP_NOp)

      RNO_POS = RNO_POS + RN2D
      RN_4S_POS = RN_4S_POS + RN4S

      RETURN

  END SUBROUTINE SPE_NPE
!=========================================================================

!-------------------------------------------------------------------------- 

  SUBROUTINE SPE_MultiLinPhion(VALUE,IND_D,IND_C,IND_N,IND_O, &
                        IND_NO,IND_T,IND_H2O, &
                        IND_NOpf,IND_Nf,IND_NMf,IND_OMf,SAP_X)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Subroutine to find linear interpolated value from SAP-Database

!     Written in Karlsruhe, April 2012, Holger Nieder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE


!!$ REAL(dp), DIMENSION(1548288) :: SAP_X !Sorted as loops: outer loop DIONP, second CONCM, innermost H2O ! op_pj_20170220
 REAL(dp), DIMENSION(:) :: SAP_X !Sorted as loops: outer loop DIONP, second CONCM, innermost H2O  ! op_pj_20170220

    INTEGER :: A,B,C,D,E,F,G,H,I,J

    REAL(dp) :: VALUE
    INTEGER :: Species

    REAL(dp) :: IND_D,IND_C,IND_N,IND_O,IND_NO,IND_T,IND_H2O
    REAL(dp) :: IND_NOpf,IND_Nf,IND_NMf,IND_OMf

    INTEGER*8 :: saparrayindex

    VALUE = 0.0

    do A=0,1
     do B=0,1
      do C=0,1
       do D=0,1
        do E=0,1
         do F=0,1
          do G=0,1
           do H=0,1
            do I=0,1
             do J=0,1

              saparrayindex = 1  &
              +(MIN(INT(IND_D)+A,4)-1) &
              +4*(MIN(INT(IND_C)+B,6)-1) &
              +4*6*(MIN(INT(IND_T)+C,7)-1) &
              +4*6*7*(MIN(INT(IND_O)+D,8)-1) &
              +4*6*7*8*(MIN(INT(IND_N)+E,6)-1) &
              +4*6*7*8*6*(MIN(INT(IND_NO)+F,8)-1) &
              +4*6*7*8*6*8*(MIN(INT(IND_NOpf)+G,3)-1) &
              +4*6*7*8*6*8*3*(MIN(INT(IND_Nf)+H,2)-1) &
              +4*6*7*8*6*8*3*2*(MIN(INT(IND_NMf)+I,2)-1) &
              +4*6*7*8*6*8*3*2*2*(MIN(INT(IND_OMf)+J,2)-1)

              VALUE = VALUE &
               +(1.0-ABS(A-IND_D+INT(IND_D))) &
               *(1.0-ABS(B-IND_C+INT(IND_C))) &
               *(1.0-ABS(C-IND_T+INT(IND_T))) &
               *(1.0-ABS(D-IND_O+INT(IND_O))) &
               *(1.0-ABS(E-IND_N+INT(IND_N))) &
               *(1.0-ABS(F-IND_NO+INT(IND_NO))) &
               *(1.0-ABS(G-IND_NOpf+INT(IND_NOpf))) &
               *(1.0-ABS(H-IND_Nf+INT(IND_Nf))) &
               *(1.0-ABS(I-IND_NMf+INT(IND_NMf))) &
               *(1.0-ABS(J-IND_OMf+INT(IND_OMf))) &
               *SAP_X(saparrayindex)

             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo

  RETURN

  END SUBROUTINE SPE_MultiLinPhion
!=========================================================================


!-------------------------------------------------------------------------- 

! op_pj_20170220+
!!$  SUBROUTINE SPE_ReadNPE(SAP_Np,SAP_NOp) !input_sap_file added 20160923
  SUBROUTINE SPE_ReadNPE
! op_pj_20170220-

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     Subroutine to read precomputed rates

!     Written in Karlsruhe, April 2012, Holger Nieder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT NONE
    
! op_pj_20170220+
!!$    REAL(dp), DIMENSION(1548288) ::  SAP_Np !DIONP-CONCM-TEMP-H2O-O-N-NO-NOpf-Nf-NMf-OMf
!!$    REAL(dp), DIMENSION(1548288) :: SAP_NOp !DIONP-CONCM-TEMP-H2O-O-N-NO-NOpf-Nf-NMf-OMf
! op_pj_20170220-

    INTEGER*8 LAUF00

    REAL :: XX !Dummy to read in unused values

    REAL :: inGLOBALTIME,inZLEV,inLAT,inLON,inRN2D,inRN4S,inRNO,inNO,inNATOM

    REAL :: inRNO2,inRNO_NEG,inRNO2_NEG,inRH_POS,inROH_POS,inROH_NEG, &
       inDIONP,inJNOi,inIONSneNO

    INTEGER :: inindDIONP,inindCONCM,inindTEMP,inindNO,inindO,inindN, &
          inindH2O, indNOpf, indNf, indNMf, indOMf

    REAL :: realNOpf, realNf, realNMf, realOMf

    INTEGER :: K

    INTEGER*8 saparrayindex

    ! op_pj_20170220+
    INTEGER, PARAMETER :: NDIMLEN = 1548288

    ALLOCATE(SAP_Np(NDIMLEN))
    ALLOCATE(SAP_NOp(NDIMLEN))
    ! op_pj_20170220-

! Read out the data from Studies. This code goes into a Subroutine, the Storage into common blocks. To be called at the beginning of the model.
    OPEN(11,FILE=input_sap_file)

! READ out FILE unit 11

!!$ DO LAUF00 =    1, 1548288 !complete array ! op_pj_20170220
    DO LAUF00 =    1, NDIMLEN !complete array ! op_pj_20170220
      SAP_Np  (LAUF00)=0.15
      SAP_NOp (LAUF00)=0.3
    ENDDO !LAUF00

!!$ DO LAUF00 =    1, 1548288 !complete array ! op_pj_20170220
    DO LAUF00 =    1, NDIMLEN !complete array ! op_pj_20170220
      READ(11,*,END=86)inindDIONP,inindCONCM,inindTEMP,inindH2O, &
                     inindO,inindN,inindNO, &
                     realNf,realNMf,realOMf,realNOpf,inRN4S,inRNO, &
                     inRNO2,inRNO_NEG,inRNO2_NEG,inRH_POS, &
                     inROH_POS,inROH_NEG,inDIONP,inJNOi,inIONSneNO

      indNf = INT(realNF) + 1
      indNMf = INT(realNMF) + 1
      indOMf = INT(realOMF) + 1
      indNOpf= INT(2.0*realNOpf) + 1

      saparrayindex = 1 &
        + (inindDIONP-1) + 4*(inindCONCM-1) + 4*6*(inindTEMP-1) &
        + 4*6*7*(inindO-1) + 4*6*7*8*(inindN-1) + 4*6*7*8*6*(inindNO-1)&
        + 4*6*7*8*6*8*(indNOpf-1) + 4*6*7*8*6*8*3*(indNf-1)&
        + 4*6*7*8*6*8*3*2*(indNMf-1) + 4*6*7*8*6*8*3*2*2*(indOMf-1)

      SAP_Np(saparrayindex) = inRN4S
      SAP_NOp(saparrayindex) = inRNO

    ENDDO !LAUF00

  86 CONTINUE

    IF(LAUF00 .LT. 200) THEN
      write(*,*)'*****WARNING: no Database for NOx production!*******'
    ENDIF

!!$ IF(LAUF00 .LT. 1548280) THEN          ! op_pj_20170220
    IF(LAUF00 .LT. (NDIMLEN-8)) THEN      ! op_pj_20170220
      write(*,*)'warning: Incomplete Database for NOx production!'
    ENDIF

    CLOSE(11)

    RETURN

  END SUBROUTINE SPE_ReadNPE

  ! =========================================================================
  SUBROUTINE SPE_IONS(nlev,time_step_len           & ! INPUT 
       ,grheight                                   & ! INPUT
       ,density,densityint                         & ! INPUT
       ,altitude                                   & ! INPUT
       ,ilat, ilon                                 & ! INPUT
       ,philat, philon                             & ! INPUT
       ,jrow, kproma                               & ! INPUT
       ,spe_method                                 & ! INPUT
       ! ka_sv_20160615+
       ,grmass,grvol                               & ! INPUT
       ! ka_sv_20160615-
!ka_sb_20160721+
       ,hn_ion_lat_bin                             & ! INPUT
       ,ikp,nkp,nlat,nlev_a                        & ! INPUT
       ,pmid                                       & ! INPUT
       ,pressure                                   & ! INPUT
       ,ionrate                                    & ! INPUT
!ka_sb_20160721-
       ,ions                                       & ! OUTPUT 
       ,maglat                                     & ! OUTPUT
       ,status                                     & ! OUTPUT
       )


    USE messy_main_constants_mem, ONLY: radius_earth, pi, R_gas & 
         ,M_air & ! mean molar mass of air
         ,N_A ! Avogadro constant [1/mol]

    IMPLICIT NONE

    ! Input Parameters
    INTEGER,                     INTENT(IN) :: nlev     ! number of levels
    REAL(dp),                    INTENT(IN) :: time_step_len ! time step length
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: grheight ! height of box [m]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: density  ! air density midlevel [g cm-3]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: densityint  ! air density interface [g cm-3]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: altitude ! altitude [km]
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: philat ! latitude philat(ilat)
    REAL(dp),    DIMENSION(:),   INTENT(IN) :: philon ! longitude philon(ilon)
    INTEGER,     DIMENSION(:,:), INTENT(IN) :: ilat  ! latitude index (1:kproma,jrow)
    INTEGER,     DIMENSION(:,:), INTENT(IN) :: ilon  ! longitude index (1:kproma,jrow)
    INTEGER,                     INTENT(IN) :: jrow
    INTEGER,                     INTENT(IN) :: kproma 
    INTEGER,                     INTENT(IN) :: spe_method
    ! ka_sv_20160615+
    REAL(dp),    DIMENSION(:,:,:), INTENT(IN) :: grmass, grvol  ! gridbox mass [kg] and volume [m³]
    ! ka_sv_20160615-
!ka_sb_20160721+
    INTEGER,                     INTENT(IN)  :: ikp,nkp,nlat,nlev_a ! necessary for spe_method 3 to select correct kp-bin
    REAL(dp), DIMENSION(:),      INTENT(IN)  :: hn_ion_lat_bin
    REAL(dp), DIMENSION(:,:,:),  INTENT(IN)  :: pmid ! & ! mid-level pressures [Pa]
    REAL(dp), DIMENSION(:),      INTENT(IN)  :: pressure   ! pressure of ionrate in external file
    REAL(dp), DIMENSION(:,:,:), INTENT(IN)   :: ionrate   ! ionrates by H. Nieder (AIMOS)
!ka_sb_20160721-

    ! Output Parameters
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: ions     ! SPE ion pair production [#/cm3/s]
    REAL(DP),  DIMENSION(kproma)  :: maglat     ! geomagnetic latitude [degrees]
    INTEGER,                    INTENT(OUT) :: status   ! error status
    
    ! Local
    REAL(DP) :: nions      ! number of produced ions
    REAL(DP) :: magplat   ! latitude of geomagnetic pole [rad]
    REAL(DP) :: magplon   ! longitude of geomagnetic pole [rad]
    REAL(DP) :: slbth      ! slab thickness [m]
    REAL(DP) :: slbeqth    ! equivalent slab thickness [m]
    ! Integral equivalent slab thickness [m]
    REAL(DP),  DIMENSION(:,:), ALLOCATABLE :: slbeqthInt  
    REAL(DP) :: rhostp     ! density at 15deg, 1Atm (-->Bethe 1951)
    REAL(DP) :: depos      ! deposited energy

    REAL(DP),  DIMENSION(:), ALLOCATABLE  :: alpha, alpha0k

    REAL(DP) :: coslat, coslon, sinlat, sinlon
    REAL(DP) :: th_u, th_l ! thickness at upper and lower interface

    ! matrices for transformation to geomagnetic latitude
    REAL(DP), DIMENSION(3,3) :: geo2mag1, geo2mag2, geo2magf 
    ! vectors to grid point in geographic and magnetic coordinates
    REAL(DP), DIMENSION(3)   :: vgeo, vmag 

!ka_sb_20160721+
    REAL(dp), DIMENSION(nlev) :: ionrate_inter ! ionrates interpolated to EMAC-pressure levels
!ka_sb_20160721-

    INTEGER :: jp,level,perionlev,k,i! counters 
    INTEGER :: ilat_a

    INTRINSIC ABS, ATAN2, COS, MATMUL, RESHAPE, SIN, SQRT, REAL, EXP, SUM

    LOGICAL, DIMENSION(nbin) ::  mask

    ! No error
    status=0

! =========================================================================
    SELECT CASE (spe_method) 

      CASE(0) ! CALCULATE IONRATES INTERNALLY (SPE_PROD)

      !  SUBROUTINE SPE_PROD 
 
     ! ****************************************************************
     ! NOx and HOx production from solar proton events
     ! ****************************************************************
 
     ALLOCATE(alpha(k_max))
     ALLOCATE(alpha0k(k_max-1))
     ALLOCATE(slbeqthInt(kproma,nlev+1))
 
     ! initialize arrays
     ions(:,:,jrow)=0._dp
 
     ! initialize pitch angles
     DO k=1, k_max
        alpha(k)=pi/2._dp/REAL(k_max,dp)*REAL(k-1,dp)
     END DO
     alpha0k=SQRT(alpha(1:k_max-1)*alpha(2:k_max))
 
     rhostp=1E-6_dp *101325._dp*M_air/(R_gas*(15._dp+273.15_dp)) ! [g cm-3]
 
    !************************************************
    ! Calculate weighting by geomagnetic latitude
    !************************************************
    ! initalize
    magplat=78.8_dp/180._dp*pi
    magplon=289.1_dp/180._dp*pi
    ! loop over all rows
    DO jp=1, kproma

       ! GEO to MAG, e.g. Hapgood (1992,1997)
       ! Initialize coordinate transformation variables and matrices
       coslat=COS(philat(ilat(jp,jrow))/180._dp*pi)
       coslon=COS(philon(ilon(jp,jrow))/180._dp*pi)
       sinlat=SIN(philat(ilat(jp,jrow))/180._dp*pi)
       sinlon=SIN(philon(ilon(jp,jrow))/180._dp*pi)
       geo2mag1=RESHAPE( (/ COS(magplat-pi/2._dp), 0._dp, -SIN(magplat-pi/2._dp),  &  ! column 1
            0._dp,             1._dp, 0._dp,                               &  ! column 2
            SIN(magplat-pi/2._dp), 0._dp, COS(magplat-pi/2._dp)/),                 &  ! column 3 
            (/3,3/) )
       geo2mag2=RESHAPE( (/ COS(magplon), -SIN(magplon), 0._dp, &  ! column 1
            SIN(magplon), COS(magplon),  0._dp,                 &  ! column 2
            0._dp        ,       0._dp,  1._dp /),              &  ! column 3 
            (/3,3/) )
       geo2magf=MATMUL(geo2mag1,geo2mag2)
       vgeo(1)=radius_earth*coslat*coslon
       vgeo(2)=radius_earth*coslat*sinlon
       vgeo(3)=radius_earth*sinlat
       vmag=MATMUL(geo2magf,vgeo)
       maglat(jp)=ATAN2(vmag(3),SQRT(vmag(1)**2+vmag(2)**2))/pi*180._dp
       ! maglon(jp)=ACOS(vmag(1)/SQRT(vmag(1)**2+vmag(2)**2))/pi*180._dp ! magn. long. not needed
       ! IF (vmag(2)<0._dp) maglon(jp)=360._dp-maglon(jp) ! magn. long. not needed
    END DO

    ! initialize upper boundary of integral equivalent slab thickness 'slbeqthInt'
    DO jp=1, kproma
       slbeqthInt(jp,1)=0._dp
    END DO 

    level_loopc01:  DO level=1, nlev
       vector_loopc01: DO jp=1, kproma
          ! Calculate slab thickness (density x slab height)
          slbth=densityint(jp,level)*grheight(jp,level) ! [g cm^-3 m]
          ! Calculate slab equivalent thickness (at STP)
          slbeqth=slbth/rhostp ! [m]
          ! Integral equivalent slab thickness at upper interface 
          slbeqthInt(jp,level+1)=slbeqthInt(jp,level)+slbeqth; ! [m]
       END DO vector_loopc01
    END DO level_loopc01

    !**************************************************
    ! Calculate ion density
    !**************************************************
    ! loop over all levels
    level_loopc02:  DO level=1, nlev
       ! loop over all rows
       vector_loopc02: DO jp=1, kproma

          !******************************************
          ! Calculate Deposited energy
          !******************************************
          depos=0._dp
          angle_loop: DO k=1, k_max-1
             th_u = slbeqthInt(jp,level)/COS(alpha0k(k))  ! eq A9 from Vitt and Jackman
             th_l = slbeqthInt(jp,level+1)/COS(alpha0k(k))
             mask = (rangestp>=th_u .AND. rangestp<=th_l) 
             ! SUM over energies E for which holds rangestp(E)>=th_u and rangestp(E)<=th_l
             depos= depos+SUM(pflxspc*energies/(grheight(jp,level)*1E2_dp) &
                     *2._dp*pi*(COS(alpha(k))-COS(alpha(k+1)))*COS(alpha0k(k)),mask)
          END DO angle_loop

          ! Only create ions if geomagnetic latitude is greater r_lat2 (55deg) 
          IF (ABS(maglat(jp))>r_lat2) THEN
             ! Scaling
             IF (ABS(maglat(jp))<r_lat1) &
                  depos=depos*(ABS(maglat(jp))-r_lat2)/(r_lat1-r_lat2)

             ! number of produced ion pairs (35eV for each ionization)
             nions=depos*1E6_dp/35._dp ! #/cm3/s


             ! add to ions stream 
             ions(jp, level,jrow)=nions
          END IF
       END DO vector_loopc02
    END DO  level_loopc02

    DEALLOCATE(alpha)
    DEALLOCATE(alpha0k)
    DEALLOCATE(slbeqthInt)

 ! END SUBROUTINE SPE_PROD
! =========================================================================

 CASE(1) ! EXTERNAL IONIZATION RATES (SPE_EXT)


    ! initialize arrays
    ions(1:kproma,:,jrow)=0._dp

    !************************************************
    ! Calculate weighting by geomagnetic latitude
    !************************************************
    ! initalize
    magplat=78.8_dp/180._dp*pi
    magplon=289.1_dp/180._dp*pi
    ! loop over all rows
    DO jp=1, kproma

       ! GEO to MAG, e.g. Hapgood (1992,1997)
       ! Initialize coordinate transformation variables and matrices
       coslat=COS(philat(ilat(jp,jrow))/180._dp*pi)
       coslon=COS(philon(ilon(jp,jrow))/180._dp*pi)
       sinlat=SIN(philat(ilat(jp,jrow))/180._dp*pi)
       sinlon=SIN(philon(ilon(jp,jrow))/180._dp*pi)
       geo2mag1=RESHAPE( (/ COS(magplat-pi/2._dp), 0._dp, -SIN(magplat-pi/2._dp),  &  ! column 1
            0._dp,             1._dp, 0._dp,                               &  ! column 2
            SIN(magplat-pi/2._dp), 0._dp, COS(magplat-pi/2._dp)/),                 &  ! column 3 
            (/3,3/) )
       geo2mag2=RESHAPE( (/ COS(magplon), -SIN(magplon), 0._dp, &  ! column 1
            SIN(magplon), COS(magplon),  0._dp,                 &  ! column 2
            0._dp        ,       0._dp,  1._dp /),              &  ! column 3 
            (/3,3/) )
       geo2magf=MATMUL(geo2mag1,geo2mag2)
       vgeo(1)=radius_earth*coslat*coslon
       vgeo(2)=radius_earth*coslat*sinlon
       vgeo(3)=radius_earth*sinlat
       vmag=MATMUL(geo2magf,vgeo)
       maglat(jp)=ATAN2(vmag(3),SQRT(vmag(1)**2+vmag(2)**2))/pi*180._dp
       ! maglon(jp)=ACOS(vmag(1)/SQRT(vmag(1)**2+vmag(2)**2))/pi*180._dp ! magn. long. not needed
       ! IF (vmag(2)<0._dp) maglon(jp)=360._dp-maglon(jp) ! magn. long. not needed
    END DO

    !**************************************************
    ! Calculate ion density and N,NO, H, OH production
    !**************************************************
    ! loop over all levels
    level_loopc12:  DO level=1, nlev
       nions=ions_ext(level) ! Externally calculated ion rates
       ! loop over all rows
       vector_loopc12: DO jp=1, kproma

          ! Only create ions if geomagnetic latitude is greater r_lat2 (55deg) 
          IF (ABS(maglat(jp))>r_lat2) THEN
             ! Scaling
             IF (ABS(maglat(jp))<r_lat1) &
                  nions=nions*(ABS(maglat(jp))-r_lat2)/(r_lat1-r_lat2)

             ! add to ions stream 
             ions(jp, level,jrow)=nions
             ! add to xnox stream
          END IF
       END DO vector_loopc12
    END DO  level_loopc12

!  END SUBROUTINE SPE_EXT
  ! =========================================================================

     CASE(2) ! EXTERNAL IONIZATION RATES FROM JAN MAIK WISSING  (SPE_AIMOS)   

      ions(1:kproma,:,jrow)=0._dp
      IF (ASSOCIATED(AIMOS_p)) ions(1:kproma,:,jrow)=ions(1:kproma,:,jrow)+AIMOS_p(1:kproma,:,jrow)*1e-6_dp
      IF (ASSOCIATED(AIMOS_e)) ions(1:kproma,:,jrow)=ions(1:kproma,:,jrow)+AIMOS_e(1:kproma,:,jrow)*1e-6_dp
      IF (ASSOCIATED(AIMOS_a)) ions(1:kproma,:,jrow)=ions(1:kproma,:,jrow)+AIMOS_a(1:kproma,:,jrow)*1e-6_dp
     
     ! END SUBROUTINE SPE_AIMOS
  ! =========================================================================
! ka_sb_20160706+
     CASE(3) ! EXTERNAL IONIZATION RATES FROM HOLGER NIEDER 
     ions(1:kproma,:,jrow)=0._dp

  DO jp=1,kproma
    ! get correct latitude bin
     DO i=1,nlat
        IF (((philat(ilat(jp,jrow)))<=hn_ion_lat_bin(i)) .AND. (philat(ilat(jp,jrow))>hn_ion_lat_bin(i+1))) THEN
          ilat_a=i
          CALL SPE_INTERP_EXT_IONRATE(pmid(jp,:,jrow), kproma,nlev &!INPUT
                                      ,ikp,ilat_a                  &!INPUT
                                      ,nkp,nlat,nlev_a,jrow        &!INPUT
                                      ,pressure                    &!INPUT
                                      ,ionrate                     &!INPUT
                                      ,ionrate_inter)               !OUTPUT

              DO level=1,nlev
                 ions(jp,level,jrow)=ions(jp,level,jrow)+ionrate_inter(level) 
              END DO
        END IF
     END DO
  END DO

! ka_sb_20160706-

  ! =========================================================================

! ka_sv_20160615+
      CASE(4) ! EXTERNAL IONIZATION RATES FROM CMIP6  (SPE_CMIP6)   

      ions(1:kproma,:,jrow)=0._dp
      IF (ASSOCIATED(CMIP6_p)) ions(1:kproma,:,jrow)=ions(1:kproma,:,jrow)+CMIP6_p(1:kproma,:,jrow)
      IF (ASSOCIATED(CMIP6_e)) ions(1:kproma,:,jrow)=ions(1:kproma,:,jrow)+CMIP6_e(1:kproma,:,jrow)
      IF (ASSOCIATED(CMIP6_g)) ions(1:kproma,:,jrow)=ions(1:kproma,:,jrow)+CMIP6_g(1:kproma,:,jrow)
      
      DO level=1,nlev
        DO jp=1,kproma
          ions(jp,level,jrow)=ions(jp,level,jrow)*grmass(jp,level,jrow)*1.e03_dp/grvol(jp,level,jrow)*1.e-6_dp
        END DO
      END DO
     
 ! =========================================================================
! ka_sv_20160615-
!ka_sv_20180508+
      CASE(5)  ! read original AIMOS textfiles
        ions(1:kproma,:,jrow)=0._dp
        ! rest done outside
        
!ka_sv_20180508-


    END SELECT    

 
    END SUBROUTINE SPE_IONS
  ! =========================================================================

    SUBROUTINE SPE_PROD_XNOX  (density                 & ! INPUT
           ,altitude                                   & ! INPUT
           ,nlev                                       & ! INPUT
           ,jrow, kproma                               & ! INPUT
           ,Nperion, NOperion                          & ! INPUT
           ,spe_method                                 & ! INPUT
           ,ions                                       & ! INPUT 
           ,maglat                                     & ! INPUT
           ,xnox, tespen, tespeno                      & ! OUTPUT 
           ,xhox, tespeh, tespeoh                      & ! OUTPUT 
           ,xhno3, tespehno3                           & ! OUTPUT 
!ka_sv_20170426+
           ,n2oprod                                    & ! INPUT
           , xn2o, tespen2o                            & ! OUTPUT
           ,prate_part                                 & !INPUT
!ka_sv_20170426-
           ,status                                     & ! OUTPUT
           )


!!    ! ****************************************************************
!!    ! NOx and HOx production from solar proton events
!!    ! ****************************************************************

    USE messy_main_constants_mem, ONLY: radius_earth, pi, R_gas & 
         ,M_air & ! mean molar mass of air
         ,N_A ! Avogadro constant [1/mol]

    IMPLICIT NONE

    ! Input Parameters
   REAL(dp),    DIMENSION(:,:), INTENT(IN) :: density  ! air density midlevel [g cm-3]
    REAL(dp),    DIMENSION(:,:), INTENT(IN) :: altitude ! altitude [km]
    INTEGER,                     INTENT(IN) :: nlev     ! number of levels
    INTEGER,                     INTENT(IN) :: jrow
    INTEGER,                     INTENT(IN) :: kproma 
    INTEGER,                     INTENT(IN) :: spe_method
    REAL(DP),    DIMENSION(kproma,nlev,jrow), INTENT(IN) :: ions     ! SPE ion pair production [#/cm3/s]
    REAL(DP),    DIMENSION(kproma), INTENT(IN) :: maglat     ! geomagnetic latitude [degrees]
    REAL(DP),    DIMENSION(kproma,nlev), INTENT(IN) :: Nperion, NOperion
!ka_sv_20170426+
    INTEGER, INTENT(IN) :: n2oprod
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: prate_part   ! auroral electrons ionization
!ka_sv_20170426-

    ! Output Parameters
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: xnox     ! SPE NOx production [g/cm3]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: tespen   ! SPE N tendency [mol/mol/s]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: tespeno  ! SPE NO tendency [mol/mol/s]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: xhox     ! SPE HOx  production [g/cm3]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: xh       ! SPE H  production [g/cm3]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: tespeh   ! SPE H tendency [mol/mol/s]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: tespeoh  ! SPE OH tendency [mol/mol/s]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: xhno3    ! SPE HNO3  production [g/cm3]
    REAL(DP),  DIMENSION(:,:,:),    POINTER :: tespehno3   ! SPE HNO3 tendency [mol/mol/s]
!ka_sv_20170426+
    REAL(dp), DIMENSION(:,:,:), POINTER :: xn2o, tespen2o
!ka_sv_20170426-
    INTEGER,                    INTENT(OUT) :: status   ! error status

    ! Local
    REAL(DP) :: N_prod
    REAL(DP) :: NO_prod
    REAL(DP) :: H_prod
    REAL(DP) :: OH_prod
    REAL(DP) :: HNO3_prod

    INTEGER :: jp,level ! counters 

    INTRINSIC ABS, REAL, EXP

    ! No error
    status=0
   
! =========================================================================
    SELECT CASE (spe_method) 

     CASE(0) ! CALCULATE IONRATES INTERNALLY (SPE_PROD)
! =========================================================================
        ! initialize arrays
        xnox(1:kproma,:,jrow)=0._dp
        xhox(1:kproma,:,jrow)=0._dp
        tespen(1:kproma,:,jrow)=0._dp
        tespeno(1:kproma,:,jrow)=0._dp
        tespeoh(1:kproma,:,jrow)=0._dp
        tespehno3(1:kproma,:,jrow)=0._dp
!ka_sv_20170426+
        xn2o(1:kproma,:,jrow)=0._dp
        tespen2o(1:kproma,:,jrow)=0._dp
!ka_sv_20170426-
    
        !******************************************************
        ! Calculate height dependent N/NO production efficiency
        !******************************************************
        !**************************************************
        ! Calculate N,NO, H, OH production
        !**************************************************
        ! loop over all levels
        level_loopc01:  DO level=1, nlev
           ! loop over all rows
           vector_loopc01: DO jp=1, kproma
    
              ! Only create ions if geomagnetic latitude is greater r_lat2 (55deg) 
              IF (ABS(maglat(jp))>r_lat2) THEN
                 ! Scaling
                 ! IF (ABS(maglat(jp))<r_lat1) 
    
                 ! altitude dependent branching ratio and production efficiency
                 N_prod =ions(jp,level,jrow) * Nperion(jp,level)
                 NO_prod=ions(jp,level,jrow) * NOperion(jp,level)
    
                 ! OH molecules per ion pair, ionization rate independent 
                 ! this is a crude approximation 
!ka_sv_20180517+
                 H_prod =max(0._dp,0.5_dp*prate_part(jp, level)*(2._dp-EXP((altitude(jp,level)-83._dp)/6._dp)))     ! prate_part instead of ions (photoionisation probably can't be done the same way as particle ionisation)
               ! will get negative in altitudes above 83 km; avoid negative production by setting it to zero (in agreement with Solomon/Jackman parameterization)
                 OH_prod=H_prod
!ka_sv_20180517-
    
                 ! HNO3 production
                 HNO3_prod=0._dp !0.5_dp*ions(jp,level,jrow)*0.4_dp * (-atan((altitude(jp,level)-48._dp)/2.5_dp)/2.8_dp+0.53_dp)
    
                 ! convert to mole fraction, volume mixing ratio [mol/mol/s]
                 N_prod =(N_prod /N_A) / (density(jp,level)/M_air)
                 NO_prod=(NO_prod/N_A) / (density(jp,level)/M_air)
                 H_prod =(H_prod /N_A) / (density(jp,level)/M_air)
                 OH_prod=(OH_prod/N_A) / (density(jp,level)/M_air)
                 HNO3_prod=(HNO3_prod/N_A) / (density(jp,level)/M_air)
    
                 ! add to xnox stream
                 xnox(jp, level,jrow)=N_prod+NO_prod
                 ! add to xhox stream
                 xhox(jp, level,jrow)=H_prod+OH_prod
                 ! add to hno3 stream
                 xhno3(jp, level,jrow)=HNO3_prod
    
                 ! add tendencies
                 tespen(jp, level, jrow)=N_prod  
                 tespeno(jp, level, jrow)=NO_prod 
                 tespeh(jp, level, jrow)=H_prod 
                 tespeoh(jp, level, jrow)=OH_prod 
                 tespehno3(jp, level, jrow)=HNO3_prod 
              END IF
           END DO vector_loopc01
        END DO  level_loopc01
    
    !  END SUBROUTINE SPE_PROD
  ! =========================================================================

     CASE(1) ! EXTERNAL IONIZATION RATES (SPE_EXT)
  ! =========================================================================
        ! initialize arrays
        xnox(1:kproma,:,jrow)=0._dp
        xhox(1:kproma,:,jrow)=0._dp
        tespen(1:kproma,:,jrow)=0._dp
        tespeno(1:kproma,:,jrow)=0._dp
        tespeoh(1:kproma,:,jrow)=0._dp
!ka_sv_20170426+
        xn2o(1:kproma,:,jrow)=0._dp
        tespen2o(1:kproma,:,jrow)=0._dp
!ka_sv_20170426-
    
        !**************************************************
        ! Calculate ion density and N,NO, H, OH production
        !**************************************************
        ! loop over all levels
        level_loopc11:  DO level=1, nlev
    
           ! loop over all rows
           vector_loopc11: DO jp=1, kproma
    !          ! Only create ions if geomagnetic latitude is greater r_lat2 (55deg) 
              IF (ABS(maglat(jp))>r_lat2) THEN
                 ! Scaling
    !             IF (ABS(maglat(jp))<r_lat1) &
    
                 N_prod=ions(jp,level,jrow)*Nperion(jp,level)
                 NO_prod=ions(jp, level,jrow)*NOperion(jp,level)
                 ! OH molecules per ion pair, ionization rate independent 
                 ! this is a crude approximation 
!ka_sv_20180517+
                 H_prod =max(0._dp,0.5_dp*prate_part(jp, level)*(2._dp-EXP((altitude(jp,level)-83._dp)/6._dp)))     ! prate_part instead of ions (photoionisation probably can't be done the same way as particle ionisation)
               ! will get negative in altitudes above 83 km; avoid negative production by setting it to zero (in agreement with Solomon/Jackman parameterization)
                 OH_prod=H_prod
!ka_sv_20180517-
    
                 ! convert to mole fraction, volume mixing ratio [mol/mol/s]
                 N_prod =(N_prod /N_A) / (density(jp,level)/M_air)
                 NO_prod=(NO_prod/N_A) / (density(jp,level)/M_air)
                 H_prod =(H_prod /N_A) / (density(jp,level)/M_air)
                 OH_prod=(OH_prod/N_A) / (density(jp,level)/M_air)
    
                 ! add to xnox stream
                 xnox(jp, level,jrow)=N_prod+NO_prod
                 ! add to xhox stream
                 xhox(jp, level,jrow)=H_prod+OH_prod
                 ! add tendencies
                 tespen(jp, level, jrow)=N_prod  
                 tespeno(jp, level, jrow)=NO_prod 
                 tespeh(jp, level, jrow)=H_prod 
                 tespeoh(jp, level, jrow)=OH_prod 
              END IF
           END DO vector_loopc11
        END DO  level_loopc11
    
    !  END SUBROUTINE SPE_EXT
  ! =========================================================================
! ka_sv_20160615+
     !CASE(2) ! EXTERNAL IONIZATION RATES FROM JAN MAIK WISSING  (SPE_AIMOS) 
    CASE(2,4) ! EXTERNAL IONIZATION RATES FROM JAN MAIK WISSING OR CMIP6 (SPE_AIMOS) 
! ka_sv_20160615-
  ! =========================================================================
    ! loop over all levels
    level_loopc21:  DO level=1, nlev
       ! loop over all rows
       vector_loopc21: DO jp=1, kproma
 
          ! altitude dependent branching ratio and production efficiency
          N_prod =ions(jp, level, jrow) * Nperion(jp,level)
          NO_prod=ions(jp, level, jrow) * NOperion(jp,level)
          ! OH molecules per ion pair, ionization rate independent 
          ! this is a crude approximation 
!ka_sv_20180517+
          H_prod =max(0._dp,0.5_dp*prate_part(jp, level)*(2._dp-EXP((altitude(jp,level)-83._dp)/6._dp)))     ! prate_part instead of ions (photoionisation probably can't be done the same way as particle ionisation)
               ! will get negative in altitudes above 83 km; avoid negative production by setting it to zero (in agreement with Solomon/Jackman parameterization)
          OH_prod=H_prod
!ka_sv_20180517-
 
          ! convert to mole fraction, volume mixing ratio [mol/mol/s]
          N_prod =(N_prod /N_A) / (density(jp,level)/M_air)
          NO_prod=(NO_prod/N_A) / (density(jp,level)/M_air)

          H_prod =(H_prod /N_A) / (density(jp,level)/M_air)
          OH_prod=(OH_prod/N_A) / (density(jp,level)/M_air)

          ! add to xnox object
          xnox(jp, level,jrow)=N_prod+NO_prod
          ! add to xhox object
          xhox(jp, level,jrow)=H_prod+OH_prod
          ! add tendencies
          tespen(jp, level, jrow)=N_prod  
          tespeno(jp, level, jrow)=NO_prod 
          tespeh(jp, level, jrow)=H_prod 
          tespeoh(jp, level, jrow)=OH_prod 
!ka_sv_20170426+
!TODO: use only electrons
          xn2o(jp,level,jrow)=0._dp
          tespen2o(jp,level,jrow)=0._dp
!ka_sv_20170426-

       END DO vector_loopc21
    END DO  level_loopc21

!ka_sb_20161118+
! ================================================================================
    
    CASE(3,5) ! EXTERNAL IONIZATION RATES FROM Holger Nieder
! =========================================================================
    ! loop over all levels
    level_loopc31:  DO level=1, nlev
       ! loop over all rows
       vector_loopc31: DO jp=1, kproma
 
          ! altitude dependent branching ratio and production efficiency
          N_prod =ions(jp, level, jrow) * Nperion(jp,level)
          NO_prod=ions(jp, level, jrow) * NOperion(jp,level)
 
          ! convert to mole fraction, volume mixing ratio [mol/mol/s]
          N_prod =(N_prod /N_A) / (density(jp,level)/M_air)
          NO_prod=(NO_prod/N_A) / (density(jp,level)/M_air)

!ka_sv_20180517+
!          H_prod =0._dp
!          OH_prod=0._dp
          H_prod =max(0._dp,0.5_dp*prate_part(jp, level)*(2._dp-EXP((altitude(jp,level)-83._dp)/6._dp)))     ! prate_part instead of ions (photoionisation probably can't be done the same way as particle ionisation)
               ! will get negative in altitudes above 83 km; avoid negative production by setting it to zero (in agreement with Solomon/Jackman parameterization)
          OH_prod=H_prod
          H_prod =(H_prod /N_A) / (density(jp,level)/M_air)
          OH_prod=(OH_prod/N_A) / (density(jp,level)/M_air)
!ka_sv_20180517-
         

          ! add to xnox object
          xnox(jp, level,jrow)=N_prod+NO_prod
          ! add tendencies
          tespen(jp, level, jrow)=N_prod  
          tespeno(jp, level, jrow)=NO_prod 
          tespeh(jp, level, jrow)=H_prod 
          tespeoh(jp, level, jrow)=OH_prod 
!ka_sv_20170426+
          xn2o(jp,level,jrow)=prate_part(jp,level)*0.0001_dp
          tespen2o(jp,level,jrow)=prate_part(jp,level)*0.0001_dp
!ka_sv_20170426-

       END DO vector_loopc31
    END DO  level_loopc31
!ka_sv_20170426+
    IF (n2oprod==0) THEN
      xn2o(1:kproma,:,jrow)=0._dp
      tespen2o(1:kproma,:,jrow)=0._dp
    END IF
!ka_sv_20170426-

!  END SUBROUTINE SPE_AIMOS
  ! =========================================================================  
!ka_sb_20161118-

    END SELECT    

END SUBROUTINE SPE_PROD_XNOX
  ! =========================================================================

  ! =========================================================================
!ka_sb_20160707+
  SUBROUTINE spe_read_ionrate(ionrate,pressure,nkp,nlat,nlev_a) 
  
  use netcdf
  USE messy_main_blather, ONLY: start_message, end_message

  IMPLICIT NONE

  !-----
  ! local variables
  !-----
  CHARACTER(len=*), PARAMETER :: substr = 'spe_read_ionrate'
  LOGICAL                     :: lex     ! file existence flag
  INTEGER                     :: fstat   ! file status
  
  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid,varid2
 
  !OUTPUT
   REAL(dp), DIMENSION(13,18,67),  INTENT(OUT)  :: ionrate
   REAL(dp), DIMENSION(67),        INTENT(OUT)  :: pressure

   INTEGER,                          INTENT(OUT)  :: nkp,nlat,nlev_a

  nkp = 13
  nlat= 18
  nlev_a=67
  
  CALL start_message(TRIM(modstr),'READ file with ionization rates from Holger Nieder', substr)
  !-----
  ! initialisation
  !-----

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  call check( nf90_open(input_spe3_file, NF90_NOWRITE, ncid) )
  write(*,*) 'File-name, input_spe3_file: ',input_spe3_file
  write(*,*) 'ncid: ',ncid

  write(*,*) 'Call check 1'

  ! Get the varid of the data variable, based on its name.
  call check( nf90_inq_varid(ncid, "ionrate", varid) )
  write(*,*) 'Call check 2: varid: ',varid
  call check( nf90_inq_varid(ncid, "pressure", varid2) )
  write(*,*) 'Call check 3:varid2: ',varid2

  ! Read the data.
  call check( nf90_get_var(ncid, varid, ionrate ))
  write(*,*) 'Call check 4'
  call check( nf90_get_var(ncid, varid2, pressure ))
  write(*,*) 'Call check 5'

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )
  write(*,*) 'Call check 6: close'

  print *,"*** SUCCESS reading example file ", input_spe3_file, "! "



  CALL end_message(TRIM(modstr),'END READ file with ionization rates from Holger Nieder', substr)

  contains
      subroutine check(status)
        integer, intent ( in) :: status
    
        if(status /= nf90_noerr) then
          print *, trim(nf90_strerror(status))
          stop "Stopped"
        end if
      end subroutine check
    
  END SUBROUTINE spe_read_ionrate
  ! =========================================================================

  ! =========================================================================

  SUBROUTINE HN_ION_PROVIDE_DATA(hn_ion_lat_bin,hn_ion_kp)
      !Calculates different values to separate different kp- and latitudinal 
      !areas of external ionization rates by H. Nieder, only needed for SPE_method=3 
       REAL(dp),    INTENT(OUT) :: hn_ion_lat_bin(19)  !OUTPUT
       REAL(dp),    INTENT(OUT) :: hn_ion_kp(14) !OUTPUT
     
       INTEGER :: i
    
       DO i=1,19
         hn_ion_lat_bin(i)=-i*10._dp+100._dp
       ! 90,80,70,...,-80,-90 (90-80N, 80-70N,...)
       END DO
       hn_ion_kp(1)=-1._dp
       hn_ion_kp(2)=0._dp
       hn_ion_kp(3)=0.5_dp
       DO i=4,14
         hn_ion_kp(i)=i*0.5_dp-1._dp
       END DO
    
    
  END SUBROUTINE HN_ION_PROVIDE_DATA
  ! =========================================================================

  ! =========================================================================

    SUBROUTINE SPE_INTERP_EXT_IONRATE(pmid, kproma,nlev &!INPUT
                                    ,ikp,ilat_a         &!INPUT
                                    ,nkp,nlat,nlev_a,jrow &!INPUT
                                    ,pressure           &!INPUT
                                    ,ionrate            &!INPUT
                                    ,ionrate_inter)      !OUTPUT 
    
      INTEGER,                     INTENT(IN)     :: kproma,nlev,ikp,ilat_a,nkp,nlat,nlev_a,jrow 
      REAL(dp), DIMENSION(nlev),   INTENT(IN)     :: pmid ! pressure needed for interpolation, since file contains pressure levels
    
      REAL(dp), DIMENSION(nlev_a), INTENT(IN)     :: pressure   ! pressure of ionrate in external file
      REAL(dp), DIMENSION(:,:,:),  INTENT(IN)     :: ionrate   ! ionrates by H. Nieder (AIMOS)
      REAL(dp), DIMENSION(nlev), INTENT(OUT)      :: ionrate_inter ! ionrates interpolated to EMAC-pressure levels
    
      !Local
      REAL(dp) :: frac
      INTEGER  :: jp,jk,l
          DO jk=nlev,1,-1
             Do l=1, nlev_a-1
                IF ((pmid(jk) .GE. pressure(l+1)) .AND. (pmid(jk) .LE. pressure(l))) THEN
    
    !  Linear interpolation between different levels
                        frac=(pmid(jk) - pressure(l)) / (pressure(l+1) - pressure(l))
                        ionrate_inter(jk)=ionrate(ikp,ilat_a,l) + frac * ((ionrate(ikp,ilat_a,l+1)-ionrate(ikp,ilat_a,l)))
                EXIT
                ELSEIF (pmid(jk)>pressure(1)) THEN
                    ionrate_inter(jk)=ionrate(ikp,ilat_a,1)
                EXIT
                ELSEIF (pmid(jk)<pressure(nlev_a)) THEN
                    ionrate_inter(jk)=ionrate(ikp,ilat_a,nlev_a)
                EXIT
                END IF
             END DO
          END DO
    
    END	SUBROUTINE SPE_INTERP_EXT_IONRATE
!ka_sb_20160707-

!ka_sb_20161114+
! =========================================================================

! =========================================================================
   !NO photoionization not implemented, needs lyman alpha

  SUBROUTINE spe_phioniz(COSSZA,RO2p,ROp4S,ROp2D,ROp2P,RN2p,RNp,RN, &
      ZO2,ZO3P,ZO1D,ZN2,concm,col_O2_box,col_O_box,col_N2_box,press,temperature,LNT,flux107_data,altitude_orig)

!  subroutine to calculate the production rates of O2+ ions, O+, N2+ and N+ out of photoionization.
!  calculates the slant column of O2 and N2 and the resulting ionization rate (resp. products) using the solar zenith angle.

    USE messy_main_constants_mem, ONLY: g, N_A, R_gas, M_O2, MO, M_N2
    implicit none

      INTEGER, PARAMETER :: specnum=17   !number of bins in reduced spectrum - multi-number-entries counted only once

      REAL(DP), DIMENSION(specnum) :: HFG
      REAL(DP), DIMENSION(specnum) :: sigO2, sigO, sigN2 !cross-sections depend on wavelength  
      REAL(DP), DIMENSION(specnum) :: c1,c2,brO2O2p,brO2Op
      REAL(DP), DIMENSION(specnum) :: brOOp4S,brOOp2D,brN2N2p
      REAL(DP), DIMENSION(specnum) :: brN2Np,brN2dis,brO2dis  
      REAL(DP), DIMENSION(specnum) :: brOOp2P
      REAL(DP), DIMENSION(specnum) :: seOOp4S
      REAL(DP), DIMENSION(specnum) :: seOOp2D, seOOp2P, seO2O2p, seO2Op, seO2dis 
      REAL(DP), DIMENSION(specnum) :: seN2N2p,seN2Np,seN2dis 
      
      INTEGER :: lauf
      INTEGER, INTENT(IN):: LNT !(nlev)

      INTEGER :: level 

      REAL(DP), INTENT(IN) :: COSSZA
      REAL(DP), DIMENSION(LNT), INTENT(OUT) :: RO2p,ROp4S,ROp2D,ROp2P,RN2p,RNp,RN
      REAL(DP) :: colO2, colO, colN2 !partial column above respective level
      REAL(DP), DIMENSION(0:LNT), INTENT(IN)  :: press 
      REAL(DP), DIMENSION(LNT), INTENT(IN)  :: ZO2,ZO3P,ZO1D,ZN2,concm !vmr (nlev)
      REAL(DP), DIMENSION(LNT) :: ZO
      REAL(DP), DIMENSION(LNT), INTENT(IN) :: col_O2_box,col_O_box,col_N2_box
      REAL(DP), DIMENSION(2), INTENT(IN) :: flux107_data
      REAL(DP) :: F107,MWF107
      REAL(DP), DIMENSION(LNT), INTENT(IN)::altitude_orig !(grheight)
      REAL(DP), DIMENSION(LNT) :: temperature,altitude
      REAL(DP) :: r1,r2, fref, Ntotal
      
      INTEGER :: rl,wl

      F107=flux107_data(1)
      MWF107=flux107_data(2)

      HFG  = (/5.01E+001,1.00E+004,2.00E+006,7.60E+006,1.66E+008, &
      4.01E+008,2.08E+009,1.72E+009,6.79E+009,2.75E+009,5.04E+009, &
      1.56E+009,3.01E+009,5.44E+008,6.06E+009,5.57E+009,6.31E+009/)
      sigO2 = (/0.0045e-18,0.034e-18,0.2251e-18,0.2101e-18, &
      0.646e-18,2.6319e-18,7.6283e-18,13.2125e-18,16.8233e-18, &
      20.3066e-18,27.0314e-18,23.5669e-18,10.498e-18,13.395e-18, &
      18.7145e-18,1.632e-18,1.15e-18/)
      sigO = (/0.0023e-18,0.017e-18,0.1125e-18,0.105e-18, &
      0.3247e-18,1.319e-18,3.7832e-18,6.0239e-18,7.7205e-18, &
      10.7175e-18,13.1253e-18,8.5159e-18,3.0031e-18,0.0e-18, &
      0.0e-18,0.0e-18,0.0/)
      sigN2 = (/0.0025e-18,0.0201e-18,0.1409e-18,1.137e-18, &
      0.3459e-18,1.5273e-18,5.0859e-18,9.9375e-18,11.7383e-18, &
      19.6514e-18,23.0931e-18,23.0346e-18,2.1434e-18,2.1775e-18, &
      2.5465e-18,0.0e-18,0.0/)
      c1 = (/0.00E+000,0.00E+000,0.00E+000,7.47E+005,6.62E+007, &
       1.66E+008,1.51E+008,3.31E+008,2.22E+009,5.47E+008,2.97E+009, & 
       6.94E+008,3.01E+009,4.19E+008,3.67E+009,4.98E+009,5.80E+009/)
      c2 = (/2.95E+002,7.60E+003,4.60E+005,9.22E+005,4.29E+006, & 
      5.68E+006,6.27E+007,9.83E+007,4.29E+007,1.08E+008,1.59E+007, & 
      8.21E+006,0.00E+000,0.00E+000,0.00E+000,0.00E+000,0.00E+000/)
      brO2O2p = (/0.0,0.0,0.0,0.0,0.108,0.347,0.553,0.624,0.649,0.759, &
      0.874,0.672,0.549,0.756,0.83,0.613,0.0/)
      brO2Op = (/1.0,1.0,1.0,1.0,0.892,0.653,0.447,0.376,0.351,0.24, &
      0.108,0.001,0.0,0.0,0.0,0.0,0.0/)
      brOOp4S = (/0.39,0.39,3.90E-001,3.90E-001,3.93E-001,3.89E-001, &
      3.67E-001,3.50E-001,3.46E-001,3.17E-001,2.98E-001,6.55E-001, &
      1.00E+000,0.00E+000,0.0,0.0,0.0/)
      brOOp2D = (/0.378,0.378,0.378,0.378,0.374,0.377,0.392,0.402, &
      0.403,0.424,0.451,0.337,0.0,0.0,0.0,0.0,0.0/)
      brOOp2P = (/0.224,0.224,0.224,0.224,0.226,0.227,0.233,0.241, &
      0.246,0.26,0.252,0.009,0.0,0.0,0.0,0.0,0.0/)
      brN2N2p = (/0.04,0.04,0.04,0.04,0.717,0.751,0.747,0.754,0.908, &
      0.996,1.0,0.679,0.0,0.0,0.0,0.0,0.0/)
      brN2Np = (/0.96,0.96,0.96,0.96,0.282,0.249,0.253,0.246,0.093, &
      0.005,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
      brN2dis = (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.32, &
      1.0,1.0,1.0,0.0,0.0/)
      brO2dis = (/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.017,0.327, &
      0.451,0.244,0.17,0.387,1.0/)

      seOOp4S = (/81.24,18.896,9.425,28.622,2.019,0.902,0.47,0.325, &
      0.209,0.084,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
      seOOp2D = (/88.526,20.691,9.365,28.199,1.962,0.853,0.418,0.253, &
      0.148,0.034,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
      seOOp2P = (/47.358,11.007,4.772,14.556,1.014,0.436,0.203,0.116, &
      0.061,0.009,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
      seO2O2p = (/134.69,32.212,13.309,39.615,2.834,1.092,0.416,0.189, &
      0.09,0.023,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
      seO2Op = (/76.136,17.944,6.981,20.338,1.437,0.521,0.163,0.052, &
      0.014,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
      seO2dis = (/87.864,20.318,17.821,56.969,4.113,2.041,1.271,0.996, &
      0.762,0.653,0.011,0.0,0.0,0.0,0.0,0.0,0.0/)
      seN2N2p = (/263.99,62.57,25.213,8.54,6.142,2.288,0.786,0.324, &
      0.169,0.031,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
      seN2Np  = (/78.674,18.31,6.948,2.295,1.647,0.571,0.146,0.037, &
      0.008,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
      seN2dis  = (/245.0,52.052,25.255,9.049,6.532,2.909,1.371,0.764, &
      0.515,0.157,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)

     level = LNT

      do lauf=1,LNT
        ZO(lauf)=ZO3P(lauf)+ZO1D(lauf)
      enddo

!     Important: altitude in cm to use matching units!
    DO lauf=1,LNT
        altitude(lauf) = altitude_orig(lauf)*100.0
    ENDDO

      r1 = 0.0138*(F107 - 71.5) + 0.005*(F107 - MWF107 + 3.9)
      r2 = 0.5943*(F107 - 71.5) + 0.381*(F107 - MWF107 + 3.9)

      do rl = 1,level
       RO2p(rl)  = 0.0
       ROp4S(rl) = 0.0
       ROp2D(rl) = 0.0
       ROp2P(rl) = 0.0
       RN2p(rl)  = 0.0_dp
       RNp(rl)   = 0.0
       RN(rl)    = 0.0
      enddo


      if(COSSZA .GT. 0.0_dp) THEN

!    calculate separately for each wavelength interval
        do wl = 1,specnum
          fref = HFG(wl) + r1*c1(wl) + r2*c2(wl)
          fref = MAX(fref,0.0_dp)

        do rl = 1,level-1
            colO2=col_O2_box(rl)
            colO=col_O_box(rl)
            colN2=col_N2_box(rl)
        
!c     FIXME: using plane atmosphere approximation - satisfying?
!ccc     Beer-Lambert; here: absorption before light meets respective layer
        fref = fref*exp((-sigO2(wl)*colO2 - sigO(wl)*colO &
                         -sigN2(wl)*colN2)/COSSZA) 
!ccc     Beer-Lambert: calculate total absorption in respective layer (per cm2), derive mean reaction rates! split up linearly in absorption parts ...

!ccc     first new colX (extinction in the box)
colO2=col_O2_box(rl+1)
colO=col_O_box(rl+1)
colN2=col_N2_box(rl+1)

        Ntotal = 1.0 - exp((-sigO2(wl)*colO2 - sigO(wl)*colO &
                            -sigN2(wl)*colN2)/COSSZA)

!     Rate for production processes from total photon absorption; sums up for all wavelength intervals
        RO2p(rl) = RO2p(rl) + 1.0/(altitude(rl))*  &
        fref*COSSZA*Ntotal*brO2O2p(wl)*(1.0+seO2O2p(wl))*             &
        (1.0-exp((-sigO2(wl)*colO2)/COSSZA))/                         &
        ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                        &
        1.0-exp((-sigO(wl)*colO)/COSSZA) +                            &
        1.0-exp((-sigN2(wl)*colN2)/COSSZA) )

!    Assumption: hv + O can produce 4S, 2D und 2P, D.I. of O2 only O+_4S 
        ROp4S(rl)  = ROp4S(rl)  + 1.0/(altitude(rl))*&
        fref*COSSZA*Ntotal*brO2O2p(wl)*(1.0+seO2O2p(wl))*             &
        (1.0-exp((-sigO2(wl)*colO2)/COSSZA))/                         &
        ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                        &
        1.0-exp((-sigO(wl)*colO)/COSSZA) +                            &
        1.0-exp((-sigN2(wl)*colN2)/COSSZA) )

        ROp4S(rl)  = ROp4S(rl)  + 1.0/(altitude(rl))*&
        fref*COSSZA*Ntotal*brOOp4S(wl)*(1.0+seOOp4S(wl))*             &
        (1.0-exp((-sigO(wl)*colO)/COSSZA))/                           &
        ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                        &
        1.0-exp((-sigO(wl)*colO)/COSSZA) +                            &
        1.0-exp((-sigN2(wl)*colN2)/COSSZA) )

        ROp2D(rl)  = ROp2D(rl)  + 1.0/(altitude(rl))*&
        fref*COSSZA*Ntotal*brOOp2D(wl)*(1.0+seOOp2D(wl))*             &
        (1.0-exp((-sigO(wl)*colO)/COSSZA))/                           &
        ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                        &
        1.0-exp((-sigO(wl)*colO)/COSSZA) +                            &
        1.0-exp((-sigN2(wl)*colN2)/COSSZA) )

        ROp2P(rl)  = ROp2P(rl)  + 1.0/(altitude(rl))*&
        fref*COSSZA*Ntotal*brOOp2P(wl)*(1.0+seOOp2P(wl))*             &
        (1.0-exp((-sigO(wl)*colO)/COSSZA))/                           &
        ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                        &
        1.0-exp((-sigO(wl)*colO)/COSSZA) +                            &
        1.0-exp((-sigN2(wl)*colN2)/COSSZA) ) 

         RN2p(rl) = RN2p(rl) + 1.0/(altitude(rl))*                    &
         fref*COSSZA*Ntotal*brN2N2p(wl)*(1.0+seN2N2p(wl))*            &
         (1.0-exp((-sigN2(wl)*colN2)/COSSZA))/                        &
         ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                       &
         1.0-exp((-sigO(wl)*colO)/COSSZA) +                           &
         1.0-exp((-sigN2(wl)*colN2)/COSSZA) )

         RNp(rl)  = RNp(rl)  + 1.0/(altitude(rl))*                    &
         fref*COSSZA*Ntotal*brN2Np(wl)*(1.0+seN2Np(wl))*              &
         (1.0-exp((-sigN2(wl)*colN2)/COSSZA))/                        &
         ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                       &
         1.0-exp((-sigO(wl)*colO)/COSSZA) +                           &
         1.0-exp((-sigN2(wl)*colN2)/COSSZA) ) 

         RN(rl)  = RN(rl)  + 1.0/(altitude(rl))*                      &
         fref*COSSZA*Ntotal*brN2Np(wl)*(1.0+seN2Np(wl))*              &
         (1.0-exp((-sigN2(wl)*colN2)/COSSZA))/                        &
         ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                       &
         1.0-exp((-sigO(wl)*colO)/COSSZA) +                           &
         1.0-exp((-sigN2(wl)*colN2)/COSSZA))

         RN(rl)  = RN(rl)  + 1.0/(altitude(rl))*                      &
         fref*COSSZA*Ntotal*brN2dis(wl)*(1.0+seN2dis(wl))*            &
         (1.0-exp((-sigN2(wl)*colN2)/COSSZA))/                        &
         ( 1.0-exp((-sigO2(wl)*colO2)/COSSZA) +                       &
         1.0-exp((-sigO(wl)*colO)/COSSZA) +                           &
         1.0-exp((-sigN2(wl)*colN2)/COSSZA))

       enddo !rl
      enddo !wl
           
           RO2p(LNT)= 0.0
           ROp4S(LNT)= 0.0
           ROp2D(LNT)= 0.0
           ROp2P(LNT)= 0.0
           RN2p(LNT)= 0.0_dp
           RNp(LNT)= 0.0
           RN(LNT)= 0.0

      endif !COSSZA gt 0.0

      DO rl = 1,level
       IF(                          &
           RO2p(rl)  .LT. 0.0 .OR.  &
           ROp4S(rl) .LT. 0.0 .OR.  &
           ROp2D(rl) .LT. 0.0 .OR.  &
           ROp2P(rl) .LT. 0.0 .OR.  &
           RN2p(rl)  .LT. 0.0 .OR.  &
           RNp(rl)   .LT. 0.0 .OR.  &
           RN(rl)    .LT. 0.0       &
           ) THEN

       write(*,*)'phioniz rate negative'
       write(*,*) RO2p(rl) , ROp4S(rl) , &
                  ROp2D(rl) , ROp2P(rl) , RN2p(rl) , RNp(rl) , RN(rl)
       write(*,*) 'r1,r2',r1,r2

       ENDIF

!    avoid rates beeing lt 0 - should not be necessary with fref GE 0
       RO2p(rl) = MAX(RO2p(rl),0.0)
       ROp4S(rl)= MAX(ROp4S(rl),0.0)
       ROp2D(rl)= MAX(ROp2D(rl),0.0)
       ROp2P(rl)= MAX(ROp2P(rl),0.0)
       RN2p(rl) = MAX(RN2p(rl),0.0_dp)
       RNp(rl)  = MAX(RNp(rl),0.0)
       RN(rl)   = MAX(RN(rl),0.0)

      enddo !rl


      END SUBROUTINE spe_phioniz

!======================================================================
!ka_sb_20161114-

  SUBROUTINE spe_read_aimos_orig(ionrate_aimos_orig,spe_aimos_dir,spe_aimos_prefix  & 
               ,year,DAYOFYEAR,spe_nlat,spe_nlon,spe_nlev,aimos_time,aimpro,aimele,aimalp,iou)
   INTEGER, INTENT(in) :: spe_nlat, spe_nlon, spe_nlev, aimos_time, year, DAYOFYEAR
   LOGICAL, INTENT(in) :: aimpro, aimele, aimalp
   CHARACTER(LEN=200), INTENT(in) :: spe_aimos_dir, spe_aimos_prefix
   INTEGER, INTENT(in) :: iou
   REAL(dp), DIMENSION(spe_nlon, spe_nlev, spe_nlat, aimos_time, 3), INTENT(out) :: ionrate_aimos_orig
   ! LOCAL
   CHARACTER(LEN=200) :: aimos_filename
   CHARACTER(LEN=4) :: stryear
   CHARACTER(LEN=3) :: strDAYOFYEAR
   INTEGER :: io_error
   CHARACTER :: chdummy
   INTEGER :: idummy
   REAL(dp) ::rdummy
   INTEGER :: i,j,k,l,m
   REAL(dp), DIMENSION(spe_nlev) :: ionrate_dummy
   
   
   write(stryear,'(i4)') year
   write(strDAYOFYEAR,'(i3)') DAYOFYEAR
   if (DAYOFYEAR .le. 99) then
     write(strDAYOFYEAR,'(i2)') DAYOFYEAR
     if (DAYOFYEAR .le. 9) then
       write(strDAYOFYEAR,'(i1)') DAYOFYEAR
     end if
   end if
   aimos_filename=TRIM(spe_aimos_dir)//TRIM(spe_aimos_prefix)//TRIM(stryear)//'_doy'//TRIM(strDAYOFYEAR)//'.dat'
   write(*,*) 'SPE: Open: ',aimos_filename
 
   open(unit=iou,file=TRIM(aimos_filename),iostat=io_error)
   if (io_error == 0) then
     do i=1,aimos_time
       do j=1,3
         do k=1,spe_nlat
           do l=1,spe_nlon
             read(iou,*) idummy,chdummy, rdummy, rdummy, rdummy, rdummy, ionrate_dummy   ! file is not strictly formatted
             do m=1,spe_nlev
               ionrate_aimos_orig(l,spe_nlev+1-m,k,i,j)=ionrate_dummy(m)/1.e06    ! change vertical direction and convert from m3 to cm3
             enddo
           enddo
         enddo
       enddo
     enddo
   else
     write(*,*) 'ERROR: Opening AIMOS file: ',TRIM(aimos_filename),' ERROR:',io_error
     stop
   end if
   close(iou)
   IF (.NOT. aimpro) ionrate_aimos_orig(:,:,:,:,1)=0._dp
   IF (.NOT. aimele) ionrate_aimos_orig(:,:,:,:,2)=0._dp
   IF (.NOT. aimalp) ionrate_aimos_orig(:,:,:,:,3)=0._dp
  
  END SUBROUTINE spe_read_aimos_orig

! ***********************************************************************
END MODULE messy_spe
! ***********************************************************************
