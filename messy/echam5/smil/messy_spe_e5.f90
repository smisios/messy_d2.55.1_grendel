! ***********************************************************************
MODULE messy_spe_e5
! ***********************************************************************

  ! MESSy-SMIL FOR SUBMODEL SPE
  !
  ! SPE = SOLAR PROTON EVENT PARAMETERIZATION
  !
  ! Authors: Andreas Baumgaertner, MPICH, Jan 2007;
  !          Stefan Versick, KIT-SCC & KIT-IMK-ASF, Jul 2016;
  !          Sabine Barthlott, KIT-IMK-ASF, Feb 2016
  ! 
  ! Jun 2014: Corrected code for AIMOS-rates, Stefan Versick (KIT)
  ! Feb 2016: Included N/NO production per ionpair based on 
  !           Nieder et al, Sabine Barthlott (KIT)
  ! Jun 2016: Included CMIP6 ionization, Stefan Versick (KIT)
  ! Jul 2016: Included AIMOS ionization rates by Holger Nieder, 
  !           Sabine Barthlott (KIT)
  ! Nov 2016: Included Photoionization, Sabine Barthlott (KIT)
  ! Dec 2017: Included possibility to use positive ions, Stefan Versick (KIT)


  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , info_bi, error_bi
  ! MESSy
  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL &
                                    , t_chaobj_cpl
  USE messy_spe

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! POINTERS TO DIAGNOSTIC CHANNEL OBJECTS:
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: ions => NULL()      ! SPE ion pair production [#/cm3/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: xnox => NULL()      ! SPE NOx production [g/cm3]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: xhox => NULL()      ! SPE HOx production [g/cm3]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: xhno3 => NULL()     ! SPE HNO3 production [g/cm3]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespen => NULL()    ! SPE N  prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespeno => NULL()   ! SPE NO prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespeh => NULL()    ! SPE H  prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespeoh => NULL()   ! SPE OH prod. tend. [mol/mol/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: tespehno3 => NULL() ! SPE OH prod. tend. [mol/mol/s]
!ka_sv_20170426+
  REAL(dp), DIMENSION(:,:,:), POINTER :: xn2o => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: tespen2o => NULL()
!ka_sv_20170426-
! ka_sv_20180518+
  REAL(dp), DIMENSION(:,:,:), POINTER :: tespeh2o => NULL()     ! SPE Destruction of H2O -> destroyed when HOx is produced
! ka_sv_20180518-
  !ionization rates due to photoionization and particles [ions/cm3/s]
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: PRate_phi_out  => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: PRate_part_out => NULL() 
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: PRate_phi_out_orig =>NULL() 
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: RO2p_out => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: ROp4S_out => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: ROp2D_out => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: ROp2P_out => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: RN2p_out =>NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: RNp_out => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: RN_out => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: RN2D_out => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER  :: RNOp_out => NULL()  !separate ionization rates due to photoionization ions/cm3/s

  ! CPL NAMELIST
  TYPE(t_chaobj_cpl) :: spe_data_int, spe_data_ion
  TYPE(t_chaobj_cpl) :: spe_kp                   !kit_sb_20160705
  TYPE(t_chaobj_cpl) :: spe_flux107, spe_cossza  !kit_sb_20161109

  ! GLOBAL PARMETERS
  ! TRACERS 
  INTEGER                            :: nlntrac_n    ! no of N  tracers
  INTEGER                            :: nlntrac_no   ! no of NO tracers
  INTEGER                            :: nlntrac_h    ! no of H  tracers
  INTEGER                            :: nlntrac_oh   ! no of OH tracers
  INTEGER                            :: nlntrac_hno3 ! no of HNO3 tracers
  INTEGER, DIMENSION(:), POINTER     :: idt_list_n
  INTEGER, DIMENSION(:), POINTER     :: idt_list_no
  INTEGER, DIMENSION(:), POINTER     :: idt_list_h
  INTEGER, DIMENSION(:), POINTER     :: idt_list_oh
  INTEGER, DIMENSION(:), POINTER     :: idt_list_hno3

  INTEGER :: idt_O3P, idt_O1D, idt_H2O, idt_N2, idt_O2 !kit_sb_20160303
!ka_sv_20170426+
  INTEGER :: idt_n2o
!ka_sv_20170426-
!ka_sv_20171213+
  INTEGER :: idt_n2d, idt_n4s, idt_n2p, idt_np, idt_o2p, idt_op, idt_nop, idt_em, idt_o4sp,idt_o2dp, idt_o2pp
!ka_sv_20171213-
  REAL(dp), DIMENSION(:), POINTER :: kp_data => NULL() !kit_sb_20160303
  !ka_sb_20161109+
  REAL(dp), DIMENSION(:), POINTER :: flux107_data => NULL() !contains F107 and MWF107 (running mean 81days)
  REAL(dp), DIMENSION(:,:), POINTER :: cossza_data => NULL()
  !ka_sb_20161109-

  ! GLOBAL PARMETERS
  CHARACTER(LEN=STRLEN_CHANNEL), PUBLIC :: AIMOS_channel = 'import_grid'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: AIMOS_p_object = 'jmionrates_p_ionrate_p'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: AIMOS_e_object = 'jmionrates_e_ionrate_e'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: AIMOS_a_object = 'jmionrates_a_ionrate_a'
  
  CHARACTER(LEN=STRLEN_CHANNEL), PUBLIC :: CMIP6_channel = 'import_grid'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: CMIP6_p_object = 'CMIP6_protons_iprp'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: CMIP6_e_object = 'CMIP6_electrons_iprm'
  CHARACTER(LEN=STRLEN_OBJECT),  PUBLIC :: CMIP6_g_object = 'CMIP6_gcr_iprg'
 
  !ka_sb_20160720+
  REAL(dp), DIMENSION(13,18,67) :: ionrate
  REAL(dp), DIMENSION(67)       :: pressure
  REAL(dp), DIMENSION(19)       :: hn_ion_lat_bin
  REAL(dp), DIMENSION(14)       :: hn_ion_kp
  INTEGER                       :: ikp,nkp,nlat,i,nlev_a
  !ka_sb_20160720- 
! op_pj_20170220+
!!$  !kit_sb_20160303+
!!$  REAL(dp), DIMENSION(1548288)  :: SAP_Np  !DIONP-CONCM-TEMP-H2O-O-N-NO-NOpf-Nf-NMf-OMf
!!$  REAL(dp), DIMENSION(1548288)  :: SAP_NOp !DIONP-CONCM-TEMP-H2O-O-N-NO-NOpf-Nf-NMf-OMf
!!$  !kit_sb_20160303-
! op_pj_20170220-

  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: ionrate_aimos_orig
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: ionrate_aimos_used => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER, PUBLIC :: ionrate_aimos_used_dec => NULL()   ! decomposed field

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: spe_initialize
  PUBLIC :: spe_init_memory
  PUBLIC :: spe_init_coupling
  PUBLIC :: spe_global_start
  PUBLIC :: spe_physc
  PUBLIC :: spe_free_memory
  !PRIVATE :: spe_read_nml_cpl
  !PRIVATE :: spe_read_data

CONTAINS

! ========================================================================
  SUBROUTINE spe_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL spe_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CTRL NAMELIST', substr)
    END IF

    CALL p_bcast(spe_method,   p_io)
    CALL p_bcast(r_lat1,       p_io)
    CALL p_bcast(r_lat2,       p_io)
    CALL p_bcast(nbin,         p_io)
    CALL p_bcast(minoffs,      p_io)
    CALL p_bcast(spect_interp_method,p_io)
    CALL p_bcast(rangerela,p_io)
    CALL p_bcast(rangerelb,p_io)
    CALL p_bcast(ion_km,p_io)
    CALL p_bcast(Nperion_km,p_io)
    CALL p_bcast(NOperion_km,p_io)
    !kit_sb_20160923+
    CALL p_bcast(npe_method,p_io)
    CALL p_bcast(input_sap_file,p_io)
    !kit_sb_20160923-
    !kit_sb_20161116+
    CALL p_bcast(switchphi,p_io)
    !kit_sb_20161116-
!ka_sv_20170426+
    CALL p_bcast(n2oprod,p_io)
!ka_sv_20170426-

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL spe_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('ERROR IN CPL NAMELIST', substr)
    END IF
    CALL p_bcast(spe_data_int%cha, p_io)
    CALL p_bcast(spe_data_int%obj, p_io)
    CALL p_bcast(spe_data_ion%cha, p_io)
    CALL p_bcast(spe_data_ion%obj, p_io)
    CALL p_bcast(AIMOS_channel, p_io)
    CALL p_bcast(AIMOS_p_object, p_io)
    CALL p_bcast(AIMOS_e_object, p_io)
    CALL p_bcast(AIMOS_a_object, p_io)
    ! ka_sb_20160705+
    CALL p_bcast(spe_kp%cha, p_io)
    CALL p_bcast(spe_kp%obj, p_io) 
    CALL p_bcast(input_spe3_file, p_io)
    ! ka_sb_20160705-
    ! ka_sb_20161109+
    CALL p_bcast(spe_flux107%cha, p_io)
    CALL p_bcast(spe_flux107%obj, p_io)
    CALL p_bcast(spe_cossza%cha, p_io)
    CALL p_bcast(spe_cossza%obj, p_io)
    ! ka_sb_20161109-
    ! ka_sv_20160615+
    CALL p_bcast(CMIP6_channel, p_io)
    CALL p_bcast(CMIP6_p_object, p_io)
    CALL p_bcast(CMIP6_e_object, p_io)
    CALL p_bcast(CMIP6_g_object, p_io)
    ! ka_sv_20160615-
    ! ka_sv_20180514+
    CALL p_bcast(spe_aimos_dir, p_io)
    CALL p_bcast(spe_aimos_prefix, p_io)
    CALL p_bcast(spe_nlat, p_io)
    CALL p_bcast(spe_nlon, p_io)
    CALL p_bcast(spe_nlev, p_io)
    CALL p_bcast(aimos_time, p_io)
    CALL p_bcast(aimpro, p_io)
    CALL p_bcast(aimele, p_io)
    CALL p_bcast(aimalp, p_io)
    CALL p_bcast(spe_hox, p_io)
    ! ka_sv_20180514-

    !ka_sb_20160707+
    SELECT CASE(spe_method) 
    CASE(0,1,2,4)
       write(*,*) 'spe_method: ',spe_method    
    CASE(3) ! External ionrates by H. Nieder
       write(*,*) 'spe_method: ',spe_method, 'AIMOS by Holger Nieder, WARNING: better only use with upper atmosphere EMAC'
       CALL SPE_READ_IONRATE(ionrate,pressure,nkp,nlat,nlev_a) !reads nc-file containing AIMOS ionrates by H.Nieder
       CALL HN_ION_PROVIDE_DATA(hn_ion_lat_bin,hn_ion_kp) !determines bin-limits (necessary for latitude & kp-values)

       ! ka_sb_20161010+
       CALL p_bcast(ionrate, p_io)
       CALL p_bcast(pressure, p_io)
       CALL p_bcast(nkp, p_io)
       CALL p_bcast(nlat, p_io) 
       CALL p_bcast(nlev_a, p_io)  
       CALL p_bcast(hn_ion_lat_bin, p_io)
       CALL p_bcast(hn_ion_kp, p_io)
       CALL p_bcast(ikp, p_io)
       ! ka_sb_20161010-
       !ka_sv_20171213+
       CALL p_bcast(calc_pos_ions, p_io)
       !ka_sv_20171213-
       !ka_sv_20180327+
       CALL p_bcast(calc_three_o_ions, p_io)
       CALL p_bcast(br_n4s, p_io)
       CALL p_bcast(br_n2d, p_io)
       !ka_sv_20180327-
    CASE(5)
       write(*,*) 'spe_method: ',spe_method, 'AIMOS by Jan-Maik Wissing'
       CALL p_bcast(calc_pos_ions, p_io)
       CALL p_bcast(calc_three_o_ions, p_io)
       CALL p_bcast(br_n4s, p_io)
       CALL p_bcast(br_n2d, p_io)
    END SELECT
    
    IF (npe_method==2) THEN
! op_pj_20170220+
!!$    CALL SPE_ReadNPE(SAP_Np,SAP_NOp)
       CALL SPE_ReadNPE
! op_pj_20170220-
    END IF
    !ka_sb_20160707-

  END SUBROUTINE spe_initialize
! ========================================================================

! ========================================================================
  SUBROUTINE spe_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    USE messy_main_channel_bi,        ONLY: GP_3D_MID
    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_init_memory'
    INTEGER :: status

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'ions', p3=ions)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'ions' &
         , 'long_name', c='SPE ion pair production')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'ions' &
         , 'units', c='ions/cm3/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'xnox', p3=xnox)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'xnox' &
         , 'long_name', c='SPE NOx production')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xnox' &
         , 'units', c='g/cm3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'xhox', p3=xhox)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'xhox' &
         , 'long_name', c='SPE HOx production')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xhox' &
         , 'units', c='g/cm3')
    CALL channel_halt(substr, status)

!ka_sv_20170426+
    CALL new_channel_object(status, modstr,'xn2o', p3=xn2o)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr,'xn2o' &
         , 'long_name',c='SPE N2O production')
    CALL channel_halt(substr,status)
    CALL new_attribute(status,modstr,'xn2o' &
         , 'units', c='g/cm3')
    CALL channel_halt(substr,status)
!ka_sv_20170426-

    CALL new_channel_object(status, modstr,'xhno3', p3=xhno3)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'xhno3' &
         , 'long_name', c='SPE HNO3 production')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'xhno3' &
         , 'units', c='g/cm3')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespen', p3=tespen)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespen' &
         , 'long_name', c='SPE N production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespen' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespeno', p3=tespeno)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespeno' &
         , 'long_name', c='SPE NO production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespeno' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespeh', p3=tespeh)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespeh' &
         , 'long_name', c='SPE NH production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespeh' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespeoh', p3=tespeoh)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespeoh' &
         , 'long_name', c='SPE OH production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespeoh' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,'tespehno3', p3=tespehno3)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespehno3' &
         , 'long_name', c='SPE HNO3 production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespehno3' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)

!ka_sv_20170426+
    CALL new_channel_object(status, modstr,'tespen2o', p3=tespen2o)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespen2o' &
         , 'long_name', c='EEP N2O production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespen2o' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)
!ka_sv_20170426-

!ka_sv_20180518+
    CALL new_channel_object(status, modstr,'tespeh2o', p3=tespeh2o)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr ,'tespeh2o' &
         , 'long_name', c='EEP H2O production tendency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'tespeh2o' &
         , 'units', c='mol/mol/s')
    CALL channel_halt(substr, status)
!ka_sv_20180518-

    !ka_sb_20170209+
    IF (npe_method == 2 .AND. (spe_method == 3 .OR. spe_method==5))  THEN

       CALL new_channel_object(status, modstr,'PRate_part', p3=PRate_part_out)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr ,'PRate_part' &
            , 'long_name', c='SPE ionisationsrate due to particle ionization')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'PRate_part' &
            , 'units', c='ions/cm3/s')
       CALL channel_halt(substr, status)
       CALL new_channel_object(status, modstr,'PRate_phi', p3=PRate_phi_out)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr ,'PRate_phi' &
            , 'long_name', c='SPE ionisationsrate due to photoionization')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'PRate_phi' &
            , 'units', c='ions/cm3/s')
       CALL channel_halt(substr, status)
       
     END IF !(npe_method==2 and spe_method==3)
    !ka_sb_20170209-

    !ka_sv_20171018+
       IF (switchphi == 1)  THEN
       
         IF (npe_method.ne.2 .AND. (spe_method.ne.3 .OR. spe_method.ne.5)) THEN
            CALL new_channel_object(status, modstr,'PRate_phi', p3=PRate_phi_out)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr ,'PRate_phi' &
               , 'long_name', c='SPE ionisationsrate due to photoionization')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'PRate_phi' &
               , 'units', c='ions/cm3/s')
            CALL channel_halt(substr, status)
         END IF
    !ka_sv_20171018-

          CALL new_channel_object(status, modstr,'RNOp_out', p3=RNOp_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'RNOp_out' &
               , 'long_name', c='RNOp_out due to photoionization')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'RNOp_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'RNp_out', p3=RNp_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'RNp_out' &
               , 'long_name', c='RNp_out due to photoionization')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'RNp_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'RN2p_out', p3=RN2p_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'RN2p_out' &
               , 'long_name', c='RN2p_out due to photoionization')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'RN2p_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'ROp2P_out', p3=ROp2P_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'ROp2P_out' &
               , 'long_name', c='ROp2P_out due to photoionization')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'ROp2P_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'ROp2D_out', p3=ROp2D_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'ROp2D_out' &
               , 'long_name', c='ROp2D_out due to photoionization')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'ROp2D_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'RO2p_out', p3=RO2p_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'RO2p_out' &
               , 'long_name', c='RO2p_out due to photoionization')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'RO2p_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'ROp4S_out', p3=ROp4S_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'ROp4S_out' &
               , 'long_name', c='ROp4S_out')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'ROp4S_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)
          
          CALL new_channel_object(status, modstr,'RN_out', p3=RN_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'RN_out' &
               , 'long_name', c='RN_out')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'RN_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)
          
          CALL new_channel_object(status, modstr,'RN2D_out', p3=RN2D_out)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr ,'RN2D_out' &
               , 'long_name', c='RN2D_out')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'RN2D_out' &
               , 'units', c='ions/cm3/s')
          CALL channel_halt(substr, status)
          
       END IF !(switchphi==1)


    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY   
       ALLOCATE(pflxspc(nbin))
       ALLOCATE(energies(nbin))
       ALLOCATE(rangestp(nbin))
    CASE(1)
       ! nothing for now
    CASE(2)
       ! nothing for now
! ka_sb_20161224+
    CASE(3)
       ! nothing for now
! ka_sb_20161224-
! ka_sv_20160615+  
    CASE(4)
       ! nothing for now
! ka_sv_20160615-
    CASE(5)
    
    END SELECT

  END SUBROUTINE spe_init_memory
! ========================================================================

! ========================================================================
  SUBROUTINE spe_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_channel,       ONLY: get_channel_object &
                                      , get_channel_object_dimvar
    USE messy_main_tracer,        ONLY: get_tracer_list &
                                      , get_tracer !ka_sb_20160224
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_ULONG
    USE messy_main_tools,         ONLY: PTR_1D_ARRAY

    IMPLICIT NONE
    INTRINSIC :: SIZE, NULL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_init_coupling'
    INTEGER :: status
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: subnames => NULL()
    TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER           :: dvs => NULL()
    CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:), POINTER  :: units => NULL()
    INTEGER :: jt
    INTEGER :: i_err !ka_sb_20160224

    CALL start_message_bi(modstr, 'LOOKING FOR TRACERS', substr)

    CALL get_tracer_list(status, GPTRSTR, 'N', idt_list_n, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_n = SIZE(idt_list_n)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_n
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... N'
          ELSE
             WRITE(*,*) ' ... N_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF

    CALL get_tracer_list(status, GPTRSTR, 'NO', idt_list_no, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_no = SIZE(idt_list_no)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_no
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... NO'
          ELSE
             WRITE(*,*) ' ... NO_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF

    CALL get_tracer_list(status, GPTRSTR, 'OH', idt_list_oh, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_oh = SIZE(idt_list_oh)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_oh
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... OH'
          ELSE
             WRITE(*,*) ' ... OH_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF

    CALL get_tracer_list(status, GPTRSTR, 'H', idt_list_h, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_h = SIZE(idt_list_h)
    IF (p_parallel_io) THEN
       DO jt=1, nlntrac_h
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... H'
          ELSE
             WRITE(*,*) ' ... H_'//TRIM(subnames(jt))
          END IF
       END DO
    END IF

!ka_sb_20160224: tracers needed for SPE_NPE (npe_method 2)
      
    CALL get_tracer(i_err, GPTRSTR, TRIM('O3P'), idx=idt_O3P)
    CALL tracer_halt(substr, i_err)

    CALL get_tracer(i_err, GPTRSTR, TRIM('O1D'), idx=idt_O1D)
    CALL tracer_halt(substr, i_err)

    CALL get_tracer(i_err, GPTRSTR, TRIM('H2O'), idx=idt_H2O)
    CALL tracer_halt(substr, i_err)

    CALL get_tracer(i_err, GPTRSTR, TRIM('N2'), idx=idt_N2)
    CALL tracer_halt(substr, i_err)

    CALL get_tracer(i_err, GPTRSTR, TRIM('O2'), idx=idt_O2)
    CALL tracer_halt(substr, i_err)

    !ka_sv_20170426+
    CALL get_tracer(i_err, GPTRSTR, TRIM('N2O'), idx=idt_n2o)
    CALL tracer_halt(substr, i_err)
    !ka_sv_20170426-

!ka_sb_20160224-

! ka_sv_20171213+
    IF (spe_method==3 .OR. spe_method==5) THEN
     IF (calc_pos_ions) THEN
      CALL get_tracer(i_err, GPTRSTR, TRIM('N2D'), idx=idt_n2d)
      CALL tracer_halt(substr, i_err)
      CALL get_tracer(i_err, GPTRSTR, TRIM('N'), idx=idt_n4s)
      CALL tracer_halt(substr, i_err)
      CALL get_tracer(i_err, GPTRSTR, TRIM('N2p'), idx=idt_n2p)
      CALL tracer_halt(substr, i_err)
      CALL get_tracer(i_err, GPTRSTR, TRIM('Np'), idx=idt_np)
      CALL tracer_halt(substr, i_err)
      CALL get_tracer(i_err, GPTRSTR, TRIM('O2p'), idx=idt_o2p)
      CALL tracer_halt(substr, i_err)
      !ka_sv_20180327+
      IF (calc_three_o_ions) THEN
        CALL get_tracer(i_err, GPTRSTR, TRIM('O4Sp'), idx=idt_o4sp)
        CALL get_tracer(i_err, GPTRSTR, TRIM('O2Dp'), idx=idt_o2dp)
        CALL get_tracer(i_err, GPTRSTR, TRIM('O2Pp'), idx=idt_o2pp)
      ELSE
        CALL get_tracer(i_err, GPTRSTR, TRIM('Op'), idx=idt_op)
      END IF
      !ka_sv_20180327-
      CALL tracer_halt(substr, i_err)
      CALL get_tracer(i_err, GPTRSTR, TRIM('NOp'), idx=idt_nop)
      CALL tracer_halt(substr, i_err)
      CALL get_tracer(i_err, GPTRSTR, TRIM('em'), idx=idt_em)
      CALL tracer_halt(substr, i_err)
     END IF ! calc_pos_ions
    END IF ! spe_method
! ka_sv_20171213-

    IF (ASSOCIATED(subnames)) DEALLOCATE(subnames)
    NULLIFY(subnames)

    CALL end_message_bi(modstr, 'LOOKING FOR TRACERS', substr)

    CALL start_message_bi(modstr, 'LOOKING FOR REQUIRED CHANNEL OBJECTS', substr)

    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY
       CALL get_channel_object(status &
            , TRIM(spe_data_int%cha), TRIM(spe_data_int%obj), p1=pflx)
       CALL channel_halt(substr, status)
       npch = SIZE(pflx)

       CALL get_channel_object_dimvar(status &
            , TRIM(spe_data_int%cha), TRIM(spe_data_int%obj) &
            , dvs, units)
       CALL channel_halt(substr, status)
       chrig => dvs(1)%ptr

    CASE(1) ! EXTERNAL IONIZATION RATES
       CALL get_channel_object(status &
            , TRIM(spe_data_ion%cha), TRIM(spe_data_ion%obj), p1=ions_ext)
       CALL channel_halt(substr, status)        
       naltitudes = SIZE(ions_ext)

    CASE(2) ! EXTERNAL IONIZATION RATES FROM AIMOS / JAN MAIK WISSING 
       CALL get_channel_object(status, TRIM(AIMOS_channel), TRIM(AIMOS_p_object) &
            , p3=AIMOS_p)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, TRIM(AIMOS_channel), TRIM(AIMOS_e_object) &
            , p3=AIMOS_e)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, TRIM(AIMOS_channel), TRIM(AIMOS_a_object) &
            , p3=AIMOS_a)
       CALL channel_halt(substr, status)

       !ka_sb_20160705+
    CASE(3)  ! use parameterization from Holger Nieder based on AIMOS rates (kp-Index needed)
       CALL get_channel_object(status &
            , TRIM(spe_kp%cha), TRIM(spe_kp%obj), p1=kp_data)
       CALL channel_halt(substr, status) 
       !ka_sb_20160705-

       ! ka_sv_20160615+
    CASE(4)  ! EXTERNAL RATES FROM CMIP6
       CALL get_channel_object(status, TRIM(CMIP6_channel) &
            , TRIM(CMIP6_p_object), p3=CMIP6_p)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, TRIM(CMIP6_channel) &
            , TRIM(CMIP6_e_object), p3=CMIP6_e)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, TRIM(CMIP6_channel) &
            , TRIM(CMIP6_g_object), p3=CMIP6_g)
       CALL channel_halt(substr, status)
       ! ka_sv_20160615+
    CASE(5)
       ! nothing
    END SELECT
    
    !ka_sv_20180515+
    IF (switchphi==1) THEN
      CALL get_channel_object(status &
            , TRIM(spe_flux107%cha), TRIM(spe_flux107%obj), p1=flux107_data)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status &
            , TRIM(spe_cossza%cha), TRIM(spe_cossza%obj), p2=cossza_data)
       CALL channel_halt(substr, status)
    END IF
    !ka_sv_20180515-

    CALL end_message_bi(modstr, 'LOOKING FOR REQUIRED CHANNEL OBJECTS', substr)

  END SUBROUTINE spe_init_coupling
! ========================================================================

! ========================================================================
  SUBROUTINE spe_global_start
  
    USE messy_main_timer, ONLY: DAY, HOUR, YEAR, DAYOFYEAR
    !USE messy_main_mpi_bi,       ONLY: p_parallel_io
    USE messy_main_tools ,    ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_global_start'
    INTEGER            :: status
    LOGICAL, SAVE :: startofday = .TRUE.
    REAL(dp), SAVE :: oldday = 0._dp
    INTEGER :: iou
!    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: ionrate_aimos_used

    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY
       CALL spe_interp
    CASE(1) ! EXTERNAL IONIZATION RATES
       ! nothing to be done currently
    CASE(2) ! EXTERNAL IONIZATION RATES FROM AIMOS
       ! nothing to be done currently
! ka_sb_20160705+
   CASE(3) ! EXTERNAL parameterization from Holger Nieder based on AIMOS rates
      ! nothing to be done currently
! ka_sb_20160705-
! ka_sv_20160615+
    CASE(4) ! EXTERNAL IONIZATION RATES FROM CMIP6
       ! nothing to be done currently   
! ka_sv_20160615-
! ka_sv_20180514+
    CASE(5) ! Read original AIMOS Data (daily files with protons, electrons and alphas included directly on model grid)
      IF (.NOT. ASSOCIATED(ionrate_aimos_used)) ALLOCATE(ionrate_aimos_used(spe_nlon,spe_nlev,spe_nlat))
      ionrate_aimos_used(:,:,:)=0._dp
      IF (abs(DAY-oldday)>0.5_dp) THEN
        startofday = .TRUE.
      END IF
      IF (startofday) THEN  ! only when a new day starts
        IF (.NOT. ALLOCATED(ionrate_aimos_orig)) THEN
          ALLOCATE(ionrate_aimos_orig(spe_nlon, spe_nlev, spe_nlat, aimos_time,3))
        END IF
        iou = find_next_free_unit(100,200)
        CALL spe_read_aimos_orig(ionrate_aimos_orig,spe_aimos_dir,spe_aimos_prefix,year,DAYOFYEAR &
                                 ,spe_nlat,spe_nlon,spe_nlev,aimos_time,aimpro,aimele,aimalp,iou)
        startofday=.FALSE.
      END IF  ! startofday
      oldday=DAY
      SELECT CASE (aimos_time)
        CASE(1)
          ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,1,:),4)
        CASE(2)
          IF (HOUR.ge.0 .AND. HOUR.le.11) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,1,:),4)
          IF (HOUR.ge.12 .AND. HOUR.le.24) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,2,:),4)
        CASE(4)
          IF (HOUR.ge.0 .AND. HOUR.le.5) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,1,:),4)
          IF (HOUR.ge.6 .AND. HOUR.le.11) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,2,:),4)
          IF (HOUR.ge.12 .AND. HOUR.le.17) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,3,:),4)
          IF (HOUR.ge.18 .AND. HOUR.le.24) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,4,:),4)
        CASE(6)
          IF (HOUR.ge.0 .AND. HOUR.le.3) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,1,:),4)
          IF (HOUR.ge.4 .AND. HOUR.le.7) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,2,:),4)
          IF (HOUR.ge.8 .AND. HOUR.le.11) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,3,:),4)
          IF (HOUR.ge.12 .AND. HOUR.le.15) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,4,:),4)
          IF (HOUR.ge.16 .AND. HOUR.le.19) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,5,:),4)
          IF (HOUR.ge.20 .AND. HOUR.le.24) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,6,:),4)
        CASE(12)
          IF (HOUR.ge.0 .AND. HOUR.le.1) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,1,:),4)
          IF (HOUR.ge.2 .AND. HOUR.le.3) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,2,:),4)
          IF (HOUR.ge.4 .AND. HOUR.le.5) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,3,:),4)
          IF (HOUR.ge.6 .AND. HOUR.le.7) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,4,:),4)
          IF (HOUR.ge.8 .AND. HOUR.le.9) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,5,:),4)
          IF (HOUR.ge.10 .AND. HOUR.le.11) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,6,:),4)
          IF (HOUR.ge.12 .AND. HOUR.le.13) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,7,:),4)
          IF (HOUR.ge.14 .AND. HOUR.le.15) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,8,:),4)
          IF (HOUR.ge.16 .AND. HOUR.le.17) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,9,:),4)
          IF (HOUR.ge.18 .AND. HOUR.le.19) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,10,:),4)
          IF (HOUR.ge.20 .AND. HOUR.le.21) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,11,:),4)
          IF (HOUR.ge.22 .AND. HOUR.le.24) ionrate_aimos_used(:,:,:)=sum(ionrate_aimos_orig(:,:,:,12,:),4)
        
        CASE DEFAULT
          WRITE(*,*) 'ERROR: aimos_time not supported'
          STOP
       END SELECT !aimos_time
! ka_sv_20180514-
    END SELECT !spe_method

  END SUBROUTINE spe_global_start
! ========================================================================

! ========================================================================
  SUBROUTINE spe_physc

    ! ECHAM5/MESSy
    USE messy_main_constants_mem, ONLY: M_air, R_gas, g, N_A
    USE messy_main_data_bi,       ONLY: pint => pressi_3d & ! level interface (above) pressure [Pa]
         , pmid => press_3d  & ! mid-level pressures [Pa]
         , press_3d &
         , t_scb             &
         , geopot_3d         &
         , temp => t_scb
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma      &
         , nlev              &
         , nproma, ngpblks
    USE messy_main_grid_def_bi,     ONLY:philat, philon    &
         , ilat, ilon        &
         , grmass, grvol
    USE messy_main_timer,         ONLY: time_step_len    
    ! ka_sb_20160303+
    USE messy_main_tracer_mem_bi, ONLY: pxtte=>qxtte, pxtm1 => qxtm1
    ! ka_sb_20160303-
    ! ka_sv_20180514+
    USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl
    ! ka_sv_20180514-

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'spe_physc'
    INTEGER                     :: status   ! error status
    REAL(DP)                    :: zpb, zpt
    INTEGER                     :: jk, jp, jt,perionlev,level
    REAL(DP)                    :: scale_height_ma = 7._dp ! Middle atmosphere scale height [km]

    ! LOCAL FIELDS
    REAL(dp),  DIMENSION(:,:), POINTER    :: grheight ! height of box [m]
    REAL(dp),  DIMENSION(:,:), POINTER    :: density  ! air density at level midpoints [g cm-3]
    REAL(dp),  DIMENSION(:,:), POINTER    :: densityint  ! air density at level interfaces [g cm-3]
    REAL(dp),  DIMENSION(:,:), POINTER    :: altitude  ! altitude [km]
! ka_sb_20160303+   
    REAL(dp),  DIMENSION(kproma)          :: maglat
    REAL(dp),  DIMENSION(:,:), POINTER    :: Nperion
    REAL(dp),  DIMENSION(:,:), POINTER    :: NOperion
    REAL(dp),  DIMENSION(kproma,nlev)     :: cair !gas density [cm-3] only needed for npe_method 2
! ka_sb_20160303-
! ka_sb_20161121+
    REAL(dp),   DIMENSION(kproma,nlev)    :: O2_t0,O3P_t0,O1D_t0,N2_t0  !vmr O2,O3,O1D,N2 at t0
    REAL(dp),   DIMENSION(kproma,nlev)    :: col_O2_box,col_O_box,col_N2_box !partial O2,O,N2 column
    REAL(dp),   DIMENSION(kproma,nlev)    :: RO2p,ROp4S,ROp2D,ROp2P,RN2p,RNp,RN,RNOp,RN2D !individual ionization rates
    REAL(dp),   DIMENSION(kproma,nlev)    :: RO2p_phi,ROp4S_phi,ROp2D_phi,ROp2P_phi,RN2p_phi,RNp_phi,RN_phi,RNOp_phi !individual ionization rates photoionization
    REAL(dp),   DIMENSION(kproma,nlev)    :: PRate_phi,PRate_part ! sum of ionization rates due to photoionization / particle ionization
! ka_sb_20161121-

    ! ALLOCATE MEMORY
    ALLOCATE(grheight(kproma,nlev))
    ALLOCATE(density(kproma,nlev))
    ALLOCATE(densityint(kproma,nlev))
    ALLOCATE(altitude(kproma,nlev))
    ALLOCATE(Nperion(kproma,nlev))
    ALLOCATE(NOperion(kproma,nlev))
    ! ka_sv_20180514+
    IF (spe_method==5) THEN
      ALLOCATE(ionrate_aimos_used_dec(nproma,nlev,ngpblks))
      CALL trp_gpdc_gpgl(-1, ionrate_aimos_used_dec, ionrate_aimos_used)
    END IF
    ! ka_sv_20180514-

    altitude(1:kproma,:)=geopot_3d(1:kproma,:,jrow)/g/1E3_dp

    ! CALCULATE DENSITY, convert from g/m3 to g/cm3
    density(1:kproma,:)=pmid(1:kproma,:,jrow)*M_air/(temp(1:kproma,:,jrow)*R_gas)*1e-6_dp
    densityint(1:kproma,:)=pint(1:kproma,1:nlev,jrow)*M_air/(temp(1:kproma,:,jrow)*R_gas)*1e-6_dp

    ! HEIGHT OF GRID-BOX [m]
    DO jk=1, nlev
       DO jp=1, kproma
          zpb = pint(jp,jk+1, jrow)
          ! ECHAM5 top layer ends at 0. Pa !!! adjust to 1. Pa
! ka_sv_20160614+
!          zpt = MAX(pint(jp, jk, jrow),1._dp)
          zpt = MAX(pint(jp, jk, jrow),pmid(jp,jk,jrow))
! ka_sv_20160614-
          grheight(jp,jk) = (1000._dp * R_gas / (M_air * g)) &
               * temp(jp,jk,jrow) * log(zpb/zpt)
       END DO
    END DO

! ka_sb_20160707+
    SELECT CASE(spe_method) 
    CASE(0,1,2,4)
    CASE(3) ! External ionrates by H. Nieder

     IF (kp_data(1)<0._dp) then
       ikp=1
     END IF
     IF (kp_data(1)>hn_ion_kp(14)) then
      ikp=14
     END IF
     DO i=2,nkp
       IF (kp_data(1)<=hn_ion_kp(i) .AND. kp_data(1)>hn_ion_kp(i-1)) THEN
         ikp=i-1
       END IF
     END DO
     END SELECT
! ka_sb_2016707-    

    CALL  SPE_IONS(nlev,time_step_len           & ! INPUT 
                  ,grheight                     & ! INPUT
                  ,density,densityint           & ! INPUT
                  ,altitude                     & ! INPUT
                  ,ilat, ilon                   & ! INPUT
                  ,philat, philon               & ! INPUT
                  ,jrow, kproma                 & ! INPUT
                  ,spe_method                   & ! INPUT
                ! ka_sv_20160615+
                  ,grmass,grvol                 & ! INPUT
                ! ka_sv_20160615-
                ! ka_sb_20160721+
                  ,hn_ion_lat_bin               & ! INPUT
                  ,ikp,nkp,nlat,nlev_a          & ! INPUT
                  ,press_3d                     &
                  ,pressure                     & ! INPUT
                  ,ionrate                      & ! INPUT
                ! ka_sb_20160721-
                  ,ions                         & ! OUTPUT   ! only ions due to particles
                  ,maglat                       & ! OUTPUT
                  ,status                       & ! OUTPUT
                         )
! ka_sv_20180514+
    IF (spe_method==5) ions(1:kproma,1:nlev,jrow)=ionrate_aimos_used_dec(1:kproma,1:nlev,jrow) ! fits in SPE_IONS
! ka_sv_20180514-

!ka_sb_20160222+
    SELECT CASE (npe_method)

    CASE (0) ! constant values
        Nperion(:,:)= 0.55_dp
        NOperion(:,:)= 0.7_dp
        !ka_sv_20171018+
        IF (switchphi == 1)  THEN
         level_loop1d0: DO jk=1,nlev
          vector_loop1d0: DO jp=1,kproma
               cair(jp,jk)  = (N_A) * pmid(jp,jk,jrow) / (R_gas*temp(jp,jk,jrow)) / 1.e06_dp ! cm^(-3)
! ka_sv_20170714+
               col_O2_box(jp,jk) = max(cair(jp,jk) * pxtm1(jp,jk,idt_O2)*grheight(jp,jk)*100,1.e-02_dp) !partial column (altitude in cm)
               col_O_box(jp,jk) = max(cair(jp,jk) * (pxtm1(jp,jk,idt_O3P)+pxtm1(jp,jk,idt_O1D))*grheight(jp,jk)*100,1.e-02_dp) !partial column (altitude in cm)
               col_N2_box(jp,jk) = max(cair(jp,jk) * pxtm1(jp,jk,idt_N2)*grheight(jp,jk)*100,1.e-02_dp) !partial column (altitude in cm)
! ka_sv_20170714-
             END DO vector_loop1d0
           END DO level_loop1d0
           DO jp=1, kproma
           !photoionization
              CALL spe_phioniz(cossza_data(jp,jrow) & 
                               ,RO2p_phi(jp,:),ROp4S_phi(jp,:),ROp2D_phi(jp,:),ROp2P_phi(jp,:),RN2p_phi(jp,:),RNp_phi(jp,:),RN_phi(jp,:)& !OUTPUT
                               ,pxtm1(jp,:,idt_O2)  & !INPUT: vmr O2 (ZO2)
                               ,pxtm1(jp,:,idt_O3P) & !INPUT: vmr O3P , sum O3P and O1D for ZO
                               ,pxtm1(jp,:,idt_O1D) & !INPUT: vmr O1D
                               ,pxtm1(jp,:,idt_N2)  & !INPUT: vmr N2 (ZN2)
                               ,cair(jp,:) & !CONCM
                               ,col_O2_box(jp,:),col_O_box(jp,:),col_N2_box(jp,:)  &  ! O2 partial columns
                               ,pint(jp,:,jrow) &
                               ,temp(jp,:,jrow) &
                               ,nlev, flux107_data(:), grheight(jp,:))
           END DO !(kproma)
          RNOp_phi = 0.0 

        ELSE 
          RO2p_phi = 0.0
          ROp4S_phi= 0.0
          ROp2D_phi= 0.0
          ROp2P_phi= 0.0
          RN2p_phi = 0.0
          RNp_phi  = 0.0
          RN_phi   = 0.0      !=RN4S
          RNOp_phi = 0.0 
       ENDIF !Photoionization
        PRate_phi(1:kproma,:)=RN_phi(1:kproma,:)
        ions(1:kproma,:,jrow)=ions(1:kproma,:,jrow)+PRate_phi(1:kproma,:) !=PRate_part+Prate_phi
       !ka_sv_20171018-
    
    CASE (1) ! using internally set values
        Nperion(1:kproma,:)=0.
        NOperion(1:kproma,:)=0.
     level_loop1a:  DO level=1, nlev
       DO perionlev=1, NMAXIONKM
          IF (ion_km(perionlev) < 0.) EXIT
          IF (altitude(1,level) >= ion_km(perionlev)) THEN
            DO jp=1,kproma
              Nperion(jp,level)=Nperion_km(perionlev)
              NOperion(jp,level)=NOperion_km(perionlev)
            END DO
          ENDIF
       END DO
      END DO level_loop1a
      !ka_sv_20171018+
      IF (switchphi == 1)  THEN
         level_loop1d1: DO jk=1,nlev
          vector_loop1d1: DO jp=1,kproma
               cair(jp,jk)  = (N_A) * pmid(jp,jk,jrow) / (R_gas*temp(jp,jk,jrow)) / 1.e06_dp ! cm^(-3)
! ka_sv_20170714+
               col_O2_box(jp,jk) = max(cair(jp,jk) * pxtm1(jp,jk,idt_O2)*grheight(jp,jk)*100,1.e-02_dp) !partial column (altitude in cm)
               col_O_box(jp,jk) = max(cair(jp,jk) * (pxtm1(jp,jk,idt_O3P)+pxtm1(jp,jk,idt_O1D))*grheight(jp,jk)*100,1.e-02_dp) !partial column (altitude in cm)
               col_N2_box(jp,jk) = max(cair(jp,jk) * pxtm1(jp,jk,idt_N2)*grheight(jp,jk)*100,1.e-02_dp) !partial column (altitude in cm)
! ka_sv_20170714-
             END DO vector_loop1d1
           END DO level_loop1d1
           DO jp=1, kproma
           !photoionization
              CALL spe_phioniz(cossza_data(jp,jrow) & 
                               ,RO2p_phi(jp,:),ROp4S_phi(jp,:),ROp2D_phi(jp,:),ROp2P_phi(jp,:),RN2p_phi(jp,:),RNp_phi(jp,:),RN_phi(jp,:)& !OUTPUT
                               ,pxtm1(jp,:,idt_O2)  & !INPUT: vmr O2 (ZO2)
                               ,pxtm1(jp,:,idt_O3P) & !INPUT: vmr O3P , sum O3P and O1D for ZO
                               ,pxtm1(jp,:,idt_O1D) & !INPUT: vmr O1D
                               ,pxtm1(jp,:,idt_N2)  & !INPUT: vmr N2 (ZN2)
                               ,cair(jp,:) & !CONCM
                               ,col_O2_box(jp,:),col_O_box(jp,:),col_N2_box(jp,:)  &  ! O2 partial columns
                               ,pint(jp,:,jrow) &
                               ,temp(jp,:,jrow) &
                               ,nlev, flux107_data(:), grheight(jp,:))
           END DO !(kproma)
          RNOp_phi = 0.0  ! ab hier nur für test
          !RO2p_phi = 0.0
          !ROp4S_phi= 0.0
          !ROp2D_phi= 0.0
          !ROp2P_phi= 0.0
          !RN2p_phi = 0.0
          !RNp_phi  = 0.0 

        ELSE 
          RO2p_phi = 0.0
          ROp4S_phi= 0.0
          ROp2D_phi= 0.0
          ROp2P_phi= 0.0
          RN2p_phi = 0.0
          RNp_phi  = 0.0
          RN_phi   = 0.0      !=RN4S
          RNOp_phi = 0.0 
       ENDIF !Photoionization
        PRate_phi(1:kproma,:)=RN_phi(1:kproma,:)
        !PRate_phi(1:kproma,:)=RO2p_phi(1:kproma,:)+ROp4S_phi(1:kproma,:)+ROp2D_phi(1:kproma,:)+ROp2P_phi(1:kproma,:)+RN2p_phi(1:kproma,:)+RNp_phi(1:kproma,:)
        ions(1:kproma,:,jrow)=ions(1:kproma,:,jrow)+PRate_phi(1:kproma,:) !=PRate_part+Prate_phi
       !ka_sv_20171018-

    CASE (2) ! NOx PRODUCTION RATES FROM HOLGER NIEDER
! calculate gas densitys
  level_loop1d: DO jk=1,nlev
    vector_loop1d: DO jp=1,kproma
        cair(jp,jk)  = (N_A) * pmid(jp,jk,jrow) / (R_gas*temp(jp,jk,jrow)) / 1.e06_dp ! cm^(-3)
! ka_sv_20170714+
        !col_O2_box(jp,jk) = cair(jp,jk) * pxtm1(jp,jk,idt_O2)*grheight(jp,jk)*100 !partial column (altitude in cm)
        !col_O_box(jp,jk) = cair(jp,jk) * (pxtm1(jp,jk,idt_O3P)+pxtm1(jp,jk,idt_O1D))*grheight(jp,jk)*100 !partial column (altitude in cm)
        !col_N2_box(jp,jk) = cair(jp,jk) * pxtm1(jp,jk,idt_N2)*grheight(jp,jk)*100 !partial column (altitude in cm)
        col_O2_box(jp,jk) = max(cair(jp,jk) * pxtm1(jp,jk,idt_O2)*grheight(jp,jk)*100.,1.e-20_dp) !partial column (altitude in cm)
        col_O_box(jp,jk) = max(cair(jp,jk) * (pxtm1(jp,jk,idt_O3P)+pxtm1(jp,jk,idt_O1D))*grheight(jp,jk)*100.,1.e-20_dp) !partial column (altitude in cm)
        col_N2_box(jp,jk) = max(cair(jp,jk) * pxtm1(jp,jk,idt_N2)*grheight(jp,jk)*100.,1.e-20_dp) !partial column (altitude in cm)
! ka_sv_20170714-
      END DO vector_loop1d
    END DO level_loop1d
!ka_sb_20161109+
! switch for photoionisation
        IF ((spe_method == 3 .OR. spe_method==5) .AND. switchphi == 1)  THEN 
           DO jp=1, kproma
           !photoionization
              CALL spe_phioniz(cossza_data(jp,jrow) & 
                               ,RO2p_phi(jp,:),ROp4S_phi(jp,:),ROp2D_phi(jp,:),ROp2P_phi(jp,:),RN2p_phi(jp,:),RNp_phi(jp,:),RN_phi(jp,:)& !OUTPUT
                               ,pxtm1(jp,:,idt_O2)  & !INPUT: vmr O2 (ZO2)
                               ,pxtm1(jp,:,idt_O3P) & !INPUT: vmr O3P , sum O3P and O1D for ZO
                               ,pxtm1(jp,:,idt_O1D) & !INPUT: vmr O1D
                               ,pxtm1(jp,:,idt_N2)  & !INPUT: vmr N2 (ZN2)
                               ,cair(jp,:) & !CONCM
                               ,col_O2_box(jp,:),col_O_box(jp,:),col_N2_box(jp,:)  &  ! O2 partial columns
                               ,pint(jp,:,jrow) &
                               ,temp(jp,:,jrow) &
                               ,nlev, flux107_data(:), grheight(jp,:))
           END DO !(kproma)
          RNOp_phi = 0.0

        ELSE 
          RO2p_phi = 0.0
          ROp4S_phi= 0.0
          ROp2D_phi= 0.0
          ROp2P_phi= 0.0
          RN2p_phi = 0.0
          RNp_phi  = 0.0
          RN_phi   = 0.0      !=RN4S
          RNOp_phi = 0.0 
       ENDIF !spe_method=3 and Photoionization   
!ka_sb_20161109-          

      DO jp=1, kproma
       DO jk=1, nlev
          CALL SPE_NPE(temp(jp,jk,jrow)              & ! INPUT
                      ,cair(jp,jk)                   & ! INPUT (Airnumberdensity)
                      ,pxtm1(jp,jk,idt_O3P)          & ! INPUT
                      ,pxtm1(jp,jk,idt_O1D)          & ! INPUT
                      ,pxtm1(jp,jk,idt_list_n(1))    & ! INPUT (old: NATOM)
                      ,pxtm1(jp,jk,idt_list_no(1))   & ! INPUT
                      ,pxtm1(jp,jk,idt_H2O)          & ! INPUT
                      ,ions(jp,jk,jrow)              & ! INPUT (old: DIONP)
                      ,NOperion(jp,jk)               & ! OUTPUT (old: RNO_POS)
                      ,Nperion(jp,jk)                & ! OUTPUT (old: RN_4S_POS)
                      ,pxtm1(jp,jk,idt_N2)           & ! INPUT
                      ,pxtm1(jp,jk,idt_O2)           & ! INPUT
                     !ka_sb_20161111+
                      ,RO2p_phi(jp,jk),ROp4S_phi(jp,jk),ROp2D_phi(jp,jk),ROp2P_phi(jp,jk),RN2p_phi(jp,jk),RNp_phi(jp,jk),RN_phi(jp,jk),RNOp_phi(jp,jk) & !INPUT & OUTPUT
                      ,RO2p(jp,jk),ROp4S(jp,jk),ROp2D(jp,jk),ROp2P(jp,jk),RN2p(jp,jk),RNp(jp,jk),RN(jp,jk),RNOp(jp,jk), RN2D(jp,jk) & !INPUT & OUTPUT
                      ,PRate_phi(jp,jk),PRate_part(jp,jk)) !OUTPUT
                     !ka_sb_20161111-
          ions(jp,jk,jrow)=ions(jp,jk,jrow)+PRate_phi(jp,jk) !=PRate_part+Prate_phi

        END DO
      END DO
      PRate_phi_out(1:kproma,:,jrow)=PRate_phi(1:kproma,:)
      PRate_part_out(1:kproma,:,jrow)=PRate_part(1:kproma,:)
      RO2p_out(1:kproma,:,jrow)=RO2p(1:kproma,:)
      ROp4S_out(1:kproma,:,jrow)=ROp4S(1:kproma,:)
      ROp2D_out(1:kproma,:,jrow)=ROp2D(1:kproma,:)
      ROp2P_out(1:kproma,:,jrow)=ROp2P(1:kproma,:)
      RN2p_out(1:kproma,:,jrow)=RN2p(1:kproma,:)
      RNp_out(1:kproma,:,jrow)=RNp(1:kproma,:)
      RN_out(1:kproma,:,jrow)=RN(1:kproma,:)
      RNOp_out(1:kproma,:,jrow)=RNOp(1:kproma,:)
      RN2D_out(1:kproma,:,jrow)=RN2D(1:kproma,:)
    END SELECT !(case npe_method)

         CALL SPE_PROD_XNOX(density                    & ! INPUT
                           ,altitude                   & ! INPUT
                           ,nlev                       & ! INPUT
                           ,jrow, kproma               & ! INPUT
                           ,Nperion, NOperion          & ! INPUT
                           ,spe_method                 & ! INPUT
                           ,ions                       & ! INPUT 
                           ,maglat                     & ! INPUT
                           ,xnox, tespen, tespeno      & ! OUTPUT 
                           ,xhox, tespeh, tespeoh      & ! OUTPUT 
                           ,xhno3, tespehno3           & ! OUTPUT 
!ka_sv_20170426+
                           ,n2oprod,xn2o,tespen2o,PRate_part & 
!ka_sv_20170426-
                           ,status                     & ! OUTPUT
                            )
!ka_sb_20160222- 
!ka_sv_20171213+
   IF (calc_pos_ions) THEN
     !pxtte(1:kproma,:,idt_n4s) = pxtte(1:kproma,:,idt_n4s) + &
     !          RN_out(1:kproma,:,jrow)/cair(1:kproma,:)
     pxtte(1:kproma,:,idt_o2p) = pxtte(1:kproma,:,idt_o2p) + &
               RO2p_out(1:kproma,:,jrow)/cair(1:kproma,:)
     !ka_sv_20180327+
     IF (calc_three_o_ions) THEN
       pxtte(1:kproma,:,idt_o4sp) = pxtte(1:kproma,:,idt_o4sp) + &
               ROp4S_out(1:kproma,:,jrow)/cair(1:kproma,:)
       pxtte(1:kproma,:,idt_o2dp) = pxtte(1:kproma,:,idt_o2dp) + &
               ROp2D_out(1:kproma,:,jrow)/cair(1:kproma,:)
       pxtte(1:kproma,:,idt_o2pp) = pxtte(1:kproma,:,idt_o2pp) + &
               ROp2P_out(1:kproma,:,jrow)/cair(1:kproma,:)
     ELSE
       pxtte(1:kproma,:,idt_op) = pxtte(1:kproma,:,idt_op) + &
               (ROp4S_out(1:kproma,:,jrow)+ROp2D_out(1:kproma,:,jrow)+ROp2P_out(1:kproma,:,jrow))/cair(1:kproma,:)
     END IF
     !ka_sv_20180327-
     pxtte(1:kproma,:,idt_n2p) = pxtte(1:kproma,:,idt_n2p) + &
               RN2p_out(1:kproma,:,jrow)/cair(1:kproma,:)
     pxtte(1:kproma,:,idt_np) = pxtte(1:kproma,:,idt_np) + &
               RNp_out(1:kproma,:,jrow)/cair(1:kproma,:)
     pxtte(1:kproma,:,idt_nop) = pxtte(1:kproma,:,idt_nop) + &
               RNOp_out(1:kproma,:,jrow)/cair(1:kproma,:)
     pxtte(1:kproma,:,idt_em) = pxtte(1:kproma,:,idt_em) + &
               RO2p_out(1:kproma,:,jrow)/cair(1:kproma,:) + &
               (ROp4S_out(1:kproma,:,jrow)+ROp2D_out(1:kproma,:,jrow)+ROp2P_out(1:kproma,:,jrow))/cair(1:kproma,:) + &
               RN2p_out(1:kproma,:,jrow)/cair(1:kproma,:) + &
               RNp_out(1:kproma,:,jrow)/cair(1:kproma,:) + &
               RNOp_out(1:kproma,:,jrow)/cair(1:kproma,:)
! ka_sv_20180323+
     pxtte(1:kproma,:,idt_n4s) = pxtte(1:kproma,:,idt_n4s) + &
               RN_out(1:kproma,:,jrow)/cair(1:kproma,:)
     pxtte(1:kproma,:,idt_n2d) = pxtte(1:kproma,:,idt_n2d) + &
               RN2D_out(1:kproma,:,jrow)/cair(1:kproma,:)
! ka_sv_20180323-


   
   ELSE ! calculate N and NO directly; no positive ions
!ka_sv_20171213-
    IF (nlntrac_n > 0) THEN
       ! ADD SPE N TO TENDENCY
       DO jt=1, nlntrac_n
          pxtte(1:kproma,:,idt_list_n(jt)) = pxtte(1:kproma,:,idt_list_n(jt)) + &
               tespen(1:kproma,:,jrow)
      !ka_sv_20171018+
!          IF (switchphi==1 .AND. npe_method.ne.2 .AND. spe_method.ne.3) THEN
!            pxtte(1:kproma,:,idt_list_n(jt))=pxtte(1:kproma,:,idt_list_n(jt))+ &
!               RN(:,:)
!          END IF
      !ka_sv_20171018-
       END DO
    END IF
    IF (nlntrac_no > 0) THEN
       ! ADD SPE NO TO TENDENCY
       DO jt=1, nlntrac_no
          pxtte(1:kproma,:,idt_list_no(jt)) = pxtte(1:kproma,:,idt_list_no(jt)) + &
               tespeno(1:kproma,:,jrow)
       END DO
    END IF
! ka_sv_20180517+
  END IF !calc_pos_ions
  IF (spe_hox) THEN
    ! net reaction is h2o -> h + oh
    ! -> it is not possible to produce more h or oh than the amount of h2o
    DO jp=1,kproma
      DO jk=1,nlev
        tespeh(jp,jk,jrow)=min(tespeh(jp,jk,jrow),pxtm1(jp,jk,idt_h2o))
        tespeoh(jp,jk,jrow)=min(tespeoh(jp,jk,jrow),pxtm1(jp,jk,idt_h2o))
      END DO
    END DO
    tespeh2o(1:kproma,:,jrow)=-tespeh(1:kproma,:,jrow)-tespeoh(1:kproma,:,jrow)
    IF (nlntrac_h > 0) THEN
       ! ADD SPE H TO TENDENCY
       DO jt=1, nlntrac_h
          pxtte(1:kproma,:,idt_list_h(jt)) = pxtte(1:kproma,:,idt_list_h(jt)) + &
               tespeh(1:kproma,:,jrow)
       END DO
    END IF
    IF (nlntrac_oh > 0) THEN
       ! ADD SPE OH TO TENDENCY
       DO jt=1, nlntrac_oh
          pxtte(1:kproma,:,idt_list_oh(jt)) = pxtte(1:kproma,:,idt_list_oh(jt)) + &
               tespeoh(1:kproma,:,jrow)
       END DO
    END IF
    ! destruction of h2o was missing before -> added here
    pxtte(1:kproma,:,idt_h2o)=pxtte(1:kproma,:,idt_h2o)+tespeh2o(1:kproma,:,jrow)
  END IF !spe_hox
! ka_sv_20180517-
!ka_sv_20170426+
   ! at the moment there is an unknown major error -> production directly disabled in the code
    IF (n2oprod==1) THEN
      !pxtte(1:kproma,:,idt_n2o)=pxtte(1:kproma,:,idt_n2o)+tespen2o(1:kproma,:,jrow)
    END IF
!ka_sv_20170426-
   

    ! DEALLOCATE MEMORY
    DEALLOCATE(grheight)
    DEALLOCATE(density)
    DEALLOCATE(densityint)
    DEALLOCATE(altitude)
!ka_sv_20161017+
    DEALLOCATE(Nperion)
    DEALLOCATE(NOperion)
!ka_sv_20161017-
!ka_sv_20180514+
    IF (ASSOCIATED(ionrate_aimos_used_dec)) DEALLOCATE(ionrate_aimos_used_dec)
!ka_sv_20180514+

  END SUBROUTINE spe_physc
! ========================================================================

! ========================================================================
  SUBROUTINE spe_free_memory

    IMPLICIT NONE
    INTRINSIC :: ALLOCATED, ASSOCIATED

    IF (ASSOCIATED(idt_list_n))  DEALLOCATE(idt_list_n)
    IF (ASSOCIATED(idt_list_no)) DEALLOCATE(idt_list_no)
    IF (ASSOCIATED(idt_list_h))  DEALLOCATE(idt_list_h)
    IF (ASSOCIATED(idt_list_oh)) DEALLOCATE(idt_list_oh)
    IF (ALLOCATED(ionrate_aimos_orig)) DEALLOCATE(ionrate_aimos_orig)
    IF (ASSOCIATED(ionrate_aimos_used)) DEALLOCATE(ionrate_aimos_used)

    CALL spe_clean ! op_pj_20170220
  
 END SUBROUTINE spe_free_memory
! ========================================================================

! ========================================================================
  SUBROUTINE spe_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ spe_data_ion, spe_data_int &
         , AIMOS_channel, AIMOS_p_object &
         , AIMOS_e_object, AIMOS_a_object &
! ka_sb_20160706+
         , spe_kp &
! ka_sb_20160706-
! ka_sb_20161109+
         , spe_flux107 &
         , spe_cossza  &
! ka_sb_20161109- 
         , CMIP6_channel, CMIP6_p_object &
         , CMIP6_e_object, CMIP6_g_object

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='spe_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    SELECT CASE (spe_method) 
    CASE(0) ! CALCULATE IONRATES INTERNALLY
       IF (TRIM(spe_data_int%cha) == '') THEN
          CALL info_bi('ERROR: empty channel name for proton flux data')
          RETURN
       ELSE
          CALL info_bi('proton flux channel :'//spe_data_int%cha)
       END IF
       IF (TRIM(spe_data_int%obj) == '') THEN
          CALL info_bi('ERROR: empty channel object name for proton flux data')
          RETURN
       ELSE
          CALL info_bi('proton flux object  :'//spe_data_int%obj)
       END IF
    CASE(1) ! EXTERNAL IONIZATION RATES
       IF (TRIM(spe_data_ion%cha) == '') THEN
          CALL info_bi('ERROR: empty channel name for ionization rates')
          RETURN
       ELSE
          CALL info_bi('ionization rates channel :'//spe_data_ion%cha)
       END IF
       IF (TRIM(spe_data_ion%obj) == '') THEN
          CALL info_bi('ERROR: empty channel object name for ionization rates')
          RETURN
       ELSE
          CALL info_bi('ionization rates object  :'//spe_data_ion%obj)
       END IF
    CASE(2) ! EXTERNAL IONIZATION RATES FROM AIMOS / JAN MAIK WISSING
       WRITE(*,*) ' Proton channel object: ',  AIMOS_channel, AIMOS_p_object
       WRITE(*,*) ' Electron channel object: ', AIMOS_channel, AIMOS_e_object
       WRITE(*,*) ' Alpha particle channel object: ', AIMOS_channel, AIMOS_a_object
! ka_sb_20160706+
    CASE(3) ! EXTERNAL IONIZATION RATES FROM AIMOS / HOLGER NIEDER
       IF (TRIM(spe_kp%cha) == '') THEN
          CALL info_bi('ERROR: empty channel name for kp')
          RETURN
       ELSE
          CALL info_bi('kp channel :'//spe_kp%cha)
       END IF
       IF (TRIM(spe_kp%obj) == '') THEN
          CALL info_bi('ERROR: empty channel object name for kp')
          RETURN
       ELSE
          CALL info_bi('kp object  :'//spe_kp%obj)
       END IF
! ka_sb_20160706-
! ka_sb_20161109+
       IF (TRIM(spe_flux107%cha) == '') THEN
          CALL info_bi('ERROR: empty channel name for flux107')
          RETURN
       ELSE
          CALL info_bi('flux107 channel :'//spe_flux107%cha)
       END IF
       IF (TRIM(spe_flux107%obj) == '') THEN
          CALL info_bi('ERROR: empty channel object name for flux107')
          RETURN
       ELSE
          CALL info_bi('flux107 object  :'//spe_flux107%obj)
       END IF
       IF (TRIM(spe_cossza%cha) == '') THEN
          CALL info_bi('ERROR: empty channel name for cossza')
          RETURN
       ELSE
          CALL info_bi('cossza channel :'//spe_cossza%cha)
       END IF
       IF (TRIM(spe_cossza%obj) == '') THEN
          CALL info_bi('ERROR: empty channel object name for flux107')
          RETURN
       ELSE
          CALL info_bi('cossza object  :'//spe_cossza%obj)
       END IF
! ka_sb_20161109-   

! ka_sv_20160615+
    CASE(4) ! EXTERNAL IONIZATION RATES FROM CMIP6
       WRITE(*,*) ' Proton channel object: ',  CMIP6_channel, CMIP6_p_object
       WRITE(*,*) ' Electron channel object: ', CMIP6_channel, CMIP6_e_object
       WRITE(*,*) ' GCR particle channel object: ', CMIP6_channel, CMIP6_g_object
! ka_sv_20160615+
! ka_sv_20180514+
     CASE(5) ! AIMOS ORIGINAL
       WRITE(*,*) ' AIMOS Directory: ', spe_aimos_dir
       WRITE(*,*) ' AIMOS fileprefix: ', spe_aimos_prefix
! ka_sv_20180514-
    END SELECT

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE spe_read_nml_cpl
! ========================================================================

! ***********************************************************************
END MODULE messy_spe_e5
! ***********************************************************************

