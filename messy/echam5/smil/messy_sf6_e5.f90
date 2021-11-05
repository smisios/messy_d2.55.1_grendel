! DESCRIPTION:
! This module is the interface (called by messy_main_control) between the 
! ECHAM5 base model and the SF6 module.
!
! AUTHOR:
! T. Reddmann, Institute for Meteorology and Climate Research, KIT, Germany
!   (2001: original code)
! O. Kirner, Steinbuch Centre for Computing, KIT, Germany
! S. Versick, Steinbuch Centre for Computing, KIT, Germany
!   (2010: Creation of this submodel) (SV,OK)
!   (2012: Conversion to EMAC2) (SV)
!   (2013: Finished version 1.0) (SV)
!   (2013: Version 1.1: climatologies added and MAOAM restructured) (SV)
!   (2014: Version 1.2: implementation in EMAC2.50; some code improvements) (SV)
!   (2015: Version 1.3: code clean up (SV)
!                       renamed from MAOAM to SF6
!   (2017: Version 2.0: changed most things to channel/object structure
!---------------------------------------------------------------------------

MODULE messy_sf6_e5

  ! ECHAM5/MESSy
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_mpi_bi,       ONLY: message
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  ! MESSy
  USE messy_SF6

  
  IMPLICIT NONE
  
  !-----------------------------------------------------------------
  ! Everything is PRIVATE, except when explicitely stated otherwise:
  !-----------------------------------------------------------------
  PRIVATE
  SAVE ! op_pj_20180315

  INTRINSIC NULL
  
  REAL(dp), DIMENSION(:,:,:), POINTER :: &
    SF6e_gas  => NULL(),&
    SF6l_gas => NULL(),&
    SF6ea_gas  => NULL(),&
    SF6la_gas => NULL(),&
    Age_of_Air_l  => NULL(), &
    Age_of_Air_la  => NULL(), &
    rk1 => NULL(), &
    rk2 => NULL(), &
    rk3 => NULL(), &
    rk4 => NULL(), &
    rk5 => NULL(), &
    rk6 => NULL(), &
    rk7 => NULL(), &
    rk8 => NULL(), &
    efield => NULL(), &
    RADO3_O3 => NULL(), &
    SF6O2_O2 => NULL(), &
    SF6O3P_O3P => NULL(), &
    SF6N2_N2 => NULL(), &
    SF6H_H => NULL(), &
    SF6HCl_HCl => NULL(), &
    SF6_pre => NULL(), &
    SF6_aoa_pre => NULL()
        
    
REAL(dp), DIMENSION(:,:,:), POINTER :: &
    lt_sf6_e => NULL()
    

  REAL(dp), DIMENSION(:,:), POINTER :: &
    cossza_2d => NULL()
  
    INTEGER :: idt_O3, idt_O1D, idt_O3P, idt_O2, idt_N2, idt_H, idt_HCl
    INTEGER :: idt_SF6e, idt_SF6ea, idt_SF6l, idt_SF6la
    
  !-----
  ! switch for tracer initialisation
  !-----
  LOGICAL :: tracer_init_required = .false.

  !-----
  ! 2d-field for tropopause index (details via namelist)
  !-----
  REAL(dp), DIMENSION(:,:), POINTER :: tp_i0 => NULL()


  TYPE(t_chaobj_cpl) :: sf6_tropop, sf6_o3, sf6_o3p, sf6_o2, sf6_n2, sf6_h, sf6_hcl, sf6_sf6, sf6_sf6_aoa
  
  PUBLIC :: SF6_initialize    
  PUBLIC :: sf6_new_tracer   
  PUBLIC :: sf6_init_memory 
  PUBLIC :: sf6_init_coupling
  PUBLIC :: sf6_init_tracer 
  PUBLIC :: sf6_physc      

CONTAINS

!=============================================================================

  SUBROUTINE sf6_initialize

    ! ECHAM5/MESSy
    USE messy_main_tools ,    ONLY: find_next_free_unit
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast, finish
    ! MESSy
    USE messy_SF6,          ONLY: sf6_read_nml_ctrl, SF6_tracer, nml_a_t0, nml_b, nml_t0, &
         tp_gap, emode, &
         autocoeff, idontcalc

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: substr = 'sf6_initialize'
    INTEGER                     :: status, iou

    !-----
    ! read namelist CTRL
    !-----
     IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL sf6_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr,'call sf6_read_nml_ctrl failed')
     END IF
     CALL p_bcast(nml_a_t0, p_io)
     CALL p_bcast(nml_b, p_io)
     CALL p_bcast(nml_t0, p_io)
     CALL p_bcast(SF6_tracer, p_io)
     CALL p_bcast(tp_gap, p_io)
     CALL p_bcast(emode, p_io)
     CALL p_bcast(autocoeff, p_io)
     CALL p_bcast(idontcalc, p_io)
     
     !-----
    ! read namelist CPL
    !-----
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL sf6_read_nml_cpl(status, iou)   ! read /CPL/
      IF (status /= 0) CALL finish(substr,'call sf6_read_nml_cpl failed')
    END IF
     CALL p_bcast(sf6_h%cha, p_io)
     CALL p_bcast(sf6_h%obj, p_io)
     CALL p_bcast(sf6_hcl%cha, p_io)
     CALL p_bcast(sf6_hcl%obj, p_io)
     CALL p_bcast(sf6_n2%cha, p_io)
     CALL p_bcast(sf6_n2%obj, p_io)
     CALL p_bcast(sf6_o2%cha, p_io)
     CALL p_bcast(sf6_o2%obj, p_io)
     CALL p_bcast(sf6_o3%cha, p_io)
     CALL p_bcast(sf6_o3%obj, p_io)
     CALL p_bcast(sf6_o3p%cha, p_io)
     CALL p_bcast(sf6_o3p%obj, p_io)
     IF (SF6_tracer.eq.1) THEN
       CALL p_bcast(sf6_sf6%cha, p_io)
       CALL p_bcast(sf6_sf6%obj, p_io)
     END IF
     CALL p_bcast(sf6_sf6_aoa%cha, p_io)
     CALL p_bcast(sf6_sf6_aoa%obj, p_io)
     CALL p_bcast(sf6_tropop%cha, p_io)
     CALL p_bcast(sf6_tropop%obj, p_io)
     
 
  END SUBROUTINE sf6_initialize

!=============================================================================
 
   SUBROUTINE sf6_new_tracer
! ---------------------------------------------------------------
     USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
     USE messy_main_tracer_tools_bi, ONLY: tracer_halt
     USE messy_main_blather_bi,      ONLY: start_message_bi
     USE messy_main_tracer,          ONLY: new_tracer, set_tracer, get_tracer &
                                          , AIR   &
                                          , AMOUNTFRACTION      &
                                          , OFF, ON                     &
                                          , R_MOLARMASS     &
                                          , R_dryreac_sf                       &
                                          , I_ADVECT, I_CONVECT, I_VDIFF   &
                                          , I_SCAV

     IMPLICIT NONE
 
     CHARACTER(LEN=*), PARAMETER :: substr='sf6_new_tracer'
     INTEGER :: i_err
     INTEGER :: status

     CALL start_message_bi(modstr,'TRACER REQUEST',substr)

       CALL new_tracer(status, GPTRSTR, 'SF6e', modstr      &
            , idx=idt_sf6e , unit='mol/mol', medium=AIR    &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6e, R_molarmass,  r=146.07_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6e, I_ADVECT, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6e, I_CONVECT,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6e, I_VDIFF,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6e, I_SCAV,   OFF)
       CALL tracer_halt(substr, status)
       
       
       CALL new_tracer(status, GPTRSTR, 'SF6ea', modstr      &
            , idx=idt_sf6ea , unit='mol/mol', medium=AIR    &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6ea, R_molarmass,  r=146.07_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6ea, I_ADVECT, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6ea, I_CONVECT,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6ea, I_VDIFF,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6ea, I_SCAV,   OFF)
       CALL tracer_halt(substr, status)

       CALL new_tracer(status, GPTRSTR, 'SF6la', modstr      &
            , idx=idt_sf6la , unit='mol/mol', medium=AIR    &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6la, R_molarmass,  r=146.07_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6la, I_ADVECT, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6la, I_CONVECT,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6la, I_VDIFF,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6la, I_SCAV,   OFF)
       CALL tracer_halt(substr, status)

       CALL new_tracer(status, GPTRSTR, 'SF6l', modstr      &
            , idx=idt_sf6l , unit='mol/mol', medium=AIR    &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6l, R_molarmass,  r=146.07_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6l, I_ADVECT, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6l, I_CONVECT,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6l, I_VDIFF,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_sf6l, I_SCAV,   OFF)
       CALL tracer_halt(substr, status)

   END SUBROUTINE sf6_new_tracer
 
!=============================================================================
 
   SUBROUTINE sf6_init_memory
! ---------------------------------------------------------------
 
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
     IMPLICIT NONE
     INTEGER :: status
     CHARACTER(LEN=*), PARAMETER           :: substr = 'sf6_init_memory'
 
 
     CALL message('sf6_init_memory','defining streams for SF6')
     ! -----
     ! Define a new output stream
     ! -----
!     CALL new_stream (maoam,modstr)
!     CALL default_stream_setting (maoam, lpost =.FALSE., &
!          contnorest=.TRUE.) ! mz_pj_20050805
 
     ! -----
     ! maoam specific molecule information
     ! -----
     CALL new_channel(status, modstr, reprid=GP_3D_MID)
     CALL channel_halt(substr, status)

     IF (SF6_tracer.eq.1) THEN

       CALL new_channel_object(status, modstr, 'SF6e' &
          , p3 = SF6e_gas)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'SF6e' &
          , 'long_name', c='SF6 vmr from external emissions')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'SF6e' &
          , 'units', c='mol/mol')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr, 'SF6ea' &
          , p3 = SF6ea_gas)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'SF6ea' &
          , 'long_name', c='SF6 vmr from external emissions with mesospheric sink')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'SF6ea' &
          , 'units', c='mol/mol')
       CALL channel_halt(substr, status)
     END IF


     CALL new_channel_object(status, modstr, 'SF6l' &
          , p3 = SF6l_gas)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'SF6l' &
          , 'long_name', c='SF6 linear tracer')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'SF6l' &
          , 'units', c='mol/mol')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'SF6la' &
          , p3 = SF6la_gas)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'SF6la' &
          , 'long_name', c='SF6 linear tracer with mesospheric sink')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'SF6la' &
          , 'units', c='mol/mol')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'Age_of_Air_l' &
          , p3 = Age_of_Air_l)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'Age_of_Air_l' &
          , 'long_name', c='Age of Air from SF6l')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'Age_of_Air_l' &
          , 'units', c='years')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'Age_of_Air_la' &
          , p3 = Age_of_Air_la)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'Age_of_Air_la' &
          , 'long_name', c='Age of Air from SF6la')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'Age_of_Air_la' &
          , 'units', c='years')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'rk1' &
          , p3 = rk1)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk1' &
          , 'long_name', c='rk1 lambda(phot, sf6 -> sfi)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk1' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'rk2' &
          , p3 = rk2)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk2' &
          , 'long_name', c='rk2 lambda(e-,   sf6 -> sf6-*)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk2' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'rk3' &
          , p3 = rk3)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk3' &
          , 'long_name', c='rk3 lambda(photo, sf6- -> sf6)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk3' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'rk4' &
          , p3 = rk4)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk4' &
          , 'long_name', c='rk4 lambda(hydro, sf6- -> hf + ..)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk4' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'rk5' &
          , p3 = rk5)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk5' &
          , 'long_name', c='rk5 lambda(M, sf-* -> sf6)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk5' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'rk6' &
          , p3 = rk6)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk6' &
          , 'long_name', c='rk6 lambda(spont, sf6-* -> sf6)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk6' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'rk7' &
          , p3 = rk7)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk7' &
          , 'long_name', c='rk7 lambda(HCl, sf6- -> sf5 + ...)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk7' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'rk8' &
          , p3 = rk8)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk8' &
          , 'long_name', c='rk8 lambda(O3, sf6- -> sf6)')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'rk8' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'lt_sf6_e' &
          , p3 = lt_sf6_e)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'lt_sf6_e' &
          , 'long_name', c='Lifetime of SF6_e')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'lt_sf6_e' &
          , 'units', c='days')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'efield' &
          , p3 = efield)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'efield' &
          , 'long_name', c='electron density')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'efield' &
          , 'units', c='#/cm^3')
     CALL channel_halt(substr, status)

  END SUBROUTINE sf6_init_memory


  SUBROUTINE sf6_init_coupling

    USE messy_main_mpi_bi,           ONLY: finish, message, p_parallel_io

    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_channel_error_bi, ONLY: channel_halt 
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: get_tracer
    USE messy_main_channel,       ONLY: get_channel_object, get_channel_info &
                                      , new_channel_object  & 
                                      , new_attribute 

    IMPLICIT NONE

    INTEGER :: status, i_err
    CHARACTER(LEN=*), PARAMETER :: substr = 'sf6_init_coupling'
    
    CALL start_message_bi(modstr,'INIT COUPLING',substr)

   IF (p_parallel_io) THEN
       WRITE(*,*) 'tropopause channel: ', TRIM(sf6_tropop%cha)
       WRITE(*,*) 'channel object containing tropopause index: ', &
            TRIM(sf6_tropop%obj)
    END IF

    CALL get_channel_object(status, &
         TRIM(sf6_tropop%cha), TRIM(sf6_tropop%obj), p2=tp_i0 )
    IF (status /= 0) &
         CALL finish(substr,'channel object for tropopause index not found')

    CALL get_channel_object(status, &
         TRIM(sf6_o3%cha), TRIM(sf6_o3%obj), p3=RADO3_O3 )
    IF (status /= 0) &
         CALL finish(substr,'channel object for O3 not found')
         
    CALL get_channel_object(status, &
         TRIM(sf6_o2%cha), TRIM(sf6_o2%obj), p3=SF6O2_O2 )
    IF (status /= 0) &
         CALL finish(substr,'channel object for O2 not found')
         
    CALL get_channel_object(status, &
         TRIM(sf6_o3p%cha), TRIM(sf6_o3p%obj), p3=SF6O3P_O3P )
    IF (status /= 0) &
         CALL finish(substr,'channel object for O3P not found')

    CALL get_channel_object(status, &
         TRIM(sf6_n2%cha), TRIM(sf6_n2%obj), p3=SF6N2_N2 )
    IF (status /= 0) &
         CALL finish(substr,'channel object for N2 not found')

    CALL get_channel_object(status, &
         TRIM(sf6_h%cha), TRIM(sf6_h%obj), p3=SF6H_H )
    IF (status /= 0) &
         CALL finish(substr,'channel object for H not found')

    CALL get_channel_object(status, &
         TRIM(sf6_hcl%cha), TRIM(sf6_hcl%obj), p3=SF6HCl_HCl )
    IF (status /= 0) &
         CALL finish(substr,'channel object for HCl not found')


    IF (SF6_tracer.eq.1) THEN
      CALL get_channel_object(status, &
         TRIM(sf6_sf6%cha), TRIM(sf6_sf6%obj), p3=SF6_pre )
      IF (status /= 0) &
         CALL finish(substr,'channel object for SF6 not found')
   ENDIF

      CALL get_channel_object(status, &
         TRIM(sf6_sf6_aoa%cha), TRIM(sf6_sf6_aoa%obj), p3=SF6_aoa_pre )
      IF (status /= 0) &
         CALL finish(substr,'channel object for SF6 aoa_pre not found')
 

    CALL get_channel_object(status, 'orbit','cossza',p2=cossza_2d)
    CALL channel_halt(substr//': channel object for cossza in orbit not found!',status)
    
    CALL end_message_bi(modstr,'INIT COUPLING',substr)

  END SUBROUTINE sf6_init_coupling
 
! !=============================================================================
 
SUBROUTINE sf6_init_tracer
 
!!   USE messy_main_tracer_bi, ONLY: main_tracer_initialize
 
   IMPLICIT NONE
 
!!   IF (tracer_init_required) CALL tracer_init(modstr)
 
END SUBROUTINE sf6_init_tracer

! !=============================================================================
 
 SUBROUTINE sf6_physc
  !--------------------------------------------------------------------------
  USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1
  USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma, nproma
  USE messy_main_grid_def_bi,     ONLY:philat_2d
  USE messy_main_data_bi,       ONLY: pmid => press_3d, temp => t_scb & ! mid-level pressures [Pa]
                                    , pint => pressi_3d  ! level interface (above) pressure [Pa]
  USE messy_main_timer,         ONLY: YEAR, DAYOFYEAR, HOUR, MINUTE, SECOND &
                                   , current_time_step, time_step_len
  USE messy_main_constants_mem, ONLY: M_air, R_gas, N_A, g
  USE messy_main_tools,         ONLY: Spline1D, Splint1D
!csv  USE messy_rad_e5,      ONLY: l_trigrad


  IMPLICIT none

  INTEGER :: jp, jk, n,  i, status
  REAL(dp), DIMENSION(16) :: hydro, zhydro
  REAL(dp), DIMENSION(2*nlev) :: sbuf
  REAL(dp), DIMENSION(nlev) :: profile

  REAL(dp), DIMENSION(nproma,nlev) :: SF6e_t0, SF6ea_t0
  REAL(dp), DIMENSION(nproma,nlev) :: SF6l_t0, SF6la_t0
  REAL(dp), DIMENSION(nproma,nlev) :: SF6e_g, SF6ea_g
  REAL(dp), DIMENSION(nproma,nlev) :: SF6l_g, SF6la_g
  REAL(dp), DIMENSION(nproma,nlev) :: Age_of_Air_l_g, Age_of_Air_la_g

  REAL(dp), DIMENSION(nproma,nlev) :: O3_t0, O_t0, O2_t0, N2_t0, H_t0, HCl_t0
  REAL(dp), DIMENSION(nproma,nlev) :: cair,c_o, c_o2, c_o3, c_n2, col_o, col_o2, col_n2, c_h, c_hcl,&
           col_tot, col_o_box, col_o2_box, col_n2_box
  REAL(dp), DIMENSION(nproma,nlev) :: grheight
  REAL(dp) :: act_time
  REAL(DP),PARAMETER               :: scale_height_ma = 7000._dp ! Middle atmosphere scale height [m]
  REAL(dp),  DIMENSION(:,:), POINTER    :: altitude  ! altitude [m]
  REAL(DP):: lamsf6 
  REAL(DP) :: SF6_lin
  REAL(DP) :: epsilon
  REAL(DP) :: end_time, end_time2, akt_vmr, old_vmr, trend
  REAL(dp), DIMENSION(nproma)      :: tp_i
  REAL(DP)                    :: zpb, zpt
  REAL(DP) :: rk(8), rk1_g(nproma,nlev), rk2_g(nproma,nlev),rk3_g(nproma,nlev),rk4_g(nproma,nlev), &
         rk5_g(nproma,nlev),rk6_g(nproma,nlev),rk7_g(nproma,nlev),rk8_g(nproma,nlev)
  REAL(dp) :: lt_sf6_e_g(nproma,nlev)

  REAL(DP) :: efield_g(nproma,nlev)
  REAL(DP) :: efield_loc
  INTEGER  :: jtpi ! op_pj_20180315

! profile from Brasseur for lev > 2; level 1 and 2 to avoid instabilities due to spline interpolation
!  hydro(1) = 1.0e-22
!  hydro(2) = 1.0e-20
!  hydro(3) = 3.2e-16
!  hydro(4) = 1.6e-14
!  hydro(5) = 4.2e-13
!  hydro(6) = 4.1e-12
!  hydro(7) = 2.2e-11
!  hydro(8) = 9.7e-11
!  hydro(9) = 4.5e-10
!  hydro(10) = 2.5e-9
!  hydro(11) = 2.0e-8
!  hydro(12) = 2.0e-7
!  hydro(13) = 4.6e-6
!  hydro(14) = 7.2e-6
!  hydro(15) = 9.0e-6
!  hydro(16) = 1.0e-5

!  zhydro(1) = 10.
!  zhydro(2) = 15000.
!  do jk=3,16
!    zhydro(jk)=(30+(jk-2)*5)*1.e3
!  end do
  
  tp_i(1:kproma) = 0.

  IF (ASSOCIATED(tp_i0)) tp_i(1:kproma) &
              = tp_i0(1:kproma,jrow)

! ACHTUNG!!! MUSS NOCH AUF SCHALTJAHRE ANGEPASST WERDEN
  act_time=YEAR+(DAYOFYEAR+(HOUR+(MINUTE+SECOND/60.)/60.)/24.)/365.

  ALLOCATE(altitude(kproma,nlev))

    ! CALCULATE ALTITUDE in m
  altitude(1:kproma,:)=-scale_height_ma*LOG(pmid(1:kproma,:,jrow)/1E5_dp)
 
! if (type_SF6_tracer.eq.1) then
!  SF6_akt=((13.053*log(exp((act_time-1984.93)/13.053) + exp(-(act_time-1984.93)/13.053)) &
!                        + (act_time-1984.93))*0.13376+0.03)/1.d12
! endif
 
! -------------------
! Linear tracer
! always switched on
! ------------------- 

  SF6_lin=(act_time-nml_t0)*nml_b+nml_a_t0
 level_loop0: DO jk=1,nlev
    vector_loop0: DO jp=1,kproma
      SF6l_t0(jp,jk) = pxtm1(jp,jk,idt_SF6l) &
                        +pxtte(jp,jk,idt_SF6l)*time_step_len
      SF6la_t0(jp,jk) = pxtm1(jp,jk,idt_SF6la) &
                        +pxtte(jp,jk,idt_SF6la)*time_step_len
 
    END DO vector_loop0
  END DO level_loop0
  
      vector_loop1: DO jp=1,kproma
       if (tp_gap.lt.tp_i(jp)) then
!!$          SF6l_g(jp,tp_i(jp)+tp_gap:nlev) = SF6_lin
!!$          SF6la_g(jp,tp_i(jp)+tp_gap:nlev) = SF6_lin
! op_pj_20180315+: Array index must be of INTEGER type
!!$          SF6l_g(jp,tp_i(jp)+tp_gap:nlev) = SF6_aoa_pre(jp,tp_i(jp)+tp_gap:nlev,jrow)
!!$          SF6la_g(jp,tp_i(jp)+tp_gap:nlev) = SF6_aoa_pre(jp,tp_i(jp)+tp_gap:nlev,jrow)
!!$          SF6l_g(jp,1:tp_i(jp)-1+tp_gap) = SF6l_t0(jp,1:tp_i(jp)-1+tp_gap)
!!$          SF6la_g(jp,1:tp_i(jp)-1+tp_gap) = SF6la_t0(jp,1:tp_i(jp)-1+tp_gap)
          jtpi = NINT(tp_i(jp))
          SF6l_g(jp,jtpi+tp_gap:nlev) = SF6_aoa_pre(jp,jtpi+tp_gap:nlev,jrow)
          SF6la_g(jp,jtpi+tp_gap:nlev) = SF6_aoa_pre(jp,jtpi+tp_gap:nlev,jrow)
          SF6l_g(jp,1:jtpi-1+tp_gap) = SF6l_t0(jp,1:jtpi-1+tp_gap)
          SF6la_g(jp,1:jtpi-1+tp_gap) = SF6la_t0(jp,1:jtpi-1+tp_gap)
! op_pj_20180315-
       else
!!$          SF6l_g(jp,nlev) = SF6_lin
!!$          SF6la_g(jp,nlev) = SF6_lin
          SF6l_g(jp,nlev) = SF6_aoa_pre(jp,nlev,jrow)
          SF6la_g(jp,nlev) = SF6_aoa_pre(jp,nlev,jrow)
          SF6l_g(jp,1:nlev-1) = SF6l_t0(jp,1:nlev-1)
          SF6la_g(jp,1:nlev-1) = SF6la_t0(jp,1:nlev-1)
       endif
    END DO vector_loop1
 
 
 level_loop1c: DO jk=1,nlev
    vector_loop1c: DO jp=1,kproma

! get current vmrs
        O3_t0(jp,jk) = RADO3_O3(jp,jk,jrow)       ! ccvo3 in KASIMA
        O_t0(jp,jk) = SF6O3P_O3P(jp,jk,jrow)
        O2_t0(jp,jk) = SF6O2_O2(jp,jk,jrow)
        N2_t0(jp,jk) = SF6N2_N2(jp,jk,jrow)
        H_t0(jp,jk) = SF6H_H(jp,jk,jrow)
        HCl_t0(jp,jk) = SF6HCl_HCl(jp,jk,jrow)
   END DO vector_loop1c
  END DO level_loop1c
! --------------------
! external SF6 tracer
! only if SF6_tracer=1
! --------------------
 
 if (SF6_tracer.eq.1) then
  level_loop0b: DO jk=1,nlev
    vector_loop0b: DO jp=1,kproma
      SF6e_t0(jp,jk) = pxtm1(jp,jk,idt_SF6e) &
                        +pxtte(jp,jk,idt_SF6e)*time_step_len
      SF6ea_t0(jp,jk) = pxtm1(jp,jk,idt_SF6ea) &
                        +pxtte(jp,jk,idt_SF6ea)*time_step_len
 
    END DO vector_loop0b
  END DO level_loop0b
  
  vector_loop1b: DO jp=1,kproma
       if (tp_gap.lt.tp_i(jp)) then
! op_pj_20180315+: Array index at (1) must be of INTEGER type
!!$          SF6e_g(jp,tp_i(jp)+tp_gap:nlev) = SF6_pre(jp,tp_i(jp)+tp_gap:nlev,jrow)
!!$          SF6ea_g(jp,tp_i(jp)+tp_gap:nlev) = SF6_pre(jp,tp_i(jp)+tp_gap:nlev,jrow)
!!$          SF6e_g(jp,1:tp_i(jp)-1+tp_gap) = SF6e_t0(jp,1:tp_i(jp)-1+tp_gap)
!!$          SF6ea_g(jp,1:tp_i(jp)-1+tp_gap) = SF6ea_t0(jp,1:tp_i(jp)-1+tp_gap)
          jtpi = NINT(tp_i(jp))
          SF6e_g(jp,jtpi+tp_gap:nlev) = SF6_pre(jp,jtpi+tp_gap:nlev,jrow)
          SF6ea_g(jp,jtpi+tp_gap:nlev) = SF6_pre(jp,jtpi+tp_gap:nlev,jrow)
          SF6e_g(jp,1:jtpi-1+tp_gap) = SF6e_t0(jp,1:jtpi-1+tp_gap)
          SF6ea_g(jp,1:jtpi-1+tp_gap) = SF6ea_t0(jp,1:jtpi-1+tp_gap)
! op_pj_20180315-
       else
          SF6e_g(jp,nlev) = SF6_pre(jp,nlev,jrow)
          SF6ea_g(jp,nlev) = SF6_pre(jp,nlev,jrow)
          SF6e_g(jp,1:nlev-1) = SF6e_t0(jp,1:nlev-1)
          SF6ea_g(jp,1:nlev-1) = SF6ea_t0(jp,1:nlev-1)
       endif
    END DO vector_loop1b
  
 END IF ! SF6_tracer=1


! calculate gas densitys
level_loop1d: DO jk=1,nlev
    vector_loop1d: DO jp=1,kproma
        cair(jp,jk)  = (N_A) * pmid(jp,jk,jrow) / (R_gas*temp(jp,jk,jrow))   ! m^(-3)
        c_o(jp,jk) = cair(jp,jk) * O_t0(jp,jk)
        c_o2(jp,jk) = cair(jp,jk) * O2_t0(jp,jk)
        c_o3(jp,jk) = cair(jp,jk) * O3_t0(jp,jk)
        c_n2(jp,jk) = cair(jp,jk) * N2_t0(jp,jk)
        c_h(jp,jk) =cair(jp,jk) * H_t0(jp,jk)
        c_hcl(jp,jk) =cair(jp,jk) * HCl_t0(jp,jk)
!calculate height of box
        zpb = pint(jp,jk+1, jrow)
          ! ECHAM5 top layer ends at 0. Pa !!! adjust to mid of uppermost level
        zpt = MAX(pint(jp, jk, jrow),pmid(jp,1,jrow))
        grheight(jp,jk) = (1000._dp * R_gas / (M_air * g)) &
               * temp(jp,jk,jrow) * log(zpb/zpt)     ! in m
        col_o_box(jp,jk)=c_o(jp,jk)*grheight(jp,jk)
        col_o2_box(jp,jk)=c_o2(jp,jk)*grheight(jp,jk)
        col_n2_box(jp,jk)=c_n2(jp,jk)*grheight(jp,jk)
        col_o(jp,jk)=0.
        col_o2(jp,jk)=0.
        col_n2(jp,jk)=0.
   END DO vector_loop1d
  END DO level_loop1d
!write(*,*) 'calculate gas densitys'

    DO jp=1,kproma
      col_o(jp,1)=col_o_box(jp,1)
      col_o2(jp,1)=col_o2_box(jp,1)
      col_n2(jp,1)=col_n2_box(jp,1)
    END DO

! calculate column above level
level_loop1e: DO jk=2,nlev
    vector_loop1e: DO jp=1,kproma
      DO i=jk-1,jk
        col_o(jp,jk)=col_o(jp,jk)+col_o_box(jp,i)
        col_o2(jp,jk)=col_o2(jp,jk)+col_o2_box(jp,i)
        col_n2(jp,jk)=col_n2(jp,jk)+col_n2_box(jp,i)
      END DO
       
   END DO vector_loop1e
  END DO level_loop1e
! convert m^2 to cm^2
  col_o=col_o*1.d-04
  col_o2=col_o2*1.d-04
  col_n2=col_n2*1.d-04

  col_tot=col_o+col_o2+col_n2
  col_tot=col_tot*1.d-17          ! magic factor from KASIMA ???

level_loop2: DO jk=1,nlev
    vector_loop2: DO jp=1,kproma


        call sf6life(lamsf6, emode, cossza_2d(jp,jrow), time_step_len, altitude(jp,jk), c_o3(jp,jk),philat_2d(jp,jrow), &
               col_o2(jp,jk), c_h(jp,jk), c_hcl(jp,jk), cair(jp,jk), rk,efield_loc, DAYOFYEAR) ! op_pj_20180418 DAYOFYEAR added
!!  rk is only stored for melec=11;
        rk1_g(jp,jk)=rk(1)
        rk2_g(jp,jk)=rk(2)
        rk3_g(jp,jk)=rk(3)
        rk4_g(jp,jk)=rk(4)
        rk5_g(jp,jk)=rk(5)
        rk6_g(jp,jk)=rk(6)
        rk7_g(jp,jk)=rk(7)
        rk8_g(jp,jk)=rk(8)
        lt_sf6_e_g(jp,jk)=lamsf6
        efield_g(jp,jk)=efield_loc

        !SF6ea_g(jp,jk)    = max(0._dp,SF6ea_g(jp,jk)-SF6ea_g(jp,jk)*lamsf6)
        SF6la_g(jp,jk)   = max(0._dp,SF6la_g(jp,jk)-SF6la_g(jp,jk)*lamsf6)
        !SF6l_g(jp,jk)   = SF6l_g(jp,jk)
        !SF6e_g(jp,jk)    = SF6l_g(jp,jk)
        

! calculate mean age of air - see messy_SF6.f90


        call calc_age(Age_of_Air_l_g(jp,jk),SF6l_g(jp,jk),nml_a_t0,nml_b,act_time,nml_t0)
        call calc_age(Age_of_Air_la_g(jp,jk),SF6la_g(jp,jk),nml_a_t0,nml_b,act_time,nml_t0)

        !pxtte(jp,jk,idt_SF6e)= pxtte(jp,jk,idt_SF6e)             & 
        !  + (SF6e_g(jp,jk)-SF6e_t0(jp,jk))                       &
        !   /time_step_len
        !pxtte(jp,jk,idt_SF6ea)= pxtte(jp,jk,idt_SF6ea)             & 
        !  + (SF6ea_g(jp,jk)-SF6ea_t0(jp,jk))                       &
        !   /time_step_len
        pxtte(jp,jk,idt_SF6l)= pxtte(jp,jk,idt_SF6l)           &
          + (SF6l_g(jp,jk)-SF6l_t0(jp,jk))                     &
           /time_step_len
        pxtte(jp,jk,idt_SF6la)= pxtte(jp,jk,idt_SF6la)           &
          + (SF6la_g(jp,jk)-SF6la_t0(jp,jk))                     &
           /time_step_len
       
!        SF6e_gas(jp,jk,jrow) = SF6e_g(jp,jk)
        SF6l_gas(jp,jk,jrow)  = SF6l_g(jp,jk)
        !SF6ea_gas(jp,jk,jrow) = SF6ea_g(jp,jk)
        SF6la_gas(jp,jk,jrow)  = SF6la_g(jp,jk)

        ! overwrite first few timesteps; specified in SF6 namelist
        if (current_time_step.le.idontcalc) then
          Age_of_Air_l_g(jp,jk) = 9999.
          Age_of_Air_la_g(jp,jk) = 9999.
        end if
           

        Age_of_Air_l(jp,jk,jrow)  = Age_of_Air_l_g(jp,jk)
        Age_of_Air_la(jp,jk,jrow)  = Age_of_Air_la_g(jp,jk)

        rk1(jp,jk,jrow)=rk1_g(jp,jk)
        rk2(jp,jk,jrow)=rk2_g(jp,jk)
        rk3(jp,jk,jrow)=rk3_g(jp,jk)
        rk4(jp,jk,jrow)=rk4_g(jp,jk)
        rk5(jp,jk,jrow)=rk5_g(jp,jk)
        rk6(jp,jk,jrow)=rk6_g(jp,jk)
        rk7(jp,jk,jrow)=rk7_g(jp,jk)
        rk8(jp,jk,jrow)=rk8_g(jp,jk)

        lt_sf6_e(jp,jk,jrow)=real(1./(lt_sf6_e_g(jp,jk)/time_step_len)/86400.)

        lt_sf6_e(jp,jk,jrow)=min(max(lt_sf6_e(jp,jk,jrow),0.d0),1.d37)

        efield(jp,jk,jrow)=efield_g(jp,jk)
        
        IF (SF6_tracer.eq.1) THEN
           SF6ea_g(jp,jk)    = max(0._dp,SF6ea_g(jp,jk)-SF6ea_g(jp,jk)*lamsf6)
           !SF6e_g(jp,jk)    = SF6l_g(jp,jk)
           pxtte(jp,jk,idt_SF6e)= pxtte(jp,jk,idt_SF6e)             & 
          + (SF6e_g(jp,jk)-SF6e_t0(jp,jk))                       &
           /time_step_len
           pxtte(jp,jk,idt_SF6ea)= pxtte(jp,jk,idt_SF6ea)             & 
          + (SF6ea_g(jp,jk)-SF6ea_t0(jp,jk))                       &
           /time_step_len
           SF6e_gas(jp,jk,jrow) = SF6e_g(jp,jk)
           SF6ea_gas(jp,jk,jrow) = SF6ea_g(jp,jk)
        END IF

    END DO vector_loop2
  END DO level_loop2

!calculations for artificial tracers
!        pxtm1(1:kproma,1,idt_trope90) = 0._dp
!   DO jk=1,nlev
!      pxtte(1:kproma,jk,idt_trope90) =  pxtte(1:kproma,jk,idt_trope90) + (pxtm1(1:kproma,jk,idt_trope90)*exp(-(time_step_len)/7776000._dp)-pxtm1(1:kproma,jk,idt_trope90))/time_step_len
!   END DO
!   DO jp=1,kproma
!     DO jk=1,tp_i(jp)-1
     ! 70 Jahre a 365 Tage
!       pxtte(jp,jk,idt_acetonitril) =  pxtte(jp,jk,idt_acetonitril) + (pxtm1(jp,jk,idt_acetonitril)*exp(-(time_step_len)/2207520000._dp)-pxtm1(jp,jk,idt_acetonitril))/time_step_len
!       pxtte(jp,jk,idt_acetonitril0) =  pxtte(jp,jk,idt_acetonitril0) + (pxtm1(jp,jk,idt_acetonitril0)*exp(-(time_step_len)/2207520000._dp)-pxtm1(jp,jk,idt_acetonitril0))/time_step_len
!     END DO
!     DO jk=tp_i(jp),nlev
     ! 228 Tage = 7.5 Monate
!      pxtte(jp,jk,idt_acetonitril) =  pxtte(jp,jk,idt_acetonitril) + (pxtm1(jp,jk,idt_acetonitril)*exp(-(time_step_len)/19699200._dp)-pxtm1(jp,jk,idt_acetonitril))/time_step_len
!       pxtte(jp,jk,idt_acetonitril0) =  pxtte(jp,jk,idt_acetonitril0) + (pxtm1(jp,jk,idt_acetonitril0)*exp(-(time_step_len)/19699200._dp)-pxtm1(jp,jk,idt_acetonitril0))/time_step_len
!     END DO
!   END DO
   
     
       

  DEALLOCATE(altitude)

END SUBROUTINE sf6_physc

!=============================================================================

!SUBROUTINE SF6_local_end

!END SUBROUTINE SF6_local_end

!=============================================================================

SUBROUTINE sf6_read_nml_cpl(status, iou)
  !------------------------------------------
  ! read coupling namelist /CPL/ from SF6.nml
  !------------------------------------------

  ! MESSy
  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  
  IMPLICIT NONE
  
  INTEGER, INTENT(out) :: status
  INTEGER, INTENT(in) :: iou
  
  CHARACTER(len=*), PARAMETER :: substr='sf6_read_nml_cpl'
  LOGICAL :: lex     ! file exists?
  INTEGER :: fstat   ! file status
  
  NAMELIST /CPL/ sf6_tropop, sf6_o3, sf6_sf6, sf6_sf6_aoa, sf6_o3p, sf6_o2, sf6_n2, sf6_h, sf6_hcl
    
  status = 1             ! initialise status flag with error code
  
  CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
  IF (.not.lex) RETURN   ! error: psc.nml does not exist
  read(iou, nml=cpl, iostat=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
  IF (fstat/=0) RETURN   ! error while reading namelist
  

  CALL read_nml_close(substr, iou, modstr)
  status = 0   ! no error
  
END SUBROUTINE sf6_read_nml_cpl

!=============================================================================

END MODULE messy_sf6_e5
