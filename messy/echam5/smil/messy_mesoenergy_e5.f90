! DESCRIPTION:
! This module is the interface (called by messy_main_control) between the 
! ECHAM5 base model and the mesoenergy module.
!
! AUTHOR:
! S. Versick, IMK-ASF, KIT, Germany

!---------------------------------------------------------------------------

MODULE messy_mesoenergy_e5

  ! ECHAM5/MESSy
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_mpi_bi,       ONLY: message, dcl
  ! MESSy
  USE messy_mesoenergy

  
  IMPLICIT NONE

  !-----------------------------------------------------------------
  ! Everything is PRIVATE, except when explicitely stated otherwise:
  !-----------------------------------------------------------------
  PRIVATE

  INTRINSIC NULL
  
  REAL(dp), DIMENSION(:,:,:), POINTER :: &
    rea_g2101  => NULL(),  &
    HR2101 => NULL(),  &
    CHEMHEATING => NULL(), &
    actloss => NULL(),  &
    oldloss => NULL(), &
    actag96 => NULL(), &
    oldag96 => NULL(), &
    actag85 => NULL(), &
    oldag85 => NULL(), &
    actag83 => NULL(), &
    oldag83 => NULL(), &
    actag62 => NULL(), &
    oldag62 => NULL(), &
    actag31 => NULL(), &
    oldag31 => NULL(), &
    actago21d => NULL(), &
    oldago21d => NULL(), &
    actago21s => NULL(), &
    oldago21s => NULL(), &
    actago1so3p => NULL(), &
    oldago1so3p => NULL(), &
    actago1so1d => NULL(), &
    oldago1so1d => NULL(), &
    actago1do3p => NULL(), &
    oldago1do3p => NULL(), &
    OHv0nd => NULL(), &
    OHv1nd => NULL(), &   !! total number density
    OHv2nd => NULL(), &
    OHv3nd => NULL(), &
    OHv4nd => NULL(), &
    OHv5nd => NULL(), &
    OHv6nd => NULL(), &
    OHv7nd => NULL(), &
    OHv8nd => NULL(), &
    OHv9nd => NULL(), &
    AG96 => NULL(), &
    AG85 => NULL(), &
    AG83 => NULL(), &
    AG62 => NULL(), &
    AG31 => NULL(), &
    AGO21d => NULL(), &
    AGO21s => NULL(), &
    AGO1SO3P => NULL(), &
    AGO1SO1D => NULL(), &
    AGO1DO3P => NULL(), &
    debug   => NULL()

  
  ! POINTER TO TEMPERATURE TENDENCY
  REAL(dp), DIMENSION(:,:,:), POINTER :: tte_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: lossg2101_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ag96_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ag85_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ag83_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ag62_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ag31_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ago21d_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ago21s_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ago1so3p_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ago1so1d_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ago1do3p_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: ohv0_ptr => NULL()

  !-----
  ! global variable for regrid events
  !-----
!!$  TYPE(RGTEVENT), DIMENSION(:), POINTER :: rgt ! op_pj_20180713 not used

  !-----
  ! switch for tracer initialisation
  !-----
  LOGICAL :: tracer_init_required = .false.


  CHARACTER (LEN=STRLEN_MEDIUM), PUBLIC        ::  &
       dummy= ''

  INTEGER :: idt_OHv0, idt_OHv1, idt_OHv2, idt_OHv3, idt_OHv4, idt_OHv5, idt_OHv6, idt_OHv7, idt_OHv8, idt_OHv9
  INTEGER :: idt_ag96, idt_ag85, idt_ag83, idt_ag62, idt_ag31, idt_ago21d, idt_ago21s, idt_ago1so3p, idt_ago1so1d, idt_ago1do3p
  INTEGER :: idt_lossg2101, idt_kJmol
  
  PUBLIC :: mesoenergy_initialize    
!!  PUBLIC :: mesoenergy_new_tracer   
  PUBLIC :: mesoenergy_init_memory 
  PUBLIC :: mesoenergy_init_coupling
!!  PUBLIC :: mesoenergy_init_tracer 
  PUBLIC :: mesoenergy_physc      
!!  PUBLIC :: idt_O2, idt_N2

CONTAINS

!=============================================================================

  SUBROUTINE mesoenergy_initialize

    ! ECHAM5/MESSy
    USE messy_main_tools ,    ONLY: find_next_free_unit
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast, finish
    ! MESSy
    USE messy_mesoenergy,          ONLY: mesoenergy_read_nml_ctrl, chemheat, rea_ent_ho3, fv, a8, a9, a10

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: substr = 'mesoenergy_initialize'
    INTEGER                     :: status, iou

    !-----
    ! read namelist CTRL
    !-----
     IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL mesoenergy_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr,'call mesoenergy_read_nml_ctrl failed')
     END IF
     CALL p_bcast(chemheat, p_io)
     CALL p_bcast(rea_ent_ho3, p_io)
     CALL p_bcast(fv, p_io)
     CALL p_bcast(a8, p_io)
     CALL p_bcast(a9, p_io)
     CALL p_bcast(a10, p_io)
    ! read namelist CPL
    !-----
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL mesoenergy_read_nml_cpl(status, iou)   ! read /CPL/
      IF (status /= 0) CALL finish(substr,'call mesoenergy_read_nml_cpl failed')
    END IF

 
  END SUBROUTINE mesoenergy_initialize

!=============================================================================
 
!!   SUBROUTINE mesoenergy_new_tracer
! ---------------------------------------------------------------
!!     USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
!!     USE messy_main_tracer_bi,     ONLY: tracer_halt
!!     USE messy_main_blather_bi,   ONLY: start_message_bi
!!     USE messy_main_tracer,        ONLY: new_tracer, set_tracer, get_tracer &
!!                                      , AIR   &
!!                                      , AMOUNTFRACTION      &
!!                                      , OFF, ON                     &
!!                                      , R_MOLARMASS     &
!!                                      , R_dryreac_sf                       &
!!                                      , I_ADVECT, I_CONVECT, I_VDIFF   &
!!                                      , I_SCAV

!!     IMPLICIT NONE
 
!!     CHARACTER(LEN=*), PARAMETER :: substr='mesoenergy_new_tracer'
!!     INTEGER :: i_err
!!     INTEGER :: status
!!     INTEGER :: idt_heat_h_o3, idt_lossg2101

!!     idt_heat_h_o3=0


!!     CALL start_message_bi(modstr,'TRACER REQUEST',substr)


!!      CALL get_tracer(i_err, GPTRSTR, 'LossG2101', idx=idt_lossg2101)
!!      CALL tracer_halt(substr, i_err)

!!   END SUBROUTINE mesoenergy_new_tracer
 
!=============================================================================
 
   SUBROUTINE mesoenergy_init_memory
! ---------------------------------------------------------------
 
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
     IMPLICIT NONE 
     INTEGER :: status
     CHARACTER(LEN=*), PARAMETER           :: substr = 'mesoenergy_init_memory'
 
 
     CALL message('mesoenergy_init_memory','defining streams for maoam')
     ! -----
     ! Define a new output stream
     ! -----
!     CALL new_stream (mesoenergy,modstr)
!     CALL default_stream_setting (mesoenergy, lpost =.FALSE., &
!          contnorest=.TRUE.) ! mz_pj_20050805
 
     ! -----
     ! mesoenergy specific molecule information
     ! -----
     CALL new_channel(status, modstr, reprid=GP_3D_MID)
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'HR2101' &
          , p3 = HR2101)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'HR2101' &
          , 'long_name', c='chemical heating rate from H + O3')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'HR2101' &
          , 'units', c='K/d')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'CHEMHEATING' &
          , p3 = CHEMHEATING)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'CHEMHEATING' &
          , 'long_name', c='chemical heating rate from KJmol')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'CHEMHEATING' &
          , 'units', c='K/d')
     CALL channel_halt(substr, status)


     CALL new_channel_object(status, modstr, 'actloss' &
          , p3 = actloss)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actloss' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actloss' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'oldloss' &
          , p3 = oldloss)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldloss' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldloss' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'actag96' &
          , p3 = actag96)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag96' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag96' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'oldag96' &
          , p3 = oldag96)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag96' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag96' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'actag85' &
          , p3 = actag85)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag85' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag85' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'oldag85' &
          , p3 = oldag85)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag85' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag85' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'actag83' &
          , p3 = actag83)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag83' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag83' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'oldag83' &
          , p3 = oldag83)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag83' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag83' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'actag62' &
          , p3 = actag62)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag62' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag62' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'oldag62' &
          , p3 = oldag62)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag62' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag62' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'actag31' &
          , p3 = actag31)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag31' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actag31' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'oldag31' &
          , p3 = oldag31)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag31' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldag31' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)


     CALL new_channel_object(status, modstr, 'actago21d' &
          , p3 = actago21d)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago21d' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago21d' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'oldago21d' &
          , p3 = oldago21d)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago21d' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago21d' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'actago21s' &
          , p3 = actago21s)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago21s' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago21s' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'oldago21s' &
          , p3 = oldago21s)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago21s' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago21s' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'actago1so3p' &
          , p3 = actago1so3p)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago1so3p' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago1so3p' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'oldago1so3p' &
          , p3 = oldago1so3p)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago1so3p' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago1so3p' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
          CALL new_channel_object(status, modstr, 'actago1so1d' &
          , p3 = actago1so1d)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago1so1d' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago1so1d' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'oldago1so1d' &
          , p3 = oldago1so1d)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago1so1d' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago1so1d' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'actago1do3p' &
          , p3 = actago1do3p)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago1do3p' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'actago1do3p' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'oldago1do3p' &
          , p3 = oldago1do3p)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago1do3p' &
          , 'long_name', c='Reaction conversion rate in the current timestep')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'oldago1do3p' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

!     CALL new_channel_object(status, modstr, 'OHv0nd' &
!          , p3 = OHv0nd)
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'OHv0nd' &
!          , 'long_name', c='vmr of OH in vibrational state v=0')
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'OHv0nd' &
!          , 'units', c='#/m-3')
!     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv1nd' &
          , p3 = OHv1nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv1nd' &
          , 'long_name', c='vmr of OH in vibrational state v=1')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv1nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv2nd' &
          , p3 = OHv2nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv2nd' &
          , 'long_name', c='vmr of OH in vibrational state v=2')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv2nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv3nd' &
          , p3 = OHv3nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv3nd' &
          , 'long_name', c='vmr of OH in vibrational state v=3')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv3nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv4nd' &
          , p3 = OHv4nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv4nd' &
          , 'long_name', c='vmr of OH in vibrational state v=4')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv4nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv5nd' &
          , p3 = OHv5nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv5nd' &
          , 'long_name', c='vmr of OH in vibrational state v=5')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv5nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv6nd' &
          , p3 = OHv6nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv6nd' &
          , 'long_name', c='vmr of OH in vibrational state v=6')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv6nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv7nd' &
          , p3 = OHv7nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv7nd' &
          , 'long_name', c='vmr of OH in vibrational state v=7')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv7nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv8nd' &
          , p3 = OHv8nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv8nd' &
          , 'long_name', c='vmr of OH in vibrational state v=8')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv8nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'OHv9nd' &
          , p3 = OHv9nd)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv9nd' &
          , 'long_name', c='vmr of OH in vibrational state v=9')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'OHv9nd' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'AG96' &
          , p3 = AG96)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG96' &
          , 'long_name', c='Airglow from OH v=9 to v=6; 1382nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG96' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'AG85' &
          , p3 = AG85)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG85' &
          , 'long_name', c='Airglow from OH v=8 to v=5; 1290nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG85' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'AG83' &
          , p3 = AG83)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG83' &
          , 'long_name', c='Airglow from OH v=8 to v=3; 728nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG83' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'AG62' &
          , p3 = AG62)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG62' &
          , 'long_name', c='Airglow from OH v=6 to v=2; 834nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG62' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'AG31' &
          , p3 = AG31)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG31' &
          , 'long_name', c='Airglow from OH v=3 to v=1; 1540nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AG31' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'AGO21d' &
          , p3 = AGO21d)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO21d' &
          , 'long_name', c='Airglow from O2(1D); 728nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO21d' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'AGO21s' &
          , p3 = AGO21s)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO21s' &
          , 'long_name', c='Airglow from O2(1S); 728nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO21s' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'AGO1SO3P' &
          , p3 = AGO1SO3P)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO1SO3P' &
          , 'long_name', c='Airglow from O(1S) to O(3P); 297nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO1SO3P' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'AGO1SO1D' &
          , p3 = AGO1SO1D)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO1SO1D' &
          , 'long_name', c='Airglow from O(1S) to O(1D); 557nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO1SO1D' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'AGO1DO3P' &
          , p3 = AGO1DO3P)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO1DO3P' &
          , 'long_name', c='Airglow from O(1D) to O(3P); 630nm')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'AGO1DO3P' &
          , 'units', c='photons/cm3/s')
     CALL channel_halt(substr, status)     

     CALL new_channel_object(status, modstr, 'debug' &
          , p3 = debug)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'debug' &
          , 'long_name', c='Just for debugging')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'debug' &
          , 'units', c='#/m-3')
     CALL channel_halt(substr, status)

  END SUBROUTINE mesoenergy_init_memory


  SUBROUTINE mesoenergy_init_coupling

    USE messy_main_mpi_bi,           ONLY: finish, message, p_parallel_io

    USE messy_main_blather_bi,        ONLY: start_message_bi, end_message_bi
    USE messy_main_tracer_mem_bi,     ONLY: GPTRSTR
    USE messy_main_channel_error_bi,  ONLY: channel_halt 
    USE messy_main_tracer_tools_bi,   ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: get_tracer
    USE messy_main_channel,       ONLY: get_channel_object, get_channel_info &
                                      , new_channel_object  & 
                                      , new_attribute 

    IMPLICIT NONE

    INTEGER :: status, i_err
    CHARACTER(LEN=*), PARAMETER :: substr = 'mesoenergy_init_coupling'

    CALL start_message_bi(modstr,'INIT COUPLING',substr)

    CALL get_channel_object(status, 'scnbuf', 'tte', p3=tte_ptr)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'tracer_gp', 'LossG2101', p3=lossg2101_ptr)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'tracer_gp', 'Airglow96', p3=ag96_ptr)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'tracer_gp', 'Airglow85', p3=ag85_ptr)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'tracer_gp', 'Airglow83', p3=ag83_ptr)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'tracer_gp', 'Airglow62', p3=ag62_ptr)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'tracer_gp', 'AirglowO21d', p3=ago21d_ptr)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'tracer_gp', 'AirglowO21sO2', p3=ago21s_ptr)
    CALL channel_halt(substr, status)

    CALL get_channel_object(status, 'tracer_gp', 'Airglow31', p3=ag31_ptr)
    CALL channel_halt(substr, status)
    
    CALL get_channel_object(status, 'tracer_gp', 'AirglowO1SO3P', p3=ago1so3p_ptr)
    CALL channel_halt(substr, status)
    
    CALL get_channel_object(status, 'tracer_gp', 'AirglowO1SO1D', p3=ago1so1d_ptr)
    CALL channel_halt(substr, status)
    
    CALL get_channel_object(status, 'tracer_gp', 'AirglowO1DO3P', p3=ago1do3p_ptr)
    CALL channel_halt(substr, status)

      CALL get_tracer(i_err, GPTRSTR, 'LossG2101', idx=idt_lossg2101)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'kJmol', idx=idt_kJmol)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'Airglow96', idx=idt_ag96)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'Airglow85', idx=idt_ag85)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'Airglow83', idx=idt_ag83)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'Airglow62', idx=idt_ag62)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'Airglow31', idx=idt_ag31)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'AirglowO21d', idx=idt_ago21d)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'AirglowO21sO2', idx=idt_ago21s)
      CALL tracer_halt(substr, i_err)
      
      CALL get_tracer(i_err, GPTRSTR, 'AirglowO1SO3P', idx=idt_ago1so3p)
      CALL tracer_halt(substr, i_err)
      
      CALL get_tracer(i_err, GPTRSTR, 'AirglowO1SO1D', idx=idt_ago1so1d)
      CALL tracer_halt(substr, i_err)
      
      CALL get_tracer(i_err, GPTRSTR, 'AirglowO1DO3P', idx=idt_ago1do3p)
      CALL tracer_halt(substr, i_err)

!      CALL get_tracer(i_err, GPTRSTR, 'OHv0', idx=idt_OHv0)
!      CALL tracer_halt(substr, i_err) 

      CALL get_tracer(i_err, GPTRSTR, 'OHv1', idx=idt_OHv1)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'OHv2', idx=idt_OHv2)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'OHv3', idx=idt_OHv3)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'OHv4', idx=idt_OHv4)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'OHv5', idx=idt_OHv5)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'OHv6', idx=idt_OHv6)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'OHv7', idx=idt_OHv7)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'OHv8', idx=idt_OHv8)
      CALL tracer_halt(substr, i_err)

      CALL get_tracer(i_err, GPTRSTR, 'OHv9', idx=idt_OHv9)
      CALL tracer_halt(substr, i_err)

    CALL end_message_bi(modstr,'INIT COUPLING',substr)

  END SUBROUTINE mesoenergy_init_coupling
 
! !=============================================================================
 
!!SUBROUTINE mesoenergy_init_tracer
 
!!   USE messy_main_tracer_bi, ONLY: main_tracer_initialize
 
!!   IMPLICIT NONE
 
!!   IF (tracer_init_required) CALL tracer_init(modstr)
 
!!END SUBROUTINE mesoenergy_init_tracer

! !=============================================================================
 
 SUBROUTINE mesoenergy_physc
  !--------------------------------------------------------------------------
  USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1
  USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma, nproma
  USE messy_main_grid_def_bi,     ONLY: grmass, grvol &
                                      , philat_2d
  USE messy_main_data_bi,       ONLY: pmid => press_3d, temp => t_scb & ! mid-level pressures [Pa]
                                    , pint => pressi_3d  ! level interface (above) pressure [Pa]
  USE messy_main_timer,         ONLY: time_step_len  
  ! op_pj_20180713+
  USE messy_main_tools,         ONLY: vmr_to_nd
  ! op_pj_20180713-
  USE messy_main_constants_mem, ONLY: M_air, R_gas, N_A, g, k_b
! op_pj_20180713: illegal use of other submodel; this needs to be replaced
!                 by channel object (SMIL) and formal parameter (to SMCL))
  USE messy_edith_msis,          ONLY: cpvf


  IMPLICIT none

  INTEGER :: jp, jk, n,  i, status

  REAL(dp), DIMENSION(nproma,nlev) :: HR2101_t0
  REAL(dp), DIMENSION(nproma,nlev) :: HR2101_g, lossg2101_g, ag96_g, ag85_g, ag83_g, ag62_g, ag31_g, ago21d_g, ago21s_g
  REAL(dp), DIMENSION(nproma,nlev) :: ago1so3p_t0, ago1so3p_g, ago1so1d_t0, ago1so1d_g, ago1do3p_t0, ago1do3p_g
  REAL(dp), DIMENSION(nproma,nlev) :: kJmol_t0, CHEMHEATING_g
  REAL(dp), DIMENSION(nproma,nlev) :: OHv0_t0, OHv1_t0, OHv2_t0, OHv3_t0, OHv4_t0, OHv5_t0, OHv6_t0, OHv7_t0, OHv8_t0, OHv9_t0
  REAL(dp), DIMENSION(nproma,nlev) :: ag96_t0, ag85_t0, ag83_t0, ag62_t0, ag31_t0, ago21d_t0, ago21s_t0
  REAL(dp), DIMENSION(nproma,nlev) :: lossg2101_t0
  REAL(dp), DIMENSION(nproma,nlev) :: cair
  REAL(dp), DIMENSION(nproma,nlev) :: cp


do jk=1,nlev
    DO jp=1, nproma
      cp(jp,jk)=cpvf(jk)
    ENDDO
  ENDDO


level_loop_tracer: DO jk=1,nlev
    vector_loop_tracer: DO jp=1,kproma
        lossg2101_t0(jp,jk) = pxtm1(jp,jk,idt_lossg2101)
        OHv1_t0(jp,jk) = pxtm1(jp,jk,idt_OHv1)
        OHv2_t0(jp,jk) = pxtm1(jp,jk,idt_OHv2)
        OHv3_t0(jp,jk) = pxtm1(jp,jk,idt_OHv3)
        OHv4_t0(jp,jk) = pxtm1(jp,jk,idt_OHv4)
        OHv5_t0(jp,jk) = pxtm1(jp,jk,idt_OHv5)
        OHv6_t0(jp,jk) = pxtm1(jp,jk,idt_OHv6)
        OHv7_t0(jp,jk) = pxtm1(jp,jk,idt_OHv7)
        OHv8_t0(jp,jk) = pxtm1(jp,jk,idt_OHv8)
        OHv9_t0(jp,jk) = pxtm1(jp,jk,idt_OHv9)
        ag96_t0(jp,jk) = pxtm1(jp,jk,idt_ag96)
        ag85_t0(jp,jk) = pxtm1(jp,jk,idt_ag85)
        ag83_t0(jp,jk) = pxtm1(jp,jk,idt_ag83)
        ag62_t0(jp,jk) = pxtm1(jp,jk,idt_ag62)
        ag31_t0(jp,jk) = pxtm1(jp,jk,idt_ag31)
        ago21d_t0(jp,jk) = pxtm1(jp,jk,idt_ago21d)
        ago21s_t0(jp,jk) = pxtm1(jp,jk,idt_ago21s)
        ago1so3p_t0(jp,jk) = pxtm1(jp,jk,idt_ago1so3p)
        ago1so1d_t0(jp,jk) = pxtm1(jp,jk,idt_ago1so1d)
        ago1do3p_t0(jp,jk) = pxtm1(jp,jk,idt_ago1do3p)
        
        kJmol_t0(jp,jk) = pxtm1(jp,jk,idt_kJmol)*1000._dp          ! conversion to J/mol
   END DO vector_loop_tracer
  END DO level_loop_tracer

  level_loop1c: DO jk=1,nlev
    vector_loop1c: DO jp=1,kproma
        call vmr_to_nd(lossg2101_g(jp,jk),lossg2101_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ag96_g(jp,jk),ag96_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ag85_g(jp,jk),ag85_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ag83_g(jp,jk),ag83_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ag62_g(jp,jk),ag62_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ag31_g(jp,jk),ag31_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ago21d_g(jp,jk),ago21d_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ago21s_g(jp,jk),ago21s_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ago1so3p_g(jp,jk),ago1so3p_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ago1so1d_g(jp,jk),ago1so1d_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(ago1do3p_g(jp,jk),ago1do3p_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv1nd(jp,jk,jrow),OHv1_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv2nd(jp,jk,jrow),OHv2_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv3nd(jp,jk,jrow),OHv3_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv4nd(jp,jk,jrow),OHv4_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv5nd(jp,jk,jrow),OHv5_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv6nd(jp,jk,jrow),OHv6_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv7nd(jp,jk,jrow),OHv7_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv8nd(jp,jk,jrow),OHv8_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(OHv9nd(jp,jk,jrow),OHv9_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(cair(jp,jk),1._dp,pmid(jp,jk,jrow),temp(jp,jk,jrow))
        call vmr_to_nd(CHEMHEATING_g(jp,jk),kJmol_t0(jp,jk),pmid(jp,jk,jrow),temp(jp,jk,jrow))
   END DO vector_loop1c
  END DO level_loop1c

  level_loop1b: DO jk=1,nlev
    vector_loop1b: DO jp=1,kproma
     actloss(jp,jk,jrow)=max((lossg2101_g(jp,jk))*grvol(jp,jk,jrow),1.e-30)
     CHEMHEATING(jp,jk,jrow)=max((CHEMHEATING_g(jp,jk))*grvol(jp,jk,jrow),1.e-30)
     actag96(jp,jk,jrow)=max((ag96_g(jp,jk)),1.e-30)
     actag85(jp,jk,jrow)=max((ag85_g(jp,jk)),1.e-30)
     actag83(jp,jk,jrow)=max((ag83_g(jp,jk)),1.e-30)
     actag62(jp,jk,jrow)=max((ag62_g(jp,jk)),1.e-30)
     actag31(jp,jk,jrow)=max((ag31_g(jp,jk)),1.e-30)
     actago21d(jp,jk,jrow)=max((ago21d_g(jp,jk)),1.e-30)
     actago21s(jp,jk,jrow)=max((ago21s_g(jp,jk)),1.e-30)
     actago1so3p(jp,jk,jrow)=max((ago1so3p_g(jp,jk)),1.e-30)
     actago1so1d(jp,jk,jrow)=max((ago1so1d_g(jp,jk)),1.e-30)
     actago1do3p(jp,jk,jrow)=max((ago1do3p_g(jp,jk)),1.e-30)
     
   END DO vector_loop1b
  END DO level_loop1b


   IF (chemheat) THEN
level_loop0: DO jk=1,nlev
    vector_loop0: DO jp=1,kproma

        HR2101_g(jp,jk) = 196825.475_dp*lossg2101_g(jp,jk)/N_A/cp(jp,jk)/grmass(jp,jk,jrow)*grvol(jp,jk,jrow)

    END DO vector_loop0
  END DO level_loop0


   level_loop1d: DO jk=1,nlev
     vector_loop1d: DO jp=1,kproma
           tte_ptr(jp,jk,jrow) = tte_ptr(jp,jk,jrow) &
         + HR2101_g(jp,jk)/time_step_len  ! -> tendency [K/s]
!           tte_ptr(jp,jk,jrow) = tte_ptr(jp,jk,jrow) &
!         + CHEMHEATING_g(jp,jk)/time_step_len
     END DO vector_loop1d
   END DO level_loop1d
  ENDIF



      level_loop1e: DO jk=1,nlev
    vector_loop1e: DO jp=1,kproma
          HR2101(jp,jk,jrow)=HR2101_g(jp,jk)/time_step_len*86400.
          CHEMHEATING(jp,jk,jrow)=CHEMHEATING_g(jp,jk)/time_step_len*86400.
           AG96(jp,jk,jrow) = actag96(jp,jk,jrow)/1.e6/time_step_len
           AG85(jp,jk,jrow) = actag85(jp,jk,jrow)/1.e6/time_step_len
           AG83(jp,jk,jrow) = actag83(jp,jk,jrow)/1.e6/time_step_len
           AG62(jp,jk,jrow) = actag62(jp,jk,jrow)/1.e6/time_step_len
           AG31(jp,jk,jrow) = actag31(jp,jk,jrow)/1.e6/time_step_len
           AGO21d(jp,jk,jrow) = actago21d(jp,jk,jrow)/1.e6/time_step_len
           AGO21s(jp,jk,jrow) = actago21s(jp,jk,jrow)/1.e6/time_step_len
           AGO1SO3P(jp,jk,jrow) = actago1so3p(jp,jk,jrow)/1.e6/time_step_len
           AGO1SO1D(jp,jk,jrow) = actago1so1d(jp,jk,jrow)/1.e6/time_step_len
           AGO1DO3P(jp,jk,jrow) = actago1do3p(jp,jk,jrow)/1.e6/time_step_len
        pxtm1(jp,jk,idt_lossg2101) = 0.d0
        pxtm1(jp,jk,idt_ag96) = 0.d0
        pxtm1(jp,jk,idt_ag85) = 0.d0
        pxtm1(jp,jk,idt_ag83) = 0.d0
        pxtm1(jp,jk,idt_ag62) = 0.d0
        pxtm1(jp,jk,idt_ag31) = 0.d0
        pxtm1(jp,jk,idt_ago21d) = 0.d0
        pxtm1(jp,jk,idt_ago21s) = 0.d0
        pxtm1(jp,jk,idt_ago1so3p) = 0.d0
        pxtm1(jp,jk,idt_ago1so1d) = 0.d0
        pxtm1(jp,jk,idt_ago1do3p) = 0.d0
        pxtm1(jp,jk,idt_kJmol) = 0.d0
    END DO vector_loop1e
   END DO level_loop1e

!! ------------------------------------------
!! Calculations of OH(v)
!! ------------------------------------------

   

END SUBROUTINE mesoenergy_physc

!=============================================================================

!SUBROUTINE mesoenergy_local_end

!END SUBROUTINE mesoenergy_local_end

!=============================================================================

SUBROUTINE mesoenergy_read_nml_cpl(status, iou)
  !------------------------------------------
  ! read coupling namelist /CPL/ from mesoenergy.nml
  !------------------------------------------

  ! MESSy
  USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
  
  IMPLICIT NONE
  
  INTEGER, INTENT(out) :: status
  INTEGER, INTENT(in) :: iou
  
  CHARACTER(len=*), PARAMETER :: substr='mesoenergy_read_nml_cpl'
  LOGICAL :: lex     ! file exists?
  INTEGER :: fstat   ! file status
  
  NAMELIST /CPL/ dummy 
    
  status = 1             ! initialise status flag with error code
  
  CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
  IF (.not.lex) RETURN   ! error: psc.nml does not exist
  read(iou, nml=cpl, iostat=fstat)
  CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
  IF (fstat/=0) RETURN   ! error while reading namelist
  

  CALL read_nml_close(substr, iou, modstr)
  status = 0   ! no error
  
END SUBROUTINE mesoenergy_read_nml_cpl

!=============================================================================

END MODULE messy_mesoenergy_e5
