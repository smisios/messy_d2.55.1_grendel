! **********************************************************************
!
!ka_av_20151229+
! original by Alexey Vlasov (2014-2016)
! adapted to MESSy-Standard by Stefan Versick (2018)
!
! **********************************************************************

! **********************************************************************
MODULE messy_edith_e5
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi, warning_bi

  ! SMCL
  USE messy_edith
  USE messy_main_mpi_bi,       ONLY: message
  USE messy_main_channel,       ONLY: t_chaobj_cpl
#ifdef MESSYTENDENCY  
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,  mtend_register,  &
                                      mtend_get_start_l, mtend_id_t,      &
                                      mtend_id_q,        mtend_id_xl,     &
                                      mtend_id_xi,       mtend_id_u,      &
                                      mtend_id_v,        mtend_id_tracer, &
                                      mtend_add_l
#endif

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE ! op_pj_20180713

  REAL(dp), DIMENSION(:,:,:), POINTER :: &
     amu, cp, g, alt
     
  REAL(dp), DIMENSION(:,:,:), POINTER :: &
     tteion, ttevdiffmol, qtevdiffmol, utevdiffmol, vtevdiffmol, &
     fricvdiffmol, uteiondrag, vteiondrag
     
  REAL(dp), DIMENSION(:,:,:), POINTER :: tte_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: um1_ptr => NULL()
  
  REAL(dp), DIMENSION(:,:), POINTER :: philat_ptr => NULL()
     
  INTEGER :: idt_O2, idt_N2, idt_O3P, idt_H, idt_myNOX, idt_CO2, idt_O1D, idt_O1S, idt_CO
  
  ! ka_sv_20150521+
  REAL(dp), DIMENSION(:,:,:), POINTER :: &
    NOcooling  => NULL(),  &
    CO2NIR     => NULL(),  &
    CO2NLTE    => NULL(),  &
    CO2nd => NULL(),   &
    N2nd => NULL(),   &
    O2nd  => NULL(),   &
    O3nd  => NULL(),   &
    Ond   => NULL(),   &
! ka_sv_20180702+
    JO2    => NULL(),  &
    JCO2    => NULL(),  &
    JO2t    => NULL(),  &
    JCO2t    => NULL(), &
    moldiffcoef => NULL()
! ka_sv_20180702-
! ka_sv_20170420+
!    debug1 => NULL()
!    debug2 => NULL(), &
!    debug3 => NULL(), &
!    debug4 => NULL(), &
!    debug5 => NULL()
! ka_sv_20170420+
    
  REAL(dp), DIMENSION(:,:,:), POINTER :: no_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: o_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: co2_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: o3_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: o2_ptr => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: n2_ptr => NULL()
  
  REAL(dp), DIMENSION(:,:), POINTER :: cossza_2d => NULL()

  TYPE(t_chaobj_cpl) :: edith_no, edith_co2, edith_o, edith_o3, edith_o2, edith_n2
  TYPE(t_chaobj_cpl) :: edith_kp                   !kit_sb_20180419
    
#ifdef MESSYTENDENCY
  INTEGER                             :: my_handle
#endif

  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_e5.f90
  PUBLIC :: edith_initialize    ! initialize submodel
  PUBLIC :: edith_init_memory
  PUBLIC :: edith_init_coupling ! set pointers for coupling to BM and other SMs
  PUBLIC :: edith_physc         ! entry point in time loop (current vector)


CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE edith_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'edith_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

   
    ! READ CTRL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find free I/O unit
       CALL edith_read_nml_ctrl(status, iou)  ! read CTRL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
    END IF
    ! BROADCAST CTRL namleist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(iamu, p_io)
    CALL p_bcast(calc_iondrag, p_io)
    CALL p_bcast(passive_nox, p_io)
    CALL p_bcast(use_chem, p_io)
    CALL p_bcast(prandtlnumber, p_io)
    CALL p_bcast(ed_solvar, p_io)
    !ka_sv_20180702+
    CALL p_bcast(edith_phot, p_io)
    CALL p_bcast(ed_o2_o3p, p_io)
    CALL p_bcast(ed_o2_o1d, p_io)
    CALL p_bcast(ed_o2_o1s, p_io)
    !ka_sv_20180702-

    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL edith_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF
    ! BROADCAST CPL namleist entries from I/O-PE to ALL OTHER PEs
    
    CALL p_bcast(edith_co2%cha, p_io)
    CALL p_bcast(edith_co2%obj, p_io)
    CALL p_bcast(edith_no%cha, p_io)
    CALL p_bcast(edith_no%obj, p_io)
    CALL p_bcast(edith_o2%cha, p_io)
    CALL p_bcast(edith_o2%obj, p_io)
    CALL p_bcast(edith_o3%cha, p_io)
    CALL p_bcast(edith_o3%obj, p_io)
    CALL p_bcast(edith_o%cha, p_io)
    CALL p_bcast(edith_o%obj, p_io)
    CALL p_bcast(edith_n2%cha, p_io)
    CALL p_bcast(edith_n2%obj, p_io)
!ka_sb_20180419+
    CALL p_bcast(edith_kp%cha, p_io)
    CALL p_bcast(edith_kp%obj, p_io) 
!ka_sb_20180419-

    ! ### PERFORM INITIAL SETUP (CALL RESPECTIVE SMCL ROUTINE(S)) HERE
    
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE edith_initialize
  ! ====================================================================

! ====================================================================
SUBROUTINE edith_new_tracer

     ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, set_tracer &
                                      , R_molarmass ! ,ON, OFF
    USE messy_main_constants_mem, ONLY: MO, MN

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'edith_new_tracer'
    INTEGER                     :: status
!    INTEGER                     :: idt_myNOX  ! tracer index (identifier)

    CALL start_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output
    IF (passive_nox.ne.0) THEN
    ! ### define new tracers here
      CALL new_tracer(status, GPTRSTR, 'myNOX', modstr, idx=idt_myNOX)
      CALL tracer_halt(substr, status)   ! terminate if error
      CALL set_tracer(status, GPTRSTR, idt_myNOX &
           , R_molarmass, r=1.0_dp*MO+1.0_dp*MN)
      CALL tracer_halt(substr, status)   ! terminate if error
    END IF

    CALL end_message_bi(modstr,'TRACER DEFINITION',substr)  ! log-output

  END SUBROUTINE edith_new_tracer
! ====================================================================


   SUBROUTINE edith_init_memory
! ---------------------------------------------------------------
 
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

     IMPLICIT NONE
     INTEGER :: status
     CHARACTER(LEN=*), PARAMETER           :: substr = 'edith_init_memory'
 
#ifdef MESSYTENDENCY
    CALL mtend_register (my_handle,mtend_id_t)
    CALL mtend_register (my_handle,mtend_id_q)
    CALL mtend_register (my_handle,mtend_id_xl)
    CALL mtend_register (my_handle,mtend_id_xi)
    CALL mtend_register (my_handle,mtend_id_u)
    CALL mtend_register (my_handle,mtend_id_v)
    CALL mtend_register (my_handle,mtend_id_tracer)
#endif
 
     CALL message('edith_init_memory','defining streams for EDITh')
 
     ! -----
     ! ubcnox specific molecule information
     ! -----
     CALL new_channel(status, modstr, reprid=GP_3D_MID)
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'amu' &
          , p3 = amu)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'amu' &
          , 'long_name', c='Molecular weight of air')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'amu' &
          , 'units', c='g/mol')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'cp' &
          , p3 = cp)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'cp' &
          , 'long_name', c='specific heat at constant pressure')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'cp' &
          , 'units', c='J/kg/K')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'g' &
          , p3 = g)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'g' &
          , 'long_name', c='gravitational acceleration')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'g' &
          , 'units', c='m/s2')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'alt' &
          , p3 = alt)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'alt' &
          , 'long_name', c='altitude based on geopotential')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'alt' &
          , 'units', c='m')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'tteion' &
          , p3 = tteion)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'tteion' &
          , 'long_name', c='temperature tendency due to iondrag')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'tteion' &
          , 'units', c='K/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'ttevdiffmol' &
          , p3 = ttevdiffmol)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'ttevdiffmol' &
          , 'long_name', c='temperature tendency due to molecular diffusion')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'ttevdiffmol' &
          , 'units', c='K/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'qtevdiffmol' &
          , p3 = qtevdiffmol)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'qtevdiffmol' &
          , 'long_name', c='moisture tendency due to molecular diffusion')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'qtevdiffmol' &
          , 'units', c='kg/kg/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'utevdiffmol' &
          , p3 = utevdiffmol)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'utevdiffmol' &
          , 'long_name', c='zonal wind tendency due to molecular diffusion')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'utevdiffmol' &
          , 'units', c='m/s/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'vtevdiffmol' &
          , p3 = vtevdiffmol)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'vtevdiffmol' &
          , 'long_name', c='meridional wind tendency due to molecular diffusion')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'vtevdiffmol' &
          , 'units', c='m/s/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'fricvdiffmol' &
          , p3 = fricvdiffmol)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'fricvdiffmol' &
          , 'long_name', c='temperature tendency due to frictional heating')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'fricvdiffmol' &
          , 'units', c='K/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'uteiondrag' &
          , p3 = uteiondrag)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'uteiondrag' &
          , 'long_name', c='zonal wind tendency due to iondrag')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'uteiondrag' &
          , 'units', c='m/s/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'vteiondrag' &
          , p3 = vteiondrag)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'vteiondrag' &
          , 'long_name', c='meridional wind tendency due to iondrag')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'vteiondrag' &
          , 'units', c='m/s/s')
     CALL channel_halt(substr, status)
     
          CALL new_channel_object(status, modstr, 'NOcooling' &
          , p3 = NOcooling ,reprid=GP_3D_MID)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'NOcooling' &
          , 'long_name', c='Cooling rate by NO')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'NOcooling' &
          , 'units', c='K/d')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'CO2NIR' &
          , p3 = CO2NIR ,reprid=GP_3D_MID)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'CO2NIR' &
          , 'long_name', c='Cooling rate by CO2 NIR')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'CO2NIR' &
          , 'units', c='K/d')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'CO2NLTE' &
          , p3 = CO2NLTE ,reprid=GP_3D_MID)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'CO2NLTE' &
          , 'long_name', c='Cooling rate by CO2 NLTE')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'CO2NLTE' &
          , 'units', c='K/d')
     CALL channel_halt(substr, status)

!ka_sv_20180702+     
     CALL new_channel_object(status, modstr, 'JO2' &
          , p3 = JO2 ,reprid=GP_3D_MID)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'JO2' &
          , 'long_name', c='SRC photolysis rate of O2')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'JO2' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'JCO2' &
          , p3 = JCO2 ,reprid=GP_3D_MID)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'JCO2' &
          , 'long_name', c='SRC photolysis rate of CO2')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'JCO2' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'JO2t' &
          , p3 = JO2t ,reprid=GP_3D_MID)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'JO2t' &
          , 'long_name', c='SRC photolysis rate of O2: tendency')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'JO2t' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)

     CALL new_channel_object(status, modstr, 'JCO2t' &
          , p3 = JCO2t ,reprid=GP_3D_MID)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'JCO2t' &
          , 'long_name', c='SRC photolysis rate of CO2: tendency')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'JCO2t' &
          , 'units', c='1/s')
     CALL channel_halt(substr, status)
     
     CALL new_channel_object(status, modstr, 'moldiffcoef' &
          , p3 = moldiffcoef ,reprid=GP_3D_MID)
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'moldiffcoef' &
          , 'long_name', c='Molecular diffusion coefficient')
     CALL channel_halt(substr, status)
     CALL new_attribute(status, modstr, 'moldiffcoef' &
          , 'units', c='m2s-1')
     CALL channel_halt(substr, status)
!ka_sv_20180702-
          
!     CALL new_channel_object(status, modstr, 'debug1' &
!          , p3 = debug1 ,reprid=GP_3D_MID)
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug1' &
!          , 'long_name', c='debug1')
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug1' &
!               , 'units', c='1')
!     CALL channel_halt(substr, status)

!     CALL new_channel_object(status, modstr, 'debug2' &
!          , p3 = debug2 ,reprid=GP_3D_MID)
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug2' &
!          , 'long_name', c='debug2')
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug2' &
!               , 'units', c='1')
!     CALL channel_halt(substr, status)
     
!          CALL new_channel_object(status, modstr, 'debug3' &
!          , p3 = debug3 ,reprid=GP_3D_MID)
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug3' &
!          , 'long_name', c='debug3')
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug3' &
!               , 'units', c='1')
!     CALL channel_halt(substr, status)
     
!          CALL new_channel_object(status, modstr, 'debug4' &
!          , p3 = debug4 ,reprid=GP_3D_MID)
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug4' &
!          , 'long_name', c='debug4')
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug4' &
!               , 'units', c='1')
!     CALL channel_halt(substr, status)
     
!               CALL new_channel_object(status, modstr, 'debug5' &
!          , p3 = debug5 ,reprid=GP_3D_MID)
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug5' &
!          , 'long_name', c='debug5')
!     CALL channel_halt(substr, status)
!     CALL new_attribute(status, modstr, 'debug5' &
!               , 'units', c='1')
!     CALL channel_halt(substr, status)


  END SUBROUTINE edith_init_memory


  ! ====================================================================
  SUBROUTINE edith_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the 
    ! basemodel and to other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object
    !
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_tracer,           ONLY: get_tracer
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_mpi_bi,           ONLY: finish
    
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'edith_init_coupling'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output
    
    CALL get_channel_object(status, 'scnbuf', 'tte', p3=tte_ptr)
    CALL channel_halt(substr, status)
    
    !CALL get_channel_object(status, 'grid_def', 'philat_2d', p2=philat_ptr)
    !CALL channel_halt(substr, status)
    
    IF (use_chem) THEN
      CALL get_tracer(status, GPTRSTR, TRIM('O2'), idx=idt_O2)
      CALL tracer_halt(substr, status)
      
            CALL get_tracer(status, GPTRSTR, TRIM('N2'), idx=idt_N2)
      CALL tracer_halt(substr, status)
      
            CALL get_tracer(status, GPTRSTR, TRIM('O3P'), idx=idt_O3P)
      CALL tracer_halt(substr, status)
      
            CALL get_tracer(status, GPTRSTR, TRIM('H'), idx=idt_H)
      CALL tracer_halt(substr, status)
      
      CALL get_tracer(status, GPTRSTR, TRIM('CO2'), idx=idt_CO2)
      CALL tracer_halt(substr, status)
    END IF
    
    IF (edith_phot) THEN
      CALL get_tracer(status, GPTRSTR, TRIM('O1D'), idx=idt_O1D)
      CALL tracer_halt(substr, status)
      
!      CALL get_tracer(status, GPTRSTR, TRIM('O1S'), idx=idt_O1S)
!      CALL tracer_halt(substr, status)
      
      CALL get_tracer(status, GPTRSTR, TRIM('CO'), idx=idt_CO)
      CALL tracer_halt(substr, status)

    END IF

    IF (passive_nox.ne.0) THEN
      CALL get_tracer(status, GPTRSTR, TRIM('myNOX'), idx=idt_myNOX)
      CALL tracer_halt(substr, status)    ! terminate on error
    END IF

!ka_sb_20180507+
    IF (ed_solvar.eq.4) THEN  ! kp-Index needed for variable solar conditions
       CALL get_channel_object(status &
            , TRIM(edith_kp%cha), TRIM(edith_kp%obj), p1=kp_data)
       CALL channel_halt(substr, status) 
    END IF   
!ka_sb_20180507-

    
    CALL get_channel_object(status, TRIM(edith_no%cha), TRIM(edith_no%obj), p3=no_ptr)
    IF (status /= 0) &
       CALL finish(substr,'channel object not found NO')
       
    CALL get_channel_object(status, TRIM(edith_o%cha), TRIM(edith_o%obj), p3=o_ptr)
    IF (status /= 0) &
       CALL finish(substr,'channel object not found O')
       
    CALL get_channel_object(status, TRIM(edith_co2%cha), TRIM(edith_co2%obj), p3=co2_ptr)
    IF (status /= 0) &
       CALL finish(substr,'channel object not found CO2')

    CALL get_channel_object(status, TRIM(edith_o2%cha), TRIM(edith_o2%obj), p3=o2_ptr)
    IF (status /= 0) &
       CALL finish(substr,'channel object not found O2')

    CALL get_channel_object(status, TRIM(edith_o3%cha), TRIM(edith_o3%obj), p3=o3_ptr)
    IF (status /= 0) &
       CALL finish(substr,'channel object not found O3')

    CALL get_channel_object(status, TRIM(edith_n2%cha), TRIM(edith_n2%obj), p3=n2_ptr)
    IF (status /= 0) &
       CALL finish(substr,'channel object not found N2')
       
    CALL get_channel_object(status, 'orbit','cossza',p2=cossza_2d)
    CALL channel_halt(substr//': channel object for cossza in orbit not found!',status)

    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output

  END SUBROUTINE edith_init_coupling
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE edith_physc

    ! ------------------------------------------------------------------
    ! This subroutine is called within the time loop.
    ! It constitutes the main entry point for additional processes 
    ! or diagnostics.
    ! Here, only the current vector of the grid-point-fields is
    ! accessible.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_timer,         ONLY: lstart, time_step_len
    USE messy_main_grid_def_bi,     ONLY:  philat_2d
    USE messy_main_grid_def_mem_bi, ONLY: kproma,nproma, jrow, nlev
    USE messy_main_data_bi,       ONLY: geopot_3d &
!ka_av_20150104+
                                      , apm1,aphm1,xlm1, xim1 &
                                      ! molecular diffusion coefficient
! op_pj_20180713+: Note: molco is on messy_main_data_bi only declared as 
!                        pointer, but no memory is ever associated to it.
!                        Below values are assigend, this should cause a 
!                        segmentation fault. But it is never used further.
!                        Thus, it has been removed!
!!$                                      , molco &
! op_pj_20180713-
                                      ! ion drag
!                                      ionvom,ionvol &
                                      , geopot_3d &
                                      , vol_3d, vom_3d       &
                                      , xtte, qte_3d
! op_pj_20180713+: illegal USE of other submodel
!!$USE messy_e5vdiff,            ONLY: cvdifts
USE messy_main_constants_mem, ONLY: cvdifts
! op_pj_20180713-

!USE messy_cloud_ori,      ONLY: vtmpc1
!! variables needed for molecular viscousity (temporary here)
! op_pj_20180713+: illegal, obsolete USE of basemodel; replaced all *%ldc
!!$USE mo_decomposition,             ONLY: ldc => local_decomposition
! op_pj_20180713-
!USE mo_tracer,                    ONLY: ntrac
USE messy_main_tracer_mem_bi,            ONLY: ntrac_gp    
! op_pj_20180713+: illegal USE of basemodel; an even better solution would
!                  be to set local pointers to channel object memory
!                  (with get_channel_object)
!!$USE mo_memory_g1a,                ONLY: tm1, qm1, xtm1
!!$USE mo_memory_g2a,                ONLY: vm1, um1
!!$!USE mo_scan_buffer,               ONLY: vo, vol, vom, tte, qte!, xtte   ! get if from messy_main_data_bi
USE messy_main_data_bi, ONLY: tm1, qm1, vm1, um1 
! op_pj_20180713-
  USE messy_main_tracer_mem_bi, ONLY: pxtm1 => qxtm1
!  USE mo_memory_g3b,   ONLY: amatsave, bmatsave, alsave, co2colxsave, co2xsave
  USE messy_main_tracer_mem_bi, ONLY: ti_gp, I_vdiff, ON
  USE messy_main_constants_mem, ONLY: argas=>R_gas, RTD
  USE messy_main_tracer,        ONLY: R_molarmass
  
!!USE mo_gaussgrid,      ONLY: gl_twomu ! op_pj_20180713 not used
!ka_sv_20180702+
! op_pj_20180713+: illegal use of other submodel; moved to messy_main_tools
!!$USE messy_mesoenergy, ONLY: vmr_to_nd, nd_to_vmr
USE messy_main_tools, ONLY: vmr_to_nd, nd_to_vmr
! op_pj_20180713-
!ka_sv_20180702-
  
 IMPLICIT NONE
 REAL(dp)    :: zgmurhoh(nproma,nlev),zdiffmol(nproma,nlev)
! ion drag 
 REAL(dp) ::  zvom(nproma,nlev),zvol(nproma,nlev)
 REAL(dp) :: amu1(nproma,nlev),cp1(nproma,nlev),grav1(nproma,nlev)

  REAL(dp) :: ztvm1(kproma, nlev)
  INTEGER  :: kbdim
  INTEGER  :: jk, jl, jt
!    IMPLICIT NONE
!ka_av_20150104-
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'edith_physc'
    INTEGER                     :: status
    
    REAL(dp) :: twoatomicvmr,oneatomicvmr
    REAL(dp),  DIMENSION(:,:), POINTER    :: altitude  ! altitude [km]
    
    REAL(dp), DIMENSION(nproma,nlev) :: NO_t0, O_t0, CO2_t0, O2_t0, O3_t0, N2_t0

! variables from vdiff_mol
REAL(dp)    :: zalpha,zmalpha,ztmst,zbet(nproma,nlev),zgam(nproma,nlev),zpr,zinvpr
REAL(dp)    :: zmah(nproma,nlev),zrhok(nproma,nlev,ntrac_gp)
REAL(dp)    :: zhh(nproma,nlev),zrdp(nproma,nlev),zr(nproma,nlev),zrt(nproma,nlev,ntrac_gp)
REAL(dp)    :: za(nproma,nlev),zb(nproma,nlev),zc(nproma,nlev)
REAL(dp)    :: zla(nproma,nlev),zlb(nproma,nlev),zrstar
REAL(dp)    :: zlc(nproma,nlev), zut(nproma,nlev), zutt(nproma,nlev,ntrac_gp)
REAL(dp)    :: zrdz, zdudz, zdvdz, zcoef
REAL(dp)    :: zmm(ntrac_gp),zexpan(ntrac_gp)
REAL(dp)    :: zsumx(nproma,nlev),zuttsum(nproma,nlev)
REAL(dp)    :: zrhowd(nproma,nlev+1,ntrac_gp)
REAL(dp)    :: zgvh(nproma,nlev),zrhoh(nproma,nlev),zth(nproma,nlev),ztvh(nproma,nlev)
REAL(dp)    :: zpxtm1(nproma,nlev,ntrac_gp)

INTEGER, DIMENSION(ntrac_gp) :: trac_id_gp
REAL(dp) :: ptte(nproma,nlev), pvol(nproma,nlev), pvom(nproma,nlev)
REAL(dp) :: pqte(nproma,nlev), pxtte(nproma,nlev,ntrac_gp)
INTEGER :: jkk

!ka_sv_20180702+
REAL(dp) :: scolo2(nproma,nlev)
LOGICAL :: sill(nproma,nlev)
!REAL(dp) :: j_src(nproma,nlev,2)
REAL(dp) :: j_src(kproma,nlev,2)
REAL(dp) :: n_o2, n_co2     ! number densities for photolysis
! vertically flipped variables
REAL(dp) :: o2flip(nproma,nlev), tmflip(nproma,nlev), pressflip(nlev), j_src_flip(kproma,nlev,2)
!ka_sv_20180702-
!ka_sv_20180705+
REAL(dp) :: amatsave(nproma,43,9)
REAL(dp) :: bmatsave(nproma,43,9)
REAL(dp) :: alsave(nproma,17)
REAL(dp) :: co2colxsave(nproma,59)
REAL(dp) :: co2xsave(nproma,59)
!ka_sv_20180705-

!------------------------------------
ALLOCATE(altitude(kproma,nlev))

       ztvm1 = 0._dp
       ztvm1(:,:) = tm1(:,:,jrow) 
!*                        &
!                          (1.0_dp+vtmpc1*qm1(:,:,jrow) -          &
!                          (xlm1(:,:,jrow)+xim1(:,:,jrow)))
       kbdim = nproma

!qm1(:,1:20)=0._dp

altitude(1:kproma,:)=geopot_3d(1:kproma,:,jrow)/9.806_dp
     
   IF (use_chem) THEN
     DO jk=1,nlev
       DO jl=1,nproma
         twoatomicvmr=pxtm1(jl,jk,idt_O2)+pxtm1(jl,jk,idt_N2)
         oneatomicvmr=pxtm1(jl,jk,idt_O3P)+pxtm1(jl,jk,idt_H)
         CALL edith_calc_local_amu(iamu,apm1(jl,jk),amu1(jl,jk))
         CALL edith_calc_local_cp(cp1(jl,jk),amu1(jl,jk),twoatomicvmr,oneatomicvmr)
         CALL edith_calc_local_g(grav1(jl,jk),philat_2d(jl,jrow),altitude(jl,jk))
       END DO
     END DO
   ELSE
     DO jk=1,nlev
       DO jl=1,nproma
         CALL edith_calc_local_amu(iamu,apm1(jl,jk),amu1(jl,jk))
         cp1(jl,jk)=1004._dp
         CALL edith_calc_local_g(grav1(jl,jk),philat_2d(jl,jrow),altitude(jl,jk))
       END DO
     END DO
   END IF

! ka_sv_20180122+
! split old vdiff_mol to 5 seperate routines

! general settings
!-----------------------------
  zrstar=argas*1.e3_dp        ! universal gas constant 8314 J/kmol/K
  zpr=prandtlnumber         ! Standard number: 0.72
  zinvpr=1._dp/zpr
  zalpha=cvdifts
  zmalpha=1._dp-zalpha
  ztmst=time_step_len

  ! preliminary define thermal diffusion factor zexpan
  ! zexpan = -0.38 for H and He, and zexpan = 0 for all other species (to be adapted)
  zexpan=0._dp

  DO jt=1,ntrac_gp
    trac_id_gp(jt) = jt
    zmm(jt)=max(1._dp,ti_gp(trac_id_gp(jt))%tp%meta%cask_r(R_molarmass))
  ENDDO

  IF (ntrac_gp .GT. 0) THEN

  ! Integrate preceding processes
    DO jt = 1, ntrac_gp
      DO jk = 1, nlev
        DO jl = 1, kproma
           zpxtm1(jl,jk,jt) = pxtm1(jl,jk,jt) + ztmst*xtte(jl,jk,jt)
        END DO
      END DO
    END DO
  ! Make sure that mixing ratios are non-negative
    DO jt = 1, ntrac_gp
      DO jk = 1, nlev
        DO jl = 1, kproma
           zpxtm1(jl,jk,jt) = MAX(zpxtm1(jl,jk,jt),1.e-30_dp)
        END DO
      END DO
    END DO

!** IF INTERACTIVE CHEMISTRY (INCLUDING ALL MAJOR SPECIES)

    zsumx=0.0_dp
    DO jt=1,ntrac_gp
      IF (ti_gp(jt)%tp%meta%cask_i(I_vdiff) == ON) THEN
        zsumx=zsumx+zpxtm1(:,:,jt)
      ENDIF
    ENDDO

  ENDIF   ! ktrac .GT. 0

!** SET LAYERS' THICKNESS, (DYNAMIC VISCOSITY)*(pg)/(RT) 
!   AND (DENSITY * TRACER DIFFUSION COEFF)

  DO jk=2,nlev
    jkk=nlev-jk+1
    zhh     (1:kproma,jk)=1._dp/(apm1(1:kproma,jk-1)-apm1(1:kproma,jk))
    zrdp    (1:kproma,jk)=ztmst*grav1(1:kproma,jkk)/(aphm1(1:kproma,jk+1)-aphm1(1:kproma,jk))
    zgvh    (1:kproma,jk)=0.5_dp * ( grav1  (1:kproma,jkk) + grav1  (1:kproma,jkk+1) )
    zmah    (1:kproma,jk)=0.5_dp * ( amu1 (1:kproma,jk) + amu1 (1:kproma,jk-1) )
    zth     (1:kproma,jk)=0.5_dp * ( tm1 (1:kproma,jk,jrow) + tm1 (1:kproma,jk-1,jrow) )
    ztvh    (:,jk)=0.5_dp * ( ztvm1(:,jk) + ztvm1(:,jk-1) )
! ka_av_20160329+
! a workaround against corrupted virtual temperature
! SV: not needed anymore
!    ztvh = zth
! ka_av_20160329-
    zrhoh   (1:kproma,jk)=aphm1(1:kproma,jk)*&
          zmah(:,jk)&
          /(zrstar*ztvh(:,jk))
    zgmurhoh(1:kproma,jk)=1.87E-5_dp * (zth(1:kproma,jk)/273.15_dp)**.69_dp &
                          * zgvh(1:kproma,jk) * zrhoh(1:kproma,jk)
    zdiffmol(1:kproma,jk)=1.87E-5 * (zth(:,jk)/273.15)**.69/zrhoh(:,jk)
    DO jt=1,ntrac_gp
      zrhok(1:kproma,jk,jt)=4.17E-6_dp * SQRT(zth(1:kproma,jk)/273.15_dp)  &
                    * SQRT(zmah(1:kproma,jk)+zmah(1:kproma,jk)*zmah(1:kproma,jk)/zmm(jt))
    ENDDO
  ENDDO

  zrdp(1:kproma,1)=ztmst*grav1(1:kproma,nlev)/(aphm1(1:kproma,2)-aphm1(1:kproma,1))
  
  DO jk=2,nlev-1
    zla(1:kproma,jk)=-zrdp(1:kproma,jk)*(zinvpr*zgmurhoh(1:kproma,jk  )*zhh(1:kproma,jk  )) 
    zlb(1:kproma,jk)= zrdp(1:kproma,jk)*(zinvpr*zgmurhoh(1:kproma,jk  )*zhh(1:kproma,jk  ) + &
                           zinvpr*zgmurhoh(1:kproma,jk+1)*zhh(1:kproma,jk+1))
    zlc(1:kproma,jk)=-zrdp(1:kproma,jk)*(zinvpr*zgmurhoh(1:kproma,jk+1)*zhh(1:kproma,jk+1))
  ENDDO

  zlb(1:kproma,1)= zrdp(1:kproma,1)*zinvpr*zgmurhoh(1:kproma,2)*zhh(1:kproma,2)
  zlc(1:kproma,1)=-zrdp(1:kproma,1)*zinvpr*zgmurhoh(1:kproma,2)*zhh(1:kproma,2)

  zla(1:kproma,nlev)=-zrdp(1:kproma,nlev)*zinvpr*zgmurhoh(1:kproma,nlev)*zhh(1:kproma,nlev)
  zlb(1:kproma,nlev)= zrdp(1:kproma,nlev)*zinvpr*zgmurhoh(1:kproma,nlev)*zhh(1:kproma,nlev)

!-----------------------------
! call temperature part of vdiff_mol
  CALL edith_vdiff_mol_temp(kproma, kbdim, nlev, nlev-1 &
                    , tm1(:,:,jrow), time_step_len &
                    , zalpha, zmalpha &
                    , zla, zlb, zlc &
                    , ptte, status)
  IF (status /= 0) CALL error_bi('EDITh ERROR: Error in molecular diffusion','edith_vdiff_mol_temp')
  ttevdiffmol(:,:,jrow) = ptte
                    
!------------------------------
! call humidity part of vdiff_mol

  CALL edith_vdiff_mol_moist( kproma, kbdim, nlev, nlev-1 &
                    , qm1(:,:,jrow), time_step_len &
                    , zalpha, zmalpha &
                    , zla, zlb, zlc &
                    , pqte, status)
  IF (status /= 0) CALL error_bi('EDITh ERROR: Error in molecular diffusion','edith_vdiff_mol_moist')
  qtevdiffmol(:,:,jrow) = pqte
  
!------------------------------
! call tracer part of vdiff_mol
  IF (ntrac_gp>0) THEN     ! only if we have tracers
    DO jt=1,ntrac_gp
      IF (ti_gp(jt)%tp%meta%cask_i(I_vdiff) == ON) THEN
        CALL edith_vdiff_mol_trac( kproma, kbdim, nlev, nlev-1 &
                    , tm1(1:kproma,:,jrow), zpxtm1(1:kproma,:,jt), time_step_len &
                    , zalpha, zmalpha &
                    , zmah(1:kproma,:), zrhok(1:kproma,:,jt), zgvh(1:kproma,:), zrhoh(1:kproma,:) &
                    , zth(1:kproma,:), zrstar, zmm(jt), zexpan(jt), zhh(1:kproma,:), zrdp(1:kproma,:) &
                    , zutt(1:kproma,:,jt), status)
      END IF   ! vdiff==ON
      IF (status /= 0) CALL error_bi('EDITh ERROR: Error in molecular diffusion','edith_vdiff_mol_trac')
    END DO ! all tracers
    !scale overall mixing ratio to 1
    !   Sum of mixing ratios is normalized because the sum is not constrained.
    !   As long as N2 is not a tracer the original sum zsumx has to be used for this.
    zuttsum=0.0_dp
    DO jt=1,ntrac_gp
      IF (ti_gp(jt)%tp%meta%cask_i(I_vdiff) == ON) THEN
        zuttsum(1:kproma,:)=zuttsum(1:kproma,:)+zutt(1:kproma,:,jt)
        
      ENDIF
      
    ENDDO

    DO jt=1,ntrac_gp
      IF (ti_gp(jt)%tp%meta%cask_i(I_vdiff) == ON) THEN
        zutt(1:kproma,:,jt)=zutt(1:kproma,:,jt)/zuttsum(1:kproma,:)    ! <- standard
        !zutt(1:kproma,:,jt)=zutt(1:kproma,:,jt)*zsumx(1:kproma,:)/zuttsum(1:kproma,:)    ! alternative way
      !ENDIF
      
      pxtte(1:kproma,:,jt)=(zutt(1:kproma,:,jt)-zpxtm1(1:kproma,:,jt))/time_step_len
      xtte(1:kproma,:,jt)=xtte(1:kproma,:,jt)+pxtte(1:kproma,:,jt)
      ENDIF
      
    ENDDO
    
    
  END IF  ! #tracers>0

!------------------------------
! call u-wind part of vdiff_mol

  CALL edith_vdiff_mol_wind( kproma, kbdim, nlev, nlev-1 &
                    , um1(:,:,jrow), vm1(:,:,jrow), time_step_len &
                    , zalpha, zmalpha &
                    , zrdp, zgmurhoh, zhh &
                    , utevdiffmol(:,:,jrow), vtevdiffmol(:,:,jrow), status)
  vom_3d(1:kproma,1:nlev,jrow) = vom_3d(1:kproma,1:nlev,jrow) &
          + utevdiffmol(:,:,jrow)
  vol_3d(1:kproma,1:nlev,jrow) = vol_3d(1:kproma,1:nlev,jrow) &
          + vtevdiffmol(:,:,jrow)

!-------------------------------
! call frictional heating in vdiff_mol

  CALL edith_vdiff_mol_fric( kproma, kbdim, nlev, nlev-1 &
                    , apm1, um1(:,:,jrow), vm1(:,:,jrow) &
                    , zgmurhoh, grav1, cp1 &
                    , fricvdiffmol(:,:,jrow), status)
                    

! ka_av_12082015+
! is this really correct???
! op_pj_20180713: see my comment above
!!$     molco(:,:,jrow)= zgmurhoh(:,:)
! ka_av_12082015-
! ka_sv_20180704+
! ka_sv_20200131+
! output of wrong variable as diffusion coefficient
!     moldiffcoef(:,:,jrow)=zgmurhoh(:,:)
     moldiffcoef(:,:,jrow)=zdiffmol(:,:)
! ka_sv_20180704+

     IF (status /= 0) CALL error_bi('Error in molecular diffusion',substr)
!ka_av_20150825+
! ka_sv_20180122-

!!!!!!! ION DRAG FROM HAMMONIA
 IF (calc_iondrag.eq.1) THEN
            
  CALL edith_iondrag(jrow, kproma, kbdim, nlev, um1(1:kproma,:,jrow), &
                vm1(1:kproma,:,jrow), qm1(1:kproma,:,jrow) &
               ,geopot_3d(1:kproma,:,jrow), cp1, tteion(1:kproma,:,jrow), &
              zvom,zvol, &
              philat_2d(1:kproma,jrow),time_step_len &
!ka_sb_20180419+
              ,kp_data(1))
!ka_sb_20180419- 
   uteiondrag(1:kproma,1:nlev,jrow)=zvom
   vteiondrag(1:kproma,1:nlev,jrow)=zvol
   vom_3d(1:kproma,1:nlev,jrow) = vom_3d(1:kproma,1:nlev,jrow) &
          + uteiondrag(1:kproma,1:nlev,jrow)
   vol_3d(1:kproma,1:nlev,jrow) = vol_3d(1:kproma,1:nlev,jrow) &
          + vteiondrag(1:kproma,1:nlev,jrow)
 END IF
 
 !! -----------------------------
 !!  radiative cooling by NO
 
 NO_t0(:,:) = no_ptr(:,:,jrow)
 O_t0(:,:) = o_ptr(:,:,jrow)
 
 DO jk=1,nlev
   DO jl=1,kproma
     CALL edith_calc_nocooling(NO_t0(jl,jk),O_t0(jl,jk),apm1(jl,jk),tm1(jl,jk,jrow),cp1(jl,jk),NOcooling(jl,jk,jrow))
   END DO
 END DO
 
 !! -----------------------------
 !!  radiative cooling by CO2
 
 CO2_t0=co2_ptr(:,:,jrow)
 O3_t0=o3_ptr(:,:,jrow)
 O2_t0=o2_ptr(:,:,jrow)
 N2_t0=n2_ptr(:,:,jrow)
 
 CALL edith_calc_co2cooling(nproma,kbdim,nlev,CO2_t0, O3_t0, O2_t0, O_t0, N2_t0, amu1, cp1, &
           apm1,tm1(:,:,jrow),cossza_2d(:,jrow), &
           CO2nlte(:,:,jrow), CO2nir(:,:,jrow), &
           amatsave(:,:,:), &
           bmatsave(:,:,:),alsave(:,:),co2colxsave(:,:), &
           co2xsave(:,:) )
 
 
 amu(1:kproma,:,jrow) = amu1(1:kproma,:)
 cp(1:kproma,:,jrow) = cp1(1:kproma,:)
 g(1:kproma,:,jrow) = grav1(1:kproma,:)
 alt(1:kproma,:,jrow) = altitude(1:kproma,:)

  DEALLOCATE(altitude)
  
  tte_ptr(:,:,jrow) = tte_ptr(:,:,jrow) &
         + tteion(:,:,jrow) + ttevdiffmol(:,:,jrow) + fricvdiffmol(:,:,jrow) &
         + NOcooling(:,:,jrow) + CO2nlte(:,:,jrow) + CO2nir(:,:,jrow)
         
    
  qte_3d(1:kproma,:,jrow)=qte_3d(1:kproma,:,jrow)+pqte(1:kproma,:)
  
! ka_sv_20180702+  
  ! calculate Schuman-Runge-Continuum Photolysis of CO2 and O2
  IF (edith_phot) THEN
    ! flip vertical dimension -> no need to change routines from KASIMA
    DO jk=1,nlev
      o2flip(:,nlev-jk+1)=O2_t0(:,jk)
      tmflip(:,nlev-jk+1)=tm1(:,jk,jrow)
      pressflip(nlev-jk+1)=apm1(1,jk)
    END DO
    
    CALL edith_rcolsph( o2flip(1:kproma,1:nlev), nlev, acos(cossza_2d(1:kproma,jrow))*RTD, kproma, scolo2(1:kproma,1:nlev),&
                         sill(1:kproma,1:nlev), pressflip(1:nlev) )
    ! scolo2 and sill are from edith_rcolsph -> level 1 is at the bottom -> no need to flip dimension in the next call
    CALL edith_photsrc( scolo2(1:kproma,1:nlev), sill(1:kproma,1:nlev), tmflip(1:kproma,1:nlev), nlev, kproma, &
                   j_src_flip(1:kproma,1:nlev,:), 2, 0)
    DO jk=1,nlev
      j_src(1:kproma,jk,:)=j_src_flip(1:kproma,nlev-jk+1,:)
    END DO
    ! j_src is flipped into EMAC grid direction

    ! save photolysis rates in output
    JO2(:,:,jrow) = j_src(:,:,1)
    JCO2(:,:,jrow) = j_src(:,:,2)
    
    ! calculate tendencies
    DO jk=1,nlev
      DO jl=1,kproma
        ! calculate number densities
        call vmr_to_nd(n_o2,max(O2_t0(jl,jk),1.e-30),apm1(jl,jk),tm1(jl,jk,jrow))
        call vmr_to_nd(n_co2,max(CO2_t0(jl,jk),1.e-30),apm1(jl,jk),tm1(jl,jk,jrow))
        ! calculate tendency in molecules per second and convert it to mols per second
        call nd_to_vmr(JO2t(jl,jk,jrow),n_o2*JO2(jl,jk,jrow),apm1(jl,jk),tm1(jl,jk,jrow))
        call nd_to_vmr(JCO2t(jl,jk,jrow),n_co2*JCO2(jl,jk,jrow),apm1(jl,jk),tm1(jl,jk,jrow))
      END DO
    END DO
    
    ! tendency should be negative -> O2 and CO2 are destroyed
    JO2t(:,:,jrow)=-JO2t(:,:,jrow)
    JCO2t(:,:,jrow)=-JCO2t(:,:,jrow)
    
    ! add tendencies
    pxtte(1:kproma,:,idt_co2) = pxtte(1:kproma,:,idt_co2) + &
               JCO2t(1:kproma,:,jrow)
    pxtte(1:kproma,:,idt_o2) = pxtte(1:kproma,:,idt_o2) + &
               JO2t(1:kproma,:,jrow)
    pxtte(1:kproma,:,idt_o3p) = pxtte(1:kproma,:,idt_o3p) - &
               JO2t(1:kproma,:,jrow)*ed_o2_o3p - &
               JCO2t(1:kproma,:,jrow)
    pxtte(1:kproma,:,idt_o1d) = pxtte(1:kproma,:,idt_o1d) - &
               JO2t(1:kproma,:,jrow)*ed_o2_o1d
!    pxtte(1:kproma,:,idt_o1s) = pxtte(1:kproma,:,idt_o1s) - &
!               JO2t(1:kproma,:,jrow)*ed_o2_o1s     ! not included yet
    pxtte(1:kproma,:,idt_co) = pxtte(1:kproma,:,idt_co) - &
               JCO2t(1:kproma,:,jrow)
  END IF   ! edith_phot
  
! ka_sv_20180702-
  

  END SUBROUTINE edith_physc
  ! ====================================================================


  ! ====================================================================
  SUBROUTINE edith_read_nml_cpl(status, iou)
   
    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
    
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    INTEGER :: dummy

    NAMELIST /CPL/ edith_no, edith_co2, edith_o, edith_o3, edith_o2, edith_n2 &
!ka_sb_20180419+
                   , edith_kp
!ka_sb_20180419-
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='edith_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

! ka_sb_20180419+
!    IF (ed_solvar.eq.4) THEN ! variable solar conditions, kp needed
         IF (TRIM(edith_kp%cha) == '') THEN
            CALL info_bi('ERROR: empty channel name for kp')
            RETURN
         ELSE
            CALL info_bi('kp channel :'//edith_kp%cha)
         END IF
         IF (TRIM(edith_kp%obj) == '') THEN
            CALL info_bi('ERROR: empty channel object name for kp')
            RETURN
         ELSE
            CALL info_bi('kp object  :'//edith_kp%obj)
         END IF
!     END IF !ed_solvar=4
! ka_sb_20180419-



    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE edith_read_nml_cpl
  ! ====================================================================

! ------------------------------------------------------------

! TODO: SUBROUTINE HAS TO BE MOVED TO SMCL LAYER



! ------------------------------------------------------------

! TODO: SUBROUTINE HAS TO BE MOVED TO SMCL LAYER



! **********************************************************************
END MODULE messy_edith_e5
! **********************************************************************
