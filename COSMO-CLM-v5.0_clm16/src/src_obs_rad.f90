!+ Module for generating feedback files for satellite radiances.
!------------------------------------------------------------------------------

MODULE src_obs_rad

!------------------------------------------------------------------------------
!
! Description:
!   This module contains data and routines for calculating first guess for
!   satellite observations.
!
! Method:
!
! Current Code Owner: DWD, Andreas Messer
!  phone:  +49  69  8062 2438
!  fax:    +49  69  8062 3721
!  email:  andreas.messer@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! V4_26        2012/12/06 Andreas Messer
!  Initial Release
! V4_27        2013/03/19 Ulrich Schaettler
!  Introduced conditional compilation for Nudging to separate nudging parts
!  from SYNSAT parts.
! V4_28        2013/07/12 Ulrich Schaettler
!  Changed pointer assignment to standard assignment, because variable is
!   no pointer any more (SX accepted that gracefully, gfortran reports error)
! V4_29        2013/10/04 Astrid Kerkweg, Ulrich Schaettler
!  Unification of MESSy interfaces and COSMO Tracer structure
! V4_30        2013-11-08 Ulrich Schaettler
!  Changed intent attribute of variable nlev to INOUT in SR prepare_rttov_input
!  Changed intent attribute of variable comm to INOUT in SR p_bcast_rad_set
! V5_00_clm9   2016/05/11 H.-J. Panitz, IMK/KIT
!  Implemented F2003 IOMSG-mechanism for better namelist error messages. (UB)
!   adapted from COSMO_5.1
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE iso_fortran_env, ONLY: stderr => error_unit, stdout => output_unit

USE data_parameters, ONLY: iintegers,ireals

USE data_constants, ONLY: B1, B2W, B3, B4W, rdv, o_m_rdv

USE data_parallel, ONLY: &
  nproc,            & ! total number of processors
  num_compute,      & ! number of compute pe's
  nboundlines,      & ! no of overlapping boundary lines of the subdomains
  isubpos,          & ! positions of subdomains in total domain
  my_world_id,      & ! rank of this subdomain in the global    communicator
  my_cart_id,       & ! rank of this subdomain in the cartesian communicator
  icomm_world,      & ! global communicator
  icomm_cart,       & ! cartesian communicator
  imp_character,    & ! character type for MPI
  imp_reals,        & ! character type for MPI
  imp_integers,     & ! integer type for MPI
  imp_logical         ! logical type for MPI

USE data_modelconfig, ONLY: &
  ie,             & ! no of gridpoints in zonal direction (local)
  je,             & ! no of gridpoints in meridional direction (local)
  ke,             & ! no of gridpoints in vertical direction
  ie_tot,         & ! no of gridpoints in zonal direction
  je_tot,         & ! no of gridpoints in meridional direction
  ke_tot,         & ! no of gridpoints in vertical direction
  pollon,         & ! logitude of the rotated north pole
  pollat,         & ! latitude of the rotated north pole
  polgam,         & ! angle between the north poles of the system
  dlon,           & ! grid point distance in zonal direction
  dlat,           & ! grid point distance in meridonal direction
  startlat_tot,   & ! rot latitude lower left gridpoint 
  startlon_tot,   & ! rot longitude lower left gridpoint
  degrad,         & ! factor for transforming degree to rad
  raddeg,         & ! factor for transforming rad to degree
  dt,             & ! Length of Timestep
  idt_qv, idt_qc, idt_qi, idt_qs, idt_qg

#ifdef NUDGING
USE data_nudge_all, ONLY: &
  nolbc,            & ! no of grid rows at lateral boundaries where obs are neg.
  doromx            ! station heights
#endif

USE data_runcontrol, ONLY: &
  nstart,           & ! First timestep of forecast
  nstop,            & ! Last timestep of forecast
  ntstep,           & ! Last timestep of forecast
  nnow,             & ! field index for last timestep
  itype_calendar,   & ! for specifying the calendar used
  lseaice,          & ! 
  hstart,           & ! 
  hstop,            & ! 
  leps,             & ! 
  iepsmem,          & ! 
  nvers               ! 

USE data_io, ONLY: &
  ydate_ini,           & ! Start date of forecast
  yncglob_source,      & ! Start date of forecast
  yncglob_institution, & ! Start date of forecast
  root                   !

#ifdef NUDGING
USE data_obs_record, ONLY: &
  imdi,             & ! missing data indicator
  fdoro               ! scaling factor to vertical distances btw model
#endif

USE data_fields, ONLY: &
  p0,               & ! presure
  pp,               & ! presure
  ps,               & ! surface presure
  t,                & ! temperature
  t_g,              & ! surface temperature
  t_2m,             & ! 2m temperature
! qv,               & ! humidity
! qc,               & ! humidity
! qi,               & ! humidity
! qs,               & ! humidity
! qg,               & ! humidity
  qv_2m,            & ! 2m humidity
  u_10m,            & ! 10m wind u-component
  v_10m,            & ! 10m wind v-component
  rlat,             & ! rotated lon lat
  hsurf,            & ! surface height
  hhl,              & ! geommetricla height of half model levels
  soiltyp,          & ! soil type
  lseamask,         & ! land/sea mask
  fr_land,          & ! geommetricla height of half model levels
  sun_el,           & ! sun
  h_ice,            & ! 
  clc_sgs,          & ! 
  clw_con,          & ! 
  clc_con,          & ! 
  rho                 ! 

USE data_satellites, ONLY: &
  extrp_const,      &
  extrp_lin,        &
  extrp_clim,       &
  p_top,            & ! Top level pressure
  t_top,            & ! Top level temperature
  q_top,            & ! Top level humidity
  rcnw,             & ! Top level humidity
  iceshape,         & ! Top level humidity
  iwc2effdiam,      & ! Top level humidity
  lcon_clw,         & ! Top level humidity
  num_sensors         ! Number of sensors in use for synsat

#ifdef NUDGING
USE data_obs_lib_cosmo, ONLY: &
  rmdi,  &
  rmdich
#endif

USE src_tracer,       ONLY :  trcr_get, trcr_errorstr
USE data_tracer,      ONLY :  T_ERR_NOTFOUND

USE environment, ONLY: &
  model_abort         ! aborts the program

#ifdef NUDGING
USE src_obs_cdfin_util, ONLY: &
  obs_assign_gridpt   ,& ! Assign gridpoint to obs according lon/lat
  obs_assign_sort_node,& ! Assign node to gridpoint
  get_global_surface     ! get global surface
#endif

USE parallel_utilities, ONLY: &  
  distribute_values,& ! Distribute values to all PEs
  gather_field,     & ! 
  gather_values,    & ! 
  global_values       ! 

USE utilities, ONLY: &
  diff_minutes,     & ! Calculate time difference in minutes
  rla2rlarot,       & ! Convert lambda to rotated latlon
  phi2phirot          ! Convert phi to rotated latlon

#ifdef RTTOV10
USE mo_rad, ONLY: &
  t_rad_set,             & ! Type of an radiance set 
  t_radv,                & ! Type of a set of radiances
  rad_set,               & ! Array of radiance sets 
  read_tovs_obs_chan_nml,& ! Subroutine to read "TOVS_OBS_CHAN" Namelist
  read_satpp_feedbk     ,& ! Subroutine to read a satpp file
  link_rad              ,& ! Combine Information from Namelists with Satpp
  sat_id2name,           & ! Get the satellite name for a given bufr id
  destruct,              & ! Release memory  
  assignment (=),        & !
  USE_PASSIVE              ! Passive bit

USE mo_rttov_ifc, ONLY:  &
  rttov_ifc_version,     & !
  rttov_fill_input,      & !
  rttov_set_input,       & !
  rttov_k_ifc,           & !
  rttov_ifc_errmsg         !
#endif

#ifdef NUDGING
USE mo_fdbk, ONLY: &
  t_fdbk,                & ! Type of feedback file
  setup_fdbk,            & ! setup a fdbk files variables
  create_fdbk,           & ! create the fdbk file
  close_fdbk,            & ! create the fdbk file
  open_fdbk_read,        & ! create the fdbk file
  open_fdbk_write,       & ! create the fdbk file
  add_verification         ! create the fdbk file

USE mo_fdbk_tables, ONLY: &
  vt_firstguess,          &
  vt_analysis,            &
  rc_ass,                 &
  OT_RAD,                 &
  VN_RAWBT,               &
  ST_SEA,                 &
  ST_ICE,                 &
  ST_LAND,                &
  ST_HIGHLAND,            &
  ST_MISMATCH,            &
  MS_LAND,                &
  MS_SEA,                 &
  MS_ICE,                 &
  MS_NO_ICE,              &
  MS_SNOW,                &
  MS_NO_SNOW,             &
  ST_ACCEPTED,            &
  ST_PASSIVE,             &
  ST_PAS_REJ,             &
  ST_DISMISS,            &
  FL_TIME,                &
  FL_AREA,                &
  FL_SURF,                &
  FL_OBSTYPE,                &
  FL_OPERATOR

USE mo_fdbk_cosmo, ONLY: &
  t_account,             &
  t_acc_header,          &
  t_acc_body,            &
  write_report

USE mo_netcdf_param  
#endif

IMPLICIT NONE

PRIVATE

#ifdef NUDGING
PUBLIC :: input_obs_radctl
PUBLIC :: input_obs_satpp
PUBLIC :: calc_obs_satpp
#endif
PUBLIC :: t_sensor
PUBLIC :: nsensors
PUBLIC :: sensors
PUBLIC :: prepare_rttov_input

!==============================================================================
! MPI Include
!------------------------------------------------------------------------------
#include "mpif.h"

!==============================================================================
! Module Parameters
!------------------------------------------------------------------------------

  CHARACTER(LEN=*), PARAMETER :: &
    yurad = 'YURAD'                ! iradiance obs diagnostic file
  INTEGER, PARAMETER ::                &
    max_satpp_files         = 16    ! maximum number of satpp files
  REAL (KIND=ireals), PARAMETER :: &
    tmax = 399.0_ireals,           & ! [K]
    tmin = 91.0_ireals,            & ! [K]
    qmax = 0.372_ireals,           & ! [kg/kg]
    qmin = TINY(0._ireals)            ! [kg/kg]

#ifdef NUDGING
  ! The following array maps a check to its assigned state  
  INTEGER (KIND=iintegers), PARAMETER :: &
    map_flg2state(0:21) = (/ST_PASSIVE,  & ! FL_OBSTYPE
                            ST_DISMISS,  & ! FL_BLACKLIST
                            ST_PAS_REJ,  & ! FL_SUSP_LOCT
                            ST_DISMISS,  & ! FL_TIME
                            ST_DISMISS,  & ! FL_AREA
                            ST_PASSIVE,  & ! FL_HEIGHT
                            ST_PASSIVE,  & ! FL_SURF
                            ST_PASSIVE,  & ! FL_CLOUD
                            ST_PAS_REJ,  & ! FL_PRACTICE
                            ST_DISMISS,  & ! FL_DATASET
                            ST_PASSIVE,  & ! FL_REDUNDANT
                            ST_PAS_REJ,  & ! FL_FLIGHTTRACK
                            ST_DISMISS,  & ! FL_MERGE
                            ST_DISMISS,  & ! FL_THIN
                            ST_PAS_REJ,  & ! FL_RULE
                            ST_PAS_REJ,  & ! FL_OBS_ERR
                            ST_PASSIVE,  & ! FL_GROSS
                            ST_PASSIVE,  & ! FL_NO_BIASCOR
                            ST_PASSIVE,  & ! FL_FG
                            ST_PAS_REJ,  & ! FL_NO_OBS
                            ST_PAS_REJ,  & ! FL_OPERATOR
                            ST_DISMISS/)   ! FL_FG_LBC

  ! Category and subcategory
  INTEGER (KIND=iintegers), PARAMETER :: &
    fdbk_rad_cat    = 21, & ! value to use for category field in fdbk file
    fdbk_rad_subcat = -1    ! value to use fur sub_cat.. field in fdbk file


  INTEGER (KIND=iintegers), PARAMETER :: &
    fdbk_height_method = 2  ! Which method to use for calculating the height
                            ! of an obs:
                            ! 1 = max of temp jacobian
                            ! 2 = mean of jacobian weighted by temp/hum
                            ! 3 = like (2) but squared sensitivities
  REAL (KIND=ireals), PARAMETER :: &
    fdbk_e_bg_t  = 0.5_ireals, & ! Typical background temp error
    fdbk_e_bg_rh = 0.1_ireals    ! Typical background hum error
#endif

!==============================================================================
! Module Types
!------------------------------------------------------------------------------

#ifdef RTTOV10
  TYPE t_rad_container
    TYPE (t_radv) :: &
      radv                 ! The radv container 
    INTEGER   (KIND=iintegers), ALLOCATABLE  :: &
      ntstep(:),    & ! time step when to calculate the radiance
      iob_tot(:),   & ! zonal index (total domain) of assigned gridpt
      job_tot(:),   & ! meridional index (total domain) of assigned gridpt
      rttov_chan(:)   ! rs%chan channels mapped to loaded rttov channels
    REAL      (KIND=ireals), POINTER :: &
      plevel(:,:) => NULL() ! height level assigned to one channel and fov
  END TYPE t_rad_container
#endif      

  TYPE t_sensor
    INTEGER (KIND=iintegers) :: &
      id(3)                 ! RTTOV instrument triplet (platform,sat,instr)
    INTEGER (KIND=iintegers), POINTER :: &
      channels(:) => NULL() ! RTTOV channels for this instrument
    LOGICAL :: &
      addcloud              ! Use cloud information for this sensor
  END TYPE

  TYPE t_time
    INTEGER (KIND=iintegers) :: &
      year,  &
      month, &
      day,   &
      hour,  &
      minute
  END TYPE

!==============================================================================
! Module Variables
!------------------------------------------------------------------------------

 INTEGER  (KIND=iintegers), SAVE :: &
    nsensors                          ! Number of sensors
  TYPE (t_sensor), SAVE, ALLOCATABLE, TARGET :: &
    sensors(:)                        ! Array of sensors
  character(len=128), save ::         &
    input_path                        ! path to input files
  character(len=128), save ::         &
    satpp_files(max_satpp_files)      ! satpp files to read in

#ifdef RTTOV10    
  INTEGER  (KIND=iintegers), SAVE :: &
    nurad                              ! radiance obs information file

  INTEGER  (KIND=iintegers), SAVE, ALLOCATABLE :: &
      nobspe_gath(:,:),             & ! number of obs gathered (pe,radv) (local)
      nobspe_tot(:,:),              & ! number of total obs   (pe,radv) (global)
      fdbk_offset(:)                  ! offset within feedback fie (pe 0 only)
  INTEGER  (KIND=iintegers),SAVE ::   &
    usesurftype = 1                   ! Which surfacetype to use: 0,1 = cosmo,satpp
  TYPE (t_rad_container), SAVE, ALLOCATABLE, TARGET :: &
    obs_rad(:)                        ! Array of radiances for each satpp (local)

#ifdef NUDGING
  TYPE (t_fdbk), SAVE, ALLOCATABLE :: &
    fdbk(:)                           ! array of fdbk files
  TYPE (t_time),SAVE :: &
    reftime   ! reference time used in varous places
#endif

#endif

!==============================================================================
! Module Interfaces
!------------------------------------------------------------------------------

#ifdef NUDGING
#ifdef RTTOV10
  INTERFACE construct
    MODULE PROCEDURE construct_rad_container
  END INTERFACE construct

  INTERFACE destruct
    MODULE PROCEDURE destruct_rad_container
  END INTERFACE destruct  
  INTERFACE p_bcast
    MODULE PROCEDURE p_bcast_rad_set
  END INTERFACE
#endif  
#endif

!==============================================================================

CONTAINS

#ifdef NUDGING
#ifdef RTTOV10
!==============================================================================
!+ Procedure in "src_obs_sat" for constructing the t_obs_container type
!------------------------------------------------------------------------------
SUBROUTINE construct_rad_container(c,inobs, ierrstat, radv)
!------------------------------------------------------------------------------
!
! Description:
!   This subroutine allocates the memory for rad derived type
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)       ::        &
    inobs           ! number of observables
  TYPE(t_radv),   INTENT (IN),OPTIONAL            ::        &
    radv            ! radv to use as template
  TYPE(t_rad_container),   INTENT (OUT)         ::        &
    c               ! the container to construct
  INTEGER   (KIND=iintegers),   INTENT (OUT)      ::        &
    ierrstat        ! error status variable
! Local variables
  INTEGER   (KIND=iintegers) :: &
    n_instr, n_chan, n_lev, n_sens
    
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE construct_rad_container
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Allocate Memory
!-------------------------------------------------------------------------------

  ALLOCATE(c%ntstep(inobs),     &
           c%iob_tot(inobs),    &
           c%job_tot(inobs),    &
           c%radv%date(inobs),  &
           c%radv%time(inobs),  &
           c%radv%dlat(inobs),  &
           c%radv%dlon(inobs),  &
           c%radv%fov (inobs),  &
           c%radv%stzen (inobs),&
           c%radv%stazi (inobs),&
           stat = ierrstat)

  IF (PRESENT(radv) .AND. ierrstat == 0) THEN
    n_instr = radv%i%n_instr
    n_chan  = radv%i%n_chan
    n_sens  = radv%i%n_sens
    n_lev   = radv%n_lev

    ALLOCATE(c%radv%landfr(n_sens,inobs), &
             c%radv%stype(n_sens,inobs),  &
             c%radv%shgt(n_sens,inobs),   &
             c%radv%valid(n_chan,inobs),   &
             c%radv%bt_obs(n_chan,inobs),  &
             stat = ierrstat)

    IF (ASSOCIATED(radv%obsnum) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%obsnum(inobs),stat=ierrstat)
    ENDIF

    IF (ASSOCIATED(radv%date_d) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%date_d(inobs),stat=ierrstat)
    ENDIF
  
    IF (ASSOCIATED(radv%time_d) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%time_d(inobs),stat=ierrstat)
    ENDIF  

    IF (ASSOCIATED(radv%scanl) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%scanl(inobs),stat=ierrstat)
    ENDIF
  
    IF (ASSOCIATED(radv%sunzen) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%sunzen(inobs),stat=ierrstat)
    ENDIF  

    IF (ASSOCIATED(radv%sunazi) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%sunazi(inobs),stat=ierrstat)
    ENDIF  

    IF (ASSOCIATED(radv%center) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%center(inobs),stat=ierrstat)
    ENDIF
  
    IF (ASSOCIATED(radv%subcenter) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%subcenter(inobs),stat=ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%mdlsfc) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%mdlsfc(inobs),stat=ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%state) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%state(n_chan,inobs),stat=ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%flags) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%flags(n_chan,inobs),stat=ierrstat)
    ENDIF
    
!    IF (ASSOCIATED(radv%check) .AND. ierrstat == 0) THEN
!      ALLOCATE(c%radv%check(n_chan,inobs),stat=ierrstat)
!    ENDIF

    IF (ASSOCIATED(radv%bt_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%bt_fg(n_chan,inobs),stat=ierrstat)
    ENDIF

    IF (ASSOCIATED(radv%bt_bcor) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%bt_bcor(n_chan,inobs),stat=ierrstat)
    ENDIF

    IF (ASSOCIATED(radv%bcor_) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%bcor_(n_chan,inobs),stat=ierrstat)
    ENDIF  
    
    IF (ASSOCIATED(radv%p) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%p(inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%h_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%h_fg(n_instr,inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%t_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%t_fg(n_lev,inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%q_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%q_fg(n_lev,inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%t2m) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%t2m(inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%q2m) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%q2m(inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%ps_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%ps_fg(inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%ts_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%ts_fg(inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%u10_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%u10_fg(inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%v10_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%v10_fg(inobs),stat = ierrstat)
    ENDIF
    
    IF (ASSOCIATED(radv%v10_abs_fg) .AND. ierrstat == 0) THEN
      ALLOCATE(c%radv%v10_abs_fg(inobs),stat = ierrstat)
    ENDIF
  ENDIF
!-------------------------------------------------------------------------------
!- Section 2: Assign values
!-------------------------------------------------------------------------------
  IF (ierrstat == 0) THEN
    c%radv%n_rec      = inobs

    IF (PRESENT(radv)) THEN
      c%radv%filename   = radv%filename
      c%radv%n_lev      = radv%n_lev
      c%radv%model_date = radv%model_date
      c%radv%i          = radv%i
    ENDIF  
  ENDIF
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
END SUBROUTINE construct_rad_container

!==============================================================================
!+ Procedure in "src_obs_sat" for destructing the t_obs_container type
!------------------------------------------------------------------------------
SUBROUTINE destruct_rad_container(c)
!------------------------------------------------------------------------------
!
! Description:
!   This subroutine allocates the memory for rad derived type
!------------------------------------------------------------------------------

! Parameters
  TYPE(t_rad_container),   INTENT (INOUT)         ::        &
    c             ! the container to destruct

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE destruct_rad_container
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Release Memory
!-------------------------------------------------------------------------------
  IF (ALLOCATED(c%ntstep))       DEALLOCATE(c%ntstep)
  IF (ALLOCATED(c%iob_tot))      DEALLOCATE(c%iob_tot)
  IF (ALLOCATED(c%job_tot))      DEALLOCATE(c%job_tot)
  IF (ALLOCATED(c%rttov_chan))   DEALLOCATE(c%rttov_chan)
  IF (ASSOCIATED(c%plevel))       DEALLOCATE(c%plevel)

  CALL destruct(c%radv)
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
END SUBROUTINE destruct_rad_container

!==============================================================================
!+ Procedure in "src_obs_sat" for distributing the rad container
!------------------------------------------------------------------------------

SUBROUTINE gather_obs_rad_container (c_loc, start_loc, end_loc, &
                                     c_tot, start_tot, end_tot, &
                                     ireceiver, ierrstat,       &
                                     numreceived)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine collects radiances in one place
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)   :: &
    start_loc,    & ! First Record to distribute
    end_loc,      & ! Last Record        
    start_tot,    & ! First Record to receive
    ireceiver       ! Node which shall receive the data
  TYPE(t_rad_container), INTENT(IN), TARGET :: &
    c_loc     ! Input radiances
  TYPE(t_rad_container), INTENT(INOUT) :: &
    c_tot     ! Collected radiances
  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable
  INTEGER   (KIND=iintegers),   INTENT (INOUT), OPTIONAL :: &
    numreceived(num_compute)      ! number of records received  
  INTEGER   (KIND=iintegers),   INTENT (OUT)   :: &
    end_tot         ! Last Record received 

! Local variables
  INTEGER   (KIND=iintegers) :: &
    inumelems(num_compute),    &  ! Number of obs per node
    ioffsets(num_compute),     &  ! Offsets in target vector
    icnt                            ! Loop counter
  CHARACTER(LEN=256) :: &
    ymsg            ! error message

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE gather_obs_rad_container
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Prepare 
!-------------------------------------------------------------------------------
  inumelems = 0
  inumelems(my_cart_id+1) = end_loc - start_loc + 1
  
  CALL gather_values(inumelems(my_cart_id+1), inumelems(:), 1, &
                     num_compute, imp_integers, ireceiver, icomm_cart, &
                     ymsg, ierrstat)
  IF (ierrstat /= 0) CALL model_abort(my_cart_id,99903, 'gather_values', &
                                   'gather_obs_rad_container', ierrstat)

  IF (PRESENT(numreceived) .AND. my_cart_id == ireceiver) THEN
    numreceived = inumelems
  ENDIF  

!  IF (my_cart_id == ireceiver) THEN
!    WRITE (stderr,*) 'ID',my_cart_id,': ', ireceiver, ':', inumelems(:,1)
!  ENDIF  

  ioffsets(1) = 0
  DO icnt = 2,num_compute
    ioffsets(icnt) = ioffsets(icnt-1) + inumelems(icnt-1)
  ENDDO

  end_tot = start_tot + SUM(inumelems) - 1

  CALL gatherv_integers(c_loc % ntstep , c_tot % ntstep)
  CALL gatherv_integers(c_loc % iob_tot, c_tot % iob_tot)
  CALL gatherv_integers(c_loc % job_tot, c_tot % job_tot)

  CALL gatherv_reals2dp(c_loc %plevel, c_tot%plevel)
  CALL gatherv_integersp(c_loc%radv%obsnum, c_tot%radv%obsnum)
  CALL gatherv_integersp(c_loc%radv%date,   c_tot%radv%date)
  CALL gatherv_integersp(c_loc%radv%time,   c_tot%radv%time)
  CALL gatherv_integersp(c_loc%radv%date_d, c_tot%radv%date_d)
  CALL gatherv_integersp(c_loc%radv%time_d, c_tot%radv%time_d)
  CALL gatherv_realsp   (c_loc%radv%dlat,   c_tot%radv%dlat)
  CALL gatherv_realsp   (c_loc%radv%dlon,   c_tot%radv%dlon)
  CALL gatherv_integersp(c_loc%radv%fov,    c_tot%radv%fov)
  CALL gatherv_integersp(c_loc%radv%scanl,  c_tot%radv%scanl)
  CALL gatherv_realsp(c_loc%radv%stzen, c_tot%radv%stzen)
  CALL gatherv_realsp(c_loc%radv%stazi, c_tot%radv%stazi)
  CALL gatherv_realsp(c_loc%radv%sunzen, c_tot%radv%sunzen)
  CALL gatherv_realsp(c_loc%radv%sunazi, c_tot%radv%sunazi)
  CALL gatherv_reals2dp(   c_loc%radv%landfr,   c_tot%radv%landfr)
  CALL gatherv_integers2dp(c_loc%radv%stype,    c_tot%radv%stype)
  CALL gatherv_reals2dp(   c_loc%radv%shgt,     c_tot%radv%shgt)
  CALL gatherv_integersp(c_loc%radv%center,     c_tot%radv%center)
  CALL gatherv_integersp(c_loc%radv%subcenter,  c_tot%radv%subcenter)
  CALL gatherv_integersp(c_loc%radv%mdlsfc,     c_tot%radv%mdlsfc)
  CALL gatherv_integers2dp(c_loc%radv%state,    c_tot%radv%state)
  CALL gatherv_integers2dp(c_loc%radv%flags,    c_tot%radv%flags)
!  CALL gatherv_integers2dp(c_loc%radv%check,    c_tot%radv%check)
  CALL gatherv_logical2dp( c_loc%radv%valid,    c_tot%radv%valid)
  CALL gatherv_reals2dp(   c_loc%radv%bt_obs,   c_tot%radv%bt_obs)
  CALL gatherv_reals2dp(c_loc%radv%bt_fg, c_tot%radv%bt_fg)
  CALL gatherv_reals2dp(c_loc%radv%bt_bcor, c_tot%radv%bt_bcor)
  CALL gatherv_reals2dp(c_loc%radv%bcor_, c_tot%radv%bcor_)
  CALL gatherv_realsp(c_loc%radv%p, c_tot%radv%p)
  CALL gatherv_reals2dp(c_loc%radv%h_fg, c_tot%radv%h_fg)
  CALL gatherv_reals2dp(c_loc%radv%t_fg, c_tot%radv%t_fg)
  CALL gatherv_reals2dp(c_loc%radv%q_fg, c_tot%radv%q_fg)
  CALL gatherv_realsp(c_loc%radv%t2m, c_tot%radv%t2m)
  CALL gatherv_realsp(c_loc%radv%q2m, c_tot%radv%q2m)
  CALL gatherv_realsp(c_loc%radv%ps_fg, c_tot%radv%ps_fg)
  CALL gatherv_realsp(c_loc%radv%ts_fg, c_tot%radv%ts_fg)
  CALL gatherv_realsp(c_loc%radv%u10_fg, c_tot%radv%u10_fg)
  CALL gatherv_realsp(c_loc%radv%v10_fg, c_tot%radv%v10_fg)
  CALL gatherv_realsp(c_loc%radv%v10_abs_fg, c_tot%radv%v10_abs_fg)

CONTAINS

  SUBROUTINE gatherv_integers(loc,tot)
  INTEGER (KIND=iintegers), INTENT(IN) :: loc(:)
  INTEGER (KIND=iintegers), INTENT(INOUT) :: tot(:)
    CALL MPI_GATHERV(loc(start_loc), inumelems(my_cart_id+1), imp_integers, &
                     tot(start_tot), inumelems(:), ioffsets(:), imp_integers, &
                     ireceiver, icomm_cart, ierrstat)
    IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
                                        'gather_obs_rad_container', ierrstat)
  END SUBROUTINE

  SUBROUTINE gatherv_integers2d(loc,tot)
  INTEGER (KIND=iintegers), INTENT(IN) :: loc(:,:)
  INTEGER (KIND=iintegers), INTENT(INOUT) :: tot(:,:)
    INTEGER (KIND=iintegers) :: m
    m = SIZE(loc,1)

    CALL MPI_GATHERV(loc(1,start_loc), m * inumelems(my_cart_id+1), imp_integers, &
                     tot(1,start_tot), m * inumelems(:), m * ioffsets(:), imp_integers, &
                     ireceiver, icomm_cart, ierrstat)
    IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
                                        'gather_obs_rad_container', ierrstat)
  END SUBROUTINE

  SUBROUTINE gatherv_integersp(loc,tot)
  INTEGER (KIND=iintegers), INTENT(IN), POINTER :: loc(:)
  INTEGER (KIND=iintegers), INTENT(INOUT), POINTER :: tot(:)
    IF (ASSOCIATED(loc)) THEN
      CALL MPI_GATHERV(loc(start_loc), inumelems(my_cart_id+1), imp_integers, &
                       tot(start_tot), inumelems(:), ioffsets(:), imp_integers, &
                       ireceiver, icomm_cart, ierrstat)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
                                          'gather_obs_rad_container', ierrstat)
    ENDIF  
  END SUBROUTINE
  
  SUBROUTINE gatherv_integers2dp(loc,tot)
  INTEGER (KIND=iintegers), INTENT(IN), POINTER :: loc(:,:)
  INTEGER (KIND=iintegers), INTENT(INOUT), POINTER :: tot(:,:)
    INTEGER (KIND=iintegers) :: m
    IF (ASSOCIATED(loc)) THEN
      m = SIZE(loc,1)

      CALL MPI_GATHERV(loc(1,start_loc), m * inumelems(my_cart_id+1), imp_integers, &
                       tot(1,start_tot), m * inumelems(:), m * ioffsets(:), imp_integers, &
                       ireceiver, icomm_cart, ierrstat)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
                                          'gather_obs_rad_container', ierrstat)
    ENDIF  
  END SUBROUTINE

  SUBROUTINE gatherv_realsp(loc,tot)
  REAL (KIND=ireals), INTENT(IN), POINTER    :: loc(:)
  REAL (KIND=ireals), INTENT(INOUT), POINTER :: tot(:)
    IF (ASSOCIATED(loc)) THEN
      CALL MPI_GATHERV(loc(start_loc), inumelems(my_cart_id+1), imp_reals, &
                       tot(start_tot), inumelems(:), ioffsets(:), imp_reals, &
                       ireceiver, icomm_cart, ierrstat)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
                                          'gather_obs_rad_container', ierrstat)
    ENDIF  
  END SUBROUTINE

  SUBROUTINE gatherv_reals2dp(loc,tot)
  REAL (KIND=ireals), INTENT(IN), POINTER  :: loc(:,:)
  REAL (KIND=ireals), INTENT(INOUT), POINTER :: tot(:,:)
    INTEGER (KIND=iintegers) :: m
    IF (ASSOCIATED(loc)) THEN
      m = SIZE(loc,1)

      CALL MPI_GATHERV(loc(1,start_loc), m * inumelems(my_cart_id+1), imp_reals, &
                       tot(1,start_tot), m*inumelems(:), m*ioffsets(:), imp_reals, &
                       ireceiver, icomm_cart, ierrstat)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
                                          'gather_obs_rad_container', ierrstat)
    ENDIF  
  END SUBROUTINE

  SUBROUTINE gatherv_logical2dp(loc,tot)
  LOGICAL , INTENT(IN), POINTER :: loc(:,:)
  LOGICAL , INTENT(INOUT), POINTER       :: tot(:,:)
    INTEGER (KIND=iintegers) :: m
    IF (ASSOCIATED(loc)) THEN
      m = SIZE(loc,1)

      CALL MPI_GATHERV(loc(1,start_loc), m * inumelems(my_cart_id+1), imp_logical, &
                       tot(1,start_tot), m * inumelems(:), m * ioffsets(:), imp_logical, &
                       ireceiver, icomm_cart, ierrstat)
      IF (ierrstat /= 0) CALL model_abort(my_cart_id,99904, 'MPI_GATHERV', &
                                          'gather_obs_rad_container', ierrstat)
    ENDIF
  END SUBROUTINE

END SUBROUTINE gather_obs_rad_container 


!==============================================================================
!+ Procedure in "src_obs_sat" for broadcasting the t_rad_set type
!------------------------------------------------------------------------------

SUBROUTINE p_bcast_rad_set(buffer,source,comm)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine broadcasts a t_Rad_set structure to all PEs
!
!------------------------------------------------------------------------------
! Parameters
  TYPE(t_rad_set),   INTENT(INOUT) :: buffer(:)
  INTEGER,           INTENT(IN)    :: source
  INTEGER, OPTIONAL, INTENT(INOUT) :: comm

! Local variables
  INTEGER :: iradset, real_type
  INTEGER :: lcom, errorcode
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE p_bcast_rad_set
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Broadcast derived type 
!-------------------------------------------------------------------------------
 
  lcom = MPI_COMM_WORLD ; IF(PRESENT(comm)) lcom = comm

  CALL MPI_Bcast(buffer,SIZE(buffer) * SIZE(TRANSFER(buffer(1),(/' '/))), & 
                 MPI_BYTE, source,lcom,errorcode)

  IF (errorcode /= MPI_SUCCESS) THEN
    PRINT *, 'MPI ERROR in MPI_Bcast: ', errorcode
    STOP 'MPI ERROR'
  ENDIF

!-------------------------------------------------------------------------------
!- Section 2: Broadcast all allocated arrays 
!-------------------------------------------------------------------------------
  DO iradset=1,SIZE(buffer)
    IF (ASSOCIATED(buffer(iradset)%chan)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%chan(buffer(iradset)%n_chan))

      CALL MPI_Bcast(buffer(iradset)%chan, buffer(iradset)%n_chan, &
                     imp_integers, source, lcom,errorcode)
    ENDIF

    IF (ASSOCIATED(buffer(iradset)%ichan)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%ichan(buffer(iradset)%n_chan))

      CALL MPI_Bcast(buffer(iradset)%ichan, buffer(iradset)%n_chan, &
                     imp_integers, source, lcom,errorcode)
    ENDIF

    IF (ASSOCIATED(buffer(iradset)%band)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%band(buffer(iradset)%n_chan))

        CALL MPI_Bcast(buffer(iradset)%band, buffer(iradset)%n_chan, &
                      imp_integers, source, lcom,errorcode)
      ENDIF

    IF (ASSOCIATED(buffer(iradset)%flag)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%flag(buffer(iradset)%n_chan))

      CALL MPI_Bcast(buffer(iradset)%flag, buffer(iradset)%n_chan, &
                     imp_integers, source, lcom,errorcode)
      ENDIF

    IF (ASSOCIATED(buffer(iradset)%var)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%var(buffer(iradset)%n_chan))

      IF (KIND(buffer(1)%var) == 4) THEN
        real_type = MPI_REAL4
      ELSE
        real_type = MPI_DOUBLE_PRECISION
      ENDIF

      CALL MPI_Bcast(buffer(iradset)%var, buffer(iradset)%n_chan, &
                     real_type, source, lcom,errorcode)
    ENDIF
   
    IF (associated(buffer(iradset)%sensor_instr)) THEN
      IF (my_world_id /= source) &
        ALLOCATE(buffer(iradset)%sensor_instr(buffer(iradset)%n_sens))
      
      CALL MPI_Bcast(buffer(iradset)%sensor_instr, buffer(iradset)%n_sens, &
                     imp_integers, source, lcom,errorcode)
    ENDIF
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE p_bcast_rad_set
!==============================================================================

!==============================================================================
!+ Procedure in "src_obs_sat" for updating the state of an observable
!------------------------------------------------------------------------------
SUBROUTINE change_use_rpt_rad(radv,ichan,irpt,flag)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine updates the state of an report in a radv structure
!
!------------------------------------------------------------------------------

! Parameters
  TYPE(t_radv), INTENT(INOUT) :: &
    radv ! input radv structure
  INTEGER (KIND=iintegers), INTENT(IN) :: &
    ichan, & ! channel
    irpt,  & ! report index
    flag     ! the flag of the check triggered
! Local Variabes
  INTEGER (KIND=iintegers) :: &
    state   ! the new state
    
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE change_use_rpt_rad
!-------------------------------------------------------------------------------
  state = map_flg2state(flag)

  IF (radv%state(ichan,irpt) < state) THEN
    radv%state(ichan,irpt) = state
  ENDIF

  radv%flags(ichan,irpt) = IBSET(radv%flags(ichan,irpt),flag)
!-------------------------------------------------------------------------------
!- End SUBROUTINE change_use_rpt_rad
!-------------------------------------------------------------------------------

END SUBROUTINE change_use_rpt_rad
!==============================================================================
#endif

!==============================================================================
!+ Procedure in "src_obs_sat" for the input of NAMELIST obs_radctl
!------------------------------------------------------------------------------

SUBROUTINE input_obs_radctl (nuspecif, nuin, ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the NAMELIST file for satellite
!   observations processing. This file contains multiple NAMELIST groups 
!   TOVS_OBS_CHAN describing a single satellite dataset.
!
! Method:
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    nuspecif,     & ! Unit number for protocolling the task
    nuin            ! Unit number for Namelist INPUT file

  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable
! Local variables
#ifdef RTTOV10
  TYPE      (t_rad_set) :: &
    tmp             ! Temporary buffer for an rad set

  INTEGER   (KIND=iintegers) :: &
    iradset,      & ! counting variable for radiance set
    istat,        & ! result of calls to mo_rad subroutines
    i               ! couting variable

!HJP Begin 2016-05-11
  CHARACTER(LEN=250)         :: iomsg_str
!HJP End   2016-05-11

! Definition of the namelist group
  NAMELIST /TOVS_OBS/ input_path, satpp_files, usesurftype
#endif  
    
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE input_obs_radctl
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Initialize variables
!-------------------------------------------------------------------------------
  input_path(:) = ''
  satpp_files(:) = ''

!-------------------------------------------------------------------------------
!- Section 2: Read in namelist group TOVS_OBS and all following groups i
!-            TOVS_OBS_CHANs
!-------------------------------------------------------------------------------
#ifdef RTTOV10
IF (my_world_id == 0) THEN
  iradset  = 0
  istat    = 0

!HJP Begin 2016-05-11
! read(nuin,nml=TOVS_OBS, iostat=istat)
  iomsg_str(:) = ' '
  READ (nuin, nml=TOVS_OBS, IOSTAT=istat, IOMSG=iomsg_str)

!HJP END   2016-05-11

  IF (istat /= 0) THEN
    ierrstat = 99001
    PRINT *, ' ERROR    *** Reading namelist group TOVS_OBS failed. *** '
!HJP Begin 2016-05-11
    WRITE (*,'(A,A)') 'Namelist-ERROR LMGRID: ', TRIM(iomsg_str)
!HJP END   2016-05-11
    RETURN
  ENDIF
  
  DO WHILE (istat == 0)
    call read_tovs_obs_chan_nml(nuin, tmp, istat, iradset + 1)

    IF (istat /= 0) THEN ! TODO Check if last namelist group or error in namelist parameters
      EXIT
    ENDIF  

    IF (iradset + 1 > SIZE(rad_set)) THEN
      ierrstat = 99002
      PRINT *, ' ERROR    *** Maximum number of radiance sets reached. *** '
      RETURN
    ENDIF  

    iradset = iradset + 1
    rad_set(iradset) = tmp 
  ENDDO
ENDIF

!-------------------------------------------------------------------------------
!- Section 3: Evaluate parameters
!-------------------------------------------------------------------------------

i = LEN_TRIM(input_path)
IF (i > 0 .AND. input_path(i:i) /= '/') THEN
    input_path = TRIM(input_path) // '/'
ENDIF    
!-------------------------------------------------------------------------------
!- Section 4: Distribute variables to all nodes
!-------------------------------------------------------------------------------

IF (nproc > 1) THEN
  call distribute_values(input_path, len(input_path), 0, &
                         imp_character, icomm_world)
  call distribute_values(usesurftype, 1, 0, imp_integers, icomm_world)
  DO i = 1,size(satpp_files)
    call distribute_values(satpp_files(i), len(satpp_files(i)), 0, &
                           imp_character, icomm_world)
  ENDDO

  call p_bcast(rad_set, 0,icomm_world)
ENDIF


!-------------------------------------------------------------------------------
!- Section 4: Output of the namelist variables
!-------------------------------------------------------------------------------

IF (my_world_id == 0) THEN
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(A28)') '0  NAMELIST obs_radctl'
  WRITE (nuspecif, '(A28)') '   -------------------'
  WRITE (nuspecif, '(A2)')  '  '

  WRITE (nuspecif, '(T7,A,T21,A,T39,A)')  'Variable', 'Value', 'Format'
  WRITE (nuspecif, '(T7,A,T21,A,T39,A)')  'input_path', input_path, 'A'
  
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A)')  'Input Files:'
  WRITE (nuspecif, '(A2)')  '  '
  DO i = 1,COUNT(satpp_files /= '')
    WRITE (nuspecif, '(T7,A,T21,A)') 'SATPP', satpp_files(i)
  ENDDO
  
  WRITE (nuspecif, '(A2)')  '  '
  WRITE (nuspecif, '(T7,A)')  'TOVS_OBS_CHAN'
  WRITE (nuspecif, '(A2)')  '  '

  WRITE (nuspecif, '(T7,A,T21,A,T39,A)')  'Variable', 'Value', 'Format'
  DO i = 1,iradset
    WRITE (nuspecif, '(T7,A,T21,A,T39,A)')  'c_satellite', &
                                            sat_id2name(rad_set(i)%satid), &
                                            'A'
  ENDDO
ENDIF
#endif

ierrstat = 0

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE input_obs_radctl
!==============================================================================

!==============================================================================
!+ Procedure in "src_obs_sat" for the input of satpp files
!------------------------------------------------------------------------------

SUBROUTINE input_obs_satpp (ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine organizes the input of the satpp files
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable
    
! Local variables
  INTEGER   (KIND=iintegers) :: &
    nsatpp                                      ! number of satpp files
#ifdef RTTOV10    
  CHARACTER(LEN=256) :: &
    ymsg,         & ! error message
    path            ! path to input file

  CHARACTER(LEN=64) :: yversion         
  CHARACTER(LEN=45) :: yfdbk_descript
  CHARACTER(LEN=12) :: yveri_ref_datim

  INTEGER   (KIND=iintegers) :: &
    istat,                                    & ! status variable
    isatpp,                                   & ! couting variable for satpp files
    nrec,                                     & ! number of satpp records to read in
    irec,                                     & ! counting variable for record
    inobs,                                    & ! counting variable for obs
    ioffs,                                    & ! 
    ichan, imissing,                          & ! counting variable for channel
    iinstr,                                   & ! counting variable for channel
    icnt, istart,iend,                        & ! counting variable
    idiffm,                                   & ! time difference [minutes]
    imoyy, imomm, imodd, imohh,               & ! Start date of forecast               
    isthght,                                  & ! Station Height
    inobs_node_loc(num_compute),              & ! number of obs per node (local)
    inobs_dropped_area,                       & ! no of obs dropped for area 
    inobs_dropped_tstart,                     & ! no of obs dropped for too early 
    inobs_dropped_tend,                       & ! no of obs dropped for too late 
    iveri_ref_date,                           & ! verification reference date
    iveri_ref_time,                           & ! verification reference date
    iveri_start,                              & ! verification reference date
    iveri_end,                                & ! verification reference date
    domain(3),                                & ! 
    iveri_run_type,                           & ! 
    iveri_run_class,                          & ! 
    iveri_exp_id,                             & ! 
    iveri_forecast_time,                      & ! 
    iveri_epsmem,                             & ! 
    varid                                       ! 
  INTEGER (KIND=iintegers), ALLOCATABLE ::&  
    ifcaststep(:),                            & ! timesteps assigned to obs
    iob_tot(:),                               & ! zonal idx of gridpt 
    job_tot(:),                               & ! meridional idx of gridpt
    iindex(:),                                & ! index vector of rads to use
    isort_index(:),                           & ! index vector of rads to use
    inode_index(:)                              ! index vector of rads to use
  INTEGER (KIND=iintegers), POINTER ::&  
    idummy(:)                                   ! Array of channels
  REAL      (KIND=ireals)     ::              &
    rlon,rlat,                                & ! Rotated Latitude/Longitude
    zio_tot,zjo_tot,                          & ! obs location in grid points
    rdummy                                      ! dummy variable
  REAL :: & 
    pole(2),lower_left(2),upper_right(2),     & ! 
    resolution(2)                               ! 
  LOGICAL :: &
    valid,                                    & ! result of satpp reading routine
    lsea,                                     & ! If the obs is above sea
    ldummy
  TYPE (t_radv), TARGET  :: &
    radv            ! Container for radiances from satpp
  TYPE (t_rad_container) :: &
    obs_rad_loc     ! Container for radiances for COSMO
  TYPE (t_rad_set), POINTER  :: &
    rs              
  TYPE (t_radv), POINTER :: &
    r               ! Container for radiances from satpp
  TYPE (t_sensor), POINTER  :: &
    sensor
#endif   
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE input_obs_satpp
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Initialize variables 
!-------------------------------------------------------------------------------
  nsatpp   = COUNT(satpp_files /= '')
  nsensors = 0

#ifdef RTTOV10
  ALLOCATE (obs_rad(nsatpp),                &
            nobspe_gath(num_compute,nsatpp),&
            nobspe_tot(num_compute,nsatpp), &
            stat = istat)

  IF (my_world_id == 0 .AND. istat == 0) THEN
    ALLOCATE(fdbk(nsatpp),        & 
             fdbk_offset(nsatpp), &
              stat = istat)
  ENDIF

  IF (istat /= 0) THEN
    ierrstat = 99101
    PRINT *, ' ERROR   *** Allocating memory failed. ***'
    RETURN
  ENDIF  
 
  READ (ydate_ini,'(I4,I2,I2,I2)') imoyy, imomm, imodd, imohh
 

  IF (my_world_id == 0) THEN
    OPEN(nurad, FILE=yurad, FORM='FORMATTED', STATUS='UNKNOWN', &
                POSITION='APPEND', IOSTAT=istat)

    IF (istat /= 0) THEN
      ierrstat = 99101
      PRINT *, ' ERROR   *** OPENING OF FILE yurad FAILED. ***'
      RETURN
    ENDIF  
    
    WRITE (nurad,*) '  '
    WRITE (nurad,'(T2,A)') 'SATPP Input'
    WRITE (nurad,'(T2,A)') '-----------'
  ENDIF  
  
  DO isatpp=1,nsatpp
!-------------------------------------------------------------------------------
!- Section 2: Read in the satpp file
!-------------------------------------------------------------------------------
    path = TRIM(input_path) // TRIM(satpp_files(isatpp))

    ! find out about the size of the satpp
    call read_satpp_feedbk(TRIM(path),radv,valid,lprint=.False., lread=.False.)

    IF (.NOT. valid) THEN
      ierrstat = 99102
      PRINT *, ' ERROR   *** Inquiring length of satpp failed. ***'
      EXIT
    ENDIF

    ! The satpp file will be read in parallel on all processors
    ! Therefore each processor has a part of the file in its memory
    nrec   = (radv%n_rec + nproc - 1) / nproc 

    call read_satpp_feedbk(TRIM(path),radv,valid,            &
                           istart = nrec * my_world_id,         &
                           iend   = nrec * (my_world_id+1) - 1)

    IF (.NOT. valid) THEN
      ierrstat = 99103
      WRITE (stderr,*) ' ERROR   *** Reading satpp failed. ***'
      EXIT
    ENDIF

    nrec = radv%n_rec
    CALL global_values(nrec,1,'SUM',imp_integers,icomm_world,0,ymsg,istat)

    IF (my_world_id == 0) THEN
      WRITE (nurad,*) '  '
      WRITE (nurad,'(T4,A,T14,A)') 'Radiance File', radv % filename
      WRITE (nurad,'(T6,A,T30,I15)') 'Number Of Records: ', nrec
    ENDIF  

!-------------------------------------------------------------------------------
!- Section 3: Link Radiance data to namelist information 
!-------------------------------------------------------------------------------
    CALL link_rad(radv,istat)

    IF (istat /= 0) THEN
      ierrstat = 99103
      WRITE (stderr,*) &
        ' ERROR   *** Associating namelist input with satpp failed ***'
      EXIT
    ENDIF
  
    rs => radv%i
!-------------------------------------------------------------------------------
!- Section 3: Prepare additional fields 
!-------------------------------------------------------------------------------
    ALLOCATE(radv%state(rs%n_chan,radv%n_rec), &
             radv%flags(rs%n_chan,radv%n_rec), &
!             radv%check(rs%n_chan,radv%n_rec), &
             stat = istat)

    IF (istat /= 0) THEN
      ierrstat = 99101
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      EXIT
    ENDIF  

    DO irec=1,radv%n_rec
      DO ichan=1,rs%n_chan
        radv%state(ichan,irec) = ST_ACCEPTED
        radv%flags(ichan,irec) = 0

        ! set initial state
        IF (.NOT. radv%valid(ichan,irec)) THEN
          CALL change_use_rpt_rad(radv,ichan,irec,FL_OBSTYPE)
        ELSE IF (BTEST(rs%flag(ichan), USE_PASSIVE)) THEN
          CALL change_use_rpt_rad(radv,ichan,irec,FL_OBSTYPE)
        ENDIF
      ENDDO
    ENDDO  
    
!-------------------------------------------------------------------------------
!- Section 3: Evaluate obs location and time 
!-------------------------------------------------------------------------------
    ALLOCATE(ifcaststep(radv%n_rec), iob_tot(radv%n_rec),          &
             job_tot(radv%n_rec), iindex(radv%n_rec+1),            &
             isort_index(radv%n_rec+1), inode_index(radv%n_rec+1), &
             stat = istat)

    IF (istat /= 0) THEN
      ierrstat = 99101
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      EXIT
    ENDIF  

    ! The following is really weird, we should get rid of this
    ! slow function, somehow
    DO irec = 1,radv%n_rec
      CALL diff_minutes(imoyy,imomm,imodd,imohh,0,                 &
                        radv % date(irec)     / 10000     , &
                        MOD(radv % date(irec) / 100  ,100), &
                        MOD(radv % date(irec)        ,100), &
                        radv % time(irec)     / 10000     , &
                        MOD(radv % time(irec) / 100  ,100), &
                        itype_calendar, idiffm,             &
                        istat                                        )

      IF (istat /= 0) THEN
        ifcaststep(irec) = -1
      ELSE  
        ifcaststep(irec) = int(real(idiffm) * 60.0_ireals / dt)
      ENDIF
    ENDDO
    
    DO irec = 1,radv%n_rec
      rlon = rla2rlarot(radv%dlat(irec), radv%dlon(irec), &
                        pollat, pollon, polgam)

      rlat = phi2phirot(radv%dlat(irec), radv%dlon(irec), &
                        pollat, pollon)

      zio_tot = 1._ireals + (rlon - startlon_tot) / dlon
      zjo_tot = 1._ireals + (rlat - startlat_tot) / dlat
     
      ! the following was stolen from func obs_assign_gridpt
      rdummy = MAX(nolbc, 1+nboundlines) + 1.E-8_ireals

      IF (zio_tot > rdummy .AND. zio_tot < ie_tot + 1.0_ireals - rdummy .AND. &
          zjo_tot > rdummy .AND. zjo_tot < je_tot + 1.0_ireals - rdummy) THEN
          
        iob_tot(irec) = NINT(zio_tot)
        job_tot(irec) = NINT(zjo_tot)

        iob_tot(irec) = MAX(iob_tot(irec), 1 + nolbc)
        iob_tot(irec) = MIN(iob_tot(irec), ie_tot - nolbc)
        job_tot(irec) = MAX(job_tot(irec), 1 + nolbc)
        job_tot(irec) = MIN(job_tot(irec), je_tot - nolbc)
      ELSE
        iob_tot(irec) = 0
        job_tot(irec) = 0
      ENDIF  
    ENDDO  

!-------------------------------------------------------------------------------
!- Section 5: Filter by timestep and location 
!-------------------------------------------------------------------------------
    inobs_dropped_area   = 0
    inobs_dropped_tstart = 0
    inobs_dropped_tend   = 0
      
    inobs = 0
    DO irec = 1,radv%n_rec
      IF (iob_tot(irec) <= 0 .OR.  job_tot(irec) <= 0) THEN
        DO ichan = 1,rs%n_chan
          CALL change_use_rpt_rad(radv,ichan,irec,FL_AREA)
        ENDDO
        inobs_dropped_area = inobs_dropped_area + 1
      ELSE IF (ifcaststep(irec) < nstart) THEN
        DO ichan = 1,rs%n_chan
          CALL change_use_rpt_rad(radv,ichan,irec,FL_TIME)
        ENDDO
        inobs_dropped_tstart = inobs_dropped_tstart + 1
      ELSE IF (ifcaststep(irec) > nstop) THEN
        DO ichan = 1,rs%n_chan
          CALL change_use_rpt_rad(radv,ichan,irec,FL_TIME)
        ENDDO
        inobs_dropped_tend = inobs_dropped_tend + 1
      ENDIF

      IF (ANY(radv%state(:,irec) < ST_DISMISS)) THEN
        inobs = inobs + 1    
        iindex(inobs) = irec
      ENDIF    
    ENDDO

    CALL global_values(inobs_dropped_area,1,'SUM',imp_integers,&
                       icomm_world,0,ymsg,istat)
    CALL global_values(inobs_dropped_tstart,1,'SUM',imp_integers,&
                       icomm_world,0,ymsg,istat)
    CALL global_values(inobs_dropped_tend,1,'SUM',imp_integers,&
                       icomm_world,0,ymsg,istat)

    IF (my_world_id == 0) THEN
      WRITE (nurad,'(T6,A,T30,I15)') '# Obs in wrong area:', inobs_dropped_area
      WRITE (nurad,'(T6,A,T30,I15)') '# Obs too early:',     inobs_dropped_tstart
      WRITE (nurad,'(T6,A,T30,I15)') '# Obs too late:',      inobs_dropped_tend
    ENDIF  

!-------------------------------------------------------------------------------
!- Section 4: Assign selected obs to nodes 
!-------------------------------------------------------------------------------
    IF (inobs > 0) THEN
      CALL obs_assign_sort_node(SIZE(iob_tot), inobs, iindex(1:inobs+1), &
                                iob_tot,job_tot, num_compute,       &
                                isubpos, nboundlines, my_cart_id,       &
                                inode_index(1:inobs+1),isort_index(1:inobs+1))
      DO icnt=1,num_compute
        inobs_node_loc(icnt) = COUNT(inode_index(1:inobs) == (icnt - 1))
      ENDDO
    ELSE
      inobs_node_loc(:) = 0
    ENDIF  

    IF (SUM(inobs_node_loc) /= inobs) THEN
      ierrstat = 99102
      WRITE (stderr,*) ' ERROR   *** Inconsistent obs numbers. ***'
      EXIT
    ENDIF

    nobspe_tot(:,isatpp) = inobs_node_loc(:)
    
    CALL global_values(nobspe_tot(:,isatpp),SIZE(nobspe_tot,1),'SUM',imp_integers,&
                       icomm_world,0,ymsg,istat)
    IF (istat /= 0) CALL model_abort(my_cart_id,99901, 'global_values', &
                                     'input_obs_satpp', istat)

    CALL distribute_values(nobspe_tot(:,isatpp),SIZE(nobspe_tot,1), 0, imp_integers, &
                       icomm_world,istat)
    IF (istat /= 0) CALL model_abort(my_cart_id,99902, 'distribute_values', &
                                     'input_obs_satpp', istat)
    
    IF (my_world_id == 0) THEN
      WRITE (ymsg,*) '(T6,A,T30,',num_compute,'(1X,I4))'
      WRITE (nurad,ymsg) 'Node:', (/(icnt,icnt=0,num_compute-1)/)
      WRITE (nurad,ymsg) 'Assigned Obs:', nobspe_tot(:,isatpp)
    ENDIF  

!-------------------------------------------------------------------------------
!- Section 5: Create obs container with obs to be distributed from this node
!-            to other nodes (ordered by nodeid) 
!-------------------------------------------------------------------------------
    call construct(obs_rad_loc, inobs, istat, radv) 

    IF (istat /= 0) THEN
      ierrstat = 99102
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      EXIT
    ENDIF 

    obs_rad_loc%ntstep(:)  = ifcaststep(isort_index(1:inobs))
    obs_rad_loc%iob_tot(:) = iob_tot(isort_index(1:inobs))
    obs_rad_loc%job_tot(:) = job_tot(isort_index(1:inobs))

    IF (ASSOCIATED(radv%obsnum)) THEN
      obs_rad_loc%radv%obsnum(:) = radv%obsnum(isort_index(1:inobs))
    ENDIF
    
    obs_rad_loc%radv%date(:) = radv%date(isort_index(1:inobs))
    obs_rad_loc%radv%time(:) = radv%time(isort_index(1:inobs))

    IF (ASSOCIATED(radv%date_d)) THEN
      obs_rad_loc%radv%date_d = radv%date_d(isort_index(1:inobs))
    ENDIF
  
    IF (ASSOCIATED(radv%time_d)) THEN
      obs_rad_loc%radv%time_d = radv%time_d(isort_index(1:inobs))
    ENDIF  

    obs_rad_loc%radv%dlat(:) = radv%dlat(isort_index(1:inobs))
    obs_rad_loc%radv%dlon(:) = radv%dlon(isort_index(1:inobs))
    obs_rad_loc%radv%fov (:) = radv%fov (isort_index(1:inobs))

    IF (ASSOCIATED(radv%scanl)) THEN
      obs_rad_loc%radv%scanl (:) = radv%scanl (isort_index(1:inobs))
    ENDIF
  
    obs_rad_loc%radv%stzen (:) = radv%stzen (isort_index(1:inobs))
    obs_rad_loc%radv%stazi (:) = radv%stazi (isort_index(1:inobs))

    IF (ASSOCIATED(radv%sunzen)) THEN
      obs_rad_loc%radv%sunzen(:) = radv%sunzen(isort_index(1:inobs))
    ENDIF  

    IF (ASSOCIATED(radv%sunazi)) THEN
      obs_rad_loc%radv%sunazi(:) = radv%sunazi(isort_index(1:inobs))
    ENDIF  

    obs_rad_loc%radv%landfr(:,:) = radv%landfr(:,isort_index(1:inobs))
    obs_rad_loc%radv%stype(:,:)  = radv%stype(:,isort_index(1:inobs))
    obs_rad_loc%radv%shgt(:,:)   = radv%shgt(:,isort_index(1:inobs))

    IF (ASSOCIATED(radv%center)) THEN
      obs_rad_loc%radv%center(:) = radv%center(isort_index(1:inobs))
    ENDIF
  
    IF (ASSOCIATED(radv%subcenter)) THEN
      obs_rad_loc%radv%subcenter(:) = radv%subcenter(isort_index(1:inobs))
    ENDIF

    obs_rad_loc%radv%valid(:,:)  = radv%valid(:,isort_index(1:inobs))
    obs_rad_loc%radv%state(:,:)  = radv%state(:,isort_index(1:inobs))
    obs_rad_loc%radv%flags(:,:)  = radv%flags(:,isort_index(1:inobs))
    obs_rad_loc%radv%bt_obs(:,:) = radv%bt_obs(:,isort_index(1:inobs))

!-------------------------------------------------------------------------------
!- Section 6: Merge Obs from all nodes at assigned node 
!-------------------------------------------------------------------------------
    call construct(obs_rad(isatpp), nobspe_tot(my_cart_id+1,isatpp), istat, obs_rad_loc%radv)

    IF (istat /= 0) THEN
      ierrstat = 99102
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      EXIT
    ENDIF 
    
    ioffs = 1
    DO icnt = 1, num_compute
      call gather_obs_rad_container(                              &
        obs_rad_loc,     ioffs, ioffs + inobs_node_loc(icnt) - 1, &
        obs_rad(isatpp),     1,                           irec, &
        icnt-1,istat, nobspe_gath(:,isatpp)) 
      ioffs = ioffs + inobs_node_loc(icnt)
    ENDDO

!-------------------------------------------------------------------------------
!- Section 7: Clean up
!-------------------------------------------------------------------------------
!CDIR NOIEXPAND
    call destruct(obs_rad_loc)
    
!CDIR NOIEXPAND
    call destruct(radv)

    DEALLOCATE(ifcaststep, iob_tot, job_tot, iindex, isort_index, inode_index)
   
!-------------------------------------------------------------------------------
!- Section 8: Allocate  & Prepare Fields fg calculations
!-------------------------------------------------------------------------------
    r     => obs_rad(isatpp)%radv
    rs    => r%i
    inobs =  r%n_rec

    ALLOCATE(r%mdlsfc (inobs),                        &
             r%bt_fg  (rs%n_chan,inobs),              &
             r%bt_bcor(rs%n_chan,inobs),              &
             r%bcor_  (rs%n_chan,inobs),              &
             r%p      (inobs),                        &
             r%h_fg   (rs%n_instr,inobs),             &
             r%t_fg   (ke + 1,inobs),                 &
             r%q_fg   (ke + 1,inobs),                 &
             r%t2m    (inobs),                        &
             r%q2m    (inobs),                        &
             r%ps_fg  (inobs),                        &
             r%ts_fg  (inobs),                        &
             r%u10_fg (inobs),                        &
             r%v10_fg (inobs),                        &
             obs_rad(isatpp)%plevel(rs%n_chan,inobs), &
             stat = istat)
    
    IF (istat /= 0) THEN
      ierrstat = 99102
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      EXIT
    ENDIF 

    r%n_lev = ke + 1

    r%mdlsfc(:)    = 0
    r%bt_fg(:,:)   = rmdi
    r%bt_bcor(:,:) = rmdi
    r%bcor_(:,:)   = 0._ireals
    r%p(:)         = rmdi
    r%h_fg(:,:)    = rmdi
    r%t_fg(:,:)    = rmdi
    r%q_fg(:,:)    = rmdi
    r%t2m(:)       = rmdi
    r%q2m(:)       = rmdi
    r%ps_fg(:)     = rmdi
    r%ts_fg(:)     = rmdi
    r%u10_fg(:)    = rmdi
    r%v10_fg(:)    = rmdi

    obs_rad(isatpp)%plevel(:,:) = rmdi

    IF (.NOT. ASSOCIATED(r%sunzen)) THEN
      ALLOCATE(r%sunzen(inobs),stat=istat)

      IF (istat /= 0) THEN
        ierrstat = 99102
        PRINT *, ' ERROR   *** Allocating memory failed. ***'
        EXIT
      ENDIF 

      r%sunzen(:) = rmdi
    END IF  
  ENDDO
  
!-------------------------------------------------------------------------------
!- Section 9: Compute RTTOV instrument information
!-------------------------------------------------------------------------------
  ! This may reserve more space than necessary if same instruments
  ! are input. Anyway, I dont see much benefit from saving some bytes
  ! to double looping over all the stuff again
  ALLOCATE (sensors(SUM(obs_rad(:)%radv%i%n_instr)), stat=istat)

  IF (istat /= 0) THEN
    ierrstat = 99102
    PRINT *, ' ERROR   *** Allocating memory failed. ***'
    RETURN
  ENDIF 
  
  !-------------------------------------------------------------------------------
  !+ Section 9.1: Merge channel information for all input files 
  !-------------------------------------------------------------------------------
  DO isatpp=1,nsatpp
    r  => obs_rad(isatpp)%radv
    rs => r%i
   
    DO iinstr = 1, rs % n_instr
      sensors(nsensors+1)%id = &
                   (/ rs%platform   , &
                      rs%rttov_satid, &
                      rs%instr(iinstr)/)
                   
      ! Search already defined sensors for sensor
      DO icnt=1,nsensors
        IF (ALL(sensors(icnt)%id(:) == sensors(nsensors+1)%id(:))) THEN
          rs%rttov_indx(iinstr) = icnt
          EXIT
        ENDIF
      ENDDO

      ! sensor not yet in use by any other radset, store it
      IF (rs%rttov_indx(iinstr) == -1) THEN
        nsensors = nsensors + 1
        rs%rttov_indx(iinstr) = nsensors
      ENDIF
      
      sensor => sensors(rs%rttov_indx(iinstr))

      ! check if the sensor entry misses some of the channels in this rad set
      IF (ASSOCIATED(sensor%channels)) THEN
        imissing = 0
        DO icnt=rs%o_ch_i(iinstr) + 1,rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)
          IF (.NOT. ANY(sensor%channels == rs%chan(icnt))) THEN
            imissing = imissing+1
          ENDIF
        ENDDO
      ELSE
        imissing = rs%n_ch_i(iinstr)
      ENDIF  

      IF (imissing > 0) THEN
        istart = MINVAL(rs%chan(rs%o_ch_i(iinstr) + 1:&
                                rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)))
        iend   = MAXVAL(rs%chan(rs%o_ch_i(iinstr) + 1:&
                                rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)))

        IF (ASSOCIATED(sensor%channels)) THEN
          ALLOCATE(idummy(imissing + SIZE(sensor%channels)))
          istart = MIN(istart,MINVAL(sensor%channels))
          iend   = MAX(iend,  MAXVAL(sensor%channels))
        ELSE
          ALLOCATE(idummy(imissing))
        ENDIF  

        ! merge existing sensor channels with radset's channels 
        ! sensors are automatically sorted  
        icnt = 0
        DO ichan=istart,iend
          IF (ASSOCIATED(sensor%channels)) THEN
            IF (ANY(ichan == sensor%channels)) THEN
              icnt = icnt+1
              idummy(icnt)= ichan
              CYCLE
            ENDIF
          ENDIF  
            
          IF (ANY(ichan == rs%chan(rs%o_ch_i(iinstr) + 1:                      &
                                   rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)))) THEN
              icnt = icnt+1
              idummy(icnt)= ichan
              CYCLE
          ENDIF    
        ENDDO

        ! interchange arrays
        IF (ASSOCIATED(sensor%channels)) DEALLOCATE(sensor%channels)
        sensor%channels => idummy
      ENDIF
    ENDDO
  ENDDO
 
  !-------------------------------------------------------------------------------
  !+ Section 9.2: map input file channels to rttov channel indices
  !-------------------------------------------------------------------------------
  ! The following is needed, as we only load these channels in rttov we
  ! really need. Therefore the channel indices of rttov will change e.g.
  ! input data channels 4-7 -> rttov channels 1-3
  DO isatpp=1,nsatpp
    r  => obs_rad(isatpp)%radv
    rs => r%i
    
    ALLOCATE(obs_rad(isatpp)%rttov_chan(rs%n_chan), stat = istat)
    
    IF (istat /= 0) THEN
      ierrstat = 99102
      PRINT *, ' ERROR   *** Allocating memory failed. ***'
      RETURN
    ENDIF 
    
    DO iinstr = 1, rs % n_instr
      sensor => sensors(rs%rttov_indx(iinstr))
     
      ! TODO: optimize the following, we know that sensor%channels is
      !       an ordered array
      DO icnt=rs%o_ch_i(iinstr) + 1,rs%o_ch_i(iinstr) + rs%n_ch_i(iinstr)
        obs_rad(isatpp)%rttov_chan(icnt) = &
          COUNT(sensor%channels(:) <= rs%chan(icnt))
      ENDDO
    ENDDO  
  ENDDO

  WHERE (sensors(:)%id(3) == 0) ! hirs has clouds
    sensors(:)%addcloud = .True.
  ELSEWHERE
    sensors(:)%addcloud = .False.
  END WHERE  
   
!-------------------------------------------------------------------------------
!- Section 10: Prepare some globally used variables
!-------------------------------------------------------------------------------
  READ (ydate_ini,'(I4,I2,I2,I2)') reftime%year, reftime%month, reftime%day, &
                                   reftime%hour
  
  reftime%minute = 0
!-------------------------------------------------------------------------------
!- Section 10: Prepare fdbk file output (pe 0 only
!-------------------------------------------------------------------------------
  IF (my_world_id == 0) THEN
    fdbk_offset(:) = 0
    yversion = yncglob_source 
    IF (yversion(1:5) == 'COSMO') yversion = yversion(7:LEN_TRIM(yversion)-6)

    iveri_ref_date = reftime%year*1000 + reftime%month*100 + reftime%day
    iveri_ref_time = reftime%hour * 100 + reftime%minute 
    
    iveri_start = hstart * 60.0_ireals
    iveri_end   = hstop  * 60.0_ireals

    resolution  = (/dlat, dlon/)
    domain      = (/ie_tot, je_tot, ke_tot/)
    pole        = (/startlat_tot, startlon_tot/)
    lower_left  = pole
    upper_right = pole + resolution * (/je_tot-1,ie_tot-1/)

    WRITE (yveri_ref_datim,'(I8.8,I4.4)') iveri_ref_date, iveri_ref_time

    IF (nvers >= 0 .AND. nvers < 1048576) THEN
      iveri_exp_id    = MOD(nvers,16384)
      iveri_run_class = nvers/16384
    ELSE
      iveri_exp_id    = imdi
      iveri_run_class = imdi
    ENDIF

    IF (leps) THEN
      WRITE (yfdbk_descript,'(A,1X,I3)') 'ensemble forecast member', &
                                        iepsmem
      iveri_epsmem = iepsmem                                  
    ELSE
      yfdbk_descript = 'deterministic forecast'
      iveri_epsmem   = -1                                  

      IF ( iveri_run_class == rc_ass .OR. iveri_end <= 300) THEN
        iveri_run_type = vt_firstguess
      ELSE
        iveri_run_type = vt_analysis
      ENDIF  
    ENDIF  

    iveri_forecast_time = INT(hstop * 100)
    
    DO isatpp=1,nsatpp
      rs   => obs_rad(isatpp)%radv%i
      inobs =  SUM(nobspe_tot(:,isatpp))

      IF (inobs <= 0) CYCLE

      CALL setup_fdbk(fdbk(isatpp)%nc, latex = .FALSE.)

      DO iinstr=1,rs%n_instr
        WRITE (ymsg((iinstr-1)*3 + 1:),'(I2.2,A1)') rs%instr(iinstr), '-'
      ENDDO
      
      WRITE (path,'(A,A,I3.3,A,I2.2,A,A,A)')          &
        TRIM(root%ydir),'/monRAD', rs%satid,'_instr', rs%grid, &
        '_', ymsg(1:rs%n_instr*3-1),'.nc'

      CALL create_fdbk(fdbk(isatpp),                 &
                       TRIM(path),                   &
                       'COSMO',                      &
                       yversion, yncglob_institution, &
                       inobs, inobs * rs%n_chan,     &
                       iveri_ref_date,               &
                       iveri_ref_time,               &
                       iveri_start,                  &
                       iveri_end,                    &
                       resolution,                   &
                       domain,                       &
                       yfdbk_descript,               &
                       yveri_ref_datim,              &
                       pole        = pole,           &
                       lower_left  = lower_left,     &
                       upper_right = upper_right,    &
                       opt         ='RAD'          )

      CALL add_verification(fdbk(isatpp),            &
                            'COSMO',                 &
                            iveri_run_type,          &
                            iveri_run_class,         &
                            yveri_ref_datim,         &
                            iveri_forecast_time,     &
                            resolution,              &
                            domain,                  &
                            'first guess',           &
                            iveri_epsmem,            &
                            iveri_exp_id,            &
                            varid                    )
    ENDDO
  ENDIF
!-------------------------------------------------------------------------------
!- Section 11: Writeout some nice information
!-------------------------------------------------------------------------------
  IF (my_world_id == 0) THEN
    WRITE (nurad,*) '  '
    WRITE (nurad,'(T4,A,T20,A)') 'Sensor Triplet' , 'Channels'
    WRITE (nurad,*) '  '

    DO icnt=1,nsensors
      ymsg = ''

      DO ichan = 1,SIZE(sensors(icnt)%channels)
        IF (LEN_TRIM(ymsg) > 50) THEN
          WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) '...'
          EXIT
        ENDIF
        
        IF (ichan == 1) THEN
          WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) sensors(icnt)%channels(ichan)
        ELSEIF (ichan == SIZE(sensors(icnt)%channels)) THEN
          IF (sensors(icnt)%channels(ichan) - 1 == &
              sensors(icnt)%channels(ichan-1)      ) THEN
            WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) &
                   '-',sensors(icnt)%channels(ichan)
          ENDIF         
        ELSE
          IF (sensors(icnt)%channels(ichan) - 1 /= &
              sensors(icnt)%channels(ichan-1)      ) THEN
             
            IF (sensors(icnt)%channels(ichan-1) - 1 == &
                sensors(icnt)%channels(ichan-2)        ) THEN
              WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) &
                     '-',sensors(icnt)%channels(ichan-1)
            ENDIF         
          
            WRITE (ymsg(LEN_TRIM(ymsg)+1:),*) &
                  ', ', sensors(icnt)%channels(ichan)
          ENDIF
        ENDIF
      ENDDO

      WRITE (nurad,'(T4,3(1x,I2),T20,A)') sensors(icnt)%id(:), TRIM(ymsg)
    ENDDO
  ENDIF
!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
#endif
  ierrstat = 0

END SUBROUTINE input_obs_satpp

!==============================================================================
!+ Procedure in "src_obs_sat" for calculating the radiances
!------------------------------------------------------------------------------

SUBROUTINE calc_obs_satpp (ierrstat)

!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the radiances
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (INOUT)   ::        &
    ierrstat        ! error status variable

! Local variables
#ifdef RTTOV10
  INTEGER   (KIND=iintegers) :: &
    istat,                       & ! Status variable
    inobs,                       & ! Number of observables to calc fg for
    iradv,                       & ! Counter for radv (satpp input file)
    iobs,                        & ! Counter for obs
    ichan,                       & ! Counter for channels
    dummy,                       & ! surface type
    isurfsens,                   & ! grid instrument index
    icnt,iinstr,                 & ! Counting variables
    ipstart,ipend,               & ! Array slice offsets
    icstart,icend,               & ! Array slice offsets
    ncalcs,                      & ! number of calculations for rttov
    nlev                           ! number of levels for rttov calcs
  INTEGER   (KIND=iintegers), ALLOCATABLE :: &
    iob_loc(:),               & ! array of zonal gridpoint indices (local)
    job_loc(:),               & ! array of meridional gridpoint indices (local)
    profs(:),                 & ! array of profile indizes
    chans(:),                 & ! array of channel indizes
    idg(:),                   &
    ish(:),                   &
    wtype(:),                 &
    stype(:)
  REAL      (KIND=ireals), ALLOCATABLE :: &
    emiss(:),                 &
    emiss_k(:),               &
    temp_k(:,:,:),            &
    humi_k(:,:,:),            &
    t2m_k (:,:),              &
    q2m_k (:,:),              &
    stemp_k (:,:),            &
    pres(:,:),                &
    cfrac(:,:,:),             &
    cloud(:,:,:)
  REAL      , ALLOCATABLE, TARGET :: &
    report_veridata(:,:,:)     
  REAL      (KIND=ireals), POINTER :: &
    h_ice_dummy(:,:,:),       & !
    qg_dummy   (:,:,:)          !
  CHARACTER(LEN=256) :: &
    ymsg            ! error message
  TYPE(t_radv), POINTER :: &
    r  
  TYPE(t_rad_set), POINTER :: &
    rs  
  TYPE(t_sensor), POINTER :: &
    sensor  
  TYPE(t_rad_container), TARGET :: &
    obs
  TYPE(t_account), ALLOCATABLE :: &
    reports(:)
  TYPE(t_acc_body), ALLOCATABLE, TARGET :: &
    report_body(:,:)
#endif    

! Tracer pointers
REAL (KIND=ireals), POINTER :: &
  qv  (:,:,:) => NULL(),    &   ! QV at nnow
  qc  (:,:,:) => NULL(),    &   ! QC at nnow
  qg  (:,:,:) => NULL(),    &   ! QG at nnow
  qi  (:,:,:) => NULL(),    &   ! QI at nnow
  qs  (:,:,:) => NULL()         ! QS at nnow

CHARACTER (LEN=50)  :: yerror       ! error message
CHARACTER (LEN=25)  :: yzroutine

!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE calc_obs_satpp
!-------------------------------------------------------------------------------

  ierrstat = 0
  yerror   = '   '
  yzroutine= 'calc_obs_satpp'

  ! retrieve the required microphysics tracers
  CALL trcr_get(ierrstat, idt_qv, ptr_tlev = nnow, ptr=qv)
  IF (ierrstat /= 0) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF
  CALL trcr_get(ierrstat, idt_qc, ptr_tlev=nnow, ptr=qc)
  IF (ierrstat /= 0) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF
  CALL trcr_get(ierrstat, idt_qg, ptr_tlev=nnow, ptr=qg)
  IF (ierrstat /= 0 .AND. ierrstat /= T_ERR_NOTFOUND) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF
  CALL trcr_get(ierrstat, idt_qi, ptr_tlev=nnow, ptr=qi)
  IF (ierrstat /= 0 .AND. ierrstat /= T_ERR_NOTFOUND) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF
  CALL trcr_get(ierrstat, idt_qs, ptr_tlev=nnow, ptr=qs)
  IF (ierrstat /= 0 .AND. ierrstat /= T_ERR_NOTFOUND) THEN
    yerror = trcr_errorstr(ierrstat)
    CALL model_abort(my_cart_id, ierrstat, yerror, yzroutine)
  ENDIF

#ifdef RTTOV10
  IF (ALLOCATED(h_ice)) THEN
    h_ice_dummy => h_ice(:,:,:)
  ELSE
    ALLOCATE(h_ice_dummy(0,0,nnow))
  ENDIF
  IF (ASSOCIATED(qg)) THEN
    qg_dummy => qg(:,:,:)
  ELSE
    ALLOCATE(qg_dummy(0,0,0))
  ENDIF
 
  DO iradv=1,SIZE(obs_rad)
    r => obs_rad(iradv)%radv
    rs => r%i
    
!-------------------------------------------------------------------------------
!- Section 1: Prepare obs to be computede
!-------------------------------------------------------------------------------
    ipstart = COUNT((obs_rad(iradv) % ntstep(:) - ntstep) < 0 ) + 1
    ipend   = COUNT((obs_rad(iradv) % ntstep(:) - ntstep) <= 0)
    
    inobs = ipend - ipstart + 1
    nlev = r%n_lev

    ! find the index of the grid instrument
    isurfsens = 0
    
    IF (rs%flag_instr >= 0) THEN
      DO icnt=1,rs%n_instr
        IF (rs%flag_instr == rs%instr(icnt)) isurfsens = icnt
      END DO  
    ELSE
      DO icnt=1,rs%n_instr
        IF (rs%grid == rs%instr(icnt)) isurfsens = icnt
      END DO  
    ENDIF  

    IF (isurfsens <= 0) THEN
        ierrstat = 99107
        PRINT *, ' ERROR   *** invalid flag_instr in namelist. ***'
        EXIT
    ENDIF    
      
     
    IF (inobs > 0) THEN
      IF (ANY(obs_rad(iradv)%ntstep(ipstart:ipend) /= ntstep)) THEN
        ierrstat = 99104
        PRINT *, ' ERROR   *** bad order of profiles. ***'
        EXIT
      ENDIF


      ALLOCATE(iob_loc(inobs),    &
               job_loc(inobs),    &
               pres(nlev,ipstart:ipend),         &
               wtype(inobs),                     &
               stype(inobs),                     &
               cloud(6,nlev,inobs),              &
               cfrac(6,nlev,inobs),              &
               idg(inobs),                       &
               ish(inobs),                       &
               stat = istat)

      IF (istat /= 0) THEN
          ierrstat = 99102
          PRINT *, ' ERROR   *** Allocating memory failed. ***'
          EXIT
      ENDIF    

      DO iobs=ipstart,ipend
        icnt = iobs - ipstart + 1
        
        iob_loc(icnt) = obs_rad(iradv)%iob_tot(iobs) - &
                        isubpos(my_cart_id,1) + 1 + nboundlines
        job_loc(icnt) = obs_rad(iradv)%job_tot(iobs) - &
                        isubpos(my_cart_id,2) + 1 + nboundlines

        ! fall back to comso values if no satpp value available
        IF (r%sunzen(iobs) < rmdich) THEN
          r%sunzen(iobs) = sun_el(iob_loc(icnt),job_loc(icnt))
        ENDIF  
      ENDDO

!-------------------------------------------------------------------------------
!- Section 2: Prepare RTTOV profile data from model statee
!-------------------------------------------------------------------------------


      CALL prepare_rttov_input(iob_loc, job_loc, extrp_const, nlev,        &
                               pp(:,:,:,nnow),ps(:,:,nnow),                &
                               t(:,:,:,nnow),t_g(:,:,nnow),                &
                               qv(:,:,:),qc(:,:,:),                        &
                               qi(:,:,:),qs(:,:,:),                        &
                               qg_dummy(:,:,:),h_ice_dummy(:,:,nnow),      &
                               pres(:,ipstart:ipend),                      &
                               r%t_fg(:,ipstart:ipend),                    &
                               r%q_fg(:,ipstart:ipend),                    &
                               r%t2m(ipstart:ipend),                       &
                               r%q2m(ipstart:ipend),                       &
                               r%ps_fg(ipstart:ipend),                     &
                               r%h_fg(1,ipstart:ipend),                    &
                               r%u10_fg(ipstart:ipend),                    &
                               r%v10_fg(ipstart:ipend),                    &
                               r%ts_fg(ipstart:ipend),                     &
                               stype(:),                                   &
                               wtype,                                      &
                               cloud, cfrac,idg,ish,                       &
                               ierrstat = istat) 

      IF (istat /= 0) THEN
          ierrstat = 99104
          PRINT *, ' ERROR   *** prepare rttov_input failed. ***'
          EXIT
      ENDIF    
 
      ! short infor on formats:
      ! surftype satpp: sea,mixed,land (0,1,2)
      ! surftype rttov: land,sea,sea-ice (0,1,2)
      ! we use the satpp encoding for radv instances here!!!
      
      IF (usesurftype == 1) THEN
        DO iinstr=1,rs%n_instr
          IF (rs%flag_instr == -2) THEN
            r%h_fg(iinstr,ipstart:ipend) = MAXVAL(r%shgt(:,ipstart:ipend),1)
          ELSE
            r%h_fg(iinstr,ipstart:ipend) = r%shgt(isurfsens,ipstart:ipend)
          ENDIF
        END DO                   

        DO iobs=ipstart,ipend
          IF (rs%flag_instr == -2) THEN
            dummy = MAXVAL(r%stype(:,iobs))
          ELSE
            dummy = r%stype(isurfsens,iobs)
          ENDIF 

          IF (dummy == 2) THEN
           ! if satpp says land(2) force land only
            r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_LAND)
            stype(iobs) = 0
          ELSE 
            IF (dummy == 0 .OR. dummy == 1) THEN
              r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_SEA)
            ENDIF
            
            IF (dummy == 1) r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_LAND)

            IF (stype(iobs) == 2) THEN
              r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_ICE)
            ELSE  
              r%mdlsfc(iobs) = IBSET(r%mdlsfc(iobs),MS_NO_ICE)
            END IF
            ! satpp says land/mixed, what now?
          END IF  
          
          IF (rs%flag_instr == -3) THEN
            IF (ANY(r%stype(:,iobs) /= r%stype(1,iobs))) THEN
              DO ichan = 1, rs%n_chan
                CALL change_use_rpt_rad(r,ichan,iobs,FL_SURF)
              ENDDO  
            ENDIF
          ENDIF
        END DO  
      ELSE
        ! use surface from cosmo
        ! prepare_rttov_input generates kilometer heights as
        ! required by rttov, but h_fg in radv is in meters
        r%h_fg(1,ipstart:ipend) = r%h_fg(1,ipstart:ipend) * 1000.0

!CDIR NOLOOPCHG      
        DO iinstr=2,rs%n_instr
          r%h_fg(iinstr,ipstart:ipend) = r%h_fg(1,ipstart:ipend)
        END DO  

        WHERE (stype(:) == 0)
          r%mdlsfc(:) = IBSET(r%mdlsfc(:),MS_LAND)
          r%stype(:,iobs) = 2
        ELSEWHERE (stype(:) == 1)
          r%mdlsfc(:) = IBSET(r%mdlsfc(:),MS_SEA)
          r%mdlsfc(:) = IBSET(r%mdlsfc(:),MS_NO_ICE)
          r%stype(:,iobs) = 0
        ELSEWHERE (stype(:) == 2)
          r%mdlsfc(:) = IBSET(r%mdlsfc(:),MS_SEA)
          r%mdlsfc(:) = IBSET(r%mdlsfc(:),MS_ICE)
          r%stype(:,iobs) = 0
        END WHERE  
      END IF
      

      istat = rttov_fill_input(                     &
        pres(1:nlev,ipstart:ipend),                 &
        r%t_fg(1:nlev,ipstart:ipend),               &
        r%q_fg(1:nlev,ipstart:ipend),               &
        r%t2m(ipstart:ipend),                       &
        r%q2m(ipstart:ipend),                       &
        r%ps_fg(ipstart:ipend),                     &
        r%h_fg(1,ipstart:ipend) * 0.001,            &
        r%u10_fg(ipstart:ipend),                    &
        r%v10_fg(ipstart:ipend),                    &
        r%ts_fg(ipstart:ipend),                     &
        stype(:),                                   &
        r%dlat (ipstart:ipend),                     &
        r%stzen(ipstart:ipend),                     &
        r%sunzen(ipstart:ipend),                    &
        satazim       = r%stazi(ipstart:ipend),     &
        cloud         = cloud(:,1:nlev-1,:),        &
        cfrac         = cfrac(:,1:nlev-1,:),        &
        idg           = idg(:),                     &
        ish           = ish(:),                     &
        watertype     = wtype(:),                   & 
        addsolar      = .FALSE.,                    &
        addrefrac     = .TRUE.,                     &
        addinterp     = .TRUE.,                     &
        ivect         = 1,                          &
        rttov9_compat = .FALSE.)
            
      DEALLOCATE(wtype,stype,cloud,cfrac,idg,ish,iob_loc,job_loc)

      IF (istat /= 0) THEN
          ierrstat = 99104
          PRINT *, ' ERROR   *** rttov_fill_input failed. ***'
          EXIT
      ENDIF    
      
!-------------------------------------------------------------------------------
!- Section 3: Compute First guess
!-------------------------------------------------------------------------------
      ALLOCATE(profs(inobs * MAXVAL(rs%n_ch_i)),              &
                chans(inobs * MAXVAL(rs%n_ch_i)),             &
                emiss(inobs * MAXVAL(rs%n_ch_i)),             &
                emiss_k(inobs * MAXVAL(rs%n_ch_i)),           &
                temp_k(nlev,MAXVAL(rs%n_ch_i),ipstart:ipend), &
                humi_k(nlev,MAXVAL(rs%n_ch_i),ipstart:ipend), &
                t2m_k (MAXVAL(rs%n_ch_i),inobs),              &
                q2m_k (MAXVAL(rs%n_ch_i),inobs),              &
                stemp_k (MAXVAL(rs%n_ch_i),inobs),            &
                stat = istat)

      IF (istat /= 0) THEN
          ierrstat = 99102
          PRINT *, ' ERROR   *** Allocating memory failed. ***'
          EXIT
      ENDIF    

      DO iinstr = 1, rs%n_instr
        sensor => sensors(rs%rttov_indx(iinstr))
            
        DO iobs = 1, inobs
          icstart = 1 + (iobs - 1) * rs%n_ch_i(iinstr)
          icend   =     (iobs)     * rs%n_ch_i(iinstr)
          profs(icstart:icend) = iobs
          chans(icstart:icend) = obs_rad(iradv)%rttov_chan(&
                                       1+rs%o_ch_i(iinstr):rs%n_ch_i(iinstr)+rs%o_ch_i(iinstr))
        ENDDO

        ncalcs = inobs * rs%n_ch_i(iinstr)

        icstart = 1                 + rs%o_ch_i(iinstr)
        icend   = rs%n_ch_i(iinstr) + rs%o_ch_i(iinstr)
            
        emiss(1:ncalcs)   = 0.0_ireals
        emiss_k(1:ncalcs) = 0.0_ireals

        istat = rttov_set_input(                            &
          hsurf=r%h_fg(iinstr,ipstart:ipend)*0.001          )
        
        ! We use the k matrix rttov interface, as we have to calculate
        ! the most sensitive height as well (for LETKF)
        istat = rttov_k_ifc(                            &
                  rs%rttov_indx(iinstr) + num_sensors,         &
                  profs(1:ncalcs),                             &
                  chans(1:ncalcs),                             &
                  emiss(1:ncalcs),                             &
                  emiss_k(1:ncalcs),                           &
                  temp_k (:,1:rs%n_ch_i(iinstr),ipstart:ipend),&
                  humi_k (:,1:rs%n_ch_i(iinstr),ipstart:ipend),&
                  t2m_k  (  1:rs%n_ch_i(iinstr),:),            &
                  q2m_k  (  1:rs%n_ch_i(iinstr),:),            &
                  stemp_k(  1:rs%n_ch_i(iinstr),:),            &
                  t_b = r%bt_fg(icstart:icend, ipstart:ipend))



        IF (istat == 0) THEN
          DO iobs=ipstart,ipend
            SELECT CASE(fdbk_height_method)
              CASE (1)
                obs_rad(iradv)%plevel(icstart:icend,iobs) = 100.0_ireals * &
                  pres(MAXLOC(ABS(temp_k(:,1:rs%n_ch_i(iinstr),iobs)),1),iobs)
              CASE (2)
                IF (ANY( SUM(temp_k(:,1:rs%n_ch_i(iinstr),iobs) + &
                             humi_k(:,1:rs%n_ch_i(iinstr),iobs) ,1) &
                    == 0)) THEN
                    PRINT *, ' WARNING   *** Normalization not possible. ***'
                ELSE
                  obs_rad(iradv)%plevel(icstart:icend,iobs) = 100.0_ireals * &
                    SUM((ABS(temp_k(:,1:rs%n_ch_i(iinstr),iobs) * fdbk_e_bg_t) + &
                       ABS(humi_k(:,1:rs%n_ch_i(iinstr),iobs) * fdbk_e_bg_rh)) * &
                      SPREAD(pres(:,iobs),2,rs%n_ch_i(iinstr)),1) / &
                    SUM((ABS(temp_k(:,1:rs%n_ch_i(iinstr),iobs) * fdbk_e_bg_t) + &
                       ABS(humi_k(:,1:rs%n_ch_i(iinstr),iobs) * fdbk_e_bg_rh)), 1)
                ENDIF    
              CASE (3)
                IF (ANY( SUM(temp_k(:,1:rs%n_ch_i(iinstr),iobs) + &
                             humi_k(:,1:rs%n_ch_i(iinstr),iobs) ,1) &
                    == 0)) THEN
                    PRINT *, ' WARNING   *** Normalization not possible. ***'
                ELSE
                  obs_rad(iradv)%plevel(icstart:icend,iobs) = 100.0_ireals * &
                    SUM((SQUARE(temp_k(:,1:rs%n_ch_i(iinstr),iobs) * fdbk_e_bg_t) + &
                       SQUARE(humi_k(:,1:rs%n_ch_i(iinstr),iobs) * fdbk_e_bg_rh)) * &
                      SPREAD(pres(:,iobs),2,rs%n_ch_i(iinstr)),1) / &
                    SUM((SQUARE(temp_k(:,1:rs%n_ch_i(iinstr),iobs) * fdbk_e_bg_t) + &
                       SQUARE(humi_k(:,1:rs%n_ch_i(iinstr),iobs) * fdbk_e_bg_rh)), 1)
                ENDIF    
            END SELECT

            ! If the highest level has the highest sensitivity
            ! we peak most likely outside the cosmo levels
            ! so set the plevel to a very small pressure
            WHERE (MAXLOC(ABS(temp_k(:,1:rs%n_ch_i(iinstr),iobs)),1) == 1)
              obs_rad(iradv)%plevel(icstart:icend,iobs) = 0.000001
            END WHERE
          END DO
        ELSE
          WRITE (stderr,*) ' WARNING   *** rttov_k_ifc failed. ***'
          DO iobs = ipstart,ipend
            DO ichan = icstart,icend
              CALL change_use_rpt_rad(r,ichan,iobs,FL_OPERATOR)
              r%bt_fg(ichan, iobs) = rmdi
            END DO
          END DO  
        ENDIF  
      ENDDO

      DEALLOCATE(pres, profs, chans, emiss, emiss_k, temp_k, humi_k, t2m_k,q2m_k, &
                 stemp_k)
      
!-------------------------------------------------------------------------------
!- Section 4: Calculate Bias correction and bias corrected obs
!-------------------------------------------------------------------------------

      ! TODO: check if the following is correct
      WHERE ( r%state(:,ipstart:ipend) < ST_DISMISS)
        r%bt_bcor(:,ipstart:ipend) =  r%bt_obs(:,ipstart:ipend) - &
                                      r%bcor_(:,ipstart:ipend)
      END WHERE    
    ENDIF 
    
!-------------------------------------------------------------------------------
!- Section 4: Collect Results at PE 0
!-------------------------------------------------------------------------------
    CALL global_values(inobs,1,'SUM',imp_integers,&
                       icomm_world,0,ymsg,istat)

    IF (my_world_id /= 0) inobs = 0

    CALL construct(obs,inobs,istat,obs_rad(iradv)%radv)

    IF (istat == 0) ALLOCATE(obs%plevel(rs%n_chan,inobs), stat=istat)
   
    IF (istat == 0) THEN
      CALL gather_obs_rad_container(obs_rad(iradv), ipstart, ipend, &
                                    obs,           1,       ncalcs, &
                                    0,istat ) 
    ELSE
       ierrstat = 99140
       PRINT *, ' ERROR   *** failed to allocate memory ***'
    ENDIF

!-------------------------------------------------------------------------------
!- Section 5: Append data to feed back file (PE 0 only)
!-------------------------------------------------------------------------------
    IF (my_world_id == 0 .AND. inobs > 0 .AND. istat == 0) THEN
      r  => obs%radv
      rs => r%i

      IF (inobs /= ncalcs ) THEN
        ierrstat = 99104
        PRINT *, ' ERROR   *** obs number mismatch ***'
        EXIT
      ENDIF

      ALLOCATE(reports(inobs),                     &
               report_veridata(1,rs%n_chan,inobs), &
               report_body(rs%n_chan,inobs),       &
               stat = istat)

      IF (istat /= 0) THEN
        ierrstat = 99105
        PRINT *, ' ERROR   *** allocating memory failed. ***'
        EXIT
      ENDIF

  !-----------------------------------------------------------------------------
  !+ Section 5.1: Populate report structure with our values 
  !-----------------------------------------------------------------------------
      DO iobs=1,inobs
        reports(iobs)%offset        = (fdbk_offset(iradv) + iobs-1) * rs%n_chan + 1 
        reports(iobs)%header%i_body = (fdbk_offset(iradv) + iobs-1) * rs%n_chan + 1

        ! Dont allocate each pointer on its own, this
        ! causes memory fragmentation and tons of allocs (slow on sx9)
        ! use a container alloc instead and just point to the memory
        reports(iobs)%body => report_body(:,iobs)
        DO icnt=1,rs%n_chan
!US          reports(iobs)%body(icnt)%veri_data => report_veridata(1:1,icnt,iobs)
          reports(iobs)%body(icnt)%veri_data = report_veridata(1:1,icnt,iobs)
        END DO
       
      
       
        CALL diff_minutes(reftime%year,reftime%month,reftime%day,         &
                          reftime%hour, reftime%minute,                   &
                          r%date(iobs)     / 10000     ,                  &
                          MOD(r%date(iobs) / 100  ,100),                  &
                          MOD(r%date(iobs)        ,100),                  &
                          r%time(iobs)     / 10000     ,                  &
                          MOD(r%time(iobs) / 100  ,100),                  &
                          itype_calendar, reports(iobs)%header%time,      &
                          istat                                   )

        CALL diff_minutes(reftime%year,reftime%month,reftime%day,          &
                          reftime%hour, reftime%minute,                    &
                          r%date_d(iobs)     / 10000     ,                 &
                          MOD(r%date_d(iobs) / 100  ,100),                 &
                          MOD(r%date_d(iobs)        ,100),                 &
                          r%time_d(iobs)     / 10000     ,                 &
                          MOD(r%time_d(iobs) / 100  ,100),                 &
                          itype_calendar, reports(iobs)%header%time_dbase, &
                          istat                                   )
        
        report_body(:,iobs)%e_o = SQRT(rs%var(:))
      END DO  

      DO iinstr=1,rs%n_instr
        DO ichan=1+rs%o_ch_i(iinstr),rs%n_ch_i(iinstr)+ rs%o_ch_i(iinstr)
          report_body(ichan,:)%level     = rs%chan(ichan)
          report_body(ichan,:)%level_sig = rs%instr_wmo(iinstr)
        END DO
      END DO  
 
      report_body(:,:)%varno     = VN_RAWBT
      report_body(:,:)%obs       = r%bt_bcor(:,:)
      report_body(:,:)%bcor      = - r%bcor_  (:,:)
      report_body(:,:)%level_typ = 0
      report_body(:,:)%qual      = imdi
      report_body(:,:)%plevel    = obs%plevel(:,:) !TODO

      WHERE (r%state(:,:) < ST_DISMISS)
        report_veridata(1,:,:) = r%bt_fg(:,:)
      ELSEWHERE
        report_veridata(1,:,:) = REAL(rmdi)
      END WHERE  
        
!      report_body(:,:)%check   = r%check(:,:)
      report_body(:,:)%check   = imdi
      report_body(:,:)%state   = r%state(:,:)
      report_body(:,:)%flags   = r%flags(:,:) 
       
      ! the folling three are combined ones of the aboves?!
      ! how combine them? TODO

      reports(:)%header%r_state = MINVAL(r%state(:,:),1)
      reports(:)%header%r_flags = 0
      DO ichan=1,rs%n_chan
        reports(:)%header%r_flags = IOR(r%state(ichan,:),        &
                                        reports(:)%header%r_flags)
      END DO  
      reports(:)%header%r_check = imdi

      reports(:)%len                  = rs%n_chan
      reports(:)%header%l_body        = rs%n_chan
      reports(:)%header%n_level       = rs%n_chan
      reports(:)%header%data_category = fdbk_rad_cat
      reports(:)%header%sub_category  = fdbk_rad_subcat
      reports(:)%header%obstype       = OT_RAD
      reports(:)%header%codetype      = 0  ! TODO: check this value
      reports(:)%header%ident         = rs%satid
      reports(:)%header%statid        = sat_id2name(rs%satid)
      reports(:)%header%lat           = r%dlat
      reports(:)%header%lon           = r%dlon
      reports(:)%header%time_nomi     = reports(:)%header%time
      reports(:)%header%z_station     = -999
      reports(:)%header%z_modsurf     = r%h_fg(isurfsens,:) ! TODO: instr. dep. height
      reports(:)%header%sta_corr      = 0
      reports(:)%header%index_x       = obs%iob_tot(:)
      reports(:)%header%index_y       = obs%job_tot(:)

      DO iinstr=1,rs%n_instr
        IF (rs%grid == rs%instr(iinstr)) reports(:)%header%instype = iinstr
      END DO  

      reports(:)%header%phase         = r%fov(:)
      reports(:)%header%surftype      = 0

      ! if fl_surf flag ist set, then for all obs of one report
      WHERE (BTEST(r%flags(1,:),FL_SURF))
        reports(:)%header%surftype = IBSET(0,ST_MISMATCH)
      ELSEWHERE (r%h_fg(isurfsens,:) > 1000.)
        reports(:)%header%surftype = IBSET(0,ST_HIGHLAND)
      ELSEWHERE (BTEST(r%mdlsfc(:),MS_LAND))
        reports(:)%header%surftype = IBSET(0,ST_LAND)
      ELSEWHERE (BTEST(r%mdlsfc(:),MS_ICE))
        reports(:)%header%surftype = IBSET(0,ST_ICE)
      ELSEWHERE (BTEST(r%mdlsfc(:),MS_SEA))
        reports(:)%header%surftype = IBSET(0,ST_SEA)
      END WHERE
        
      reports(:)%header%source        = iradv
      reports(:)%header%subset        = r%fov(:)
      reports(:)%header%sat_zenit     = r%stzen(:)
      reports(:)%header%mdlsfc        = r%mdlsfc(:)
      reports(:)%header%sun_zenit     = r%sunzen

      IF (ASSOCIATED(r%center)) THEN
        reports(:)%header%center = r%center(:)
      ELSE
        reports(:)%header%center = -1
      ENDIF  

      IF (ASSOCIATED(r%subcenter)) THEN
        reports(:)%header%sub_center = r%subcenter(:)
      ELSE
        reports(:)%header%sub_center = -1
      ENDIF  

      IF (ASSOCIATED(r%scanl)) THEN
        reports(:)%header%record = r%scanl(:) 
      ELSE 
        reports(:)%header%record = imdi 
      ENDIF  
      
      IF (ASSOCIATED(r%flg_prc)) THEN
        reports(:)%header%flg_1dvar     = r%flg_prc(:)
      ELSE  
        reports(:)%header%flg_1dvar     = imdi
      ENDIF

      IF (ASSOCIATED(r%flg_cld)) THEN
        reports(:)%header%flg_cld       = r%flg_cld(:)
      ELSE  
        reports(:)%header%flg_cld       = imdi
      ENDIF
      
      IF (ASSOCIATED(r%obsnum)) THEN
        reports(:)%header%obs_id        = r%obsnum(:)
      ELSE
        reports(:)%header%obs_id        = imdi
      ENDIF

     
  !-----------------------------------------------------------------------------
  !+ Section 5.2: Write data to file 
  !-----------------------------------------------------------------------------
      CALL write_report(fdbk(iradv), reports, inobs, inobs*rs%n_chan,    &
                    fdbk_offset(iradv)+1,fdbk_offset(iradv)*rs%n_chan+1, &
                    imdi, REAL(rmdich), istat,ymsg)
      
      DEALLOCATE(reports, report_veridata, report_body)
      fdbk_offset(iradv) = fdbk_offset(iradv) + inobs
      
      IF (istat /= NF_NOERR) THEN
        ierrstat = 99105
        WRITE (stderr,*) ' ERROR   *** write_report error ',istat,':',TRIM(ymsg),'***'
      ENDIF
  !-----------------------------------------------------------------------------
  !+ Section 5.3: Close file if done 
  !-----------------------------------------------------------------------------
      ! close file if all obs have been written
      IF (fdbk_offset(iradv) == SUM(nobspe_tot(:,iradv))) THEN
        CALL close_fdbk(fdbk(iradv))
      ENDIF  
    ENDIF  

    CALL destruct(obs)
    
    IF (ierrstat /= 0) EXIT
  ENDDO

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
  IF (.NOT.ALLOCATED(h_ice)) DEALLOCATE(h_ice_dummy)
  IF (.NOT.ASSOCIATED  (qg)) DEALLOCATE(qg_dummy)
#endif

CONTAINS
  ELEMENTAL REAL(KIND=ireals) FUNCTION square(val)
  REAL (KIND=ireals),INTENT(IN) :: val
    square = val * val
  END FUNCTION
END SUBROUTINE calc_obs_satpp

#endif

!==============================================================================
!+ Procedure in "src_obs_sat" for preventing supersaturation
!------------------------------------------------------------------------------
ELEMENTAL REAL(KIND=ireals) FUNCTION q_sat(p, t)
  REAL(KIND=ireals), intent(in) :: p, t
  REAL(KIND=ireals) :: p_sat

  p_sat = B1 * exp(B2w * (t - B3) / (t - B4w))
  q_sat = Rdv * p_sat / max((p - O_m_rdv * p_sat), 1.0_ireals)

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------
END FUNCTION q_sat


!==============================================================================
!+ Procedure in "src_obs_sat" for preparing rttov input fields
!------------------------------------------------------------------------------

SUBROUTINE prepare_rttov_input(iob_loc, job_loc, extrp_type, nlev,      &
                               pp, ps,t, t_g, qv, qc, qi, qs, qg, h_ice,&
                               pres, temp, humidity,                    &
                               t2m, hum2m, psurf, s_hgt, u10m, v10m,t_s,&
                               stype,wtype,                             &
                               cloud, cfrac,idg,ish,                    &
                               ierrstat)


!------------------------------------------------------------------------------
!
! Description:
!   This subroutine calculates the radiances
!
!------------------------------------------------------------------------------

! Parameters
  INTEGER   (KIND=iintegers),   INTENT (IN)      ::        &
    iob_loc(:),   & ! zonal indices of profile
    job_loc(:),   & ! meridional indices of profile
    extrp_type      ! Extrapolation type to use
  REAL   (KIND=ireals),   INTENT (IN)   ::        &
    pp(:,:,:),    & ! pressure state
    ps(:,:),      & ! surface pressure
    t(:,:,:),     & ! temperature state
    t_g(:,:),     & ! surface temperature state
    qv(:,:,:),    & ! humidity state
    qc(:,:,:),    & !
    qi(:,:,:),    & !
    qs(:,:,:),    & ! 
    h_ice(:,:),   & ! 
    qg(:,:,:)       ! 
  INTEGER   (KIND=iintegers),   INTENT (INOUT) :: &
    nlev              ! number of levels
  INTEGER   (KIND=iintegers),   INTENT (OUT) :: &
    stype(:),       & ! surface type
    wtype(:)          ! water type of grounde
  REAL   (KIND=ireals),   INTENT (OUT) :: &
    pres(:,:),      & ! presure profiles
    temp(:,:),      & ! temperature profiles
    humidity(:,:),  & ! humidity profiles
    t2m(:),         & ! 2m temperature
    hum2m(:),       & ! 2m humidity
    psurf(:),       & ! surface pressure
    s_hgt(:),       & ! surface height
    u10m(:),        & ! 10m u wind-component
    v10m(:),        & ! 10m u wind-component
    t_s(:)            ! surface temperature
  REAL   (KIND=ireals),   INTENT (OUT), OPTIONAL :: &
    cloud(:,:,:),   & ! cloud water/ice content
    cfrac(:,:,:)      ! cloud fraction
  INTEGER   (KIND=iintegers),  INTENT (OUT), OPTIONAL :: &
    idg(:),     &     ! ice water cloud scheme
    ish(:)            ! ice crystal shape
  INTEGER   (KIND=iintegers),   INTENT (OUT), OPTIONAL ::   &
    ierrstat        ! error status
! Constants
  REAL   (KIND=ireals), PARAMETER :: &
    psmax = 1099.9_ireals ! Max surface pressure [hPa]
! Local variables  
  INTEGER (KIND=iintegers) :: &
    i,j,k,          & ! indices for gridpoint, level
    iprof,          & ! indice for profile
    nprof,          & ! number of profiles
    n_add             ! additional profiles on model top
  REAL (KIND=ireals) :: &
    avg_p_1,        & ! Average presure on highest model level
    avg_p_2,        & ! Average presure on 2nd-highest model level
    clch,           & ! 
    clwh,           & ! 
    qsh,            & ! 
    qih,            & ! 
    w,              & ! 
    w_con,          & ! 
    clw_conh          ! 
  LOGICAL :: &
    lflag             !
    
    
!-------------------------------------------------------------------------------
!- End of header -
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Begin SUBROUTINE prepare_rttov_input
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!- Section 1: Initialize some variables
!-------------------------------------------------------------------------------

  n_add = 0
  nprof = SIZE(iob_loc)

  ierrstat = 0

!-------------------------------------------------------------------------------
!- Section 1.1: Calculate number of extra levels to add to the top
!-------------------------------------------------------------------------------

  avg_p_1 = 0
  avg_p_2 = 1
 
  IF ( extrp_type /= extrp_const) THEN
    DO iprof = 1,nprof
      i = iob_loc(iprof)
      j = job_loc(iprof)

      avg_p_1 = avg_p_1 + p0(i,j,1) + pp(i,j,1)
      avg_p_2 = avg_p_2 + p0(i,j,2) + pp(i,j,2)
    ENDDO  

    avg_p_1 = avg_p_1 / REAL(SIZE(iob_loc))
    avg_p_2 = avg_p_2 / REAL(SIZE(iob_loc))

    IF (avg_p_1 < p_top) THEN
      n_add = MAX(1,NINT(ABS((p_top - avg_p_1)/(avg_p_1-avg_p_2))))
    ENDIF
  ENDIF  

  IF (nlev < n_add + ke + 1) THEN
    ierrstat = 99101
    PRINT *, ' ERROR   *** to less levels available. ***'
    RETURN
  ENDIF
  
  ! additional levels on top plus model levels plus 1 surface level
  nlev = n_add + ke + 1

!-------------------------------------------------------------------------------
!- Section 2.1: Fill atmospheric fields
!-------------------------------------------------------------------------------
  DO k=1, ke
    DO iprof = 1,nprof
      i = iob_loc(iprof)
      j = job_loc(iprof)

      pres    (k+n_add,iprof) = (p0(i,j,k) + pp(i,j,k)) / 100.0_ireals
      humidity(k+n_add,iprof) = qv(i,j,k)
      temp    (k+n_add,iprof) = t(i,j,k)
    ENDDO  
  ENDDO

  DO k = 1,n_add
    DO iprof = 1,nprof
      pres(k,iprof) = (p_top + (k-1) * (pres(n_add+1,iprof) - p_top) &
                      / REAL(n_add)) / 100.0_ireals

      IF (extrp_type == extrp_lin) THEN
        temp(k,iprof) = temp(n_add + 1, iprof) + (n_add - k + 1) * &
                        (temp(n_add+1, iprof) - temp(n_add+2, iprof))
        humidity(k,iprof) = humidity(n_add + 1, iprof) + (n_add - k + 1) * &
                            (humidity(n_add+1, iprof) - humidity(n_add+2, iprof))
      ELSEIF (extrp_type == extrp_clim) THEN
        temp(k,iprof) = temp(n_add + 1, iprof) + (n_add - k + 1) * &
                        (t_top - temp(n_add+1, iprof)) / REAL(n_add)
        
        humidity(k,iprof) = humidity(n_add + 1, iprof) + (n_add - k + 1) * &
                            (q_top - humidity(n_add+1, iprof)) / REAL(n_add)
      ELSE                      
        temp(k,iprof)     = temp(n_add + 1, iprof)
        humidity(k,iprof) = humidity(n_add + 1, iprof)
      ENDIF                  

      temp(k,iprof) = MIN(tmax,MAX(tmin,temp(k,iprof)))
      humidity(k,iprof) = MAX(qmin,MIN(qmax,                               &
                                       q_sat(pres(k,iprof), temp(k,iprof)),&
                                       humidity(k,iprof)))
    ENDDO  
  ENDDO

!-------------------------------------------------------------------------------
!- Section 2.2: Fill Surface information
!-------------------------------------------------------------------------------

  lflag = ALL(SHAPE(h_ice) > 0)

  DO iprof = 1,nprof
    i = iob_loc(iprof)
    j = job_loc(iprof)
    
    psurf(iprof) = ps(i,j) / 100.0_ireals
    psurf(iprof) = MIN(psurf(iprof), psmax)

    pres(nlev,iprof)     = psurf(iprof)
    humidity(nlev,iprof) = qv_2m(i,j)
    temp(nlev,iprof)     = t_2m(i,j)

    t2m(iprof)   = t_2m(i,j)
    hum2m(iprof) = MIN(qv_2m(i,j), q_sat(psurf(iprof), t2m(iprof)))

    s_hgt(iprof)         = hsurf(i,j) * 0.001_ireals
    u10m(iprof)          = u_10m(i,j)
    v10m(iprof)          = v_10m(i,j)
    t_s(iprof)           = t_g(i,j)

    IF (soiltyp(i,j) == 9) THEN
      stype(iprof) = 1
      IF (lflag) THEN
        IF (h_ice(i,j) > 0.01_ireals) stype(iprof) = 2
      ENDIF  
    ELSE IF (soiltyp(i,j) == 10) THEN   
      stype(iprof) = 2
    ELSE  
      stype(iprof) = 0
    ENDIF

    IF (lseaice .AND. stype(iprof) /= 0) THEN
      IF (lseamask(i,j)) THEN
        wtype(iprof) = 1
      ELSE
        wtype(iprof) = 0
      ENDIF
    ELSE  
      wtype(iprof) = 0
    ENDIF  
  ENDDO

!-------------------------------------------------------------------------------
!- Section 2.3: Convert Humidity
!-------------------------------------------------------------------------------
  humidity(1:nlev,:) = humidity(1:nlev,:) / (1.0_ireals - humidity(1:nlev,:)) &
                       * rcnw
  hum2m(:) = hum2m(:) / (1.0_ireals - hum2m(:)) * rcnw

!-------------------------------------------------------------------------------
!- Section 2.4: Compute Cloud Cover an Contents
!-------------------------------------------------------------------------------
  IF (present(cloud)) THEN
!CDIR COLLAPSE  
    cloud(:,:,:) = 0._ireals
!CDIR COLLAPSE  
    cfrac(:,:,:) = 0._ireals

    idg  (:)     = iwc2effdiam
    ish  (:)     = iceshape

    lflag = ALL(SHAPE(qg) > 0)

    DO k=1, ke
      DO iprof = 1,nprof
        i = iob_loc(iprof)
        j = job_loc(iprof)

        qsh = qs(i,j,k)

        IF (lflag) qsh = qsh + qg(i,j,k)

        IF (qsh > 1.0E-7_ireals) THEN
          w = 1.0_ireals
        ELSE
          w = clc_sgs(i,j,k)
        ENDIF

        IF (clw_con(i,j,k) > 0._ireals) THEN
          w_con = clc_con(i,j,k) * (1.0_ireals - w)
        ELSE
          w_con = 0._ireals
        ENDIF

        clch = MAX(0.0_ireals, MIN(1.0_ireals,w + w_con))

        IF (clch > 0._ireals) THEN
          ! Liquid Wate content
          clwh = qc(i,j,k) * 1000._ireals * rho(i,j,k) 

          ! Convective Clouds
          IF (lcon_clw .AND. clc_con(i,j,k) > 0._ireals) THEN
            clw_conh = clw_con(i,j,k) * 1000._ireals * rho(i,j,k)
          ELSE
            clw_conh = 0._ireals 
          ENDIF  

          ! Ice content
          IF(qi(i,j,k) > 1.E-7) THEN
            qih = (qi(i,j,k) + qsh) * 1000._ireals * rho(i,j,k)
          ELSE
            qih = (            qsh) * 1000._ireals * rho(i,j,k)
          ENDIF  

          ! More than 10g/m^3 ice/snow is physically not reasonable
          ! and might cause a crash of RTTOV
          qih = min(qih,10._ireals)

#ifdef RTTOV10
          IF (rttov_ifc_version >= 10) THEN
            ! RTTOV10 can handle multiple cloud types at the same time
            cloud(6,k+n_add,iprof) = qih
            cloud(3,k+n_add,iprof) = clw_conh
            cloud(1,k+n_add,iprof) = clwh
            ! the cloud fraction must be specified once only - not type
            ! specific
            cfrac(1,k+n_add,iprof) = clch
#else
          IF (.false.) THEN
#endif          
          ELSEIF ( (clwh + clw_conh) <= 0._ireals .AND. qih > 0._ireals) THEN
            cloud(6,k+n_add,iprof) = qih
            cfrac(6,k+n_add,iprof) = clch
          ELSEIF ( (clwh + qih) <= 0._ireals   .AND. clw_conh > 0._ireals) THEN
            cloud(3,k+n_add,iprof) = clw_conh
            cfrac(3,k+n_add,iprof) = clch
          ELSEIF ( (clwh + qih + clw_conh) > 0._ireals) THEN
            cloud(1,k+n_add,iprof) = ((clwh + qih) * w + clw_conh * w_con) / (w + w_con)
            cfrac(1,k+n_add,iprof) = clch
          ENDIF
        ENDIF
      ENDDO  
    ENDDO
  ENDIF

!------------------------------------------------------------------------------
!- End of the Subroutine
!------------------------------------------------------------------------------

END SUBROUTINE prepare_rttov_input
!------------------------------------------------------------------------------
! End of module src_obs_rad
!------------------------------------------------------------------------------

END MODULE src_obs_rad
