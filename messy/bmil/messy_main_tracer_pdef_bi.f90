! Note: not all MBMs use CHANNEL
#if defined(ECHAM5)
#define MESSYCHANNELS
#endif
#ifdef COSMO
#define MESSYCHANNELS
#endif
#ifdef MBM_MPIOM
#define MESSYCHANNELS
#endif
#ifdef BLANK
#define MESSYCHANNELS
#endif
#ifdef MBM_CLAMS
#define MESSYCHANNELS
#endif
#ifdef CESM1
#define MESSYCHANNELS
#endif
#ifdef VERTICO
#define MESSYCHANNELS
#endif
#ifdef ICON
#define MESSYCHANNELS
#endif

!***********************************************************************
MODULE messy_main_tracer_pdef_bi
!***********************************************************************

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_tracer,        ONLY: NMAXSETID
  USE messy_main_tracer_pdef
#ifdef MESSYCHANNELS
  USE messy_main_channel_mem,   ONLY: dom_curid
#endif

  IMPLICIT NONE
  PRIVATE

#ifndef MESSYCHANNELS
  INTEGER, PARAMETER :: dom_curid = 1
#endif

! SET TIMELEVEL over which PDEF integrates the tracer mass
#if defined(COSMO) || defined(ICON)
  LOGICAL, PARAMETER :: int_newtl = .TRUE.
#else
  LOGICAL, PARAMETER :: int_newtl = .FALSE.
#endif

  ! ROUTINES
  PUBLIC :: main_tracer_pdef_initialize
  PUBLIC :: main_tracer_pdef_init_mem
  PUBLIC :: main_tracer_pdef_init_cpl
  PUBLIC :: main_tracer_pdef_global_start
  PUBLIC :: main_tracer_pdef_global_end
  PUBLIC :: main_tracer_pdef_free_mem

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_initialize

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i
#if defined(ECHAM5) || defined(BLANK) || defined(MBM_CLAMS) || defined(VERTICO) || defined(CESM1) || defined(MBM_MPIOM)
    CHARACTER(LEN=*), PARAMETER :: trigstr = 'next'
#endif
#if defined(COSMO) || defined(ICON)
    CHARACTER(LEN=*), PARAMETER :: trigstr = 'present'
#endif

    ! INITIALIZE CTRL_PDEF
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL tracer_pdef_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(L_DIAGOUT, p_io)
    CALL p_bcast(L_TASK_AGGREGATE, p_io)
    DO i=1, NMAXSETID
       CALL p_bcast(TPD_DEFAULT(i)%set, p_io)
       CALL p_bcast(TPD_DEFAULT(i)%name, p_io)
       CALL p_bcast(TPD_DEFAULT(i)%subname, p_io)
       CALL p_bcast(TPD_DEFAULT(i)%lswitch, p_io)
       CALL p_bcast(TPD_DEFAULT(i)%rtol, p_io)
    END DO
    DO i=1, NMAXTPDEF
       CALL p_bcast(TPD(i)%set, p_io)
       CALL p_bcast(TPD(i)%name, p_io)
       CALL p_bcast(TPD(i)%subname, p_io)
       CALL p_bcast(TPD(i)%lswitch, p_io)
       CALL p_bcast(TPD(i)%rtol, p_io)
    END DO

  END SUBROUTINE main_tracer_pdef_initialize
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_init_mem

#ifdef MESSYCHANNELS
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    ! + if task-aggregation is OFF
    USE messy_main_channel_repr,       ONLY: new_representation         &
                                           , set_representation_decomp  &
                                           , IRANK, PIOTYPE_COL
    USE messy_main_channel_repr,       ONLY: REPR_UNDEF
    USE messy_main_channel_dimensions, ONLY: new_dimension, DIMID_UNDEF &
                                           , get_dimension_info
    USE messy_main_channel_bi,         ONLY: DC_AG
#ifdef ICON
    USE messy_main_channel_mem,        ONLY: dom_unbound
#endif
    USE messy_main_channel_bi,         ONLY: SCALAR
#endif
    ! - if task-aggregation is OFF
    USE messy_main_mpi_bi,           ONLY: p_nprocs
    USE messy_main_blather_bi,       ONLY: error_bi, info_bi
    USE messy_main_constants_mem,    ONLY: M_air
    USE messy_main_tracer,           ONLY: NSETID, TRSET, R_molarmass &
                                         , param2string &
                                         , AIR, AEROSOL, OCEAN, CLOUD

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_init_mem'
    INTEGER                               :: jt
    INTEGER                               :: n
    INTEGER                               :: ntrac
    INTEGER                               :: rid
#ifdef MESSYCHANNELS
    INTEGER                               :: status
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: mem
    ! + if task-aggregation is OFF
    INTEGER, SAVE  :: REPR_TRACER_PDEF
    INTEGER        :: DIMID_NCPUS
    INTEGER        :: dimlen
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    ! - if task-aggregation is OFF
#endif

    CALL start_message_bi(submodstr, 'INITIALIZE MEMORY', substr)

    CALL tracpdef_initmem(p_nprocs)

#ifdef MESSYCHANNELS
    task_aggr: IF (.NOT. L_TASK_AGGREGATE) THEN
       ! shared between TRACER_PDEF, QTIMER AND RND ...
       CALL get_dimension_info(status, 'NCPUS', DIMID_NCPUS, dimlen)
       IF (status == 0) THEN
          IF (dimlen /= p_nprocs) &
               CALL error_bi('dimension NCPUS exists with wrong length',substr)
       ELSE
          CALL new_dimension(status, DIMID_NCPUS, 'NCPUS', p_nprocs)
          CALL channel_halt(substr, status)
       END IF

       ! NEW REPRESENTATION
       ! NOTE: DC_AG requires p_pe at 4th rank ...
       CALL new_representation(status, REPR_TRACER_PDEF, 'REPR_TRACER_PDEF' &
            , rank = 1, link = '---x', dctype = DC_AG             &
            , dimension_ids = (/ DIMID_NCPUS /)                   &
            , ldimlen       = (/ 1 /)                             &
            , axis          = '---N'                              &
#ifdef ICON
            , dom_id        = dom_unbound                         &
#endif
            )
       CALL channel_halt(substr, status)

       nseg = 1
       ALLOCATE(start(nseg,IRANK))
       ALLOCATE(cnt(nseg,IRANK))
       ALLOCATE(meml(nseg,IRANK))
       ALLOCATE(memu(nseg,IRANK))

       start(:,:) = 1
       cnt(:,:)   = 1
       meml(:,:)  = 1
       memu(:,:)  = 1

       start(:,1) = 1
       cnt(:,4)   = 1 ! p_nprocs
       meml(:,1)  = 1
       memu(:,4)  = 1 ! p_nprocs

       CALL set_representation_decomp(status, REPR_TRACER_PDEF &
            , start, cnt, memu, meml, piotype=PIOTYPE_COL)
       CALL channel_halt(substr, status)

       DEALLOCATE(start) ; NULLIFY(start)
       DEALLOCATE(cnt)   ; NULLIFY(cnt)
       DEALLOCATE(meml)  ; NULLIFY(meml)
       DEALLOCATE(memu)  ; NULLIFY(memu)
    END IF task_aggr
#endif

    set_loop: DO n=1, NSETID

       ntrac  = TRSET(n)%ntrac
       IF (ntrac == 0) CYCLE

#ifdef MESSYCHANNELS
       ! CHANNEL AND OBJECTS
       CALL info_bi(' new channel: '//submodstr//'_'//TRIM(TRSET(n)%name) &
            , substr)
       IF (L_TASK_AGGREGATE) THEN
          rid = SCALAR
       ELSE
          rid = REPR_TRACER_PDEF
       END IF
       CALL new_channel(status, submodstr//'_'//TRIM(TRSET(n)%name) &
            , reprid=rid)
       CALL channel_halt(substr, status)
#endif

       tracer_loop: DO jt=1, ntrac

          ! check, if tracer has properties
          IF (TRSET(n)%ti(jt)%tp%meta%cask_r(R_molarmass) <= 0.0_dp) THEN
             CALL error_bi('TRACER '//&
                  TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)//&
                  ' IN TRACER SET '//TRIM(TRSET(n)%name)//&
                  ' HAS MOLAR MASS <= 0.0' , substr)
          ENDIF

          SELECT CASE(TRSET(n)%ti(jt)%tp%ident%medium)
             !
          CASE(AIR, AEROSOL, CLOUD)
             !
             ! WHICH UNIT ?
             SELECT CASE(TRIM(TRSET(n)%ti(jt)%tp%ident%unit))
             CASE('kg/kg','kg kg-1')
                XWRK(n)%unit(jt) = 'kg'
                XWRK(n)%scalf(jt) = 1.0_DP
             CASE('mol/mol','mol mol-1')  ! -> kg/kg
                XWRK(n)%unit(jt) = 'kg'
                XWRK(n)%scalf(jt) = &
                     TRSET(n)%ti(jt)%tp%meta%cask_r(R_molarmass) / M_air
             ! mz_sg_20160122+: added units for passive prod./loss tracers
             CASE('mol/mol/s','mol mol-1 s-1')  ! -> kg/kg/s
                XWRK(n)%unit(jt) = 'kg/s'
                XWRK(n)%scalf(jt) = &
                     TRSET(n)%ti(jt)%tp%meta%cask_r(R_molarmass) / M_air
             ! mz_sg_20160122-
             CASE('1/mol','mol-1')    ! -> 1/kg
                XWRK(n)%unit(jt) = '1'
                XWRK(n)%scalf(jt) = 1.0_DP / M_air
             CASE DEFAULT
                XWRK(n)%unit(jt) = &
                     'kg*('//TRIM(TRSET(n)%ti(jt)%tp%ident%unit)//')'
                XWRK(n)%scalf(jt) = 1.0_DP
                XWRK(n)%lok(jt) = .FALSE.
             END SELECT
             !
          CASE(OCEAN)
             !
             SELECT CASE(TRIM(TRSET(n)%ti(jt)%tp%ident%unit))
             CASE('kmol/m^3')
                XWRK(n)%unit(jt) = 'kg'
                XWRK(n)%scalf(jt) = &
                     TRSET(n)%ti(jt)%tp%meta%cask_r(R_molarmass)
             CASE DEFAULT
                XWRK(n)%unit(jt) = &
                     'kg*('//TRIM(TRSET(n)%ti(jt)%tp%ident%unit)//')'
                XWRK(n)%scalf(jt) = 1.0_DP
                XWRK(n)%lok(jt) = .FALSE.
             END SELECT
             !
          CASE DEFAULT
             !
             CALL error_bi('integration not implemented for medium '//&
                param2string(TRSET(n)%ti(jt)%tp%ident%medium, 'medium'),substr)
             !
          END SELECT

#ifdef MESSYCHANNELS
          mem => XWRK(n)%mem_mass_p(jt,:,:,:,:)
          CALL new_channel_object(status, submodstr//'_'//TRIM(TRSET(n)%name) &
               , 'MP_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)  &
               , mem=mem)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)     &
               , 'MP_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'long_name', c='positive mass of '&
               &//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname))
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)     &
               , 'MP_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'units', c=TRIM(XWRK(n)%unit(jt)))
          CALL channel_halt(substr, status)

          mem => XWRK(n)%mem_mass_n(jt,:,:,:,:)
          CALL new_channel_object(status, submodstr//'_'//TRIM(TRSET(n)%name) &
               , 'MN_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)  &
               , mem=mem)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)     &
               , 'MN_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'long_name', c='negative mass of '&
               &//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname))
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)  &
               , 'MN_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'units', c=TRIM(XWRK(n)%unit(jt)))
          CALL channel_halt(substr, status)

          mem => XWRK(n)%mem_mass_c(jt,:,:,:,:)
          CALL new_channel_object(status, submodstr//'_'//TRIM(TRSET(n)%name) &
               , 'MC_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)  &
               , mem=mem)
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)  &
               , 'MC_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'long_name', c='mass correction of '&
               &//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname))
          CALL channel_halt(substr, status)
          !
          CALL new_attribute(status, submodstr//'_'//TRIM(TRSET(n)%name)  &
               , 'MC_'//TRIM(TRSET(n)%ti(jt)%tp%ident%fullname) &
               , 'units', c=TRIM(XWRK(n)%unit(jt)))
          CALL channel_halt(substr, status)
#endif

       END DO tracer_loop

    END DO set_loop

    CALL end_message_bi(submodstr, 'INITIALIZE MEMORY',substr)

  END SUBROUTINE main_tracer_pdef_init_mem
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_init_cpl

    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: error_bi, warning_bi
#ifdef ICON
    USE messy_main_bmluse_bi,     ONLY: p_patch, nproma, grf_bdywidth_c &
                                      , get_indices_c, min_rlcell_int
    USE messy_main_tracer_mem_bi, ONLY: L_GPTRSTR
    USE messy_main_channel_mem,   ONLY: n_dom
#endif
#ifdef COSMO
    USE messy_main_grid_def_mem_bi,ONLY: ie, je, ke, istartpar,iendpar &
                                      , jstartpar, jendpar
#endif
#if defined(ICON) || defined(COSMO)
    USE messy_main_tracer,        ONLY: STRLEN_TRSET
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_init_cpl'
#if defined(ICON) || defined(COSMO)
    CHARACTER(LEN=STRLEN_TRSET) :: setname = ' '
    INTEGER :: status
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: mask
#endif
#if defined(ICON)
    INTEGER :: jg, jb
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx
#endif

    CALL start_message_bi(submodstr &
         , 'SET TRACER SWITCHES/TOLERANCES', substr)

    CALL tracpdef_settings(p_parallel_io)
    IF (.NOT. L_TASK_AGGREGATE) &
         CALL warning_bi('APPLICATION OF TOLERANCE CRITERIA NOT POSSIBLE' &
         , substr)

    CALL end_message_bi(submodstr &
         , 'SET TRACER SWITCHES/TOLERANCES', substr)

    CALL start_message_bi(submodstr, 'DEFINE MASKS', substr)

#ifdef ICON
    DO  jg=1, n_dom
       setname =L_GPTRSTR(jg)
       ALLOCATE(mask(nproma, p_patch(jg)%nlev, p_patch(jg)%nblks_c))
       mask = 0._dp
       ! computations on prognostic points
       i_startblk = p_patch(jg)%cells%start_block(grf_bdywidth_c+1)
       i_endblk   = p_patch(jg)%cells%end_block(min_rlcell_int-1)

       DO jb = i_startblk, i_endblk

          CALL get_indices_c(p_patch(jg), jb, i_startblk, i_endblk, &
               & i_startidx, i_endidx, grf_bdywidth_c+1, min_rlcell_int-1)

          mask(i_startidx:i_endidx,:,jb) = 1.0_dp

       END DO

       CALL tracpdef_mask(setname, mask, status)
       CALL mask_error(status, setname)
       DEALLOCATE(MASK)
    END DO
#endif

#ifdef COSMO
    ALLOCATE(mask(ie,je,ke))
    mask(:,:,:) = 0._dp
    mask(istartpar:iendpar, jstartpar:jendpar,:) = 1._dp
    setname = 'gp'
    CALL tracpdef_mask(setname, mask, status)
    CALL mask_error(status, setname)
    DEALLOCATE(MASK)
#endif

    CALL end_message_bi(submodstr, 'DEFINE MASKS', substr)

CONTAINS

  SUBROUTINE mask_error(status, setname)

    INTEGER, INTENT(IN) :: status
    CHARACTER(LEN=*)    :: setname

    SELECT CASE (status)
    CASE (0)
       ! ALL IS FINE
    CASE(2)
       CALL error_bi( 'tracer set '//TRIM(setname)//' not found'&
              , 'tracpdef_mask')
    CASE(3)
       CALL warning_bi( 'tracer number in tracer set '//TRIM(setname)//&
            &'''is zero', 'tracpdef_mask')
    CASE(4)
       CALL error_bi( 'DIMENSIONS OF MASK PROVIDED for TRACER-SET '&
            &//TRIM(setname)//' DO NOT FIT', 'tracpdef_mask')
    CASE DEFAULT
       CALL error_bi( 'unspecified error', 'tracpdef_mask')
    END SELECT

  END SUBROUTINE mask_error

  END SUBROUTINE main_tracer_pdef_init_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_global_start

#if defined(MESSYCHANNELS)
    USE messy_main_channel,          ONLY: get_channel_output, STRLEN_CHANNEL
    USE messy_main_channel_mem,      ONLY: l_dom
    USE messy_main_channel_error_bi, ONLY: channel_halt
    !
    USE messy_main_tracer,           ONLY: NSETID, TRSET
    USE messy_main_tools,            ONLY: split_name_domain
#endif

    IMPLICIT NONE

#if defined(MESSYCHANNELS)
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_global_start'
    INTEGER :: status
    INTEGER :: ntrac
    LOGICAL :: lout, lstat
    INTEGER :: n
    CHARACTER(LEN=STRLEN_CHANNEL) :: zname
    INTEGER :: zdomid
    LOGICAL :: l_this_domain

    set_loop: DO n=1, NSETID

       ntrac  = TRSET(n)%ntrac
       IF (ntrac == 0) CYCLE

       IF (l_dom) THEN
          ! additional condition for domain-aware basemodels
          CALL split_name_domain(status, &
               TRIM(TRSET(n)%name), zname, zdomid)
          l_this_domain = (zdomid == dom_curid)
       ELSE
          l_this_domain = .TRUE.
       END IF

       IF (l_this_domain) THEN
          CALL get_channel_output(status, submodstr//'_'//TRIM(TRSET(n)%name) &
               , lout, lstat)
          CALL channel_halt(substr, status)

          lnow(n) = (lout .OR. lstat)

          ! Note: Do NOT set lnow(n) to .FALSE. in ELSE-branch!
          !       This will deactivate calculation of other tracer set(s)!
       END IF

    END DO set_loop
#else
    lnow(:) = .TRUE.
#endif

  END SUBROUTINE main_tracer_pdef_global_start
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_global_end

    USE messy_main_mpi_bi,           ONLY: p_pe &
                                         , p_parallel_io
    ! + for task aggregation ON
    USE messy_main_mpi_bi,           ONLY: p_sum
    ! - for task aggregation ON
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_timer,            ONLY: time_step_len
    USE messy_main_tracer,           ONLY: NSETID, TRSET

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_pdef_global_end'
    INTEGER :: status
    INTEGER :: n
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zsum

    ! 1st step: sum on each PE

    CALL tracpdef_integrate(status, 1, time_step_len, p_pe, lnewtl=int_newtl)

    ! status always 0 ...

    ! BROADCAST RESULT TO ALL PEs
    DO n=1, NSETID
       IF (.NOT. lnow(n)) CYCLE
       IF (TRSET(n)%ntrac == 0) CYCLE
       IF (L_TASK_AGGREGATE) THEN
          ! Note: Use p_sum here, which is based on
          ! MPI_ALLREDUCE(...MPI_SUM...),
          !       ALTHOUGH this potentially causes parallel decomposition
          !       dependent results (for the sake of an increased performance!).
          !       For the purely diagnostic quantitiy 'global tracer mass',
          !       which has no feedback to the dynamics, this approach might be
          !       justifiable.
          ALLOCATE(zsum(SIZE(XWRK(n)%mass_pe,1),SIZE(XWRK(n)%mass_pe,2)))
          zsum(:,:) = p_sum(XWRK(n)%mass_pe(:,:,p_pe))
          XWRK(n)%mass_n(:) = zsum(:,1)
          XWRK(n)%mass_p(:) = zsum(:,2)
          XWRK(n)%mass_c(:) = zsum(:,3)
          DEALLOCATE(zsum)
       ELSE
          XWRK(n)%mass_n(:) = XWRK(n)%mass_pe(:,1,p_pe)
          XWRK(n)%mass_p(:) = XWRK(n)%mass_pe(:,2,p_pe)
          XWRK(n)%mass_c(:) = XWRK(n)%mass_pe(:,3,p_pe)
       END IF
    END DO

    IF (L_TASK_AGGREGATE) THEN
       ! application of global tolerance criteria not possible without
       ! task aggregation, because global mass is not available
       CALL tracpdef_integrate(status, 2, time_step_len, p_pe)
       IF (status /= 0) THEN
          CALL error_bi('negative mass of tracer exceeds tolerance', substr)
       END IF
    END IF

    ! diagnostic output
    CALL tracpdef_print(p_parallel_io)

  END SUBROUTINE main_tracer_pdef_global_end
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_pdef_free_mem

    IMPLICIT NONE

    CALL tracpdef_freemem

  END SUBROUTINE main_tracer_pdef_free_mem
! ----------------------------------------------------------------------

!***********************************************************************
END MODULE messy_main_tracer_pdef_bi
!***********************************************************************
