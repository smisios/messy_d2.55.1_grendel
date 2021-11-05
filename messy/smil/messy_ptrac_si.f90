! **********************************************************************
!
! MODULE FOR TESTING MASS CONSERVATION AND MONOTONICITY OF
! EULERIAN TRANSPORT ALGORITHMS WITH PROGNOSTIC TRACERS
!
! MESSy SMIL
!
! Author : Patrick Joeckel, MPICH, July  2002
!
! References:
!
! **********************************************************************
#include "messy_main_ppd_bi.inc"

#if defined(ECHAM5) || defined(COSMO) || defined(BLANK) || defined(CESM1) || defined(VERTICO) || defined(ICON)
#define MESSYCHANNELS
#endif

! **********************************************************************
MODULE messy_ptrac_si
! **********************************************************************

  ! BMIL
  USE messy_main_blather_bi,     ONLY: start_message_bi, end_message_bi
  USE messy_main_constants_mem,  ONLY: STRLEN_MEDIUM, STRLEN_XLONG
  USE messy_main_tracer,         ONLY: T_TRACPROP_IO, AIR, AMOUNTFRACTION &
                                     , SINGLE, STRLEN_TRSET
  USE messy_ptrac

  IMPLICIT NONE
  PRIVATE
  INTRINSIC :: NULL

  ! PUBLIC INTERFACE ROUTINES
  PUBLIC :: ptrac_initialize  ! initialize
  PUBLIC :: ptrac_new_tracer  ! define new tracers
  PUBLIC :: ptrac_init_memory ! allocate memory
  !PRIVATE :: ptrac_init_memory_gp
  !PRIVATE :: ptrac_init_memory_lg
  !PRIVATE :: ptrac_init_memory_cl
  PUBLIC :: ptrac_vdiff            ! mz_pj_20081030 ICON TEST PHASE

  ! PRIVATE HELPER ROUTINES
  !PRIVATE :: ptrac_read_nml_cpl   ! initialize /CPL/-namelist )

  TYPE T_TRAC_IO
     CHARACTER(LEN=10*STRLEN_TRSET+10) :: trset  = ''
     CHARACTER(LEN=10*STRLEN_XLONG+10) :: trlist = ''
     CHARACTER(LEN=STRLEN_MEDIUM)      :: unit   = ''
     INTEGER                           :: medium   = AIR
     INTEGER                           :: quantity = AMOUNTFRACTION
     INTEGER                           :: type     = SINGLE
     ! FOR AERSOL TRACERS
     CHARACTER(LEN=STRLEN_MEDIUM)      :: aermod = modstr
     INTEGER                           :: mode  = 0
     REAL(DP)                          :: radi  = 0.0_dp
     REAL(DP)                          :: sigma = 1.00_dp
     REAL(DP)                          :: dens  = 0.0_dp
  END type T_TRAC_IO

  ! WORSPACE TO SAVE INFORMATION FROM INITIALIZE FOR INIT_MEMORY
  TYPE T_XTRAC
     CHARACTER(LEN=2*STRLEN_TRSET+1),  DIMENSION(:), POINTER :: trs  => NULL()
     CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: trnb => NULL()
     CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: trns => NULL()
  END type T_XTRAC

  INTEGER, PARAMETER :: NMAXTRAC = 1000
  TYPE(T_TRAC_IO), DIMENSION(NMAXTRAC), SAVE :: TRAC
  TYPE(T_XTRAC),   DIMENSION(NMAXTRAC), SAVE :: XTRAC

  INTEGER, PARAMETER :: NMAXTRACPROP = 1000
  TYPE(T_TRACPROP_IO), DIMENSION(NMAXTRACPROP), SAVE :: TPROP

CONTAINS

! ------------------------------------------------------------------------
  SUBROUTINE  ptrac_initialize

    ! BMIL
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ptrac_initialize'
    INTEGER         :: iou    ! I/O unit
    INTEGER         :: status ! error status
    INTEGER         :: i

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL ptrac_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
       !
    END IF

    DO i=1, NMAXTRAC
       CALL p_bcast(TRAC(i)%trset, p_io)
       CALL p_bcast(TRAC(i)%trlist, p_io)
       CALL p_bcast(TRAC(i)%unit, p_io)
       CALL p_bcast(TRAC(i)%medium, p_io)
       CALL p_bcast(TRAC(i)%quantity, p_io)
       CALL p_bcast(TRAC(i)%type, p_io)
       !
       CALL p_bcast(TRAC(i)%aermod, p_io)
       CALL p_bcast(TRAC(i)%mode, p_io)
       CALL p_bcast(TRAC(i)%radi, p_io)
       CALL p_bcast(TRAC(i)%sigma, p_io)
       CALL p_bcast(TRAC(i)%dens, p_io)
    END DO

    DO i=1, NMAXTRACPROP
       CALL p_bcast(TPROP(i)%trset, p_io)
       CALL p_bcast(TPROP(i)%trlist, p_io)
       CALL p_bcast(TPROP(i)%caskname, p_io)
       CALL p_bcast(TPROP(i)%cont, p_io)
    END DO

  END SUBROUTINE ptrac_initialize
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE ptrac_new_tracer

    ! BMIL
    USE messy_main_blather_bi,       ONLY: warning_bi, error_bi
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    ! SMCL
    USE messy_main_tracer,        ONLY: new_tracer, get_tracer_set          &
                                      , full2base_sub, get_tracer, TR_EXIST &
                                      , STRLEN_TRSET, set_tracer_properties &
                                      , AEROSOL, set_tracer, I_AEROSOL_MODE &
                                      , R_AEROSOL_DENSITY, S_AEROSOL_MODEL  &
                                      , get_tracer_set_info
    USE messy_main_tools,         ONLY: strcrack
#if defined(ICON)
    USE messy_main_tools,         ONLY: int2str
    USE messy_main_channel_mem,   ONLY: dom_current
#endif

    IMPLICIT NONE
    INTRINSIC :: TRIM, ASSOCIATED
#if defined(ICON)
    INTRINSIC :: INDEX
#endif

    ! LOCAL
    INTEGER :: status, statmp
    CHARACTER(LEN=*), PARAMETER :: substr = 'ptrac_new_tracer'
    CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: trn => NULL()
    CHARACTER(LEN=2*STRLEN_TRSET+1),  DIMENSION(:), POINTER :: trs => NULL()
    INTEGER  :: i, n, j, m, k
    CHARACTER(LEN=STRLEN_MEDIUM)     :: basename    = '' ! name of tracer
    CHARACTER(LEN=STRLEN_MEDIUM)     :: subname     = '' ! OPTIONAL subname
    CHARACTER(LEN=STRLEN_MEDIUM)     :: sm = ''
    INTEGER                          :: idt
    LOGICAL                          :: l_enabled
    CHARACTER(LEN=2*STRLEN_TRSET+1)  :: trsname
#if defined(ICON)
    CHARACTER(LEN=2)                 :: istr = '  '
    INTEGER                          :: idx
    INTEGER                          :: jg
#endif

    CALL start_message_bi(modstr, 'TRACER REQUEST', substr)

    ! TRACERS
    lentry: DO i=1, NMAXTRAC

       IF (TRIM(TRAC(i)%trset)    == '') CYCLE
       IF (TRIM(TRAC(i)%trlist)   == '') CYCLE

       ! split tracer sets into individual names
       CALL strcrack(TRIM(TRAC(i)%trset), ';', XTRAC(i)%trs, m)

       ! split list of tracers into individual names
       CALL strcrack(TRIM(TRAC(i)%trlist), ';', trn, n)
       ALLOCATE(XTRAC(i)%trnb(n))
       ALLOCATE(XTRAC(i)%trns(n))

       lset: DO k=1, m

          trsname = TRIM(XTRAC(i)%trs(k))
#if defined(ICON)
          ! check for specific patch
          CALL int2str(istr, dom_current)
          idx = INDEX(trsname,'_D')
          IF (idx==0) THEN
             ! if patch is unspecified, use it as wildcard for all patches;
             ! -> need to append current patch number
             trsname = TRIM(trsname)//'_D'//istr
          ELSE
             ! get specific patch; cycle for all others
             READ(trsname(idx+2:idx+3),*) jg
             IF (jg /= dom_current) CYCLE
          ENDIF
#endif

!!$       CALL get_tracer_set(status, TRIM(XTRAC(i)%trs(k)))
          CALL get_tracer_set(status, TRIM(trsname))
          IF (status /= 0) THEN 
             CALL warning_bi(&
!!$               'TRACER SET '//TRIM(XTRAC(i)%trs(k))//' NOT AVAILABLE!' &
                  'TRACER SET '//TRIM(trsname)//' NOT AVAILABLE!' &
                  , substr)
             CYCLE
          END IF

          ltrac: DO j=1, n
             ! convert fullname into basename and subname ...
             CALL full2base_sub(status, trn(j), basename, subname)
             CALL tracer_halt(substr, status)
             ! ... and save it for later ...
             XTRAC(i)%trnb(j) = TRIM(basename)
             XTRAC(i)%trns(j) = TRIM(subname)

!!$          CALL new_tracer(status, TRIM(XTRAC(i)%trs(k)) &
             CALL new_tracer(status, TRIM(trsname)    &
                  , TRIM(basename), modstr        &
                  , idx = idt                     &
                  , subname = TRIM(subname)       &
                  , unit = TRIM(TRAC(i)%unit)     &
                  , medium = TRAC(i)%medium       &
                  , quantity = TRAC(i)%quantity   &
                  , type = TRAC(i)%type )

             IF (status == TR_EXIST) THEN
!!$             CALL get_tracer(statmp, TRIM(XTRAC(i)%trs(k)), TRIM(basename) &
                CALL get_tracer(statmp, TRIM(trsname), TRIM(basename) &
                     , subname = TRIM(subname), submodel = sm)
                CALL tracer_halt(substr, statmp)
                IF (TRIM(sm) /= modstr) THEN
                   CALL tracer_halt(substr, status)
                END IF
             END IF

             ! special for aerosol
             IF (TRAC(i)%medium == AEROSOL) THEN
!!$             CALL set_tracer(status, TRIM(XTRAC(i)%trs(k)), idt &
                CALL set_tracer(status, TRIM(trsname), idt &
                     , S_AEROSOL_MODEL, TRAC(i)%aermod)
                CALL tracer_halt(substr, status)
!!$             CALL set_tracer(status, TRIM(XTRAC(i)%trs(k)), idt &
                CALL set_tracer(status, TRIM(trsname), idt &
                     , I_AEROSOL_MODE, TRAC(i)%mode)
                CALL tracer_halt(substr, status)
!!$             CALL set_tracer(status, TRIM(XTRAC(i)%trs(k)), idt &
                CALL set_tracer(status, TRIM(trsname), idt &
                     , R_AEROSOL_DENSITY, TRAC(i)%dens)
                CALL tracer_halt(substr, status)
             ENDIF
          END DO ltrac

       END DO lset

       IF (ASSOCIATED(trn)) THEN
          DEALLOCATE(trn)
          NULLIFY(trn)
       END IF

    END DO lentry

    ! TRACER PROPERTIES
    entry_loop: DO i=1, NMAXTRACPROP

       IF (TRIM(TPROP(i)%trset)    == '') CYCLE
       IF (TRIM(TPROP(i)%trlist)   == '') CYCLE
       IF (TRIM(TPROP(i)%caskname) == '') CYCLE
       IF (TRIM(TPROP(i)%cont)     == '') CYCLE

       ! split tracer sets into individual names
       CALL strcrack(TRIM(TPROP(i)%trset), ';', trs, m)

       ! split list of tracers into individual names
       CALL strcrack(TRIM(TPROP(i)%trlist), ';', trn, n)

       set_loop: DO k=1, m

          trsname = trs(k)
#if defined(ICON)
          ! check for specific patch
          CALL int2str(istr, dom_current)
          idx = INDEX(trsname,'_D')
          IF (idx==0) THEN
             ! if patch is unspecied, use it as wildcard for all patches;
             ! -> need to append current patch number
             trsname = TRIM(trsname)//'_D'//istr
          ELSE
             ! get specific patch; cycle for all others
             READ(trsname(idx+2:idx+3),*) jg
             IF (jg /= dom_current) CYCLE
          ENDIF    
#endif
          CALL get_tracer_set(status, TRIM(trsname))
          IF (status /= 0) THEN 
             CALL warning_bi(&
                  'TRACER SET '//TRIM(trsname)//' NOT AVAILABLE!', substr)
             CYCLE
          END IF

          CALL get_tracer_set_info(status, TRIM(trsname), l_enabled = l_enabled)
          IF (.NOT. l_enabled) THEN 
             CALL warning_bi(&
                  'TRACER SET '//TRIM(trsname)//' NOT ENABLED!', substr)
             CYCLE
          END IF

          tracer_loop: DO j=1, n
             ! convert fullname into basename and subname
             CALL full2base_sub(status, TRIM(trn(j)), basename, subname)
             CALL tracer_halt(substr, status)

             ! check if tracer exists and is from ptrac
             ! Note: you can still overwrite tracer properties in tracer.nml
             CALL get_tracer(statmp, TRIM(trsname), TRIM(basename) &
                  , subname = TRIM(subname), submodel = sm)
             CALL tracer_halt(substr, statmp)
             IF (TRIM(sm) /= modstr) THEN
                CALL error_bi(&
                     'ATTEMPT TO MODIFY TRACER '//TRIM(trn(j))//&
                     &' FROM SUBMODEL '//TRIM(sm), substr)
             END IF
          END DO tracer_loop

       END DO set_loop

       IF (ASSOCIATED(trn)) THEN
          DEALLOCATE(trn)
          NULLIFY(trn)
       END IF

       IF (ASSOCIATED(trs)) THEN
          DEALLOCATE(trs)
          NULLIFY(trs)
       END IF

    END DO entry_loop

#ifndef ICON
    CALL set_tracer_properties(status, TPROP, lprint=p_parallel_io)
#else
    CALL set_tracer_properties(status, TPROP, lprint=p_parallel_io &
         , dom_id = dom_current)
#endif
    CALL tracer_halt(substr, status)

    CALL end_message_bi(modstr, 'TRACER REQUEST', substr)

  END SUBROUTINE ptrac_new_tracer
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE ptrac_init_memory

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED
    
    ! LOCAL
    INTEGER :: i

    CALL ptrac_init_memory_gp
    CALL ptrac_init_memory_lg
    CALL ptrac_init_memory_cl

    DO i=1, NMAXTRAC
       IF (ASSOCIATED(XTRAC(i)%trnb)) THEN
          DEALLOCATE(XTRAC(i)%trnb) ; NULLIFY(XTRAC(i)%trnb)
       END IF
       IF (ASSOCIATED(XTRAC(i)%trns)) THEN
          DEALLOCATE(XTRAC(i)%trns) ; NULLIFY(XTRAC(i)%trns)
       END IF
    END DO

  END SUBROUTINE ptrac_init_memory
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE ptrac_init_memory_gp

#ifdef MESSYCHANNELS
    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, ngpblks

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DC_BC, DC_GP &
                                         , DIMID_LON, DIMID_LEV, DIMID_LAT
#if ! (defined(BLANK) || defined(VERTICO) || defined(ICON))
    USE messy_main_channel_bi,       ONLY: gp_nseg, gp_start, gp_cnt &
                                         , gp_meml, gp_memu
#endif
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt

    ! SMCL
    USE messy_main_channel_dimensions,  ONLY: new_dimension
    USE messy_main_channel_repr,        ONLY: new_representation, AUTO &
                                            , set_representation_decomp  &
                                            , IRANK, PIOTYPE_COL &
                                            , repr_def_axes
    USE messy_main_channel,             ONLY: new_channel, new_channel_object &
                                            , new_attribute
    USE messy_main_tracer,              ONLY: get_tracer, AEROSOL &
                                            , I_AEROSOL_MODE

    IMPLICIT NONE
    INTRINSIC :: TRIM, SIZE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'ptrac_init_memory_gp'
    LOGICAL                          :: lchannel
    !
    ! CHANNEL OBJECTS FOR AEROSOL
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: radius_ptr => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: density_ptr => NULL()
    REAL(DP), DIMENSION(:),       POINTER :: sigma_ptr => NULL()
    !
    ! CHANNEL MANAGEMENT
    INTEGER :: DIMID_NMODE   ! DIMENSION ID: NUMBER OF (PSEUDO-) AEROSOL MODES
    INTEGER :: REPR_PTRAC_4D ! REPRESENTATION ID
    INTEGER :: REPR_PTRAC_1D ! REPRESENTATION ID
    !
    INTEGER                          :: status
    INTEGER                          :: i, m, n, j, k, jt
    CHARACTER(LEN=STRLEN_MEDIUM)     :: sm = ''
    INTEGER                          :: med
    INTEGER                          :: mode
    !
    ! NUMBER OF MODES
    INTEGER                      :: maxmode, jm
    CHARACTER(LEN=STRLEN_MEDIUM) :: reprname
    !
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    ! CHANNEL REQUIRED ?
    lchannel = .FALSE.
    maxmode = 0
    DO i=1, NMAXTRAC
       IF (TRIM(TRAC(i)%trset)    == '') CYCLE
       IF (TRIM(TRAC(i)%trlist)   == '') CYCLE
       m = SIZE(XTRAC(i)%trs)
       n = SIZE(XTRAC(i)%trnb)
       DO k=1, m
          IF (TRIM(XTRAC(i)%trs(k)) /= GPTRSTR) CYCLE
          DO j=1, n
             CALL get_tracer(status, GPTRSTR, TRIM(XTRAC(i)%trnb(j)) &
                  , subname = TRIM(XTRAC(i)%trns(j)), idx=jt &
                  , submodel=sm, medium=med)
             CALL tracer_halt(substr,status)
             IF ( (TRIM(sm) == modstr) .AND. (med == AEROSOL) ) THEN
                lchannel = .TRUE.
                CALL get_tracer(status, GPTRSTR, jt, I_AEROSOL_MODE, i=mode)
                CALL tracer_halt(substr,status)
                IF (mode > maxmode) maxmode = mode
             END IF
          END DO
       END DO
    END DO

    IF (.NOT. lchannel) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL new_channel(status, modstr//'_'//GPTRSTR)
    CALL channel_halt(substr, status)

    aerosol_cond: IF (maxmode > 0) THEN

       ! ADDITIONAL DIMENSION
       CALL new_dimension(status, DIMID_NMODE &
            , 'PTRAC_NMODE'//'_'//GPTRSTR, maxmode)
       CALL channel_halt(substr, status)

       ! NEW REPRESENTATIONS
       reprname = ''
       WRITE(reprname,'(a14,i3.3,a1,a2)') &
            'REPR_PTRAC_4D_',maxmode,'_',GPTRSTR
!!$#if defined(ECHAM5) || defined(CESM1)
!!$       CALL new_representation(status, REPR_PTRAC_4D, TRIM(reprname)  &
!!$            , rank = 4, link = 'xxxx', dctype = DC_GP                 &
!!$            , dimension_ids = (/ DIMID_LON, DIMID_LEV                 &
!!$            ,                    DIMID_NMODE, DIMID_LAT /)            &
!!$            , ldimlen       = (/ nproma, AUTO, AUTO, ngpblks   /)     &
!!$            , output_order  = (/ 3,1,4,2 /)                           &
!!$            , axis = 'XZNY'                                           &
!!$            )
!!$       CALL channel_halt(substr, status)
!!$#endif
!!$#ifdef COSMO
!!$       CALL new_representation(status, REPR_PTRAC_4D, TRIM(reprname)  &
!!$            , rank = 4, link = 'xxxx', dctype = DC_GP                 &
!!$            , dimension_ids = (/ DIMID_LON, DIMID_LAT                 &
!!$            ,                    DIMID_NMODE, DIMID_LEV /)            &
!!$            , ldimlen       = (/ nproma, ngpblks, AUTO, AUTO   /)     &
!!$            , output_order  = (/ 3,1,2,4 /)                           &
!!$            , axis = 'XYNZ'                                           &
!!$            )
!!$       CALL channel_halt(substr, status)
!!$#endif
       CALL new_representation(status, REPR_PTRAC_4D, TRIM(reprname)  &
            , rank = 4, link = 'xxxx', dctype = DC_GP                 &
            , dimension_ids = (/ &
            _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_NMODE) /) &
            , ldimlen       = (/ &
            _RI_XYZN_(nproma, ngpblks, AUTO, AUTO) /)   &
            , output_order  = (/ _IN_XYZN_, _IX_XYZN_   &        ! E: 3,1,4,2
                               , _IY_XYZN_, _IZ_XYZN_ /)       & ! C: 3,1,2,4
            , axis = repr_def_axes(_RI_XYZN_('X','Y','Z','N')) &
            )
       CALL channel_halt(substr, status)
#if ! (defined(BLANK) || defined(VERTICO) || defined(ICON))
       nseg = gp_nseg
       ALLOCATE(start(nseg,IRANK))
       ALLOCATE(cnt(nseg,IRANK))
       ALLOCATE(meml(nseg,IRANK))
       ALLOCATE(memu(nseg,IRANK))

       start(:,:) = gp_start(:,:)
       cnt(:,:) = gp_cnt(:,:)
       meml(:,:) = gp_meml(:,:)
       memu(:,:) = gp_memu(:,:)

!!$       cnt(:,3) = maxmode
!!$       memu(:,3) = maxmode
       cnt(:,_IN_XYZN_) = maxmode
       memu(:,_IN_XYZN_) = maxmode

       CALL set_representation_decomp(status, REPR_PTRAC_4D &
            , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
       CALL channel_halt(substr, status)
       
       DEALLOCATE(start) ; NULLIFY(start)
       DEALLOCATE(cnt)   ; NULLIFY(cnt)
       DEALLOCATE(meml)  ; NULLIFY(meml)
       DEALLOCATE(memu)  ; NULLIFY(memu)
#endif       
       ! ----------------------------------------------------------
       
       reprname = ''
       WRITE(reprname,'(a14,i3.3,a1,a2)') &
            'REPR_PTRAC_1D_',maxmode,'_',GPTRSTR
       CALL new_representation(status, REPR_PTRAC_1D, TRIM(reprname) &
            , rank = 1, link = 'x---', dctype = DC_BC                &
            , dimension_ids = (/ DIMID_NMODE /)                      &
            , ldimlen       = (/ AUTO   /)                           &
            , axis = 'N---'                                          &
            )
       CALL channel_halt(substr, status)
#if ! (defined(BLANK) || defined(VERTICO) || defined(ICON)) 
       nseg = 1
       ALLOCATE(start(nseg,IRANK))
       ALLOCATE(cnt(nseg,IRANK))
       ALLOCATE(meml(nseg,IRANK))
       ALLOCATE(memu(nseg,IRANK))

       start(:,:) = 1
       cnt(:,:) = 1
       meml(:,:) = 1
       memu(:,:) = 1

       start(:,1) = 1
       cnt(:,1)   = maxmode
       meml(:,1)  = 1
       memu(:,1)  = maxmode

       CALL set_representation_decomp(status, REPR_PTRAC_1D &
            , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
       CALL channel_halt(substr, status)
       
       DEALLOCATE(start) ; NULLIFY(start)
       DEALLOCATE(cnt)   ; NULLIFY(cnt)
       DEALLOCATE(meml)  ; NULLIFY(meml)
       DEALLOCATE(memu)  ; NULLIFY(memu)
#endif
       ! ----------------------------------------------------------

       ! NEW OBJECTS
       IF (p_parallel_io) WRITE(*,*) ' ... wetradius'
       CALL new_channel_object(status, modstr//'_'//GPTRSTR &
            , 'wetradius'                     &
            , p4 = radius_ptr                 &
            , reprid = REPR_PTRAC_4D          &
            , lrestreq = .FALSE. )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr//'_'//GPTRSTR, 'wetradius'   &
            , 'long_name', c='const. ambient aerosol radius')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//GPTRSTR, 'wetradius'   &
            , 'units', c='m')
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... densaer'
       CALL new_channel_object(status, modstr//'_'//GPTRSTR &
            , 'densaer'                       &
            , p4 = density_ptr                &
            , reprid = REPR_PTRAC_4D          &
            , lrestreq = .FALSE. )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr//'_'//GPTRSTR, 'densaer'   &
            , 'long_name', c='aerosol density')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//GPTRSTR, 'densaer'   &
            , 'units', c='kg/m3')
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... sigma'
       CALL new_channel_object(status, modstr//'_'//GPTRSTR &
            , 'sigma'                         &
            , p1 = sigma_ptr                  &
            , reprid = REPR_PTRAC_1D          &
            , lrestreq = .FALSE. )
       sigma_ptr(:) = 1.00_dp
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr//'_'//GPTRSTR, 'sigma'         &
            , 'long_name', c='standard deviation of aerosol mode')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//GPTRSTR, 'sigma'   &
            , 'units', c=' ')
       CALL channel_halt(substr, status)
       
       ! INITIALIZE RADIUS AND SIGMA
       ! -> FIRST TRACER OF MODE DETERMINES RADIUS, SIGMA, and DENSITY
       mode_loop: DO jm = 1, maxmode
          DO i=1, NMAXTRAC
             IF (TRIM(TRAC(i)%trset)    == '') CYCLE
             IF (TRIM(TRAC(i)%trlist)   == '') CYCLE
             m = SIZE(XTRAC(i)%trs)
             n = SIZE(XTRAC(i)%trnb)
             DO k=1, m
                IF (TRIM(XTRAC(i)%trs(k)) /= GPTRSTR) CYCLE
                DO j=1, n
                   CALL get_tracer(status, GPTRSTR, TRIM(XTRAC(i)%trnb(j)) &
                        , subname = TRIM(XTRAC(i)%trns(j)), idx=jt &
                        , submodel=sm, medium=med)
                   CALL tracer_halt(substr,status)
                   IF (TRIM(sm) /= modstr) CYCLE
                   IF (med /= AEROSOL) CYCLE
                   CALL get_tracer(status, GPTRSTR, jt, I_AEROSOL_MODE, i=mode)
                   CALL tracer_halt(substr,status)
                   IF (mode == jm) THEN
                      radius_ptr(_RI_XYZN_(:,:,:,jm))  = TRAC(i)%radi
                      sigma_ptr(jm)         = TRAC(i)%sigma
                      density_ptr(_RI_XYZN_(:,:,:,jm)) = (2._dp/3._dp) * TRAC(i)%dens
                      EXIT
                   END IF
                END DO
             END DO
          END DO
       END DO mode_loop

    END IF aerosol_cond

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)
#endif

  END SUBROUTINE ptrac_init_memory_gp
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE ptrac_init_memory_lg

#if defined(ECHAM5)
#ifdef MESSYCHANNELS
    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_pe

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DC_IX, DC_BC
    USE messy_main_tracer_mem_bi,    ONLY: LGTRSTR
    USE messy_main_transform_bi,     ONLY: get_dc_index
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt

    ! SMCL
    USE messy_main_channel_dimensions,  ONLY: new_dimension, get_dimension &
                                            , t_dimension
    USE messy_main_channel_repr,        ONLY: new_representation, AUTO &
                                            , set_representation_decomp  &
                                            , IRANK, PIOTYPE_COL
    USE messy_main_channel,             ONLY: new_channel, new_channel_object &
                                            , new_attribute
    USE messy_main_tracer,              ONLY: get_tracer, AEROSOL &
                                            , I_AEROSOL_MODE

    IMPLICIT NONE
    INTRINSIC :: TRIM, SIZE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'ptrac_init_memory_lg'
    LOGICAL                          :: lchannel
    !
    ! CHANNEL OBJECTS FOR AEROSOL
    REAL(DP), DIMENSION(:,:), POINTER :: radius_ptr => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: density_ptr => NULL()
    REAL(DP), DIMENSION(:),   POINTER :: sigma_ptr => NULL()
    !
    ! CHANNEL MANAGEMENT
    TYPE(t_dimension), POINTER       :: dim => NULL()
    ! INDEX BOUNDARIES NGCELL
    INTEGER, DIMENSION(:,:),   POINTER :: IDX => NULL()
    INTEGER :: NGCELL, NCELL ! NUMBER OF CELLS
    INTEGER :: DIMID_NGCELL  ! DIMENSION ID: NUMBER OF CELLS
    INTEGER :: DIMID_NMODE   ! DIMENSION ID: NUMBER OF (PSEUDO-) AEROSOL MODES
    INTEGER :: REPR_PTRAC_2D ! REPRESENTATION ID
    INTEGER :: REPR_PTRAC_1D ! REPRESENTATION ID
    !
    INTEGER                          :: status
    INTEGER                          :: i, m, n, j, k, jt
    CHARACTER(LEN=STRLEN_MEDIUM)     :: sm = ''
    INTEGER                          :: med
    INTEGER                          :: mode
    !
    ! NUMBER OF MODES
    INTEGER                      :: maxmode, jm
    CHARACTER(LEN=STRLEN_MEDIUM) :: reprname
    !
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    ! LG ACTIVE
    CALL get_dimension(status, 'NCELL', dim)
    IF (status /= 0) RETURN
    !CALL channel_halt(substr, status)

    ! INIT DIMENSION INFO
    DIMID_NGCELL = dim%id
    NGCELL       = dim%len
    ! SET THE NUMBER OF CELLS PER PE
    CALL get_dc_index(NGCELL, IDX)
    NCELL = IDX(p_pe,2) - IDX(p_pe,1) + 1

    ! CHANNEL REQUIRED ?
    lchannel = .FALSE.
    maxmode = 0
    DO i=1, NMAXTRAC
       IF (TRIM(TRAC(i)%trset)    == '') CYCLE
       IF (TRIM(TRAC(i)%trlist)   == '') CYCLE
       m = SIZE(XTRAC(i)%trs)
       n = SIZE(XTRAC(i)%trnb)
       DO k=1, m
          IF (TRIM(XTRAC(i)%trs(k)) /= LGTRSTR) CYCLE
          DO j=1, n
             CALL get_tracer(status, LGTRSTR, TRIM(XTRAC(i)%trnb(j)) &
                  , subname = TRIM(XTRAC(i)%trns(j)), idx=jt &
                  , submodel=sm, medium=med)
             CALL tracer_halt(substr,status)
             IF ( (TRIM(sm) == modstr) .AND. (med == AEROSOL) ) THEN
                lchannel = .TRUE.
                CALL get_tracer(status, LGTRSTR, jt, I_AEROSOL_MODE, i=mode)
                CALL tracer_halt(substr,status)
                IF (mode > maxmode) maxmode = mode
             END IF
          END DO
       END DO
    END DO

    IF (.NOT. lchannel) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL new_channel(status, modstr//'_'//LGTRSTR)
    CALL channel_halt(substr, status)

    aerosol_cond: IF (maxmode > 0) THEN

       ! ADDITIONAL DIMENSION
       CALL new_dimension(status, DIMID_NMODE &
            , 'PTRAC_NMODE'//'_'//LGTRSTR, maxmode)
       CALL channel_halt(substr, status)

       ! NEW REPRESENTATIONS
       reprname = ''
       WRITE(reprname,'(a14,i3.3,a1,a2)') &
            'REPR_PTRAC_2D_',maxmode,'_',LGTRSTR
       CALL new_representation(status, REPR_PTRAC_2D, TRIM(reprname)  &
            , rank = 2, link = 'xx--', dctype = DC_IX                 &
            , dimension_ids = (/ DIMID_NGCELL, DIMID_NMODE /)         &
            , ldimlen       = (/ NCELL, AUTO /)                       &
            )
       CALL channel_halt(substr, status)

       nseg = 1
       ALLOCATE(start(nseg,IRANK))
       ALLOCATE(cnt(nseg,IRANK))
       ALLOCATE(meml(nseg,IRANK))
       ALLOCATE(memu(nseg,IRANK))

       start(:,:) = 1
       cnt(:,:) = 1
       meml(:,:) = 1
       memu(:,:) = 1

       start(:,1) = IDX(p_pe,1)
       cnt(:,1)   = NCELL
       meml(:,1)  = 1
       memu(:,1)  = NCELL

       start(:,2) = 1
       cnt(:,2)   = maxmode
       meml(:,2)  = 1
       memu(:,2)  = maxmode

       CALL set_representation_decomp(status, REPR_PTRAC_2D &
            , start, cnt, memu, meml, piotype=PIOTYPE_COL)
       CALL channel_halt(substr, status)

       DEALLOCATE(start) ; NULLIFY(start)
       DEALLOCATE(cnt)   ; NULLIFY(cnt)
       DEALLOCATE(meml)  ; NULLIFY(meml)
       DEALLOCATE(memu)  ; NULLIFY(memu)

       ! ----------------------------------------------------------
       
       reprname = ''
       WRITE(reprname,'(a14,i3.3,a1,a2)') &
            'REPR_PTRAC_1D_',maxmode,'_',LGTRSTR
       CALL new_representation(status, REPR_PTRAC_1D, TRIM(reprname) &
            , rank = 1, link = 'x---', dctype = DC_BC                &
            , dimension_ids = (/ DIMID_NMODE /)                      &
            , ldimlen       = (/ AUTO   /)                           &
            )
       CALL channel_halt(substr, status)

       nseg = 1
       ALLOCATE(start(nseg,IRANK))
       ALLOCATE(cnt(nseg,IRANK))
       ALLOCATE(meml(nseg,IRANK))
       ALLOCATE(memu(nseg,IRANK))

       start(:,:) = 1
       cnt(:,:) = 1
       meml(:,:) = 1
       memu(:,:) = 1

       start(:,1) = 1
       cnt(:,1)   = maxmode
       meml(:,1)  = 1
       memu(:,1)  = maxmode

       CALL set_representation_decomp(status, REPR_PTRAC_1D &
            , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
       CALL channel_halt(substr, status)
       
       DEALLOCATE(start) ; NULLIFY(start)
       DEALLOCATE(cnt)   ; NULLIFY(cnt)
       DEALLOCATE(meml)  ; NULLIFY(meml)
       DEALLOCATE(memu)  ; NULLIFY(memu)

       ! ----------------------------------------------------------

       ! NEW OBJECTS
       IF (p_parallel_io) WRITE(*,*) ' ... wetradius'
       CALL new_channel_object(status, modstr//'_'//LGTRSTR &
            , 'wetradius'                     &
            , p2 = radius_ptr                 &
            , reprid = REPR_PTRAC_2D          &
            , lrestreq = .FALSE. )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr//'_'//LGTRSTR, 'wetradius'   &
            , 'long_name', c='const. ambient aerosol radius')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//LGTRSTR, 'wetradius'   &
            , 'units', c='m')
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... densaer'
       CALL new_channel_object(status, modstr//'_'//LGTRSTR &
            , 'densaer'                       &
            , p2 = density_ptr                &
            , reprid = REPR_PTRAC_2D          &
            , lrestreq = .FALSE. )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr//'_'//LGTRSTR, 'densaer'   &
            , 'long_name', c='aerosol density')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//LGTRSTR, 'densaer'   &
            , 'units', c='kg/m3')
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... sigma'
       CALL new_channel_object(status, modstr//'_'//LGTRSTR &
            , 'sigma'                         &
            , p1 = sigma_ptr                  &
            , reprid = REPR_PTRAC_1D          &
            , lrestreq = .FALSE. )
       CALL channel_halt(substr, status)
       sigma_ptr(:) = 1.00_dp
       !
       CALL new_attribute(status, modstr//'_'//LGTRSTR, 'sigma'         &
            , 'long_name', c='standard deviation of aerosol mode')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//LGTRSTR, 'sigma'   &
            , 'units', c=' ')
       CALL channel_halt(substr, status)
       
       ! INITIALIZE RADIUS AND SIGMA
       ! -> FIRST TRACER OF MODE DETERMINES RADIUS, SIGMA, and DENSITY
       mode_loop: DO jm = 1, maxmode
          DO i=1, NMAXTRAC
             IF (TRIM(TRAC(i)%trset)    == '') CYCLE
             IF (TRIM(TRAC(i)%trlist)   == '') CYCLE
             m = SIZE(XTRAC(i)%trs)
             n = SIZE(XTRAC(i)%trnb)
             DO k=1, m
                IF (TRIM(XTRAC(i)%trs(k)) /= LGTRSTR) CYCLE
                DO j=1, n
                   CALL get_tracer(status, LGTRSTR, TRIM(XTRAC(i)%trnb(j)) &
                        , subname = TRIM(XTRAC(i)%trns(j)), idx=jt &
                        , submodel=sm, medium=med)
                   CALL tracer_halt(substr,status)
                   IF (TRIM(sm) /= modstr) CYCLE
                   IF (med /= AEROSOL) CYCLE
                   CALL get_tracer(status, LGTRSTR, jt, I_AEROSOL_MODE, i=mode)
                   CALL tracer_halt(substr,status)
                   IF (mode == jm) THEN
                      radius_ptr(:,jm)  = TRAC(i)%radi
                      sigma_ptr(jm)     = TRAC(i)%sigma
                      density_ptr(:,jm) = (2._dp/3._dp) * TRAC(i)%dens
                      EXIT
                   END IF
                END DO
             END DO
          END DO
       END DO mode_loop

    END IF aerosol_cond

    ! CLEAN
    DEALLOCATE(IDX); NULLIFY(IDX)

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)
#endif
#endif
  END SUBROUTINE ptrac_init_memory_lg
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE ptrac_init_memory_cl

#if defined(ECHAM5)
#ifdef MESSYCHANNELS
    ! BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_pe

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DC_IX, DC_BC
    USE messy_main_tracer_mem_bi,    ONLY: CLTRSTR
    USE messy_main_transform_bi,     ONLY: get_dc_index
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt

    ! SMCL
    USE messy_main_channel_dimensions,  ONLY: new_dimension, get_dimension &
                                            , t_dimension
    USE messy_main_channel_repr,        ONLY: new_representation, AUTO &
                                            , set_representation_decomp  &
                                            , IRANK, PIOTYPE_COL
    USE messy_main_channel,             ONLY: new_channel, new_channel_object &
                                            , new_attribute
    USE messy_main_tracer,              ONLY: get_tracer, AEROSOL &
                                            , I_AEROSOL_MODE

    IMPLICIT NONE
    INTRINSIC :: TRIM, SIZE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'ptrac_init_memory_cl'
    LOGICAL                          :: lchannel
    !
    ! CHANNEL OBJECTS FOR AEROSOL
    REAL(DP), DIMENSION(:,:), POINTER :: radius_ptr => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: density_ptr => NULL()
    REAL(DP), DIMENSION(:),   POINTER :: sigma_ptr => NULL()
    !
    ! CHANNEL MANAGEMENT
    TYPE(t_dimension), POINTER       :: dim => NULL()
    ! INDEX BOUNDARIES NGCELL
    INTEGER, DIMENSION(:,:),   POINTER :: IDX => NULL()
    INTEGER :: dnparts_max, dnparts ! NUMBER OF CELLS
    INTEGER :: DIMID_dnparts_max  ! DIMENSION ID: NUMBER OF CELLS
    INTEGER :: DIMID_NMODE   ! DIMENSION ID: NUMBER OF (PSEUDO-) AEROSOL MODES
    INTEGER :: REPR_PTRAC_2D ! REPRESENTATION ID
    INTEGER :: REPR_PTRAC_1D ! REPRESENTATION ID
    !
    INTEGER                          :: status
    INTEGER                          :: i, m, n, j, k, jt
    CHARACTER(LEN=STRLEN_MEDIUM)     :: sm = ''
    INTEGER                          :: med
    INTEGER                          :: mode
    !
    ! NUMBER OF MODES
    INTEGER                      :: maxmode, jm
    CHARACTER(LEN=STRLEN_MEDIUM) :: reprname
    !
    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    ! LG ACTIVE
    CALL get_dimension(status, 'dnparts_max', dim)
    IF (status /= 0) RETURN
    CALL channel_halt(substr, status)

    ! INIT DIMENSION INFO
    DIMID_dnparts_max = dim%id
    dnparts_max       = dim%len
    ! SET THE NUMBER OF CELLS PER PE
    CALL get_dc_index(dnparts_max, IDX)
    dnparts = IDX(p_pe,2) - IDX(p_pe,1) + 1

    ! CHANNEL REQUIRED ?
    lchannel = .FALSE.
    maxmode = 0
    DO i=1, NMAXTRAC
       IF (TRIM(TRAC(i)%trset)    == '') CYCLE
       IF (TRIM(TRAC(i)%trlist)   == '') CYCLE
       m = SIZE(XTRAC(i)%trs)
       n = SIZE(XTRAC(i)%trnb)
       DO k=1, m
          IF (TRIM(XTRAC(i)%trs(k)) /= CLTRSTR) CYCLE
          DO j=1, n
             CALL get_tracer(status, CLTRSTR, TRIM(XTRAC(i)%trnb(j)) &
                  , subname = TRIM(XTRAC(i)%trns(j)), idx=jt &
                  , submodel=sm, medium=med)
             CALL tracer_halt(substr,status)
             IF ( (TRIM(sm) == modstr) .AND. (med == AEROSOL) ) THEN
                lchannel = .TRUE.
                CALL get_tracer(status, CLTRSTR, jt, I_AEROSOL_MODE, i=mode)
                CALL tracer_halt(substr,status)
                IF (mode > maxmode) maxmode = mode
             END IF
          END DO
       END DO
    END DO

    IF (.NOT. lchannel) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    CALL new_channel(status, modstr//'_'//CLTRSTR)
    CALL channel_halt(substr, status)

    aerosol_cond: IF (maxmode > 0) THEN

       ! ADDITIONAL DIMENSION
       CALL new_dimension(status, DIMID_NMODE &
            , 'PTRAC_NMODE'//'_'//CLTRSTR, maxmode)
       CALL channel_halt(substr, status)

       ! NEW REPRESENTATIONS
       reprname = ''
       WRITE(reprname,'(a14,i3.3,a1,a2)') &
            'REPR_PTRAC_2D_',maxmode,'_',CLTRSTR
       CALL new_representation(status, REPR_PTRAC_2D, TRIM(reprname)  &
            , rank = 2, link = 'xx--', dctype = DC_IX                 &
            , dimension_ids = (/ DIMID_dnparts_max, DIMID_NMODE /)         &
            , ldimlen       = (/ dnparts, AUTO /)                       &
            )
       CALL channel_halt(substr, status)

       nseg = 1
       ALLOCATE(start(nseg,IRANK))
       ALLOCATE(cnt(nseg,IRANK))
       ALLOCATE(meml(nseg,IRANK))
       ALLOCATE(memu(nseg,IRANK))

       start(:,:) = 1
       cnt(:,:) = 1
       meml(:,:) = 1
       memu(:,:) = 1

       start(:,1) = IDX(p_pe,1)
       cnt(:,1)   = dnparts
       meml(:,1)  = 1
       memu(:,1)  = dnparts

       start(:,2) = 1
       cnt(:,2)   = maxmode
       meml(:,2)  = 1
       memu(:,2)  = maxmode

       CALL set_representation_decomp(status, REPR_PTRAC_2D &
            , start, cnt, memu, meml, piotype=PIOTYPE_COL)
       CALL channel_halt(substr, status)

       DEALLOCATE(start) ; NULLIFY(start)
       DEALLOCATE(cnt)   ; NULLIFY(cnt)
       DEALLOCATE(meml)  ; NULLIFY(meml)
       DEALLOCATE(memu)  ; NULLIFY(memu)

       ! ----------------------------------------------------------
       
       reprname = ''
       WRITE(reprname,'(a14,i3.3,a1,a2)') &
            'REPR_PTRAC_1D_',maxmode,'_',CLTRSTR
       CALL new_representation(status, REPR_PTRAC_1D, TRIM(reprname) &
            , rank = 1, link = 'x---', dctype = DC_BC                &
            , dimension_ids = (/ DIMID_NMODE /)                      &
            , ldimlen       = (/ AUTO   /)                           &
            )
       CALL channel_halt(substr, status)

       nseg = 1
       ALLOCATE(start(nseg,IRANK))
       ALLOCATE(cnt(nseg,IRANK))
       ALLOCATE(meml(nseg,IRANK))
       ALLOCATE(memu(nseg,IRANK))

       start(:,:) = 1
       cnt(:,:) = 1
       meml(:,:) = 1
       memu(:,:) = 1

       start(:,1) = 1
       cnt(:,1)   = maxmode
       meml(:,1)  = 1
       memu(:,1)  = maxmode

       CALL set_representation_decomp(status, REPR_PTRAC_1D &
            , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
       CALL channel_halt(substr, status)
       
       DEALLOCATE(start) ; NULLIFY(start)
       DEALLOCATE(cnt)   ; NULLIFY(cnt)
       DEALLOCATE(meml)  ; NULLIFY(meml)
       DEALLOCATE(memu)  ; NULLIFY(memu)

       ! ----------------------------------------------------------

       ! NEW OBJECTS
       IF (p_parallel_io) WRITE(*,*) ' ... wetradius'
       CALL new_channel_object(status, modstr//'_'//CLTRSTR &
            , 'wetradius'                     &
            , p2 = radius_ptr                 &
            , reprid = REPR_PTRAC_2D          &
            , lrestreq = .FALSE. )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr//'_'//CLTRSTR, 'wetradius'   &
            , 'long_name', c='const. ambient aerosol radius')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//CLTRSTR, 'wetradius'   &
            , 'units', c='m')
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... densaer'
       CALL new_channel_object(status, modstr//'_'//CLTRSTR &
            , 'densaer'                       &
            , p2 = density_ptr                &
            , reprid = REPR_PTRAC_2D          &
            , lrestreq = .FALSE. )
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr//'_'//CLTRSTR, 'densaer'   &
            , 'long_name', c='aerosol density')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//CLTRSTR, 'densaer'   &
            , 'units', c='kg/m3')
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... sigma'
       CALL new_channel_object(status, modstr//'_'//CLTRSTR &
            , 'sigma'                         &
            , p1 = sigma_ptr                  &
            , reprid = REPR_PTRAC_1D          &
            , lrestreq = .FALSE. )
       CALL channel_halt(substr, status)
       sigma_ptr(:) = 1.00_dp ! op_pj_20171123
       !
       CALL new_attribute(status, modstr//'_'//CLTRSTR, 'sigma'         &
            , 'long_name', c='standard deviation of aerosol mode')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_'//CLTRSTR, 'sigma'   &
            , 'units', c=' ')
       CALL channel_halt(substr, status)
       
       ! INITIALIZE RADIUS AND SIGMA
       ! -> FIRST TRACER OF MODE DETERMINES RADIUS, SIGMA, and DENSITY
       mode_loop: DO jm = 1, maxmode
          DO i=1, NMAXTRAC
             IF (TRIM(TRAC(i)%trset)    == '') CYCLE
             IF (TRIM(TRAC(i)%trlist)   == '') CYCLE
             m = SIZE(XTRAC(i)%trs)
             n = SIZE(XTRAC(i)%trnb)
             DO k=1, m
                IF (TRIM(XTRAC(i)%trs(k)) /= CLTRSTR) CYCLE
                DO j=1, n
                   CALL get_tracer(status, CLTRSTR, TRIM(XTRAC(i)%trnb(j)) &
                        , subname = TRIM(XTRAC(i)%trns(j)), idx=jt &
                        , submodel=sm, medium=med)
                   CALL tracer_halt(substr,status)
                   IF (TRIM(sm) /= modstr) CYCLE
                   IF (med /= AEROSOL) CYCLE
                   CALL get_tracer(status, CLTRSTR, jt, I_AEROSOL_MODE, i=mode)
                   CALL tracer_halt(substr,status)
                   IF (mode == jm) THEN
                      radius_ptr(:,jm)  = TRAC(i)%radi
                      sigma_ptr(jm)     = TRAC(i)%sigma
                      density_ptr(:,jm) = (2._dp/3._dp) * TRAC(i)%dens
                      EXIT
                   END IF
                END DO
             END DO
          END DO
       END DO mode_loop

    END IF aerosol_cond

    ! CLEAN
    DEALLOCATE(IDX); NULLIFY(IDX)

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)
#endif
#endif
  END SUBROUTINE ptrac_init_memory_cl
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
SUBROUTINE ptrac_vdiff

!!$#ifdef ICON
!!$  USE messy_main_tracer_mem_bi, ONLY: xtte
!!$  USE messy_main_data_bi,       ONLY: jrow
!!$
!!$  IMPLICIT NONE
!!$
!!$  INTEGER :: znlev
!!$
!!$  IF (jrow == 1) THEN
!!$     znlev = SIZE(xtte,2)
!!$
!!$     xtte(:,znlev,jrow,6) = 1.0E-12_dp
!!$  END IF
!!$#endif

  END SUBROUTINE ptrac_vdiff
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE ptrac_read_nml_cpl(status, iou)

    ! SMCL
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'ptrac_read_nml_cpl'

    NAMELIST /CPL/ TRAC, TPROP

    ! LOCAL
    LOGICAL :: lex      ! file exists ?
    INTEGER :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE ptrac_read_nml_cpl
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_ptrac_si
! **********************************************************************
