! **********************************************************************
PROGRAM TRACER
! **********************************************************************

  ! BMIL
  USE messy_main_grid_def_mem_bi, ONLY: ngpblks, jrow
  USE messy_main_tracer_tools_bi, ONLY: tracer_halt
  USE messy_main_tracer_bi

  ! SMIL
  USE messy_ptrac_si
  USE messy_othersm_si

  IMPLICIT NONE

  INTEGER, PARAMETER :: n_time_steps = 3
  INTEGER            :: nt

  ! ====================================================================
  ! INITIALISATION PHASE
  ! ====================================================================
  
  ! INITIALIZE ---------------------------------------------------------
  CALL bml_initialize                 ! BML: initialize base model
  CALL main_tracer_initialize         ! BMIL: read namelist(s)
  CALL ptrac_initialize               ! SMIL: initalize submodel
  CALL othersm_initialize             ! SMIL: initalize submodel
  !# ADD OTHER SUBMODELS HERE

  ! NEW_TRACER ---------------------------------------------------------
  CALL main_tracer_new_tracer(1)      ! BMIL: define new tracer set(s)
  CALL ptrac_new_tracer               ! SMIL: add tracers to set(s)
  !# ADD OTHER SUBMODELS HERE
  CALL main_tracer_new_tracer(2)      ! BMIL: define tracer families
  CALL main_tracer_new_tracer(3)      ! BMIL: diagnostic output

  ! INIT_MEMORY --------------------------------------------------------
  CALL main_tracer_init_memory(1)     ! BMIL: setup tracer memory
  CALL ptrac_init_memory              ! SMIL: setup submodel memory
  !# ADD OTHER SUBMODELS HERE
  !
  ! a. ASSOCIATE TRACER MEMORY WITH MESSY DATA TRANSFER / EXPORT INTERFACE
  !    (NOTE: THIS IS NOT IMPLEMENTED IN THIS SIMPLE BOX MODEL)
  ! b. SET META INFORMATION OF FAMILY MEMBERS TO 'FRACTION'
  !    (NOTE: THIS IS FOR A PROPER ADVECTION INITIALISATION ...
  !           AND NOT REALLY REQUIRED IN THIS SIMPLE BOX MODEL)
  CALL main_tracer_init_memory(2)     ! BMIL: ...

  ! INIT_COUPLING -----------------------------------------------------
  ! RESET META INFORMATION OF FAMILY MEMBERS TO 'TRACER'
  ! (NOTE: THIS IS USED AFTER A PROPER ADVECTION INITIALISATION ...
  !           AND NOT REALLY REQUIRED IN THIS BOX MODEL)
  CALL main_tracer_init_coupling      ! BMIL: ...
  CALL othersm_init_coupling          ! SMIL: coupling between SMs
  !# ADD OTHER SUBMODELS HERE
  !

  ! INIT_TRACER -------------------------------------------------------
  ! TRACER INITIALISATION
  ! a. CHECK TRACER INITIALISATION FROM RESTART FILES
  !    (NOTE: THIS IS NOT IMPLEMENTED IN THIS SIMPLE BOX MODEL)
  CALL main_tracer_init_tracer(1)     ! BMIL: ...
  CALL othersm_init_tracer            ! SMIL: initialize tracers
  !# ADD OTHER SUBMODELS HERE
  !
  ! b. CHECK TRACER INITIALISATION FROM SUBMODELS AND SET
  !    TO CONSTANT IF NOT YET INITIALIZED
  CALL main_tracer_init_tracer(2)     ! BMIL: ...
  ! c. DIAGNOSTIC OUTPUT OF TRACER INITIALISATION
  CALL main_tracer_init_tracer(3)     ! BMIL: diagn. out. of tracer init.

  ! ====================================================================
  ! TIME INTEGRATION PHASE
  ! ====================================================================
  time_loop: DO nt=1, n_time_steps

     CALL main_tracer_global_start(1)    ! BMIL
     CALL othersm_global_start           ! SMIL: (here 'diffusion')
     !# ADD OTHER SUBMODELS HERE
     !# NOTE: IN THE SUBMODEL SMIL, YOU CAN APPLY THE
     !#       CONVERSION ROUTINE SEQUENCE
     !#            CALL main_tracer_fconv_glb('sum',substr,S1TRSTR) ! TYPE-2
     !#            CALL main_tracer_fconv_glb('t2f',substr,S1TRSTR) ! TYPE-1
     !#            ... access families ...
     !#            CALL main_tracer_fconv_glb('f2t',substr,S1TRSTR) ! TYPE-1
     !#            CALL main_tracer_fconv_glb('rsc',substr,S1TRSTR) ! TYPE-2
     !
     CALL main_tracer_global_start(2)    ! BMIL

     CALL main_tracer_beforeadv           ! BML

     !# ADD ADVECTION HERE

     CALL main_tracer_afteradv           ! BML

     dim3_loop: DO jrow=1, ngpblks

        ! ### CALL BASE MODEL PHYSICS HERE
        ! e.g., change pressi_3d, gboxarea_2d, ...

        CALL main_tracer_local_start

        CALL othersm_physc           ! SMIL: (here: 'decay')
        !# ADD OTHER SUBMODELS HERE, WHICH CHANGE TRACERS
        !# NOTE: IN THE SUBMODEL SMIL, YOU CAN APPLY THE
        !#       CONVERSION ROUTINE SEQUENCE
        !#            CALL main_tracer_fconv_loc('sum',substr,S1TRSTR) ! TYPE-2
        !#            CALL main_tracer_fconv_loc('t2f',substr,S1TRSTR) ! TYPE-1
        !#            ... access families ...
        !#            CALL main_tracer_fconv_loc('f2t',substr,S1TRSTR) ! TYPE-1
        !#            CALL main_tracer_fconv_loc('rsc',substr,S1TRSTR) ! TYPE-2

     END DO dim3_loop

     CALL main_tracer_global_end
     !# ADD OTHER SUBMODELS HERE

     CALL bml_output               ! BML: output routine(s)

  END DO time_loop

  !CALL test_chemprop ! mz_rs_20160521

  ! ====================================================================
  ! FINALIZING PHASE
  ! ====================================================================
  CALL main_tracer_free_memory

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE bml_initialize

    USE messy_main_constants_mem,   ONLY: dp, g
    USE messy_main_data_bi,         ONLY: pressi_3d
    USE messy_main_grid_def_bi,     ONLY: gboxarea_2d, grmass
    USE messy_main_grid_def_mem_bi, ONLY: nlev

    IMPLICIT NONE

    ! surface area of grid boxes
    gboxarea_2d(:,:) = 1.0E6_dp       ! 1 km^2

    ! pressure
    pressi_3d(:,1,:) =   1000.0_dp ! Pa
    pressi_3d(:,2,:) =  50000.0_dp ! Pa
    pressi_3d(:,3,:) = 100000.0_dp ! Pa

    ! grid box airmass
    grmass(:,:,:) = ((pressi_3d(:,2:nlev+1,:) - pressi_3d(:,1:nlev,:)) / g) &
         * SPREAD(gboxarea_2d,2,nlev)

  END SUBROUTINE bml_initialize
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE bml_output

    IMPLICIT NONE

    WRITE(*,*) 'END OF TIME STEP : ',nt

  END SUBROUTINE bml_output
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  ! mz_rs_20160521+
  SUBROUTINE test_chemprop
    USE messy_main_tracer
    USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM

    INTEGER :: specnum
    INTEGER :: status
    INTEGER :: mydata_int
    REAL(DP) :: mydata_real
    CHARACTER(LEN=STRLEN_MEDIUM) :: mydata_char
    INTEGER :: idx ! op_pj_20160823
    CHARACTER(LEN=*), PARAMETER :: substr = 'test_chemprop'

    PRINT *
    PRINT *, "TEST_CHEMPROP"
    PRINT *, "Use chemprop directly, e.g. for loop:"
    DO specnum = 5,8
      WRITE(*,'(A,A10,A,I3,A,F10.4,A,A)') &
        " kppname: ",     TRIM(chemprop(specnum)%kppname), &
        " I_charge: ",    chemprop(specnum)%cask_i(I_charge), &
        " R_MOLARMASS: ", chemprop(specnum)%cask_r(R_MOLARMASS), &
        " S_casrn: ",     chemprop(specnum)%cask_s(S_casrn)
    ENDDO

    PRINT *
    PRINT *, "Use get_chemprop:"
    status = get_chemprop("NH3", I_charge, mydata_int)
    IF (status /=0) STOP 'ERROR I_charge'
    status = get_chemprop("NH3", R_MOLARMASS, mydata_real)
    IF (status /=0) STOP 'ERROR R_MOLARMASS'
    status = get_chemprop("NH3", S_casrn, mydata_char)
    IF (status /=0) STOP 'ERROR S_casrn'
    WRITE(*,'(A,A10,A,I3,A,F10.4,A,A)') &
      " kppname: ",     "NH3", &
      " I_charge: ",    mydata_int, &
      " R_MOLARMASS: ", mydata_real, &
      " S_casrn: ",     mydata_char

    ! op_pj_20160823+
    PRINT *
    PRINT *, "Use get_tracer:"
    CALL get_tracer(status, 's1', 'NH3', idx=idx)
    CALL tracer_halt(substr,status)
    CALL get_tracer(status, 's1', idx, I_charge, mydata_int)
    CALL tracer_halt(substr,status)
    CALL get_tracer(status, 's1', idx, R_molarmass, mydata_real)
    CALL tracer_halt(substr,status)
    CALL get_tracer(status, 's1', idx, S_casrn, mydata_char)
    CALL tracer_halt(substr,status)
    WRITE(*,'(A,A10,A,I3,A,F10.4,A,A)') &
      " TRACER:  ",     "NH3", &
      " I_charge: ",    mydata_int, &
      " R_MOLARMASS: ", mydata_real, &
      " S_casrn: ",     mydata_char
    ! op_pj_20160823-

  END SUBROUTINE test_chemprop
  ! mz_rs_20160521-
  ! --------------------------------------------------------------------

! **********************************************************************
END PROGRAM TRACER
! **********************************************************************
