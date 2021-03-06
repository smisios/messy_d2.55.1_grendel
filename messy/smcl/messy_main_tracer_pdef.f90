! ************************************************************************
MODULE messy_main_tracer_pdef
! ************************************************************************

  ! MODULE FOR FORCING POSITIVE TRACERS (MESSy-SMCL)
  !
  ! Authors:
  !    Patrick Joeckel, MPICH, March 2007
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM
  USE messy_main_tracer,        ONLY: NMAXSETID, NSETID, TRSET, STRLEN_TRSET &
                                    , modstr, TRRANK

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: NULL

  PUBLIC :: DP

  PUBLIC :: tracpdef_initmem
  PUBLIC :: tracpdef_settings
  PUBLIC :: tracpdef_airmass
  PUBLIC :: tracpdef_mask
  PUBLIC :: tracpdef_integrate
  PUBLIC :: tracpdef_freemem
  PUBLIC :: tracpdef_print
  !
  PUBLIC :: tracer_pdef_read_nml_ctrl

  INTERFACE tracpdef_airmass
     MODULE PROCEDURE tracpdef_airmass_const
     MODULE PROCEDURE tracpdef_airmass_3d
  END INTERFACE

  ! ----------- <

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodstr = 'tracer_pdef'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: submodver = '1.3'

  TYPE IO_TRACPDEF
     CHARACTER(len=STRLEN_TRSET)  :: set=''      ! name of tracer set
     CHARACTER(len=STRLEN_MEDIUM) :: name=''     ! name of tracer family
     CHARACTER(len=STRLEN_MEDIUM) :: subname=''  ! subname of tracer family
     LOGICAL, DIMENSION(2)        :: lswitch = (/ .FALSE., .FALSE. /)
     REAL(DP)                     :: rtol = 1.0_dp
  END TYPE IO_TRACPDEF
  PUBLIC :: IO_TRACPDEF

  ! MAX NUMBER OF NAMELIST ENTRIES
  INTEGER, PARAMETER, PUBLIC :: NMAXTPDEF = 500 * NMAXSETID  ! 500 per set

  TYPE T_SET_WORKSPACE
     LOGICAL,  DIMENSION(:,:),       POINTER :: xlswitch => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: xrtol    => NULL()
     !
     LOGICAL                                 :: ldiagonly
     !
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: mem_mass_p => NULL()
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: mem_mass_n => NULL()
     REAL(DP), DIMENSION(:,:,:,:,:), POINTER :: mem_mass_c => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: mass_p => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: mass_n => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: mass_c => NULL()
     ! ...
     REAL(DP), DIMENSION(:,:,:),     POINTER :: airmass => NULL()
     ! ...
     REAL(DP), DIMENSION(:,:,:),     POINTER :: mask    => NULL()
     ! ...
     REAL(DP), DIMENSION(:,:,:),     POINTER :: mass_pe => NULL()
     LOGICAL,  DIMENSION(:),         POINTER :: lok     => NULL()
     REAL(DP), DIMENSION(:),         POINTER :: scalf   => NULL()
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: unit => NULL()
  END TYPE T_SET_WORKSPACE
  PUBLIC :: T_SET_WORKSPACE

  ! CTRL-NAMELIST
  LOGICAL,                                 PUBLIC :: L_DIAGOUT = .FALSE.
  TYPE(IO_TRACPDEF), DIMENSION(NMAXSETID), PUBLIC :: TPD_DEFAULT
  TYPE(IO_TRACPDEF), DIMENSION(NMAXTPDEF), PUBLIC :: TPD
  LOGICAL,                                 PUBLIC :: L_TASK_AGGREGATE = .TRUE.

  ! GLOBAL VARIABLES
  TYPE(T_SET_WORKSPACE), DIMENSION(:), POINTER, PUBLIC :: XWRK => NULL()

  ! TRIGGER CALCULATION IN THIS CALL?
  LOGICAL, DIMENSION(:), ALLOCATABLE, PUBLIC :: lnow

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_initmem(nprocs)

    IMPLICIT NONE
    INTRINSIC :: SIZE, ASSOCIATED

    ! I/O
    INTEGER, INTENT(IN) :: nprocs

    ! LOCAL
    INTEGER :: n, s1, s2, sx, ntrac

    ! INIT
    ALLOCATE(XWRK(NSETID))
    ALLOCATE(lnow(NSETID))
    lnow(:) = .TRUE.

    set_loop: DO n=1, NSETID

       ! NUMBER OF TRACERS IN THIS SET
       ntrac = TRSET(n)%ntrac

       ! IF ONE OF TENDENCY OR VALUE OF T-1 IS NOT PRESENT
       ! MASS CAN ONLY BE DIAGNOSED BUT NOT CORRECTED
       XWRK(n)%ldiagonly = (.NOT. ASSOCIATED(TRSET(n)%xtte)) .OR. &
            (.NOT. ASSOCIATED(TRSET(n)%xtm1))

       ALLOCATE(XWRK(n)%xlswitch(ntrac,2))
       XWRK(n)%xlswitch(:,:) = .FALSE.       ! default

       ALLOCATE(XWRK(n)%xrtol(ntrac))
       XWRK(n)%xrtol(:) = 1.0_dp             ! default
       !
       ALLOCATE(XWRK(n)%mem_mass_n(ntrac,1,1,1,1))
       ALLOCATE(XWRK(n)%mem_mass_p(ntrac,1,1,1,1))
       ALLOCATE(XWRK(n)%mem_mass_c(ntrac,1,1,1,1))
       !
       XWRK(n)%mass_n => XWRK(n)%mem_mass_n(:,1,1,1,1)
       XWRK(n)%mass_p => XWRK(n)%mem_mass_p(:,1,1,1,1)
       XWRK(n)%mass_c => XWRK(n)%mem_mass_c(:,1,1,1,1)

       XWRK(n)%mass_n(:) = 0.0_dp
       XWRK(n)%mass_p(:) = 0.0_dp
       XWRK(n)%mass_c(:) = 0.0_dp

       s1 = SIZE(TRSET(n)%xt, 1)
       s2 = SIZE(TRSET(n)%xt, 2)

       SELECT CASE(TRRANK)
       CASE(3)
          !s3 = SIZE(TRSET(n)%xt, 3) ! = ntrac
          sx = SIZE(TRSET(n)%xt, 4)
          !s5 = SIZE(TRSET(n)%xt, 5) ! = 1
       CASE(4)
          sx = SIZE(TRSET(n)%xt, 3)
          !s4 = SIZE(TRSET(n)%xt, 4) ! = ntrac
          !s5 = SIZE(TRSET(n)%xt, 5) ! = 1
       END SELECT
       ALLOCATE(XWRK(n)%airmass(s1, s2, sx))

       XWRK(n)%airmass(:,:,:) = 0.0_dp

       ALLOCATE(XWRK(n)%mask(s1, s2, sx))

       XWRK(n)%mask(:,:,:) = 1.0_dp

       ALLOCATE(XWRK(n)%mass_pe(ntrac, 3, 0:nprocs-1))
       XWRK(n)%mass_pe(:,:,:) = 0.0_DP

       ALLOCATE(XWRK(n)%lok(ntrac))
       XWRK(n)%lok(:) = .TRUE.

       ALLOCATE(XWRK(n)%scalf(ntrac))
       XWRK(n)%scalf(:) = 0.0_DP

       ALLOCATE(XWRK(n)%unit(ntrac))
       XWRK(n)%unit(:) = ''

    END DO set_loop

  END SUBROUTINE tracpdef_initmem
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_settings(ldiagout)

    USE messy_main_tracer, ONLY: get_tracer, STRLEN_TRSET
    USE messy_main_tools,  ONLY: split_name_domain

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    LOGICAL, INTENT(IN) :: ldiagout

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracpdef_settings'
    INTEGER :: i, ierr, jt
    INTEGER :: n, nn
    INTEGER :: ntrac
    CHARACTER(LEN=STRLEN_TRSET) :: zname
    INTEGER :: domain
    INTEGER :: status

    set_loop: DO n=1, NSETID

       ntrac = TRSET(n)%ntrac

       IF (ntrac == 0) CYCLE

       ! INIT: SET DEFAULT
       default_loop: DO nn=1, NMAXSETID
          IF (TRIM(TPD_DEFAULT(nn)%set) == '') CYCLE
          CALL split_name_domain(status, TRIM(TRSET(n)%name) &
               , zname, domain)
          IF ( (TRIM(TPD_DEFAULT(nn)%set) /= TRIM(TRSET(n)%name)) .AND. &
               (TRIM(TPD_DEFAULT(nn)%set) /= TRIM(zname)) ) CYCLE

          XWRK(n)%xlswitch(:,1) = TPD_DEFAULT(nn)%lswitch(1)
          XWRK(n)%xlswitch(:,2) = TPD_DEFAULT(nn)%lswitch(2)
          XWRK(n)%xrtol(:)      = TPD_DEFAULT(nn)%rtol
          IF (ldiagout) THEN
             WRITE(*,*) 'DEFAULT FOR SET '//TRIM(TRSET(n)%name)//': ' &
                  ,TPD_DEFAULT(nn)%lswitch(:), TPD_DEFAULT(nn)%rtol
          END IF
       END DO default_loop

       ! LOOP AND SET SPECIAL
       DO i=1, NMAXTPDEF
          IF (TRIM(TPD(i)%set)  == '') CYCLE
          IF (TRIM(TPD(i)%name) == '') CYCLE
          CALL split_name_domain(status, TRIM(TRSET(n)%name) &
               , zname, domain)
          IF ( (TRIM(TPD(i)%set) /= TRIM(TRSET(n)%name)) .AND. &
               (TRIM(TPD(i)%set) /= TRIM(zname)) ) CYCLE

!!$       CALL get_tracer(ierr, TRIM(TPD(i)%set), TRIM(TPD(i)%name) &
          CALL get_tracer(ierr, TRIM(TRSET(n)%name), TRIM(TPD(i)%name) &
               , subname=TRIM(TPD(i)%subname), idx=jt)
          IF (ierr == 0) THEN
             XWRK(n)%xlswitch(jt,:) = TPD(i)%lswitch(:)
             XWRK(n)%xrtol(jt)      = TPD(i)%rtol
          ELSE
             IF (ldiagout) THEN
                IF (TRIM(TPD(i)%subname) == '') THEN
                   WRITE(*,*) 'TRACER '''//TRIM(TPD(i)%name)//&
!!$                     &''' not found in set '//TRIM(TPD(i)%set)//&
                        &''' not found in set '//TRIM(TRSET(n)%name)//&
                        &' ... skipping ...'
                ELSE
                   WRITE(*,*) 'TRACER '''//TRIM(TPD(i)%name)//&
                        &'_'//TRIM(TPD(i)%subname)//&
!!$                     &''' not found in set '//TRIM(TPD(i)%set)//&
                        &''' not found in set '//TRIM(TRSET(n)%name)//&
                        &' ... skipping ...'
                END IF
             END IF
          END IF
       END DO

       IF (ldiagout) THEN
          WRITE(*,*) '======================================================='
          WRITE(*,*) 'TRACER SET: ',TRIM(TRSET(n)%name)
          WRITE(*,*) '-------------------------------------------------------'
          DO jt=1, ntrac
             WRITE(*,*) jt, TRSET(n)%ti(jt)%tp%ident%fullname &
                  , XWRK(n)%xlswitch(jt,:), XWRK(n)%xrtol(jt)
          END DO
          WRITE(*,*) '======================================================='
       END IF

    END DO set_loop

  END SUBROUTINE tracpdef_settings
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_airmass_3d(setname, airmass)

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*),           INTENT(IN)  :: setname
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: airmass

    ! LOCAL
    LOGICAL :: l_found
    INTEGER :: n

    l_found = .FALSE.
    set_loop: DO n=1, NSETID
       IF (TRIM(TRSET(n)%name) == TRIM(setname)) THEN
          l_found = .TRUE.
          EXIT
       END IF
    END DO set_loop

    IF (.NOT. l_found) RETURN

    IF (TRSET(n)%ntrac == 0) RETURN

    XWRK(n)%airmass(:,:,:) = airmass(:,:,:)

  END SUBROUTINE tracpdef_airmass_3d
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_airmass_const(setname, airmass)

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: setname
    REAL(DP),         INTENT(IN)  :: airmass

    ! LOCAL
    LOGICAL :: l_found
    INTEGER :: n

    l_found = .FALSE.
    set_loop: DO n=1, NSETID
       IF (TRIM(TRSET(n)%name) == TRIM(setname)) THEN
          l_found = .TRUE.
          EXIT
       END IF
    END DO set_loop

    IF (.NOT. l_found) RETURN

    IF (TRSET(n)%ntrac == 0) RETURN

    XWRK(n)%airmass(:,:,:) = airmass

  END SUBROUTINE tracpdef_airmass_const
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_mask(setname, mask, status)

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! I/O
    CHARACTER(LEN=*),           INTENT(IN)  :: setname
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: mask
    INTEGER,                    INTENT(OUT) :: status

    ! LOCAL
    LOGICAL :: l_found
    INTEGER :: n

    status = 1

    l_found = .FALSE.
    set_loop: DO n=1, NSETID
       IF (TRIM(TRSET(n)%name) == TRIM(setname)) THEN
          l_found = .TRUE.
          EXIT
       END IF
    END DO set_loop

    IF (.NOT. l_found) THEN
       status = 2
       RETURN
    END IF

    IF (TRSET(n)%ntrac == 0) THEN
       status = 3
       RETURN
    END IF

    IF (SIZE(XWRK(n)%mask,1) /= SIZE(mask,1) .OR. &
         (SIZE(XWRK(n)%mask,2) /= SIZE(mask,2)).OR. &
         (SIZE(XWRK(n)%mask,3) /= SIZE(mask,3)) ) THEN
       status = 4
    END IF

    XWRK(n)%mask(:,:,:) = mask(:,:,:)

    status = 0

  END SUBROUTINE tracpdef_mask
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_integrate(status, flag, time_step_len, p_pe, lnewtl)

    USE messy_main_constants_mem, ONLY: FLAGGED_BAD, TINY_DP

    IMPLICIT NONE

    ! I/O
    INTEGER,   INTENT(OUT)      :: status
    INTEGER,   INTENT(IN)       :: flag
    REAL(DP),  INTENT(IN)       :: time_step_len
    INTEGER,   INTENT(IN)       :: p_pe
    LOGICAL,   INTENT(IN), OPTIONAL :: lnewtl

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracpdef_integrate'
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zxt
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zxtte
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zxtp1
    REAL(DP), DIMENSION(:),     ALLOCATABLE :: ratio
    REAL(DP)                                :: div
    INTEGER                                 :: jt
    INTEGER                                 :: n
    INTEGER                                 :: s1, s2, sx
    INTEGER                                 :: ntrac
    ! um_ak_20110330 LOGICAL                :: l_2tls = .FALSE. !um_ak_20090317
    LOGICAL                                 :: l_newtl = .FALSE.! um_ak_20110330

    INTRINSIC :: MERGE, SUM, MAX, ABS, TINY, TRIM, SIZE

    status = 0

    IF (PRESENT(lnewtl)) l_newtl = lnewtl ! um_ak_20110330

    SELECT CASE(flag)

    CASE(1)

       set_loop1: DO n=1, NSETID

          IF (.NOT. lnow(n)) CYCLE

          ntrac = TRSET(n)%ntrac
          IF (ntrac == 0) CYCLE

          SELECT CASE(TRRANK)
          CASE(3)
             s1 = SIZE(TRSET(n)%xt(:,:,:,:,:), 1)
             s2 = SIZE(TRSET(n)%xt(:,:,:,:,:), 2)
             !s3 = SIZE(TRSET(n)%xt(:,:,:,:,:), 3) ! = ntrac
             sx = SIZE(TRSET(n)%xt(:,:,:,:,:), 4)
             !s5 = SIZE(TRSET(n)%xt(:,:,:,:,:), 5) ! = 1
          CASE(4)
             s1 = SIZE(TRSET(n)%xt(:,:,:,:,:), 1)
             s2 = SIZE(TRSET(n)%xt(:,:,:,:,:), 2)
             sx = SIZE(TRSET(n)%xt(:,:,:,:,:), 3)
             !s4 = SIZE(TRSET(n)%xt(:,:,:,:,:), 4) ! = ntrac
             !s5 = SIZE(TRSET(n)%xt(:,:,:,:,:), 5) ! = 1
          END SELECT

          ! INIT
          ALLOCATE(zxt(s1,s2,sx))
          ALLOCATE(zxtte(s1,s2,sx))
          ALLOCATE(zxtp1(s1,s2,sx))

          diag_only: IF (XWRK(n)%ldiagonly) THEN

             tracer_loop1: DO jt=1, ntrac

                IF (.NOT.XWRK(n)%lok(jt)) CYCLE

                ! INTEGRATE MASS (NEGATIVE)
                SELECT CASE(TRRANK)
                CASE(3)
                   zxt(:,:,:) = &
                        MERGE(TRSET(n)%xt(:,:,jt,:,1), &
                        0.0_DP, TRSET(n)%xt(:,:,jt,:,1)<0.0_DP)
                CASE(4)
                   zxt(:,:,:) = &
                        MERGE(TRSET(n)%xt(:,:,:,jt,1), &
                        0.0_DP, TRSET(n)%xt(:,:,:,jt,1)<0.0_DP)
                END SELECT

                zxt(:,:,:) = XWRK(n)%mask(:,:,:) * &
                     MERGE(0.0_dp, zxt(:,:,:),     &
                     ABS(zxt(:,:,:)-FLAGGED_BAD)<TINY_DP)
                XWRK(n)%mass_pe(jt, 1, p_pe) = &
                     SUM( zxt(:,:,:) *          &
                     XWRK(n)%airmass(:,:,:) *   &
                     XWRK(n)%mask(:,:,:)) *     &
                     XWRK(n)%scalf(jt)

                ! INTEGRATE MASS (POSITIVE)
                SELECT CASE(TRRANK)
                CASE(3)
                   zxt(:,:,:) = &
                        MERGE(TRSET(n)%xt(:,:,jt,:,1), &
                        0.0_DP, TRSET(n)%xt(:,:,jt,:,1)>0.0_DP)
                CASE(4)
                   zxt(:,:,:) = &
                        MERGE(TRSET(n)%xt(:,:,:,jt,1), &
                        0.0_DP, TRSET(n)%xt(:,:,:,jt,1)>0.0_DP)
                END SELECT

                zxt(:,:,:) =  XWRK(n)%mask(:,:,:) * &
                     MERGE(0.0_dp, zxt(:,:,:),      &
                     ABS(zxt(:,:,:)-FLAGGED_BAD)<TINY_DP)
                XWRK(n)%mass_pe(jt, 2, p_pe) =              &
                     SUM( zxt(:,:,:) *          &
                     XWRK(n)%airmass(:,:,:)) *     &
                     XWRK(n)%scalf(jt)

                ! NO CORRECTION POSSIBLE
                XWRK(n)%mass_pe(jt, 3, p_pe) = 0.0_DP

             END DO tracer_loop1

          ELSE

             tracer_loop2: DO jt=1, ntrac

                IF (.NOT.XWRK(n)%lok(jt)) CYCLE

                ! um_ak_20110330 if2tls: IF (.NOT. l_2tls) THEN ! um_ak_20090317
                ifnewtl: IF (.NOT. l_newtl) THEN ! um_ak_20110330
                   ! ACTUAL VALUE AT t+1
                   SELECT CASE(TRRANK)
                   CASE(3)
                      zxtp1(:,:,:) = &
                           TRSET(n)%xtm1(:,:,jt,:,1) + &
                           TRSET(n)%xtte(:,:,jt,:,1)*time_step_len
                   CASE(4)
                      zxtp1(:,:,:) = &
                           TRSET(n)%xtm1(:,:,:,jt,1) + &
                           TRSET(n)%xtte(:,:,:,jt,1)*time_step_len
                   END SELECT

                   ! PSEUDO-TENDENCY TO FORCE RESULTS AT t+1 >= 0.0
                   zxtte(:,:,:) = &
                        MAX(-zxtp1(:,:,:)/time_step_len, 0.0_DP)

                   ! INTEGRATE MASS (NEGATIVE)
                   zxt(:,:,:) = XWRK(n)%mask(:,:,:)  *    &
                        MERGE(zxtp1(:,:,:), 0.0_DP, zxtp1(:,:,:)<0.0_DP)
                   XWRK(n)%mass_pe(jt, 1, p_pe) =              &
                        SUM( zxt(:,:,:) *          &
                        XWRK(n)%airmass(:,:,:)) *  &
                        XWRK(n)%scalf(jt)
                   ! INTEGRATE MASS (POSITIVE)
                   zxt(:,:,:) = XWRK(n)%mask(:,:,:)  *    &
                         MERGE(zxtp1(:,:,:), 0.0_DP, zxtp1(:,:,:)>0.0_DP)
                   XWRK(n)%mass_pe(jt, 2, p_pe) =              &
                        SUM( zxt(:,:,:) *          &
                        XWRK(n)%airmass(:,:,:)) *  &
                        XWRK(n)%scalf(jt)

                   ! RESET NEGATIVES ON REQUEST
                   XWRK(n)%mass_pe(jt, 3, p_pe) = 0.0_DP
                   IF (XWRK(n)%xlswitch(jt,1)) THEN
                      SELECT CASE(TRRANK)
                      CASE(3)
                         TRSET(n)%xtte(:,:,jt,:,1) =      &
                              TRSET(n)%xtte(:,:,jt,:,1) + &
                              zxtte(:,:,:)
                      CASE(4)
                         TRSET(n)%xtte(:,:,:,jt,1) =      &
                              TRSET(n)%xtte(:,:,:,jt,1) + &
                              zxtte(:,:,:)
                      END SELECT
                      XWRK(n)%mass_pe(jt, 3, p_pe) =              &
                           SUM(zxtte(:,:,:) *         &
                           XWRK(n)%airmass(:,:,:)  *  &
                           XWRK(n)%mask(:,:,:) ) *    &
                           XWRK(n)%scalf(jt) * time_step_len
                   END IF
                ELSE ! l_newtl ! um_ak_20090317+

                   ! ACTUAL VALUE AT t+1
                   SELECT CASE(TRRANK)
                   CASE(3)
                      zxtp1(:,:,:) =  TRSET(n)%xt(:,:,jt,:,1)
                   CASE(4)
                      zxtp1(:,:,:) =  TRSET(n)%xt(:,:,:,jt,1)
                   END SELECT

                   ! INTEGRATE MASS (NEGATIVE)
                   zxt(:,:,:) =  XWRK(n)%mask(:,:,:) * &
                       MERGE(zxtp1(:,:,:), 0.0_DP,zxtp1(:,:,:)<0.0_DP)

                   XWRK(n)%mass_pe(jt, 1, p_pe) =  &
                        SUM( zxt(:,:,:) *          &
                        XWRK(n)%airmass(:,:,:) ) * &
                        XWRK(n)%scalf(jt)

                   ! INTEGRATE MASS (POSITIVE)
                   zxt(:,:,:) = XWRK(n)%mask(:,:,:) *&
                        MERGE(zxtp1(:,:,:), 0.0_DP,zxtp1(:,:,:)>0.0_DP)

                   XWRK(n)%mass_pe(jt, 2, p_pe) = &
                        SUM( zxt(:,:,:) *          &
                        XWRK(n)%airmass(:,:,:) ) * &
                        XWRK(n)%scalf(jt)

                   ! RESET NEGATIVES ON REQUEST
                   XWRK(n)%mass_pe(jt, 3, p_pe) = 0.0_DP
                   IF (XWRK(n)%xlswitch(jt,1)) THEN
                      SELECT CASE(TRRANK)
                      CASE(3)
                         TRSET(n)%xt(:,:,jt,:,1) = &
                              MERGE(TRSET(n)%xt(:,:,jt,:,1) &
                              , 0.0_DP,TRSET(n)%xt(:,:,jt,:,1)>0.0_DP)
                      CASE(4)
                         TRSET(n)%xt(:,:,:,jt,1) = &
                              MERGE(TRSET(n)%xt(:,:,:,jt,1) &
                              , 0.0_DP,TRSET(n)%xt(:,:,:,jt,1)>0.0_DP)
                      END SELECT
                      XWRK(n)%mass_pe(jt, 3, p_pe) =              &
                           XWRK(n)%mass_pe(jt, 1, p_pe)
                   END IF

                ENDIF ifnewtl ! um_ak_20090317-
             END DO tracer_loop2

          END IF diag_only

          ! CLEAN UP
          DEALLOCATE(zxt)
          DEALLOCATE(zxtte)
          DEALLOCATE(zxtp1)

       END DO set_loop1

    CASE(2)

       set_loop2: DO n=1, NSETID

          IF (.NOT. lnow(n)) CYCLE

          ntrac = TRSET(n)%ntrac

          IF (ntrac == 0) CYCLE

          ! INIT
          ALLOCATE(ratio(ntrac))
          ratio(:) = 0.0_DP

! op_pj_20120924+
!          Note: This integration is now performed in the BMIL routine
!                main_tracer_pdef_global_end (messy_main_tracer_pdef_bi.f90)
!                with the MPI_ALLREDUCE(...MPI_SUM...) based function p_sum
!                for performance reasons. See note there.
!!$          ! INTEGRATE MASSES
!!$          XWRK(n)%mass_n(:)  = SUM(XWRK(n)%mass_pe(:,1,:),2) ! SUM OVER p_pe
!!$          XWRK(n)%mass_p(:)  = SUM(XWRK(n)%mass_pe(:,2,:),2) ! SUM OVER p_pe
!!$          XWRK(n)%mass_c(:)  = SUM(XWRK(n)%mass_pe(:,3,:),2) ! SUM OVER p_pe
! op_pj_20120924-

          ! CALCULATE RATIO
          DO jt=1, ntrac
             IF (.NOT.XWRK(n)%lok(jt)) THEN
                ratio(jt) = -777.0_DP
                CYCLE
             END IF
             !
             ! CONSERVED MASS
             div = ABS(XWRK(n)%mass_n(jt) + XWRK(n)%mass_p(jt))
             IF (div > TINY(0._DP)) THEN
                ratio(jt) = ABS(XWRK(n)%mass_n(jt))/div
             ELSE
                ratio(jt) = -999.0_DP
             END IF
          END DO

          ! APPLY TOLERANCE CRITERIUM
          tracer_loop3: DO jt=1, ntrac
             IF (.NOT.XWRK(n)%xlswitch(jt,2)) CYCLE
             IF (.NOT.XWRK(n)%lok(jt)) CYCLE

             IF (ratio(jt) > XWRK(n)%xrtol(jt)) THEN

                WRITE(*,*) substr,' :NEGATIVE MASS OF TRACER '''//&
                     &TRIM(TRSET(n)%ti(jt)%tp%ident%fullname)//&
                     &''' IN SET '''//TRIM(TRSET(n)%name)//&
                     &''' EXCEEDS TOLERANCE: '            &
                     , XWRK(n)%mass_n(jt),XWRK(n)%mass_p(jt),ratio(jt)
                status = 1 ! terminate simulation
                RETURN
             END IF
          END DO tracer_loop3

          DEALLOCATE(ratio)

       END DO set_loop2

    END SELECT

  END SUBROUTINE tracpdef_integrate
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_freemem

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: n

    set_loop: DO n=1, NSETID

       DEALLOCATE(XWRK(n)%xlswitch)
       DEALLOCATE(XWRK(n)%xrtol)

       DEALLOCATE(XWRK(n)%mem_mass_n)
       DEALLOCATE(XWRK(n)%mem_mass_p)
       DEALLOCATE(XWRK(n)%mem_mass_c)

       DEALLOCATE(XWRK(n)%airmass)

       DEALLOCATE(XWRK(n)%mass_pe)
       DEALLOCATE(XWRK(n)%lok)
       DEALLOCATE(XWRK(n)%scalf)
       DEALLOCATE(XWRK(n)%unit)

    END DO set_loop

    DEALLOCATE(XWRK)
    DEALLOCATE(lnow)

  END SUBROUTINE tracpdef_freemem
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracpdef_print(ldiagout)

    IMPLICIT NONE
    INTRINSIC :: SUM, TRIM

    ! I/O
    LOGICAL, INTENT(IN) :: ldiagout

    ! LOCAL
    INTEGER :: n, ntrac, jt

    IF (.NOT. L_DIAGOUT) RETURN
    IF (.NOT. ldiagout)  RETURN

    set_loop: DO n=1, NSETID

       ntrac = TRSET(n)%ntrac
       IF (.NOT. lnow(n)) CYCLE
       IF (ntrac == 0) CYCLE

       WRITE(*,*) '========================================================'
       WRITE(*,*) 'TRACER SET: '//TRIM(TRSET(n)%name)
       WRITE(*,*) 'AIR MASS  : ',SUM(XWRK(n)%airmass)
       WRITE(*,*) '--------------------------------------------------------'
       WRITE(*,'(a15,1x,3(a12,1x))') 'TRACER', 'MP', 'MN', 'MC'
       WRITE(*,*) '--------------------------------------------------------'

       tracer_loop1: DO jt=1, ntrac

          IF (.NOT.XWRK(n)%lok(jt)) CYCLE

          WRITE(*,'(a15,1x,3(e12.4,1x),a)') &
               TRSET(n)%ti(jt)%tp%ident%fullname  &
               , XWRK(n)%mass_p(jt), XWRK(n)%mass_n(jt), XWRK(n)%mass_c(jt) &
               , TRIM(XWRK(n)%unit(jt))

       END DO tracer_loop1

       WRITE(*,*) '========================================================'

    END DO set_loop

  END SUBROUTINE tracpdef_print
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracer_pdef_read_nml_ctrl(status, iou)

    ! Author: Patrick Joeckel, MPICH, Aug 2005

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'tracer_pdef_read_nml_ctrl'

    NAMELIST /CTRL_PDEF/ L_DIAGOUT, TPD_DEFAULT, TPD, L_TASK_AGGREGATE

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CTRL_PDEF', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL_PDEF, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL_PDEF', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE tracer_pdef_read_nml_ctrl
! ----------------------------------------------------------------------

! ************************************************************************
END MODULE messy_main_tracer_pdef
! ************************************************************************
