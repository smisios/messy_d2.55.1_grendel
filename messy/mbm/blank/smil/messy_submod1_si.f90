! **********************************************************************
MODULE messy_submod1_si
! **********************************************************************
! This submodel is based on the submodel of the CHANNEL box model
! Authors: Patrick Joeckel, DLR, Oberpfaffenhofen (original code)
!          Astrid Kerkweg, UNI-MZ, Mainz (adapted to example submodel)
! **********************************************************************

  ! SMCL
  USE messy_submod1

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  ! MODULE VARIABLES
  REAL(DP), DIMENSION(:),     POINTER :: state
  REAL(DP), DIMENSION(:,:,:), POINTER :: f01 => NULL()
  REAl(DP),                   POINTER :: s01 => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: fbnd    => NULL()  ! with boundaries
  INTEGER,  DIMENSION(4)              :: nbnd

  PUBLIC :: submod1_init_memory
  PUBLIC :: submod1_global_end
  PUBLIC :: submod1_free_memory

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE submod1_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi, ONLY: SCALAR, GP_3D_MID, ARRAY &
                                   , GP_3D_MID_BND 
    USE messy_main_channel,    ONLY: new_channel, new_channel_object &
                                   , new_attribute, get_channel_object_info

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'submod1_init_memory'
    INTEGER :: status

    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'state', p1=state, reprid=ARRAY &
         , lrestreq = .TRUE.)  ! save for restart
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'state' &
         , 'long_name', c = 'random state vector')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'f01', p3=f01, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'f01' &
         , 'long_name', c = 'random field')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 's01', p0=s01, reprid=SCALAR)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 's01' &
         , 'long_name', c = 'average of random field')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'fbnd', p3=fbnd &
         , reprid=GP_3D_MID_BND)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'fbnd' &
         , 'long_name', c = 'a field with a boundary')
    CALL channel_halt(substr, status)

    CALL get_channel_object_info(status, modstr, 'fbnd', nbounds=nbnd)
    CALL channel_halt(substr, status)

  END SUBROUTINE submod1_init_memory
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod1_global_end

    ! BML/BMIL
    USE messy_main_timer,           ONLY: lresume

    IMPLICIT NONE
    INTRINSIC :: SUM, UBOUND, REAL, SIZE

    ! LOCAL
    INTEGER :: n1, n2, n3
    INTEGER :: i1, i2, i3
    REAL(DP), DIMENSION(:), POINTER :: harvest => NULL()
    LOGICAL, SAVE :: lfirst = .TRUE.
    ! mz_ab_20090921+
    INTEGER, DIMENSION(3) :: s, e
    ! mz_ab_20090921-

    IF (lfirst) THEN
       CALL RANINI(state, lresume)  ! set initial value
       lfirst = .FALSE.             ! only once
    END IF

    n1=SIZE(f01,1)
    n2=SIZE(f01,2)
    n3=SIZE(f01,3)

    ALLOCATE(harvest(n1*n2*n3))
    CALL RANINI(state, .FALSE., harvest)
    DO i1=1, n1
       DO i2=1, n2
          DO i3=1, n3
             f01(i1,i2,i3) = harvest(i1*i2*i3)
          END DO
       END DO
    END DO

    s01 = SUM(f01)/REAL(n1*n2*n3,dp)

    ! mz_ab_20090921+
    s(:) = 1 - nbnd(1:3)
    e(:) = UBOUND(fbnd) - nbnd(1:3)
    CALL populate_bnds(fbnd, s, e)   
    ! mz_ab_20090921-

  END SUBROUTINE submod1_global_end
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod1_free_memory

    IMPLICIT NONE

    CALL RANCLEAN

  END SUBROUTINE submod1_free_memory
  ! --------------------------------------------------------------------

! **********************************************************************
END MODULE messy_submod1_si
! **********************************************************************
