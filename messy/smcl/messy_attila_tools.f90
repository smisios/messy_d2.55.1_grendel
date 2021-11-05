! ***************************************************************************
MODULE messy_attila_tools
! ***************************************************************************

  ! MODULE FOR TRANSFORMATIONS/TRANSPOSITIONS BETWEEN GRIDPOINT SPACE
  ! AND LAGRANGIAN SPACE
  !
  ! ALL TRANSFORMATIONS/TRANSPOSITIONS HERE
  ! ARE CONSISTENT ON EACH CPU SEPARATELY
  !
  ! NOTE:
  ! ECHAM5-DECOMPOSITION CONFORM ROUTINES ARE IN
  ! messy_attila_tools_e5.f90
  !
  ! Authors: Patrick Joeckel and Michael Traub, MPICH, Oct 2003
  !

!#define GENERIC

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: ABS, INDEX, PRESENT, SIZE, SUM, INT &
       , SQRT, REAL, TINY, NINT

  PUBLIC :: DP

  ! ----------- <

  ! order of (x,y,z) in index position vector NPOS
  ! ATTILA:
  !         NPOS(3, :) -> LEV (z)
  !         NPOS(4, :) -> LON (x)
  !         NPOS(5, :) -> LAT (y)
  !
  !         GP(LON, LEV, LAT)
  !
  ! INDEX OF xyz IN NPOS
  INTEGER, DIMENSION(3), PARAMETER :: lgindx = (/ 4, 5, 3 /) ! NPOS
  ! ORDER OF x=1, y=2, z=3 IN ATTILA GP-FIELDS (PMBOX, NGCELL) -> xzy
  INTEGER, DIMENSION(3), PARAMETER :: ovec_a = (/ 1, 3, 2 /) ! PMBOX
  REAL(DP),              PARAMETER, PUBLIC :: EPS = TINY(0.0_DP) ! zero

  PUBLIC :: gp2lg             ! GRIDPOINT -> LAGRANGE
  PUBLIC :: lg2gp             ! LAGRANGE  -> GRIDPOINT
  PUBLIC :: gpsfemis2lgemis   ! GRIDPOINT SURFACE EMISSION -> LAGRANGIAN EMISS.

#ifdef GENERIC

  ! generic routines for arbitrary order; double indexing

  INTERFACE gpsfemis2lgemis
     MODULE PROCEDURE gpsfemis2lgemis_gen ! generic
  END INTERFACE

  INTERFACE gp2lg
     !MODULE PROCEDURE gp2lg_xyzn_gen  ! generic (NOT IMPLEMENTED !!!)
     MODULE PROCEDURE gp2lg_xyz_gen  ! generic
     MODULE PROCEDURE gp2lg_xy_gen   ! generic
  END INTERFACE

  INTERFACE lg2gp
     MODULE PROCEDURE lg2gp_xyz_gen  ! generic
     !MODULE PROCEDURE lg2gp_xyzn_gen ! generic (NOT IMPLEMENTED !!!)
  END INTERFACE

#else

  ! faster, but less flexible; no double indexing

  INTERFACE gpsfemis2lgemis
     MODULE PROCEDURE gpsfemis2lgemis_xy
  END INTERFACE

  INTERFACE gp2lg
     MODULE PROCEDURE gp2lg_xzny
     MODULE PROCEDURE gp2lg_xzy
     MODULE PROCEDURE gp2lg_xy
  END INTERFACE

  INTERFACE lg2gp
     MODULE PROCEDURE lg2gp_xzy
     MODULE PROCEDURE lg2gp_xzny
  END INTERFACE

#endif

  ! FOR lg2gp
  INTEGER, PARAMETER, PUBLIC :: LG2GP_SUM = 1
  INTEGER, PARAMETER, PUBLIC :: LG2GP_AVE = 2
  INTEGER, PARAMETER, PUBLIC :: LG2GP_STD = 3
  INTEGER, PARAMETER, PUBLIC :: LG2GP_AVEGT0 = 4 ! average all elements > 0

CONTAINS

#ifdef GENERIC
  ! =======================================================================
  SUBROUTINE gp2lg_xyz_gen(gp, lg, gpr, lmcons, gorder, status)
    
    ! 1) NON-MASS CONSERVING TRANSFORMATION (DEFAULT):
    !    (lmass_cons = .false.).AND.(.NOT.PRESENT(gpr))
    !    THE LAGRANGIAN CELLS RECEIVE THE VALUE OF THE RESPECTIVE
    !    GRIDPOINT BOX, INDEPENDENT OF THE NUMBER OF LAGRANGIAN CELLS
    !    PER GRIDPOINT BOX.
    !    (temperature, rate-coeff., NON-transported offline
    !    tracers, initialisation of tracers with given volume mixing ratio,
    !    etc.)
    !
    ! 2) LOCAL (IN TIME) MASS CONSERVING TRANSFORMATION:
    !    (lmass_cons = .true.).AND.(.NOT.PRESENT(gpr))
    !    THE VALUE OF THE GRIDPOINT BOX IS MASS-WEIGHTED DISTRIBUTED
    !    AMONG ALL LAGRANGIAN CELLS IN THE RESPECTIVE GRIDPOINT BOX.
    !    INFORMATION IN GRIDPOINT BOXES, WHERE CURRENTLY NO LAGRANGIAN
    !    BOX IS LOCATED, IS LOST.
    !
    ! 3) TOTAL MASS CONSERVING TRANSFORMATION:
    !    (lmass_cons = .true.).AND.(PRESENT(gpr))
    !    THE VALUE OF THE GRIDPOINT BOX IS MASS-WEIGHTED DISTRIBUTED
    !    AMONG ALL LAGRANGIAN CELLS IN THE RESPECTIVE GRIDPOINT BOX.
    !    INFORAMTION IN GRIDPOINT BOXES, WHERE CURRENTLY NO LAGRANGIAN
    !    CELL IS LOCATED, IS SAVED/ACCUMULATED UNTIL THE NEXT LAGRANGIAN
    !    CELL COMES ACROSS (gpr).
    !
    !
    !    NOTES: 
    !    SINCE gpr IS USED TO CONSERVE THE TOTAL MASS, [gpr] SHOULD BE IN KG.
    !    HOWEVER, THE TIMESTEP LENGTH, THE RATIO OF MOLAR MASSES (AIR/TRACER),
    !    ETC. ARE SCALING CONSTANTS (AS LONG AS THE TIME STEP DOES NOT VARY
    !    DURING THE INTEGRATION), AND ARE THEREFORE OMITTED.
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, PMBOX, AMCELL, GNCB

    IMPLICIT NONE

    ! I/O
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:,:), INTENT(IN)              :: gp
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(OUT)             :: lg
    ! grid point field with REST (if no CELL is in GRID-BOX)
    REAL(dp), DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: gpr
    ! for mass conserving transformation
    LOGICAL,                INTENT(IN),    OPTIONAL :: lmcons
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=3),       INTENT(IN),    OPTIONAL :: gorder
    ! error status
    INTEGER,                INTENT(OUT),   OPTIONAL :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gp2lg_xyz_gen'
    INTEGER  :: ovec(3)   ! order of (x,y,z) in global gridpoint field
    INTEGER  :: jn
    INTEGER  :: i1_a, i2_a, i3_a    ! INDEX IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2, i3          ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    LOGICAL  :: zlmrest
    REAL(dp) :: weight
    INTEGER  :: n(3), j(3), jx, jy, jz

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    zlmrest = .FALSE.
    IF (PRESENT(gpr)) THEN
       IF (SIZE(gpr) > 0) zlmrest = .TRUE.
    END IF
    !
    lg(:) = 0.0_dp

    ovec   = (/ 1, 2, 3 /)             ! DEFAULT xyz FOR GP
    IF (PRESENT(gorder)) THEN
       IF ((INDEX(gorder, 'x')) == 0) RETURN
       IF ((INDEX(gorder, 'y')) == 0) RETURN
       IF ((INDEX(gorder, 'z')) == 0) RETURN
       ovec(INDEX(gorder, 'x')) = 1
       ovec(INDEX(gorder, 'y')) = 2
       ovec(INDEX(gorder, 'z')) = 3
       IF (SUM(ovec) /= 6) RETURN
    END IF
    !
    n(1) = SIZE(gp, INDEX(gorder, 'x'))
    n(2) = SIZE(gp, INDEX(gorder, 'y'))
    n(3) = SIZE(gp, INDEX(gorder, 'z'))

    cell_loop: DO jn = 1, NCELL
       i1 = NPOS(lgindx(ovec(1)), jn)
       i2 = NPOS(lgindx(ovec(2)), jn)
       i3 = NPOS(lgindx(ovec(3)), jn)
       !
       i1_a = NPOS(lgindx(ovec_a(1)), jn)
       i2_a = NPOS(lgindx(ovec_a(2)), jn)
       i3_a = NPOS(lgindx(ovec_a(3)), jn)
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (NINT(GNCB(i1_a, i2_a, i3_a)) == 0) THEN
          WRITE(*,*) substr,' WARNING: GNCB=0 for ',&
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = ( (gp(i1, i2, i3)*PMBOX(i1_a, i2_a, i3_a)/AMCELL)   &
                  +   gpr(i1, i2, i3) ) / GNCB(i1_a, i2_a, i3_a)
          ELSE
             lg(jn) = gp(i1, i2, i3)*PMBOX(i1_a, i2_a, i3_a) &
                  /(GNCB(i1_a, i2_a, i3_a)*AMCELL)
          END IF
       ELSE
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = (gp(i1, i2, i3) + gpr(i1, i2, i3))
          ELSE
             lg(jn) = gp(i1, i2, i3)
          END IF
       END IF
    END DO cell_loop
    
    IF (zlmrest) THEN
       ! SAVE GRIDPOINT-REST (WHERE NO CELL IS LOCATED)
       DO jx=1, n(1)           ! x
          DO jy=1, n(2)        ! y
             DO jz=1, n(3)     ! z
                j(:) = (/jx, jy, jz/)
                i1   = j(ovec(1))
                i2   = j(ovec(2))
                i3   = j(ovec(3))
                i1_a = j(ovec_a(1))
                i2_a = j(ovec_a(2))
                i3_a = j(ovec_a(3))
                !
                IF (zlmcons) THEN
                   weight = PMBOX(i1_a, i2_a, i3_a)/AMCELL
                ELSE
                   weight = 1.0_dp
                END IF
                !
                IF (NINT(GNCB(i1_a,i2_a,i3_a)) == 0) THEN
                   gpr(i1,i2,i3) = gpr(i1,i2,i3) + gp(i1,i2,i3) * weight
                ELSE
                   gpr(i1,i2,i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
    END IF

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gp2lg_xyz_gen
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gp2lg_xy_gen(gp, lg, gpr, lmcons, klev, gorder, llev, status)

    ! DESCRIPTION: see gp2lg_xyz
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, PMBOX, AMCELL, GNCB, NLEV

    IMPLICIT NONE

    ! I/O
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:),   INTENT(IN)              :: gp
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(OUT)             :: lg
    ! grid point field with REST MASS (if no CELL is in GRID-BOX)
    REAL(dp), DIMENSION(:,:),   INTENT(INOUT), OPTIONAL :: gpr
    ! for mass conserving transformation
    LOGICAL,                INTENT(IN),    OPTIONAL :: lmcons
    ! gridpoint level
    INTEGER,                INTENT(IN),    OPTIONAL :: klev
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=2),       INTENT(IN),    OPTIONAL :: gorder
    ! TRUE IF LEVEL=ZKLEV
    LOGICAL, DIMENSION(:),   POINTER,      OPTIONAL :: llev
    ! error status
    INTEGER,                INTENT(OUT),   OPTIONAL :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gp2lg_xy_gen'
    INTEGER  :: ovec(2)   ! order of (x,y) in global gridpoint field
    INTEGER  :: jn
    INTEGER  :: i1_a, i2_a, i3_a    ! INDEX IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2              ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    LOGICAL  :: zlmrest
    REAL(dp) :: weight
    INTEGER  :: n(2), j(3), jx, jy
    INTEGER  :: zklev

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    zlmrest = .FALSE.
    IF (PRESENT(gpr)) THEN
       IF (SIZE(gpr) > 0) zlmrest = .TRUE.
    END IF
    !
    IF (PRESENT(klev)) THEN
       zklev = klev
    ELSE
       zklev = NLEV  ! DEAFAULT: surface level
    END IF
    ! 
    lg(:) = 0.0_dp
    !
    IF (PRESENT(llev)) THEN
       ALLOCATE(llev(NCELL))
       llev = .FALSE.
    END IF


    ovec   = (/ 1, 2 /)             ! DEFAULT xy FOR GP
    IF (PRESENT(gorder)) THEN
       IF ((INDEX(gorder, 'x')) == 0) RETURN
       IF ((INDEX(gorder, 'y')) == 0) RETURN
       ovec(INDEX(gorder, 'x')) = 1
       ovec(INDEX(gorder, 'y')) = 2
       IF (SUM(ovec) /= 3) RETURN
    END IF
    !
    n(1) = SIZE(gp, INDEX(gorder, 'x'))
    n(2) = SIZE(gp, INDEX(gorder, 'y'))


    cell_loop: DO jn = 1, NCELL
       IF (INT(NPOS(3, jn)) /= zklev) CYCLE
       !
       i1 = NPOS(lgindx(ovec(1)), jn)
       i2 = NPOS(lgindx(ovec(2)), jn)
       !
       i1_a = NPOS(lgindx(ovec_a(1)), jn)
       i2_a = NPOS(lgindx(ovec_a(2)), jn)   ! = zklev
       i3_a = NPOS(lgindx(ovec_a(3)), jn)
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (NINT(GNCB(i1_a, i2_a, i3_a)) == 0) THEN
          WRITE(*,*) substr,' WARNING: GNCB=0 for ', &
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = ( (gp(i1, i2)*PMBOX(i1_a, i2_a, i3_a)/AMCELL)   &
                  +   gpr(i1, i2) ) / GNCB(i1_a, i2_a, i3_a)
          ELSE
             lg(jn) = gp(i1, i2)*PMBOX(i1_a, i2_a, i3_a) &
                  /(GNCB(i1_a, i2_a, i3_a)*AMCELL)
          END IF
       ELSE
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = (gp(i1, i2) + gpr(i1, i2))
          ELSE
             lg(jn) = gp(i1, i2)
          END IF
       END IF

       ! IF THE LEVEL IS ZKLEV SET THE LOGICAL FIELD TO TRUE
       IF (PRESENT(llev)) THEN
          llev(jn) =.TRUE.
       END IF

    END DO cell_loop
    
    IF (zlmrest) THEN
       ! SAVE GRIDPOINT-REST (WHERE NO CELL IS LOCATED)
       DO jx=1, n(1)           ! x
          DO jy=1, n(2)        ! y
             j(:) = (/ jx, jy, zklev/)
             i1   = j(ovec(1))
             i2   = j(ovec(2))
             i1_a = j(ovec_a(1))
             i2_a = j(ovec_a(2))
             i3_a = j(ovec_a(3))
             !
             IF (zlmcons) THEN
                weight = PMBOX(i1_a, i2_a, i3_a)/AMCELL
             ELSE
                weight = 1.0_dp
             END IF
             !
             IF (NINT(GNCB(i1_a,i2_a,i3_a)) == 0) THEN
                gpr(i1,i2) = gpr(i1,i2) + gp(i1,i2) * weight
             ELSE
                gpr(i1, i2) = 0.0_dp
             ENDIF
          END DO
       END DO
    END IF

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gp2lg_xy_gen
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE lg2gp_xyz_gen(lg, gp, method, lmcons, gorder  &
       , ltm1, status, fngp)
    
    ! 1) NON-MASS CONSERVING TRANSFORMATION (DEFAULT):
    !    (lmass_cons = .false.)
    !    THE GRIDPOINT BOX RECEIVES THE 
    !      - SUM
    !      - UNWEIGHTED AVERAGE (NUMBER OF CELLS PER BOX)
    !      - STANDARD DEVIATION 
    !    OF ALL LAGRANGIAN CELLS IN THE RESPECTIVE GRIDPOINT BOX
    !    (temperature, rate-coeff., NON-transported offline
    !    tracers, initialisation of tracers with given volume mixing ratio,
    !    etc.)
    !
    ! 2) MASS CONSERVING TRANSFORMATION:
    !    (lmass_cons = .true.)
    !    THE MASS OF ALL LAGRANGIAN CELLS IS RELEASED INTO THE
    !    HOSTING GRIDPOINT BOX
    !     - SUM
    !     - MASS WEIGHTED AVERAGE
    !     - STANDARD DEVIATION
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, NPOSM1, PMBOX, AMCELL&
                              , NCB, NCBM1

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! I/O
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(IN)             :: lg
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT)            :: gp
    ! method
    INTEGER,                    INTENT(IN)             :: method
    ! for mass conserving transformation
    LOGICAL,                    INTENT(IN),   OPTIONAL :: lmcons
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=3),           INTENT(IN),   OPTIONAL :: gorder
    ! use positions of CELLS at t-1 ?
    LOGICAL,                    INTENT(IN),   OPTIONAL :: ltm1
    ! error status
    INTEGER,                    INTENT(OUT),  OPTIONAL :: status
    ! number of 'events'
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT),  OPTIONAL :: fngp

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lg2gp_xyz_gen'
    INTEGER  :: ovec(3)   ! order of (x,y,z) in global gridpoint field
    INTEGER  :: jn
    INTEGER  :: i1_a, i2_a, i3_a    ! INDEX IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2, i3          ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    INTEGER  :: s1, s2, s3
    LOGICAL  :: zltm1
    REAL(dp) :: weight
    REAL(dp),  DIMENSION(:,:,:), ALLOCATABLE :: zgph, num, zfngp
    REAL(dp),  DIMENSION(:,:),   POINTER     :: ZNPOS => NULL()
    INTEGER,   DIMENSION(:,:,:), POINTER     :: ZNCB  => NULL()

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    s1 = SIZE(gp,1)
    s2 = SIZE(gp,2)
    s3 = SIZE(gp,3)
    !
    gp(:,:,:) = 0.0_dp
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    IF (PRESENT(ltm1)) THEN
       zltm1 = ltm1
    ELSE
       zltm1 = .FALSE.    ! DEFAULT: POSITIONS at t
    END IF
    !
    IF (zltm1) THEN
       ZNPOS => NPOSM1
       ZNCB  => NCBM1
    ELSE
       ZNPOS => NPOS
       ZNCB  => NCB
    END IF

    ovec   = (/ 1, 2, 3 /)             ! DEFAULT xyz FOR GP
    IF (PRESENT(gorder)) THEN
       IF ((INDEX(gorder, 'x')) == 0) RETURN
       IF ((INDEX(gorder, 'y')) == 0) RETURN
       IF ((INDEX(gorder, 'z')) == 0) RETURN
       ovec(INDEX(gorder, 'x')) = 1
       ovec(INDEX(gorder, 'y')) = 2
       ovec(INDEX(gorder, 'z')) = 3
       IF (SUM(ovec) /= 6) RETURN
    END IF

    IF (method == LG2GP_STD) THEN
       ALLOCATE(zgph(s1,s2,s3))
       zgph(:,:,:) = 0.0_dp
       ALLOCATE(num(s1,s2,s3))
       num(:,:,:) = 0.0_dp
    END IF

    IF (method == LG2GP_AVEGT0) THEN
       ALLOCATE(zfngp(s1,s2,s3))
       zfngp(:,:,:) = 0.0_dp
    END IF

    ! MASS CONSERVING TRANSFORMATION
    cell_loop: DO jn = 1, NCELL
       i1 = ZNPOS(lgindx(ovec(1)), jn)
       i2 = ZNPOS(lgindx(ovec(2)), jn)
       i3 = ZNPOS(lgindx(ovec(3)), jn)
       !
       i1_a = ZNPOS(lgindx(ovec_a(1)), jn)
       i2_a = ZNPOS(lgindx(ovec_a(2)), jn)
       i3_a = ZNPOS(lgindx(ovec_a(3)), jn)
       !
       ! NOTE: NCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (ZNCB(i1_a, i2_a, i3_a) == 0) THEN
          WRITE(*,*) substr,' WARNING: NCB=0 for ', &
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          weight = AMCELL/PMBOX(i1_a, i2_a, i3_a)
       ELSE
          weight = 1.0_dp
       END IF
       !
       SELECT CASE(method)
       CASE(LG2GP_SUM)
          gp(i1, i2, i3) = gp(i1, i2, i3) + lg(jn) * weight
       CASE(LG2GP_AVE)
          gp(i1, i2, i3) = gp(i1, i2, i3) + (lg(jn) * weight) &
               / REAL(ZNCB(i1_a, i2_a, i3_a),dp)
       CASE(LG2GP_STD)
          zgph(i1, i2, i3) = zgph(i1, i2, i3) + (lg(jn) * weight) &
               / REAL(ZNCB(i1_a, i2_a, i3_a),dp)
          gp(i1, i2, i3) = gp(i1, i2, i3) + ((lg(jn) * weight)**2) &
               / REAL(ZNCB(i1_a, i2_a, i3_a),dp)
          num(i1, i2, i3) = REAL(ZNCB(i1_a, i2_a, i3_a),dp)
       CASE(LG2GP_AVEGT0)
          IF (lg(jn) > 0.0_dp) THEN
             gp(i1, i2, i3) = gp(i1, i2, i3) + (lg(jn) * weight)
             zfngp(i1, i2, i3) = zfngp(i1, i2, i3) + 1.0_dp
          END IF
       CASE DEFAULT
          STATUS = 2
          RETURN     ! ERROR: UNKNOWN METHOD
       END SELECT
    END DO cell_loop

    IF (method == LG2GP_STD) THEN
       ! STANDARD DEVIATION (<.> = average = SUM(.)/N)
       ! SIGMA = SQRT( SUM((xi-<x>)**2)/(N-1) )
       !       = SQRT(| N<x2> - N<x>2 | / (N-1))
       !       = SQRT(| <x2> - <x>2 |) * SQRT(N/(N-1))
       DO i1=1, s1
          DO i2=1, s2
             DO i3=1, s3
                IF (ABS(num(i1,i2,i3)-1.0_dp) > EPS) THEN
                   num(i1,i2,i3) = num(i1,i2,i3) / (num(i1,i2,i3) - 1.0_dp) 
                ELSE
                   num(i1,i2,i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
       !
       gp(:,:,:) = SQRT( ABS(gp(:,:,:)-zgph(:,:,:)**2) * num(:,:,:) )
    END IF

    IF (method == LG2GP_AVEGT0) THEN
       DO i1=1, s1
          DO i2=1, s2
             DO i3=1, s3
                IF (zfngp(i1, i2, i3) > 0.0_dp) THEN
                   gp(i1, i2, i3) = gp(i1, i2, i3) / zfngp(i1, i2, i3)
                ELSE
                   gp(i1, i2, i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
       IF (PRESENT(fngp)) fngp(:,:,:) = zfngp(:,:,:)
    END IF

    ! CLEAN MEMORY
    IF (ALLOCATED(zgph))  DEALLOCATE(zgph)
    IF (ALLOCATED(num))   DEALLOCATE(num)
    IF (ALLOCATED(zfngp)) DEALLOCATE(zfngp)

    IF (PRESENT(status)) status = 0

  END SUBROUTINE lg2gp_xyz_gen
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gpsfemis2lgemis_gen(gp, lg, method, gpr, lmcons, gorder, status)

    ! DESCRIPTION: convert gridpoint emission flux to lagrangian emission flux
    !
    ! METHODS:
    !      1: SURFACE LAYER + REST
    !      2: LOWEST CELL IN BOUNDARY LAYER + REST
    !      3: ALL CELLS in BOUNDARY LAYER (weighted) + REST
    !      4: ALL CELLS in BOUNDARY LAYER (unweighted) + REST
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, PMBOX, AMCELL, GNCB, NLEV, KHPBL

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! I/O
    ! global grid point emission flux 
    REAL(dp), DIMENSION(:,:),   INTENT(IN)              :: gp
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(OUT)             :: lg
    ! METHOD
    INTEGER,                    INTENT(IN)              :: method
    ! grid point field with REST MASS (if no CELL is in GRID-BOX)
    REAL(dp), DIMENSION(:,:),   INTENT(INOUT), OPTIONAL :: gpr
    ! for mass conserving transformation
    LOGICAL,                    INTENT(IN),    OPTIONAL :: lmcons
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=2),           INTENT(IN),    OPTIONAL :: gorder
    ! error status
    INTEGER,                    INTENT(OUT),   OPTIONAL :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gpsfemis2lgemis_gen'
    INTEGER  :: ovec(2)   ! order of (x,y) in global gridpoint field
    INTEGER  :: jn, jk
    INTEGER  :: i1_a, i2_a, i3_a    ! INDICES IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2              ! INDICES IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    LOGICAL  :: zlmrest
    INTEGER  :: n(2), j(3), jx, jy
    INTEGER  :: ilev, flev
    REAL(dp) :: weight
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: gpn

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    lg(:) = 0.0_dp
    ! 
    zlmrest = .FALSE.
    IF (PRESENT(gpr)) THEN
       IF (SIZE(gpr) > 0) zlmrest = .TRUE.
    END IF

    ! NUMBER OF CELLS FOR EMISSION
    IF (zlmrest) THEN
       ALLOCATE(gpn(SIZE(gpr,1),SIZE(gpr,2)))
       gpn(:,:) = 0
    END IF

    ovec   = (/ 1, 2 /)             ! DEFAULT xy FOR GP
    IF (PRESENT(gorder)) THEN
       IF ((INDEX(gorder, 'x')) == 0) RETURN
       IF ((INDEX(gorder, 'y')) == 0) RETURN
       ovec(INDEX(gorder, 'x')) = 1
       ovec(INDEX(gorder, 'y')) = 2
       IF (SUM(ovec) /= 3) RETURN
    END IF
    !
    n(1) = SIZE(gp, INDEX(gorder, 'x'))
    n(2) = SIZE(gp, INDEX(gorder, 'y'))

    cell_loop: DO jn = 1, NCELL
       !
       i1 = NPOS(lgindx(ovec(1)), jn)
       i2 = NPOS(lgindx(ovec(2)), jn)
       !
       i1_a = NPOS(lgindx(ovec_a(1)), jn)  ! longitude index
       i2_a = NPOS(lgindx(ovec_a(2)), jn)  ! level index
       i3_a = NPOS(lgindx(ovec_a(3)), jn)  ! latitude index
       !
       ilev = INT(NPOS(3, jn))  ! SHOULD BE = i2_a ???
       !
       ! KHPBL is first layer in free troposphere
       flev = KHPBL(i1_a, i3_a) + 1
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (NINT(GNCB(i1_a, i2_a, i3_a)) == 0) THEN
          WRITE(*,*) substr,' WARNING: GNCB=0 for ', &
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       SELECT CASE(METHOD)
       CASE(1)
          ! PUT EMISSION EXCLUSIVELY IN SURFACE LAYER CELL(S)
          IF (ilev /= nlev) CYCLE
          ! NOTE: ALL CELLS IN ONE GRIDPOINT BOX
          weight = 1.0_dp / GNCB(i1_a, i2_a, i3_a)
          !
       CASE(2)
          ! DISTRIBUTE EMISSION AMONG LOWEST BOUNDARY LAYER CELL(S)
          ! NOTES:
          IF (ilev < flev) CYCLE
          !
          !  - EMISSION ONLY, IF NO CELL BELOW THIS LAYER
          IF (ilev < nlev) THEN
             IF (ABS(SUM(GNCB(i1_a,ilev+1:nlev,i3_a))) > EPS) CYCLE
          END IF
          ! NOTE: ALL CELLS ARE IN ONE GRIDPOINT BOX
          weight = 1.0_dp / GNCB(i1_a, i2_a, i3_a)
          !
       CASE(3)
          ! DISTRIBUTE EMISSION AMONG ALL CELLS IN BOUNDARY LAYER
          ! (WEIGHT WITH DISTANCE FROM SURFACE => SOME VERT. DIFF.)
          IF (ilev < flev) CYCLE
          ! NOTE: CELLS ARE LOCATED IN COLUMN BETWEEN flev AND nlev,
          !       HOWEVER, WEIGHTED MASS PER LAYER IS DISTRIBUTED 
          !       EQUALLY AMONG CELLS OF THIS LAYER
          weight = 0.0_dp
          ! SUM ALL CELLS IN BL AND WEIGHT WITH LEVEL
          DO jk=flev,nlev
             weight = weight + REAL(jk,dp)
          END DO
          weight = REAL(ilev,dp) / weight  ! normed weight for layer ilev
          !
! op_sb_20140312+
!!$          weight = weight / GNCB(i1_a, jk, i3_a)
          weight = weight / GNCB(i1_a, ilev, i3_a)
! op_sb_20140312-
          !
       CASE(4)
          ! DISTRIBUTE EMISSION AMONG ALL CELLS IN BOUNDARY LAYER
          ! (UNWEIGHTED => MAX. VERTICAL DIFFUSION)
          IF (ilev < flev) CYCLE
          ! NOTE: CELLS ARE LOCATED IN COLUMN BETWEEN flev AND nlev
          weight = 1.0_dp / SUM(GNCB(i1_a,flev:nlev,i3_a))
          !
       CASE DEFAULT
          ! UNKNOWN METHOD
          IF (PRESENT(status)) status = 2
          RETURN
          !
       END SELECT
       !
       IF (zlmcons) THEN   ! MASS CONSERVING TRANSFORMATION
          ! NOTE: GRIDPOINT EMISSION IS AT LEVEL nlev
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = ( (gp(i1, i2)*PMBOX(i1_a, nlev, i3_a)/AMCELL) &
                    +   gpr(i1, i2) ) * weight
          ELSE
             lg(jn) =   (gp(i1, i2)*PMBOX(i1_a, nlev, i3_a)/AMCELL) &
                    * weight
          END IF
       ELSE               ! NON-MASS CONSERVING TRANSFORMATION
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = (gp(i1, i2) + gpr(i1, i2))
          ELSE
             lg(jn) = gp(i1, i2)
          END IF
       END IF             ! MASS CONSERVING ?
       
    END DO cell_loop

    ! RE-SET REST MASS FLUX
    IF (zlmrest) THEN
       gridbox_loop: DO jx=1, n(1)  ! x
          DO jy=1, n(2)             ! y
             j(:) = (/ jx, jy, nlev/)
             i1   = j(ovec(1))
             i2   = j(ovec(2))
             i1_a = j(ovec_a(1))
             i2_a = j(ovec_a(2))    ! = nlev
             i3_a = j(ovec_a(3))
             !
             ! UPDATE NUMBER OF CELLS FOR EMISSION
             ! KHPBL is first layer in free troposphere
             flev = KHPBL(i1_a, i3_a) + 1             
             !
             SELECT CASE(METHOD)
             CASE(1)
                !
                ! PUT EMISSION EXCLUSIVELY IN SURFACE LAYER CELL(S)
                gpn(i1, i2) = NINT(GNCB(i1_a, i2_a, i3_a))
                !
             CASE(2)
                !
                ! DISTRIBUTE EMISIION AMONG LOWEST BOUNDARY LAYER CELL(S)
                DO jk=nlev, flev, -1
                   gpn(i1, i2) = NINT(GNCB(i1_a, jk, i3_a))
                   IF (gpn(i1, i2) > 0 ) EXIT
                END DO
                !
             CASE(3, 4)
                !
                ! DISTRIBUTE EMISSION AMONG ALL CELLS IN BOUNDARY LAYER
                DO jk=nlev, flev, -1
                   gpn(i1, i2) = gpn(i1, i2) + NINT(GNCB(i1_a, jk, i3_a))
                END DO
                !
             CASE DEFAULT
                ! UNKNOWN METHOD
                IF (PRESENT(status)) status = 2
                RETURN
                !
             END SELECT
             !
             IF (gpn(i1, i2) == 0) THEN
                ! NOTE: GRIDPOINT EMISSION IS AT LEVEL nlev
                IF (zlmcons) THEN
                   gpr(i1,i2) = gpr(i1,i2) &
                        + gp(i1,i2)*PMBOX(i1_a, nlev, i3_a)/AMCELL
                ELSE
                   gpr(i1,i2) = gpr(i1,i2) + gp(i1,i2)
                END IF
             ELSE
                gpr(i1, i2) = 0.0_dp
             END IF
             !
          END DO
       END DO gridbox_loop
    END IF

    ! CLEAN MEMORY
    IF (ALLOCATED(gpn)) DEALLOCATE(gpn)

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gpsfemis2lgemis_gen
  ! =======================================================================

  ! =======================================================================
  ! #######################################################################
#else
  ! #######################################################################
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gp2lg_xzy(gp, lg, gpr, lmcons, gorder, status)
    
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, PMBOX, AMCELL, GNCB

    IMPLICIT NONE

    ! I/O
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:,:), INTENT(IN)              :: gp
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(OUT)             :: lg
    ! grid point field with REST (if no CELL is in GRID-BOX)
    REAL(dp), DIMENSION(:,:,:), INTENT(INOUT), OPTIONAL :: gpr
    ! for mass conserving transformation
    LOGICAL,                    INTENT(IN),    OPTIONAL :: lmcons
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=3),           INTENT(IN),    OPTIONAL :: gorder
    ! error status
    INTEGER,                    INTENT(OUT),   OPTIONAL :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gp2lg_xzy'
    INTEGER  :: jn
    INTEGER  :: i1_a, i2_a, i3_a    ! INDEX IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2, i3          ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    LOGICAL  :: zlmrest
    REAL(dp) :: weight
    INTEGER  :: n1, n2, n3

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(gorder)) THEN
       IF (gorder /= 'xzy') THEN
          WRITE(*,*) substr, &
               ': WRONG ORDER! PLEASE RECOMPILE WITH ''#define GENERIC'''
          RETURN
       END IF
    END IF
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    zlmrest = .FALSE.
    IF (PRESENT(gpr)) THEN
       IF (SIZE(gpr) > 0) zlmrest = .TRUE.
    END IF
    !
    lg(:) = 0.0_dp

    cell_loop: DO jn = 1, NCELL
       i1 = NPOS(4, jn) ! x
       i2 = NPOS(3, jn) ! z
       i3 = NPOS(5, jn) ! y
       !
       i1_a = NPOS(4, jn) ! x
       i2_a = NPOS(3, jn) ! z
       i3_a = NPOS(5, jn) ! y
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (NINT(GNCB(i1_a, i2_a, i3_a)) == 0) THEN
          WRITE(*,*) substr,' WARNING: GNCB=0 for ',&
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = ( (gp(i1, i2, i3)*PMBOX(i1_a, i2_a, i3_a)/AMCELL)   &
                  +   gpr(i1, i2, i3) ) / GNCB(i1_a, i2_a, i3_a)
          ELSE
             lg(jn) = gp(i1, i2, i3)*PMBOX(i1_a, i2_a, i3_a) &
                  /(GNCB(i1_a, i2_a, i3_a)*AMCELL)
          END IF
       ELSE
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = (gp(i1, i2, i3) + gpr(i1, i2, i3))
          ELSE
             lg(jn) = gp(i1, i2, i3)
          END IF
       END IF
    END DO cell_loop
    
    IF (zlmrest) THEN
       n1 = SIZE(gp, 1) ! x
       n2 = SIZE(gp, 2) ! z
       n3 = SIZE(gp, 3) ! y
       ! SAVE GRIDPOINT-REST (WHERE NO CELL IS LOCATED)
       DO i1=1, n1           ! x
          DO i2=1, n2        ! z
             DO i3=1, n3     ! y
                !
                IF (zlmcons) THEN
                   weight = PMBOX(i1, i2, i3)/AMCELL  ! xzy
                ELSE
                   weight = 1.0_dp
                END IF
                !
                IF (NINT(GNCB(i1,i2,i3)) == 0) THEN
                   gpr(i1,i2,i3) = gpr(i1,i2,i3) + gp(i1,i2,i3) * weight
                ELSE
                   gpr(i1,i2,i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
    END IF

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gp2lg_xzy
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gp2lg_xzny(gp, lg, gpr, lmcons, gorder, status)
    
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, PMBOX, AMCELL, GNCB

    IMPLICIT NONE

    ! I/O
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:,:,:), INTENT(IN)              :: gp
    ! lagrange vetor
    REAL(dp), DIMENSION(:,:),     INTENT(OUT)             :: lg
    ! grid point field with REST (if no CELL is in GRID-BOX)
    REAL(dp), DIMENSION(:,:,:,:), INTENT(INOUT), OPTIONAL :: gpr
    ! for mass conserving transformation
    LOGICAL,                      INTENT(IN),    OPTIONAL :: lmcons
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=4),             INTENT(IN),    OPTIONAL :: gorder
    ! error status
    INTEGER,                      INTENT(OUT),   OPTIONAL :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gp2lg_xzny'
    INTEGER  :: jn
    INTEGER  :: i1_a, i2_a, i3_a    ! INDEX IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2, i3          ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    LOGICAL  :: zlmrest
    REAL(dp) :: weight
    INTEGER  :: n1, n2, n3

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(gorder)) THEN
       IF (gorder /= 'xzny') THEN
          WRITE(*,*) substr, &
               ': WRONG ORDER! NOT IMPLEMENTED AS GENERIC!'
          RETURN
       END IF
    END IF
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    zlmrest = .FALSE.
    IF (PRESENT(gpr)) THEN
       IF (SIZE(gpr) > 0) zlmrest = .TRUE.
    END IF
    !
    lg(:,:) = 0.0_dp

    cell_loop: DO jn = 1, NCELL
       i1 = NPOS(4, jn) ! x
       i2 = NPOS(3, jn) ! z
       i3 = NPOS(5, jn) ! y
       !
       i1_a = NPOS(4, jn) ! x
       i2_a = NPOS(3, jn) ! z
       i3_a = NPOS(5, jn) ! y
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (NINT(GNCB(i1_a, i2_a, i3_a)) == 0) THEN
          WRITE(*,*) substr,' WARNING: GNCB=0 for ',&
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn,:) = ( (gp(i1, i2, :, i3)*PMBOX(i1_a, i2_a, i3_a)/AMCELL) &
                  +   gpr(i1, i2, :, i3) ) / GNCB(i1_a, i2_a, i3_a)
          ELSE
             lg(jn,:) = gp(i1, i2, :, i3)*PMBOX(i1_a, i2_a, i3_a) &
                  /(GNCB(i1_a, i2_a, i3_a)*AMCELL)
          END IF
       ELSE
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn,:) = (gp(i1, i2, :, i3) + gpr(i1, i2, :, i3))
          ELSE
             lg(jn,:) = gp(i1, i2, :, i3)
          END IF
       END IF
    END DO cell_loop
    
    IF (zlmrest) THEN
       n1 = SIZE(gp, 1) ! x
       n2 = SIZE(gp, 2) ! z
       n3 = SIZE(gp, 3) ! y
       ! SAVE GRIDPOINT-REST (WHERE NO CELL IS LOCATED)
       DO i1=1, n1           ! x
          DO i2=1, n2        ! z
             DO i3=1, n3     ! y
                !
                IF (zlmcons) THEN
                   weight = PMBOX(i1, i2, i3)/AMCELL  ! xzy
                ELSE
                   weight = 1.0_dp
                END IF
                !
                IF (NINT(GNCB(i1,i2,i3)) == 0) THEN
                   gpr(i1,i2,:,i3) = gpr(i1,i2,:,i3) + gp(i1,i2,:,i3) * weight
                ELSE
                   gpr(i1,i2,:,i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
    END IF

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gp2lg_xzny
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gp2lg_xy(gp, lg, gpr, lmcons, klev, gorder, llev, status)

    ! DESCRIPTION: see gp2lg_xyz
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, PMBOX, AMCELL, GNCB, NLEV

    IMPLICIT NONE

    ! I/O
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:),   INTENT(IN)              :: gp
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(OUT)             :: lg
    ! grid point field with REST MASS (if no CELL is in GRID-BOX)
    REAL(dp), DIMENSION(:,:),   INTENT(INOUT), OPTIONAL :: gpr
    ! for mass conserving transformation
    LOGICAL,                    INTENT(IN),    OPTIONAL :: lmcons
    ! gridpoint level
    INTEGER,                    INTENT(IN),    OPTIONAL :: klev
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=2),           INTENT(IN),    OPTIONAL :: gorder
    ! TRUE IF LEVEL=ZKLEV
    LOGICAL, DIMENSION(:),      POINTER,      OPTIONAL :: llev
    ! error status
    INTEGER,                    INTENT(OUT),   OPTIONAL :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gp2lg_xy'
    INTEGER  :: jn
    INTEGER  :: i1_a, i2_a     ! INDEX IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2         ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    LOGICAL  :: zlmrest
    REAL(dp) :: weight
    INTEGER  :: n1, n2
    INTEGER  :: zklev

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(gorder)) THEN
       IF (gorder /= 'xy') THEN
          WRITE(*,*) substr, &
               ': WRONG ORDER! PLEASE RECOMPILE WITH ''#define GENERIC'''
          RETURN
       END IF
    END IF
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    zlmrest = .FALSE.
    IF (PRESENT(gpr)) THEN
       IF (SIZE(gpr) > 0) zlmrest = .TRUE.
    END IF
    !
    IF (PRESENT(klev)) THEN
       zklev = klev
    ELSE
       zklev = NLEV  ! DEAFAULT: surface level
    END IF
    ! 
    lg(:) = 0.0_dp
    !
    IF (PRESENT(llev)) THEN
       ALLOCATE(llev(NCELL))
       llev = .FALSE.
    END IF

    cell_loop: DO jn = 1, NCELL
       IF (INT(NPOS(3, jn)) /= zklev) CYCLE  ! z
       !
       i1 = NPOS(4, jn) ! x
       i2 = NPOS(5, jn) ! y
       !
       i1_a = NPOS(4, jn) ! x
       i2_a = NPOS(5, jn) ! y
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, zklev, i2_a) !!!
       IF (NINT(GNCB(i1_a, zklev, i2_a)) == 0) THEN
          WRITE(*,*) substr,' WARNING: GNCB=0 for ', &
               jn, i1_a, zklev, i2_a
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = ( (gp(i1, i2)*PMBOX(i1_a, zklev, i2_a)/AMCELL)   &
                  +   gpr(i1, i2) ) / GNCB(i1_a, zklev, i2_a)
          ELSE
             lg(jn) = gp(i1, i2)*PMBOX(i1_a, zklev, i2_a) &
                  /(GNCB(i1_a, zklev, i2_a)*AMCELL)
          END IF
       ELSE
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = (gp(i1, i2) + gpr(i1, i2))
          ELSE
             lg(jn) = gp(i1, i2)
          END IF
       END IF

       ! IF THE LEVEL IS ZKLEV SET THE LOGICAL FIELD TO TRUE
       IF (PRESENT(llev)) THEN
          llev(jn) =.TRUE.
       END IF

    END DO cell_loop
    
    IF (zlmrest) THEN
       n1 = SIZE(gp, 1) ! x
       n2 = SIZE(gp, 2) ! y
       ! SAVE GRIDPOINT-REST (WHERE NO CELL IS LOCATED)
       DO i1=1, n1           ! x
          DO i2=1, n2        ! y
             !
             IF (zlmcons) THEN
                weight = PMBOX(i1, zklev, i2)/AMCELL ! xzy
             ELSE
                weight = 1.0_dp
             END IF
             !
             IF (NINT(GNCB(i1,zklev,i2)) == 0) THEN
                gpr(i1,i2) = gpr(i1,i2) + gp(i1,i2) * weight
             ELSE
                gpr(i1, i2) = 0.0_dp
             ENDIF
          END DO
       END DO
    END IF

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gp2lg_xy
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE lg2gp_xzy(lg, gp, method, lmcons, gorder  &
       , ltm1, status, fngp)
    
    ! 1) NON-MASS CONSERVING TRANSFORMATION (DEFAULT):
    !    (lmass_cons = .false.)
    !    THE GRIDPOINT BOX RECEIVES THE 
    !      - SUM
    !      - UNWEIGHTED AVERAGE (NUMBER OF CELLS PER BOX)
    !      - STANDARD DEVIATION 
    !    OF ALL LAGRANGIAN CELLS IN THE RESPECTIVE GRIDPOINT BOX
    !    (temperature, rate-coeff., NON-transported offline
    !    tracers, initialisation of tracers with given volume mixing ratio,
    !    etc.)
    !
    ! 2) MASS CONSERVING TRANSFORMATION:
    !    (lmass_cons = .true.)
    !    THE MASS OF ALL LAGRANGIAN CELLS IS RELEASED INTO THE
    !    HOSTING GRIDPOINT BOX
    !     - SUM
    !     - MASS WEIGHTED AVERAGE
    !     - STANDARD DEVIATION
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, NPOSM1, PMBOX, AMCELL&
                              , NCB, NCBM1

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! I/O
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(IN)             :: lg
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT)            :: gp
    ! method
    INTEGER,                    INTENT(IN)             :: method
    ! for mass conserving transformation
    LOGICAL,                    INTENT(IN),   OPTIONAL :: lmcons
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=3),           INTENT(IN),   OPTIONAL :: gorder
    ! use positions of CELLS at t-1 ?
    LOGICAL,                    INTENT(IN),   OPTIONAL :: ltm1
    ! error status
    INTEGER,                    INTENT(OUT),  OPTIONAL :: status
    ! number of 'events'
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT),  OPTIONAL :: fngp

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lg2gp_xzy'
    INTEGER  :: jn
    INTEGER  :: i1_a, i2_a, i3_a    ! INDEX IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2, i3          ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    INTEGER  :: s1, s2, s3
    LOGICAL  :: zltm1
    REAL(dp) :: weight
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zgph, num, zfngp
    REAL(dp), DIMENSION(:,:),   POINTER     :: ZNPOS => NULL()
    INTEGER,  DIMENSION(:,:,:), POINTER     :: ZNCB  => NULL()

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(gorder)) THEN
       IF (gorder /= 'xzy') THEN
          WRITE(*,*) substr, &
               ': WRONG ORDER! PLEASE RECOMPILE WITH ''#define GENERIC'''
          RETURN
       END IF
    END IF
    !
    s1 = SIZE(gp,1)
    s2 = SIZE(gp,2)
    s3 = SIZE(gp,3)
    !
    gp(:,:,:) = 0.0_dp
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    IF (PRESENT(ltm1)) THEN
       zltm1 = ltm1
    ELSE
       zltm1 = .FALSE.    ! DEFAULT: POSITIONS at t
    END IF
    !
    IF (zltm1) THEN
       ZNPOS => NPOSM1
       ZNCB  => NCBM1
    ELSE
       ZNPOS => NPOS
       ZNCB  => NCB
    END IF

    IF (method == LG2GP_STD) THEN
       ALLOCATE(zgph(s1,s2,s3))
       zgph(:,:,:) = 0.0_dp
       ALLOCATE(num(s1,s2,s3))
       num(:,:,:) = 0.0_dp
    END IF

    IF (method == LG2GP_AVEGT0) THEN
       ALLOCATE(zfngp(s1,s2,s3))
       zfngp(:,:,:) = 0.0_dp
    END IF

    ! MASS CONSERVING TRANSFORMATION
    cell_loop: DO jn = 1, NCELL
       i1 = ZNPOS(4, jn) ! x
       i2 = ZNPOS(3, jn) ! z
       i3 = ZNPOS(5, jn) ! y
       !
       i1_a = ZNPOS(4, jn) ! x
       i2_a = ZNPOS(3, jn) ! z
       i3_a = ZNPOS(5, jn) ! y
       !
       ! NOTE: NCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (ZNCB(i1_a, i2_a, i3_a) == 0) THEN
          WRITE(*,*) substr,' WARNING: NCB=0 for ', &
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          weight = AMCELL/PMBOX(i1_a, i2_a, i3_a)
       ELSE
          weight = 1.0_dp
       END IF
       !
       SELECT CASE(method)
       CASE(LG2GP_SUM)
          gp(i1, i2, i3) = gp(i1, i2, i3) + lg(jn) * weight
       CASE(LG2GP_AVE)
          gp(i1, i2, i3) = gp(i1, i2, i3) + (lg(jn) * weight) &
               / REAL(ZNCB(i1_a, i2_a, i3_a), dp)
       CASE(LG2GP_STD)
          zgph(i1, i2, i3) = zgph(i1, i2, i3) + (lg(jn) * weight) &
               / REAL(ZNCB(i1_a, i2_a, i3_a), dp)
          gp(i1, i2, i3) = gp(i1, i2, i3) + ((lg(jn) * weight)**2) &
               / REAL(ZNCB(i1_a, i2_a, i3_a), dp)
          num(i1, i2, i3) = REAL(ZNCB(i1_a, i2_a, i3_a), dp)
       CASE(LG2GP_AVEGT0)
          IF (lg(jn) > 0.0_dp) THEN
             gp(i1, i2, i3) = gp(i1, i2, i3) + (lg(jn) * weight)
             zfngp(i1, i2, i3) = zfngp(i1, i2, i3) + 1.0_dp
          END IF
       CASE DEFAULT
          STATUS = 2
          RETURN     ! ERROR: UNKNOWN METHOD
       END SELECT
    END DO cell_loop

    IF (method == LG2GP_STD) THEN
       ! STANDARD DEVIATION (<.> = average = SUM(.)/N)
       ! SIGMA = SQRT( SUM((xi-<x>)**2)/(N-1) )
       !       = SQRT(| N<x2> - N<x>2 | / (N-1))
       !       = SQRT(| <x2> - <x>2 |) * SQRT(N/(N-1))
       DO i1=1, s1
          DO i2=1, s2
             DO i3=1, s3
                IF (ABS(num(i1,i2,i3)-1.0_dp) > EPS) THEN
                   num(i1,i2,i3) = num(i1,i2,i3) / (num(i1,i2,i3) - 1.0_dp) 
                ELSE
                   num(i1,i2,i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
       !
       gp(:,:,:) = SQRT( ABS(gp(:,:,:)-zgph(:,:,:)**2) * num(:,:,:) )
    END IF

    IF (method == LG2GP_AVEGT0) THEN
       DO i1=1, s1
          DO i2=1, s2
             DO i3=1, s3
                IF (zfngp(i1, i2, i3) > 0.0_dp) THEN
                   gp(i1, i2, i3) = gp(i1, i2, i3) / zfngp(i1, i2, i3)
                ELSE
                   gp(i1, i2, i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
       IF (PRESENT(fngp)) fngp(:,:,:) = zfngp(:,:,:)
    END IF

    ! CLEAN MEMORY
    IF (ALLOCATED(zgph))  DEALLOCATE(zgph)
    IF (ALLOCATED(num))   DEALLOCATE(num)
    IF (ALLOCATED(zfngp)) DEALLOCATE(zfngp)

    IF (PRESENT(status)) status = 0

  END SUBROUTINE lg2gp_xzy
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE lg2gp_xzny(lg, gp, method, lmcons, gorder  &
       , ltm1, status, fngp)
    
    ! 1) NON-MASS CONSERVING TRANSFORMATION (DEFAULT):
    !    (lmass_cons = .false.)
    !    THE GRIDPOINT BOX RECEIVES THE 
    !      - SUM
    !      - UNWEIGHTED AVERAGE (NUMBER OF CELLS PER BOX)
    !      - STANDARD DEVIATION 
    !    OF ALL LAGRANGIAN CELLS IN THE RESPECTIVE GRIDPOINT BOX
    !    (temperature, rate-coeff., NON-transported offline
    !    tracers, initialisation of tracers with given volume mixing ratio,
    !    etc.)
    !
    ! 2) MASS CONSERVING TRANSFORMATION:
    !    (lmass_cons = .true.)
    !    THE MASS OF ALL LAGRANGIAN CELLS IS RELEASED INTO THE
    !    HOSTING GRIDPOINT BOX
    !     - SUM
    !     - MASS WEIGHTED AVERAGE
    !     - STANDARD DEVIATION
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, NPOSM1, PMBOX, AMCELL&
                              , NCB, NCBM1

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! I/O
    ! lagrange vetor
    REAL(dp), DIMENSION(:,:),     INTENT(IN)             :: lg
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:,:,:), INTENT(OUT)            :: gp
    ! method
    INTEGER,                      INTENT(IN)             :: method
    ! for mass conserving transformation
    LOGICAL,                      INTENT(IN),   OPTIONAL :: lmcons
    ! order of x,y,z,n - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=4),             INTENT(IN),   OPTIONAL :: gorder
    ! use positions of CELLS at t-1 ?
    LOGICAL,                      INTENT(IN),   OPTIONAL :: ltm1
    ! error status
    INTEGER,                      INTENT(OUT),  OPTIONAL :: status
    ! number of 'events'
    REAL(dp), DIMENSION(:,:,:,:), INTENT(OUT),  OPTIONAL :: fngp

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lg2gp_xzny'
    INTEGER  :: jn
    INTEGER  :: i1_a, i2_a, i3_a    ! INDEX IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2, i3          ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    INTEGER  :: s1, s2, s3, s4, i
    LOGICAL  :: zltm1
    REAL(dp) :: weight
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: zgph
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE :: num
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: zfngp
    REAL(dp), DIMENSION(:,:),     POINTER     :: ZNPOS => NULL()
    INTEGER,  DIMENSION(:,:,:),   POINTER     :: ZNCB  => NULL()

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(gorder)) THEN
       IF (gorder /= 'xzny') THEN
          WRITE(*,*) substr, &
               ': WRONG ORDER! NOT IMPLEMENTED AS GENERIC!'
          RETURN
       END IF
    END IF
    !
    s1 = SIZE(gp,1)
    s2 = SIZE(gp,2)
    s3 = SIZE(gp,3)
    s4 = SIZE(gp,4)
    !
    gp(:,:,:,:) = 0.0_dp
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    IF (PRESENT(ltm1)) THEN
       zltm1 = ltm1
    ELSE
       zltm1 = .FALSE.    ! DEFAULT: POSITIONS at t
    END IF
    !
    IF (zltm1) THEN
       ZNPOS => NPOSM1
       ZNCB  => NCBM1
    ELSE
       ZNPOS => NPOS
       ZNCB  => NCB
    END IF

    IF (method == LG2GP_STD) THEN
       ALLOCATE(zgph(s1,s2,s3,s4))
       zgph(:,:,:,:) = 0.0_dp
       ALLOCATE(num(s1,s2,s4))
       num(:,:,:) = 0.0_dp
    END IF

    IF (method == LG2GP_AVEGT0) THEN
       ALLOCATE(zfngp(s1,s2,s3,s4))
       zfngp(:,:,:,:) = 0.0_dp
    END IF

    ! MASS CONSERVING TRANSFORMATION
    cell_loop: DO jn = 1, NCELL
       i1 = ZNPOS(4, jn) ! x
       i2 = ZNPOS(3, jn) ! z
       i3 = ZNPOS(5, jn) ! y
       !
       i1_a = ZNPOS(4, jn) ! x
       i2_a = ZNPOS(3, jn) ! z
       i3_a = ZNPOS(5, jn) ! y
       !
       ! NOTE: NCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (ZNCB(i1_a, i2_a, i3_a) == 0) THEN
          WRITE(*,*) substr,' WARNING: NCB=0 for ', &
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          weight = AMCELL/PMBOX(i1_a, i2_a, i3_a)
       ELSE
          weight = 1.0_dp
       END IF
       !
       SELECT CASE(method)
       CASE(LG2GP_SUM)
          gp(i1, i2, :, i3) = gp(i1, i2, :, i3) + lg(jn,:) * weight
       CASE(LG2GP_AVE)
          gp(i1, i2, :, i3) = gp(i1, i2, :, i3) + (lg(jn,:) * weight) &
               / REAL(ZNCB(i1_a, i2_a, i3_a), dp)
       CASE(LG2GP_STD)
          zgph(i1, i2, :, i3) = zgph(i1, i2, :, i3) + (lg(jn,:) * weight) &
               / REAL(ZNCB(i1_a, i2_a, i3_a), dp)
          gp(i1, i2, :, i3) = gp(i1, i2, :, i3) + ((lg(jn,:) * weight)**2) &
               / REAL(ZNCB(i1_a, i2_a, i3_a), dp)
          num(i1, i2, i3) = REAL(ZNCB(i1_a, i2_a, i3_a), dp)
       CASE(LG2GP_AVEGT0)
          DO i=1, s3
             IF (lg(jn,i) > 0.0_dp) THEN
                gp(i1, i2, i, i3) = gp(i1, i2, i, i3) + (lg(jn,i) * weight)
                zfngp(i1, i2, i, i3) = zfngp(i1, i2, i, i3) + 1.0_dp
             END IF
          END DO
       CASE DEFAULT
          STATUS = 2
          RETURN     ! ERROR: UNKNOWN METHOD
       END SELECT
    END DO cell_loop

    IF (method == LG2GP_STD) THEN
       ! STANDARD DEVIATION (<.> = average = SUM(.)/N)
       ! SIGMA = SQRT( SUM((xi-<x>)**2)/(N-1) )
       !       = SQRT(| N<x2> - N<x>2 | / (N-1))
       !       = SQRT(| <x2> - <x>2 |) * SQRT(N/(N-1))
       DO i1=1, s1
          DO i2=1, s2
             DO i3=1, s4
                IF (ABS(num(i1,i2,i3)-1.0_dp) > EPS) THEN
                   num(i1,i2,i3) = num(i1,i2,i3) / (num(i1,i2,i3) - 1.0_dp) 
                ELSE
                   num(i1,i2,i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
       !
       DO i=1, s3
          gp(:,:,i,:) = SQRT( ABS(gp(:,:,i,:)-zgph(:,:,i,:)**2) * num(:,:,:) )
       END DO
    END IF

    IF (method == LG2GP_AVEGT0) THEN
       DO i1=1, s1
          DO i2=1, s2
             DO i=1, s3
                DO i3=1, s4
                   IF (zfngp(i1, i2, i, i3) > 0.0_dp) THEN
                      gp(i1, i2, i, i3) = gp(i1, i2, i, i3) &
                           / zfngp(i1, i2, i, i3)
                   ELSE
                      gp(i1, i2, i, i3) = 0.0_dp
                   END IF
                END DO
             END DO
          END DO
       END DO
       IF (PRESENT(fngp)) fngp(:,:,:,:) = zfngp(:,:,:,:)
    END IF

    ! CLEAN MEMORY
    IF (ALLOCATED(zgph)) DEALLOCATE(zgph)
    IF (ALLOCATED(num))  DEALLOCATE(num)
    IF (ALLOCATED(zfngp)) DEALLOCATE(zfngp)

    IF (PRESENT(status)) status = 0

  END SUBROUTINE lg2gp_xzny
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gpsfemis2lgemis_xy(gp, lg, method, gpr, lmcons, gorder, status)

    ! DESCRIPTION: convert gridpoint emission flux to lagrangian emission flux
    !
    ! METHODS:
    !      1: SURFACE LAYER + REST
    !      2: LOWEST CELL IN BOUNDARY LAYER + REST
    !      3: ALL CELLS in BOUNDARY LAYER (weighted) + REST
    !      4: ALL CELLS in BOUNDARY LAYER (unweighted) + REST
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_attila, ONLY: NCELL, NPOS, PMBOX, AMCELL, GNCB, NLEV, KHPBL

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! I/O
    ! global grid point emission flux 
    REAL(dp), DIMENSION(:,:),   INTENT(IN)              :: gp
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(OUT)             :: lg
    ! METHOD
    INTEGER,                    INTENT(IN)              :: method
    ! grid point field with REST MASS (if no CELL is in GRID-BOX)
    REAL(dp), DIMENSION(:,:),   INTENT(INOUT), OPTIONAL :: gpr
    ! for mass conserving transformation
    LOGICAL,                    INTENT(IN),    OPTIONAL :: lmcons
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=2),           INTENT(IN),    OPTIONAL :: gorder
    ! error status
    INTEGER,                    INTENT(OUT),   OPTIONAL :: status

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'gpsfemis2lgemis_xy'
    INTEGER  :: jn, jk
    INTEGER  :: i1_a, i2_a, i3_a    ! INDICES IN ATTILA GP-FIELDS
    INTEGER  :: i1, i2              ! INDICES IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    LOGICAL  :: zlmrest
    INTEGER  :: ilev, flev
    REAL(dp) :: weight
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: gpn
    INTEGER  :: n1, n2

    ! INIT
    IF (PRESENT(status)) status = 1  ! ERROR
    !
    IF (PRESENT(gorder)) THEN
       IF (gorder /= 'xy') THEN
          WRITE(*,*) substr, &
               ': WRONG ORDER! PLEASE RECOMPILE WITH ''#define GENERIC'''
          RETURN
       END IF
    END IF
    !
    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    lg(:) = 0.0_dp
    ! 
    zlmrest = .FALSE.
    IF (PRESENT(gpr)) THEN
       IF (SIZE(gpr) > 0) zlmrest = .TRUE.
    END IF

    ! NUMBER OF CELLS FOR EMISSION
    IF (zlmrest) THEN
       ALLOCATE(gpn(SIZE(gpr,1),SIZE(gpr,2)))
       gpn(:,:) = 0
    END IF

    cell_loop: DO jn = 1, NCELL
       !
       i1 = NPOS(4, jn)  ! x
       i2 = NPOS(5, jn)  ! y
       !
       i1_a = NPOS(4, jn)  ! longitude index
       i2_a = NPOS(3, jn)  ! level index
       i3_a = NPOS(5, jn)  ! latitude index
       !
       ilev = INT(NPOS(3, jn))  ! SHOULD BE = i2_a ???
       !
       ! KHPBL is first layer in free troposphere
       flev = KHPBL(i1_a, i3_a) + 1
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_a, i2_a, i3_a) !!!
       IF (NINT(GNCB(i1_a, i2_a, i3_a)) == 0) THEN
          WRITE(*,*) substr,' WARNING: GNCB=0 for ', &
               jn, i1_a, i2_a, i3_a
          CYCLE
       END IF
       !
       SELECT CASE(METHOD)
       CASE(1)
          ! PUT EMISSION EXCLUSIVELY IN SURFACE LAYER CELL(S)
          IF (ilev /= nlev) CYCLE
          ! NOTE: ALL CELLS IN ONE GRIDPOINT BOX
          weight = 1.0_dp / GNCB(i1_a, i2_a, i3_a)
          !
       CASE(2)
          ! DISTRIBUTE EMISSION AMONG LOWEST BOUNDARY LAYER CELL(S)
          ! NOTES:
          IF (ilev < flev) CYCLE
          !
          !  - EMISSION ONLY, IF NO CELL BELOW THIS LAYER
          IF (ilev < nlev) THEN
             IF (ABS(SUM(GNCB(i1_a,ilev+1:nlev,i3_a))) > EPS) CYCLE
          END IF
          ! NOTE: ALL CELLS ARE IN ONE GRIDPOINT BOX
          weight = 1.0_dp / GNCB(i1_a, i2_a, i3_a)
          !
       CASE(3)
          ! DISTRIBUTE EMISSION AMONG ALL CELLS IN BOUNDARY LAYER
          ! (WEIGHT WITH DISTANCE FROM SURFACE => SOME VERT. DIFF.)
          IF (ilev < flev) CYCLE
          ! NOTE: CELLS ARE LOCATED IN COLUMN BETWEEN flev AND nlev,
          !       HOWEVER, WEIGHTED MASS PER LAYER IS DISTRIBUTED 
          !       EQUALLY AMONG CELLS OF THIS LAYER
          weight = 0.0_dp
          ! SUM ALL CELLS IN BL AND WEIGHT WITH LEVEL
          DO jk=flev,nlev
             weight = weight + REAL(jk, dp)
          END DO
          weight = REAL(ilev, dp) / weight  ! normed weight for layer ilev
          !
! op_sb_20140312+
!!$          weight = weight / GNCB(i1_a, jk, i3_a)
          weight = weight / GNCB(i1_a, ilev, i3_a)
! op_sb_20140312-
          !
       CASE(4)
          ! DISTRIBUTE EMISSION AMONG ALL CELLS IN BOUNDARY LAYER
          ! (UNWEIGHTED => MAX. VERTICAL DIFFUSION)
          IF (ilev < flev) CYCLE
          ! NOTE: CELLS ARE LOCATED IN COLUMN BETWEEN flev AND nlev
          weight = 1.0_dp / SUM(GNCB(i1_a,flev:nlev,i3_a))
          !
       CASE DEFAULT
          ! UNKNOWN METHOD
          IF (PRESENT(status)) status = 2
          RETURN
          !
       END SELECT
       !
       IF (zlmcons) THEN   ! MASS CONSERVING TRANSFORMATION
          ! NOTE: GRIDPOINT EMISSION IS AT LEVEL nlev
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = ( (gp(i1, i2)*PMBOX(i1_a, nlev, i3_a)/AMCELL) &
                    +   gpr(i1, i2) ) * weight
          ELSE
             lg(jn) =   (gp(i1, i2)*PMBOX(i1_a, nlev, i3_a)/AMCELL) &
                    * weight
          END IF
       ELSE               ! NON-MASS CONSERVING TRANSFORMATION
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = (gp(i1, i2) + gpr(i1, i2))
          ELSE
             lg(jn) = gp(i1, i2)
          END IF
       END IF             ! MASS CONSERVING ?
       
    END DO cell_loop

    ! RE-SET REST MASS FLUX
    IF (zlmrest) THEN
       n1 = SIZE(gp, 1) ! x
       n2 = SIZE(gp, 2) ! y
       gridbox_loop: DO i1 = 1, n1  ! x
          DO i2 = 1, n2             ! y
             i1_a = i1
             i3_a = i2
             !
             ! UPDATE NUMBER OF CELLS FOR EMISSION
             ! KHPBL is first layer in free troposphere
             flev = KHPBL(i1_a, i3_a) + 1             
             !
             SELECT CASE(METHOD)
             CASE(1)
                !
                ! PUT EMISSION EXCLUSIVELY IN SURFACE LAYER CELL(S)
                gpn(i1, i2) = NINT(GNCB(i1_a, nlev, i3_a))
                !
             CASE(2)
                !
                ! DISTRIBUTE EMISIION AMONG LOWEST BOUNDARY LAYER CELL(S)
                DO jk=nlev, flev, -1
                   gpn(i1, i2) = NINT(GNCB(i1_a, jk, i3_a))
                   IF (gpn(i1, i2) > 0 ) EXIT
                END DO
                !
             CASE(3, 4)
                !
                ! DISTRIBUTE EMISSION AMONG ALL CELLS IN BOUNDARY LAYER
                DO jk=nlev, flev, -1
                   gpn(i1, i2) = gpn(i1, i2) + NINT(GNCB(i1_a, jk, i3_a))
                END DO
                !
             CASE DEFAULT
                ! UNKNOWN METHOD
                IF (PRESENT(status)) status = 2
                RETURN
                !
             END SELECT
             !
             IF (gpn(i1, i2) == 0) THEN
                ! NOTE: GRIDPOINT EMISSION IS AT LEVEL nlev
                IF (zlmcons) THEN
                   gpr(i1,i2) = gpr(i1,i2) &
                        + gp(i1,i2)*PMBOX(i1_a, nlev, i3_a)/AMCELL
                ELSE
                   gpr(i1,i2) = gpr(i1,i2) + gp(i1,i2)
                END IF
             ELSE
                gpr(i1, i2) = 0.0_dp
             END IF
             !
          END DO
       END DO gridbox_loop
    END IF

    ! CLEAN MEMORY
    IF (ALLOCATED(gpn)) DEALLOCATE(gpn)

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gpsfemis2lgemis_xy
  ! =======================================================================

#endif
! ***************************************************************************
END MODULE messy_attila_tools
! ***************************************************************************
