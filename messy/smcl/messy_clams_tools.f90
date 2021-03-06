! ***************************************************************************
MODULE messy_clams_tools
! ***************************************************************************

  ! MODULE FOR TRANSFORMATIONS/TRANSPOSITIONS BETWEEN GRIDPOINT SPACE
  ! AND LAGRANGIAN SPACE
  ! 
  ! see below info to a mass conserving transformation !!!!
  !
  ! ALL TRANSFORMATIONS/TRANSPOSITIONS HERE
  ! ARE CONSISTENT ON EACH CPU SEPARATELY
  !
  ! Authors: Routines originally for ATTILA:
  !          Patrick Joeckel and Michael Traub, MPICH, Oct 2003
  !          Adopted to ClaMS:
  !          Sabine Brinkop, DLR July 2019

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP
  USE messy_clams_global,       ONLY: dnparts, cmcell 
  USE messy_clams,              ONLY: nlev=>nlev_i  ! ECHAM5-levels 

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: ABS, INDEX, PRESENT, SIZE, SUM, INT &
       , SQRT, REAL, TINY, NINT

  PUBLIC :: DP

  ! ----------- <

  ! index position vectors:
  ! 
  !         LEV (z)
  !         LON (x)
  !         LAT (y)
  !
  !         GP(LON, LEV, LAT)
  !
  ! INDEX OF xyz IN NPOS
  INTEGER, DIMENSION(3), PARAMETER :: lgindx = (/ 1,3,2/) ! POS

  ! ORDER OF x=1, y=2, z=3 IN GP-FIELDS  -> xzy
  INTEGER, DIMENSION(3), PARAMETER :: ovec_c = (/ 1, 3, 2 /) 
  REAL(DP),              PARAMETER, PUBLIC :: EPS = TINY(0.0_DP) ! zero

  PUBLIC :: cl2gp             ! CLAMS Lagrange  -> GRIDPOINT
  PUBLIC :: gpsfemis2clemis   ! GRIDPOINT SURFACE EMISSION -> LAGRANGIAN EMISS
  PUBLIC :: gp2cl             ! GRIDPOINT -> LAGRANGE

  INTERFACE gpsfemis2clemis
     MODULE PROCEDURE gpsfemis2clemis_xy
  END INTERFACE

  INTERFACE gp2cl
     MODULE PROCEDURE gp2cl_xzny
     MODULE PROCEDURE gp2cl_xzy
     MODULE PROCEDURE gp2cl_xy
  END INTERFACE

  INTERFACE cl2gp
     MODULE PROCEDURE cl2gp_xzy
     MODULE PROCEDURE cl2gp_xzny
  END INTERFACE

  ! FOR cl2gp
  INTEGER, PARAMETER, PUBLIC :: CL2GP_SUM = 1
  INTEGER, PARAMETER, PUBLIC :: CL2GP_AVE = 2
  INTEGER, PARAMETER, PUBLIC :: CL2GP_STD = 3
  INTEGER, PARAMETER, PUBLIC :: CL2GP_AVEGT0 = 4 ! average all elements > 0

CONTAINS
     
  ! =======================================================================
  SUBROUTINE cl2gp_xzy(lg, gp, method, gorder, lmcons  &
       , status, fngp)
    
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
    ! NOTE NOTE NOTE NOTE NOTE NOTE **********************************
    ! 2) MASS CONSERVING TRANSFORMATION:
    !    Please note: CLaMS parcels have no mass, therefore
    !    the CLaMS parcel mass is set to 1.
    !    (lmass_cons = .true.)
    !    THE MASS OF ALL LAGRANGIAN CELLS IS RELEASED INTO THE
    !    HOSTING GRIDPOINT BOX
    !     - SUM
    !     - MASS WEIGHTED AVERAGE
    !     - STANDARD DEVIATION
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_clams, ONLY: POS, PC_G, pmbox
    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! I/O
    ! lagrange vetor
    REAL(dp), DIMENSION(:),     INTENT(IN)             :: lg
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT)            :: gp
    ! method
    INTEGER,                    INTENT(IN)             :: method
    LOGICAL,                    INTENT(IN), OPTIONAL   :: lmcons
    ! order of x,y,z - indices in gridpoint fields gp, gpr
    CHARACTER(LEN=3),           INTENT(IN),   OPTIONAL :: gorder
    ! error status
    INTEGER,                    INTENT(OUT),  OPTIONAL :: status
    ! number of 'events'
    REAL(dp), DIMENSION(:,:,:), INTENT(OUT),  OPTIONAL :: fngp

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cl2gp_xzy'
    INTEGER  :: jn
    INTEGER  :: i1_c, i2_c, i3_c    ! INDEX IN CLAMS GP-FIELDS
    INTEGER  :: i1, i2, i3          ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    INTEGER  :: s1, s2, s3
    REAL(dp) :: weight
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zgph, num, zfngp
    REAL(dp), DIMENSION(:,:),   POINTER     :: ZNPOS => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER     :: ZNCB  => NULL()

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

    IF (PRESENT(lmcons)) THEN
       zlmcons = lmcons
    ELSE
       zlmcons = .FALSE. ! DEFAULT: NON-MASS CONSERVING TRAFO
    ENDIF
    !
    ZNPOS => POS
    ZNCB  => PC_G
  
    IF (method == CL2GP_STD) THEN
       ALLOCATE(zgph(s1,s2,s3))
       zgph(:,:,:) = 0.0_dp
       ALLOCATE(num(s1,s2,s3))
       num(:,:,:) = 0.0_dp
    END IF

    IF (method == CL2GP_AVEGT0) THEN
       ALLOCATE(zfngp(s1,s2,s3))
       zfngp(:,:,:) = 0.0_dp
    END IF

    ! MASS CONSERVING TRANSFORMATION
    cell_loop: DO jn = 1, dnparts
       i1 = ZNPOS(jn,1) ! x
       i2 = ZNPOS(jn,2) ! z
       i3 = ZNPOS(jn,3) ! y
       !
       i1_c = ZNPOS(jn,1) ! x
       i2_c = ZNPOS(jn,2) ! z
       i3_c = ZNPOS(jn,3) ! y
       !
       ! NOTE: NCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_c, i2_c, i3_c) !!!
       IF (ZNCB(i1_c, i2_c, i3_c) == 0) THEN
          WRITE(*,*) substr,' WARNING: PC_G=0 for ', &
               jn, i1_c, i2_c, i3_c
          CYCLE
       END IF
       IF (zlmcons) THEN
          weight = CMCELL/PMBOX(i1_c, i2_c, i3_c)
       ELSE
          weight = 1.0_dp
       END IF
       !
       SELECT CASE(method)
       CASE(CL2GP_SUM)
          gp(i1, i2, i3) = gp(i1, i2, i3) + lg(jn)* weight
       CASE(CL2GP_AVE)
          gp(i1, i2, i3) = gp(i1, i2, i3) + (lg(jn)* weight) &
               / REAL(ZNCB(i1_c, i2_c, i3_c), dp)
       CASE(CL2GP_STD)
          zgph(i1, i2, i3) = zgph(i1, i2, i3) + (lg(jn)* weight) &
               / REAL(ZNCB(i1_c, i2_c, i3_c), dp)
          gp(i1, i2, i3) = gp(i1, i2, i3) + ((lg(jn)* weight)**2) &
               / REAL(ZNCB(i1_c, i2_c, i3_c), dp)
          num(i1, i2, i3) = REAL(ZNCB(i1_c, i2_c, i3_c), dp)
       CASE(CL2GP_AVEGT0)
          IF (lg(jn) > 0.0_dp) THEN
             gp(i1, i2, i3) = gp(i1, i2, i3) + (lg(jn)* weight)
             zfngp(i1, i2, i3) = zfngp(i1, i2, i3) + 1.0_dp
          END IF
       CASE DEFAULT
          STATUS = 2
          RETURN     ! ERROR: UNKNOWN METHOD
       END SELECT
    END DO cell_loop

    IF (method == CL2GP_STD) THEN
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

    IF (method == CL2GP_AVEGT0) THEN
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

  END SUBROUTINE cl2gp_xzy
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE cl2gp_xzny(lg, gp, method, gorder, lmcons, status, fngp)
    
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

    USE messy_clams, ONLY: POS, PC_G, pmbox

    IMPLICIT NONE

    INTRINSIC :: ALLOCATED

    ! I/O
    ! lagrange vetor
    REAL(dp), DIMENSION(:,:),     INTENT(IN)             :: lg
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(dp), DIMENSION(:,:,:,:), INTENT(OUT)            :: gp
    ! method
    INTEGER,                      INTENT(IN)             :: method
    CHARACTER(LEN=4),             INTENT(IN),   OPTIONAL :: gorder
    ! error status
    INTEGER,                      INTENT(OUT),  OPTIONAL :: status
    LOGICAL,                      INTENT(IN),   OPTIONAL :: lmcons
    ! number of 'events'
    REAL(dp), DIMENSION(:,:,:,:), INTENT(OUT),  OPTIONAL :: fngp

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cl2gp_xzny'
    INTEGER  :: jn
    INTEGER  :: i1_c, i2_c, i3_c    ! INDEX IN CLMAS GP-FIELDS
    INTEGER  :: i1, i2, i3          ! INDEX IN EXTERNAL GP-FIELDS
    LOGICAL  :: zlmcons
    INTEGER  :: s1, s2, s3, s4, i
    REAL(dp) :: weight
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: zgph
    REAL(dp), DIMENSION(:,:,:),   ALLOCATABLE :: num
    REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: zfngp
    REAL(dp), DIMENSION(:,:),     POINTER     :: ZNPOS => NULL()
    REAL(dp), DIMENSION(:,:,:),   POINTER     :: ZNCB  => NULL()

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
    ZNPOS => POS
    ZNCB  => PC_G
 
    IF (method == CL2GP_STD) THEN
       ALLOCATE(zgph(s1,s2,s3,s4))
       zgph(:,:,:,:) = 0.0_dp
       ALLOCATE(num(s1,s2,s4))
       num(:,:,:) = 0.0_dp
    END IF

    IF (method == CL2GP_AVEGT0) THEN
       ALLOCATE(zfngp(s1,s2,s3,s4))
       zfngp(:,:,:,:) = 0.0_dp
    END IF

    ! MASS CONSERVING TRANSFORMATION
    cell_loop: DO jn = 1, dnparts

       if (ZNPOS(jn,1) .lt. 0 .or. ZNPOS(jn,2) .lt. 0 .or. &
             ZNPOS(jn,3) .lt. 0 ) then
    !  
       else
       i1 = ZNPOS(jn,1) ! x
       i2 = ZNPOS(jn,2) ! z
       i3 = ZNPOS(jn,3) ! y
       !
       i1_c = ZNPOS(jn,1) ! x
       i2_c = ZNPOS(jn,2) ! z
       i3_c = ZNPOS(jn,3) ! y
       !
       ! NOTE: NCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_c, i2_c, i3_c) !!!
       IF (ZNCB(i1_c, i2_c, i3_c) == 0.) THEN
          WRITE(*,*) substr,' WARNING: PC_G=0 for ', &
               jn, i1_c, i2_c, i3_c
          CYCLE
       END IF

       IF (zlmcons) THEN
          weight = CMCELL/PMBOX(i1_c, i2_c, i3_c)
       ELSE
          weight = 1.0_dp
       END IF
       
       !
       SELECT CASE(method)
       CASE(CL2GP_SUM)
          gp(i1, i2, :, i3) = gp(i1, i2, :, i3) + lg(jn,:) * weight
       CASE(CL2GP_AVE)
          gp(i1, i2, :, i3) = gp(i1, i2, :, i3) + (lg(jn,:)  * weight) &
               / ZNCB(i1_c, i2_c, i3_c)
       CASE(CL2GP_STD)
          zgph(i1, i2, :, i3) = zgph(i1, i2, :, i3) + (lg(jn,:) * weight) &
               / ZNCB(i1_c, i2_c, i3_c)
          gp(i1, i2, :, i3) = gp(i1, i2, :, i3) + (lg(jn,:) * weight)**2 &
               / ZNCB(i1_c, i2_c, i3_c)
          num(i1, i2, i3) = ZNCB(i1_c, i2_c, i3_c)
       CASE(CL2GP_AVEGT0)
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

     ENDIF
    END DO cell_loop

    IF (method == CL2GP_STD) THEN
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

    IF (method == CL2GP_AVEGT0) THEN
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

  END SUBROUTINE cl2gp_xzny

  ! =======================================================================

 SUBROUTINE gpsfemis2clemis_xy(gp, lg, method, gpr, lmcons, gorder, status)

    ! DESCRIPTION: convert gridpoint emission flux to lagrangian emission flux
    !
    ! METHODS:
    !      1: SURFACE LAYER + REST
    !      2: LOWEST CELL IN BOUNDARY LAYER + REST
    !      3: ALL CELLS in BOUNDARY LAYER (weighted) + REST
    !      4: ALL CELLS in BOUNDARY LAYER (unweighted) + REST
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_clams_global,       ONLY: dnparts 
    USE messy_clams,              ONLY: POS, pc_g, khpbl, pmbox

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
    CHARACTER(LEN=*), PARAMETER :: substr = 'gpsfemis2clemis_xy'
    INTEGER  :: jn, jk
    INTEGER  :: i1_c, i2_c, i3_c    ! INDICES IN CLaMS GP-FIELDS
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
               ': WRONG ORDER!'
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
    cell_loop_emis: DO jn = 1, dnparts
       !
       i1_c = INT(POS(jn,1))  ! longitude index
       i2_c = INT(POS(jn,2))  ! level index
       i3_c = INT(POS(jn,3))  ! latitude index
       !
       ilev = INT(POS(jn,2)) 
       !
       ! KHPBL is first layer in free troposphere
       if (pos(jn,1) .ge. 1. .and. pos(jn,2) .ge. 1.      &
                                 .and. pos(jn,3) .ge. 1.) then ! op_sb_20190802
        flev = INT(KHPBL(i1_c, i3_c)) + 1
    
       !
       ! NOTE: pc_g MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_c, i2_c, i3_c) !!!
       IF (NINT(pc_g(i1_c, i2_c, i3_c)) == 0) THEN
          WRITE(*,*) substr,' WARNING: pc_g=0 for ', &
               jn, i1_c, i2_c, i3_c
          CYCLE
       END IF
       
       !
       SELECT CASE(METHOD)
       CASE(1)
          ! PUT EMISSION EXCLUSIVELY IN SURFACE LAYER CELL(S)
          IF (ilev /= nlev) CYCLE
          ! NOTE: ALL CELLS IN ONE GRIDPOINT BOX
          weight = 1.0_dp / pc_g(i1_c, i2_c, i3_c)
          !
       CASE(2)
          ! DISTRIBUTE EMISSION AMONG LOWEST BOUNDARY LAYER CELL(S)
          ! NOTES:
          IF (ilev < flev) CYCLE
          !
          !  - EMISSION ONLY, IF NO CELL BELOW THIS LAYER
          IF (ilev < nlev) THEN
             IF (ABS(SUM(pc_g(i1_c,ilev+1:nlev,i3_c))) > EPS) CYCLE
          END IF
          ! NOTE: ALL CELLS ARE IN ONE GRIDPOINT BOX
          weight = 1.0_dp / pc_g(i1_c, i2_c, i3_c)
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
          weight = weight / pc_g(i1_c, ilev, i3_c)
          !
       CASE(4)
          ! DISTRIBUTE EMISSION AMONG ALL CELLS IN BOUNDARY LAYER
          ! (UNWEIGHTED => MAX. VERTICAL DIFFUSION)
          IF (ilev < flev) CYCLE
          ! NOTE: CELLS ARE LOCATED IN COLUMN BETWEEN flev AND nlev
          weight = 1.0_dp / SUM(pc_g(i1_c,flev:nlev,i3_c))
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
             lg(jn) = ( (gp(i1_c, i3_c)*PMBOX(i1_c, nlev, i3_c)/CMCELL) &
                    +   gpr(i1_c, i3_c) ) * weight
             lg(jn) = ( (gp(i1_c, i3_c)) + gpr(i1_c, i3_c) ) * weight
          ELSE
               lg(jn) =   (gp(i1_c, i3_c)*PMBOX(i1_c, nlev, i3_c)/CMCELL) &
                   * weight
            
          END IF
       ELSE               ! NON-MASS CONSERVING TRANSFORMATION
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = (gp(i1_c, i3_c) + gpr(i1_c, i3_c))
          ELSE
             lg(jn) = gp(i1_c, i3_c)
          END IF
       END IF             ! MASS CONSERVING ?
       endif   !  pos > 0?
    END DO cell_loop_emis

    ! RE-SET REST MASS FLUX
    IF (zlmrest) THEN
       n1 = SIZE(gp, 1) ! x
       n2 = SIZE(gp, 2) ! y
       gridbox_loop: DO i1 = 1, n1  ! x
          DO i2 = 1, n2             ! y
             i1_c = i1
             i3_c = i2
             !
             ! UPDATE NUMBER OF CELLS FOR EMISSION
             ! KHPBL is first layer in free troposphere
             flev = KHPBL(i1_c, i3_c) + 1             
             !
             SELECT CASE(METHOD)
             CASE(1)
                !
                ! PUT EMISSION EXCLUSIVELY IN SURFACE LAYER CELL(S)
                gpn(i1, i2) = NINT(pc_g(i1_c, nlev, i3_c))
                !
             CASE(2)
                !
                ! DISTRIBUTE EMISIION AMONG LOWEST BOUNDARY LAYER CELL(S)
                DO jk=nlev, flev, -1
                   gpn(i1, i2) = NINT(pc_g(i1_c, jk, i3_c))
                   IF (gpn(i1, i2) > 0 ) EXIT
                END DO
                !
             CASE(3, 4)
                !
                ! DISTRIBUTE EMISSION AMONG ALL CELLS IN BOUNDARY LAYER
                DO jk=nlev, flev, -1
                   gpn(i1, i2) = gpn(i1, i2) + NINT(pc_g(i1_c, jk, i3_c))
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
                        + gp(i1,i2)*PMBOX(i1_c, nlev, i3_c)/CMCELL
                
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

  END SUBROUTINE gpsfemis2clemis_xy
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gp2cl_xzy(gp, lg, gpr, lmcons, gorder, status)
    
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_clams, ONLY: POS, pc_g, pmbox

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
    CHARACTER(LEN=*), PARAMETER :: substr = 'gp2cl_xzy'
    INTEGER  :: jn
    INTEGER  :: i1_c, i2_c, i3_c    ! INDEX IN ATTILA GP-FIELDS
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
               ': WRONG ORDER!'
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

    cell_loop: DO jn = 1, dnparts
       i1 = INT(POS(jn,1)) ! x
       i2 = INT(POS(jn,2)) ! z
       i3 = INT(POS(jn,3)) ! y
       !
       i1_c = INT(POS(jn,1)) ! x
       i2_c = INT(POS(jn,2)) ! z
       i3_c = INT(POS(jn,3)) ! y
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_c, i2_c, i3_c) !!!
       IF (NINT(pc_g(i1_c, i2_c, i3_c)) == 0) THEN
          WRITE(*,*) substr,' WARNING: pc_g=0 for ',&
               jn, i1_c, i2_c, i3_c
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = ( (gp(i1, i2, i3)*PMBOX(i1_c, i2_c, i3_c)/CMCELL)   &
                  +   gpr(i1, i2, i3) ) / pc_g(i1_c, i2_c, i3_c)
          ELSE
             lg(jn) = gp(i1, i2, i3)*PMBOX(i1_c, i2_c, i3_c) &
                  /(pc_g(i1_c, i2_c, i3_c)*CMCELL)
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
                   weight = PMBOX(i1, i2, i3)/CMCELL  ! xzy
                ELSE
                   weight = 1.0_dp
                END IF
                !
                IF (NINT(pc_g(i1,i2,i3)) == 0) THEN
                   gpr(i1,i2,i3) = gpr(i1,i2,i3) + gp(i1,i2,i3) * weight
                ELSE
                   gpr(i1,i2,i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
    END IF

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gp2cl_xzy
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gp2cl_xzny(gp, lg, gpr, lmcons, gorder, status)
    
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_clams, ONLY: POS, pc_g, pmbox
  
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
    CHARACTER(LEN=*), PARAMETER :: substr = 'gp2cl_xzny'
    INTEGER  :: jn
    INTEGER  :: i1_c, i2_c, i3_c    ! INDEX IN ATTILA GP-FIELDS
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
               ': WRONG ORDER!'
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

    cell_loop: DO jn = 1, dnparts
       i1 = INT(POS(jn,1)) ! x
       i2 = INT(POS(jn,2)) ! z
       i3 = INT(POS(jn,3)) ! y
       !
       i1_c = INT(POS(jn,1)) ! x
       i2_c = INT(POS(jn,2)) ! z
       i3_c = INT(POS(jn,3)) ! y
       !
       ! NOTE: GNCB MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_c, i2_c, i3_c) !!!
       IF (NINT(pc_g(i1_c, i2_c, i3_c)) == 0) THEN
          WRITE(*,*) substr,' WARNING: pc_g=0 for ',&
               jn, i1_c, i2_c, i3_c
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn,:) = ( (gp(i1, i2, :, i3)*PMBOX(i1_c, i2_c, i3_c)/CMCELL) &
                  +   gpr(i1, i2, :, i3) ) / pc_g(i1_c, i2_c, i3_c)
          ELSE
             lg(jn,:) = gp(i1, i2, :, i3)*PMBOX(i1_c, i2_c, i3_c) &
                  /(pc_g(i1_c, i2_c, i3_c)*CMCELL)
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
                   weight = PMBOX(i1, i2, i3)/CMCELL  ! xzy
                ELSE
                   weight = 1.0_dp
                END IF
                !
                IF (NINT(pc_g(i1,i2,i3)) == 0) THEN
                   gpr(i1,i2,:,i3) = gpr(i1,i2,:,i3) + gp(i1,i2,:,i3) * weight
                ELSE
                   gpr(i1,i2,:,i3) = 0.0_dp
                END IF
             END DO
          END DO
       END DO
    END IF

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gp2cl_xzny
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE gp2cl_xy(gp, lg, gpr, lmcons, klev, gorder, llev, status)

    ! DESCRIPTION: see gp2cl_xyz
    !
    ! Author : Patrick Joeckel, MPICH, Oct  2003

    USE messy_clams, ONLY: POS, pc_g, pmbox
  
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
    CHARACTER(LEN=*), PARAMETER :: substr = 'gp2cl_xy'
    INTEGER  :: jn
    INTEGER  :: i1_c, i2_c     ! INDEX IN CLaMS GP-FIELDS
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
               ': WRONG ORDER!'
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
       ALLOCATE(llev(dnparts))
       llev = .FALSE.
    END IF

    cell_loop: DO jn = 1, dnparts
       IF (INT(POS(jn,2)) /= zklev) CYCLE  ! z
       !
       i1 = INT(POS(jn,1)) ! x
       i2 = INT(POS(jn,3)) ! y
       !
       i1_c = INT(POS(jn,1)) ! x
       i2_c = INT(POS(jn,3)) ! y
       !
       ! NOTE: pc_g MUST NOT BE ZERO, SINCE AT LEAST LAGRANGIAN
       !       BOX WITH INDEX jn IS IN GRIDPOINT BOX WITH
       !       INDICES (i1_c, zklev, i3_c) !!!
       IF (NINT(pc_g(i1_c, zklev, i2_c)) == 0) THEN
          WRITE(*,*) substr,' WARNING: pc_g=0 for ', &
               jn, i1_c, zklev, i2_c
          CYCLE
       END IF
       !
       IF (zlmcons) THEN
          IF (zlmrest) THEN
             ! ADD SAVED GRIDPOINT-REST
             lg(jn) = ( (gp(i1, i2)*PMBOX(i1_c, zklev, i2_c)/CMCELL)   &
                  +   gpr(i1, i2) ) / pc_g(i1_c, zklev, i2_c)

          ELSE
             lg(jn) = gp(i1, i2)*PMBOX(i1_c, zklev, i2_c) &
                  /(pc_g(i1_c, zklev, i2_c)*CMCELL)
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
                weight = PMBOX(i1, zklev, i2)/CMCELL ! xzy
             ELSE
                weight = 1.0_dp
             END IF
             !
             IF (NINT(pc_g(i1,zklev,i2)) == 0) THEN
                gpr(i1,i2) = gpr(i1,i2) + gp(i1,i2) * weight
             ELSE
                gpr(i1, i2) = 0.0_dp
             ENDIF
          END DO
       END DO
    END IF

    IF (PRESENT(status)) status = 0

  END SUBROUTINE gp2cl_xy
 
! ***************************************************************************
END MODULE messy_clams_tools
! ***************************************************************************
