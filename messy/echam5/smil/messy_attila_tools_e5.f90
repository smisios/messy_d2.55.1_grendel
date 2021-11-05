! ***************************************************************************
MODULE messy_attila_tools_e5
! ***************************************************************************

  ! MODULE FOR TRANSFORMATIONS/TRANSPOSITIONS BETWEEN GRIDPOINT SPACE
  ! AND LAGRANGIAN SPACE, BOTH IN ECHAM5 DECOMPOSITION
  !
  ! Author: Patrick Joeckel, MPICH, Oct 2003
  !

  USE messy_attila_tools,       ONLY: LG2GP_SUM,  LG2GP_AVE,  LG2GP_STD &
                                    , LG2GP_AVEGT0
  USE messy_attila,             ONLY: GNCB, DP
  
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: GNCB
  PUBLIC :: gp2lg_e5
  PUBLIC :: lg2gp_e5
  PUBLIC :: LG2GP_SUM,  LG2GP_AVE,  LG2GP_STD, LG2GP_AVEGT0
  PUBLIC :: gpsfemis2lgemis_e5
  PUBLIC :: lggpte2lgte_e5

  INTERFACE gp2lg_e5
     MODULE PROCEDURE gp2lg_xy_e5
     MODULE PROCEDURE gp2lg_xyz_e5
     MODULE PROCEDURE gp2lg_xyzn_e5
  END INTERFACE

  INTERFACE lg2gp_e5
     MODULE PROCEDURE lg2gp_xyz_e5
     MODULE PROCEDURE lg2gp_xyzn_e5
  END INTERFACE

  INTERFACE lggpte2lgte_e5
     MODULE PROCEDURE lggpte2lgte_xyz_e5
     MODULE PROCEDURE lggpte2lgte_xyzn_e5
  END INTERFACE

  INTRINSIC :: PRESENT, ABS, MERGE, ASSOCIATED

CONTAINS

!----------------------------------------------------------------------------
  SUBROUTINE gp2lg_xy_e5(gpl, lgl, gprl, lmcons, klev, llev)

    USE messy_attila_tools,      ONLY: gp2lg
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl, M_LOC
    USE messy_main_mpi_bi,       ONLY: finish

    IMPLICIT NONE

    ! I/O
    ! local GRIDPOINT FIELD (channel)
    REAL(DP), DIMENSION(:,:), POINTER              :: gpl  ! INTENT(IN)
    ! local LAGRANGIAN FIELD (channel)
    REAL(DP), DIMENSION(:),   POINTER              :: lgl  ! INTENT(OUT)
    ! local GRIDPOINT FIELD (channel) WITH REST OF gpr
    REAL(DP), DIMENSION(:,:), POINTER,    OPTIONAL :: gprl ! INTENT(INOUT)
    ! SWITCH FOR MASS CONSERVATION
    LOGICAL,              INTENT(IN), OPTIONAL :: lmcons
    ! WHICH LEVEL (DEFAULT: SURFACE)
    INTEGER,              INTENT(IN), OPTIONAL :: klev
    ! TRUE IF BOX IS AT THIS LEVEL
    LOGICAL, DIMENSION(:), POINTER,   OPTIONAL :: llev

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'gp2lg_xy_e5'
    REAL(DP), DIMENSION(:,:), POINTER  :: gpg  => NULL() ! globalized field
    REAL(DP), DIMENSION(:,:), POINTER  :: gprg => NULL() ! globalized rest
    INTEGER                        :: status

    ! INIT
    IF (PRESENT(llev)) THEN
       IF (ASSOCIATED(llev)) THEN
          DEALLOCATE(llev)
          NULLIFY(llev)
       END IF
    END IF

    ! GLOBAL FIELD ON ALL CPUs
    CALL trp_gpdc_gpgl(1, gpl, gpg)
    ! GLOBAL REST ON ALL CPUs
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(1, gprl, gprg)

    ! TRANSFER FIELD INTO LAGRANGIAN SPACE (MASS CONSERVING) AND SAVE REST
    IF (PRESENT(gprl)) THEN
       CALL gp2lg(gpg, lgl, gpr=gprg &
            , lmcons=lmcons, klev=klev, gorder='xy', llev=llev, status=status)
    ELSE
       CALL gp2lg(gpg, lgl           &
            , lmcons=lmcons, klev=klev, gorder='xy', llev=llev, status=status)
    END IF
    IF (status /= 0) &
         CALL finish(substr,' Error in gp2lg!')

    ! GET REST FROM ALL CPUs; PUT BACK INTO CHANNEL
    ! NOTE: 
    !   - AVE, SINCE ALL CPUs calculate (gp2lg) the same REST with GNCB !!!
    !   - LOC even better
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(-1, gprl, gprg, M_LOC)

    ! CLEAN UP
    IF (ASSOCIATED(gpg))  DEALLOCATE(gpg)
    IF (ASSOCIATED(gprg)) DEALLOCATE(gprg)

  END SUBROUTINE gp2lg_xy_e5
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  SUBROUTINE gp2lg_xyz_e5(gpl, lgl, gprl, lmcons)

    USE messy_attila_tools,      ONLY: gp2lg
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl, M_LOC
    USE messy_main_mpi_bi,       ONLY: finish

    IMPLICIT NONE

    ! I/O
    ! local GRIDPOINT FIELD (channel)
    REAL(DP), DIMENSION(:,:,:), POINTER              :: gpl  ! INTENT(IN)
    ! local LAGRANGIAN FIELD (channel)
    REAL(DP), DIMENSION(:),     POINTER              :: lgl  ! INTENT(OUT)
    ! local GRIDPOINT FIELD (channel) WITH REST OF gpr
    REAL(DP), DIMENSION(:,:,:), POINTER,    OPTIONAL :: gprl ! INTENT(INOUT)
    ! SWITCH FOR MASS CONSERVATION
    LOGICAL,                INTENT(IN), OPTIONAL :: lmcons

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER      :: substr = 'gp2lg_xyz_e5'
    REAL(DP), DIMENSION(:,:,:), POINTER  :: gpg  => NULL() ! globalized field
    REAL(DP), DIMENSION(:,:,:), POINTER  :: gprg => NULL() ! globalized rest
    INTEGER                          :: status

    ! GLOBAL FIELD ON ALL CPUs
    CALL trp_gpdc_gpgl(1, gpl, gpg)
    ! GLOBAL REST ON ALL CPUs
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(1, gprl, gprg)

    ! TRANSFER FIELD INTO LAGRANGIAN SPACE (MASS CONSERVING) AND SAVE REST
    IF (PRESENT(gprl)) THEN
       CALL gp2lg(gpg, lgl, gpr=gprg &
            , lmcons=lmcons, gorder='xzy', status=status)
    ELSE
       CALL gp2lg(gpg, lgl           &
            , lmcons=lmcons, gorder='xzy', status=status)
    END IF
    IF (status /= 0) &
         CALL finish(substr,' Error in gp2lg!')

    ! GET REST FROM ALL CPUs; PUT BACK INTO CHANNEL
    ! NOTE:
    !   - AVE, SINCE ALL CPUs calculate (gp2lg) the same REST with GNCB !!!
    !   - LOC, even better
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(-1, gprl, gprg, M_LOC)

    ! CLEAN UP
    IF (ASSOCIATED(gpg))  DEALLOCATE(gpg)
    IF (ASSOCIATED(gprg)) DEALLOCATE(gprg)

  END SUBROUTINE gp2lg_xyz_e5
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  SUBROUTINE gp2lg_xyzn_e5(gpl, lgl, gprl, lmcons)

    USE messy_attila_tools,      ONLY: gp2lg
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl, M_LOC
    USE messy_main_mpi_bi,       ONLY: finish

    IMPLICIT NONE

    ! I/O
    ! local GRIDPOINT FIELD (channel)
    REAL(DP), DIMENSION(:,:,:,:), POINTER              :: gpl  ! INTENT(IN)
    ! local LAGRANGIAN FIELD (channel)
    REAL(DP), DIMENSION(:,:),     POINTER              :: lgl  ! INTENT(OUT)
    ! local GRIDPOINT FIELD (channel) WITH REST OF gpr
    REAL(DP), DIMENSION(:,:,:,:), POINTER,    OPTIONAL :: gprl ! INTENT(INOUT)
    ! SWITCH FOR MASS CONSERVATION
    LOGICAL,                  INTENT(IN), OPTIONAL :: lmcons

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER        :: substr = 'gp2lg_xyzn_e5'
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: gpg  => NULL() ! globalized field
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: gprg => NULL() ! globalized rest
    INTEGER                            :: status

    ! GLOBAL FIELD ON ALL CPUs
    CALL trp_gpdc_gpgl(1, gpl, gpg)
    ! GLOBAL REST ON ALL CPUs
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(1, gprl, gprg)

    ! TRANSFER FIELD INTO LAGRANGIAN SPACE (MASS CONSERVING) AND SAVE REST
    IF (PRESENT(gprl)) THEN
       CALL gp2lg(gpg, lgl, gpr=gprg &
            , lmcons=lmcons, gorder='xzny', status=status)
    ELSE
       CALL gp2lg(gpg, lgl           &
            , lmcons=lmcons, gorder='xzny', status=status)
    END IF
    IF (status /= 0) &
         CALL finish(substr,' Error in gp2lg!')

    ! GET REST FROM ALL CPUs; PUT BACK INTO CHANNEL
    ! NOTE:
    !   - AVE, SINCE ALL CPUs calculate (gp2lg) the same REST with GNCB !!!
    !   - LOC, even better
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(-1, gprl, gprg, M_LOC)

    ! CLEAN UP
    IF (ASSOCIATED(gpg))  DEALLOCATE(gpg)
    IF (ASSOCIATED(gprg)) DEALLOCATE(gprg)

  END SUBROUTINE gp2lg_xyzn_e5
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  SUBROUTINE lg2gp_xyz_e5(lgl, gpl, method, lmcons, ltm1 &
       , fill_value, fill_field)

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, ngpblks, nlon, ngl 
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl, M_SUM
    USE messy_main_mpi_bi,       ONLY: finish, message
    USE messy_attila_mem_e5,     ONLY: SGNCB, SGNCBM1
    ! MESSy
    USE messy_main_timer,    ONLY: lstart !, lresume
    USE messy_attila_tools,  ONLY: lg2gp, EPS &
                                 , LG2GP_SUM, LG2GP_AVE, LG2GP_STD &
                                 , LG2GP_AVEGT0

    IMPLICIT NONE

    ! I/O
    ! lagrange vetor
    REAL(DP), DIMENSION(:),     INTENT(IN)              :: lgl
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(DP), DIMENSION(:,:,:), POINTER                 :: gpl
    ! method
    INTEGER,                INTENT(IN)              :: method
    ! for mass conserving transformation
    LOGICAL,                INTENT(IN),    OPTIONAL :: lmcons
    ! use positions of CELLS at t-1 ?
    LOGICAL,                INTENT(IN),    OPTIONAL :: ltm1
    ! optional fill-value where no CELL is located
    REAL(DP),                   INTENT(IN),    OPTIONAL :: fill_value
    ! optional GP-field to fill in where no CELL is located
    REAL(DP), DIMENSION(:,:,:), INTENT(IN),    OPTIONAL :: fill_field

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lg2gp_xyz_e5'
    INTEGER :: status
    REAL(DP), DIMENSION(:,:,:), POINTER :: gp     => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: zgpl   => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: ZSGNCB => NULL()
    ! number of 'events'
    REAL(DP), DIMENSION(:,:,:), POINTER :: zfngp  => NULL()
    INTEGER                         :: i1,i2,i3

    IF (PRESENT(ltm1)) THEN
       IF ((ltm1).AND.lstart) THEN
          CALL message(substr,' WARNING: t-1 not defined at t=0;'&
               &//' transformation with ltm1=.true. not possible!')
          RETURN ! t-1 not defined at t=0
       END IF
       !
       IF (ltm1) THEN
          ZSGNCB => SGNCBM1
       ELSE
          ZSGNCB => SGNCB 
       END IF
    ELSE
       ZSGNCB => SGNCB 
    END IF

    ALLOCATE(gp(nlon, nlev, ngl))

    SELECT CASE(method)
    CASE(LG2GP_AVEGT0)
       ALLOCATE(zfngp(nlon, nlev, ngl))
       CALL lg2gp(lgl, gp, method, lmcons=lmcons   &
            , gorder='xzy', ltm1=ltm1, status=status, fngp=zfngp)
    CASE DEFAULT
       CALL lg2gp(lgl, gp, LG2GP_SUM, lmcons=lmcons   &
            , gorder='xzy', ltm1=ltm1, status=status)
    END SELECT
    IF (status /= 0) &
         CALL finish(substr,' Error in lg2gp!')

    ! NOTE:
    !
    ! THE RESULTS OF lg2gp ARE CONSISTENT FOR ALL CELLS ON 1 CPU !!!
    !
    ! (z)lmcons            .FALSE.   | .TRUE.
    ! => WEIGHT OF x     UNWEIGHTED  | MASS WEIGHTED (AMCELL/PMBOX)
    !         -----------------------------------------------------
    ! gp sum | SUM          SUM(x)
    ! gp ave | AVERAGE      <x> = SUM(x)/NCB
    ! gp std | STD.DEV      SQRT(| <x2> - <x>2 |) * SQRT(NCB/(NCB-1))
    !
    ! THESE 'PER CPU'-RESULTS NEED TO BE TRANSPOSED/TRANSFORMED 
    ! BACK INTO CONSISTENT FIELDS ON ALL CPUs
    !
    ! THE FOLLOWING CALCULATIONS ARE  INDEPENDENT OF MASS CONSERVATION
    ! ((z)lmcons), SINCE THE WEIGHT AMCELL/PMBOX IS CONSTANT FOR
    ! ALL CELLS IN ONE BOX


    SELECT CASE(METHOD)
    CASE(LG2GP_SUM)
       ! SUM
       ! SUMMATION OVER ALL CELLS -> SUMMATION OVER ALL CPUs (CPU_SUM)
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       !
    CASE(LG2GP_AVE)
       ! AVERAGE
       ! AVERAGE OVER ALL CELLS -> METHOD 1: CPU_SUM(gp_sum)/GNCB
       !                           METHOD 2: CPU_SUM(gp_ave*NCB)/GNCB
       !            (NCB not available in CHANNEL DECOMP.
       !             gp_ave need not to be calculated -> USE METHOD 1 )       
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       ! (S)GNCB(i1,i2,i3) = 0 -> NO CELL IN BOX i1,i2,i3
       !                       -> SUM/AVE MUST BE ZERO
       !                       -> 0/0 -> 0/1
       gpl(:,:,:) = gpl(:,:,:) / &
            MERGE(ZSGNCB(:,:,:), 1.0_dp, ABS(ZSGNCB(:,:,:)) > EPS)
       !
    CASE(LG2GP_STD)
       ! STANDARD DEVIATION
       ! STANDARD DEVIATION OVER ALL CELLS -> 
       ! SIGMA = SQRT(| <x2> - <x>2 |) * SQRT(GNCB/(GNCB-1))
       !       = SQRT(CPU_SUM((x - <x>)**2)/(GNCB-1))
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       ALLOCATE(zgpl(nproma, nlev, ngpblks))
       zgpl(:,:,:) = gpl(:,:,:) / &
            MERGE(ZSGNCB(:,:,:), 1.0_dp, ABS(ZSGNCB(:,:,:)) > EPS)
       gpl(:,:,:) = (gpl(:,:,:)-zgpl(:,:,:))**2  &
            / MERGE(ZSGNCB(:,:,:)-1.0_dp, 1.0_dp, ABS(ZSGNCB(:,:,:)-1.0_dp) &
            > EPS)
       !
    CASE(LG2GP_AVEGT0)
       !
       ! AVERAGE OVER ALL ELEMETNS > 0
       ! ... SUM OVER ALL CELLS (> 0) AND DIVIDE BY NUMBER OF CELLS (> 0)
       gp(:,:,:) = gp(:,:,:) * zfngp(:,:,:)        ! CPU average -> CPU sum
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)      ! sum over all CPUs
       ALLOCATE(zgpl(nproma, nlev, ngpblks))          ! # of el. > 0 ...
       CALL trp_gpdc_gpgl(-1, zgpl, zfngp, M_SUM)  ! ... sum over all CPUs
       DO i1=1, nproma
          DO i2=1, nlev
             DO i3=1,ngpblks
                IF (zgpl(i1,i2,i3) > 0.0) THEN
                   gpl(i1,i2,i3) = gpl(i1,i2,i3) / zgpl(i1,i2,i3)
                ELSE
                   gpl(i1,i2,i3) = 0.0
                END IF
             END DO
          END DO
       END DO
       ! use "number of elements > 0" for optional filling below
       ZSGNCB => zgpl
       !
    CASE DEFAULT
       CALL finish(substr, 'invalid METHOD parameter')
    END SELECT

    ! OPTIONAL FILLING
    IF ( PRESENT(fill_value) .AND. PRESENT(fill_field) ) &
         CALL finish(substr, 'fill_value AND fill_field are ambiguous')
    IF (PRESENT(fill_value)) THEN
       gpl(:,:,:) = MERGE(gpl(:,:,:), fill_value, ABS(ZSGNCB(:,:,:)) > EPS)
    END IF
    IF (PRESENT(fill_field)) THEn
       gpl(:,:,:) = MERGE(gpl(:,:,:), fill_field(:,:,:), &
            ABS(ZSGNCB(:,:,:)) > EPS)
    END IF

    ! CLEAN UP
    IF (ASSOCIATED(gp))    DEALLOCATE(gp)
    IF (ASSOCIATED(zgpl))  DEALLOCATE(zgpl)
    IF (ASSOCIATED(zfngp)) DEALLOCATE(zfngp)

  END SUBROUTINE lg2gp_xyz_e5
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  SUBROUTINE lg2gp_xyzn_e5(lgl, gpl, method, lmcons, ltm1 &
       , fill_value, fill_field)

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, ngpblks, nlon, ngl  
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl, M_SUM
    USE messy_main_mpi_bi,       ONLY: finish, message
    USE messy_attila_mem_e5,     ONLY: SGNCB, SGNCBM1
    ! MESSy
    USE messy_main_timer,    ONLY: lstart !, lresume
    USE messy_attila_tools,  ONLY: lg2gp, EPS &
                                 , LG2GP_SUM, LG2GP_AVE, LG2GP_STD &
                                 , LG2GP_AVEGT0

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    ! lagrange vetor
    REAL(DP), DIMENSION(:,:),     INTENT(IN)              :: lgl
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(DP), DIMENSION(:,:,:,:), POINTER                 :: gpl
    ! method
    INTEGER,                  INTENT(IN)              :: method
    ! for mass conserving transformation
    LOGICAL,                  INTENT(IN),    OPTIONAL :: lmcons
    ! use positions of CELLS at t-1 ?
    LOGICAL,                  INTENT(IN),    OPTIONAL :: ltm1
    ! optional fill-value where no CELL is located
    REAL(DP),                     INTENT(IN),    OPTIONAL :: fill_value
    ! optional GP-field to fill in where no CELL is located
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN),    OPTIONAL :: fill_field

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lg2gp_xyzn_e5'
    INTEGER :: status
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: gp     => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zgpl   => NULL()
    REAL(DP), DIMENSION(:,:,:),   POINTER :: ZSGNCB => NULL()
    INTEGER :: s3, i
    ! number of 'events'
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zfngp  => NULL()
    INTEGER                           :: i1,i2,i3

    s3 = SIZE(lgl,2)

    IF (PRESENT(ltm1)) THEN
       IF ((ltm1).AND.lstart) THEN
          CALL message(substr,' WARNING: t-1 not defined at t=0;'&
               &//' transformation with ltm1=.true. not possible!')
          RETURN ! t-1 not defined at t=0
       END IF
       !
       IF (ltm1) THEN
          ZSGNCB => SGNCBM1
       ELSE
          ZSGNCB => SGNCB 
       END IF
    ELSE
       ZSGNCB => SGNCB 
    END IF

    ALLOCATE(gp(nlon, nlev, s3, ngl))

    SELECT CASE(method)
    CASE(LG2GP_AVEGT0)
       ALLOCATE(zfngp(nlon, nlev, s3, ngl))
       CALL lg2gp(lgl, gp, method, lmcons=lmcons   &
            , gorder='xzny', ltm1=ltm1, status=status, fngp=zfngp)
    CASE DEFAULT
       CALL lg2gp(lgl, gp, LG2GP_SUM, lmcons=lmcons   &
            , gorder='xzny', ltm1=ltm1, status=status)
    END SELECT
    IF (status /= 0) &
         CALL finish(substr,' Error in lg2gp!')

    ! NOTE:
    !
    ! THE RESULTS OF lg2gp ARE CONSISTENT FOR ALL CELLS ON 1 CPU !!!
    !
    ! (z)lmcons            .FALSE.   | .TRUE.
    ! => WEIGHT OF x     UNWEIGHTED  | MASS WEIGHTED (AMCELL/PMBOX)
    !         -----------------------------------------------------
    ! gp sum | SUM          SUM(x)
    ! gp ave | AVERAGE      <x> = SUM(x)/NCB
    ! gp std | STD.DEV      SQRT(| <x2> - <x>2 |) * SQRT(NCB/(NCB-1))
    !
    ! THESE 'PER CPU'-RESULTS NEED TO BE TRANSPOSED/TRANSFORMED 
    ! BACK INTO CONSISTENT FIELDS ON ALL CPUs
    !
    ! THE FOLLOWING CALCULATIONS ARE  INDEPENDENT OF MASS CONSERVATION
    ! ((z)lmcons), SINCE THE WEIGHT AMCELL/PMBOX IS CONSTANT FOR
    ! ALL CELLS IN ONE BOX


    SELECT CASE(METHOD)
    CASE(LG2GP_SUM)
       ! SUM
       ! SUMMATION OVER ALL CELLS -> SUMMATION OVER ALL CPUs (CPU_SUM)
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       !
    CASE(LG2GP_AVE)
       ! AVERAGE
       ! AVERAGE OVER ALL CELLS -> METHOD 1: CPU_SUM(gp_sum)/GNCB
       !                           METHOD 2: CPU_SUM(gp_ave*NCB)/GNCB
       !            (NCB not available in CHANNEL DECOMP.
       !             gp_ave need not to be calculated -> USE METHOD 1 )       
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       ! (S)GNCB(i1,i2,i3) = 0 -> NO CELL IN BOX i1,i2,i3
       !                       -> SUM/AVE MUST BE ZERO
       !                       -> 0/0 -> 0/1
       DO i=1, s3
          gpl(:,:,i,:) = gpl(:,:,i,:) / &
            MERGE(ZSGNCB(:,:,:), 1.0_dp, ABS(ZSGNCB(:,:,:)) > EPS)
       END DO
       !
    CASE(LG2GP_STD)
       ! STANDARD DEVIATION
       ! STANDARD DEVIATION OVER ALL CELLS -> 
       ! SIGMA = SQRT(| <x2> - <x>2 |) * SQRT(GNCB/(GNCB-1))
       !       = SQRT(CPU_SUM((x - <x>)**2)/(GNCB-1))
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       ALLOCATE(zgpl(nproma, nlev, s3, ngpblks))
       DO i=1, s3
          zgpl(:,:,i,:) = gpl(:,:,i,:) / &
            MERGE(ZSGNCB(:,:,:), 1.0_dp, ABS(ZSGNCB(:,:,:)) > EPS)
          gpl(:,:,i,:) = (gpl(:,:,i,:)-zgpl(:,:,i,:))**2  &
               / MERGE(ZSGNCB(:,:,:)-1.0_dp, 1.0_dp, ABS(ZSGNCB(:,:,:)-1.0_dp) &
               > EPS)
       END DO
       !
    CASE(LG2GP_AVEGT0)
       !
       ! AVERAGE OVER ALL ELEMETNS > 0
       ! ... SUM OVER ALL CELLS (> 0) AND DIVIDE BY NUMBER OF CELLS (> 0)
       gp(:,:,:,:) = gp(:,:,:,:) * zfngp(:,:,:,:)  ! CPU average -> CPU sum
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)      ! sum over all CPUs
       ALLOCATE(zgpl(nproma, nlev, s3, ngpblks))      ! # of el. > 0 ...
       CALL trp_gpdc_gpgl(-1, zgpl, zfngp, M_SUM)  ! ... sum over all CPUs
       DO i1=1, nproma
          DO i2=1, nlev
             DO i=1, s3
                DO i3=1,ngpblks
                   IF (zgpl(i1,i2,i,i3) > 0.0) THEN
                      gpl(i1,i2,i,i3) = gpl(i1,i2,i,i3) / zgpl(i1,i2,i,i3)
                   ELSE
                      gpl(i1,i2,i,i3) = 0.0
                   END IF
                END DO
             END DO
          END DO
       END DO
       !
    CASE DEFAULT
       CALL finish(substr, 'invalid METHOD parameter')
    END SELECT

    ! OPTIONAL FILLING
    IF ( PRESENT(fill_value) .AND. PRESENT(fill_field) ) &
         CALL finish(substr, 'fill_value AND fill_field are ambiguous')
    IF (PRESENT(fill_value)) THEN
       SELECT CASE(METHOD)
       CASE(LG2GP_AVEGT0)
          DO i=1, s3
             gpl(:,:,i,:) = MERGE(gpl(:,:,i,:), fill_value &
                  , ABS(zgpl(:,:,i,:)) > EPS)
          END DO
       CASE DEFAULT
          DO i=1, s3
             gpl(:,:,i,:) = MERGE(gpl(:,:,i,:), fill_value &
                  , ABS(ZSGNCB(:,:,:)) > EPS)
          END DO
       END SELECT
    END IF
    IF (PRESENT(fill_field)) THEN
       SELECT CASE(METHOD)
       CASE(LG2GP_AVEGT0)
          DO i=1, s3
             gpl(:,:,i,:) = MERGE(gpl(:,:,i,:), fill_field(:,:,i,:), &
                  ABS(zgpl(:,:,i,:)) > EPS)
          END DO
       CASE DEFAULT
          DO i=1, s3
             gpl(:,:,i,:) = MERGE(gpl(:,:,i,:), fill_field(:,:,i,:), &
                  ABS(ZSGNCB(:,:,:)) > EPS)
          END DO
       END SELECT
    END IF

    ! CLEAN UP
    IF (ASSOCIATED(gp))   DEALLOCATE(gp)
    IF (ASSOCIATED(zgpl)) DEALLOCATE(zgpl)
    IF (ASSOCIATED(zfngp)) DEALLOCATE(zfngp)

  END SUBROUTINE lg2gp_xyzn_e5
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  SUBROUTINE gpsfemis2lgemis_e5(gpl, lgl, method, gprl, lmcons)
    
    USE messy_attila_tools,      ONLY: gpsfemis2lgemis
    USE messy_main_transform_bi, ONLY: trp_gpdc_gpgl, M_LOC
    USE messy_main_mpi_bi,       ONLY: finish

    IMPLICIT NONE

    ! I/O
    ! local GRIDPOINT FIELD (channel)
    REAL(DP), DIMENSION(:,:), POINTER              :: gpl  ! INTENT(IN)
    ! local LAGRANGIAN FIELD (channel)
    REAL(DP), DIMENSION(:),   POINTER              :: lgl  ! INTENT(OUT)
    ! METHOD
    INTEGER,                        INTENT(IN) :: method
    ! local GRIDPOINT FIELD (channel) WITH REST OF gpr
    REAL(DP), DIMENSION(:,:), POINTER,    OPTIONAL :: gprl ! INTENT(INOUT)
    ! SWITCH FOR MASS CONSERVATION
    LOGICAL,              INTENT(IN), OPTIONAL :: lmcons

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER    :: substr = 'gpsfemis2lgemis_e5'
    REAL(DP), DIMENSION(:,:), POINTER  :: gpg  => NULL() ! globalized field
    REAL(DP), DIMENSION(:,:), POINTER  :: gprg => NULL() ! globalized rest
    INTEGER                        :: status

    ! GLOBAL FIELD ON ALL CPUs
    CALL trp_gpdc_gpgl(1, gpl, gpg)
    ! GLOBAL REST ON ALL CPUs
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(1, gprl, gprg)

    ! TRANSFER FIELD INTO LAGRANGIAN SPACE (MASS CONSERVING) AND SAVE REST
    IF (PRESENT(gprl)) THEN
       CALL gpsfemis2lgemis(gpg, lgl, method, gpr=gprg   &
            , lmcons=lmcons, gorder='xy', status=status)
    ELSE
       CALL gpsfemis2lgemis(gpg, lgl, method             &
            , lmcons=lmcons, gorder='xy', status=status)
    END IF
    IF (status /= 0) THEN
       IF (status == 2) &
            CALL finish(substr, 'Unknown method !')
       CALL finish(substr,' Error in gpsfemis2lgemis!')
      END IF

    ! GET REST FROM ALL CPUs; PUT BACK INTO CHANNEL
    ! NOTE: Since all CPUs calculate (gpsfemis2lgemis) 
    ! the same REST with GNCB, the decomposed array can
    ! simply be 'picked' out locally (M_LOC) ...
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(-1, gprl, gprg, M_LOC)

    ! CLEAN UP
    IF (ASSOCIATED(gpg))  DEALLOCATE(gpg)
    IF (ASSOCIATED(gprg)) DEALLOCATE(gprg)

  END SUBROUTINE gpsfemis2lgemis_e5
!----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE lggpte2lgte_xyz_e5(gp, gpte, lg, lgte)

    USE messy_attila_mem_e5,     ONLY: SGNCB
    USE messy_attila,            ONLY: NCELL

    IMPLICIT NONE

    ! NOTE:
    ! THIS IS FOR CALCULATION OF PSEUDO-LG TENDENCIES BY PROCESSES IN GP
    ! REPRESENTATION:
    !
    ! positive definite !
    !
    ! SEQUENCE:
    !    CALL lg2gp_e5(lg, gp, LG2GP_AVE, fill_value=0._dp)
    !    CALL ... process(gp -> gpte) ...
    !    CALL lggpte2lgte_e5(gp, gpte, lg, lgte)
    !    lgte is the resulting LG process tendency

    ! I/O
    ! GP field transformed (AVE) from lg (below)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: gp
    ! GP tendency (arbitrary process)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: gpte
    ! original LG field corresponding to gp (above)
    REAL(DP), DIMENSION(:),     INTENT(IN)  :: lg
    ! resulting LG tendency
    REAL(DP), DIMENSION(:),     POINTER     :: lgte       ! INTENT(OUT)

    ! LOCAL
    REAL(DP), DIMENSION(:,:,:), POINTER  :: gpsc => NULL()
    INTEGER :: n1, n2, n3
    INTEGER :: i1, i2, i3
    INTEGER :: jc
    
    ! INIT
    n1 = SIZE(gp,1)
    n2 = SIZE(gp,2)
    n3 = SIZE(gp,3)

    ! MEMORY
    ALLOCATE(gpsc(n1,n2,n3))
    IF (ASSOCIATED(lgte)) THEN
       DEALLOCATE(lgte)
       NULLIFY(lgte)
    END IF
    ALLOCATE(lgte(NCELL))

    DO i1=1,n1
       DO i2=1,n2
          DO i3=1,n3
             IF ( gpte(i1,i2,i3) > TINY(0.) ) THEN
                ! tendency > 0._dp   === source
                IF (NINT(SGNCB(i1,i2,i3)) >= 1) THEN
                   gpsc(i1,i2,i3) = &
                        gpte(i1,i2,i3) / SGNCB(i1,i2,i3)
                ELSE
                   gpsc(i1,i2,i3) = 0.
                ENDIF
             ELSE
                ! tendency < 0._dp   === sink
                IF (gp(i1,i2,i3) > TINY(0.) ) THEN
                   gpsc(i1,i2,i3) =    &
                        gpte(i1,i2,i3) /  &
                        gp(i1,i2,i3) 
                ELSE
                   gpsc(i1,i2,i3) = 0.
                ENDIF
             END IF
          ENDDO
       ENDDO
    ENDDO

    CALL gp2lg_e5(gpsc, lgte)

    DO jc=1, NCELL
       IF (lgte(jc) < TINY(0.)) lgte(jc) = lgte(jc) * lg(jc) 
    ENDDO

    ! CLEAN UP
    DEALLOCATE(gpsc)

  END SUBROUTINE lggpte2lgte_xyz_e5
!----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  SUBROUTINE lggpte2lgte_xyzn_e5(gp, gpte, lg, lgte)

    USE messy_attila_mem_e5,     ONLY: SGNCB
    USE messy_attila,            ONLY: NCELL

    IMPLICIT NONE

    ! I/O
    ! GP field transformed (AVE) from lg (below)
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN)  :: gp
    ! GP tendency (arbitrary process)
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN)  :: gpte
    ! original LG field corresponding to gp (above)
    REAL(DP), DIMENSION(:,:),     INTENT(IN)  :: lg
    ! resulting LG tendency
    REAL(DP), DIMENSION(:,:),     POINTER     :: lgte  ! INTENT(OUT)

    ! LOCAL
    REAL(DP), DIMENSION(:,:,:,:), POINTER  :: gpsc => NULL()
    INTEGER :: n1, n2, n3, n4
    INTEGER :: i1, i2, i3, i4
    INTEGER :: jc
    
    ! INIT
    n1 = SIZE(gp,1)
    n2 = SIZE(gp,2)
    n3 = SIZE(gp,3)
    n4 = SIZE(gp,4)

    ! MEMORY
    ALLOCATE(gpsc(n1,n2,n3,n4))
    IF (ASSOCIATED(lgte)) THEN
       DEALLOCATE(lgte)
       NULLIFY(lgte)
    END IF
    ALLOCATE(lgte(NCELL,n3))

    DO i1=1,n1
       DO i2=1,n2
          DO i3=1,n3
             DO i4=1,n4
                IF ( gpte(i1,i2,i3,i4) > TINY(0.) ) THEN
                   ! tendency > 0._dp   === source
                   IF (NINT(SGNCB(i1,i2,i4)) >= 1) THEN
                      gpsc(i1,i2,i3,i4) = &
                           gpte(i1,i2,i3,i4) / SGNCB(i1,i2,i4)
                   ELSE
                      gpsc(i1,i2,i3,i4) = 0.
                   ENDIF
                ELSE
                   ! tendency < 0._dp   === sink
                   IF (gp(i1,i2,i3,i4) > TINY(0.) ) THEN
                      gpsc(i1,i2,i3,i4) =    &
                           gpte(i1,i2,i3,i4) /  &
                           gp(i1,i2,i3,i4) 
                   ELSE
                      gpsc(i1,i2,i3,i4) = 0.
                   ENDIF
                END IF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    CALL gp2lg_e5(gpsc, lgte)

    DO jc=1, NCELL
       DO i3 = 1, n3
          IF (lgte(jc,i3) < TINY(0.)) lgte(jc,i3) = lgte(jc,i3) * lg(jc,i3) 
       ENDDO
    END DO

    ! CLEAN UP
    DEALLOCATE(gpsc)

  END SUBROUTINE lggpte2lgte_xyzn_e5
!----------------------------------------------------------------------------

! ***************************************************************************
END MODULE messy_attila_tools_e5
! ***************************************************************************
