! ***************************************************************************
MODULE messy_clams_tools_e5
! ***************************************************************************

  ! MODULE FOR TRANSFORMATIONS/TRANSPOSITIONS BETWEEN GRIDPOINT SPACE
  ! AND LAGRANGIAN SPACE, BOTH IN ECHAM5 DECOMPOSITION
  !
  ! Author: Patrick Joeckel, MPICH, Oct 2003
  !         Sabine Brinkop, DLR, Jul 2019 (module applied for CLaMS)

  USE messy_clams_tools,           ONLY: CL2GP_SUM,  CL2GP_AVE,  CL2GP_STD &
                                       , CL2GP_AVEGT0
  USE messy_clams_global,          ONLY: dnparts_max
!  USE messy_clams,             ONLY: pc_g

  USE messy_main_channel_error_bi, ONLY: channel_halt
  USE messy_main_channel,          ONLY: get_channel_object
  USE messy_main_constants_mem,    ONLY: DP


  IMPLICIT NONE
  PRIVATE

!  PUBLIC :: pc_g
  PUBLIC :: gp2cl_e5
  PUBLIC :: cl2gp_e5
  PUBLIC :: CL2GP_SUM,  CL2GP_AVE,  CL2GP_STD, CL2GP_AVEGT0
  PUBLIC :: gpsfemis2clemis_e5
  PUBLIC :: dp

  INTERFACE gp2cl_e5
     MODULE PROCEDURE gp2cl_xy_e5
     MODULE PROCEDURE gp2cl_xyz_e5
     MODULE PROCEDURE gp2cl_xyzn_e5
  END INTERFACE

  INTERFACE cl2gp_e5
     MODULE PROCEDURE cl2gp_xyz_e5
     MODULE PROCEDURE cl2gp_xyzn_e5
  END INTERFACE

  INTRINSIC :: PRESENT, ABS, MERGE, ASSOCIATED

  ! number of clams parcel per echam grid box global field in decomposition
  REAL(dp), DIMENSION(:,:,:), POINTER :: spc_g => NULL()

CONTAINS
!----------------------------------------------------------------------------
  SUBROUTINE gp2cl_xy_e5(gpl, lgl, gprl, lmcons, klev, llev)

    USE messy_clams_tools,       ONLY: gp2cl
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
    CHARACTER(LEN=*), PARAMETER    :: substr = 'gp2cl_xy_e5'
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
       CALL gp2cl(gpg, lgl, gpr=gprg &
            , lmcons=lmcons, klev=klev, gorder='xy', llev=llev, status=status)
    ELSE
       CALL gp2cl(gpg, lgl           &
            , lmcons=lmcons, klev=klev, gorder='xy', llev=llev, status=status)
    END IF
    IF (status /= 0) &
         CALL finish(substr,' Error in gp2cl!')

    ! GET REST FROM ALL CPUs; PUT BACK INTO CHANNEL
    ! NOTE: 
    !   - AVE, SINCE ALL CPUs calculate (gp2cl) the same REST with GNCB !!!
    !   - LOC even better
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(-1, gprl, gprg, M_LOC)

    ! CLEAN UP
    IF (ASSOCIATED(gpg))  DEALLOCATE(gpg)
    IF (ASSOCIATED(gprg)) DEALLOCATE(gprg)

  END SUBROUTINE gp2cl_xy_e5
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE gp2cl_xyz_e5(gpl, lgl, gprl, lmcons)

    USE messy_clams_tools,      ONLY: gp2cl
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
    CHARACTER(LEN=*), PARAMETER      :: substr = 'gp2cl_xyz_e5'
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
       CALL gp2cl(gpg, lgl, gpr=gprg &
            , lmcons=lmcons, gorder='xzy', status=status)
    ELSE
       CALL gp2cl(gpg, lgl           &
            , lmcons=lmcons, gorder='xzy', status=status)
    END IF
    IF (status /= 0) &
         CALL finish(substr,' Error in gp2cl!')

    ! GET REST FROM ALL CPUs; PUT BACK INTO CHANNEL
    ! NOTE:
    !   - AVE, SINCE ALL CPUs calculate (gp2cl) the same REST with GNCB !!!
    !   - LOC, even better
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(-1, gprl, gprg, M_LOC)

    ! CLEAN UP
    IF (ASSOCIATED(gpg))  DEALLOCATE(gpg)
    IF (ASSOCIATED(gprg)) DEALLOCATE(gprg)

  END SUBROUTINE gp2cl_xyz_e5
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
  SUBROUTINE gp2cl_xyzn_e5(gpl, lgl, gprl, lmcons)

    USE messy_clams_tools,       ONLY: gp2cl
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
    CHARACTER(LEN=*), PARAMETER        :: substr = 'gp2cl_xyzn_e5'
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
       CALL gp2cl(gpg, lgl, gpr=gprg &
            , lmcons=lmcons, gorder='xzny', status=status)
    ELSE
       CALL gp2cl(gpg, lgl           &
            , lmcons=lmcons, gorder='xzny', status=status)
    END IF
    IF (status /= 0) &
         CALL finish(substr,' Error in gp2cl!')

    ! GET REST FROM ALL CPUs; PUT BACK INTO CHANNEL
    ! NOTE:
    !   - AVE, SINCE ALL CPUs calculate (gp2cl) the same REST with GNCB !!!
    !   - LOC, even better
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(-1, gprl, gprg, M_LOC)

    ! CLEAN UP
    IF (ASSOCIATED(gpg))  DEALLOCATE(gpg)
    IF (ASSOCIATED(gprg)) DEALLOCATE(gprg)

  END SUBROUTINE gp2cl_xyzn_e5
!----------------------------------------------------------------------------
!--------------------------------------------------------------------------
  SUBROUTINE cl2gp_xyz_e5(lgl, gpl, method,  &
        fill_value, fill_field)

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, ngpblks, nlon, ngl 
    USE messy_main_transform_bi,    ONLY: trp_gpdc_gpgl, M_SUM
    USE messy_main_mpi_bi,          ONLY: finish, message
    ! MESSy
    USE messy_main_timer,   ONLY: lstart !, lresume
    USE messy_clams_tools,  ONLY:  cl2gp, EPS &
                                 , CL2GP_SUM, CL2GP_AVE, CL2GP_STD &
                                 , CL2GP_AVEGT0

    IMPLICIT NONE

    ! I/O
    ! lagrange vetor
    REAL(DP), DIMENSION(:),     INTENT(IN)              :: lgl
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(DP), DIMENSION(:,:,:), POINTER                 :: gpl
    ! method
    INTEGER,                INTENT(IN)              :: method
    ! optional fill-value where no CELL is located
    REAL(DP),                   INTENT(IN),    OPTIONAL :: fill_value
    ! optional GP-field to fill in where no CELL is located
    REAL(DP), DIMENSION(:,:,:), INTENT(IN),    OPTIONAL :: fill_field

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cl2gp_xyz_e5'
    INTEGER :: status
    REAL(DP), DIMENSION(:,:,:), POINTER :: gp     => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: zgpl   => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: ZSGNCB => NULL()
    ! number of 'events'
    REAL(DP), DIMENSION(:,:,:), POINTER :: zfngp  => NULL()
    INTEGER                         :: i1,i2,i3

    CALL get_channel_object(status, 'clams', 'SPC_G', p3=spc_g)
    CALL channel_halt(substr, status)
    
    ! parcel number per ECHAM grid box
    ZSGNCB => spc_g

    ALLOCATE(gp(nlon, nlev, ngl))

    SELECT CASE(method)
    CASE(CL2GP_AVEGT0)
       ALLOCATE(zfngp(nlon, nlev, ngl))
       CALL cl2gp(lgl, gp, method,  &
            gorder='xzy', status=status, fngp=zfngp)
    CASE DEFAULT
       CALL cl2gp(lgl, gp, CL2GP_SUM,   &
            gorder='xzy', status=status)
    END SELECT
    IF (status /= 0) &
         CALL finish(substr,' Error in cl2gp!')

    ! NOTE:
    !
    ! THE RESULTS OF cl2gp ARE CONSISTENT FOR ALL CELLS ON 1 CPU !!!
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
    CASE(CL2GP_SUM)
       ! SUM
       ! SUMMATION OVER ALL CELLS -> SUMMATION OVER ALL CPUs (CPU_SUM)
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       !
    CASE(CL2GP_AVE)
       ! AVERAGE
       ! AVERAGE OVER ALL CELLS -> METHOD 1: CPU_SUM(gp_sum)/GNCB
       !                           METHOD 2: CPU_SUM(gp_ave*NCB)/GNCB
       !            (NCB not available in CHANNEL DECOMP.
       !             gp_ave need not to be calculated -> USE METHOD 1 )       
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       ! (S)PC_G(i1,i2,i3) = 0 -> NO CELL IN BOX i1,i2,i3
       !                       -> SUM/AVE MUST BE ZERO
       !                       -> 0/0 -> 0/1
       gpl(:,:,:) = gpl(:,:,:) / &
            MERGE(ZSGNCB(:,:,:), 1.0_dp, ABS(ZSGNCB(:,:,:)) > EPS)
       !
    CASE(CL2GP_STD)
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
    CASE(CL2GP_AVEGT0)
       !
       ! AVERAGE OVER ALL ELEMENTS > 0
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

  END SUBROUTINE cl2gp_xyz_e5
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
  SUBROUTINE cl2gp_xyzn_e5(lgl, gpl, method,  &
       fill_value, fill_field)

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, ngpblks, nlon, ngl  
    USE messy_main_transform_bi,    ONLY: trp_gpdc_gpgl, M_SUM
    USE messy_main_mpi_bi,          ONLY: finish, message
    ! MESSy
    USE messy_main_timer,        ONLY: lstart !, lresume
    USE messy_clams_tools,       ONLY: cl2gp, EPS &
                                 , CL2GP_SUM, CL2GP_AVE, CL2GP_STD &
                                 , CL2GP_AVEGT0

    IMPLICIT NONE

    INTRINSIC :: SIZE

    ! I/O
    ! lagrange vetor
    REAL(DP), DIMENSION(:,:),     INTENT(IN)              :: lgl
    ! global grid point field ([mol/mol] or [kg/kg])
    REAL(DP), DIMENSION(:,:,:,:), POINTER                 :: gpl
    ! method
    INTEGER,                  INTENT(IN)              :: method
    ! optional fill-value where no CELL is located
    REAL(DP),                     INTENT(IN),    OPTIONAL :: fill_value
    ! optional GP-field to fill in where no CELL is located
    REAL(DP), DIMENSION(:,:,:,:), INTENT(IN),    OPTIONAL :: fill_field

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'cl2gp_xyzn_e5'
    INTEGER :: status
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: gp     => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zgpl   => NULL()
    REAL(DP), DIMENSION(:,:,:),   POINTER :: ZSGNCB => NULL()
    INTEGER :: s3, i
    ! number of 'events'
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zfngp  => NULL()
    INTEGER                           :: i1,i2,i3

    s3 = SIZE(lgl,2)

    CALL get_channel_object(status, 'clams', 'SPC_G', p3=spc_g)
    CALL channel_halt(substr, status)

    ZSGNCB => spc_g

    ALLOCATE(gp(nlon, nlev, s3, ngl))

    SELECT CASE(method)
    CASE(CL2GP_AVEGT0)
       ALLOCATE(zfngp(nlon, nlev, s3, ngl))
       CALL cl2gp(lgl, gp, method,  &
            gorder='xzny', status=status, fngp=zfngp)
    CASE DEFAULT
       CALL cl2gp(lgl, gp, CL2GP_SUM,   &
            gorder='xzny', status=status)
    END SELECT
    IF (status /= 0) &
         CALL finish(substr,' Error in cl2gp!')

    ! NOTE:
    !
    ! THE RESULTS OF cl2gp ARE CONSISTENT FOR ALL CELLS ON 1 CPU !!!
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
    CASE(CL2GP_SUM)
       ! SUM
       ! SUMMATION OVER ALL CELLS -> SUMMATION OVER ALL CPUs (CPU_SUM)
       CALL trp_gpdc_gpgl(-1, gpl, gp, M_SUM)
       !
    CASE(CL2GP_AVE)
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
    CASE(CL2GP_STD)
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
    CASE(CL2GP_AVEGT0)
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
       CASE(CL2GP_AVEGT0)
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
       CASE(CL2GP_AVEGT0)
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

  END SUBROUTINE cl2gp_xyzn_e5
!----------------------------------------------------------------------------

 SUBROUTINE gpsfemis2clemis_e5(gpl, lgl, method, gprl, lmcons)
    
    USE messy_clams_tools,       ONLY: gpsfemis2clemis
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
    CHARACTER(LEN=*), PARAMETER    :: substr = 'gpsfemis2clemis_e5'
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
       CALL gpsfemis2clemis(gpg, lgl, method, gpr=gprg   &
            , lmcons=lmcons, gorder='xy', status=status)
    ELSE
       CALL gpsfemis2clemis(gpg, lgl, method             &
            , lmcons=lmcons, gorder='xy', status=status)  
    END IF
    IF (status /= 0) THEN
       IF (status == 2) &
            CALL finish(substr, 'Unknown method !')
       CALL finish(substr,' Error in gpsfemis2clemis!')
      END IF

    ! GET REST FROM ALL CPUs; PUT BACK INTO CHANNEL
    ! NOTE: Since all CPUs calculate (gpsfemis2clemis) 
    ! the same REST with , the decomposed array can
    ! simply be 'picked' out locally (M_LOC) ...
    IF (PRESENT(gprl)) &
         CALL trp_gpdc_gpgl(-1, gprl, gprg, M_LOC)

    ! CLEAN UP
    IF (ASSOCIATED(gpg))  DEALLOCATE(gpg)
    IF (ASSOCIATED(gprg)) DEALLOCATE(gprg)

  END SUBROUTINE gpsfemis2clemis_e5
!----------------------------------------------------------------------------

! ***************************************************************************
END MODULE messy_clams_tools_e5
! ***************************************************************************
