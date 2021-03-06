! ******************************************************************
MODULE messy_main_grid_trafo_nrgd
! ******************************************************************
  ! ------------------------------------------------------------------
  ! Author: Patrick Joeckel, MPICH, Mainz, June 2002
  !         Astrid  Kerkweg, UniMz, Mainz, 2012-2013
  !         -> re-structured/expanded for more general GRID application
  ! ------------------------------------------------------------------
  
  USE messy_main_constants_mem,        ONLY: SP, DP, I4, I8
  USE messy_main_grid_trafo_nrgd_base, ONLY: t_axis, PRINT_AXIS, INIT_AXIS &
                                           , NREGRID, NREGRID_STAT
  USE messy_main_grid_trafo,    ONLY: RGEMPTY, PACK_GEOHYBGRID_NCVAR &
                                    , SORT_GEOHYBGRID_NCVAR          &
                                    , SORT_GEOHYBGRID                &
                                    , BALANCE_GEOHYBGRID_NCVAR       &
                                    , BALANCE_GEOHYBGRID_TIME        &
                                    , CHECK_NCVAR_ON_GEOHYBGRID      &
                                    , SWITCH_GEOHYBGRID              &
                                    , COMPLETE_GEOHYBGRID            &
                                    , BALANCE_GEOHYBGRID             &
                                    , CHECK_GEOHYBGRID, H2PSIG       &
                                    , RGTstr, GTRF_NRGD              &
                                    , RG_IDX, RG_IXF, RG_INT, RG_EXT
  USE messy_main_grid_netcdf,   ONLY: t_ncatt, t_ncvar, t_narray, GRD_MAXSTRLEN&
                                    , NULL_DIMID, NULL_VARID                   &
                                    , VTYPE_DOUBLE, VTYPE_INT,  VTYPE_REAL     &
                                    , VTYPE_BYTE,   VTYPE_CHAR, VTYPE_UNDEF    &
                                    , ERRMSG, RGMLE, RGMLEC, RGMLI,  RGMLIC    &
                                    , RGMLW, RGMLWC, RGMLVM, RGMLVL            &
                                    , RGMLVLC, RGMLVMC, RGMSG, MSGMODE         &
                                    , MSGMODE_S, MSGMODE_E, MSGMODE_VL         &
                                    , MSGMODE_W, MSGMODE_VM, MSGMODE_I         &
                                    , nf90_inq_libvers, nf90_global, nf90_char &
                                    , INIT_NARRAY, SCALE_NARRAY, CAT_NARRAY    &
                                    , COPY_NARRAY, COPY_NCDIM                  &
                                    , IMPORT_NCVAR, EXPORT_NCVAR, SCAN_NCVAR   &
                                    , RENAME_NCVAR, INIT_NCVAR, COPY_NCVAR     &
                                    , MAXFRAC2IDX_NCVAR, IDX2FRAC_NCVAR        &
                                    , IMPORT_NCATT, EXPORT_NCATT, ADD_NCATT    &
                                    , INIT_NCATT, string                       &
                                    , QDEF_NCVAR, QCMP_NCDIM                   &
                                    , EXPORT_NCATT, COPY_NCATT, PRINT_NCVAR
  USE messy_main_grid,          ONLY: t_geohybgrid                         &
                                    , COPY_GEOHYBGRID, INIT_GEOHYBGRID     &
                                    , IMPORT_GEOHYBGRID, EXPORT_GEOHYBGRID &
                                    , PRINT_GEOHYBGRID
  USE messy_main_tools,         ONLY: INT2STR

  IMPLICIT NONE

  INTRINSIC :: TRIM, SIZE, PRESENT, ASSOCIATED, REAL, COS, IAND

  PRIVATE   :: TRIM, SIZE, PRESENT, ASSOCIATED, REAL, COS, IAND

  REAL (DP), PARAMETER :: PI = 3.141592653589_DP 

  CHARACTER(*), PARAMETER :: NCREGRIDVERS = '1.5b'

  INTERFACE REGRID_CONTROL
     MODULE PROCEDURE REGRID_CONTROL
  END INTERFACE

  PUBLIC :: REGRID_CONTROL

CONTAINS
  
  ! ------------------------------------------------------------------
  SUBROUTINE REGRID_CONTROL( grid_in, grid_out, tvar, var        &
       , RG_TYPE, lint                &
       , lrgx, lrgy, lrgz             &
       , lfirsto                      &
       , lpresaxis                    & 
       , lwork                        &
       , lstatout                     &
       , grid_conv)

    IMPLICIT NONE

    TYPE (t_geohybgrid), INTENT(IN)          :: grid_in     ! input  grid info
    TYPE (t_geohybgrid), INTENT(IN)          :: grid_out    ! output grid info
    TYPE (t_ncvar), DIMENSION(:), POINTER    :: tvar  ! list of input  variables
    TYPE (t_ncvar), DIMENSION(:), POINTER    :: var   ! list of output variables
    INTEGER,  INTENT(IN)                     :: RG_TYPE(:)  ! regrid type
    LOGICAL,  INTENT(IN)                     :: lint        ! input time ?
    LOGICAL          , INTENT(IN) , OPTIONAL :: lrgx     ! regrid along 'lon'
    LOGICAL          , INTENT(IN) , OPTIONAL :: lrgy     ! regrid along 'lat'
    LOGICAL          , INTENT(IN) , OPTIONAL :: lrgz     ! regrid along 'lev'
    LOGICAL          , INTENT(IN) , OPTIONAL :: lfirsto  ! first output step
    LOGICAL          , INTENT(IN),  OPTIONAL :: lpresaxis ! pressure axis regrid
    LOGICAL          , INTENT(IN),  OPTIONAL :: lwork     ! i_am_worker = T ?
    ! disable ncregrid statistics independent from MSG_MODE;
    ! in ncregrid statistic output is always (!) switched on, if any messages
    ! shall be written; for 2-way coupling we need the error messages, but
    ! the statistics makes the output file big and unreadable)
    LOGICAL            , INTENT(IN),  OPTIONAL :: lstatout
    TYPE (t_geohybgrid), INTENT(OUT), OPTIONAL :: grid_conv

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'REGRID_CONTROL'
    LOGICAL :: lx    ! regrid along x
    LOGICAL :: ly    ! regrid along y
    LOGICAL :: lz    ! regrid along z
    INTEGER :: nvars ! number of variables
    INTEGER :: mvars ! 'effective' number of variables
    INTEGER :: status
    TYPE (t_ncvar), DIMENSION(:), POINTER  :: xivar ! list of variables (input)
    TYPE (t_ncvar), DIMENSION(:), POINTER  :: xovar ! list of variables (output)
    TYPE (t_ncvar)                         :: qvari ! temporal variable
    TYPE (t_ncvar)                         :: svari ! temporal variable (sorted)
    TYPE (t_ncvar)                         :: pvari ! temporal variable (packed)
    TYPE (t_ncvar)                         :: qvaro ! temporal variable
    TYPE (t_ncvar), DIMENSION(:), POINTER  :: svaro ! temporal variable (sorted)
    TYPE (t_ncvar), DIMENSION(:), POINTER  :: pvaro ! temporal variable (packed)
    ! NAMELIST
    TYPE (t_geohybgrid)                   :: gi  ! input grid
    TYPE (t_geohybgrid)                   :: go  ! output grid
    TYPE (t_geohybgrid)                   :: gg  ! grdfile grid
    LOGICAL, SAVE                         :: lp     ! pressure or sigma
    ! GEOHYBRID-GRIDS
    TYPE (t_geohybgrid)                   :: gis, gix ! sorted and index input-grid
    TYPE (t_ncatt), SAVE                  :: rggatt   ! RG global attribute
    LOGICAL                               :: ok       ! result OK ?
    TYPE (t_geohybgrid)                   :: ggs, ggx ! sorted and index grid
    ! AXES FOR REGRIDDING
    TYPE (t_axis),  DIMENSION(:), POINTER :: sax, dax ! source and dest. axes
    ! N-ARRAYS FOR REGRIDDING
    TYPE (t_narray),DIMENSION(:), POINTER :: nai, nao ! input/output for regridder
    ! REGRIDDING STATISTICS
    REAL (DP),      DIMENSION(:), POINTER :: sovl, dovl
    INTEGER,        DIMENSION(:), POINTER :: rcnt
    ! VARIABLES
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ivar(:)   ! input variable names
    CHARACTER(LEN=GRD_MAXSTRLEN), POINTER :: ovar(:)   ! output variable names
    REAL                        , POINTER :: scl(:)    ! scaling
    INTEGER                     , POINTER :: RGT(:)    ! regridding type
    INTEGER, DIMENSION(:)       , POINTER :: dims      ! order of grid dimensions
    INTEGER                               :: axes(3)   ! grid axes of variable
    LOGICAL                               :: llwork
    LOGICAL                               :: llstatout
    INTEGER                               :: ixf(2)

    IF (PRESENT(lwork)) THEN
       llwork = lwork
    ELSE
       llwork = .TRUE.
    END IF

    IF (PRESENT(lstatout)) THEN
       llstatout = lstatout
    ELSE
       llstatout = .TRUE.
    END IF

    ! all worker PEs compute in parallel
    CALL REGRID_CONTROL_INIT

    IF (llwork)   THEN
       CALL REGRID_CONTROL_WORK
    END IF

    RETURN

  CONTAINS

    SUBROUTINE REGRID_CONTROL_INIT

      IMPLICIT NONE

      INTEGER :: i     ! counter

      ! INIT
      IF (PRESENT(lrgx)) THEN
         lx = lrgx
      ELSE
         lx = .TRUE.               ! DEFAULT
      END IF
      !
      IF (PRESENT(lrgy)) THEN
         ly = lrgy
      ELSE
         ly = .TRUE.               ! DEFAULT
      END IF
      !
      IF (PRESENT(lrgz)) THEN
         lz = lrgz
      ELSE
         lz = .TRUE.               ! DEFAULT
      END IF
      !
      IF (PRESENT(lpresaxis)) THEN
         lp = lpresaxis
      ELSE
         lp = .FALSE.
      ENDIF

      ! INITIALIZE OUTPUT PARAMETER
      IF (ASSOCIATED(var)) THEN
         DO i=1, SIZE(var)
            CALL INIT_NCVAR(var(i))
         END DO
         DEALLOCATE(var, STAT=status)
         CALL ERRMSG(substr,status,1)
      END IF
      NULLIFY(var)

      ! INIT
      NULLIFY(dims)
      NULLIFY(ivar)
      NULLIFY(ovar)
      NULLIFY(scl)
      NULLIFY(RGT)
      NULLIFY(sax)
      NULLIFY(dax)
      NULLIFY(xivar)
      NULLIFY(xovar)
      NULLIFY(svaro)
      NULLIFY(pvaro)
      NULLIFY(nai)
      NULLIFY(nao)
      NULLIFY(sovl)
      NULLIFY(dovl)
      NULLIFY(rcnt)

      CALL RGMSG(substr, RGMLVL, &
           '-------------------------------------------------------')
      CALL RGMSG(substr, RGMLVL, &
           'START REGRIDDING PROCEDURE (NCREGRID VERSION '//&
           &TRIM(NCREGRIDVERS)//')' )
      CALL RGMSG(substr, RGMLVL, &
           ' <Author: Patrick Joeckel, MPICH, June 2002>' )
      CALL RGMSG(substr, RGMLVL, &
           ' ( COMPILED WITH netCDF Fortran90 LIBRARY' )
      CALL RGMSG(substr, RGMLVL, &
           '   VERSION '//TRIM(NF90_inq_libvers())//' )' )
      CALL RGMSG(substr, RGMLVL, &
           '-------------------------------------------------------')

      RETURN

    END SUBROUTINE REGRID_CONTROL_INIT

    SUBROUTINE REGRID_CONTROL_WORK

      IMPLICIT NONE

      INTEGER :: i, ix

      nvars = SIZE(tvar)
      ALLOCATE(xivar(nvars), STAT=status)
      ALLOCATE(RGT(nvars)) 

      CALL ERRMSG(substr,status,3)

      CALL RGMSG(substr, RGMLVL, &
           '-------------------------------------------------------')
      ! copy raw grid from namelist input
      CALL COPY_GEOHYBGRID(gi, grid_in)
      CALL COPY_GEOHYBGRID(gg, grid_out)

      ! CHECK INPUT GRID
      CALL RGMSG(substr, RGMLIC, ' ... CHECKING IN-GRID ...')
      CALL CHECK_GEOHYBGRID(gi)
      CALL RGMSG(substr, RGMLIC, ' ... O.K.')

      ! CHECK OUTPUT GRID
      IF (lx .OR. ly) THEN
         CALL RGMSG(substr, RGMLIC, ' ... SWITCHING ...')
         CALL SWITCH_GEOHYBGRID(gg, lx, ly, lz)
      ENDIF

      CALL RGMSG(substr, RGMLIC, ' ... CHECKING OUT-GRID...')
      CALL CHECK_GEOHYBGRID(gg)
      CALL RGMSG(substr, RGMLIC, ' ... O.K.')

      ! TIME BALANCING (TODO ???)
      CALL RGMSG(substr, RGMLVM, '>>> BALANCING TIME AXIS ...')
      CALL BALANCE_GEOHYBGRID_TIME(gi, gg, lint)
      CALL RGMSG(substr, RGMLVM, '<<< ... END BALANCING TIME AXIS !')

      ! SURFACE PRESSURE
      CALL RGMSG(substr, RGMLVM, '>>> BALANCING SURFACE PRESSURE ...')
      CALL BALANCE_GEOHYBGRID_PS(status, gi, gg)
      CALL RGMSG(substr, RGMLVM, '<<< ... END BALANCING SURFACE PRESSURE !')

      IF ( (.not. ly .AND. .not. lx ) .and. lz) THEN
         CALL RGMSG(substr, RGMLIC, ' ... SWITCHING OUTPUT GRID ...')
         CALL SWITCH_GEOHYBGRID(gg, lx, ly, lz)
      ENDIF

      CALL RGMSG(substr, RGMLVM, '>>> SORTING INPUT GRID ...')
      CALL SORT_GEOHYBGRID(gi, gis, gix)
      CALL RGMSG(substr, RGMLVM, '<<< ... END SORTING INPUT GRID !')

      CALL RGMSG(substr, RGMLVM, '>>> SORTING OUTPUT GRID ...')
      CALL SORT_GEOHYBGRID(gg, ggs, ggx)
      CALL RGMSG(substr, RGMLVM, '<<< ... END SORTING OUTPUT GRID !')

      CALL RGMSG(substr, RGMLVM, '>>> COMPLETING INPUT GRID ...')
      CALL COMPLETE_GEOHYBGRID(gis, gix)
      CALL RGMSG(substr, RGMLVM, '<<< ... END COMPLETING INPUT GRID !')

      CALL RGMSG(substr, RGMLVM, '>>> COMPLETING OUTPUT GRID ...')
      CALL COMPLETE_GEOHYBGRID(ggs, ggx)
      CALL RGMSG(substr, RGMLVM, '<<< ... END COMPLETING OUTPUT GRID !')

      ! BALANCE GRID
      CALL RGMSG(substr, RGMLVM, '>>> BALANCING INPUT/OUTPUT GRID ...')
      CALL BALANCE_GEOHYBGRID(gis, ggs)
      CALL BALANCE_GEOHYBGRID(gix, ggx)
      CALL RGMSG(substr, RGMLVM, '<<< END BALANCING INPUT/OUTPUT GRID !')

      ! CONVERT HYBRID GRIDS TO AXES LIST PAIRS
      CALL RGMSG(substr, RGMLVM, '>>> CONSTRUCTING REGRIDDING AXES ...')
      CALL GEOHYBGRID_AXES(gis, sax, ggs, dax, lp, ly)
      CALL RGMSG(substr, RGMLVM, '<<< END CONSTRUCTING REGRIDDING AXES !')

      ! GET INPUT VAR
      mvars = 0

      ! *********************************************

      DO i=1, nvars  ! LOOP OVER VARIABLE NAMES
         ! CHECK
         CALL RGMSG(substr, RGMLIC, '  ... CHECKING FOR GRID CONFORMITY ...')
         CALL CHECK_NCVAR_ON_GEOHYBGRID(tvar(i), gi, dims, axes, ok)

         ! NOT NEEDED HERE
         axes(:) = 0
         DEALLOCATE(dims, STAT=status)
         CALL ERRMSG(substr,status,4)
         NULLIFY(dims)
         !
         IF (ok) THEN
            ! SAVE IN LIST
            mvars = mvars + 1
            CALL COPY_NCVAR(xivar(mvars), tvar(i))
            RGT(mvars)    = RG_TYPE(i)      ! shift
            CALL  RGMSG(substr, RGMLVL,  '  ... OK ')
         ELSE
            CALL  RGMSG(substr, RGMLIC, &
                 '  ... NOT CONSISTENT WITH GRID! IGNORING !', .TRUE.)
         END IF
         ! CLEAN
      END DO        ! LOOP OVER VARIABLE NAMES

      ! *********************************************
      ! WHAT TO DO WITH IMPORTED VARIABLES
      ! ***************************************************
      ! ---------------------------------------------------
      CALL RGMSG(substr, RGMLVL, 'REGRIDDING MODE ...')
      ! ALLOCATE SPACE
      ALLOCATE(xovar(mvars), svaro(mvars), pvaro(mvars), STAT=status)
      CALL ERRMSG(substr,status,6)
      ALLOCATE(nai(mvars), STAT=status)
      CALL ERRMSG(substr,status,7)
      DO i=1, mvars             ! LOOP OVER ALL VALID VARIABLES
         CALL RGMSG(substr, RGMLVL, &
              'VARIABLE '''//TRIM(xivar(i)%name)//''' ...')
         ! CONVERT IDX fields to FRACTION
         IF ((RGT(i) == RG_IDX).OR.(RGT(i) == RG_IXF)) THEN
            ! get the number of index classes from the variable attributes,
            ! if given
            ixf(:) = 0
            DO ix=1,xivar(i)%natts
               IF (TRIM(xivar(i)%att(ix)%name) == 'mmig_ixf_range') THEN
                  ixf(1:2) = INT(xivar(i)%att(ix)%dat%vi(1:2))
               END IF
            END DO
            CALL RGMSG(substr, RGMLVMC, ' ... INDEX FRACTIONS ...')
            IF (ANY(ixf /= 0)) THEN
               CALL IDX2FRAC_NCVAR(xivar(i), qvari, ixf=ixf)
            ELSE
               CALL IDX2FRAC_NCVAR(xivar(i), qvari)
            END IF
            CALL RGMSG(substr, RGMLVMC, ' ... ->'''//TRIM(qvari%name)//'''')
         ELSE
            CALL COPY_NCVAR(qvari, xivar(i))
         END IF

         ! GET ORDER INFORMATION
         CALL RGMSG(substr, RGMLIC, ' ... CHECKING ...')
         CALL CHECK_NCVAR_ON_GEOHYBGRID(qvari, gis, dims, axes, ok)
         IF (.NOT. OK) &
              CALL  RGMSG(substr, RGMLIC, &
              '  ... NOT CONSISTENT WITH GRID GIS! IGNORING !', .TRUE.) 
         ! BALANCE OUTPUT VARIABLE
         CALL RGMSG(substr, RGMLIC, ' ... BALANCING ...')
         CALL BALANCE_GEOHYBGRID_NCVAR(qvari, axes, ggs, xovar(i))
         ! SORT VARIABLE ACCORDING TO GRID
         CALL RGMSG(substr, RGMLIC, ' ... SORTING ...')
         CALL SORT_GEOHYBGRID_NCVAR(qvari, gix, axes, svari)
         CALL BALANCE_GEOHYBGRID_NCVAR(svari, axes, ggs, svaro(i))

         ! PACK VARIABLE
         CALL RGMSG(substr, RGMLIC, ' ... PACKING ...')
         CALL PACK_GEOHYBGRID_NCVAR(svari, dims, axes, pvari)
         CALL BALANCE_GEOHYBGRID_NCVAR(pvari, axes, ggs, pvaro(i))

         CALL INIT_NARRAY(pvaro(i)%dat)
         ! FILL N-ARRAY
         CALL COPY_NARRAY(nai(i), pvari%dat)
         ! CLEAN
         CALL INIT_NCVAR(qvari)
         CALL INIT_NCVAR(svari)
         CALL INIT_NCVAR(pvari)
         DEALLOCATE(dims, STAT=status)
         CALL ERRMSG(substr,status,8)
         NULLIFY(dims)
         axes(:) = 0
         CALL RGMSG(substr, RGMLIC, ' ... DONE (VARIABLE) !')
      END DO
      !
      ! REGRID
      IF (mvars > 0) THEN

         CALL NREGRID(nai, sax, dax, nao, RGT(1:mvars), sovl, dovl, rcnt)
         ! OUTPUT STATISTICS
         IF (MSGMODE > MSGMODE_S .AND. llstatout) THEN
            CALL NREGRID_STAT(sax, dax, sovl, dovl, nai, nao, rcnt)
         END IF

         ! CLEAN
         DEALLOCATE(sovl, STAT=status)
         CALL ERRMSG(substr,status,9)
         NULLIFY(sovl)
         DEALLOCATE(dovl, STAT=status)
         CALL ERRMSG(substr,status,10)
         NULLIFY(dovl)
         DEALLOCATE(rcnt, STAT=status)
         CALL ERRMSG(substr,status,11)
         NULLIFY(rcnt)

      END IF

      DO i=1, mvars
         CALL RGMSG(substr, RGMLVL, &
              'VARIABLE '''//TRIM(xovar(i)%name)//''' ...')
         ! COPY N-ARRAY
         CALL COPY_NARRAY(pvaro(i)%dat, nao(i))   ! 'DATA'
         CALL CHECK_NCVAR_ON_GEOHYBGRID(svaro(i), ggs, dims, axes, ok)
         IF (.NOT. OK) &
              CALL  RGMSG(substr, RGMLIC, &
              '  ... NOT CONSISTENT WITH GRID GGS! IGNORING !', .TRUE.) 
         ! UNPACK VARIABLE
         CALL RGMSG(substr, RGMLIC, ' ... UN-PACKING ...')
         CALL PACK_GEOHYBGRID_NCVAR(pvaro(i), dims, axes, svaro(i), .true.)
         !  UN-SORT
         CALL RGMSG(substr, RGMLIC, ' ... UN-SORTING ...')
         CALL SORT_GEOHYBGRID_NCVAR(svaro(i), ggx, axes, qvaro, .true.)
         !  UN-IDX
         CALL INIT_NCVAR(xovar(i))

         IF (RGT(i) == RG_IDX) THEN
            CALL RGMSG(substr, RGMLVMC, ' ... MAXIMUM INDEX FRACTION ...')
            CALL MAXFRAC2IDX_NCVAR(qvaro, xovar(i))
            xovar(i)%name = TRIM(xivar(i)%name)
            CALL RGMSG(substr, RGMLVMC, &
                 ' ... ->'''//TRIM(xovar(i)%name)//'''')
         ELSE
            CALL COPY_NCVAR(xovar(i), qvaro)
         END IF

         ! CLEAN
         CALL INIT_NCVAR(pvaro(i))
         CALL INIT_NCVAR(svaro(i))
         CALL INIT_NCVAR(qvaro)
         DEALLOCATE(dims, STAT=status)
         CALL ERRMSG(substr,status,12)
         NULLIFY(dims)
         axes(:) = 0
         CALL RGMSG(substr, RGMLIC, ' ... DONE (VARIABLE) !')
      END DO

      ! CLEAN
      DO i=1, mvars
         CALL INIT_NARRAY(nai(i))
         CALL INIT_NARRAY(nao(i))
      END DO
      IF (mvars > 0) THEN
         DEALLOCATE(nai, STAT=status)
         CALL ERRMSG(substr,status,13)
         NULLIFY(nai)
         DEALLOCATE(nao, STAT=status)
         CALL ERRMSG(substr,status,14)
         NULLIFY(nao)
         DEALLOCATE(svaro, pvaro, STAT=status)
         CALL ERRMSG(substr,status,15)
         NULLIFY(svaro)
         NULLIFY(pvaro)
      END IF
      !
      ! UN-SORT COMPLETED AND BALANCED GRID
      CALL SORT_GEOHYBGRID(ggs, go, ggx, .true.)
      !
      ! OUTPUT GRID AND VARIABLES
      IF (TRIM(go%file) /= '') THEN
         go%file = TRIM(grid_out%file)    ! name of output-file
         go%t = grid_out%t                ! output time step
         CALL EXPORT_GEOHYBGRID(go)
         ! EXPORT GLOBL ATTRIBUTE, ONLY AT FIRST TIME STEP
         IF (PRESENT(lfirsto)) THEN
            IF (lfirsto) THEN

               ! GET OLD ATTRIBUTE
               CALL IMPORT_NCATT(rggatt, varid=NF90_GLOBAL   &
                    ,attname = 'RG_HISTORY'     &
                    ,file = TRIM(go%file), lnostop=.true.)
               ! CHECK IF EMPTY OR CHAR-ATT
               IF ((rggatt%dat%type == VTYPE_UNDEF).OR.  &
                    (rggatt%dat%type == VTYPE_CHAR)) THEN
                  ! APPEND NEW DATA
                  CALL CAT_NARRAY(rggatt%dat, grid_out%att%dat)
                  rggatt%xtype = NF90_CHAR
                  rggatt%len   = rggatt%dat%dim(1)
               ELSE  ! EXISTS WITH NON-CHAR TYPE
                  CALL RGMSG(substr, RGMLW, &
                       'ATTRIBUTE '''//TRIM(rggatt%name)//'''')
                  CALL RGMSG(substr, RGMLWC, &
                       'EXISTS ALREADY, BUT IS NOT OF TYPE CHAR !')
               END IF
               CALL EXPORT_NCATT(rggatt, file=TRIM(go%file)  &
                    ,clobber=.true.)
            END IF
         END IF
         ! EXPORT VARIABLES
         DO i=1, mvars
            xovar(i)%ustep = go%t
            CALL EXPORT_NCVAR(xovar(i), file=TRIM(go%file))
         END DO
         CALL INIT_NCATT(rggatt)
      END IF

      ! RETURN VALUES TO SUBROUTINE CALL
      CALL RGMSG(substr, RGMLVM, '... RETURNING REGRIDDED DATA')
      ALLOCATE(var(mvars), STAT=status)
      CALL ERRMSG(substr,status,16)
      DO i=1, mvars
         xovar(i)%ustep = go%t
         CALL COPY_NCVAR(var(i), xovar(i))
      END DO
      IF (PRESENT(grid_conv)) THEN
         CALL RGMSG(substr, RGMLVM, '... RETURNING OUTPUT GRID')
         CALL COPY_GEOHYBGRID(grid_conv, go)
      END IF

      ! CLEAN
      DO i=1, mvars
         CALL INIT_NCVAR(xovar(i))
      END DO
      DEALLOCATE(xovar, STAT=status)
      CALL ERRMSG(substr,status,17)
      NULLIFY(xovar)
      !
      CALL RGMSG(substr, RGMLVL, '... DONE (REGRIDDING MODE) !')

      ! CLEAN
      IF (ASSOCIATED(xivar)) THEN
         DO i=1, SIZE(xivar)
            CALL INIT_NCVAR(xivar(i))
         END DO
         ! DO NOT DEALLOCATE/NULLIFY 'xivar' HERE, SINCE IT
         ! HAS BEEN ALLOCATED OUTSIDE THE TIME LOOP !
      END IF

      IF (ASSOCIATED(sax)) THEN
         DO i=1, SIZE(sax)
            CALL INIT_AXIS(sax(i))
         END DO
         DEALLOCATE(sax, STAT=status)
         CALL ERRMSG(substr,status,18)
         NULLIFY(sax)
      END IF
      !
      IF (ASSOCIATED(dax)) THEN
         DO i=1, SIZE(dax)
            CALL INIT_AXIS(dax(i))
         END DO
         DEALLOCATE(dax, STAT=status)
         CALL ERRMSG(substr,status,19)
         NULLIFY(dax)
      END IF
      !
      CALL INIT_GEOHYBGRID(gis)
      CALL INIT_GEOHYBGRID(gix)
      CALL INIT_GEOHYBGRID(gg)
      CALL INIT_GEOHYBGRID(ggx)
      CALL INIT_GEOHYBGRID(ggs)
      CALL INIT_GEOHYBGRID(gi)
      CALL INIT_GEOHYBGRID(go)

      ! CLEAN
      IF (ASSOCIATED(ivar)) DEALLOCATE(ivar)
      IF (ASSOCIATED(ovar)) DEALLOCATE(ovar)
      IF (ASSOCIATED(scl))  DEALLOCATE(scl)
      IF (ASSOCIATED(RGT))  DEALLOCATE(RGT)
      NULLIFY(ivar, ovar, scl, RGT)

      IF (ASSOCIATED(xivar)) THEN
         ! CALL INIT_NCVAR(xivar(i)) HAS ALREADY BEEN DONE INSIDE
         ! (AT THE END OF) THE TIME LOOP
         DEALLOCATE(xivar, STAT=status)
         CALL ERRMSG(substr,status,20)
         NULLIFY(xivar)
      END IF

    END SUBROUTINE REGRID_CONTROL_WORK

  END SUBROUTINE REGRID_CONTROL
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE GEOHYBGRID_AXES(g, a, g2, a2, pflag, lhor)

    ! Note: NO CHECKING
    !       g, g2 must be 'ordered', 'complete', and 'consistent'

    IMPLICIT NONE

    ! I/O
    TYPE (t_geohybgrid), INTENT(IN)      :: g  ! GEO-HYBRID-GRID
    TYPE (t_axis), DIMENSION(:), POINTER :: a  ! LIST OF AXES FOR REGRIDDING
    TYPE (t_geohybgrid), INTENT(IN)     , OPTIONAL :: g2  ! GEO-HYBRID-GRID
    TYPE (t_axis), DIMENSION(:), POINTER, OPTIONAL :: a2  ! LIST OF AXES
    LOGICAL, OPTIONAL, INTENT(IN)        :: pflag  ! .true.: pressure axis
    ! .false. sigma-axis (default)
    ! .true.  horizontal interpolation requested?
    LOGICAL, OPTIONAL, INTENT(IN)        :: lhor
    ! .false. => switch off spherical geometry
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'GEOHYBGRID_AXES'
    INTEGER                        :: n      ! number of axes
    INTEGER                        :: status
    LOGICAL                        :: lp     ! local presure axis flag
    INTEGER                        :: i      ! counter
    INTEGER                        :: vtype
    INTEGER                        :: ndep_lat ! number of lat-dimension in list
    INTEGER                        :: ndep_lon ! number of lon-dimension in list
    INTEGER                        :: dpc      ! dependency counter
    TYPE (t_narray)                :: zps      ! local surface pressure
    LOGICAL                        :: llhor

    IF (PRESENT(lhor)) THEN
       llhor = lhor
    ELSE
       llhor = .TRUE.
    END IF

    ! CHECK SUBROUTINE CALL
    IF (PRESENT(g2).NEQV.PRESENT(a2)) THEN
       CALL RGMSG(substr, RGMLE, &
            '2nd GRID AND 2nd AXES LIST MUST BE PRESENT ON SUBROUTINE CALL !')
    END IF

    ! INIT
    IF (PRESENT(pflag)) THEN
       lp = pflag
    ELSE
       lp = .false.     ! DEFAULT: sigma coordinates
    END IF
    ndep_lon = 0
    ndep_lat = 0

    ! COUNT AXES
    n = 1             ! START WITH MINIMUM 1 AXIS (LAST AXIS IS 'FREE' DIM)

    IF (QDEF_NCVAR(g%loni)) n=n+1
    IF (QDEF_NCVAR(g%lati)) n=n+1
    IF (QDEF_NCVAR(g%hyai).OR.QDEF_NCVAR(g%hybi)) n=n+1

    ! ALLOCATE SPACE
    ALLOCATE(a(n), STAT=status)
    CALL ERRMSG(substr,status,1)
    DO i=1, n
       CALL INIT_AXIS(a(i))
    END DO

    IF (PRESENT(g2)) THEN
       ALLOCATE(a2(n), STAT=status)
       CALL ERRMSG(substr,status,2)
       DO i=1, n
          CALL INIT_AXIS(a2(i))
       END DO
    END IF

    ! ASSIGN DATA
    n = 0

    CALL RGMSG(substr, RGMLI, 'AXIS CONSTRUCTION ...')

    ! NOTE: THE ORDER OF TESTING THE FOLLOWING DIMENSIONS
    !       (HERE: LON, LAT, LEV) SHOULD NOT BE CHANGED !
    !       IT HAS TO BE CONSISTENT WITH OTHER ROUTINES

    ! LONGITUDE
    IF (QDEF_NCVAR(g%loni)) THEN
       CALL RGMSG(substr, RGMLIC, ' ... ADDING LONGITUDE AXIS ...')
       n=n+1
       a(n)%lm     = g%lonc     ! LONGITUDE IS MODULO AXIS
       a(n)%ndp    = 1          ! LONGITUDE IS ...
       ALLOCATE(a(n)%dep(1), STAT=status)
       CALL ERRMSG(substr,status,3)
       a(n)%dep(1) = n          ! ... INDEPENDENT
       CALL COPY_NARRAY(a(n)%dat, g%loni%dat)
       ndep_lon    = n
       IF (PRESENT(g2)) THEN
          IF (QDEF_NCVAR(g2%loni)) THEN
             CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
             a2(n)%lm  = g2%lonc
             a2(n)%ndp = 1
             ALLOCATE(a2(n)%dep(1), STAT=status)
             CALL ERRMSG(substr,status,4)
             a2(n)%dep(1) = n
             CALL COPY_NARRAY(a2(n)%dat, g2%loni%dat)
          ELSE !UNDEFINED => INVARIANT
             CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
          END IF
       END IF
    END IF

    ! LATITUDE
    IF (QDEF_NCVAR(g%lati)) THEN
       CALL RGMSG(substr, RGMLIC, ' ... ADDING LATITUDE AXIS ...')
       n=n+1
       a(n)%lm     = .false.    ! LATITUDE IS NON-MODULO AXIS
       a(n)%ndp    = 1          ! LATITUDE IS ...
       ALLOCATE(a(n)%dep(1), STAT=status)
       CALL ERRMSG(substr,status,5)
       a(n)%dep(1) = n          ! ... INDEPENDENT
       CALL COPY_NARRAY(a(n)%dat, g%lati%dat)
       ndep_lat    = n
       IF (PRESENT(g2)) THEN
          IF (QDEF_NCVAR(g2%lati)) THEN
             CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
             a2(n)%lm  = .false.
             a2(n)%ndp = 1
             ALLOCATE(a2(n)%dep(1), STAT=status)
             CALL ERRMSG(substr,status,6)
             a2(n)%dep(1) = n
             CALL COPY_NARRAY(a2(n)%dat, g2%lati%dat)
          ELSE !UNDEFINED => INVARIANT
             CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
          END IF
       END IF
       !

       IF (llhor) THEN
          ! TAKE INTO ACCOUNT SPHERICAL GEOMETRY ...
          vtype = a(n)%dat%type
          SELECT CASE(vtype)
          CASE(VTYPE_REAL)
             a(n)%dat%vr = COS(((a(n)%dat%vr - REAL(90., SP))/ &
                  REAL(180., SP))*REAL(PI,SP))
          CASE(VTYPE_DOUBLE)
             a(n)%dat%vd = COS(((a(n)%dat%vd - REAL(90., DP))/ &
                  REAL(180., DP))*PI)
          CASE(VTYPE_INT)
             CALL RGMSG(substr, RGMLE, &
                  'LATITUDE AXIS OF TYPE INTEGER IS NOT SUPPORTED !')
          CASE(VTYPE_BYTE)
             CALL RGMSG(substr, RGMLE, &
                  'LATITUDE AXIS OF TYPE BYTE IS NOT SUPPORTED !')
          CASE(VTYPE_CHAR)
             CALL RGMSG(substr, RGMLE, &
                  'LATITUDE AXIS OF TYPE CHAR IS NOT SUPPORTED !')
          CASE(VTYPE_UNDEF)
             CALL RGMSG(substr, RGMLE, &
                  'LATITUDE AXIS IS UNDEFINED !')
          CASE DEFAULT
             CALL RGMSG(substr, RGMLE, &
                  'TYPE OF LATITUDE AXIS IS NOT RECOGNIZED !')
          END SELECT
          !
          ! ... ALSO FOR 2nd GRID
          IF (PRESENT(g2)) THEN
             IF (QDEF_NCVAR(g2%lati)) THEN
                vtype = a2(n)%dat%type
                SELECT CASE(vtype)
                CASE(VTYPE_REAL)
                   a2(n)%dat%vr = COS(((a2(n)%dat%vr - REAL(90., SP))/ &
                        REAL(180., SP))*REAL(PI,SP))
                CASE(VTYPE_DOUBLE)
                   a2(n)%dat%vd = COS(((a2(n)%dat%vd - REAL(90., DP))/ &
                        REAL(180., DP))*PI)
                CASE(VTYPE_INT)
                   CALL RGMSG(substr, RGMLE, &
                        'LATITUDE AXIS OF TYPE INTEGER IS NOT SUPPORTED !')
                CASE(VTYPE_BYTE)
                   CALL RGMSG(substr, RGMLE, &
                        'LATITUDE AXIS OF TYPE BYTE IS NOT SUPPORTED !')
                CASE(VTYPE_CHAR)
                   CALL RGMSG(substr, RGMLE, &
                        'LATITUDE AXIS OF TYPE CHAR IS NOT SUPPORTED !')
                CASE(VTYPE_UNDEF)
                   CALL RGMSG(substr, RGMLE, &
                        'LATITUDE AXIS IS UNDEFINED !')
                CASE DEFAULT
                   CALL RGMSG(substr, RGMLE, &
                        'TYPE OF LATITUDE AXIS IS NOT RECOGNIZED !')
                END SELECT
             END IF
          END IF
       END IF
    ENDIF

    ! LEVELS
    IF (QDEF_NCVAR(g%hyai).OR.QDEF_NCVAR(g%hybi)) THEN
       CALL RGMSG(substr, RGMLIC, ' ... ADDING VERTICAL AXIS ...')
       n=n+1
       a(n)%lm     = .false.     ! VERTICAL AXIS IS NON-MODULO AXIS

       ! CALCULATE AXIS DATA
       IF (lp) THEN
          CALL RGMSG(substr, RGMLIC, '     -> PRESSURE AXIS ...')
       ELSE
          CALL RGMSG(substr, RGMLIC, '     -> DIMENSIONLESS AXIS ...')
       END IF
       CALL PS2PS(g%ps, zps)
       CALL H2PSIG(a(n)%dat,g%hyai%dat,g%hybi%dat,zps,g%p0%dat,lp)
       CALL INIT_NARRAY(zps)

       ! SET DEPENDENCIES
       a(n)%ndp = a(n)%dat%n
       ALLOCATE(a(n)%dep(a(n)%ndp), STAT=status)
       CALL ERRMSG(substr,status,7)
       a(n)%dep(:) = 0
       dpc = 1
       a(n)%dep(dpc) = n
       !
       IF (a(n)%ndp > 1) THEN
          DO i=1, g%ps%ndims
             ! CHECK DIM LENGTH AND NAME/ID
             IF (QCMP_NCDIM(g%ps%dim(i), g%lonm%dim(1)) > 1) THEN ! LONGITUDE
                dpc = dpc + 1
                a(n)%dep(dpc) = ndep_lon
                CALL RGMSG(substr, RGMLIC, '     -> DEPENDING ON LONGITUDE ...')
             END IF
             IF (QCMP_NCDIM(g%ps%dim(i), g%latm%dim(1)) > 1) THEN ! LATITUDE
                dpc = dpc + 1
                a(n)%dep(dpc) = ndep_lat
                CALL RGMSG(substr, RGMLIC, '     -> DEPENDING ON LATITUDE ...')
             END IF
          END DO

          IF ((dpc > 1).AND.(dpc /= a(n)%ndp)) THEN
             CALL RGMSG(substr, RGMLE, 'DEPENDENCY MISMATCH !')
          END IF
       END IF

       ! 2nd GRID
       IF (PRESENT(g2)) THEN
          IF (QDEF_NCVAR(g2%hyai).OR.QDEF_NCVAR(g2%hybi)) THEN
             CALL RGMSG(substr, RGMLIC, ' ... DIFFERENT FOR OUTPUT GRID ...')
             a2(n)%lm     = .false.     ! VERTICAL AXIS IS NON-MODULO AXIS

             ! CALCULATE AXIS DATA
             CALL PS2PS(g2%ps, zps)
             CALL H2PSIG(a2(n)%dat,g2%hyai%dat,g2%hybi%dat, &
                  zps,g2%p0%dat,lp)
             CALL INIT_NARRAY(zps)

             ! SET DEPENDENCIES
             a2(n)%ndp = a2(n)%dat%n
             ALLOCATE(a2(n)%dep(a2(n)%ndp), STAT=status)
             CALL ERRMSG(substr,status,8)
             a2(n)%dep(:) = 0
             dpc = 1
             a2(n)%dep(dpc) = n
             !
             IF (a2(n)%ndp > 1) THEN
                DO i=1, g2%ps%ndims
                   ! CHECK DIM LENGTH AND NAME/ID
                   IF (QDEF_NCVAR(g2%lonm)) THEN
                      IF (QCMP_NCDIM(g2%ps%dim(i), g2%lonm%dim(1)) > 1) THEN
                         ! LONGITUDE
                         dpc = dpc + 1
                         a2(n)%dep(dpc) = ndep_lon
                         CALL RGMSG(substr, RGMLIC, &
                              '     -> DEPENDING ON LONGITUDE ...')
                      END IF
                   END IF
                   IF (QDEF_NCVAR(g2%latm)) THEN
                      IF (QCMP_NCDIM(g2%ps%dim(i), g2%latm%dim(1)) > 1) THEN
                         ! LATITUDE
                         dpc = dpc + 1
                         a2(n)%dep(dpc) = ndep_lat
                         CALL RGMSG(substr, RGMLIC, &
                              '     -> DEPENDING ON LATITUDE ...')
                      END IF
                   END IF
                END DO
                IF ((dpc > 1).AND.(dpc /= a2(n)%ndp)) THEN
                   CALL RGMSG(substr, RGMLE, 'DEPENDENCY MISMATCH !')
                END IF
             END IF
          ELSE ! UNDEFINED => INVARIANT
             CALL RGMSG(substr, RGMLIC, ' ... INVARIANT ...')
          END IF
       END IF ! 2nd GRID

    END IF ! LEVELS

    CALL RGMSG(substr, RGMLIC, '... END AXES CONSTRUCTION !')

  END SUBROUTINE GEOHYBGRID_AXES
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE PS2PS(var, na)

    IMPLICIT NONE

    ! I/O
    TYPE (t_ncvar) , INTENT(IN)  :: var
    TYPE (t_narray), INTENT(OUT) :: na

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'PS2PS'
    INTEGER :: i
    INTEGER :: uidpos
    INTEGER, DIMENSION(:), POINTER :: vec
    INTEGER :: dc
    INTEGER :: status
    INTEGER :: vtype

    NULLIFY(vec)

    uidpos = 0
    DO i=1, var%ndims
       IF (var%dim(i)%fuid) THEN
          uidpos = i
          EXIT
       END IF
    END DO

    IF (uidpos == 0) THEN  ! no unlimited ID
       CALL COPY_NARRAY(na, var%dat)
    ELSE                   ! 'remove' unlimited ID
       IF ((var%dim(uidpos)%len /= 1).OR.(var%dat%dim(uidpos) /= 1)) THEN
          CALL RGMSG(substr, RGMLE, &
               'DIMENSION LENGTH OF UNLIMITED DIMENSION MUST BE 1 !')
       END IF
       ALLOCATE(vec(var%ndims-1), STAT=status)
       CALL ERRMSG(substr,status,1)
       dc = 0
       DO i=1, var%ndims
          IF (.NOT.var%dim(i)%fuid) THEN
             dc = dc + 1
             vec(dc) = var%dim(i)%len
          END IF
       END DO
       vtype = var%dat%type
       CALL INIT_NARRAY(na, var%ndims-1, vec, vtype)
       SELECT CASE(vtype)
       CASE(VTYPE_REAL)
          na%vr(:) = var%dat%vr(:)
       CASE(VTYPE_DOUBLE)
          na%vd(:) = var%dat%vd(:)
       CASE(VTYPE_INT)
          na%vi(:) = var%dat%vi(:)
       CASE(VTYPE_BYTE)
          na%vb(:) = var%dat%vb(:)
       CASE(VTYPE_CHAR)
          na%vc(:) = var%dat%vc(:)
       CASE(VTYPE_UNDEF)
          CALL RGMSG(substr, RGMLE, 'N-ARRAY IS UNDEFINED !')
       CASE DEFAULT
          CALL RGMSG(substr, RGMLE, 'N-ARRAY IS UNRECOGNIZED !')
       END SELECT
       DEALLOCATE(vec,STAT=status)
       CALL ERRMSG(substr,status,2)
       NULLIFY(vec)
    END IF

  END SUBROUTINE PS2PS
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE BALANCE_GEOHYBGRID_PS(status, gi, go)

    IMPLICIT NONE

    ! I/O
    INTEGER,                    INTENT(OUT)   :: status
    TYPE (t_geohybgrid),        INTENT(INOUT) :: gi, go

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'BALANCE_GEOHYBGRID_PS'
    LOGICAL :: err
    LOGICAL :: linv    ! PRE-REGRID AVAILABLE go%ps

    err = (.NOT.QDEF_NCVAR(gi%ps)).AND.(.NOT.QDEF_NCVAR(go%ps))
    err = err .AND. ((QDEF_NCVAR(gi%hybi).OR.QDEF_NCVAR(gi%hybm)) .OR. &
         (QDEF_NCVAR(go%hybi).OR.QDEF_NCVAR(go%hybm)))

    IF (err) THEN
       CALL RGMSG(substr, RGMLE, &
            'HYBRID-B-COEFFICIENTS NEED I_PS AND/OR G_PS IN NAMELIST !')
    END IF

    err = (.NOT.QDEF_NCVAR(gi%p0)).AND.(.NOT.QDEF_NCVAR(go%p0))
    err = err .AND. ((QDEF_NCVAR(gi%hyai).OR.QDEF_NCVAR(gi%hyam)) .OR. &
         (QDEF_NCVAR(go%hyai).OR.QDEF_NCVAR(go%hyam)))

    IF (err) THEN
       CALL RGMSG(substr, RGMLE, &
            'HYBRID-A-COEFFICIENTS NEED I_P0 AND/OR G_P0 IN NAMELIST !')
    END IF

    ! RETURN, IF 2D
    err = (.NOT.QDEF_NCVAR(gi%p0)).AND.(.NOT.QDEF_NCVAR(gi%ps)).AND.     &
         (.NOT.QDEF_NCVAR(gi%hyai)).AND.(.NOT.QDEF_NCVAR(gi%hybi)).AND. &
         (.NOT.QDEF_NCVAR(gi%hyam)).AND.(.NOT.QDEF_NCVAR(gi%hybm))
    IF (err) THEN
       CALL RGMSG(substr, RGMLI, &
            'INPUT GRID IS 2-D! NO SURFACE PRESSURE REGRIDDING REQUIRED!')
       RETURN
    END IF

    ! BALANCE PO IN ANY CASE
    IF (QDEF_NCVAR(gi%p0).AND.(.NOT.QDEF_NCVAR(go%p0))) THEN
       CALL COPY_NCVAR(go%p0, gi%p0)
    ELSE
       IF (QDEF_NCVAR(go%p0).AND.(.NOT.QDEF_NCVAR(gi%p0))) THEN
          CALL COPY_NCVAR(gi%p0, go%p0)
       END IF
    END IF

    ! NOW PS ...

    ! CASE A: gi%ps AND go%ps BOTH AVAILABLE
    ! ADJUST TIME AXIS
    ! NOTE: THE TIME BALANCING (INCLUDING THAT FOR THE SURFACE PRESSURE
    !       TIME DIMENSION) IS PERFORMED IN
    !       SUBROUTINE BALANCE_GEOHYBGRID_TIME !

    IF (QDEF_NCVAR(gi%ps).AND.(QDEF_NCVAR(go%ps))) THEN
       RETURN
    END IF

    ! CASE B: gi%ps AND go%ps BOTH UNAVAILABLE
    ! NOTHING TO DO !

    IF (.NOT.QDEF_NCVAR(gi%ps).AND.(.NOT.QDEF_NCVAR(go%ps))) THEN
       RETURN
    END IF

    ! SURFACE PRESSURE NEEDS TO BE PRE-REGRIDDED, IF NOT AVAILABLE
    ! CASE C: go%ps AVAILABLE, BUT ON WRONG HORIZONTAL GRID

    linv = QDEF_NCVAR(go%ps).AND.( &
         ((.NOT.QDEF_NCVAR(go%lonm)).AND.(.NOT.QDEF_NCVAR(go%loni))).OR. &
         ((.NOT.QDEF_NCVAR(go%latm)).AND.(.NOT.QDEF_NCVAR(go%lati))) )

    IF (QDEF_NCVAR(go%ps).AND. linv ) THEN
       CALL RGMSG(substr, RGMLE,  &
            'REGRIDDING 3-D DISTRIBUTIONS', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'ONTO A DESTINATION SURFACE PRESSURE COORDINATE', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'USING (AN) INVARIANT HORIZONTAL DIMENSION(S)', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'IS NOT POSSIBLE DUE TO A LACK OF INFORMATION', .false.)
       CALL RGMSG(substr, RGMLEC, &
            '(HORIZONTAL DESTINATION GRID)', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'FOR PRE-REGRIDDING THE DESTINATION SURFACE PRESSURE !', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'PLEASE PERFORM 2-D PRE-REGRIDDING OF SURFACE PRESSURE', .false.)
       CALL RGMSG(substr, RGMLEC, &
            'IN SEPARATE STEP! (NRGD)')
    END IF

    ! CASE D: gi%ps XOR go%ps NOT AVAILABLE
    IF (QDEF_NCVAR(gi%ps).AND.(.NOT.QDEF_NCVAR(go%ps))) THEN
       CALL REGRID_GEOHYBGRID_PS(gi, go)
    ELSE
       IF (QDEF_NCVAR(go%ps).AND.(.NOT.QDEF_NCVAR(gi%ps))) THEN
          CALL REGRID_GEOHYBGRID_PS(go, gi)
       END IF
    END IF

    status = 0

  END SUBROUTINE BALANCE_GEOHYBGRID_PS
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE REGRID_GEOHYBGRID_PS(gi, go)

    IMPLICIT NONE

    ! I/O
    TYPE (t_geohybgrid),      INTENT(IN)    :: gi
    TYPE (t_geohybgrid),      INTENT(INOUT) :: go

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER            :: substr = 'REGRID_GEOHYBGRID_PS'
    TYPE (t_geohybgrid)                    :: gih, goh, gihs, gohs
    TYPE (t_geohybgrid)                    :: gix, gox
    INTEGER, DIMENSION(:),         POINTER :: dims    ! order of grid dimensions
    INTEGER                                :: axes(3) ! grid axes of variable
    TYPE (t_axis),   DIMENSION(:), POINTER :: ai, ao  ! source and dest. axes
    TYPE (t_ncvar)                         :: psi, pso
    TYPE (t_ncvar)                         :: psis, psos
    TYPE (t_ncvar)                         :: psip, psop
    TYPE (t_narray), DIMENSION(:), POINTER :: nao       ! re-gridder I/O
    REAL (DP),       DIMENSION(:), POINTER :: sovl, dovl
    INTEGER,         DIMENSION(:), POINTER :: rcnt
    LOGICAL                                :: ok
    INTEGER                                :: i
    INTEGER                                :: status

    status = 3999

    ! INIT
    NULLIFY(dims)
    NULLIFY(ai)
    NULLIFY(ao)
    NULLIFY(nao)
    NULLIFY(sovl)
    NULLIFY(dovl)
    NULLIFY(rcnt)

    ! CREATE 2-D GRIDS
    CALL INIT_GEOHYBGRID(gih)
    CALL INIT_GEOHYBGRID(goh)
    !
    gih%file = TRIM(gi%file)
    gih%t    = gi%t
    !
    goh%file = TRIM(go%file)
    goh%t    = go%t
    !
    CALL COPY_NCVAR(gih%lonm, gi%lonm)
    CALL COPY_NCVAR(gih%loni, gi%loni)
    !
    CALL COPY_NCVAR(gih%latm, gi%latm)
    CALL COPY_NCVAR(gih%lati, gi%lati)
    !
    CALL COPY_NCVAR(goh%lonm, go%lonm)
    CALL COPY_NCVAR(goh%loni, go%loni)
    !
    CALL COPY_NCVAR(goh%latm, go%latm)
    CALL COPY_NCVAR(goh%lati, go%lati)
    !
    CALL COPY_NCVAR(gih%timem, gi%timem)
    CALL COPY_NCVAR(gih%timei, gi%timei)
    !
    CALL COPY_NCVAR(goh%timem, go%timem)
    CALL COPY_NCVAR(goh%timei, go%timei)

    CALL SORT_GEOHYBGRID(gih, gihs, gix)
    CALL SORT_GEOHYBGRID(goh, gohs, gox)

    CALL COMPLETE_GEOHYBGRID(gihs, gix, go%ranges)
    CALL COMPLETE_GEOHYBGRID(gohs, gox, go%ranges)

    CALL BALANCE_GEOHYBGRID(gihs, gohs)
    CALL BALANCE_GEOHYBGRID(gix, gox)

    CALL GEOHYBGRID_AXES(gihs, ai, gohs, ao, .false.)

    CALL CHECK_NCVAR_ON_GEOHYBGRID(gi%ps, gihs, dims, axes, ok)
    IF (.NOT.ok) THEN
       CALL RGMSG(substr, RGMLE, 'PS NOT GRID CONFORM !')
    END IF

    CALL COPY_NCVAR(psi, gi%ps)
    CALL SORT_GEOHYBGRID_NCVAR(psi, gix, axes, psis)
    CALL BALANCE_GEOHYBGRID_NCVAR(psis, axes, gohs, psos)

    CALL PACK_GEOHYBGRID_NCVAR(psis, dims, axes, psip)
    CALL BALANCE_GEOHYBGRID_NCVAR(psip, axes, gohs, psop)

    DEALLOCATE(dims, STAT=status)
    CALL ERRMSG(substr,status,1)
    NULLIFY(dims)

    CALL NREGRID( (/ psip%dat /), ai, ao, nao, (/ RG_INT /), sovl, dovl, rcnt)
    IF (IAND(MSGMODE, MSGMODE_VM) == MSGMODE_VM) THEN
       CALL NREGRID_STAT(ai, ao, sovl, dovl, (/ psip%dat /), nao, rcnt)
    END IF

    CALL COPY_NARRAY(psop%dat, nao(1))
    CALL CHECK_NCVAR_ON_GEOHYBGRID(psos, gohs, dims, axes, ok)

    IF (.NOT.ok) THEN
       CALL RGMSG(substr, RGMLE, 'PS NOT GRID CONFORM !')
    END IF

    CALL PACK_GEOHYBGRID_NCVAR(psop, dims, axes, psos, .true.)

    CALL SORT_GEOHYBGRID_NCVAR(psos, gox, axes, pso, .true.)

    CALL COPY_NCVAR(go%ps, pso)

    ! CLEAN UP
    DEALLOCATE(dims, STAT=status)
    CALL ERRMSG(substr,status,2)
    NULLIFY(dims)
    DEALLOCATE(sovl, dovl, rcnt, STAT=status)
    CALL ERRMSG(substr,status,3)
    NULLIFY(sovl, dovl, rcnt)
    CALL INIT_NCVAR(psi)
    CALL INIT_NCVAR(pso)
    CALL INIT_NCVAR(psis)
    CALL INIT_NCVAR(psos)
    CALL INIT_NCVAR(psip)
    CALL INIT_NCVAR(psop)
    CALL INIT_GEOHYBGRID(gih)
    CALL INIT_GEOHYBGRID(goh)
    CALL INIT_GEOHYBGRID(gihs)
    CALL INIT_GEOHYBGRID(gohs)
    CALL INIT_GEOHYBGRID(gix)
    CALL INIT_GEOHYBGRID(gox)
    DO i=1, SIZE(nao)
       CALL INIT_NARRAY(nao(i))
    END DO
    DEALLOCATE(nao, STAT=status)
    CALL ERRMSG(substr,status,4)
    NULLIFY(nao)
    DO i=1, SIZE(ai)
       CALL INIT_AXIS(ai(i))
       CALL INIT_AXIS(ao(i))
    END DO
    DEALLOCATE(ai, ao, STAT=status)
    CALL ERRMSG(substr,status,5)
    NULLIFY(ai,ao)

  END SUBROUTINE REGRID_GEOHYBGRID_PS
  ! ------------------------------------------------------------------

! ******************************************************************
END MODULE messy_main_grid_trafo_nrgd
! ******************************************************************
