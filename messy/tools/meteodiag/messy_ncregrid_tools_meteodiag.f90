! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_TOOLS_METEODIAG
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, October 2002
! ******************************************************************

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  ! USER INTERFACE
  ! - DATA IMPORT
  PUBLIC  :: RGTOOL_METEODIAG_READ_NCVAR  ! READ ONE FIELD WITH NCREGRID
  PUBLIC  :: RGTOOL_METEODIAG_SAVE_NCVAR  ! SAVE ONE FIELD WITH NCREGRID


CONTAINS

! ####### PUBLIC ROUTINES ##########################################

! ------------------------------------------------------------------
SUBROUTINE RGTOOL_METEODIAG_READ_NCVAR(modstr, vname, t, dat       &
                               , lrg, lrgx, lrgy, lrgz, lok, grid &  ! OPTIONAL
                               , hyam, hybm, p0, ps         &  ! OPTIONAL
                               , hyai, hybi                 &  ! OPTIONAL
                               , latm, lonm, lati, loni     &  ! OPTIONAL
                               )

  ! PERFORMS ONE NCREGRID STEP FOR ONE FIELD
  ! AND RETURNS DATA AS 4D ARRAY (x,z,n,y),
  ! AND OPTIONALLY GRID-STRUCTURE AS ARRAYS
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, October 2002

  ! MESSY
  USE messy_main_tools,      ONLY: find_next_free_unit
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf, ONLY: ncvar, init_ncvar
  USE messy_ncregrid_geohyb, ONLY: geohybgrid, init_geohybgrid, COPY_GEOHYBGRID
  USE messy_ncregrid_tools,  ONLY: rgtool_read_ncvar, rgtool_convert &
                                 , rgtool_g2c

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, PRESENT, SIZE, TRIM

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)             :: modstr    ! calling module
  CHARACTER(LEN=*), INTENT(IN)             :: vname     ! name of variable
  INTEGER,          INTENT(IN)             :: t         ! netCDF time step
  REAL(dp), DIMENSION(:,:,:,:), POINTER    :: dat       ! (local) data field
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrg       ! regrid really ?
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgx      ! regrid in x
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgy      ! regrid in y
  LOGICAL,          INTENT(IN),   OPTIONAL :: lrgz      ! regrid in z
  LOGICAL,          INTENT(OUT),  OPTIONAL :: lok       ! OK?
  TYPE(geohybgrid),  OPTIONAL              :: grid      ! grid-structure
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyam  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybm  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: p0    ! reference pressure
  REAL(dp), DIMENSION(:,:), POINTER,  OPTIONAL :: ps    ! (local) surf. press.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hyai  ! hybrid-A-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: hybi  ! hybrid-B-coeff.
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: latm  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lonm  ! longitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: lati  ! latitude
  REAL(dp), DIMENSION(:),   POINTER,  OPTIONAL :: loni  ! longitude
  TYPE(ncvar)                        :: var       ! nc-variable

  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_METEODIAG_READ_NCVAR'
  INTEGER                            :: iou       ! logical I/O unit
!!$  INTEGER                            :: status    ! memory status
!!$  TYPE(ncvar)                        :: var       ! nc-variable
  TYPE(geohybgrid)                   :: lgrid      ! grid-structure
  !
  LOGICAL                            :: llrg, llrgx, llrgy, llrgz, llok
  INTEGER                            :: klat
  INTEGER                            :: klon
  INTEGER                            :: klev
  INTEGER                            :: kparam

  ! INITIALIZE
  ! I/O
  IF (ASSOCIATED(dat)) THEN
     DEALLOCATE(dat)
  END IF
  NULLIFY(dat)
  !
  IF (PRESENT(ps)) THEN
     IF (ASSOCIATED(ps)) THEN
        DEALLOCATE(ps)
     END IF
     NULLIFY(ps)
  END IF
  !
  IF (PRESENT(hyam)) THEN
     IF (ASSOCIATED(hyam)) THEN
        DEALLOCATE(hyam)
     END IF
     NULLIFY(hyam)
  END IF
  !
  IF (PRESENT(hybm)) THEN
     IF (ASSOCIATED(hybm)) THEN
        DEALLOCATE(hybm)
     END IF
     NULLIFY(hybm)
  END IF
  !
  IF (PRESENT(hyai)) THEN
     IF (ASSOCIATED(hyai)) THEN
        DEALLOCATE(hyai)
     END IF
     NULLIFY(hyai)
  END IF
  !
  IF (PRESENT(hybi)) THEN
     IF (ASSOCIATED(hybi)) THEN
        DEALLOCATE(hybi)
     END IF
     NULLIFY(hybi)
  END IF
  !
  IF (PRESENT(p0)) THEN
     IF (ASSOCIATED(p0)) THEN
        DEALLOCATE(p0)
     END IF
     NULLIFY(p0)
  END IF
  !
  IF (PRESENT(latm)) THEN
     IF (ASSOCIATED(latm)) THEN
        DEALLOCATE(latm)
     END IF
     NULLIFY(latm)
  END IF
  !
  IF (PRESENT(lonm)) THEN
     IF (ASSOCIATED(lonm)) THEN
        DEALLOCATE(lonm)
     END IF
     NULLIFY(lonm)
  END IF
  !
  IF (PRESENT(lati)) THEN
     IF (ASSOCIATED(lati)) THEN
        DEALLOCATE(lati)
     END IF
     NULLIFY(lati)
  END IF
  !
  IF (PRESENT(loni)) THEN
     IF (ASSOCIATED(loni)) THEN
        DEALLOCATE(loni)
     END IF
     NULLIFY(loni)
  END IF
  !
  ! LOCAL
  IF (PRESENT(lrg)) THEN
     llrg = lrg
  ELSE
     llrg = .true.
  END IF
  IF (PRESENT(lrgx)) THEN
     llrgx = lrgx
  ELSE
     llrgx = .true.
  END IF
  IF (PRESENT(lrgy)) THEN
     llrgy = lrgy
  ELSE
     llrgy = .true.
  END IF
  IF (PRESENT(lrgz)) THEN
     llrgz = lrgz
  ELSE
     llrgz = .true.
  END IF
  !
  iou = find_next_free_unit(100,200)
  CALL RGTOOL_READ_NCVAR(iou, TRIM(modstr)//'.nml', vname, t, var &
       ,lgrid, llrg, llrgx, llrgy, llrgz, llok)
  IF (llok) THEN
     CALL RGTOOL_CONVERT(var, dat, lgrid, order='xzny')
     klon = SIZE(dat, 1)
     klev = SIZE(dat, 2)
     kparam = SIZE(dat, 3)
     klat = SIZE(dat, 4)
     CALL RGTOOL_G2C(lgrid, hyam, hybm, p0, ps      &
          ,hyai, hybi                    &
          ,latm, lonm                    &
          ,lati, loni, klat, klon, klev)
  END IF
  ! Delete variable
  CALL INIT_NCVAR(var)
  IF (PRESENT(grid)) CALL COPY_GEOHYBGRID(grid,lgrid)
  IF (PRESENT(lok)) THEN
     lok = llok
  END IF


END SUBROUTINE RGTOOL_METEODIAG_READ_NCVAR
! ------------------------------------------------------------------

! ------------------------------------------------------------------
SUBROUTINE RGTOOL_METEODIAG_SAVE_NCVAR(file, vname, t, dat, grid) !       &

  ! DATA FROM 4D ARRAY (x,z,n,y) IS SAVED TO FILE
  !
  ! Author: Andreas Baumgaertner, MPICH, Mainz, July 2010

  ! MESSY
  USE messy_main_tools,      ONLY: find_next_free_unit
  ! NCREGRID
  USE messy_ncregrid_base
  USE messy_ncregrid_netcdf, ONLY: ncvar, init_ncvar, EXPORT_NCDIM, EXPORT_NCVAR
  USE messy_ncregrid_geohyb, ONLY: geohybgrid, init_geohybgrid, EXPORT_GEOHYBGRID
  USE messy_ncregrid_tools,  ONLY: rgtool_read_ncvar,rgtool_convert_dat2var &
                                 , rgtool_g2c

  IMPLICIT NONE

  INTRINSIC :: ASSOCIATED, PRESENT, SIZE, TRIM

  ! I/O
  CHARACTER(LEN=*), INTENT(IN)             :: file      ! filename
  CHARACTER(LEN=*), INTENT(IN)             :: vname     ! name of variable
  INTEGER,          INTENT(IN)             :: t         ! netCDF time step
  REAL(dp), DIMENSION(:,:), POINTER        :: dat       ! (local) data field
  TYPE(geohybgrid), INTENT(INOUT)          :: grid      ! grid-structure


  ! LOCAL
  CHARACTER(LEN=*), PARAMETER :: substr = 'RGTOOL_E5_SAVE_NCVAR'
  INTEGER                                  :: iou       ! logical I/O unit
  REAL(dp), DIMENSION(:,:,:,:), POINTER    :: ldat      ! (local) data field
  TYPE(ncvar)                              :: var       ! nc-variable
  !

  grid%file = file
  grid%t    = t
  CALL EXPORT_GEOHYBGRID(grid)

  ALLOCATE(ldat(SIZE(dat,1),SIZE(dat,2),1,1))
  ldat(:,:,1,1) = dat
  CALL RGTOOL_CONVERT_DAT2VAR(var, ldat, vname, grid, order='xyzn')
  DEALLOCATE(ldat)

  
  CALL EXPORT_NCVAR(var, file=file)

END SUBROUTINE RGTOOL_METEODIAG_SAVE_NCVAR
! ------------------------------------------------------------------


! ******************************************************************
! ------------------------------------------------------------------
END MODULE MESSY_NCREGRID_TOOLS_METEODIAG
! ------------------------------------------------------------------
! ******************************************************************
