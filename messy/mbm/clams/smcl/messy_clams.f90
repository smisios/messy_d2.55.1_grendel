!**************************************************************************
 MODULE messy_clams
!**************************************************************************
! CLaMS MODULE 
!----------- 

  USE messy_main_constants_mem, ONLY: DP
  USE messy_clams_global,       ONLY: dnparts_max     ! op_sb_20190806

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: clams_read_nml

!----------- 
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clams'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

  ! op_sb_20190819+
  REAL(DP), DIMENSION(:,:),   POINTER, PUBLIC :: POS   => NULL() ! parcel positions
  ! number of clams parcel per echam5 grid box
  REAL(DP), DIMENSION(:,:,:), POINTER, PUBLIC :: pc_g  => NULL() 
  ! index of pblh (global field)
  REAL(DP), DIMENSION(:,:),   POINTER,PUBLIC  :: khpbl  => NULL()
  ! ECHAM5 DRY grid mass
  REAL(dp), DIMENSION(:,:,:), POINTER,PUBLIC  :: grmass => NULL() 
  REAL(dp), DIMENSION(:,:,:), POINTER,PUBLIC  :: pmbox  => NULL() 
  INTEGER, PUBLIC :: nlev_i
  ! op_sb_20190819-
!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clams_read_nml(status, iou)


    USE messy_clams_global, ONLY: username, initfile, first_initfile, &
                                  met_dir, met_freq, met_prefix, &
                                  theta_dir, theta_prefix, &
                                  rres, maxspec, rres_shuffle, buffersize, & ! op_sb_20190806
                                  init_h2o_emac, lperpetuum, rank, ldiagout, timestep_clamsout, &
                                  levdotname, corrfile, use_3d, nparams, paramnames, &
                                  clams_gridding, clams_grid_verbose, n_cltr, cl_grid_tracers
    USE messy_main_tools,   ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou

    CHARACTER(LEN=*), PARAMETER :: substr='clams_read_nml'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    LOGICAL              :: l_print  ! write control output

    NAMELIST /CTRL/ initfile, first_initfile, met_dir, met_freq, met_prefix, &
                    theta_dir, theta_prefix, &
                    rres, rres_shuffle, username, buffersize, &
                    init_h2o_emac, lperpetuum, &
                    ldiagout, timestep_clamsout, &
                    levdotname, corrfile, use_3d, nparams, paramnames, &
                    clams_gridding, clams_grid_verbose, &
                    n_cltr, cl_grid_tracers

    status = 1 !ERROR

    if (rank==0 .and. ldiagout) then
       l_print = .true.
    else
       l_print = .false.
    endif

    ! Read namelist variables:

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, l_print)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, l_print)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    CALL read_nml_close(substr, iou, modstr, l_print)

    status = 0 !NO ERROR

  END SUBROUTINE clams_read_nml
!--------------------------------------------------------------------

!**************************************************************************
 END MODULE messy_clams
!**************************************************************************
