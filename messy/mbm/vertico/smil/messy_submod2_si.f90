! NOTE: The preprocessor directives (_RI_*_) (<R>ank <I>dentifier) are used
!       to flip the meaning of ranks for the different spatial dimensions for
!       different basemodels, such as X (index jp), Y (index jrow), 
!       Z (jk) and N,M (jt, e.g., number of tracers, mode number) 
!       and C for full rank (colon, :) 
!       and I for special index and K for vertical index.
!
#ifdef ECHAM5
#define _RI_JZ_   1:nlev,jrow
#define _RI_YZ_   1:nlev,1:ngpblks
#endif
#ifdef COSMO
#define _RI_JZ_   jrow,1:nlev
#define _RI_YZ_   1:je,1:nlev
#endif
#if defined(BLANK) || defined(VERTICO)
#define _RI_JZ_   jrow,1:nlev
#define _RI_YZ_   1:nlat,1:nlev
#endif
! **********************************************************************
MODULE messy_submod2_si
! **********************************************************************
! This submodel is based on the submodel of the CHANNEL box model
! Authors: Patrick Joeckel, DLR, Oberpfaffenhofen (original code)
!          Astrid Kerkweg, UNI-MZ, Mainz (adapted to example submodel)
! **********************************************************************
  ! SMCL
  USE messy_submod2

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  ! MODULE VARIABLES
  REAL(DP), DIMENSION(:,:,:), POINTER :: f01      => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: f02physc => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: f02gend  => NULL()

  PUBLIC :: submod2_init_memory
  PUBLIC :: submod2_init_coupling
  PUBLIC :: submod2_physc
  PUBLIC :: submod2_global_end

CONTAINS

  ! --------------------------------------------------------------------
  SUBROUTINE submod2_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'submod2_init_memory'
    INTEGER :: status

    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt(substr, status)


    CALL new_channel_object(status, modstr, 'f02physc', p3=f02physc &
         , reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'f02physc' &
         , 'long_name', c = 'random field from submodel1')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'f02gend', p3=f02gend &
         , reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'f02gend' &
         , 'long_name', c = 'random field from submodel1')
    CALL channel_halt(substr, status)

  END SUBROUTINE submod2_init_memory
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod2_init_coupling

    ! BML/BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'submod2_init_coupling'
    INTEGER :: status

    CALL get_channel_object(status,'submod1','f01', p3=f01)
    CALL channel_halt(substr, status)

  END SUBROUTINE submod2_init_coupling
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod2_physc

    USE messy_main_grid_def_mem_bi, ONLY: nlon, jrow, nlev

    IMPLICIT NONE
    
    !f02physc(1:nlon,jrow, 1:nlev) = f01(1:nlon,jrow, 1:nlev)
    f02physc(1:nlon,_RI_JZ_) = f01(1:nlon,_RI_JZ_)

  END SUBROUTINE submod2_physc
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE submod2_global_end

    USE messy_main_grid_def_mem_bi, ONLY: nlon, nlat, nlev

    IMPLICIT NONE
    
    f02gend(1:nlon,_RI_YZ_) = f01(1:nlon,_RI_YZ_)

  END SUBROUTINE submod2_global_end
  ! --------------------------------------------------------------------

! **********************************************************************
END MODULE messy_submod2_si
! **********************************************************************
