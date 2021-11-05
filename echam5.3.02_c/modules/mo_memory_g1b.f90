MODULE mo_memory_g1b

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
  USE mo_memory_base, ONLY: delete_stream, add_stream_element
#ifdef MESSY
  ! mz_pj_20040326+
  USE messy_main_tracer_mem_bi, ONLY: xtf, xtf_a
  ! mz_pj_20040326-
#endif

  IMPLICIT NONE
  SAVE   ! mz_jb_20040830
  PRIVATE

  PUBLIC :: construct_g1b ! construct the g1b table
  PUBLIC :: destruct_g1b  ! destruct  the g1b table

  PUBLIC :: g1b           ! the g1b table
#ifdef MESSY
  PUBLIC :: xtf, xtf_a    ! mz_pj_20040329
#endif

  ! declaration of predefined fields within this module 

  REAL(dp), POINTER, PUBLIC :: vof(:,:,:)
  REAL(dp), POINTER, PUBLIC :: df(:,:,:)
  REAL(dp), POINTER, PUBLIC :: tf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: alpsf(:,:)
  REAL(dp), POINTER, PUBLIC :: qf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xlf(:,:,:)
  REAL(dp), POINTER, PUBLIC :: xif(:,:,:)
#ifndef MESSY
  REAL(dp), POINTER, PUBLIC :: xtf(:,:,:,:)
#endif

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: g1b

CONTAINS

#ifndef MESSY
  SUBROUTINE construct_g1b (lnlon, lnlev, lntrac, lngl, &
                             nlon,  nlev,  ntrac,  ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lntrac, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ntrac,  ngl
#else
  SUBROUTINE construct_g1b (lnlon, lnlev, lngl, &
                             nlon,  nlev,  ngl) 

    INTEGER, INTENT (in) :: lnlon, lnlev, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ngl
#endif

    INTEGER :: nlevp1
    INTEGER :: dim1(3), dim1p(3)
    INTEGER :: dim2(2), dim2p(2)
#ifndef MESSY
    INTEGER :: dim3(4), dim3p(4)
#endif

    ! construct the g1b table
    !
    ! all information specific to this table is set in this subroutine

    nlevp1 = nlev  + 1

    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    dim1p = (/ lnlon, lnlev, lngl  /)
    dim1  = (/  nlon,  nlev,  ngl  /)

    dim2p = (/ lnlon, lngl  /)
    dim2  = (/  nlon,  ngl  /)

#ifndef MESSY
    dim3p = (/ lnlon, lnlev,  lntrac, lngl /)
    dim3  = (/  nlon,  nlev,   ntrac,  ngl /)
#endif

    CALL add_stream_element (g1b, 'vof',   vof,   dim1p, dim1)
    CALL add_stream_element (g1b, 'df',    df,    dim1p, dim1)
    CALL add_stream_element (g1b, 'tf',    tf,    dim1p, dim1)
    CALL add_stream_element (g1b, 'alpsf', alpsf, dim2p, dim2)
    CALL add_stream_element (g1b, 'qf',    qf,    dim1p, dim1)
    CALL add_stream_element (g1b, 'xlf',   xlf,   dim1p, dim1)
    CALL add_stream_element (g1b, 'xif',   xif,   dim1p, dim1)
#ifndef MESSY
    CALL add_stream_element (g1b, 'xtf',   xtf,   dim3p, dim3)
#endif

  END SUBROUTINE construct_g1b

  SUBROUTINE destruct_g1b

    CALL delete_stream (g1b)

  END SUBROUTINE destruct_g1b

END MODULE mo_memory_g1b
