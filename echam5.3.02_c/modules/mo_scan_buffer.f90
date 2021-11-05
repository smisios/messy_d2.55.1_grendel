#ifndef MESSY

MODULE mo_scan_buffer
  !
  ! Authors:
  !  ?                               original source and changes
  !  A. Rhodin,      DWD, June 2002, new subroutine cleanup_scanbuf
  !

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc=>local_decomposition

  IMPLICIT NONE

  !                                        ! set in  > used in

  REAL(dp), ALLOCATABLE :: dtm    (:,:,:)  ! ffti    ! dyn
  REAL(dp), ALLOCATABLE :: dtl    (:,:,:)  ! ffti    ! dyn
  REAL(dp), ALLOCATABLE :: dalpsl (:,:)    ! ffti    ! dyn
  REAL(dp), ALLOCATABLE :: dalpsm (:,:)    ! ffti    ! dyn
  REAL(dp), ALLOCATABLE :: dm     (:,:,:)  ! scan1sl ! fftd

  REAL(dp), ALLOCATABLE :: vo     (:,:,:)  ! ffti    ! dyn,scan1sl,statd,tf2n
  REAL(dp), ALLOCATABLE :: d      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: t      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: alps   (:,:)    ! ffti
  REAL(dp), ALLOCATABLE :: u      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: dudl   (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: v      (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: dvdl   (:,:,:)  ! ffti
  REAL(dp), ALLOCATABLE :: vol    (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE :: vom    (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE :: rh     (:,:,:)  ! fftd
  REAL(dp), ALLOCATABLE :: qte    (:,:,:)
  REAL(dp), ALLOCATABLE :: xlte   (:,:,:)
  REAL(dp), ALLOCATABLE :: xite   (:,:,:)
  REAL(dp), ALLOCATABLE :: xtte   (:,:,:,:)
  REAL(dp), ALLOCATABLE :: tte    (:,:,:)
  REAL(dp), ALLOCATABLE :: alpste (:,:)
  REAL(dp), ALLOCATABLE :: u0     (:,:)              ! fftd
  REAL(dp), ALLOCATABLE :: du0    (:,:)              ! fftd
  REAL(dp), ALLOCATABLE :: ul     (:,:)              ! fftd
  REAL(dp), ALLOCATABLE :: alnpr  (:,:,:)
  REAL(dp), ALLOCATABLE :: alpha  (:,:,:)
  REAL(dp), ALLOCATABLE :: vervel (:,:,:)

  LOGICAL, SAVE, PRIVATE :: lnot_used   = .TRUE.

CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE m_bufscan

    USE mo_tracer,        ONLY: ntrac

    INTEGER :: ngpblks, nlev, nproma, ngl

    IF (lnot_used) THEN

       ngl     = ldc% nglat
       ngpblks = ldc% ngpblks
       nlev    = ldc% nlev
       nproma  = ldc% nproma
       ! zero for test_scan_buffer
       ALLOCATE (dtm    (nproma,nlev,ngpblks))       ;dtm    = 0.0_dp
       ALLOCATE (dtl    (nproma,nlev,ngpblks))       ;dtl    = 0.0_dp
       ALLOCATE (dalpsl (nproma,ngpblks))            ;dalpsl = 0.0_dp
       ALLOCATE (dalpsm (nproma,ngpblks))            ;dalpsm = 0.0_dp
       ALLOCATE (dm     (nproma,nlev,ngpblks))       ;dm     = 0.0_dp

       ALLOCATE (vo     (nproma,nlev,ngpblks))       ;vo     = 0.0_dp
       ALLOCATE (d      (nproma,nlev,ngpblks))       ;d      = 0.0_dp
       ALLOCATE (t      (nproma,nlev,ngpblks))       ;t      = 0.0_dp
       ALLOCATE (alps   (nproma,ngpblks))            ;alps   = 0.0_dp
       ALLOCATE (u      (nproma,nlev,ngpblks))       ;u      = 0.0_dp
       ALLOCATE (dudl   (nproma,nlev,ngpblks))       ;dudl   = 0.0_dp
       ALLOCATE (v      (nproma,nlev,ngpblks))       ;v      = 0.0_dp
       ALLOCATE (dvdl   (nproma,nlev,ngpblks))       ;dvdl   = 0.0_dp
       ALLOCATE (vol    (nproma,nlev,ngpblks))       ;vol    = 0.0_dp
       ALLOCATE (vom    (nproma,nlev,ngpblks))       ;vom    = 0.0_dp
       ALLOCATE (rh     (nproma,nlev,ngpblks))       ;rh     = 0.0_dp
       ALLOCATE (qte    (nproma,nlev,ngpblks))       ;qte    = 0.0_dp
       ALLOCATE (xlte   (nproma,nlev,ngpblks))       ;xlte   = 0.0_dp
       ALLOCATE (xite   (nproma,nlev,ngpblks))       ;xite   = 0.0_dp
       ALLOCATE (xtte   (nproma,nlev,ntrac,ngpblks)) ;xtte   = 0.0_dp
       ALLOCATE (tte    (nproma,nlev,ngpblks))       ;tte    = 0.0_dp
       ALLOCATE (alpste (nproma,ngpblks))            ;alpste = 0.0_dp
       ALLOCATE (u0     (nlev,ngl))                  ;u0     = 0.0_dp
       ALLOCATE (du0    (nlev,ngl))                  ;du0    = 0.0_dp
       ALLOCATE (ul     (nlev,ngl))                  ;ul     = 0.0_dp
       ALLOCATE (alnpr  (nproma,nlev,ngpblks))       ;alnpr  = 0.0_dp
       ALLOCATE (alpha  (nproma,nlev,ngpblks))       ;alpha  = 0.0_dp
       ALLOCATE (vervel (nproma,nlev,ngpblks))       ;vervel = 0.0_dp

       lnot_used = .FALSE.

    ENDIF

  END SUBROUTINE m_bufscan
  !------------------------------------------------------------------------------
  SUBROUTINE cleanup_scanbuffer
    !------------------------------------
    ! deallocate variables in this module
    !------------------------------------

    IF (.NOT. lnot_used) THEN

       DEALLOCATE (dtm    )
       DEALLOCATE (dtl    )
       DEALLOCATE (dalpsl )
       DEALLOCATE (dalpsm )
       DEALLOCATE (dm     )
       DEALLOCATE (vo     )
       DEALLOCATE (d      )
       DEALLOCATE (t      )
       DEALLOCATE (alps   )
       DEALLOCATE (u      )
       DEALLOCATE (dudl   )
       DEALLOCATE (v      )
       DEALLOCATE (dvdl   )
       DEALLOCATE (vol    )
       DEALLOCATE (vom    )
       DEALLOCATE (rh     )
       DEALLOCATE (qte    )
       DEALLOCATE (xlte   )
       DEALLOCATE (xite   )
       DEALLOCATE (xtte   )
       DEALLOCATE (tte    )
       DEALLOCATE (alpste )
       DEALLOCATE (u0     )
       DEALLOCATE (du0    )
       DEALLOCATE (ul     )
       DEALLOCATE (alnpr  )
       DEALLOCATE (alpha  )
       DEALLOCATE (vervel )

       lnot_used   = .TRUE.

    END IF

  END SUBROUTINE cleanup_scanbuffer
  !----------------------------------------------------------------------------

END MODULE mo_scan_buffer

#else

MODULE mo_scan_buffer
  !
  ! Authors:
  !  ?                               original source and changes
  !  A. Rhodin,      DWD, June 2002, new subroutine cleanup_scanbuf
  !  P. Joeckel,     MPICH, Nov. 2004, converted to stream 

  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_tracer_mem_bi, ONLY: xtte, xtte_a, &
                                      xtte_c ! op_sb_20191007

  ! NOTE: xtte=>xtte, xtte_a=>xtte_a MUST be renamed and USEd,
  !       since ECHAM5 USEs the tracer memory from here;
  !       A local pointer-assignment in m_bufscan is not possible, since the
  !       memory for xtte and xtte_a is ALLOCATED AFTER m_bufscan
  !       is called:
  !         init_memory (mo_memory_streams.f90)
  !          |--> m_bufscan (mo_scan_buffer.f90)
  !          |--> messy_init_memory (messy_main_control_e5.f90)
  !               |--> main_tracer_init_memory (messy_main_tracer_bi.f90)
  !                    |--> ...
  !                          |--> setup_tracer_set
  !                          |--> get_tracer_set
  !

  IMPLICIT NONE
  SAVE

  !                                                             ! set in  ! used in
  REAL(DP), POINTER, DIMENSION(:,:,:)  :: dtm    => NULL()  ! ffti    ! dyn
  REAL(DP), POINTER, DIMENSION(:,:,:)  :: dtl    => NULL()  ! ffti    ! dyn
  REAL(DP), POINTER, DIMENSION(:,:)    :: dalpsl => NULL()  ! ffti    ! dyn
  REAL(DP), POINTER, DIMENSION(:,:)    :: dalpsm => NULL()  ! ffti    ! dyn
  REAL(DP), POINTER, DIMENSION(:,:,:)  :: dm     => NULL()  ! scan1sl ! fftd
  REAL(DP), POINTER, DIMENSION(:,:,:)  :: vo     => NULL()  ! ffti    ! dyn,scan1sl,statd,tf2n
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: d      => NULL()  ! ffti
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: t      => NULL()  ! ffti
  REAL(DP), POINTER, DIMENSION(:,:)     :: alps   => NULL()  ! ffti
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: u      => NULL()  ! ffti
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: dudl   => NULL()  ! ffti
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: v      => NULL()  ! ffti
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: dvdl   => NULL()  ! ffti
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: vol    => NULL()  ! fftd
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: vom    => NULL()  ! fftd
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: rh     => NULL()  ! fftd
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: qte    => NULL()
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: xlte   => NULL()
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: xite   => NULL()
!!$  REAL(DP), POINTER, DIMENSION(:,:,:,:) :: xtte   => NULL()
!!$  REAL(dp), POINTER, DIMENSION(:,:,:,:) :: xtte_a => NULL()
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: tte    => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)     :: alpste => NULL()
  REAL(DP), POINTER, DIMENSION(:,:)     :: u0     => NULL() ! fftd
  REAL(DP), POINTER, DIMENSION(:,:)     :: du0    => NULL() ! fftd
  REAL(DP), POINTER, DIMENSION(:,:)     :: ul     => NULL() ! fftd
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: alnpr  => NULL()
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: alpha  => NULL()
  REAL(DP), POINTER, DIMENSION(:,:,:)   :: vervel => NULL()

CONTAINS
  !------------------------------------------------------------------------------
  SUBROUTINE m_bufscan

    USE mo_decomposition, ONLY: dcl => local_decomposition
    USE mo_memory_base,   ONLY: t_stream, new_stream, add_stream_element &
                              , GRIDPOINT, default_stream_setting 

    USE mo_exception,     ONLY: finish
    
    ! LOCAL
    INTEGER                               :: ngpblks, nlev, nproma, ngl, nlat
    TYPE(t_stream),               POINTER :: stream_ptr
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: p4
    LOGICAL, PARAMETER                    :: lverbose = .FALSE.

    ngl     = dcl% nglat
    ngpblks = dcl% ngpblks
    nlev    = dcl% nlev
    nproma  = dcl% nproma
    nlat    = dcl% nlat
    
    ! CREAT NEW STREAM 
    CALL new_stream(stream_ptr, 'scnbuf', lpost=.false., lrerun=.false.)
    CALL default_stream_setting(stream_ptr, repr=GRIDPOINT, &
         lpost=.false., lrerun=.false.)
    
    ! ADD STREAM ELEMENTS (3D) nproma x nlev x ngpblks
    CALL add_stream_element(stream_ptr, 'dtm',  dtm  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'dtl',  dtl  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'dm',   dm  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'vo',   vo  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'd',    d  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 't',    t  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'u',    u  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'dudl', dudl  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'v',    v  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'dvdl', dvdl  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'vol',  vol  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'vom',  vom  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'rh',   rh  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'qte',  qte  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'xlte', xlte  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'xite', xite  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'tte',  tte  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'alnpr',alnpr  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'alpha',alpha  &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'vervel',vervel  &
         , verbose=lverbose)
    
    ! ADD STREAM ELEMENTS (2D) nproma x ngpblks
    CALL add_stream_element(stream_ptr, 'dalpsl', dalpsl &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'dalpsm', dalpsm &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'alps',   alps   &
         , verbose=lverbose)
    CALL add_stream_element(stream_ptr, 'alpste', alpste &
         , verbose=lverbose)
    
    ALLOCATE (u0     (nlev,ngl))                  ;u0     = 0.0_dp
    ALLOCATE (du0    (nlev,ngl))                  ;du0    = 0.0_dp
    ALLOCATE (ul     (nlev,ngl))                  ;ul     = 0.0_dp

  END SUBROUTINE m_bufscan
  !-------------------------------------------------------------------------
  SUBROUTINE cleanup_scanbuffer

    IMPLICIT NONE
    
    NULLIFY (dtm    )
    NULLIFY (dtl    )
    NULLIFY (dalpsl )
    NULLIFY (dalpsm )
    NULLIFY (dm     )
    NULLIFY (vo     )
    NULLIFY (d      )
    NULLIFY (t      )
    NULLIFY (alps   )
    NULLIFY (u      )
    NULLIFY (dudl   )
    NULLIFY (v      )
    NULLIFY (dvdl   )
    NULLIFY (vol    )
    NULLIFY (vom    )
    NULLIFY (rh     )
    NULLIFY (qte    )
    NULLIFY (xlte   )
    NULLIFY (xite   )
!!$    NULLIFY (xtte   )
!!$    NULLIFY (xtte_a )
    NULLIFY (tte    )
    NULLIFY (alpste )
    DEALLOCATE(u0); NULLIFY (u0     )
    DEALLOCATE(du0); NULLIFY (du0    )
    DEALLOCATE(ul); NULLIFY (ul     )
    NULLIFY (alnpr  )
    NULLIFY (alpha  )
    NULLIFY (vervel )
    
    NULLIFY (vo)
    NULLIFY (d)
    NULLIFY (t)
    NULLIFY (alps)
    NULLIFY (u)
    NULLIFY (dudl)
    NULLIFY (v)
    NULLIFY (dvdl)
    NULLIFY (vol)
    NULLIFY (vom)
    NULLIFY (rh)
    NULLIFY (qte)
    NULLIFY (xlte)
    NULLIFY (xite)
    NULLIFY (xtte)
    NULLIFY (xtte_a)
    NULLIFY (xtte_c) ! op_sb_20191007
    NULLIFY (tte)
    NULLIFY (alpste)
    NULLIFY (alnpr)
    NULLIFY (alpha)
    NULLIFY (vervel)
    
  END SUBROUTINE cleanup_scanbuffer
  !----------------------------------------------------------------------------

END MODULE mo_scan_buffer

#endif
