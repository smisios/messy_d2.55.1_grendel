MODULE mo_memory_gl
!------------------------------------------------------------------------------
  !
  ! Modules used
  !

  USE mo_kind,        ONLY: dp
  USE mo_linked_list, ONLY: t_stream
#ifndef MESSY
  USE mo_memory_base, ONLY: new_stream, delete_stream, add_stream_element, &
                            default_stream_setting, add_stream_reference
#else
  USE mo_memory_base, ONLY: new_stream, delete_stream, add_stream_element, &
                            default_stream_setting
#endif
  USE mo_netCDF,      ONLY: max_dim_name
#ifndef MESSY
  USE mo_tracdef,     ONLY: trlist
#else
  ! mz_pj_20040325+
  USE messy_main_tracer_mem_bi, ONLY: xt, xt_a ,pxt, pxt_a, &
                                      xt_c, pxt_c ! op_sb_20191007
  ! mz_pj_20040325-
#endif
  USE mo_filename,    ONLY: trac_filetype

  IMPLICIT NONE
!------------------------------------------------------------------------------
  !
  ! Public entities
  !

  PRIVATE
#ifdef MESSY
  ! mz_pj_20040325+
  PUBLIC :: xt, xt_a, xt_c ! op_sb_20191007
  ! mz_pj_20040325-
#endif

  PUBLIC :: construct_gl ! construct the gl table
  PUBLIC :: destruct_gl  ! destruct  the gl table

  PUBLIC :: gl           ! the gl table
#ifndef MESSY
  PUBLIC :: tracer       ! the tracer table
#endif
!------------------------------------------------------------------------------

  ! declaration of predefined fields within this module 

  REAL(dp), SAVE, POINTER, PUBLIC :: q(:,:,:)    ! mz_jb_20031205: SAVE added
  REAL(dp), SAVE, POINTER, PUBLIC :: xl(:,:,:)   ! mz_jb_20031205: SAVE added
  REAL(dp), SAVE, POINTER, PUBLIC :: xi(:,:,:)   ! mz_jb_20031205: SAVE added
#ifndef MESSY
  REAL(dp), SAVE, POINTER, PUBLIC :: xt(:,:,:,:) ! mz_ak_20030611: SAVE added
#endif
  REAL(dp), POINTER, PUBLIC :: lammp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: phimp(:,:,:)
  REAL(dp), POINTER, PUBLIC :: sigmp(:,:,:)

  ! declaration of table with 3d-field entries

  TYPE (t_stream), POINTER :: gl
#ifndef MESSY
  TYPE (t_stream), POINTER :: tracer
#endif

  ! private storage for tracer fields

#ifndef MESSY
  REAL(dp), POINTER :: pxt (:,:,:,:,:)
#endif

CONTAINS
!------------------------------------------------------------------------------
#ifndef MESSY
  SUBROUTINE construct_gl (lnlon, lnlev, lntrac, lngl, &
                            nlon,  nlev,  ntrac,  ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lntrac, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ntrac,  ngl
#else
  SUBROUTINE construct_gl (lnlon, lnlev, lngl, &
                            nlon,  nlev, ngl)

    INTEGER, INTENT (in) :: lnlon, lnlev, lngl
    INTEGER, INTENT (in) ::  nlon,  nlev,  ngl
#endif

    INTEGER                      :: dim1(3), dim1p(3)
#ifndef MESSY
    INTEGER                      :: dim2(4), dim2p(4)
    CHARACTER (len=max_dim_name) :: dim1n(3), dim2n(4)
    INTEGER                      :: i
    REAL(dp), POINTER            :: p3(:,:,:), p4(:,:,:,:)
#else
    CHARACTER (len=max_dim_name) :: dim1n(3)
#endif

    ! construct the gl table
    !
    ! all information specific to this table is set in this subroutine


    ! overwrite default entries for the predefined fields
    ! allocate the predefined fields

    ! assign pointers

    dim1p = (/  lnlon,  lnlev, lngl  /)
    dim1  = (/   nlon,   nlev,  ngl  /)
    dim1n = (/   "lon ","lev ","lat "/)

#ifndef MESSY
    dim2p = (/ lnlon,    lnlev,   lntrac,  lngl    /)
    dim2  = (/  nlon,     nlev,    ntrac,   ngl    /)
    dim2n = (/  "lon   ","lev   ","ntrac ","lat   "/)
#endif

    CALL default_stream_setting (gl ,dimnames=dim1n ,lrerun=.TRUE. ,table=128)

    CALL add_stream_element (gl,'q',    q,    code=133 ,longname='specific humidity'        ,units='kg/kg')
    CALL add_stream_element (gl,'xl',   xl,   code=153 ,longname='cloud water'              ,units='kg/kg')
    CALL add_stream_element (gl,'xi',   xi,   code=154 ,longname='cloud ice'                ,units='kg/kg')
    CALL add_stream_element (gl,'lammp',lammp,lpost=.FALSE., &
         lrerun=.FALSE.,contnorest=.TRUE.) 
    CALL add_stream_element (gl,'phimp',phimp,lpost=.FALSE., &
         lrerun=.FALSE.,contnorest=.TRUE.) 
    CALL add_stream_element (gl,'sigmp',sigmp,lpost=.FALSE., &
         lrerun=.FALSE.,contnorest=.TRUE.) 

#ifndef MESSY
    !
    ! Special handling for tracers
    !

    IF (ntrac > 0) THEN
      !
      ! Allocate a 5d-array (with dummy index 5) for tracers. This array 
      ! is referenced by the 4-d array 'XT' and by the 3-d arrays of 
      ! individual tracers.
      !
      ALLOCATE (pxt(lnlon, lnlev, lntrac, lngl, 1))
      !
      ! Set meta information on XT array.
      ! Set restart flag, obtain reference to memory info entry.
      !
      p4 => pxt(:,:,:,:,1)
      CALL add_stream_element (gl, 'xt', xt, dim2p, dim2,       &
                               dimnames   = dim2n,              &
                               lrerun     = trlist% oldrestart, &
                               contnorest = .TRUE.,             &   
                               mem_info   = trlist% mixt,       &
                               lpost      = .FALSE.,            &
                               p4=p4)
      !
      ! setup a new output stream
      !
      CALL new_stream (tracer ,'tracer',trac_filetype)
      !
      ! add entries for geopotential, log surface pressure, grid-box area
      !
      CALL add_stream_reference (tracer, 'geosp'   ,'g3b'   ,lpost=.TRUE.)
      CALL add_stream_reference (tracer, 'lsp'     ,'sp'    ,lpost=.TRUE.)
      CALL add_stream_reference (tracer, 'aps'     ,'g3b'   ,lpost=.TRUE.)    
      CALL add_stream_reference (tracer, 'gboxarea','geoloc',lpost=.TRUE.)
      !
      ! provide additional meta-information for individual tracers.
      !
      DO i = 1, lntrac
        p4 => pxt(:,:,i,:,:)

        CALL add_stream_element (tracer, trlist% ti(i)% fullname, p3,      &
                                 dim1p, dim1,                              &
                                 dimnames   = dim1n,                       &
                                 units      = trlist% ti(i)% units,        &
                                 longname   = trlist% ti(i)% longname,     &
                                 lrerun     = trlist% ti(i)% nrerun==1,    &
                                 contnorest = .TRUE.,                      &
                                 mem_info   = trlist% mi(i)% xt,           &
                                 lpost      = trlist% ti(i)% nwrite==1,    &
                                 code       = trlist% ti(i)% code,         &
                                 table      = trlist% ti(i)% table,        &
                                 bits       = trlist% ti(i)% gribbits,     &
                                 tracidx    = i,                           &
                                 p4         = p4)
      END DO
    ELSE
      NULLIFY (tracer)
      NULLIFY (pxt)
      CALL add_stream_element (gl,'xt', xt, dim2p, dim2, &
                               lpost   =.FALSE.,         &
                               lrerun  =.FALSE.,         &
                               mem_info=trlist% mixt)
    END IF
#endif

  END SUBROUTINE construct_gl
!------------------------------------------------------------------------------
  SUBROUTINE destruct_gl

    CALL delete_stream (gl)
#ifndef MESSY
    CALL delete_stream (tracer)

    IF(ASSOCIATED (pxt)) DEALLOCATE (pxt)
#endif

  END SUBROUTINE destruct_gl
!------------------------------------------------------------------------------
END MODULE mo_memory_gl
