MODULE mo_memory_streams

! Description:
!
! Initialisation and deinitialisation routine for the ECHAM memory buffer
!
! Method:
!
! Initialise:
! Call 'new_stream' for each each memory buffer to create it.
! Afterwards call the specific construct routines to fill each buffer
!
! Authors: A. Rhodin, MPI, May 2001
!          L. Kornblueh, MPI, April 2003, added port test
! 
  USE mo_doctor,        ONLY: nout
  USE mo_memory_base,   ONLY: new_stream, print_memory_use, print_sinfo, &
                              set_stream, suffix, delete_streams
#ifndef MESSY
  USE mo_sub_nml,       ONLY: set_stream_element_nml, set_stream_nml
#endif
  USE mo_memory_sp,     ONLY: sp  ,construct_sp  ,destruct_sp
  USE mo_memory_ls,     ONLY: ls  ,construct_ls  ,destruct_ls
  USE mo_memory_gl,     ONLY: gl  ,construct_gl  ,destruct_gl
  USE mo_memory_f,      ONLY: f   ,construct_f   ,destruct_f
  USE mo_memory_g1a,    ONLY: g1a ,construct_g1a ,destruct_g1a
  USE mo_memory_g1b,    ONLY: g1b ,construct_g1b ,destruct_g1b
  USE mo_memory_g2a,    ONLY: g2a ,construct_g2a ,destruct_g2a
  USE mo_memory_g2b,    ONLY: g2b ,construct_g2b ,destruct_g2b
  USE mo_memory_g3a,    ONLY:      construct_g3a
  USE mo_memory_g3b,    ONLY: g3b ,construct_g3b ,destruct_g3b
  USE mo_buffer_fft,    ONLY:      construct_fft ,destruct_fft
#ifndef MESSY
  USE mo_tracer,        ONLY: ntrac
#endif
  USE mo_control,       ONLY: ngl, nhgl, nlev, nlon, nmp1, nnp1, nsp,   &
                              ldebugmem ,loldrerun,                     &
                              nhf1 ,nhg1 ,nhg2 ,nhg3 ,nhgl1
  USE mo_decomposition, ONLY: lc => local_decomposition
  USE mo_mpi,           ONLY: p_pe, p_io
  USE mo_filename,      ONLY: standard_rerun_file, out_expname, &
                              ftyp => out_filetype
  USE mo_port_test,     ONLY: init_port_test
#ifdef MESSY
  USE mo_scan_buffer,   ONLY: m_bufscan
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init_memory
  PUBLIC :: free_memory
!==============================================================================
CONTAINS
!==============================================================================
  SUBROUTINE init_memory
    !--------------------------------
    ! declare standard output streams
    !--------------------------------
    CALL new_stream (sp  ,'sp' ,ftyp,lrerun=.FALSE.,lpost=.TRUE. ,post_suf='')
    CALL new_stream (ls  ,'ls' ,ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='')
    CALL new_stream (f   ,'f'  ,ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='')
    CALL new_stream (gl  ,'gl' ,ftyp,lrerun=.FALSE.,lpost=.TRUE. ,post_suf='')
    CALL new_stream (g1a ,'g1a',ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='')
    CALL new_stream (g1b ,'g1b',ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='')
    CALL new_stream (g2a ,'g2a',ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='')
    CALL new_stream (g2b ,'g2b',ftyp,lrerun=.FALSE.,lpost=.FALSE.,post_suf='')
    CALL new_stream (g3b ,'g3b',ftyp,lrerun=.FALSE.,lpost=.TRUE. ,post_suf='')
    !-------------------------------------------------------------
    ! set restart flags and restart file suffixes (old convention)
    !-------------------------------------------------------------
    IF(loldrerun) THEN
      standard_rerun_file='unit'
      CALL set_stream (f   ,lrerun=.TRUE.  ,rest_suf= suffix(nhf1 )) ! 31
      CALL set_stream (g1a ,lrerun=.TRUE.  ,rest_suf= suffix(nhg1 )) ! 35
      CALL set_stream (g2a ,lrerun=.TRUE.  ,rest_suf= suffix(nhg2 )) ! 36
      CALL set_stream (g3b ,lrerun=.TRUE.  ,rest_suf= suffix(nhg3 )) ! 37
      CALL set_stream (gl  ,lrerun=.TRUE.  ,rest_suf= suffix(nhgl1)) ! 32
    !-------------------------------------------------------------
    ! set restart flags and restart file suffixes (new convention)
    !-------------------------------------------------------------
    ELSE
      standard_rerun_file='rerun_'//out_expname
      CALL set_stream (f   ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 31
      CALL set_stream (g1a ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 35
      CALL set_stream (g2a ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 36
      CALL set_stream (g3b ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 37
      CALL set_stream (gl  ,lrerun=.TRUE.  ,rest_suf= '_echam') ! 32
    ENDIF
    !----------------------------------
    ! declare fields within each buffer
    !----------------------------------
    CALL construct_sp  (nlev,       lc%nsnm0, lc%snsp, &
                        nlev,       nnp1,     nsp)
    CALL construct_ls  (lc%nllev, lc%nllevp1, lc%nlnm0, lc%lnsp, &
                        nlev,                 nnp1,     nsp)
    CALL construct_f   (lc%nllev, lc%nllevp1, lc%nlm,  lc%nlat/2, &  
                        nlev,                 nmp1,    nhgl)
#ifndef MESSY
    CALL construct_g1b (lc%nproma,   lc%nlev, ntrac, lc%ngpblks, &
                        nlon,       nlev,    ntrac, ngl)
#else
       CALL construct_g1b (lc%nproma,   lc%nlev, lc%ngpblks, &
            nlon,       nlev, ngl)
#endif
    CALL construct_g2a (lc%nproma,   lc%nlev, lc%ngpblks, &
                        nlon,       nlev,    ngl)
    CALL construct_g2b (lc%nproma,   lc%nlev, lc%ngpblks, &
                        nlon,       nlev,    ngl)
    CALL construct_g3b

    CALL construct_g3a

#ifndef MESSY
    CALL construct_gl  (lc%nproma,   lc%nlev, ntrac, lc%ngpblks, &
                        nlon,       nlev,    ntrac, ngl)
    CALL construct_g1a (lc%nproma,   lc%nlev, ntrac, lc%ngpblks, &
                        nlon,       nlev,    ntrac, ngl)
#else
    CALL construct_gl  (lc%nproma,   lc%nlev, lc%ngpblks, &
         nlon,       nlev, ngl)

    CALL construct_g1a (lc%nproma,   lc%nlev, lc%ngpblks, &
         nlon,       nlev, ngl)
#endif
    CALL construct_fft (lc)

#ifndef MESSY
    CALL call_init_submodel_memory
#else
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL m_bufscan
    CALL messy_init_memory
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
    CALL init_port_test

#ifndef MESSY
    CALL set_stream_element_nml
    CALL set_stream_nml
#endif

    IF (ldebugmem) THEN
      WRITE (nout, '(/,a)') ' Global memory buffers:'
      WRITE (nout, '(a)') '    sp-buffer:  '
      CALL print_sinfo(sp)
      WRITE (nout, '(a)') '    ls-buffer:  '
      CALL print_sinfo(ls)
      WRITE (nout, '(a)') '     f-buffer:  '
      CALL print_sinfo(f)
      WRITE (nout, '(a)') '    gl-buffer:  '
      CALL print_sinfo(gl)
      WRITE (nout, '(a)') '   g1a-buffer:  '
      CALL print_sinfo(g1a)
      WRITE (nout, '(a)') '   g1b-buffer:  '
      CALL print_sinfo(g1b)
      WRITE (nout, '(a)') '   g2a-buffer:  '
      CALL print_sinfo(g2a)
      WRITE (nout, '(a)') '   g2b-buffer:  '
      CALL print_sinfo(g2b)
      WRITE (nout, '(a)') '   g3b-buffer:  '
      CALL print_sinfo(g3b)
    END IF

    IF (p_pe == p_io) THEN
      WRITE (nout, '(/,a)') ' Global memory buffers:'
      WRITE (nout, '(a)', ADVANCE='NO') '    sp-buffer:  '
      CALL print_memory_use(sp)
      WRITE (nout, '(a)', ADVANCE='NO') '     f-buffer:  '
      CALL print_memory_use(f)
      WRITE (nout, '(a)', ADVANCE='NO') '    gl-buffer:  '
      CALL print_memory_use(gl)
      WRITE (nout, '(a)', ADVANCE='NO') '   g1a-buffer:  '
      CALL print_memory_use(g1a)
      WRITE (nout, '(a)', ADVANCE='NO') '   g1b-buffer:  '
      CALL print_memory_use(g1b)
      WRITE (nout, '(a)', ADVANCE='NO') '   g2a-buffer:  '
      CALL print_memory_use(g2a)
      WRITE (nout, '(a)', ADVANCE='NO') '   g2b-buffer:  '
      CALL print_memory_use(g2b)
      WRITE (nout, '(a)', ADVANCE='NO') '   g3b-buffer:  '
      CALL print_memory_use(g3b)
    END IF

  END SUBROUTINE init_memory
!------------------------------------------------------------------------------
  !
  ! deallocate all memory allocated in the ECHAM memory buffer
  !
  SUBROUTINE free_memory

#ifdef MESSY
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL messy_free_memory             ! mz_pj_20030424
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif
    CALL destruct_sp
    CALL destruct_ls
    CALL destruct_f
    CALL destruct_gl
    CALL destruct_g1a
    CALL destruct_g1b
    CALL destruct_g2a
    CALL destruct_g2b
    CALL destruct_g3b
    CALL destruct_fft

#ifndef MESSY
    CALL call_free_submodel_memory
#endif

    CALL delete_streams

  END SUBROUTINE free_memory
!------------------------------------------------------------------------------
END MODULE mo_memory_streams
