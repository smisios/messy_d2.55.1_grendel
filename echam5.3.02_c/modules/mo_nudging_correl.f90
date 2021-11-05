MODULE MO_NUDGING_CORREL

#ifdef MESSY

  IMPLICIT NONE

  PRIVATE
!EOX
  ! !PUBLIC MEMBER FUNCTIONS:

  PUBLIC :: Nudg_Correl

CONTAINS

  SUBROUTINE Nudg_Correl

    ! !DESCRIPTION: 
    ! correlations between the nudging term and the tendency diagnostic terms

    ! !USES:
    USE mo_kind,                ONLY: dp
    USE mo_nudging_buffer,      ONLY: sdtac, svtac, sttac, sstac, nio_index
    USE mo_nudging_constants,   ONLY: ndunit
    USE mo_control,             ONLY: nlev,nsp
#ifdef OBSOLETE
    USE mo_post,                ONLY: lppt, lppp, lppvo, lppd
#endif
    USE mo_diag_tendency,       ONLY: &
         ndvor, nddiv, ndtem, ndprs, &! number of diagnostic terms
         pdvor, pddiv, pdtem, pdprs, &! accumulated terms of equations
         dio_index
    USE mo_mpi,            ONLY: p_parallel_io
    USE mo_exception,      ONLY: finish, message
    USE mo_time_conversion,ONLY: print_date, TC_PRN_NATIVE
    USE mo_time_control,   ONLY: next_date, get_interval_seconds, ev_putdata
    USE mo_spectral,       ONLY: corrsp
    USE mo_transpose,      ONLY: gather_sp
    USE mo_decomposition,  ONLY: gdc=>global_decomposition
    ! !INPUT PARAMETERS: 

!EOP

    INTEGER       :: jk, jt, isec1, isec2
    REAL(kind=dp) :: cc_tem(ndtem), cc_vor(ndvor), cc_div(nddiv), cc_prs(ndprs)
    REAL(kind=dp) :: cc_tem_l(nlev, ndtem), cc_vor_l(nlev, ndvor), &
                     cc_div_l(nlev, nddiv), cc_prs_l(nlev, ndprs)
    CHARACTER(len=256) :: mess
    REAL(dp), POINTER :: wrk_ptr(:,:,:,:) => NULL()
    REAL(dp), POINTER :: nudg_ptr(:,:,:)  => NULL()
    REAL(dp), POINTER :: p2a(:,:)         => NULL()
    REAL(dp), POINTER :: p2b(:,:)         => NULL()
    REAL(dp), POINTER :: p3a(:,:,:)       => NULL()
    REAL(dp), POINTER :: p3b(:,:,:)       => NULL()

    INTRINSIC trim

    isec1 = get_interval_seconds(ev_putdata(nio_index))
    isec2 = get_interval_seconds(ev_putdata(dio_index))

    CALL print_date(next_date,TC_PRN_NATIVE,mess=mess)
    WRITE(ndunit,'(a,a)') '###  ',TRIM(mess)

    IF (isec1 /= isec2) THEN
      CALL message('Nudg_Correl', &
        'accumulation interval mismatch (Nudging, TDiag) - NO correlation performed!')
      WRITE(ndunit,'(a)') &
        '### accumulation interval mismatch -- no correlation performed'
      RETURN
    END IF
!--------------------------------------------------------------------------
#ifdef OBSOLETE
    IF (lppd) THEN    ! divergence equation
  
      if (p_parallel_io) then
        ALLOCATE(wrk_ptr(nlev,2,nsp, NDDIV))
        ALLOCATE(nudg_ptr(nlev,2,nsp))
      else
        NULLIFY(wrk_ptr)
        NULLIFY(nudg_ptr)
      endif

      DO jt=1,NDDIV
        IF (ASSOCIATED(wrk_ptr)) then
          p3a => wrk_ptr(:,:,:,jt)
        ELSE
          NULLIFY(p3a)
        END IF
        p3b => pddiv(:,:,:,jt)
        call gather_sp(p3a, p3b, gdc) 
      END DO
      CALL gather_sp(nudg_ptr, sdtac, gdc)

      if (p_parallel_io) then
        WRITE(ndunit,'(a,20(1x,a7))') 'DIV_    ',&
          'dynAPCG', 'VerAdv ', 'VerDiff', 'GWDrag ', &
          'CuCall ', 'TimFilt', 'SemiImp', 'HorDiff', &
          'SumTend'
        cc_div(:)=0._dp
        DO jt=1,NDDIV
          cc_div(jt) = corrsp(nudg_ptr,wrk_ptr(:,:,:,jt))
        ENDDO
        WRITE(ndunit,'(a,20(1x,f7.3))') 'DIV_ALL ',cc_div(:)

        DO jk=1,nlev
          cc_div(:)=0._dp
          DO jt=1,NDDIV
            p2a => nudg_ptr(jk,:,:)
            p2b => wrk_ptr(jk,:,:,jt)
            cc_div(jt) = corrsp(p2a,p2b)
          END DO
          WRITE(ndunit,'(a,i3,20(1x,f7.3))') 'DIV_L',jk,cc_div(:)
        END DO
      ENDIF
      
      IF (ASSOCIATED(wrk_ptr)) DEALLOCATE(wrk_ptr)
      IF (ASSOCIATED(nudg_ptr)) DEALLOCATE(nudg_ptr)
    END IF
!-----------------------------------------------------------------------------
    IF (lppvo) THEN     ! vorticity equation
      
      if (p_parallel_io) then
        ALLOCATE(wrk_ptr(nlev,2,nsp, NDVOR))
        ALLOCATE(nudg_ptr(nlev,2,nsp))
      else
        NULLIFY(wrk_ptr)
        NULLIFY(nudg_ptr)
      endif

      DO jt=1,NDVOR
        IF (ASSOCIATED(wrk_ptr)) then
          p3a => wrk_ptr(:,:,:,jt)
        ELSE
          NULLIFY(p3a)
        END IF
        p3b => pdvor(:,:,:,jt)
        call gather_sp(p3a, p3b, gdc) 
      END DO
 
      CALL gather_sp(nudg_ptr, svtac, gdc)

      if (p_parallel_io) then
        WRITE(ndunit,'(a,20(1x,a7))') 'VOR_    ',&
           'AdPrCor', 'VerAdv ', 'VerDiff', 'GWDrag ', &
           'CuCall ', 'TimFilt', 'SemiImp', 'HorDiff', &
           'SumTend'
        cc_vor(:)=0._dp

        DO jt=1,NDVOR
          cc_vor(jt) = corrsp(nudg_ptr,wrk_ptr(:,:,:,jt))
        END DO
        WRITE(ndunit,'(a,20(1x,f7.3))') 'VOR_ALL ',cc_vor(:)

        DO jk=1,nlev
          DO jt=1,NDVOR
            p2a => nudg_ptr(jk,:,:)
            p2b => wrk_ptr(jk,:,:,jt)
            cc_vor(jt) = corrsp(p2a,p2b)
          END DO
          WRITE(ndunit,'(a,i3,20(1x,f7.3))') 'VOR_L',jk,cc_vor(:)
        END DO
      ENDIF

      IF (ASSOCIATED(wrk_ptr)) DEALLOCATE(wrk_ptr)
      IF (ASSOCIATED(nudg_ptr)) DEALLOCATE(nudg_ptr)     

    END IF
!--------------------------------------------------------------------------------
    IF (lppt) THEN    ! temperature equation

      if (p_parallel_io) then
        ALLOCATE(wrk_ptr(nlev,2,nsp, NDTEM))
        ALLOCATE(nudg_ptr(nlev,2,nsp))
      else
        NULLIFY(wrk_ptr)
        NULLIFY(nudg_ptr)
      endif

      DO jt=1,NDTEM
        IF (ASSOCIATED(wrk_ptr)) then
          p3a => wrk_ptr(:,:,:,jt)
        ELSE
          NULLIFY(p3a)
        END IF
        p3b => pdtem(:,:,:,jt)
        call gather_sp(p3a, p3b, gdc) 
      END DO
      CALL gather_sp(nudg_ptr, sttac, gdc)

      if (p_parallel_io) then
        WRITE(ndunit,'(a,20(1x,a7))') 'TEM_    ',&
          'HorAdv ', 'VerAdv ', 'EnConv ', 'Radheat', &
          'VerDiff', 'GWDrag ', 'CuCall ', 'Convect', &
          'TimFilt', 'SemiImp', 'HorDiff', 'RadLong', &
          'RadSola', 'SumTend'
        cc_tem(:)=0._dp
        DO jt=1,NDTEM
          cc_tem(jt) = corrsp(nudg_ptr,wrk_ptr(:,:,:,jt))
        END DO
        WRITE(ndunit,'(a,20(1x,f7.3))') 'TEM_ALL ',cc_tem(:)

        DO jk=1,nlev
        cc_tem(:)=0._dp
          DO jt=1,NDTEM
            p2a => nudg_ptr(jk,:,:)
            p2b => wrk_ptr(jk,:,:,jt)
            cc_tem(jt) = corrsp(p2a,p2b)
          END DO
          WRITE(ndunit,'(a,i3,20(1x,f7.3))') 'TEM_L',jk,cc_tem(:)
        END DO
      ENDIF
    
      IF (ASSOCIATED(wrk_ptr)) DEALLOCATE(wrk_ptr)
      IF (ASSOCIATED(nudg_ptr)) DEALLOCATE(nudg_ptr)
    END IF
!-------------------------------------------------------------------------
    IF (lppp) THEN    ! pressure equation

      if (p_parallel_io) then
        ALLOCATE(wrk_ptr(1,2,nsp, NDTEM))
!        ALLOCATE(nudg_ptr(nlev,2,nsp))
        ALLOCATE(nudg_ptr(1,2,nsp))
      else
        NULLIFY(wrk_ptr)
        NULLIFY(nudg_ptr)
      endif

      DO jt=1,NDPRS
        IF (ASSOCIATED(wrk_ptr)) then
          p3a => wrk_ptr(:,:,:,jt)
        ELSE
          NULLIFY(p3a)
        END IF
        p3b => pdprs(:,:,:,jt)
        call gather_sp(p3a, p3b, gdc) 
      END DO
      CALL gather_sp(nudg_ptr, sstac, gdc)

      if (p_parallel_io) then
        WRITE(ndunit,'(a,20(1x,a7))') 'PRS_    ',&
          'ConvInt', 'TimFilt', 'SemiImp', 'SumTend'
        p2a => nudg_ptr(1,:,:)
        DO jt=1,NDPRS
          p2b => wrk_ptr(1,:,:,jt) ! mz_jb_20041214
          cc_prs(jt) = corrsp(p2a,p2b)
        END DO
        WRITE(ndunit,'(a,i3,20(1x,f7.3))') 'PRS_L',0,cc_prs(:)
      END IF

      IF (ASSOCIATED(wrk_ptr)) DEALLOCATE(wrk_ptr)
      IF (ASSOCIATED(nudg_ptr)) DEALLOCATE(nudg_ptr)
    END IF
#endif

  END SUBROUTINE Nudg_Correl

!=============================================================================

#endif

END MODULE MO_NUDGING_CORREL
