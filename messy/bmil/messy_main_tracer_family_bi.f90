#include "messy_main_ppd_bi.inc"

!***********************************************************************
MODULE messy_main_tracer_family_bi
!***********************************************************************

  ! MODULE FOR TRACER FAMILIES (MESSy-SMIL)
  !
  ! Authors: Patrick Joeckel,  MPICH, Jan 2004, June 2007
  !          Astrid Kerkweg,   MPICH, May 2004, June 2007
  !          Joachim Buchholz, MPICH, November 2004
  !
#ifdef MBM_TRACER
  USE messy_main_tracer_mem_bi, ONLY: S1TRSTR
#else
  USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
#endif
#if defined(ICON)
  USE messy_main_tracer_mem_bi, ONLY: L_GPTRSTR
#endif

  USE messy_main_tracer,        ONLY: NSETID, TRSET
  USE messy_main_tracer_family

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: main_tracer_family_initialize      ! read namelist and initialize
  PUBLIC :: main_tracer_family_new_tracer      ! request family-tracers
  PUBLIC :: main_tracer_family_init_mem        ! attributes: family + fract.
  PUBLIC :: main_tracer_family_init_cpl        ! attributes: tracers
  !
  ! FOR USE IN BML, BEFORE/AFTER ADVECTION
  PUBLIC :: main_tracer_family_beforeadv
  !
  PUBLIC :: main_tracer_family_afteradv
  !
  PUBLIC :: main_tracer_family_free_mem        ! free memory
  !
  ! SPECIFIC FOR TYPE-1 (TRANSPORT FAMILIES)
  PUBLIC :: tracfamily_1_f2t
  PUBLIC :: tracfamily_1_t2f
  !PRIVATE :: tracfamily_meta_1_t2f

  ! SPECIFIC FOR TYPE-2 (SACLING, TAGGING)
  PUBLIC :: tracfamily_2_rsc    ! re-scale
  PUBLIC :: tracfamily_2_sum    ! summation

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_family_initialize

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi, error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, j

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL tracer_family_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    CALL start_message_bi(submodstr, 'INITIALISATION', substr)

    IF (p_parallel_io) THEN
       CALL tracfamily_init(status)
    END IF
    CALL p_bcast(status, p_io)
    IF (status /= 0) CALL error_bi('tracfamily_init reported an error', substr)

    ! BROADCAST RESULTS
    CALL p_bcast(l_verbose,   p_io)
    CALL p_bcast(i_diag_pe,   p_io)
    CALL p_bcast(i_diag_jrow, p_io)
    CALL p_bcast(NTF, p_io)
    DO i=1, NMAXTFAM
       CALL p_bcast(XTF(i)%IO%set, p_io)
       CALL p_bcast(XTF(i)%IO%type, p_io)
       CALL p_bcast(XTF(i)%IO%l_rescale, p_io)
       CALL p_bcast(XTF(i)%IO%name, p_io)
       CALL p_bcast(XTF(i)%IO%subname, p_io)
       CALL p_bcast(XTF(i)%fidt, p_io)
       CALL p_bcast(XTF(i)%nt, p_io)
       DO j=1, NMAXTRAC
          CALL p_bcast(XTF(i)%idt(j), p_io)
          CALL p_bcast(XTF(i)%IO%tracer(j), p_io)
          CALL p_bcast(XTF(i)%weight(j), p_io)
       END DO
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) '----> ',NTF,' TRACER FAMILY/IES REQUESTED !'
    END IF

    CALL end_message_bi(submodstr, 'INITIALISATION', substr)

  END SUBROUTINE main_tracer_family_initialize
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_family_new_tracer

    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                      , error_bi
    USE messy_main_tracer,        ONLY: tracer_error_str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'main_tracer_family_new_tracer'
    INTEGER                      :: status

    CALL start_message_bi(submodstr, 'REQUEST TRACERS', substr)

    CALL tracfamily_newtrac(status, p_parallel_io)
    IF (status /=0) CALL error_bi(tracer_error_str(status), substr)

    CALL end_message_bi(submodstr, 'REQUEST TRACERS', substr)

  END SUBROUTINE main_tracer_family_new_tracer
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_family_init_mem

    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_tracer,        ONLY: tracer_error_str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_init_mem'
    INTEGER :: status
#ifdef ICON
    INTEGER :: jg
#endif

    ! flag=1: reset meta information of family-members
    ! flag=2: set meta information of family-members to fraction

#if defined(ECHAM5) || defined(COSMO) || defined(CESM1)
    ! set meta information of family-members to fraction
    CALL tracfamily_meta(status, 2, substr, GPTRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)

!!$    CALL tracfamily_meta(status, 2, substr, LGTRSTR, p_parallel_io)
!!$    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)
#endif
#if defined(ICON)
    DO jg=1, SIZE(L_GPTRSTR)
       ! set meta information of family-members to fraction
       CALL tracfamily_meta(status, 2, substr, L_GPTRSTR(jg), p_parallel_io)
       IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)
    END DO
#endif
#ifdef MBM_TRACER
    ! set meta information of family-members to fraction
    CALL tracfamily_meta(status, 2, substr, S1TRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi( tracer_error_str(status), substr)
#endif

    ! ### ADD MORE TRACER SETS HERE

  END SUBROUTINE main_tracer_family_init_mem
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE main_tracer_family_init_cpl

    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_tracer,        ONLY: tracer_error_str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_init_cpl'
    INTEGER :: status
#ifdef ICON
    INTEGER :: jg
#endif

    ! reset meta information of family-members

#if defined(ECHAM5) || defined(CESM1) || defined(COSMO)  || defined(MESSYDWARF)
    CALL tracfamily_meta(status, 1, substr, GPTRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

!!$    CALL tracfamily_meta(status, 1, substr, LGTRSTR, p_parallel_io)
!!$    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)
#endif

#ifdef MBM_TRACER
    CALL tracfamily_meta(status, 1, substr, S1TRSTR, p_parallel_io)
    IF (status /= 0) CALL error_bi(tracer_error_str(status), substr)
#endif

#ifdef ICON
    DO jg=1, SIZE(L_GPTRSTR)
       CALL tracfamily_meta(status, 1, substr, L_GPTRSTR(jg), p_parallel_io)
       IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)
    END DO
#endif

    ! NOW STATUS IS: TRACERS ACTIVE
    CALL tracfamily_initmode(_IY_XYZN_)

    ! ### ADD MORE TRACER SETS HERE

  END SUBROUTINE main_tracer_family_init_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
 SUBROUTINE main_tracer_family_beforeadv

   USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma, npromz
   USE messy_main_timer,           ONLY: time_step_len
   USE messy_main_mpi_bi,          ONLY: p_pe
   USE messy_main_blather_bi,      ONLY: error_bi
   USE messy_main_tracer,          ONLY: tracer_error_str

   IMPLICIT NONE

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_beforeadv'
   INTEGER :: jjrow, jp, status

   ! NOTE: CALL SUMMATION OF TYPE-2 TRACERS BEFORE TRANSFORMATION OF
   !       TPYE-1 FAMILIES (T2F), I.E. BEFORE THEY ARE
   !       POTENTIALLY CONVERTED TO FRACTIONS (AS MEMBER OF FAMILY TYPE-1);
   !       TYPE-2 FAMILY MEMBERS WITHOUT RESCALING OPTION ARE ALLOWED TO
   !       BE MEMBER OF A TYPE-1 FAMILY

#if defined(ECHAM5) || defined(COSMO) || defined(ICON)
   DO jjrow=1, ngpblks

      CALL tracfamily_2_sum(GPTRSTR, jjrow)

      IF (jjrow == ngpblks) THEN
         jp = npromz
      ELSE
         jp = nproma
      END IF

      CALL tracfamily_1_t2f(status, substr, p_pe, GPTRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

   END DO

#endif
#if defined(CESM1)
   DO jjrow=1, ngpblks

      CALL tracfamily_2_sum(GPTRSTR, jjrow)

      jp = npromz(jjrow)

      CALL tracfamily_1_t2f(status, substr, p_pe, GPTRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

   END DO

#endif
#ifdef MBM_TRACER
   DO jjrow=1, ngpblks

      CALL tracfamily_2_sum(S1TRSTR, jjrow)

      IF (jjrow == ngpblks) THEN
         jp = npromz
      ELSE
         jp = nproma
      END IF

      CALL tracfamily_1_t2f(status, substr, p_pe, S1TRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

      ! ### ADD MORE TRACER SETS HERE ...

   END DO

   ! ### .. OR HERE

#endif

 END SUBROUTINE main_tracer_family_beforeadv
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
 SUBROUTINE main_tracer_family_afteradv

   USE messy_main_grid_def_mem_bi, ONLY: ngpblks, nproma, npromz
   USE messy_main_timer,           ONLY: time_step_len
   USE messy_main_mpi_bi,          ONLY: p_pe
   USE messy_main_blather_bi,      ONLY: error_bi
   USE messy_main_tracer,          ONLY: tracer_error_str

   IMPLICIT NONE

   ! LOCAL
   CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_family_afteradv'
   INTEGER :: jjrow, jp, status

   ! NOTE: CALL RE-SCALING OF TYPE-2 TRACERS AFTER TRANSFORMATION OF
   !       TPYE-1 FAMILIES (F2T), I.E. AFTER THEY ARE
   !       POTENTIALLY CONVERTED BACK FROM FRACTIONS TO TRACERS
   !       (AS MEMBER OF FAMILY TYPE-1);
   !       TYPE-2 FAMILY MEMBERS WITHOUT RESCALING OPTION ARE ALLOWED TO
   !       BE MEMBER OF A TYPE-1 FAMILY

#if defined(ECHAM5) || defined(COSMO) || defined(ICON)
   DO jjrow=1, ngpblks
      IF (jjrow == ngpblks) THEN
         jp = npromz
      ELSE
         jp = nproma
      END IF

      CALL tracfamily_1_f2t(status, substr, p_pe, GPTRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

      CALL tracfamily_2_rsc(GPTRSTR, time_step_len, jjrow)
   END DO

#endif
#if defined(CESM1)
   DO jjrow=1, ngpblks
      jp = npromz(jjrow)

      CALL tracfamily_1_f2t(status, substr, p_pe, GPTRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

      CALL tracfamily_2_rsc(GPTRSTR, time_step_len, jjrow)
   END DO

#endif
#ifdef MBM_TRACER
   DO jjrow=1, ngpblks
      IF (jjrow == ngpblks) THEN
         jp = npromz
      ELSE
         jp = nproma
      END IF

      CALL tracfamily_1_f2t(status, substr, p_pe, S1TRSTR, &
           time_step_len, jjrow, jp)
      IF (status /= 0) CALL error_bi(tracer_error_str(status),substr)

      CALL tracfamily_2_rsc(S1TRSTR, time_step_len, jjrow)
   END DO
#endif

   ! ### ADD MORE TRACER SETS HERE

 END SUBROUTINE main_tracer_family_afteradv
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
 SUBROUTINE main_tracer_family_free_mem

   IMPLICIT NONE

   CALL tracfamily_freemem

 END SUBROUTINE main_tracer_family_free_mem
! ----------------------------------------------------------------------

!***********************************************************************
!***********************************************************************
! SUBROUTINES (moved from SMCL because of rank order dependence)
!***********************************************************************
!***********************************************************************

 ! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_1_f2t(status, callstr, p_pe, setname, &
       ztmst, jjrow, ksize)

    ! CONVERT FAMILIES TO TRACERS

    USE messy_main_tracer,        ONLY: t_trinfo_tp, STRLEN_TRSET
    USE messy_main_tools,         ONLY: split_name_domain

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, MAXVAL, MINVAL, PRESENT, SIZE, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: callstr
    INTEGER,          INTENT(IN)            :: p_pe
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    REAL(DP),         INTENT(IN)            :: ztmst
    INTEGER,          INTENT(IN)            :: jjrow
    INTEGER,          INTENT(IN), OPTIONAL  :: ksize

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER                :: substr = 'tracfamily_1_f2t'
    TYPE(t_trinfo_tp),   DIMENSION(:), POINTER :: ti     => NULL()
    REAL(DP), DIMENSION(:,:,:,:),      POINTER :: zxt    => NULL()
    REAL(DP), DIMENSION(:,:,:,:),      POINTER :: zxtte  => NULL()
    REAL(DP), DIMENSION(:,:,:,:),      POINTER :: zxtm1  => NULL()
    INTEGER :: i, j
    INTEGER :: n
    INTEGER :: kproma
    INTEGER :: idt, idt2
    CHARACTER(LEN=STRLEN_TRSET) :: zname
    INTEGER                     :: domain

    status = 0

    set_loop: DO n=1, NSETID

       ! ONLY THIS SET ...
       IF (TRIM(TRSET(n)%name) /= TRIM(setname)) CYCLE

       ! NO EMPTY SETS ...
       IF (TRSET(n)%ntrac == 0) CYCLE

       ti => TRSET(n)%ti

       ! NO INCOMPLETE SETS W.R.T. TIMEFILTER
       ! (SKIP SILENTLY)
       IF (ASSOCIATED(TRSET(n)%xt)) THEN
          zxt => TRSET(n)%xt(:,:,:,:,1)
       ELSE
          CYCLE
       END IF

       IF (ASSOCIATED(TRSET(n)%xtte)) THEN
          zxtte => TRSET(n)%xtte(:,:,:,:,1)
       ELSE
          CYCLE
       END IF

       IF (ASSOCIATED(TRSET(n)%xtm1)) THEN
          zxtm1 => TRSET(n)%xtm1(:,:,:,:,1)
       ELSE
          CYCLE
       END IF

       IF (l_verbose .AND.(jjrow==i_diag_jrow).AND.(p_pe==i_diag_pe)) THEN
          WRITE(*,*) '======================================================'
          WRITE(*,*) substr,' (p_pe=',p_pe,'; jrow = ',jjrow,')'
          WRITE(*,*) ' ... CALLED FROM ',TRIM(callstr), &
               ' FOR SET ',TRIM(setname)
       END IF

       ! CHECK STATUS
       IF (IS_TRACER(n)%ptr(jjrow)) THEN
          status = 2000 ! TRACER MODE IS ALREADY ACTIVE
          RETURN
       END IF

       IF (PRESENT(ksize)) THEN
          kproma = ksize
       ELSE
          kproma = SIZE(zxt,1)
       END IF

       family_loop: DO i=1, NTF

          ! ONLY TRACER-FAMILIES FROM THIS SET ...
          CALL split_name_domain(status, TRIM(setname) &
               , zname, domain)
          status = 0
          IF ( (TRIM(XTF(i)%IO%set) /= TRIM(setname))  .AND. &
               (TRIM(XTF(i)%IO%set) /= TRIM(zname)) ) CYCLE

          ! ONLY FAMILIES ...
          IF (XTF(i)%fidt == 0) CYCLE

          ! CONVERSION ONLY FOR TRANSPORT-FAMILIES ...
          IF (XTF(i)%IO%type /= FTYPE_TRAN) CYCLE


          IF (l_verbose .AND.(jjrow==i_diag_jrow).AND.(p_pe==i_diag_pe)) THEN
             !
             ! NOTE: THE META-INFORMATION OF THE TRACERS IS RE-SET
             !       jjrow; THUS FOR jrow > 1 THE UNIT IS
             !       (see ti(XTF(i)%idt(j))%tp = XTF(i)%ti(j)) FOR EVERY
             !       'mol/mol' IN THE DIAGNOSTIC OUTPUT
             !
             WRITE(*,*) &
                  '------------------------------------------------------'
             WRITE(*,*) '### TRACERS AND TENDENCIES BEFORE f2t: '//&
                  &'(p_pe=',p_pe,'; jrow = ',jjrow,')'
             ! FAMILY
             WRITE(*,*) &
                  ti(XTF(i)%fidt)%tp%ident%fullname,        &
                  ' ('//ti(XTF(i)%fidt)%tp%ident%unit//')'

             idt2 = XTF(i)%fidt
             WRITE(*,'(6(e12.4,1x))') &
                  MINVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MINVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MINVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt2) ))

             ! TRACERS
             DO j=1, XTF(i)%nt
                IF (XTF(i)%idt(j) == 0) CYCLE
                WRITE(*,*) &
                     ti(XTF(i)%idt(j))%tp%ident%fullname,         &
                     ' ('//ti(XTF(i)%idt(j))%tp%ident%unit//')'
                idt = XTF(i)%idt(j)
                WRITE(*,'(6(e12.4,1x))') &
                     MINVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt) )),   &
                     MAXVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt) )),   &
                     MINVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt) )), &
                     MAXVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt) )), &
                     MINVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt) )), &
                     MAXVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt) ))
             END DO
          END IF

          tracer_loop: DO j=1, XTF(i)%nt

             IF (XTF(i)%idt(j) == 0) CYCLE
             idt  = XTF(i)%idt(j)
             idt2 = XTF(i)%fidt
             zxt( _RI_XYZN_(:,jjrow,:,idt) ) =    &
                  zxt( _RI_XYZN_(:,jjrow,:,idt) ) &
                  * zxt(_RI_XYZN_(:,jjrow,:,idt2) ) / XTF(i)%weight(j)

             zxtm1( _RI_XYZN_(:,jjrow,:,idt) ) =    &
                  zxtm1( _RI_XYZN_(:,jjrow,:,idt) ) &
                  * zxtm1( _RI_XYZN_(:,jjrow,:,idt2) ) / XTF(i)%weight(j)

             zxtte( _RI_XYZN_(:,jjrow,:,idt) ) =        &
                  ( zxtte( _RI_XYZN_(:,jjrow,:,idt) ) * &
                  (zxtm1( _RI_XYZN_(:,jjrow,:,idt2) )   &
                  + zxtte( _RI_XYZN_(:,jjrow,:,idt2) )  &
                  * ztmst)/XTF(i)%weight(j) &
                  - zxtm1( _RI_XYZN_(:,jjrow,:,idt) ) ) / ztmst

             ! RESET SAVED TRACER INFORMATION
             ti(XTF(i)%idt(j))%tp = XTF(i)%ti(j)

          END DO tracer_loop

          IF (l_verbose .AND.(jjrow==i_diag_jrow).AND.(p_pe==i_diag_pe)) THEN
             WRITE(*,*) '### TRACERS AND TENDENCIES AFTER f2t: '//&
                  &'(p_pe=',p_pe,'; jrow = ',jjrow,')'
             ! FAMILY
             WRITE(*,*) &
                  ti(XTF(i)%fidt)%tp%ident%fullname,        &
                  ' ('//ti(XTF(i)%fidt)%tp%ident%unit//')'
             idt2 = XTF(i)%fidt
             WRITE(*,'(6(e12.4,1x))') &
                  MINVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MINVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MINVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt2) ))

             ! TRACERS
             DO j=1, XTF(i)%nt
                IF (XTF(i)%idt(j) == 0) CYCLE
                WRITE(*,*) &
                     ti(XTF(i)%idt(j))%tp%ident%fullname,         &
                     ' ('//ti(XTF(i)%idt(j))%tp%ident%unit//')'
                idt  = XTF(i)%idt(j)
                WRITE(*,'(6(e12.4,1x))') &
                     MINVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt) )),   &
                     MAXVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt) )),   &
                     MINVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt) )), &
                     MAXVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt) )), &
                     MINVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt) )), &
                     MAXVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt) ))
             END DO
             WRITE(*,*) '-----------------------------------------------------'
          END IF

       END DO family_loop

       IF (l_verbose .AND.(jjrow==i_diag_jrow).AND.(p_pe==i_diag_pe)) &
            WRITE(*,*) '======================================================'

       ! TRACERS ACTIVE
       IS_TRACER(n)%ptr(jjrow) = .TRUE.

    END DO set_loop

  END SUBROUTINE tracfamily_1_f2t
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_1_t2f(status, callstr, p_pe, setname, &
       ztmst, jjrow, ksize, l_frac)

    ! CONVERT TRACERS TO FAMILIES

    USE messy_main_tracer,        ONLY: t_trinfo_tp, STRLEN_TRSET
    USE messy_main_constants_mem, ONLY: TINY_DP
    USE messy_main_tools,         ONLY: split_name_domain

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, MAXVAL, MINVAL, PRESENT, SIZE, TRIM

    ! I/O
    INTEGER,          INTENT(OUT)           :: status
    CHARACTER(LEN=*), INTENT(IN)            :: callstr
    INTEGER,          INTENT(IN)            :: p_pe
    CHARACTER(LEN=*), INTENT(IN)            :: setname
    REAL(DP),         INTENT(IN)            :: ztmst
    INTEGER,          INTENT(IN)            :: jjrow
    INTEGER,          INTENT(IN), OPTIONAL  :: ksize
    LOGICAL,          INTENT(IN), OPTIONAL  :: l_frac

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER                :: substr = 'tracfamily_1_t2f'
    TYPE(t_trinfo_tp),   DIMENSION(:), POINTER :: ti     => NULL()
    REAL(DP), DIMENSION(:,:,:,:),      POINTER :: zxt    => NULL()
    REAL(DP), DIMENSION(:,:,:,:),      POINTER :: zxtte  => NULL()
    REAL(DP), DIMENSION(:,:,:,:),      POINTER :: zxtm1  => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: zptr_xm1 => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: zptr_xte => NULL()
    REAL(DP)                            :: zf
    INTEGER :: i, j, jj1, jj2
    INTEGER :: n
    INTEGER :: kproma
    LOGICAL :: zl_frac
    INTEGER :: idt, idt2
    INTEGER :: idx
    CHARACTER(LEN=STRLEN_TRSET) :: zname
    INTEGER                     :: domain

    status = 0
    IF (PRESENT(l_frac)) THEN
       ! CALCULATE FAMILIES ONLY, IF FALSE
       zl_frac = l_frac
    ELSE
       ! CALCULATE FAMILIES AND TRACER FRACTIONS
       zl_frac = .TRUE. ! default
    END IF

    set_loop: DO n=1, NSETID

       ! ONLY THIS SET ...
       IF (TRIM(TRSET(n)%name) /= TRIM(setname)) CYCLE

       ! NO EMPTY SETS ...
       IF (TRSET(n)%ntrac == 0) CYCLE

       ti => TRSET(n)%ti

       ! NO INCOMPLETE SETS W.R.T. TIMEFILTER
       ! (SKIP SILENTLY)
       IF (ASSOCIATED(TRSET(n)%xt)) THEN
          zxt => TRSET(n)%xt(:,:,:,:,1)
       ELSE
          CYCLE
       END IF

       IF (ASSOCIATED(TRSET(n)%xtte)) THEN
          zxtte => TRSET(n)%xtte(:,:,:,:,1)
       ELSE
          CYCLE
       END IF

       IF (ASSOCIATED(TRSET(n)%xtm1)) THEN
          zxtm1 => TRSET(n)%xtm1(:,:,:,:,1)
       ELSE
          CYCLE
       END IF

       IF (l_verbose .AND.(jjrow==i_diag_jrow).AND.(p_pe==i_diag_pe)) THEN
          WRITE(*,*) '======================================================'
          WRITE(*,*) substr,' (p_pe=',p_pe,'; jrow = ',jjrow,')'
          WRITE(*,*) ' ... CALLED FROM ',TRIM(callstr), &
               ' FOR SET ',TRIM(setname)
       END IF

       ! CHECK STATUS
       IF (.NOT.IS_TRACER(n)%ptr(jjrow)) THEN
          status = 2001    ! FAMILY MODE IS ALREADY ACTIVE
          RETURN
       END IF

       IF (PRESENT(ksize)) THEN
          kproma = ksize
       ELSE
          kproma = SIZE(zxt,1)
       END IF

       family_loop: DO i = 1, NTF

          ! ONLY TRACER-FAMILIES FROM THIS SET ...
          CALL split_name_domain(status, TRIM(setname) &
               , zname, domain)
          status = 0
          IF ( (TRIM(XTF(i)%IO%set) /= TRIM(setname)) .AND. &
               (TRIM(XTF(i)%IO%set) /= TRIM(zname)) ) CYCLE

          ! ONLY FAMILIES ...
          IF (XTF(i)%fidt == 0) CYCLE

          ! CONVERSION ONLY FOR TRANSPORT-FAMILIES
          IF (XTF(i)%IO%type /= FTYPE_TRAN) CYCLE

          IF (l_verbose .AND.(jjrow==i_diag_jrow).AND.(p_pe==i_diag_pe)) THEN
             !
             ! NOTE: THE META-INFORMATION OF THE TRACERS IS RE-SET
             !       (see CALL tracfamily_meta_1_t2f below) FOR EVERY
             !       jjrow; THUS FOR jrow > 1 THE UNIT IS
             !       'frac. of ...' IN THE DIAGNOSTIC OUTPUT
             !
             WRITE(*,*) '-------------------------------------------------'
             WRITE(*,*) '### TRACERS AND TENDENCIES BEFORE t2f: '//&
                  &'(p_pe=',p_pe,'; jrow = ',jjrow,')'
             ! FAMILY
             WRITE(*,*) &
                  ti(XTF(i)%fidt)%tp%ident%fullname,        &
                  ' ('//ti(XTF(i)%fidt)%tp%ident%unit//')'
             idt2 = XTF(i)%fidt
             WRITE(*,'(6(e12.4,1x))') &
                  MINVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MINVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MINVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt2) ))

             ! TRACERS
             DO j=1, XTF(i)%nt
                IF (XTF(i)%idt(j) == 0) CYCLE
                WRITE(*,*) &
                     ti(XTF(i)%idt(j))%tp%ident%fullname,        &
                     ' ('//ti(XTF(i)%idt(j))%tp%ident%unit//')'
                idt = XTF(i)%idt(j)
                WRITE(*,'(6(e12.4,1x))') &
                     MINVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MAXVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MINVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MAXVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MINVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MAXVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt) ))
             END DO
          END IF

          ! RE-INITIALIZE FAMILY TRACER WITH ZERO

          idt2 = XTF(i)%fidt
          zxt  (_RI_XYZN_(:,jjrow,:,idt2) ) = 0.0_DP
          zxtm1(_RI_XYZN_(:,jjrow,:,idt2) ) = 0.0_DP
          zxtte(_RI_XYZN_(:,jjrow,:,idt2) ) = 0.0_DP

          ! FIRST LOOP: SUMMATION (FAMILY)
          tracer_loop1: DO j=1, XTF(i)%nt
             IF (XTF(i)%idt(j) == 0) CYCLE

             idt = XTF(i)%idt(j)
             idt2 = XTF(i)%fidt

             zxt( _RI_XYZN_(:,jjrow,:,idt2) ) =    &
                  zxt( _RI_XYZN_(:,jjrow,:,idt2) ) &
                  + zxt(_RI_XYZN_(:,jjrow,:,idt) ) * XTF(i)%weight(j)

             zxtm1( _RI_XYZN_(:,jjrow,:,idt2) ) =    &
                  zxtm1( _RI_XYZN_(:,jjrow,:,idt2) ) &
                  + zxtm1( _RI_XYZN_(:,jjrow,:,idt) ) * XTF(i)%weight(j)

             zxtte( _RI_XYZN_(:,jjrow,:,idt2) ) =    &
                  zxtte( _RI_XYZN_(:,jjrow,:,idt2) ) &
                  + zxtte( _RI_XYZN_(:,jjrow,:,idt) ) * XTF(i)%weight(j)

          END DO tracer_loop1

          fractions: IF (zl_frac) THEN

             zptr_xm1 => zxtm1( _RI_XYZN_(:,jjrow,:,:) )
             zptr_xte => zxtte( _RI_XYZN_(:,jjrow,:,:) )

             ! SECOND LOOP: FRACTION (TRACERS AND TENDENCIES)
             tracer_loop2: DO j=1, XTF(i)%nt
                IF (XTF(i)%idt(j) == 0) CYCLE

                idt  = XTF(i)%idt(j)
                idt2 = XTF(i)%fidt
                !
                ! XT
                DO jj2=1, SIZE(zxt,_IZ_XYZN_)
                   idx = jj2
                   DO jj1=1, SIZE(zxt,_IX_XYZN_)
                      IF (ABS(zxt( _RI_XYZN_(jj1,jjrow,idx,idt2) )) &
                           > TINY_DP) THEN
                         zxt( _RI_XYZN_(jj1,jjrow,idx,idt) ) =    &
                              zxt( _RI_XYZN_(jj1,jjrow,idx,idt) ) &
                              *XTF(i)%weight(j)                      &
                              / zxt( _RI_XYZN_(jj1,jjrow,idx,idt2) )
                      ELSE
                         zxt( _RI_XYZN_(jj1,jjrow,idx,idt) ) = 0.0_DP
                      END IF
                   END DO
                END DO

                ! XTTE
                DO jj2=1, SIZE(zptr_xte,_IZ_X_ZN_)
                   idx = jj2
                   DO jj1=1, SIZE(zptr_xte,_IX_X_ZN_)
                      zf = zxtm1( _RI_XYZN_(jj1,jjrow,idx,idt2) )   &
                           + zxtte( _RI_XYZN_(jj1,jjrow,idx,idt2) ) &
                           * ztmst
                      IF (ABS(zf) > TINY_DP) THEN
                         zptr_xte( _RI_X_ZN_(jj1,idx,idt) ) =          &
                              ( zxtm1( _RI_XYZN_(jj1,jjrow,idx,idt) )  &
                              + zxtte( _RI_XYZN_(jj1,jjrow,idx,idt) )  &
                              * ztmst ) &
                              *XTF(i)%weight(j)                           &
                              / zf
                      ELSE
                         zptr_xte( _RI_X_ZN_(jj1,idx,idt) ) = 0.0_DP
                      END IF
                   END DO
                END DO

                ! XTM1
                DO jj1=1, SIZE(zptr_xm1,_IX_X_ZN_)
                   DO jj2=1, SIZE(zptr_xm1,_IZ_X_ZN_)
                      idx = jj2
                      IF (ABS(zxtm1( _RI_XYZN_(jj1,jjrow,idx,idt2) )) &
                           > TINY_DP) THEN
                         zptr_xm1( _RI_X_ZN_(jj1,idx,idt) ) =       &
                              zxtm1( _RI_XYZN_(jj1,jjrow,idx,idt) ) &
                              *XTF(i)%weight(j)                  &
                              / zxtm1( _RI_XYZN_(jj1,jjrow,idx,idt2) )
                      ELSE
                         zptr_xm1( _RI_X_ZN_(jj1,idx,idt) ) = 0.0_DP
                      END IF
                   END DO
                END DO

                ! RESET TRACER INFO (AFTER LAST ROW)
                CALL tracfamily_meta_1_t2f(ti(XTF(i)%idt(j))%tp, &
                     XTF(i)%IO%name, XTF(i)%IO%subname )

             END DO tracer_loop2

          END IF fractions

          IF (l_verbose .AND.(jjrow==i_diag_jrow).AND.(p_pe==i_diag_pe)) THEN
             WRITE(*,*) '### TRACERS AND TENDENCIES AFTER t2f: '//&
                  &'(p_pe=',p_pe,'; jrow = ',jjrow,')'
             ! FAMILY
             WRITE(*,*) &
                  ti(XTF(i)%fidt)%tp%ident%fullname,        &
                  ' ('//ti(XTF(i)%fidt)%tp%ident%unit//')'
             idt2 = XTF(i)%fidt
             WRITE(*,'(6(e12.4,1x))') &
                  MINVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MINVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MINVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt2) )),  &
                  MAXVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt2) ))

             ! TRACERS
             DO j=1, XTF(i)%nt
                IF (XTF(i)%idt(j) == 0) CYCLE
                WRITE(*,*) &
                     ti(XTF(i)%idt(j))%tp%ident%fullname,        &
                     ' ('//ti(XTF(i)%idt(j))%tp%ident%unit//')'
                idt = XTF(i)%idt(j)
                WRITE(*,'(6(e12.4,1x))') &
                     MINVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MAXVAL(zxt( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MINVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MAXVAL(zxtm1( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MINVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt) )),  &
                     MAXVAL(zxtte( _RI_XYZN_(1:kproma,jjrow,:,idt) ))
             END DO
             WRITE(*,*) '-------------------------------------------------'
          END IF

       END DO family_loop

       IF (l_verbose .AND.(jjrow==i_diag_jrow).AND.(p_pe==i_diag_pe)) &
            WRITE(*,*) '======================================================'

       IF (zl_frac) THEN
          ! THE STATUS IS NOW: FAMILIES ACTIVE
          IS_TRACER(n)%ptr(jjrow) = .FALSE.
       END IF

    END DO set_loop

  END SUBROUTINE tracfamily_1_t2f
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_2_rsc(setname, ztmst, jjrow)

    USE messy_main_tracer,        ONLY: STRLEN_TRSET
    USE messy_main_tools,         ONLY: split_name_domain

    IMPLICIT NONE

    INTRINSIC :: ABS, ASSOCIATED, SIZE, TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: setname
    REAL(DP),         INTENT(IN) :: ztmst
    INTEGER,          INTENT(IN) :: jjrow

    ! LOCAL
    REAL(dp), PARAMETER                   :: SMALL = 1.e-20_dp
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zxt    => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zxtte  => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zxtm1  => NULL()
    INTEGER                               :: i, j, n, n1, n2
    REAL(DP), DIMENSION(:,:),     POINTER :: a => NULL()
    LOGICAL                               :: lcomplete
    INTEGER                               :: idt, idt2
    CHARACTER(LEN=STRLEN_TRSET)           :: zname
    INTEGER                               :: domain
    INTEGER                               :: status

    set_loop: DO n=1, NSETID

       ! ONLY THIS SET ...
       IF (TRIM(TRSET(n)%name) /= TRIM(setname)) CYCLE

       ! NO EMPTY SETS ...
       IF (TRSET(n)%ntrac == 0) CYCLE

       ! NO INCOMPLETE SETS W.R.T. TIMEFILTER
       ! (SKIP SILENTLY)
       IF (ASSOCIATED(TRSET(n)%xt)) THEN
          zxt => TRSET(n)%xt(:,:,:,:,1)
       ELSE
          CYCLE
       END IF

       IF (ASSOCIATED(TRSET(n)%xtte)) THEN
          zxtte => TRSET(n)%xtte(:,:,:,:,1)
       END IF

       IF (ASSOCIATED(TRSET(n)%xtm1)) THEN
          zxtm1 => TRSET(n)%xtm1(:,:,:,:,1)
       END IF

       lcomplete = ASSOCIATED(zxtte) .AND. ASSOCIATED(zxtm1)

       n1 = SIZE(zxt,1)
       n2 = SIZE(zxt,2)

       ALLOCATE(a(n1, n2))

       family_loop: DO i = 1, NTF

          CALL split_name_domain(status, TRIM(TRSET(n)%name) &
               , zname, domain)
          IF ( (TRIM(XTF(i)%IO%set) /= TRIM(TRSET(n)%name)) .AND. &
               (TRIM(XTF(i)%IO%set) /= TRIM(zname)) ) CYCLE

          IF (XTF(i)%fidt == 0) CYCLE

          ! RE-SCALING ONLY FOR SCALING-FAMILIES
          IF (XTF(i)%IO%type /= FTYPE_SCAL) CYCLE

          ! RE-SCALING SWITCHED OFF
          IF (.NOT. XTF(i)%IO%l_rescale) CYCLE

          complete: IF (lcomplete) THEN

             ! CALCULATE SUM OF TRACERS
             a(:,:) = 0.0_DP
             tracer_loop1: DO j=1, XTF(i)%nt
                IF (XTF(i)%idt(j) == 0) CYCLE
                idt = XTF(i)%idt(j)
                a(:,:) = a(:,:) + (zxtm1(_RI_XYZN_(:,jjrow,:,idt) ) &
                     + zxtte(_RI_XYZN_(:,jjrow,:,idt) ) * ztmst) &
                     * XTF(i)%weight(j)
             END DO tracer_loop1

             idt2 = XTF(i)%fidt
             WHERE(ABS(a(:,:)) >= SMALL)
                a(:,:) = &
                     (zxtm1(_RI_XYZN_(:,jjrow,:,idt2) ) +        &
                     zxtte(_RI_XYZN_(:,jjrow,:,idt2) ) * ztmst) &
                     / a(:,:)
             ELSEWHERE
                a(:,:) = 1.0_DP
             ENDWHERE

             tracer_loop2: DO j=1, XTF(i)%nt
                idt = XTF(i)%idt(j)
                zxtte(_RI_XYZN_(:,jjrow,:,idt) ) = &
                     zxtte(_RI_XYZN_(:,jjrow,:,idt) ) * a(:,:) +  &
                     zxtm1(_RI_XYZN_(:,jjrow,:,idt) ) * (a(:,:) - &
                     1.0_DP) &
                     / ztmst
             END DO tracer_loop2

          ELSE

             a(:,:) = 0.0_DP

             tracer_loop3: DO j=1, XTF(i)%nt
                IF (XTF(i)%idt(j) == 0) CYCLE
                idt = XTF(i)%idt(j)
                a(:,:) = a(:,:) + (zxt(_RI_XYZN_(:,jjrow,:,idt) ) &
                     * XTF(i)%weight(j))
             END DO tracer_loop3

             idt2 = XTF(i)%fidt
             WHERE(ABS(a(:,:)) >= SMALL)
                a(:,:) = zxt(_RI_XYZN_(:,jjrow,:,idt2) ) / a(:,:)
             ELSEWHERE
                a(:,:) = 1.0_DP
             ENDWHERE

             tracer_loop4: DO j=1, XTF(i)%nt
                idt = XTF(i)%idt(j)
                zxt(_RI_XYZN_(:,jjrow,:,idt) ) = &
                     zxt(_RI_XYZN_(:,jjrow,:,idt) ) * a(:,:)
             END DO tracer_loop4

          END IF complete

       END DO family_loop

       DEALLOCATE(a)

    END DO set_loop

  END SUBROUTINE tracfamily_2_rsc
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE tracfamily_2_sum(setname, jjrow)

    USE messy_main_tracer,        ONLY: STRLEN_TRSET
    USE messy_main_tools,         ONLY: split_name_domain

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: setname
    INTEGER,          INTENT(IN) :: jjrow

    ! LOCAL
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zxt    => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zxtte  => NULL()
    REAL(DP), DIMENSION(:,:,:,:), POINTER :: zxtm1  => NULL()
    INTEGER                               :: i, j, n
    LOGICAL                               :: lcomplete
    INTEGER                               :: idt, idt2
    CHARACTER(LEN=STRLEN_TRSET)           :: zname
    INTEGER                               :: domain
    INTEGER                               :: status

    set_loop: DO n=1, NSETID

       ! ONLY THIS SET ...
       IF (TRIM(TRSET(n)%name) /= TRIM(setname)) CYCLE

       ! NO EMPTY SETS ...
       IF (TRSET(n)%ntrac == 0) CYCLE

       ! NO INCOMPLETE SETS W.R.T. TIMEFILTER
       ! (SKIP SILENTLY)
       IF (ASSOCIATED(TRSET(n)%xt)) THEN
          zxt => TRSET(n)%xt(:,:,:,:,1)
       ELSE
          CYCLE
       END IF

       IF (ASSOCIATED(TRSET(n)%xtte)) THEN
          zxtte => TRSET(n)%xtte(:,:,:,:,1)
       END IF

       IF (ASSOCIATED(TRSET(n)%xtm1)) THEN
          zxtm1 => TRSET(n)%xtm1(:,:,:,:,1)
       END IF

       lcomplete = ASSOCIATED(zxtte) .AND. ASSOCIATED(zxtm1)

       family_loop: DO i = 1, NTF

          CALL split_name_domain(status, TRIM(TRSET(n)%name) &
               , zname, domain)
          IF ( (TRIM(XTF(i)%IO%set) /= TRIM(TRSET(n)%name)) .AND. &
               (TRIM(XTF(i)%IO%set) /= TRIM(zname) ) ) CYCLE

          IF (XTF(i)%fidt == 0) CYCLE

          ! SUMMATION ONLY FOR SCALING-FAMILIES
          IF (XTF(i)%IO%type /= FTYPE_SCAL) CYCLE

          ! RE-INITIALIZE FAMILY TRACER WITH ZERO
          idt2 = XTF(i)%fidt
          zxt  (_RI_XYZN_(:,jjrow,:,idt2) ) = 0.0_DP
          IF (lcomplete) zxtm1(_RI_XYZN_(:,jjrow,:,idt2) ) = 0.0_DP
          IF (lcomplete) zxtte(_RI_XYZN_(:,jjrow,:,idt2) ) = 0.0_DP

          ! SUMMATION (FAMILY)
          tracer_loop_3: DO j=1, XTF(i)%nt
             IF (XTF(i)%idt(j) == 0) CYCLE
             idt  = XTF(i)%idt(j)
             idt2 = XTF(i)%fidt

             zxt(_RI_XYZN_(:,jjrow,:,idt2) ) = &
                  zxt  (_RI_XYZN_(:,jjrow,:,idt2) ) &
                  + zxt  (_RI_XYZN_(:,jjrow,:,idt) ) * XTF(i)%weight(j)

             IF (.NOT. lcomplete) CYCLE

             zxtm1(_RI_XYZN_(:,jjrow,:,idt2) ) = &
                  zxtm1(_RI_XYZN_(:,jjrow,:,idt2) ) &
                  + zxtm1(_RI_XYZN_(:,jjrow,:,idt) ) * XTF(i)%weight(j)

             zxtte(_RI_XYZN_(:,jjrow,:,idt2) ) = &
                  zxtte(_RI_XYZN_(:,jjrow,:,idt2) ) &
                  + zxtte(_RI_XYZN_(:,jjrow,:,idt) ) * XTF(i)%weight(j)
          END DO tracer_loop_3

       END DO family_loop

    END DO set_loop

  END SUBROUTINE tracfamily_2_sum
! ----------------------------------------------------------------------

!***********************************************************************
END MODULE messy_main_tracer_family_bi
!***********************************************************************
