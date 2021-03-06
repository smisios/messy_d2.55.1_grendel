! **********************************************************************
!
! SUBMODEL: DIUMOD
!
! THIS SUBMODEL IS USED TO APPROXIMATE THE DIURNAL ABS. SIN-VARIATION 
! OF A TRACER FROM AN AVERAGE DAILY VALUE
!
! Author : Sergey Gromov, MPI-C, November 2014
!
! References:
!
! * SO FAR NONE
!
! **********************************************************************

!> \mainpage Description
!!
!!   DIUMOD (<i>diu</i>rnal <i>mod</i>ulation) sub-model is made to simulate solar
!!   zenith angle (SZA) correlated variations in tracer abundance around a given average
!!   value.
!!
!!   The initially intended application of this sub-model is to better represent (mimic) the
!!   diurnal variations in OH (or any short-lived photolysis-correlated tracer) whose
!!   average daily (or montly) value is prescribed. The latter is conventionally applied
!!   in global models as a constant throughout a day.
!!
!!   The calculations involved here approximate the solar declination to derive the (cosine of) SZA and
!!   use local solar time to gauge the modulation amplitude around a given baseline. The formulation follows.
!!
!!   The modulation factor \f$\mu(\eta)\f$ is derived using modulation function \f$f(\eta)\f$ 
!!   (using the local solar day fraction \f$\eta\f$ of 0 at noon, \f$\pm0.5\f$ at night):
!!
!!   \f$\\
!!   f(\eta) = c + a\cos( {2\pi\eta}),\\
!!   c = {{( {ma + mi} )} \mathord{/
!!    {\vphantom {{( {ma + mi} )} 2}}\\
!!    \kern-\nulldelimiterspace} 2}\\
!!   a = {{( {ma - mi} )} \mathord{/
!!    {\vphantom {{( {ma - mi} )} 2}}\\
!!    \kern-\nulldelimiterspace} 2}\\
!!   \f$
!!
!!   where \f$a\f$ is amplitude sought to fit \f$f(\eta)\f$ about the average value \f$c\f$ (or \f${f(\eta)}\f$),
!!   with \f$mi\f$ and \f$ma\f$ being the minimum and maximum cosine of SZA at a givel latitude, respectively.
!!
!!   Except for the cases of polar day/night (i.e., when SZA does not traverse zero), the solution to the
!!   above equation gives the times of sun-set/rise is:
!!
!!   \f$\\
!!   f({{\eta _r}}) = 0:\;\;\;{\eta _r} =  \pm {{\arccos ({ - {c \mathord{/
!!    {\vphantom {c a}}
!!    \kern-\nulldelimiterspace} a}})} \mathord{/
!!    {\vphantom {{\arccos ({ - {c \mathord{/
!!    {\vphantom {c a}}
!!    \kern-\nulldelimiterspace} a}})} {2\pi }}}
!!    \kern-\nulldelimiterspace} {2\pi }}\\
!!   \f$
!!
!!   and allows to derive the daily average value \f$\bar f\f$ of
!!
!!   \f$\\
!!   \bar f = \int\limits_{ - {\eta _r}}^{{\eta _r}} {f(\eta ) = 2c{\eta _r} + {{a\sin ({2\pi {\eta _r}})} \mathord{/
!!    {\vphantom {{a\sin ({2\pi {\eta _r}})} \pi }}
!!    \kern-\nulldelimiterspace} \pi }} \\
!!   \f$
!!
!!   Finally, the modulation factor is calculated by normalising the modulation function:
!!
!!   \f$
!!   \mu (\eta ) = {{f(\eta )} \mathord{/
!!    {\vphantom {{f(\eta )} {\bar f}}}
!!    \kern-\nulldelimiterspace} {\bar f}}
!!   \f$
!!
!!   Integrating \f$\mu(\eta)\f$ over a period of one day yields unity. In the case of polar day/night,
!!   \f$c\f$ is set to unity and \f$a\f$ is scaled to vary similar proportional to the cos(SZA) variations.
!!
!!   Further calculus involves derivation of 
!!   - approximate solar declination (using day of year \f${N_{day}}\f$) \n
!!   \f$\\
!!   {\delta _ \odot } = \frac{{ - 23.44^\circ }}{{180^\circ  \cdot \pi }}\cos \left( {2\pi \frac{{{N_{day}} + 10}}{{365}}} \right)\\
!!   \f$
!!   - local solar day fraction \n
!!   \f$\\
!!   \eta  = \left[ {(hou{r_{{\rm{UTC}}}} + 12) + \frac{{24}}{{360}}\lambda } \right]\bmod 24
!!   \f$
!!   - cosine of SZA \n
!!   \f$\\
!!   \cos(\theta) = \sin \varphi \sin {\delta _ \odot } + \cos \varphi \cos {\delta _ \odot }\cos \eta \\
!!   \f$
!!
!!   where \f$\phi\f$ and \f$\lambda\f$ denote local latitude and longitude, respectively.
!!
!!   DIUMOD requires an input field to modulate and the destination tracer field (replaced
!!   every time step) to output to, defined by the time of the coupling.
!!   A working configuration (MESSy v2.50) was tested using the average values from a 3D field
!!   read by IMPORT submodel with the output to a tracer created by TREXP submodel.
!!

!> \brief DIUMOD submodel interface module.

!> \authors Sergey Gromov, MPI-C, 2014-2015
!>   - original DIUMOD si code

!> \version 1.0
!>   - initial release
!>   - successful test with MESSy 2.50 (OHVAR experiment) using imported OH field and TREXP-created OH tracers
!>   - code cleanup for MESSy 2.55 release (checks with GNU/NAG compilers)

!> \todo
!>   - optional parameter for phase shift

! **********************************************************************
MODULE messy_diumod_si
  !
  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                     mtend_get_start_l,      &
                                     mtend_add_l,            &
                                     mtend_register,         &
                                     mtend_id_tracer
#endif
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  USE messy_main_constants_mem, ONLY: DP, TRACNAMELEN => STRLEN_MEDIUM
  ! MESSy
  USE messy_diumod

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  !> structural type for input data record
  TYPE T_DIUMOD_IO

     CHARACTER(LEN=STRLEN_CHANNEL) :: cha_in = ''       !< channel name of the input field to modulate
     CHARACTER(LEN=STRLEN_OBJECT)  :: obj_in = ''       !< object name of the input field to modulate
     CHARACTER(LEN=TRACNAMELEN)    :: tra_out = ''      !< name of the output (modulated) tracer
     REAL(dp)                      :: sca_fac = 1.0_dp  !< constant scaling factor
  END TYPE T_DIUMOD_IO

  TYPE T_DIUMOD
     TYPE(T_DIUMOD_IO)            :: io
     LOGICAL                      :: ok
     REAL(DP), DIMENSION(:,:,:), POINTER :: dat_in  => NULL()  !< pointer to the input field
     INTEGER                             :: idt_out = 0        !< output tracer index
  END TYPE T_DIUMOD

  INTEGER, PARAMETER              :: NMAXDIUMOD = 50                   !< maximum no. of input records (modulation jobs)
  TYPE(T_DIUMOD_IO), DIMENSION(NMAXDIUMOD), SAVE :: DIUMOD             !< CTRL namelist to read in
  TYPE(T_DIUMOD),    DIMENSION(:), POINTER, SAVE :: XDIUMOD => NULL()  !< data to work with (read from nml)
  INTEGER,                                  SAVE :: NDIUMOD            !< no. of modulation jobs

! op_pj_20160510+
#ifdef MESSYTENDENCY
  INTEGER, SAVE :: my_handle
#endif
! op_pj_20160510-

  PUBLIC :: diumod_initialize
  PUBLIC :: diumod_init_memory    ! op_pj_20160510
  PUBLIC :: diumod_init_coupling
  PUBLIC :: diumod_local_start
  !PUBLIC :: diumod_local_end
  PUBLIC :: diumod_free_memory
  PRIVATE :: diumod_read_nml_cpl

CONTAINS

! ======================================================================
! PUBLIC ROUTINES
! ======================================================================

!> \brief Initialisation
!> \details Initialises data structures, reads input (NML)

! ----------------------------------------------------------------------
  SUBROUTINE diumod_initialize

    ! DIUMOD MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,ONLY: error_bi
    ! MESSy
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'diumod_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i

    ALLOCATE(XDIUMOD(NMAXDIUMOD))

    ! initialize CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL diumod_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('', substr)
    END IF

    CALL start_message_bi(modstr, 'INITIALISATION', substr)

    IF (p_parallel_io) THEN
       NDIUMOD = 0
       DO i=1, NMAXDIUMOD

          IF (TRIM(DIUMOD(i)%obj_in) == '') CYCLE

          ! active entry
          NDIUMOD = NDIUMOD + 1
          XDIUMOD(NDIUMOD)%io = DIUMOD(i)
          WRITE(*,*) 'in -> mod -> out: ', TRIM(XDIUMOD(NDIUMOD)%io%obj_in), ' (',      &
                                           TRIM(XDIUMOD(NDIUMOD)%io%cha_in), ') -> [',  &
                                           TRIM(XDIUMOD(NDIUMOD)%io%tra_out), ']'
       END DO ! NMAXDIUMOD
    END IF

    CALL p_bcast(NDIUMOD, p_io)

    ! broadcast all results
    DO i=1, NDIUMOD
       ! I/O user interface
       CALL p_bcast(XDIUMOD(i)%io%cha_in,  p_io)
       CALL p_bcast(XDIUMOD(i)%io%obj_in,  p_io)
       CALL p_bcast(XDIUMOD(i)%io%tra_out, p_io)
       CALL p_bcast(XDIUMOD(i)%io%sca_fac, p_io)
       !
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NDIUMOD,'DIUMOD OBJECTS INFO INITIALIZED'
    END IF

! op_pj_20160510+
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif
! op_pj_20160510-

    CALL end_message_bi(modstr,'INITIALISATION ',substr)

  END SUBROUTINE diumod_initialize

! op_pj_20160510+
! ---------------------------------------------------------------------
  SUBROUTINE diumod_init_memory

    IMPLICIT NONE

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
#endif

  END SUBROUTINE diumod_init_memory
! op_pj_20160510-
! ---------------------------------------------------------------------

!> \brief Coupling
!> \details Couples DIUMOD with CHANNEL/TRACER data structures

  SUBROUTINE diumod_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io!, p_io, p_bcast
    USE messy_main_blather_bi,       ONLY: error_bi!, info_bi, warning_bi
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR!, LGTRSTR

    ! MESSy
    USE messy_main_tracer,           ONLY: get_tracer
    USE messy_main_channel,          ONLY: get_channel_object!, get_channel_object_info

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'diumod_init_coupling'
    INTEGER                      :: status
    INTEGER                      :: i
    ! REPRESENTATION IDs
!   INTEGER                      :: rid_in, rid_out

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    a_loop: DO i=1, NDIUMOD

      XDIUMOD(i)%ok = .FALSE.

      ! input field
      NULLIFY(XDIUMOD(i)%dat_in)
      CALL get_channel_object(status    &
           , TRIM(XDIUMOD(i)%io%cha_in) &
           , TRIM(XDIUMOD(i)%io%obj_in)  &
           , p3=XDIUMOD(i)%dat_in)

      IF (status /= 0) THEN
         CALL error_bi('  '//substr//'  input: '//&
              TRIM(XDIUMOD(i)%io%obj_in)//' ('//&
              TRIM(XDIUMOD(i)%io%cha_in)//') - not found' , substr)
         NULLIFY(XDIUMOD(i)%dat_in)
         CYCLE
      END IF

      ! output tracer and tendency
      CALL get_tracer(status, GPTRSTR, TRIM(XDIUMOD(i)%io%tra_out), &
           idx = XDIUMOD(i)%idt_out )
      IF (status /= 0) THEN
         CALL error_bi('  output tracer: '//TRIM(XDIUMOD(i)%io%tra_out)//&
                 &' - not found', substr)
         XDIUMOD(i)%idt_out = 0
         CYCLE
      END IF

      IF (p_parallel_io) THEN
        WRITE(*,*) 'in -> out [TN]: ', TRIM(XDIUMOD(i)%io%obj_in),  ' (',     &
                                       TRIM(XDIUMOD(i)%io%cha_in),  ') -> ',  &
                                       TRIM(XDIUMOD(i)%io%tra_out), ' [',     &
                                       XDIUMOD(i)%idt_out,']'
      END IF

!     CALL get_channel_object_info(status, & 
!             TRIM(XDIUMOD(i)%io%cha_in),   &
!             TRIM(XDIUMOD(i)%io%obj_in),   &
!             reprid=rid_in)
!     CALL get_channel_object_info(status, &
!             TRIM(XDIUMOD(i)%io%cha_out),  &
!             TRIM(XDIUMOD(i)%io%obj_out),  &
!             reprid=rid_out)
!
!     IF ((rid_in /= GP_3D_MID) .OR. (rid_out /= GP_3D_MID)) THEN
!        CALL info_bi( '    '//&
!             &'representations are not supported ... skipping', ' ')
!        NULLIFY(XDIUMOD(i)%dat_in)
!        NULLIFY(XDIUMOD(i)%dat_out)
!        NULLIFY(XDIUMOD(i)%dat_mod)
!        CYCLE
!     END IF

     XDIUMOD(i)%ok = .TRUE.

    END DO a_loop ! NDIUMOD

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE diumod_init_coupling

! ---------------------------------------------------------------------

!> \brief Calculations
!> \details Calculates modulation depending on the solar zenith angle and local time and 
!> stores modulated input fields in specified tracers

  SUBROUTINE diumod_local_start

#ifdef DIUMOD_DEBUG
    USE messy_main_mpi_bi,        ONLY: p_parallel_io!, p_io, p_bcast
#endif
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: qxtte, qxtm1
#endif
    USE messy_main_timer,         ONLY: DAYOFYEAR, HOUR, MINUTE, SECOND, &
                                        time_step_len!, lstart, delta_time
    USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma
    USE messy_main_grid_def_bi,     ONLY:  &
#ifdef DIUMOD_DEBUG
                                        philat_2d, &
#endif
                                        philon_2d, &
                                        sinlat_2d, coslat_2d
    USE messy_main_constants_mem, ONLY: pi

    IMPLICIT NONE

    INTRINSIC MODULO, SIN, COS, ACOS

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'diumod_local_start'
    INTEGER  :: i
    INTEGER  :: jp, jk, jt
#ifndef MESSYTENDENCY
    REAL(DP) :: xten              !< additional nudging tendency
#else
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vstart
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vtend
#endif

    REAL(dp)                              :: dec              !< apx. solar declination
    REAL(DP), DIMENSION(:),   POINTER     :: lst => NULL()    !< local solar time
    REAL(DP), DIMENSION(:),   POINTER     :: a => NULL()      !< amplitude
    REAL(DP), DIMENSION(:),   POINTER     :: c => NULL()      !< baseline 
    REAL(DP), DIMENSION(:),   POINTER     :: h => NULL()      !< interm. calc.
    REAL(DP), DIMENSION(:),   POINTER     :: m => NULL()      !< modulator

#ifdef MESSYTENDENCY
    ALLOCATE(vstart(kproma,nlev))
    ALLOCATE(vtend(kproma,nlev))
#endif
    ALLOCATE(lst(kproma))
    ALLOCATE(a(kproma))
    ALLOCATE(c(kproma))
    ALLOCATE(h(kproma))
    ALLOCATE(m(kproma))

!> calculating parameters

   lst(1:kproma) = MODULO( &
        ( ( REAL(HOUR,dp)+REAL(MINUTE,dp)/60._dp+REAL(SECOND,dp)/3600._dp + 12._dp ) + &
               24._dp/360._dp*philon_2d(1:kproma,jrow) ), 24._dp )/12._dp*pi

   dec = ( -23.44_dp * COS( 2*pi * (REAL(DAYOFYEAR,dp)+10._dp)/365._dp ) ) / 180._dp * pi

   a(1:kproma) = coslat_2d(1:kproma,jrow) * COS(dec)
   c(1:kproma) = -sinlat_2d(1:kproma,jrow) * SIN(dec)

   WHERE ( c .lt. a ) c = -c

   WHERE ( c .ge. a )
     !> polar day/night conditions
     m = 1._dp + a/c * COS( lst )
   ELSEWHERE
     !> all other cases
     h = ACOS( -c/a )
     m = ( c + a*COS( lst ) ) / ( ( c*h + a*SIN(h) ) / pi)
   ENDWHERE

!   m(:) = 1._dp

!   WHERE (c(1:kproma) .gt. a(1:kproma))
!     m(1:kproma) = 1._dp + a(1:kproma)/c(1:kproma)*COS( lst(1:kproma) )
!   ELSEWHERE
!     h(1:kproma) = ACOS( -c(1:kproma)/a(1:kproma) )
!     m(1:kproma) = ( c(1:kproma) + a(1:kproma)*COS( lst(1:kproma) ) ) / &
!                   ( ( c(1:kproma)*h(1:kproma) + a(1:kproma)*SIN(h(1:kproma)) ) / pi)
!
!   ENDWHERE

#ifdef DIUMOD_DEBUG
   IF (p_parallel_io) THEN
      jp = 1
      WRITE(*,*) substr,' diag ( @PIO @jp = ',jp,' ) ---------------------'
      WRITE(*,*) '     DOY: ',DAYOFYEAR
      WRITE(*,*) ' lat/lon: ',philat_2d(jp,jrow),' / ',philon_2d(jp,jrow)
      WRITE(*,*) '     dec: ',dec
      WRITE(*,*) '    hour: ',(REAL(HOUR,dp)+REAL(MINUTE,dp)/60._dp+REAL(SECOND,dp)/3600._dp)
      WRITE(*,*) '  lst(h): ',lst(jp)/pi*12_dp
      WRITE(*,*) '       a: ',a(jp)
      WRITE(*,*) '       c: ',c(jp)
      WRITE(*,*) '       h: ',h(jp)
      WRITE(*,*) '       m: ',m(jp)
   END IF
#endif


    !> looping over DIUMOD entries
    mod_loop: DO i=1, NDIUMOD
      IF (XDIUMOD(i)%ok) THEN      !< if active
        jt = XDIUMOD(i)%idt_out    !< output tracer idx
#ifdef MESSYTENDENCY
        ! SET START VALUE AND INIT TENDENCY
!        CALL mtend_get_start_l(mtend_id_tracer, v0=vstart, idt=jt)
        CALL mtend_get_start_l(jt, v0=vstart)
        vtend(:,:) = 0.0_dp
#endif
        level_loop: DO jk=1, nlev  ! LOOP OVER LEVELS
           vector_loop: DO jp=1, kproma
#ifndef MESSYTENDENCY
             xten = -1._dp * &
               ( ( qxtm1(jp,jk,jt)+qxtte(jp,jk,jt)*time_step_len ) - &
                 ( XDIUMOD(i)%dat_in(jp,jk,jrow)*XDIUMOD(i)%io%sca_fac*m(jp) ) ) / time_step_len
             qxtte(jp,jk,jt) = qxtte(jp,jk,jt) + xten
#else
! op_pj_20160510+: unclear, if this fixis correct
!!$             vtend(jp,jk) = -1._DP * ( vstart(jp,jk) - &
!!$               XDIUMOD(i)%field(jp,jk,jt))
             vtend(jp,jk) = -1._DP * ( ( vstart(jp,jk) - &
                 ( XDIUMOD(i)%dat_in(jp,jk,jrow)*XDIUMOD(i)%io%sca_fac*m(jp) ) ) ) / time_step_len
! op_pj_20160510-
#endif
           END DO vector_loop
        END DO level_loop
#ifdef MESSYTENDENCY
!        CALL mtend_add_l(my_handle, mtend_id_tracer, px=vtend, idt=jt)
        CALL mtend_add_l(my_handle, jt, px=vtend)
#endif
      END IF

    END DO mod_loop

    !> CLEAN UP
#ifdef MESSYTENDENCY
    DEALLOCATE(vstart)
    DEALLOCATE(vtend)
#endif
    DEALLOCATE(lst)
    DEALLOCATE(a)
    DEALLOCATE(c)
    DEALLOCATE(h)
    DEALLOCATE(m)

  END SUBROUTINE diumod_local_start


! ---------------------------------------------------------------------
!> \brief De-initialisation
!> \details Frees up memory used for modulation jobs

  SUBROUTINE diumod_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    IF (ASSOCIATED(XDIUMOD)) DEALLOCATE(XDIUMOD)
    NULLIFY(XDIUMOD)

  END SUBROUTINE diumod_free_memory
! ---------------------------------------------------------------------


! ======================================================================
! PRIVATE ROUTINES
! ======================================================================

! ----------------------------------------------------------------------
!> \brief Reads namelist
!> \details Checks for presence and reads DIUMOD namelist

  SUBROUTINE diumod_read_nml_cpl(status, iou)

    ! DIUMOD MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling'
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2007

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     !< error status
    INTEGER, INTENT(IN)  :: iou        !< I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'diumod_read_nml_cpl'

    NAMELIST /CPL/ DIUMOD

    LOGICAL :: lex   !< file exists?
    INTEGER :: fstat !< file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not. lex) RETURN   ! namelist file (<modstr>.nml) not available

    READ (iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! ...done, without error


  END SUBROUTINE diumod_read_nml_cpl


! **********************************************************************
END MODULE messy_diumod_si
! **********************************************************************
