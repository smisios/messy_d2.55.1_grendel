#include "messy_main_ppd_bi.inc"

MODULE MESSY_AEROPT_SI

!  AUTHORS:   - Holger Tost, MPI Chemie, Mainz
!             - CESM1 interface (_c1) by Andreas Baumgertner
!             - merged _e5 and _c1 into _si by Patrick Joeckel

#if defined(ECHAM5) || defined (CESM1) || defined (MBM_RAD)

  ! ECHAM5/MESSy
  USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks
#ifdef CESM1
  USE messy_main_grid_def_mem_bi, ONLY: npromz
#endif
  USE messy_main_mpi_bi,        ONLY: p_bcast,p_io
  USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, ntrac => ntrac_gp, ti => ti_gp, &
                                      t_trinfo_tp
  ! MESSy
  USE messy_main_tracer,        ONLY: get_tracer

  USE messy_aeropt_mem

  IMPLICIT NONE

  ! number of calls for the optical properties
  INTEGER,SAVE :: n_calls = 0
  INTEGER,SAVE :: n_lut   = 0
  INTEGER,SAVE :: n_inp   = 0
  INTEGER,SAVE :: n_merge = 0
  LOGICAL,SAVE :: lconv_req = .FALSE.
  REAL(dp), DIMENSION(:), POINTER, SAVE :: lw_bands => NULL()
  REAL(dp), DIMENSION(:), POINTER, SAVE :: sw_bands => NULL()
  REAL(dp), POINTER, SAVE :: plw, psw

  TYPE t_optset_io
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel_aot_lw   = ""
     CHARACTER(LEN=STRLEN_OBJECT)  :: object_aot_lw    = ""
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel_aot_sw   = ""
     CHARACTER(LEN=STRLEN_OBJECT)  :: object_aot_sw    = ""
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel_gamma    = ""
     CHARACTER(LEN=STRLEN_OBJECT)  :: object_gamma     = ""
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel_omega    = ""
     CHARACTER(LEN=STRLEN_OBJECT)  :: object_omega     = ""
     CHARACTER(LEN=STRLEN_MEDIUM)  :: merge_name1      = ""
     CHARACTER(LEN=STRLEN_MEDIUM)  :: merge_name2      = ""
  END type t_optset_io

  TYPE t_optset
     TYPE(t_optset_io)          :: io
     TYPE(t_aero_set),  POINTER :: as
  END type t_optset

  TYPE(t_optset), DIMENSION(:), POINTER :: optset
 
CONTAINS


!===============================================================================

  SUBROUTINE AEROPT_INITIALIZE

    USE messy_main_tools,        ONLY: find_next_free_unit
    USE messy_main_blather_bi,   ONLY: start_message_bi, end_message_bi, &
                                       error_bi, warning_bi
    USE messy_main_tools,        ONLY: match_wild
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io,     &
                                       p_bcast
    USE messy_main_grid_def_bi,  ONLY: cetah
    USE messy_aeropt,            ONLY: aeropt_read_nml_ctrl
    USE messy_aeropt_tanre,      ONLY: su_aero_tanre, tanre_init_su
    USE messy_main_import_bi,    ONLY: import_grid_submodstr

    IMPLICIT NONE

    ! LOCAL
    ! name of subroutine
    CHARACTER(LEN=*), PARAMETER     :: substr = 'aeropt_initialize' 
    INTEGER                         :: iou    ! I/O unit
    INTEGER                         :: status ! error status
    TYPE nml_set_struct
      CHARACTER(LEN=STRLEN_MEDIUM)  :: name             = ""
      LOGICAL                       :: tracer           =.TRUE.
      CHARACTER(LEN=STRLEN_MEDIUM)  :: tracer_set       = ""
      CHARACTER(LEN=256)            :: exclude_str_full = ""
      CHARACTER(LEN=256)            :: exclude_str_spec = ""
      LOGICAL                       :: lcalc_seasalt    =.FALSE.
      LOGICAL                       :: lextmix          =.FALSE.
      CHARACTER(LEN=STRLEN_MEDIUM)  :: aermodel         = ""
      CHARACTER(LEN=STRLEN_MEDIUM)  :: wetradius        = ""
      CHARACTER(LEN=STRLEN_MEDIUM)  :: dryradius        = ""
      CHARACTER(LEN=STRLEN_MEDIUM)  :: aernumber        = ""
      INTEGER                       :: lut_number       = 0
      LOGICAL                       :: lcalc_jval       = .FALSE.
      REAL(dp), DIMENSION(max_diag_wavelens) :: diag_wavelength = 0.0_dp
    END TYPE nml_set_struct
   
    TYPE(nml_set_struct), DIMENSION(MAX_AEROPT_CALLS) :: read_aero_sets
   
    INTEGER                   :: count, i
    CHARACTER(LEN=2)          :: idx_nr
    LOGICAL                   :: found

    TYPE import_struct
      CHARACTER(LEN=STRLEN_MEDIUM)  :: name             = ""
      CHARACTER(LEN=STRLEN_CHANNEL) :: channel_aot_lw   = ""
      CHARACTER(LEN=STRLEN_OBJECT)  :: object_aot_lw    = ""
      CHARACTER(LEN=STRLEN_CHANNEL) :: channel_aot_sw   = ""
      CHARACTER(LEN=STRLEN_OBJECT)  :: object_aot_sw    = ""
      CHARACTER(LEN=STRLEN_CHANNEL) :: channel_gamma    = ""
      CHARACTER(LEN=STRLEN_OBJECT)  :: object_gamma     = ""
      CHARACTER(LEN=STRLEN_CHANNEL) :: channel_omega    = ""
      CHARACTER(LEN=STRLEN_OBJECT)  :: object_omega     = ""
      INTEGER                       :: inp_dim          = 4
    END TYPE import_struct
    TYPE(import_struct), DIMENSION(MAX_AEROPT_CALLS) :: read_input_sets

    TYPE merge_struct
       CHARACTER(LEN=STRLEN_MEDIUM)  :: name             = ""
       CHARACTER(LEN=STRLEN_MEDIUM)  :: merge_name1      = ""
       REAL(dp)                      :: weight1          = 1.0
       CHARACTER(LEN=STRLEN_MEDIUM)  :: merge_name2      = ""
       REAL(dp)                      :: weight2          = 1.0
       REAL(dp)                      :: pressure_bound1  = 0.0
       REAL(dp)                      :: pressure_bound2  = 0.0
    END type merge_struct
    TYPE(merge_struct), DIMENSION(MAX_AEROPT_CALLS) :: read_merge_sets

    NAMELIST /CPL/      read_aero_sets
    NAMELIST /CPL_IMP/  read_input_sets
    NAMELIST /CPL_MERGE/read_merge_sets
    

    CALL start_message_bi(modstr,'INITIALIZATION',substr)

    !--- Read namelist and control variables:

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      ! *** CALL CORE ROUTINE:
      CALL aeropt_read_nml_ctrl(status, iou)
      IF (status /= 0)  CALL error_bi('error in aeropt_read_nml_ctrl', substr)

      !--- Read CPL namelist

      iou = find_next_free_unit(100,200)
      CALL aeropt_read_nml_cpl(status, iou)
      IF (status /= 0) CALL error_bi('error in aeropt_read_nml_cpl', substr)
      
      count = 0
      DO i=1,MAX_AEROPT_CALLS
        IF ( TRIM(read_aero_sets(i)%name) == "" ) CYCLE
        count = count + 1
      END DO
      n_calls = count

      count = 0
      DO i=1,MAX_AEROPT_CALLS
        IF ( read_lut_sets(i)%lut_number == 0 ) CYCLE
        count = count + 1
      END DO
      n_lut = count

      count = 0
      DO i=1,MAX_AEROPT_CALLS
        IF ( TRIM(read_input_sets(i)%name) == "" ) CYCLE
        count = count + 1
      END DO
      n_inp = count

      count = 0
      DO i=1,MAX_AEROPT_CALLS
         IF ( TRIM(read_merge_sets(i)%name) == "" ) CYCLE
         count = count + 1
      END DO
      n_merge = count
    END IF

    CALL p_bcast(n_calls, p_io)
    CALL p_bcast(n_lut, p_io)
    CALL p_bcast(n_inp, p_io)
    CALL p_bcast(n_merge, p_io)

    ALLOCATE(optset(n_calls + n_inp + n_merge))
    DO i=1, n_calls + n_inp + n_merge
       ALLOCATE(optset(i)%as)
    END DO

    ALLOCATE(tab_set(n_lut))

    IF (p_parallel_io ) THEN
      count = 0
      DO i=1,MAX_AEROPT_CALLS
        IF ( TRIM(read_aero_sets(i)%name) == "" ) CYCLE
        count = count + 1

        optset(count)%as%aero_set_name       = read_aero_sets(i)%name
        optset(count)%as%l_tracer            = read_aero_sets(i)%tracer
        optset(count)%as%tracerset           = read_aero_sets(i)%tracer_set
        optset(count)%as%exclude_string_full = read_aero_sets(i)%exclude_str_full
        optset(count)%as%exclude_string_spec = read_aero_sets(i)%exclude_str_spec
        optset(count)%as%lut_number          = read_aero_sets(i)%lut_number
        optset(count)%as%aermodelname        = read_aero_sets(i)%aermodel
        optset(count)%as%lcalc_seasalt       = read_aero_sets(i)%lcalc_seasalt
        optset(count)%as%lcalc_jval          = read_aero_sets(i)%lcalc_jval
        optset(count)%as%l_extmixt           = read_aero_sets(i)%lextmix
        optset(count)%as%diag_wavelen(:)     = read_aero_sets(i)%diag_wavelength(:)
      END DO
    ENDIF

    DO i=1,n_calls
      call p_bcast(optset(i)%as%aero_set_name       , p_io)
      call p_bcast(optset(i)%as%l_tracer            , p_io)
      call p_bcast(optset(i)%as%tracerset           , p_io)
      call p_bcast(optset(i)%as%exclude_string_full , p_io)
      call p_bcast(optset(i)%as%exclude_string_spec , p_io)
      call p_bcast(optset(i)%as%lut_number          , p_io)
      call p_bcast(optset(i)%as%aermodelname        , p_io)
      call p_bcast(optset(i)%as%lcalc_seasalt       , p_io)
      call p_bcast(optset(i)%as%lcalc_jval          , p_io)
      call p_bcast(optset(i)%as%l_extmixt           , p_io)
      call p_bcast(optset(i)%as%diag_wavelen(:)     , p_io)
    END DO

    IF (p_parallel_io ) THEN
      count = 0
      DO i=1,MAX_AEROPT_CALLS
        IF ( read_lut_sets(i)%lut_number == 0 ) CYCLE
        count = count + 1
        tab_set(count)%table_number    = read_lut_sets(i)%lut_number
        tab_set(count)%rad_sw_filename = read_lut_sets(i)%sw_filename
        tab_set(count)%rad_lw_filename = read_lut_sets(i)%lw_filename
      END DO
    ENDIF

    DO i=1,n_lut
      call p_bcast(tab_set(i)%table_number    , p_io)
      call p_bcast(tab_set(i)%rad_sw_filename , p_io)
      call p_bcast(tab_set(i)%rad_lw_filename , p_io)
    END DO

    IF (p_parallel_io ) THEN
      count = n_calls
      DO i=1,MAX_AEROPT_CALLS
        IF ( TRIM(read_input_sets(i)%name) == "" ) CYCLE
        count = count + 1   
        optset(count)%as%aero_set_name       = read_input_sets(i)%name
        optset(count)%io%channel_aot_lw      = read_input_sets(i)%channel_aot_lw
        optset(count)%io%object_aot_lw       = read_input_sets(i)%object_aot_lw
        optset(count)%io%channel_aot_sw      = read_input_sets(i)%channel_aot_sw
        optset(count)%io%object_aot_sw       = read_input_sets(i)%object_aot_sw
        optset(count)%io%channel_gamma       = read_input_sets(i)%channel_gamma
        optset(count)%io%object_gamma        = read_input_sets(i)%object_gamma
        optset(count)%io%channel_omega       = read_input_sets(i)%channel_omega
        optset(count)%io%object_omega        = read_input_sets(i)%object_omega
        optset(count)%as%inp_dim             = read_input_sets(i)%inp_dim
      END DO
    END IF

    DO i=n_calls+1,n_calls+n_inp
      call p_bcast(optset(i)%as%aero_set_name , p_io)
      call p_bcast(optset(i)%io%channel_aot_lw, p_io)
      call p_bcast(optset(i)%io%object_aot_lw , p_io)
      call p_bcast(optset(i)%io%channel_aot_sw, p_io)
      call p_bcast(optset(i)%io%object_aot_sw , p_io)
      call p_bcast(optset(i)%io%channel_gamma , p_io)
      call p_bcast(optset(i)%io%object_gamma  , p_io)
      call p_bcast(optset(i)%io%channel_omega , p_io)
      call p_bcast(optset(i)%io%object_omega  , p_io)
      call p_bcast(optset(i)%as%inp_dim       , p_io)
    END DO

    IF (p_parallel_io ) THEN
       count = n_calls + n_inp
       DO i=1,MAX_AEROPT_CALLS
          IF ( TRIM(read_merge_sets(i)%name) == "" ) CYCLE
          count = count + 1
          optset(count)%as%aero_set_name       = read_merge_sets(i)%name
          optset(count)%io%merge_name1         = read_merge_sets(i)%merge_name1
          optset(count)%io%merge_name2         = read_merge_sets(i)%merge_name2
          optset(count)%as%weight1             = REAL(read_merge_sets(i)%weight1,dp)
          optset(count)%as%weight2             = REAL(read_merge_sets(i)%weight2,dp)
          optset(count)%as%press_bound1        = REAL(read_merge_sets(i)%pressure_bound1,dp)
          optset(count)%as%press_bound2        = REAL(read_merge_sets(i)%pressure_bound2,dp)
       ENDDO
    ENDIF

    DO i=n_calls+n_inp+1, n_calls+n_inp+n_merge
      call p_bcast(optset(i)%as%aero_set_name , p_io)
      call p_bcast(optset(i)%io%merge_name1   , p_io)
      call p_bcast(optset(i)%io%merge_name2   , p_io)
      call p_bcast(optset(i)%as%weight1       , p_io)
      call p_bcast(optset(i)%as%weight2       , p_io)
      call p_bcast(optset(i)%as%press_bound1  , p_io)
      call p_bcast(optset(i)%as%press_bound2  , p_io)
    END DO

    DO i=1,n_calls+n_inp+n_merge
       optset(i)%as%lconv_sw = ( TRIM(optset(i)%io%channel_aot_sw) == &
            TRIM(import_grid_submodstr) )
       !
       IF (TRIM(optset(i)%io%channel_aot_sw) /= '') THEN
          IF (optset(i)%as%lconv_sw) THEN
             CALL warning_bi('object '//TRIM(optset(i)%io%object_aot_sw)//&
                  ' imported from channel '//&
                  TRIM(optset(i)%io%channel_aot_sw)//'! Assuming unit [1/m],'//&
                  ' therefore scaling accordingly ...' &
                  , substr)
          ELSE
             CALL warning_bi('object '//TRIM(optset(i)%io%object_aot_sw)//&
                  ' provided by channel '//&
                  TRIM(optset(i)%io%channel_aot_sw)//'! Assuming unit [1/box],'//&
                  ' therefore no scaling ...' &
                  , substr)
          END IF
       END IF

       optset(i)%as%lconv_lw = (TRIM(optset(i)%io%channel_aot_lw) == &
            TRIM(import_grid_submodstr) )
       !
       IF (TRIM(optset(i)%io%channel_aot_lw) /= '') THEN
          IF (optset(i)%as%lconv_lw) THEN
             CALL warning_bi('object '//TRIM(optset(i)%io%object_aot_lw)//&
                  ' imported from channel '//&
                  TRIM(optset(i)%io%channel_aot_lw)//'! Assuming unit [1/m],'//&
                  ' therefore scaling accordingly ...' &
                  , substr)
          ELSE
             CALL warning_bi('object '//TRIM(optset(i)%io%object_aot_lw)//&
                  ' provided by channel '//&
                  TRIM(optset(i)%io%channel_aot_lw)//'! Assuming unit [1/box],'//&
                  ' therefore no scaling ...' &
                  , substr)
          END IF
       END IF

       lconv_req = lconv_req .OR. optset(i)%as%lconv_sw .OR. optset(i)%as%lconv_lw
    END DO

    DO i=1,n_calls
       optset(i)%as%tanre=.FALSE.
      IF ( .NOT. (optset(i)%as%l_tracer) )  THEN
        IF (MATCH_WILD( 'TANRE*', TRIM(optset(i)%as%aero_set_name))) THEN
          optset(i)%as%tanre = .TRUE.
          IF (.NOT. tanre_init_su) &
            CALL su_aero_tanre(cetah)
        ENDIF
      ENDIF
    ENDDO

    DO i=1,n_calls
       IF (i < 10) then
          write (idx_nr,'(A,I1)') "0",i
        ELSE
          write (idx_nr,'(I2)') i
        END IF
      IF (optset(i)%as%TANRE) CYCLE
      found = .FALSE.
      DO count=1,n_lut
        IF (optset(i)%as%lut_number == tab_set(count)%table_number) &
          FOUND = .TRUE.
      END DO
      
      IF (.NOT. FOUND) CALL error_bi(&
        'Table number of optset '//idx_nr//' does not exist!', substr)
    END DO

    

    CALL end_message_bi(modstr,'INITIALIZATION',substr)

!--------------------------------------------------------------------------------
    CONTAINS

!--------------------------------------------------------------------------------
      SUBROUTINE aeropt_read_nml_cpl(status, iou)

        ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES

        USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

        IMPLICIT NONE
        SAVE
        
        ! I/O
        INTEGER, INTENT(OUT) :: status ! error status
        INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
        
        !--- Local variables:

        CHARACTER(LEN=*), PARAMETER :: substr = 'aeropt_read_nml_cpl'
        LOGICAL                     :: lex          ! file exists ?
        INTEGER                     :: fstat        ! file status

        ! INITIALIZE
        
        status = 1 ! ERROR

        ! INITIALIZE GLOBAL CONTROL VARIABLES
        ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

        !--- 1) Read CPL namelist:

        CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
        IF (.not.lex) RETURN    ! <modstr>.nml does not exist

        READ(iou, NML=CPL, IOSTAT=fstat)
        CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
        IF (fstat /= 0) RETURN  ! error while reading namelist
    

        CALL read_nml_close(substr, iou, modstr)

        !--- 2) Read CPL_IMP namelist:

        CALL read_nml_open(lex, substr, iou+1, 'CPL_IMP', modstr)
        IF (.not.lex) RETURN    ! <modstr>.nml does not exist

        READ(iou+1, NML=CPL_IMP, IOSTAT=fstat)
        CALL read_nml_check(fstat, substr, iou+1, 'CPL_IMP', modstr)
        IF (fstat /= 0) RETURN  ! error while reading namelist
    
        !--- 3) Read CPL_MERGE namelist:

        CALL read_nml_open(lex, substr, iou+1, 'CPL_MERGE', modstr)
        IF (.not.lex) RETURN    ! <modstr>.nml does not exist

        READ(iou+1, NML=CPL_MERGE, IOSTAT=fstat)
        CALL read_nml_check(fstat, substr, iou+1, 'CPL_MERGE', modstr)
        IF (fstat /= 0) RETURN  ! error while reading namelist
        
        CALL read_nml_close(substr, iou+1, modstr)

        status = 0  ! no ERROR


      END SUBROUTINE aeropt_read_nml_cpl
!-------------------------------------------------------------------------------
  END SUBROUTINE AEROPT_INITIALIZE
!===============================================================================

  SUBROUTINE AEROPT_INIT_MEMORY

    ! ECHAM5/MESSy
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
 

    CHARACTER(LEN=*), PARAMETER :: substr ='aeropt_init_memory'

    CALL start_message_bi(modstr,'MEMORY INITIALIZATION',substr)


    CALL end_message_bi(modstr,'MEMORY INITIALIZATION',substr)

  END SUBROUTINE AEROPT_INIT_MEMORY

!===============================================================================

  SUBROUTINE AEROPT_INIT_COUPLING

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, dcl, reorder
#ifdef ECHAM5
    USE messy_main_mpi_bi,           ONLY: dcl
    USE messy_main_grid_def_bi,      ONLY: coslon, sinlon, gl_twomu
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
#endif
#ifdef CESM1
    USE messy_main_grid_def_mem_bi, ONLY: nlev
    USE messy_main_grid_def_bi,     ONLY: coslon_2d, sinlon_2d, sinlat_2d
#endif
#ifdef MBM_RAD
    USE messy_main_grid_def_bi,     ONLY: coslon, sinlon, gl_twomu
    USE messy_main_grid_def_mem_bi, ONLY:nlev &
                                      , nlon, nlat                     &
                                      , bi_glon => glon                &
                                      , bi_glat => glat
#endif
    USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                        error_bi, info_bi
    USE messy_main_channel,       ONLY: get_channel_object, get_channel_info, &
                                        new_channel_object, new_attribute,    &
                                        new_channel
    USE messy_main_channel_error_bi,   ONLY: channel_halt 
    USE messy_main_channel_bi,    ONLY: GP_3D_MID, DC_GP  &
                                      , DIMID_LON, DIMID_LAT, DIMID_LEV &
                                      , gp_nseg, gp_start, gp_cnt &
                                      , gp_meml, gp_memu &
                                      , GP_2D_HORIZONTAL
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO &
                                           , set_representation_decomp &
                                           , IRANK, PIOTYPE_COL &
                                           , REPR_DEF_AXES
                                        
    ! MESSy
    USE messy_main_tools,         ONLY: strcrack, match_wild
    USE messy_main_constants_mem, ONLY: M_cl => MCl, rho_h2o, m_h2o, m_na => MNa
    USE messy_main_tracer,        ONLY: r_molarmass, r_aerosol_density, &
                                        i_aerosol_mode, aerosol, family
    
    ! MESSy/AEROPT
    USE messy_aeropt_mem,         ONLY: nmod, aerspec, num_diag_elements,     &
                                        max_diag_wavelens,                    &
                                        diag_bc, diag_oc,                     &
                                        diag_du, diag_ss, diag_waso, diag_h2o,&
                                        diag_tot, &
                                        ri_h2o, ri_waso, ri_bc, ri_ss, ri_oc, &
                                        ri_du, nsol
    USE messy_aeropt,             ONLY: aeropt_initialise
    USE messy_aeropt_sets,        ONLY: map_optset_struct
    USE messy_aeropt_input,       ONLY: aeropt_check_lut, aeropt_finalise_lut
    USE messy_aeropt_tanre,       ONLY: zaes_x, zael_x, zaeu_x, zaed_x,       &
                                        tanre_init, naer
#if defined(ECHAM5) || defined(MBM_RAD)
    USE messy_aeropt_tanre,       ONLY: init_aero_tanre
#endif
#ifdef CESM1
     USE messy_aeropt_tanre,      ONLY: init_aero_tanre_dc
#endif


#ifdef ECHAM5
    REAL(dp):: laes(dcl%nglon, dcl%nglat)
    REAL(dp):: lael(dcl%nglon, dcl%nglat)
    REAL(dp):: laeu(dcl%nglon, dcl%nglat)
    REAL(dp):: laed(dcl%nglon, dcl%nglat)
#endif
#ifdef MBM_RAD
    REAL(dp):: laes(nlon, nlat)
    REAL(dp):: lael(nlon, nlat)
    REAL(dp):: laeu(nlon, nlat)
    REAL(dp):: laed(nlon, nlat)
#endif

    INTEGER          :: status
    CHARACTER(LEN=2) :: nm, inp_set_nr               
    CHARACTER(LEN=3) :: nmo                             
    CHARACTER(LEN=4) :: wl    
    CHARACTER(LEN=*), PARAMETER :: substr ='aeropt_init_coupling'
                                                        
    INTEGER  :: jm, n, i, ji, jt, j, k
#ifdef CESM1
    INTEGER  ::  zjrow, zkproma
#endif
    INTEGER  :: counter, counter2, dummy
    CHARACTER(LEN=26), POINTER     :: strname_full(:) => NULL()
    CHARACTER(LEN=26), POINTER     :: strname_spec(:) => NULL()
    CHARACTER(LEN=26), POINTER     :: strname     (:) => NULL()
    INTEGER  :: n_out_full, n_out_spec
    LOGICAL  :: no_count
    LOGICAL  :: Na_found  = .FALSE.
    LOGICAL  :: Cl_found  = .FALSE.
    LOGICAL  :: H2O_found = .FALSE.

  ! auxiliary pointer for data managment
    REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  p4 => NULL()
    REAL(dp), DIMENSION(:,:,:),     POINTER ::  p3 => NULL()

    INTEGER                               :: DIMID_NSW
    INTEGER                               :: DIMID_JPBAND
    INTEGER                               :: DIMID_JV
    INTEGER                               :: REPR_AEROPT_4D_NSW
    INTEGER                               :: REPR_AEROPT_4D_JPBAND
    INTEGER                               :: REPR_AEROPT_4D_JVAL
    CHARACTER(LEN=1)                      :: ichar

    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    LOGICAL :: cpl_rad = .FALSE.
    CHARACTER(LEN=STRLEN_CHANNEL) :: radiationmodelname=""
    CHARACTER(LEN=256)  :: cname
    CHARACTER(LEN=256)  :: cname2
    CHARACTER(LEN=2)    :: str_num
    
    REAL(dp) :: delta

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)



    IF (CPL_rad) THEN
    ! get nsw (number of shortwave wavenlengths) and 
    !     jpband (number of longwave wavelengths) 
    ! from the radiation scheme

      CALL get_channel_info(status, TRIM(radiationmodelname))
      CALL channel_halt(&
        substr//': channel for the radiation model not found: '//TRIM(radiationmodelname), &
        status)

      CALL get_channel_object(status, TRIM(radiationmodelname), 'nsw', p0=psw)
      CALL channel_halt(substr//': channel object for nsw not found!', status)
      CALL get_channel_object(status, TRIM(radiationmodelname), 'sw_bands', p1=sw_bands)
      CALL channel_halt(substr//': channel object for band_bounds (sw) not found!', status)
      CALL get_channel_object(status, TRIM(radiationmodelname), 'jpband', p0=plw)
      CALL channel_halt(substr//': channel object for jpband not found!', status)
      CALL get_channel_object(status, TRIM(radiationmodelname), 'lw_bands', p1=lw_bands)
      CALL channel_halt(substr//': channel object for band_bounds (lw) not found!', status)
      
      nsw = NINT(psw)
      jpband = NINT(plw)
    END IF

     ! 4D NSW representation / decomposition
    CALL new_dimension(status, DIMID_NSW, 'AEROPT_NSW', nsw)
    CALL channel_halt(substr, status)
    CALL new_representation(status, REPR_AEROPT_4D_NSW, &
         'REPR_AEROPT_4D_NSW', rank = 4, link = 'xxxx', dctype = DC_GP      &
         , dimension_ids = &
         (/ _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_NSW) /)        &
         , ldimlen      = (/ _RI_XYZN_(nproma , ngpblks, AUTO, AUTO) /)     &
         , output_order = (/ _IN_XYZN_, _IX_XYZN_ , _IY_XYZN_, _IZ_XYZN_ /) &
         , axis = repr_def_axes(_RI_XYZN_('X','Y','Z','N') )                &
         )
    CALL channel_halt(substr, status)

    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))
    
    start(:,:) = gp_start(:,:)
    cnt(:,:) = gp_cnt(:,:)
    meml(:,:) = gp_meml(:,:)
    memu(:,:) = gp_memu(:,:)
    
    cnt(:, _IN_XYZN_) = nsw
    memu(:, _IN_XYZN_) = nsw
    
    CALL set_representation_decomp(status, REPR_AEROPT_4D_NSW &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    ! 4D JPBAND represenation / decomposition

    CALL new_dimension(status, DIMID_JPBAND, 'AEROPT_JPBAND', jpband)
    CALL channel_halt(substr, status)
    CALL new_representation(status, REPR_AEROPT_4D_JPBAND, &
         'REPR_AEROPT_4D_JPBAND', rank = 4, link = 'xxxx', dctype = DC_GP     &
         , dimension_ids =  &
         (/ _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_JPBAND) /)       &
         , ldimlen       = (/ _RI_XYZN_(nproma , ngpblks, AUTO, AUTO) /)      &
         , output_order  = (/ _IN_XYZN_, _IX_XYZN_ , _IY_XYZN_, _IZ_XYZN_ /)  &
         , axis = repr_def_axes(_RI_XYZN_('X','Y','Z','N'))                   &
         )
    CALL channel_halt(substr, status)

    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))
    
    start(:,:) = gp_start(:,:)
    cnt(:,:) = gp_cnt(:,:)
    meml(:,:) = gp_meml(:,:)
    memu(:,:) = gp_memu(:,:)
    
    cnt(:,_IN_XYZN_) = jpband
    memu(:,_IN_XYZN_) = jpband
    
    CALL set_representation_decomp(status, REPR_AEROPT_4D_JPBAND &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu) 

    ! 4D JVAL representation / decomposition
    CALL new_dimension(status, DIMID_JV, 'AEROPT_JV', n_jv)
    CALL channel_halt(substr, status)
    CALL new_representation(status, REPR_AEROPT_4D_JVAL, &
         'REPR_AEROPT_4D_JVAL', rank = 4, link = 'xxxx', dctype = DC_GP      &
         , dimension_ids =  &
         (/  _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_JV) /)         &
         , ldimlen       = (/ _RI_XYZN_(nproma , ngpblks, AUTO, AUTO) /)     &
         , output_order  = (/ _IN_XYZN_, _IX_XYZN_ , _IY_XYZN_, _IZ_XYZN_ /) &
         , axis = repr_def_axes(_RI_XYZN_('X','Y','Z','N'))                  &
         )
    CALL channel_halt(substr, status)

    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))
    
    start(:,:) = gp_start(:,:)
    cnt(:,:) = gp_cnt(:,:)
    meml(:,:) = gp_meml(:,:)
    memu(:,:) = gp_memu(:,:)
    
    cnt(:,3) = n_jv
    memu(:,3) = n_jv
    
    CALL set_representation_decomp(status, REPR_AEROPT_4D_JVAL &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    DO i=1,n_calls
      ! loop over the aeropt_calls to get the information from the desired
      ! aerosol submodel(s)
      IF ( optset(i)%as%TANRE ) CYCLE
      CALL get_channel_info(status, TRIM(optset(i)%as%aermodelname))
      CALL channel_halt(&
        substr//': channel for the aerosol model not found: '//TRIM(optset(i)%as%aermodelname), &
        status)
      CALL get_channel_object(status, TRIM(optset(i)%as%aermodelname), 'sigma', p1=optset(i)%as%sigma)
      CALL channel_halt(substr//': channel object for sigma not found!', status)

      CALL get_channel_object(status, TRIM(optset(i)%as%aermodelname), 'wetradius', &
        p4=optset(i)%as%wetradius)
      CALL channel_halt(substr//': channel object for wet radius not found!', status)
      CALL get_channel_object(status, TRIM(optset(i)%as%aermodelname), 'dryradius', &
        p4=optset(i)%as%dryradius)
      CALL channel_halt(substr//': channel object for dry radius not found!', status)
      CALL get_channel_object(status, TRIM(optset(i)%as%aermodelname), 'anumber', &
        p4=optset(i)%as%aernumber)
      IF (status /= 0) &
        CALL info_bi(' channel object for aerosol number not found!'//&
        ' Checking for tracer for aerosol number concentrations!', substr)
      
      optset(i)%as%nmod = size(optset(i)%as%sigma)
      IF (optset(i)%as%nmod .GT. max_modes) CALL error_bi(&
        'Number of modes too large! Recompile with increased max_modes.' &
        , substr)

      call strcrack(optset(i)%as%aermodelname,"_",strname,dummy)

      SELECT CASE (TRIM(strname(1)))

      CASE('m7', 'M7', 'gmxe', 'GMXE', 'GMXe')
        optset(i)%as%NS = 1
        optset(i)%as%KS = 2
        optset(i)%as%AS = 3
        optset(i)%as%CS = 4

        DO jm = 1, 3
          optset(i)%as%radmod(jm+1) = jm
          optset(i)%as%radmod(jm+4) = jm
        END DO

        optset(i)%as%nsol = 4

      CASE('MADE', 'made')
        optset(i)%as%NS = 0
        optset(i)%as%KS = 1
        optset(i)%as%AS = 2
        optset(i)%as%CS = 3

        DO jm = 1, 3
          optset(i)%as%radmod(jm) = jm
        END DO

        optset(i)%as%nsol = 3

      CASE('MADE3', 'made3')
        optset(i)%as%NS = 0
        optset(i)%as%KS = 1
        optset(i)%as%AS = 4
        optset(i)%as%CS = 7

        DO jm = 1, 3
          optset(i)%as%radmod(jm)   = 1
          optset(i)%as%radmod(jm+3) = 2
          optset(i)%as%radmod(jm+6) = 3
        END DO

        optset(i)%as%nsol = 3

      END SELECT

    END DO
    DO i=1,n_calls + n_inp + n_merge

      IF ((I > n_calls) .AND. (I <= n_calls + n_inp)) CYCLE

      cname=modstr//"_"//optset(i)%as%aero_set_name
      CALL new_channel(status, TRIM(cname), reprid=GP_3D_MID)
      CALL channel_halt(substr, status)

      IF ( optset(i)%as%TANRE ) THEN
        CALL new_channel_object(status, TRIM(cname), 'zaes' &
          , p2=zaes_x, reprid=GP_2D_HORIZONTAL)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'zaes' &
          ,'units', c=' ')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, TRIM(cname), 'zael' &
          , p2=zael_x, reprid=GP_2D_HORIZONTAL)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'zael' &
          ,'units', c=' ')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, TRIM(cname), 'zaeu' &
          , p2=zaeu_x, reprid=GP_2D_HORIZONTAL)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'zaeu' &
          ,'units', c=' ')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, TRIM(cname), 'zaed' &
          , p2=zaed_x, reprid=GP_2D_HORIZONTAL)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'zaed' &
          ,'units', c=' ')
        CALL channel_halt(substr, status)
        
        IF (.NOT. tanre_init) THEN
#ifdef ECHAM5
          CALL init_aero_tanre(dcl%nglon, dcl%nglat,     &
            dcl%glat, dcl%glon,       & 
            coslon, sinlon, gl_twomu, &
            laes, lael, laeu, laed)
          CALL reorder(zaes_x, laes)
          CALL reorder(zael_x, lael)
          CALL reorder(zaeu_x, laeu)
          CALL reorder(zaed_x, laed)
#endif

#ifdef MBM_RAD
          CALL init_aero_tanre(nlon, nlat, &
            bi_glat, bi_glon,              & 
            coslon, sinlon, gl_twomu,      &
            laes, lael, laeu, laed)
          CALL reorder(zaes_x, laes)
          CALL reorder(zael_x, lael)
          CALL reorder(zaeu_x, laeu)
          CALL reorder(zaed_x, laed)
#endif

#ifdef CESM1
          DO zjrow =1,ngpblks
             zkproma = npromz(zjrow)
             CALL init_aero_tanre_dc(zkproma,     &
                  coslon_2d(:,zjrow), sinlon_2d(:,zjrow), sinlat_2d(:,zjrow), &
                  zaes_x(1:zkproma,zjrow), zael_x(1:zkproma,zjrow), &
                  zaeu_x(1:zkproma,zjrow), zaed_x(1:zkproma,zjrow))
          END DO
#endif
        ENDIF
        
      END IF

      ALLOCATE(optset(i)%as%aot_lw_5d(_RI_XYZN_(nproma,ngpblks,nlev,jpband),1))
      ALLOCATE(optset(i)%as%aot_sw_5d(_RI_XYZN_(nproma,ngpblks,nlev,nsw),1))
      ALLOCATE(optset(i)%as%omega_sw_5d(_RI_XYZN_(nproma,ngpblks,nlev,nsw),1))
      ALLOCATE(optset(i)%as%gamma_sw_5d(_RI_XYZN_(nproma,ngpblks,nlev,nsw),1))

      IF (optset(i)%as%lcalc_jval) THEN
         write(*,*) "Jval setup for wavelength"
        ! determine the number of calculations for JVAL
        optset(i)%as%n_jv_bands = n_jv_calc
        DO k=1,n_jv_calc
           IF (l_mean_wavelength) THEN
              ! calculate mean wavelength from the boundaries
              delta = jval_wavelen_u(k) - jval_wavelen_l(k)
              jval_wavelen(k) = jval_wavelen_l(k) + 1./2. * delta
           ELSE  
              ! specific based on landgraf&crutzen
              jval_wavelen(k) = jval_wavelen_lc(k)
           ENDIF
        END DO

        ALLOCATE(optset(i)%as%jv_asca_5d(_RI_XYZN_(nproma,ngpblks,nlev,n_jv),1))
        ALLOCATE(optset(i)%as%jv_aabs_5d(_RI_XYZN_(nproma,ngpblks,nlev,n_jv),1))
        ALLOCATE(optset(i)%as%jv_ga_5d(_RI_XYZN_(nproma,ngpblks,nlev,n_jv),1))
        IF (i < 10) then
          write (str_num,'(A,I1)') "0",i
        ELSE
          write (str_num,'(I2)') i
        END IF
        cname2=modstr//"_jval_"//str_num
        CALL new_channel(status, TRIM(cname2), reprid=GP_3D_MID)
        CALL channel_halt(substr, status)

        p4 => optset(i)%as%jv_asca_5d(:,:,:,:,1)
        CALL new_channel_object(status, TRIM(cname2), 'aer_asca', &
          p4=optset(i)%as%jv_asca, lrestreq=.FALSE., &
          reprid = REPR_AEROPT_4D_JVAL, mem=p4)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname2), 'aer_asca', &
          'long_name', c='aerosol scattering extinction for JVAL (per layer)' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname2), 'aer_asca', &
          'units', c=' - ' )
        CALL channel_halt(substr, status)

        p4 => optset(i)%as%jv_aabs_5d(:,:,:,:,1)
        CALL new_channel_object(status, TRIM(cname2), 'aer_aabs', &
          p4=optset(i)%as%jv_aabs, lrestreq=.FALSE., &
          reprid = REPR_AEROPT_4D_JVAL, mem=p4)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname2), 'aer_aabs', &
          'long_name', c='aerosol absorbing extinction for JVAL (per layer)' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname2), 'aer_aabs', &
          'units', c=' - ' )
        CALL channel_halt(substr, status)

        p4 => optset(i)%as%jv_ga_5d(:,:,:,:,1)
        CALL new_channel_object(status, TRIM(cname2), 'aer_ga', &
          p4=optset(i)%as%jv_ga, lrestreq=.FALSE., &
          reprid = REPR_AEROPT_4D_JVAL, mem=p4)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname2), 'aer_ga', &
          'long_name', c='aerosol asymmetry factor for JVAL (per layer)' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname2), 'aer_ga', &
          'units', c=' - ' )
        CALL channel_halt(substr, status)

        DO ji=1,n_jv
          IF (ji < 10) then
            write (str_num,'(A,I1)') "0",ji
          ELSE
            write (str_num,'(I2)') ji
          END IF
          p4 => optset(i)%as%jv_asca_5d(_RI_XYZN_(:,:,:,ji),:)
          CALL new_channel_object(status, TRIM(cname2), 'aer_asca_b'//str_num, &
            mem=p4, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, TRIM(cname2), 'aer_asca_b'//str_num, &
            'long_name', c='aerosol scattering extinction for JVAL (per layer) for band '//str_num )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, TRIM(cname2), 'aer_asca_b'//str_num, &
            'units', c=' - ' )
          CALL channel_halt(substr, status)
          
          p4 => optset(i)%as%jv_aabs_5d(_RI_XYZN_(:,:,:,ji),:)
          CALL new_channel_object(status, TRIM(cname2), 'aer_aabs_b'//str_num, &
            mem=p4, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, TRIM(cname2), 'aer_aabs_b'//str_num, &
            'long_name', c='aerosol absorbing extinction for JVAL (per layer) for band '//str_num )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, TRIM(cname2), 'aer_aabs_b'//str_num, &
            'units', c=' - ' )
          CALL channel_halt(substr, status)

          p4 => optset(i)%as%jv_ga_5d(_RI_XYZN_(:,:,:,ji),:)
          CALL new_channel_object(status, TRIM(cname2), 'aer_ga_b'//str_num, &
            mem=p4, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, TRIM(cname2), 'aer_ga_b'//str_num, &
            'long_name', c='aerosol asymmetry factor for JVAL (per layer) for band '//str_num )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, TRIM(cname2), 'aer_ga_b'//str_num, &
            'units', c=' - ' )
          CALL channel_halt(substr, status)
        END DO
      
      ELSE
        optset(i)%as%n_jv_bands = 0
      ENDIF
      p4 => optset(i)%as%aot_lw_5d(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(cname), 'aot_lw', &
        p4=optset(i)%as%aot_lw, lrestreq=.FALSE., &
        reprid = REPR_AEROPT_4D_JPBAND, mem=p4)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(cname), 'aot_lw', &
        'long_name', c='longwave aerosol optical thickness' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(cname), 'aot_lw', &
        'units', c=' - ' )
      CALL channel_halt(substr, status)

      DO ji=1,jpband
        IF (ji < 10) then
          write (str_num,'(A,I1)') "0",ji
        ELSE
          write (str_num,'(I2)') ji
        END IF
        p4 => optset(i)%as%aot_lw_5d(_RI_XYZN_(:,:,:,ji),:)
        CALL new_channel_object(status, TRIM(cname), 'aot_lw_B'//str_num, &
          mem=p4, lrestreq=.TRUE.)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'aot_lw_B'//str_num, &
          'long_name', &
          c='longwave aerosol optical thickness for Band '//str_num )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'aot_lw_B'//str_num, &
          'units', c=' - ' )
        CALL channel_halt(substr, status)
      END DO

      p4 => optset(i)%as%aot_sw_5d(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(cname), 'aot_sw', &
        p4=optset(i)%as%aot_sw, lrestreq=.FALSE., &
        reprid = REPR_AEROPT_4D_NSW, mem=p4)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(cname), 'aot_sw', &
        'long_name', c='shortwave aerosol optical thickness' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(cname), 'aot_sw', &
        'units', c=' - ' )
      CALL channel_halt(substr, status)
      
      p4 => optset(i)%as%omega_sw_5d(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(cname), 'omega_sw', &
        p4=optset(i)%as%omega_sw, lrestreq=.FALSE., &
        reprid = REPR_AEROPT_4D_NSW, mem=p4)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(cname), 'omega_sw', &
        'long_name', c='shortwave aerosol single scattering albedo' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(cname), 'omega_sw', &
        'units', c=' - ' )
      CALL channel_halt(substr, status)

      p4 => optset(i)%as%gamma_sw_5d(:,:,:,:,1)
      CALL new_channel_object(status, TRIM(cname), 'gamma_sw', &
        p4=optset(i)%as%gamma_sw, lrestreq=.FALSE., &
        reprid = REPR_AEROPT_4D_NSW, mem=p4)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(cname), 'gamma_sw', &
        'long_name', c='shortwave aerosol asymmetry factor')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, TRIM(cname), 'gamma_sw', &
        'units', c=' - ' )
      CALL channel_halt(substr, status)
      
      DO ji=1,nsw
        IF (ji < 10) then
          write (str_num,'(A,I1)') "0",ji
        ELSE
          write (str_num,'(I2)') ji
        END IF

        p4 => optset(i)%as%aot_sw_5d(_RI_XYZN_(:,:,:,ji),:)
        CALL new_channel_object(status, TRIM(cname), 'aot_sw_B'//str_num, &
          mem=p4, lrestreq=.TRUE.)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'aot_sw_B'//str_num, &
          'long_name', &
          c='shortwave aerosol optical thickness for Band '//str_num )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'aot_sw_B'//str_num, &
          'units', c=' - ' )
        CALL channel_halt(substr, status)

        p4 => optset(i)%as%omega_sw_5d(_RI_XYZN_(:,:,:,ji),:)
        CALL new_channel_object(status, TRIM(cname), 'omega_sw_B'//str_num, &
          mem=p4, lrestreq=.TRUE.)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'omega_sw_B'//str_num, &
          'long_name', &
          c='shortwave aerosol single scattering albedo for Band '//str_num )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'omega_sw_B'//str_num, &
          'units', c=' - ' )
        CALL channel_halt(substr, status)

        p4 => optset(i)%as%gamma_sw_5d(_RI_XYZN_(:,:,:,ji),:)
        CALL new_channel_object(status, TRIM(cname), 'gamma_sw_B'//str_num, &
          mem=p4, lrestreq=.TRUE.)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'gamma_sw_B'//str_num, &
          'long_name', &
          c='shortwave aerosol asymmetry factor for Band '//str_num )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, TRIM(cname), 'gamma_sw_B'//str_num, &
          'units', c=' - ' )
        CALL channel_halt(substr, status)
      END DO

      n = 0
      DO jm=1, max_diag_wavelens
        IF (optset(i)%as%diag_wavelen(jm) <= 0.0_dp) CYCLE
        n = n + 1
      END DO
      optset(i)%as%num_opt_wavelens = n

      IF (optset(i)%as%num_opt_wavelens > 0) THEN

        nmod => optset(i)%as%nmod

        ALLOCATE(optset(i)%as%aot_opt(num_diag_elements,max_diag_wavelens,nmod+1))
        ALLOCATE(optset(i)%as%extcoeff_opt(num_diag_elements,max_diag_wavelens,nmod+1))

        DO jm=1,nmod + 1
          IF (jm == nmod + 1) THEN
            nmo = 'TOT'
          ELSE
            write (nm, '(I2)') INT(jm)
            nm = adjustl(nm)
            IF (jm < 10) THEN
              nmo = 'M0'//TRIM(nm)
            ELSE
              nmo = 'M'//nm
            ENDIF


          ENDIF

          DO n = 1, optset(i)%as%num_opt_wavelens

            IF (optset(i)%as%diag_wavelen(n) < 1.) THEN
              write (wl, '(I4)') INT(optset(i)%as%diag_wavelen(n) * 1000.0 + 0.5)
            ELSE
              write (wl, "(I2, 'um')") INT(optset(i)%as%diag_wavelen(n) + 0.5)
            ENDIF

            wl = adjustl(wl)

            CALL new_channel_object(status,TRIM(cname),           &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_bc',    &
              p3 = optset(i)%as%aot_opt(diag_bc,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_bc',    &
              'long_name',                                        &
              c='AOT (BC), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_bc',    &
              'units', c=' - ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_oc',    &
              p3 = optset(i)%as%aot_opt(diag_oc,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_oc',    &
              'long_name',                                        &
              c='AOT (OC), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_oc',    &
              'units', c=' - ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_du',    &
              p3 = optset(i)%as%aot_opt(diag_du,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_du',    &
              'long_name',                                        &
              c='AOT (DU), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_du',    &
              'units', c=' - ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_waso',  &
              p3 = optset(i)%as%aot_opt(diag_waso,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_waso',  &
              'long_name',                                        &
              c='AOT (WASO), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_waso',  &
              'units', c=' - ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_h2o',   &
              p3 = optset(i)%as%aot_opt(diag_h2o,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_h2o',   &
              'long_name',                                        &
              c='AOT (H2O), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_h2o',   &
              'units', c=' - ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_ss',    &
              p3 = optset(i)%as%aot_opt(diag_ss,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
            'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_ss',      &
            'long_name',                                          &
            c='AOT (SS), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_ss',    &
              'units', c=' - ' )
            CALL channel_halt(substr, status)
          
            CALL new_channel_object(status,TRIM(cname),           &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_total', &
              p3 = optset(i)%as%aot_opt(diag_tot,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_total', &
              'long_name',                                        &
              c='AOT (total), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'aot_opt_' //TRIM(nmo)//'_'// trim(wl) // '_total', &
              'units', c=' - ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_bc',    &
             p3 = optset(i)%as%extcoeff_opt(diag_bc,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_bc',    &
              'long_name',                                        &
              c='ExtCoeff (BC), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
             'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_bc',    &
              'units', c=' 1/m ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_oc',    &
              p3 = optset(i)%as%extcoeff_opt(diag_oc,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_oc',    &
              'long_name',                                        &
              c='ExtCoeff (OC), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_oc',    &
              'units', c=' 1/m ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_du',    &
              p3 = optset(i)%as%extcoeff_opt(diag_du,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_du',    &
              'long_name',                                        &
              c='ExtCoeff (DU), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_du',    &
              'units', c=' 1/m ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_waso',    &
              p3 = optset(i)%as%extcoeff_opt(diag_waso,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_waso',    &
              'long_name',                                        &
              c='ExtCoeff (WASO), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_waso',    &
              'units', c=' 1/m ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_h2o',    &
              p3 = optset(i)%as%extcoeff_opt(diag_h2o,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_h2o',    &
              'long_name',                                        &
              c='ExtCoeff (H2O), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_h2o',    &
              'units', c=' 1/m ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_ss',    &
              p3 = optset(i)%as%extcoeff_opt(diag_ss,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_ss',    &
              'long_name',                                        &
              c='ExtCoeff (SS), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_ss',    &
              'units', c=' 1/m ' )
            CALL channel_halt(substr, status)

            CALL new_channel_object(status,TRIM(cname),           &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_total',    &
              p3 = optset(i)%as%extcoeff_opt(diag_tot,n,jm)%ptr )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_total',    &
              'long_name',                                        &
              c='ExtCoeff (total), '//trim(wl)//' nm ('//trim(nmo)//')')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, TRIM(cname),               &
              'extcoeff_opt_' //TRIM(nmo)//'_'// trim(wl) // '_total',    &
              'units', c=' 1/m ' )
            CALL channel_halt(substr, status)

            do ji=1,num_diag_elements
              optset(i)%as%aot_opt(ji,n,jm)%ptr(:,:,:) = 0._dp ! initialization
              optset(i)%as%extcoeff_opt(ji,n,jm)%ptr(:,:,:) = 0._dp ! initialization
            enddo
          END DO ! n-loop (optional wavelengths)
        END DO
      END IF ! IF (num_opt_wavelens > 0)

    END DO ! set loop


    ! loop for import of climatologies read via an arbitrary channel object
    ! e.g. from import
    DO i=n_calls + 1, n_calls + n_inp
       write (inp_set_nr,'(I2)') i

       ! check for 3D or 4D objects

       SELECT CASE (optset(i)%as%inp_dim)
       CASE(3)
          
          ! DETERMINE number of wavelengthbands for longwave

          optset(i)%as%inp_nlw = 0
          DO j=1,99
             write (nm, '(I2)') INT(j)
             nm = adjustl(nm)
             IF (j < 10) THEN
                nm = '0'//TRIM(nm)
             ENDIF
             CALL get_channel_object(status, TRIM(optset(i)%io%channel_aot_lw), &
               TRIM(optset(i)%io%object_aot_lw)//'_B'//nm, p3=p3)
           !  print*, "object name: ",TRIM(optset(i)%io%channel_aot_lw), &
           !           TRIM(optset(i)%io%object_aot_lw)//'_B'//nm, "toll",status
             IF (status == 0) THEN 
                !CALL info_bi( 'channel / object : '//&
                !     &TRIM(optset(i)%io%channel_aot_lw)//' / '//&
                !     &TRIM(optset(i)%io%object_aot_lw)//'_B'//nm, substr)
                optset(i)%as%inp_nlw = j
             ELSE
                EXIT
             END IF
          END DO
          !print*, " lw input bands ", optset(i)%as%inp_nlw
          ALLOCATE(optset(i)%as%inp_lw(optset(i)%as%inp_nlw))

          ! DETERMINE number of wavelengthbands for shortwave

          optset(i)%as%inp_nsw = 0
          DO j=1,500
             write (nm, '(I2)') INT(j)
             nm = adjustl(nm)
             IF (j < 10) THEN
                nm = '0'//TRIM(nm)
             ENDIF
             CALL get_channel_object(status, TRIM(optset(i)%io%channel_aot_sw), &
               TRIM(optset(i)%io%object_aot_sw)//'_B'//nm, p3=p3)
             IF (status == 0) THEN 
                !CALL info_bi( 'channel / object : '//&
                !     &TRIM(optset(i)%io%channel_aot_sw)//' / '//&
                !     &TRIM(optset(i)%io%object_aot_sw)//'_B'//nm, substr)
                optset(i)%as%inp_nsw = j
             ELSE
                EXIT
             END IF
          END DO
          !print*, " sw input bands ", optset(i)%as%inp_nsw
          ALLOCATE(optset(i)%as%inp_sw(optset(i)%as%inp_nsw))

          ! ASSOCIATE channel objects

          DO j=1,optset(i)%as%inp_nsw
             write (nm, '(I2)') INT(j)
             nm = adjustl(nm)
             IF (j < 10) THEN
                nmo = '0'//TRIM(nm)
             ELSE
                nmo = nm
             ENDIF
             CALL get_channel_object(status, TRIM(optset(i)%io%channel_aot_sw), &
                  TRIM(optset(i)%io%object_aot_sw)//'_B'//nmo, &
                  p3=optset(i)%as%inp_sw(j)%extinct)
             CALL channel_halt(&
                  substr//': channel object for aot_sw not found for set '//inp_set_nr//' Band '//nmo//' !', status)
             CALL get_channel_object(status, TRIM(optset(i)%io%channel_gamma), &
                  TRIM(optset(i)%io%object_gamma)//'_B'//nmo, &
                  p3=optset(i)%as%inp_sw(j)%gamma)
             CALL channel_halt(&
                  substr//': channel object for omega not found for set '//inp_set_nr//' Band '//nmo//' !', status)
             CALL get_channel_object(status, TRIM(optset(i)%io%channel_omega), &
                  TRIM(optset(i)%io%object_omega)//'_B'//nmo, &
                  p3=optset(i)%as%inp_sw(j)%omega)
             CALL channel_halt(&
                  substr//': channel object for gamma not found for set '//inp_set_nr//' Band '//nmo//' !', status)
          END DO
          DO j=1,optset(i)%as%inp_nlw
             write (nm, '(I2)') INT(j)
             nm = adjustl(nm)
             IF (j < 10) THEN
                nmo = '0'//TRIM(nm)
             ELSE
                nmo = nm
             ENDIF
             CALL get_channel_object(status, TRIM(optset(i)%io%channel_aot_lw), &
                  TRIM(optset(i)%io%object_aot_lw)//'_B'//nmo, &
                  p3=optset(i)%as%inp_lw(j)%extinct)
             CALL channel_halt(&
                  substr//': channel object for aot_lw not found for set '//inp_set_nr//' Band '//nmo//' !', status)
          END DO

       CASE(4)
          CALL get_channel_object(status, TRIM(optset(i)%io%channel_aot_lw), &
               TRIM(optset(i)%io%object_aot_lw), p4=optset(i)%as%aot_lw)
          CALL channel_halt(&
               substr//': channel object for aot_lw not found for set '//inp_set_nr//' !', status)
          CALL get_channel_object(status, TRIM(optset(i)%io%channel_aot_sw), &
               TRIM(optset(i)%io%object_aot_sw), p4=optset(i)%as%aot_sw)
          CALL channel_halt(&
               substr//': channel object for aot_sw not found for set '//inp_set_nr//' !', status)
          CALL get_channel_object(status, TRIM(optset(i)%io%channel_gamma), &
               TRIM(optset(i)%io%object_gamma), p4=optset(i)%as%gamma_sw)
          CALL channel_halt(&
               substr//': channel object for gamma not found for set '//inp_set_nr//' !', status)
          CALL get_channel_object(status, TRIM(optset(i)%io%channel_omega), &
               TRIM(optset(i)%io%object_omega), p4=optset(i)%as%omega_sw)
          CALL channel_halt(&
               substr//': channel object for omega not found for set '//inp_set_nr//' !', status)
       END SELECT
    END DO

    CALL info_bi('index searching for merged data sets', substr)
    DO i=n_calls + n_inp + 1,n_calls + n_inp + n_merge
       DO j=1,n_calls + n_inp + n_merge
          IF (optset(i)%io%merge_name1 == optset(j)%as%aero_set_name) &
               optset(i)%as%merge_idx1 = j
          IF (optset(i)%io%merge_name2 == optset(j)%as%aero_set_name) &
               optset(i)%as%merge_idx2 = j
       END DO
       IF (optset(i)%as%merge_idx1 == 0) &
            call error_bi('set 1 for merging not found', substr)
       IF (optset(i)%as%merge_idx2 == 0) &
            call error_bi('set 2 for merging not found', substr)
    END DO

    ! 
    CALL info_bi('creating lookup - tables', substr)
    DO i=1,n_lut      
      CALL aeropt_initialise(p_parallel_IO, i)
    END DO

    CALL info_bi('checking lookup - tables', substr)
    DO i=1,n_calls
      IF (optset(i)%as%TANRE) CYCLE
      CALL aeropt_finalise_lut(optset(i)%as, p_parallel_io, status)
      IF (status /= 0) CALL error_bi(&
        'Lookup Table diagnostic wavelengths entries do not match the requirements!',&
        substr)
      CALL aeropt_check_lut(optset(i)%as, p_parallel_io)
    END DO


    ! coupling to tracers
    ! looped over all aeropt_calls
    DO i=1,n_calls
      ! create exclude-list for species in all modes (tracer - basename)
      call strcrack(optset(i)%as%exclude_string_full, ';', strname_full, n_out_full)
    
    ! create exclude-list for single species (tracer - fullname)
      call strcrack(optset(i)%as%exclude_string_spec, ';', strname_spec, n_out_spec)

      IF (optset(i)%as%TANRE) THEN
        ALLOCATE(optset(i)%as%spec_exclude(naer))
        optset(i)%as%spec_exclude(:) = .FALSE.
        DO jt=1,naer
          write (ichar,'(I1)') jt
          DO ji=1,n_out_full
            IF (TRIM(strname_full(ji)) == TRIM(ichar)) THEN
              optset(i)%as%spec_exclude(jt) = .TRUE.
            ENDIF
          END DO
        END DO
      ENDIF
  
      IF (optset(i)%as%TANRE) CYCLE
          
      nmod => optset(i)%as%nmod
      ALLOCATE(optset(i)%as%aerspec_e5(nmod))
      call map_optset_struct(optset(i)%as)

      H2O_found   = .FALSE.
      DO jm = 1,nmod
        counter = 0
        Do jt=1,ntrac
          no_count = .false.
          IF (ti(jt)%tp%ident%medium/=AEROSOL) CYCLE
          IF (ti(jt)%tp%ident%quantity /= 1) CYCLE
          IF (ti(jt)%tp%meta%cask_i(I_aerosol_mode) /= jm) CYCLE
          IF (ti(jt)%tp%ident%type==FAMILY) CYCLE
               

          DO ji=1,n_out_full
            IF (ti(jt)%tp%ident%basename==TRIM(strname_full(ji))) THEN
              no_count = .true.
              EXIT
            END IF
          ENDDO
          DO ji=1,n_out_spec
            IF (ti(jt)%tp%ident%fullname==TRIM(strname_spec(ji))) THEN
              no_count = .true.
              EXIT
            END IF
          ENDDO
          IF (no_count) CYCLE
          counter = counter + 1

          ! check if H2O tracer exists in this mode
          IF (ti(jt)%tp%ident%basename=="H2O" .or. &
            ti(jt)%tp%ident%basename=="h2o") THEN
            H2O_found   =.true.
          END IF
        ENDDO

        aerspec(jm)%count0 = counter

        ! if no H2O tracer for this mode has been found
        IF (.NOT.H2O_found) THEN
          IF (jm <= nsol) THEN
            aerspec(jm)%l_calc_water = .true.
            counter = counter + 1
          ENDIF
        ENDIF

        ! calculate sea salt from Na+ and Cl-
        IF (lcalc_seasalt) THEN
          Na_found    = .false.
          Cl_found    = .false.
          aerspec(jm)%Na_idx  = 0
          aerspec(jm)%Cl_idx  = 0
          aerspec(jm)%SS_idx  = 0
          aerspec(jm)%Na_idx2 = 0
          aerspec(jm)%Cl_idx2 = 0
          aerspec(jm)%seasalt_idx2 = 0
          aerspec(jm)%l_calc_seasalt = .TRUE.
          DO jt=1,ntrac
            no_count = .false.
            IF (ti(jt)%tp%ident%medium/=AEROSOL) CYCLE
            IF (ti(jt)%tp%ident%quantity /= 1) CYCLE
            IF (ti(jt)%tp%meta%cask_i(I_aerosol_mode) /= jm) CYCLE
            IF (ti(jt)%tp%ident%type==FAMILY) CYCLE
        
            ! check if one of the components is excluded from the calculations
            DO ji=1,n_out_full
              IF (ti(jt)%tp%ident%basename==TRIM(strname_full(ji))) THEN
                no_count = .true.
                EXIT
              END IF
            ENDDO
            DO ji=1,n_out_spec
              IF (ti(jt)%tp%ident%fullname==TRIM(strname_spec(ji))) THEN
                no_count = .true.
                EXIT
              END IF
            ENDDO
            IF (no_count) CYCLE

          ! check if Na or Nap exists in this mode
            IF (ti(jt)%tp%ident%basename=="Na" .or. &
              ti(jt)%tp%ident%basename=="Nap") THEN
              Na_found   =.true.
              aerspec(jm)%Na_idx = jt
            ENDIF
            ! check if Cl or Clm exists in this mode
            IF (ti(jt)%tp%ident%basename=="Cl" .or. &
              ti(jt)%tp%ident%basename=="Clm") THEN
              Cl_found   =.true.
              aerspec(jm)%Cl_idx = jt
            ENDIF
         
            ! check if a seasalt tracer already exists 
            ! and is not excluded from the calculations

            IF (ti(jt)%tp%ident%basename=="SS") &
              aerspec(jm)%ss_idx=jt
          END DO
        
          ! if both compounds are found and no other problem occured, 
          ! increase the counter by one, for an artificial seasalt species
        
          IF ( (aerspec(jm)%Cl_idx /= 0) .AND. (aerspec(jm)%Na_idx /= 0) &
            .AND. (aerspec(jm)%ss_idx == 0) ) &
            counter = counter + 1
        END IF


        ! define number of considered species
        aerspec(jm)%count = counter
      
        ! save tracer number index for each mode
        Do jt=1,ntrac
          IF (ti(jt)%tp%ident%medium/=AEROSOL) CYCLE
          IF (ti(jt)%tp%meta%cask_i(I_aerosol_mode) /= jm) CYCLE
          IF (ti(jt)%tp%ident%quantity == 2) THEN
            aerspec(jm)%idx_num = jt
            EXIT
          ENDIF
        ENDDO


        ! allocate size of working array for considered compounds
        ALLOCATE(optset(i)%as%aerspec_e5(jm)%paerml(nproma,nlev,aerspec(jm)%count))
        !               of tracer indices for working array
        ALLOCATE(optset(i)%as%aerspec_e5(jm)%idx_trac(aerspec(jm)%count))
        !               of radiation indices for working array
        ALLOCATE(optset(i)%as%aerspec_e5(jm)%rad_spec(aerspec(jm)%count))
        !               of molarmass for working array
        ALLOCATE(optset(i)%as%aerspec_e5(jm)%molarmass(aerspec(jm)%count))
        !               of aerosol densities for working array
        ALLOCATE(optset(i)%as%aerspec_e5(jm)%density(aerspec(jm)%count))
        !               of compound names in the working array
        ALLOCATE(optset(i)%as%aerspec_e5(jm)%name(aerspec(jm)%count))

        aerspec(jm)%ram2cmr = 1._dp / EXP(1.5_dp*(LOG(sigma(jm))**2))

        counter2 = 0
        DO jt=1,ntrac
          no_count = .false.
          IF (ti(jt)%tp%ident%medium/=AEROSOL) CYCLE
          IF (ti(jt)%tp%ident%quantity /= 1) CYCLE
          IF (ti(jt)%tp%meta%cask_i(I_aerosol_mode) /= jm) CYCLE
          IF (ti(jt)%tp%ident%type==FAMILY) CYCLE
        
          ! check if one of the components is excluded from the calculations
          DO ji=1,n_out_full
            IF (ti(jt)%tp%ident%basename==TRIM(strname_full(ji))) THEN
              no_count = .true.
              EXIT
            END IF
          ENDDO
          DO ji=1,n_out_spec
            IF (ti(jt)%tp%ident%fullname==TRIM(strname_spec(ji))) THEN
              no_count = .true.
              EXIT
            END IF
          ENDDO
          IF (no_count) CYCLE
          counter2 = counter2 + 1
          aerspec(jm)%idx_trac(counter2)  = jt
          aerspec(jm)%molarmass(counter2) = ti(jt)%tp%meta%cask_r(R_molarmass)
          aerspec(jm)%name(counter2)      = TRIM(ti(jt)%tp%ident%basename)
          aerspec(jm)%density(counter2)   = &
            ti(jt)%tp%meta%cask_r(R_AEROSOL_DENSITY)
        END DO
        IF (.NOT.H2O_found) THEN
          IF (jm <= nsol) THEN
            counter2 = counter2 + 1
            aerspec(jm)%idx_trac(counter2)  = 0
            aerspec(jm)%molarmass(counter2) = M_H2O
            aerspec(jm)%name(counter2)      = "H2O"
            aerspec(jm)%density(counter2)   = rho_h2o
            aerspec(jm)%h2o_idx = counter2
          ENDIF
        ENDIF

        IF (lcalc_seasalt) THEN
          IF ( (aerspec(jm)%Cl_idx == 0) .OR. (aerspec(jm)%Na_idx == 0) ) &
            aerspec(jm)%l_calc_seasalt = .false.          
          IF ( (aerspec(jm)%l_calc_seasalt) .AND. &
            (aerspec(jm)%ss_idx == 0) ) THEN
            counter2 = counter2 + 1
            aerspec(jm)%idx_trac(counter2)  = 0
            aerspec(jm)%molarmass(counter2) = M_Na + M_Cl
            aerspec(jm)%name(counter2)      = "SS"
            aerspec(jm)%density(counter2)   = 2.2e3_dp
            aerspec(jm)%seasalt_idx2 = counter2
          ENDIF
          DO jt=1,aerspec(jm)%count
            IF (aerspec(jm)%name(jt) == "Na" .OR. &
              aerspec(jm)%name(jt) == "Nap") aerspec(jm)%Na_idx2 = jt
            IF (aerspec(jm)%name(jt) == "Cl" .OR. &
              aerspec(jm)%name(jt) == "Clm") aerspec(jm)%Cl_idx2 = jt
            IF (aerspec(jm)%name(jt) == "SS" ) aerspec(jm)%seasalt_idx2 = jt
          ENDDO
        ELSE
          aerspec(jm)%l_calc_seasalt = .false.
        ENDIF

        DO jt=1,aerspec(jm)%count
          SELECT CASE  ( TRIM(aerspec(jm)%name(jt)) )
          CASE ('SS', 'ss', 'Ss')
            aerspec(jm)%rad_spec(jt) = ri_ss
          CASE ('DU', 'du', 'Du')
            aerspec(jm)%rad_spec(jt) = ri_du
          CASE ('OC', 'oc', 'Oc', 'POM', 'POM2', 'POMphob', 'POMphil', 'SOA' )
            aerspec(jm)%rad_spec(jt) = ri_oc
          CASE ('BC', 'bc', 'Bc', 'EC', 'ec', 'Ec', 'BC2', 'BCphob', 'BCphil', 'BCtag')
            aerspec(jm)%rad_spec(jt) = ri_bc
          CASE ('H2O', 'h2o')
            aerspec(jm)%rad_spec(jt) = ri_h2o
          CASE default
            aerspec(jm)%rad_spec(jt) = ri_waso
            
          END SELECT
          IF (MATCH_WILD( 'WSOC*', TRIM(aerspec(jm)%name(jt)) )) THEN
            aerspec(jm)%rad_spec(jt) = ri_oc 
          ENDIF

! needed for oracle SOA produced
          IF (MATCH_WILD( 'fPOA*', TRIM(aerspec(jm)%name(jt)) )) THEN
            aerspec(jm)%rad_spec(jt) = ri_oc 
          ENDIF
          IF (MATCH_WILD( 'fSOA*', TRIM(aerspec(jm)%name(jt)) )) THEN
            aerspec(jm)%rad_spec(jt) = ri_oc 
          ENDIF
          IF (MATCH_WILD( 'bbPOA*', TRIM(aerspec(jm)%name(jt)) )) THEN
            aerspec(jm)%rad_spec(jt) = ri_oc 
          ENDIF
          IF (MATCH_WILD( 'bbSOA*', TRIM(aerspec(jm)%name(jt)) )) THEN
            aerspec(jm)%rad_spec(jt) = ri_oc 
          ENDIF
          IF (MATCH_WILD( 'SOAv*', TRIM(aerspec(jm)%name(jt)) )) THEN
            aerspec(jm)%rad_spec(jt) = ri_oc 
          ENDIF
        ENDDO
      END DO

    


    END DO ! loop of aeropt calls




    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

  END SUBROUTINE AEROPT_INIT_COUPLING

!===============================================================================

  SUBROUTINE AEROPT_RADIATION

    CALL AEROPT_DRIVER
    CALL AEROPT_MERGE

  END SUBROUTINE AEROPT_RADIATION

!===============================================================================

  SUBROUTINE AEROPT_DRIVER

    USE messy_main_grid_def_mem_bi, ONLY: nlev, nlevp1, kproma, jrow
    USE messy_main_grid_def_bi,     ONLY: grmass, grvol   
    USE messy_main_data_bi,         ONLY: aphm1, apm1, tm1,              &
                                          slm, icecov, seacov,           &
                                          tslm1, tsi, tsw      
    USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte,   pxtm1 => qxtm1

    USE messy_main_constants_mem,  ONLY: M_air, M_h2o, rho_h2o, pi

    USE messy_aeropt_mem,          ONLY: nmod, aerspec, nsw, jpband,    &
                                         num_diag_elements,             &
                                         max_diag_wavelens,             &
                                         num_opt_wavelens
    USE messy_aeropt,              ONLY: calc_aero_properties, min_val, &
                                         min_val_omega,                 &
                                         vrev2
    USE messy_aeropt_tanre,        ONLY: naer, calc_aero_tanre,         &
                                         lw_optics_tanre, sw_optics_tanre
    USE messy_aeropt_sets,         ONLY: map_optset_struct, map_lut_struct
    USE messy_main_timer,          ONLY: time_step_len

    REAL(dp), DIMENSION(:,:,:),     POINTER :: prwet    => NULL()
    REAL(dp), DIMENSION(:,:,:),     POINTER :: paernl   => NULL()
    REAL(dp), DIMENSION(:,:,:),     POINTER :: paernl2  => NULL()
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: zaot_opt => NULL()
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: zextcoeff_opt => NULL()    
    
    REAL(dp) :: zaot_sw(kproma,nlev,nsw)
    REAL(dp) :: zaot_lw(kproma,nlev,jpband)
    REAL(dp) :: zomega_sw(kproma,nlev,nsw)
    REAL(dp) :: zgamma_sw(kproma,nlev,nsw)
    
    REAL(dp), DIMENSION(kproma,nlev,n_jv)   :: z_jv_asca, z_jv_aabs, z_jv_ga

    REAL(dp) :: lrhoa(kproma,nlev)

   ! loop indices
    INTEGER  ::  jl, jk, jm, jt, jj, ji, i, j
    ! auxilliary indices
    INTEGER  :: idx1
    REAL(dp) :: mass_w, mol_wat, volume, delta
    INTEGER  :: isw, ipband ! op_pj_20130429

    ! "aerosol" field
    REAL(dp), DIMENSION(:,:,:), POINTER :: zsaer => NULL()
    ! layer thickness
    REAL(dp), DIMENSION(:,:),   POINTER :: zdp   => NULL()
    ! air temperature at full levels
    REAL(DP), DIMENSION(:,:),   POINTER :: ptf   => NULL()
    ! air temperature at half levels
    REAL(dp), DIMENSION(:,:),   POINTER :: zth   => NULL()

    DO jk= 1,nlev
      DO jl = 1,kproma
        lrhoa(jl,jk) = grmass(_RI_XYZ__(jl,jrow,jk)) / grvol(_RI_XYZ__(jl,jrow,jk))
      ENDDO
    ENDDO

    DO i=1,n_calls
      CALL MAP_OPTSET_STRUCT(optset(i)%as)
      IF ( optset(i)%as%TANRE ) THEN

        ALLOCATE(zsaer(kproma,nlev,naer))
        ALLOCATE(zdp(kproma,nlev))
        ALLOCATE(ptf(kproma,nlev))
        ALLOCATE(zth(kproma,nlevp1))

        ptf(1:kproma,1:nlev) = tm1(_RI_XYZ__(1:kproma,jrow,1:nlev))

        ! layer pressure thickness
        zdp(1:kproma,:)=aphm1(1:kproma,2:nlev+1)-aphm1(1:kproma,1:nlev)
        ! temperature at half levels
        DO jk=2,nlev
          DO jl = 1, kproma
            zth(jl,jk) = ( ptf(jl,jk-1)*apm1(jl,jk-1)*(apm1(jl,jk)-aphm1(jl,jk)  ) &
                       + ptf(jl,jk) * apm1(jl,jk) * (aphm1(jl,jk)-apm1(jl,jk-1)) ) &
                       / (aphm1(jl,jk) * (apm1(jl,jk)-apm1(jl,jk-1)) )
          END DO
        END DO
        DO jl = 1, kproma
          zth(jl,nlevp1) = (slm(jl,jrow)    * tslm1(jl,jrow)**4._dp     &
                         +  icecov(jl,jrow) * tsi(jl,jrow)**4._dp       &
                         +  seacov(jl,jrow) * tsw(jl,jrow)**4._dp )**0.25_dp
          zth(jl,1) = ptf(jl,1)-apm1(jl,1)*(ptf(jl,1)-zth(jl,2)) &
                    / (apm1(jl,1)-aphm1(jl,2))
        END DO
      
        zsaer(1:kproma,:,:)=EPSILON(1._dp)
        ! function for the calculation of the "aerosol" field
        zsaer(1:kproma,:,1:naer)=calc_aero_tanre(jrow,kproma,nlev,zdp &
             ,apm1(1:kproma,:),zth)
 
        ! CHECK for exclusion of species from the TANRE climatology, e.g.
        ! 4 = volcanic aerosol (eliminated via the namelist exclude_string)
        DO j=1,naer
          IF ( optset(i)%as%SPEC_EXCLUDE(j) ) zsaer(1:kproma,:,j) = 0._dp
        END DO


        ! aerosol values are at least epsilon
        zsaer(1:kproma,:,:)=MAX(zsaer(1:kproma,:,:),EPSILON(1._dp))

        CALL lw_optics_tanre(aot_lw(_RI_XYZN_(1:kproma,jrow,:,1:jpband)), nlev, kproma, &
                             zsaer(1:kproma,:,:) )
        CALL sw_optics_tanre(aot_sw(_RI_XYZN_(1:kproma,jrow,:,1:nsw)),  &
                             omega_sw(_RI_XYZN_(1:kproma,jrow,:,1:nsw)), &
                             gamma_sw(_RI_XYZN_(1:kproma,jrow,:,1:nsw)), nlev, kproma,  &
                             zsaer(1:kproma,:,:) )
        DEALLOCATE(zsaer); NULLIFY(ZSAER)
        DEALLOCATE(zdp);   NULLIFY(ZDP)
        DEALLOCATE(ptf);   NULLIFY(PTF)
        DEALLOCATE(zth);   NULLIFY(ZTH)

        IF (lcalc_jval) THEN
          DO j=1,n_jv
            ! scattering extinction
            jv_asca(_RI_XYZN_(1:kproma,jrow,:,j)) = &
              aot_sw(_RI_XYZN_(1:kproma,jrow,:,1)) * omega_sw(_RI_XYZN_(1:kproma,jrow,:,1))
            ! absorbing extinction
            jv_aabs(_RI_XYZN_(1:kproma,jrow,:,j)) = &
              aot_sw(_RI_XYZN_(1:kproma,jrow,:,1)) * (1._dp - omega_sw(_RI_XYZN_(1:kproma,jrow,:,1)))
            ! asymmetry factor
            jv_ga(_RI_XYZN_(1:kproma,jrow,:,j)) = gamma_sw(_RI_XYZN_(1:kproma,jrow,:,1))
          END DO
        END IF

        ! FLIP VERTICAL DIRECTION INTO NATURAL GP ODER (TOP to BOTTOM)
        DO isw=1, nsw
           aot_sw(_RI_XYZN_(1:kproma,jrow,:,isw)) = &
                vrev2(aot_sw(_RI_XYZN_(1:kproma,jrow,:,isw)))
           omega_sw(_RI_XYZN_(1:kproma,jrow,:,isw)) = &
                vrev2(omega_sw(_RI_XYZN_(1:kproma,jrow,:,isw)))
           gamma_sw(_RI_XYZN_(1:kproma,jrow,:,isw)) = &
                vrev2(gamma_sw(_RI_XYZN_(1:kproma,jrow,:,isw)))
        END DO
        DO ipband = 1, jpband
           aot_lw(_RI_XYZN_(1:kproma,jrow,:,ipband)) = &
                vrev2(aot_lw(_RI_XYZN_(1:kproma,jrow,:,ipband)))
        END DO

      ELSE

        ALLOCATE(prwet(kproma,nlev,nmod))
        ALLOCATE(paernl(kproma,nlev,nmod))
        ALLOCATE(paernl2(kproma,nlev,nmod))
        ALLOCATE(zaot_opt(kproma,nlev,num_diag_elements,max_diag_wavelens,nmod+1))
        ALLOCATE(zextcoeff_opt(kproma,nlev,num_diag_elements,max_diag_wavelens,nmod+1))

        CALL MAP_LUT_STRUCT(lut_number)

        DO jm=1,nmod
          aerspec(jm)%paerml(:,:,:) = 0.0_dp
          DO jt=1,aerspec(jm)%count0
            IDX1 = aerspec(jm)%idx_trac(jt)

            DO jk=1,nlev
              DO jl=1,kproma
            ! unit conversion to ??mol / m^3 
                aerspec(jm)%paerml(jl,jk,jt) = MAX(0._dp, 1.e9_dp / M_air * &
                                               lrhoa(jl,jk) *               &
                                             ( pxtm1(jl,jk,IDX1) +          &
                                               pxtte(jl,jk,IDX1) * time_step_len) )
              ENDDO
            ENDDO
          ENDDO
          IDX1 = aerspec(jm)%idx_num
          DO jk=1,nlev
            DO jl=1,kproma
              ! number in 1/mol
              paernl2(jl,jk,jm) = MAX(0._dp, ( pxtm1(jl,jk,IDX1) +      &
                                  pxtte(jl,jk,IDX1) * time_step_len ) )
              ! number in 1/cm^3
              paernl(jl,jk,jm)  = paernl2(jl,jk,jm) * lrhoa(jl,jk) *    &
                                  1.e-3_dp / M_air
              ! wet radius in cm
              prwet(jl,jk,jm)   = MAX(0._dp, wetradius(_RI_XYZN_(jl,jrow,jk,jm)) * 1.e2_dp)
            ENDDO
          ENDDO
      
          ! special treatment for seasalt, being calculated from NaCl
          IF (aerspec(jm)%l_calc_seasalt) THEN
            DO jk=1,nlev
              DO jl=1,kproma
                delta = MIN(aerspec(jm)%paerml(jl,jk,aerspec(jm)%Cl_idx2),   &
                  aerspec(jm)%paerml(jl,jk,aerspec(jm)%Na_idx2) )
                
                aerspec(jm)%paerml(jl,jk,aerspec(jm)%seasalt_idx2) =   &
                  aerspec(jm)%paerml(jl,jk,aerspec(jm)%seasalt_idx2) + delta
                aerspec(jm)%paerml(jl,jk,aerspec(jm)%Cl_idx2) =        &
                  aerspec(jm)%paerml(jl,jk,aerspec(jm)%Cl_idx2)      - delta
                aerspec(jm)%paerml(jl,jk,aerspec(jm)%Na_idx2) =        &
                  aerspec(jm)%paerml(jl,jk,aerspec(jm)%Na_idx2)      - delta
              ENDDO
            ENDDO
          END IF
          
          ! special treatment for aerosol water, being calculated from 
          ! differences of dry and wet volume
          IF (aerspec(jm)%l_calc_water) THEN
            DO jk=1,nlev
              DO jl=1,kproma
                ! volume of water of a single particle [m^3]
                volume = 4._dp/3._dp * pi * (wetradius(_RI_XYZN_(jl,jrow,jk,jm))**3._dp - &
                                             dryradius(_RI_XYZN_(jl,jrow,jk,jm))**3._dp) 
                ! mass of water of a single particle [kg]
                ! rho_H2O = 999.97_dp    ! density of H2O [kg/m3]
                mass_w = volume * rho_H2O
                ! moles of water of a single particle [??mol]
                mol_wat = mass_w * 1.e9_dp / M_H2O
                ! total moles of water of all particles [??mol/m^3]
                ! water = mol_wat [??mol] * paernl [1/cm^3] * 1e6
                aerspec(jm)%paerml(jl,jk,aerspec(jm)%h2o_idx) = MAX(0._dp, &
                  mol_wat * paernl(jl,jk,jm) * 1.e6_dp)
                
              END DO
            END DO
          ENDIF
          
        ENDDO
        zaot_opt(:,:,:,:,:) = 0._dp
        zextcoeff_opt(:,:,:,:,:) = 0._dp
        zaot_sw(:,:,:)      = 0._dp
        zaot_lw(:,:,:)      = 0._dp
        zomega_sw(:,:,:)    = 0._dp
        zgamma_sw(:,:,:)    = 0._dp

        z_jv_asca(:,:,:)    = 0._dp
        z_jv_aabs(:,:,:)    = 0._dp
        z_jv_ga(:,:,:)      = 0._dp
        
        CALL calc_aero_properties(kproma, nlev, paernl, paernl2, prwet,       &
                                  aphm1(1:kproma,:),                          &
                                  zaot_opt,                                   &
                                  zextcoeff_opt,                              &
                                  zaot_sw, zaot_lw,                           &
                                  zomega_sw, zgamma_sw,                       &
                                  z_jv_asca, z_jv_aabs, z_jv_ga  )

        DO jk=1,nlev
          DO jl=1,kproma
            aot_sw(_RI_XYZN_(jl,jrow,jk,1:nsw))    = MAX(zaot_sw(jl,jk,1:nsw),MIN_VAL)
            aot_lw(_RI_XYZN_(jl,jrow,jk,1:jpband)) = MAX(zaot_lw(jl,jk,1:jpband), MIN_VAL)
            omega_sw(_RI_XYZN_(jl,jrow,jk,1:nsw))  = MIN(zomega_sw(jl,jk,1:nsw),     &
                                              MIN_VAL_OMEGA)
            gamma_sw(_RI_XYZN_(jl,jrow,jk,1:nsw))  = MAX(zgamma_sw(jl,jk,1:nsw), MIN_VAL)

            do ji=1,num_diag_elements
              do jj = 1,num_opt_wavelens
                do jm=1,nmod+1
                  aot_opt(ji,jj,jm)%ptr(_RI_XYZ__(jl,jrow,jk)) = &
                    zaot_opt(jl,jk,ji,jj,jm)
                  extcoeff_opt(ji,jj,jm)%ptr(_RI_XYZ__(jl,jrow,jk)) = &
                    zextcoeff_opt(jl,jk,ji,jj,jm)  
                end do
              enddo
            end do
            
          END DO
        END DO
        
        IF (n_jv_bands > 0) THEN
          DO jk=1,nlev
            DO jl=1,kproma
              jv_asca(_RI_XYZN_(jl,jrow,jk,1:n_jv))  = MAX(z_jv_asca(jl,jk,1:n_jv), MIN_VAL)
              jv_aabs(_RI_XYZN_(jl,jrow,jk,1:n_jv))  = MAX(z_jv_aabs(jl,jk,1:n_jv), MIN_VAL)
              jv_ga(_RI_XYZN_(jl,jrow,jk,1:n_jv))    = MAX(z_jv_ga(jl,jk,1:n_jv), MIN_VAL)
            END DO
          END DO
        END IF

        DEALLOCATE(prwet);    NULLIFY(prwet)
        DEALLOCATE(paernl);   NULLIFY(paernl)
        DEALLOCATE(paernl2);  NULLIFY(paernl2)
        DEALLOCATE(zaot_opt); NULLIFY(zaot_opt)
        DEALLOCATE(zextcoeff_opt); NULLIFY(zextcoeff_opt) 

      END IF! TANRE OR CALCULATION
    END DO ! aeropt_calls

    
  END SUBROUTINE AEROPT_DRIVER

!===============================================================================

  SUBROUTINE AEROPT_MERGE
    
    USE messy_main_grid_def_mem_bi, ONLY: nlev, kproma, jrow
    USE messy_main_grid_def_bi,     ONLY: deltaz
    USE messy_main_data_bi,         ONLY: apm1

    USE messy_aeropt,              ONLY: min_val, min_val_omega
    USE messy_main_constants_mem,  ONLY: g

    IMPLICIT NONE
    INTRINSIC :: MIN, MAX
    
    INTEGER  :: i, jl, jk, ji, jidx, idx
    REAL(dp) :: dummy1_sw(kproma,nsw,2), dummy2_sw(kproma,nsw,2), &
                dummy3_sw(kproma,nsw,2)
    REAL(dp) :: dummy_lw(kproma,jpband,2)

    INTEGER  :: kidx(kproma,2), delta_kidx(kproma)
    REAL(dp) :: k_weight(kproma,nlev,2)
    REAL(dp) :: zheight(kproma,nlev)

    ! If there is any set which needs conversion from 1/m to 1/box
    ! the conversion factors are calculated here
    IF ( lconv_req ) THEN
       zheight(:,:) = 0._dp
       DO jk=2,nlev
          zheight(1:kproma,jk) =deltaz(_RI_XYZ__(1:kproma,jrow,jk))
       END DO
    ENDIF
      
    merge_loop: DO i= n_calls + n_inp + 1, n_calls + n_inp + n_merge


       ! DETERMINE LEVEL INDICES OF PRESSURE BOUNDARIES (DEPENDENT ON
       ! HORIZONTAL POSITION)
       kidx(:,:) = 1
       !
       vecor_loop: DO jl=1,kproma
          ! NOTE: %press_bound1 >= %press_bound2 assumed 
          !       (bottom to top)
          DO jk=1,nlev   
             kidx(jl,1) = jk
             IF (apm1(jl,jk) > optset(i)%as%press_bound1) EXIT
          END DO

          DO jk=1,nlev   
             kidx(jl,2) = jk
             IF (apm1(jl,jk) > optset(i)%as%press_bound2) EXIT
          END DO
          
          DELTA_kidx(jl) = ABS(kidx(jl,1) - kidx(jl,2))
          
          k_weight(jl,:,:) = 1._dp
          
          IF (kidx(jl,1) == kidx(jl,2)) THEN
             k_weight(jl,kidx(jl,1):nlev,1) = 1.0_dp
             IF (kidx(jl,1) > 1) THEN
                k_weight(jl,1:kidx(jl,1)-1,1)  = 0.0_dp
                k_weight(jl,1:kidx(jl,1)-1,2)  = 1.0_dp
             END IF
             k_weight(jl,kidx(jl,1):nlev,2) = 0.0_dp
          ELSE
             DO jk=1,nlev
                IF (jk >= kidx(jl,1)) THEN
                   k_weight(jl,jk,1) = 1.0_dp
                   k_weight(jl,jk,2) = 0.0_dp
                ELSE IF (jk < kidx(jl,2)) THEN
                   k_weight(jl,jk,1) = 0.0_dp
                   k_weight(jl,jk,2) = 1.0_dp
                ELSE
                   k_weight(jl,jk,1) =  (REAL(jk,dp) - &
                        REAL(kidx(jl,2),dp))/REAL(DELTA_kidx(jl),dp)
                   k_weight(jl,jk,2) = 1._dp - k_weight(jl,jk,1)
                END IF
             END DO
          END IF

       END DO vecor_loop

       ! APPLY WEIGHTS FOR MERGING BOTH COMPOENENTS
       level_loop: DO jk=1,nlev

          component_loop: DO jidx =1, 2

             ! INDICES OF COMPONENTS TO BE MERGED
             IF (jidx == 1) idx = optset(i)%as%merge_idx1
             IF (jidx == 2) idx = optset(i)%as%merge_idx2

             SELECT CASE(OPTSET(idx)%AS%INP_DIM)
             CASE(4)
                DO jl=1,kproma

                   dummy1_sw(jl,:,jidx) = optset(idx)%as%aot_sw(_RI_XYZN_(jl,jrow,jk,:))   &
                        * k_weight(jl,jk,jidx)
                   dummy2_sw(jl,:,jidx) = optset(idx)%as%omega_sw(_RI_XYZN_(jl,jrow,jk,:)) &
                        * k_weight(jl,jk,jidx)
                   dummy3_sw(jl,:,jidx) = optset(idx)%as%gamma_sw(_RI_XYZN_(jl,jrow,jk,:)) &
                        * k_weight(jl,jk,jidx)
                   dummy_lw(jl,:,jidx) = optset(idx)%as%aot_lw(_RI_XYZN_(jl,jrow,jk,:)) &
                     * k_weight(jl,jk,jidx)
                END DO
             CASE(3)
                DO ji = 1,nsw
                   DO jl=1,kproma
                      dummy1_sw(jl,ji,jidx) = &
                           optset(idx)%as%inp_sw(ji)%extinct(_RI_XYZ__(jl,jrow,jk)) &
                           * k_weight(jl,jk,jidx)
                      dummy2_sw(jl,ji,jidx) = &
                           optset(idx)%as%inp_sw(ji)%omega(_RI_XYZ__(jl,jrow,jk)) &
                           * k_weight(jl,jk,jidx)
                      dummy3_sw(jl,ji,jidx) = &
                           optset(idx)%as%inp_sw(ji)%gamma(_RI_XYZ__(jl,jrow,jk)) &
                           * k_weight(jl,jk,jidx)
                   END DO
                END DO
                DO ji = 1,jpband
                   DO jl=1,kproma
                      dummy_lw(jl,ji,jidx) = &
                           optset(idx)%as%inp_lw(ji)%extinct(_RI_XYZ__(jl,jrow,jk)) &
                           * k_weight(jl,jk,jidx)
                   END DO
                END DO
             END SELECT

             ! conversion from 1/m to 1/box
             IF (optset(idx)%as%lconv_sw) THEN
                DO jl=1,kproma
                   dummy1_sw(jl,:,jidx) = dummy1_sw(jl,:,jidx) * zheight(jl,jk)
                END DO
             ENDIF
             IF (optset(idx)%as%lconv_lw) THEN
                DO jl=1,kproma
                   dummy_lw(jl,:,jidx) = dummy_lw(jl,:,jidx) * zheight(jl,jk)
                END DO
             END IF
             
          END DO component_loop

          ! APPLY ADDITIONAL WEIGHTS FROM NAMELIST
          DO jl=1,kproma 
             optset(i)%as%aot_sw(_RI_XYZN_(jl,jrow,jk,:)) = MAX(min_val,         &
                  dummy1_sw(jl,:,1) * optset(i)%as%weight1 +          &
                  dummy1_sw(jl,:,2) * optset(i)%as%weight2 )
             optset(i)%as%omega_sw(_RI_XYZN_(jl,jrow,jk,:)) = MIN(min_val_omega, &
                  dummy2_sw(jl,:,1) * optset(i)%as%weight1 +          &
                  dummy2_sw(jl,:,2) * optset(i)%as%weight2 )
             optset(i)%as%gamma_sw(_RI_XYZN_(jl,jrow,jk,:)) = MAX(min_val,       &
                  dummy3_sw(jl,:,1) * optset(i)%as%weight1 +          &
                  dummy3_sw(jl,:,2) * optset(i)%as%weight2 )
             optset(i)%as%aot_lw(_RI_XYZN_(jl,jrow,jk,:)) = MAX(min_val,         &
                  dummy_lw(jl,:,1) * optset(i)%as%weight1 +           &
                  dummy_lw(jl,:,2) * optset(i)%as%weight2 )
          END DO

       END DO level_loop
       
    END DO merge_loop

  END SUBROUTINE AEROPT_MERGE

!===============================================================================


  SUBROUTINE AEROPT_FREE_MEMORY

    USE messy_aeropt_tanre, ONLY: cleanup_aero_tanre

    INTEGER :: i

    do i=1,n_calls
      IF (ASSOCIATED(optset(i)%as%aot_opt)) DEALLOCATE(optset(i)%as%aot_opt)
      IF (ASSOCIATED(optset(i)%as%extcoeff_opt)) DEALLOCATE(optset(i)%as%extcoeff_opt)
    END do
    CALL cleanup_aero_tanre
    DO i = 1,SIZE(optset)
       DEALLOCATE(optset(i)%as)
       NULLIFY(optset(i)%as)
    END DO
  
  END SUBROUTINE AEROPT_FREE_MEMORY


!===============================================================================

#endif

END MODULE MESSY_AEROPT_SI
