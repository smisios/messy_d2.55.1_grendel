!**********************************************************************
MODULE messy_clamsrdfrc_si

#if defined(ECHAM5)
!**********************************************************************
!  Submodel interface for clamsrdfrc (CLaMS radiative forcing in EMAC)
!  In this version, only H2O forcing is completely implemented!!! 
!**********************************************************************
 
  USE messy_clamsrdfrc
  USE messy_main_timer_event, ONLY: time_event, io_time_event 

  IMPLICIT NONE
  INTRINSIC :: NULL

  ! Global variables in RDFRC:
! op_pj_20160714+
!!$  TYPE(time_event)   :: rdfrcevent
!!$  TYPE(io_time_event):: io_rdfrcevent
  TYPE(time_event), SAVE    :: rdfrcevent
  TYPE(io_time_event), SAVE :: io_rdfrcevent
! op_pj_20160714-

  ! MODULE VARIABLES:
  REAL(PREC), DIMENSION(:), POINTER :: LAT_D        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_D        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_D        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: O3_CLAMS_D   => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: H2O_CLAMS_D  => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LAT_G        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_G        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_G        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: O3_CLAMS_G   => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: H2O_CLAMS_G  => NULL()

  REAL(PREC), DIMENSION(:,:,:), POINTER :: E5_ZETA_D => NULL()
  REAL(PREC), DIMENSION(:,:,:), POINTER :: E5_ZETA_G => NULL()

  REAL(PREC), DIMENSION(:,:,:), POINTER :: O3_ECHAM_D  => NULL()
  REAL(PREC), DIMENSION(:,:,:), POINTER :: O3_ECHAM_G  => NULL()
  REAL(PREC), DIMENSION(:,:,:), POINTER :: H2O_ECHAM_D => NULL()
  REAL(PREC), DIMENSION(:,:,:), POINTER :: H2O_ECHAM_G => NULL()

  ! Part of ECHAM fields for interpolation
  REAL(PREC), DIMENSION(:,:,:), ALLOCATABLE :: O3_ECHAM_P  
  REAL(PREC), DIMENSION(:,:,:), ALLOCATABLE :: H2O_ECHAM_P 
  REAL(PREC), DIMENSION(:,:,:), POINTER :: ZETA_ECHAM_P => NULL()

  PUBLIC :: clamsrdfrc_initialize
  PUBLIC :: clamsrdfrc_init_memory
  PUBLIC :: clamsrdfrc_init_coupling
  PUBLIC :: clamsrdfrc_global_end

!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsrdfrc_initialize

    USE messy_main_tools,        ONLY: find_next_free_unit
    USE messy_main_timer,        ONLY: delta_time
    USE messy_main_timer_bi,     ONLY: timer_event_init
    USE messy_clamsrdfrc_tools,  ONLY: nc_read_ap_s_info
    USE messy_clams_global,      ONLY: initfile
    USE messy_clamsrdfrc_global, ONLY: lev_window 

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsrdfrc_initialize'
    INTEGER :: status, iou

    ! Read namelist variables:
    iou = find_next_free_unit(100,200)

    ! Read namelist and set default values:
    CALL clamsrdfrc_read_nml(status, iou)

    CALL nc_read_ap_s_info(initfile, lev_window)

    ! Define forcing time event:
    io_rdfrcevent%counter = timestep_rdfrc   
    io_rdfrcevent%unit = 'hours'
    io_rdfrcevent%adjustment = 'exact'
    io_rdfrcevent%offset = -delta_time
    CALL timer_event_init (rdfrcevent, io_rdfrcevent, 'RDFRC_Event', 'present')

  END SUBROUTINE clamsrdfrc_initialize
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsrdfrc_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: new_channel, new_attribute,&
                                           new_channel_object

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsrdfrc_init_memory'
    INTEGER       :: status
 
    ! Define channel RDFRC
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt(substr, status)

  END SUBROUTINE clamsrdfrc_init_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsrdfrc_init_coupling

    ! BMIL
    USE messy_main_channel_error_bi,    ONLY: channel_halt
    ! SMCL
    USE messy_main_channel,             ONLY: get_channel_object

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsrdfrc_init_coupling'
    INTEGER :: status

    ! CLaMS positions
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=LEV_D)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'E5_ZETA', p3=E5_ZETA_D)
    CALL channel_halt(substr, status)

    IF (use_o3forc) THEN
       CALL get_channel_object(status, 'clams', 'O3', p1=O3_CLAMS_D)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'rad4all', 'O3_pre', p3=O3_ECHAM_D)
       CALL channel_halt(substr, status)
    END IF
    IF (use_h2oforc) THEN
       CALL get_channel_object(status, 'clams', 'H2O', p1=H2O_CLAMS_D)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'rad4all', 'H2O_pre', p3=H2O_ECHAM_D)
       CALL channel_halt(substr, status)
    END IF

  END SUBROUTINE clamsrdfrc_init_coupling
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamsrdfrc_global_end

    USE messy_main_timer_bi,         ONLY: event_state
    USE messy_main_timer,            ONLY: current_date
    USE messy_main_transform_bi,     ONLY: trp_gpdc_gpgl, gather_glix
    USE messy_clamsbmix_global,      ONLY: lev_in_down
    USE messy_clams_global,          ONLY: nparts_max, dnparts, dnparts_max &
                                         , nz, ny, nx, mdi, eps
    USE messy_main_mpi_bi,           ONLY: p_pe, p_barrier, p_bcast
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: set_channel_output

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamsrdfrc_global_end'

    LOGICAL :: lrdfrcevent
    INTEGER :: ilev, ilat, ilon
    INTEGER :: bound_up, bound_down 
    INTEGER :: status,i
    LOGICAL :: check_mdi

    lrdfrcevent = event_state(rdfrcevent, current_date)

    IF (lrdfrcevent) THEN

       IF (.NOT.use_o3forc .AND. .NOT.use_h2oforc ) THEN
          WRITE(*,*) 'RDFRC: No forcing selected.'
          RETURN
       END IF

       WRITE(*,*) 'Active RDFRC Event'

       ! Collect global ECHAM field:
       IF (use_o3forc)  CALL trp_gpdc_gpgl (1, O3_ECHAM_D, O3_ECHAM_G)
       IF (use_h2oforc) CALL trp_gpdc_gpgl (1,H2O_ECHAM_D,H2O_ECHAM_G)
       CALL trp_gpdc_gpgl (1,  E5_ZETA_D,  E5_ZETA_G)

       ! Collect global CLaMS airparcel positions:
       ALLOCATE (LAT_G(nparts_max))
       ALLOCATE (LON_G(nparts_max))
       ALLOCATE (LEV_G(nparts_max))
       CALL gather_glix(LAT_G,LAT_D)
       CALL gather_glix(LON_G,LON_D)
       CALL gather_glix(LEV_G,LEV_D)

       ! Find ECHAM levels for interpolation: (E5_ZETA(nx, nlev, nlat))
       ! Upper boundary:
       bound_up = 1
       loop_up: DO ilev = 1, nz
          DO ilon = 1, nx
             DO ilat = 1, ny
                IF (E5_ZETA_G(ilon,ilev,ilat).LT.bound_up_zeta .AND. &
                     .NOT. ABS((E5_ZETA_G(ilon,ilev,ilat)-mdi)/mdi)<=eps) THEN
                   bound_up = ilev
                   EXIT loop_up
                END IF
             END DO
          END DO
          bound_up = ilev
       END DO loop_up
       ! Lower boundary:
       bound_down = nz
       loop_down: DO ilev = nz, 1, -1
          DO ilon = 1, nx
             DO ilat = 1, ny
                IF (E5_ZETA_G(ilon,ilev,ilat).GT.lev_in_down .AND. &
                     .NOT. ABS((E5_ZETA_G(ilon,ilev,ilat)-mdi)/mdi)<=eps) THEN
                   bound_down = ilev
                   EXIT loop_down
                END IF
             END DO
          END DO
       END DO loop_down

       ALLOCATE(ZETA_ECHAM_P(nx, bound_down-bound_up+1, ny))
       ZETA_ECHAM_P => E5_ZETA_G(:,bound_up:bound_down,:)
       ! O3 radiative forcing:
       IF (use_o3forc) THEN
          ALLOCATE (O3_CLAMS_G(nparts_max))
          CALL gather_glix(O3_CLAMS_G, O3_CLAMS_D)
          IF (p_pe==0) THEN
             ALLOCATE(O3_ECHAM_P(nx, bound_down-bound_up+1, ny))
             O3_ECHAM_P = O3_ECHAM_G(:,bound_up:bound_down,:)
             CALL rdfrc(LAT_G, LON_G, LEV_G, O3_CLAMS_G, O3_ECHAM_P,&
                  & ZETA_ECHAM_P,'O3')
             DEALLOCATE (O3_CLAMS_G)
          END IF
          CALL p_barrier()
       END IF


       ! Water vapor radiative forcing:
       IF (use_h2oforc) THEN

          check_mdi = .FALSE.

          ALLOCATE (H2O_CLAMS_G(nparts_max))
          CALL gather_glix(H2O_CLAMS_G,H2O_CLAMS_D)
          IF (p_pe==0) THEN
             ALLOCATE(H2O_ECHAM_P(nx, bound_down-bound_up+1, ny))
             H2O_ECHAM_P = H2O_ECHAM_G(:,bound_up:bound_down,:)
             CALL rdfrc(LAT_G, LON_G, LEV_G, H2O_CLAMS_G, H2O_ECHAM_P, ZETA_ECHAM_P,'H2O')
             DEALLOCATE (H2O_CLAMS_G)

             ! CLaMS Stratospheric water vapor to EMAC in specified regions:
             ! Above tropical level:
             DO ilev = bound_up, lev_bound_trop
                DO ilon = 1, nx
                   DO ilat = 1, ny
                      IF (ABS(H2O_ECHAM_P(ilon,ilev-bound_up+1,ilat)-mdi)/mdi<=eps) THEN
                         H2O_ECHAM_G(ilon,ilev,ilat) = H2O_ECHAM_P(ilon,ilev-bound_up+1,ilat)
                      ELSE
                         check_mdi = .TRUE.
                      END IF
                   END DO
                END DO
             END DO
             
             ! Add extratropical lower startosphere (works for latitudes in
             ! descending order):
             DO ilev = lev_bound_trop+1, lev_bound_extr 
                DO ilon = 1, nx
                   DO ilat = 1,bound_trop_n-1 ! Northern hemisphere
                      IF (ABS(H2O_ECHAM_P(ilon,ilev-bound_up+1,ilat)-mdi)/mdi&
                           &<=eps) THEN
                         H2O_ECHAM_G(ilon,ilev,ilat) = H2O_ECHAM_P(ilon,ilev&
                              &-bound_up+1,ilat)
                         ELSE
                            check_mdi = .TRUE.
                      END IF
                   END DO
                   DO ilat = bound_trop_s+1,ny ! Southern hemisphere
                      IF (ABS(H2O_ECHAM_P(ilon,ilev-bound_up+1,ilat)-mdi)/mdi&
                           &<=eps) THEN  
                         H2O_ECHAM_G(ilon,ilev,ilat) = H2O_ECHAM_P(ilon,ilev&
                              &-bound_up+1,ilat)  
                         ELSE
                            check_mdi = .TRUE.
                      END IF
                   END DO
                END DO
                
             END DO
             
             ! Check for missing vals in H2O_ECHAM_G:
             IF (check_mdi) THEN
                WRITE(*,*) 'Warning in RDFRC: Missing values occured in H2O&
                     & field and are not replaced in EMAC water vapor!'
             END IF

          END IF
          CALL p_barrier()
       END IF
       
       IF (use_o3forc)  CALL p_bcast( O3_ECHAM_G,0)
       IF (use_h2oforc) CALL p_bcast(H2O_ECHAM_G,0)

       ! von teilfeld zurueck auf ganzes ECHAM-feld(wie in tracer_init)
       IF (use_o3forc)  THEN
          CALL trp_gpdc_gpgl (-1, O3_ECHAM_D, O3_ECHAM_G)
       END IF
       IF (use_h2oforc) THEN
          CALL trp_gpdc_gpgl (-1,H2O_ECHAM_D,H2O_ECHAM_G)
       END IF
       
       CALL set_channel_output(status, 'clamsrdfrc', .TRUE.)
       CALL channel_halt(substr, status)
       CALL messy_write_output

       ! Deallocate arrays:
       NULLIFY(ZETA_ECHAM_P)
       IF (p_pe==0) THEN
          IF (use_o3forc)  DEALLOCATE(O3_ECHAM_P)
          IF (use_h2oforc) DEALLOCATE(H2O_ECHAM_P)
          DEALLOCATE (LAT_G)
          DEALLOCATE (LON_G)
          DEALLOCATE (LEV_G)
       END IF

    END IF

  END SUBROUTINE clamsrdfrc_global_end
!--------------------------------------------------------------------

!--------------------------------------------------------------------
#endif
END MODULE messy_clamsrdfrc_si
