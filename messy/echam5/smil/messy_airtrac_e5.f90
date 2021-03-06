! **********************************************************************
! INTERFACE TO ECHAM5 FOR
! AIRTRAC module to track air traffic emissions and calculate
! simple chemistry on trajectories
! Author : Christine Froemming, DLR, 2010/2011
!          Sabine Brinkop & Patrick Joeckel, DLR, 2014
!
! **********************************************************************
!
! Reference:
!
! * Frömming, C., Grewe, V., Brinkop, S., and Jöckel, P.: Documentation of the
!   EMAC submodels AIRTRAC 1.0 and CONTRAIL 1.0, supplementary material of
!   Grewe et al., 2014a, Geoscientific Model Development, 7, 175–201,
!   https://doi.org/10.5194/gmd-7-175-2014, 2014.
! **********************************************************************

! **********************************************************************
MODULE messy_airtrac_e5
! **********************************************************************

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi,error_bi
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_3D_ARRAY &
                                    , PTR_0D_ARRAY
  USE messy_main_constants_mem, ONLY: MN, MO, MH, MC, M_H2O
  USE messy_airtrac

  IMPLICIT NONE
  PRIVATE
  INTRINSIC :: NULL, REAL

  ! BACKGROUND AND DIAGNOSTIC (PROD/LOSS) TRACER IN GRIDPOINT SPACE
  ! REQUIRED FOR PERTURBATION CHEMISTRY
  ! (YOU MIGHT WANT TO ADD MORE TRACERS, PROD-/LOSSRATES)
  ! NUMBER
  INTEGER, PARAMETER:: ngptrac = 22
  ! NAMES
  CHARACTER (LEN=*), DIMENSION(ngptrac), PARAMETER :: gpname=(/ &
       'O3      ', 'NO      ', 'NO2     ', 'ProdO3N ', 'LossO3N ', 'LossO3Y ',&
       'LossNOx ', 'LossHNO3', 'HNO3    ', 'OH      ', 'HO2     ', 'ProdOH1 ',&
       'ProdOH2 ', 'ProdOH3 ', 'LossOH1 ', 'LossOH2 ', 'LossOH3 ', 'LossOH4 ',&
       'LossOH5 ', 'ProdHO21', 'LossHO21', 'LossHO22'/)
  ! INDICES
  INTEGER, DIMENSION(ngptrac):: gpidx = 0
  ! FOR CONVERSION FROM GP TO LG
  TYPE(PTR_1D_ARRAY), DIMENSION(ngptrac), SAVE :: gp2lgtptr
  ! (incl. tendencies)
  TYPE(PTR_1D_ARRAY), DIMENSION(ngptrac), SAVE :: gp2lgtptr_te

  ! channel object for HNO3 tendency from SCAV ...
!qqq implemenmt alternative WITH TENDENCY
  REAL(dp), DIMENSION(:,:,:), POINTER :: hno3_tte_scav    => NULL()
  ! .. converted from GP to LG
  REAL(dp), DIMENSION(:),     POINTER :: hno3_tte_scav_lg => NULL()

  ! LAGRANGIAN PROGNOSTIC TRACERS
  ! NUMBER
  INTEGER, PARAMETER:: nlgtrac=13
  ! NAMES
  CHARACTER (LEN=*), DIMENSION(nlgtrac), PARAMETER :: lgname=(/ &
        'airNOx    ', 'airO3     ', 'airProdO3N', 'airLossO3N', 'airLossO3Y' &
       ,'airLossNOx', 'airHNO3   ', 'airOH     ', 'airHO2    ', 'airCH4    ' &
       ,'airH2O    ', 'airScvHNO3', 'airLosHNO3' &
       /)
  ! MOLAR MASSES
  REAL(DP), DIMENSION(nlgtrac):: lg_molarmass=(/ &
       MN+MO,       3.0_dp*MO, 3.0_dp*MO,    3.0_dp*MO,    3.0_dp*MO,   &
       MN+MO, MH+MN+3.0_dp*MO,     MO+MH, MH+2.0_dp*MO, MC+4.0_dp*MH,   &
       M_H2O, MH+MN+3.0_dp*MO, MH+MN+3.0_dp*MO &
       /)
  ! INDICES
  INTEGER, DIMENSION(:,:),POINTER:: lgidx=> NULL()
  ! POINTER TO TENDENCIES
  TYPE(PTR_1D_ARRAY), DIMENSION(:,:), POINTER, SAVE :: lgtptr_te => NULL()

  ! channel objects for coupling to cloud and convect
  REAL(dp), DIMENSION(:,:,:), POINTER :: cv_rform    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: cv_sform    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: cv_cover    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: cl_rform    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: cl_sform    => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: cl_aclc     => NULL()
  ! LG
  REAL(dp), DIMENSION(:),     POINTER :: rain_sum_lg    => NULL()

  ! POINTER FOR CHANNEL OBJECT WITH NUMBER OF EMISSION POINTS (AS REAL)
  REAL(dp), POINTER :: r_emis_points => NULL()

  ! -------------------------------------------------------------

  ! PUBLIC SMIL ROUTINES
  PUBLIC :: airtrac_initialize       ! global initialisation of module
  PUBLIC :: airtrac_new_tracer       ! define new tracers
  PUBLIC :: airtrac_init_memory      ! allocate memory and define stream
  PUBLIC :: airtrac_init_coupling    ! Initialize
  PUBLIC :: airtrac_global_end
  PUBLIC :: airtrac_free_memory

  ! -------------------------------------------------------------

CONTAINS

! ************************************************************************
! PUBLIC ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ------------------------------------------------------------------------
  SUBROUTINE  airtrac_initialize

    ! AIRTRAC MODULE ROUTINE (ECHAM-5 INTERFACE)
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! Author: Christine Froemming, DLR, August 2010

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,        ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtrac_initialize'
    INTEGER         :: iou    ! I/O unit
    INTEGER         :: status ! error status

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL airtrac_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi('', substr)
    END IF
    !
    CALL p_bcast(n_emis_points, p_io)

  END SUBROUTINE airtrac_initialize
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE airtrac_new_tracer

    ! AIRTRAC MODULE ROUTINE (SMIL)
    ! defines AIRTRAC specific Lagrangian tracers
    ! Author: Christine Froemming, DLR, 2010/2011

    USE messy_main_mpi_bi,          ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,   ONLY: LGTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tools,           ONLY: int2str
    USE messy_main_tracer,          ONLY: new_tracer, set_tracer       &
                                        , ON, OFF, R_molarmass, I_scav &
                                        , I_drydep, I_sedi             &
                                        , I_mix ,I_advect              &
                                        , I_convect, I_vdiff

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! LOCAL
    INTEGER :: status, jt,je
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtrac_new_tracer'
    CHARACTER(LEN=3):: istr=''

    CALL start_message_bi(modstr, 'DEFINE LG TRACERS', substr)

    ALLOCATE (lgidx(n_emis_points,nlgtrac))
    lgidx(:,:)=0

    DO je=1, n_emis_points
       CALL int2str(istr,je)

       DO jt=1, nlgtrac
          IF (p_parallel_io) &
               WRITE(*,*) 'creating lg tracer ', TRIM(lgname(jt))//'_'//istr

          CALL new_tracer(status, LGTRSTR, TRIM(lgname(jt)),modstr      &
               , subname=istr, idx=lgidx(je,jt) , unit='mol/mol')
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, lgidx(je,jt) &
               , R_molarmass, r=lg_molarmass(jt))
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, lgidx(je,jt), I_advect, i=OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, lgidx(je,jt), I_vdiff, i=OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, lgidx(je,jt), I_convect, i=OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, lgidx(je,jt), I_scav, i=OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, lgidx(je,jt), I_drydep, i=OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, lgidx(je,jt), I_sedi, i=OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, lgidx(je,jt), I_mix, i=ON)
          CALL tracer_halt(substr, status)
       END DO
    END DO

    CALL end_message_bi(modstr, 'DEFINE LG TRACERS', substr)

  END SUBROUTINE airtrac_new_tracer
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE airtrac_init_memory

    ! AIRTRAC MODULE ROUTINE (SMIL)
    ! defines AIRTRAC specific channel objects and allocates memory for
    ! global fields
    ! Author: Christine Froemming, DLR, 2010/2011

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_tools,            ONLY: int2str
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: LG_ATTILA, SCALAR
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtrac_init_memory'
    CHARACTER(LEN=3)            :: istr=''
    INTEGER                     :: jt,je
    INTEGER                     :: status

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ALLOCATE (lgtptr_te(n_emis_points, nlgtrac))
    DO je=1, n_emis_points
       DO jt=1, nlgtrac
          lgtptr_te(je,jt)%ptr=>NULL()
       END DO
    END DO

    ! DEFINE NEW CHANNEL
    CALL new_channel(status, modstr,lrestreq = .TRUE.)
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) 'creating channel object r_emis_points'
    CALL new_channel_object(status, modstr, 'r_emis_points' &
         , p0=r_emis_points, reprid=SCALAR)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'r_emis_points'  &
         , 'units', c='' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'r_emis_points' &
         , 'long_name', c='number of emission points ')
    CALL channel_halt(substr, status)
    r_emis_points=REAL(n_emis_points)
    IF (p_parallel_io) &
         WRITE(*,*) 'r_emis_points = ',r_emis_points &
         , 'n_emis_points = ',n_emis_points

    DO jt=1,ngptrac
       IF (p_parallel_io) WRITE(*,*) 'creating channel object ' &
            , TRIM(gpname(jt))
       CALL new_channel_object(status, modstr, TRIM(gpname(jt)) &
            , p1 = gp2lgtptr(jt)%ptr, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(gpname(jt))  &
            , 'units', c='mol/mol' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(gpname(jt)) &
            , 'long_name', c='transformed gp tracer '//TRIM(gpname(jt)))
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) 'creating channel object ' &
            , TRIM(gpname(jt))//'_te'
       CALL new_channel_object(status, modstr, TRIM(gpname(jt))//'_te' &
            , p1 = gp2lgtptr_te(jt)%ptr, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(gpname(jt))//'_te' &
            , 'units', c='mol/mol/s' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, TRIM(gpname(jt))//'_te' &
            , 'long_name' &
            , c='transformed gp tracer tendency'//TRIM(gpname(jt)))
       CALL channel_halt(substr, status)
    END DO

    ! new channel object for transformed HNO3 tendency from SCAV
    IF (p_parallel_io) WRITE(*,*) 'creating channel object hno3_tte_scav_lg'
    CALL new_channel_object(status, modstr, 'hno3_tte_scav_lg' &
         , p1 = hno3_tte_scav_lg, reprid=LG_ATTILA)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hno3_tte_scav_lg' &
         , 'units', c='mol/mol/s' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hno3_tte_scav_lg' &
         , 'long_name', c='transformed gp hno3 scavenging tendency')
    CALL channel_halt(substr, status)


    IF (p_parallel_io) WRITE(*,*) 'creating channel object rain_sum_lg'
    CALL new_channel_object(status, modstr, 'rain_sum_lg' &
         , p1 = rain_sum_lg, reprid=LG_ATTILA)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rain_sum_lg' &
         , 'units', c='mol/mol' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rain_sum_lg' &
         , 'long_name', c='transformed gp rain sum from convect and cloud')
    CALL channel_halt(substr, status)

    DO je=1,n_emis_points
       CALL int2str(istr,je)
       DO jt=1,nlgtrac
          IF (p_parallel_io) WRITE(*,*) 'creating lg channel object ' &
               , TRIM(lgname(jt))//'_'//istr//'_te'
          CALL new_channel_object(status, modstr &
               , TRIM(lgname(jt))//'_'//istr//'_te' &
               , p1 = lgtptr_te(je,jt)%ptr, reprid=LG_ATTILA)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr &
               , TRIM(lgname(jt))//'_'//istr//'_te' &
               , 'units', c='mol/mol/s')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr &
               , TRIM(lgname(jt))//'_'//istr//'_te' &
               , 'long_name' &
               , c='lg tracer tendency '//TRIM(lgname(jt))//'_'//istr)
          CALL channel_halt(substr, status)
       END DO
    END DO

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE airtrac_init_memory
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE airtrac_init_coupling

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer,           ONLY: get_tracer
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_blather_bi,       ONLY: info_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'airtrac_init_coupling'
    INTEGER                     :: status
    INTEGER :: jt
    CHARACTER(LEN=3) :: istr = ''  ! len=3: 001-999

    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output

    ! GET INDICES OF GP TRACERS
    DO jt=1, ngptrac
       IF (p_parallel_io) WRITE(*,*) 'looking for gp tracer ',TRIM(gpname(jt))
       CALL get_tracer(status, GPTRSTR, TRIM(gpname(jt)), idx=gpidx(jt))
       CALL TRACER_HALT(substr,status)
    END DO

    ! SET POINTER TO object HNO3_scte in channel scav_gp
    CALL info_bi('looking for channel scav_gp / object HNO3_scte',' ')
    CALL get_channel_object(status, 'scav_gp', 'HNO3_scte' &
        , p3=hno3_tte_scav)
    CALL channel_halt(substr, status)


    IF (p_parallel_io) &
         WRITE(*,*) 'Checking for parameters from convect channel ...'
    !
    CALL get_channel_object(status, 'convect', 'cv_rform', p3=cv_rform)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'convect', 'cv_sform', p3=cv_sform)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'convect', 'cv_cover', p3=cv_cover)
    CALL channel_halt(substr, status)

    IF (p_parallel_io) &
         WRITE(*,*) 'Checking for parameters from cloud channel ...'
    !
    CALL get_channel_object(status, 'cloud', 'rain_form', p3=cl_rform)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'cloud', 'snow_form', p3=cl_sform)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'cloud', 'aclc'     , p3=cl_aclc)
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr,'COUPLING',substr)

  END SUBROUTINE airtrac_init_coupling
! ------------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE airtrac_global_end

      ! ECHAM5/MESSy
      USE messy_main_tracer_mem_bi,   ONLY: qxtte_a &
                                          , qxtm1_a &
                                          , xtte    &
                                          , xtm1    &
                                          , NCELL
      USE messy_attila_tools_e5,      ONLY: gp2lg_e5
      USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks, nlev
      USE messy_main_timer,           ONLY: ztmst=>time_step_len
      USE messy_main_data_bi,         ONLY: qm1_3d, qte_3d
      USE messy_main_constants_mem,   ONLY: M_H2O, M_AIR

      IMPLICIT NONE
      INTRINSIC :: MAX

      ! LOCAL
      REAL(DP), PARAMETER :: scvmr = M_air/M_H2O
      REAL(DP), DIMENSION(:,:,:), POINTER     :: zgp => NULL ()
      REAL(DP), DIMENSION(:,:,:), POINTER     :: zgp_te => NULL ()
      REAL(DP), DIMENSION(:,:,:), POINTER     :: rain_sum_gp => NULL ()
      REAL(DP), DIMENSION(:,:,:), POINTER     :: zh2o_gp => NULL ()
      REAL(DP), DIMENSION(:), POINTER         :: zh2o_lg => NULL ()
      TYPE(PTR_1D_ARRAY), DIMENSION(nlgtrac)  :: zlg
      INTEGER                                 :: jt, je

      ! CALCULATE START VALUES FOR GP BACKGROUND TRACERS AND TRANSFORM TO LG
      ALLOCATE(zgp(   nproma, nlev, ngpblks))
      ALLOCATE(zgp_te(nproma, nlev, ngpblks))
      zgp(:,:,:)    = 0.0
      zgp_te(:,:,:) = 0.0
      DO jt=1, ngptrac
         zgp(:,:,:)    = xtm1(:,:,gpidx(jt),:) + xtte(:,:,gpidx(jt),:)*ztmst
         zgp_te(:,:,:) = xtte(:,:,gpidx(jt),:)
         CALL gp2lg_e5(zgp,    gp2lgtptr(jt)%ptr)
         CALL gp2lg_e5(zgp_te, gp2lgtptr_te(jt)%ptr)
      END DO
      DEALLOCATE(zgp);    NULLIFY(zgp)
      DEALLOCATE(zgp_te); NULLIFY(zgp_te)

      ALLOCATE(hno3_tte_scav_lg(NCELL))
      ALLOCATE(zh2o_lg(NCELL))
      ALLOCATE(rain_sum_lg(NCELL))

      ! TRANSFORM HNO3 tendency from SCAV to LG
      CALL gp2lg_e5(hno3_tte_scav, hno3_tte_scav_lg)

      ! Transform background water vapour
      ALLOCATE(zh2o_gp(nproma, nlev, ngpblks))
      zh2o_gp(:,:,:)    =  qm1_3d(:,:,:) + qte_3d(:,:,:)*ztmst
      zh2o_gp(:,:,:)    =  scvmr * (zh2o_gp(:,:,:)/(1.0_DP-zh2o_gp(:,:,:)))
      CALL gp2lg_e5(zh2o_gp,zh2o_lg)
      DEALLOCATE(zh2o_gp)

      ALLOCATE(rain_sum_gp(nproma, nlev, ngpblks))
      rain_sum_gp(:,:,:) = ((cv_rform(:,:,:) + cv_sform(:,:,:)) &
           * cv_cover(:,:,:)                        &
           + (cl_rform(:,:,:) + cl_sform(:,:,:))    &
           * cl_aclc(:,:,:)) * M_AIR / M_H2O
      rain_sum_gp(:,:,:) = MAX(0._dp, rain_sum_gp(:,:,:))
      !
      CALL gp2lg_e5(rain_sum_gp, rain_sum_lg)
      DEALLOCATE(rain_sum_gp)

      ! MEMORY FOR START VALUES OF LAGRANGIAN TRACERS
      DO jt=1, nlgtrac
         ALLOCATE(zlg(jt)%ptr(NCELL))
      END DO

      emis_pts: DO je=1, n_emis_points

         ! CALCULATE START VALUES OF LAGRANGIAN TRACERS
         DO jt=1, nlgtrac
            zlg(jt)%ptr(:)= qxtm1_a(:,lgidx(je,jt)) &
                 + qxtte_a(:,lgidx(je,jt))*ztmst
         END DO

         CALL airtrac_integrate( ztmst    &        ! time step length
              , gp2lgtptr(1)%ptr(:)       &        ! O3 background
              , gp2lgtptr(2)%ptr(:)       &        ! NO background
              , gp2lgtptr(3)%ptr(:)       &        ! NO2 background
              , gp2lgtptr_te(4)%ptr(:)    &        ! ProdO3N tendency
              , gp2lgtptr_te(5)%ptr(:)    &        ! LossO3N tendency
              , gp2lgtptr_te(6)%ptr(:)    &        ! LossO3Y tendency
              , gp2lgtptr_te(7)%ptr(:)    &        ! LossNOx tendency
              , gp2lgtptr_te(8)%ptr(:)    &        ! LossHNO3 tendency
              , gp2lgtptr(9)%ptr(:)       &        ! HNO3 background
              , hno3_tte_scav_lg(:)       &        ! HNO3 scavenging tendency
              , gp2lgtptr(10)%ptr(:)      &        ! OH
              , gp2lgtptr(11)%ptr(:)      &        ! HO2
              , gp2lgtptr_te(12)%ptr(:)   &        ! prodoh1
              , gp2lgtptr_te(13)%ptr(:)   &        ! prodoh2
              , gp2lgtptr_te(14)%ptr(:)   &        ! prodoh3
              , gp2lgtptr_te(15)%ptr(:)   &        ! lossoh1
              , gp2lgtptr_te(16)%ptr(:)   &        ! lossoh2
              , gp2lgtptr_te(17)%ptr(:)   &        ! lossoh3
              , gp2lgtptr_te(18)%ptr(:),gp2lgtptr_te(19)%ptr(:) & ! lossoh4,5
              , gp2lgtptr_te(20)%ptr(:),gp2lgtptr_te(21)%ptr(:) & ! prodho21,2
              , gp2lgtptr_te(22)%ptr(:)  &         ! lossho22
              , zlg(1)%ptr(:),zlg(2)%ptr(:) &    ! airNOx,airO3
              , zlg(3)%ptr(:),zlg(4)%ptr(:) &    ! airProdO3N,airLossO3N
              , zlg(5)%ptr(:),zlg(6)%ptr(:) &    ! airLossO3Y,airLossNOx
              , zlg(7)%ptr(:),zlg(8)%ptr(:) &    ! airHNO3,airOH
              , zlg(9)%ptr(:)               &    ! airHO2
              , lgtptr_te(je,1)%ptr(:),lgtptr_te(je,2)%ptr(:)   & ! airNOx_te,airO3_te
              , lgtptr_te(je,3)%ptr(:),lgtptr_te(je,4)%ptr(:)   & ! airProdO3N_te,airLossO3N_te
              , lgtptr_te(je,5)%ptr(:),lgtptr_te(je,6)%ptr(:)   & ! airLossO3Y_te,airLossNOx_te
              , lgtptr_te(je,7)%ptr(:),lgtptr_te(je,8)%ptr(:)   & ! airHNO3_te,airOH_te
              , lgtptr_te(je,9)%ptr(:),lgtptr_te(je,10)%ptr(:)  & ! airHO2_te,airCH4_te
              , lgtptr_te(je,12)%ptr(:),lgtptr_te(je,13)%ptr(:) & ! hno3scav_te, hno3loss_te
              , zh2o_lg(:),zlg(11)%ptr(:)                          &  ! H2O background,airH2O
              , rain_sum_lg(:),lgtptr_te(je,11)%ptr(:)          &  ! rain sum,H2O tendency from integrate
              )


         ! UPDATE LAGRANGIAN TRACER TENDENCIES
         DO jt=1,nlgtrac
            qxtte_a(:,lgidx(je,jt)) = qxtte_a(:,lgidx(je,jt)) &
                 + lgtptr_te(je,jt)%ptr(:)
         END DO

      END DO emis_pts

      ! FREE TEMPORARY MEMORY
      DO jt=1, nlgtrac
         DEALLOCATE(zlg(jt)%ptr); NULLIFY(zlg(jt)%ptr)
      END DO

      IF (ASSOCIATED (zh2o_lg)) then
         DEALLOCATE(zh2o_lg); NULLIFY(zh2o_lg)
      ENDIF

      IF (ASSOCIATED (hno3_tte_scav_lg)) then
         DEALLOCATE(hno3_tte_scav_lg); NULLIFY(hno3_tte_scav_lg)
      ENDIF

      IF (ASSOCIATED (rain_sum_lg)) then
         DEALLOCATE(rain_sum_lg); NULLIFY(rain_sum_lg)
      ENDIF

  END SUBROUTINE airtrac_global_end
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE airtrac_free_memory

    IMPLICIT NONE

    INTEGER :: je, jt

    DEALLOCATE (lgidx)

    DO je=1,n_emis_points
       DO jt=1,nlgtrac
          lgtptr_te(je,jt)%ptr=>NULL()
       END DO
    END DO
    DEALLOCATE (lgtptr_te)

  END SUBROUTINE airtrac_free_memory
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_airtrac_e5
! **********************************************************************
