#include "messy_main_ppd_bi.inc"

MODULE MESSY_CAT_E5

  ! This module should allow the calculation of CAT indices
  ! and associated turbulence
  ! a) diagnostically
  ! b) allow the formulation of a CAT driven turbulent mixing scheme
  ! Author: H. Tost
  ! JGU Mainz, Aug, 2014

  USE messy_cat 
  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi, warning_bi

  USE messy_main_constants_mem, ONLY: dp

  USE messy_main_data_bi,       ONLY: tpoteq => tpoteq2_3d
  
  IMPLICIT NONE

  REAL(dp), DIMENSION(:,:,:), POINTER :: dudx => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: dvdx => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: dudy => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: dvdy => NULL()
  ! REAL(dp), DIMENSION(:,:,:), POINTER :: zup  => NULL() !!!
  ! REAL(dp), DIMENSION(:,:,:), POINTER :: zvp  => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: zup_cos  => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: zvp_cos  => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: v_total  => NULL() !!!

!!! careful, this is the absolute of the horizontal wind for vertical shear
  REAL(dp), DIMENSION(:,:,:), POINTER :: dvdz_down   => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: dvdz_cen    => NULL() !!!

  REAL(dp), DIMENSION(:,:,:), POINTER :: DVSI_down   => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: DVSI_cen    => NULL() !!!

  TYPE div_term
     REAL(dp), DIMENSION(:,:,:), POINTER :: d_term   => NULL()
  END type div_term
  TYPE(div_term), DIMENSION(:), POINTER :: time_cnt  => NULL()
  INTEGER  :: nti
 
  REAL(dp), DIMENSION(:,:,:), POINTER :: DVT       => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: DDVSI     => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: DDVSI_MOD => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: D_DDVSI   => NULL() !!!


  REAL(dp), DIMENSION(:,:,:), POINTER :: dpotdz   => NULL() !!!
  REAL(dp), DIMENSION(:,:,:), POINTER :: N2       => NULL() !!!

  INTEGER :: up_idx, lo_idx

CONTAINS
!--------------------------------------------------------------------------------
  SUBROUTINE CAT_INITIALIZE

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast,  &
                                     finish, message 
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    INTEGER   :: iou      ! I/O unit
    INTEGER   :: status   ! status
    CHARACTER(LEN=*), PARAMETER::substr='cat_initialize'
    CHARACTER(LEN=4)  :: str_num

    ! Initialize main-ctrl
    status = 1
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CAT CORE ROUTINE:
       CALL cat_read_nml_ctrl(iou, status)
       if (status/=0) CALL FINISH('CAT INIT')
    END IF
   
    CALL MESSAGE (substr, 'CAT active')
    
    CALL P_BCAST(cat_param, p_io)
    IF (cat_param == 1)   &
      CALL MESSAGE (substr, 'Ellrod & Knox TI calculation')
    str_num=""
    CALL P_BCAST(nhours, p_io)
    write (str_num,'(F4.1)') nhours
    CALL MESSAGE (substr, str_num//' hr forecast time selected')
    str_num=""
    CALL P_BCAST(DIV_C, p_io)
    write (str_num,'(F4.1)') DIV_C
    CALL MESSAGE (substr, 'Constant of DVT scaling = '//str_num)
    CALL P_BCAST(l_tracmix, p_io)
    IF (L_tracmix) CALL MESSAGE (substr, 'Tracer mixing by CAT enabled')
    CALL P_BCAST(USE_DDVSI, p_io)
    IF (USE_DDVSI) CALL MESSAGE (substr, &
         'Mixing coefficient determined by inmodified turbulence index')

    ! INITIALIZE COUPLING-CONTROL
!!$    IF (p_parallel_io) THEN
!!$       iou = find_next_free_unit(100,200)
!!$       CALL cat_read_nml_cpl(status, iou)
!!$       IF (status /= 0) CALL finish('cat COUPLING INIT')
!!$    END IF

  END SUBROUTINE CAT_INITIALIZE
!--------------------------------------------------------------------------------
  SUBROUTINE cat_read_nml_cpl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE
!    NAMELIST /CPL/

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='cat_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
 
!!$    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
!!$    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
!!$
!!$    READ(iou, NML=CPL, IOSTAT=fstat)
!!$    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
!!$    IF (fstat /= 0) RETURN  ! error while reading namelist
!!$  
!!$    CALL read_nml_close(substr, iou, modstr)
!!$    status = 0 ! NO ERROR

  END SUBROUTINE cat_read_nml_cpl
!--------------------------------------------------------------------------------

  SUBROUTINE CAT_INIT_MEMORY

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_timer,            ONLY: delta_time 


    CHARACTER(LEN=*), PARAMETER::substr='cat_init_memory'
    INTEGER :: status, i
    CHARACTER(LEN=2) :: str_num

    CALL start_message_bi(modstr,'MEMORY INITIALIZATION', substr)

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dudx', p3=dudx)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dudx' &
         , 'long_name', c='zonal gradient of u wind')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dudx', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dudy', p3=dudy)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dudy' &
         , 'long_name', c='meridional gradient of u wind')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dudy', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dvdx', p3=dvdx)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dvdx' &
         , 'long_name', c='zonal gradient of v wind')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dvdx', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dvdy', p3=dvdy)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dvdy' &
         , 'long_name', c='meridional gradient of v wind')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dvdy', 'units', c='1/s')
    CALL channel_halt(substr, status)

    !CALL new_channel_object(status, modstr, 'zup', p3=zup) !!!
    !CALL channel_halt(substr, status)
    !CALL new_attribute(status, modstr, 'zup' &
    !     , 'long_name', c='u wind')
    !CALL channel_halt(substr, status)
    !CALL new_attribute(status, modstr, 'zup', 'units', c='m/s')
    !CALL channel_halt(substr, status)

    !CALL new_channel_object(status, modstr, 'zvp', p3=zvp) !!!
    !CALL channel_halt(substr, status)
    !CALL new_attribute(status, modstr, 'zvp' &
    !     , 'long_name', c='v wind')
    !CALL channel_halt(substr, status)
    !CALL new_attribute(status, modstr, 'zvp', 'units', c='m/s')
    !CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'zup_cos', p3=zup_cos) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zup_cos' &
         , 'long_name', c='u wind * cos(lat)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zup_cos', 'units', c='m/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'zvp_cos', p3=zvp_cos) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zvp_cos' &
         , 'long_name', c='v wind * cos(lat)')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'zvp_cos', 'units', c='m/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'v_total', p3=v_total) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v_total' &
         , 'long_name', c='total wind vector')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'v_total', 'units', c='m/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vertical_shear_down', p3=dvdz_down) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vertical_shear_down' &
         , 'long_name', c='vertical shear of the 2D horizontal wind')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vertical_shear_down', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'vertical_shear_cen', p3=dvdz_cen) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vertical_shear_cen' &
         , 'long_name', c='vertical shear of the 2D horizontal wind')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'vertical_shear_cen', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DVSI_down', p3=DVSI_down) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DVSI_down' &
         , 'long_name', c='Ellrod turbulence index')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DVSI_down', 'units', c='1/s^2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DVSI_cen', p3=DVSI_cen) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DVSI_cen' &
         , 'long_name', c='Ellrod turbulence index')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DVSI_cen', 'units', c='1/s^2')
    CALL channel_halt(substr, status)

 !   nhours = 6._dp  ! be careful with choice
    nti = 3600._dp * nhours / delta_time

    ALLOCATE(time_cnt(nti))
    DO i=1,nti
        IF (i < 10) then
           write (str_num,'(A,I1)') "0",i
        ELSE
           write (str_num,'(I2)') i
        END IF
       CALL new_channel_object(status, modstr, 'div_t'//str_num, p3=time_cnt(i)%d_term)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'div_t'//str_num &
            , 'long_name', c='part of divergence trend term')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'div_t'//str_num, 'units', c='1/s')
       CALL channel_halt(substr, status)
    END DO
 
    CALL new_channel_object(status, modstr, 'DVT', p3=DVT) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DVT' &
         , 'long_name', c='divergence trend term')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DVT', 'units', c='1/s')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DDVSI', p3=DDVSI) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DDVSI' &
         , 'long_name', c='full Ellrod turbulence index')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DDVSI', 'units', c='1/s^2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'DDVSI_MOD', p3=DDVSI_MOD) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DDVSI_MOD' &
         , 'long_name', c='modified Ellrod turbulence index')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'DDVSI_MOD', 'units', c='1/s^2')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'D_DDVSI', p3=D_DDVSI) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'D_DDVSI' &
         , 'long_name', c='mixing coefficient')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'D_DDVSI', 'units', c='-')
    CALL channel_halt(substr, status)
    D_DDVSI = 0._dp


!!$    CALL new_channel_object(status, modstr, 'tpoteq', p3=tpoteq) !!!
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'tpoteq' &
!!$         , 'long_name', c='eq. pot. temperature')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr, 'tpoteq', 'units', c='K')
!!$    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'dpotdz', p3=dpotdz) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dpotdz' &
         , 'long_name', c='deriv. of tpoteq')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dpotdz', 'units', c='K/m')
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'N2', p3=N2) !!!
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'N2' &
         , 'long_name', c='Brunt Vaisala Frequency')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'N2', 'units', c='1/s^2')
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr,'MEMORY INITIALIZATION', substr)

  END SUBROUTINE CAT_INIT_MEMORY
!--------------------------------------------------------------------------------
  SUBROUTINE CAT_INIT_COUPLING

    USE messy_main_grid_def_mem_bi, ONLY: nlev, vct, nvclev

    USE messy_main_mpi_bi,     ONLY: message 
    USE messy_main_tools,      ONLY: int2str
    REAL(dp) :: press_help(nlev)
    CHARACTER(LEN=*), PARAMETER::substr='cat_init_coupling'
    CHARACTER(LEN=3)  :: str_num1, str_num2
    INTEGER           :: jk

    IF (.NOT. l_tracmix) RETURN

    DO jk=1,nlev
       press_help(jk) = vct(jk) + vct(jk+nvclev) * 101325._dp
    ENDDO
    up_idx = 1
    DO jk=1,nlev
       IF (press_help(jk) > 7000._dp)  EXIT
       up_idx = jk
    END DO
    lo_idx = nlev
    DO jk=nlev,up_idx,-1
       IF (press_help(jk) < 50000._dp)  EXIT
       lo_idx = jk
    END DO
    CALL int2str(str_num1,up_idx)
    CALL int2str(str_num2,lo_idx)
    CALL MESSAGE(substr, &
         'Tracer mixing only active between levels '//&
         str_num1//' and '//str_num2) 

  END SUBROUTINE CAT_INIT_COUPLING
!--------------------------------------------------------------------------------
  SUBROUTINE CAT_GLOBAL_START

    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks
    USE messy_main_grid_def_bi,     ONLY: coslat_2d
    USE messy_main_data_bi,         ONLY: u => u_scb,             &
                                          v => v_scb,             &
                                          vom => vom_scb,         &
                                          vol => vol_scb,         &
                                          geopot => geopot_3d
    USE messy_main_constants_mem, ONLY: radius_earth, pi, g
    USE messy_main_timer,         ONLY: ztmst => time_step_len

    USE messy_main_mpi_bi,        ONLY: p_pe, dcl, dcg
   
    USE mo_transpose,             ONLY: tr_gp_ffsl


    REAL(dp) :: zup1(nproma,nlev,ngpblks)
    REAL(dp) :: zvp1(nproma,nlev,ngpblks)

    REAL(dp), DIMENSION(:,:,:), POINTER :: ub => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: vb => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: dudxb => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: dvdxb => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: dudyb => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: dvdyb => NULL()

    REAL(dp), DIMENSION(:,:,:), POINTER :: u2 => NULL()
    REAL(dp), DIMENSION(:,:,:), POINTER :: v2 => NULL()

    INTEGER :: jfirst, jlast, plon, plat, klev, plev, kstart, kstep
    INTEGER :: kproma, zjrow, nlat
    INTEGER :: j, jk, jl, jm, i, jt
    INTEGER, PARAMETER :: ng = 1
    INTEGER, PARAMETER :: nq = 1

    REAL(dp) :: cor_lat(nproma)
    REAL(dp) :: dx, dy
    REAL(dp) :: pasl, kappa

    REAL(dp), DIMENSION(nproma,nlev) :: dz !!!

    REAL(dp), PARAMETER :: N2_lim = 3.e-4_dp !!!
    REAL(dp), PARAMETER :: t_norm = (86400._dp/4.e6_dp) !!!
    REAL(dp) :: D_norm
    REAL(dp) :: N2_lt(nproma,nlev) !!!
    

    jfirst = dcl%ffsl%lats
    jlast  = dcl%ffsl%latn

    plon  = dcl%nlon
    plat  = dcl%nglat
    plev  = dcl%nlev
    jm    = dcl%nlat
    kstart = dcl%set_b
    kstep  = dcl%nprocb
    klev   = (plev-kstart) / kstep + 1

    ALLOCATE(ub(plon, jfirst:jlast, klev))
    ALLOCATE(vb(plon, jfirst:jlast, klev))
    ALLOCATE(dudxb(plon, jfirst:jlast, klev))
    ALLOCATE(dvdxb(plon, jfirst:jlast, klev))
    ALLOCATE(dudyb(plon, jfirst:jlast, klev))
    ALLOCATE(dvdyb(plon, jfirst:jlast, klev))

    ALLOCATE(u2(plon, jfirst-ng:jlast+ng, klev))
    ALLOCATE(v2(plon, jfirst-ng:jlast+ng, klev))

    dudx = 0._dp
    dudy = 0._dp
    dvdx = 0._dp
    dvdy = 0._dp
    u2   = 0._dp
    v2   = 0._dp
    ! zup  = 0._dp !!!
    ! zvp  = 0._dp !!!
    zup_cos  = 0._dp !!!
    zvp_cos  = 0._dp !!!
    v_total  = 0._dp !!!
    dz = 0._dp !!!
    D_norm = ztmst / t_norm !!!

    ! collect horizontal winds in local decomposition
    DO zjrow=1, dcl%ngpblks
       IF ( zjrow == dcl%ngpblks ) THEN
          kproma = dcl%npromz
       ELSE
          kproma = dcl%nproma
       END IF
       DO jl=1,kproma
          cor_lat(jl) = 1._dp / coslat_2d(jl,zjrow)
       END DO
       DO jk=1,nlev
          DO jl=1,kproma
             zup1(jl,jk,zjrow)=u(jl,jk,zjrow) +vom(jl,jk,zjrow)*ztmst
             zvp1(jl,jk,zjrow)=v(jl,jk,zjrow) +vol(jl,jk,zjrow)*ztmst

      !       zup(jl,jk,zjrow) = zup1(jl,jk,zjrow) !!!
      !       zvp(jl,jk,zjrow) = zvp1(jl,jk,zjrow) !!!

             ! correct for cos(latitude)
!             zup1(jl,jk,zjrow)= zup1(jl,jk,zjrow) * cor_lat(jl)
!             zvp1(jl,jk,zjrow)= zvp1(jl,jk,zjrow) * cor_lat(jl)

             zup1(jl,jk,zjrow)= zup1(jl,jk,zjrow) * coslat_2D(jl,zjrow) !!!
             zvp1(jl,jk,zjrow)= zvp1(jl,jk,zjrow) * coslat_2D(jl,zjrow) !!!

             zup_cos(jl,jk,zjrow) = zup1(jl,jk,zjrow) !!!
             zvp_cos(jl,jk,zjrow) = zvp1(jl,jk,zjrow) !!!
          END DO
       END DO
    END DO

    ! transform wind with tr_gp_ffsl to ffsl representation

    CALL tr_gp_ffsl (dcl ,1 ,zup1, ub )
    CALL tr_gp_ffsl (dcl ,1 ,zvp1, vb )

    ! collect neighbouring "ghost" latitudes for determination of the meridional gradient

    DO jk=1,klev
       DO j=jfirst,jlast
          DO jl=1,plon
             u2(jl,j,jk) = ub(jl,j,jk)
             v2(jl,j,jk) = vb(jl,j,jk)
          ENDDO
       ENDDO
    ENDDO

#ifndef NOMPI
    CALL ghost_info(u2, plon, jm, klev, nq, jfirst, jlast, ng, ng)
    CALL ghost_info(v2, plon, jm, klev, nq, jfirst, jlast, ng, ng)
#endif

    ! determine zonal gradients
    ! Note that the dependency on cos(latitude) of dx is corrected 
    !      later after retransformation to GP
    dx = radius_earth * 2._dp * pi / plon

    DO j=jfirst,jlast
       DO i=2,plon-1
          dudxb(i,j,1:klev) = &
               (u2(i+1,j,1:klev) - u2(i-1,j,1:klev) ) / (2. * dx)
          dvdxb(i,j,1:klev) = &
               (v2(i+1,j,1:klev) - v2(i-1,j,1:klev) ) / (2. * dx)
       END DO
       dudxb(1,j,1:klev) = &
            (u2(2,j,1:klev) - u2(plon,j,1:klev) ) / (2. * dx)
       dvdxb(1,j,1:klev) = &
            (v2(2,j,1:klev) - v2(plon,j,1:klev) ) / (2. * dx)
       dudxb(plon,j,1:klev) = &
            (u2(1,j,1:klev) - u2(plon-1,j,1:klev) ) / (2. * dx)
       dvdxb(plon,j,1:klev) = &
            (v2(1,j,1:klev) - v2(plon-1,j,1:klev) ) / (2. * dx)
    END DO

    ! determine meridional gradients

    dy = radius_earth * pi / jm

    DO j=jfirst,jlast
       DO i=1,plon
          dudyb(i,j,1:klev) = &
               (u2(i,j+1,1:klev) - u2(i,j-1,1:klev) ) / (2. * dy)
          dvdyb(i,j,1:klev) = &
               (v2(i,j+1,1:klev) - v2(i,j-1,1:klev) ) / (2. * dy)
       END DO
    END DO
   
    ! retransform gradients from ffsl to GP

    CALL tr_gp_ffsl (dcl ,-1 ,dudx       ,dudxb          )
    CALL tr_gp_ffsl (dcl ,-1 ,dvdx       ,dvdxb          )
    CALL tr_gp_ffsl (dcl ,-1 ,dudy       ,dudyb          )
    CALL tr_gp_ffsl (dcl ,-1 ,dvdy       ,dvdyb          )

    ! correct dx for cosine(latitude)
    ! moved outside transformation for an easier access to coslat_2d
    DO zjrow=1, dcl%ngpblks
       IF ( zjrow == dcl%ngpblks ) THEN
          kproma = dcl%npromz
       ELSE
          kproma = dcl%nproma
       END IF
       DO jk=1,nlev
          DO jl=1,kproma
             dudx(jl,jk,zjrow)=dudx(jl,jk,zjrow) / coslat_2d(jl,zjrow)
             dvdx(jl,jk,zjrow)=dvdx(jl,jk,zjrow) / coslat_2d(jl,zjrow)
          END DO
       END DO
    END DO

    DEALLOCATE(ub)   ; NULLIFY(ub)
    DEALLOCATE(vb)   ; NULLIFY(vb)
    DEALLOCATE(dudxb); NULLIFY(dudxb)
    DEALLOCATE(dvdxb); NULLIFY(dvdxb)
    DEALLOCATE(dudyb); NULLIFY(dudyb)
    DEALLOCATE(dvdyb); NULLIFY(dvdyb)
    DEALLOCATE(u2)   ; NULLIFY(u2)
    DEALLOCATE(v2)   ; NULLIFY(v2)

    ! perform calculations for indices in GP

    
    ! determine total wind vector for determination of the vertical wind shear   !!!
    DO zjrow=1, dcl%ngpblks
       IF ( zjrow == dcl%ngpblks ) THEN
          kproma = dcl%npromz
       ELSE
          kproma = dcl%nproma
       END IF

       DO jk=1,nlev
          DO jl=1,kproma
             v_total(jl,jk,zjrow) = sqrt(zup1(jl,jk,zjrow)*zup1(jl,jk,zjrow) + &
	                                 zvp1(jl,jk,zjrow)*zvp1(jl,jk,zjrow) )

          ENDDO
       ENDDO

       DO jk=1,nlev-1
          DO jl=1,kproma
             ! determine layer thickness dz
             dz(jl,jk) = (geopot(jl,jk,zjrow) - geopot(jl,jk+1,zjrow)) / g
          ENDDO
       ENDDO

       DO jk=2,nlev-1
          DO jl=1,kproma
     ! determine vertical wind shear
     ! downwards oriented 
             dvdz_down(jl,jk,zjrow) = (v_total(jl,jk,zjrow) - &
                                       v_total(jl,jk+1,zjrow)) / dz(jl,jk)
     ! centered 
             dvdz_cen(jl,jk,zjrow) = (v_total(jl,jk-1,zjrow) - &
                                      v_total(jl,jk+1,zjrow)) / &
                                    ( dz(jl,jk) + dz(jl,jk-1) )

     ! determine static stability, centered
             dpotdz(jl,jk,zjrow) = (tpoteq(jl,jk-1,zjrow) - &
                                     tpoteq(jl,jk+1,zjrow)) / &
                                   ( dz(jl,jk) + dz(jl,jk-1) )

          END DO
      END DO

     ! calculate (original) Ellrod Index
     DO jk=2,nlev-1
	DO jl=1,kproma

	DVSI_down(jl,jk,zjrow) = sqrt(                                      &
                           (dudx(jl,jk,zjrow) - dvdy(jl,jk,zjrow))  * &
                           (dudx(jl,jk,zjrow) - dvdy(jl,jk,zjrow))  + &
                           (dvdx(jl,jk,zjrow) + dudy(jl,jk,zjrow))  * &
                           (dvdx(jl,jk,zjrow) + dudy(jl,jk,zjrow))) * &
                            dvdz_down(jl,jk,zjrow)  

	DVSI_cen(jl,jk,zjrow) = sqrt(                                      &
                           (dudx(jl,jk,zjrow) - dvdy(jl,jk,zjrow))  * &
                           (dudx(jl,jk,zjrow) - dvdy(jl,jk,zjrow))  + &
                           (dvdx(jl,jk,zjrow) + dudy(jl,jk,zjrow))  * &
                           (dvdx(jl,jk,zjrow) + dudy(jl,jk,zjrow))) * &
                            dvdz_cen(jl,jk,zjrow)

        ENDDO
      ENDDO

      ! calculate divergence trend term
      DO i=2,nti
         time_cnt(i-1)%d_term(:,:,zjrow) = time_cnt(i)%d_term(:,:,zjrow)
      ENDDO

      DO jk=1,nlev
         DO jl=1,kproma
            time_cnt(nti)%d_term(jl,jk,zjrow) = &
                            dudx(jl,jk,zjrow) + dvdy(jl,jk,zjrow)

            DVT(jl,jk,zjrow) = DIV_C * &
                 (time_cnt(nti)%d_term(jl,jk,zjrow) - time_cnt(1)%d_term(jl,jk,zjrow)) 

            ! set negative DVT-values to zero
            IF (DVT(jl,jk,zjrow) < 0._dp) THEN
               DVT(jl,jk,zjrow) = 0._dp
            END IF

         END DO
      END DO

      ! calculate D-DVSI and Brunt-Vaisala-Frequency
      DO jk=2,nlev-1
         DO jl=1,kproma
            DDVSI(jl,jk,zjrow) = DVSI_cen(jl,jk,zjrow) + DVT(jl,jk,zjrow)
            N2(jl,jk,zjrow)    = g / tpoteq(jl,jk,zjrow) * dpotdz(jl,jk,zjrow)
         ENDDO
      ENDDO


      ! calculate with N^2 modified D-DVSI
      DO jk=2,nlev-1
         DO jl=1,kproma
            ! use N2-values lower N2_lim(3.e-4) only
            IF (N2(jl,jk,zjrow) < N2_lim) THEN
                N2_lt(jl,jk) = N2(jl,jk,zjrow)
                DDVSI_MOD(jl,jk,zjrow) = (N2_lim-N2_lt(jl,jk))/N2_lim * &
                                         DDVSI(jl,jk,zjrow)
             ELSE
                DDVSI_MOD(jl,jk,zjrow) = 0._dp
             END IF
         
         ENDDO
      ENDDO

     

      ! exchange of tracer between different layers
      IF (USE_DDVSI) THEN
         DO jk=1,nlev-1
            DO jl=1,kproma
               ! interpolation for DDVSI at box-boundaries
               D_DDVSI(jl,jk,zjrow) = D_norm * &
                    (DDVSI(jl,jk,zjrow) + DDVSI(jl,jk+1,zjrow))/2._dp
               D_DDVSI(jl,jk,zjrow) = MAX(0._dp,D_DDVSI(jl,jk,zjrow))
            END DO
         END DO
      ELSE         
         DO jk=1,nlev-1
            DO jl=1,kproma
               D_DDVSI(jl,jk,zjrow) = D_norm * &
                    (DDVSI_MOD(jl,jk,zjrow) + DDVSI_MOD(jl,jk+1,zjrow))/2._dp
               D_DDVSI(jl,jk,zjrow) = MAX(0._dp,D_DDVSI(jl,jk,zjrow))
            END DO
         END DO
      END IF

   END DO  ! from loop over zjrow=1, dcl%ngpblks

!----------------------------------------------------------------------------
  CONTAINS

    SUBROUTINE ghost_info(q, im, jm, km, nq, jfirst, jlast, nd, ng)

      USE messy_main_mpi_bi, ONLY: p_isend, p_irecv, p_wait
    
      ! No buffer needed.
    
      INTEGER, INTENT(in):: im, jm, km, nq, jfirst, jlast
      INTEGER, INTENT(in):: nd         ! ghost dimension of q
      INTEGER, INTENT(in):: ng         ! zones to be ghosted
                                     ! nd may not be equal to ng if update =.F.
      REAL(dp) :: q(im,jfirst-ng:jlast+ng,km,nq)
      INTEGER :: qsize

!#ifdef OLD_GHOST_UPDATE

      INTEGER :: i, k

      qsize = im*ng

      DO i=1,nq
         DO k=1,km
            IF ( jfirst > 1 ) CALL p_isend(q(1,jfirst    ,k,i),dcl%ffsl%pe_s,1234,qsize)
            IF ( jlast < jm ) CALL p_isend(q(1,jlast-ng+1,k,i),dcl%ffsl%pe_n,9876,qsize)
         ENDDO
      ENDDO
      
      DO i=1,nq
         DO k=1,km
            IF ( jfirst > 1 ) CALL p_irecv(q(1,jfirst-ng,k,i),dcl%ffsl%pe_s,9876,qsize)
            IF ( jlast < jm ) CALL p_irecv(q(1,jlast+1,  k,i),dcl%ffsl%pe_n,1234,qsize)
         ENDDO
      ENDDO
      
      CALL p_wait

!!$!#else
!!$
!!$    REAL(wp) :: q_send1(im,ng,km,nq)
!!$    REAL(wp) :: q_send2(im,ng,km,nq) 
!!$    REAL(wp) :: q_recv1(im,ng,km,nq) 
!!$    REAL(wp) :: q_recv2(im,ng,km,nq)
!!$
!!$    qsize = im*ng*km*nq
!!$
!!$    IF ( jfirst > 1 ) CALL  p_irecv(q_recv1(1,1,1,1),dcl%ffsl%pe_s,9876,qsize) 
!!$    IF ( jlast < jm ) CALL  p_irecv(q_recv2(1,1,1,1),dcl%ffsl%pe_n,1234,qsize)
!!$
!!$    IF ( jfirst > 1 ) THEN
!!$      q_send1(:,1:ng,:,:) = q(:,jfirst:jfirst+ng-1,:,:)
!!$      CALL  p_isend(q_send1(1,1,1,1),dcl%ffsl%pe_s,1234,qsize)
!!$    END IF
!!$
!!$    IF ( jlast < jm ) THEN
!!$      q_send2(:,1:ng,:,:) = q(:,jlast-ng+1:jlast,:,:)
!!$      CALL  p_isend(q_send2(1,1,1,1),dcl%ffsl%pe_n,9876,qsize)
!!$    END IF
!!$
!!$    CALL p_wait
!!$
!!$    IF ( jfirst > 1 )  q(:,jfirst-ng:jfirst-1,:,:) = q_recv1(:,1:ng,:,:)
!!$    IF ( jlast < jm )  q(:,jlast+1:jlast+ng,:,:)   = q_recv2(:,1:ng,:,:)
!!$
!!$#endif

    END SUBROUTINE ghost_info

!--------------------------------------------------------------------------------    

  END SUBROUTINE CAT_GLOBAL_START
!--------------------------------------------------------------------------------

  SUBROUTINE CAT_PHYSC
    
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1, &
                                         ti_gp, ntrac => ntrac_gp
    USE messy_main_grid_def_mem_bi, ONLY: kproma,nlev,jrow
    USE messy_main_data_bi,         ONLY: pressi=>pressi_3d
    USE messy_main_constants_mem,   ONLY: g
    USE messy_main_timer,           ONLY: lstart, ztmst => time_step_len

    INTEGER  :: jl, jk, jt
    REAL(dp) :: xtp1(kproma,nlev,ntrac)
    REAL(dp) :: xtte_cat(kproma,nlev,ntrac)
    REAL(dp) :: inv_ztmst
    REAL(dp) :: delp(kproma,nlev), mean_delp(kproma,nlev)
    REAL(dp) :: sum_bef(kproma), sum_af(kproma)
    
 ! apply tracermixing by CAT with exchange coefficient determined in CAT_GLOBAL_START
    IF (.NOT. l_tracmix) RETURN
    IF (LSTART) RETURN
    xtte_cat = 0._dp
    inv_ztmst = 1._dp/ztmst

    do jk=up_idx, lo_idx+1
       do jl=1,kproma
          delp(jl,jk) = pressi(jl,jk+1,jrow) - pressi(jl,jk,jrow)
       enddo
    enddo
    do jk=up_idx, lo_idx+1
       do jl=1,kproma
          mean_delp(jl,jk) = 0.5_dp * (delp(jl,jk+1) + delp(jl,jk))
       enddo
    enddo

    do jt=1,ntrac
       do jk=up_idx,lo_idx+1
          do jl=1,kproma
               ! updated tracer mixing ratio
             xtp1(jl,jk,jt) = pxtm1(jl,jk,jt) + &
                  pxtte(jl,jk,jt) * ztmst
          end do
         end do

         ! sum before
         sum_bef(:) = 0._dp
         do jk=up_idx,lo_idx+1
            do jl=1,kproma
               sum_bef(jl) = sum_bef(jl) + &
                    xtp1(jl,jk,jt) * delp(jl,jk)/g
            enddo
         enddo


         do jk=up_idx,lo_idx
            DO jl=1,kproma
               ! mixing of tracers due to CAT
               ! mass conserving and homogeneous tracer stays homogeneous
               xtte_cat(jl,jk,jt)   = xtte_cat(jl,jk,jt) + D_DDVSI(jl,jk,jrow) * &
                                    ( (xtp1(jl,jk+1,jt) - xtp1(jl,jk,jt)) * mean_delp(jl,jk)/delp(jl,jk) ) * inv_ztmst
               xtte_cat(jl,jk+1,jt) = xtte_cat(jl,jk+1,jt) + D_DDVSI(jl,jk,jrow) * &
                                    ( (xtp1(jl,jk,jt) - xtp1(jl,jk+1,jt))*mean_delp(jl,jk)/delp(jl,jk+1) ) * inv_ztmst
            end do
         end do

         ! sum_after
         sum_af(:) = 0._dp
         do jk=up_idx,lo_idx+1
            do jl=1,kproma
               sum_af(jl) = sum_af(jl) + &
                    (xtp1(jl,jk,jt)+xtte_cat(jl,jk,jt)*ztmst) * delp(jl,jk)/g
            enddo
         enddo
         do jl=1,kproma
            If (ABS(sum_af(jl)-sum_bef(jl)) > 1.e-15_dp) THEN
               print*, "problem with cat mixing!", sum_af(jl), sum_bef(jl), jt, jl, up_idx, lo_idx+1,&
               " tracer causing problem: before: ", xtp1(jl,up_idx:lo_idx+1,jt), " after: ", &
                    xtp1(jl,up_idx:lo_idx+1,jt) + xtte_cat(jl,up_idx:lo_idx+1,jt) * ztmst, &
                    "tendency: ", xtte_cat(jl,up_idx:lo_idx+1,jt), " delp: ",delp(jl,up_idx:lo_idx+1)
            ENDIF
         END do
         do jk=up_idx,lo_idx+1
            do jl=1,kproma
               ! updated tracer mixing ratio tendency
               pxtte(jl,jk,jt) = pxtte(jl,jk,jt)        + & 
                                 xtte_cat(jl,jk,jt) 
            end do
         end do
      end do
         
            
  END SUBROUTINE CAT_PHYSC

!--------------------------------------------------------------------------------
  SUBROUTINE CAT_FREE_MEMORY
  END SUBROUTINE CAT_FREE_MEMORY
!--------------------------------------------------------------------------------
END MODULE MESSY_CAT_E5
