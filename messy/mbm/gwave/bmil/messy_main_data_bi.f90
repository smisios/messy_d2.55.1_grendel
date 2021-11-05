MODULE messy_main_data_bi

  USE messy_main_grid_def_mem_bi, ONLY: nlevp1, nlev
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE

  PUBLIC 

  CHARACTER(LEN=*), PARAMETER :: modstr = 'GWAVE'
  CHARACTER(LEN=*), PARAMETER :: modver = '1.0'

  LOGICAL, PARAMETER :: l2tls = .FALSE.      

  REAL(dp), DIMENSION(nlev) :: &
       tm1_init, um1_init, vm1_init, geopot_init, aphm1_init


! ncdump -v vct_a,vct_b /mpcdata/projects/modeldata/DATA/ECHAM5/echam5.3.02/init/T42/T42L39MA_19780101_spec.nc
! ncks -v hyam,hybm,um1,vm1,tm1,geopot -d time,1,1 -d lat,40,40 -d lon,0,0 test01_________19780101_0500_ECHAM5.nc var.nc; ncdump var.nc
!#include "data_L39.inc"
#include "data_L90.inc"
!#include "data_L50_Medvedev.inc"

  REAL(dp), DIMENSION(:,:,:), POINTER :: &
       tm1  => NULL(), &
       um1  => NULL(), &
       vm1  => NULL(), &
       press_3d => NULL(), &
       pressi_3d => NULL(), &
       geopot_3d => NULL(), &
       tte_3d   => NULL(), &
       tte_scb  => NULL(), &
       vol_scb  => NULL(), &
       vom_scb  => NULL()

  REAL(dp), DIMENSION(:,:,:), POINTER :: &
       vom_3d  => NULL(), &
       vol_3d  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: &
       tm1_2d  => NULL(), &
       um1_2d  => NULL(), &
       vm1_2d  => NULL(), &
       press => NULL(), &
       pressi => NULL(), &
       geopot => NULL()
  
  REAL(dp), DIMENSION(1,nlevp1) :: aphm1
  REAL(dp), DIMENSION(1,nlev)  :: apm1
      

  ! not used
  REAL(dp), DIMENSION(1,1) :: aprflux, dalpslm1, dalpsmm1, alpsm1
  REAL(dp), DIMENSION(1,1,1) :: vom1, dm1  &
       ,dtlm1, dtmm1, dudlm1, dvdlm1
  
 PUBLIC :: main_data_init_memory

CONTAINS

  !***************************************************************************
  SUBROUTINE main_data_init_memory

    ! MESSy/BMIL
    USE messy_main_grid_def_bi,      ONLY: ceta, hyam, hybm, philat
    USE messy_main_grid_def_mem_bi,  ONLY: apzero
    USE messy_main_channel_error_bi, ONLY: channel_halt
    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_channel_object_reference    &
                                      , new_attribute
    USE messy_main_channel_repr,  ONLY: get_representation_id

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_data_init_memory'
    INTEGER :: status
    INTEGER :: reprid,  repr_int3d   
    INTEGER :: k

    CALL get_representation_id(status, 'GP_3D_MID', reprid)
    CALL channel_halt(substr, status)
    CALL get_representation_id(status, 'GP_3D_INT', repr_int3d)
    CALL channel_halt(substr, status)

    ! create new channel
    CALL new_channel (status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tm1', &
         p3=tm1, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr,  'um1', &
         p3=um1, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr,  'vm1', &
         p3=vm1, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr,  'geopot', &
         p3=geopot_3d, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr,  'tte', &
         p3=tte_scb, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr,  'vom', &
         p3=vom_scb, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    vom_3d => vom_scb
    CALL new_channel_object(status, modstr,  'vol', &
         p3=vol_scb, reprid=reprid, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    vol_3d => vol_scb

    ! pressure at middle of box ("full level pressure")
    CALL new_channel_object(status, modstr,  'press', &
         p3=press_3d, reprid=reprid)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'long_name', c='pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'press', 'units', c='Pa')
    CALL channel_halt(substr, status)
    ! ------------------------------------------------------------------------
    ! pressure at level interfaces ("half level pressure")
    CALL new_channel_object(status, modstr,  'pressi', &
         p3=pressi_3d, reprid=repr_int3d, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi', &
         'long_name', c='interface pressure')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'pressi', 'units', c='Pa')
    CALL channel_halt(substr, status)


    apm1(1,:)=hyam+hybm*apzero
    aphm1(1,1)=0.
    do k=1,nlev
       aphm1(1,k+1)=2* apm1(1,k)-aphm1(1,k)
    enddo
!!$    aphm1(1,1:nlev)=aphm1_init

    ceta=hyam/apzero+hybm

    tm1(1,:,1)=tm1_init
    um1(1,:,1)=um1_init/cos(philat(1)/57.)
    vm1(1,:,1)=vm1_init/cos(philat(1)/57.)
    geopot_3d(1,:,1)=geopot_init

    press_3d(1,:,1)=apm1(1,:)
    pressi_3d(1,:,1)= aphm1(1,:)
    
    tm1_2d=>tm1(:,:,1)
    um1_2d=>um1(:,:,1)
    vm1_2d=>vm1(:,:,1)
    geopot=>geopot_3d(:,:,1)

    press=>press_3d(:,:,1)
    pressi=>pressi_3d(:,:,1)

    tte_3d => tte_scb(:,:,1:1)
     
  END SUBROUTINE main_data_init_memory

END MODULE messy_main_data_bi
