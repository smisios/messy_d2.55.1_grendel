MODULE mo_varlist

#ifndef NO_NEW_IO

  USE mo_kind,     ONLY : dp
  USE mo_param1,   ONLY : ke
  USE mo_commo1,   ONLY : uko,vke,uoo,voe,tiestw
  ! USE mo_commo1,   ONLY : zo,tho,uko,vke,uoo,voe,sao,po,wo,rhoo,sictho &
  !                        ,sicomo,sicuo,sicve,txo,tye,alatv     &
  !                        ,alonv,eminpo,fswr,ftdew,z1o,depto,dlxp      &
  !                        ,dlyp,tafo,avo,dvo,sicsno,alatu,alatu,fu10   &
  !                        ,bolx,boly,fclou,weto,dduo,dlxu,dlyu    &
  !                        ,ddue,dlxv,dlyv,alat,alon,ddpo,deuto,amsue   &
  !                        ,amsuo,deute,wgo,fprec,hibete,hibeto    &
  !                        ,hibzete,hibzeto,rivrun &
  !                        , tiestw


!  USE mo_mean,     ONLY: wtmix,tmceo,tmcdo,rinu
!  USE mo_commoau2, ONLY: prech,qseo,qlao,qswo,qlwo
   USE mo_grid, ONLY: thkcello
!  USE mo_diagnosis, ONLY: psitro,zmld,sictru,sictrv,flum           &
!!       ,sss,sst,sst_sqr,zo_sqr                                    &
!       ,pem!,amld
   USE mo_eddydiag
   USE mo_basin_masks, ONLY: rbek
   USE mo_swr_absorption, ONLY : swsum,swrab,heatabs,swr_frac

  USE mo_parallel, ONLY : stop_all
  !  USE mo_restart, ONLY: rdt
  USE mo_util_string
  USE mo_range_map, range_map => list_type
  USE mo_grid, ONLY: get_level_index_by_depth

#include "pointer_intent_macro.inc"

  IMPLICIT NONE
  PRIVATE
  SAVE

#ifndef NOCDI
  INCLUDE 'cdi.inc'
#endif

  ! Grid types
  character(*), parameter, public :: grid_p = 'p'
  character(*), parameter, public :: grid_p_minus = 'p-'
  character(*), parameter, public :: grid_u = 'u'
  character(*), parameter, public :: grid_u_plus = 'u+'
  character(*), parameter, public :: grid_v = 'v'
  character(*), parameter, public :: grid_v_plus = 'v+'
  character(*), parameter, public :: grid_s = 's'
  character(*), parameter, public :: grid_s_minus = 's-'
  character(*), parameter, public :: grid_single_point = 'g'
  character(*), parameter, public :: grid_zonal_mean = 'zm'

  ! Z axis types
  character(*), parameter :: z_center = 'c'
  character(*), parameter :: z_interface = 'i'
  character(*), parameter :: z_surface = 's'
  character(*), parameter :: z_sediment = 'sed'
  character(*), parameter :: z_ground = 'g'

  ! Reduction operator types.
  character(*), parameter :: reduce_integral = 'i'

  TYPE, PUBLIC :: varlist_element
     CHARACTER(len=128)   :: name
     CHARACTER(len=128)   :: std_name
     CHARACTER(len=32)    :: unit
     INTEGER              :: icode
     REAL(dp), POINTER    :: scalar => NULL()
     REAL(dp), POINTER    :: array_1d(:) => NULL()
     REAL(dp), POINTER    :: array_2d(:,:) => NULL()
     REAL(dp), POINTER    :: array_3d(:,:,:) => NULL()
     INTEGER              :: i_shape(3)
     INTEGER              :: varid
     INTEGER              :: gridid = -1
     CHARACTER(len=5)     :: grid
     INTEGER              :: zaxisid = -1
     CHARACTER(len=7)     :: zaxis
     INTEGER              :: top_level = -1
     INTEGER              :: top_level_index = -1
     INTEGER              :: level = -1
     INTEGER              :: level_index = -1
     INTEGER              :: mask_level_index = -1
     LOGICAL              :: column_integrate = .false.
     LOGICAL              :: surface_integrate = .false.
     REAL(dp)             :: factor = 1._dp
     LOGICAL              :: staggered = .false.
     LOGICAL              :: registered = .false.
  END TYPE varlist_element

  TYPE, PUBLIC :: varlist
     TYPE (varlist_element)  :: vardata
     TYPE (varlist), POINTER :: next ! pointer field
  END TYPE varlist


  INTERFACE new_var
     MODULE PROCEDURE new_var_r
     MODULE PROCEDURE new_var_1d
     MODULE PROCEDURE new_var_2d
     MODULE PROCEDURE new_var_2d_subarr
     MODULE PROCEDURE new_var_3d
     MODULE PROCEDURE new_var_3d_subarr
  END INTERFACE


  TYPE (varlist), POINTER, PUBLIC :: ocean_varlist

  TYPE(range_map) :: zaxis_ids
  INTEGER :: i_zaxisid = -1, g_zaxisid = -1

  INTEGER, PUBLIC :: TsSecOffset, TsCodeperSec
  INTEGER, PUBLIC :: TsRegOffset, TsCodeperReg

  ! Public module procedures must be declared here!

  PUBLIC :: get_var_by_code
  PUBLIC :: print_varlist
  PUBLIC :: build_ocean_varlist
  PUBLIC :: new_var
  PUBLIC :: generate_gridid
  PUBLIC :: generate_zaxisid
  PUBLIC :: varlist_add
  PUBLIC :: varlist_get_element_ref
  PUBLIC :: varlist_is_code_registered
  PUBLIC :: varlist_element_to_string
  PUBLIC :: varlist_post_process_ocean_data

CONTAINS

  SUBROUTINE varlist_add(this,vardata)
    TYPE (varlist), POINTER :: this
    TYPE (varlist_element), INTENT(in) :: vardata
    TYPE (varlist), POINTER :: current

    ALLOCATE(current) ! create new node
    current%vardata=vardata
    CALL element_set_level_info(current%vardata)
    current%next => this
    this => current ! update head of list

  END SUBROUTINE varlist_add

  FUNCTION varlist_get_element_ref(this, code)
    TYPE (varlist_element), POINTER :: varlist_get_element_ref
    TYPE (varlist), POINTERINTENT(IN) :: this
    INTEGER,INTENT(IN) :: code

    TYPE (varlist), POINTER :: current

    current => this
    DO WHILE (ASSOCIATED(current))
       IF ( current%vardata%icode == code) THEN
          varlist_get_element_ref  => current%vardata
          RETURN
       ENDIF
       current => current%next ! make current alias of next varnode
    END DO

    varlist_get_element_ref => NULL()

  END FUNCTION varlist_get_element_ref

  FUNCTION get_var_by_code(this,code)
    TYPE (varlist_element) :: get_var_by_code
    TYPE (varlist), POINTERINTENT(IN) :: this
    INTEGER,INTENT(IN) :: code

    TYPE (varlist_element), POINTER :: element_pointer
    CHARACTER(100) :: message

    element_pointer => varlist_get_element_ref(this, code)

    if(associated(element_pointer)) then
       get_var_by_code = element_pointer
    else
       WRITE(message,*)" STOP: invalid code number in IOCTL", code
       CALL stop_all(message)
    end if

  END FUNCTION get_var_by_code

  function varlist_is_code_registered(this, code)
    TYPE(varlist), POINTERINTENT(IN) :: this
    integer, intent(in) :: code
    logical :: varlist_is_code_registered

    type(varlist_element) :: item

    item = get_var_by_code(this, code)
    varlist_is_code_registered = item%registered

  end function varlist_is_code_registered

  !>
  !! Format list for output.
  !!
  function varlist_element_to_string(this, verbose) result(the_string)

    type(varlist_element), intent(in) :: this
    logical, intent(in), optional :: verbose
    character(2048) :: the_string

    logical :: do_verbose
    integer :: pos, length
    character(100) :: buffer

    ! Default to non-verbose output.
    do_verbose = merge(verbose, .false., present(verbose))

    the_string = ''
    pos = 1

    ! prefix
    if(do_verbose) then
        the_string(pos:pos+14) = 'varlist_element'
        pos = pos + 15
    end if
    the_string(pos:pos) = '('
    pos = pos + 1

    ! %name
    if(do_verbose) then
        the_string(pos:pos+5) = 'name: '
        pos = pos + 6
    end if
    length = len_trim(this%name)
    the_string(pos:pos) = "'"
    the_string(pos+1:pos+length) = this%name(1:length)
    the_string(pos+length+1:pos+length+1) = "'"
    pos = pos + length + 2

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %std_name
    if(do_verbose) then
        the_string(pos:pos+9) = 'std_name: '
        pos = pos + 10
    end if
    length = len_trim(this%std_name)
    the_string(pos:pos) = "'"
    the_string(pos+1:pos+length) = this%std_name(1:length)
    the_string(pos+length+1:pos+length+1) = "'"
    pos = pos + length + 2

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %unit
    if(do_verbose) then
        the_string(pos:pos+5) = 'unit: '
        pos = pos + 6
    end if
    length = len_trim(this%unit)
    the_string(pos:pos) = "'"
    the_string(pos+1:pos+length) = this%unit(1:length)
    the_string(pos+length+1:pos+length+1) = "'"
    pos = pos + length + 2

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %icode
    if(do_verbose) then
        the_string(pos:pos+6) = 'icode: '
        pos = pos + 7
    end if
    write(buffer, '(I0)') this%icode
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %scalar, %array_1d, %array_2d, %array_3d
    if(do_verbose) then
        the_string(pos:pos+11) = 'associated: '
        pos = pos + 12
    end if
    if(associated(this%scalar)) then
       the_string(pos:pos+5) = 'scalar'
       pos = pos + 6
    end if
    if(associated(this%array_1d)) then
       the_string(pos:pos+7) = 'array_1d'
       pos = pos + 8
    end if
    if(associated(this%array_2d)) then
       the_string(pos:pos+7) = 'array_2d'
       pos = pos + 8
    end if
    if(associated(this%array_3d)) then
       the_string(pos:pos+7) = 'array_3d'
       pos = pos + 8
    end if

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %i_shape
    if(do_verbose) then
        the_string(pos:pos+8) = 'i_shape: '
        pos = pos + 9
    end if
    write(buffer, '("[", I0, ", ", I0, ", ", I0, "]")') this%i_shape
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %varid
    if(do_verbose) then
        the_string(pos:pos+6) = 'varid: '
        pos = pos + 7
    end if
    write(buffer, '(I0)') this%varid
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %gridid
    if(do_verbose) then
        the_string(pos:pos+7) = 'gridid: '
        pos = pos + 8
    end if
    write(buffer, '(I0)') this%gridid
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %grid
    if(do_verbose) then
        the_string(pos:pos+5) = 'grid: '
        pos = pos + 6
    end if
    length = len_trim(this%grid)
    the_string(pos:pos) = "'"
    pos = pos + 1
    the_string(pos:pos+length-1) = this%grid(1:length)
    pos = pos + length
    the_string(pos:pos) = "'"
    pos = pos + 1

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %zaxisid
    if(do_verbose) then
        the_string(pos:pos+8) = 'zaxisid: '
        pos = pos + 9
    end if
    write(buffer, '(I0)') this%zaxisid
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %zaxis
    if(do_verbose) then
        the_string(pos:pos+6) = 'zaxis: '
        pos = pos + 7
    end if
    length = len_trim(this%zaxis)
    the_string(pos:pos) = "'"
    pos = pos + 1
    the_string(pos:pos+length-1) = this%zaxis(1:length)
    pos = pos + length
    the_string(pos:pos) = "'"
    pos = pos + 1

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %top_level_index
    if(do_verbose) then
        the_string(pos:pos+10) = 'top_level: '
        pos = pos + 11
    end if
    write(buffer, '(I0)') this%top_level
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %top_level_index
    if(do_verbose) then
        the_string(pos:pos+16) = 'top_level_index: '
        pos = pos + 17
    end if
    write(buffer, '(I0)') this%top_level_index
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %level
    if(do_verbose) then
        the_string(pos:pos+6) = 'level: '
        pos = pos + 7
    end if
    write(buffer, '(I0)') this%level
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %level_index
    if(do_verbose) then
        the_string(pos:pos+12) = 'level_index: '
        pos = pos + 13
    end if
    write(buffer, '(I0)') this%level_index
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %mask_level_index
    if(do_verbose) then
        the_string(pos:pos+17) = 'mask_level_index: '
        pos = pos + 18
    end if
    write(buffer, '(I0)') this%mask_level_index
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %column_integrate
    if(do_verbose) then
        the_string(pos:pos+17) = 'column_integrate: '
        pos = pos + 18
    end if
    buffer = logical_to_string(this%column_integrate)
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %surface_integrate
    if(do_verbose) then
        the_string(pos:pos+18) = 'surface_integrate: '
        pos = pos + 19
    end if
    buffer = logical_to_string(this%surface_integrate)
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %factor
    if(do_verbose) then
        the_string(pos:pos+7) = 'factor: '
        pos = pos + 8
    end if
    write(buffer, '(EN11.3)') this%factor
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %staggered
    if(do_verbose) then
        the_string(pos:pos+10) = 'staggered: '
        pos = pos + 11
    end if
    buffer = logical_to_string(this%staggered)
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! separator
    the_string(pos:pos+1) = ', '
    pos = pos + 2

    ! %registered
    if(do_verbose) then
        the_string(pos:pos+11) = 'registered: '
        pos = pos + 12
    end if
    buffer = logical_to_string(this%registered)
    length = len_trim(buffer)
    the_string(pos:pos+length-1) = buffer(1:length)
    pos = pos + length

    ! suffix
    the_string(pos:pos) = ')'

  end function varlist_element_to_string



  SUBROUTINE print_varlist(filename,list)

    USE mo_parallel,ONLY : p_pe,p_io
    USE mo_io_config, ONLY: next_free_unit

    TYPE (varlist), POINTER :: list, current
    INTEGER :: id
    CHARACTER(*) :: filename

    IF (p_pe == p_io) THEN
       id=next_free_unit()
       OPEN(id,file=TRIM(filename),status='unknown',                  &
            access='sequential',form='formatted')

       ! transverse the list and print the values
       current => list ! make current as alias of list
       DO WHILE ( ASSOCIATED(current) )

          WRITE(id, '("&PARAMETER")')
          WRITE(id, '("  CODE=", I0)') current%vardata%icode
          WRITE(id, '("  NAME=", A)') TRIM(current%vardata%name)
          WRITE(id, '("  LONG_NAME=""", A, """")') &
               capitalize(translate(TRIM(current%vardata%std_name), '_', ' '))
          WRITE(id, '("  UNITS=""", A, """")') TRIM(current%vardata%unit)
          WRITE(id, '("/")')

          current => current%next ! make current alias of next varnode
       END DO

       CLOSE(id)
    ENDIF

  END SUBROUTINE print_varlist



  SUBROUTINE generate_gridid(this)

#ifndef NOCDI
    USE mo_parallel, ONLY : p_io, p_pe, gather, stop_all
    USE mo_param1, ONLY : ie_g, je_g
    USE mo_commo1, ONLY : alat_g, alon_g, alatu, alonu, alatv, alonv, alatpsi_g, alonpsi_g
    REAL(wp), ALLOCATABLE :: alonu_g(:,:)
    REAL(wp), ALLOCATABLE :: alatu_g(:,:)
    REAL(wp), ALLOCATABLE :: alonv_g(:,:)
    REAL(wp), ALLOCATABLE :: alatv_g(:,:)
    REAL(wp) :: zonal_lat(180),zonal_lon(1)
#endif

    TYPE (varlist), POINTERINTENT(IN) :: this
    TYPE (varlist), POINTER :: current

#ifndef NOCDI

    INTEGER :: p_gridid,u_gridid,v_gridid,s_gridid,g_gridid,zm_gridid
    INTEGER :: i

    ! transverse the list

    !    IF ( p_pe == p_io ) WRITE(0,*) 'generate horizontal gridid '

    IF ( p_pe == p_io ) THEN

       ALLOCATE(alonu_g(ie_g,je_g))
       ALLOCATE(alatu_g(ie_g,je_g))
       ALLOCATE(alonv_g(ie_g,je_g))
       ALLOCATE(alatv_g(ie_g,je_g))

    ELSE

       ALLOCATE(alonu_g(0,0))
       ALLOCATE(alatu_g(0,0))
       ALLOCATE(alonv_g(0,0))
       ALLOCATE(alatv_g(0,0))

    ENDIF

    CALL gather(alonu,alonu_g,p_io)
    CALL gather(alatu,alatu_g,p_io)
    CALL gather(alonv,alonv_g,p_io)
    CALL gather(alatv,alatv_g,p_io)


    IF ( p_pe == p_io ) THEN

       p_gridid = gridcreate(grid_curvilinear, (ie_g*je_g))
       CALL griddefxsize(p_gridid, ie_g)
       CALL griddefysize(p_gridid, je_g)
       CALL griddefxvals(p_gridid, alon_g)
       CALL griddefyvals(p_gridid, alat_g)

       u_gridid = gridcreate(grid_curvilinear, (ie_g*je_g))
       CALL griddefxsize(u_gridid, ie_g)
       CALL griddefysize(u_gridid, je_g)
       CALL griddefxvals(u_gridid, alonu_g)
       CALL griddefyvals(u_gridid, alatu_g)

       v_gridid = gridcreate(grid_curvilinear, (ie_g*je_g))
       CALL griddefxsize(v_gridid, ie_g)
       CALL griddefysize(v_gridid, je_g)
       CALL griddefxvals(v_gridid, alonv_g)
       CALL griddefyvals(v_gridid, alatv_g)

       s_gridid = gridcreate(grid_curvilinear, (ie_g*je_g))
       CALL griddefxsize(s_gridid, ie_g)
       CALL griddefysize(s_gridid, je_g)
       CALL griddefxvals(s_gridid, alonpsi_g)
       CALL griddefyvals(s_gridid, alatpsi_g)

!      grid for zonal mean sections (S->N)
       zonal_lat(1)=-89.5_wp
       DO i=1,179
         zonal_lat(i+1) = zonal_lat(i) + 1._wp
       ENDDO
       zonal_lon(1)=0.0_dp

       zm_gridid = gridcreate(grid_lonlat,180)
       CALL griddefxsize(zm_gridid, 1)
       CALL griddefysize(zm_gridid, 180)
       CALL griddefxvals(zm_gridid, zonal_lon)
       CALL griddefyvals(zm_gridid, zonal_lat)


       g_gridid = gridcreate(grid_lonlat,1)
       CALL griddefxsize(g_gridid, 1)
       CALL griddefysize(g_gridid, 1)
       CALL griddefxvals(g_gridid, zonal_lon(1) )
       CALL griddefyvals(g_gridid, zonal_lon(1) )


       current => this
       DO WHILE ( ASSOCIATED(current) )

          SELECT CASE (current%vardata%grid)

          CASE (grid_p, grid_p_minus)

             current%vardata%gridid = p_gridid

          CASE (grid_u, grid_u_plus)

             current%vardata%gridid = u_gridid

          CASE (grid_v, grid_v_plus)

             current%vardata%gridid = v_gridid

          CASE (grid_s, grid_s_minus)

             current%vardata%gridid = s_gridid

          CASE (grid_single_point)

             current%vardata%gridid  = g_gridid

          CASE (grid_zonal_mean)

             current%vardata%gridid  = zm_gridid

          CASE default !! unsupported value

             CALL STOP_ALL(' generate_gridid : argument unsupported ')


          END SELECT

          current => current%next ! make current alias of next varnode

       END DO

    ENDIF
    DEALLOCATE (alonu_g,alatu_g,alonv_g,alatv_g)


#endif

  END SUBROUTINE generate_gridid

  !>
  !! Create Z axis information for output.
  !!
  SUBROUTINE generate_zaxisid(this)

#ifndef NOCDI

    USE mo_parallel, ONLY : p_io, p_pe
    USE mo_commo1, ONLY : tiestu

#endif


    TYPE (varlist), POINTERINTENT(IN) :: this
    TYPE (varlist), POINTER :: current
    TYPE (varlist_element), POINTER :: var

    INTEGER :: id
    ! @todo This should not be necessary
    ! @todo but for now unallocated variables also need a grid...
    INTEGER :: s_zaxisid

#ifndef NOCDI

    IF ( p_pe == p_io ) THEN

       ! @todo This should not be necessary
       ! @todo but for now unallocated variables also need a grid...
       IF(.NOT. range_map_get(zaxis_ids, 0, 0, id)) THEN
          id = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, 1)
          CALL zaxisDefLevels(id, (/ 0.0_dp /))
          ! Don't forget to store new z axis for future use.
          CALL range_map_add(zaxis_ids, 0, 0, id)
       END IF
       s_zaxisid = id

       ! Traverse the given variable definition list.
       current => this
       DO WHILE ( ASSOCIATED(current) )
          var => current%vardata

          ! Skip unallocated variables.
          if(var%i_shape(3) > 0) then

             IF(var%zaxis == z_center) THEN

                if(.not. range_map_get(zaxis_ids, int(tiestu(1)), &
                                       int(tiestu(var%i_shape(3))), id)) then
                   id = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, var%i_shape(3))
                   call zaxisDefLevels(id, tiestu(1:var%i_shape(3)))
                   call range_map_add(zaxis_ids, int(tiestu(1)), &
                                      int(tiestu(var%i_shape(3))), id)
                end if
                var%zaxisid = id

             ELSE IF(var%zaxis == z_interface) THEN

                if(i_zaxisid == -1) then
                   i_zaxisid = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, ke+1)
                   call zaxisDefLevels(i_zaxisid, tiestw(1:ke+1))
                end if
                var%zaxisid = i_zaxisid

             ELSE IF(var%level >= 0 .or. var%top_level >= 0) THEN
                ! With the 's=nnn' specification,
                ! we create a surface at a given level below surface.

                ! If we do no have a z axis at the given level, create one.
                IF(.NOT. range_map_get(zaxis_ids, &
                                       var%top_level, var%level, id)) THEN
                   if( var%top_level == var%level ) then
                      id = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, 1)
                      CALL zaxisDefLevels(id, (/ real(var%level, dp) /))
                   else
                      id = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, var%i_shape(3))
                      CALL zaxisDefLevels(id, &
                           tiestu(var%top_level_index:var%level_index))
                   end if
                   ! Don't forget to store new z axis for future use.
                   CALL range_map_add(zaxis_ids, &
                                      var%top_level, var%level, id)
                END IF
                var%zaxisid = id

             ELSE IF(var%zaxis == z_sediment) THEN

                ! Z axis for sediment layers is set separately

             ELSE IF(var%zaxis == z_ground) THEN

                if(g_zaxisid == -1) then
                   g_zaxisid = zaxisCreate(ZAXIS_DEPTH_BELOW_SEA, 1)
                   call zaxisDefLevels(g_zaxisid, (/ 0.0_dp /))
                end if
                var%zaxisid  = g_zaxisid

             ELSE !! unsupported value

                write(0, '("generate_zaxisid: ", A)') &
                     trim(varlist_element_to_string(var, .true.))
                CALL STOP_ALL('generate_zaxisid : argument unsupported ')


             END IF

          else
             ! @todo This should not be necessary
             ! @todo but for now unallocated variables also need a grid...
             var%zaxisid = s_zaxisid
          end if ! ishape(3) > 0

          current => current%next ! make current alias of next varnode
       END DO ! WHILE ( ASSOCIATED(current) )

    ENDIF ! p_pe == p_io

#endif

  END SUBROUTINE generate_zaxisid



  !>
  !! Compute level information from data fields.
  !!
  !! Z axis    2D                         3D
  !! --------- -----------------------    ------------------------------
  !! 's' means
  !! 's=0'     var, mask(0)               var(0), mask(0)
  !! 's=100'   var, mask(100)             var(100), mask(100)
  !! 'si' means
  !! 'si>0'    error                      sum(var(:)*thick(:)), mask(0)
  !! 'si>100'  error                      sum(var(100:)*thick(100:), mask(100)
  !! 'si<100'  error                      sum(var(:100)*thick(:100)), mask(0)

  subroutine element_set_level_info(this)

     type(varlist_element), intent(inout) :: this

     integer :: pos, level, level_index
     character(10) :: reduce_op
     character(10) :: select_op

     ! Evaluate Z axis information.

     ! Check prefix for surface specification, e.g. 's'.
     ! Note that this might also be 'sed' for sediment,
     ! that is dealt with later.
     if(this%zaxis(:len(z_surface)) == z_surface) then

        level = 0
        pos = len(z_surface)
        reduce_op = ''
        select_op = ''

        ! Check for reduction operator, e.g. 'si'
        ! Currently only integral is supported.
        if(this%zaxis(pos+1:pos+len(reduce_integral)) == reduce_integral) then
            reduce_op = reduce_integral
            pos = pos + len(reduce_integral)
        end if
        ! Check for selection operator.
        ! Currently only single level (=) or levels from surface (<) supported.
        if(scan(this%zaxis(pos+1:pos+1), '<=>') > 0) then
            select_op = this%zaxis(pos+1:pos+1)
            pos = pos + 1
            if(verify(this%zaxis(pos+1:), '0123456789 ') == 0) then
                read(this%zaxis(pos+1:), '(I5)') level
                pos = len(this%zaxis)
            else
                ! No valid selection operator, back to square one.
                pos = 1
            end if
        end if

        ! Check if we are finished.
        ! Other characters mean we have to skip to next check.
        if(this%zaxis(pos+1:) == '') then

            level_index = get_level_index_by_depth(real(level, dp))

            this%top_level_index = 1
            this%level_index = this%i_shape(3)

            if(this%i_shape(3) == 1) then
                ! Scalar, 1D, or 2D case.
                ! Specified level is used for masking only.
                this%top_level = level
                this%level = level
                this%mask_level_index = level_index
            else if(this%i_shape(3) > 1) then
                ! 3D case. Check for reduction/selection.
                ! Set default for selection wrt reduction
                if(select_op == '') then
                    if(reduce_op == reduce_integral) then
                        select_op = '>'
                    else
                        select_op = '='
                    end if
                end if
                if(select_op == '>' .or. select_op == '=') then
                    this%top_level = level
                    this%top_level_index = level_index
                end if
                if(select_op == '<' .or. select_op == '=') then
                    this%level = level
                    this%level_index = level_index
                end if
                this%i_shape(3) = this%top_level_index - this%level_index + 1
                ! Check for reduction operation.
                if(reduce_op == reduce_integral) then
                    this%column_integrate = .true.
                    ! Reduce 3D field
                    this%i_shape(3) = 1
                    ! Collapse to one level
                    if(select_op == '>') then
                        this%level = level
                    else if(select_op == '<') then
                        this%top_level = 0
                        this%level = 0
                    end if
                    ! Masking is defined by top level
                    this%mask_level_index = this%top_level_index
                end if
            end if ! i_shape(3) > 1

        end if ! this%zaxis(pos+1:) == '' .and. all(this%i_shape(1:2) > 1)

     end if ! z_surface

     if(this%zaxis == z_sediment) then
        ! Sediment layers are available for every depth,
        ! so we may only use the surface mask.
        this%mask_level_index = 1
     end if

     ! Evaluate grid information.
     if( this%grid == grid_single_point .and. &
         (this%i_shape(2)>1 .or. this%i_shape(1)>1) ) then
        ! Reduce 2D or 1D surface to single point
        ! or a 3D volume to a single line.
        this%surface_integrate = .true.
        this%i_shape(1) = 1
        this%i_shape(2) = 1
     end if

  end subroutine element_set_level_info


  FUNCTION new_var_r(name, std_name, unit, icode, field, grid, zaxis, factor, &
                     staggered)

    TYPE (varlist_element) :: new_var_r

    CHARACTER(len=*),       INTENT(in)    :: name
    CHARACTER(len=*),       INTENT(in)    :: std_name
    CHARACTER(len=*),       INTENT(in)    :: unit
    INTEGER,                INTENT(in)    :: icode
    REAL(dp),TARGET,        INTENT(in)    :: field
    CHARACTER(len=*),       INTENT(in)    :: grid
    CHARACTER(len=*),       INTENT(in)    :: zaxis
    real(dp), optional, intent(in) :: factor
    logical, optional, intent(in) :: staggered

    new_var_r%name     = name
    new_var_r%std_name = std_name
    new_var_r%unit     = unit
    new_var_r%icode    = icode
    new_var_r%grid     = grid
    new_var_r%zaxis    = zaxis
    if(present(factor)) new_var_r%factor = factor
    if(present(staggered)) new_var_r%staggered = staggered

    new_var_r%scalar => field
    new_var_r%i_shape(:) = 1

    new_var_r%array_1d => NULL()
    new_var_r%array_2d => NULL()
    new_var_r%array_3d => NULL()

  END FUNCTION new_var_r

  FUNCTION new_var_1d(name, std_name, unit, icode, field, grid, zaxis, factor, &
                      staggered)

    TYPE (varlist_element) :: new_var_1d

    CHARACTER(len=*),       INTENT(in)    :: name
    CHARACTER(len=*),       INTENT(in)    :: std_name
    CHARACTER(len=*),       INTENT(in)    :: unit
    INTEGER,                INTENT(in)    :: icode
    REAL(dp),ALLOCATABLE,TARGET,INTENT(in) :: field(:)
    CHARACTER(len=*),       INTENT(in)    :: grid
    CHARACTER(len=*),       INTENT(in)    :: zaxis
    real(dp), optional, intent(in) :: factor
    logical, optional, intent(in) :: staggered

    new_var_1d%name     = name
    new_var_1d%std_name = std_name
    new_var_1d%unit     = unit
    new_var_1d%icode    = icode
    new_var_1d%grid     = grid
    new_var_1d%zaxis    = zaxis
    if(present(factor)) new_var_1d%factor = factor
    if(present(staggered)) new_var_1d%staggered = staggered

    IF (ALLOCATED(field)) THEN
       new_var_1d%array_1d => field
       new_var_1d%i_shape(1) = SIZE(field,1)
       new_var_1d%i_shape(2:3) = 1
    ELSE
       new_var_1d%array_1d => NULL()
       new_var_1d%i_shape(:) = 0
    ENDIF

    new_var_1d%array_2d => NULL()
    new_var_1d%array_3d => NULL()
    new_var_1d%scalar => NULL()

  END FUNCTION new_var_1d

  FUNCTION new_var_2d(name, std_name, unit, icode, field, grid, zaxis, factor, &
                      staggered)

    TYPE (varlist_element) :: new_var_2d

    CHARACTER(len=*),       INTENT(in)    :: name
    CHARACTER(len=*),       INTENT(in)    :: std_name
    CHARACTER(len=*),       INTENT(in)    :: unit
    INTEGER,                INTENT(in)    :: icode
    REAL(dp),ALLOCATABLE,TARGET,INTENT(in) :: field(:,:)
    CHARACTER(len=*),       INTENT(in)    :: grid
    CHARACTER(len=*),       INTENT(in)    :: zaxis
    real(dp), optional, intent(in) :: factor
    logical, optional, intent(in) :: staggered

    new_var_2d%name     = name
    new_var_2d%std_name = std_name
    new_var_2d%unit     = unit
    new_var_2d%icode    = icode
    new_var_2d%grid     = grid
    new_var_2d%zaxis    = zaxis
    if(present(factor)) new_var_2d%factor = factor
    if(present(staggered)) new_var_2d%staggered = staggered

    IF (ALLOCATED(field)) THEN
       new_var_2d%array_2d => field
       new_var_2d%i_shape(1) = SIZE(field,1)
       new_var_2d%i_shape(2) = SIZE(field,2)
       new_var_2d%i_shape(3) = 1
    ELSE
       new_var_2d%array_2d => NULL()
       new_var_2d%i_shape(:) = 0
    ENDIF

    new_var_2d%array_3d => NULL()
    new_var_2d%array_1d => NULL()
    new_var_2d%scalar => NULL()

  END FUNCTION new_var_2d

  FUNCTION new_var_2d_subarr(name, std_name, unit, icode, field, itracer, &
       grid, zaxis, factor, staggered)

    TYPE (varlist_element) :: new_var_2d_subarr

    CHARACTER(len=*),       INTENT(in)    :: name
    CHARACTER(len=*),       INTENT(in)    :: std_name
    CHARACTER(len=*),       INTENT(in)    :: unit
    INTEGER,                INTENT(in)    :: icode
    INTEGER,                INTENT(in)    :: itracer
    REAL(dp),ALLOCATABLE,TARGET,INTENT(in) :: field(:,:,:)
    CHARACTER(len=*),       INTENT(in)    :: grid
    CHARACTER(len=*),       INTENT(in)    :: zaxis
    real(dp), optional, intent(in) :: factor
    logical, optional, intent(in) :: staggered

    new_var_2d_subarr%name     = name
    new_var_2d_subarr%std_name = std_name
    new_var_2d_subarr%unit     = unit
    new_var_2d_subarr%icode    = icode
    new_var_2d_subarr%grid     = grid
    new_var_2d_subarr%zaxis    = zaxis
    if(present(factor)) new_var_2d_subarr%factor = factor
    if(present(staggered)) new_var_2d_subarr%staggered = staggered

    IF (ALLOCATED(field) .AND. itracer /= 0 ) THEN
       new_var_2d_subarr%array_2d => field(:,:,itracer)
       new_var_2d_subarr%i_shape(1) = SIZE(field,1)
       new_var_2d_subarr%i_shape(2) = SIZE(field,2)
       new_var_2d_subarr%i_shape(3) = 1
    ELSE
       new_var_2d_subarr%array_2d => NULL()
       new_var_2d_subarr%i_shape(:) = 0
    ENDIF

    new_var_2d_subarr%array_3d => NULL()
    new_var_2d_subarr%array_1d => NULL()
    new_var_2d_subarr%scalar => NULL()

  END FUNCTION new_var_2d_subarr

  FUNCTION new_var_3d(name, std_name, unit, icode, field, grid, zaxis, factor, &
                      staggered)

    TYPE (varlist_element) :: new_var_3d

    CHARACTER(len=*),       INTENT(in)    :: name
    CHARACTER(len=*),       INTENT(in)    :: std_name
    CHARACTER(len=*),       INTENT(in)    :: unit
    INTEGER,                INTENT(in)    :: icode
    REAL(dp),ALLOCATABLE,TARGET,INTENT(in) :: field(:,:,:)
    CHARACTER(len=*),       INTENT(in)    :: grid
    CHARACTER(len=*),       INTENT(in)    :: zaxis
    real(dp), optional, intent(in) :: factor
    logical, optional, intent(in) :: staggered

    new_var_3d%name     = name
    new_var_3d%std_name = std_name
    new_var_3d%unit     = unit
    new_var_3d%icode    = icode
    new_var_3d%grid     = grid
    new_var_3d%zaxis    = zaxis
    if(present(factor)) new_var_3d%factor = factor
    if(present(staggered)) new_var_3d%staggered = staggered

    IF (ALLOCATED(field)) THEN
       new_var_3d%array_3d => field
       new_var_3d%i_shape(1) = SIZE(field,1)
       new_var_3d%i_shape(2) = SIZE(field,2)
       new_var_3d%i_shape(3) = SIZE(field,3)
    ELSE
       new_var_3d%array_3d => NULL()
       new_var_3d%i_shape(:) = 0

    ENDIF

    new_var_3d%array_2d => NULL()
    new_var_3d%array_1d => NULL()
    new_var_3d%scalar => NULL()

  END FUNCTION new_var_3d

  FUNCTION new_var_3d_subarr(name, std_name, unit, icode, field, itracer, &
       grid, zaxis, factor, staggered)

    TYPE (varlist_element) :: new_var_3d_subarr

    CHARACTER(len=*),       INTENT(in)    :: name
    CHARACTER(len=*),       INTENT(in)    :: std_name
    CHARACTER(len=*),       INTENT(in)    :: unit
    INTEGER,                INTENT(in)    :: icode
    INTEGER,                INTENT(in)    :: itracer
    REAL(dp),ALLOCATABLE,TARGET,INTENT(in) :: field(:,:,:,:)
    CHARACTER(len=*),       INTENT(in)    :: grid
    CHARACTER(len=*),       INTENT(in)    :: zaxis
    real(dp), optional, intent(in) :: factor
    logical, optional, intent(in) :: staggered

    new_var_3d_subarr%name     = name
    new_var_3d_subarr%std_name = std_name
    new_var_3d_subarr%unit     = unit
    new_var_3d_subarr%icode    = icode
    new_var_3d_subarr%grid     = grid
    new_var_3d_subarr%zaxis    = zaxis
    if(present(factor)) new_var_3d_subarr%factor = factor
    if(present(staggered)) new_var_3d_subarr%staggered = staggered

    IF (ALLOCATED(field) .AND. itracer /= 0 ) THEN
       new_var_3d_subarr%array_3d => field(:,:,:,itracer)
       new_var_3d_subarr%i_shape(1) = SIZE(field,1)
       new_var_3d_subarr%i_shape(2) = SIZE(field,2)
       new_var_3d_subarr%i_shape(3) = SIZE(field,3)
    ELSE
       new_var_3d_subarr%array_3d => NULL()
       new_var_3d_subarr%i_shape(:) = 0
    ENDIF

    new_var_3d_subarr%array_2d => NULL()
    new_var_3d_subarr%array_1d => NULL()
    new_var_3d_subarr%scalar => NULL()

  END FUNCTION new_var_3d_subarr

  SUBROUTINE build_ocean_varlist

    USE mo_parallel, ONLY : p_pe, p_io
    USE mo_commo1
    USE mo_diagnosis
    USE mo_fluxes1

    IMPLICIT NONE

    INTEGER n,ncode

    ! build up the list

    NULLIFY(ocean_varlist) ! initially nullify list (empty)

    CALL varlist_add(ocean_varlist,new_var('zo','sea_surface_height_above_sea_level','m',1,zo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('tho','sea_water_potential_temperature','C',2,tho,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('uko','sea_water_x_velocity','m s-1',3,uko,'u','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('vke','sea_water_y_velocity','m s-1',4,vke,'v','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('sao','sea_water_salinity','psu',5,sao,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('po','sea_water_pressure','pa',6,po,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('wo','upward_sea_water_velocity','m s-1',7,wo,'p','i'))
    CALL varlist_add(ocean_varlist,new_var('rhoo','sea_water_density','kg m-3',8,rhoo,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('uaccel','sea_water_x_acceleration','m s-2',9,uaccel,'u','c'))
    CALL varlist_add(ocean_varlist,new_var('vaccel','sea_water_y_acceleration','m s-2',10,vaccel,'v','c'))
    CALL varlist_add(ocean_varlist,new_var('zo_sqr','square_of_sea_surface_height_above_sea_level','m2',11,zo_sqr,'p','s')) ! cmip5
    CALL varlist_add(ocean_varlist,new_var('sst','sea_surface_temperature','K',12,sst,'p','s')) ! cmip5
    CALL varlist_add(ocean_varlist,new_var('sictho','sea_ice_thickness','m',13,sictho,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('sst_sqr','square_of_sea_surface_temperature','K2',14,sst_sqr,'p','s')) ! cmip5
    CALL varlist_add(ocean_varlist,new_var('sicomo','sea_ice_area_fraction','1',15,sicomo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('sss','sea_surface_salinity','C',16,sss,'p','s')) ! cmip5
    ! cmip5
    CALL varlist_add(ocean_varlist,new_var('bottom_pressure','sea_water_pressure_at_sea_floor','dbar',17,bottom_pressure,'p','g'))
    CALL varlist_add(ocean_varlist,new_var('rhopoto','sea_water_potential_density','kg m-3',18,rhopoto,'p','c')) ! cmip5

    CALL varlist_add(ocean_varlist,new_var('wmo','upward_ocean_mass_transport','kg s-1',21,wmo,'p','i')) ! cmip5
    CALL varlist_add(ocean_varlist,new_var('wmosq','square_of_upward_ocean_mass_transport','kg2 s-2',22,wmosq,'p','i')) ! cmip5
    CALL varlist_add(ocean_varlist,new_var('umo','ocean_mass_x_transport','kg s-1',23,umo,'u','c',staggered=.true.)) ! cmip5
    CALL varlist_add(ocean_varlist,new_var('vmo','ocean_mass_y_transport','kg s-1',24,vmo,'v','c',staggered=.true.)) ! cmip5

    CALL varlist_add(ocean_varlist,new_var('psitro','ocean_barotropic_mass_streamfunction','kg s-1',27,psitro,'s','s'))

    CALL varlist_add(ocean_varlist,new_var('sicuo','sea_ice_x_velocity','m s-1',35,sicuo,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('sicve','sea_ice_y_velocity','m s-1',36,sicve,'v','s'))
    ! cmip5 ; same as sicuo, but masked if sea ice thickness < 1e-3 m
    CALL varlist_add(ocean_varlist,new_var('usi','zonal_sea_ice_velocity','m s-1',37,usi,'u','s'))
    ! cmip5 ; same as sicve, but masked if sea ice thickness < 1e-3 m
    CALL varlist_add(ocean_varlist,new_var('vsi','meridional_sea_ice_velocity','m s-1',38,vsi,'v','s'))

    CALL varlist_add(ocean_varlist,new_var('tauwatu','Eastward_Ocean_Stress_On_Sea_Ice','N m-2',50,tauwatu,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('tauwatv','Northward_Ocean_Stress_On_Sea_Ice','N m-2',51,tauwatv,'v','s'))
    CALL varlist_add(ocean_varlist,new_var('txo','downward_eastward_momentum_flux_in_air','pa kg-1 m-3',52,txo,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('tye','downward_northward_momentum_flux_in_air','pa kg-1 m-3',53,tye,'v','s'))
    CALL varlist_add(ocean_varlist,new_var('alatv','latitude_at_v_vector_point','degree',56,alatv,'v','s'))
    CALL varlist_add(ocean_varlist,new_var('alonv','longitude_at_v_vector_point','degree',57,alonv,'v','s'))
    CALL varlist_add(ocean_varlist,new_var('prech','water_flux_in_ocean_without_flux_correction','m s-1',65,prech,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('eminpo','water_flux_correction','m s-1',67,eminpo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('flum','surface_net_downward_heat_flux_where_sea','w m-2',70,flum,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('pem','water_flux_into_ocean','m s-1',79,pem,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('fswr','downwelling_shortwave_flux_in_air','w m-2',80,fswr,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('ftdew','dew_point_temperature','k',81,ftdew,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('z1o','sea_surface_height_above_sea_level_change','m',82,z1o,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('depto','depth_at_pressure_point','m',84,depto,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('dlxp','grid_x_distance_at_pressure_point','m',85,dlxp,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('dlyp','grid_y_distance_at_pressure_point','m',86,dlyp,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('tafo','air_temperature','c',92,tafo,'p','s'))

!   (zm) zonal mean sections from south to north
    CALL varlist_add(ocean_varlist,new_var('global_hfl','global_ocean_implied_heattransport'                           &
         ,'W',93,global_hfl,'zm','s'))
    CALL varlist_add(ocean_varlist,new_var('atlantic_hfl','atlantic_ocean_implied_heattransport'                       &
         ,'W',94,atlantic_hfl,'zm','s'))
    CALL varlist_add(ocean_varlist,new_var('indopacific_hfl','indopacific_ocean_implied_heattransport'                 &
         ,'W',95,indopacific_hfl,'zm','s'))

    CALL varlist_add(ocean_varlist,new_var('global_wfl','global_ocean_implied_freshwater_transport'                           &
         ,'m3',96,global_wfl,'zm','s'))
    CALL varlist_add(ocean_varlist,new_var('atlantic_wfl','atlantic_ocean_implied_freshwater_transport'                       &
         ,'m3',97,atlantic_wfl,'zm','s'))
    CALL varlist_add(ocean_varlist,new_var('indopacific_wfl','indopacific_ocean_implied_freshwater_transport'                 &
         ,'m3',98,indopacific_wfl,'zm','s'))

    CALL varlist_add(ocean_varlist,new_var('global_moc','global_ocean_meridional_overturning_streamfunction'           &
         ,'kg s-1',100,global_moc,'zm','i'))
    CALL varlist_add(ocean_varlist,new_var('atlantic_moc','atlantic_ocean_meridional_overturning_streamfunction'       &
         ,'kg s-1',101,atlantic_moc,'zm','i'))
    CALL varlist_add(ocean_varlist,new_var('indopacific_moc','indopacific_ocean_meridional_overturning_streamfunction' &
         ,'kg s-1',102,indopacific_moc,'zm','i'))

    CALL varlist_add(ocean_varlist,new_var('rivrun','runoff_flux',' m s-1',105,rivrun,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('avo','sea_water_vertical_viscosity','m2 s-1',110,avo,'p','i'))
    CALL varlist_add(ocean_varlist,new_var('dvo','sea_water_vertical_diffusivity','m2 s-1',111,dvo,'p','i'))
    !  CALL varlist_add(ocean_varlist,new_var('sea_water_vertical_diffusivity_due_to_wind','m2 s-1',112,wtmix,'p'))

    CALL varlist_add(ocean_varlist,new_var('uih0','product_of_sea_ice_x_velocity_and_sea_ice_thickness' &
         ,'m2 s-1',113,uih0,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('vih0','product_of_sea_ice_y_velocity_and_sea_ice_thickness' &
         ,'m2 s-1',114,vih0,'v','s'))

    CALL varlist_add(ocean_varlist,new_var('uic0','product_of_sea_ice_x_velocity_and_sea_ice_concentration' &
         ,'m s-1',115,uic0,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('vic0','product_of_sea_ice_y_velocity_and_sea_ice_concentration' &
         ,'m s-1',116,vic0,'v','s'))

    CALL varlist_add(ocean_varlist,new_var('uisn0','product_of_sea_ice_x_velocity_and_snow_thickness' &
         ,'m2 s-1',117,uisn0,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('visn0','product_of_sea_ice_y_velocity_and_snow_thickness' &
         ,'m2 s-1',118,visn0,'v','s'))

    CALL varlist_add(ocean_varlist,new_var('tt0','square_of_sea_water_potential_temperature' &
         ,'C2',119,tt0,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('ss0','square_of_sea_water_salinity' &
         ,'psu2',120,tt0,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('pp0','square_of_sea_water_pressure' &
         ,'pa2',121,pp0,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('rr0','square_of_sea_water_density' &
         ,'kg2 m-6',122,rr0,'p','c'))

!   123
!   ..
!   133

    CALL varlist_add(ocean_varlist,new_var('swsum','total_penetrating_fraction_of_solar_flux','frac',134,swsum,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('swrab','penetrating_fraction_of_solar_flux_with_depth','frac',135,swrab,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('heatabs','total_subsurface_heating_due_to_absorption','J m-2',136,heatabs,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('swr_frac','available_light_fraction','frac',137,swr_frac,'p','c'))

    CALL varlist_add(ocean_varlist,new_var('transix','sea_ice_x_transport','kg s-1',138,transix,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('transiy','sea_ice_y_transport','kg s-1',139,transiy,'v','s'))
    CALL varlist_add(ocean_varlist,new_var('sicsno','surface_snow_thickness_where_sea_ice','m',141,sicsno,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('sictru','sea_ice_x_transport','kg s-1',142,sictru,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('sictrv','sea_ice_y_transport','kg s-1',143,sictrv,'v','s'))
    CALL varlist_add(ocean_varlist,new_var('qseo','surface_downward_sensible_heat_flux_where_sea','w m-2',146,qseo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('qlao','surface_downward_latent_heat_flux_where_sea','w m-2',147,qlao,'p','s'))
    !  CALL varlist_add(ocean_varlist,new_var('surface_downward_net_minus_shortwave_flux_where_sea','w m-2',148,qnet)
    CALL varlist_add(ocean_varlist,new_var('alatu','latitude_at_u_vector_point','degree',154,alatu,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('alonu','longitude_at_u_vector_point','degree',155,alonu,'u','s'))
    !  CALL varlist_add(ocean_varlist,new_var('tmceo', &
    !      'tendency_of_ocean_potential_energy_content_due_to_diffusion', &
    !      'w m-2',157,tmceo,'p','s'))
    !  CALL varlist_add(ocean_varlist,new_var('tmcdo', &
    !      'tendency_of_ocean_potential_energy_content_due_to_convection', &
    !      'w m-2',158,tmcdo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('bolx','gm_x_coefficient ','m2 s-1',159,bolx,'u+','s'))
    CALL varlist_add(ocean_varlist,new_var('boly','gm_y_coefficient','m2 s-1',160,boly,'v+','s'))
    CALL varlist_add(ocean_varlist,new_var('fclou','cloud_area_fraction','1',164,fclou,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('fu10','10m_wind_velocity','m s-1',171,fu10,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('weto','sea_binary_mask_at_pressure_point','1',172,weto,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('dlxpsi','grid_x_distance_at_psi_point','m',174,dlxpsi,'s','s'))
    CALL varlist_add(ocean_varlist,new_var('dlypsi','grid_y_distance_at_psi_point','m',175,dlypsi,'s','s'))

    CALL varlist_add(ocean_varlist,new_var('qswo','surface_net_downward_shortwave_flux_where_sea'                    &
         ,'w m-2',176,qswo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('qlwo','surface_net_downward_longwave_flux_where_sea'                     &
         ,'w m-2',177,qlwo,'p','s'))

!178
!..
!180

    CALL varlist_add(ocean_varlist,new_var('zmld_sqr','square_of_ocean_mixed_layer_thickness_defined_by_sigma_t'     &
         ,'m2',181,zmld_sqr,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('amld','maxium_ocean_mixed_layer_thickness','m',182,amld,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('zmld','ocean_mixed_layer_thickness_defined_by_sigma_t'                   &
         ,'m',183,zmld,'p','s'))

    CALL varlist_add(ocean_varlist,new_var('dduo','ocean_level_thickness_at_u_vector_point','m',184,dduo,'u+','c'))
    CALL varlist_add(ocean_varlist,new_var('dlxu','grid_x_distance_at_u_vector_point','m',185,dlxu,'u+','s'))
    CALL varlist_add(ocean_varlist,new_var('dlyu','grid_y_distance_at_u_vector_point','m',186,dlyu,'u+','s'))
     CALL varlist_add(ocean_varlist,new_var('ddue','ocean_level_thickness_at_v_vector_point','m',187,ddue,'v+','c'))
    CALL varlist_add(ocean_varlist,new_var('dlxv','grid_x_distance_at_v_vector_point','m',188,dlxv,'v+','s'))
    CALL varlist_add(ocean_varlist,new_var('dlyv','grid_y_distance_at_v_vector_point','m',189,dlyv,'v+','s'))
    CALL varlist_add(ocean_varlist,new_var('alat','latitude_at_pressure_point','degree',190,alat,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('alon','longitude_at_pressure_point','degree',191,alon,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('ddpo','ocean_level_thickness_at_pressure_point','m',192,ddpo,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('deuto','depth_at_u_vector_point','m',193,deuto,'u+','s'))
    CALL varlist_add(ocean_varlist,new_var('amsue','sea_binary_mask_at_v_vector_point','1',194,amsue,'v+','c'))
    CALL varlist_add(ocean_varlist,new_var('amsuo','sea_binary_mask_at_u_vector_point','1',195,amsuo,'u+','c'))
    CALL varlist_add(ocean_varlist,new_var('deute','depth_at_v_vector_point','m',196,deute,'v+','s'))
    CALL varlist_add(ocean_varlist,new_var('thkcello','cell_thickness','m',197,thkcello,'p','c'))
    CALL varlist_add(ocean_varlist,new_var('wgo','gm_upward_sea_water_velocity','m s-1',207,wgo,'p','i'))
    CALL varlist_add(ocean_varlist,new_var('fprec','precipitation_flux','m s-1',213,fprec,'p','s'))
    !  CALL varlist_add(ocean_varlist,new_var('richardson_number','1',214,rinu,'p'))

    CALL varlist_add(ocean_varlist,new_var('aofltxwo','zonal wind stress on water','N m-2',215,aofltxwo,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('aofltywe','meridional wind stress on water','N m-2',216,aofltywe,'v','s'))

    CALL varlist_add(ocean_varlist,new_var('aofltxio','zonal wind stress on ice','N m-2',217,aofltxio,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('aofltyie','meridional wind stress on ice','N m-2',218,aofltyie,'v','s'))
    CALL varlist_add(ocean_varlist,new_var('aoflfrwo','liquid freshwater flux','m s-1',219,aoflfrwo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('aoflfrio','solid freshwater flux over ice','m s-1',220,aoflfrio,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('aoflrhio','residual heat flux used to melt snow/ice','W m-2',221,aoflrhio,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('aoflnhwo','net heat flux over water','W m-2',222,aoflnhwo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('aoflshwo','downwelling solar radiation','W m-2',223,aoflshwo,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('aoflchio','conductive heat flux through ice','W m-2',224,aoflchio,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('aoflwsvo','10m_wind_velocity','m s-1',225,aoflwsvo,'p','s'))


    CALL varlist_add(ocean_varlist,new_var('strairx','zonal_wind_stress_on_ice','N m-2',226,strairx,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('strairy','meridional_wind_stress_on_ice','N m-2',227,strairy,'v','s'))
    CALL varlist_add(ocean_varlist,new_var('strocx','zonal_ocean_stress_on_ice','N m-2',228,strocx,'u','s'))
    CALL varlist_add(ocean_varlist,new_var('strocy','meridional_ocean_stress_on_ice','N m-2',229,strocy,'v','s'))

    CALL varlist_add(ocean_varlist,new_var('ut0','product_of_sea_water_x_velocity_and_sea_water_potential_temperature' &
         ,'C m s-1',230,ut0,'u','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('us0','product_of_sea_water_x_velocity_and_sea_water_salinity'              &
         ,'psu m s-1',231,us0,'u','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('uu0','square_of_sea_water_x_velocity','m2 s-2',232,uu0,'u','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('uv0','product_of_sea_water_x_velocity_and_sea_water_y_velocity'            &
         ,'m2 s-2',233,uv0,'u','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('uw0','product_of_sea_water_x_velocity_and_upward_sea_water_velocity'       &
         ,'m2 s-2',234,uw0,'u','c',staggered=.true.))

    CALL varlist_add(ocean_varlist,new_var('vt0','product_of_sea_water_y_velocity_and_sea_water_potential_temperature' &
         ,'C m s-1',235,vt0,'v','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('vs0','product_of_sea_water_y_velocity_and_sea_water_salinity'              &
         ,'psu m s-1',236,vs0,'v','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('vv0','square_of_sea_water_y_velocity','m2 s-2',237,vv0,'v','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('vu0','product_of_sea_water_y_velocity_and_sea_water_x_velocity'            &
         ,'m2 s-2',238,vu0,'v','c',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('vw0','product_of_sea_water_y_velocity_and_upward_sea_water_velocity'       &
         ,'m2 s-2',239,vw0,'v','c',staggered=.true.))

    CALL varlist_add(ocean_varlist,new_var('wt0','product_of_upward_sea_water_velocity_and_sea_water_potential_temperature' &
         ,'C m s-1',240,wt0,'s','i',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('ws0','product_of_upward_sea_water_velocity_and_sea_water_salinity'         &
         ,'psu m s-1',241,ws0,'s','i',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('ww0','square_of_upward_sea_water_velocity','m2 s-2',242,ww0,'s','i',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('wu0','product_of_upward_sea_water_velocity_and_sea_water_x_velocity'       &
         ,'m2 s-2',243,wu0,'s','i',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('wv0','product_of_upward_sea_water_velocity_and_sea_water_y_velocity'       &
         ,'m2 s-2',244,wv0,'s','i',staggered=.true.))

    CALL varlist_add(ocean_varlist,new_var('u100','sea_water_x_velocity_at_100m','m s-1',245,uko,'u','s=100',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('v100','sea_water_y_velocity_at_100m','m s-1',246,vke,'v','s=100',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('t100','sea_water_potential_temperature_at_100m','C',247,tho,'p','s=100'))
    CALL varlist_add(ocean_varlist,new_var('s100','sea_water_salinity_at_100m','psu',248,sao,'p','s=100'))

    CALL varlist_add(ocean_varlist,new_var('u2000','sea_water_x_velocity_at_2000m','m s-1',249,uko,'u','s=2000',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('v2000','sea_water_y_velocity_at_2000m','m s-1',250,vke,'v','s=2000',staggered=.true.))
    CALL varlist_add(ocean_varlist,new_var('t2000','sea_water_potential_temperature_at_2000m','C',251,tho,'p','s=2000'))
    CALL varlist_add(ocean_varlist,new_var('s2000','sea_water_salinity_at_2000m','psu',252,sao,'p','s=2000'))

    CALL varlist_add(ocean_varlist,new_var('hibete','hibete','fixme',501,hibete,'s','s'))
    CALL varlist_add(ocean_varlist,new_var('hibeto','hibeto','fixme',502,hibeto,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('hibzete','hibzete','fixme',503,hibzete,'s','s'))
    CALL varlist_add(ocean_varlist,new_var('hibzeto','hibzeto','fixme',504,hibzeto,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('tice','ice_surface_temperature','C',99,tice,'p','s'))
    !  CALL varlist_add(ocean_varlist,new_var('rdt','rdt','1',999,rdt,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('rbek','map_of_the_diagnostic_regions','1',511,rbek,'p','s'))
    CALL varlist_add(ocean_varlist,new_var('secmap','map_of_the_diagnostic_sections','1',512,secmap,'p','s'))


    CALL varlist_add(ocean_varlist,new_var('total_volume_of_liquid_sea_water'                                   &
         ,'total_volume_of_liquid_sea_water','m3',513,global_volume,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('total_mass_of_sea_water'                                            &
         ,'total_mass_of_sea_water','kg',514,global_mass,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('global_average_sea_water_potential_temperature'                     &
         ,'global_average_sea_water_potential_temperature','K',515,global_mean_temperature,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('global_average_sea_water_salinity'                                  &
         ,'global_average_sea_water_salinity','psu',516,global_mean_salinity,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('global_salt_content'                                  &
         ,'global_salt_content','kg',517,global_salt_content,'g','g'))

    ! to be implemented
    ! global_average_sea_level_change

    CALL varlist_add(ocean_varlist,new_var('gmsl_st'                                  &
         ,'global_average_steric_sea_level_change','m',518,gmsl_st,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('gmsl_eu'                                  &
         ,'global_average_eustatic_sea_level_change','m',519,gmsl_eu,'g','g'))

    ! global_average_thermosteric_sea_level_change (maybe not - requires recalculation of density !!)


    TsSecOffset=610
    TsCodeperSec=10

    DO n=1,SIZE(section)
      ncode=TsSecOffset+(n-1)*TsCodeperSec
      CALL varlist_add_section(ncode,n)
    ENDDO

    TsRegOffset=820
    TsCodeperReg=16

    DO n=greenland_iceland_norwegian_sea,global_ocean
      ncode=TsRegOffset+(n-1)*TsCodeperReg
      CALL varlist_add_region(ncode,n)
    ENDDO


    !  WRITE(0,*)'call print_varlist'
    IF ( p_pe == p_io ) CALL print_varlist('mpiom.partab',ocean_varlist)

    CALL generate_gridid(ocean_varlist)

    CALL generate_zaxisid(ocean_varlist)

  END SUBROUTINE build_ocean_varlist


  !>
  !! Distribute restart data to mirroring variables or temporaries.
  !!
  SUBROUTINE varlist_post_process_ocean_data
    uoo = uko
    voe = vke
  END SUBROUTINE varlist_post_process_ocean_data


  SUBROUTINE varlist_add_section(ncode,n)

    USE mo_diagnosis, ONLY : section

    INTEGER :: ncode,n

    CALL varlist_add(ocean_varlist,new_var('NetHeatTransport_'//section(n)%SecName       &
         ,'ocean_heat_transport_across_line_'//section(n)%SecName                        &
         ,'W m3 kg-1',ncode,section(n)%NetHeatTransport,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('NetSaltTransport_'//section(n)%SecName       &
         ,'ocean_salt_transport_across_line_'//section(n)%SecName                        &
         ,'psu m3 s-1',ncode+1,section(n)%NetSaltTransport,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('NetWaterTransport_'//section(n)%SecName     &
         ,'ocean_water_transport_across_line_'//section(n)%SecName                       &
         ,'kg s-1',ncode+2,section(n)%NetWaterTransport,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('SiceTransport_'//section(n)%SecName          &
         ,'ocean_seaice_transport_across_line_'//section(n)%SecName                      &
         ,'kg s-1',ncode+3,section(n)%SiceTransport,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('layer2transport_'//section(n)%SecName        &
           ,'layer2_transport_across_line_'//section(n)%SecName                          &
           ,'kg s-1',ncode+4,section(n)%layer2transport,'g','g'))

    CALL varlist_add(ocean_varlist,new_var('layer1transport_'//section(n)%SecName        &
           ,'layer1_transport_across_line_'//section(n)%SecName                          &
           ,'kg s-1',ncode+5,section(n)%layer1transport,'g','g'))


  END SUBROUTINE varlist_add_section

  SUBROUTINE varlist_add_region(ncode,n)

    USE mo_diagnosis, ONLY : region

    INTEGER :: ncode,n

      CALL varlist_add(ocean_varlist,new_var('eisab_'//region(n)%RegionName      &
           ,'sea_ice_extent_in_region_'//region(n)%RegionName                    &
           ,'m2',ncode,region(n)%eisab,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('eiscb_'//region(n)%RegionName      &
           ,'sea_ice_volume_in_region_'//region(n)%RegionName                    &
           ,'m3',ncode+1,region(n)%eiscb,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('hflb_'//region(n)%RegionName       &
           ,'heat_uptake_in_region_'//region(n)%RegionName                       &
           ,'W',ncode+2,region(n)%hflb,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('wflb_'//region(n)%RegionName       &
           ,'freshwater_uptake_in_region_'//region(n)%RegionName                 &
           ,'m3 s-1',ncode+3,region(n)%wflb,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('tem_0m_'//region(n)%RegionName     &
           ,'sst_in_region_'//region(n)%RegionName                               &
           ,'C',ncode+4,region(n)%tem_0m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('salt_0m_'//region(n)%RegionName    &
           ,'sss_in_region_'//region(n)%RegionName                               &
           ,'psu',ncode+5,region(n)%salt_0m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('u2pv2_0m_'//region(n)%RegionName   &
           ,'kinetic_energy_in_region_'//region(n)%RegionName                    &
           ,'m2 s-2',ncode+6,region(n)%u2pv2_0m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('tem_200m_'//region(n)%RegionName   &
           ,'temperature_200m_in_region_'//region(n)%RegionName                  &
           ,'C',ncode+7,region(n)%tem_200m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('salt_200m_'//region(n)%RegionName  &
           ,'salinity_200m_in_region_'//region(n)%RegionName                     &
           ,'psu',ncode+8,region(n)%salt_200m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('u2pv2_200m_'//region(n)%RegionName &
           ,'kinetic_energy_in_region_'//region(n)%RegionName                    &
           ,'m2 s-2',ncode+9,region(n)%u2pv2_200m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('tem_800m_'//region(n)%RegionName   &
           ,'temperature_800m_in_region_'//region(n)%RegionName                  &
           ,'C',ncode+10,region(n)%tem_800m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('salt_800m_'//region(n)%RegionName  &
           ,'salinity_800m_in_region_'//region(n)%RegionName                     &
           ,'psu',ncode+11,region(n)%salt_800m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('u2pv2_800m_'//region(n)%RegionName &
           ,'kinetic_energy_in_region_'//region(n)%RegionName                    &
           ,'m2 s-2',ncode+12,region(n)%u2pv2_800m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('tem_2000m_'//region(n)%RegionName  &
           ,'temperature_2000m_in_region_'//region(n)%RegionName                 &
           ,'C',ncode+13,region(n)%tem_2000m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('salt_2000m_'//region(n)%RegionName &
           ,'salinity_2000m_in_region_'//region(n)%RegionName                    &
           ,'psu',ncode+14,region(n)%salt_2000m,'g','g'))

      CALL varlist_add(ocean_varlist,new_var('u2pv2_2000m_'//region(n)%RegionName &
           ,'kinetic_energy_in_region_'//region(n)%RegionName                     &
           ,'m2 s-2',ncode+15,region(n)%u2pv2_2000m,'g','g'))


  END SUBROUTINE varlist_add_region

#endif/*ndef NO_NEW_IO */

END MODULE mo_varlist
