MODULE MO_CDI

#ifndef NOCDI
  USE mo_kind,   ONLY: wp
  USE mo_mpi
  USE mo_parallel

  USE mo_param1, ONLY: ie,je,ie_g,je_g,ke,kep 
  USE mo_commo1, ONLY: tiestu,tiestw,uko,vke,tho,sao,zo,z1o   &
       ,hibete,hibeto,hibzete,hibzeto,depto,dvo,avo,wo        &
       ,sictho,sicomo,sicsno,sicuo,sicve

  IMPLICIT NONE

  INCLUDE 'cdi.inc'

  TYPE variable_list 
     CHARACTER(len=64)  :: name
     CHARACTER(len=32)  :: unit
     INTEGER            :: code
     INTEGER            :: dim
     REAL(wp), POINTER  :: array(:,:,:)
     INTEGER            :: varID
  END TYPE variable_list

  TYPE (variable_list), PUBLIC   :: variables(19)
  !  TYPE (variable_list), POINTER :: variables


  INTERFACE new_var
     MODULE PROCEDURE new_var_2d
     MODULE PROCEDURE new_var_3d
  END INTERFACE 

CONTAINS


  SUBROUTINE new_var_2d(variable, name, unit, code, dim, field)

    USE mo_param1, ONLY: ie,je,ie_g,je_g 
    USE mo_parallel

    TYPE (variable_list), INTENT(inout) :: variable 
    CHARACTER(len=*),     INTENT(in)    :: name
    CHARACTER(len=*),     INTENT(in)    :: unit
    INTEGER,              INTENT(in)    :: code
    INTEGER,              INTENT(in)    :: dim
    REAL(wp),             INTENT(in)    :: field(ie,je)
    REAL(wp)                            :: field_g(ie,je)     
    REAL(wp),TARGET                     :: field2_g(ie_g,je_g,1)     

    call gather_arr(field,field_g,0)


    field2_g(:,:,1) = field_g(:,:) 

    variable%name  =  name
    variable%unit  =  unit
    variable%code  =  code
    variable%dim   =  dim    
    variable%array => field2_g

  END SUBROUTINE new_var_2d


 SUBROUTINE new_var_3d(variable, name ,unit, code, dim, field)

    USE mo_param1, ONLY: ie,je,ie_g,je_g 
    USE mo_parallel
    TYPE (variable_list), INTENT(inout) :: variable 
    CHARACTER(len=*),     INTENT(in)    :: name
    CHARACTER(len=*),     INTENT(in)    :: unit
    INTEGER,              INTENT(in)    :: code
    INTEGER,              INTENT(in)    :: dim
    REAL(wp),             INTENT(in)    :: field(:,:,:)     
    REAL(wp), ALLOCATABLE,TARGET        :: field_g(:,:,:)     


    allocate(field_g(ie_g,je_g,dim))

    call gather_arr(field,field_g,0)

    variable%name  =  name
    variable%unit  =  unit
    variable%code  =  code
    variable%dim  =  dim
    variable%array => field_g

  END SUBROUTINE new_var_3d



  SUBROUTINE create_restart_varlist

  IF (p_pe==p_io)THEN
     CALL new_var(variables(1),'sea_water_x_velocity','ms-1',3,ke,uko)
     CALL new_var(variables(2),'sea_water_y_velocity','ms-1',4,ke,vke)
     CALL new_var(variables(3),'sea_water_potential_temperature','C',2,ke,tho)
     CALL new_var(variables(4),'sea_water_salinity','psu',5,ke,sao)
     CALL new_var(variables(5),'sea_surface_height_above_sea_level','m',1,1,zo)
     CALL new_var(variables(6),'sea_surface_height_above_sea_level_change','m',82,1,z1o)
     CALL new_var(variables(7),'sea_ice_thickness','m',13,1,sictho)
     CALL new_var(variables(8),'sea_ice_concentration','1',15,1,sicomo)
     CALL new_var(variables(9),'sea_ice_x_velocity','m',35,1,sicuo)
     CALL new_var(variables(10),'sea_ice_y_velocity','1',36,1,sicve)
     CALL new_var(variables(11),'snow_thickness','1',141,1,sicsno)
     CALL new_var(variables(12),'hibete','??',501,1,hibete)
     CALL new_var(variables(13),'hibeto','??',502,1,hibeto)
     CALL new_var(variables(14),'hibzete','??',503,1,hibzete)
     CALL new_var(variables(15),'hibzeto','??',504,1,hibzeto)
     CALL new_var(variables(16),'dvo','??',111,kep,dvo)
     CALL new_var(variables(17),'avo','??',110,kep,avo)
     CALL new_var(variables(18),'wo','ms-1',7,kep,wo)
     CALL new_var(variables(19),'depto','m',84,1,depto)     
  ENDIF



END SUBROUTINE create_restart_varlist

#endif /%NOCDI%/

END MODULE MO_CDI
