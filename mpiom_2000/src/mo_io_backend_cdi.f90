!>
!! Short description goes here, don't forget period at end.
!!
!! @todo add module description for mo_file_list_cdi.f90
!!
MODULE mo_io_backend_cdi

  USE mo_kind
  USE mo_util_string

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vlist_new
  PUBLIC :: vlist_get_taxis
  PUBLIC :: vlist_var_new

  INCLUDE 'cdi.inc'

CONTAINS

  !>
  !! Create a new variable list and initialize time axis information.
  !!
  INTEGER FUNCTION vlist_new()
    vlist_new = vlistCreate()
!    CALL vlistDefTaxis(vlist_new, taxisCreate(TAXIS_ABSOLUTE))
    CALL vlistDefTaxis(vlist_new, taxisCreate(TAXIS_RELATIVE))
  END FUNCTION vlist_new

  !>
  !! Get time axis ID from given variable list.
  !!
  INTEGER FUNCTION vlist_get_taxis(this)
    INTEGER, INTENT(in) :: this
    vlist_get_taxis = vlistInqTaxis(this)
  END FUNCTION vlist_get_taxis

  !>
  !! Create a new variable for a given variable list.
  !!
  !! This requires a couple of settings, concerning alphanumerical ID (name),
  !! naming according to a set standard, human readable description (long name),
  !! unit of measure, numerical ID according to set standard, representation of
  !! missing values, and numerical datatype representation.
  !!
  INTEGER FUNCTION vlist_var_new(vlist, grid, zaxis, name, std_name, unit, &
       code, missval, high_precision)

    INTEGER, INTENT(in) :: vlist, grid, zaxis
    CHARACTER(*), INTENT(in) :: name
    CHARACTER(*), INTENT(in) :: std_name
    CHARACTER(*), INTENT(in) :: unit
    INTEGER, INTENT(in) :: code
    REAL(dp), INTENT(in) :: missval
    LOGICAL, INTENT(in) :: high_precision

    vlist_var_new = vlistDefVar(vlist, grid, zaxis, TIME_VARIABLE)
    CALL vlistDefVarName(vlist, vlist_var_new, name)
    CALL vlistDefVarStdname(vlist, vlist_var_new, &
         translate(tolower(std_name), ' ', '_'))
    CALL vlistDefVarLongname(vlist, vlist_var_new, &
         capitalize(translate(TRIM(std_name), '_', ' ')))
    CALL vlistDefVarUnits(vlist, vlist_var_new, unit)
    CALL vlistDefVarCode(vlist, vlist_var_new, code)
    CALL vlistDefVarMissval(vlist, vlist_var_new, missval)
    CALL vlistDefVarDatatype(vlist, vlist_var_new, &
         MERGE(DATATYPE_FLT64, DATATYPE_FLT32, high_precision))
  END FUNCTION vlist_var_new

END MODULE mo_io_backend_cdi
