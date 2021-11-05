module mo_iso_c_kinds
#ifdef HAVE_CONFIG_H
#include "config.h"
#ifdef HAVE_FORTRAN_MODULE_ISO_C_BINDING
  use iso_c_binding
#endif
#ifdef C_INT64_T
  integer, parameter :: c_int64_t = C_INT64_T
#endif
#ifdef C_FLOAT
  integer, parameter :: c_float = C_FLOAT
#endif
#else
  USE mo_kind, ONLY : i8, sp
  integer, parameter :: c_int64_t = i8
  integer, parameter :: c_float = sp
#endif
  public c_float, c_int64_t
end module mo_iso_c_kinds
