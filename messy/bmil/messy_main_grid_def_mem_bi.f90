MODULE messy_main_grid_def_mem_bi

  USE messy_main_constants_mem, ONLY: dp

#ifdef ECHAM5
#include "messy_main_grid_def_mem_echam5.inc"
#endif
#ifdef CESM1
#include "messy_main_grid_def_mem_cesm1.inc"
#endif
#ifdef COSMO
#include "messy_main_grid_def_mem_cosmo.inc"
#endif
#ifdef ICON
#include "messy_main_grid_def_mem_icon.inc"
#endif

#if defined(BLANK) || defined(MESSYDWARF)
#if defined(MBM_RAD)
#include "messy_main_grid_def_mem_rad.inc"
#else
#if defined(MBM_DISSOC)
#include "messy_main_grid_def_mem_dissoc.inc"
#else
#if defined(MBM_QBO)
#include "messy_main_grid_def_mem_qbo.inc"
#else
#include "messy_main_grid_def_mem_blank.inc"
#endif
#endif
#endif
#endif
#if defined(MBM_CLAMS)
#include "messy_main_grid_def_mem_clams.inc"
#endif
#if defined(MBM_MPIOM)
#include "messy_main_grid_def_mem_mpiom.inc"
#endif
#ifdef VERTICO
#include "messy_main_grid_def_mem_vertico.inc"
#endif
 ! additional info see grid_bi:
 INTEGER, ALLOCATABLE, DIMENSION(:) :: BASEGRID_ID

END MODULE messy_main_grid_def_mem_bi
