#include "messy_main_ppd_bi.inc"

!*****************************************************************************
MODULE messy_main_grid_def_bi
!*****************************************************************************
!
!*****************************************************************************
!
!*****************************************************************************
!
#if defined(ECHAM5)
#include "messy_main_grid_def_echam5.inc"
#endif

!==============================================================================

#ifdef COSMO
#include "messy_main_grid_def_cosmo.inc"
#endif

!==============================================================================

#ifdef CESM1
#include "messy_main_grid_def_cesm1.inc"
#endif

!==============================================================================

#ifdef ICON
#include "messy_main_grid_def_icon.inc"
#endif

!==============================================================================


#if defined(MBM_MPIOM)
#include "messy_main_grid_def_mpiom.inc"
#endif

!==============================================================================

#if defined(BLANK) || defined(MESSYDWARF)
#if defined(MBM_RAD)
#include "messy_main_grid_def_rad.inc"
#else
#if defined(MBM_DISSOC)
#include "messy_main_grid_def_dissoc.inc"
#else
#if defined(MBM_QBO)
#include "messy_main_grid_def_qbo.inc"
#else
#include "messy_main_grid_def_blank.inc"
#endif
#endif
#endif
#endif

!==============================================================================

#ifdef MBM_CLAMS
#include "messy_main_grid_def_clams.inc"
#endif

!==============================================================================

#if defined(VERTICO)
#include "messy_main_grid_def_vertico.inc"
#endif

!==============================================================================

!*****************************************************************************
END MODULE messy_main_grid_def_bi
!*****************************************************************************
