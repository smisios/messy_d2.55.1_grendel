! -*- f90 -*-
#ifdef NEVER
!
! NOTE:
! The pre-processor macros, named _RI_ (for <R>ank <I>dentifier), are
! used to flip the meaning of ranks (e.g. the different spatial dimensions)
! for the different basemodels. X, Y and Z denote the spatial dimensions,
! and N a free parameter dimension, e.g. for the number of tracers,
! the number of (aerosol) modes, etc.
!
#endif

! ----------------------------------------------------------------------
! SETTINGS FOR VARIOUS BASEMODELS
! ----------------------------------------------------------------------
! SO FAR FOR ALL MESSy BASEMODES
#define _RI_XY___(X,Y) X,Y
#define _IX_XY___ 1
#define _IY_XY___ 2

! ORDER: h,v,[n],h (TRRANK=3)  OR h,v,h,[n] (TRRANK=4)
#if defined(ECHAM5) || defined (MBM_TRACER) || defined (CESM1) || defined (ICON)
#define _XRI_hvh_
#endif

! ORDER: h,h,[n],v (TRRANK=3)  OR h,h,v,[n] (TRRANK=4)
#if defined (COSMO) || defined (BLANK) || defined(VERTICO) || defined (MBM_CLAMS) || defined(MBM_MPIOM) || defined(MESSYDWARF)
#define _XRI_hhv_
#endif

#if defined (ICON) || defined(COSMO) || defined(MESSYDWARF) || defined (BLANK)
#define _RI_TRRANK 4
#else
#define _RI_TRRANK 3
#endif

! ----------------------------------------------------------------------
! GENERIC SETTINGS
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ORDER: h,v,[n],h OR h,v,h,[n]
#if defined(_XRI_hvh_)
! ----------------------------------------------------------------------

! op_pj_20170228: macros required for TTRANK=3 (messy_main_tracer)
#if _RI_TRRANK == 3

#define _RI_XYZN_(X,Y,Z,N) X,Z,N,Y
#define _IX_XYZN_ 1
#define _IY_XYZN_ 4
#define _IZ_XYZN_ 2
#define _IN_XYZN_ 3

#define _RI_XYZ__(X,Y,Z) X,Z,Y
#define _IX_XYZ__ 1
#define _IY_XYZ__ 3
#define _IZ_XYZ__ 2

  ! classical RI to pass one `jrow` to SMCL:
#define _RI_X_ZN_(X,Z,N) X,Z,N
#define _IX_X_ZN_ 1
#define _IZ_X_ZN_ 2
#define _IN_X_ZN_ 3

  ! only ONEMIS, BIOBURN, GMXE, VEG:
#define _RI_XY_N_(X,Y,N) X,N,Y
#define _IX_XY_N_ 1
#define _IY_XY_N_ 3
#define _IN_XY_N_ 2
#endif

! op_pj_20170228: macros required for TTRANK=4 (messy_main_tracer)
#if _RI_TRRANK == 4

#define _RI_XYZN_(X,Y,Z,N) X,Z,Y,N
#define _IX_XYZN_ 1
#define _IY_XYZN_ 3
#define _IZ_XYZN_ 2
#define _IN_XYZN_ 4

#define _RI_XYZ__(X,Y,Z) X,Z,Y
#define _IX_XYZ__ 1
#define _IY_XYZ__ 3
#define _IZ_XYZ__ 2

#define _RI_X_ZN_(X,Z,N) X,Z,N
#define _IX_X_ZN_ 1
#define _IZ_X_ZN_ 2
#define _IN_X_ZN_ 3

#define _RI_XY_N_(X,Y,N) X,Y,N
#define _IX_XY_N_ 1
#define _IY_XY_N_ 2
#define _IN_XY_N_ 3

#endif
!
! ----------------------------------------------------------------------
#endif
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! ORDER: h,h,[n],v OR h,h,v,[n]
#if defined(_XRI_hhv_)
! ----------------------------------------------------------------------

! op_pj_20170228: macros required for TTRANK=3 (messy_main_tracer)
#if _RI_TRRANK == 3

#define _RI_XYZN_(X,Y,Z,N) X,Y,N,Z
#define _IX_XYZN_ 1
#define _IY_XYZN_ 2
#define _IZ_XYZN_ 4
#define _IN_XYZN_ 3

#define _RI_XYZ__(X,Y,Z) X,Y,Z
#define _IX_XYZ__ 1
#define _IY_XYZ__ 2
#define _IZ_XYZ__ 3

#define _RI_X_ZN_(X,Z,N) X,N,Z
#define _IX_X_ZN_ 1
#define _IZ_X_ZN_ 3
#define _IN_X_ZN_ 2

#define _RI_XY_N_(X,Y,N) X,Y,N
#define _IX_XY_N_ 1
#define _IY_XY_N_ 2
#define _IN_XY_N_ 3

#endif

! op_pj_20170228: macros required for TTRANK=4 (messy_main_tracer)
#if _RI_TRRANK == 4
#define _RI_XYZN_(X,Y,Z,N) X,Y,Z,N
#define _IX_XYZN_ 1
#define _IY_XYZN_ 2
#define _IZ_XYZN_ 3
#define _IN_XYZN_ 4

#define _RI_XYZ__(X,Y,Z) X,Y,Z
#define _IX_XYZ__ 1
#define _IY_XYZ__ 2
#define _IZ_XYZ__ 3

#define _RI_X_ZN_(X,Z,N) X,Z,N
#define _IX_X_ZN_ 1
#define _IZ_X_ZN_ 2
#define _IN_X_ZN_ 3

#define _RI_XY_N_(X,Y,N) X,Y,N
#define _IX_XY_N_ 1
#define _IY_XY_N_ 2
#define _IN_XY_N_ 3

#endif
!
! ----------------------------------------------------------------------
#endif
! ----------------------------------------------------------------------
