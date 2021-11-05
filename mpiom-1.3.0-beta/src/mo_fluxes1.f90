      MODULE MO_FLUXES1
      USE MO_PARAM1
      IMPLICIT NONE
! ---------------------------------------------------------------------
!
!*    *COMMON* *FLUXES1* - Atmospheric fluxes on the ocean grid.
!
!     S. Legutke          *DKRZ*           30.05.97
!sv   S. Venzke           *MPI*            21.07.99
!sv   restricted to odd fields for c-Grid and added wind stress velocity
!
!*    VARIABLE        TYPE       PURPOSE.
!     --------        ----       --------
!     *AOFLTXW*       *REAL*    Zonal wind stress on water
!     *AOFLTXI*       *REAL*    Zonal wind stress on snow/ice
!     *AOFLTYW*       *REAL*    Meridional wind stress on water
!     *AOFLTYI*       *REAL*    Meridional wind stress on snow/ice
!     *AOFLFRW*       *REAL*    liquid freshwater flux (over water and ice)
!     *AOFLFRI*       *REAL*    solid freshwater flux (over ice only)
!     *AOFLRHI*       *REAL*    residual heat flux used to melt snow/ice
!     *AOFLNHW*       *REAL*    net heat flux over water
!     *AOFLSHW*       *REAL*    downwelling solar radiation
!     *AOFLCHI*       *REAL*    conductive heat flux through ice
!sv   *AOFLWSV        *REAL*    wind stress velocity
!
!*    PURPOSE
!     -------
!     The ocean can be driven either by fluxes (gpp option FLUXES, used
!     with the coupled model) or by surface variables and fluxes. 
!     In the first case the forcing fields reside in COMMON block FLUXES, 
!     in the second case, they are in COMMON block OCEVAL.
!    Note: if (Prec.-evap.(over ice))*comp. is downward and the atmos.
!     temperature at the blending height is above 0 (i.e. it is raining), 
!     then this is added to AOFLFRW before being passed to the ocean 
!     This saves transfer of 1 field. AOFLFRW also includes runoff.
!    Note: downwelling solar radiation is multiplied by surface albedo in
!     the atmosphere model already. This is not the case with the other version.
!    Note: all fluxes besides freshwater are per m^2, and not per grid cell, i.e.
!     they are not multiplied with ice compactness.
!    Note: the relaxation heat flux is declared in COMMON FLUXES2.h
!
! ----------------------------------------------------------------------
!
#ifdef __coupled
      REAL, POINTER :: &
     &AOFLTXWO(:,:),AOFLTYWE(:,:),                                  &
     &AOFLTXIO(:,:),AOFLTYIE(:,:),                                  &
     &AOFLFRIO(:,:),AOFLFRWO(:,:),                                  &
     &AOFLRHIO(:,:),AOFLCHIO(:,:),                                  &
     &AOFLNHWO(:,:),AOFLSHWO(:,:),                                  &
     &AOFLWSVO(:,:),AOFLDHWO(:,:)
#endif
      CONTAINS

      SUBROUTINE alloc_mem_fluxes1
!
! Allocate memory for arrays used in coupling, which are assigned in mo_couple
! nsk, 03.09.2004

#ifdef __coupled
        ALLOCATE(AOFLTXWO(IE,JE),AOFLTYWE(IE,JE),                                  &
     &AOFLTXIO(IE,JE),AOFLTYIE(IE,JE),                                  &
     &AOFLFRIO(IE,JE),AOFLFRWO(IE,JE),                                  &
     &AOFLRHIO(IE,JE),AOFLCHIO(IE,JE),                                  &
     &AOFLNHWO(IE,JE),AOFLSHWO(IE,JE),                                  &
     &AOFLWSVO(IE,JE),AOFLDHWO(IE,JE))

#endif
      END SUBROUTINE alloc_mem_fluxes1

      END MODULE MO_FLUXES1
