      SUBROUTINE SPINUP_BGC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp)

!**********************************************************************
!
!**** *SPINUP_BGC* - relaxation of the alkalinity during
!                    bgc spinup.
!     P.Wetzel,        *MPI, HH*    09.04.02
!
!     Modified
!     --------
!
!     Purpose
!     -------
!     Calculate global mean atmospheric pCO2 and adjust the alkalinity 
!     during bgc spinup state.
!
!     Method
!     -------
!
!**   Interface.
!     ----------
!
!     *CALL*       *SPINUP_BGC(kpie,kpje,kpke,pddpo,pdlxp,pdlyp)*
!
!     called by ini_bgc.f90
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *pdlxp*   - size of scalar grid cell (1st dimension) [m].
!     *REAL*    *pdlyp*   - size of scalar grid cell (2nd dimension) [m].
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_control_bgc
      use mo_param1_bgc 
      use mo_parallel
      
      implicit none

      INTEGER :: kpie,kpje,kpke,i,j,k,l
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL ::    ztotarea, area
      REAL ::    zatc, meanatc, alkcor

! Total area of grid

      ztotarea=0.      
      DO j=2,kpje-1
      DO i=2,kpie-1
         ztotarea = ztotarea + pdlxp(i,j)*pdlyp(i,j)
      ENDDO
      ENDDO

      call global_sum(ztotarea)


! Mean atmosperic CO2

      zatc  = 0.
      DO j=2,kpje-1
      DO i=2,kpie-1
         area  = pdlxp(i,j)*pdlyp(i,j)
         zatc  = zatc  + atm(i,j,iatmco2) *area
      ENDDO
      ENDDO 

      call global_sum(zatc)
      meanatc = zatc/ztotarea

      write(0,*)'ATTN: spinup bgc',zatc,ztotarea,meanatc


!  Correction factor for total alkalinity      
 
      alkcor = 1. + 0.01*(meanatc/278. -1) ! 0.003
! from c4mip:                    286.2

      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        IF(pddpo(i,j,k).GT.0.5) THEN
          ocetra(i,j,k,ialkali)=ocetra(i,j,k,ialkali)*alkcor
      ENDIF
      ENDDO
      ENDDO
      ENDDO

      WRITE(io_stdo_bgc,*) ' '
      WRITE(io_stdo_bgc,*) 'Correction factor for total alkalinity'
      WRITE(io_stdo_bgc,*) '-------------------------------------------'
      WRITE(io_stdo_bgc,*) 'global mean atm. CO2[ppm]: ',meanatc
      WRITE(io_stdo_bgc,*) 'Alk.correction factor: ', alkcor
      

      RETURN
      END
