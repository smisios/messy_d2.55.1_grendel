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
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_control_bgc
      USE mo_param1_bgc
      USE mo_parallel

      IMPLICIT NONE

      INTEGER :: kpie,kpje,kpke,i,j,k
      REAL(wp) :: pddpo(kpie,kpje,kpke)
      REAL(wp) :: pdlxp(kpie,kpje)
      REAL(wp) :: pdlyp(kpie,kpje)
      REAL(wp) :: meanatc,alkcor

      meanatc = 0._wp
      CALL global_mean_2d(atm(:,:,iatmco2),meanatc)

      WRITE(0,*) 'ATTN: spinup bgc'

      !  Correction factor for total alkalinity
      IF (ABS(fspinbgc - 1._wp) < EPSILON(1._wp)) THEN
         ! Compute factor from ratio of actual and target atmospheric CO2 concentration
         alkcor = 1._wp + 0.01_wp*(meanatc/278._wp -1._wp) ! 0.003
         ! from c4mip:                    286.2
      ELSE
         ! Use factor provided in namelist
         alkcor = fspinbgc
         write(0,*) 'ATTN: using alkalinity correction from namelist', alkcor
      END IF

      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        IF (pddpo(i,j,k) > 0.5_wp) THEN
          ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali)*alkcor
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
