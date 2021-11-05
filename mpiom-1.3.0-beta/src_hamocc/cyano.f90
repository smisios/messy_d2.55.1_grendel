       SUBROUTINE CYANO(kpie,kpje,kpke,pddpo)
!**********************************************************************
!
!**** *CYANO* -  .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,             *MPI-MaD, HH*    10.04.01
!     - included : surface reduction of gaseous nitrogen
!
!     Purpose
!     -------
!     Nitrate reduction by cyano bacteria (2NO3 + O2 => N2O + O2).
!
!     Method:
!     ------
!
!     *CALL*       *CYANO(kpie,kpje,kpke,pddpo)*
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
!     .
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      USE mo_control_bgc
      use mo_param1_bgc 


implicit none

      INTEGER :: kpie,kpje,kpke, i,j,k
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: oldocetra,contppm
      
      contppm=1./0.35e-3    ! only used if DIFFAT
!
!  N-fixation by cyano bacteria from atmospheric nitrogen 
!  if nitrate is below redfield ratio wrt phosphate (this is not a surface flux!)
!  N2 = nitrogen, NO3 = nitrate
!  from tech report: 
!  uptake of atmospheric nitrogen and its immediate release as nitrate by diazotrophs
!
      DO j=1,kpje
      DO i=1,kpie
        IF(pddpo(i,j,1).GT.0.5) THEN
          IF(ocetra(i,j,1,iano3).LT.(rnit*ocetra(i,j,1,iphosph))) THEN

            oldocetra = ocetra(i,j,1,iano3)
            ocetra(i,j,1,iano3)= ocetra(i,j,1,iano3)*(1-bluefix)       &
     &                          +bluefix*rnit*ocetra(i,j,1,iphosph)

            ocetra(i,j,1,igasnit)=ocetra(i,j,1,igasnit)-                          &
     &      (ocetra(i,j,1,iano3)-oldocetra)*.5          ! *(1./2.) half of the nitrogen from gaseous pool (?)

            ocetra(i,j,1,ioxygen)=ocetra(i,j,1,ioxygen)-                          &
     &      (ocetra(i,j,1,iano3)-oldocetra)*1.5              ! 1.5 (?) should be 172./122. redfield ratios O2/C , PO4/C

!#ifdef DIFFAT
!            atm(i,j,iatmn2)=atm(i,j,iatmn2) -                          &
!     &      (ocetra(i,j,1,iano3)-oldocetra)/2.*contppm*pddpo(i,j,1)
!#endif
            ENDIF  
         ENDIF  
      ENDDO
      ENDDO

      RETURN
      END
