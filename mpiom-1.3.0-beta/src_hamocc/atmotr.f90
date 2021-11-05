      SUBROUTINE ATMOTR(kpie,kpje,kpke,kplmon,pdlxp,pdlyp)
!*******************************************************************
!
!**** *ATMOTR* - calculate the atmospheric diffusion
!
!     Maier-Reimer,          *MPI-Met, HH*
!
!     P.Wetzel,              *MPI-Met, HH*    27.03.02
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!     - calculate the atmospheric diffusion for surface CO2 exchange.
!
!     Method
!     -------
!     -
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *CALL*       *ATMOTR(kpie,kpje,kpke)*
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_control_bgc
      USE mo_bgcmean
      use mo_param1_bgc 

      use mo_parallel
      implicit none
      
      INTEGER :: kpie,kpje,kpke,kplmon
      INTEGER :: i,j,k,l,is
      REAL, DIMENSION(kpie,kpje) :: ucor, vcor
      REAL, DIMENSION(kpie,kpje) :: pdlxp,pdlyp
      REAL :: ppm2con
 
! distribution of emissions for DIFFAT, don't want that for coupled run
#ifndef __cpl_co2
#ifdef PANTHROPOCO2 

! ppm2con: atmospheric weight: ~10000kg/m^2, avrg. ~29 g/mol
! --> 350 kmol/m^2 --> 1ppm ~ 0.35e-3 kmol/m^2
      ppm2con=0.35e-3
      
! North America:      
        atm(25,35,iatmco2) = atm(25,35,iatmco2)                         &
     &      +0.27*ems_per_step/(pdlxp(25,35)*pdlyp(25,35)*ppm2con)
! Asia:      
        atm(95,45,iatmco2) = atm(95,45,iatmco2)                         &
     &      +0.42*ems_per_step/(pdlxp(95,45)*pdlyp(95,45)*ppm2con)
! Europe:      
        atm(75,35,iatmco2) = atm(75,35,iatmco2)                         &
     &      +0.17*ems_per_step/(pdlxp(75,35)*pdlyp(75,35)*ppm2con)      
! South America:      
        atm(45,65,iatmco2) = atm(45,65,iatmco2)                         &
     &      +0.07*ems_per_step/(pdlxp(45,65)*pdlyp(45,65)*ppm2con)
! Africa:      
        atm(75,70,iatmco2) = atm(75,70,iatmco2)                         &
     &      +0.05*ems_per_step/(pdlxp(75,70)*pdlyp(75,70)*ppm2con)
! Australia:      
        atm(113,71,iatmco2) = atm(113,71,iatmco2)                       &
     &      +0.02*ems_per_step/(pdlxp(113,71)*pdlyp(113,71)*ppm2con)
#endif
#endif 
      

      do j=1,kpje
      do i=1,kpie
        ucor(i,j)=0.
        vcor(i,j)=0.
      enddo
      enddo
      
!      do l=1,natm ! outer loop over all tracers
          
!     call perio2(123,atm(1,1,l))
      call bounds_exch('p',atm)
      
      do l=1,3     ! 3 hard coded !!!!!!!!!!!!!!!!!!!!!!!!!  
      do is=1,4 ! speed up
      
      do j=2,kpje-1
      do i=2,kpie-1
!        ucor(i,j)=0.01*(atm(i+1,j)-atm(i,j))*dlxu(i,j)*dlyu(i,j)
!        vcor(i,j)=0.01*(atm(i,j)-atm(i,j+1))*dlxv(i,j)*dlyv(i,j)*atdifv(i,j)
        ucor(i,j)=0.2*(atm(i+1,j,l)-atm(i,j,l))*pdlxp(i,j)*pdlyp(i,j)
        vcor(i,j)=0.2*(atm(i,j,l)-atm(i,j+1,l))*pdlxp(i,j)*pdlyp(i,j)*atdifv(i,j)
      enddo
      enddo
!     call perio2(124,ucor)
!     call perio2(125,vcor)
      call bounds_exch('u+',ucor)
      call bounds_exch('v+',vcor)
      
      do j=2,kpje
      do i=2,kpie
        atm(i,j,l)=atm(i,j,l)+(ucor(i,j)-ucor(i-1,j)                &
     &           +vcor(i,j-1)-vcor(i,j))/(pdlxp(i,j)*pdlyp(i,j))
      enddo
      enddo
!     call perio2(126,atm(1,1,l))

      enddo  ! speed up
      enddo  ! outer loop over all tracers
      call bounds_exch('p',atm)

#ifdef DIFFAT
      do j=1,kpje
      do i=1,kpie
        bgcm2d(i,j,jatmco2)=bgcm2d(i,j,jatmco2)+atm(i,j,iatmco2)
      bgcm2d(i,j,jatmo2) =bgcm2d(i,j,jatmo2) +atm(i,j,iatmo2)
      bgcm2d(i,j,jatmn2) =bgcm2d(i,j,jatmn2) +atm(i,j,iatmn2)
      enddo
      enddo
#endif
      
      RETURN
      END
