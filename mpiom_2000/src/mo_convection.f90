!> Calulation of ocean convection.
!! @author: Aiko Voigt, Helmuth Haak, Max Planck Institute for Meteorology
!> @date Last modified on 09.03.2009 by Helmuth Haak
!> @todo Documentation, Generalisation for arbitrary tracer, check nlopps
!> This module holds subroutines for ocean convection parametrization to
!! reduce static instability.
!! Included subroutines are:
!! convection : Parametrization schemes for convective adjustment.
!! NLOPPS     : Plume convection model
!!
!! Corresponding to the name list parameter 'ioconv' (1,2,3,4) different schemes are used:
!! ioconv = 1 : enhancement of vertical eddy diffusivity (default)
!!              This is sofar the only schme that works with HAMOCC tracers.
!! ioconv = 2 : complete mixing of pot. temperature and salinity
!! ioconv = 3 : interchange pot. temperature and salinity
!! ioconv = 4 : plume convection model following a modified version of the OOPS model

 MODULE mo_convection

  USE mo_param1
  USE mo_param3      ! This is not actually used here but only in subroutine adisitj which is called from octher
                     ! and is supposed to be expanded inline. The purpose of this dependency is to force compilation
                     ! of mo_param3.f90 before compilation of mo_octher.f90, because mo_param3.mod is needed for
                     ! the inline expansion of adisitj.

  USE mo_boundsexch, ONLY : bounds_exch

  IMPLICIT NONE

 CONTAINS

  !>
  !! Parametrization schemes for convective adjustment.
  !!
  !! Four exclusive choices by specifying the namelist paramter ioconv are possible
  !! ioconv = 1 : enhancement of vertical eddy diffusivity and change of
  !! (default)    vertical eddy viscosity in case of unstable stratification in subroutine calc_rinum
  !!              immediately called after subroutine convection, no change of pot. temperature tho
  !!              or salinity sao in subroutine convection
  !!              (corresponds to former CPP flag NURDIF)
  !! ioconv = 2 : complete mixing of pot. temperature tho and salinity sao of upper and lower box
  !!              in case of unstable stratification following Bryan, K., 1969:
  !!              A Numerical Method for the Study of the Circulation of
  !!              the World Ocean, Journal of Computational Physics, 4(3):347-376.
  !!              (corresponds to former CPP flag NURMISCH)
  !! ioconv = 3 : interchange pot. temperature and salinity of upper and lower box
  !!              in case of unstable stratification obeying conservation of heat and salt, should
  !!              result in larger penetrative depth compaired to complete mixing of ioconv=2
  !!              (corresponds to former CPP flag UMKLAP)
  !! ioconv = 4 : plume convection model following a modified version of the OOPS model described in
  !!              Paluskiewicz, T. and R. D. Romea, 1997: A one-dimensional model for
  !!              the parameterization of deep convection in the ocean, Dynamics of
  !!              Atmospheres and Oceans, 26, 95-130,
  !!              invoked by call of subroutine nloops in subroutine convection
  !!              a one-to-one documentation for nlopps is still missing
  !!              (corresponds to former CPP flag PLUME)
  !!
  !! rewrite by Aiko Voigt, Max Planck Institute for Meteorology (Jan 13, 2009)

  SUBROUTINE convection

    USE mo_planetary_constants, ONLY: g, inv_rhoref_water  &
           ,rhoicwa,rhosnwa
    USE mo_mean, only : kcondep
    USE mo_commo1, ONLY : dt, ioconv, vk1e, sao, &
         tho, preffw, stabio, po, tiestu, tiestw, s1o, dz, lweto, weto, kcondep, &
         ddpo, zo, sictho, sicsno, dzw, rhoo
    USE mo_mpi, only : p_abort

  !! local scalars

    INTEGER i         !< loop variable in W-E direction
    INTEGER j         !< loop variable in N-S direction
    INTEGER k         !< loop variable for vertical levels (k=1 top, k=ke bottom)

    REAL(wp) :: shelp(ie,je),thelp(ie,je),rhelp(ie,je)

    REAL(wp) stabio1      !< vertical density gradient
    REAL(wp) disti        !< vertical distance of two pressure points, i.e. 1/dz(k)
    REAL(wp) dtts         !< time step paramter for plume convection subroutine nlopps
    REAL(wp) tupper       !< pot. temperature of upper box
    REAL(wp) tlower       !< pot. temperature of lower box
    REAL(wp) supper       !< salinity of upper box
    REAL(wp) slower       !< salinity of lower box
    REAL(wp) ddhelp       !< thickness of upper box including sea-level for top box
    REAL(wp) sq           !< depth-weighted average of salinity of upper and lower box
    REAL(wp) tq           !< depth-weighted average of pot. temperature of upper and lower box
    REAL(wp) switc        !< switch to detect static instability, i.e. 0=stable, 1=unstable

  !! local arrays

    REAL(wp) rhuppo(ie,je) !< in-situ density of the upper box
    REAL(wp) zsurf(ke)     !< mask for top layer, i.e. zsurf=1 for k=1 and 0 otherwise
    REAL(wp) :: dzdzo !< weighting factor for pressure interpolation from u level to w level
    REAL(wp) :: dzdzu !< weighting factor for pressure interpolation from u level to w level

!$omp parallel private(i,j,k, disti,stabio1,tq,sq,switc,    &
!$omp                  tupper, tlower,                    &
!$omp                  supper, slower, ddhelp,dzdzo,dzdzu)



  !! mask for top layer
    zsurf(1) = 1._wp
    DO k=2,ke
       zsurf(k) = 0._wp
    ENDDO


!$omp do
    DO j=1,je   !! big j loop till end of subroutine convection

    !! calculate in-situ density, baroclinic pressure po, and vertical density
    !! gradient stabio

    !! top layer k =1
       DO i=1,ie
          vk1e(i, j, 1) = 0._wp
          shelp(i,j)=sao(i,j,1)
          thelp(i,j)=tho(i,j,1)
       ENDDO
    !! pot. temperature to in-situ temperature thelp
    !! preff is reference pressure set in beleg.f90 according to hydrostatic
       CALL adisitj(thelp,shelp,preffw(1),j)
    !! gives in-situ density rhelp
       CALL rho1j(thelp,shelp,preffw(1),rhelp,j)

       rhoo(:,j,1)=rhelp(:,j)   ! save density for later use

       DO i=1,ie
          stabio(i,j,1) = 0._wp
          po(i,j,1)     = g*tiestu(1)*inv_rhoref_water*rhelp(i,j)  !! pressure for top layer according to hydrostatic equation

          s1o(i,j,1)    = rhelp(i,j)
       ENDDO

       !uwe  include downward propagation of tho and sao in land

       DO k=2,ke
          DO i=1,ie
             IF(.NOT. lweto(i,j,k)) THEN
                tho(i,j,k)=tho(i,j,k-1)
                sao(i,j,k)=sao(i,j,k-1)
             ENDIF
          ENDDO
       ENDDO

    !! layers below top: k=2,..,ke
       DO k=2,ke

          disti = 1._wp / dz(k)
          DO i=1,ie
             shelp(i,j)=sao(i,j,k)
             thelp(i,j)=tho(i,j,k)
          ENDDO
          CALL adisitj(thelp,shelp,preffw(k),j)      !! pot. temperature to in-situ temperature thelp
          CALL rho1j(thelp,shelp,preffw(k),rhelp,j)  !! gives in-situ density rhelp

          rhoo(:,j,k)=rhelp(:,j)   ! save density for later use

          !! calculate in-situ density rhuppo of the above layer k-1
          DO i=1,ie
             shelp(i,j)=sao(i,j,k-1)
             thelp(i,j)=tho(i,j,k-1)
          ENDDO
          CALL adisitj(thelp,shelp,preffw(k),j)
          CALL rho1j(thelp,shelp,preffw(k),rhuppo,j) !! rhuppo is in-situ density of the above layer k-1

          dzdzo=(tiestw(k)-tiestu(k-1))*inv_rhoref_water
          dzdzu=(tiestu(k)-tiestw(k))*inv_rhoref_water
          ! in-situ density at the bottom of the grid box (w level)
          rhoo(:, j, k-1) = 0.5_wp * (rhoo(:, j, k-1) + rhuppo(:, j))

          DO i=1,ie

             s1o(i,j,k)=rhelp(i,j)

             !! calculate vertical density stabio gradient between upper and lower box
             stabio1 = disti * ( rhelp(i,j) - rhuppo(i,j) )  !! vertical density gradient 1/delta_z * (rho(k)-rho(k-1))
                                                             !! stabio1 < 0 instable stratification
                                                             !! stabio1 > 0 stable stratification
             stabio(i,j,k)=MAX(stabio1, 0._wp)                   !! set negative values to zero (why???)
                                                             !! stabio  = 0 instable stratification
                                                             !! stabio1 > 0 stable stratification

             !! po = hydrostatic pressure / ref_water with rhoref_water = 1025 kgm**-3
             !! integrated over depth according to hydrostatic equation,
             !! use the distance weighted average of upper and lower densities
             !! po is defined on the level of the horizontal velocities (u,v)

             po(i,j,k) = po(i,j,k-1) + g*(dzdzu*rhelp(i,j)+dzdzo*rhuppo(i,j))

             !! switch for static instability, only used for ioconv=1,2,3
             !! switc = 0 : static stability, upper box less dense than lower box
             !! switc = 1 : static instability, upper box denser than lower box

             switc = MAX(0._wp, -stabio1/(1.e-11_wp + ABS(stabio1)))*weto(i,j,k)

             !! select convection parametrization
             !! could be done with select case statement if subroutine is rearranged

             IF (ioconv .eq. 1) THEN       !! enhanced vertical eddy diffusivity

                vk1e(i, j, k) = 0._wp
                IF(kcondep(i,j).EQ.k-1) kcondep(i,j) = kcondep(i,j)+NINT(switc)

             ELSE IF (ioconv .eq. 2) THEN      !! complete mixing of upper and lower box

                vk1e(i,j,k) = switc

                if (switc .gt. 0.5_wp .and. lweto(i,j,k)) then

                   ddhelp = ddpo(i,j,k-1) !! thickness of upper box
                   !! thickness of upper box including sea-level for top box
                   IF (k .EQ. 2) THEN
                     ddhelp = ddhelp + (zo(i,j) &
                          - sictho(i,j) * rhoicwa-sicsno(i,j)*rhosnwa)
                   END IF
                   tq = (ddhelp*tho(i,j,k-1)+ddpo(i,j,k)*tho(i,j,k)) / (ddpo(i,j,k)+ddhelp)
                   sq = (ddhelp*sao(i,j,k-1)+ddpo(i,j,k)*sao(i,j,k)) / (ddpo(i,j,k)+ddhelp)

                   tho(i,j,k-1) = tq
                   tho(i,j,k)   = tq
                   sao(i,j,k-1) = sq
                   sao(i,j,k)   = sq

                   IF(kcondep(i,j).EQ.k-1) kcondep(i,j) = kcondep(i,j)+NINT(switc)

                endif

             ELSE IF (ioconv .eq. 3) THEN     !! interchange pot. temperature and salinity of upper and lower box
                                         !! plus correction term for conservation of heat and salt

                vk1e(i,j,k) = switc

                IF (switc .GE. 0.5_wp) THEN   !! static instability
                    tupper = tho(i,j,k-1)
                    tlower = tho(i,j,k)
                    supper = sao(i,j,k-1)
                    slower = sao(i,j,k)
                    ddhelp = ddpo(i,j,k-1)+zsurf(k-1)*(zo(i,j)-sictho(i,j)*rhoicwa &
                             -sicsno(i,j)*rhosnwa)

                    !! set thinner box to values of thicker box
                    !! set thicker box to values of thinner box including correction term
                    !! for heat and salt conservation
                    IF(ddpo(i,j,k).GT.ddhelp) THEN   !! lower box thicker than upper box
                       tho(i,j,k-1)=tlower
                       sao(i,j,k-1)=slower
                       tho(i,j,k)=tlower+(tupper-tlower)*(ddhelp/ddpo(i,j,k))
                       sao(i,j,k)=slower+(supper-slower)*(ddhelp/ddpo(i,j,k))
                    ELSE                             !! lower box thinner than upper box
                       tho(i,j,k)=tupper
                       sao(i,j,k)=supper
                       tho(i,j,k-1)=tupper+(tlower-tupper)*(ddpo(i,j,k)/ddhelp)
                       sao(i,j,k-1)=supper+(slower-supper)*(ddpo(i,j,k)/ddhelp)
                    ENDIF

                 ENDIF

                 IF(kcondep(i,j).EQ.k-1) kcondep(i,j) = kcondep(i,j)+NINT(switc)

             ENDIF     ! end of if test for convection paramterizations

             !! nothing to be done here for ioconv=4, nlopps is called after end of j loop
             !! attention: should we set vk1e in case of ioconv=4 to some value,
             !! e.g. switc as done for ioconv=2,3 or 0 as done for ioconv=1

          ENDDO
       ENDDO

       DO i=1,ie
         shelp(i,j)=sao(i,j,ke)
         thelp(i,j)=tho(i,j,ke)
       ENDDO

       CALL adisitj(thelp,shelp,preffw(kep),j)      !! pot. temperature to in-situ temperature thelp
       CALL rho1j(thelp,shelp,preffw(kep),rhuppo,j)  !! gives in-situ density rhelp

       rhoo(:, j, ke) = 0.5_wp * (rhoo(:, j, ke) + rhuppo(:, j))
       !
    ENDDO ! j-loop
!$omp end do
    CALL bounds_exch(1,'p',rhoo)

!! call plume convection subroutine nlopps if ioconv set to 4
!$omp single
    IF (ioconv.eq.4) THEN
        dtts=dt              !! set nlopps time step to MPIOM time step
        CALL nlopps(tho,sao,dzw,ddpo,preffw,kcondep,dtts)
    ENDIF

    IF (ioconv.lt. 1 .or. ioconv .gt. 4 ) THEN
       WRITE(0,*) 'unsupported ioconv=',ioconv
       CALL P_ABORT
    ENDIF

!$omp end single

!$omp end parallel

  END SUBROUTINE convection


      SUBROUTINE NLOPPS(TA,SA,DZ,DP,PRES,KDEPCON,DTTS)
!---------------------------------------------------------
!
!     NLOPPS:   MODIFIED FROM LOPPS BY E. SKYLLINGSTAD AND T. PALUSZKIEWICZ
!
!     VERSION: DECEMBER 11, 1996
!
!     NLOPPS:  THIS VERSION OF LOPPS IS SIGNIFICANTLY DIFFERENT FROM
!     THE ORIGINAL CODE DEVELOPED BY R. ROMEA AND T. PALUSKIEWICZ.  THE
!     CODE USES A FLUX CONSTRAINT TO CONTROL THE CHANGE IN T AND S AT
!     EACH GRID LEVEL.  FIRST, A PLUME PROFILE OF T,S, AND W ARE
!     DETERMINED USING THE STANDARD PLUME MODEL, BUT WITH A DETRAINING
!     MASS INSTEAD OF ENTRAINING.  THUS, THE T AND S PLUME
!     CHARACTERISTICS STILL CHANGE, BUT THE PLUME CONTRACTS IN SIZE
!     RATHER THAN EXPANDING ALA CLASSICAL ENTRAINING PLUMES.  THIS
!     IS HEURISTICALLY MORE IN LINE WITH LARGE EDDY SIMULATION RESULTS.
!     AT EACH GRID LEVEL, THE CONVERGENCE OF PLUME VELOCITY DETERMINES
!     THE FLUX OF T AND S, WHICH IS CONSERVED BY USING AN UPSTREAM
!     ADVECTION.  THE VERTICAL VELOCITY IS BALANCED SO THAT THE AREA
!     WEIGHTED UPWARD VELOCITY EQUALS THE AREA WEIGHTED DOWNDRAFT
!     VELOCITY, ENSURING MASS CONSERVATION. THE PRESENT IMPLEMENTATION
!     ADJUSTS THE PLUME FOR A TIME PERIOD EQUAL TO THE TIME FOR 1/2 OF
!     THE MASS OF THE FASTEST MOVING LEVEL TO MOVE DOWNWARD.  AS A
!     CONSEQUENCE, THE MODEL DOES NOT COMPLETELY ADJUST THE PROFILE AT
!     EACH MODEL TIME STEP, BUT PROVIDES A SMOOTH ADJUSTMENT OVER TIME.
!
!
!---------------------------------------------------------
!
      USE MO_PARAM1
      USE mo_planetary_constants, ONLY: g
!
      REAL(wp) TA(IE,JE,KE),SA(IE,JE,KE),DP(IE,JE,KE)
      REAL(wp) DZ(KE),PRES(KE)
      REAL(wp) THELP(KE),SHELP(KE),THELP1(KE),SHELP1(KE)
      REAL(wp) THELP2(KE),SHELP2(KE)
!
      REAL(wp) TTEMP(KE),STEMP(KE),TAA(KE),SAA(KE)
      REAL(wp) TDA(KE),SDA(KE),MDA(KE)!, WDA(KE)
      REAL(wp) AD(KE),SD(KE),TD(KE),WD(KE),MD(KE)
      REAL(wp) DE(KE),DD(KE)
      REAL(wp) PLUMEENTRAINMENT(KE)
      REAL(wp) GRIDTHICKNESS(KE)
      REAL(wp) EPPS
      REAL(wp) WET

!
      REAL(wp) WSQR,RADIUS
      REAL(wp) SMIX,THMIX
      REAL(wp) D1,D2
      REAL(wp) DZ1,DZ2
      REAL(wp) STARTINGFLUX,OLDFLUX,NEWFLUX,ENTRAINRATE
      REAL(wp) DTTS,DT
      INTEGER NTIME,NN,KMX,IC
      INTEGER KDEPCON(IE,JE) !SJM DEPTH OF CONVECTION

      INTEGER i         !< loop variable in W-E direction
      INTEGER j         !< loop variable in N-S direction
      INTEGER k         !< loop variable for vertical levels (k=1 top, k=ke bottom)
      INTEGER KMAX
      INTEGER MAXDEPTH
      INTEGER K2
!
!
! INPUT THE VARIABLES THROUGH A COMMON
!
!
!      LOGICAL DEBUG,DONE,PROBLEM
      INTEGER, PARAMETER :: max_abe_iterations=1

!SJ***SENS.PCN: CHANGE THE PLUME RADIUS INTO 700 M
!      PARAMETER ( PLUMERADIUS          =  700.D0   )
      REAL(wp), PARAMETER :: plumeradius = 500.e0_wp
      REAL(wp), PARAMETER :: stability_threshold = -1.e-4_wp
      REAL(wp), PARAMETER :: fractional_area = .1E0_wp
      REAL(wp), PARAMETER :: vertical_velocity = .03E0_wp
      REAL(wp), PARAMETER :: entrainment_rate = -.05E0_wp
      REAL(wp), PARAMETER :: e2 = 2.E0_wp * entrainment_rate
!
!
!-----MAY WANT TO SETUP AN OPTION TO GET THIS ONLY ON FIRST CALL
!     OTHERWISE IT IS REPETIVE
!     GRIDDZ IS INITIALIZE BY CALL TO SETUPGRID
!
!      DTTS=2400.
!SJM   DTTS = 72000.
!      DTTS=1920.   !!!CSJM
!
        DO I=1,IE
        DO J=1,JE
        KDEPCON(I,J)=0
        ENDDO
        ENDDO

        DO K=1,KE
           GRIDTHICKNESS(K) = DZ(K)
        ENDDO
!
!
! MODIFIED TO LOOP OVER SLAB
!
      DO 10 J=1,JE
!
      DO 100 I=1,IE
!
      KMAX=0
      DO K=1,KE
         epps = 1.E-30_wp
         wet = MAX(0._wp, dp(i, j, k)/(dp(i, j, k) - epps))
         IF (WET .NE. 0._wp) KMAX=K
      ENDDO
!
      IF(KMAX.LE.1) GOTO 100
!
      DO K=1,KMAX
         STEMP(K)=SA(I,J,K)
         TTEMP(K)=TA(I,J,K)
!       IF(I.EQ.5.AND.J.EQ.35)PRINT*,'OLD=',STEMP(K),K
         SHELP1(K)=STEMP(K)
         THELP1(K)=TTEMP(K)
      ENDDO
!
      DO K=1,KMAX-1
! INITIALIZE THE PLUME T,S,DENSITY, AND W VELOCITY
!
          SD(K)=STEMP(K)
          TD(K)=TTEMP(K)
!
          SHELP(K)=STEMP(K)
          THELP(K)=TTEMP(K)
!
          CALL ADISIT1(THELP1(K),SHELP1(K),PRES(K))
          CALL RHO2(THELP1(K),SHELP1(K),PRES(K),DD(K))
          DE(K)=DD(K)
!
          WD(K)=VERTICAL_VELOCITY
! GUESS AT INITIAL TOP GRID CELL VERTICAL VELOCITY
!
!          WD(K) = 0.03
! THESE ESTIMATES OF INITIAL PLUME VELOCITY BASED ON PLUME SIZE AND
! TOP GRID CELL WATER MASS
!          WD(K) = 0.5*DZ(K)/(DTTS*FRACTIONAL_AREA)
!          WD(K) = 0.5*DZ(K)/DTTS
!
          WSQR=WD(K)*WD(K)
          plumeentrainment(k) = 0.0_wp
!
          RADIUS=PLUMERADIUS
          STARTINGFLUX=RADIUS*RADIUS*WD(K)*DD(K)
          OLDFLUX=STARTINGFLUX
!
          DZ2=GRIDTHICKNESS(K)
!
          DO K2=K,KMAX-1
!  CALCULATE DENSITY FOR UPPER LAYER
            CALL ADISIT1(THELP(K2),SHELP(K2),PRES(K2+1))
            CALL RHO2(THELP(K2),SHELP(K2),PRES(K2+1),D1)
!  CALCULATE DENSITY FOR LOWER LAYER
            CALL ADISIT1(THELP1(K2+1),SHELP1(K2+1),PRES(K2+1))
            CALL RHO2(THELP1(K2+1),SHELP1(K2+1),PRES(K2+1),D2)
!
            DE(K2+1)=D2
!
! TO START DOWNWARD, PARCEL HAS TO INITIALLY BE HEAVIER THAN ENVIRONMENT
! BUT AFTER IT HAS STARTED MOVING, WE CONTINUE PLUME UNTIL PLUME TKE OR
! FLUX GOES NEGATIVE
!
!SJ***SENS.CONV_NS:
!          IF(J.GE.34) STABILITY_THRESHOLD = -1.E-5
!
            IF (D2-D1 .LT. STABILITY_THRESHOLD.OR.K2.NE.K) THEN
                 DZ1=DZ2
                 DZ2=GRIDTHICKNESS(K2+1)
!
! DEFINE MASS FLUX ACCORDING TO EQ. 4 FROM PAPER
                 newflux = oldflux + e2 * radius * wd(k2) * dd(k2) * 0.50_wp * &
     &              (DZ1+DZ2)
!
                 PLUMEENTRAINMENT(K2+1) = NEWFLUX/STARTINGFLUX
!
!SJ***SENS.CONV_NS: MODIFIED BY SJKIM
                 IF (newflux .LT. 1000.0_wp) THEN
                     MAXDEPTH = K2
                     IF(MAXDEPTH.EQ.K) GOTO 1000
                     GOTO 1
                 ENDIF
!
! ENTRAINMENT RATE IS BASICALLY A SCALED MASS FLUX DM/M
!
                 ENTRAINRATE = (NEWFLUX - OLDFLUX)/NEWFLUX
                 OLDFLUX = NEWFLUX
!
!
! MIX VAR'S ARE THE AVERAGE ENVIRONMENTAL VALUES OVER THE TWO GRID LEVELS
!
                 SMIX=(DZ1*STEMP(K2)+DZ2*STEMP(K2+1))/(DZ1+DZ2)
                 THMIX=(DZ1*TTEMP(K2)+DZ2*TTEMP(K2+1))/(DZ1+DZ2)
!
! FIRST COMPUTE THE NEW SALINITY AND TEMPERATURE FOR THIS LEVEL
! USING EQUATIONS 3.6 AND 3.7 FROM THE PAPER
!
!
!
                  SD(K2+1)=SD(K2) - ENTRAINRATE*(SMIX - SD(K2))
                  TD(K2+1)=TD(K2) - ENTRAINRATE*(THMIX - TD(K2))
!
        IF (sd(k2+1) .LE. 0._wp) PRINT *, i, j, k, k2, entrainrate, newflux
!
                  SHELP2(K2+1)=SD(K2+1)
                  THELP2(K2+1)=TD(K2+1)
!
!
! COMPUTE THE DENSITY AT THIS LEVEL FOR THE BUOYANCY TERM IN THE
! VERTICAL K.E. EQUATION
!
!
           CALL ADISIT1(THELP2(K2+1),SHELP2(K2+1),PRES(K2+1))
           CALL RHO2(THELP2(K2+1),SHELP2(K2+1),PRES(K2+1),DD(K2+1))
!
! NEXT, SOLVE FOR THE VERTICAL VELOCITY K.E. USING COMBINED EQ. 4
! AND EQ 5 FROM THE PAPER
!
!
                 WSQR = WSQR - WSQR*ABS(ENTRAINRATE)+ g *             &
     &             (DZ1*(DD(K2)-DE(K2))/DE(K2)                          &
     &             +DZ2*(DD(K2+1)-DE(K2+1))/DE(K2+1))
!
! IF NEGATIVE K.E. THEN PLUME HAS REACHED MAX DEPTH, GET OUT OF LOOP
!
                 IF (wsqr .LT. 0.0_wp) THEN
                     MAXDEPTH = K2
                     IF(MAXDEPTH.EQ.K) GOTO 1000
                     GOTO 1
                 ENDIF
                 WD(K2+1)=SQRT(WSQR)
!
! COMPUTE A NEW RADIUS BASED ON THE NEW MASS FLUX AT THIS GRID LEVEL
                 RADIUS=SQRT(NEWFLUX/(WD(K2)*DD(K2)))
              ELSE
                 MAXDEPTH=K2
                 IF(MAXDEPTH.EQ.K) GOTO 1000
                 GOTO 1
              ENDIF
          ENDDO
!
! PLUME HAS REACHED THE BOTTOM
!
          MAXDEPTH=KMAX
!
 1         CONTINUE
!
          AD(K)=FRACTIONAL_AREA
          IC=0
!
! START ITERATION ON FRACTIONAL AREA, NOT USED IN OGCM IMPLEMENTATION
!
!
!
          DO IC=1,MAX_ABE_ITERATIONS
!
!
! NEXT COMPUTE THE MASS FLUX BETWEEN EACH GRID BOX USING THE ENTRAINMENT
!
             MD(K)=WD(K)*AD(K)
!
             DO K2=K+1,MAXDEPTH
               MD(K2)=MD(K)*PLUMEENTRAINMENT(K2)
             ENDDO
!
! NOW MOVE ON TO CALCULATE NEW TEMPERATURE USING FLUX FROM
! TD, SD, WD, TA, SA, AND WE. VALUES FOR THESE VARIABLES ARE AT
! CENTER OF GRID CELL, USE WEIGHTED AVERAGE TO GET BOUNDARY VALUES
!
! USE A TIMESTEP LIMITED BY THE GCM MODEL TIMESTEP AND THE MAXIMUM PLUME
! VELOCITY (CFL CRITERIA)
!
!
! CALCULATE THE WEIGHTED WD, TD, AND SD
!
             DT = DTTS
             DO K2=K,MAXDEPTH-1
                DT = MIN(DT,DZ(K2)/WD(K2))
!
! TIME INTEGRATION WILL BE INTEGER NUMBER OF STEPS TO GET ONE
! GCM TIME STEP
!
                NTIME = NINT(0.5*REAL(INT(DTTS/DT)))
                IF(NTIME.EQ.0) THEN
                   NTIME = 1
                ENDIF
!
! MAKE SURE AREA WEIGHTED VERTICAL VELOCITIES MATCH; IN OTHER WORDS
! MAKE SURE MASS IN EQUALS MASS OUT AT THE INTERSECTION OF EACH GRID
! CELL.
!
                MDA(K2) = (MD(K2)*DZ(K2)+MD(K2+1)*DZ(K2+1))/            &
     &                    (DZ(K2)+DZ(K2+1))
!
!                WDA(K2) = (WD(K2)*DZ(K2)+WD(K2+1)*DZ(K2+1))/
!     *                    (DZ(K2)+DZ(K2+1))
!
                TDA(K2) = TD(K2)
                SDA(K2) = SD(K2)
!
                TAA(K2) = TTEMP(K2+1)
                SAA(K2) = STEMP(K2+1)
!
             ENDDO
!
             DT = MIN(DT,DTTS)
!
             TDA(MAXDEPTH) = TD(MAXDEPTH)
             SDA(MAXDEPTH) = SD(MAXDEPTH)
!
! DO TOP AND BOTTOM POINTS FIRST
!
             KMX = MAXDEPTH-1
!
             DO NN=1,NTIME
!
               TTEMP(K) =  TTEMP(K)-                                    &
     &                  (MDA(K)*(TDA(K)-TAA(K)))*DT/DZ(K)
!
               STEMP(K) =  STEMP(K)-                                    &
     &                  (MDA(K)*(SDA(K)-SAA(K)))*DT/DZ(K)
!
!
! NOW DO INNER POINTS IF THERE ARE ANY
!
               IF(MAXDEPTH-K.GT.1) THEN
                 DO K2=K+1,MAXDEPTH-1
!
                   TTEMP(K2) = TTEMP(K2) +                              &
     &              (MDA(K2-1)*(TDA(K2-1)-TAA(K2-1))-                   &
     &              MDA(K2)*(TDA(K2)-TAA(K2)))                          &
     &              *DT/DZ(K2)
!
!
                  STEMP(K2) = STEMP(K2) +                               &
     &              (MDA(K2-1)*(SDA(K2-1)-SAA(K2-1))-                   &
     &              MDA(K2)*(SDA(K2)-SAA(K2)))                          &
     &              *DT/DZ(K2)
!
                 ENDDO
               ENDIF
!
               TTEMP(KMX+1) =  TTEMP(KMX+1)+                            &
     &                  (MDA(KMX)*(TDA(KMX)-TAA(KMX)))*                 &
     &                  DT/DZ(KMX+1)
!
               STEMP(KMX+1) =  STEMP(KMX+1)+                            &
     &                  (MDA(KMX)*(SDA(KMX)-                            &
     &                  SAA(KMX)))*DT/DZ(KMX+1)
!
! SET THE ENVIRONMENTAL TEMP AND SALINITY TO EQUAL NEW FIELDS
!
                DO K2=1,MAXDEPTH-1
                  TAA(K2) = TTEMP(K2+1)
                  SAA(K2) = STEMP(K2+1)
                ENDDO
!
! END LOOP ON NUMBER OF TIME INTEGRATION STEPS
!
             ENDDO
          ENDDO
!
! ASSUME THAT IT CONVERGED, SO UPDATE THE TA AND SA WITH NEW FIELDS
!
          DO K2=K,MAXDEPTH
            SA(I,J,K2) = STEMP(K2)
            TA(I,J,K2) = TTEMP(K2)
          ENDDO
!SJM
      KDEPCON(I,J) = MAXDEPTH
!
! JUMP HERE IF K = MAXDEPTH OR IF LEVEL NOT UNSTABLE, GO TO NEXT
! PROFILE POINT
!
 1000     CONTINUE
!
! END LOOP ON K, MOVE ON TO NEXT POSSIBLE PLUME
!
      ENDDO
!
! I LOOP
!
 100  CONTINUE
!
! J LOOP
  10  CONTINUE
      RETURN
      END SUBROUTINE NLOPPS

 END MODULE mo_convection
