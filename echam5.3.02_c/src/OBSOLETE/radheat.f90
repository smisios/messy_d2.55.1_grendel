SUBROUTINE radheat (kproma, kbdim, klev,  klevp1                       &
                  , krow                                               &
                  , pi0                                                &
                  , ptm1,         pqm1                                 &
                  , ptrsof,       ptrsol                               &
                  , pemtef,       pemter                               &
                  , pemtef0,      ptrsof0                              &
                  , ptrfll,       ptrflw,      ptrfli                  &
                  , psofll,       psoflw,      psofli                  &
                  , ptrfllac,     ptrflwac,    ptrfliac                &
                  , psofllac,     psoflwac,    psofliac                &
                  , psrad0,       psrads                               &
                  , psradl,       psrafl                               &
                  , psrad0u,      psradsu                              &
                  , psraf0,       psrafs                               &
                  , psrad0d                                            &
                  , ptrad0,       ptrads                               &
                  , ptradl,       ptrafl                               &
                  , ptraf0,       ptrafs                               &
                  , ptradsu                                            &
                  , ptslm1,       ptsi,        ptsw                    &
                  , palbedo                                            &
                  , palsol,       palsow,      palsoi                  &
                  , paphm1,       papm1                                &
                  , ptslnew,      ptte                                 &
                  , pfrl,         pfrw,        pfri   )
!
!
!
!**** *RADHEAT* - COMPUTES TEMPERATURE CHANGES DUE TO RADIATION.
!
!
!     SUBJECT.
!     --------
!
!          THIS ROUTINE COMPUTES THE TENDENCIES OF THE ATMOSPHERE'S
!     TEMPERATURE DUE TO THE EFFECTS OF LONG WAVE AND SHORT WAVE
!     RADIATION. THE COMPUTATION IS DONE ON THE T-1 TIME LEVEL USING
!     VALUES OF ATMOSPHERIC TRANSMISIVITIES AND EMISSIVITIES THAT HAVE
!     BEEN STORED AT THE LAST FULL RADIATION TIME STEP. THE SURFACE
!     SOLAR FLUX LATER TO BE USED IN THE SOIL PROCESS CALCULATIONS IS
!     ALSO STORED.
!
!**   INTERFACE.
!     ----------
!
!          *RADHEAT* IS CALLED FROM *PHYSC*.
!
!     INPUT ARGUMENTS.
!     ----- ---------
!
!
!     OUTPUT ARGUMENTS.
!     ------ ---------
!
!
!     METHOD.
!     -------
!
!     PRODUCT OF SOLAR
!     INFLUX BY TRANSMISSIVITIES LEADS TO SOLAR FLUXES. THEN THE
!     TEMPERATURES ARE INTERPOLATED/EXTRAPOLATED TO THE LAYER BOUNDARIES
!     (AT THE BOTTOM ONE TAKES THE SURFACE TEMPERATURE) AND A PRODUCT BY
!     EMISSIVITIES OF SIGMA*T**4 GIVES THERMAL FLUXES. THE TWO FLUXES
!     ARE ADDED AND DIVERGENCES COMPUTED TO GIVE HEATING RATES.
!
!     EXTERNALS.
!     ----------
!
!          *SOLANG*.
!
!     AUTHOR.
!     ------
!
!     U. SCHLESE    DKRZ-HAMBURG    JUNE 1995

!     Modifications
!     U. Schlese, December 1999:  version for coupling
!     U. Schlese, July 2000, *solang* removed, weighted surface fluxes
!     I. Kirchner, May 2002, tendency diagnose bugfix for surface fluxes
!
USE mo_kind,              ONLY: dp
USE mo_control,           ONLY: ltdiag
USE mo_constants,         ONLY: g, cpd, vtmpc2, stbo
USE mo_radiation,         ONLY: cemiss
USE mo_diag_tendency,     ONLY: pdiga
USE mo_time_control,      ONLY: delta_time
!
IMPLICIT NONE
!
INTEGER :: kproma,kbdim,klev,klevp1,krow

REAL(dp) ::                                                            &
        pi0(kbdim),                                                    &
        ptm1(kbdim,klev),      pqm1(kbdim,klev),                       &
        ptrsol(kbdim,klevp1),  ptrsof(kbdim,2),                        &
        pemter(kbdim,klevp1),  pemtef(kbdim,2),                        &
        ptrsof0(kbdim,klevp1), pemtef0(kbdim,klevp1),                  &
        ptrfll(kbdim),         ptrflw(kbdim),         ptrfli(kbdim),   &
        psofll(kbdim),         psoflw(kbdim),         psofli(kbdim),   &
        ptrfllac(kbdim),       ptrflwac(kbdim),       ptrfliac(kbdim), &
        psofllac(kbdim),       psoflwac(kbdim),       psofliac(kbdim), &
        psrad0(kbdim),         psrads(kbdim),                          &
        psradl(kbdim),         psrafl(kbdim),                          &
        psrad0u(kbdim),        psradsu(kbdim),                         &
        psraf0(kbdim),         psrafs(kbdim),                          &
        psrad0d(kbdim),                                                &
        ptrad0(kbdim),         ptrads(kbdim),                          &
        ptradl(kbdim),         ptrafl(kbdim),                          &
        ptraf0(kbdim),         ptrafs(kbdim),                          &
        ptradsu(kbdim),                                                &
        ptslm1(kbdim),         ptsi(kbdim),           ptsw(kbdim),     &
        palbedo(kbdim),                                                &
        palsol(kbdim),         palsow(kbdim),         palsoi(kbdim),   &
        paphm1(kbdim,klevp1),  papm1(kbdim,klev),                      &
        ptslnew(kbdim),        ptte(kbdim,klev),                       &
        pfrl(kbdim),           pfrw(kbdim),           pfri(kbdim)
!
!
!    Local arrays
!
REAL(dp) :: zti(kbdim,klevp1),                                         &
            zflxs(kbdim,klevp1),  zflxt(kbdim,klevp1),                 &
            zflxs0(kbdim,klevp1), zflxt0(kbdim,klevp1),                &
            ztsnew(kbdim),        zteffl4(kbdim)
REAL(dp) :: ztrps(kbdim),  ztrpt(kbdim),  ztrpss(kbdim), ztrpts(kbdim)
REAL(dp) :: zdtime, zcons3, zdtdt, zfltop, zflbot, ztrdown, ztrdown0,  &
            zflts, zfltt, zffact, zsr0u, zsrsu, zdp1, zdp2, ztrsu
INTEGER :: jrow, jk, jl
!
! ----------------------------------------------------------------------
!
!*     1.   COMPUTATIONAL CONSTANTS.
!           ------------- ----------
!
100 CONTINUE
!
!
  zcons3=g/cpd
  zdtime = delta_time
!
!
  jrow = krow
!
!     ------------------------------------------------------------------
!
!*         3.     TEMPERATURES AT LAYERS' BOUDARIES.
!                 ------------ -- ------- ----------
!
300  CONTINUE
!
!*         3.1     INTERPOLATION PROPER.
!
310  CONTINUE
     DO 312 jk=2,klev
        DO 311 jl=1,kproma
           zti(jl,jk)=(ptm1(jl,jk-1)*papm1(jl,jk-1)                    &
                     *(papm1(jl,jk)-paphm1(jl,jk))                     &
                     +ptm1(jl,jk)*papm1(jl,jk)                         &
                     *(paphm1(jl,jk)-papm1(jl,jk-1)))                  &
                     /(paphm1(jl,jk)*(papm1(jl,jk)-papm1(jl,jk-1)))
311     END DO
312  END DO
!
!*        3.2     SURFACE AND TOP OF ATMOSPHERE TEMPERATURE.
!
320  CONTINUE
     DO 321 jl=1,kproma
!  -  fractional surface coverage:
        zti(jl,klevp1)=(pfrl(jl)*ptslm1(jl)**4                         &
                       +pfri(jl)*ptsi(jl)**4                           &
                       +pfrw(jl)*ptsw(jl)**4)**0.25_dp
        zteffl4(jl)=ptslm1(jl)**3*(4._dp*ptslnew(jl)-3._dp*ptslm1(jl))
        ztsnew(jl)=(pfrl(jl)*zteffl4(jl)                               &
                   +pfri(jl)*ptsi(jl)**4                               &
                   +pfrw(jl)*ptsw(jl)**4)**0.25_dp
        zti(jl,1)=ptm1(jl,1)-papm1(jl,1)*(ptm1(jl,1)-zti(jl,2))        &
                  /(papm1(jl,1)-paphm1(jl,2))
321  END DO
!
!     ------------------------------------------------------------------
!
!*         4.    UPDATE FLUXES AND COMPUTE HEATING RATES.
!                ------ ------ --- ------- ------- -----
!
!    4.1 Fluxes at top of the atmosphere
!
     DO 401 jl=1,kproma
       zflxs(jl,1)=pi0(jl)*ptrsol(jl,1)
       zflxt(jl,1)=pemter(jl,1)
       zflxs0(jl,1)=pi0(jl)*ptrsof0(jl,1)
       zflxt0(jl,1)=pemtef0(jl,1)
401  END DO
!
!
!     4.2  Fluxes and heating rates except for lowest layer
!
      DO 403 jk=1,klev-1
        DO 402 jl=1,kproma
          zfltop=zflxs(jl,jk)+zflxt(jl,jk)
          zflxs(jl,jk+1)=pi0(jl)*ptrsol(jl,jk+1)
          zflxt(jl,jk+1)=pemter(jl,jk+1)
          zflbot=zflxs(jl,jk+1)+zflxt(jl,jk+1)
          zdtdt=-zcons3*(zflbot-zfltop)/((paphm1(jl,jk+1)              &
                -paphm1(jl,jk))*(1._dp+vtmpc2*pqm1(jl,jk)))
          ptte(jl,jk)=ptte(jl,jk)+zdtdt
!
          zflxs0(jl,jk+1)=pi0(jl)*ptrsof0(jl,jk+1)
          zflxt0(jl,jk+1)=pemtef0(jl,jk+1)

402     END DO
403  END DO
!
!
!     4.3  Lowest layer
!
     DO 404 jl=1,kproma      
       ztrdown=pemter(jl,klevp1)+cemiss*stbo*zti(jl,klevp1)**4
       zflxt(jl,klevp1)=ztrdown-cemiss*stbo*ztsnew(jl)**4
       ztrdown0=pemtef0(jl,klevp1)+cemiss*stbo*zti(jl,klevp1)**4
       zflxt0(jl,klevp1)=ztrdown0-cemiss*stbo*ztsnew(jl)**4
       ptrfll(jl)=ztrdown-cemiss*stbo*zteffl4(jl)
       ptrflw(jl)=ztrdown-cemiss*stbo*ptsw(jl)**4
       ptrfli(jl)=ztrdown-cemiss*stbo*ptsi(jl)**4
       zflxs(jl,klevp1)=pi0(jl)*ptrsol(jl,klevp1)
       zflxs0(jl,klevp1)=pi0(jl)*ptrsof0(jl,klevp1)
       psofll(jl)=(1._dp-palsol(jl))*zflxs(jl,klevp1)                  &
                                                   /(1._dp-palbedo(jl))
       psoflw(jl)=(1._dp-palsow(jl))*zflxs(jl,klevp1)                  &
                                                   /(1._dp-palbedo(jl))
       psofli(jl)=(1._dp-palsoi(jl))*zflxs(jl,klevp1)                  &
                                                   /(1._dp-palbedo(jl))
       zfltop=zflxs(jl,klev)+zflxt(jl,klev)
       zflbot=zflxs(jl,klevp1)+zflxt(jl,klevp1)
       zdtdt=-zcons3*(zflbot-zfltop)/((paphm1(jl,klevp1)               &
                 -paphm1(jl,klev))*(1._dp+vtmpc2*pqm1(jl,klev)))
       ptte(jl,klev)=ptte(jl,klev)+zdtdt
404  END DO
!
     IF (ltdiag) THEN
       ! tendency diagnostics
       DO jk = 1, klev
         DO jl = 1,kproma
           zflts = zflxs(jl,jk+1)-zflxs(jl,jk)
           zfltt = zflxt(jl,jk+1)-zflxt(jl,jk)
           zffact = - zcons3/((paphm1(jl,jk+1)-paphm1(jl,jk)) &
                * (1._dp+vtmpc2*pqm1(jl,jk)) )
           pdiga(jl,jk,24,jrow) = pdiga(jl,jk,24,jrow) + zfltt*zffact
           pdiga(jl,jk,25,jrow) = pdiga(jl,jk,25,jrow) + zflts*zffact
         END DO
       END DO
     END IF
 
!     ------------------------------------------------------------------
!
!*         5.     Diagnostics of top and surface fluxes
!
!
!
     DO 510 jl = 1, kproma
       psrad0(jl) = psrad0(jl) + zdtime*zflxs(jl,1)
       ptrad0(jl) = ptrad0(jl) + zdtime*zflxt(jl,1)
       zsr0u = zflxs(jl,1) - pi0(jl)
       psrad0u(jl) = psrad0u(jl) + zdtime*zsr0u
       psrads(jl) = psrads(jl) + zdtime*zflxs(jl,klevp1)
       ptrads(jl) = ptrads(jl) + zdtime*zflxt(jl,klevp1)
       zsrsu = -zflxs(jl,klevp1)*(1._dp/(1._dp-palbedo(jl))-1.)
       psradsu(jl) = psradsu(jl) + zdtime*zsrsu
       psrad0d(jl) = psrad0d(jl) + zdtime*pi0(jl)
       psraf0(jl) = psraf0(jl) + zdtime*pi0(jl)*ptrsof(jl,1)
       psrafs(jl) = psrafs(jl) + zdtime*pi0(jl)*ptrsof(jl,2)
510  END DO
!
!
! Diagnostics of fluxes at 200mb
!
  DO 524 jk=1,klev
     DO 523 jl=1,kproma
        IF (paphm1(jl,jk) .LE. 20000._dp .AND.                         &
            paphm1(jl,jk+1) .GE. 20000._dp) THEN
            zdp1=paphm1(jl,jk)-paphm1(jl,jk+1)
            zdp2=paphm1(jl,jk)-20000._dp
            ztrps(jl)=zflxs(jl,jk)-(zflxs(jl,jk)-zflxs(jl,jk+1))       &
                        *(zdp2/zdp1)
            ztrpt(jl)=zflxt(jl,jk)-(zflxt(jl,jk)-zflxt(jl,jk+1))       &
                        *(zdp2/zdp1)
            ztrpss(jl)=zflxs0(jl,jk)-(zflxs0(jl,jk)-zflxs0(jl,jk+1))   &
                        *(zdp2/zdp1)
            ztrpts(jl)=zflxt0(jl,jk)-(zflxt0(jl,jk)-zflxt0(jl,jk+1))   &
                        *(zdp2/zdp1)
        END IF
523 END DO
524 END DO
  DO 528 jl=1,kproma
     psradl(jl)  =psradl(jl) + ztrps(jl)*zdtime
     ptradl(jl)  =ptradl(jl) + ztrpt(jl)*zdtime
     psrafl(jl)  =psrafl(jl) + ztrpss(jl)*zdtime
     ptrafl(jl)  =ptrafl(jl) + ztrpts(jl)*zdtime
528 END DO
!
!
     DO 520 jl=1,kproma
       ztrsu = -cemiss*stbo*ztsnew(jl)**4+(cemiss-1._dp)               &
               *(stbo*zti(jl,klevp1)**4+pemter(jl,klevp1)/cemiss)
       ptradsu(jl)=ptradsu(jl)+zdtime*ztrsu
       ptraf0(jl) =ptraf0(jl) +zdtime*pemtef(jl,1)
       ptrafs(jl) =ptrafs(jl) +zdtime*(pemtef(jl,2)+cemiss*stbo        &
                                  *(zti(jl,klevp1)**4-ztsnew(jl)**4))
520  END DO
!
!        Area weighted surface fluxes over land/water/ice
!
     DO 530 jl=1,kproma
       psofllac(jl)=psofllac(jl)+zdtime*psofll(jl)*pfrl(jl)
       psoflwac(jl)=psoflwac(jl)+zdtime*psoflw(jl)*pfrw(jl)
       psofliac(jl)=psofliac(jl)+zdtime*psofli(jl)*pfri(jl)
       ptrfllac(jl)=ptrfllac(jl)+zdtime*ptrfll(jl)*pfrl(jl)
       ptrflwac(jl)=ptrflwac(jl)+zdtime*ptrflw(jl)*pfrw(jl)
       ptrfliac(jl)=ptrfliac(jl)+zdtime*ptrfli(jl)*pfri(jl)
530  END DO     
!
  RETURN
END SUBROUTINE radheat
