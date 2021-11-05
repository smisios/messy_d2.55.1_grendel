!=======================================================================
!
! *** BLOCK DATA BLKPAR
! *** THIS SUBROUTINE PROVIDES INITIAL (DEFAULT) VALUES TO PROGRAM
!     PARAMETERS VIA DATA STATEMENTS
!
! *** WRITTEN BY ATHANASIOS NENES
! *** MODIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES
!
!=======================================================================
!

      BLOCK DATA BLKPAR
!
      INCLUDE 'messy_cloud_parameters_uaf.inc'
!
      DATA AMA    /29d-3/               ! Air molecular weight
      DATA GRAV   /9.81d0/              ! g constant
      DATA RGAS   /8.31d0/              ! Universal gas constant
      DATA Dw     /2.75d-10/            ! Water Molecule Diameter
      DATA AMW    /18d-3/               ! Water molecular weight
      DATA DENW   /1d3/                 ! Water density
      DATA DHV    /2.25d6/              ! Water enthalpy of vaporization
      DATA CPAIR  /1.0061d3/            ! Air Cp

!     Data for FHH exponent calculation
   
      DATA D11   /-0.1907/
      DATA D12   /-1.6929/
      DATA D13   /1.4963/
      DATA D14   /-0.5644/ 
      DATA D15   /0.0711/
      ! for C2
      DATA D21   /-3.9310/
      DATA D22   /7.0906/
      DATA D23   /-5.3436/
      DATA D24   /1.8025/ 
      DATA D25   /-0.2131/
      ! for C3
      DATA D31   /8.4825/
      DATA D32   /-14.9297/
      DATA D33   /11.4552/
      DATA D34   /-3.9115/ 
      DATA D35   /0.4647/
      ! for C4
      DATA D41   /-5.1774/
      DATA D42   /8.8725/
      DATA D43   /-6.8527/
      DATA D44   /2.3514/ 
      DATA D45   /-0.2799/
!
      DATA MAXIT   /100/                 ! Max iterations for solution
      DATA EPSC     /1d-6/                ! Convergence criterion
!
      DATA PI      /3.1415927d0/         ! Some constants
      DATA ZERO    /0d0/
      DATA GREAT   /1D30/
      DATA SQ2PI   /2.5066282746d0/
!
      DATA CCNSPST /.FALSE./             ! Internal consistency check
!
! *** END OF BLOCK DATA SUBPROGRAM *************************************
!
      END
!=======================================================================
!=======================================================================
!
! *** SUBROUTINE CCNSPEC
! *** THIS SUBROUTINE CALCULATES THE CCN SPECTRUM OF THE AEROSOL USING
!     THE APPROPRIATE FORM OF KOHLER THEORY
!
! *** ORIGINALLY WRITTEN BY ATHANASIOS NENES FOR ONLY KOHLER PARTICLES
! *** MODIFIED BY PRASHANT KUMAR AND ATHANSIOS NENES TO INCLUDE 
! *** ACTIVATION BY FHH PARTICLES
!
!=======================================================================
!
      SUBROUTINE CCNSPEC (TPI,DPGI,SIGI,TPARC,PPARC,NMODES, &
     &                    AKKI,ei,A,B,SG)
!
      INCLUDE 'messy_cloud_parameters_uaf.inc'
      DOUBLE PRECISION DPGI(NMODES),SIGI(NMODES),TPI(NMODES), &
     & AKKI(NMODES),TP(NSMX),AKK(NSMX),SG(NSMX),ei(NSMX),A,B,TPARC,PPARC
!
      DOUBLE PRECISION Dpcm

      NMD  = NMODES                ! Save aerosol params in COMMON
      DO I=1,NMD

!     For soluble aerosol -- Kohler modes

        IF (MODE(I).EQ.1) THEN
          DPG(I) = MAX(DPGI(I),1.d-9)
          SIG(I) = SIGI(I)
          AKK(I) = MAX(AKKI(I),1.d-3)
          TP(I)  = TPI(I)

!     For insoluble aerosol -- FHH modes (absorption activation)

        ELSEIF (MODE(I).EQ.2) THEN
          DPG(I) = MAX(DPGI(I),1.d-7)/ei(I)**(1./3.)
          SIG(I) = SIGI(I)
          TP(I)  = TPI(I)
          AKK(I) = MAX(AKKI(I),1.d-3)
        ENDIF
      ENDDO
!C
      TEMPER = TPARC                                ! Save parcel props in COMMON
      PRES   = PPARC
      CALL PROPS                                    ! Thermophysical properties
      AKOH   = 4D0*AMW*SURT/RGAS/TEMPER/DENW        ! Kelvin parameter
!C
      DO K=1,NMD
        IF (MODE(K).EQ.1) THEN                    ! Kohler modes
           PAR1   = 4D0/27D0/AKK(K)/DPG(K)**3         
           PAR2   = SQRT(PAR1*AKOH**3)
           SG(K)  = EXP(PAR2) - 1D0                     
        ELSEIF (MODE(K).EQ.2) THEN                ! FHH modes
           CALL DpcFHH(DPG(K),TPARC,AKK(K),ei(K),A,B,Dpcm)
           Dpc(K) = Dpcm
        SG(K)  = (AKOH/Dpc(K))+ &
     &  (-(1.-ei(K))*DPG(K)**3.*AKK(K)/(Dpc(K)**3.-ei(K)*DPG(K)**3.))+ &
     &  (-A*(((Dpc(K)-ei(K)**(1./3.)*DPG(K))/(2*Dw))**(-B)))
        ENDIF
      ENDDO   
!C
!C *** END OF SUBROUTINE CCNSPEC ****************************************
!C
      RETURN
      END
!C=======================================================================
!C=======================================================================
!C
!C *** SUBROUTINE DpcFHH
!C *** THIS SUBROUTINE CALCULATES THE CRITICAL PARTICLE DIAMETER
!C     ACCORDING TO THE FHH ADSOSPRTION ISOTHERM THEORY.
!C
!C *** WRITTEN BY PRASHANT KUMAR AND ATHANASIOS NENES
!C *** UPDATED BY VLASSIS KARYDIS AND ATHANASIOS NENES
!C=======================================================================
!C
      SUBROUTINE DpcFHH(Ddry,TPARC,dhk,dei,A,B,Dc)
!C
      Include 'messy_cloud_parameters_uaf.inc'
      DOUBLE PRECISION Ddry,Dpcm,Dpcl,Dpcu,&
     &FDpcl,FDpcu,FDpcm,Dc,A,B,dei,dhk
       INTEGER N
       PARAMETER (ITMAX=100)
               
      TEMPER = TPARC
      CALL PROPS


            Dpcl = Ddry         !Lower Limit
            Dpcu = 1000.*Ddry   !Upper Limit
         N=0
 100      N=N+1
         IF(N.gt.ITMAX) then 
!        write(6,*) 'Dpcm, ITMAX too small'
!        write(6,*) 'dry',Ddry,'kappa',dhk,'ei',dei
         Dpcm=Ddry*1.5
         goto 200
         end if
         FDpcl = -AKOH/Dpcl**2. + & 
     &   3.*(1.-dei)*Ddry**3.*dhk*Dpcl**2./(Dpcl**3.-dei*Ddry**3.)**2. + & 
     &   (A*B/(2.*Dw))*(((Dpcl-dei**(1./3.)*Ddry)/ &
     &   (2.*Dw))**(-B-1.))


         FDpcu = -AKOH/Dpcu**2. + & 
     &   3.*(1.-dei)*Ddry**3.*dhk*Dpcu**2./(Dpcu**3.-dei*Ddry**3.)**2. + &
     &   (A*B/(2.*Dw))*(((Dpcu-dei**(1./3.)*Ddry)/                       & 
     &   (2.*Dw))**(-B-1.))


            Dpcm = (Dpcu+Dpcl)/2.


         FDpcm = -AKOH/Dpcm**2. + & 
     &   3.*(1.-dei)*Ddry**3.*dhk*Dpcm**2./(Dpcm**3.-dei*Ddry**3.)**2. + &
     &   (A*B/(2.*Dw))*(((Dpcm-dei**(1./3.)*Ddry)/                       &
     &   (2.*Dw))**(-B-1.))

            If ((FDpcl*FDpcm).Le.0) Then

               If (ABS(FDpcm).Le.10e-8) Then
                  Goto 200
               Else
                   Dpcl = Dpcl
                   Dpcu = Dpcm
                   goto 100
               End if

            Else If ((FDpcl*FDpcm).GE.0) Then

                If (ABS(FDpcm).Le.10e-8) Then
                   Goto 200
                Else
                    Dpcl = Dpcm
                    Dpcu = Dpcu
                    goto 100
                End if

            Else If ((FDpcl*FDpcm).Eq.0) Then
                    Goto 200
            End if

200   Dc = Dpcm
      
      RETURN
      END
!C *** END OF SUBROUTINE DpcFHH ***************************************
!C=======================================================================
!C=======================================================================
!C
!C *** SUBROUTINE PDFACTIV
!C *** THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING
!C     TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
!C     MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
!C     and Nenes (2004). THIS ROUTINE CALCULATES FOR A PDF OF
!C     UPDRAFT VELOCITIES.
!C
!C *** WRITTEN BY ATHANASIOS NENES
!C
!C=======================================================================
!C
      SUBROUTINE PDFACTIV (WPARC,TP,AKK,ei,A,B,ACCOM,SG,SIGW,&
     & TPARC,PPARC,NACT,SMAX,NDM)
!C
      INCLUDE 'messy_cloud_parameters_uaf.inc'
      DOUBLE PRECISION NACT, NACTI, A, B,ACCOM,&
     & TP(NSMX),AKK(NSMX),ei(NSMX),SG(NSMX),NDM(NSMX),NDMI(NSMX)
      REAL             PDF
!C
!C *** Case where updraft is very small
!C
      IF (WPARC.LE.1d-6) THEN
         SMAX  = 0d0
         NACT  = 0d0
         DO j=1,NMD
         NDM(J) = 0d0
         END DO
         RETURN
      ENDIF
!C
!C *** Single updraft case
!C
      IF (SIGW.LT.1e-10) THEN
         CALL ACTIVATE (WPARC,TP,AKK,ei,A,B,ACCOM,SG,NACT,SMAX,NDM)
!C
!C *** PDF of updrafts
!C
      ELSE
         NACT  = ZERO
         SMAX  = ZERO
         PLIMT = 1e-3                                                   ! Probability of High Updraft limit
         PROBI = SQRT(-2.0*LOG(PLIMT*SIGW*SQ2PI))
         WHI   = WPARC + SIGW*PROBI                                     ! Upper updrft limit
         WLO   = 0.05                                                   ! Low updrft limit
         SCAL  = 0.5*(WHI-WLO)                                          ! Scaling for updrafts
         DO I=1,Npgauss
            WPI  = WLO + SCAL*(1.0-XGS(i))                              ! Updraft
         CALL ACTIVATE (WPI,TP,AKK,ei,A,B,ACCOM,SG,NACTI,SMAXI,NDMI)           ! # of drops
            PDF  = (1.0/SQ2PI/SIGW)*EXP(-0.5*((WPI-WPARC)/SIGW)**2)     ! Prob. of updrafts
            NACT = NACT + WGS(i)*(PDF*NACTI)                            ! Integral for drops
            SMAX = SMAX + WGS(i)*(PDF*SMAXI)                            ! Integral for Smax
            IF (PDF.LT.PLIMT) GOTO 100
         ENDDO
 100     NACT = NACT*SCAL                                               ! Scale Integrals
         SMAX = SMAX*SCAL
        DO j=1,NMD
        NDM(J)=NDMI(J)*SCAL
        END DO
      ENDIF
!C
      RETURN
!C
!C *** END OF SUBROUTINE PDFACTIV ***************************************
!C
      END
!C=======================================================================
!C=======================================================================
!C
!C *** SUBROUTINE WACTIV
!C *** THIS SUBROUTINE CALCULATES THE UPDRAFT NECESSARY TO ACHIEVE A DROP
!C     CONCENTRATION.
!C
!C *** WRITTEN BY ATHANASIOS NENES
!C
!C=======================================================================
!C
      SUBROUTINE WACTIV (NACT, WPARC, TP,SMAX)
!C
      INCLUDE 'messy_cloud_parameters_uaf.inc'
      DOUBLE PRECISION NACT, NACT1, NACT2, NACT3
      DOUBLE PRECISION TP(NSMX),SG(NSMX),NDM(NSMX)
      DOUBLE PRECISION AKK(NSMX),ei(NSMX)
!C
!C *** INITIAL VALUES FOR BISECTION **************************************
!C
      X1   = 1e-3          ! Low value of updraft
      CALL ACTIVATE (X1, TP,AKK,ei,A,B,ACCOM,SG,NACT1,SMAX,NDM)
      Y1   = NACT1/NACT - 1d0
!C
      X2   = 20           ! High value of updraft
!     CALL ACTIVATE (X2,TP, NACT2,AKK,ei,ACCOM,A,B,SG,SMAX,NDM)
! arguments seem to be mixed up - compare to next call below
      CALL ACTIVATE (X2,TP,AKK,ei,A,B,ACCOM,SG,NACT2,SMAX,NDM)
      !Y2   = NACT2/NACT - 1d0
!C
!C *** PERFORM BISECTION ************************************************
!C
20    DO 30 I=1,MAXIT
         X3   = 0.5*(X1+X2)
         CALL ACTIVATE (X3,TP,AKK,ei,A,B,ACCOM,SG,NACT3, SMAX,NDM)
         Y3   = NACT3/NACT - 1d0
!C
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            !Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
!C
         IF (ABS(X2-X1) .LE. EPSC*X1) GOTO 40
         NITER = I
!C
30    CONTINUE
!
! *** CONVERGED ; RETURN ***********************************************
!
40    X3    = 0.5*(X1+X2)
      CALL ACTIVATE (X3,TP,AKK,ei,A,B,ACCOM,SG,NACT3, SMAX,NDM)
      Y3    = NACT3/NACT - 1d0
      WPARC = X3
!C
      RETURN
!C
!C *** END OF SUBROUTINE WACTIVE ****************************************
!C
      END
!C=======================================================================
!C=======================================================================
!C
!C *** SUBROUTINE ACTIVATE
!C *** THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING
!C     TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
!C     MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
!C     and Nenes (in preparation).
!C
!C *** WRITTEN BY ATHANASIOS NENES FOR KOHLER PARTICLES
!C *** MODIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES TO INCLUDE FHH 
!C     PARTICLES 
!C
!C=======================================================================
!C
      SUBROUTINE ACTIVATE (WPARC,TP,AKK,ei,A,B,ACCOM,SG,NDRPL,SMAX,NDM)
      INCLUDE 'messy_cloud_parameters_uaf.inc'
      DOUBLE PRECISION NDRPL, WPARCEL,A,B,ACCOM,BET2,BETA
      DOUBLE PRECISION TP(NSMX),AKK(NSMX),ei(NSMX),SG(NSMX),NDM(NSMX)
      DOUBLE PRECISION C1, C2, C3, C4, X_FHH(NSMX)
!C
!C *** Setup common block variables
!C
      PRESA = PRES/1.013d5                  ! Pressure (Pa)
      DV    = (0.211d0/PRESA)*(TEMPER/273d0)**1.94
      DV    = DV*1d-4                       ! Water vapor diffusivity in air
      DBIG  = 5.0d-6
      DLOW  = 0.207683*((ACCOM)**(-0.33048))
      DLOW  = DLOW*1d-6
!C
!C Compute an average diffusivity Dv as a function of ACCOM
!C
      COEF  = ((2*PI*AMW/(RGAS*TEMPER))**0.5)
      DV    = (DV/(DBIG-DLOW))*((DBIG-DLOW)-(2*DV/ACCOM)*COEF*&
     &        (DLOG((DBIG+(2*DV/ACCOM)*COEF)/(DLOW+(2*DV/ACCOM)*&
     &        COEF))))                      ! Non-continuum effects

      WPARCEL = WPARC
!
! *** Setup constants
!
      ALFA = GRAV*AMW*DHV/CPAIR/RGAS/TEMPER/TEMPER-GRAV*AMA/RGAS/TEMPER
      BET1 = PRES*AMA/PSAT/AMW + AMW*DHV*DHV/CPAIR/RGAS/TEMPER/TEMPER
      BET2 = RGAS*TEMPER*DENW/PSAT/DV/AMW/4d0 +&
     &       DHV*DENW/4d0/AKA/TEMPER*(DHV*AMW/RGAS/TEMPER - 1d0)
      BETA = 0.5d0*PI*BET1*DENW/BET2/ALFA/WPARC/DAIR
      CF1  = 0.5*(((1/BET2)/(ALFA*WPARC))**0.5)
      CF2  = AKOH/3d0
!
!C     DETERMINATION OF EXPONENT FOR FHH PARTICLES
!
      C1     = (D11)+(D12/A)+(D13/(A*A))+(D14/(A*A*A))+(D15/(A*A*A*A))
      C2     = (D21)+(D22/A)+(D23/(A*A))+(D24/(A*A*A))+(D25/(A*A*A*A))
      C3     = (D31)+(D32/A)+(D33/(A*A))+(D34/(A*A*A))+(D35/(A*A*A*A))
      C4     = (D41)+(D42/A)+(D43/(A*A))+(D44/(A*A*A))+(D45/(A*A*A*A))
      do J=1,NMD 
       X_FHH(J)  = (C1) + (C2/B) + (C3/(B*B)) + (C4/(B*B*B))
       X_FHH(J) = X_FHH(J)*exp(log(-1.5/X_FHH(J))*(1.-ei(J))** &
     & (0.1693*exp(-0.988*AKK(J)))) !unified FHH theory --- VAK
      end do
!
! *** INITIAL VALUES FOR BISECTION *************************************
!     
      X1   = 1.0d-5   ! Min cloud supersaturation -> 0
      CALL SINTEGRAL (X1,NDRPL,WPARCEL,TP,X_FHH,BET2,SG,&
     & SINTEG1,SINTEG2,SINTEG3,NDM)
      Y1   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X1 - 1d0
!     
      X2   = 0.1d0      ! MAX cloud supersaturation = 10%
      CALL SINTEGRAL (X2,NDRPL,WPARCEL,TP,X_FHH,BET2,SG,&
     & SINTEG1,SINTEG2,SINTEG3,NDM)
      !Y2   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X2 - 1d0
!
! *** PERFORM BISECTION ************************************************
!
20    DO 30 I=1,MAXIT
         X3   = 0.5*(X1+X2)
         CALL SINTEGRAL (X3,NDRPL,WPARCEL,TP,X_FHH,BET2,SG,&
     &   SINTEG1,SINTEG2,SINTEG3,NDM)
         Y3 = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X3 - 1d0
!
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
             !Y2    = Y3
             X2    = X3
         ELSE
             Y1    = Y3
             X1    = X3
         ENDIF
!
         IF (ABS(X2-X1) .LE. EPSC*X1) GOTO 40
             NITER = I

30    CONTINUE

! *** CONVERGED ; RETURN ***********************************************
40    X3   = 0.5*(X1+X2)
!
      CALL SINTEGRAL (X3,NDRPL,WPARCEL,TP,X_FHH,BET2,SG,&
     & SINTEG1,SINTEG2,SINTEG3,NDM)
      Y3   = (SINTEG1*CF1+SINTEG2*CF2+SINTEG3*CF1)*BETA*X3 - 1d0
      
      SMAX = X3

      RETURN
!C
!C *** END OF SUBROUTINE ACTIVATE ***************************************
!C
      END
!C=======================================================================
!C=======================================================================
!C
!C *** SUBROUTINE SINTEGRAL
!C *** THIS SUBROUTINE CALCULATES THE CONDENSATION INTEGRALS, ACCORDING
!C     TO THE POPULATION SPLITTING ALGORITHM AND THE SUBSEQUENT VERSIONS:
!C
!C       - Nenes and Seinfeld (2003)       Population Splitting
!C       - Fountoukis and Nenes (2004)     Modal formulation
!C       - Barahona and Nenes (2010)       Approach for large CCN
!C       - Morales and Nenes (2014)        Population Splitting revised
!C
!C *** WRITTEN BY ATHANASIOS NENES for Kohler Particles
!C *** MODFIFIED BY PRASHANT KUMAR AND ATHANASIOS NENES TO INCLUDE FHH
!C     PARTICLES
!C=======================================================================
!C
      SUBROUTINE SINTEGRAL (SPAR, SUMMA, WPARCEL, TP, XFHH, BET2, SG,&
     &                      SUM, SUMMAT, SUMFHH,NDM)
!C
      INCLUDE 'messy_cloud_parameters_uaf.inc'
      DOUBLE PRECISION SUM, SUMMAT, SUMMA, Nd(NSMX),WPARCEL,TP(NSMX),&
     &                 INTEG1(NSMX),INTEG2(NSMX),SG(NSMX),BET2&
     &                 ,SUMFHH,INTEG1F(NSMX),NdF(NSMX),XFHH(NSMX)&
     &                 ,NDM(NSMX)

      REAL             ERF1,ERF2,ERF3,ERF4,ERF5,ERF6,ERF4F,ERF5F,ERF66F
      REAL             ORISM1, ORISM2, ORISM3, ORISM4, ORISM5,ORISM6
      REAL             intaux1p1, intaux1p2, DLGSP1,DLGSP2
      REAL             scrit
!C
      REAL             ORISM1F, ORISM2F, ORISM3F, ORISM4F, ORISM5F,&
     &                 ORISM6F, ORISM7F, ORISM8F, ORISM9F, &
     &                 ORISM66F
!      
      SQTWO  = SQRT(2d0)
!C
!C ** Population Splitting -- Modified by Ricardo Morales 2013

      DESCR  = 1d0 - (16d0/9d0)*ALFA*WPARCEL*BET2*(AKOH/SPAR**2)**2
      IF (DESCR.LE.0d0) THEN
         CRIT2  = .TRUE.             
         scrit  = ((16d0/9d0)*ALFA*WPARCEL*BET2*(AKOH**2))**(0.25d0)    ! Scrit - (only for DELTA < 0 )
         RATIO  = (2.0d7/3.0)*AKOH*(SPAR**(-0.3824)-scrit**(-0.3824))   ! Computing sp1 and sp2 (sp1 = sp2)
         RATIO  = 1/SQTWO + RATIO
         IF (RATIO.GT.1.0) RATIO = 1.0
         SSPLT2 = SPAR*RATIO
      ELSE
         CRIT2  = .FALSE.
         SSPLT1 = 0.5d0*(1d0-SQRT(DESCR))                               ! min root --> sp1
         SSPLT2 = 0.5d0*(1d0+SQRT(DESCR))                               ! max root --> sp2
         SSPLT1 = SQRT(SSPLT1)*SPAR                                     ! Multiply ratios with Smax
         SSPLT2 = SQRT(SSPLT2)*SPAR
      ENDIF
!C
      SSPLT = SSPLT2  ! Store Ssplit in COMMON
!C
!C *** Computing the condensation integrals I1 and I2
!C
      SUM       = 0.0d0   !Contribution of integral 1 for Kohler 
      SUMMAT    = 0.0d0   !Contribution of integral 2 for kohler
      SUMMA     = 0.0d0   !Variable that stores all droplets
      SUMFHH    = 0.0d0   !Contribution of FHH integral
!C
      DO J = 1, NMD
!C
      IF (MODE(J).EQ.1) THEN          ! Kohler modes
!C
        DLGSG  = DLOG(SIG(J))                            !ln(sigmai)
        DLGSP  = DLOG(SG(J)/SPAR)                        !ln(sg/smax)
        DLGSP2 = DLOG(SG(J)/SSPLT2)                      !ln(sg/sp2)
!C
        ORISM1 = 2.d0*DLGSP2/(3.d0*SQTWO*DLGSG)          ! u(sp2)
        ORISM2 = ORISM1 - 3.d0*DLGSG/(2.d0*SQTWO)        ! u(sp2)-3ln(sigmai)/(2sqrt(2)
        ORISM5 = 2.d0*DLGSP/(3.d0*SQTWO*DLGSG)           ! u(smax)
        ORISM3 = ORISM5 - 3.d0*DLGSG/(2.d0*SQTWO)        ! u(smax)-3ln(sigmai)/(2sqrt(2)
        DEQ    = AKOH*2d0/SG(j)/3d0/SQRT(3d0)            ! Dp0 = Dpc/sqrt(3) - Equilibrium diameter

        ERF2   = erfp(ORISM2)
        ERF3   = erfp(ORISM3)
 
        INTEG2(J) = (EXP(9D0/8D0*DLGSG*DLGSG)*TP(J)/SG(J))*&
     &              (ERF2 - ERF3)                          ! I2(sp2,smax)

        IF (CRIT2) THEN     

          ORISM6 = (SQTWO*DLGSP2/3d0/DLGSG)-(1.5d0*DLGSG/SQTWO)
          ERF6   = erfp(ORISM6)

          INTEG1(J) = 0.0d0
          DW3       = TP(j)*DEQ*EXP(9D0/8D0*DLGSG*DLGSG)*&   ! 'inertially' limited particles
     &           (1d0-ERF6)*((BET2*ALFA*WPARCEL)**0.5d0)

        ELSE
 
          EKTH    = EXP(9D0/2d0*DLGSG*DLGSG)
          DLGSP1  = DLOG(SG(J)/SSPLT1)                      ! ln(sg/sp1)
          ORISM4  = ORISM1 + 3.d0*DLGSG/SQTWO               ! u(sp2) + 3ln(sigmai)/sqrt(2)
          ERF1    = erfp(ORISM1)
          ERF4    = erfp(ORISM4)

          intaux1p2 =  TP(J)*SPAR*((1-ERF1) -&
     &              0.5d0*((SG(J)/SPAR)**2)*EKTH*(1-ERF4))  ! I1(0,sp2)

          ORISM1  = 2.d0*DLGSP1/(3.d0*SQTWO*DLGSG)          ! u(sp1)
          ORISM4  = ORISM1 + 3.d0*DLGSG/SQTWO               ! u(sp1) + 3ln(sigmai)/sqrt(2)
          ORISM6  = (SQTWO*DLGSP1/3d0/DLGSG)-(1.5d0*DLGSG/SQTWO)

          ERF1 = erfp(ORISM1)
          ERF4 = erfp(ORISM4)
          ERF6 = erfp(ORISM6)

          intaux1p1 = TP(J)*SPAR*((1-ERF1) -&
     &              0.5d0*((SG(J)/SPAR)**2)*EKTH*(1-ERF4))    ! I1(0,sp1)

          INTEG1(J) = (intaux1p2-intaux1p1)                   ! I1(sp1,sp2) = I1(0,sp2) - I1(0,sp1)
!
          DW3 = TP(j)*DEQ*EXP(9D0/8D0*DLGSG*DLGSG)*&          ! 'inertially' limited particles.
     &       (1d0-ERF6)*((BET2*ALFA*WPARCEL)**0.5d0)
 
        ENDIF

!C *** Calculate number of Drops

        ERF5     = erfp(ORISM5)
! 
        Nd(J)    = (TP(J)/2.0)*(1.0-ERF5)
        SUM      = SUM    + INTEG1(J) + DW3           !SUM OF INTEGRAL 1 FOR KOHLER
        SUMMAT   = SUMMAT + INTEG2(J)                 !SUM OF INTEGRAL 2 FOR KOHLER
        SUMMA    = SUMMA  + Nd(J)                     !SUM OF ACTIVATED KOHLER PARTICLES
        NDM(J)   = Nd(J)
!C
      ELSEIF (MODE(J).EQ.2) THEN                      ! FHH modes
!C      
        DLGSGF  = DLOG(SIG(J))                        ! ln(sigma,i)
        DLGSPF  = DLOG(SG(J)/SPAR)                    ! ln(sg/smax)
        ORISM1F = (SG(J)*SG(J))/(SPAR*SPAR)           ! (sg/smax)^2
        ORISM2F = EXP(2D0*XFHH(J)*XFHH(J)*DLGSGF*DLGSGF)    ! exp(term)
        ORISM3F = SQTWO*XFHH(J)*DLGSGF                   ! sqrt(2).x.ln(sigma,i)
        ORISM4F = DLGSPF/(-1*ORISM3F)                 ! Umax
        ORISM5F = ORISM3F - ORISM4F
        ERF5F   = erfp(ORISM5F)
        ORISM6F = ERF5F
        ORISM7F = ORISM6F + 1
        ORISM8F = 0.5*ORISM1F*ORISM2F*ORISM7F
        ERF4F   = erfp(ORISM4F)
        ORISM9F = ORISM8F + ERF4F - 1

        INTEG1F(J) =-1*TP(J)*SPAR*ORISM9F

!
!      ! Correction for kinetically limited large particles (D. Barahona et al. ACPD 2009)
!
      DLGSPF     = DLOG(SG(j)/SSPLT)
      ORISM66F  = (-DLGSPF/XFHH(J)/SQTWO/DLGSGF)-(-XFHH(J)*DLGSGF/SQTWO)
         DEQF=Dpc(j)
      ERF66F=erfp(ORISM66F)
      DW3F=TP(j)*DEQF*EXP(XFHH(J)**2*DLGSGF*DLGSGF)* &
     &    (1d0-ERF66F)*((BET2*ALFA*WPARCEL)**0.5d0)
!C
!C *** Calculate number of drops activated by FHH theory
!C
        ERF4F   = erfp(ORISM4F)

        NdF(J)  = (TP(J)/2.0)*(1-ERF4F)
        SUMFHH  = SUMFHH + INTEG1F(J) + DW3F         !Sum of Integral 1 for FHH + GCCN correction
        SUMMA   = SUMMA + NdF(J)              !Sum of ACTIVATED Kohler + FHH particles
        NDM(J)  = NdF(J)
      ENDIF

      ENDDO
      RETURN
!C
      END
!C=======================================================================
!C=======================================================================
!C
!C *** SUBROUTINE PROPS
!C *** THIS SUBROUTINE CALCULATES THE THERMOPHYSICAL PROPERTIES
!C
!C *** WRITTEN BY ATHANASIOS NENES
!C
!C=======================================================================
!C
      SUBROUTINE PROPS
      INCLUDE 'messy_cloud_parameters_uaf.inc'
      REAL  VPRES, SFT
!C
      !PRESA = PRES/1.013d5                  ! Pressure (Pa)
      DAIR  = PRES*AMA/RGAS/TEMPER          ! Air density
      AKA   = (4.39+0.071*TEMPER)*1d-3      ! Air thermal conductivity
      PSAT  = VPRES(SNGL(TEMPER))*(1e5/1.0d3)! Saturation vapor pressure
      SURT  = SFT(SNGL(TEMPER))             ! Surface Tension for water (J m-2)
!C
      RETURN
!C
!C *** END OF SUBROUTINE PROPS ******************************************
!C
      END
!C=======================================================================
!C=======================================================================
!C
!C *** FUNCTION VPRES
!C *** THIS FUNCTION CALCULATES SATURATED WATER VAPOUR PRESSURE AS A
!C     FUNCTION OF TEMPERATURE. VALID FOR TEMPERATURES BETWEEN -50 AND
!C     50 C.
!C
!C========================= ARGUMENTS / USAGE ===========================
!C
!C  INPUT:
!C     [T]
!C     REAL variable.
!C     Ambient temperature expressed in Kelvin.
!C  OUTPUT:
!C     [VPRES]
!C     REAL variable.
!C     Saturated vapor pressure expressed in mbar.
!C
!C=======================================================================
!C
      REAL FUNCTION VPRES (T)
      REAL A(0:6), T
      DATA A/6.107799610E+0, 4.436518521E-1, 1.428945805E-2,&
     &       2.650648471E-4, 3.031240396E-6, 2.034080948E-8,&
     &       6.136820929E-11/

      TTEMP = T-273
      VPRES = A(6)*TTEMP
      DO I=5,1,-1
         VPRES = (VPRES + A(I))*TTEMP
      ENDDO
      VPRES = VPRES + A(0)
      RETURN
      END
!C=======================================================================
!C=======================================================================
!C
!C *** FUNCTION SFT
!C *** THIS FUNCTION CALCULATES WATER SURFACE TENSION AS A
!C     FUNCTION OF TEMPERATURE. VALID FOR TEMPERATURES BETWEEN -40 AND
!C     40 C.
!C
!C ======================== ARGUMENTS / USAGE ===========================
!C
!C  INPUT:
!C     [T]
!C     REAL variable.
!C     Ambient temperature expressed in Kelvin.
!C
!C  OUTPUT:
!C     [SFT]
!C     REAL variable.
!C     Surface Tension expressed in J m-2.
!C
!C=======================================================================
!C
      REAL FUNCTION SFT (T)
      REAL T
!C
      TPARS = T-273
      SFT   = 0.0761-1.55e-4*TPARS
!C
      RETURN
      END
!C=======================================================================
!C ***********************************************************************
!C
      SUBROUTINE GAULEG (X,W,N)
!C
!C Calculation of points and weights for N point GAUSS integration
!C ***********************************************************************
      DIMENSION X(N), W(N)
      PARAMETER (EPSC=1.E-6)
      PARAMETER (X1=-1.0, X2=1.0)
!C
!C Calculation
!C
      M=(N+1)/2
      XM=0.5d0*(X2+X1)
      XL=0.5d0*(X2-X1)
      DO 12 I=1,M
        Z=COS(3.141592654d0*(I-.25d0)/(N+.5d0))
1       CONTINUE
          P1=1.d0
          P2=0.d0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.d0*J-1.)*Z*P2-(J-1.d0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.d0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPSC)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.d0*XL/((1.d0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END

!C=======================================================================
!C
!C *** REAL FUNCTION erfp
!C *** THIS SUBROUTINE CALCULATES THE ERROR FUNCTION USING A
!C *** POLYNOMIAL APPROXIMATION
!C
!C=======================================================================
!C
      REAL*8 FUNCTION erfp(x)
        REAL :: x
        REAL*8 :: AA(4), axx, y
        DATA AA /0.278393d0,0.230389d0,0.000972d0,0.078108d0/
        
        y = dabs(dble(x))
        axx = 1.d0 + y*(AA(1)+y*(AA(2)+y*(AA(3)+y*AA(4))))
        axx = axx*axx
        axx = axx*axx
        axx = 1.d0 - (1.d0/axx)
        if(x.le.0.) then
          erfp = -axx
        else
          erfp = axx
        endif
      RETURN
      END FUNCTION
!C
!C=======================================================================
!C
!C *** REAL FUNCTION ERF
!C *** THIS SUBROUTINE CALCULATES THE ERROR FUNCTION
!C
!C *** OBTAINED FROM NUMERICAL RECIPIES
!C
!C=======================================================================
!C
      SUBROUTINE CALCERF(X,ERF)
      IF(X.LT.0.)THEN
        ERF=-GAMMP(.5,X**2)
      ELSE
        ERF=GAMMP(.5,X**2)
      ENDIF
      RETURN
      END
!C
!C=======================================================================
!C
      FUNCTION GAMMLN_K(XX)
!C
!C=======================================================================
!C
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,&
     &    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN_K=TMP+LOG(STP*SER)
      RETURN
      END
!C
!C=======================================================================
!C

      FUNCTION GAMMP(A,X)
!C
!C=======================================================================
!C
      IF(X.LT.0..OR.A.LE.0.)PAUSE
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END


!C
!C=======================================================================
!C

      SUBROUTINE GCF(GAMMCF,A,X,GLN)
!C
!C=======================================================================
!C
      PARAMETER (ITMAX=100,EPSC=3.E-7)
      GLN=GAMMLN_K(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPSC)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*LOG(X)-GLN)*G
      RETURN
      END


!C
!C=======================================================================
!C
      SUBROUTINE GSER(GAMSER,A,X,GLN)
!C
!C=======================================================================
!C
      PARAMETER (ITMAX=100,EPSC=3.E-7)
      GLN=GAMMLN_K(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)PAUSE
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPSC)GO TO 1
11    CONTINUE
      PAUSE 'A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END

!C=======================================================================
!C=======================================================================

