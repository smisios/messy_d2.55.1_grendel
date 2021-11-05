C=======================================================================
C
C *** ISORROPIA CODE
C *** PROGRAM MAIN
C *** THIS IS THE MAIN ROUTINE OF THE PROGRAM. IT READS INPUT FROM THE
C     USER, THEN CALLS SUBROUTINE ISOROPIA
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      PROGRAM MAIN
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INCLUDE 'main.inc'
      INTEGER   ERRSTKI
      LOGICAL   DEXS,       IEXS,       EOF,     STKOFL
      CHARACTER ERRMSGI*40, VERSION*15, SCASE*15
      DOUBLE PRECISION LMASS, LMOLE, GAS, MOLAL, IONIC
      DIMENSION AERLIQ(15),  AERSLD(19),  CNTRL(2), OTHER(9), MOLAL(10),
     &          ERRSTKI(25), ERRMSGI(25), GASAQ(3), W(8),     GAS(3)
C
C *** OBTAIN INFORMATION ABOUT ISORROPIA *******************************
C
      CALL ISORINF (VERSION, NCMP, NIONS, NGASAQ, NSOL, NERR, TINY, GRT)
C
C *** PRINT LOGO *******************************************************
C
      WRITE(*,*)
      WRITE(*,*) '-----------------------------------------------------'
      WRITE(*,*) ' ISORROPIA v',VERSION
      WRITE(*,*) ' Copyright 1996-2008'
      WRITE(*,*) ' University of Miami, Carnegie Mellon University'
      WRITE(*,*) ' Georgia Institute of Technology'
      WRITE(*,*) ' Written by Athanasios Nenes (nenes@eas.gatech.edu)'
      WRITE(*,*) ' and Christos Fountoukis',
     &           ' (Christos.Fountoukis@chbe.gatech.edu)'
      WRITE(*,*) '-----------------------------------------------------'
      WRITE(*,*)
C
C *** INITIALIZE VARIABLES *********************************************
C 
      CNFFIL = 'isrpia.cnf'     ! Name of configuration file
      CALL CHRBLN (CNFFIL, ICE) ! Position of last non-blank character
C
C *** INPUT CONFIGURATION IF CNF FILE EXISTS *****************************
C
      CALL INPCNF (IERR)        ! Get configuration from file
C
      IF (IERR.EQ.0) THEN
         WRITE(*,3000) 'Parameters read from file [',CNFFIL(1:ICE),']'
      ELSE
         WRITE(*,3000) 'Configuration file [',CNFFIL(1:ICE),'] ',
     &                 'not found, using defaults'
      ENDIF
C
C *** OBTAIN RUNFILE NAME - SETUP I/O FILES ******************************
C
6     WRITE (*,*) ' '
      WRITE (*,*) 'File name with runs [Enter=screen input]: '
      READ  (*,'(A)') INPFIL
C      INPFIL='nonurban2.INP'
      CALL CHRBLN(INPFIL,IPE)
      IF (IPE.EQ.1 .AND. INPFIL(1:1).EQ.' ') THEN
         INPFIL = 'SCREEN'
      ENDIF
      CALL CHRBLN(INPFIL,IPE)
      IPD = INDEX(INPFIL,'.')
      IF (IPD.EQ.0) IPD = IPE+1
      IF (INPFIL(1:6).EQ.'SCREEN' .OR. INPFIL(1:6).EQ.'screen') THEN
         IFIL = .FALSE.
      ELSE
         IFIL = .TRUE.        
      ENDIF
      OUTFIL = INPFIL(1:IPD-1)
      DATFIL = OUTFIL
      CALL APPENDEXT (OUTFIL, '.txt', .TRUE.)
      CALL APPENDEXT (DATFIL, '.dat', .TRUE.)
C
C *** INITIALIZE INPUT; DETRMINE WHERE INPUT IS COMING FROM ************
C
      INQUIRE (FILE=INPFIL, EXIST=IEXS)   ! Does INP file exist?
      IF (IFIL) THEN
         IF (IEXS) THEN                   ! INP files exists
            II = 14
            OPEN (UNIT=II, FILE=INPFIL, STATUS='OLD')
4           READ (II,*,ERR=4) INUNIT
5           READ (II,*,ERR=5) IPROBI, METSTBLI
         ELSE                             ! INP file doesnt exist
            IFIL = .FALSE.
            WRITE (*,3000) 'Input file [',INPFIL(1:IPE),'] doesn''t ',
     &                     'exist; try again.'
            GOTO 6
         ENDIF
      ENDIF
C
C *** GENERATE INITIAL OUTPUT TO FILES *********************************
C
      IF (OFIL) THEN                       ! REPORT FILE
         IO = 11
         OPEN (IO, FILE=OUTFIL, STATUS='UNKNOWN')
         CALL PUSHEND(IO)
      ELSE
         IO = 0   
      ENDIF
C
      IF (DFIL) THEN                       ! DATA FILE
         ID = 12
         INQUIRE (FILE=DATFIL, EXIST=DEXS)
         OPEN (UNIT=ID, FILE=DATFIL, STATUS='UNKNOWN', RECL=4096)
         CALL PUSHEND (ID)
         IF (.NOT.DEXS) WRITE (ID,3100)
      ENDIF
C
      IUSR = 99                            ! USER OUTPUT
      CALL USROUT(1)
C
C *** INPUT DATA (LOOP) ************************************************
C
10    CALL INPDAT (EOF)
      IF (EOF) GOTO 50    ! EOF encountered
C
C *** WRITE INPUT TO REPORT FILE ***************************************
C
      WRITE (IO,3000) '*** [ INPUT ] ***'
      WRITE (IO,3000) ' '
      IF (IPROBI.EQ.0) THEN
         WRITE (IO,1500) 'TOTAL:',  '(ug/m3 air)','(umol/m3 air)'
      ELSE
         WRITE (IO,1500) 'AEROSOL:','(ug/m3 air)','(umol/m3 air)'
      ENDIF
      WRITE (IO,2000) '[Na   ] ', WI(1)*23.0*1.D6, WI(1)*1.D6
      WRITE (IO,2000) '[H2SO4] ', WI(2)*98.0*1.D6, WI(2)*1.D6
      WRITE (IO,2000) '[NH3  ] ', WI(3)*17.0*1.D6, WI(3)*1.D6
      WRITE (IO,2000) '[HNO3 ] ', WI(4)*63.0*1.D6, WI(4)*1.D6
      WRITE (IO,2000) '[HCL  ] ', WI(5)*36.5*1.D6, WI(5)*1.D6
      WRITE (IO,2000) '[Ca   ] ', WI(6)*40.1*1.D6, WI(6)*1.D6
      WRITE (IO,2000) '[K    ] ', WI(7)*39.1*1.D6, WI(7)*1.D6
      WRITE (IO,2000) '[Mg   ] ', WI(8)*24.3*1.D6, WI(8)*1.D6
      WRITE (IO,3000) ' '
      WRITE (IO,2800) 'TEMPERATURE (K) [', TEMPI,     ']'
      WRITE (IO,2800) 'REL. HUMIDITY   [', RHI*100.0, '] %'
      WRITE (IO,3000) ' '
C
C *** CALL ISORROPIA **************************************************
C
      CNTRL(1) = IPROBI     ! 0=FORWARD PROBLEM, 1=REVERSE PROBLEM
      CNTRL(2) = METSTBLI   ! 0=SOLID+LIQUID AEROSOL, 1=METASTABLE
      print*, "CNTRL", CNTRL
C
      CALL ISOROPIA (WI, RHI, TEMPI,  CNTRL, 
     &               W,  GAS, AERLIQ, AERSLD, SCASE, OTHER)
C
      CALL ISERRINF (ERRSTKI, ERRMSGI, NOFER, STKOFL) ! Get error status
C
      IF (NOFER.GT.0) THEN
         IF(ERRSTKI(NOFER).GT.1000) GOTO 40           ! FATAL ERROR, ABORT.
      ENDIF
C
C *** SAVE RESULTS TO WORK ARRAYS (units = mole/m3, kg/m3 for water) ****
C
      GNH3 = GAS(1)                  ! Gaseous aerosol species
      GHNO3= GAS(2)
      GHCL = GAS(3)
C
      DO 11 I=1,7                ! Liquid aerosol species
         MOLAL(I) = AERLIQ(I)
  11  CONTINUE
      MOLAL(8) = AERLIQ(13)
      MOLAL(9) = AERLIQ(14)
      MOLAL(10)= AERLIQ(15)
      WATER = 18.0D-3*AERLIQ(7+1)
      DO 12 I=1,NGASAQ
         GASAQ(I) = AERLIQ(7+1+I)
  12  CONTINUE
      COH = AERLIQ(12)
C
      CNANO3  = AERSLD(01)            ! Solid aerosol species
      CNH4NO3 = AERSLD(02) 
      CNACL   = AERSLD(03)
      CNH4CL  = AERSLD(04)
      CNA2SO4 = AERSLD(05)
      CNH42S4 = AERSLD(06)
      CNAHSO4 = AERSLD(07)
      CNH4HS4 = AERSLD(08)
      CLC     = AERSLD(09)
      CCASO4  = AERSLD(10)
      CCANO32 = AERSLD(11)
      CCACL2  = AERSLD(12)
      CK2SO4  = AERSLD(13)
      CKHSO4  = AERSLD(14)
      CKNO3   = AERSLD(15)
      CKCL    = AERSLD(16)
      CMGSO4  = AERSLD(17)
      CMGNO32 = AERSLD(18)
      CMGCL2  = AERSLD(19)
C
      SULRAT  = OTHER(2)
      SULRATW = OTHER(3)
      SODRAT  = OTHER(4)
      IONIC   = OTHER(5)
      SO4RAT  = OTHER(7)
      CRNARAT = OTHER(8)
      CRRAT   = OTHER(9)
C
C *** PERFORM OTHER TYPES OF CALCULATIONS *****************************
C
      GMASS = GNH3*17.0 + GHNO3*63.0 + GHCL*36.5
      GMOLE = GNH3 + GHNO3 + GHCL
C
      LMASS = WATER*1.D3    + MOLAL(1)      + MOLAL(2)*23.0 + 
     &        MOLAL(3)*18.0 + MOLAL(4)*35.5 + MOLAL(5)*96.0 + 
     &        MOLAL(6)*97.0 + MOLAL(7)*62.0 + MOLAL(8)*40.1 +
     &        MOLAL(9)*39.1 + MOLAL(10)*24.3+ GASAQ(1)*17.0 +
     &        GASAQ(2)*36.5 + GASAQ(3)*63.0
      LMOLE = WATER*1.D3/18.0 + MOLAL(1) + MOLAL(2) + MOLAL(3) + 
     &        MOLAL(4)        + MOLAL(5) + MOLAL(6) + MOLAL(7) +
     &        MOLAL(8)        + MOLAL(9) + MOLAL(10)+ GASAQ(1) +
     &        GASAQ(2)        + GASAQ(3)
C
      SMASS = CNH42S4*132.0+CNH4HS4*115.0+CLC   *247.0+CNH4NO3*80.0+
     &        CNH4CL * 53.5+CNA2SO4*142.0+CNANO3* 85.0+CNACL  *58.5+
     &        CNAHSO4*120.0+CCASO4*136.1 +CCANO32*164.0+CCACL2*111.0+
     &        CK2SO4*174.2 +CKHSO4*136.1 +CKNO3*101.1 +CKCL*74.5   +
     &        CMGSO4*120.3 +CMGNO32*148.3 +CMGCL2*95.2

      SMOLE = CNH42S4+CNH4HS4+CLC+CNH4NO3+CNH4CL+CNA2SO4+CNANO3+CNACL+
     &        CNAHSO4+CCASO4+CCANO32+CCACL2+CK2SO4+CKHSO4+CKNO3+
     &        CKCL+CMGSO4+CMGNO32+CMGCL2
C
C *** WRITE DATA POINT IN WORKSHEET FILE ********
C
C UGR/M3 FOR ALL SPECIES
C
      IF (DFIL) THEN
         WRITE (ID,3200) W(1)*23.D6, W(2)*98.D6,  W(3)*17.D6, 
     &                W(4)*63.D6, W(5)*36.5D6, W(6)*40.1D6,
     &                W(7)*39.1D6, W(8)*24.3D6, RHI, TEMPI,
     &                GNH3 *17.D6, GHCL *36.5D6, GHNO3*63.D6, 
     &                CNACL*58.5D6, 
     &                CNANO3*85.D6, CNA2SO4*142.D6, CNAHSO4*120.D6,
     &                CNH4CL*53.5D6, CNH4NO3*80.D6, CNH42S4*132.D6, 
     &                CNH4HS4*115.D6, CLC*247.D6,    CCASO4*136.1D6,
     &                CCANO32*164.0D6,CCACL2*111.0D6,CK2SO4*174.2D6,
     &                CKHSO4*136.1D6,CKNO3*101.1D6, CKCL*74.5D6,
     &                CMGSO4*120.3D6,CMGNO32*148.3D6,CMGCL2*95.2D6,
     &                MOLAL(1)*1.D6, MOLAL(2)*23.D6, MOLAL(3)*18.D6, 
     &                MOLAL(4)*35.5D6, MOLAL(5)*96.D6,
     &                MOLAL(6)*97.D6, MOLAL(7)*62.D6, MOLAL(8)*40.1D6,
     &                MOLAL(9)*39.1D6, MOLAL(10)*24.3D6,
     &                (W(3)-GNH3)*17.D6,
     &                (W(5)-GHCL)*36.5D6,(W(4)-GHNO3)*63.D6,
     &                WATER*1.D9,LMASS*1.D6, SMASS*1.D6,' "',SCASE,'"'
      ENDIF
C
C *** WRITE REPORT  *****************************
C
      CALL CHRBLN (SCASE, IEND)
      WRITE (IO,3000) '*** [ SOLUTION ] ***'
      WRITE (IO,3000) ' '
      WRITE (IO,3000) 'GENERAL INFO:'
      IF (IPROBI.EQ.0) THEN
         WRITE (IO,3000) 'PROBLEM TYPE   [FOREWARD]'
      ELSE
         WRITE (IO,3000) 'PROBLEM TYPE   [REVERSE ]'
      ENDIF
      IF (METSTBLI.EQ.1) THEN
         WRITE (IO,3000) 'AEROSOL STATE  [LIQUID ONLY (METASTABLE)]'
      ELSE 
         WRITE (IO,3000) 'AEROSOL STATE  [SOLID + LIQUID POSSIBLE ]'
      ENDIF
      IF (IACALCI.EQ.0) THEN
         WRITE (IO,3000) 'ACTIVITY COEFS [KUSSIK-MEISSNER; FULL ',
     &                   'CALCULATIONS]'
      ELSE
         WRITE (IO,3000) 'ACTIVITY COEFS [KUSSIK-MEISSNER; TABULATED ',
     &                   'COEFFICIENTS]'
      ENDIF
      IF (NADJI.EQ.1) THEN
         WRITE (IO,3000) 'MASS CONSERV.  [FORCED TO MACHINE PRECISION]'
      ELSE
         WRITE (IO,3000) 'MASS CONSERV.  [NORMAL SOLUTION PROCEDURE]'
      ENDIF
      WRITE (IO,3000) 'REGIME CASE          [',SCASE(1:IEND),']'
      WRITE (IO,2500) 'SO4 RATIO            [',SULRAT,  ']'
      WRITE (IO,2500) 'Na  RATIO            [',SODRAT,  ']'
      WRITE (IO,2500) 'CRUSTAL RATIO        [',CRRAT, ']'
      WRITE (IO,2500) 'CRUSTAL & Na RATIO   [',CRNARAT,']'
      WRITE (IO,2500) 'SO4 RATIO W/CRUSTALS [',SO4RAT, ']'
      WRITE (IO,2500) 'SO4 POOR RATIO       [', SULRATW, ']'
      WRITE (IO,3000) ' '
      WRITE (IO,3000) 'COMPOSITION:'
      WRITE (IO,1500) 'GAS+AEROSOL:', '(ug/m3 air)', '(umol/m3 air)'
      WRITE (IO,2000) '[Na   ]', W(1)*23.0D6, W(1)*1.0D6
      WRITE (IO,2000) '[H2SO4]', W(2)*98.0D6, W(2)*1.0D6
      WRITE (IO,2000) '[NH3  ]', W(3)*17.0D6, W(3)*1.0D6
      WRITE (IO,2000) '[HNO3 ]', W(4)*63.0D6, W(4)*1.0D6
      WRITE (IO,2000) '[HCL  ]', W(5)*36.5D6, W(5)*1.0D6
      WRITE (IO,2000) '[Ca   ]', W(6)*40.1D6, W(6)*1.0D6
      WRITE (IO,2000) '[K    ]', W(7)*39.1D6, W(7)*1.0D6
      WRITE (IO,2000) '[Mg   ]', W(8)*24.3D6, W(8)*1.0D6
      WRITE (IO,3000) ' '
C
C *** GAS PHASE ************************
C
      IF (GNH3+GHNO3+GHCL.GT.0.0) THEN
         WRITE (IO,1500) 'GAS:','(ug/m3 air)','(umol/m3 air)','(% mass)'
     &                  ,'(% mole)'
         WRITE (IO,2000) '[NH3  ]', GNH3 *17.0*1.D6, GNH3 *1.D6, 
     &                    GNH3 *17.0/GMASS*1.D2, GNH3 /GMOLE*1.D2
         WRITE (IO,2000) '[HNO3 ]', GHNO3*63.0*1.D6, GHNO3*1.D6, 
     &                    GHNO3*63.0/GMASS*1.D2, GHNO3/GMOLE*1.D2
         WRITE (IO,2000) '[HCL  ]', GHCL *36.5*1.D6, GHCL *1.D6, 
     &                    GHCL *36.5/GMASS*1.D2, GHCL /GMOLE*1.D2
      ELSE
         WRITE (IO,3000) 'NOTHING REMAINS IN GAS PHASE'
      ENDIF
      WRITE (IO,3000) ' '
C
C *** AEROSOL PHASE (GENERAL) **********
C
      WRITE (IO,1500) 'AEROSOL MASS:', '(ug/m3 air)', '(% of total)'
      WRITE (IO,2700) '[TOTAL ]', (SMASS+LMASS)*1.D6
      IF (WATER.GT.TINY) WRITE (IO,2700) '[LIQUID]', LMASS*1.D6, 
     &                                    LMASS/(LMASS+SMASS)*1.D2
      IF (SMASS.GT.TINY) WRITE (IO,2700) '[SOLID ]', SMASS*1.D6, 
     &                                    SMASS/(LMASS+SMASS)*1.D2
      WRITE (IO,3000) ' '
C
C *** SOLID AEROSOL PHASE **************
C
      IF (SMASS .GT. TINY) THEN
         WRITE (IO,1500) 'SOLID AEROSOL:','(ug/m3 air)','(umol/m3 air)',
     &                   '(% mass)','(% mole)'
         WRITE (IO,2000) '[NaNO3    ]', CNANO3 * 85.0*1.D6,CNANO3 *1.D6,
     &                    CNANO3 * 85.0/SMASS*1.D2, CNANO3 /SMOLE*1.D2
         WRITE (IO,2000) '[Na2SO4   ]', CNA2SO4*142.0*1.D6,CNA2SO4*1.D6,
     &                    CNA2SO4*142.0/SMASS*1.D2, CNA2SO4/SMOLE*1.D2
         WRITE (IO,2000) '[NaHSO4   ]', CNAHSO4*120.0*1.D6,CNAHSO4*1.D6,
     &                    CNAHSO4*120.0/SMASS*1.D2, CNAHSO4/SMOLE*1.D2
         WRITE (IO,2000) '[NaCL     ]', CNACL  * 58.5*1.D6,CNACL  *1.D6,
     &                    CNACL  * 58.5/SMASS*1.D2, CNACL  /SMOLE*1.D2
         WRITE (IO,2000) '[NH4CL    ]', CNH4CL * 53.5*1.D6,CNH4CL *1.D6,
     &                    CNH4CL * 53.5/SMASS*1.D2, CNH4CL /SMOLE*1.D2
         WRITE (IO,2000) '[NH4NO3   ]', CNH4NO3* 80.0*1.D6,CNH4NO3*1.D6,
     &                    CNH4NO3* 80.0/SMASS*1.D2, CNH4NO3/SMOLE*1.D2
         WRITE (IO,2000) '[(NH4)2SO4]', CNH42S4*132.0*1.D6,CNH42S4*1.D6,
     &                    CNH42S4*132.0/SMASS*1.D2, CNH42S4/SMOLE*1.D2
         WRITE (IO,2000) '[NH4HSO4  ]', CNH4HS4*115.0*1.D6,CNH4HS4*1.D6,
     &                    CNH4HS4*115.0/SMASS*1.D2, CNH4HS4/SMOLE*1.D2
         WRITE (IO,2000) '[LC       ]', CLC    *247.0*1.D6,CLC    *1.D6,
     &                    CLC    *247.0/SMASS*1.D2, CLC    /SMOLE*1.D2
         WRITE (IO,2000) '[CaSO4    ]', CCASO4 *136.1*1.D6,CCASO4 *1.D6,
     &                    CCASO4 *136.1/SMASS*1.D2, CCASO4 /SMOLE*1.D2
         WRITE (IO,2000) '[Ca(NO3)2 ]', CCANO32*164.0*1.D6,CCANO32*1.D6,
     &                    CCANO32*164.0/SMASS*1.D2, CCANO32 /SMOLE*1.D2
         WRITE (IO,2000) '[CaCl2    ]', CCACL2 *111.0*1.D6,CCACL2 *1.D6,
     &                    CCACL2 *111.0/SMASS*1.D2, CCACL2 /SMOLE*1.D2
         WRITE (IO,2000) '[K2SO4    ]', CK2SO4 *174.2*1.D6,CK2SO4 *1.D6,
     &                    CK2SO4 *174.2/SMASS*1.D2, CK2SO4 /SMOLE*1.D2
         WRITE (IO,2000) '[KHSO4    ]', CKHSO4 *136.1*1.D6,CKHSO4 *1.D6,
     &                    CKHSO4 *136.1/SMASS*1.D2, CKHSO4 /SMOLE*1.D2
         WRITE (IO,2000) '[KNO3     ]', CKNO3  *101.1*1.D6,CKNO3  *1.D6,
     &                    CKNO3  *101.1/SMASS*1.D2, CKNO3  /SMOLE*1.D2
         WRITE (IO,2000) '[KCl      ]', CKCL   *74.5 *1.D6,CKCL   *1.D6,
     &                    CKCL   *74.5 /SMASS*1.D2, CKCL   /SMOLE*1.D2
         WRITE (IO,2000) '[MgSO4    ]', CMGSO4 *120.3*1.D6,CMGSO4 *1.D6,
     &                    CMGSO4 *120.3/SMASS*1.D2, CMGSO4 /SMOLE*1.D2
         WRITE (IO,2000) '[MgNO32   ]', CMGNO32*148.3*1.D6,CMGNO32*1.D6,
     &                    CMGNO32*148.3/SMASS*1.D2, CMGNO32 /SMOLE*1.D2
         WRITE (IO,2000) '[MgCl2    ]', CMGCL2 *95.3 *1.D6,CMGCL2 *1.D6,
     &                    CMGCL2 *95.3 /SMASS*1.D2, CMGCL2 /SMOLE*1.D2
      ELSE
         WRITE (IO,3000) 'NO SOLID AEROSOL PHASE'
      ENDIF
      WRITE (IO,3000) ' '
C
C *** LIQUID AEROSOL PHASE **************
C
      IF (WATER.GT.TINY) THEN
         WRITE (IO,1500)'LIQUID AEROSOL:','(ug/m3 air)','(umol/m3 air)',
     &                   '(% mass)', '(% mole)'
         WRITE (IO,2000)'[WATER ]', WATER*1.D9,         WATER*1.D9/18.,
     &                   WATER  * 1.D3/LMASS*1.D2, WATER/18.0/LMOLE*1.D5
         WRITE (IO,2000)'[H+    ]', MOLAL(1)*1.D6,      MOLAL(1)*1.D6,
     &                   MOLAL(1)*1.0 /LMASS*1.D2, MOLAL(1)  /LMOLE*1.D2
         WRITE (IO,2000)'[Na+   ]', MOLAL(2)*23.0*1.D6, MOLAL(2)*1.D6,
     &                   MOLAL(2)*23.0/LMASS*1.D2, MOLAL(2)  /LMOLE*1.D2
         WRITE (IO,2000)'[NH4+  ]', MOLAL(3)*18.0*1.D6, MOLAL(3)*1.D6,
     &                   MOLAL(3)*18.0/LMASS*1.D2, MOLAL(3)  /LMOLE*1.D2
         WRITE (IO,2000)'[Cl-   ]', MOLAL(4)*35.5*1.D6, MOLAL(4)*1.D6,
     &                   MOLAL(4)*35.5/LMASS*1.D2, MOLAL(4)  /LMOLE*1.D2
         WRITE (IO,2000)'[NO3-  ]', MOLAL(7)*62.0*1.D6, MOLAL(7)*1.D6,
     &                   MOLAL(7)*62.0/LMASS*1.D2, MOLAL(7)  /LMOLE*1.D2
         WRITE (IO,2000)'[SO4-- ]', MOLAL(5)*96.0*1.D6, MOLAL(5)*1.D6,
     &                   MOLAL(5)*96.0/LMASS*1.D2, MOLAL(5)  /LMOLE*1.D2
         WRITE (IO,2000)'[HSO4- ]', MOLAL(6)*97.0*1.D6, MOLAL(6)*1.D6,
     &                   MOLAL(6)*97.0/LMASS*1.D2, MOLAL(6)  /LMOLE*1.D2
         WRITE (IO,2000)'[Ca   ]', MOLAL(8)*40.1*1.D6, MOLAL(8)*1.D6,
     &                   MOLAL(8)*40.1/LMASS*1.D2, MOLAL(8)  /LMOLE*1.D2
         WRITE (IO,2000)'[K    ]', MOLAL(9)*39.1*1.D6, MOLAL(9)*1.D6,
     &                   MOLAL(9)*39.1/LMASS*1.D2, MOLAL(9)  /LMOLE*1.D2
         WRITE (IO,2000)'[Mg   ]',MOLAL(10)*24.3*1.D6,MOLAL(10)*1.D6,
     &                   MOLAL(10)*24.3/LMASS*1.D2, MOLAL(10)/LMOLE*1.D2
         WRITE (IO,2000)'[NH3aq ]', GASAQ(1)*17.0*1.D6, GASAQ(1)*1.D6,
     &                   GASAQ(1)*17.0/LMASS*1.D2, GASAQ(1)  /LMOLE*1.D2
         WRITE (IO,2000)'[HClaq ]', GASAQ(2)*36.5*1.D6, GASAQ(2)*1.D6,
     &                   GASAQ(2)*36.5/LMASS*1.D2, GASAQ(2)  /LMOLE*1.D2
         WRITE (IO,2000)'[HNO3aq]', GASAQ(3)*63.0*1.D6, GASAQ(3)*1.D6,
     &                   GASAQ(3)*63.0/LMASS*1.D2, GASAQ(3)  /LMOLE*1.D2
C         WRITE (IO,2000)'[H2SO4 ]', CH2SO4  *98.0*1.D6, CH2SO4  *1.D6
         IF (MOLAL(1).GT.0.0) 
     &   WRITE (IO,2000) '[pH    ]',-LOG10(MOLAL(1)/WATER)
         WRITE (IO,2000) '[IONIC STRENGTH]', IONIC
      ELSE
         WRITE (IO,3000) 'NO LIQUID AEROSOL PHASE'
      ENDIF
      WRITE (IO,3000) ' '
C
C *** MASS BALANCE **********************
C
      WRITE (IO,3000) '*** [ MASS BALANCE (% ERROR) ] ***'
      WRITE (IO,3000) ' '
C
      IF (W(1).GT.TINY) THEN
         ERRNA  = 2.D0*CNA2SO4  + CNAHSO4 + CNANO3 + CNACL + MOLAL(2) -
     &            W(1)
         WRITE (IO,2500) 'Na  TOTAL  [', ERRNA/W(1)*100.0D0, ']'
      ENDIF
C
      IF (W(2).GT.TINY) THEN
         ERRSO4 = 2.D0*CLC  + CNH42S4  + CNH4HS4 +  MOLAL(5) + MOLAL(6)+
     &            CNA2SO4  + CNAHSO4  + CCASO4 + CMGSO4 + CK2SO4 +
     &            CKHSO4   - W(2)
         WRITE (IO,2500) 'SO4 TOTAL  [', ERRSO4/W(2)*100.0D0, ']'
      ENDIF
C
      IF (W(3).GT.TINY) THEN
         ERRNH4 = 3.D0*CLC  + 2.0*CNH42S4 + CNH4HS4 + CNH4NO3 + CNH4CL +
     &            MOLAL(3)  + GASAQ(1)    + GNH3    - W(3)
         WRITE (IO,2500) 'NH4 TOTAL  [', ERRNH4/W(3)*100.0D0, ']'
      ENDIF
C
      IF (W(4).GT.TINY) THEN
         ERRNO3 = CNH4NO3 + CNANO3 + MOLAL(7) + GASAQ(3) + GHNO3 +
     &            CKNO3   + 2D0*CMGNO32 + 2D0*CCANO32 - W(4)
         WRITE (IO,2500) 'NO3 TOTAL  [', ERRNO3/W(4)*100.0D0, ']'
      ENDIF
C
      IF (W(5).GT.TINY) THEN
         ERRCL  = CNH4CL + CNACL + MOLAL(4) + GASAQ(2) + GHCL +
     &            CKCL   + 2D0*CMGCL2+ 2D0*CCACL2 - W(5)
         WRITE (IO,2500) 'CL  TOTAL  [', ERRCL/W(5)*100.0D0, ']'
      ENDIF
C
      IF (W(6).GT.TINY) THEN
         ERRCA  = CCASO4 + CCANO32 + CCACL2 + MOLAL(8) - W(6)
         WRITE (IO,2500) 'CA  TOTAL  [', ERRCA/W(6)*100.0D0, ']'
      ENDIF
C
      IF (W(7).GT.TINY) THEN
         ERRK   = CKCL + CKNO3 + 2D0*CK2SO4 + CKHSO4 + MOLAL(9) -
     &            W(7)
         WRITE (IO,2500) 'K   TOTAL  [', ERRK/W(7)*100.0D0, ']'
      ENDIF
C
      IF (W(8).GT.TINY) THEN
         ERRMG  = CMGSO4 + CMGCL2 + CMGNO32 + MOLAL(10) - W(8)
         WRITE (IO,2500) 'MG  TOTAL  [', ERRMG/W(8)*100.0D0, ']'
      ENDIF
C      WRITE (IO,2600) 'ACTIVITY CALLED  [', ICLACT, '] TIMES'
      WRITE (IO,3000) ' '
C
C *** CHARGE BALANCE **********************
C
      IF (WATER.GT.TINY) THEN
         CHPL  = MOLAL(1) + MOLAL(2) + MOLAL(3) + 2.D0*MOLAL(8) +
     &           MOLAL(9) + 2.D0*MOLAL(10)
         CHNG  = MOLAL(4) + 2.D0*MOLAL(5) + MOLAL(6) + MOLAL(7) + COH
         WRITE (IO,3000) '*** [ CHARGE BALANCE (umole/m3) ] ***'
         WRITE (IO,3000) ' '
         WRITE (IO,2500) 'TOTAL POSITIVE [', CHPL*1.D6, ']'
         WRITE (IO,2500) 'TOTAL NEGATIVE [', CHNG*1.D6, ']'
         WRITE (IO,3000) ' '
      ENDIF
C
C *** ERROR STATUS REPORT ***********************
C
40    WRITE (IO,3000) '*** [ ERROR MESSAGES ] ***'
      IF (NOFER.EQ.0) THEN
         CALL ERRSTAT (IO, 0, ' ')
      ELSE
         DO 60 I=1,NOFER
            WRITE (IO,3000) ' '
            CALL ERRSTAT (IO, ERRSTKI(I), ERRMSGI(I))
60       CONTINUE
      ENDIF
      WRITE (IO,3000) ' '
C
C *** WRITE DIVIDING BAR IN RESULT FILE *******************************
C
      IF (OFIL) THEN
      WRITE (IO,3000) '================================================'
      WRITE (IO,3000) ' '
      ENDIF
C
C *** USER OUTPUT ROUTINE *********************************************
C
      CALL USROUT(2)
C
C *** REREAD DATA POINT? **********************************************
C
      IF (IFIL .AND. .NOT.EOF) GOTO 10
C
C *** INFORM USER OF OUTPUT AND STOP *******************************
C
50    CALL CHRBLN (OUTFIL,IOE)
      CALL CHRBLN (DATFIL,IDE)
      IF (IFIL) CLOSE (UNIT=II, STATUS='KEEP')
      IF (OFIL) THEN
         CLOSE (IO, STATUS='KEEP')
         WRITE(*,*) 'Results saved in file [',OUTFIL(1:IOE),']'
      ENDIF
      IF (DFIL) THEN
         CLOSE (IO, STATUS='KEEP')
         IF (NOFER.GT.0) THEN
            IF (ERRSTKI(NOFER).LT.1000)
     &      WRITE(*,*) 'Data    saved in file [',DATFIL(1:IDE),']'
         ELSE
            WRITE(*,*) 'Data    saved in file [',DATFIL(1:IDE),']'
         ENDIF
      ENDIF
      CALL USROUT (3)
      STOP
C
C *** FORMAT STATEMENTS ************************************************
C
 1500 FORMAT(1X,A,T19,A:T34,A:T50,A:T62,A)
 2000 FORMAT(1X,A,T19,1PE10.3:T35,1PE10.3:T50,0PF7.3:T62,0PF7.3)
 2500 FORMAT(1X,A,1PE10.3:A)
 2600 FORMAT(1X,A,I4:A)
 2700 FORMAT(1X,A,T19,1PE10.3:T36,0PF7.3)
 2800 FORMAT(1X,A,F6.2:A)
 3000 FORMAT(1X,A:A:A:A)
 3100 FORMAT(1X,'"NATOT" "SO4TOT" "NH4TOT" "NO3TOT" "CLTOT" "CATOT" ',
     &          '"KTOT" "MGTOT" "RH" ',
     &          '"TEMP" "GNH3" "GHCL" "GHNO3" "CNACL" "CNANO3" ',
     &          '"CNA2SO4" "CNAHSO4" "CNH4CL" "CNH4NO3" "CNH42S4" ',
     &          '"CNH4HS4" "CLC" "CCASO4" "CCANO32" "CCACL2" ',
     &          '"CK2SO4" "CKHSO4" "CKNO3" "CKCL" "CMGSO4" ',
     &          '"CMGNO32" "CMGCL2" "HLIQ" "NALIQ" "NH4LIQ" "CLLIQ" ',
     &          '"SO4LIQ" "HSO4LIQ" "NO3LIQ" "CaLIQ" "KLIQ" "MgLIQ" ',
     &          '"NH4AER" "CLAER" ',
     &          '"NO3AER" "WATER" "LMASS" "SMASS" "CASE"')
 3200 FORMAT(48(1X,1PE14.7):1X,A:A:A:A)
C
C *** END OF MAIN PROGRAM **********************************************
C
      END



C
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE INPCNF
C *** THIS SUBROUTINE IS CALLED BY MAIN, AND CONFIGURES THE DRIVER. 
C     CONFIGURATION IS FIRST READ FROM A FILE. IF THIS FILE DOES NOT
C     EXIST, THEN DEFAULTS ARE ASSIGNED.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE INPCNF (IERR)
      INCLUDE 'isrpia.inc'
      INCLUDE 'main.inc'
      LOGICAL CNFX
C
C *** DEFAULT PARAMETERS ***********************************************
C
      DFLTU  = 1.0D0         ! Default units   (0:umol/m3, 1:ug/m3)
      DFLTP  = 0.0D0         ! Default problem (0:forward, 1:reverse)
      DFLTC  = 1.0D0         ! Default concentrations (umole/m3)
      DFLTR  = 0.9D0         ! Default RH
      DFLTT  = 298.0D0       ! Default T
C
C *** READ FROM FILE ***************************************************
C
      INQUIRE (FILE=CNFFIL, EXIST=CNFX)   ! Does CNF file exist?
10    IF (CNFX) THEN
         IC=13
         OPEN (IC, FILE=CNFFIL, STATUS='OLD')
         READ (IC,'(/)')
         READ (IC,  *  ) OFIL
         READ (IC,'( )')
         READ (IC,  *  ) DFIL
         READ (IC,'(//)')
         READ (IC,  *  ) EPS
         READ (IC,'( )')
         READ (IC,  *  ) MAXIT
         READ (IC,'( )')
         READ (IC,  *  ) NSWEEP
         READ (IC,'( )')
         READ (IC,  *  ) EPSACT
         READ (IC,'( )')
         READ (IC,  *  ) NDIV 
         READ (IC,'( )')
         READ (IC,  *  ) IACALCI
         READ (IC,'( )')
         READ (IC,  *  ) NADJI
         CLOSE(IC, STATUS='KEEP')
         IERR = 0                   ! SETUP ISORROPIA
C         
         CALL SETPARM (-1, IACALCI, EPS, MAXIT, NSWEEP, EPSACT, NDIV,
     &                 NADJI)
C
C *** FILE NOT AVAILABLE, DEFAULTS *************************************
C
      ELSE
         IFIL   = .FALSE.       ! .T. input from 'INPFIL', .F. from screen
         OFIL   = .TRUE.        ! .T. report in 'OUTFIL', .F. on screen
         DFIL   = .TRUE.        ! .T. create data file, .F. dont create
         INPFIL = 'input.dat'   ! Name of input file, if IFIL=.T.
         OUTFIL = 'output.txt'  ! Name of results file, if OFIL=.T.
         DATFIL = 'output.dat'  ! Name of data file, if DFIL=.T.
         IERR   = 1             ! configuration file not found ; use default
      ENDIF
      print*,"setparm",iacalci,eps,maxit, nsweep, epsact, ndiv, nadj
C
      RETURN
C
C *** END OF SUBROUTINE INPCNF *****************************************
C
      END




C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE INPDAT
C *** THIS SUBROUTINE IS CALLED BY MAIN, AND CONFIGURES THE DRIVER. 
C     CONFIGURATION IS FIRST READ FROM A FILE. IF THIS FILE DOES NOT
C     EXIST, THEN DEFAULTS ARE ASSIGNED.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE INPDAT (EOF)
      INCLUDE 'isrpia.inc'
      INCLUDE 'main.inc'
      LOGICAL   EOF
      CHARACTER CUNIT*14, PFX*4
      PARAMETER (TWO=2.0)
C
C *** BY DEFAULT, SET EOF FLAG TO .FALSE. ******************************
C
      EOF = .FALSE.
C
C *** READ FROM FILE ***************************************************
C
      IF (IFIL) THEN    
200      READ (II,*,END=100,ERR=200) (WI(I),I=1,NCOMP), RHI, TEMPI
      ELSE
C
C *** READ FROM CONSOLE ************************************************
C 
1        CALL INPTD (DD, DFLTU,'Input units (0=umol/m3 air, '//
     &                         '1=ug/m3 air)', '(F2.0)', IERR)
         IF (IERR.NE.0) GOTO 1
         INUNIT= DD
C    
         IF (INUNIT.EQ.0) THEN
            CUNIT='(umol/m3 air)?'
         ELSE
            CUNIT='(ug/m3 air)?'
         ENDIF
C
2        CALL INPTD (DD, DFLTP,'Problem type? (0=forward, 1=reverse)',
     &                         '(F2.0)', IERR)
         IF (IERR.NE.0) GOTO 2
         IPROBI= DD
C
         IF (IPROBI.EQ.0) THEN
            PFX = 'tot '
         ELSE
            PFX = 'aer '       
         ENDIF
C
3        CALL INPTD (DD, DFLTP,'Aerosol state? (0=solid+liquid '//
     &                         'possible, 1=metastable)','(F2.0)',IERR)
         IF (IERR.NE.0) GOTO 3
         METSTBLI= DD
C    
11       CALL INPTD (WI(1),DFLTC,'Na'//PFX//' '//CUNIT,'(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 11
C
12       CALL INPTD (WI(2),DFLTC,'SO4'//PFX//CUNIT,    '(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 12
C
13       CALL INPTD (WI(3),DFLTC,'NH3'//PFX//CUNIT,    '(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 13
C
14       CALL INPTD (WI(4),DFLTC,'NO3'//PFX//CUNIT,    '(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 14
C
15       CALL INPTD (WI(5),DFLTC,'Cl'//PFX//' '//CUNIT,'(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 15
C
151      CALL INPTD (WI(6),DFLTC,'Ca'//PFX//' '//CUNIT,'(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 151
C
152      CALL INPTD (WI(7),DFLTC,'K'//PFX//' '//CUNIT,'(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 152
C
153      CALL INPTD (WI(8),DFLTC,'Mg'//PFX//' '//CUNIT,'(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 153
C
17      CALL INPTD (RHI,  DFLTR,'RH (0-1)?',          '(F5.3)',IERR)
         IF (IERR.NE.0) GOTO 17
C
18       CALL INPTD (TEMPI,DFLTT,'Temperature (K)?',   '(F6.1)',IERR)
         IF (IERR.NE.0) GOTO 18
      ENDIF
C
C *** CONVERT TO moles/m3 *********************************************
C
      IF (INUNIT.EQ.0) THEN        ! User gave umoles/m3
         DO 20 I=1,NCOMP
            WI(I) = MAX(WI(I)*1D-6,ZERO)
20       CONTINUE
      ELSE                         ! User gave ug/m3
         DO 30 I=1,NCOMP
            WI(I) = MAX(WI(I)/WMW(I)*1D-6,ZERO)
30       CONTINUE
      ENDIF
C
C *** DATA POINT READ, EOF NOT ENCOUNTERED *****************************
C
40    RETURN
C
C *** EOF ENCOUNTERED IN INPUT FILE, SET EOF=.TRUE. AND RETURN *********
C
100   EOF = .TRUE.
      GOTO 40
C
C *** END OF SUBROUTINE INPDAT *****************************************
C
      END
C=======================================================================
C
C *** ISORROPIA CODE
C *** SUBROUTINE USROUT
C *** THIS SUBROUTINE IS USED TO PROVIDE SPECIAL PRINTOUTS.
C
C *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
C *** GEORGIA INSTITUTE OF TECHNOLOGY
C *** WRITTEN BY ATHANASIOS NENES
C *** UPDATED BY CHRISTOS FOUNTOUKIS
C
C=======================================================================
C
      SUBROUTINE USROUT (ICODE)
      INCLUDE 'isrpia.inc'
      INCLUDE 'main.inc'
C
C =======================================================================
C                            INITIAL OUTPUT
C =======================================================================
C
      IF (ICODE.EQ.1) THEN
C
C =======================================================================
C                         OUTPUT AFTER EACH RUN
C =======================================================================
C
      ELSE IF (ICODE.EQ.2) THEN 
C
C =======================================================================
C                             FINAL OUTPUT
C =======================================================================
C
      ELSE IF (ICODE.EQ.3) THEN 
      ENDIF
C
C *** END OF USROUT SUBROUTINE *****************************************
C
      RETURN
      END
