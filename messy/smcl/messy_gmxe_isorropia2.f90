MODULE MESSY_GMXE_ISORROPIA2

  USE messy_main_constants_mem,          ONLY: dp, M_H2O
  USE messy_gmxe_mem,                    ONLY: td, ngas, nanions, ncations
  USE messy_gmxe_kappa,                  ONLY: salt, nss

  IMPLICIT NONE

  SAVE



! for output from isorropia
  INTEGER :: iso_l2cation(15)
  INTEGER :: iso_s2cation(19,2)
  INTEGER :: iso_l2anion(15)
  INTEGER :: iso_s2anion(19,2)
  INTEGER :: iso_gas(3)

  INTEGER :: iso_salt2kappa(23)
  INTEGER :: hso4_idx, hp_idx
  REAL(dp), PARAMETER :: TINYCO=tiny(0._dp)

CONTAINS
!----------------------------------------------------------------------------

  SUBROUTINE init_species_isorropia(idx)
    
    INTEGER  :: jc
    INTEGER  :: idx

    IF (IDX == 1) THEN

      DO jc=1,ngas
        SELECT CASE (TRIM(td%gas(jc)%name))
        CASE('NH3')
          iso_gas(1) = jc
        CASE('HNO3')
          iso_gas(2) = jc
        CASE('HCl')
          iso_gas(3) = jc
        END SELECT
      END DO
      
      iso_l2anion(:)   = 0
      iso_s2anion(:,1) = 0
      iso_s2anion(:,2) = 1
      DO jc=1,nanions
        SELECT CASE (TRIM(td%anion(jc)%name))
        CASE('SO4mm')
          td%anion(jc)%iso_idx = 2
          iso_l2anion(5)    = jc         
          iso_s2anion(5,1)  = jc         
          iso_s2anion(6,1)  = jc         
          iso_s2anion(9,1)  = jc         
          iso_s2anion(10,1) = jc         
          iso_s2anion(13,1) = jc
          iso_s2anion(17,1) = jc         
        CASE('HSO4m')
          td%anion(jc)%iso_idx = 2
          hso4_idx          = jc
          iso_l2anion(6)    = jc
          iso_s2anion(7,1)  = jc         
          iso_s2anion(8,1)  = jc         
          iso_s2anion(9,1)  = jc         
          iso_s2anion(14,1) = jc         
        CASE('Clm')
          td%anion(jc)%iso_idx = 5
          iso_l2anion(4)    = jc
          iso_l2anion(10)   = jc
          iso_s2anion(3,1)  = jc         
          iso_s2anion(4,1)  = jc         
          iso_s2anion(12,1) = jc         
          iso_s2anion(12,2) = 2         
          iso_s2anion(16,1) = jc         
          iso_s2anion(19,1) = jc         
          iso_s2anion(19,2) = 2         
        CASE('NO3m')
          td%anion(jc)%iso_idx = 4
          iso_l2anion(7)    = jc
          iso_l2anion(11)   = jc
          iso_s2anion(1,1)  = jc         
          iso_s2anion(2,1)  = jc         
          iso_s2anion(11,1) = jc         
          iso_s2anion(11,2) = 2         
          iso_s2anion(15,1) = jc         
          iso_s2anion(18,1) = jc         
          iso_s2anion(18,2) = 2  
        END SELECT
      END DO
      
      iso_l2cation(:)   = 0
      iso_s2cation(:,1) = 0
      iso_s2cation(:,2) = 1
      DO jc=1,ncations
        SELECT CASE (TRIM(td%cation(jc)%name))
        CASE('Nap')
          td%cation(jc)%iso_idx = 1
          iso_l2cation(2)   = jc
          iso_s2cation(1,1) = jc
          iso_s2cation(3,1) = jc
          iso_s2cation(5,1) = jc
          iso_s2cation(5,2) = 2
          iso_s2cation(7,1) = jc
        CASE('NH4p')
          td%cation(jc)%iso_idx = 3
          iso_l2cation(3)   = jc
          iso_l2cation(9)   = jc
          iso_s2cation(2,1) = jc
          iso_s2cation(4,1) = jc
          iso_s2cation(6,1) = jc
          iso_s2cation(6,2) = 2
          iso_s2cation(8,1) = jc
          iso_s2cation(9,1) = jc
          iso_s2cation(9,2) = 3
        CASE('Mgpp')
          td%cation(jc)%iso_idx = 8
          iso_l2cation(15)   = jc
          iso_s2cation(17,1) = jc
          iso_s2cation(18,1) = jc
          iso_s2cation(19,1) = jc
        CASE('Capp')
          td%cation(jc)%iso_idx = 6
          iso_l2cation(13)   = jc
          iso_s2cation(10,1) = jc
          iso_s2cation(11,1) = jc
          iso_s2cation(12,1) = jc
        CASE('Kp')
          td%cation(jc)%iso_idx = 7
          iso_l2cation(14)   = jc
          iso_s2cation(13,1) = jc
          iso_s2cation(13,2) = 2
          iso_s2cation(14,1) = jc
          iso_s2cation(15,1) = jc
          iso_s2cation(16,1) = jc
        CASE('Hp')
          td%HP_idx          = jc
          hp_idx             = jc
        END SELECT
      END DO
      
    END IF

    IF (IDX == 2) THEN
      
      DO jc = 1, nss
        salt(jc)%iso_idx = 0
        SELECT CASE (TRIM(salt(jc)%name))
        CASE('Nap+Clm')
          salt(jc)%iso_idx = 1
        CASE('Nap+SO4mm')
          salt(jc)%iso_idx = 2
        CASE('Nap+NO3m')
          salt(jc)%iso_idx = 3
        CASE('NH4p+SO4mm')
          salt(jc)%iso_idx = 4
        CASE('NH4p+NO3m')
          salt(jc)%iso_idx = 5
        CASE('NH4p+Clm')
          salt(jc)%iso_idx = 6
        CASE('H2SO4')
          salt(jc)%iso_idx = 7
        CASE('NH4p+HSO4m')
          salt(jc)%iso_idx = 9
        CASE('Nap+HSO4m')
          salt(jc)%iso_idx = 12
        CASE('LC')
          salt(jc)%iso_idx = 13
        CASE('Capp+SO4mm')
          salt(jc)%iso_idx = 14
        CASE('Capp+NO3m')
          salt(jc)%iso_idx = 15
        CASE('Capp+Clm')
          salt(jc)%iso_idx = 16
        CASE('Kp+SO4mm')
          salt(jc)%iso_idx = 17
        CASE('Kp+HSO4m')
          salt(jc)%iso_idx = 18
        CASE('Kp+NO3m')
          salt(jc)%iso_idx = 19
        CASE('Kp+Clm')
          salt(jc)%iso_idx = 20
        CASE('Mgpp+SO4mm')
          salt(jc)%iso_idx = 21
        CASE('Mgpp+NO3m')
          salt(jc)%iso_idx = 22
        CASE('Mgpp+Clm')
          salt(jc)%iso_idx = 23
        END SELECT
      END DO
    
!!$      DO jc=1, nss
!!$        print*, "salt_indices", salt(jc)%name, salt(jc)%iso_idx
!!$      END DO
    END IF
  END SUBROUTINE init_species_isorropia

!----------------------------------------------------------------------------

  SUBROUTINE ISORROPIA_INTERFACE(xmions, xpions, ygases, &
                               RHI, TEMPI, h2o_aer, frac_h2so4_sulf, &
                               jl, jk, jm, ldry, lhyster )

      INTEGER  :: jl, jk, jm
      REAL(dp) :: xmions(2,nanions)
      REAL(dp) :: xpions(2,ncations)
      REAL(dp) :: xmions0(2,nanions)
      REAL(dp) :: xpions0(2,ncations)
      REAL(dp) :: ygases(ngas)
      REAL(dp) :: RHI, TEMPI
      REAL(dp) :: h2o_aer, diff
      REAL(dp) :: frac_h2so4_sulf               ! fraction of sulphate from H2SO4
                                                ! might be relevant for AEROPT
      
      REAL(dp) :: WI(8), WLI(8)
      INTEGER  :: I, jc,jcomp,kcomp
      
      REAL(dp) :: W(8)
      REAL(dp) :: GAS(3)
      REAL(dp) :: AERLIQ(15),  AERSLD(19), OTHER(9), CNTRL(2), SALTS(23)
      CHARACTER(LEN=15) :: SCASE
      REAL(dp) :: delta_hp, delta_hso4, hso4_pre

      REAL(dp) :: SUM_BEFORE, SUM_AFTER

      CHARACTER(LEN=40) :: ERRMSG(25)
      INTEGER           :: ERRSTK(25), NOFER
      LOGICAL           :: STKOFL, ldry, lhyster

      EXTERNAL ISOROPIA

!      IF (jk /= 19 ) RETURN

      WI(:)      = 0._dp
      W(:)       = 0._dp
      GAS(:)     = 0._dp
      AERLIQ(:)  = 0._dp
      AERSLD(:)  = 0._dp
      h2o_aer    = 0._dp

      delta_hso4 = 0._dp
      delta_hp   = 0._dp
      SALTS(:)   = 0._dp

      NOFER      = 0

      SUM_BEFORE = 0._dp
      SUM_AFTER  = 0._dp
!!$      DO jc=1,ncations
!!$        sum_before = sum_before + xpions(1,jc) + xpions(2,jc)
!!$      END DO
!!$      DO jc=1,nanions
!!$        sum_before = sum_before + xmions(1,jc) + xmions(2,jc)
!!$      END DO

!!$      DO jc=1,ncations
!!$        sum_before = sum_before + xpions(1,jc)*td%cation(jc)%molmass + xpions(2,jc)*td%cation(jc)%molmass
!!$      END DO
!!$      DO jc=1,nanions
!!$        sum_before = sum_before + xmions(1,jc)*td%anion(jc)%molmass + xmions(2,jc)*td%anion(jc)%molmass
!!$      END DO
      xmions0(:,:) = xmions(:,:)
      xpions0(:,:) = xpions(:,:)
      IF (hso4_idx /= 0) THEN
        hso4_pre = xmions(1,hso4_idx) + xmions(2,hso4_idx)
      END IF
      DO jc = 1,nanions
        jcomp = td%anion(jc)%iso_idx
        IF (jcomp == 0) CYCLE
        WI(jcomp) = WI(jcomp) + xmions(1,jc)
        xmions(1,jc) = 0._dp
      END DO
      DO jc = 1,ncations
        jcomp = td%cation(jc)%iso_idx
        IF (jcomp == 0) CYCLE
        WI(jcomp) = WI(jcomp) + xpions(1,jc)
        xpions(1,jc) = 0._dp
      END DO

!!$      print*, "before iso anions" 
!!$      DO jc=1,nanions 
!!$        print*, "Anion: ", jc,xmions(1,jc), td%anion(jc)%name, td%anion(jc)%iso_idx
!!$      END DO
!!$      print*, "before iso cations" 
!!$      DO jc=1,ncations 
!!$        print*, "cation: ", jc,xpions(1,jc), td%cation(jc)%name, td%cation(jc)%iso_idx
!!$      END DO

!!$      ! Na+ ions
!!$      WI(1) = xpions(1,3)
!!$      ! S(VI)(2)- ions
!!$      WI(2) = (xmions(1,2) + xmions(1,3))
!!$      ! NH4+ ions
!!$      WI(3) = xpions(1,2)
!!$      ! NO3- ions
!!$      WI(4) = xmions(1,4)
!!$      ! Cl- ions
!!$      WI(5) = xmions(1,5)
!!$      ! Ca++ ions
!!$      WI(6) = xpions(1,5)
!!$      ! K+ ions
!!$      WI(7) = xpions(1,4)
!!$      ! Mg++ ions
!!$      WI(8) = xpions(1,6)

      WLI(:) = WI(:)
     
     

!!$  CNTRL(1) = IPROBI     ! 0=FORWARD PROBLEM, 1=REVERSE PROBLEM
!!$  CNTRL(2) = METSTBLI   ! 0=SOLID+LIQUID AEROSOL, 1=METASTABLE
      CNTRL(1) = 0._dp     ! 0=FORWARD PROBLEM, 1=REVERSE PROBLEM
      CNTRL(2) = 0._dp     ! 0=SOLID+LIQUID AEROSOL, 1=METASTABLE
      IF (.NOT.lhyster) CNTRL(2) = 1._dp

!!$        print*, "starting isorropia2 :", jl,jk,WI, RHI, TEMPI, CNTRL, &
!!$          "xp:1 ",xpions(1,:),"xp:2 ",xpions(2,:),&
!!$          "xm:1 ",xmions(1,:),"xm:2 ",xmions(2,:)
      DO I=1,8
        WI(I) = MAX(0._dp, WI(I))
      ENDDO

      IF (MAXVAL(WI(:)) > 0._dp) THEN
        CALL ISOROPIA (WI, RHI, TEMPI, CNTRL, W,  GAS, AERLIQ, AERSLD, &
          SCASE, OTHER, SALTS)
        CALL ISERRINF (ERRSTK, ERRMSG, NOFER, STKOFL)

!      IF (NOFER > 0) print*, "Warning: error in ISORROPIA2!", NOFER, ERRMSG
!      print*, scase,nofer,stkofl,"stack ",errstk, errmsg
!!$      
!!$      ! *** SAVE RESULTS TO WORK ARRAYS (units = mole/m3, kg/m3 for water) ****
      ! secure positive definiteness of all compounds
        GAS(:)    = MAX(GAS(:),    0._dp)
        AERLIQ(:) = MAX(AERLIQ(:), 0._dp)
        AERSLD(:) = MAX(AERSLD(:), 0._dp)
                

      
        ! Liquid aerosol species

!!$!      xpions(1,1) = AERLIQ(1)
!!$      ! Na+
!!$      xpions(1,3) = AERLIQ(2)
!!$      ! NH4+                   undissociated NH3
!!$      xpions(1,2) = AERLIQ(3) + AERLIQ(9)
!!$      ! Cl-                    undissociated HCl
!!$      xmions(1,5) = AERLIQ(4) + AERLIQ(10)
!!$      ! NO3-                   undissociated HNO3
!!$      xmions(1,4) = AERLIQ(7) + AERLIQ(11) 
!!$      ! SO4--
!!$      xmions(1,2) = AERLIQ(5)
!!$      ! HSO4-
!!$      xmions(1,3) = AERLIQ(6)
!!$      ! Ca++
!!$      xpions(1,5) = AERLIQ(13)
!!$      ! Mg++
!!$      xpions(1,6) = AERLIQ(15)
!!$      ! K+
!!$      xpions(1,4) = AERLIQ(14)
!!$      ! OH-
!!$!      xmions(1,10) = AERLIQ(12)
      
        DO jc=1,15
          jcomp = iso_l2anion(jc)
          kcomp = iso_l2cation(jc)
          IF (jcomp /= 0) &
            xmions(1,jcomp) = xmions(1,jcomp) + AERLIQ(jc)
          IF (kcomp /= 0) &
            xpions(1,kcomp) = xpions(1,kcomp) + AERLIQ(jc)
        END DO

!!$      ! Solid aerosol species
!!$      ! Na+             NaNO3        NaCl           Na2SO4            NaHSO4
!!$      xpions(2,3) = AERSLD(01) + AERSLD(03) + 2._dp * AERSLD(05) + AERSLD(07)
!!$      ! NH4+           NH4NO3        NH4Cl          (NH4)2SO4          
!!$      xpions(2,2) = AERSLD(02) + AERSLD(04) + 2._dp * AERSLD(06) + &
!!$      !                NH4HSO4        (NH4)3H(SO4)2
!!$                      AERSLD(08) + 3._dp * AERSLD(09) 
!!$      ! Cl-            NaCl          NH4Cl           CaCl2  
!!$      xmions(2,5)  = AERSLD(03) + AERSLD(04) + 2._dp * AERSLD(12) + &
!!$      !                KCl           MgCl2
!!$                      AERSLD(16) + 2._dp * AERSLD(19)
!!$      ! NO3-           NaNO3         NH4NO3          Ca(NO3)2
!!$      xmions(2,4)  = AERSLD(01) + AERSLD(02) + 2._dp * AERSLD(11) + &
!!$      !                KNO3          Mg(NO3)2
!!$                      AERSLD(15) + 2._dp * AERSLD(18)
!!$      ! SO4--          Na2SO4        (NH4)2SO4     CaSO4
!!$      xmions(2,2)  = AERSLD(05) + AERSLD(06) + AERSLD(10) + &
!!$      !                K2SO4          MgSO4        (NH4)3H(SO4)2
!!$                      AERSLD(13) + AERSLD(17) + AERSLD(09) 
!!$      ! HSO4-          NaHSO4        NH4HSO4     KHSO4      (NH4)3H(SO4)2
!!$      xmions(2,3)  = AERSLD(07) + AERSLD(08) + AERSLD(14) + AERSLD(09) 
!!$      ! Ca++           CaCl2        Ca(NO3)2     CaSO4
!!$      xpions(2,5)  = AERSLD(12) + AERSLD(11) + AERSLD(10)
!!$      ! Mg++           MgCl2        Mg(NO3)2     MgSO4
!!$      xpions(2,6)  = AERSLD(19) + AERSLD(18) + AERSLD(17)
!!$      ! K+             KCl            KNO3            K2SO4           KHSO4
!!$      xpions(2,4)  = AERSLD(16) + AERSLD(15) + 2._dp * AERSLD(13) + AERSLD(14)

        xmions(2,:) = 0._dp
        xpions(2,:) = 0._dp
        DO jc=1,19
          jcomp = iso_s2anion(jc,1)
          kcomp = iso_s2cation(jc,1)
          IF (jcomp /= 0) &
            xmions(2,jcomp) = xmions(2,jcomp) + iso_s2anion(jc,2) * AERSLD(jc)
          IF (kcomp /= 0) &
            xpions(2,kcomp) = xpions(2,kcomp) + iso_s2cation(jc,2) * AERSLD(jc)
        END DO


      ! WATER
      !     conversion from kg/m^3 -> ug/m^3
!        h2o_aer = aerliq(8) * M_H2O * 1.e6_dp
        IF (ldry) THEN
           h2o_aer = 0._dp
        ELSE
           h2o_aer = aerliq(8) * 18._dp * 1.e6_dp
        ENDIF


!!$      ! THMSS array for kappa calculations
!!$      SALTS(:)  = MAX(SALTS(:) * 1.e6_dp, 0._dp)
!!$      ! NACL
!!$      thmss(jl,jk,jm,na,awscl) = SALTS(1)
!!$      ! NA2SO4
!!$      thmss(jl,jk,jm,na,awssu) = SALTS(2)
!!$      ! NANO3
!!$      thmss(jl,jk,jm,na,awsno) = SALTS(3)
!!$      ! (NH4)2SO4
!!$      thmss(jl,jk,jm,na,awasu) = SALTS(4)
!!$      ! NH4NO3
!!$      thmss(jl,jk,jm,na,awano) = SALTS(5)
!!$      ! NH4CL
!!$      thmss(jl,jk,jm,na,awacl) = SALTS(6)
!!$      ! H2SO4
!!$      thmss(jl,jk,jm,na,awhsa) = SALTS(7)
!!$
!!$      ! NH4HSO4
!!$      thmss(jl,jk,jm,na,awahs) = SALTS(9)
!!$      ! NAHSO4
!!$      thmss(jl,jk,jm,na,awshs) = SALTS(12)
!!$      ! LC
!!$      thmss(jl,jk,jm,na,awalc) = SALTS(13)
!!$      ! CASO4
!!$      thmss(jl,jk,jm,na,awcsu) = SALTS(14) 
!!$      ! CANO32
!!$      thmss(jl,jk,jm,na,awcno) = SALTS(15)
!!$      ! CACL2
!!$      thmss(jl,jk,jm,na,awccl) = SALTS(16)
!!$      ! K2SO4
!!$      thmss(jl,jk,jm,na,awpsu) = SALTS(17)
!!$      ! KHSO4
!!$      thmss(jl,jk,jm,na,awphs) = SALTS(18)
!!$      ! KNO3
!!$      thmss(jl,jk,jm,na,awpno) = SALTS(19)
!!$      ! KCL
!!$      thmss(jl,jk,jm,na,awpcl) = SALTS(20)
!!$      ! MGSO4
!!$      thmss(jl,jk,jm,na,awmsu) = SALTS(21)
!!$      ! MGNO32
!!$      thmss(jl,jk,jm,na,awmno) = SALTS(22)
!!$      ! MGCL2
!!$      thmss(jl,jk,jm,na,awmcl) = SALTS(23)

        DO jc=1,nss
          jcomp = salt(jc)%iso_idx
          IF (jcomp == 0) CYCLE
!          print*, "in filling salts: ", jc, salt(jc)%name, jcomp, SALTS(jcomp)
          salt(jc)%mass(jl,jk,jm) = MAX(SALTS(jcomp)*1.e6_dp,0._dp)
        END DO

      ENDIF
!      if (NOFER > 0) then
!        print*, jl,jk
!        DO i=1,nofer
!          call ERRSTAT(6,errstk(i), errmsg(i))
!        enddo
!        print*, WI(1), xpions(1,3), xpions(2,3)
!      endif
      

!      print*, ygases(4), ygases(5), ygases(15)
        
      DO jc = 1,ncations
        jcomp = td%cation(jc)%iso_idx
        IF (jcomp == 0) CYCLE
        IF (TRIM(td%cation(jc)%name)=='NH4p') THEN
          DO i=1,3
            diff = WLI(jcomp) - xpions(1,jc) - xpions(2,jc) - GAS(1)
            !        print*, "stuck in NH4+"
            IF (DIFF > 0._dp) THEN
              GAS(1) = GAS(1) + diff 
            ELSE
              IF (-1._dp * DIFF <= gas(1)) THEN
                gas(1) = gas(1) + diff
              ELSE
                diff = diff + gas(1)
                gas(1) = 0._dp
                IF (-1._dp * diff <= XPIONS(2,jc)) THEN
                  XPIONS(2,jc) = XPIONS(2,jc) + diff
                ELSE
                  diff = diff + XPIONS(2,jc)
                  XPIONS(2,jc) = 0._dp
                  XPIONS(1,jc) = XPIONS(1,jc) + diff
                ENDIF
              ENDIF
            ENDIF
            IF ( ABS( WLI(jcomp) - xpions(1,jc) - xpions(2,jc) - gas(1)) < &
              SPACING(WLI(jcomp))) EXIT
          END DO
        ELSE
          DO i=1,3
            diff = WLI(jcomp) - xpions(1,jc) - xpions(2,jc)
            IF (DIFF > 0._dp) THEN
              XPIONS(1,jc) = XPIONS(1,jc) + diff 
            ELSE
              IF (-1._dp * DIFF <= XPIONS(2,jc)) THEN
                XPIONS(2,jc) = XPIONS(2,jc) + diff
              ELSE
                diff = diff + xpions(2,jc)
                XPIONS(2,jc) = 0._dp
                XPIONS(1,jc) = XPIONS(1,jc) + diff
              ENDIF
            ENDIF
            IF ( ABS( WLI(jcomp) - xpions(1,jc) - xpions(2,jc)) &
              < SPACING(WLI(jcomp))) EXIT
          END DO
        END IF
      END DO


      DO jc = 1,nanions
        jcomp = td%anion(jc)%iso_idx
        IF (jcomp == 0) CYCLE
        SELECT CASE (TRIM(td%anion(jc)%name))
        CASE('NO3m')
        ! NO3- ions
          DO i=1,3
            diff = WLI(jcomp) - xmions(1,jc) - xmions(2,jc) - gas(2)
            IF (DIFF > 0._dp) THEN
              GAS(2) = GAS(2) + diff 
            ELSE
              IF (-1._dp * DIFF <= gas(2)) THEN
                gas(2) = gas(2) + diff
              ELSE
                diff = diff + gas(2)
                gas(2) = 0._dp
                IF (-1._dp * diff <= XMIONS(2,jc)) THEN
                  XMIONS(2,jc) = XMIONS(2,jc) + diff
                ELSE
                  diff = diff + XMIONS(2,jc)
                  XMIONS(2,jc) = 0._dp
                  XMIONS(1,jc) = XMIONS(1,jc) + diff
                ENDIF
              ENDIF
            ENDIF
            IF( ABS( WLI(jcomp) - xmions(1,jc) - xmions(2,jc) - gas(2)) < &
                SPACING(WLI(4))) EXIT
          END DO
        CASE('Clm')
        ! Cl- ions
          DO i=1,3
            diff = WLI(jcomp) - xmions(1,jc) - xmions(2,jc) - gas(3)
            IF (DIFF > 0._dp) THEN
              GAS(3) = GAS(3) + diff 
            ELSE
              IF (-1._dp * DIFF <= gas(3)) THEN
                gas(3) = gas(3) + diff
              ELSE
                diff = diff + gas(3)
                gas(3) = 0._dp
                IF (-1._dp * diff <= XMIONS(2,jc)) THEN
                  XMIONS(2,jc) = XMIONS(2,jc) + diff
                ELSE
                  diff = diff + XMIONS(2,jc)
                  XMIONS(2,jc) = 0._dp
                  XMIONS(1,jc) = XMIONS(1,jc) + diff
                ENDIF
              ENDIF
            ENDIF
            IF( ABS( WLI(jcomp) - xmions(1,jc) - xmions(2,jc) - gas(3)) < &
              SPACING(WLI(jcomp)) ) EXIT
          END DO
        CASE('SO4mm')
        ! S(VI) - ions
          DO i=1,3
            diff = WLI(jcomp) - xmions(1,jc) - xmions(1,hso4_idx) - xmions(2,jc) - xmions(2,hso4_idx)
!!$            IF (DIFF .ne. 0._dp ) &
!!$              print*, "diff: so42m 0: ", i, diff, WLI(jcomp), xmions(1,jc), xmions(1,hso4_idx), xmions(2,jc), xmions(2,hso4_idx)
            IF (DIFF > 0._dp) THEN
              xmions(1,jc) = xmions(1,jc) + diff 
            ELSE IF (DIFF < 0._dp) THEN
              IF (-1._dp * DIFF <= xmions(1,jc)) THEN
                xmions(1,jc) = xmions(1,jc) + diff
              ELSE
                diff = diff + xmions(1,jc)
                xmions(1,jc) = 0._dp
                IF (-1._dp * diff <= XMIONS(2,jc)) THEN
                  XMIONS(2,jc) = XMIONS(2,jc) + diff
                ELSE
                  diff = diff + XMIONS(2,jc)
                  XMIONS(2,jc) = 0._dp
                  IF (-1._dp * diff <= XMIONS(1,hso4_idx)) THEN
                    XMIONS(1,hso4_idx) = XMIONS(1,hso4_idx) + diff
                  ELSE
                    diff = diff + XMIONS(1,hso4_idx)
                    XMIONS(1,hso4_idx) = 0._dp
                    XMIONS(2,hso4_idx) = XMIONS(2,hso4_idx) + diff
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
!!$            print*, "diff: so42m: ", i, diff, WLI(jcomp), xmions(1,jc), xmions(1,hso4_idx), xmions(2,jc), xmions(2,hso4_idx)

            IF (ABS(WLI(jcomp) - xmions(1,jc) - xmions(1,hso4_idx) - xmions(2,jc) &
              - xmions(2,hso4_idx)) < SPACING(WLI(jcomp)) ) EXIT
          END DO

          frac_h2so4_sulf = 0._dp
          IF ( (xmions(1,jc) + xmions(2,jc)) > 0._dp ) &
            frac_h2so4_sulf = MIN(1._dp, salts(7) / (xmions(1,jc) + xmions(2,jc)))

        END SELECT
      END DO

        ! Gaseous aerosol species
!!$      !NH3
!!$      ygases(15) = GAS(1)
!!$      !HNO3
!!$      yGASES(4)  = GAS(2)
!!$      !HCl
!!$      ygases(5)  = GAS(3)

        ! including conversion to umol/m3
      DO jc=1,3
        jcomp = iso_gas(jc)
        IF (jcomp /= 0) &
          ygases(jcomp) = GAS(jc) * 1.e6_dp
      END DO
      

!      END IF
      ! H+
      ! change in HSO4   HSO4m_pre   HSO4m_l       HSO4m_s
      delta_hso4     =   hso4_pre - xmions(1,hso4_idx) - xmions(2,hso4_idx)
      ! change in H+  
      delta_hp   = gas(1) - GAS(2) - gas(3) + delta_hso4 
!      print*, "change in hp",xpions(1,1), delta_hp,ygases(15) ,yGASES(4) ,ygases(5), delta_hso4,xpions(1,1) + delta_hp
      xpions(1,hp_idx) = xpions(1,hp_idx) + delta_hp
!      print*, WLI(2), xmions(1,2), xmions(2,2), xmions(1,3), xmions(2,3), xmions(3,2),xmions(3,3)
 
!!$! Mass check
!!$      sum_after = 0._dp
!!$      DO jc=1,ncations
!!$        sum_after = sum_after + xpions(1,jc)*td%cation(jc)%molmass + xpions(2,jc)*td%cation(jc)%molmass
!!$      END DO
!!$      DO jc=1,nanions
!!$        sum_after = sum_after + xmions(1,jc)*td%anion(jc)%molmass + xmions(2,jc)*td%anion(jc)%molmass
!!$      END DO
!!$      do jc=1,3
!!$        jcomp = iso_gas(jc)
!!$        IF (jcomp /= 0) &
!!$          sum_after = sum_after + gas(jc)*td%gas(jcomp)%molmass
!!$      end do
!!$      IF (ABS(SUM_BEFORE - SUM_AFTER) > 1.e-10_dp * sum_after) &
!!$        print*, sum_before, sum_after, SUM_BEFORE - SUM_AFTER, &
!!$        "xmions0: ","l: ", xmions0(1,:), "s: ", xmions0(2,:), &
!!$        "xpions0: ","l: ", xpions0(1,:), "s: ", xpions0(2,:), &
!!$        "GAS: ", GAS, &
!!$        "xmions: ","l: ", xmions(1,:), "s: ", xmions(2,:), &
!!$        "xpions: ","l: ", xpions(1,:), "s: ", xpions(2,:), &
!!$        "hp: ", delta_hp, delta_hso4

      do i=1,nanions
        xmions(1,i)  =  xmions(1,i) * 1.e6_dp
        xmions(2,i)  =  xmions(2,i) * 1.e6_dp
      enddo
      do i=1,ncations
        xpions(1,i)  =  xpions(1,i) * 1.e6_dp
        xpions(2,i)  =  xpions(2,i) * 1.e6_dp
      enddo
          

    END SUBROUTINE ISORROPIA_INTERFACE

!------------------------------------------------------------------------------------


END MODULE MESSY_GMXE_ISORROPIA2
