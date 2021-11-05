!> Tidal Module
!> Real time computation of luni-solar gravitational potential
!> @author Maik Thomas, Ann-Margrit Hellmich
!> @date Last modified 9.9.2009 Malte Mueller
!! Included subroutines are:
!! alloc_mem_tidal   -- Allocation of memory
!! foreph_ini        -- Initialization of tidal module
!!    eph            -- determination of the phase with respect to 1.1.2000
!! foreph            -- main routine for calculation of realtime lunisolar
!!                      gravitational potential
!!    ephvsop87      -- ephemeredes of sun and moon
!!        sidt2      -- computes sidereal time
!!        obliq      -- calculation of obliquity of ecliptic
!!        Sun_n      -- calculation of position of sun
!!          anomaly  -- calculation of true anomaly AT and eccentric anomaly AE
!!        moon      -- calculation of position of moon
!!        aufb2
!!    eqecl          -- conversion of ecliptic into equatorial coordinates
!! tipouv            -- horizontal gradients of gravitational potential





MODULE mo_tidal
  USE mo_kind, ONLY: wp, dp, sp, i4
      USE mo_param1

      IMPLICIT NONE
      !Earth Tides ( maik thomas, emr  pers. comm )
      REAL(wp), PARAMETER :: sal=0.69_wp

      INTEGER :: mmccdt      ! Time

      REAL(wp), ALLOCATABLE :: acc(:,:),ass(:,:),acs(:,:)
      REAL(wp), ALLOCATABLE :: tipoto(:,:)
      REAL(wp), ALLOCATABLE :: colato(:,:),silato(:,:)
      REAL(wp), ALLOCATABLE :: colono(:,:),silono(:,:)

      ! coordinates of the reference tidal station diagnostic
      INTEGER, DIMENSION(191) :: icoord,jcoord

      CONTAINS

      SUBROUTINE alloc_mem_tidal
      ! Allocation of Memory

      ALLOCATE(acc(ie,je),ass(ie,je),acs(ie,je))
      ALLOCATE(tipoto(ie,je))
      ALLOCATE(colato(ie,je),silato(ie,je))
      ALLOCATE(colono(ie,je),silono(ie,je))

      END SUBROUTINE alloc_mem_tidal

      SUBROUTINE foreph_ini
      ! Initialization of tidal module
      ! Determination of Julian Day of first time step
      ! Projection of mpiom grid on tidal module internal coordinates
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS

      IMPLICIT NONE

      INTEGER :: i, j, jcc, moph

      mmccdt = 0; jcc = 0; moph = 0
      CALL eph(lyear1,lmont1,jcc,moph)
      ! FIXME: why not trot from mo_planetary_constants?
      mmccdt = (jcc + moph) * NINT(86400._wp / dt)

      IF (icontro.NE.0) THEN
      WRITE(0,*)'tidal: phase relative to 2000 :'    &
      ,'year= ',lyear1, 'month= ',lmont1, 'yearoff= ',jcc,' monoff= ',moph ,'mmccdt= ',mmccdt
      ENDIF

      DO J=1,JE
         DO I=1,IE
            colato(i,j)=cos(giph(2*i,2*j))
            silato(i,j)=sin(giph(2*i,2*j))
            colono(I,J)=COS(gila(2*I,2*J))
            silono(I,J)=SIN(gila(2*I,2*J))
            acc(i,j) = colono(I,J)*colono(I,J)
            ass(i,j) = silono(I,J)*silono(I,J)
            acs(i,j) = colono(I,J)*silono(I,J)
         ENDDO
      ENDDO

      END SUBROUTINE foreph_ini

      SUBROUTINE eph(jul,mon,jahrph,moph)
        USE mo_commo1, ONLY : nfixYearLen
        USE mo_model_time, ONLY: monlen
      !    Berechnet Phase [in Tagen]  bezueglich des 1.1.2000
      !       jahrph=phase zu beginn des Jahres
      !       gesamtphase=jahrph+moph
        INTEGER ::  jul,mon,jahrph,moph,jpl,jj,m

        ! Before year 2000
        IF(jul.LT.2000) THEN
           jahrph=0
           DO jj=jul,1999
              jpl=365

        IF ((MOD(jj,4).EQ.0  .AND. MOD(jj,100).NE.0)  .OR. MOD(jj,400).EQ.0  ) &
             jpl = 366

        if (nfixYearLen .eq. 365) jpl=365
        if (nfixYearLen .eq. 360) jpl=360

              jahrph=jahrph-jpl
           ENDDO
        ENDIF
        ! year 2000
        IF(jul.EQ.2000) jahrph=0
        ! After year 2000
        IF(jul.GT.2000) THEN
           jahrph=0
           DO jj=2000,jul-1
              jpl=365

              IF ((MOD(jj,4).EQ.0  .AND. MOD(jj,100).NE.0)  .OR. MOD(jj,400).EQ.0  ) &
             jpl = 366

              IF (nfixYearLen .EQ. 365) jpl=365
              IF (nfixYearLen .EQ. 360) jpl=360


              jahrph=jahrph+jpl
           ENDDO
        ENDIF

        moph=0

        IF (mon.GT.1)THEN
           DO m=1,mon-1
              moph=moph+monlen(m,jul)
           ENDDO
        ENDIF

      END SUBROUTINE eph

      SUBROUTINE foreph
!     calculates the realtime gravitational potential of sun & moon
!     input : mmccdt
!     output: dres(3,2)

      Use MO_PARAM1
      Use MO_PARALLEL
      USE mo_boundsexch, ONLY : bounds_exch
      Use MO_COMMO1
      Use MO_UNITS
      USE mo_constants, ONLY: agratorad
!HERE      USE mo_constants, ONLY: aradgra

      IMPLICIT NONE


      REAL(wp) :: dres(3,2),crim3,rkomp,erdrad,rekts,dekls
      REAL(wp) :: cris3, rektm,deklm,deklm2,dekls2,sidm,sidmq
      REAL(wp) :: rkosp,codm,codmq,sids,sidsq,cods,codsq,sidm2
      REAL(wp) :: sids2,sip2,cop2,argp,hamp,hasp,alatr
      INTEGER :: i,j

      mmccdt = mmccdt + 1

      rkomp = -4.113e-07_wp
      rkosp = 0.46051_wp * rkomp
      ! FIXME: replace with radius from mo_planetary constants
      erdrad = 6371000._wp

      CALL ephvsop87(dres)

      rekts=dres(1,1)
      dekls=dres(2,1)
      cris3=dres(3,1)

      rektm=dres(1,2)
      deklm=dres(2,2)
      crim3=dres(3,2)


      deklm2 = deklm * 2._wp
      dekls2 = dekls * 2._wp
      sidm   = Sin(deklm)
      sidmq  = sidm*sidm
      codm   = Cos(deklm)
      codmq  = codm*codm
      sids   = Sin(dekls)
      sidsq  = sids*sids
      cods   = Cos(dekls)
      codsq  = cods*cods
      sidm2=Sin(deklm2)
      sids2=Sin(dekls2)

      DO j=1,je
      DO i=1,ie

         sip2 = 2._wp * silato(i, j) * colato(i, j)
         cop2=colato(i,j)**2-silato(i,j)**2
         argp=alon(i,j)*agratorad
         alatr= alat(i,j)*agratorad
!HERE      argp=alon(i,j)*aradgra
!HERE      alatr= alat(i,j)*aradgra

         hamp = rektm + argp
         hasp = rekts + argp

         ! FIXME: move 1./3. to parameter
         tipoto(i,j) = erdrad * rkomp * crim3 &
              * (3._wp * (silato(i,j)**2 - 1._wp/3._wp) * (sidmq - 1._wp/3._wp)&
              &  + SIN(2._wp * alatr) * sidm2 * COS(hamp) &
              &  + colato(i,j)**2 * codmq * COS(2._wp * hamp))        &
              &  + erdrad * rkosp * cris3 &
              &    * (3._wp * (silato(i, j)**2 - 1._wp/3._wp) &
              &       * (sidsq - 1._wp/3._wp) &
              &       + SIN(2._wp * alatr) * sids2 * COS(hasp) &
              &       + colato(i,j)**2 * codsq * COS(2._wp * hasp))

      END DO
      END DO

      CALL bounds_exch(1,'p',tipoto,'ocpheme 1')
      CALL bounds_exch(1,'u+',dlxu,'foreph 7')
      CALL bounds_exch(1,'v+',dlyv,'foreph 8')
      CALL bounds_exch(1,'u+',amsuo,'foreph 9')
      CALL bounds_exch(1,'v+',amsue,'foreph 10')

!      if (p_pe == p_io) write(73)1,1,1,ie_g*je_g
!      call write_slice(73,tipoto)

      END SUBROUTINE foreph

      SUBROUTINE tipouv
      ! calculates the horizontal gradient of the tidal potential

      Use MO_PARAM1
      Use MO_PARALLEL
      USE mo_boundsexch, ONLY : bounds_exch
      Use MO_COMMO1
      Use MO_UNITS


      IMPLICIT NONE


      integer :: i,j,k


      DO k=1,ke
      DO j=2,je1
      DO i=2,ie1
         uoo(i,j,k) = uoo(i,j,k)+sal*(amsuo(i,j,k)*dt*                 &
     &   (tipoto(i,j)-tipoto(i+1,j))/dlxu(i,j))
         voe(i,j,k) = voe(i,j,k)+sal*(amsue(i,j,k)*dt*                 &
     &   (tipoto(i,j+1)-tipoto(i,j))/dlyv(i,j))
      END DO
      END DO
      END DO

      Call bounds_exch(1,'u',uoo,'ocpheme 2')
      Call bounds_exch(1,'v',voe,'ocpheme 3')

      RETURN
      END SUBROUTINE tipouv


      Subroutine ephvsop87(res2)
!     calcuates the ephemeredes of sun/moon
!     input: mmccdt : time step information
!     output: res(3,2) right ascension, declination and geocentric distance
!             of Sun and Moon
      Use MO_PARAM1
      Use MO_COMMO1
      Use MO_UNITS
!      IMPLICIT NONE

!  optional preparation for calculation of potential
!  according to ephaufb.f by Maik Thomas
!
      INTEGER fnut
      REAL(wp) :: pi2,pic,T,sidt,ecl,nutob,nutl,res(3,2),res2(3,2)

      ! FIXME: use pi from mo_constant here
      pi2 = dacos(-1.d0) * 2._wp
      ! FIXME: replace pic with agratorad from mo_constants
      pic = dacos(-1.d0)/DBLE(180.d0)

!     Transformation of Time into fractional julian centuries T
      ! FIXME: mo_planetary_constants trot ?
      t = (REAL(mmccdt-1, wp) * dt / 86400._wp)/36525._wp

!     corresponding sidereal time Greenwich
      Call sidt2(pic,pi2,T,sidt)

!  set fnut (perform nutation -> 1; don't -> 0)
      fnut=0

! C obliquity of the ecliptic
      CALL obliq(fnut,pic,pi2,T,ecl,nutob,nutl)

! C Sun
      CALL sun_n(fnut,pic,pi2,T,ecl,nutl,res)

! C Moon
      CALL moon(fnut,pic,pi2,T,ecl,nutl,res)

! C modifications in preparaion of calculation of potentials
      CALL aufb2(sidt,res,res2)
!
      END SUBROUTINE ephvsop87

      SUBROUTINE sidt2(pic,pi2,T,sidt)
!  convertion of Julian Centuries since J2000 T to siderical time sidt
!  according to Duffett, 1990
!
        ! FIXME: this looks precision sensitive: consider dp instead of wp
      REAL(wp) :: pic,pi2,T,sidt,JD,T2,T3
!
      T2=T*T
      T3=T2*T
!
!  Julian days
      jd = t * 36525.0_wp + 2451545.0_wp
!
!  mean siderial time of Greenwich in rad
      sidt = (280.46061837_dp + 360.98564736629_dp * (jd - 2451545.0_dp) &
           + 0.000387933_wp * T2 - T3/REAL(38710000, dp))*pic
!  between 0 and 2pi :
      IF (sidt .LT. 0._wp) THEN
        Call negangle2(pi2,sidt)
      ENDIF
      IF(sidt.Ge.pi2) THEN
        CALL langle2(pi2,sidt)
      ENDIF
!
      END SUBROUTINE sidt2

      Subroutine obliq(fnut,pic,pi2,T,ecl,nutob,nutl)
      !  calculation of obliquity of ecliptic
      !  according to Duffett, 1990

      INTEGER :: fnut
      REAL(wp) :: pic,pi2,T,ecl,nutob,nutl,A,B,C,T1,T2,T3
      REAL(wp) :: L1,L2,D1,D2,M1,M2,N1,N2
!
!  see page 57
! correction terms if requested
      If(fnut.Eq.1) Then
        t1 = t + 1._wp
        T2=T1*T1
!
        a = 100.0021358_wp * t1
        ! FIXME: consider using intrinsic AINT here
        b = 360._wp * (a - REAL(INT(a), wp))
        l1 = 279.6967_wp + 0.000303_wp * t2 + b
        l2 = 2.0_wp * l1 * pic
        a = 1336.855231_wp * t1
        ! FIXME: consider using intrinsic AINT here
        b = 360._wp * (a - REAL(INT(a), wp))
        d1 = 270.4342_wp - 0.001133_wp * t2 + b
        d2 = 2.0_wp * d1 * pic
        a = 99.99736056_wp * t1
        ! FIXME: consider using intrinsic AINT here
        b = 360._wp * (a - REAL(INT(a), wp))
        m1 = (358.4758_wp - 0.00015_wp * t2 + b) * pic
        a = 1325.552359_wp * t1
        b = 360._wp * (a - REAL(INT(a), wp))
        m2 = (296.1046_wp + 0.009192_wp * t2 + b) * pic
        a = 5.372616667_wp * t1
        b = 360._wp * (a - REAL(INT(a), wp))
        n1 = (259.1833_wp + 0.002078_wp * t2 - b) * pic
        n2 = 2._wp * n1
! correction term for nutation in longitude
        nutl = ((-17.2327_wp - 0.01737_wp * t1) * SIN(n1)               &
             + (-1.2729_wp - 0.00013_wp * t1) * SIN(l2) + 0.2088_wp * SIN(n2) &
             - 0.2037_wp * SIN(d2) + (0.1261_wp - 0.00031_wp * t1) * SIN(m1)  &
             + 0.0675_wp * SIN(m2) - (0.0497_wp - 0.00012_wp * t1) * SIN(l2 + m1) &
             - 0.0342_wp * SIN(d2 - n1) - 0.0261_wp * SIN(d2 + m2)      &
             + 0.0214_wp * SIN(l2 - m1) - 0.0149_wp * SIN(l2 - d2 + m2) &
             + 0.0124_wp * SIN(l2 - n1) + 0.0114_wp * SIN(d2 - m2)) &
             / 3600.0_dp * pic
! correction term for nutation in obliquity of the ecliptic
        nutob = ((9.21_wp + 0.00091_wp * t1) * COS(n1)                  &
             + (0.5522_wp - 0.00029_wp * t1) * COS(l2) - 0.0904_wp * COS(n2) &
             + 0.0884_wp * COS(d2) + 0.0216_wp * COS(l2 + m1)           &
             + 0.0183_wp * COS(d2 - n1) + 0.0113_wp * COS(d2 + m2)      &
             - 0.0093_wp * COS(l2 - m1) - 0.0066_wp * COS(l2 - n1)) &
             / 3600.0_dp * pic
!
      Else
        nutob = 0.0_wp
        nutl = 0.0_wp
      Endif
! obliquity of the ecliptic
      T2=T*T
      T3=T2*T
      c = 46.815_wp * t + 0.0006_wp * t2 - 0.00181_wp * t3
      ecl = (23.43929167_wp - c/3600.0_dp) * pic + nutob
!
      END SUBROUTINE obliq

      SUBROUTINE Sun_n(fnut,pic,pi2,T,ecl,nutl,res)
!  calculation of position of the Sun
!  according to Duffett, 1990
!
      Integer fnut
      REAL(wp) :: pic,pi2,T,T1,T2,T3,ecl,A,B,nutl
      REAL(wp) :: L,M1,EC,AT,AE
      REAL(wp) :: A1,B1,C1,D1,E1,H1,D2,D3,S1,S2,S3,SW,X1,X2,res(3,2)
!
!  see page 116
      t1 = t + 1._wp
      T2=T1*T1
      T3=T2*T1
!
      a = 100.0021359_wp * t1
      ! FIXME: consider using intrinsic AINT here
      b = 360._wp * (a - REAL(INT(a), wp))
      l = (279.69668_wp + 0.0003025_wp * t2 + b) * pic
      a = 99.99736042_wp * t1
      ! FIXME: consider using intrinsic AINT here
      b = 360._wp * (a - REAL(INT(a), wp))
      m1 = (358.47583_wp - 0.00015_wp * t2 + 0.0000033_wp * t3 + b) * pic
      ec = 0.01675104_wp - 0.0000418_wp * t1 - 0.000000126_wp * t2
!
!  true and eccentric anomaly in rad
      Call anomaly(pi2,M1,EC,AT,AE)
!
!  various arguments in rad
      a = 62.55209472_wp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360._wp * (a - REAL(INT(a), wp))
      a1 = (153.23_wp + b) * pic
      a = 125.1041894_wp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360._wp * (a - REAL(INT(a), wp))
      b1 = (216.57_wp + b) * pic
      a = 91.56766028_wp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360._wp * (a - REAL(INT(a), wp))
      c1 = (312.69_wp + b) * pic
      a = 1236.853095_wp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360._wp * (a - REAL(INT(a), wp))
      d1 = (350.74_wp - 0.00144_wp * t2 + b) * pic
      e1 = (231.19_wp + 20.2_wp * t1) * pic
      a = 183.1353208_wp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360._wp * (a - REAL(INT(a), wp))
      h1 = (353.4_wp + b) * pic
!
      d2 = (0.00134_wp * COS(a1) + 0.00154_wp * COS(b1) + 0.002_wp * COS(c1) &
           + 0.00179_wp * SIN(d1) + 0.00178_wp * SIN(e1)) * pic
      d3 = 0.00000543_wp * SIN(a1) + 0.00001575_wp * SIN(b1) &
           + 0.00001627_wp * SIN(c1) + 0.00003076_wp * COS(d1) &
           + 0.00000927_wp * SIN(h1)
!
!  geocentric ecliptic coordinates of the Sun
      S1=AT+L-M1+D2
      If(fnut.Eq.1) Then
        S1=S1+nutl
      Endif
      IF (s1 .LT. 0._wp) CALL negangle2(pi2,S1)
      If(S1.Ge.pi2) Call langle2(pi2,S1)
      s2 = 0.0_wp
      s3 = 1.0000002_wp * (1.0_wp - ec * COS(ae)) + d3
!
!  geocentric equatorial coordinates of the Sun
      SW  =  -1._wp
      Call eqecl(pic,pi2,S1,S2,X1,X2,ecl,SW)
      res(1,1)=X1
      res(2,1)=X2
      res(3,1)=S3
!
      End Subroutine sun_n

      Subroutine moon(fnut,pic,pi2,T,ecl,nutl,res)
!  calculation of position of the Moon
!  according to Duffett, 1990
!
      Integer fnut
      REAL(wp) :: pic,pi2,T,T1,T2,T3,ecl,nutl,A,B,C,SW
      REAL(wp) :: Q,M1,M2,M3,M4,M5,M6,ML,MS,MD,ME,MF,NA,S1,S2,S3,S4,E,E2
      REAL(wp) :: L,G,W1,W2,PM,MO1,MO2,MO3,X1,X2,res(3,2)
!
!  see page 157
      t1 = t + 1._wp
      T2=T1*T1
      T3=T2*T1
!
      q = t1 * 36525._wp
      m1 = q/27.32158213_wp
      !FIXME: consider using intrinsic AINT here
      m1 = 360._wp * (m1 - REAL(INT(m1), wp))
      m2 = q / 365.2596407_wp
      !FIXME: consider using intrinsic AINT here
      m2 = 360._wp * (m2 - REAL(INT(m2), wp))
      m3 = q / 27.55455094_wp
      !FIXME: consider using intrinsic AINT here
      m3 = 360._wp * (m3 - REAL(INT(m3), wp))
      m4 = q / 29.53058868_wp
      !FIXME: consider using intrinsic AINT here
      m4 = 360._wp * (m4 - REAL(INT(m4), wp))
      m5 = q / 27.21222039_wp
      !FIXME: consider using intrinsic AINT here
      m5 = 360._wp * (m5 - REAL(INT(m5), wp))
      m6 = q / 6798.363307_wp
      !FIXME: consider using intrinsic AINT here
      m6 = 360._wp * (m6 - REAL(INT(m6), wp))
      ml = 270.434164_wp + m1 - 0.001133_wp * t2 + 0.0000019_wp * t3
      ms = 358.475833_wp + m2 - 0.00015_wp * t2 + 0.0000033_wp * t3
      md = 296.104608_wp + m3 + 0.009192_wp * t2 + 0.0000144_wp * t3
      me = 350.737486_wp + m4 - 0.001436_wp * t2 + 0.0000019_wp * t3
      mf = 11.250889_wp + m5 - 0.003211_wp * t2 - 0.0000003_wp * t3
      na = (259.183275_wp - m6 + 0.002078_wp * t2 + 0.0000022_wp * t3) * pic
      s2 = SIN(na)
      a = (51.2_wp + 20.2_wp * t1) * pic
      s1 = SIN(a)
      b = (346.56_wp + 132.87_wp * t1 - 0.0091731_wp * t2) * pic
      s3 = 0.003964_wp * SIN(b)
      c = na + (275.05_wp - 2.3_wp * t1) * pic
      s4 = SIN(c)
      ml = (ml + 0.000233_wp * s1 + s3 + 0.001964_wp * s2) * pic
      ms = (ms - 0.001778_wp * s1) * pic
      md = (md + 0.000817_wp * s1 + s3 + 0.002541_wp * s2) * pic
      mf = (mf + s3 - 0.024691_wp * s2 - 0.004328_wp * s4) * pic
      me = (me + 0.002011_wp * s1 + s3 + 0.001964_wp * s2) * pic
      e = 1._wp - 0.002495_wp * t1 + 0.00000752_wp * t2
      e2 = e * e

!  ecliptic longitude MO1
      l = 6.28875_wp * SIN(md) + 1.274018_wp * SIN(2._wp * me - md)        &
           + 0.658309_wp * SIN(2._wp * me) + 0.213616_wp * SIN(2._wp * md) &
           - e * 0.185596_wp * SIN(ms) - 0.114336_wp * SIN(2._wp * mf)     &
           + 0.058793_wp * SIN(2._wp * (me - md))                          &
           + 0.057212_wp * e * SIN(2._wp * me - ms - md) + 0.05332_wp      &
           & * SIN(2._wp * me + md)                                        &
           + 0.045874_wp * e * SIN(2._wp * me - ms)                        &
           + 0.041024_wp * e * SIN(md - ms)                                &
           - 0.034718_wp * SIN(me) - e * 0.030465_wp * SIN(md + ms)        &
           + 0.015326_wp * SIN(2._wp * (me - mf))                          &
           - 0.012528_wp * SIN(2._wp * mf + md)                            &
           - 0.01098_wp * SIN(2._wp * mf - md)                             &
           + 0.010674_wp * SIN(4._wp * me - md)                            &
           + 0.010034_wp * SIN(3._wp * md)                                 &
           + 0.008548_wp * SIN(4._wp * me - 2._wp * md)                    &
           - e * 0.00791_wp * SIN(ms - md + 2._wp * me)                    &
           - e * 0.006783_wp * SIN(2._wp * me + ms)                        &
           + 0.005162_wp * SIN(md - me) + e * 0.005_wp * SIN(me + ms)      &
           + 0.003862_wp * SIN(4._wp * me)                                 &
           + e * 0.004049_wp * SIN(md - ms + 2._wp * me)                   &
           + 0.003996_wp * SIN(2._wp * (md + me))                          &
           + 0.003665_wp * SIN(2._wp * me - 3._wp * md)                    &
           + e * 0.002695_wp * SIN(2._wp * md - ms)                        &
           + 0.002602_wp * SIN(md - 2._wp * (mf + me))                     &
           + e * 0.002396_wp * SIN(2._wp * (me - md) - ms)                 &
           - 0.002349_wp * SIN(me + md)                                    &
           + e2 * 0.002249_wp * SIN(2._wp * (me - ms))                     &
           - e * 0.002125_wp * SIN(ms + 2._wp * md)                        &
           - e2 * 0.002079_wp * SIN(2._wp * ms)                            &
           + e2 * 0.002059_wp * SIN(2._wp * (me - ms) - md)                &
           - 0.001773_wp * SIN(2._wp * (me - mf) + md)                     &
           - 0.001595_wp * SIN(2._wp * (me + mf))                          &
           + e * 0.00122_wp * SIN(4._wp * me - ms - md)                    &
           - 0.00111_wp * SIN(2._wp * (md + mf))                           &
           + 0.000892_wp * SIN(md - 3._wp * me)                            &
           - e * 0.000811_wp * SIN(ms + md + 2._wp * me)                   &
           + e * 0.000761_wp * SIN(4._wp * me - ms - 2._wp * md)           &
           + e2 * 0.000704_wp * SIN(md - 2._wp * (ms + me))                &
           + e * 0.000693_wp * SIN(ms - 2._wp * (md - me))                 &
           + e * 0.000598_wp * SIN(2._wp * (me - mf) - ms)                 &
           + 0.00055_wp * SIN(md + 4._wp * me)                             &
           + 0.000538_wp * SIN(4._wp * md)                                 &
           + e * 0.000521_wp * SIN(4._wp * me - ms)                        &
           + 0.000486_wp * SIN(2._wp * md - me)                            &
           + e2 * 0.000717_wp * SIN(md - 2._wp * ms)
      MO1=ML+L*pic
      If(fnut.Eq.1) Then
        MO1=MO1+nutl
      Endif
      IF (MO1 .LT. 0._wp) CALL negangle2(pi2,MO1)
      IF(MO1 .GE. pi2) CALL langle2(pi2,MO1)

!  ecliptic latitude MO2
      g = 5.128189_wp * SIN(mf) + 0.280606_wp * SIN(md + mf)                 &
           + 0.277693_wp * SIN(md - mf) + 0.173238_wp * SIN(2._wp * me - mf) &
           + 0.055413_wp * SIN(2._wp * me + mf - md)                         &
           + 0.046272_wp * SIN(2._wp * me - mf - md)                         &
           + 0.032573_wp * SIN(2._wp * me + mf)                              &
           + 0.017198_wp * SIN(2._wp * md + mf)                              &
           + 0.009267_wp * SIN(2._wp * me - mf + md)                         &
           + 0.008823_wp * SIN(2._wp * md - mf)                              &
           + e * 0.008247_wp * SIN(2._wp * me - ms - mf)                     &
           + 0.004323_wp * SIN(2._wp * (me + md) - mf)                       &
           + 0.0042_wp * SIN(2._wp * me + md + mf)                           &
           + e * 0.003372_wp * SIN(mf - ms - 2._wp * me)                     &
           + e * 0.002472_wp * SIN(2._wp * me - md + mf - ms)                &
           + e * 0.002222_wp * SIN(2._wp * me + mf - ms)                     &
           + e * 0.002072_wp * SIN(2._wp * me - md - mf - ms)                &
           + e * 0.001877_wp * SIN(mf - ms + md)                             &
           + 0.001828_wp * SIN(4._wp * me - md - mf)                         &
           - e * 0.001803_wp * SIN(ms + mf) - 0.00175_wp * SIN(3._wp * mf)   &
           + e * 0.00157_wp * SIN(md - mf - ms) - 0.001487_wp * SIN(me + mf) &
           - e * 0.001481_wp * SIN(mf + ms + md)                             &
           + e * 0.001417_wp * SIN(mf - ms - md)                             &
           + e * 0.00135_wp * SIN(mf - ms) + 0.00133_wp * SIN(mf - me)       &
           + 0.001106_wp * SIN(mf + 3._wp * md)                              &
           + 0.00102_wp * SIN(4._wp * me - mf)                               &
           + 0.000833_wp * SIN(mf + 4._wp * me - md)                         &
           + 0.000781_wp * SIN(md - 3._wp * mf)                              &
           + 0.00067_wp * SIN(mf + 3._wp * me - 2._wp * md)                  &
           + 0.000606_wp * SIN(2._wp * me - 3._wp * mf)                      &
           + 0.000597_wp * SIN(2._wp * (me + md) - mf)                       &
           + e * 0.000492_wp * SIN(2._wp * me + md - ms - mf)                &
           + 0.00045_wp * SIN(2._wp * (md - me) - mf)                        &
           + 0.000439_wp * SIN(3._wp * me - mf)                              &
           + 0.000423_wp * SIN(mf + 2._wp * (me + md))                       &
           + 0.000422_wp * SIN(2._wp * me - 3._wp * md - mf)                 &
           - e * 0.000367_wp * SIN(mf + ms + 2._wp * me - md)                &
           - e * 0.000353_wp * SIN(mf + ms + 2._wp * me)                     &
           + 0.000331_wp * SIN(mf + 4._wp * me)                              &
           + e * 0.000317_wp * SIN(2._wp * me + md - ms + mf)                &
           + e2 * 0.000306_wp * SIN(2._wp * (me - ms) - mf)                  &
           - 0.000283_wp *SIN(md + 3._wp * mf)
      w1 = 0.0004664_wp * COS(na)
      w2 = 0.0000754_wp * COS(c)
      mo2 = g * pic * (1.0_wp - w1 - w2)

!  horizontal parallax PM
      pm = 0.950724_wp + 0.051818_wp * COS(md)                             &
           + 0.009531_wp * COS(2._wp * me - md)                            &
           + 0.007843_wp * COS(2._wp * me) + 0.002824_wp * COS(2._wp * md) &
           + 0.000857_wp * COS(2._wp * me + md)                            &
           + e * 0.000533_wp * COS(2._wp * me - ms)                        &
           + e * 0.000401_wp * COS(2._wp * me - md - ms)                   &
           + e * 0.00032_wp * COS(md - ms) - 0.000271_wp * COS(me)         &
           - e * 0.000264_wp * COS(md + ms)                                &
           - 0.000198_wp * COS(2._wp * mf - md)                            &
           + 0.000173_wp * COS(3._wp * md)                                 &
           + 0.000167_wp * COS(4._wp * me - md)- e * 0.000111_wp * COS(ms) &
           + 0.000103_wp * COS(4._wp * me - 2._wp * md)                    &
           - 0.000084_wp * COS(2._wp * md - 2._wp * me)                    &
           - e * 0.000083_wp * COS(2._wp * me + ms)                        &
           + 0.000079_wp * COS(2._wp * me + 2._wp * md)                    &
           + 0.000072_wp * COS(4._wp * me)                                 &
           + e * 0.000064_wp * COS(2._wp * me - ms + md)                   &
           - e * 0.000063_wp * COS(2._wp * me + ms - md)                   &
           + e * 0.000041_wp * COS(ms + me)                                &
           + e * 0.000035_wp * COS(2._wp * md - ms)                        &
           - 0.000033_wp * COS(3._wp * md - 2._wp * me)                    &
           - 0.00003_wp * COS(md + me)                                     &
           - 0.000029_wp * COS(2._wp * (mf - me))                          &
           - e * 0.000029_wp * COS(2._wp * md + ms)                        &
           + e2 * 0.000026_wp * COS(2._wp * (me - ms))                     &
           - 0.000023_wp * COS(2._wp * (mf - me) + md)                     &
           + e * 0.000019_wp * COS(4._wp * me - md - ms)
      PM=PM*pic

!  geocentric distance MO3 in km
      ! FIXME: how is this related to radius in mo_planetary_constants?
      mo3 = 6378.14_wp/SIN(pm)

!  geocentric equatorial coordinates of the Moon
      sw = -1._wp
      Call eqecl(pic,pi2,MO1,MO2,X1,X2,ecl,SW)
      res(1,2)=X1
      res(2,2)=X2
      res(3,2)=MO3
!
      END SUBROUTINE moon

      SUBROUTINE anomaly(pi2,AM,EC,AT,AE)
!  calculation of true anomaly AT and eccentric anomaly AE
!  given mean anomaly AM and eccentricity EC
!  according to Duffett, 1990
!
      REAL(wp) :: pi2,AM,EC,AT,AE,M,D,A
!
!  see page 113
      ! FIXME: consider using intrinsic AINT here
      m = am - pi2 * REAL(INT(am / pi2), wp)
      AE=M
  1   D=AE-(EC*Sin(AE))-M
      IF (ABS(d) .GE. 0.000006_wp) THEN
        d = d/(1.0_wp - ec * COS(ae))
        AE=AE-D
        Goto 1
      Else
        a = SQRT((1.0_wp + ec)/(1.0_wp - ec))*TAN(ae/2.0_wp)
        at = 2.0_wp * ATAN(a)
      Endif
!
      END SUBROUTINE anomaly

      SUBROUTINE eqecl(pic,pi2,X,Y,P,Q,ecl,SW)
!  conversion of ecliptic into equatorial coordinates
!  according to Duffett, 1990
!
! if SW=+1: equatorial (X,Y..alpha,delta) to ecliptic (P,Q..lambda,beta)
! if SW=-1: equatorial (X,Y..lambda,beta) to ecliptic (P,Q..alpha,delta)
!
      REAL(wp) :: pic,pi2,ecl,P,Q,X,Y,SW
!
!  see page 62
      p = ATAN2((SIN(x) * COS(ecl) + TAN(y) * SIN(ecl) * sw), COS(X))
      IF (p .LT. 0._wp) CALL negangle2(pi2,p)
      IF (P .GE. pi2) CALL langle2(pi2,P)
      Q=Asin(Sin(Y)*Cos(ecl)-Cos(Y)*Sin(ecl)*Sin(X)*SW)
!
      END SUBROUTINE eqecl

      SUBROUTINE aufb2(sidt,res,res2)
!
! modifications according to "ephaufb.f" by Maik Thomas
! (rekt(rad)->sid.time.green.-r.asc.; dekl(rad); cri3->(a/r)^3
! for Sun and Moon)
!
      REAL(wp) :: sidt,h(3)
      REAL(wp) :: res(3,2),res2(3,2)
!
      res2(1,1)=sidt-res(1,1)
      res2(1,2)=sidt-res(1,2)
      res2(2,1)=res(2,1)
      res2(2,2)=res(2,2)
      h(1) = 1._wp / res(3,1)
      h(2) = 384400._wp / res(3,2)
      res2(3,1)=h(1)*h(1)*h(1)
      res2(3,2)=h(2)*h(2)*h(2)
!
      Return
      END SUBROUTINE aufb2

      SUBROUTINE negangle2(pi2,x)
!  transformation of negative angles
!  to angles in interval [0°;360°)
!
      Logical endvar
      REAL(wp) :: pi2,x
      ! FIXME: replace the following bad loop with DO WHILE
      endvar=.False.
    1 If(.Not.endvar) Then
        x=x+pi2
        endvar = (x .GE. 0._wp)
        Goto 1
      Endif
!
      Return
      END  SUBROUTINE  negangle2

      SUBROUTINE langle2(pi2,x)
!  transformation of large angles
!  to angles in interval [0°;360°)
!
      Logical endvar
      REAL(wp) :: pi2,x
!
      endvar=.False.
    1 If(.Not.endvar) Then
        x=x-pi2
        endvar=(x.Lt.pi2)
        Goto 1
      Endif
!
      Return
      END SUBROUTINE    langle2

      SUBROUTINE init_tide_timeseries
        ! Initialization of 191 stations, for which sea level
        ! will stored for each time step
        ! Information about stations:
        ! 1-102 Le Provost open ocean tide gauge stations
        ! 103-108   Stations in Gulf of Maine
        ! 109-118   MERICA Station in Hudson Bay and Strait   ( see Saucier et al. 2004)
        ! 119-124   Yellow and East China Sea (see Kang et al. 2002)   IC AH WI DH CJ SG
        ! 125-191   Global map of coastal tide stations

        USE mo_units, ONLY : io_ou_lpv
        USE mo_io_config, ONLY: next_free_unit
        USE mo_mpi,    ONLY :p_pe,p_io
        USE mo_commo1
        USE mo_grid, ONLY : p_suchij

        IMPLICIT NONE
        REAL(wp),DIMENSION(102)              ::  loncoord1,latcoord1
        REAL(wp),DIMENSION(89)               ::  loncoord2,latcoord2
        REAL(wp) :: dist
        INTEGER :: ij

        ! FIXME: whoever wrote this should be ashamed: this belongs
        ! in a data file, it's not code!
        ! 102 Le Provost open ocean tide gauge data
        DATA latcoord1  /60.200_wp, 57.150_wp, 54.033_wp, 53.600_wp, 53.517_wp, 44.483_wp, 41.417_wp, 40.300_wp, &
             39.467_wp, 37.800_wp, 37.150_wp, 33.983_wp, 33.917_wp, 32.633_wp, 32.367_wp, 28.450_wp, &
             28.233_wp, 26.583_wp, 26.450_wp, 24.767_wp, 14.917_wp, 14.700_wp,  7.000_wp,  0.933_wp, &
             0.233_wp,  0.017_wp, -0.017_wp, -1.417_wp, -3.833_wp, -7.917_wp,-15.917_wp,-17.067_wp, &
             -18.050_wp,-20.500_wp,-38.517_wp,-53.533_wp,-54.283_wp,-54.517_wp,-56.483_wp,-56.700_wp, &
             -60.050_wp,-61.467_wp,  4.233_wp,  4.190_wp, -0.687_wp, -4.617_wp, -9.400_wp,-10.417_wp, &
             -10.433_wp,-12.117_wp,-12.783_wp,-19.668_wp,-20.150_wp,-28.317_wp,-37.017_wp,-37.883_wp, &
             -50.033_wp,-53.000_wp,-60.017_wp, 59.067_wp, 56.133_wp, 53.417_wp, 53.317_wp, 52.833_wp, &
             49.583_wp, 46.767_wp, 43.333_wp, 40.683_wp, 36.500_wp, 33.467_wp, 32.000_wp, 28.217_wp, &
             27.083_wp, 25.000_wp, 24.300_wp, 23.867_wp, 21.533_wp, 19.283_wp, 18.933_wp, 18.717_wp, &
             16.733_wp, 15.227_wp,  9.500_wp,  7.107_wp,  6.983_wp,  1.983_wp,  1.362_wp, -0.433_wp, &
             -2.017_wp, -2.867_wp, -8.525_wp, -8.917_wp,-14.283_wp,-17.517_wp,-17.750_wp,-21.198_wp, &
                    -23.117_wp,-27.150_wp,-29.067_wp,-31.533_wp,-33.617_wp,-54.500_wp /

        DATA loncoord1 / -28.767_wp, -10.083_wp, -52.783_wp, -13.850_wp, -25.100_wp, -40.500_wp, -27.950_wp, -15.050_wp, &
             -31.117_wp, -67.983_wp, -20.083_wp, -29.400_wp, -41.183_wp, -16.917_wp, -64.700_wp, -76.800_wp, &
             -67.533_wp, -43.967_wp, -69.317_wp, -89.650_wp, -23.500_wp, -48.833_wp, -51.550_wp, -29.283_wp, &
             -41.217_wp, -20.017_wp,  -9.983_wp,   5.617_wp, -32.400_wp, -14.417_wp,  -5.700_wp, -13.667_wp, &
             -36.133_wp, -29.333_wp, -11.150_wp, -57.017_wp, -36.500_wp,   3.333_wp, -62.983_wp, -52.533_wp, &
             -47.083_wp, -61.283_wp,  52.867_wp,  73.527_wp,  73.152_wp,  55.450_wp,  46.200_wp, 105.667_wp, &
             56.667_wp,  96.883_wp,  45.250_wp,  63.418_wp,  57.483_wp,  66.833_wp, 132.017_wp,  77.583_wp, &
             132.149_wp,  73.400_wp, 132.117_wp,-175.133_wp,-144.367_wp, 154.267_wp,-135.633_wp, 173.183_wp, &
             -132.783_wp,-130.817_wp,-160.067_wp,-169.350_wp,-163.967_wp,-134.183_wp, 149.800_wp,-177.367_wp, &
             142.183_wp,-133.917_wp, 153.967_wp,-166.289_wp,-109.933_wp, 166.617_wp,-152.483_wp,-111.017_wp, &
             -169.533_wp, 145.742_wp, 138.100_wp, 171.373_wp, 158.233_wp,-157.467_wp, 172.929_wp, -90.283_wp, &
             147.267_wp, -95.017_wp, 179.207_wp,-140.067_wp,-170.683_wp,-149.500_wp, 168.317_wp,-159.770_wp, &
             -134.950_wp,-109.433_wp, 167.933_wp, 159.067_wp, -78.833_wp, 159.000_wp /

        ! -----    Stations Hudson Bay and Gulf of Maine  ------- !
        ! 1-6             St.Helena Boston Portland Eastport St.John Halifax
        ! 7-14    MERICA       2        4        6      7     8     9      10    25    ( see Saucier et al. 2004)
        ! 15-16    Leaf Basin and  Hudson Strait
        ! 17-22    Yellow and East China Sea (see Kang et al. 2002)   IC AH WI DH CJ SG
        ! 23-86    Global map of tide stations

        DATA latcoord2 / -15.5_wp, 42.4_wp, 43.8_wp, 44.8_wp, 45.2_wp, &
             44.7_wp, 60.5_wp, 60.75_wp, 60.75_wp, 64.75_wp, 65.5_wp, &
             63.75_wp, 62.75_wp, 64.75_wp, 58.75_wp, 62.0_wp, 37.47_wp, &
             35.67_wp, 35.62_wp, 34.67_wp, 33.95_wp, 33.23_wp, &
             46.23_wp, 44.67_wp, 45.25_wp, 44.54_wp, 43.39_wp, 42.21_wp, &
             41.3_wp,  41.21_wp, 39.21_wp, 38.58_wp, 34.14_wp, 32.47_wp, &
             32.02_wp, 30.24_wp, 24.5_wp,  30.24_wp,  9.21_wp, 10.23_wp, &
             -25.01_wp, -34.11_wp, -33.54_wp, -29.15_wp, 38.42_wp, 42.14_wp, &
             43.22_wp, 48.23_wp, 50.1_wp, 51.7_wp, 55.37_wp, 53.63_wp, &
             51.12_wp, 58.0_wp, 60.15_wp, -23.39_wp, -2.12_wp, &
             1.5_wp, 3.54_wp, 8.58_wp, 16.5_wp, 31.51_wp, 32.52_wp, &
             32.71_wp, 32.72_wp, 37.81_wp, 41.45_wp, 44.38_wp, 46.13_wp, &
             48.22_wp, 48.42_wp, 54.32_wp, 55.33_wp, 57.03_wp,  60.34_wp, &
             51.52_wp,  19.73_wp, 21.31_wp, 21.96_wp, 19.28_wp, 16.74_wp, &
             28.22_wp, 42.97_wp, 39.07_wp, 33.47_wp, 31.57_wp, 32.73_wp, &
             -33.51_wp , 52.76_wp /

        DATA loncoord2 / -5.5_wp, -71.0_wp, -70.0_wp, -67.0_wp, -66.0_wp, &
             -63.5_wp, -82.1_wp, -87.0_wp,-92.0_wp, -87.0_wp, -81.9_wp, &
             -80.0_wp,-77.1_wp, -81.0_wp, -69.83_wp, -71.58_wp, 126.59_wp, &
             126.13_wp, 126.30_wp, 125.43_wp, 126.30_wp, 126.57_wp, &
             -63.12_wp, -63.52_wp, -66.06_wp, -66.59_wp, -70.15_wp, &
             -71.03_wp, -71.2_wp,-72.05_wp,-74.25_wp, -74.58_wp, -77.57_wp, &
             -79.56_wp, -80.54_wp, -81.26_wp, -81.8_wp, -87.13_wp, &
             -79.55_wp, -75.32_wp, -47.56_wp, 18.16_wp, 18.26_wp, &
             16.52_wp, -9.25_wp, -8.44_wp, -8.24_wp, -4.3_wp, -5.54_wp, &
             -5.01_wp, -7.33_wp, 0.19_wp, 1.32_wp, 7.57_wp, -1.14_wp, &
             -70.24_wp, -80.55_wp, -78.44_wp, -77.06_wp, -79.34_wp, &
             -99.55_wp, -116.38_wp, -117.15_wp, -117.17_wp, -117.17_wp, &
             -122.47_wp, -124.11_wp, -124.03_wp, -123.46_wp, &
             -124.37_wp, -123.37_wp, -130.17_wp, -131.63_wp, -135.21_wp, &
             -145.45_wp,-176.38_wp, -155.06_wp,-157.87_wp,  -159.36_wp, &
             -166.62_wp, -169.53_wp, -177.37_wp, 144.38_wp, 141.72_wp, &
             135.78_wp, 131.42_wp , 129.87_wp, 151.14_wp , 1.85_wp/



        IF (p_pe==p_io)     THEN
          PRINT *,"Init tide stations"
        ENDIF

        DO ij = 1, 102
          CALL p_suchij(latcoord1(ij), loncoord1(ij), 1, icoord(ij), jcoord(ij), dist , 1._wp)
          !if (p_pe==p_io)     then
          !   print *,"Station",ij,  icoord(ij),jcoord(ij)
          !   print *,"Coordi1",ij,  latcoord1(ij),loncoord1(ij)
          !   print *,"Coordi2",ij,  180/3.14_wp*giph_g(2*icoord(ij),2*jcoord(ij)),180/3.14_wp*gila_g(2*icoord(ij),2*jcoord(ij))
          !endif
        ENDDO

        DO ij = 1, 89
          CALL p_suchij(latcoord2(ij), loncoord2(ij), 1, icoord(ij+102), jcoord(ij+102), dist , 1._wp)
          !   if (p_pe==p_io) write(0,*) ij,icoord(ij+102),jcoord(ij+102)
          !if (p_pe==p_io)     then
          !   WRITE(0,*) "Station",ij+102,  icoord(ij+102),jcoord(ij+102)
          !   WRITE(0,*) "Coordi1",ij,  latcoord2(ij),loncoord2(ij)
          !   WRITE(0,*) "Coordi2",ij, &
          !        180/3.14_wp*giph_g(2*icoord(ij+102),2*jcoord(ij+102)),180/3.14_wp*gila_g(2*icoord(ij+102),2*jcoord(ij+102))
          !endif
        ENDDO
        !        if (p_pe==p_io) print *, "MERICA Station 7: Hardwired: 507,17"
        !        icoord(102+12)=507; jcoord(102+12)=17

        PRINT *,"Open File tidestation.dat:"

        io_ou_lpv = next_free_unit()

        OPEN(io_ou_lpv,FILE='tidestation_mpiom',form='unformatted')


      END SUBROUTINE init_tide_timeseries


      SUBROUTINE p_tide_timeseries

        USE mo_mpi,    ONLY : p_pe, p_io
        USE mo_parallel, ONLY : p_ioff, p_joff, global_sum, have_g_js
        USE mo_commo1, ONLY  : zo,lbounds_exch_tp,lyear,lmonts,lday
        USE mo_units, ONLY : io_ou_lpv

        IMPLICIT NONE

        INTEGER                          :: ij,ii,jj,jb,idate

        INTEGER(i4)                 :: j1,j2,j3,j4
        REAL(dp), DIMENSION(191)              :: sl108

        jb = 2
        IF ( lbounds_exch_tp .AND. have_g_js ) jb = 3

        DO ij=1,191
          ii=icoord(ij)-p_ioff
          jj=jcoord(ij)-p_joff
          IF (ii .GE. 2 .AND. ii .LE. ie-1 &
               .AND. jj .GE. jb .AND. jj .LE. je - 1) THEN
            sl108(ij) = zo(ii,jj)
          ELSE
            sl108(ij) = 0.0_dp
         ENDIF
       ENDDO

       CALL global_sum(sl108)

       IF(p_pe==p_io) THEN
         !    print *,"write stations ..."
          idate=lyear*10000+lmonts*100+lday
          j1=INT(idate,i4)
          j2=198
          j3=0
          j4=191

          WRITE(io_ou_lpv) j1,j2,j3,j4
          WRITE(io_ou_lpv) REAL(sl108,sp)

!         WRITE(io_ou_lpv,'(102f9.5)') (sl108(ij),ij=1,102)
!         WRITE(io_ou_lpv,'(89f9.5)') (sl108(ij),ij=103,191)

       ENDIF


     END SUBROUTINE p_tide_timeseries

   END MODULE mo_tidal


