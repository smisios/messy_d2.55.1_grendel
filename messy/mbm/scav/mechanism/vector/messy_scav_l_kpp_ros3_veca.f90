      SUBROUTINE ROS3(N,TIN,Tnext,Hmin,Hmax,Hstart, &
                         y,AbsTol,RelTol)
                         
      USE messy_scav_l_kpp_s_mem

!       L-stable Rosenbrock 3(2), with 
!     strongly A-stable embedded formula for error control.  
!
!     All the arguments aggree with the KPP syntax.
!
!  INPUT ARGUMENTS:
!     y = Vector of (NVAR) concentrations, contains the
!         initial values on input
!     [T, Tnext] = the integration interval
!     Hmin, Hmax = lower and upper bounds for the selected step-size.
!          Note that for Step = Hmin the current computed
!          solution is unconditionally accepted by the error
!          control mechanism.
!     AbsTol, RelTol = (NVAR) dimensional vectors of 
!          componentwise absolute and relative tolerances.
!     FUN = name of routine of derivatives. KPP syntax.
!          See the header below.
!     JAC_SP = name of routine that computes the Jacobian, in
!          sparse format. KPP syntax. See the header below.
!     Info(1) = 1  for  autonomous   system
!             = 0  for nonautonomous system 
!
!  OUTPUT ARGUMENTS:
!     y = the values of concentrations at Tend.
!     T = equals Tend on output.
!     Info(2) = # of FUN calls.
!     Info(3) = # of JAC_SP calls.
!     Info(4) = # of accepted steps.
!     Info(5) = # of rejected steps.
!    
!     Adrian Sandu, April 1996
!     The Center for Global and Regional Environmental Research

      IMPLICIT NONE !mz_rs_20030728

      ! mz_rs_20030728+
      REAL(dp) dround, gam, c21, c31, c32, b1, b2, b3, d1, d2, d3, a21
      REAL(dp) a31, a32, alpha2, alpha3, g1, g2, g3, tau(nvect), x1(nvect), x2(nvect), x3(nvect), ytol(nvect)
      INTEGER    ier
      ! mz_rs_20030728-
      REAL(dp)     K1(nvect,NVAR), K2(nvect,NVAR), K3(nvect,NVAR)
      REAL(dp)     F1(nvect,NVAR), JAC(nvect,LU_NONZERO_V+1) ! mz_rs_20020604: 1 added
      REAL(dp)     JAC_LOC(nvect,LU_NONZERO_V+1)       ! mz_bs_20051905
      REAL(dp)     Hmin,Hmax,Hstart,ghinv(nvect),uround
      REAL(dp)     y(nvect,NVAR), ynew(nvect,NVAR)
      REAL(dp)     AbsTol(NVAR), RelTol(NVAR)
      REAL(dp)     TIN,T(nvect), Tnext, H(nvect), Hold(nvect), Tplus(nvect)
      REAL(dp)     ERR(nvect), factor, facmax
      INTEGER    n,i,j
      INTEGER  :: ncfn(nvect), njac(nvect), Naccept(nvect), Nreject(nvect)
      LOGICAL    IsReject(nvect) 
      ! EXTERNAL   FUN, JAC_SP
      INTEGER      JL
      INTEGER  ::  igoto
! iwstop - number of boxes with unfinished calculations
! iwork  - integer array containing indices of unfinished boxes
      INTEGER  :: iwstop,iwork(nvect)


!     Initialization of counters, etc.
      uround = 1.e-15
      dround = SQRT(uround)
      DO jl=1,nvect
        H(jl) = MAX(1.e-8_dp, Hstart) ! mz_rs_20040830: dp added
        T(jl) = TIN
        Tplus(jl) = T(jl)
        IsReject(jl) = .false.
      ENDDO
      gam=   .43586652150845899941601945119356e+00
      c21=  -.10156171083877702091975600115545e+01
      c31=   .40759956452537699824805835358067e+01
      c32=   .92076794298330791242156818474003e+01
       b1=   .10000000000000000000000000000000e+01
       b2=   .61697947043828245592553615689730e+01
       b3=  -.42772256543218573326238373806514e+00
       d1=   .50000000000000000000000000000000e+00
       d2=  -.29079558716805469821718236208017e+01
       d3=   .22354069897811569627360909276199e+00
       a21 = 1.e0
       a31 = 1.e0
       a32 = 0.e0
       alpha2 = gam
       alpha3 = gam
       g1=   .43586652150845899941601945119356e+00
       g2=   .24291996454816804366592249683314e+00
       g3=   .21851380027664058511513169485832e+01
       Naccept(:)  = 0
       Nreject(:)  = 0
       Ncfn(:)     = 0
       Njac(:)     = 0
       iwork(:)    = 0

! === Starting the time loop ===      
 10    continue  

       iwstop=0
       DO jl=1,nvect
         IF( idone(jl) == 1 ) THEN
           iwstop=iwstop+1
           iwork(iwstop)=jl
         ENDIF
       ENDDO

       call JAC_SP(iwork, iwstop, y, JAC_LOC)  ! mz_bs_20051905

!CDIR NODEP       
       DO i=1,iwstop
! achtung hier minimalen zeitschritt zulassen
! achtung keine endlosschleife bauen
! falls box schon fertig, soll mit hmin weiter gerechnet werden bis alle boxen
! fertig sind
         jl=iwork(i)
         Tplus(jl) = T(jl) + H(jl)
         Njac(jl) = Njac(jl) + 1 
         if ( (Tplus(jl) .gt. Tnext) .and. (T(jl) .lt. Tnext) ) then
           H(jl) = Tnext - T(jl)
           Tplus(jl) = Tnext
           H(jl) = MAX(Hmin,H(jl))
         end if
         gHinv(jl) = -1.0e0/(gam*H(jl))
       ENDDO

!CDIR NODEP       
       DO j=1,LU_NONZERO_V
         DO i=1,iwstop
           jl=iwork(i)
           JAC(jl,j) = -JAC_LOC(jl,j) 
         ENDDO
       ENDDO

!CDIR NODEP       
       DO i=1,iwstop
         jl=iwork(i)
! loop unrolling has to be done here
