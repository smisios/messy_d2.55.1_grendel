       ENDDO

       call DECOMP (iwork, iwstop, JAC)

       call FUN(iwork, iwstop, y, F1)
!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         ncfn(jl) = ncfn(jl) + 1 
       ENDDO
! ====== NONAUTONOMOUS CASE ===============
         ! mz_rs_20040830: dp added to next line:
!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         tau(jl) = sign(dround*max( 1.0e-6_dp, abs(T(jl)) ), T(jl))
       ENDDO

       call FUN(iwork, iwstop, y, K2)
!CDIR NODEP
      DO i=1,iwstop
         jl=iwork(i)
         ncfn(jl) = ncfn(jl) + 1 
       ENDDO

       DO j = 1,NVAR
!CDIR NODEP
         DO i=1,iwstop
           jl=iwork(i)
           K3(jl,j) = ( K2(jl,j)-F1(jl,j) )/tau(jl)
         ENDDO
       ENDDO
 
! ----- STAGE 1 (NONAUTONOMOUS) -----
!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         x1(jl) = g1*H(jl)
       ENDDO

       DO j = 1,NVAR
!CDIR NODEP
         DO i=1,iwstop
           jl=iwork(i)
           K1(jl,j) =  F1(jl,j) + x1(jl)*K3(jl,j)
         ENDDO
       ENDDO

       call SOLVE (iwork, iwstop, JAC, K1)
      
! ----- STAGE 2 (NONAUTONOMOUS) -----
       DO j = 1,NVAR
!CDIR NODEP
         DO i=1,iwstop
           jl=iwork(i)
           ynew(jl,j) = y(jl,j) + K1(jl,j) 
         ENDDO
       ENDDO

       CALL FUN(iwork, iwstop, ynew, F1)
!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         ncfn(jl) = ncfn(jl) + 1 
       ENDDO

!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         x1(jl) = c21/H(jl)
         x2(jl) = g2*H(jl)
       ENDDO
 
       DO j = 1,NVAR
!CDIR NODEP
         DO i=1,iwstop
           jl=iwork(i)
           K2(jl,j) = F1(jl,j) + x1(jl)*K1(jl,j) + x2(jl)*K3(jl,j)
         ENDDO
       ENDDO

       CALL SOLVE (iwork, iwstop, JAC, K2)
       
! ----- STAGE 3  (NONAUTONOMOUS) -----
!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         x1(jl) = c31/H(jl)
         x2(jl) = c32/H(jl)
         x3(jl) = g3*H(jl)
       ENDDO

       DO j = 1,NVAR
!CDIR NODEP
         DO i=1,iwstop
           jl=iwork(i)
           K3(jl,j) = F1(jl,j) + x1(jl)*K1(jl,j) + x2(jl)*K2(jl,j) + x3(jl)*K3(jl,j)
         ENDDO
       ENDDO

       call SOLVE (iwork, iwstop, JAC, K3)

! ---- The Solution ---

       DO j = 1,NVAR
!CDIR NODEP
         DO i=1,iwstop
           jl=iwork(i)
           ynew(jl,j) = y(jl,j) + b1*K1(jl,j) + b2*K2(jl,j) + b3*K3(jl,j) 
         ENDDO
       ENDDO

! ====== Error estimation ========

       DO jl=1,nvect
         ERR(jl)=0.e0
       ENDDO

       DO j=1,NVAR
!CDIR NODEP
         DO i=1,iwstop
           jl=iwork(i)
           ytol(jl) = AbsTol(j) + RelTol(j)*abs(ynew(jl,j))
           ERR(jl)=ERR(jl)+((d1*K1(jl,j)+d2*K2(jl,j)+d3*K3(jl,j))/ytol(jl))**2
         ENDDO
       ENDDO

!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         ERR(jl) = MAX( uround, SQRT( ERR(jl)/float(NVAR) ) )
       ENDDO

! ======= Choose the stepsize ===============================

!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         factor = 0.9/ERR(jl)**(1.e0/3.e0)
         if (IsReject(jl)) then
            facmax=1.0
         else
            facmax=10.0
         end if 
         factor = MAX( 1.0e-1_dp, min(factor,facmax) ) ! mz_rs_20040830: dp added
         Hold(jl) = H(jl)
         H(jl) = min( Hmax, MAX(Hmin,factor*H(jl)) )    
       ENDDO

! ======= Rejected/Accepted Step ============================

!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         if ( (ERR(jl).gt.1).and.(Hold(jl).gt.Hmin) ) then
           IsReject(jl) = .true.
           Nreject(jl)  = Nreject(jl)+1
           IF ((Nreject(jl) >INFO(jl,5)) .AND. (INFO(jl,5) > 0)) THEN
             Info(jl,2) = Ncfn(jl)
             Info(jl,3) = Njac(jl)
             Info(jl,4) = Naccept(jl)
             Info(jl,5) = Nreject(jl)
             idone(jl) = 0
          ENDIF
         else
           IsReject(jl) = .false.
