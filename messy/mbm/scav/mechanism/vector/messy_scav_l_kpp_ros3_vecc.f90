           T(jl) = Tplus(jl)  
           Naccept(jl) = Naccept(jl)+1   
         end if
       ENDDO

! ======= End of the time loop ===============================
! das folgende statement wird komplizierter
! a) es muss gewaehrleistet sein, dass der zeitschritt d.i. H immer mindestens Hmin ist.
! b) das go to darf nur nach einem globalen check (ueber alle nvect) erfolgen
!CDIR NODEP
       DO i=1,iwstop
         jl=iwork(i)
         if ( T(jl) .lt. Tnext .and. idone(jl) == 1) THEN
           idone(jl) = 1
         else
           idone(jl) = 0
           H(jl) = Hmin
         ENDIF
       ENDDO
       
       igoto = sum(idone)

       if ( igoto .ge. 1) go to 10

       Info(1:nvect,2) = Ncfn(1:nvect)
       Info(1:nvect,3) = Njac(1:nvect)
       Info(1:nvect,4) = Naccept(1:nvect)
       Info(1:nvect,5) = Nreject(1:nvect)
      
      END SUBROUTINE ROS3
