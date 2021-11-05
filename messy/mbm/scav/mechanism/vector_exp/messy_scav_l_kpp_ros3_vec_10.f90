         end if
       ENDDO

! ======= End of the time loop ===============================
! das folgende statement wird komplizierter
! a) es muss gewaehrleistet sein, dass der zeitschritt d.i. H immer mindestens Hmin ist.
! b) das go to darf nur nach einem globalen check (ueber alle nvect) erfolgen
       igoto = 0
       DO jl=1,iwstop
         if ( T_loc(jl) .lt. Tnext ) THEN
           idone_loc(jl) = 1
           igoto = 1
         else
           idone_loc(jl) = 0
           H_loc(jl) = Hmin
         ENDIF
       ENDDO
       
! let's globalize
       DO i=1,iwstop
         jl=iwork(i)
