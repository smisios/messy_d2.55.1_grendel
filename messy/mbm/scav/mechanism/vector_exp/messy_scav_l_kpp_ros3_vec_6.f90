       ENDDO

!CDIR NOIEXPAND
       CALL SOLVE (iwstop, JAC_loc, K2)
       
! ----- STAGE 3  (AUTONOMOUS) -----
       DO jl=1,iwstop
         x1(jl) = c31/H_LOC(jl)
         x2(jl) = c32/H_LOC(jl)
       ENDDO

! loop unrolling has to be done here ros3 17
       DO jl=1,iwstop
