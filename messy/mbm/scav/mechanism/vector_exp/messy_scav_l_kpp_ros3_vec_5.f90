       ENDDO

!CDIR NOIEXPAND
       CALL FUN(iwstop, ynew, F1, RCONST_LOC, F_LOC)

       DO jl=1,iwstop
         x1(jl) = c21/H_LOC(jl)
       ENDDO
 
! loop unrolling has to be done here ros3 16
       DO jl=1,iwstop
