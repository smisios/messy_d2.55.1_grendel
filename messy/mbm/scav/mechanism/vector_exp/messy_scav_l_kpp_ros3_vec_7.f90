       ENDDO

!CDIR NOIEXPAND
       call SOLVE (iwstop, JAC_loc, K3)

! ---- The Solution ---

! loop unrolling has to be done here ros3 18
       DO jl=1,iwstop
