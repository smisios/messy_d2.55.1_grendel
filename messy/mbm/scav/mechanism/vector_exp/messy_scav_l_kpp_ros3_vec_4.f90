       ENDDO

!CDIR NOIEXPAND
       call DECOMP (iwstop, JAC_LOC)

!CDIR NOIEXPAND
       call FUN(iwstop, y_loc, K1, RCONST_loc, F_loc)

! ====== AUTONOMOUS CASE ===============

!CDIR NOIEXPAND
       call SOLVE (iwstop, JAC_loc, K1)


! ----- STAGE 2 (AUTONOMOUS) -----
       DO jl=1,iwstop
