       ENDDO

!CDIR NOIEXPAND
       call JAC_SP(iwstop, y_loc, JAC_LOC, RCONST_LOC, F_LOC)  ! mz_bs_20051905

       DO jl=1,iwstop
! achtung hier minimalen zeitschritt zulassen
! achtung keine endlosschleife bauen
! falls box schon fertig, soll mit hmin weiter gerechnet werden bis alle boxen
! fertig sind
         Tplus(jl) = T_loc(jl) + H_loc(jl)
         if ( (Tplus(jl) .gt. Tnext) .and. (T_loc(jl) .lt. Tnext) ) then
           H_loc(jl) = Tnext - T_loc(jl)
           Tplus(jl) = Tnext
           H_loc(jl) = MAX(Hmin,H_loc(jl))
         end if
         gHinv(jl) = -1.0e0/(gam*H_loc(jl))
       ENDDO

! loop unrolling has to be done here ros3 11
       DO jl=1,iwstop
