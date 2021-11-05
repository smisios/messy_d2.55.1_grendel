       ENDDO

       DO jl=1,iwstop
         ERR(jl) = MAX( uround, SQRT( ERR(jl)/float(NVAR) ) )
       ENDDO

! ======= Choose the stepsize ===============================

       DO jl=1,iwstop
         factor = 0.9/ERR(jl)**(1.e0/3.e0)
         if (IsReject_loc(jl)) then
            facmax=1.0
         else
            facmax=10.0
         end if
         factor = MAX( 1.0e-1_dp, min(factor,facmax) ) ! mz_rs_20040830: dp added
         Hold(jl) = H_loc(jl)
         H_loc(jl) = min( Hmax, MAX(Hmin,factor*H_loc(jl)) )
       ENDDO

! ======= Rejected/Accepted Step ============================

! loop unrolling has to be done here ros3 20
       DO jl=1,iwstop
         if ( (ERR(jl).gt.1).and.(Hold(jl).gt.Hmin) ) then
           IsReject_loc(jl) = .true.
         else
           IsReject_loc(jl) = .false.
