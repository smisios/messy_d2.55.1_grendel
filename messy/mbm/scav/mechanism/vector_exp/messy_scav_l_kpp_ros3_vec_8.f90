       ENDDO

! ====== Error estimation ========

       DO jl=1,iwstop
         ERR(jl)=0.e0
       ENDDO

! loop unrolling has to be done here ros3 19
       DO jl=1,iwstop
