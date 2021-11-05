
         H(jl)=H_loc(i)
         T(jl)=T_loc(i)
         IsReject(jl) = IsReject_loc(i)
         idone(jl) = idone_loc(i)
       ENDDO

       if ( igoto .ge. 1) go to 10

      END SUBROUTINE ROS3
