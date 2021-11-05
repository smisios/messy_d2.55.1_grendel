  ENDDO
  !***********************************************************************
  !rvk call to old routine with indirect adressing 
  CALL DECOMP (J_chem)
  ! new calls with direct addressing (needs to be preprocessed by write_decomp2d.f)
  ! uncomment if wanted
  ! additional input: sxmin = sxmin = sqrt(smallest representable number > 0
  !     sxmin = 1.e-19 ! 4.e-154 for CRAY or double precision

  !     call DECOMP_D (J_chem,ier,sxmin)  ! box model call
  !     call DECOMP_2D (J_chem,ier,sxmin)  ! 2D call
  !***********************************************************************

  CALL F_VAR(y,rad,k1)

  CALL SOLVE (J_chem,k1)

  DO i = 1,nvar
    DO j = 1,nvect
      w1(j,i) = MAX(0.0_dp, y(j,i) + dtin * k1(j,i) ) ! mz_rs_20040830: dp added
    ENDDO
  ENDDO

  CALL F_VAR(w1,rad,k2)

  DO i = 1,nvar
    DO j = 1,nvect
      k2(j,i) = k2(j,i) - 2.0 * k1(j,i)
    ENDDO
  ENDDO
  CALL SOLVE (J_chem,k2)
  DO i = 1,nvar
    DO j = 1,nvect
      ! mz_rs_20040830: dp added to next line:
      y(j,i) = MAX(0.0_dp,y(j,i) + 1.5*dtin * k1(j,i) + 0.5 *dtin*k2(j,i) )
    ENDDO
  ENDDO
END SUBROUTINE ROS2
