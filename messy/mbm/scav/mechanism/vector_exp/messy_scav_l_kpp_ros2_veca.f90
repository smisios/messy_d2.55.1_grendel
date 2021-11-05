SUBROUTINE ROS2(y,  t,  dtin)

  USE messy_scav_l_vec_kpp_s_mem

  ! Solve: d/dt c = f(t,c)
  !
  ! written by Edwin Spee, CWI, Amsterdam.
  ! email: Edwin.Spee@cwi.nl
  ! if non-existent, try http://edwin-spee.mypage.org/
  ! last update: July 28, 1997
  !
  ! on input:
  ! y      vector containing variable species at time=t
  ! fix    fixed species
  ! rconst reaction constants
  ! t      start time
  ! dtin   step size 
  !
  ! on output:
  ! y      vector containing variable species at time = tend
  ! other parameters are not changed
  !
  ! Uses datastructure and subroutine JACVAR_SP,DECOMP,F_VAR,SOLVE
  ! from KPP
  !
  ! Integration method for Ros2:
  ! C_{n+1} = C_n + 3/2 dt k_1 + 1/2 dt k_2
  ! k_1 = S f(t_n, C_n)
  ! k_2 = S [f(t_{n+1},C_n + dt k_1) - 2 k_1]
  !
  ! where g = 1.0 + sqrt(0.5),
  ! S = (I - g dt J ) ^ {-1}
  ! with J the Jacobian
  !
  ! mz_rs_20010423: real changed to REAL(dp)

  IMPLICIT NONE !mz_rs_20030728

  REAL(dp) y(nvect,nvar), rad(nvect,nrad), k1(nvect,nvar), k2(nvect,nvar) & ! mz_rs_20020623: rad added
    ,w1(nvect,nvar) &
    ,g,dtin,t &
    ,J_chem(nvect,LU_NONZERO_V+1) &
    ,sxmin
  INTEGER ier,it,nt,j,i

  g = 1.0 + SQRT(0.5)
  CALL JACVAR_SP(y,rad,J_chem)
  DO j=1,LU_NONZERO_V
    DO i=1,nvect
    J_chem(i,j) = -g*dtin*J_chem(i,j)
    ENDDO
  ENDDO
  !rvk optionally cut this out and replace it by directly addressed statements

  DO i=1,nvect
