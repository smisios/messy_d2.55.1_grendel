!> contains grid specific routines that are independent of global
!> state or paralllelization
MODULE mo_grid_elementals
  USE mo_constants, ONLY: agratorad
  USE mo_kind,      ONLY: dp, wp

  TYPE grid_dist_2d
    SEQUENCE
    REAL(dp) :: dist
    INTEGER  :: i, j, pe
  END TYPE grid_dist_2d

  !> subroutine to determine i and j indices of nearest point
  !> from latitude and longitude
  !>
  !> @param ie   dimension of grid searched in i direction
  !> @param je   dimension of grid searched in j direction
  !> @param alat latitude (deg.)
  !> @param alon longitude (deg.)
  !> @wetref     if land shall be considered too : .false.
  !>              ocean points only                 .true.
  !>
  !> @param ipos  i index is stored in this variable
  !> @param jpos  j index is stored in this variable
  !> @param dist  distance from point in m
  !> @param tp_grid correct upper bound for tripolar grid halo
CONTAINS
  SUBROUTINE suchij_2d(isize, jsize, alat, alon, ipos, jpos, dist, wetref, &
       weto, alat_tab, alon_tab, tp_grid)
    IMPLICIT NONE
    REAL(wp), INTENT(in) :: alat, alon
    LOGICAL, INTENT(in) :: wetref, tp_grid
    INTEGER, INTENT(in) :: isize, jsize
    INTEGER, INTENT(out) :: ipos, jpos
    REAL(wp), INTENT(out) :: dist
    REAL(wp), INTENT(in) :: alon_tab(isize, jsize), alat_tab(isize, jsize), &
         weto(isize, jsize)
    REAL(wp) :: ala, alam, alo, aphi, ax, ay, az, xx, yy, zz, d(isize)
    INTEGER :: i, j, jb

    aphi = agratorad * alat
    alam = agratorad * alon

    xx = COS(alam) * COS(aphi)
    yy = SIN(alam) * COS(aphi)
    zz = SIN(aphi)

    dist = 1.e20_wp
    ipos = -1
    jpos = -1

    ! start searches in line 3 on tripolar grids
    jb = MERGE(2, 3, .not. tp_grid)
    !
!CDIR NOLOOPCHG
    DO j = jb, jsize-1
!CDIR NOVECTOR
      DO i = 2, isize-1
        IF(.NOT. wetref .OR. weto(i,j) .GT. 0.5_wp)THEN
          alo = agratorad * alon_tab(i,j)
          ala = agratorad * alat_tab(i,j)
          ax = COS(ala) * COS(alo)
          ay = COS(ala) * SIN(alo)
          az = SIN(ala)
          d(i)=(xx-ax)**2+(yy-ay)**2+(zz-az)**2
        END IF
      END DO
    
      DO i = 2, isize-1
        IF(.NOT. wetref .OR. weto(i,j) .GT. 0.5_wp)THEN
          IF (d(i)<dist) THEN
            dist = d(i)
            ipos = i
            jpos = j
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    dist=SQRT(dist)*6350000._wp
    !       write(6, *)'in suchij ende: ', ipos, jpos, weto(ipos, jpos)
  END SUBROUTINE suchij_2d
END MODULE mo_grid_elementals
