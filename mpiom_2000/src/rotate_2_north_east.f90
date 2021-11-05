!----------------------------------------------------------------------
!
!     rotation of vectors: in ocean models with rotated grids velocity
!     vectors are given in the direction of grid lines and rows. they
!     have to be rotated in latitudinal and longitudinal direction.
!
!     note: this routine assumes positive meridional flow for a flow
!           from grid point(i,j) to grid point(i,j+1) and positive
!           zonal flow for a flow from grid point(i,j) to point(i+1,j).
!           this is not the case for mpiom!
!
!           if this routine is used to rotate data of mpiom, the
!           logical change_sign_v needs to be true.
!j. jungclaus: 22.01.04:
!note here for the coupling fields u-i,v_j are on the non-verlapping
! (ie-2) grid, furthermore, the velocity fields were previously
! interpolated onto the scalar points !
!
!h.haak: 07.10.2005 vectorisation and omp directives
!----------------------------------------------------------------------
SUBROUTINE rotate_2_north_east(u_i,v_j,ix,iy)
  !
  USE mo_kind, ONLY: wp
  USE mo_param1, ONLY: ie_g, je_g
  USE mo_commo1,   ONLY: alon_g,alat_g,wetol1_g
  USE mo_constants, ONLY: api, agratorad
  USE mo_units, ONLY: io_stdout

  IMPLICIT NONE
  INTEGER, INTENT(in) :: ix, iy
  REAL(wp), INTENT(inout) :: u_i(ix,iy),v_j(ix,iy)! vector component in i-direction

  !-----------------------------------------------------------------------
  !     local variables
  !-----------------------------------------------------------------------

  REAL(wp) lat(ie_g-2,je_g),lon(ie_g-2,je_g)     ! latitudes and longitudes
  REAL(wp) u_lon(ie_g-2,je_g),v_lat(ie_g-2,je_g) ! vector component in logitudinal direction

  REAL(wp) dlat_i, dlat_j,dlon_i,dlon_j,dist_i,dist_j
  REAL(wp) lat_factor

  INTEGER i,j,ip1,im1,jp1,jm1

  LOGICAL change_sign_u,change_sign_v
  !-----------------------------------------------------------------------!
  !     specification whether change in sign is needed for the input arrays
  change_sign_u=.FALSE.
  change_sign_v=.TRUE.

  !     transformation to radians
  !     -------------------------
  IF (ix .NE. ie_g-2) WRITE(io_stdout,*) 'alarm boundary in rotate!'
!$OMP PARALLEL &
!$OMP PRIVATE(i,j,ip1,jp1,im1,jm1,dlat_i,dlat_j,dlon_i,dlon_j,dist_i,dist_j)

  !$OMP DO
  DO j=1,iy
     DO i=1,ix
        lat(i,j)=alat_g(i+1,j)*agratorad
        lon(i,j)=alon_g(i+1,j)*agratorad
     END DO
  END DO
  !$OMP END DO

  !     initialization
  !     --------------
  !$OMP DO
  DO i = 1, ix
     DO j = 1, iy
        v_lat(i,j) = 0.0_wp
        u_lon(i,j) = 0.0_wp
     END DO
  END DO
  !$OMP END DO

  IF (change_sign_u) THEN
     !$OMP DO
     DO i = 1, ix
        DO j = 1, iy
           u_i(i, j) = u_i(i, j) * (-1.0_wp)
        END DO
     END DO
     !$OMP END DO
  ENDIF


  IF (change_sign_v) THEN
     !$OMP DO
     DO i = 1, ix
        DO j = 1, iy
           v_j(i, j) = v_j(i, j) * (-1.0_wp)
        END DO
     END DO
     !$OMP END DO
  ENDIF


  !     rotation
  !     --------
  !$OMP DO
  DO i = 1, ix
     DO j = 1, iy

        ip1 = i + 1
        im1 = i - 1
        jp1 = j + 1
        jm1 = j - 1

        IF (ip1 > ix) ip1 = ip1 - ix ! the 0-meridian
        IF (im1 < 1 ) im1 = ix
        IF (jp1 > iy) THEN           ! treatment of the last..
           jp1 = j
        ENDIF
        IF (jm1 < 1 ) THEN ! .. and the fist grid-row
           jm1 = j
        ENDIF

        !                 difference in latitudes
        dlat_i = lat(ip1,j) - lat(im1,j)
        dlat_j = lat(i,jp1) - lat(i,jm1)

        !                 difference in longitudes
        dlon_i = lon(ip1,j) - lon(im1,j)
        IF (dlon_i >   api)  dlon_i = dlon_i - (2._wp * api)
        IF (dlon_i < (-api)) dlon_i = dlon_i + (2._wp * api)
        dlon_j = lon(i,jp1) - lon(i,jm1)
        IF (dlon_j >   api)  dlon_j = dlon_j - (2._wp * api)
        IF (dlon_j < (-api)) dlon_j = dlon_j + (2._wp * api)

        lat_factor = COS(lat(i,j))
        dlon_i = dlon_i * lat_factor
        dlon_j = dlon_j * lat_factor

        !                 projection by scalar product
        !                 ----------------------------
        u_lon(i,j) = u_i(i,j)*dlon_i + v_j(i,j)*dlat_i
        v_lat(i,j) = u_i(i,j)*dlon_j + v_j(i,j)*dlat_j

        dist_i = SQRT(dlon_i**2+dlat_i**2)
        dist_j = SQRT(dlon_j**2+dlat_j**2)
        IF (dist_i /= 0._wp .AND. dist_j /= 0._wp) THEN
           u_lon(i,j) = u_lon(i,j)/dist_i
           v_lat(i,j) = v_lat(i,j)/dist_j
        ELSE
           u_lon(i, j) = 0.0_wp
           v_lat(i, j) = 0.0_wp
        ENDIF

        !                  absold = sqrt(u_i(i,j)**2 + v_j(i,j)**2)
        !                  absnew = sqrt(u_lon(i,j)**2 + v_lat(i,j)**2)
        !                  print*, absold, absnew

        !                 test orthogonality
        !                 ------------------
        !                  if ((dlon_i*dlon_j+dlat_j*dlat_i) > 0.1) then
        !                     write(io_stdout,*) 'orthogonal? ', i, j,      &
        !                         (dlon_i*dlon_j+dlat_j*dlat_i)
        !                  endif

     END DO
  END DO
  !$OMP END DO


  ! write back to input field!
  !$OMP DO
  DO i=1,ix
     DO j=1,iy
        u_i(i,j)=u_lon(i,j)*wetol1_g(i+1,j)
        v_j(i,j)=v_lat(i,j)*wetol1_g(i+1,j)
     END DO
  END DO

  !$OMP END DO

  !$OMP END PARALLEL


END SUBROUTINE rotate_2_north_east
