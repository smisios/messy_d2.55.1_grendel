MODULE mo_rotation

  USE mo_kind, ONLY: i4, dp
  USE mo_commo1, ONLY: alonu, alatu, alonv, alatv, gila, giph
  USE mo_param1, ONLY: ie, je, ito, jto
  USE mo_constants, ONLY: agratorad

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rotate2_ini, rotate2_u, rotate2_v

  REAL(dp), ALLOCATABLE :: ru_poli(:,:,:), ru_polj(:,:,:), &
                           rv_poli(:,:,:), rv_polj(:,:,:)
                           

CONTAINS

  SUBROUTINE rotate2_ini

    REAL(dp)  :: gixyz(ito,jto,3), r_i(3), r_j(3), r_iabs, r_jabs

    INTEGER(i4)  :: i, j

    ALLOCATE(ru_poli(ie,je,3), ru_polj(ie,je,3), &
             rv_poli(ie,je,3), rv_polj(ie,je,3))

    DO j=1,jto
      DO i=1,ito
        CALL polxyz(1._dp, gila(i,j),  giph(i,j),  gixyz(i,j,:))
      ENDDO
    ENDDO

    DO j=2,je-1
      DO i=2,ie-1

        r_i(:)=gixyz(2*i+1,2*j,:)-gixyz(2*i,2*j,:)
        r_iabs=SQRT(SUM(r_i(:)**2))
        r_i(:)=r_i(:)/r_iabs

        CALL dek2geo(r_i, ru_poli(i,j,:), (/alonu(i,j)*agratorad,alatu(i,j)*agratorad/))


        r_j(:)=gixyz(2*i+1,2*j-1,:)-gixyz(2*i+1,2*j,:)
        r_jabs=SQRT(SUM(r_j(:)**2))
        r_j(:)=r_j(:)/r_jabs

        CALL dek2geo(r_j, ru_polj(i,j,:), (/alonu(i,j)*agratorad,alatu(i,j)*agratorad/))

        r_i(:)=gixyz(2*i,2*j+1,:)-gixyz(2*i-1,2*j+1,:)
        r_iabs=SQRT(SUM(r_i(:)**2))
        r_i(:)=r_i(:)/r_iabs

        CALL dek2geo(r_i, rv_poli(i,j,:), (/alonv(i,j)*agratorad,alatv(i,j)*agratorad/))

        r_j(:)=gixyz(2*i,2*j,:)-gixyz(2*i,2*j+1,:)
        r_jabs=SQRT(SUM(r_j(:)**2))
        r_j(:)=r_j(:)/r_jabs

        CALL dek2geo(r_j, rv_polj(i,j,:), (/alonv(i,j)*agratorad,alatv(i,j)*agratorad/))

      ENDDO
    ENDDO

  END SUBROUTINE rotate2_ini

!-------------------------------------------------------------------------------

  SUBROUTINE rotate2_u(fieldu, fieldv, ie, je)

    INTEGER(i4), INTENT(IN)     :: ie, je
    REAL(dp),    INTENT(IN)     :: fieldv(ie,je)
    REAL(dp),    INTENT(INOUT)  :: fieldu(ie,je)
    REAL(dp)                    :: fieldu_rot, fieldv_rot, absin, absrot

    INTEGER(i4)  :: i, j

    DO j=2,je-1
      DO i=2,ie-1

        absin = SQRT(fieldu(i,j)**2 + fieldv(i,j)**2)

        IF (absin .GT. 0.0_dp) THEN

          fieldu_rot = fieldu(i,j)*ru_poli(i,j,1) + fieldv(i,j)*ru_poli(i,j,2)
          fieldv_rot = fieldu(i,j)*ru_polj(i,j,1) + fieldv(i,j)*ru_polj(i,j,2)

          !H       BETRAG NACH DER DREHUNG
          absrot = SQRT(fieldu_rot**2+fieldv_rot**2)

          !H       QUOTIENT ALS KORREKTURFAKTOR

          fieldu(i,j) = fieldu_rot * absin / absrot

        ELSE

          fieldu(i,j)=0._dp
        ENDIF

      ENDDO
    ENDDO

    fieldu(:,je)=0._dp

  END SUBROUTINE rotate2_u

!-------------------------------------------------------------------------------

  SUBROUTINE rotate2_v(fieldu, fieldv, ie, je)

    INTEGER(i4), INTENT(IN)     :: ie, je
    REAL(dp),    INTENT(INOUT)  :: fieldv(ie,je)
    REAL(dp),    INTENT(IN)     :: fieldu(ie,je)
    REAL(dp)                    :: fieldu_rot, fieldv_rot, absin, absrot

    INTEGER(i4)  :: i, j

    DO j=2,je-1
      DO i=2,ie-1

        absin = SQRT(fieldu(i,j)**2 + fieldv(i,j)**2)

        IF (absin .GT. 0.0_dp) THEN

          fieldu_rot = fieldu(i,j)*rv_poli(i,j,1) + fieldv(i,j)*rv_poli(i,j,2)
          fieldv_rot = fieldu(i,j)*rv_polj(i,j,1) + fieldv(i,j)*rv_polj(i,j,2)

          !H       BETRAG NACH DER DREHUNG
          absrot = SQRT(fieldu_rot**2+fieldv_rot**2)

          !H       QUOTIENT ALS KORREKTURFAKTOR

          fieldv(i,j) = fieldv_rot * absin / absrot

        ELSE

          fieldv(i,j)=0._dp
        ENDIF

      ENDDO
    ENDDO

    fieldv(:,je)=0._dp

  END SUBROUTINE rotate2_v

!-------------------------------------------------------------------------------

  SUBROUTINE xyzpol(xyz,r,alam,phi)

    REAL(dp), INTENT(IN)  :: xyz(3)
    REAL(dp), INTENT(OUT) :: r, phi, alam
  
    r = SQRT(xyz(1)*xyz(1)+xyz(2)*xyz(2)+xyz(3)*xyz(3))
    phi = ASIN(xyz(3)/r)
    alam = ATAN2(xyz(1),xyz(2))

  END SUBROUTINE xyzpol

!-------------------------------------------------------------------------------

  SUBROUTINE polxyz(r,alam,phi,xyz)

    REAL(dp), INTENT(IN)   :: r, alam, phi
    REAL(dp), INTENT(OUT)  :: xyz(3)
    REAL(dp)               :: zc

    xyz(3) = r * SIN(phi)
    zc = r * COS(phi)
    xyz(1) = zc * COS(alam)
    xyz(2) = zc * SIN(alam)

  END SUBROUTINE polxyz

!-------------------------------------------------------------------------------

  SUBROUTINE geo2dek(vc,vg,p)

  ! convert geographical vector components to cartesian

    REAL(dp), INTENT(OUT) :: vc(3)     ! x,y,z vector components
    REAL(dp), INTENT(IN)  :: vg(3)     ! lam, phi, r vector components
    REAL(dp), INTENT(IN)  :: p(2)      ! lam, phi (in radians)
    REAL(dp)              :: lam, phi

    lam=p(1)
    phi=p(2)
    vc(1)= vg(3) * COS(phi) * COS(lam) - vg(1) * SIN(lam) &
         - vg(2) * SIN(phi) * COS(lam)
    vc(2)= vg(3) * COS(phi) * SIN(lam) + vg(1) * COS(lam) &
         - vg(2) * SIN(phi) * SIN(lam)
    vc(3)= vg(3) * SIN(phi) + vg(2) * COS(phi)

  END SUBROUTINE geo2dek

!-------------------------------------------------------------------------------

  SUBROUTINE dek2geo(vc,vg,p)

  ! convert cartesian vector components to geographical

    REAL(dp), INTENT(IN)    :: vc(3)     ! x,y,z vector components
    REAL(dp), INTENT(OUT)   :: vg(3)     ! lam, phi, r vector components
    REAL(dp), INTENT(IN)    :: p(2)      ! lam, phi (in radians)
    REAL(dp)                :: lam, phi

    lam = p(1)
    phi = p(2)
    vg(1) =-vc(1) * SIN(lam) + vc(2) * COS(lam)
    vg(2) = vc(3) * COS(phi) - vc(1) * SIN(phi) * COS(lam) &
          - vc(2) * SIN(phi) * SIN(lam)
    vg(3) = vc(1) * COS(phi) * COS(lam) + vc(2) * COS(phi) * SIN(lam) &
          + vc(3) * SIN(phi)

  END SUBROUTINE dek2geo

END MODULE mo_rotation

