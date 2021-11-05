!>
!! @ingroup common
!!
MODULE mo_commobbl

  USE mo_kind, ONLY: wp
  USE mo_param1, ONLY: ie, ie1, je, je1, ke
  USE MO_COMMO1, ONLY: kbot, rhoo, dduo, ddue, dlxv, dlyu, uko, vke, almzer, &
       depto, weto, dlxu, dlyv
  USE mo_planetary_constants, ONLY: g
  USE mo_boundsexch, ONLY: bounds_exch
#ifdef DEBUG
  USE mo_mpi, ONLY: p_sum
  USE mo_units, ONLY: io_stdout
#endif
  IMPLICIT NONE

  REAL(wp), ALLOCATABLE:: alpx(:,:),alpy(:,:)                           &
       ,ubbl(:,:),vbbl(:,:)
  ! , bblflx(:,:)
  INTEGER, ALLOCATABLE :: kupubbl(:,:),kdwubbl(:,:)                 &
       ,kupvbbl(:,:),kdwvbbl(:,:)


CONTAINS

  SUBROUTINE alloc_mem_commobbl

    ALLOCATE(alpx(ie,je),alpy(ie,je),                   &
         ubbl(ie,je),vbbl(ie,je),                             &
         kupubbl(ie,je),kdwubbl(ie,je),                       &
         kupvbbl(ie,je),kdwvbbl(ie,je) )
    ! if(!bbl_diag) allocate(bblflx(ie,je))
                alpx(:,:) = 0.0_wp
                alpy(:,:) = 0.0_wp
                !bblflx(:,:) = 0.0_wp
                ubbl(:,:) = 0.0_wp
                vbbl(:,:) = 0.0_wp
                kupubbl(:,:)=0
                kdwubbl(:,:)=0
                kupvbbl(:,:)=0
                kdwvbbl(:,:)=0

  END SUBROUTINE alloc_mem_commobbl

  !>
  !! Calculate sites of possible slope convection prior to time-stepping.
  !!
  !! @author jhj
  !! @date Sept. 2, 1999
  !!
  SUBROUTINE findalfa
    REAL(wp) :: alp
    INTEGER :: i, j

    alpx(:,:) = 0._wp
    alpy(:,:) = 0._wp

    !:: sweep in x dir
    DO j=2,je1
      DO i=2,ie1
        alp= weto(i,j,1)*weto(i+1,j,1)*                                 &
             &        (depto(i+1,j)-depto(i,j))/dlxu(i,j)
        IF (ABS(alp) .GT. 1.e-8_wp) THEN
          alpx(i,j)=alp
        END IF
      END DO
    END DO
#ifdef DEBUG
    WRITE(io_stdout,*)'max number of slope sites x = ', &
         p_sum(COUNT(alpx > 1.e-8_wp))
#endif
    !:: sweep in y dir
    DO j=2,je1
      DO i=2,ie1
        alp= weto(i,j,1)*weto(i,j+1,1)*                                 &
             &        (depto(i,j)-depto(i,j+1))/dlyv(i,j)
        IF (ABS(alp) .GT. 1.e-8_wp) THEN
          alpy(i,j)=alp
        END IF
      END DO
    END DO
#ifdef DEBUG
    WRITE(io_stdout,*)'max number of slope sites y = ', &
         p_sum(COUNT(alpy > 1.e-8_wp))
#endif

    CALL bounds_exch(1,'u',alpx,'findalpha 1')
    CALL bounds_exch(1,'v',alpy,'findalpha 2')

    !jj initialize bbl thickness ans flux fields (diagnostic)
    ! DO j=1,je
    !   DO i=1,ie
    !     bblflx(i,j) = 0.0_wp
    !   END DO
    ! END DO
  END SUBROUTINE findalfa

  !> calculate bbl transport before use in tracer advection
  !! JHJ FEB.24. 2000
  SUBROUTINE slopetrans
    INTEGER :: i,j
    ! LINEAR FRICTION COEFF FOR CAMPIN&GOOSSE FLUX
    ! REAL(wp), PARAMETER :: ymue = 1.e-6_wp, flxfac = g / (rhoref_water * ymue)

    CALL sweep(alpx, kdwubbl, kupubbl, uko, 1, 0, dduo, ubbl, dlyu, 'u')
    CALL sweep(alpy, kdwvbbl, kupvbbl, vke, 0, 1, ddue, vbbl, dlxv, 'v')

    !:: SAVE ABS OF TRANSPORT FOR DIAGNOSTIC


    ! DO J=1,JE
    !    DO I=1,IE
    !       BBLFLX(I,J)=SQRT(UBBL(I,J)**2+VBBL(I,J)**2)*1.E-6_wp
    !    ENDDO
    ! ENDDO
  END SUBROUTINE slopetrans


  SUBROUTINE sweep(alpha, kdw_bbl, kup_bbl, velocity, &
       imask, jmask, level_thickness, bbl_transport, grid_dist, grid_type)
    REAL(wp), INTENT(in) :: alpha(ie, je), velocity(ie, je, ke)
    INTEGER, INTENT(inout) :: kdw_bbl(ie, je), kup_bbl(ie, je)
    INTEGER, INTENT(in) :: imask, jmask ! these are supposed to be 0 or 1
    CHARACTER(len=*), INTENT(in) :: grid_type
    REAL(wp), INTENT(in) :: level_thickness(ie, je, ke), grid_dist(ie, je)
    REAL(wp), INTENT(inout) :: bbl_transport(ie, je)

    REAL(wp) :: bottom_velocity
    ! MAX BBL LAYER THICKNESS
    REAL(wp), PARAMETER :: bblmax = 500._wp
    REAL(wp) :: aux(ie,je)
    REAL(wp) :: alp, ddro1, dzwbbl
    INTEGER :: i, j, ido, jdo, ish, jsh, k, kdw, kup, jivz
    DO j = 2,je-1
      DO i = 2, ie-1

        alp = alpha(i, j)
        !:: SET ZONAL INDEX FOR SHELF (SH) AND DEEP OCEAN (DO)
        jivz = NINT(alp / (almzer + ABS(alp)))

        !:: IF SIGN(SLOPE) IS POSITIVE (I.E. DEP(I+1)> DEP(I))
        !:: THEN THE FLOW IS FROM ISH=I TO IDO=I+1 AND VICE VERSA
        !:: SET MERIDIONAL INDEX FOR SHELF (SH) AND DEEP OCEAN (DO)
        !:: IF SIGN(SLOPE) IS POSITIVE (I.E. DEP(J)> DEP(J+1))
        !:: THEN THE FLOW IS FROM JSH=J+1 TO JDO=J AND VICE VERSA
        ido = i + MAX(0,  jivz) * imask
        ish = i + MAX(0, -jivz) * imask
        jdo = j + MAX(0, -jivz) * jmask
        jsh = j + MAX(0,  jivz) * jmask
        kup = kbot(ish, jsh)
        IF (KUP .LT. 1) THEN
          ddro1 = 0._wp
          bbl_transport(i, j) = 0.0_wp
          kup_bbl(i, j) = 0
          kdw_bbl(i, j) = 0
        ELSE
          ddro1 = rhoo(ish, jsh, kup) - rhoo(ido, jdo, kup)
          !:: CHECK CONDITION FOR DOWNSLOPE CONV
          IF (ddro1 .GT. 1.E-8_wp) THEN
            !:: SEARCH LEVEL OF EQUAL IN-SITU DENSITY
            kdw = kup
            DO WHILE (kdw < kbot(ido,jdo))
              IF (rhoo(ish, jsh, kdw + 1) .LE. rhoo(ido, jdo, kdw + 1)) EXIT
              kdw = kdw + 1
            END DO
            !:: DO SLOPE CONV ONLY IF THERE IS A LESS DENSE DEEPER
            !:: LAYER IN THE NEIGHBORING BOX, OR, IF BOTH CELLS ARE
            !:: IN THE BOTTOMMOST LAYER BUT HAVE DIFFERENT DEPTHS:
            IF (kdw .GT. kup) THEN
              !:: SAVE UP AND DOWNLOCATION
              kup_bbl(i, j) = kup
              kdw_bbl(i, j) = kdw
              !:SET THICKNESS OF DOWNFLOW TO BML THICKNESS OF BBL
              !:: OR TO THE TOTAL THICKNESS
              dzwbbl = MIN(bblmax, level_thickness(i, j, kup))
              !:: USE BOTTOM VELOCITY AT POINT I,J ONLY IF
              !   IT IS DIRECTED DOWNSLOPE
              bottom_velocity = velocity(i,j, kup) * dzwbbl
              IF (alp * bottom_velocity .GT. 0.0_wp) THEN
                bbl_transport(i, j) = bottom_velocity * grid_dist(i, j)
              ELSE
                bbl_transport(i, j) = 0.0_wp
              END IF
            ELSE
              bbl_transport(i, j) = 0.0_wp
              kup_bbl(i, j) = 0
              kdw_bbl(i, j) = 0
            END IF
          ELSE
            bbl_transport(i, j) = 0.0_wp
            kup_bbl(i, j) = 0
            kdw_bbl(i, j) = 0
          END IF
        END IF

      END DO
    END DO
!!$OMP SINGLE
    CALL bounds_exch(1, grid_type, bbl_transport, 'mo_commobbl 1')
!!$OMP END SINGLE
    AUX(:,:) = REAL(kup_bbl(:,:), wp)
!!$OMP SINGLE
    CALL bounds_exch(1, grid_type//'+', aux, 'mo_commobbl 2')
!!$OMP END SINGLE
    kup_bbl(:,:) = INT(AUX(:,:))

    aux(:,:) = REAL(kdw_bbl(:,:), wp)
!!$OMP SINGLE
    CALL bounds_exch(1, grid_type//'+', aux,'mo_commobbl 3')
!!$OMP END SINGLE

    kdw_bbl(:,:) = INT(aux(:,:))
  END SUBROUTINE sweep


END MODULE mo_commobbl
