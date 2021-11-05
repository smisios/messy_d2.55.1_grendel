  SUBROUTINE beleg
!> subroutine beleg
!> computes or checks values of some layer depth configuration fields
!> initialize some model fields

!----------------------------------------------------------------
!
!     BBBBB   EEEEEE  L       EEEEEE  GGGGG
!     B    B  E       L       E       G
!     BBBBB   EEEEE   L       EEEEE   GGGGG
!     B    B  E       L       E       G    G
!     BBBBB   EEEEEE  LLLLLL  EEEEEE  GGGGGG
!
!
!----------------------------------------------------------------

      USE mo_kind, ONLY: wp
      USE mo_param1
      USE mo_mpi
      USE mo_parallel
      USE mo_commo1
      USE mo_units
      USE mo_planetary_constants, ONLY : rhoref_water, slpref

      IMPLICIT NONE

      REAL(wp) :: rho
      INTEGER  :: i, j, k

!
!> Layer depth configuration (layer thicknesses dzw starting from the top
!> are read in via namelist ocedzw in main program mpiom) 
!>
!>    tiestu(k) : depth of u-point in layer k
!>    tiestw(k) : depth of w-point of the upper boundary of layer k

      tiestu(1) = 0.5_wp * dzw(1)
      tiestw(1) = 0.0_wp

      DO k=1,ke
        dwi(k)       = 1._wp / dzw(k)
        tiestw(k+1)  = tiestw(k) + dzw(k)
      END DO
!
!> Calculation of vector point distances dz(1,..,kep)
!
      DO k=2,ke
        tiestu(k) = 0.5_wp * ( tiestw(k+1) + tiestw(k) )
      END DO
      tiestu(kep) = 9000._wp

      dz(1) = tiestu(1)
      DO k=2,kep
        dz(k) = tiestu(k) - tiestu(k-1)
      END DO

      DO k=1,ke
        WRITE(io_stdout,4546) k, tiestw(k), dzw(k)
        WRITE(io_stdout,4547) k, tiestu(k), dz(k)
      END DO

      WRITE(io_stdout,4546) kep, tiestw(k), zero
      WRITE(io_stdout,4547) kep, tiestu(k), dz(k)

4546  FORMAT(' LAYER ',i2,' W-POINT DEPTH ',f6.0,30x,' THICKNESS : ',   &
     &f6.1)
4547  FORMAT(' LAYER ',i2,' U-POINT DEPTH ',15x,f6.0,' DISTANCE  : ',   &
     &f6.1)

!
!  End of layer depth configuration
!

      WRITE(io_stdout,*)' REFERENCE STRATIFICATION : '

      preffw(:) = 0.0001_wp * rhoref_water * tiestw(:)
      preff(1:ke) = 0.0001_wp * rhoref_water * tiestu(1:ke)

      DO k=1,ke

        di(k)    = 1._wp / dz(k)

        WRITE(io_stdout,'(a,i3,a,f12.4,a,3f10.5)') &
              ' LAYER ', k,' RHO : ', rho(saf(k),taf(k),preff(k)), &
              ' S,T,P ', saf(k), taf(k), preff(k)

        WRITE(io_stdout,*) ie_g, je_g, ke

        DO j=1,je
          DO i=1,ie
            stabio(i,j,k)=zero

            ddue(i,j,k)=zero
            dduo(i,j,k)=zero

            ddpo(i,j,k)=zero
            dpio(i,j,k)=zero

            weto(i,j,k)  = zero
            amsue(i,j,k) = zero
            amsuo(i,j,k) = zero

            po(i,j,k)=zero

            dvo(i,j,k)=zero

            vke(i,j,k) = zero
            uko(i,j,k) = zero

            uk1o(i,j,k)=zero
            vk1e(i,j,k)=zero

            wo(i,j,k)=zero

            t1o(i,j,k)=zero

            s1o(i,j,k)=zero

!> Temperature and salinity fields  set to reference stratification
!>   to start from a different stratification specify
!>   SBR status0 and CALL status0
            tho(i,j,k)=taf(k)
            sao(i,j,k)=saf(k)

          END DO
        END DO
      END DO

      DO i=1,ilt
        b(i)=zero
        x(i)=zero
      END DO
!
! For layer kep
!
      di(kep) = 1._wp / dz(ke)

      DO j=1,je
        DO i=1,ie
          wo(i,j,kep)=zero
          dvo(i,j,kep)=zero
        END DO
      END DO
!
      DO j=1,je
        DO i=1,ie
          dvo(i,j,1)=zero
 
          eminpo(i,j)=zero

          v1e(i,j)=zero
          u1o(i,j)=zero

          z1o(i,j)=zero
          zo(i,j)=zero
        END DO
      END DO

      fslp(:,:)=slpref

  END SUBROUTINE beleg
