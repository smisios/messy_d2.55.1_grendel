      SUBROUTINE OCSCHEP
!
!     HORIZONTAL DIFFUSION OF MOMENTUM
!
!     BY ERNST MAIER-REIMER
!     MODIFIED BY UWE MIKOLAJEWICZ 2/99
!           MAKE AUSF 2D, ADAPT TO SPATIALLY VARIABLE GRID
!           AUS TIME CONSTANT
!
      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      USE mo_boundsexch, ONLY : bounds_exch

      IMPLICIT NONE

      INTEGER i,j,k
      REAL(wp) P11(IE,JE),P22(IE,JE),P12(IE,JE),AUSF(IE,JE),AUSFPSI(IE,JE)
      REAL(wp) DLXXYY, UX, UY, VX, VY

!
!UWE  AUS DIMENSIONLESS MOMENTUM DIFFUSION COEFFICIENT
!
!     AUSF     2D-MOMENTUM DIFFUSION COEFFICENT P-POINTS      [M2/S]
!     AUSFPSI  THE SAME FOR PSI-POINTS
!
!     SPATIAL DEPENDENCE ACCORDING TO BRYAN ET AL. 1975
!
      DO J=1,JE
        DO I=1,IE
          p11(i, j) = 0._wp
          p12(i, j) = 0._wp
          p22(i, j) = 0._wp
          ausfpsi(i, j) = 0._wp
          DLXXYY=MAX(DLXP(I,J),DLYP(I,J))
          AUSF(I,J) = AUS * DLXXYY**2
        END DO
      END DO
!
      DO J=2,JE-1
       DO I=2,IE-1
        ausfpsi(i, j) = 0.25_wp * (ausf(i, j) + ausf(i+1, j)            &
     &                     +AUSF(I,J+1)+AUSF(I+1,J+1))
       ENDDO
      ENDDO
!
      CALL bounds_exch(1,'s',AUSFPSI,'ocschep 1')
!
      DO K=1,KE
        DO J=2,JE1
          DO I=2,IE1
! ZUSATZTERM KANN EVTL. AUF 0 GESETZT WERDEN, DURCHSTROEME!!!!
!
#ifdef SCHEP
            UK1O(I,J,K)=UKO(I,J,K)*AMSUO(I,J,K)                         &
                 &      -2._wp * (uko(i, j-1, k) + uko(i, j+1, k))      &
                 &        * (1._wp - amsuo(i, j, k))
#else
            UK1O(I,J,K)=UKO(I,J,K)*AMSUO(I,J,K)
#endif /*SCHEP*/
#ifdef SCHEP
            VK1E(I,J,K)=VKE(I,J,K)*AMSUE(I,J,K)                         &
                 &      -2._wp * (vke(i-1, j, k) + vke(i+1, j, k)) &
                 &        * (1._wp - amsue(i, j, k))
#else
            VK1E(I,J,K)=VKE(I,J,K)*AMSUE(I,J,K)
#endif /*SCHEP*/
          END DO
        END DO
      END DO
!
      CALL bounds_exch(1,'u',UK1O,'ocschep 2')
      CALL bounds_exch(1,'v',VK1E,'ocschep 3')
!
      DO K=1,KE
        DO J=2,JE1
          DO I=2,IE1
            UX=(DLYU(I,J)*DDUO(I,J,K)*UKO(I,J,K)                              &
     &     -DLYU(I-1,J)*DDUO(I-1,J,K)*UKO(I-1,J,K))                     &
     &      *DPIO(I,J,K)/(DLXP(I,J)*DLYP(I,J))
!     DV/DY AUF P-PUNKT
            VY=(DDUE(I,J-1,K)*DLXV(I,J-1)*VKE(I,J-1,K)                        &
     &    -DDUE(I,J,K)*DLXV(I,J)*VKE(I,J,K))                            &
     &     *DPIO(I,J,K)/(DLXP(I,J)*DLYP(I,J))

!     DV/DX AUF PSI-PUNKT
            vx = (vk1e(i+1, j, k) - vk1e(i, j, k)) &
                 * 2._wp / (dlxu(i, j) + dlxu(i, j+1))
!     DU/DY AUF PSI-PUNKT
            uy = (uk1o(i, j, k) - uk1o(i, j+1, k)) &
                 * 2._wp / (dlxv(i, j) + dlxv(i+1, j))
            P11(I,J)=AUSF(I,J)*DDPO(I,J,K)*(UX-VY)
            P22(I,J)=-P11(I,J)
            P12(I,J)=AUSFPSI(I,J)*DDPSIO(I,J,K)*(UY+VX)
          END DO
        END DO


      CALL bounds_exch(1,'p',P11,'ocschep 4')
      CALL bounds_exch(1,'p',P22,'ocschep 5')
      CALL bounds_exch(1,'s',P12,'ocschep 6')
!
      DO J=2,JE1
        DO I=2,IE1
!
!replace dwio
!     UKO(I,J,K)=UOO(I,J,K)+                                            &
!    & DWIO(I,J,K)*(DTDXUO(I,J)*(P11(I+1,J)-P11(I,J))                   &
!    & +(DPYO(I,J-1)*P12(I,J-1)-DPYO(I,J)*P12(I,J)))
          UKO(I,J,K)=UOO(I,J,K)+                                            &
     & (amsuo(i,j,k)/(almzer+dduo(i,j,k)))                              &
     &     *(dtdxuo(i,j)*(p11(i+1,j)-p11(i,j))                          &
     & +(dpyo(i,j-1)*p12(i,j-1)-dpyo(i,j)*p12(i,j)))

!
!replace dwie
!     VKE(I,J,K)=VOE(I,J,K)                                             &
!    & + DWIE(I,J,K)*((DTDXUE(I,J)*P12(I,J)-DTDXUE(I-1,J)*P12(I-1,J))   &
!    & + DPYE(I,J)*(P22(I,J)-P22(I,J+1)))
          VKE(I,J,K)=VOE(I,J,K)                                             &
     & + (amsue(i,j,k)/(almzer+ddue(i,j,k)))                            &
     &           *((dtdxue(i,j)*p12(i,j)-dtdxue(i-1,j)*p12(i-1,j))      &
     & + dpye(i,j)*(p22(i,j)-p22(i,j+1)))

        END DO
      END DO
!
      DO J=2,JE1
        DO I=2,IE1
          UOO(I,J,K)=UKO(I,J,K)
          VOE(I,J,K)=VKE(I,J,K)
        END DO
      END DO

     ENDDO

      CALL bounds_exch(1,'u',UKO,'ocschep 8')
      CALL bounds_exch(1,'v',VKE,'ocschep 9')
      CALL bounds_exch(1,'u',UOO,'ocschep 10')
      CALL bounds_exch(1,'v',VOE,'ocschep 11')

      RETURN
      END
