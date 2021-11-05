      SUBROUTINE compute_mld(kpie,kpje,kpke,ptiestw,pmld)
#ifdef PGBC
      use mo_bgcmean
      USE mo_dynamic

      implicit none
      INTEGER :: i,j,k,kpie,kpje,kpke

      REAL(wp) :: ptiestw(kpke+1)
      REAL(wp) :: pmld(kpie,kpje)

      DO j = 1, kpje
      DO i = 1, kpie
          nbgc_mld(i,j) = 1
      ENDDO
      ENDDO

      DO k = 2, kpke
      DO j = 1, kpje
      DO i = 1, kpie
        if (pmld(i,j).gt.0.5) then
          if (ptiestw(k).le.pmld(i,j)) nbgc_mld(i,j) = k - 1
        endif
      ENDDO
      ENDDO
      ENDDO

      DO j = 1, kpje
      DO i = 1, kpie
!        if (pmld(i,j).gt.0.5) then
!          WRITE(*,*)'mld : ',i,j,pmld(i,j),ptiestw(nbgc_mld(i,j)+1),nbgc_mld(i,j)
!        endif
        bgc_zmld(i,j) = bgc_zmld(i,j) + pmld(i,j)
        bgc_nmld(i,j) = bgc_nmld(i,j) + nbgc_mld(i,j)
      ENDDO
      ENDDO


      RETURN
#endif
      END
