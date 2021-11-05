SUBROUTINE octdiff_trf(trf)


  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_levitus
  USE mo_commoau1
  USE mo_units
  USE mo_octdiff
  USE mo_mpi


  !     COMPUTES DIFFUSION OF TEMPERATURE trf AND SALINITY SAO
  !
  !     UWE MIKOLAJEWICZ 12/99
  !
  !     VERTICAL DIFFUSION IMPLICITE
  !     HORIZONTAL DIFFUSION BOTH HARMONIC AND BIHARMONIC, BOTH EXPLICITE
  !                                AH00         AULAPTS
  !     ACTUAL COEFFICIENTS SPATIALLY VARIABLE WITH RESOLUTION
  !                                AHX,AHY      AULX,AULY
  !
  !     VARIABLE THICKNESS OF SURFACE LAYER NOW INCLUDED
  !
  !     Changes R. Johanni, 2003-11-10:
  !     DSLOPX, DSLOPY, DVSLOP were calculated but nowhere used -> removed
  !
  !     Changes R. Smith, 2004-09-27
  !     tracer-independent matrices calculated separately first for use
  !     with multiple tracers

  IMPLICIT NONE

  INTEGER i, j, k, l
  REAL fakres

  REAL zzsurf(ke)
  REAL tv_diff, sv_diff, dlxyp, dlxyp1
  REAL tflux(ie,je,ke), tfluz, tfluy(ie,je,ke)
  REAL trf(ie,je,ke), ten(ie,je,ke)

!#ifdef bounds_exch_save
  CALL bounds_exch('p',trf,'octdiff_trf 1')
!#endif

  DO k=1,ke
     zzsurf(k)=0.
  ENDDO
  zzsurf(1)=1.

!$OMP PARALLEL PRIVATE(i,j,k,l,tv_diff,sv_diff,dlxyp,dlxyp1,tfluz)

  DO k = 1, ke
!$OMP DO
     DO j = 1, je
        DO i = 1, ie
           t1o(i,j,k)=trf(i,j,k)*vol_term(i,j,k)
        END DO
     END DO
!$OMP END DO
  END DO

  IF(ah00.GT.almzer)THEN


     DO k = 1, ke
!$OMP DO
        DO j = 1, je
           DO i = 1, ie
              tflux(i,j,k)=0.0
           END DO
        END DO
!$OMP END DO
     END DO

     !X DIRECTION - tracer dependent

#ifdef ISOPYK

     DO k=1,ke
!$OMP DO
        DO j=1,je
           DO i=1,ie1
              IF(amsuo(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i+1,j,k)-trf(i,j,k))/dlxu(i,j)
                 sv_diff=(sao(i+1,j,k)-sao(i,j,k))/dlxu(i,j)
                 dlxyp=dlxp(i,j)*dlyp(i,j)
                 dlxyp1=dlxp(i+1,j)*dlyp(i+1,j)
                 tfluz=0
                 !
                 !       triangle left,upw
                 !
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_lo(i,j,k)*(tv_diff  &
                      +sloplox(i,j,k)*                                &
                      (trf(i,j,ko(j,k))-trf(i,j,k))/dz(k))
                 tfluz=zcoeff_lox(i,j,k)*(sloplox(i,j,k)*             &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*                  &
                      dlxyp/dz(k)+1.*dlxyp*tv_diff)

                 t1o(i,j,ko(j,k))=t1o(i,j,ko(j,k))-tfluz
                 t1o(i,j,k)=t1o(i,j,k)+tfluz
                 !
                 !       triangle left down
                 !
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_lu(i,j,k)*(tv_diff  &
                      +sloplux(i,j,k)*                                &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))/dz(ku(j,k)))            
                 tfluz=zcoeff_lux(i,j,k)*(sloplux(i,j,k)*             &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*                  &
                      dlxyp/dz(ku(j,k))+1.*dlxyp*tv_diff)

                 t1o(i,j,k)=t1o(i,j,k)-tfluz
                 t1o(i,j,ku(j,k))=t1o(i,j,ku(j,k))+tfluz

              ENDIF
           END DO
        END DO
!$OMP END DO
     END DO


     DO k=1,ke
!$OMP DO
        DO j=1,je
           DO i=1,ie1
              IF(amsuo(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i+1,j,k)-trf(i,j,k))/dlxu(i,j)
                 dlxyp=dlxp(i,j)*dlyp(i,j)
                 dlxyp1=dlxp(i+1,j)*dlyp(i+1,j)
                 tfluz=0

                 !
                 !       triangle right,upw
                 !
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_ro(i,j,k)*(tv_diff &
                      +sloprox(i,j,k)*                               &
                      (trf(i+1,j,ko(j,k))-trf(i+1,j,k))/dz(k))
                 tfluz=zcoeff_rox(i,j,k)*(sloprox(i,j,k)*            &
                      (trf(i+1,j,ko(j,k))-trf(i+1,j,k))*             &
                      dlxyp1/dz(k)+1.*dlxyp1*tv_diff)

                 t1o(i+1,j,ko(j,k))=t1o(i+1,j,ko(j,k))-tfluz
                 t1o(i+1,j,k)=t1o(i+1,j,k)+tfluz

                 !
                 !       triangle right down
                 !
                 tflux(i,j,k)=tflux(i,j,k)+xcoeff_ru(i,j,k)*(tv_diff &
                      +sloprux(i,j,k)*                               &
                      (trf(i+1,j,k)-trf(i+1,j,ku(j,k)))/dz(ku(j,k)))
                 tfluz=zcoeff_rux(i,j,k)*(sloprux(i,j,k)*            &
                      (trf(i+1,j,k)-trf(i+1,j,ku(j,k)))*             &
                      dlxyp1/dz(ku(j,k))+1.*dlxyp1*tv_diff)

                 t1o(i+1,j,k)=t1o(i+1,j,k)-tfluz
                 t1o(i+1,j,ku(j,k))=t1o(i+1,j,ku(j,k))+tfluz

              ENDIF
           END DO
        END DO
!$OMP END DO
     END DO

#else

     DO k=1,ke
!$OMP DO
        DO j=1,je
           DO i=1,ie1
              IF(amsuo(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i+1,j,k)-trf(i,j,k))/dlxu(i,j)
                 tflux(i,j,k)=tv_diff*xflux(i,j,k)
              ENDIF
           END DO
        END DO
!$OMP END DO
     END DO

#endif
!$OMP SINGLE
     call bounds_exch('u',tflux,'octdiff_trf 2 ')
!$OMP END SINGLE

!!$!reformulated to avoid bank conflicts on sx6 
!!$!!CDIR NOLOOPCHG
!!$     DO k = 1, ke
!!$!$OMP DO
!!$        DO j = 1, je
!!$           DO i = 1, ie1
!!$              t1o(i,j,k)   = t1o(i,j,k)   + tflux(i,j,k)
!!$              t1o(i+1,j,k) = t1o(i+1,j,k) - tflux(i,j,k)
!!$           END DO
!!$        END DO
!!$!$OMP END DO
!!$     END DO

!$OMP DO
     DO k = 1, ke
        DO j = 1, je
!hh not needed due to following bounds_exch
!              t1o(1,j,k)   = t1o(1,j,k)   + tflux(1,j,k)
!              t1o(ie,j,k)  = t1o(ie,j,k)  - tflux(ie1,j,k)
           DO i = 2, ie1
              t1o(i,j,k)   = t1o(i,j,k)   + tflux(i,j,k) - tflux(i-1,j,k)
           END DO
        END DO
     END DO 
!$OMP END DO

!$OMP SINGLE
     call bounds_exch('p',t1o,'octdiff_trf 3 ')
!$OMP END SINGLE

!$OMP DO
     DO k = 1, ke
        DO j = 1, je
           DO i = 1, ie
              tfluy(i,j,k) = 0.0
           END DO
        END DO
     END DO
!$OMP END DO

     !Y DIRECTION - tracer dependent values

#ifdef ISOPYK

     DO k=1,ke
!$OMP DO
        DO j=1,je1
           DO i=1,ie
              IF(amsue(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))/dlyv(i,j)
                 dlxyp=dlxp(i,j)*dlyp(i,j)
                 dlxyp1=dlxp(i,j+1)*dlyp(i,j+1)
                 tfluz=0
                 !
                 !       triangle left,upw
                 !
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_lo(i,j,k)*(tv_diff &
                      +sloploy(i,j,k)*                               &
                      (trf(i,j,ko(j,k))-trf(i,j,k))/dz(k))
                 tfluz=zcoeff_loy(i,j,k)*(sloploy(i,j,k)*            &
                      (trf(i,j,ko(j,k))-trf(i,j,k))*                 &
                      dlxyp/dz(k)+1.*dlxyp*tv_diff)

                 t1o(i,j,ko(j,k))=t1o(i,j,ko(j,k))-tfluz
                 t1o(i,j,k)=t1o(i,j,k)+tfluz
                 !
                 !       triangle left down
                 !
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_lu(i,j,k)*(tv_diff &
                      +slopluy(i,j,k)*                               &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))/dz(ku(j,k)))
                 tfluz=zcoeff_luy(i,j,k)*(slopluy(i,j,k)*            &
                      (trf(i,j,k)-trf(i,j,ku(j,k)))*                 &
                      dlxyp/dz(ku(j,k))+1.*dlxyp*tv_diff)

                 t1o(i,j,k)=t1o(i,j,k)-tfluz
                 t1o(i,j,ku(j,k))=t1o(i,j,ku(j,k))+tfluz

              END IF
           END DO
        END DO
!$OMP END DO
     END DO



     DO k=1,ke
!$OMP DO
        DO j=1,je1
           DO i=1,ie
              IF(amsue(i,j,k).GT.0.)THEN
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))/dlyv(i,j)
                 dlxyp=dlxp(i,j)*dlyp(i,j)
                 dlxyp1=dlxp(i,j+1)*dlyp(i,j+1)
                 tfluz=0
                 !
                 !       triangle right,upw
                 !
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_ro(i,j,k)*(tv_diff &
                      +sloproy(i,j,k)*                               &
                      (trf(i,j+1,ko(j,k))-trf(i,j+1,k))/dz(k))
                 tfluz=zcoeff_roy(i,j,k)*(sloproy(i,j,k)*            &
                      (trf(i,j+1,ko(j,k))-trf(i,j+1,k))*             &
                      dlxyp1/dz(k)+1.*dlxyp1*tv_diff)

                 t1o(i,j+1,ko(j,k))=t1o(i,j+1,ko(j,k))-tfluz
                 t1o(i,j+1,k)=t1o(i,j+1,k)+tfluz
                 !
                 !       triangle right down
                 !
                 tfluy(i,j,k)=tfluy(i,j,k)+ycoeff_ru(i,j,k)*(tv_diff &
                      +slopruy(i,j,k)*                               &
                      (trf(i,j+1,k)-trf(i,j+1,ku(j,k)))/dz(ku(j,k)))
                 tfluz=zcoeff_ruy(i,j,k)*(slopruy(i,j,k)*            &
                      (trf(i,j+1,k)-trf(i,j+1,ku(j,k)))*             &
                      dlxyp1/dz(ku(j,k))+1.*dlxyp1*tv_diff)

                 t1o(i,j+1,k)=t1o(i,j+1,k)-tfluz
                 t1o(i,j+1,ku(j,k))=t1o(i,j+1,ku(j,k))+tfluz

              END IF
           END DO
        END DO
!$OMP END DO
     END DO

#else
      tfluy=0.
     DO k=1,ke
!$OMP DO
        DO j=1,je1
           DO i=1,ie                
              IF(amsue(i,j,k).GT.0.)THEN        
                 tv_diff=(trf(i,j+1,k)-trf(i,j,k))/dlyv(i,j)
                 tfluy(i,j,k)=yflux(i,j,k)*tv_diff
              ENDIF
           END DO
        END DO
!$OMP END DO
     END DO

#endif

!        SUMMATION
!$OMP SINGLE
     CALL bounds_exch('v',tfluy,'octdiff_trf 4')
!$OMP END SINGLE

     DO k=1,ke
!$OMP DO
        DO j=2,je1
           DO i=1,ie 
              t1o(i,j,k)=t1o(i,j,k)+tfluy(i,j,k)-tfluy(i,j-1,k)
!              t1o(i,j+1,k)=t1o(i,j+1,k)-tfluy(i,j,k)
           END DO
        END DO
!$OMP END DO
     END DO

!$OMP SINGLE
     CALL bounds_exch('p',T1O,'octdiff trf 5')
!$OMP END SINGLE


     DO k=1,ke
!$OMP DO
        DO j=2,je1
           DO i=1,ie
              ten(i,j,k)=t1o(i,j,k)/vol_term(i,j,k)-trf(i,j,k)
           END DO
        END DO
!$OMP END DO
     END DO

  ELSE ! ah00.gt.almzer

    DO k=1,ke
!$OMP DO
      DO j=1,je
        DO i=1,ie
          ten(i,j,k)=0.
        END DO
      END DO
!$OMP END DO
   END DO 

  END IF



  ! VERTICAL DIFFUSION -tracer dependent


   DO k=1,ke
!$OMP DO
     DO j=1,je 
       DO i=1,ie 
         t1o(i,j,k)=trf(i,j,k)
       END DO
     END DO
!$OMP END DO 
   END DO


  !
  !     IMPLICIT VERTICAL DIFFUSION
  !



  DO k=2,ke
!$OMP DO
    DO j=2,je-1 
      DO i=1,ie 
        t1o(i,j,k)=t1o(i,j,k)-tridsy(i,j,k-1,1)*t1o(i,j,k-1)
      END DO
    END DO  
!$OMP END DO  
  END DO



!$OMP DO
  DO j = 2,je-1
    DO i = 1,ie
      trf(i,j,ke)=t1o(i,j,ke)/tridsy(i,j,ke,2)
    END DO
  END DO
!$OMP END DO



  DO k = ke-1, 1, -1      
!$OMP DO
    DO j = 2, je-1
      DO i = 1,ie 
        trf(i,j,k)=(t1o(i,j,k)-tridsy(i,j,k,3)*trf(i,j,k+1))              &
                  /tridsy(i,j,k,2)

      END DO
    END DO
!$OMP END DO
  END DO


!$OMP SINGLE
  CALL bounds_exch('p',trf,'octdiff_trf 6')
!$OMP END SINGLE





  DO k=1,ke
!$OMP DO
    DO j = 2,je1
      DO i = 1,ie
        trf(i,j,k)=trf(i,j,k)+ten(i,j,k)
      END DO
    END DO
!$OMP END DO
  END DO

!$OMP SINGLE
  CALL bounds_exch('p',trf,'octdiff_trf 7')
!$OMP END SINGLE


  !
  ! BIHARMONIC tracer DIFFUSION - tracer dependent
  !
  IF(aulapts.GT.almzer)THEN

    DO k=1,ke
!$OMP DO
      DO j=2,je1
        DO i=2,ie1
          t1o(i,j,k)=weto(i,j,k)*(                                           &
               (weto(i-1,j,k)*(trf(i-1,j,k)-trf(i,j,k))/dlxu(i-1,j)          &
               +weto(i+1,j,k)*(trf(i+1,j,k)-trf(i,j,k))/dlxu(i,j))/dlxp(i,j) &
               +(weto(i,j-1,k)*(trf(i,j-1,k)-trf(i,j,k))/dlyv(i,j-1)         &
               +weto(i,j+1,k)*(trf(i,j+1,k)-trf(i,j,k))/dlyv(i,j))/dlyp(i,j))
        END DO
      END DO
!$OMP END DO
    END DO


!$OMP SINGLE
    CALL bounds_exch('p',t1o,'octdiff_trf 8')
!$OMP END SINGLE


    DO k=1,ke
!$OMP DO
      DO j=2,je1
        DO i=2,ie1
          trf(i,j,k)=trf(i,j,k)                                              &
               - (weto(i,j,k)/((ddpo(i,j,k)+almzer                           &
               +zzsurf(k)*(zo(i,j)-rhoicwa*sictho(i,j)-rhosnwa*sicsno(i,j))) &
               *dlxp(i,j)*dlyp(i,j)))*dzw(k)*(                               &
               weto(i-1,j,k)*aulx(i-1,j)*dlyu(i-1,j)*                        &
               (t1o(i-1,j,k)-t1o(i,j,k))/dlxu(i-1,j)                         &
               +weto(i+1,j,k)*aulx(i,j)*dlyu(i,j)*                           &
               (t1o(i+1,j,k)-t1o(i,j,k))/dlxu(i,j)                           &
               +weto(i,j-1,k)*auly(i,j-1)*dlxv(i,j-1)*                       &
               (t1o(i,j-1,k)-t1o(i,j,k))/dlyv(i,j-1)                         &
               +weto(i,j+1,k)*auly(i,j)*dlxv(i,j)*                           &
               (t1o(i,j+1,k)-t1o(i,j,k))/dlyv(i,j))
        END DO
      END DO
!$OMP END DO
    END DO


!$OMP SINGLE
    CALL bounds_exch('p',trf,'octdiff_trf 9')
!$OMP END SINGLE


  ENDIF ! aulapts > alzmer


!$OMP END PARALLEL


  
END SUBROUTINE octdiff_trf

