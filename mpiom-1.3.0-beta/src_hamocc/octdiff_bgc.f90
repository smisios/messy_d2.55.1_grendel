      SUBROUTINE OCTDIFF_BGC
! *****************************************************************************
! 
!     Computes diffusion of ocean tracer
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
!    Stephanie Legutke   MPI, M&D             &9.7.01
!    - no use of dry values in vertical diffusion
!      (no dependence on dry values: checked)
!      DVO was not 0 at interfaces to dry cells; this
!      made the extrapolation of bottom values downward necessary;
!      this is no longer necessary with the new formulation;
!    - no use of dry values in horizontal (isopycnal) diffusion
!      (no dependence on dry values: checked)
!    - loop over imal eliminated for gbc tracer
!    - no use of array T1O in vertical diffusion (results identical)
!    - no use of array THTEN
!    - loop imal moved to main (where OCTDIFF is called)
!
!    P. Wetzel, E Maier-Reimer, MPI-MET     01.08.03
!    - no isopycnal diffusion
!    - TRIDSY(ie,ke,3) --> tridsy(ie,je,ke,3)
!    - computational efficiency & omp
! *****************************************************************************

      USE mo_carbch
      USE mo_control_bgc
      
      USE MO_COMMO1
      USE MO_COMMOAU1
      USE MO_PARAM1_BGC

      use mo_parallel

implicit none

      INTEGER :: i,j,k,l,ll

      REAL   AHX(IE,JE),AHY(IE,JE),AHXD(IE,JE),AHYD(IE,JE)             &
     &      ,ZZSURF(KE)    
     
      REAL :: tflux, tfluy
      REAL :: tridsy(ie,je,ke,3)
      REAL :: TMPVAR(ie)
  

      DO K=1,KE
       ZZSURF(K)=0.
      ENDDO
      ZZSURF(1)=1.

!$OMP parallel DO private (i,j)
      DO J=1,JE
       DO I=1,IE
        AHX(I,J)=AH00*MAX(DLXU(I,J),DLYU(I,J))
        AHY(I,J)=AH00*MAX(DLYU(I,J),DLYV(I,J))
        AHXD(I,J)=AHX(I,J)*DT*DLYU(I,J)
        AHYD(I,J)=AHY(I,J)*DT*DLXV(I,J)
       ENDDO
      ENDDO


! transformation to upperdiagonal (LU1) form only once
!$OMP parallel DO private (i,j,k)
      DO j=1,je
      K=1
      DO  I=1,IE
        TRIDSY(I,J,K,1) = - DT*DVO(I,J,K)*WETO(I,J,K)*DI(K)              &
     &   /(DDPO(I,J,K)+ZZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)            &
     &                      -RHOSNWA*SICSNO(I,J))+ALMZER)
        TRIDSY(I,J,K,3) = - DT*DVO(I,J,K+1)*WETO(I,J,K) * DI(K+1)        &
     &   /(DDPO(I,J,K)+ZZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)            &
     &                        -RHOSNWA*SICSNO(I,J))+ALMZER)
        TRIDSY(I,J,K,2) = 1. - TRIDSY(I,J,K,1) - TRIDSY(I,J,K,3)

      ENDDO

      DO  K=2,KE
      DO  I=1,IE
        TRIDSY(I,J,K,1) = - DT*DVO(I,J,K)  *WETO(I,J,K)*DI(K)  *DPIO(I,J,K)
        TRIDSY(I,J,K,3) = - DT*DVO(I,J,K+1)*WETO(I,J,K)*DI(K+1)*DPIO(I,J,K)
        TRIDSY(I,J,K,2) = 1. - TRIDSY(I,J,K,1) - TRIDSY(I,J,K,3)

      ENDDO
      ENDDO

      DO  K=2,KE
      DO  I=1,IE
        TRIDSY(I,J,K-1,1) = TRIDSY(I,J,K,1) / TRIDSY(I,J,K-1,2)
        TRIDSY(I,J,K,2)   = TRIDSY(I,J,K,2)                              &
     &    - TRIDSY(I,J,K-1,3) * TRIDSY(I,J,K,1) / TRIDSY(I,J,K-1,2)
      ENDDO
      ENDDO
      enddo
! endblock LU 1

!
!  looping over  marine bgc tracer ...
!
!$OMP parallel DO private (i,j,k,l,T1O,TFLUX,TFLUY)
      DO l=1,nocetra

!  ... harmonic (isopycnal/horrizontal) diffusion

      IF(AH00.GT.ALMZER)THEN


!  Cell tracer content accounting for surface elevation and ice draft
!$OMP parallel DO private (i,j,k)
      DO J=1,JE
      DO K=1,KE
      DO I=1,IE
         T1O(I,J,K)=OCETRA(I,J,K,l)*DLXP(I,J)*DLYP(I,J)                &
     &  *(ZZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)-RHOSNWA*SICSNO(I,J))  &
     &                +DDPO(I,J,K))*WETO(I,J,K)
      ENDDO
      ENDDO
      ENDDO

!
!     DIFFUSION IN X-DIRECTION
!

!$OMP parallel DO private (i,j,k,TMPVAR,TFLUX)
      DO J=2,JE-1
      DO K=1,KE
      DO I=1,IE-1
        IF(AMSUO(I,J,K).GT.0.)THEN

        TFLUX=AHXD(I,J)*DDUO(I,J,K)*                                   &
     &            (OCETRA(I+1,J,K,l)-OCETRA(I,J,K,l))/DLXU(I,J)
         T1O(I  ,J,K)=T1O(I  ,J,K)+TFLUX
       TMPVAR(I+1) =TFLUX
         T1O(I+1,J,K)=T1O(I+1,J,K)-TMPVAR(I+1)

        ENDIF
      ENDDO
      ENDDO
      ENDDO

      CALL bounds_exch(T1O)


!
!     DIFFUSION IN Y-DIRECTION
!

! RJ: The following loop may NOT be parallelized with OMP in this form!!!!!
!!!!!$OMP parallel DO private (i,j,k)

        DO J=1,JE1
        DO K=1,KE
        DO I=2,IE1
        IF(AMSUE(I,J,K).GT.0.)THEN

         TFLUY=AHYD(I,J)*DDUE(I,J,K)*                                  &
     &            (OCETRA(I,J+1,K,l)-OCETRA(I,J,K,l))/DLYV(I,J)
         T1O(I,J  ,K)=T1O(I,J  ,K)+TFLUY
         T1O(I,J+1,K)=T1O(I,J+1,K)-TFLUY


         ENDIF

        ENDDO
       ENDDO
      ENDDO

      CALL bounds_exch(T1O)

!$OMP parallel DO private (i,j,k)

      DO J=1,JE
      DO K=1,KE
      DO I=1,IE
         T1O(I,J,K)=WETO(I,J,K)*(T1O(I,J,K)/(DLXP(I,J)*DLYP(I,J)       &
     &   *( ZZSURF(K)*(ZO(I,J)-RHOICWA*SICTHO(I,J)-RHOSNWA*SICSNO(I,J))&
     &     +DDPO(I,J,K)+(1.-WETO(I,J,K)))) - OCETRA(I,J,K,l))
      ENDDO
      ENDDO
      ENDDO


      ENDIF
! End of harmonic diffusion ...


!
!  Vertical diffusion (implicite)
!

!$OMP parallel DO private (i,j,k,ll)
      do j=1,je
      DO  K=2,KE
      DO  I=1,IE
        OCETRA(I,J,K,l)                                                &
     &= OCETRA(I,J,K,l) - TRIDSY(I,J,K-1,1) * OCETRA(I,J,K-1,l)
      ENDDO
      ENDDO

      K=KE
      DO  I=1,IE
        OCETRA(I,J,K,l) = OCETRA(I,J,K,l) / TRIDSY(I,J,K,2)
      ENDDO
      
      DO  K=1,KE1
        Ll=KE-K
        DO  I=1,IE
          OCETRA(I,J,Ll,l)=(OCETRA(I,J,Ll,l)-TRIDSY(I,J,Ll,3)            &
     &                    *OCETRA(I,J,Ll+1,l))/TRIDSY(I,J,Ll,2)
      ENDDO
      ENDDO

      enddo




!:: Update with fluxes from isopycnal/horizontal harmonic diffusion

!$OMP parallel DO private (i,j,k)

       DO j=1,je
        DO k=1,ke
         DO i=1,ie
          OCETRA(I,J,K,l)=OCETRA(I,J,K,l)+T1O(I,J,K)
          ocetra(I,J,K,l)=max(0.,ocetra(I,J,K,l))
         ENDDO
        ENDDO
       ENDDO


!##! $OMP END PARALLEL


      ENDDO  ! End of loop on tracer


      IF( kchck .EQ. 1) THEN 
         WRITE(io_stdo_bgc,*)                                          &
     &   'after octdiff_bgc tracer 4, layer 5 , step ',ldt
         CALL EXTR                                                     &
     &(ie,je,ocetra(1,1,5,4),ddpo(1,1,5),99999.00,io_stdo_bgc)
     
       ENDIF
       
      RETURN

      END
