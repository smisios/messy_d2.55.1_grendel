MODULE mo_contro

    USE MO_PARAM1
    USE MO_MPI
    USE MO_PARALLEL
    USE MO_COMMO1
    USE MO_COMMOAU1
    USE mo_planetary_constants, ONLY : rhoicwa, rhosnwa
    USE MO_UNITS

    IMPLICIT NONE

    REAL(wp) :: gl_saltcont_old
    REAL(wp) :: gl_saltcont

  CONTAINS

    SUBROUTINE CONTRO(L)

      INTEGER :: l,i,j,k,jb
      REAL(wp) :: delt_salt,sum1
      REAL(wp) :: sums(IE,JE)


      IF (icontro.NE.0) THEN

        !     Summing up in the following way should give identical results
        !     independent of the number of processors
        !     (at least if optimization is turned off)

        jb=MERGE(3,2,( lbounds_exch_tp .AND. have_g_js ))

        IF (L.EQ.-999) THEN

          sums(:,:) = 0._wp

          DO K=1,KE
            DO J=jb,JE
              DO I=2,IE-1
                sums(i,j)=sums(i,j)+DLXP(I,J)*DLYP(I,J)*DDPO(I,J,K)*SAO(I,J,K)
              ENDDO
            ENDDO
          ENDDO

          DO J=jb,JE
            DO I=2,IE-1
              sums(i,j)=sums(i,j)+DLXP(I,J)*DLYP(I,J)*WETO(I,J,1)*(SAO(I,J,1) &
                   *(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)    &
                   +SICE*SICTHO(I,J)*RHOICWA)
            ENDDO
          ENDDO

          gl_saltcont_old = global_array_sum(sums)

        ELSE

          gl_saltcont_old=gl_saltcont

        ENDIF

        sums(:,:) = 0._wp
        DO K=1,KE
          DO J=jb,JE
            DO I=2,IE-1
              sums(i,j)=sums(i,j)+DLXP(I,J)*DLYP(I,J)*DDPO(I,J,K)*SAO(I,J,K)
            ENDDO
          ENDDO
        ENDDO

        DO J=jb,JE
          DO I=2,IE-1
            sums(i,j)=sums(i,j)+DLXP(I,J)*DLYP(I,J)*WETO(I,J,1)*(SAO(I,J,1) &
                 *(ZO(I,J)-SICTHO(I,J)*RHOICWA-SICSNO(I,J)*RHOSNWA)    &
                 +SICE*SICTHO(I,J)*RHOICWA)
          ENDDO
        ENDDO

        gl_saltcont = global_array_sum(sums)

        delt_salt = (gl_saltcont / gl_saltcont_old) - 1._wp

        IF (p_pe==p_io) THEN
          WRITE(IO_STDOUT,*)'SALZCHECK:',L,gl_saltcont,delt_salt
          WRITE(0,*)'SALZCHECK:',L,gl_saltcont,delt_salt
        ENDIF

!!$        sums(:,:) = 0.
!!$
!!$        DO J=jb,JE
!!$          DO I=2,IE-1
!!$            sums(i,j)=sums(i,j)+DLXP(I,J)*DLYP(I,J)*WETO(I,J,1)*SICTHO(I,J)
!!$          ENDDO
!!$        ENDDO
!!$
!!$        SUM1 = global_array_sum(sums)
!!$
!!$
!!$        IF (p_pe==p_io) THEN
!!$          WRITE(IO_STDOUT,*)'ICECHECK: ',SUM1,L
!!$          WRITE(0,*)'ICECHECK: ',SUM1,L
!!$        ENDIF


        IF (icontro.EQ.999) THEN
          CALL WRTE_DEBUG(L)
        ENDIF

      ENDIF

    END SUBROUTINE CONTRO

  END MODULE mo_contro
