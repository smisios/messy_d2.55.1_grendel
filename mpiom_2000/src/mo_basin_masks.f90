! handles ocean basin diagnostics
MODULE mo_basin_masks
  USE mo_kind, ONLY : dp, wp
  USE mo_param1, ONLY: ie, je, ie_g, je_g
  USE mo_commo1, ONLY: iwetol1_g, lwetol1_g, lbounds_exch_tp, &
       di, dz, dzw, tiestu
  USE mo_mpi, ONLY: p_pe, p_io
  USE mo_parallel, ONLY: p_ioff, p_joff, p_bcast
  USE mo_units, ONLY: io_stdout
  USE mo_io_config, ONLY: next_free_unit
  IMPLICIT NONE
  INTEGER, ALLOCATABLE :: ibek(:,:), ibek_g(:,:)
  REAL(wp), ALLOCATABLE,TARGET :: rbek(:,:)
CONTAINS
  SUBROUTINE load_basin_masks(istart)
    INTEGER, INTENT(in) :: istart
    CHARACTER*5 CIJJU
    INTEGER ijju, icou, jmend, j, i, jmmm, jmm, jmanf
    INTEGER io_in_bgin
    io_in_bgin = next_free_unit()
    IF (p_pe == p_io) OPEN(IO_IN_BGIN,FILE='BEK',FORM='FORMATTED', &
         ACTION='readwrite')
    !
    ALLOCATE(ibek(ie, je), ibek_g(ie_g, je_g))
    ALLOCATE(rbek(ie, je))
    !
    DO J=1,JE_G
      DO I=1,IE_G
        IBEK_G(I,J) = 9 * iwetol1_g(i, j)
      END DO
    END DO
    JMMM = (JE_G - 1)/120
    IF (ISTART .GT. 0 ) THEN
      IF (p_pe == p_io) THEN
        DO JMM=0,JMMM
          JMANF=1+JMM*120
          JMEND=MIN((JMM+1)*120,JE_G)
          DO I=2,IE_G-1
            IF ( lbounds_exch_tp ) THEN
              READ(IO_IN_BGIN,'(a5,1X,120I1)')CIJJU,                        &
                   &                    (IBEK_G(I,J),J=JMEND,JMANF,-1)
            ELSE
              READ(IO_IN_BGIN,'(I3,2X,120I1)')IJJU,                         &
                   &                    (IBEK_G(I,J),J=JMEND,JMANF,-1)
            END IF
          END DO
        END DO
      END IF
      CALL p_bcast(IBEK_G,p_io)
    END IF !ISTART
!
    WRITE(IO_STDOUT,*)'read IBEK successful'
    DO J=1,JE_G
      IBEK_G(1, J) = IBEK_G(IE_G-1, J)
      IBEK_G(IE_G, J) = IBEK_G(2, J)
    END DO

    DO J=2,JE_G-1
      DO I=2,IE_G-1
        IBEK_G(I,J)=IBEK_G(I,J) * iwetol1_g(i,j)
        IF (lwetol1_g(i,j) .AND. IBEK_G(I,J) .EQ. 0) THEN
          IBEK_G(I,J) = MAX(IBEK_G(I+1,J),IBEK_G(I-1,J),IBEK_G(I,J+1),    &
               IBEK_G(I,J-1))
        END IF
      END DO
    END DO

    IF (p_pe == p_io) THEN
      REWIND(IO_IN_BGIN)
      DO JMM=0,JMMM
        JMANF=1+JMM*120
        JMEND=MIN((JMM+1)*120,JE_G)
        WRITE(IO_STDOUT,*)'IBEK, JM ',JMM,JMEND,JMANF
        DO I=2,IE_G-1
          IF ( lbounds_exch_tp ) THEN
            WRITE(IO_IN_BGIN,'(I5,1X,120I1)')I,                           &
                 (IBEK_G(I,J),J=JMEND,JMANF,-1)
          ELSE
            WRITE(IO_IN_BGIN,'(I3,2X,120I1)')I,                           &
                 (IBEK_G(I,J),J=JMEND,JMANF,-1)
          END IF
        END DO
      END DO
      CLOSE(IO_IN_BGIN)
    END IF

    IBEK(:,:) = IBEK_G(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)

    rbek(:,:)=REAL(ibek(:,:),dp)

    IF (p_pe == p_io) THEN
      WRITE(IO_STDOUT,*) 'DZ ', DZ,DI,DZW,TIESTU
      WRITE(IO_STDOUT,*)'WETO:'
      DO JMM=0,JMMM
        JMANF=1+JMM*120
        JMEND=MIN((JMM+1)*120,JE_G)
        WRITE(IO_STDOUT,*)'JMM ',JMM,JMANF,JMEND
        WRITE(IO_STDOUT,6061)0,(MOD(J,10),J=JMEND,JMANF,-1)
        WRITE(IO_STDOUT,*)'        J  <=== '
        !         DO I=1,IE_G
        !            WRITE(IO_STDOUT,6061)I,(NINT(WETOL1_G(I,J))*MOD(J,10)       &
        !                 -10*NINT(WETOL1_G(I,J)-1.),J=JMEND,JMANF,-1)
        !         END DO
      END DO

      DO J=1,JE_G
        ICOU=0
        DO I=2,IE_G-1
          IF (lwetol1_g(i, j)) ICOU = ICOU + 1
        END DO
        WRITE(IO_STDOUT,*)'ZAEHL : ', J, ICOU
      END DO
6061  FORMAT(I4,1X,120I1)
    END IF
  END SUBROUTINE load_basin_masks
END MODULE mo_basin_masks
