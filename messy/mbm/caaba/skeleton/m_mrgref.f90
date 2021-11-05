! slightly modified version of mrgref.f90 from ORDERPACK 2.0 at
! http://www.fortran-2000.com/rank/
! The web page says:
! "Users can freely download ORDERPACK 2.0 from this site."

MODULE m_mrgref
  USE messy_mecca_kpp_precision, ONLY: DP

  PUBLIC :: mrgref
  PRIVATE :: R_mrgref, I_mrgref, D_mrgref
  INTERFACE mrgref
    MODULE PROCEDURE d_mrgref, r_mrgref, i_mrgref
  END INTERFACE mrgref
CONTAINS

  SUBROUTINE D_mrgref (XVALT, IRNGT)
    !   Ranks array XVALT into index array IRNGT, using merge-sort
    ! __________________________________________________________
    !   This version is not optimized for performance, and is thus
    !   not as difficult to read as some other ones.
    !   Michel Olagnon - April 2000
    ! __________________________________________________________
    ! __________________________________________________________
    REAL(DP), DIMENSION (:), INTENT (In) :: XVALT
    INTEGER, DIMENSION (:), INTENT (Out) :: IRNGT
    ! __________________________________________________________
    !
    INTEGER, DIMENSION (:), ALLOCATABLE :: JWRKT
    INTEGER :: LMTNA, LMTNC
    INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    !
    NVAL = MIN (SIZE(XVALT), SIZE(IRNGT))
    IF (NVAL <= 0) THEN
      RETURN
    END IF
    !
    !  Fill-in the index array, creating ordered couples
    !
    DO IIND = 2, NVAL, 2
      IF (XVALT(IIND-1) < XVALT(IIND)) THEN
        IRNGT (IIND-1) = IIND - 1
        IRNGT (IIND) = IIND
      ELSE
        IRNGT (IIND-1) = IIND
        IRNGT (IIND) = IIND - 1
      END IF
    END DO
    IF (MODULO (NVAL, 2) /= 0) THEN
      IRNGT (NVAL) = NVAL
    END IF
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    ALLOCATE (JWRKT(1:NVAL))
    LMTNC = 2
    LMTNA = 2
    !
    !  Iteration. Each time, the length of the ordered subsets
    !  is doubled.
    !
    DO
      IF (LMTNA >= NVAL) EXIT
      IWRKF = 0
      LMTNC = 2 * LMTNC
      IWRK = 0
      !
      !   Loop on merges of A and B into C
      !
      DO
        IINDA = IWRKF
        IWRKD = IWRKF + 1
        IWRKF = IINDA + LMTNC
        JINDA = IINDA + LMTNA
        IF (IWRKF >= NVAL) THEN
          IF (JINDA >= NVAL) EXIT
          IWRKF = NVAL
        END IF
        IINDB = JINDA
        !
        !   Shortcut for the case when the max of A is smaller
        !   than the min of B (no need to do anything)
        !
        IF (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) THEN
          IWRK = IWRKF
          CYCLE
        END IF
        !
        !  One steps in the C subset, that we create in the final rank array
        !
        DO
          IF (IWRK >= IWRKF) THEN
            !
            !  Make a copy of the rank array for next iteration
            !
            IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
            EXIT
          END IF
          !
          IWRK = IWRK + 1
          !
          !  We still have unprocessed values in both A and B
          !
          IF (IINDA < JINDA) THEN
            IF (IINDB < IWRKF) THEN
              IF (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
                & THEN
                IINDB = IINDB + 1
                JWRKT (IWRK) = IRNGT (IINDB)
              ELSE
                IINDA = IINDA + 1
                JWRKT (IWRK) = IRNGT (IINDA)
              END IF
            ELSE
              !
              !  Only A still with unprocessed values
              !
              IINDA = IINDA + 1
              JWRKT (IWRK) = IRNGT (IINDA)
            END IF
          ELSE
            !
            !  Only B still with unprocessed values
            !
            IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
            IWRK = IWRKF
            EXIT
          END IF
          !
        END DO
      END DO
      !
      !  The Cs become As and Bs
      !
      LMTNA = 2 * LMTNA
    END DO
    !
    !  Clean up
    !
    DEALLOCATE (JWRKT)
    RETURN
    !
  END SUBROUTINE D_mrgref

  SUBROUTINE R_mrgref (XVALT, IRNGT)
    !   Ranks array XVALT into index array IRNGT, using merge-sort
    ! __________________________________________________________
    !   This version is not optimized for performance, and is thus
    !   not as difficult to read as some other ones.
    !   Michel Olagnon - April 2000
    ! __________________________________________________________
    ! _________________________________________________________
    REAL, DIMENSION (:), INTENT (In) :: XVALT
    INTEGER, DIMENSION (:), INTENT (Out) :: IRNGT
    ! __________________________________________________________
    !
    INTEGER, DIMENSION (:), ALLOCATABLE :: JWRKT
    INTEGER :: LMTNA, LMTNC
    INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    !
    NVAL = MIN (SIZE(XVALT), SIZE(IRNGT))
    IF (NVAL <= 0) THEN
      RETURN
    END IF
    !
    !  Fill-in the index array, creating ordered couples
    !
    DO IIND = 2, NVAL, 2
      IF (XVALT(IIND-1) < XVALT(IIND)) THEN
        IRNGT (IIND-1) = IIND - 1
        IRNGT (IIND) = IIND
      ELSE
        IRNGT (IIND-1) = IIND
        IRNGT (IIND) = IIND - 1
      END IF
    END DO
    IF (MODULO (NVAL, 2) /= 0) THEN
      IRNGT (NVAL) = NVAL
    END IF
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    ALLOCATE (JWRKT(1:NVAL))
    LMTNC = 2
    LMTNA = 2
    !
    !  Iteration. Each time, the length of the ordered subsets
    !  is doubled.
    !
    DO
      IF (LMTNA >= NVAL) EXIT
      IWRKF = 0
      LMTNC = 2 * LMTNC
      IWRK = 0
      !
      !   Loop on merges of A and B into C
      !
      DO
        IINDA = IWRKF
        IWRKD = IWRKF + 1
        IWRKF = IINDA + LMTNC
        JINDA = IINDA + LMTNA
        IF (IWRKF >= NVAL) THEN
          IF (JINDA >= NVAL) EXIT
          IWRKF = NVAL
        END IF
        IINDB = JINDA
        !
        !   Shortcut for the case when the max of A is smaller
        !   than the min of B (no need to do anything)
        !
        IF (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) THEN
          IWRK = IWRKF
          CYCLE
        END IF
        !
        !  One steps in the C subset, that we create in the final rank array
        !
        DO
          IF (IWRK >= IWRKF) THEN
            !
            !  Make a copy of the rank array for next iteration
            !
            IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
            EXIT
          END IF
          !
          IWRK = IWRK + 1
          !
          !  We still have unprocessed values in both A and B
          !
          IF (IINDA < JINDA) THEN
            IF (IINDB < IWRKF) THEN
              IF (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
                & THEN
                IINDB = IINDB + 1
                JWRKT (IWRK) = IRNGT (IINDB)
              ELSE
                IINDA = IINDA + 1
                JWRKT (IWRK) = IRNGT (IINDA)
              END IF
            ELSE
              !
              !  Only A still with unprocessed values
              !
              IINDA = IINDA + 1
              JWRKT (IWRK) = IRNGT (IINDA)
            END IF
          ELSE
            !
            !  Only B still with unprocessed values
            !
            IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
            IWRK = IWRKF
            EXIT
          END IF
          !
        END DO
      END DO
      !
      !  The Cs become As and Bs
      !
      LMTNA = 2 * LMTNA
    END DO
    !
    !  Clean up
    !
    DEALLOCATE (JWRKT)
    RETURN
    !
  END SUBROUTINE R_mrgref
  SUBROUTINE I_mrgref (XVALT, IRNGT)
    !   Ranks array XVALT into index array IRNGT, using merge-sort
    ! __________________________________________________________
    !   This version is not optimized for performance, and is thus
    !   not as difficult to read as some other ones.
    !   Michel Olagnon - April 2000
    ! __________________________________________________________
    ! __________________________________________________________
    INTEGER, DIMENSION (:), INTENT (In)  :: XVALT
    INTEGER, DIMENSION (:), INTENT (Out) :: IRNGT
    ! __________________________________________________________
    !
    INTEGER, DIMENSION (:), ALLOCATABLE :: JWRKT
    INTEGER :: LMTNA, LMTNC
    INTEGER :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
    !
    NVAL = MIN (SIZE(XVALT), SIZE(IRNGT))
    IF (NVAL <= 0) THEN
      RETURN
    END IF
    !
    !  Fill-in the index array, creating ordered couples
    !
    DO IIND = 2, NVAL, 2
      IF (XVALT(IIND-1) < XVALT(IIND)) THEN
        IRNGT (IIND-1) = IIND - 1
        IRNGT (IIND) = IIND
      ELSE
        IRNGT (IIND-1) = IIND
        IRNGT (IIND) = IIND - 1
      END IF
    END DO
    IF (MODULO (NVAL, 2) /= 0) THEN
      IRNGT (NVAL) = NVAL
    END IF
    !
    !  We will now have ordered subsets A - B - A - B - ...
    !  and merge A and B couples into     C   -   C   - ...
    !
    ALLOCATE (JWRKT(1:NVAL))
    LMTNC = 2
    LMTNA = 2
    !
    !  Iteration. Each time, the length of the ordered subsets
    !  is doubled.
    !
    DO
      IF (LMTNA >= NVAL) EXIT
      IWRKF = 0
      LMTNC = 2 * LMTNC
      IWRK = 0
      !
      !   Loop on merges of A and B into C
      !
      DO
        IINDA = IWRKF
        IWRKD = IWRKF + 1
        IWRKF = IINDA + LMTNC
        JINDA = IINDA + LMTNA
        IF (IWRKF >= NVAL) THEN
          IF (JINDA >= NVAL) EXIT
          IWRKF = NVAL
        END IF
        IINDB = JINDA
        !
        !   Shortcut for the case when the max of A is smaller
        !   than the min of B (no need to do anything)
        !
        IF (XVALT(IRNGT(JINDA)) <= XVALT(IRNGT(JINDA+1))) THEN
          IWRK = IWRKF
          CYCLE
        END IF
        !
        !  One steps in the C subset, that we create in the final rank array
        !
        DO
          IF (IWRK >= IWRKF) THEN
            !
            !  Make a copy of the rank array for next iteration
            !
            IRNGT (IWRKD:IWRKF) = JWRKT (IWRKD:IWRKF)
            EXIT
          END IF
          !
          IWRK = IWRK + 1
          !
          !  We still have unprocessed values in both A and B
          !
          IF (IINDA < JINDA) THEN
            IF (IINDB < IWRKF) THEN
              IF (XVALT(IRNGT(IINDA+1)) > XVALT(IRNGT(IINDB+1))) &
                & THEN
                IINDB = IINDB + 1
                JWRKT (IWRK) = IRNGT (IINDB)
              ELSE
                IINDA = IINDA + 1
                JWRKT (IWRK) = IRNGT (IINDA)
              END IF
            ELSE
              !
              !  Only A still with unprocessed values
              !
              IINDA = IINDA + 1
              JWRKT (IWRK) = IRNGT (IINDA)
            END IF
          ELSE
            !
            !  Only B still with unprocessed values
            !
            IRNGT (IWRKD:IINDB) = JWRKT (IWRKD:IINDB)
            IWRK = IWRKF
            EXIT
          END IF
          !
        END DO
      END DO
      !
      !  The Cs become As and Bs
      !
      LMTNA = 2 * LMTNA
    END DO
    !
    !  Clean up
    !
    DEALLOCATE (JWRKT)
    RETURN
    !
  END SUBROUTINE I_mrgref
END MODULE m_mrgref
