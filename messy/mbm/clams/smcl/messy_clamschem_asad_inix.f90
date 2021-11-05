Module messy_clamschem_asad_inix

contains

!
!     inix  - initialise indice arrays
!
!     Glenn Carver           Centre for Atmospheric Science,
!                            University of Cambridge.
!
!
! Purpose: To set up the arrays holding the indicies of reactions,
!     species etc used in the main parts of the code e.g. prls.
!
!
!          Called from ASAD_CINIT
!
!
!     Interface
!     ---------
!     On exit from this routine, the following arrays will have
!     been set: ngrp, nprdx3, nprdx2, nprdx1, nlall, nlstst, nlf,
!     nldepd, nldepw and nlemit.
!
!     IMPORTANT! This routine MUST be called after inrats since it
!     expects the nspi array to be initialised.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INIX

        USE messy_clamschem_asad_mod
        USE messy_clamschem_asad_dummy, ONLY: ereport
        
        IMPLICIT NONE

!       Local variables

        INTEGER :: j                    ! Loop variable
        INTEGER :: jr                   ! Loop variable
        INTEGER :: js                   ! Loop variable
        INTEGER :: jp                   ! Loop variable
        INTEGER :: iprd
        INTEGER :: igrp3
        INTEGER :: iloss
        INTEGER :: ngmax
        INTEGER :: is
        INTEGER :: ispec
        INTEGER :: jdry
        INTEGER :: jwet
        INTEGER :: jpnpx3

        INTEGER :: igrp(3)
        INTEGER :: idepd(jpspec)
        INTEGER :: idepw(jpspec)

        CHARACTER (LEN=2)  :: itype
        CHARACTER (LEN=72) :: cmessage  ! Error message

        LOGICAL :: gdry
        LOGICAL :: gwet


!       1.  Initialise the types of specie arrays.
!           ---------- --- ----- -- ------ -------

        jpnpx3 = (jpnr/(3*3))+3*3
        nstst = 0
        nf = 0
        DO js = 1, jpspec
          itype = ctype(js)
          nlall(js) = js
          IF (itype==jpfm .OR. itype==jpif .OR. itype==jpna) THEN
            nstst = nstst + 1
            nlstst(nstst) = js
          END IF
          IF ( itype==jpfm .OR. itype==jpif .OR. itype==jpsp ) THEN
            nf = nf + 1
            nlf(nf) = js
          END IF
        END DO

!       2.  Fill in the 3, 2 & 1 term product index arrays.
!           ---- -- --- -- - - - ---- ------- ----- -------

        DO js = 1, jpspec

!         2.1  For each species, scan the reactions for groups
!              of 3, 2 & 1 product terms and loss terms.

          iprd  = 0
          igrp3 = 0
          DO jr = 1, jpnr
            IF ( nfrpx(jr) == 0 ) THEN

!             If reaction does NOT have fractional products.

              DO jp = 3, jpmsp
                IF ( nspi(jr,jp) == js ) THEN
                  iprd = iprd + 1
                  igrp(iprd) = jr
                  IF ( iprd == 3 ) THEN
                    iprd = 0
                    igrp3 = igrp3 + 1
                    DO j = 1, 3
                      nprdx3(j,igrp3,js) = igrp(j)
                    END DO
                  END IF
                END IF
              END DO
            END IF
          END DO
          ngrp(js,3) = igrp3
          IF ( iprd == 2 ) THEN
            ngrp(js,2) = 1
            nprdx2(1,js) = igrp(1)
            nprdx2(2,js) = igrp(2)
          ELSE IF ( iprd == 1 ) THEN
            ngrp(js,1) = 1
            nprdx1(js) = igrp(1)
          END IF

          iloss = 0
          igrp3 = 0
          DO jr = 1, jpnr
            DO jp = 1, 2
              IF ( nspi(jr,jp) == js ) THEN
                iloss = iloss + 1
                igrp(iloss) = jr
                IF ( iloss == 3 ) THEN
                  iloss = 0
                  igrp3 = igrp3 + 1
                  DO j = 1, 3
                    nprdx3(j,igrp3,jpspec+js) = igrp(j)
                  END DO
                END IF
              END IF
            END DO
          END DO

          ngrp(jpspec+js,3) = igrp3
          IF ( iloss == 2 ) THEN
            ngrp(js+jpspec,2) = 1
            nprdx2(1,jpspec+js) = igrp(1)
            nprdx2(2,jpspec+js) = igrp(2)
          ELSE IF ( iloss == 1 ) THEN
            ngrp(js+jpspec,1) = 1
            nprdx1(jpspec+js) = igrp(1)
          END IF

        END DO

!       2.2  Check array size is ok.

        ngmax = ngrp(1,3)
        DO js = 2, jpspec
          ngmax = max( ngmax, ngrp(js,3) )
        END DO
        IF ( ngmax > jpnpx3 ) THEN
           WRITE (6,*) '** ASAD ERROR: The parameter jpnpx3 is too',   &
         ' small. Please change it not less than ',ngmax
          cmessage = ' Parameter jpnpx3 is too small'

          CALL EREPORT('ASAD_INIX',ngmax,cmessage)
        END IF

!       2.3  Build table for fractional products

        nnfrp = 0
        DO jr = 1, jpnr
          IF ( nfrpx(jr) /= 0 ) THEN
            DO jp = 3, jpspb
              is = nspi(jr,jp)
              IF ( is /= 0 ) THEN
                nnfrp = nnfrp + 1
                IF ( nnfrp > jpfrpx ) THEN
                  WRITE (6,*) '*ASAD ERROR: The parameter jpfrpx is',  &
                 ' too small. Please increase it.'
                  cmessage=' Parameter jpfrpx is too small for no. '// &
                            'of fractional products'

                  CALL EREPORT('ASAD_INIX',nnfrp,cmessage)
                END IF
                ntabfp(nnfrp,1) = is
                ntabfp(nnfrp,2) = jr
                ntabfp(nnfrp,3) = nfrpx(jr) + (jp-3)
              END IF
            END DO
          END IF
        END DO

!       3.   Set up lists for deposition and emissions.
!            --- -- ----- --- ---------- --- ----------

        ndepd = 0
        ndepw = 0
        nemit = 0
        DO  js = 1, jpspec
          IF ( ldepd(js) ) THEN
            ndepd = ndepd + 1
            nldepd(ndepd) = js
          END IF
          IF ( ldepw(js) ) THEN
            ndepw = ndepw + 1
            nldepw(ndepw) = js
          END IF
          IF ( lemit(js) ) THEN
            nemit = nemit + 1
            nlemit(nemit) = js
          END IF
        END DO

!       3.1.  Now form the list used in prls.f

        DO js = 1, jpspec
          idepd(js) = nldepd(js)
          idepw(js) = nldepw(js)
        END DO

!       First look for species with both dry and wet on.

        nldepx(1) = 7
        nldepx(2) = nldepx(1) - 1

        DO ispec = 1,jpspec
          gdry = .false.
          gwet = .false.

          DO j = 1, jpspec
            IF ( idepd(j) == ispec ) THEN
              gdry = .true.
              jdry = j
            END IF
            IF ( idepw(j) == ispec ) THEN
              gwet = .true.
              jwet = j
            END IF
          END DO

          IF ( gdry .AND. gwet ) THEN
            nldepx(2)         = nldepx(2) + 1
            nldepx(nldepx(2)) = ispec
            idepd(jdry)       = 0
            idepw(jwet)       = 0
          END IF
        END DO   ! ispec loop

!       Now get the remaining species with dry dep. only

        nldepx(3) = nldepx(2) + 1
        nldepx(4) = nldepx(3) - 1
        DO j = 1, jpspec
          IF ( idepd(j) /= 0 ) THEN
            nldepx(4) = nldepx(4) + 1
            nldepx(nldepx(4)) = idepd(j)
          END IF
        END DO

!       Remaining species with wet deposition on.

        nldepx(5) = nldepx(4) + 1
        nldepx(6) = nldepx(5) - 1
        DO j = 1, jpspec
          IF ( idepw(j) /= 0 ) THEN
            nldepx(6) = nldepx(6) + 1
            nldepx(nldepx(6)) = idepw(j)
          END IF
        END DO

        RETURN
        END SUBROUTINE ASAD_INIX

      End Module messy_clamschem_asad_inix
