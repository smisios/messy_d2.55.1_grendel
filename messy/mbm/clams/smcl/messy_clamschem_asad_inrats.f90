Module messy_clamschem_asad_inrats

contains

!
!**** *inrats* - to read species and reaction data 
!
!     Glenn Carver & Paul Brown, Centre for Atmospheric Science,
!                                University of Cambridge.
!
!
! Purpose: Reads species and reaction data. Combines reactions into one
!          array and reorders them to put single reactant reactions first
!          to improve code in the prls routine.
!
!
!          Called from ASAD_CINIT
!
!
!     Interface
!     ---------
!     Reads the information:
!                  - Chosen chemistry , contains information
!                    on species involved in the chemistry and
!                    the families/tracers to which they belong.
!                  - Data for bimolecular reactions.
!                  - Data for trimolecular reactions.
!                  - Data for photolysis reactions.
!                  - Data for heterogeneous reactions.
!     from chem module
!
!     Method
!     ------
!     The file specifies the species types using 2 letter
!     codes for easier reading.
!
!             ctype         Meaning
!             'FM'          Family member
!             'FT'          Tracer but will be put into a family
!                           if lifetime becomes short.
!             'TR'          Tracer, advected by calling model.
!             'SS'          Steady state species.
!             'CT'          Constant species.
!
!     Local variables
!     ---------------
!     iadv         Counter of number of model tracers/families.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INRATS

        USE messy_clamschem_asad_mod
        USE messy_clamschem_asad_mod_clams, ONLY: mype, prstatus_oper, printstatus, &
                                                  n_chem_diags, &
                                                  jpctr, jpspec, jpbk, jptk, jphk, &
                                                  jppj, jpnr, &
                                                  lhet
        USE messy_clamschem_defs_mod
        USE messy_clams_global,         ONLY: prec
        USE messy_clamschem_global,     ONLY: ipart
        USE messy_clamschem_asad_dummy, ONLY: ereport

        IMPLICIT NONE

!       Local variables

        INTEGER :: errcode                ! Variable passed to ereport

        INTEGER :: ispb(jpbk+1,jpspb)
        INTEGER :: ispt(jptk+1,jpspt)
        INTEGER :: ispj(jppj+1,jpspj)
        INTEGER :: isph(jphk+1,jpsph)
        INTEGER :: ifrpbx(jpbk+1)
        INTEGER :: ifrpjx(jppj+1)
        INTEGER :: ifrptx(jptk+1)
        INTEGER :: ifrphx(jphk+1)
        INTEGER :: ilmin(jpspec)
        INTEGER :: ilmaj(jpspec)
        INTEGER :: ixdumt(jptk+1)
        INTEGER :: ixdumj(jppj+1)
        INTEGER :: ixdumh(jphk+1)
        INTEGER :: ierror
        INTEGER :: iadv                  ! Counter for advected species
        INTEGER :: inadv                 ! Counter for non-advected species
        INTEGER :: imajor                ! Counter
        INTEGER :: iminor                ! Counter
        INTEGER :: ix                    ! Counter
        INTEGER :: icount                ! Counter
        INTEGER :: ifam                  ! Index
        INTEGER :: imaj                  ! Index
        INTEGER :: iflag                 ! Used to test family order
        INTEGER :: idummy                ! Dummy variable
        INTEGER :: istat                 ! Tag for communication
        INTEGER :: j                     ! Loop variable
        INTEGER :: jf                    ! Loop variable
        INTEGER :: jadv                  ! Loop variable
        INTEGER :: jb                    ! Loop variable
        INTEGER :: jctr                  ! Loop variable
        INTEGER :: jh                    ! Loop variable
        INTEGER :: jj                    ! Loop variable
        INTEGER :: jp                    ! Loop variable
        INTEGER :: jr                    ! Loop variable
        INTEGER :: js                    ! Loop variable
        INTEGER :: jspb                  ! Loop variable
        INTEGER :: jsph                  ! Loop variable
        INTEGER :: jspj                  ! Loop variable
        INTEGER :: jspt                  ! Loop variable
        INTEGER :: jt                    ! Loop variable
        INTEGER :: k                     ! Loop variable
        INTEGER :: ind                   ! Loop index

        REAL(PREC) :: zdumt(1)                 ! Dummy array for trimol file
        REAL(PREC) :: zdumj(1)                 ! Dummy array for photol file
        REAL(PREC) :: zdumh(1)                 ! Dummy array for heter file

        CHARACTER (LEN=10) :: cmntb(jpbk+1)
        CHARACTER (LEN=10) :: cmntt(jptk+1)
        CHARACTER (LEN=10) :: cmntj(jppj+1)
        CHARACTER (LEN=10) :: cmnth(jphk+1)
        CHARACTER (LEN=10), PARAMETER :: nullx='          '
        CHARACTER (LEN=72) :: cmessage        ! Error message

        LOGICAL :: gtype3
        LOGICAL :: gdebug
        LOGICAL :: L_exist
        LOGICAL :: L_fa


!       1.  Read chosen chemistry file and determine chemistry
!           ---- ------ --------- ---- --- --------- ---------

! Initialise local counters, see asad_cinit for (e.g.) spb and frpb 
        ispb(:,:) = 0
        ifrpbx(:) = 0
        ispt(:,:) = 0
        ifrptx(:) = 0
        ispj(:,:) = 0
        ifrpjx(:) = 0
        isph(:,:) = 0
        ifrphx(:) = 0

!       1.1.1  Set types and find which species are in families.

! Replace read statement by module:

        IF (size(chch_defs) /= jpspec) THEN
          errcode=1
          cmessage=' jpspec and chch_defs are inconsistent'
          WRITE(6,*) cmessage

          CALL EREPORT('ASAD_INRATS',errcode,cmessage)
        ENDIF
        DO k=1,jpspec
          speci(k)  = chch_defs(k)%speci
          nodd(k)   = chch_defs(k)%nodd
          ctype(k)  = chch_defs(k)%ctype
          family(k) = chch_defs(k)%family
        ENDDO

        iadv = 0
        inadv = 0
        ntrf = 0
        nnaf = 0
        ntr3 = 0
        DO js = 1, jpspec
          gtype3 = .false.
          IF ( ctype(js) /= jpfm .AND. ctype(js) /= jpif )              &
            family(js)='          '
          IF ( ctype(js) == jpsp ) THEN               ! tracers
            iadv = iadv + 1
            ntrf = ntrf + 1
            IF ( iadv > jpctr ) THEN
              WRITE (6,*) '** ASAD ERROR in subroutine inrats'
              WRITE (6,*) '** Parameter jpctr is too low; found',iadv
              WRITE (6,*) '** tracers so far with ',jpspec-js
              WRITE (6,*) '** species to check.'
              cmessage = 'ASAD ERROR: jpctr is too low'

              CALL EREPORT('ASAD_INRATS',iadv,cmessage)
            ENDIF
            advt(iadv)  = speci(js)
            nltrf(ntrf) = iadv
          ELSE IF( ctype(js) == jpna )THEN        ! non-advected species 
            inadv = inadv + 1 
            nnaf = nnaf + 1 
            IF ( inadv > n_chem_diags ) THEN 
              cmessage = 'ASAD ERROR: Too many non-advected tracers' 
 
              CALL EREPORT('ASAD_INRATS',iadv,cmessage) 
            ENDIF 
            nadvt(inadv)  = speci(js) 
          ELSE IF( ctype(js) == jpfm .OR. ctype(js) == jpif )THEN
! ju_nt_20140205             
!!$            cmessage = ' Family chemistry not available in this version'
!!$            CALL EREPORT('ASAD_INRATS',js,cmessage)
            IF ( ctype(js) == jpif ) THEN
              iadv = iadv + 1
              ntr3 = ntr3 + 1
              IF ( iadv  >   jpctr ) THEN
                WRITE (6,*) '** ASAD ERROR in subroutine inrats'
                WRITE (6,*) '** Parameter jpctr is too low; found',iadv
                WRITE (6,*) '** tracers so far with ',jpspec-js
                WRITE (6,*) '** species to check.'
                cmessage = 'ERROR in jpctr'

                CALL EREPORT('ASAD_INRATS',js,cmessage)
              ENDIF
              advt(iadv)  = speci(js)
              nltr3(ntr3) = iadv
            ENDIF
            l_fa=.true.
            DO jadv = 1, iadv
              IF ( family(js) == advt(jadv) ) THEN
                l_fa=.false.
                EXIT
              ENDIF
            ENDDO
            IF (l_fa) THEN
              iadv = iadv + 1
              ntrf = ntrf + 1
              IF ( iadv  >   jpctr ) then
                WRITE (6,*) '***** ASAD ERROR in subroutine inrats'
                WRITE (6,*) '***** Param jpctr is too low; found',iadv
                WRITE (6,*) '***** tracers so far with ',jpspec-js
                WRITE (6,*) '***** species to check.'
                cmessage = 'INRATS ERROR : jpctr is too low'

                CALL EREPORT('ASAD_INRATS',iadv,cmessage)
              ENDIF
              advt(iadv) = family(js)
              nltrf(ntrf) = iadv
            ENDIF      ! l_fa
          ENDIF
        ENDDO

!       1.2 Find major species of families

        DO jadv = 1, iadv
          DO js = 1, jpspec
            IF( family(js) == advt(jadv) .and. ctype(js) /= jpif )     &
              majors(jadv) = js
            IF( speci(js) == advt(jadv) ) majors(jadv) = js
          ENDDO
        ENDDO

!       1.3 Allocate families to species

        DO js = 1, jpspec
          moffam(js) = 0
          madvtr(js) = 0
          DO jadv = 1, iadv
            IF (family(js) == advt(jadv) ) moffam(js) = jadv
            IF (speci(js)  == advt(jadv) ) madvtr(js) = jadv
            IF (family(js) == advt(jadv) .AND. js > majors(jadv)) THEN
              WRITE (6,*) '** ASAD ERROR: '
              WRITE (6,*) 'RE-ORDER SPECIES FILE SO THAT THE MAJOR '
              WRITE (6,*) 'SPECIES OF A FAMILY OCCURS AFTER THE OTHERS'
              cmessage = 'INRATS ERROR : Order of species is incorrect'

              CALL EREPORT('ASAD_INRATS',jadv,cmessage)
            ENDIF
          ENDDO
        ENDDO

!       1.4  Build the list of major and minor species

        nlmajmin(1) = 5
        imajor      = 0
        iminor      = 0
        DO js = 1, jpspec
          ifam = moffam(js)
          imaj = 0
          IF ( ifam /= 0 ) THEN
           imaj = majors(ifam)
           IF ( imaj /= js ) THEN
              iminor = iminor + 1
              ilmin(iminor) = js
            ELSE
              imajor = imajor + 1
              ilmaj(imajor) = js
            ENDIF
          ENDIF
        ENDDO
        nlmajmin(2) = nlmajmin(1) + imajor - 1
        nlmajmin(3) = nlmajmin(1) + imajor
        nlmajmin(4) = nlmajmin(3) + iminor - 1
        DO j = nlmajmin(1), nlmajmin(2)
          nlmajmin(j) = ilmaj(j-nlmajmin(1)+1)
        ENDDO
        DO j = nlmajmin(3), nlmajmin(4)
          nlmajmin(j) = ilmin(j-nlmajmin(3)+1)
        ENDDO

        IF ( iadv /= jpctr ) then
          WRITE (6,*) '** ASAD ERROR: Number of advected tracers',     &
           ' specified in chch_defs does not match jpctr'
          WRITE (6,*) 'Found ',iadv,' but expected ',jpctr
          cmessage = 'INRATS ERROR : iadv and jpctr do not match'

          CALL EREPORT('ASAD_INRATS',iadv,cmessage)
        ENDIF

!       2.  Write details of chemistry selection to log file
!           ----- ------- -- --------- --------- -- --- ----

        IF (mype == 0 .AND. printstatus >= prstatus_oper .and. ipart==1) THEN
          WRITE(6,*)
          WRITE(6,*)'  ***  CHEMISTRY INFORMATION  ***'
          WRITE(6,*)
          WRITE(6,*)'ASAD IS TREATING ADVECTED TRACERS IN THE ORDER:'
          WRITE(6,*)
          WRITE(6,'(5(2x,i2,1x,a10))')(jctr,advt(jctr), jctr= 1,jpctr)
          WRITE(6,*)
          WRITE(6,*)'ASAD IS TREATING NON-ADVECTED TRACERS IN THE ORDER:' 
          WRITE(6,*) 
          WRITE(6,'(5(2x,i2,1x,a10))')(jctr,nadvt(jctr), jctr= 1,nnaf) 
          WRITE(6,*) 
          WRITE(6,*)'IF THE TRACERS WERE NOT INITIALISED IN THIS '
          WRITE(6,*)'ORDER THEN THE MODEL RESULTS ARE WORTHLESS '
          WRITE(6,*)

          iflag = 0
          DO jctr = 1, jpctr
            IF (advt(jctr) /= speci(majors(jctr)) ) THEN
              IF (iflag == 0 )then
                WRITE(6,*)'THE MAJOR MEMBER OF EACH OF THE FAMILIES'
                WRITE(6,*)'IS GIVEN BELOW. IF THIS IS NOT ACCEPTABLE,'
                WRITE(6,*)'THEN YOU MUST REORDER THE SPECIES IN'//      &
                          ' CHCH_DEFS'
                WRITE(6,*)'SO THE MAJOR SPECIES FOLLOWS THE OTHERS.'
                WRITE(6,*)
                iflag = 1
              ENDIF
              WRITE(6,'(a10,1x,a10)') advt(jctr), speci(majors(jctr))
            ENDIF
          ENDDO
        ENDIF     ! End of IF mype statement

!       3.  Bimolecular ratefile
!           ----------- --------

!       Get bimolecular rates from module

        IF (size(ratb_defs) /= jpbk) THEN 
          errcode=1
          cmessage='size of ratb_defs is inconsistent with jpbk'

          CALL EREPORT('ASAD_INRATS',errcode,cmessage)
        END IF 
        icount=1
        DO k=1,jpbk
          spb(k,1) = ratb_defs(k)%react1
          spb(k,2) = ratb_defs(k)%react2
          spb(k,3) = ratb_defs(k)%prod1
          spb(k,4) = ratb_defs(k)%prod2
          spb(k,5) = ratb_defs(k)%prod3
          spb(k,6) = ratb_defs(k)%prod4
          ab(k,1)  = ratb_defs(k)%K0
          ab(k,2)  = ratb_defs(k)%alpha
          ab(k,3)  = ratb_defs(k)%beta
          IF (ratb_defs(k)%pyield1 > 1e-18) THEN
            ifrpbx(k)     = icount
            frpb(icount)  = ratb_defs(k)%pyield1
            IF (spb(k,4) /= nullx) frpb(icount+1) = ratb_defs(k)%pyield2
            IF (spb(k,5) /= nullx) frpb(icount+2) = ratb_defs(k)%pyield3
            IF (spb(k,6) /= nullx) frpb(icount+3) = ratb_defs(k)%pyield4
            icount = icount + 4
          ENDIF
        ENDDO

        DO jb = 1, jpbk
          DO js = 1, jpspec
            DO jspb = 1, jpspb
              IF ( speci(js) == spb(jb,jspb) ) ispb(jb,jspb) = js
            ENDDO
          ENDDO
        ENDDO

! Load in bimol fractional prod coefs to frpx array
        DO jf = 1, jpfrpb
          frpx(jf) = frpb(jf)
        END DO

!       4.  Trimolecular ratefile
!           ------------ --------

!       Get trimolecular rates from module for UM version

        IF (size(ratt_defs) /= jptk) THEN
          errcode=1
          cmessage='size of ratt_defs is inconsistent with jptk'

          CALL EREPORT('ASAD_INRATS',errcode,cmessage)
        ENDIF
        icount=1
        DO k=1,jptk
          spt(k,1) = ratt_defs(k)%react1
          spt(k,2) = ratt_defs(k)%react2
          spt(k,3) = ratt_defs(k)%prod1
          spt(k,4) = ratt_defs(k)%prod2
          at(k,1)  = ratt_defs(k)%F
          at(k,2)  = ratt_defs(k)%K1
          at(k,3)  = ratt_defs(k)%alpha1
          at(k,4)  = ratt_defs(k)%beta1
          at(k,5)  = ratt_defs(k)%K2
          at(k,6)  = ratt_defs(k)%alpha2
          at(k,7)  = ratt_defs(k)%beta2
          IF (ratt_defs(k)%pyield1 > 1e-18) THEN
            ifrptx(k)      = icount
            frpt(icount)   = ratt_defs(k)%pyield1
            IF (spt(k,4) /= nullx) frpt(icount+1) = ratt_defs(k)%pyield2
            icount = icount + 2
          ENDIF
        ENDDO

        DO jt = 1, jptk
          DO js = 1, jpspec
            DO jspt = 1, jpspt
              IF (speci(js) == spt(jt,jspt) ) ispt(jt,jspt) = js
            ENDDO
          ENDDO
        ENDDO

! Load in trimol fractional prod coefs to frpx array
        DO jf = 1, jpfrpt
          ind = jf + jpfrpb
          frpx(ind) = frpt(jf)
        END DO


!       5.  Photolysis ratefile
!           ---------- --------

!       use module to get spj

        IF (size(ratj_defs) /= jppj) THEN
          errcode=1
          cmessage='size of ratj_defs is not equal to jppj'

          CALL EREPORT('ASAD_INRATS',errcode,cmessage)
        ENDIF
        icount=1
        DO k=1,jppj
          spj(k,1) = ratj_defs(k)%react1
          spj(k,2) = ratj_defs(k)%react2
          spj(k,3) = ratj_defs(k)%prod1
          spj(k,4) = ratj_defs(k)%prod2
          spj(k,5) = ratj_defs(k)%prod3
          spj(k,6) = ratj_defs(k)%prod4
          IF (ratj_defs(k)%pyield1 > 1e-18) THEN
            ifrpjx(k)     = icount
            frpj(icount)  = ratj_defs(k)%pyield1
            IF (spj(k,4) /= nullx) frpj(icount+1) = ratj_defs(k)%pyield2
            IF (spj(k,5) /= nullx) frpj(icount+2) = ratj_defs(k)%pyield3
            IF (spj(k,6) /= nullx) frpj(icount+3) = ratj_defs(k)%pyield4
            icount = icount + 4
          ENDIF
        ENDDO

        DO jj = 1, jppj
          DO js = 1, jpspec
            DO jspj = 1, jpspj
              IF (speci(js) == spj(jj,jspj) ) ispj(jj,jspj) = js
            ENDDO
          ENDDO
        ENDDO

! Load in photol fractional prod coefs to frpx array
        DO jf = 1, jpfrpj
          ind = jf + jpfrpb + jpfrpt
          frpx(ind) = frpj(jf)
        END DO

!       6.  Heterogeneous ratefile
!           ------------- --------

!!!!!
! ju_nt_20140210
!$$        IF (jphk > 0) THEN
        IF (jphk > 0 .and. lhet) THEN

!         use module to get sph

          IF (size(rath_defs) /= jphk) THEN
            errcode=1
            cmessage='size of rath_defs is not equal to jphk'

            CALL EREPORT('ASAD_INRATS',errcode,cmessage)
          ENDIF
          icount=1
          DO k=1,jphk
            sph(k,1) = rath_defs(k)%react1
            sph(k,2) = rath_defs(k)%react2
            sph(k,3) = rath_defs(k)%prod1
            sph(k,4) = rath_defs(k)%prod2
            sph(k,5) = rath_defs(k)%prod3
            sph(k,6) = rath_defs(k)%prod4
            IF (rath_defs(k)%pyield1 > 1e-18) THEN
              ifrphx(k)     = icount
              frph(icount)  = rath_defs(k)%pyield1
              IF (sph(k,4) /= nullx) frph(icount+1)=rath_defs(k)%pyield2
              IF (sph(k,5) /= nullx) frph(icount+2)=rath_defs(k)%pyield3
              IF (sph(k,6) /= nullx) frph(icount+3)=rath_defs(k)%pyield4
              icount = icount + 4
            ENDIF
          ENDDO

          DO jh = 1, jphk
            DO js = 1, jpspec
              DO jsph = 1, jpsph
                IF (speci(js) == sph(jh,jsph) ) isph(jh,jsph) = js
              ENDDO
            ENDDO
          ENDDO

! Load in het fractional prod coefs to frpx array
          DO jf = 1, jpfrph
            ind = jf + jpfrpb + jpfrpt + jpfrpj
            frpx(ind) = frph(jf)
          END DO

        END IF       ! jphk > 0

!       7.  Reorder reactions, putting single reactants first.
!           ------- ---------- ------- ------ --------- ------

        nuni = 0

!       7.1  Single reactants; scan ratefiles in turn.

!
        DO jr = 1, jpbk
          IF (ispb(jr,2) == 0 ) THEN
            nuni = nuni + 1
            nbrkx(jr) = nuni
            DO jp = 1, jpspb
              nspi(nuni,jp) = ispb(jr,jp)
            ENDDO
            IF ( ifrpbx(jr) /= 0 ) nfrpx(nuni) = ifrpbx(jr)
          ENDIF
        ENDDO

        DO jr = 1, jptk
          IF ( ispt(jr,2) == 0 ) THEN
            nuni = nuni + 1
            ntrkx(jr) = nuni
            DO jp = 1, jpspt
              nspi(nuni,jp) = ispt(jr,jp)
            ENDDO
            IF (ifrptx(jr) /= 0) nfrpx(nuni) = ifrptx(jr)+jpfrpb
          ENDIF
        ENDDO

        DO jr = 1, jppj
          IF ( ispj(jr,2) == 0 ) THEN
            nuni = nuni + 1
            nprkx(jr) = nuni
            DO jp = 1, jpspj
              nspi(nuni,jp) = ispj(jr,jp)
            ENDDO
            IF (ifrpjx(jr) /= 0) nfrpx(nuni) = ifrpjx(jr)+jpfrpb+jpfrpt
          ENDIF
        ENDDO

!!!!!
! ju_nt_20140210
!$$        IF ( jphk > 0 ) THEN
        IF ( jphk > 0 .and. lhet) THEN
          DO jr = 1, jphk
            IF ( isph(jr,2) == 0 ) THEN
              nuni = nuni + 1
              nhrkx(jr) = nuni
              DO jp = 1, jpsph
                nspi(nuni,jp) = isph(jr,jp)
              ENDDO
              IF (ifrphx(jr) /= 0) nfrpx(nuni) = ifrphx(jr)+jpfrpb+         &
                                                 jpfrpt+jpfrpj
            ENDIF
          ENDDO
        ENDIF

!       7.2  Two reactants; copy remaining reactions

        ix = nuni
        DO jr = 1, jpbk
          IF ( ispb(jr,2) /= 0 ) THEN
            ix = ix + 1
            nbrkx(jr) = ix
            DO jp = 1, jpspb
              nspi(ix,jp) = ispb(jr,jp)
            ENDDO

            IF ( ifrpbx(jr) /= 0 ) nfrpx(ix) = ifrpbx(jr)
          ENDIF
        ENDDO

        DO jr = 1, jptk
          IF ( ispt(jr,2) /= 0 ) THEN
            ix = ix + 1
            ntrkx(jr) = ix
            DO jp = 1, jpspt
              nspi(ix,jp) = ispt(jr,jp)
            ENDDO
            IF (ifrptx(jr) /= 0) nfrpx(ix) = ifrptx(jr)+jpfrpb
          ENDIF
        ENDDO

        DO jr = 1, jppj
          IF (ispj(jr,2) /= 0 ) THEN
            ix = ix + 1
            nprkx(jr) = ix
            DO jp = 1, jpspj
              nspi(ix,jp) = ispj(jr,jp)
            ENDDO
            IF (ifrpjx(jr) /= 0) nfrpx(ix) = ifrpjx(jr)+jpfrpb+jpfrpt
          ENDIF
        ENDDO

!!!!!
! ju_nt_20140210
!$$        IF ( jphk > 0 ) THEN
        IF ( jphk > 0 .and. lhet) THEN
          DO jr = 1, jphk
            IF ( isph(jr,2) /= 0 ) THEN
              ix = ix + 1
              nhrkx(jr) = ix
              DO jp = 1, jpsph
                nspi(ix,jp) = isph(jr,jp)
              ENDDO
              IF (ifrphx(jr) /= 0) nfrpx(ix) = ifrphx(jr)+jpfrpb+       &
                                               jpfrpt+jpfrpj
            ENDIF
          ENDDO
        ENDIF

!!!!!
! ju_nt_20140210
!$$        IF ( ix /= jpnr ) THEN
        IF ( (lhet .and. ix /= jpnr) .or. &
             (.not. lhet .and. ix /= jpnr-jphk) ) THEN
          WRITE (6,*) '*** INTERNAL ASAD ERROR: Number of reactions',   &
                      ' placed in nspi array does not equal jpnr. '
          WRITE (6,*) '                         Check that reaction',   &
                      ' files and value of jpnr in chem_data module are' 
          WRITE (6,*) '                         consistent. Found: ',   &
                      ix,' jpnr: ',jpnr
          cmessage = 'No of rxns in nspi array is not equal to jpnr'

          CALL EREPORT('ASAD_INRATS',ix,cmessage)
        ENDIF

        RETURN
        END SUBROUTINE ASAD_INRATS


      End Module messy_clamschem_asad_inrats
