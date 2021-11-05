Module messy_clamschem_asad_cinit

contains

!**** *asad_cinit* - chemistry initialization routine.
!
!     Glenn D. Carver            Centre for Atmospheric Science,
!                                University of Cambridge.
!
!
!     slightly modified, Rolf Mueller
!     call to inhet commented out, is called from main routine
!
! Purpose: To initialize variables used in the chemistry
!
!          Called from chem before call to  ASAD_CDRIVE
!
!
!     Method
!     ------
!     Input arguments are checked and copied to common. Species
!     and reaction data are read in. Other variables used in the
!     chemistry are initialised.
!
!     -- Calculation of peps
!
!     At several places in the ASAD code, we needed a min. value
!     to use to guard against zero divides. We have to compute
!     this to allow for all possible computer hardwares and precisions
!     that ASAD might be run at.
!
!     Externals
!     ---------
!     inrats      - Reads species and reaction data.
!     inphot      - Initialises photolysis scheme.
!     inwdep      - Initialises wet deposition data
!                   (user supplied routine)
!     inddep      - Initialises dry deposition data
!                   (user supplied routine)
!     inemit      - Initialises emission data
!                   (user supplied routine).
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_CINIT(p_field)

        USE messy_clamschem_asad_mod
        USE messy_clamschem_asad_dummy, ONLY: ereport, &
                                              inwdep, inddep
        USE messy_clams_global,     ONLY: prec, rank
        USE messy_clamschem_global, ONLY: ipart
        USE messy_clamschem_inemit, ONLY: inemit
        USE messy_clamschem_inphot, ONLY: inphot
        USE messy_clamschem_asad_inijac,    ONLY: ASAD_INIJAC
        USE messy_clamschem_asad_inix,      ONLY: ASAD_INIX
        USE messy_clamschem_asad_inrats,    ONLY: ASAD_INRATS
        USE messy_clamschem_asad_setsteady, ONLY: ASAD_SETSTEADY
        USE messy_clamschem_asad_mod_clams, ONLY: jpnr, jpctr, jpspec, &
                                       printstatus, prstatus_oper, &
                                       lphotol

        IMPLICIT NONE


        INTEGER, INTENT(IN) :: p_field

!       Local variables

        INTEGER :: j                      ! Loop variable
        INTEGER :: jc                     ! Loop variable
        INTEGER :: jf                     ! Loop variable
        INTEGER :: jg                     ! Loop variable
        INTEGER :: jl                     ! Loop variable
        INTEGER :: jp                     ! Loop variable
        INTEGER :: jr                     ! Loop variable
        INTEGER :: js                     ! Loop variable
        INTEGER :: jtr                    ! Loop variable
        INTEGER :: jx                     ! Loop variable
        INTEGER :: jpnpx3                 ! Loop variable
        INTEGER :: errcode                ! Variable passed to ereport
        INTEGER, PARAMETER :: nrsteps_max = 200  ! max steps
        INTEGER, PARAMETER :: nit0_max    = 50   ! max
        INTEGER, PARAMETER :: nitfg_max   = 50   ! max
        INTEGER, PARAMETER :: nitnr_max   = 50   ! max

        REAL(PREC) :: sfmin

        CHARACTER (LEN=10), PARAMETER :: nullx='          '
        CHARACTER (LEN=72) :: cmessage    ! Error message


!       Initialisations for COMMON variables.

        ljacx = .true.
        
        jpnpx3 = (jpnr/(3*3))+3*3


!       1.  Copy input arguments to COMMON and check them.
!           ---- ----- --------- -- ------ --- ----- -----

!       Logical arguments
        DO j = 1, jpspec
          lemit(j) = .false.
        ENDDO

        IF (nrsteps < 0 .OR. nrsteps > 200) THEN
          cmessage = ' NRSTEPS IS OUT OF RANGE, RESETTING'
          write(6,*) 'NRSTEPS = ',nrsteps,' Reset to: ',nrsteps_max
          nrsteps = nrsteps_max
          errcode=-1

          CALL EREPORT('ASAD_CINIT',errcode,cmessage)
        ENDIF

        IF (nit0 < 0 .OR. nit0 > 50) THEN
          cmessage = ' NIT0 IS OUT OF RANGE, RESETTING'
          write(6,*) 'NIT0 = ',nit0,' Reset to: ',nit0_max
          nit0 = nit0_max
          errcode=-1

          CALL EREPORT('ASAD_CINIT',errcode,cmessage)
        ENDIF

        IF (nitfg < 0 .OR. nitfg > 50) THEN
          cmessage = ' NITFG IS OUT OF RANGE, RESETTING'
          write(6,*) 'NITFG = ',nitfg,' Reset to: ',nitfg_max
          nitfg = nitfg_max
          errcode=-1

          CALL EREPORT('ASAD_CINIT',errcode,cmessage)
       ENDIF

       IF (nitnr < 0 .OR. nitnr > 50 )then
         cmessage = ' NITNR IS OUT OF RANGE, RESETTING'
         write(6,*) 'NITNR = ',nitnr,' Reset to: ',nitnr_max
         nitnr = nitnr_max
         errcode=-1

         CALL EREPORT('ASAD_CINIT',errcode,cmessage)
       ENDIF

!      2.  Assign chemical timestep values for ASAD
!          ----------------------------------------

! Set up timestep counting. Interval depends on solver: IMPACT and Rosenbrock
! are run every dynamical timestep. Newton-Raphson is run every 2 / 3
! timesteps for a 30/20 minutes dynamical timestep.

        IF (method == 3) THEN          ! N-R solver
! ju_nt_20140319: use "cdt" from chem
!!$          ncsteps = 1
!!$          cdt = REAL(kcdt)
!!$          interval = NINT(cdt/dtime)
           ncsteps = int( dtime / cdt )
           if ( ncsteps*cdt /= dtime ) then
              ncsteps = ncsteps + 1
              cdt     = dtime / float(ncsteps)
           endif
           interval = 1 !???
        ELSE IF (method == 1) THEN     ! IMPACT
! use about 15 or 10 minutes, depending on dynamical timestep
! ju_nt_20140217: use "cdt" from chem
!!$          IF (dtime < tslimit) THEN
!!$            cdt = dtime
!!$            ncsteps = 1
!!$            interval = 1
!!$          ELSE
!!$            cdt = 0.5*dtime
!!$            ncsteps = 2
!!$            interval = 1
!!$          END IF
           ncsteps = int( dtime / cdt )
           if ( ncsteps*cdt /= dtime ) then
              ncsteps = ncsteps + 1
              cdt     = dtime / float(ncsteps)
           endif
           interval = 1 !???
! ju_nt_20140304: set timestep for SVODE 
        ELSE IF (method == 11) THEN  ! SVODE
           ncsteps = int( dtime / cdt )
           if ( ncsteps*cdt /= dtime ) then
              ncsteps = ncsteps + 1
              cdt     = dtime / float(ncsteps)
           endif
           interval = 1 !???
        ELSE                           ! Other solvers
          cdt = dtime
          ncsteps = 1
          interval = 1
        END IF

! ju_nt_20140403
!        IF (printstatus >= prstatus_oper) THEN
        IF (printstatus >= prstatus_oper .and. rank==0 .and. ipart==1) THEN
          WRITE(6,*) 'Interval for chemical solver set to: ', interval
          WRITE(6,*) 'Timestep for chemical solver set to: ', cdt
          WRITE(6,*) 'No steps for chemical solver set to: ', ncsteps
        END IF

        IF (ABS(cdt*ncsteps - dtime*interval) > 1e-4) THEN
          cmessage=' chemical timestep does not fit dynamical timestep'
          errcode = kcdt
          CALL EREPORT('ASAD_CINIT',errcode,cmessage)
        END IF

!       2.1  Set photolysis frequency.

        IF (kfphot < 0 .AND. abs(kfphot) > dtime) THEN
          write (6,*) '**CINIT WARNING: VALUE OF KFPHOT ',kfphot,      &
          'EXCEEDS THE MODEL TIMESTEP. ROUTINE PHOTOL WILL ONLY',      &
          ' BE CALLED ONCE.'
          nfphot = 0
        ELSEIF ( kfphot > 0 .AND. kfphot > ncsteps ) THEN
          write (6,*) '**CINIT WARNING: FREQUENCY KFPHOT ',kfphot,     &
           ' EXCEEDS THE TOTAL NUMBER OF CHEMICAL SUBSTEPS. ROUTINE ', &
           ' PHOTOL WILL BE CALLED ONCE ONLY.'
          nfphot = 0
        ELSEIF (kfphot < 0) THEN
          nfphot = int( abs(kfphot)/cdt )
        ELSE
          nfphot = kfphot
        ENDIF

!       2.2  Compute minimum safe value (see Method above)

!!!!!
! ju_nt_20140129
        !sfmin = tiny(1.0d0)
        sfmin = tiny(sfmin)
!!!!! folgende Zeile aus "alten" ICG-ASAD-Code ergaenzt: ???
!!!!! -> verursacht "Floating point underflow"-Meldung beim Compilieren !!!
        !if (huge(sfmin)*tiny(sfmin) < 1.0) sfmin=1./huge(sfmin)
        sfmin = 10.0d0**(int(log10(sfmin))+1)
        peps  = 1.0d19 * sfmin
!!!!! ACHTUNG: fuer TESTS:
!!!!!   peps = 1.0d-38
        ! jug test
        ! if (peps < 1E-40) peps = 1E-40

!       3.  Set fixed vmrs (added for CLaMS):
! ju_nt_20140219        
        fco2 = 350.0e-6
        fh2 = 5.0e-7
        fn2 = 0.78084
        fo2 = 0.20945
        fch4 = 1.9e-6

!!!!! => sub. asad_cinit_loop
!!$!       4.  Clear the species arrays
!!$
!!$        f      = 0.0
!!$        fdot   = 0.0
!!$        ej     = 0.0
!!$        linfam = .false.
!!$        linfam = .false.
!!$
!!$        y    = 0.0
!!$        ydot = 0.0
!!$        prod = 0.0
!!$        slos = 0.0
!!$        dpd  = 0.0
!!$        dpw  = 0.0
!!$        emr  = 0.0

!!!!! => sub. asad_cinit_loop
!!$!       5.   Clear the rates and index arrays.
!!$!            ----- --- ----- --- ----- -------
!!$        rk   = 0.0
!!$        prk  = 0.0
        nspi = 0

        DO js = 1, jpspec
          ngrp(js,1)          = 0
          ngrp(js,2)          = 0
          ngrp(js,3)          = 0
          nprdx2(1,js)        = 0
          nprdx2(2,js)        = 0
          nprdx1(js)          = 0
          ngrp(js+jpspec,1)   = 0
          ngrp(js+jpspec,2)   = 0
          ngrp(js+jpspec,3)   = 0
          nprdx2(1,js+jpspec) = 0
          nprdx2(2,js+jpspec) = 0
          nprdx1(js+jpspec)   = 0
          nlall(js)           = 0
          nlstst(js)          = 0
          nlf(js)             = 0
          nlmajmin(js)        = 0
          nldepd(js)          = 0
          nldepw(js)          = 0
          nlemit(js)          = 0
          nldepx(js)          = 0
        ENDDO

        DO js = 1, 2*jpspec
          DO jx = 1, jpnpx3
            nprdx3(1,jx,js) = 0
            nprdx3(2,jx,js) = 0
            nprdx3(3,jx,js) = 0
          ENDDO
        ENDDO

        nbrkx = 0
        ntrkx = 0
        nprkx = 0
        nhrkx = 0

        njacx3(1,:,:) = 0
        njacx3(2,:,:) = 0
        njacx3(3,:,:) = 0

        njcgrp(:,1) = 0
        njcgrp(:,2) = 0
        njcgrp(:,3) = 0
        njacx2(1,:) = 0
        njacx2(2,:) = 0
        njacx1(:)   = 0
        nltrf(:)    = 0
        nltr3(:)    = 0

        DO jc = 1, jpctr
          nmpjac(jc) = 0
          DO jp = 1, jppjac
            npjac1(jp,jc) = 0
          ENDDO
        ENDDO

! Initialise the character arrays
        spb(:,:) = nullx
        spt(:,:) = nullx
        spj(:,:) = nullx
        sph(:,:) = nullx

! Initialise the fractional product arrays
        frpb(:)  = 0.0
        frpt(:)  = 0.0
        frpj(:)  = 0.0
        frph(:)  = 0.0
        frpx(:)  = 0.0
        nfrpx    = 0

        ntabfp(1:jpfrpx,1) = 0
        ntabfp(1:jpfrpx,2) = 0
        ntabfp(1:jpfrpx,3) = 0
        nmzjac = 0
        nmsjac = 0
        nzjac1 = 0
        nsjac1 = 0
        ntabpd = 0
        ztabpd = 0.0
        npdfr  = 0

!       6.  Read chemistry data
!           ---- --------- ----

! DEPENDS ON: asad_inrats
        CALL asad_inrats

! Check that deposition and emission is not on for constant species
        DO js = 1, jpspec
          IF ( ldepd(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Dry deposition turned on for constant species'

            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
          IF ( ldepw(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Wet deposition turned on for constant species'

            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
          IF ( lemit(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Emission turned on for constant species'

            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
        ENDDO

        IF ( method == 3) THEN   ! For Newton-Raphson solver only
! DEPENDS ON: asad_setsteady
           CALL asad_setsteady    ! Initialize steady-state species
!!!!! ???
! ju_nt_20140319
!           CALL EREPORT('ASAD_CINIT',1,'asad_setsteady missing!!!')
        ENDIF

! ju_nt_20140226
!!$        IF ( method >= 10 ) THEN
!!$          DO j = 1, jpspec
!!$            IF ( ctype(j)  ==  jpfm .or. ctype(j)  ==  jpif ) THEN
!!$              WRITE(6,*) '*** ASAD ERROR: You cannot use families ',   &
!!$            ' with one of the stiff integrators. If method  >=  10 ',  &
!!$            ' you cannot have species specified as ',jpfm,' or ',jpif
!!$              cmessage = 'ASAD ABORTED'
!!$
!!$              CALL EREPORT('ASAD_CINIT',j,cmessage)
!!$            ENDIF
!!$          ENDDO
!!$        ENDIF

!       7.  Set up the index arrays.
!           --- -- --- ----- -------

! DEPENDS ON: asad_inix
        CALL asad_inix
! DEPENDS ON: asad_inijac
        CALL asad_inijac

!       8.  Initialise photolysis and heterogeneous chemistry
!           ---------- ---------- --- ------------- ---------

! These are dummy routines at vn7.0 of the UM
!! DEPENDS ON: asad_inphot
!        CALL asad_inphot
! ju_nt_20140131
         if (lphotol) call inphot 
!! DEPENDS ON: asad_inhet
!        CALL asad_inhet

!       9.  Read deposition and emission data
!           ---- ---------- --- -------- ----

!!!!!
! ju_nt_20140128
        IF ( ndepw /= 0 ) CALL INWDEP
        IF ( ndepd /= 0 ) CALL INDDEP
        IF ( nemit /= 0 ) CALL inemit

        RETURN
        END SUBROUTINE ASAD_CINIT




SUBROUTINE ASAD_CINIT_loop

  USE messy_clamschem_asad_mod

  IMPLICIT NONE

!       4.  Clear the species arrays

        f      = 0.0
        fdot   = 0.0
        ej     = 0.0
        linfam = .false.

        y    = 0.0
        ydot = 0.0
        prod = 0.0
        slos = 0.0
        dpd  = 0.0
        dpw  = 0.0
        emr  = 0.0


!       5.   Clear the rates and index arrays.
!            ----- --- ----- --- ----- -------
        rk   = 0.0
        prk  = 0.0
!        nspi = 0

END SUBROUTINE ASAD_CINIT_loop

END Module messy_clamschem_asad_cinit
