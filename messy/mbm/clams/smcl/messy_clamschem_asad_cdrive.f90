Module messy_clamschem_asad_cdrive

contains

!    asad_cdrive  - chemistry driver routine
!
!    Glenn D. Carver & Paul D. Brown, Centre for Atmospheric Science,
!                                     University of Cambridge.
!
!
! Purpose: Main chemistry driver routine
!
!     Chemistry driver routine. If necessary, model concentrations are
!     converted from vmr to number density. If the chemistry is to be
!     "process-split" from the transport, the chemistry tendencies are
!     integrated by the chosen method and an average chemistry tendency
!     over the model timestep is returned. If the chemistry is not to
!     be integrated, instantaneous chemistry tendencies are returned.
!     The tracer tendencies returned to the calling routine are the
!     final values held in the chemistry.
!
!     Note: It is important for conservation that the concentrations
!     passed to this routine are non-negative.
!
!     Interface: Called from main program chem
!
!
!     Arguments:
!        cdot        - Tracer tendencies due to chemistry.
!        ftr         - Tracer concentrations.
!        pp          - Pressure (Nm-2).
!        pt          - Temperature (K).
!        pq          - Water vapor field (vmr).
!        nlev        - Model level
!        dryrt       - Dry deposition rates (s-1)
!        wetrt       - Wet deposition rates (s-1)
!        n_points    - No. of points calculations be done.
!
!     Method
!     ------
!     Since the heterogeneous rates may depend on species
!     concentrations, the call to hetero is made inside the
!     loop over chemistry sub-steps.
!
!     Photolysis rates may be computed at frequency set by the
!     variable nfphot and so the call to photol is also inside
!     the chemical sub-step loop.
!
!     Externals
!     ---------
!     asad_bimol  - Calculates bimolecular rate coefficients.
!     asad_trimol - Calculates trimolecular rate coefficients.
!     asad_photol - Calculates photolysis rate coefficients.
!     asad_drydep - Calculates dry deposition rates.
!     asad_wetdep - Calculates wet deposition rates.
!     asad_emissn - Calculates emission rates.
!     asad_totnud - Calculates total number densities.
!     asad_impact - IMPACT time integration scheme.
!     asad_hetero - Calculates heterogeneous rate coefficients.
!     asad_posthet- Performs housekeeping after heterogeneous chemistry.
!     asad_ftoy   - Partitions families.
!     asad_diffun - Calculates chemistry tendencies.
!
!     Local variables
!     ---------------
!     ifam       Family index of in/out species.
!     itr        Tracer index of in/out species.
!     iodd       Number of odd atoms in in/out species.
!     gfirst     .true. when the species need to be
!                initialised on the first chemical step.
!     gphot      .true. if photol needs to be called.
!     ghet       .true. if hetero needs to be called. (J.U. Grooss, 30.11.98)
!     hetdt      timestep for heterogenous chemistry (GG, 30.8.2000, module messy_clamschem_asad_mod_clams)
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
SUBROUTINE ASAD_CDRIVE(cdot, ftr, pp, pt, nlev, n_points)

  USE messy_clamschem_asad_mod
  USE messy_clamschem_asad_dummy,     ONLY: ereport, &
                                            wetdep, drydep
  USE messy_clamschem_asad_mod_clams, ONLY: lphotol, lhet, nfhet, hetdt, mype
  USE messy_clamschem_asad_bedriv,    ONLY: ASAD_BEDRIV
  USE messy_clamschem_asad_bimol,     ONLY: ASAD_BIMOL
  USE messy_clamschem_asad_diffun,    ONLY: ASAD_DIFFUN
  USE messy_clamschem_asad_ftoy,      ONLY: ASAD_FTOY
  USE messy_clamschem_asad_impact,    ONLY: ASAD_IMPACT
  USE messy_clamschem_asad_jac,       ONLY: ASAD_JAC
  USE messy_clamschem_asad_spmjpdriv, ONLY: ASAD_SPMJPDRIV
  USE messy_clamschem_asad_totnud,    ONLY: ASAD_TOTNUD
  USE messy_clamschem_asad_trimol,    ONLY: ASAD_TRIMOL
  USE messy_clamschem_asad_vodedriv,  ONLY: VODEDRIV
  USE messy_clams_global,     ONLY: prec
! op_pj_20170110+
!!$ USE messy_clamschem_global, ONLY: asad_gfirst
  USE messy_clams_global,     ONLY: asad_gfirst
! op_pj_20170110-
  USE messy_clamschem_emissn, ONLY: emissn
  USE messy_clamschem_hetero, ONLY: hetero
  USE messy_clamschem_photol, ONLY: photol

  IMPLICIT NONE

  ! Subroutine interface
  INTEGER, INTENT(IN) :: n_points              ! No of points
  INTEGER, INTENT(IN) :: nlev                  ! Model level

  REAL(PREC), INTENT(IN) :: pp(n_points)             ! Pressure
  REAL(PREC), INTENT(IN) :: pt(n_points)             ! Temperature

  REAL(PREC), INTENT(INOUT) :: ftr(n_points,jpctr)   ! Tracer concs
  REAL(PREC), INTENT(OUT)   :: cdot(n_points,jpctr)  ! Tracer tendencies

!       Local variables

  INTEGER :: errcode               ! Variable passed to ereport

  INTEGER :: jtr                                ! Loop variable
  INTEGER :: jl                                 ! Loop variable
  INTEGER :: js                                 ! Loop variable
  INTEGER :: nl
  INTEGER :: ifam
  INTEGER :: itr
  INTEGER :: iodd
  
  LOGICAL :: gfirst
  LOGICAL :: gphot, ghet

  CHARACTER(len=72) :: cmessage          ! Error message


!     ------------------------------------------------------------------
!       1.  Initialise variables and arrays

!       1.1   Clear tendencies to avoid contributions from levels
!             on which no chemistry is performed

        DO jtr = 1, jpctr
          DO jl = 1, n_points
            cdot(jl,jtr) = 0.0
          ENDDO
        ENDDO

!       1.2  Copy pressure and temperature to common

        DO jl = 1, n_points
          p(jl) = pp(jl)
          t(jl) = pt(jl)
        ENDDO


!     ------------------------------------------------------------------
!       2.  Calculate total number densities

        CALL asad_totnud(n_points)

!     ------------------------------------------------------------------
!       3.  Read model tracer concentrations into working array,
!           and if necessary, convert vmr to number densities

        IF ( lvmr ) THEN
          do jtr = 1, jpctr
            do jl  = 1, n_points
              ftr(jl,jtr) = ftr(jl,jtr) * tnd(jl)
              f(jl,jtr)   = ftr(jl,jtr)
            ENDDO
          ENDDO
        ELSE
          do jtr = 1, jpctr
            do jl  = 1, n_points
              f(jl,jtr)   = ftr(jl,jtr)
            ENDDO
          ENDDO
        ENDIF

!     ------------------------------------------------------------------
!       4.  Calculate reaction rate coefficients
!           --------- -------- ---- ------------

        CALL asad_bimol (n_points)

        CALL asad_trimol(n_points)


! Folgendes auskommentiert:
!          CALL ASAD_HETERO(n_points, cld_f, cld_l, rc_het)

!     ------------------------------------------------------------------
!       5.  Calculate deposition and emission rates
!           --------- ---------- --- -------- -----

        ! Routines WETDEP and DRYDEP: dummy routines 
        ! (messy_clamschem_asad_dummy.f90)
        IF ( ndepw /= 0 ) CALL WETDEP
        IF ( ndepd /= 0 ) CALL DRYDEP
        ! Routine emissn from package "chem" (messy_clamschem_emissn.f90)
        IF ( nemit /= 0 ) CALL emissn

!     ------------------------------------------------------------------
!       6.  Integrate chemistry by chosen method. Otherwise,
!           simply calculate tendencies due to chemistry
!           ------ --------- ---------- --- -- ---------

        if (ncsteps <= 0) then
           print*,'=== WARNING (from cdrive): No integration will be'
           print*,'=== carried out --> check parameter ncsteps!'
           print*,'=== ncsteps = ',ncsteps
        endif
        
        gphot = .true.
        IF ( method /= 0 ) THEN

          DO jsubs = 1, ncsteps

! ju_jug_150924  determination of gfirst changed back to previous version
!                gfirst is true at the beginning of a chem timestep, when 
!                rates need to be calculated
            gfirst = jsubs  ==  1
            !gfirst = asad_gfirst   ! true at the beginning of run and after mix
            asad_gfirst = .FALSE.
            if (mype==0) write (*,*) 'cdrive: gfirst=',gfirst
            gphot  = gfirst
            IF ( nfphot /= 0 .AND. .NOT.gfirst )                       &
                          gphot = mod(jsubs-1,nfphot) == 0

!!!!!!!! Folgende Zeilen ergaenzt:
            gphot = (gphot .and. lphotol)
            ghet  = (gfirst .and. lhet)

            ! Belegung von hetdt mit dtime=cdt*ncsteps (default)
         
            hetdt = cdt*ncsteps

            if ( nfhet /= 0 .and. .not. gfirst ) then
               ghet = (mod(jsubs-1,nfhet)==0) .and. lhet
               
               ! Berechnung von hetdt, wenn nfhet /= 0
               hetdt = cdt/float(nfhet)
            endif
!!!!!!!! Ende der hinzugefuegten Zeilen

!           ---------------------------------------------------
!           NON-STIFF integrators take values in the range 1-9.
!           ---------------------------------------------------

!           6.1  IMPACT integration: first compute heterogeneous
!                and photolysis rates, species and tendencies.
!                ===============================================

           nl = n_points
            IF ( method == 1 ) THEN
              ! Routine photol from package "chem" (messy_clamschem_photol.f90)
              IF ( gphot ) CALL photol
              ! Routine hetero from package "chem" (messy_clamschem_hetero.f90)
              IF ( ghet )  CALL hetero 
              CALL asad_ftoy( gfirst, nitfg, n_points )
              CALL asad_diffun( nl )
              CALL asad_jac( n_points )
              CALL asad_impact( gfirst, n_points )



!           6.2.  Quasi-steady state scheme.
!           ================================

            ELSEIF ( method == 2 ) THEN
              cmessage='QSSA not in UM6.5 build'
              errcode=1
              CALL EREPORT('ASAD_CDRIVE',errcode,cmessage)

!           6.3   Sparse Newton-Raphson solver
!           ==================================

            ELSEIF ( method == 3 ) THEN
              ! Routine photol from package "chem" (messy_clamschem_photol.f90)
              IF ( gphot ) CALL photol
              ! Routine hetero from package "chem" (messy_clamschem_hetero.f90)
              IF ( ghet )  CALL hetero
              CALL ASAD_FTOY( gfirst, nitfg, n_points )
               
              CALL ASAD_SPMJPDRIV(nlev,n_points)

!           6.5   Backward Euler solver
!           ===========================

            ELSEIF ( method == 5 ) then
              ! Routine photol from package "chem" (messy_clamschem_photol.f90)
              IF ( gphot ) CALL photol
              ! Routine hetero from package "chem" (messy_clamschem_hetero.f90)
              IF ( ghet )  CALL hetero
              CALL ASAD_FTOY( gfirst, nitfg, n_points )
              CALL ASAD_BEDRIV(nlev,mype,n_points)

!           -------------------------------------------------
!           STIFF integrators take values in the range 10-19.
!           -------------------------------------------------

!           6.10  NAG BDF stiff integrator.

            ELSEIF ( method == 10 ) THEN
              cmessage='NAG solver not build'
              errcode=1
              CALL EREPORT('ASAD_CDRIVE',errcode,cmessage)

!           6.11  SVODE ODE stiff integrator from NETLIB.

            ELSEIF ( method == 11 ) THEN
!!!!!
               if ( gphot ) call photol
               if ( ghet ) call hetero
               call asad_ftoy( gfirst, nitfg, n_points )
               call vodedriv( 22, n_points )

!!!!!
            ELSE
               print*,'=== WARNING (from cdrive): No integration will be'
               print*,'=== carried out --> check parameter method!'
               print*,'=== method = ',method
               
            ENDIF

!           6.12  Do any final work for the heterogeneous
!                 chemistry before the end of the time loop.
!                 =========================================


          ENDDO        ! End of looping over chemical timesteps

        ELSE           ! Method is equal to zero

!       6.99  Not integrating: just compute tendencies.

          method = 0
              ! Routine photol from package "chem" (messy_clamschem_photol.f90)
              IF ( gphot ) CALL photol
              ! Routine hetero from package "chem" (messy_clamschem_hetero.f90)
              IF ( ghet )  CALL hetero
          CALL asad_ftoy( .true., nit0, n_points )
          CALL asad_diffun( nl )

        ENDIF        ! End of IF statement for method


!     ------------------------------------------------------------------
!       7.  Determine concentrations and tendencies to be returned to
!           the model depending on whether or not the chemistry has
!           been integrated  -- -- ------- -- --- --- --------- ---
!           ---- ----------

!       7.1  Obtain model family concentrations by subtracting
!            concentrations of in/out species from ASAD families

        DO js = 1, jpspec
          IF ( ctype(js) == jpif ) THEN
            ifam = moffam(js)
            itr  = madvtr(js)
            iodd = nodd(js)
            DO jl = 1, n_points
              IF ( linfam(jl,itr) ) f(jl,ifam) =                       &
                                    f(jl,ifam) - iodd*f(jl,itr)
            ENDDO
          ENDIF
        ENDDO

!       7.2  Returned values of concentration and chemical tendency

        DO jtr = 1, jpctr
          IF ( method /= 0 ) THEN
            DO jl = 1, n_points
              cdot(jl,jtr) = ( f(jl,jtr)-ftr(jl,jtr)) / (cdt*ncsteps)
              ftr(jl,jtr)  = f(jl,jtr)
            ENDDO
          ELSE
            DO jl = 1, n_points
              cdot(jl,jtr) = fdot(jl,jtr)
              ftr(jl,jtr)  = f(jl,jtr)
            ENDDO
          ENDIF
        ENDDO


!     ------------------------------------------------------------------
!       8.  If necessary, convert from number densities back to vmr
!           -- ---------- ------- ---- ------ --------- ---- -- ---

        IF ( lvmr ) THEN
          DO jtr = 1, jpctr
            DO jl = 1, n_points
              ftr(jl,jtr)  = ftr(jl,jtr)  / tnd(jl)
              cdot(jl,jtr) = cdot(jl,jtr) / tnd(jl)
            ENDDO
          ENDDO
        ENDIF

        RETURN
        END SUBROUTINE ASAD_CDRIVE


      End Module messy_clamschem_asad_cdrive
