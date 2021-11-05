!
!  Description:
!    Module defining ASAD arrays, variables, and parameters
!
!  Adapted to CLaMS:
!   Use CLaMS definitions and routines instead of UKCA modules
!   
!   Jens-Uwe Grooss, Nicole Thomas, 2014-2018
!

MODULE messy_clamschem_asad_mod

USE messy_clams_global,             ONLY: prec
USE messy_clamschem_asad_mod_clams, ONLY: jpctr, jpspec, jpdw, jpdd, jpnr, jpbk, jptk, jphk, jppj


IMPLICIT NONE
SAVE
PUBLIC

! Interface section
INTERFACE asad_mod_init
  MODULE PROCEDURE asad_mod_init
END INTERFACE asad_mod_init

INTERFACE asad_mod_final
  MODULE PROCEDURE asad_mod_final
END INTERFACE asad_mod_final

! ju_nt_20140214: asad_mod_init_loop, asad_mod_final_loop added
INTERFACE asad_mod_init_loop
  MODULE PROCEDURE asad_mod_init_loop
END INTERFACE asad_mod_init_loop
INTERFACE asad_mod_final_loop
  MODULE PROCEDURE asad_mod_final_loop
END INTERFACE asad_mod_final_loop

REAL(PREC), ALLOCATABLE :: wp(:)           ! water vapour field
REAL(PREC), ALLOCATABLE :: dpd(:,:)
REAL(PREC), ALLOCATABLE :: dpw(:,:)
REAL(PREC), ALLOCATABLE :: emr(:,:)
REAL(PREC), ALLOCATABLE :: fj(:,:,:)       ! Full jacobian
REAL(PREC), ALLOCATABLE :: qa(:,:)
REAL(PREC), ALLOCATABLE :: ratio(:,:)
REAL(PREC), ALLOCATABLE :: p(:)            ! pressure in Pa (!)
REAL(PREC), ALLOCATABLE :: t(:)            ! Temperature in K
REAL(PREC), ALLOCATABLE :: t300(:)
REAL(PREC), ALLOCATABLE :: tnd(:)          ! total number density
REAL(PREC), ALLOCATABLE :: pmintnd(:)
REAL(PREC), ALLOCATABLE :: f(:,:)          ! tracer concentrations
REAL(PREC), ALLOCATABLE :: fdot(:,:)
REAL(PREC), TARGET, ALLOCATABLE :: pd(:,:)
REAL(PREC), POINTER     :: prod(:,:) 
REAL(PREC), POINTER     :: slos(:,:) 
REAL(PREC), ALLOCATABLE :: y(:,:)
REAL(PREC), ALLOCATABLE :: ydot(:,:)
REAL(PREC), ALLOCATABLE :: ej(:,:)
REAL(PREC), ALLOCATABLE :: rk(:,:)
REAL(PREC), ALLOCATABLE :: prk(:,:)
REAL(PREC), ALLOCATABLE :: deriv(:,:,:)
REAL(PREC), ALLOCATABLE :: ab(:,:)
REAL(PREC), ALLOCATABLE :: at(:,:)
REAL(PREC), ALLOCATABLE :: aj(:,:)
REAL(PREC), ALLOCATABLE :: ah(:,:)
REAL(PREC), ALLOCATABLE :: ztabpd(:,:)
REAL(PREC), ALLOCATABLE :: spfj(:,:)       ! Sparse full Jacobian

INTEGER, ALLOCATABLE :: madvtr(:)
INTEGER, ALLOCATABLE :: majors(:)
INTEGER, ALLOCATABLE :: moffam(:)
INTEGER, ALLOCATABLE :: nodd(:)
INTEGER, ALLOCATABLE :: nltrf(:)
INTEGER, ALLOCATABLE :: nltr3(:)
INTEGER, ALLOCATABLE :: nltrim(:,:)
INTEGER, ALLOCATABLE :: nlpdv(:,:)
INTEGER, ALLOCATABLE :: nfrpx(:)     ! index to fractional product array nfrpx:
! one entry for each reaction. If zero, there are no fractional products. If nonzero,
! contains the array element in frpx for the first coefficient for that reaction.
INTEGER, ALLOCATABLE :: ntabfp(:,:)  ! Table used for indexing the fractional
!    products:  ntabfp(i,1) contains the species no.,
!               ntabfp(i,2) contains the reaction no., and
!               ntabfp(i,3) contains the array location in the frpx array.
INTEGER, ALLOCATABLE :: ntabpd(:,:)
INTEGER, ALLOCATABLE :: npdfr(:,:)
INTEGER, ALLOCATABLE :: ngrp(:,:)
INTEGER, ALLOCATABLE :: njcgrp(:,:)
INTEGER, ALLOCATABLE :: nprdx3(:,:,:)
INTEGER, ALLOCATABLE :: nprdx2(:,:)
INTEGER, ALLOCATABLE :: nprdx1(:)
INTEGER, ALLOCATABLE :: njacx3(:,:,:)
INTEGER, ALLOCATABLE :: njacx2(:,:)
INTEGER, ALLOCATABLE :: njacx1(:)
INTEGER, ALLOCATABLE :: nmpjac(:)
INTEGER, ALLOCATABLE :: npjac1(:,:)
INTEGER, ALLOCATABLE :: nbrkx(:)
INTEGER, ALLOCATABLE :: ntrkx(:)
INTEGER, ALLOCATABLE :: nprkx(:)
INTEGER, ALLOCATABLE :: nhrkx(:)
INTEGER, ALLOCATABLE :: nlall(:)
INTEGER, ALLOCATABLE :: nlstst(:)
INTEGER, ALLOCATABLE :: nlf(:)
INTEGER, ALLOCATABLE :: nlmajmin(:)
INTEGER, ALLOCATABLE :: nldepd(:)
INTEGER, ALLOCATABLE :: nldepw(:)
INTEGER, ALLOCATABLE :: nlemit(:)
INTEGER, ALLOCATABLE :: nldepx(:)
INTEGER, ALLOCATABLE :: njcoth(:,:)
INTEGER, ALLOCATABLE :: nmzjac(:)
INTEGER, ALLOCATABLE :: nzjac1(:,:)
INTEGER, ALLOCATABLE :: njcoss(:,:)
INTEGER, ALLOCATABLE :: nmsjac(:)
INTEGER, ALLOCATABLE :: nsjac1(:,:)
INTEGER, ALLOCATABLE :: nsspt(:)
INTEGER, ALLOCATABLE :: nspi(:,:)
INTEGER, ALLOCATABLE :: nsspi(:,:)
INTEGER, ALLOCATABLE :: nssi(:)
INTEGER, ALLOCATABLE :: nssrt(:)
INTEGER, ALLOCATABLE :: nssri(:,:)
INTEGER, ALLOCATABLE :: nssrx(:,:)
REAL(PREC), ALLOCATABLE :: frpb(:)         ! fractional product array (bimol)
REAL(PREC), ALLOCATABLE :: frpt(:)         ! fractional product array (trimol)
REAL(PREC), ALLOCATABLE :: frpj(:)         ! fractional product array (phot)
REAL(PREC), ALLOCATABLE :: frph(:)         ! fractional product array (het)
REAL(PREC), ALLOCATABLE :: frpx(:)         ! fractional product array (total)
! sparse algebra
INTEGER, ALLOCATABLE :: pntr(:,:)  ! Map of nonzero entries
INTEGER, ALLOCATABLE :: pntr1(:,:) ! modified map (after decomposition)
INTEGER, ALLOCATABLE :: pntr2(:,:) ! Map of nonzero entries, before reordering
INTEGER, ALLOCATABLE :: ro(:)         ! reordering of tracers to minimize fill-in

INTEGER, ALLOCATABLE :: ilcf(:)
INTEGER, ALLOCATABLE :: ilss(:)
INTEGER, ALLOCATABLE :: ilct(:)
INTEGER, ALLOCATABLE :: ilftr(:)
INTEGER, ALLOCATABLE :: ilft(:)
INTEGER, ALLOCATABLE :: ilstmin(:)

LOGICAL, ALLOCATABLE :: linfam(:,:)
LOGICAL, ALLOCATABLE :: ldepd(:)     ! T for dry deposition
LOGICAL, ALLOCATABLE :: ldepw(:)     ! T for wet deposition
LOGICAL, ALLOCATABLE :: lemit(:)     ! T for emission

LOGICAL, ALLOCATABLE :: converged(:) ! arrays for spimpmjp solver
LOGICAL, ALLOCATABLE :: positive(:)

CHARACTER(len=10), ALLOCATABLE :: advt(:)       ! advected tracers
CHARACTER(len=10), ALLOCATABLE :: nadvt(:)      ! non-advected species 
CHARACTER(len=10), ALLOCATABLE :: family(:)     ! family
CHARACTER(len=10), ALLOCATABLE :: speci(:)      ! species names
CHARACTER(len=2),  ALLOCATABLE :: ctype(:)      ! species type
CHARACTER(len=10), ALLOCATABLE :: spb(:,:)      ! species from bimolecular rates
CHARACTER(len=10), ALLOCATABLE :: spt(:,:)      ! species from termolecular rates
CHARACTER(len=10), ALLOCATABLE :: spj(:,:)      ! species from photolysis rates
CHARACTER(len=10), ALLOCATABLE :: sph(:,:)      ! species from heterogenous rates

! ju_nt_20140304: pmin and ptol changed by JUG to the following (June 2006)
! ju_jug_20150223: pmin and ptol decreased to 1.0e-15/1.0e-4, optimized for
!                  NR-solver with individual convergence check.
!                  was before 1.0e-20/1.0e-5  

REAL(PREC), PARAMETER    :: pmin = 1.0e-15
REAL(PREC), PARAMETER    :: ptol = 1.0e-4            ! tolerance for time integration
REAL(PREC), PARAMETER    :: ftol = 1.0e-3            ! tolerance in family member iteration
REAL(PREC), PARAMETER    :: fmax = 1.E+18            ! maximum concentration without divergence

INTEGER, PARAMETER :: kfphot=0
INTEGER, PARAMETER :: jpss = 16
INTEGER, PARAMETER :: jpssr = 51
INTEGER, PARAMETER :: jpeq = 2      ! dimension for dissociation arrays 

INTEGER, PARAMETER :: jpfrpd=100
INTEGER, PARAMETER :: jpab    = 3
INTEGER, PARAMETER :: jpat    = 7
INTEGER, PARAMETER :: jpaj    = 3
INTEGER, PARAMETER :: jpah    = 3
INTEGER, PARAMETER :: jpspb   = 6
INTEGER, PARAMETER :: jpspt   = 4
INTEGER, PARAMETER :: jpspj   = 6
INTEGER, PARAMETER :: jpsph   = 6
INTEGER, PARAMETER :: jpmsp   = jpspb
INTEGER, PARAMETER :: jppjac  = 10


! Fractional product parameters - these are initialsed to jp values in asad_mod_init
INTEGER :: jpfrpb
INTEGER :: jpfrpt
INTEGER :: jpfrpj
INTEGER :: jpfrph
INTEGER :: jpfrpx

INTEGER, PARAMETER :: spfjsize_max =1000 ! maximum number of
                                         ! nonzero matrix elements

CHARACTER(LEN=2),  PARAMETER :: jpfm = 'FM'     ! Family member
CHARACTER(LEN=2),  PARAMETER :: jpif = 'FT'     ! Family member depending in timestep
CHARACTER(LEN=2),  PARAMETER :: jpsp = 'TR'     ! Independent tracer
CHARACTER(LEN=2),  PARAMETER :: jpna = 'SS'     ! Steady-state species
CHARACTER(LEN=2),  PARAMETER :: jpco = 'CT'     ! Constant
CHARACTER(LEN=2),  PARAMETER :: jpcf = 'CF'     ! Constant with spatial field

LOGICAL, PARAMETER :: lvmr=.true.    ! T for volume mixing ratio
LOGICAL      :: o1d_in_ss      ! T for steady state,
LOGICAL      :: o3p_in_ss      ! these are set in routine:
LOGICAL      :: n_in_ss        ! asad_mod_init
LOGICAL      :: h_in_ss        ! 
INTEGER, PARAMETER :: nss_o1d=1      ! indicies of deriv array
INTEGER, PARAMETER :: nss_o3p=2      !    "         "
INTEGER, PARAMETER :: nss_n=3        !    "         "
INTEGER, PARAMETER :: nss_h=4        !    "         "

REAL(PREC) :: fch4,fco2,fh2,fn2,fo2
REAL(PREC)    :: cdt                       ! chemistry timestep
REAL(PREC)    :: peps                      ! 
REAL(PREC)    :: dtime                     ! model timestep
REAL(PREC), PARAMETER :: tslimit = 1200.0  ! timestep limit for some solvers

INTEGER, PARAMETER :: kcdt = 3600    ! timestep for N-R solver !! MAY NOT WANT HERE
INTEGER :: nrsteps
INTEGER :: nitnr          ! Iterations in ftoy for IMPACT solver
INTEGER :: nitfg          ! Max no of iterations in ftoy
INTEGER :: ntrf           ! Counter for tracers
INTEGER :: ntr3
INTEGER :: nnaf           ! Counter for non-advected tracers
INTEGER :: nuni
INTEGER :: nsst           ! No of steady-state species
INTEGER :: ncsteps        ! No of chemical steps
INTEGER :: nit0=20        ! ftoy iterations with method=0
INTEGER :: nfphot
INTEGER :: jsubs
INTEGER :: jlst
INTEGER, PARAMETER :: int_method_impact = 1
INTEGER, PARAMETER :: int_method_NR = 3
INTEGER, PARAMETER :: int_method_BE = 5
INTEGER, PARAMETER :: int_method_BE_explicit = 10
INTEGER :: method         ! chemistry integration method
INTEGER :: interval       ! interval in timesteps between calls to chemistry
INTEGER :: nnfrp          ! Total number of fractional products
INTEGER :: nstst          ! No of steady state species
INTEGER :: nf
INTEGER :: ndepd          ! No of dry deposited species
INTEGER :: ndepw          ! No of wet deposited species
INTEGER :: nemit          ! No of emitted species
INTEGER :: ntro3, ntroh, ntrho2, ntrno 
INTEGER :: nspo1d, nspo3p, nspo3, nspoh
INTEGER :: nspho2, nspno, nspn, nsph
INTEGER :: ih_o3, ih_h2o2, ih_so2, ih_hno3  ! index for soluble species
INTEGER :: ihso3_h2o2                       ! Index for HSO3- + H2O2(aq) reaction
INTEGER :: ihso3_o3                         ! Index for HSO3- + O3(aq) reaction
INTEGER :: iso3_o3                          ! Index for SO3-- + O3(aq) reaction
INTEGER :: iso2_oh                          ! Index for SO2 + OH reaction
INTEGER :: ih2o2_oh                         ! Index for H2O2 + OH reaction
INTEGER :: ihno3_oh                         ! Index for HNO3 + OH reaction
INTEGER :: in2o5_h                          ! Index for N2O5 => HONO2 heterog. reaction
INTEGER :: iho2_h                           ! Index for HO2 + HO2 => H2O2 heterog. "

LOGICAL :: lsvjac         ! Flag for saving jacobian if not recalculated
LOGICAL :: ljacx

CONTAINS

! ######################################################################
SUBROUTINE asad_mod_init

! To allocate and initialise ASAD arrays and variables

! ju_nt_20140324
USE messy_clamschem_defs_mod,        ONLY: chch_t, chch_defs
USE messy_clamschem_asad_dummy,      ONLY: ereport
USE messy_clamschem_asad_mod_clams,  ONLY: theta_field_size, model_levels

IMPLICIT NONE

INTEGER :: i
INTEGER :: errcode
CHARACTER(LEN=72) :: cmessage

LOGICAL, SAVE :: firstcall=.true.

! Fractional product parameters - set total number allowed to be
! total number of reactions * 2 as each fractional product has 4 potentials!!
jpfrpb  = (jpspb-2)*jpbk
jpfrpt  = (jpspt-2)*jptk   
jpfrpj  = (jpspj-2)*jppj   
jpfrph  = (jpsph-2)*jphk   
jpfrpx  = jpfrpb + jpfrpt + jpfrpj + jpfrph  ! Total FPs


IF (.NOT. ALLOCATED(madvtr)) ALLOCATE(madvtr(jpspec))
IF (.NOT. ALLOCATED(majors)) ALLOCATE(majors(jpctr))
IF (.NOT. ALLOCATED(moffam)) ALLOCATE(moffam(jpspec))
IF (.NOT. ALLOCATED(nodd)) ALLOCATE(nodd(jpspec))
IF (.NOT. ALLOCATED(nltrf)) ALLOCATE(nltrf(jpctr))
IF (.NOT. ALLOCATED(nltr3)) ALLOCATE(nltr3(jpctr))
IF (.NOT. ALLOCATED(advt)) ALLOCATE(advt(jpctr))
IF (.NOT. ALLOCATED(nadvt)) ALLOCATE(nadvt(jpspec-jpctr))
IF (.NOT. ALLOCATED(family)) ALLOCATE(family(jpspec))
IF (.NOT. ALLOCATED(speci)) ALLOCATE(speci(jpspec))
IF (.NOT. ALLOCATED(ctype)) ALLOCATE(ctype(jpspec))
IF (.NOT. ALLOCATED(nspi)) ALLOCATE(nspi(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nsspt)) ALLOCATE(nsspt(jpss))
IF (.NOT. ALLOCATED(nsspi)) ALLOCATE(nsspi(jpss,jpssr))
IF (.NOT. ALLOCATED(nssi)) ALLOCATE(nssi(jpss))
IF (.NOT. ALLOCATED(nssrt)) ALLOCATE(nssrt(jpss))
IF (.NOT. ALLOCATED(nssri)) ALLOCATE(nssri(jpss,jpssr))
IF (.NOT. ALLOCATED(nssrx)) ALLOCATE(nssrx(jpss,jpssr)) 
IF (.NOT. ALLOCATED(ldepd)) ALLOCATE(ldepd(jpspec))
IF (.NOT. ALLOCATED(ldepw)) ALLOCATE(ldepw(jpspec))
IF (.NOT. ALLOCATED(lemit)) ALLOCATE(lemit(jpspec))
IF (.NOT. ALLOCATED(nltrim)) ALLOCATE(nltrim(0:jpctr,3))
IF (.NOT. ALLOCATED(nlpdv)) ALLOCATE(nlpdv((jpspj-2)*jppj,2))
IF (.NOT. ALLOCATED(ab)) ALLOCATE(ab(jpbk+1,jpab))
IF (.NOT. ALLOCATED(at)) ALLOCATE(at(jptk+1,jpat))
IF (.NOT. ALLOCATED(aj)) ALLOCATE(aj(jppj+1,jpaj))
IF (.NOT. ALLOCATED(ah)) ALLOCATE(ah(jphk+1,jpah))
IF (.NOT. ALLOCATED(spb)) ALLOCATE(spb(jpbk+1,jpspb))
IF (.NOT. ALLOCATED(spt)) ALLOCATE(spt(jptk+1,jpspt))
IF (.NOT. ALLOCATED(spj)) ALLOCATE(spj(jppj+1,jpspj))
IF (.NOT. ALLOCATED(sph)) ALLOCATE(sph(jphk+1,jpsph))
IF (.NOT. ALLOCATED(frpb)) ALLOCATE(frpb(jpfrpb))
IF (.NOT. ALLOCATED(frpt)) ALLOCATE(frpt(jpfrpt))
IF (.NOT. ALLOCATED(frpj)) ALLOCATE(frpj(jpfrpj))
IF (.NOT. ALLOCATED(frph)) ALLOCATE(frph(jpfrph))
IF (.NOT. ALLOCATED(frpx)) ALLOCATE(frpx(jpfrpx))
IF (.NOT. ALLOCATED(ztabpd)) ALLOCATE(ztabpd(jpfrpd,2))
IF (.NOT. ALLOCATED(nfrpx)) ALLOCATE(nfrpx(jpnr))
IF (.NOT. ALLOCATED(ntabfp)) ALLOCATE(ntabfp(jpfrpx,3))
IF (.NOT. ALLOCATED(ntabpd)) ALLOCATE(ntabpd(jpfrpd,3))
IF (.NOT. ALLOCATED(npdfr)) ALLOCATE(npdfr(jpnr,2))
IF (.NOT. ALLOCATED(ngrp)) ALLOCATE(ngrp(2*jpspec,3))
IF (.NOT. ALLOCATED(njcgrp)) ALLOCATE(njcgrp(jpctr,3))
IF (.NOT. ALLOCATED(nprdx3)) ALLOCATE(nprdx3(3,(jpnr/(3*3))+3*3,2*jpspec))
IF (.NOT. ALLOCATED(nprdx2)) ALLOCATE(nprdx2(2,2*jpspec))
IF (.NOT. ALLOCATED(nprdx1)) ALLOCATE(nprdx1(2*jpspec))
IF (.NOT. ALLOCATED(njacx3)) ALLOCATE(njacx3(3,(jpnr/(3*3))+3*3,jpctr))
IF (.NOT. ALLOCATED(njacx2)) ALLOCATE(njacx2(2,jpctr))
IF (.NOT. ALLOCATED(njacx1)) ALLOCATE(njacx1(jpctr))
IF (.NOT. ALLOCATED(nmpjac)) ALLOCATE(nmpjac(jpctr))
IF (.NOT. ALLOCATED(npjac1)) ALLOCATE(npjac1(jppjac,jpctr))
IF (.NOT. ALLOCATED(nbrkx)) ALLOCATE(nbrkx(jpbk+1))
IF (.NOT. ALLOCATED(ntrkx)) ALLOCATE(ntrkx(jptk+1))
IF (.NOT. ALLOCATED(nprkx)) ALLOCATE(nprkx(jppj+1))
IF (.NOT. ALLOCATED(nhrkx)) ALLOCATE(nhrkx(jphk+1))
IF (.NOT. ALLOCATED(nlall)) ALLOCATE(nlall(jpspec))
IF (.NOT. ALLOCATED(nlstst)) ALLOCATE(nlstst(jpspec))
IF (.NOT. ALLOCATED(nlf)) ALLOCATE(nlf(jpspec))
IF (.NOT. ALLOCATED(nlmajmin)) ALLOCATE(nlmajmin(jpspec))
IF (.NOT. ALLOCATED(nldepd)) ALLOCATE(nldepd(jpspec))
IF (.NOT. ALLOCATED(nldepw)) ALLOCATE(nldepw(jpspec))
IF (.NOT. ALLOCATED(nlemit)) ALLOCATE(nlemit(jpspec))
IF (.NOT. ALLOCATED(nldepx)) ALLOCATE(nldepx(jpspec))
IF (.NOT. ALLOCATED(njcoth)) ALLOCATE(njcoth(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nmzjac)) ALLOCATE(nmzjac(jpctr))
IF (.NOT. ALLOCATED(nzjac1)) ALLOCATE(nzjac1(jpnr,jpctr))
IF (.NOT. ALLOCATED(njcoss)) ALLOCATE(njcoss(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nmsjac)) ALLOCATE(nmsjac(jpctr))
IF (.NOT. ALLOCATED(nsjac1)) ALLOCATE(nsjac1(jpnr,jpctr))

! the following had save attribs.
IF (.NOT. ALLOCATED(ilcf)) ALLOCATE(ilcf(jpspec))
IF (.NOT. ALLOCATED(ilss)) ALLOCATE(ilss(jpspec))
IF (.NOT. ALLOCATED(ilct)) ALLOCATE(ilct(jpspec))
IF (.NOT. ALLOCATED(ilftr)) ALLOCATE(ilftr(jpspec))
IF (.NOT. ALLOCATED(ilft)) ALLOCATE(ilft(jpspec))
IF (.NOT. ALLOCATED(ilstmin)) ALLOCATE(ilstmin(jpspec))

! Set integration method (1 = IMPACT; 3 = N-R solver; 5 = Backward-Euler)
! ju_nt_20140124: "method" is set in "chem"

! Initialize variables that may be changed in cinit
nrsteps = 45            ! No of N-R steps, set > 50 to debug convergence failures
nitnr   = 10            ! Iterations in ftoy
nitfg   = 10            ! Max number of iterations in ftoy

! Initialise arrays
njcoth(:,:) = 0
dtime   = 1200.0         ! model timestep, reset in chemistry_ctl

! ju_nt_20140214: moved to asad_mod_init_loop:

IF (method == int_method_NR) THEN
  IF (.NOT. ALLOCATED(pntr))  ALLOCATE(pntr(jpctr, jpctr))
  IF (.NOT. ALLOCATED(pntr1))  ALLOCATE(pntr1(jpctr, jpctr))
  IF (.NOT. ALLOCATED(pntr2))  ALLOCATE(pntr2(jpctr, jpctr))
  IF (.NOT. ALLOCATED(ro))  ALLOCATE(ro(jpctr))
END IF

! Find out which species are in steady state (for N-R solver)
o1d_in_ss = .FALSE.
o3p_in_ss = .FALSE.
n_in_ss   = .FALSE.
h_in_ss   = .FALSE.
DO i=1,jpspec
 IF (chch_defs(i)%speci=='O(1D)     ' .AND.                             &
     chch_defs(i)%ctype(1:2)==jpna) o1d_in_ss=.TRUE.
 IF (chch_defs(i)%speci=='O(3P)     ' .AND.                             &
     chch_defs(i)%ctype(1:2)==jpna) o3p_in_ss=.TRUE.
 IF (chch_defs(i)%speci=='N         ' .AND.                             &
     chch_defs(i)%ctype(1:2)==jpna) n_in_ss=.TRUE.
 IF (chch_defs(i)%speci=='H         ' .AND.                             &
     chch_defs(i)%ctype(1:2)==jpna) h_in_ss=.TRUE.
ENDDO

RETURN
END SUBROUTINE asad_mod_init

! ######################################################################
SUBROUTINE asad_mod_init_loop

! To allocate and initialise ASAD arrays and variables

USE messy_clamschem_asad_mod_clams,             ONLY: theta_field_size, model_levels

IMPLICIT NONE

LOGICAL, SAVE :: firstcall=.true.

! nullify prod and slos on firstcall to give DISSASSOCIATED attribute
IF (firstcall) THEN
   nullify(prod)
   nullify(slos)
   firstcall = .false.
end IF

! pd is a TARGET
IF (.NOT. ALLOCATED(pd)) ALLOCATE(pd(theta_field_size,2*jpspec))

IF (.NOT. ALLOCATED(wp)) ALLOCATE(wp(theta_field_size))
IF (.NOT. ALLOCATED(linfam)) ALLOCATE(linfam(theta_field_size,0:jpctr))

IF (.NOT. ALLOCATED(dpd)) ALLOCATE(dpd(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(dpw)) ALLOCATE(dpw(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(emr)) ALLOCATE(emr(theta_field_size,jpspec))

IF (.NOT. ALLOCATED(fj)) ALLOCATE(fj(theta_field_size,jpctr,jpctr))
IF (.NOT. ALLOCATED(qa)) ALLOCATE(qa(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(ratio)) ALLOCATE(ratio(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(p)) ALLOCATE(p(theta_field_size))
IF (.NOT. ALLOCATED(t)) ALLOCATE(t(theta_field_size))
IF (.NOT. ALLOCATED(t300)) ALLOCATE(t300(theta_field_size))
IF (.NOT. ALLOCATED(tnd)) ALLOCATE(tnd(theta_field_size))
IF (.NOT. ALLOCATED(pmintnd)) ALLOCATE(pmintnd(theta_field_size))
IF (.NOT. ALLOCATED(f)) ALLOCATE(f(theta_field_size,jpctr))
IF (.NOT. ALLOCATED(fdot)) ALLOCATE(fdot(theta_field_size,jpctr))
IF (.NOT. ALLOCATED(y)) ALLOCATE(y(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(ydot)) ALLOCATE(ydot(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(ej)) ALLOCATE(ej(theta_field_size,jpctr))
IF (.NOT. ALLOCATED(rk)) ALLOCATE(rk(theta_field_size,jpnr))
IF (.NOT. ALLOCATED(prk)) ALLOCATE(prk(theta_field_size,jpnr))
IF (.NOT. ALLOCATED(deriv)) ALLOCATE(deriv(theta_field_size,4,4))


! EQUIVALENCE ( pd(1,1), prod(1,1) )
! EQUIVALENCE ( pd(1,jpspec+1), slos(1,1) )
prod => pd(:,1:jpspec)
slos => pd(:,jpspec+1:2*jpspec)

deriv(:,:,:) = 1.0     ! Temp fix for deriv being uninitialised in first
                       ! solver iteration 

IF (method == int_method_NR) THEN
  IF (.NOT. ALLOCATED(spfj))  ALLOCATE(spfj(theta_field_size,spfjsize_max))
END IF

RETURN
END SUBROUTINE asad_mod_init_loop

! ######################################################################
SUBROUTINE asad_mod_final_loop

! To deallocate ASAD arrays

IMPLICIT NONE

! pd is a TARGET
IF (ALLOCATED(pd)) DEALLOCATE(pd)
nullify(prod)
nullify(slos)

!!!!!
! ju_nt_20140212
! following deallocates added:
IF (ALLOCATED(wp)) DEALLOCATE(wp)
IF (ALLOCATED(linfam)) DEALLOCATE(linfam)

IF (ALLOCATED(dpd)) DEALLOCATE(dpd)
IF (ALLOCATED(dpw)) DEALLOCATE(dpw)
IF (ALLOCATED(emr)) DEALLOCATE(emr)

IF (ALLOCATED(fj)) DEALLOCATE(fj)
IF (ALLOCATED(qa)) DEALLOCATE(qa)
IF (ALLOCATED(ratio)) DEALLOCATE(ratio)
IF (ALLOCATED(p)) DEALLOCATE(p)
IF (ALLOCATED(t)) DEALLOCATE(t)
IF (ALLOCATED(t300)) DEALLOCATE(t300)
IF (ALLOCATED(tnd)) DEALLOCATE(tnd)
IF (ALLOCATED(pmintnd)) DEALLOCATE(pmintnd)
IF (ALLOCATED(f)) DEALLOCATE(f)
IF (ALLOCATED(fdot)) DEALLOCATE(fdot)
IF (ALLOCATED(y)) DEALLOCATE(y)
IF (ALLOCATED(ydot)) DEALLOCATE(ydot)
IF (ALLOCATED(ej)) DEALLOCATE(ej)
IF (ALLOCATED(rk)) DEALLOCATE(rk)
IF (ALLOCATED(prk)) DEALLOCATE(prk)
IF (ALLOCATED(deriv)) DEALLOCATE(deriv)


IF (method == int_method_NR) THEN       ! sparse vars
   IF (ALLOCATED(spfj))   DEALLOCATE(spfj)
ENDIF

RETURN
END SUBROUTINE asad_mod_final_loop

! ######################################################################
SUBROUTINE asad_mod_final

! To deallocate ASAD arrays

IMPLICIT NONE

IF (ALLOCATED(madvtr)) DEALLOCATE(madvtr)
IF (ALLOCATED(majors)) DEALLOCATE(majors)
IF (ALLOCATED(moffam)) DEALLOCATE(moffam)
IF (ALLOCATED(nodd)) DEALLOCATE(nodd)

IF (ALLOCATED(nltrf)) DEALLOCATE(nltrf)
IF (ALLOCATED(nltr3)) DEALLOCATE(nltr3)
IF (ALLOCATED(advt)) DEALLOCATE(advt)
IF (ALLOCATED(nadvt)) DEALLOCATE(nadvt)
IF (ALLOCATED(family)) DEALLOCATE(family)
IF (ALLOCATED(speci)) DEALLOCATE(speci)
IF (ALLOCATED(ctype)) DEALLOCATE(ctype)
IF (ALLOCATED(nspi)) DEALLOCATE(nspi)
IF (ALLOCATED(nsspt)) DEALLOCATE(nsspt)
IF (ALLOCATED(nsspi)) DEALLOCATE(nsspi)
IF (ALLOCATED(nssi)) DEALLOCATE(nssi)
IF (ALLOCATED(nssrt)) DEALLOCATE(nssrt)
IF (ALLOCATED(nssri)) DEALLOCATE(nssri)
IF (ALLOCATED(nssrx)) DEALLOCATE(nssrx)
IF (ALLOCATED(ldepd)) DEALLOCATE(ldepd)
IF (ALLOCATED(ldepw)) DEALLOCATE(ldepw)
IF (ALLOCATED(lemit)) DEALLOCATE(lemit)
IF (ALLOCATED(nltrim)) DEALLOCATE(nltrim)
IF (ALLOCATED(nlpdv)) DEALLOCATE(nlpdv)
IF (ALLOCATED(ab)) DEALLOCATE(ab)
IF (ALLOCATED(at)) DEALLOCATE(at)
IF (ALLOCATED(aj)) DEALLOCATE(aj)
IF (ALLOCATED(ah)) DEALLOCATE(ah)
IF (ALLOCATED(spb)) DEALLOCATE(spb)
IF (ALLOCATED(spt)) DEALLOCATE(spt)
IF (ALLOCATED(spj)) DEALLOCATE(spj)
IF (ALLOCATED(sph)) DEALLOCATE(sph)
IF (ALLOCATED(frpb)) DEALLOCATE(frpb)
IF (ALLOCATED(frpt)) DEALLOCATE(frpt)
IF (ALLOCATED(frpj)) DEALLOCATE(frpj)
IF (ALLOCATED(frph)) DEALLOCATE(frph)
IF (ALLOCATED(frpx)) DEALLOCATE(frpx)
IF (ALLOCATED(ztabpd)) DEALLOCATE(ztabpd)
IF (ALLOCATED(nfrpx)) DEALLOCATE(nfrpx)
IF (ALLOCATED(ntabfp)) DEALLOCATE(ntabfp)
IF (ALLOCATED(ntabpd)) DEALLOCATE(ntabpd)
IF (ALLOCATED(npdfr)) DEALLOCATE(npdfr)
IF (ALLOCATED(ngrp)) DEALLOCATE(ngrp)
IF (ALLOCATED(njcgrp)) DEALLOCATE(njcgrp)
IF (ALLOCATED(nprdx3)) DEALLOCATE(nprdx3)
IF (ALLOCATED(nprdx2)) DEALLOCATE(nprdx2)
IF (ALLOCATED(nprdx1)) DEALLOCATE(nprdx1)
IF (ALLOCATED(njacx3)) DEALLOCATE(njacx3)
IF (ALLOCATED(njacx2)) DEALLOCATE(njacx2)
IF (ALLOCATED(njacx1)) DEALLOCATE(njacx1)
IF (ALLOCATED(nmpjac)) DEALLOCATE(nmpjac)
IF (ALLOCATED(npjac1)) DEALLOCATE(npjac1)
IF (ALLOCATED(nbrkx)) DEALLOCATE(nbrkx)
IF (ALLOCATED(ntrkx)) DEALLOCATE(ntrkx)
IF (ALLOCATED(nprkx)) DEALLOCATE(nprkx)
IF (ALLOCATED(nhrkx)) DEALLOCATE(nhrkx)
IF (ALLOCATED(nlall)) DEALLOCATE(nlall)
IF (ALLOCATED(nlstst)) DEALLOCATE(nlstst)
IF (ALLOCATED(nlf)) DEALLOCATE(nlf)
IF (ALLOCATED(nlmajmin)) DEALLOCATE(nlmajmin)
IF (ALLOCATED(nldepd)) DEALLOCATE(nldepd)
IF (ALLOCATED(nldepw)) DEALLOCATE(nldepw)
IF (ALLOCATED(nlemit)) DEALLOCATE(nlemit)
IF (ALLOCATED(nldepx)) DEALLOCATE(nldepx)
IF (ALLOCATED(njcoth)) DEALLOCATE(njcoth)
IF (ALLOCATED(nmzjac)) DEALLOCATE(nmzjac)
IF (ALLOCATED(nzjac1)) DEALLOCATE(nzjac1)
IF (ALLOCATED(njcoss)) DEALLOCATE(njcoss)
IF (ALLOCATED(nmsjac)) DEALLOCATE(nmsjac)
IF (ALLOCATED(nsjac1)) DEALLOCATE(nsjac1)


IF (method == int_method_NR) THEN       ! sparse vars
   IF (ALLOCATED(pntr))   DEALLOCATE(pntr)
   IF (ALLOCATED(pntr1))   DEALLOCATE(pntr1)
   IF (ALLOCATED(pntr2))   DEALLOCATE(pntr2)
   IF (ALLOCATED(ro))   DEALLOCATE(ro)
ENDIF


RETURN
END SUBROUTINE asad_mod_final

END MODULE messy_clamschem_asad_mod
