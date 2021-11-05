Module messy_clamschem_asad_vodedriv

contains

     subroutine vodedriv( kmf, n_points )
!
!     VODEDRIV  - Driver for SVODE stiff/non-stiff ODE integrator.
!
!     Glenn Carver      Centre for Atmospheric Science
!                       University of Cambridge
!
!     ASAD: vodedriv                  Version: vodedriv.f 4.1 01/15/97
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the SVODE stiff/non-stiff integrator.
!
!     Interface
!     ---------
!
!        kmf  - Method switch. Used to select stiff (mf=22) or 
!               non-stiff (Adams method, mf=10) integrator.
!
!     Called from chemistry driver routine *cdrive*.
!
!     This routine assumes that all the bi-,tri-, phot- and het-
!     reaction rates have been computed prior to this routine.
!     It is also assumed that the species array, y, has been set
!     from the array, f, passed to ASAD, and that constant species
!     have been set. This can be done by calling the routine fyinit.
!
!     Method.
!     -------
!     See the comments at the top of the svode.f file for full details
!     of the SVODE code and the reference therein. SVODE is a 
!     variable coefficient ODE solver, available from the NETLIB WWW
!     site as part of the ODE (not ODEPACK) package. It contains both
!     stiff and nonstiff solvers. The stiff solver has been found to
!     work very well and compares very well with the results from the
!     NAG D02EAF routine whilst executing in a quarter of the time for
!     our test problem. SVODE can cope with very stiff systems such as
!     when integrating O(1D).  The computational cost of SVODE is still
!     significant however.
!
!     The SVODE integrator is passed a single gridpoint at a time.
!     The tolerances set by default below have been found to work well
!     for a test problem using a detailed stratospheric chemistry 
!     scheme.
!
!     This routine will use either the stiff or non-stiff (Adams
!     method) integrators built into the SVODE code depending on the
!     value of the argument, kmf.
!
!     Note that we must reset the ISTATE argument before each call
!     since we are passing different sets of y each time (this isn't
!     true if jpnl isn't 1 but we ignore this extra cost in this case).
!
!     IMPORTANT!
!     The stiff method (mf=22) in the SVODE routine is not intended to
!     work with families or steady state species. The only species types
!     allowed with this integrator are 'TR' (tracers) and 'CT' (constant
!     species). This should be checked for in CINIT at the start of the
!     run. All types are allowed with mf=10.
!
!     IMPORTANT!
!     The nag integrator only computes one point at a time. Since the
!     inner loop in prls.f is over the spatial domain this means that
!     the code will not vectorise if this integrator is used. Users 
!     should be aware that this will incur a very significant 
!     performance penalty.
!
!     Externals
!     ---------
!
!       svode  - Variable coefficient ODE solver, single precision
!                version from NETLIB. NOTE! This routine will use
!                routines from the BLAS and LINPACK NETLIB libraries.
!       vodycn - Interface routine to ASAD routine ycn which returns
!                tendencies for svode. Declared as external.
!
!     Local variables
!     ---------------
!
!     neq  - No. of equations for SVODE to solve.
!     atol - Absolute tolerance; errors below this concentration
!            are ignored.
!     rtol - Relative tolerance; error must be less than
!            y*rtol + atol (see svode.f for more details).
!
!----------------------------------------------------------------------
!

       USE messy_clamschem_asad_mod,       ONLY: f, pmintnd, cdt, ptol, y
       USE messy_clamschem_asad_mod_clams, ONLY: jpctr, jl_current
       !USE messy_clamschem_asad_jac,       ONLY: ASAD_JAC
       USE messy_clamschem_asad_dummy,     ONLY: asad_jac_dummy
       USE messy_clamschem_asad_svode,     ONLY: SVODE
       USE messy_clams_global,             ONLY: PREC

      !external vodycn

      integer :: n_points
 
      integer :: lrw, liw 
      real(prec), dimension(:), allocatable :: rwork, zf
      integer,    dimension(:), allocatable :: iwork
      real(prec), dimension(1) :: rpar
      integer,    dimension(1) :: ipar

      real(prec) :: starti, stopi, rtol, atol
      integer :: kmf, imf, neq, itol, itask, iopt, jl, &
                 istate, js, jdummy

      lrw = 22 + 9*jpctr + 2*jpctr*jpctr
      liw = 30 + jpctr
      allocate (rwork(lrw))
      allocate (iwork(liw))
      allocate (zf(jpctr))
      rwork = 0.0
      iwork = 0


!     -----------------------------------------------------------------
!          1.  Initialise local variables.
!              ---------- ----- ----------

      if ( kmf .ne. 10 .and. kmf .ne. 22 ) then
         write (6,*) '*** ASAD ERROR: kmf in vodedriv.f must be ',  &
         'either 10 or 22. It was set to ',kmf
         stop 'ASAD aborted'
      endif
      imf = kmf

      neq = jpctr
      itol = 1
      itask = 1
      iopt = 0

!     -----------------------------------------------------------------
!        2.  Integrate gridpoints individually.
!            --------- ---------- -------------

      do 210 jl = 1, n_points
      rtol = ptol
      atol = pmintnd(jl)
      jl_current=jl

      do 220 js = 1, neq
      zf(js) = f(jl,js)
 220  continue

      istate = 1
      starti = 0.0
      stopi = cdt

! ju_nt_20160617
!!$      call svode(vodycn,neq,zf,starti,stopi,itol,rtol,atol,itask,  &
!!$                 istate,iopt,rwork,lrw,iwork,liw,jdummy,imf,  &
!!$                 rpar,ipar)
!!$  ATTENTION: CHECK PARAMETERS !!!
!!$      call svode(vodycn,neq,zf,starti,stopi,itol,(/rtol/),(/atol/),itask,  &
!!$                 istate,iopt,rwork,lrw,iwork,liw,asad_jac,imf,  &
!!$                 rpar,ipar)
      call svode(vodycn,neq,zf,starti,stopi,itol,(/rtol/),(/atol/),itask,  &
                 istate,iopt,rwork,lrw,iwork,liw,asad_jac_dummy,imf,  &
                 rpar,ipar)
     if ( istate .ne. 2 ) then
         write (6,*) '** SVODE integrator failed with ISTATE=',istate
         stop 'ASAD ABORTED'
      endif

     do 250 js = 1, neq
      f(jl,js) = zf(js)
 250  continue

 210  continue

      deallocate (rwork,iwork,zf)

!     -----------------------------------------------------------------

      return
      end  subroutine vodedriv

      subroutine vodycn(neq,t,f,dfdt,rpar,ipar)

!     VODYCN  - Interface routine between SVODE and ASAD routine YCN.

!     The ASAD routine ycn organises the computation of the rates of
!     changes for the stiff integrators in ASAD. This routine simply
!     calls ycn with the correct arguments.

!     See the source code for svode.f for full details of the 
!     arguments passed to this routine.

      USE messy_clams_global,       ONLY: prec
      USE messy_clamschem_asad_ycn, ONLY: ASAD_YCN

!!!!! ju_nt_200160615
!      dimension f(neq), dfdt(neq)
      real(prec), dimension(neq) :: f, dfdt
      real(prec) :: t
!!!!!

      call asad_ycn(t,f,dfdt)
!      call ycn(t,f,dfdt)

      return
      end subroutine vodycn

   End Module messy_clamschem_asad_vodedriv
