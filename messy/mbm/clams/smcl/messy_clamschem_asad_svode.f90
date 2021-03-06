Module messy_clamschem_asad_svode

contains

!*********************************************************
!
!  ASAD: This is the SVODE integrator complete code. Taken
!        from the ODE (not ODEPACK) subroutine package.
!        See NETLIB (http://www.hensa.ac.uk/ftp/mirrors/netlib/master/)
!        for more details.
!
!        svode.f 4.1 01/15/97
!
!*********************************************************

      SUBROUTINE SVODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK, &
                  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
                  RPAR, IPAR)
      USE messy_clams_global,        ONLY: prec
      USE messy_clamschem_asad_blas, ONLY: scopy, sscal
      EXTERNAL F, JAC
      REAL(PREC) Y, T, TOUT, RTOL, ATOL, RWORK, RPAR
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, &
              MF, IPAR
      DIMENSION Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW), &
                RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! SVODE.. Variable-coefficient Ordinary Differential Equation solver,
! with fixed-leading coefficient implementation.
! This version is in single precision.
!
! SVODE solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! SVODE is a package based on the EPISODE and EPISODEB packages, and
! on the ODEPACK user interface standard, with minor modifications.
!-----------------------------------------------------------------------
! Revision History (YYMMDD)
!   890615  Date Written
!   890922  Added interrupt/restart ability, minor changes throughout.
!   910228  Minor revisions in line format,  prologue, etc.
!   920227  Modifications by D. Pang:
!           (1) Applied subgennam to get generic intrinsic names.
!           (2) Changed intrinsic names to generic in comments.
!           (3) Added *DECK lines before each routine.
!   920721  Names of routines and labeled Common blocks changed, so as
!           to be unique in combined single/double precision code (ACH).
!   920722  Minor revisions to prologue (ACH).
!-----------------------------------------------------------------------
! References..
!
! 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE: A Variable
!    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10 (1989),
!    pp. 1038-1051.  Also, LLNL Report UCRL-98412, June 1988.
! 2. G. D. Byrne and A. C. Hindmarsh, "A Polyalgorithm for the
!    Numerical Solution of Ordinary Differential Equations,"
!    ACM Trans. Math. Software, 1 (1975), pp. 71-96.
! 3. A. C. Hindmarsh and G. D. Byrne, "EPISODE: An Effective Package
!    for the Integration of Systems of Ordinary Differential
!    Equations," LLNL Report UCID-30112, Rev. 1, April 1977.
! 4. G. D. Byrne and A. C. Hindmarsh, "EPISODEB: An Experimental
!    Package for the Integration of Systems of Ordinary Differential
!    Equations with Banded Jacobians," LLNL Report UCID-30132, April
!    1976.
! 5. A. C. Hindmarsh, "ODEPACK, a Systematized Collection of ODE
!    Solvers," in Scientific Computing, R. S. Stepleman et al., eds.,
!    North-Holland, Amsterdam, 1983, pp. 55-64.
! 6. K. R. Jackson and R. Sacks-Davis, "An Alternative Implementation
!    of Variable Step-Size Multistep Formulas for Stiff ODEs," ACM
!    Trans. Math. Software, 6 (1980), pp. 295-318.
!-----------------------------------------------------------------------
! Authors..
!
!               Peter N. Brown and Alan C. Hindmarsh
!               Computing and Mathematics Research Division, L-316
!               Lawrence Livermore National Laboratory
!               Livermore, CA 94550
! and
!               George D. Byrne
!               Exxon Research and Engineering Co.
!               Clinton Township
!               Route 22 East
!               Annandale, NJ 08801
!-----------------------------------------------------------------------
! Summary of usage.
!
! Communication between the user and the SVODE package, for normal
! situations, is summarized here.  This summary describes only a subset
! of the full set of options available.  See the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  See also the example
! problem (with program and output) following this summary.
!
! A. First provide a subroutine of the form..
!
!           SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!           REAL T, Y, YDOT, RPAR
!           DIMENSION Y(NEQ), YDOT(NEQ)
!
! which supplies the vector function f by loading YDOT(i) with f(i).
!
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  If the problem is nonstiff,
! use a method flag MF = 10.  If it is stiff, there are four standard
! choices for MF (21, 22, 24, 25), and SVODE requires the Jacobian
! matrix in some form.  In these cases (MF .gt. 0), SVODE will use a
! saved copy of the Jacobian matrix.  If this is undesirable because of
! storage limitations, set MF to the corresponding negative value
! (-21, -22, -24, -25).  (See full description of MF below.)
! The Jacobian matrix is regarded either as full (MF = 21 or 22),
! or banded (MF = 24 or 25).  In the banded case, SVODE requires two
! half-bandwidth parameters ML and MU.  These are, respectively, the
! widths of the lower and upper parts of the band, excluding the main
! diagonal.  Thus the band consists of the locations (i,j) with
! i-ML .le. j .le. i+MU, and the full bandwidth is ML+MU+1.
!
! C. If the problem is stiff, you are encouraged to supply the Jacobian
! directly (MF = 21 or 24), but if this is not feasible, SVODE will
! compute it internally by difference quotients (MF = 22 or 25).
! If you are supplying the Jacobian, provide a subroutine of the form..
!
!           SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
!           REAL T, Y, PD, RPAR
!           DIMENSION Y(NEQ), PD(NROWPD,NEQ)
!
! which supplies df/dy by loading PD as follows..
!     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
! the partial derivative of f(i) with respect to y(j).  (Ignore the
! ML and MU arguments in this case.)
!     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
! df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
! PD from the top down.
!     In either case, only nonzero elements need be loaded.
!
! D. Write a main program which calls subroutine SVODE once for
! each point at which answers are desired.  This should also provide
! for possible use of logical unit 6 for output of error messages
! by SVODE.  On the first call to SVODE, supply arguments as follows..
! F      = Name of subroutine for right-hand side vector f.
!          This name must be declared external in calling program.
! NEQ    = Number of first order ODE-s.
! Y      = Array of initial values, of length NEQ.
! T      = The initial value of the independent variable.
! TOUT   = First point where output is desired (.ne. T).
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
! RTOL   = Relative tolerance parameter (scalar).
! ATOL   = Absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution.. Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = Integer flag (input and output).  Set ISTATE = 1.
! IOPT   = 0 to indicate no optional input used.
! RWORK  = Real work array of length at least..
!             20 + 16*NEQ                      for MF = 10,
!             22 +  9*NEQ + 2*NEQ**2           for MF = 21 or 22,
!             22 + 11*NEQ + (3*ML + 2*MU)*NEQ  for MF = 24 or 25.
! LRW    = Declared length of RWORK (in user's DIMENSION statement).
! IWORK  = Integer work array of length at least..
!             30        for MF = 10,
!             30 + NEQ  for MF = 21, 22, 24, or 25.
!          If MF = 24 or 25, input in IWORK(1),IWORK(2) the lower
!          and upper half-bandwidths ML,MU.
! LIW    = Declared length of IWORK (in user's DIMENSION).
! JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).
!          If used, this name must be declared external in calling
!          program.  If not used, pass a dummy name.
! MF     = Method flag.  Standard values are..
!          10 for nonstiff (Adams) method, no Jacobian used.
!          21 for stiff (BDF) method, user-supplied full Jacobian.
!          22 for stiff method, internally generated full Jacobian.
!          24 for stiff method, user-supplied banded Jacobian.
!          25 for stiff method, internally generated banded Jacobian.
! RPAR,IPAR = user-defined real and integer arrays passed to F and JAC.
! Note that the main program must declare arrays Y, RWORK, IWORK,
! and possibly ATOL, RPAR, and IPAR.
!
! E. The output from the first call (or any call) is..
!      Y = Array of computed values of y(t) vector.
!      T = Corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if SVODE was successful, negative otherwise.
!          -1 means excess work done on this call. (Perhaps wrong MF.)
!          -2 means excess accuracy requested. (Tolerances too small.)
!          -3 means illegal input detected. (See printed message.)
!          -4 means repeated error test failures. (Check all input.)
!          -5 means repeated convergence failures. (Perhaps bad
!             Jacobian supplied or wrong choice of MF or tolerances.)
!          -6 means error weight became zero during problem. (Solution
!             component i vanished, and ATOL or ATOL(i) = 0.)
!
! F. To continue the integration after a successful return, simply
! reset TOUT and call SVODE again.  No other parameters need be reset.
!
!-----------------------------------------------------------------------
! EXAMPLE PROBLEM
!
! The following is a simple example problem, with the coding
! needed for its solution by SVODE.  The problem is from chemical
! kinetics, and consists of the following three rate equations..
!     dy1/dt = -.04*y1 + 1.e4*y2*y3
!     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
!     dy3/dt = 3.e7*y2**2
! on the interval from t = 0.0 to t = 4.e10, with initial conditions
! y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
!
! The following coding solves this problem with SVODE, using MF = 21
! and printing results at t = .4, 4., ..., 4.e10.  It uses
! ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
! y2 has much smaller values.
! At the end of the run, statistical quantities of interest are
! printed. (See optional output in the full description below.)
! To generate Fortran source code, replace C in column 1 with a blank
! in the coding below.
!
!     EXTERNAL FEX, JEX
!     REAL ATOL, RPAR, RTOL, RWORK, T, TOUT, Y
!     DIMENSION Y(3), ATOL(3), RWORK(67), IWORK(33)
!     NEQ = 3
!     Y(1) = 1.0E0
!     Y(2) = 0.0E0
!     Y(3) = 0.0E0
!     T = 0.0E0
!     TOUT = 0.4E0
!     ITOL = 2
!     RTOL = 1.E-4
!     ATOL(1) = 1.E-8
!     ATOL(2) = 1.E-14
!     ATOL(3) = 1.E-6
!     ITASK = 1
!     ISTATE = 1
!     IOPT = 0
!     LRW = 67
!     LIW = 33
!     MF = 21
!     DO 40 IOUT = 1,12
!       CALL SVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
!    1            IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
!       WRITE(6,20)T,Y(1),Y(2),Y(3)
! 20    FORMAT(' At t =',E12.4,'   y =',3E14.6)
!       IF (ISTATE .LT. 0) GO TO 80
! 40    TOUT = TOUT*10.
!     WRITE(6,60) IWORK(11),IWORK(12),IWORK(13),IWORK(19),
!    1            IWORK(20),IWORK(21),IWORK(22)
! 60  FORMAT(/' No. steps =',I4,'   No. f-s =',I4,
!    1       '   No. J-s =',I4,'   No. LU-s =',I4/
!    2       '  No. nonlinear iterations =',I4/
!    3       '  No. nonlinear convergence failures =',I4/
!    4       '  No. error test failures =',I4/)
!     STOP
! 80  WRITE(6,90)ISTATE
! 90  FORMAT(///' Error halt.. ISTATE =',I3)
!     STOP
!     END
!
!     SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
!     REAL RPAR, T, Y, YDOT
!     DIMENSION Y(NEQ), YDOT(NEQ)
!     YDOT(1) = -.04E0*Y(1) + 1.E4*Y(2)*Y(3)
!     YDOT(3) = 3.E7*Y(2)*Y(2)
!     YDOT(2) = -YDOT(1) - YDOT(3)
!     RETURN
!     END
!
!     SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
!     REAL PD, RPAR, T, Y
!     DIMENSION Y(NEQ), PD(NRPD,NEQ)
!     PD(1,1) = -.04E0
!     PD(1,2) = 1.E4*Y(3)
!     PD(1,3) = 1.E4*Y(2)
!     PD(2,1) = .04E0
!     PD(2,3) = -PD(1,3)
!     PD(3,2) = 6.E7*Y(2)
!     PD(2,2) = -PD(1,2) - PD(3,2)
!     RETURN
!     END
!
! The following output was obtained from the above program on a
! Cray-1 computer with the CFT compiler.
!
! At t =  4.0000e-01   y =  9.851680e-01  3.386314e-05  1.479817e-02
! At t =  4.0000e+00   y =  9.055255e-01  2.240539e-05  9.445214e-02
! At t =  4.0000e+01   y =  7.158108e-01  9.184883e-06  2.841800e-01
! At t =  4.0000e+02   y =  4.505032e-01  3.222940e-06  5.494936e-01
! At t =  4.0000e+03   y =  1.832053e-01  8.942690e-07  8.167938e-01
! At t =  4.0000e+04   y =  3.898560e-02  1.621875e-07  9.610142e-01
! At t =  4.0000e+05   y =  4.935882e-03  1.984013e-08  9.950641e-01
! At t =  4.0000e+06   y =  5.166183e-04  2.067528e-09  9.994834e-01
! At t =  4.0000e+07   y =  5.201214e-05  2.080593e-10  9.999480e-01
! At t =  4.0000e+08   y =  5.213149e-06  2.085271e-11  9.999948e-01
! At t =  4.0000e+09   y =  5.183495e-07  2.073399e-12  9.999995e-01
! At t =  4.0000e+10   y =  5.450996e-08  2.180399e-13  9.999999e-01
!
! No. steps = 595   No. f-s = 832   No. J-s =  13   No. LU-s = 112
!  No. nonlinear iterations = 831
!  No. nonlinear convergence failures =   0
!  No. error test failures =  22
!-----------------------------------------------------------------------
! Full description of user interface to SVODE.
!
! The user interface to SVODE consists of the following parts.
!
! i.   The call sequence to subroutine SVODE, which is a driver
!      routine for the solver.  This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions is
!        * a description of optional input available through the
!          call sequence,
!        * a description of optional output (in the work arrays), and
!        * instructions for interrupting and restarting a solution.
!
! ii.  Descriptions of other routines in the SVODE package that may be
!      (optionally) called by the user.  These provide the ability to
!      alter error message handling, save and restore the internal
!      COMMON, and obtain specified derivatives of the solution y(t).
!
! iii. Descriptions of COMMON blocks to be declared in overlay
!      or similar environments.
!
! iv.  Description of two routines in the SVODE package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! Part i.  Call Sequence.
!
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RWORK and IWORK are also used for conditional and
! optional input and optional output.  (The term output here refers
! to the return from subroutine SVODE to the user's calling program.)
!
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 in the input.
!
! The descriptions of the call arguments are as follows.
!
! F      = The name of the user-supplied subroutine defining the
!          ODE system.  The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  Subroutine F is to
!          compute the function f.  It is to have the form
!               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!               REAL T, Y, YDOT, RPAR
!               DIMENSION Y(NEQ), YDOT(NEQ)
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output.  Y and YDOT are arrays of length NEQ.
!          (In the DIMENSION statement above, NEQ  can be replaced by
!          *  to make  Y  and  YDOT  assumed size arrays.)
!          Subroutine F should not alter Y(1),...,Y(NEQ).
!          F must be declared EXTERNAL in the calling program.
!
!          Subroutine F may access user-defined real and integer
!          work arrays RPAR and IPAR, which are to be dimensioned
!          in the main program.
!
!          If quantities computed in the F routine are needed
!          externally to SVODE, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use SVINDY instead.
!
! NEQ    = The size of the ODE system (number of first order
!          ordinary differential equations).  Used only for input.
!          NEQ may not be increased during the problem, but
!          can be decreased (with ISTATE = 3 in the input).
!
! Y      = A real array for the vector of dependent variables, of
!          length NEQ or more.  Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          On the first call, Y must contain the vector of initial
!          values.  In the output, Y contains the computed solution
!          evaluated at T.  If desired, the Y array may be used
!          for other purposes between calls to the solver.
!
!          This array is passed as the Y argument in all calls to
!          F and JAC.
!
! T      = The independent variable.  In the input, T is used only on
!          the first call, as the initial point of the integration.
!          In the output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
!
! TOUT   = The next value of t at which a computed solution is desired.
!          Used only for input.
!
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should .ne. T for the next call.
!          For the initial T, an input value of TOUT .ne. T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  Integration in either direction
!          (forward or backward in t) is permitted.
!
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT .ne. T).
!          Otherwise, TOUT is required on every call.
!
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal t interval, whose endpoints are
!          TCUR - HU and TCUR.  (See optional output, below, for
!          TCUR and HU.)
!
! ITOL   = An indicator for the type of error control.  See
!          description below under ATOL.  Used only for input.
!
! RTOL   = A relative error tolerance parameter, either a scalar or
!          an array of length NEQ.  See description below under ATOL.
!          Input only.
!
! ATOL   = An absolute error tolerance parameter, either a scalar or
!          an array of length NEQ.  Input only.
!
!          The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver.  The solver will
!          control the vector e = (e(i)) of estimated local errors
!          in Y, according to an inequality of the form
!                      rms-norm of ( e(i)/EWT(i) )   .le.   1,
!          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be non-negative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!
!             ITOL    RTOL       ATOL          EWT(i)
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation.  See Part iv below.
!
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL (i.e. of EWT) should be scaled
!          down uniformly.
!
! ITASK  = An index specifying the task to be performed.
!          Input only.  ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration.  This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RWORK(1).
!
!          Note..  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT (exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at T = TOUT are returned first).
!
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!
!          In the input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK.  Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, ML, MU,
!             and any of the optional input except H0.
!             (See IWORK description for ML and MU.)
!          Note..  A preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (Such a call is sometimes useful to include
!          the initial conditions in the output.)
!          Thus the first call for which TOUT .ne. T requires
!          ISTATE = 1 in the input.
!
!          In the output, ISTATE has the following values and meanings.
!           1  means nothing was done, as TOUT was equal to T with
!              ISTATE = 1 in the input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T.  (MXSTEP is an optional input
!              and is normally 500.)  To continue, the user may
!              simply reset ISTATE to a value .gt. 1 and call again.
!              (The excess work step counter will be reset to 0.)
!              In addition, the user may increase MXSTEP to avoid
!              this error return.  (See optional input below.)
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  This was detected before
!              completing the requested task, but the integration
!              was successful as far as T.  To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3.  The optional output TOLSF may be used for this
!              purpose.  (Note.. If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  See written message for details.
!              Note..  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration.  Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!
!          Note..  Since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other input, before
!          calling the solver again.
!
! IOPT   = An integer flag to specify whether or not any optional
!          input is being used on this call.  Input only.
!          The optional input is listed separately below.
!          IOPT = 0 means no optional input is being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means optional input is being used.
!
! RWORK  = A real working array (single precision).
!          The length of RWORK must be at least
!             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          LWM = length of work space for matrix-related data..
!          LWM = 0             if MITER = 0,
!          LWM = 2*NEQ**2 + 2  if MITER = 1 or 2, and MF.gt.0,
!          LWM = NEQ**2 + 2    if MITER = 1 or 2, and MF.lt.0,
!          LWM = NEQ + 2       if MITER = 3,
!          LWM = (3*ML+2*MU+2)*NEQ + 2 if MITER = 4 or 5, and MF.gt.0,
!          LWM = (2*ML+MU+1)*NEQ + 2   if MITER = 4 or 5, and MF.lt.0.
!          (See the MF description for METH and MITER.)
!          Thus if MAXORD has its default value and NEQ is constant,
!          this length is..
!             20 + 16*NEQ                    for MF = 10,
!             22 + 16*NEQ + 2*NEQ**2         for MF = 11 or 12,
!             22 + 16*NEQ + NEQ**2           for MF = -11 or -12,
!             22 + 17*NEQ                    for MF = 13,
!             22 + 18*NEQ + (3*ML+2*MU)*NEQ  for MF = 14 or 15,
!             22 + 17*NEQ + (2*ML+MU)*NEQ    for MF = -14 or -15,
!             20 +  9*NEQ                    for MF = 20,
!             22 +  9*NEQ + 2*NEQ**2         for MF = 21 or 22,
!             22 +  9*NEQ + NEQ**2           for MF = -21 or -22,
!             22 + 10*NEQ                    for MF = 23,
!             22 + 11*NEQ + (3*ML+2*MU)*NEQ  for MF = 24 or 25.
!             22 + 10*NEQ + (2*ML+MU)*NEQ    for MF = -24 or -25.
!          The first 20 words of RWORK are reserved for conditional
!          and optional output.
!
!          The following word in RWORK is a conditional input..
!            RWORK(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot.  Required if ITASK is
!                       4 or 5, and ignored otherwise.  (See ITASK.)
!
! LRW    = The length of the array RWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! IWORK  = An integer work array.  The length of IWORK must be at least
!             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
!             30 + NEQ  otherwise (abs(MF) = 11,12,14,15,21,22,24,25).
!          The first 30 words of IWORK are reserved for conditional and
!          optional input and optional output.
!
!          The following 2 words in IWORK are conditional input..
!            IWORK(1) = ML     These are the lower and upper
!            IWORK(2) = MU     half-bandwidths, respectively, of the
!                       banded Jacobian, excluding the main diagonal.
!                       The band is defined by the matrix locations
!                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
!                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
!                       These are required if MITER is 4 or 5, and
!                       ignored otherwise.  ML and MU may in fact be
!                       the band parameters for a matrix to which
!                       df/dy is only approximately equal.
!
! LIW    = the length of the array IWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! Note..  The work arrays must not be altered between calls to SVODE
! for the same problem, except possibly for the conditional and
! optional input, and except for the last 3*NEQ words of RWORK.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside SVODE between calls, if
! desired (but not for use by F or JAC).
!
! JAC    = The name of the user-supplied routine (MITER = 1 or 4) to
!          compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  It is to have the form
!               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD,
!                               RPAR, IPAR)
!               REAL T, Y, PD, RPAR
!               DIMENSION Y(NEQ), PD(NROWPD, NEQ)
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!          PD is to be loaded with partial derivatives (elements of the
!          Jacobian matrix) in the output.  PD must be given a first
!          dimension of NROWPD.  T and Y have the same meaning as in
!          Subroutine F.  (In the DIMENSION statement above, NEQ can
!          be replaced by  *  to make Y and PD assumed size arrays.)
!               In the full matrix case (MITER = 1), ML and MU are
!          ignored, and the Jacobian is to be loaded into PD in
!          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!               In the band matrix case (MITER = 4), the elements
!          within the band are to be loaded into PD in columnwise
!          manner, with diagonal lines of df/dy loaded into the rows
!          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters. (See IWORK).
!          The locations in PD in the two triangular areas which
!          correspond to nonexistent matrix elements can be ignored
!          or loaded arbitrarily, as they are overwritten by SVODE.
!               JAC need not provide df/dy exactly.  A crude
!          approximation (possibly with a smaller bandwidth) will do.
!               In either case, PD is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Each call to JAC is preceded by a call to F with the same
!          arguments NEQ, T, and Y.  Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user COMMON block by F and not recomputed by JAC,
!          if desired.  Also, JAC may alter the Y array, if desired.
!          JAC must be declared external in the calling program.
!               Subroutine JAC may access user-defined real and integer
!          work arrays, RPAR and IPAR, whose dimensions are set by the
!          user in the main program.
!
! MF     = The method flag.  Used only for input.  The legal values of
!          MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25,
!          -11, -12, -14, -15, -21, -22, -24, -25.
!          MF is a signed two-digit integer, MF = JSV*(10*METH + MITER).
!          JSV = SIGN(MF) indicates the Jacobian-saving strategy..
!            JSV =  1 means a copy of the Jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            JSV = -1 means a copy of the Jacobian is not saved
!                     (valid only for MITER = 1, 2, 4, or 5).
!          METH indicates the basic linear multistep method..
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method..
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      full (NEQ by NEQ) Jacobian.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) full Jacobian
!                      (using NEQ extra calls to F per df/dy value).
!            MITER = 3 means chord iteration with an internally
!                      generated diagonal Jacobian approximation
!                      (using 1 extra call to F per df/dy evaluation).
!            MITER = 4 means chord iteration with a user-supplied
!                      banded Jacobian.
!            MITER = 5 means chord iteration with an internally
!                      generated banded Jacobian (using ML+MU+1 extra
!                      calls to F per df/dy evaluation).
!          If MITER = 1 or 4, the user must supply a subroutine JAC
!          (the name is arbitrary) as described above under JAC.
!          For other values of MITER, a dummy argument can be used.
!
! RPAR     User-specified array used to communicate real parameters
!          to user-supplied subroutines.  If RPAR is a vector, then
!          it must be dimensioned in the user's main program.  If it
!          is unused or it is a scalar, then it need not be
!          dimensioned.
!
! IPAR     User-specified array used to communicate integer parameter
!          to user-supplied subroutines.  The comments on dimensioning
!          RPAR apply to IPAR.
!-----------------------------------------------------------------------
! Optional Input.
!
! The following is a list of the optional input provided for in the
! call sequence.  (See also Part ii.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of this input requires IOPT = 1, and in that
! case all of this input is examined.  A value of zero for any
! of these optional input variables will cause the default value to be
! used.  Thus to use a subset of the optional input, simply preload
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! NAME    LOCATION      MEANING AND DEFAULT VALUE
!
! H0      RWORK(5)  The step size to be attempted on the first step.
!                   The default value is determined by the solver.
!
! HMAX    RWORK(6)  The maximum absolute step size allowed.
!                   The default value is infinite.
!
! HMIN    RWORK(7)  The minimum absolute step size allowed.
!                   The default value is 0.  (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
!
! MAXORD  IWORK(5)  The maximum order to be allowed.  The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
!
! MXSTEP  IWORK(6)  Maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 500.
!
! MXHNIL  IWORK(7)  Maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value.  The default value is 10.
!
!-----------------------------------------------------------------------
! Optional Output.
!
! As optional additional output from SVODE, the variables listed
! below are quantities related to the performance of SVODE
! which are available to the user.  These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of this output is defined
! on any successful return from SVODE, and on any return with
! ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, output relevant to the error will be defined,
! as noted below.
!
! NAME    LOCATION      MEANING
!
! HU      RWORK(11) The step size in t last used (successfully).
!
! HCUR    RWORK(12) The step size to be attempted on the next step.
!
! TCUR    RWORK(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
!
! TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise).  If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
!
! NST     IWORK(11) The number of steps taken for the problem so far.
!
! NFE     IWORK(12) The number of f evaluations for the problem so far.
!
! NJE     IWORK(13) The number of Jacobian evaluations so far.
!
! NQU     IWORK(14) The method order last used (successfully).
!
! NQCUR   IWORK(15) The order to be attempted on the next step.
!
! IMXER   IWORK(16) The index of the component of largest magnitude in
!                   the weighted local error vector ( e(i)/EWT(i) ),
!                   on an error return with ISTATE = -4 or -5.
!
! LENRW   IWORK(17) The length of RWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENIW   IWORK(18) The length of IWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! NLU     IWORK(19) The number of matrix LU decompositions so far.
!
! NNI     IWORK(20) The number of nonlinear (Newton) iterations so far.
!
! NCFN    IWORK(21) The number of convergence failures of the nonlinear
!                   solver so far.
!
! NETF    IWORK(22) The number of error test failures of the integrator
!                   so far.
!
! The following two arrays are segments of the RWORK array which
! may also be of interest to the user as optional output.
! For each array, the table below gives its internal name,
! its base address in RWORK, and its description.
!
! NAME    BASE ADDRESS      DESCRIPTION
!
! YH      21             The Nordsieck history array, of size NYH by
!                        (NQCUR + 1), where NYH is the initial value
!                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
!                        of YH contains HCUR**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the
!                        solution, evaluated at t = TCUR.
!
! ACOR     LENRW-NEQ+1   Array of size NEQ used for the accumulated
!                        corrections on each step, scaled in the output
!                        to represent the estimated local error in Y
!                        on the last step.  This is the vector e in
!                        the description of the error control.  It is
!                        defined only on a successful return from SVODE.
!
!-----------------------------------------------------------------------
! Interrupting and Restarting
!
! If the integration of a given problem by SVODE is to be
! interrrupted and then later continued, such as when restarting
! an interrupted run or alternating between two or more ODE problems,
! the user should save, following the return from the last SVODE call
! prior to the interruption, the contents of the call sequence
! variables and internal COMMON blocks, and later restore these
! values before the next SVODE call for that problem.  To save
! and restore the COMMON blocks, use subroutine SVSRCO, as
! described below in part ii.
!
! In addition, if non-default values for either LUN or MFLAG are
! desired, an extra call to XSETUN and/or XSETF should be made just
! before continuing the integration.  See Part ii below for details.
!
!-----------------------------------------------------------------------
! Part ii.  Other Routines Callable.
!
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with SVODE.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!
!     FORM OF CALL                  FUNCTION
!  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
!                             output of messages from SVODE, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!
!  CALL XSETF(MFLAG)          Set a flag to control the printing of
!                             messages by SVODE.
!                             MFLAG = 0 means do not print. (Danger..
!                             This risks losing valuable information.)
!                             MFLAG = 1 means print (the default).
!
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!
!  CALL SVSRCO(RSAV,ISAV,JOB) Saves and restores the contents of
!                             the internal COMMON blocks used by
!                             SVODE. (See Part iii below.)
!                             RSAV must be a real array of length 49
!                             or more, and ISAV must be an integer
!                             array of length 40 or more.
!                             JOB=1 means save COMMON into RSAV/ISAV.
!                             JOB=2 means restore COMMON from RSAV/ISAV.
!                                SVSRCO is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with SVODE.
!
!  CALL SVINDY(,,,,,)         Provide derivatives of y, of various
!        (See below.)         orders, at a specified point T, if
!                             desired.  It may be called only after
!                             a successful return from SVODE.
!
! The detailed instructions for using SVINDY are as follows.
! The form of the call is..
!
!  CALL SVINDY (T, K, RWORK(21), NYH, DKY, IFLAG)
!
! The input parameters are..
!
! T         = Value of independent variable where answers are desired
!             (normally the same as the T last returned by SVODE).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional output for TCUR and HU.)
! K         = Integer order of the derivative desired.  K must satisfy
!             0 .le. K .le. NQCUR, where NQCUR is the current order
!             (see optional output).  The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by SVODE directly.  Since NQCUR .ge. 1, the first
!             derivative dy/dt is always available with SVINDY.
! RWORK(21) = The base address of the history array YH.
! NYH       = Column length of YH, equal to the initial value of NEQ.
!
! The output parameters are..
!
! DKY       = A real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = Integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
!-----------------------------------------------------------------------
! Part iii.  COMMON Blocks.
! If SVODE is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in..
!   (1) the call sequence to SVODE,
!   (2) the two internal COMMON blocks
!         /SVOD01/  of length  81  (48 single precision words
!                         followed by 33 integer words),
!         /SVOD02/  of length  9  (1 single precision word
!                         followed by 8 integer words),
!
! If SVODE is used on a system in which the contents of internal
! COMMON blocks are not preserved between calls, the user should
! declare the above two COMMON blocks in his main program to insure
! that their contents are preserved.
!
!-----------------------------------------------------------------------
! Part iv.  Optionally Replaceable Solver Routines.
!
! Below are descriptions of two routines in the SVODE package which
! relate to the measurement of errors.  Either routine can be
! replaced by a user-supplied version, if desired.  However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note.. The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) SEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above..
!     SUBROUTINE SEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the SVODE call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by SEWSET.
!
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparison with
! errors in Y(i).  The EWT array returned by SEWSET is passed to the
! SVNORM routine (See below.), and also used by SVODE in the computation
! of the optional output IMXER, the diagonal Jacobian approximation,
! and the increments for difference quotient Jacobians.
!
! In the user-supplied version of SEWSET, it may be desirable to use
! the current values of derivatives of y.  Derivatives up to order NQ
! are available from the history array YH, described above under
! Optional Output.  In SEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of h**j/factorial(j).  On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
! can be obtained by including in SEWSET the statements..
!     REAL RVOD, H, HU
!     COMMON /SVOD01/ RVOD(48), IVOD(33)
!     COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!     NQ = IVOD(28)
!     H = RVOD(21)
! Thus, for example, the current value of dy/dt can be obtained as
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
! unnecessary when NST = 0).
!
! (b) SVNORM.
! The following is a real function routine which computes the weighted
! root-mean-square norm of a vector v..
!     D = SVNORM (N, V, W)
! where..
!   N = the length of the vector,
!   V = real array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = sqrt( (1/N) * sum(V(i)*W(i))**2 ).
! SVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by subroutine SEWSET.
!
! If the user supplies this function, it should return a non-negative
! value of SVNORM suitable for use in the error control in SVODE.
! None of the arguments should be altered by SVNORM.
! For example, a user-supplied SVNORM routine might..
!   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of Y.
!-----------------------------------------------------------------------
! Other Routines in the SVODE Package.
!
! In addition to subroutine SVODE, the SVODE package includes the
! following subroutines and function routines..
!  SVHIN     computes an approximate step size for the initial step.
!  SVINDY    computes an interpolated value of the y vector at t = TOUT.
!  SVSTEP    is the core integrator, which does one step of the
!            integration and the associated error control.
!  SVSET     sets all method coefficients and test constants.
!  SVNLSD    solves the underlying nonlinear system -- the corrector.
!  SVJAC     computes and preprocesses the Jacobian matrix J = df/dy
!            and the Newton iteration matrix P = I - (h/l1)*J.
!  SVSOL     manages solution of linear system in chord iteration.
!  SVJUST    adjusts the history array on a change of order.
!  SEWSET    sets the error weight vector EWT before each step.
!  SVNORM    computes the weighted r.m.s. norm of a vector.
!  SVSRCO    is a user-callable routines to save and restore
!            the contents of the internal COMMON blocks.
!  SACOPY    is a routine to copy one two-dimensional array to another.
!  SGEFA and SGESL   are routines from LINPACK for solving full
!            systems of linear algebraic equations.
!  SGBFA and SGBSL   are routines from LINPACK for solving banded
!            linear systems.
!  SAXPY, SSCAL, and SCOPY are basic linear algebra modules (BLAS).
!  R1MACH    sets the unit roundoff of the machine.
!  XERRWV, XSETUN, XSETF, LUNSAV, and MFLGSV handle the printing of all
!            error messages and warnings.  XERRWV is machine-dependent.
! Note..  SVNORM, R1MACH, LUNSAV, and MFLGSV are function routines.
! All the others are subroutines.
!
! The intrinsic and external routines used by the SVODE package are..
! ABS, MAX, MIN, REAL, SIGN, SQRT, and WRITE.
!
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block SVOD01 --------------------
!
      REAL(PREC) ACNRM, CCMXJ, CONP, CRATE, DRC, EL,  &
           ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
           RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
              L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
              LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
              N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
              NSLP, NYH
!
! Type declarations for labeled COMMON block SVOD02 --------------------
!
      REAL(PREC) HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
! ju_nt_20160617
!!$      EXTERNAL SVNLSD
      LOGICAL IHIT
      REAL(PREC) ATOLI, BIG, EWTI, FOUR, H0, HMAX, HMX, HUN, ONE,   &
         PT2, RH, RTOLI, SIZE, TCRIT, TNEXT, TOLSF, TP, TWO, ZERO
      INTEGER I, IER, IFLAG, IMXER, JCO, KGO, LENIW, LENJ, LENP, LENRW,  &
         LENWM, LF0, MBAND, ML, MORD, MU, MXHNL0, MXSTP0, NITER, NSLAST
! op_pj_20160830+
!!$   CHARACTER*80 MSG
      CHARACTER(LEN=80) MSG
! op_pj_20160830-
!
! Type declaration for function subroutines called ---------------------
!
! ju_nt_20160617
!!$      REAL(PREC) R1MACH, SVNORM
!
      DIMENSION MORD(2)
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to SVODE.
!-----------------------------------------------------------------------
      SAVE MORD, MXHNL0, MXSTP0
      SAVE ZERO, ONE, TWO, FOUR, PT2, HUN
!-----------------------------------------------------------------------
! The following internal COMMON blocks contain variables which are
! communicated between subroutines in the SVODE package, or which are
! to be saved between calls to SVODE.
! In each block, real variables precede integers.
! The block /SVOD01/ appears in subroutines SVODE, SVINDY, SVSTEP,
! SVSET, SVNLSD, SVJAC, SVSOL, SVJUST and SVSRCO.
! The block /SVOD02/ appears in subroutines SVODE, SVINDY, SVSTEP,
! SVNLSD, SVJAC, and SVSRCO.
!
! The variables stored in the internal COMMON blocks are as follows..
!
! ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
! CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
! CONP   = The saved value of TQ(5).
! CRATE  = Estimated corrector convergence rate constant.
! DRC    = Relative change in H*RL1 since last SVJAC call.
! EL     = Real array of integration coefficients.  See SVSET.
! ETA    = Saved tentative ratio of new to old H.
! ETAMAX = Saved maximum value of ETA to be allowed.
! H      = The step size.
! HMIN   = The minimum absolute value of the step size H to be used.
! HMXI   = Inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
! HNEW   = The step size to be attempted on the next step.
! HSCAL  = Stepsize in scaling of YH array.
! PRL1   = The saved value of RL1.
! RC     = Ratio of current H*RL1 to value on last SVJAC call.
! RL1    = The reciprocal of the coefficient EL(1).
! TAU    = Real vector of past NQ step sizes, length 13.
! TQ     = A real vector of length 5 in which SVSET stores constants
!          used for the convergence test, the error test, and the
!          selection of H at a new order.
! TN     = The independent variable, updated on each step taken.
! UROUND = The machine unit roundoff.  The smallest positive real number
!          such that  1.0 + UROUND .ne. 1.0
! ICF    = Integer flag for convergence failure in SVNLSD..
!            0 means no failures.
!            1 means convergence failure with out of date Jacobian
!                   (recoverable error).
!            2 means convergence failure with current Jacobian or
!                   singular matrix (unrecoverable error).
! INIT   = Saved integer flag indicating whether initialization of the
!          problem has been done (INIT = 1) or not.
! IPUP   = Saved flag to signal updating of Newton matrix.
! JCUR   = Output flag from SVJAC showing Jacobian status..
!            JCUR = 0 means J is not current.
!            JCUR = 1 means J is current.
! JSTART = Integer flag used as input to SVSTEP..
!            0  means perform the first step.
!            1  means take a new step continuing from the last.
!            -1 means take the next step with a new value of MAXORD,
!                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
!          On return, SVSTEP sets JSTART = 1.
! JSV    = Integer flag for Jacobian saving, = sign(MF).
! KFLAG  = A completion code from SVSTEP with the following meanings..
!               0      the step was succesful.
!              -1      the requested error could not be achieved.
!              -2      corrector convergence could not be achieved.
!              -3, -4  fatal error in VNLS (can not occur here).
! KUTH   = Input flag to SVSTEP showing whether H was reduced by the
!          driver.  KUTH = 1 if H was reduced, = 0 otherwise.
! L      = Integer variable, NQ + 1, current order plus one.
! LMAX   = MAXORD + 1 (used for dimensioning).
! LOCJS  = A pointer to the saved Jacobian, whose storage starts at
!          WM(LOCJS), if JSV = 1.
! LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers
!          to segments of RWORK and IWORK.
! MAXORD = The maximum order of integration method to be allowed.
! METH/MITER = The method flags.  See MF.
! MSBJ   = The maximum number of steps between J evaluations, = 50.
! MXHNIL = Saved value of optional input MXHNIL.
! MXSTEP = Saved value of optional input MXSTEP.
! N      = The number of first-order ODEs, = NEQ.
! NEWH   = Saved integer to flag change of H.
! NEWQ   = The method order to be used on the next step.
! NHNIL  = Saved counter for occurrences of T + H = T.
! NQ     = Integer variable, the current integration method order.
! NQNYH  = Saved value of NQ*NYH.
! NQWAIT = A counter controlling the frequency of order changes.
!          An order change is about to be considered if NQWAIT = 1.
! NSLJ   = The number of steps taken as of the last Jacobian update.
! NSLP   = Saved value of NST as of last Newton matrix update.
! NYH    = Saved value of the initial value of NEQ.
! HU     = The step size in t last used.
! NCFN   = Number of nonlinear convergence failures so far.
! NETF   = The number of error test failures of the integrator so far.
! NFE    = The number of f evaluations for the problem so far.
! NJE    = The number of Jacobian evaluations so far.
! NLU    = The number of matrix LU decompositions so far.
! NNI    = Number of nonlinear iterations so far.
! NQU    = The method order last used.
! NST    = The number of steps taken for the problem so far.
!-----------------------------------------------------------------------
      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),  &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,  &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
                      NSLP, NYH
      COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

      DATA  MORD(1) /12/, MORD(2) /5/, MXSTP0 /500/, MXHNL0 /10/
       DATA ZERO /0.0E0/, ONE /1.0E0/, TWO /2.0E0/, FOUR /4.0E0/,  &
          PT2 /0.2E0/, HUN /100.0E0/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .gt. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
!-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .NE. 1) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all input and various initializations.
!
! First check legality of the non-optional input NEQ, ITOL, IOPT,
! MF, ML, and MU.
!-----------------------------------------------------------------------
 20   IF (NEQ .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ .GT. N) GO TO 605
 25   N = NEQ
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      JSV = SIGN(1,MF)
      MF = ABS(MF)
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0 .OR. MITER .GT. 5) GO TO 608
      IF (MITER .LE. 3) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
! Next process and check the optional input. ---------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = ZERO
      HMXI = ZERO
      HMIN = ZERO
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. ZERO) GO TO 615
      HMXI = ZERO
      IF (HMAX .GT. ZERO) HMXI = ONE/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. ZERO) GO TO 616
!-----------------------------------------------------------------------
! Set work array pointers and check lengths LRW and LIW.
! Pointers to segments of RWORK and IWORK are named by prefixing L to
! the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
! Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
! Within WM, LOCJS is the location of the saved Jacobian (JSV .gt. 0).
!-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      JCO = MAX(0,JSV)
      IF (MITER .EQ. 0) LENWM = 0
      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
        LENWM = 2 + (1 + JCO)*N*N
        LOCJS = N*N + 3
      ENDIF
      IF (MITER .EQ. 3) LENWM = 2 + N
      IF (MITER .EQ. 4 .OR. MITER .EQ. 5) THEN
        MBAND = ML + MU + 1
        LENP = (MBAND + ML)*N
        LENJ = MBAND*N
        LENWM = 2 + LENP + JCO*LENJ
        LOCJS = LENP + 3
        ENDIF
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LACOR = LSAVF + N
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 30 + N
      IF (MITER .EQ. 0 .OR. MITER .EQ. 3) LENIW = 30
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
! Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. ZERO) GO TO 619
        IF (ATOLI .LT. ZERO) GO TO 620
 70     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
! If ISTATE = 3, set flag to signal parameter changes to SVSTEP. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 90
! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      CALL SCOPY (N, RWORK(LWM), 1, RWORK(LSAVF), 1)
! Reload WM(1) = RWORK(LWM), since LWM may have changed. ---------------
 90   IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
!-----------------------------------------------------------------------
! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! and the calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.
!-----------------------------------------------------------------------
 100  UROUND = R1MACH(4)
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
      IF (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO)  &
         H0 = TCRIT - T
 110  JSTART = 0
      IF (MITER .GT. 0) RWORK(LWM) = SQRT(UROUND)
      CCMXJ = PT2
      MSBJ = 50
      NHNIL = 0
      NST = 0
      NJE = 0
      NNI = 0
      NCFN = 0
      NETF = 0
      NLU = 0
      NSLJ = 0
      NSLAST = 0
      HU = ZERO
      NQU = 0
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (N, T, Y, RWORK(LF0), RPAR, IPAR)
      NFE = 1
! Load the initial value vector in YH. ---------------------------------
      CALL SCOPY (N, Y, 1, RWORK(LYH), 1)
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = ONE
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 621
 120    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
      IF (H0 .NE. ZERO) GO TO 180
! Call SVHIN to set initial step size H0 to be attempted. --------------
      CALL SVHIN (N, T, RWORK(LYH), RWORK(LF0), F, RPAR, IPAR, TOUT,  &
         UROUND, RWORK(LEWT), ITOL, ATOL, Y, RWORK(LACOR), H0,  &
         NITER, IER)
      NFE = NFE + NITER
      IF (IER .NE. 0) GO TO 622
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. ONE) H0 = H0/RH
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      CALL SSCAL (N, H0, RWORK(LF0), 1)
      GO TO 270
!-----------------------------------------------------------------------
! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
 200  NSLAST = NST
      KUTH = 0
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL SVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(ONE + HUN*UROUND)
      IF ((TP - TOUT)*H .GT. ZERO) GO TO 623
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. ZERO) GO TO 625
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 245
      CALL SVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
!-----------------------------------------------------------------------
! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator SVSTEP.
!
! This is a looping point for the integration steps.
!
! First check for too many steps being taken, update EWT (if not at
! start of problem), check for too much accuracy being requested, and
! check for H below the roundoff level in T.
!-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL SEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 510
 260    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*SVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. ONE) GO TO 280
      TOLSF = TOLSF*TWO
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'SVODE--  Warning..internal T (=R1) and H (=R2) are'
      CALL XERRWV (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWV (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      (H = step size). solver will continue anyway'
      CALL XERRWV (MSG, 50, 101, 1, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'SVODE--  Above warning has been issued I1 times.  '
      CALL XERRWV (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      it will not be issued again for this problem'
      CALL XERRWV (MSG, 50, 102, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
 290  CONTINUE
!-----------------------------------------------------------------------
! CALL SVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR,
!              WM, IWM, F, JAC, F, SVNLSD, RPAR, IPAR)
!-----------------------------------------------------------------------
      CALL SVSTEP (Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),  &
         RWORK(LSAVF), Y, RWORK(LACOR), RWORK(LWM), IWORK(LIWM),  &
         F, JAC, F, SVNLSD, RPAR, IPAR)
      KGO = 1 - KFLAG
! Branch on KFLAG.  Note..In this version, KFLAG can not be set to -3.
!  KFLAG .eq. 0,   -1,  -2
      GO TO (300, 530, 540), KGO
!-----------------------------------------------------------------------
! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0).  Test for stop conditions.
!-----------------------------------------------------------------------
 300  INIT = 1
      KUTH = 0
      GO TO (310, 400, 330, 340, 350), ITASK
! ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL SVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. ZERO) GO TO 400
      GO TO 250
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. ZERO) GO TO 345
      CALL SVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
      GO TO 250
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
!-----------------------------------------------------------------------
! Block G.
! The following block handles all successful returns from SVODE.
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional output is loaded into the work
! arrays before returning.
!-----------------------------------------------------------------------
 400  CONTINUE
      CALL SCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = HNEW
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NEWQ
      IWORK(19) = NLU
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      RETURN
!-----------------------------------------------------------------------
! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input.  First the error message routine is called.
! if there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH, T is set to TN, and the illegal input
! The optional output is loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'SVODE--  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWV (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWV (MSG, 50, 201, 1, 1, MXSTEP, 0, 1, TN, ZERO)
      ISTATE = -1
      GO TO 580
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'SVODE--  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWV (MSG, 50, 202, 1, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
! Too much accuracy requested for machine precision. -------------------
 520  MSG = 'SVODE--  At T (=R1), too much accuracy requested  '
      CALL XERRWV (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      for precision of machine..  see TOLSF (=R2) '
      CALL XERRWV (MSG, 50, 203, 1, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'SVODE--  At T(=R1) and step size H(=R2), the error'
      CALL XERRWV (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      test failed repeatedly or with abs(H) = HMIN'
      CALL XERRWV (MSG, 50, 204, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
! KFLAG = -2.  Convergence failed repeatedly or with abs(H) = HMIN. ----
 540  MSG = 'SVODE--  At T (=R1) and step size H (=R2), the    '
      CALL XERRWV (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWV (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      or with abs(H) = HMIN   '
      CALL XERRWV (MSG, 30, 205, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -5
! Compute IMXER if relevant. -------------------------------------------
 560  BIG = ZERO
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
! Set Y vector, T, and optional output. --------------------------------
 580  CONTINUE
      CALL SCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = NLU
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      RETURN
!-----------------------------------------------------------------------
! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! First the error message routine is called.   If the illegal input
! is a negative ISTATE, the run is aborted (apparent infinite loop).
!-----------------------------------------------------------------------
 601  MSG = 'SVODE--  ISTATE (=I1) illegal '
      CALL XERRWV (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'SVODE--  ITASK (=I1) illegal  '
      CALL XERRWV (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
      GO TO 700
 603  MSG='SVODE--  ISTATE (=I1) .gt. 1 but SVODE not initialized      '
      CALL XERRWV (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      GO TO 700
 604  MSG = 'SVODE--  NEQ (=I1) .lt. 1     '
      CALL XERRWV (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
      GO TO 700
 605  MSG = 'SVODE--  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWV (MSG, 50, 5, 1, 2, N, NEQ, 0, ZERO, ZERO)
      GO TO 700
 606  MSG = 'SVODE--  ITOL (=I1) illegal   '
      CALL XERRWV (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
      GO TO 700
 607  MSG = 'SVODE--  IOPT (=I1) illegal   '
      CALL XERRWV (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
      GO TO 700
 608  MSG = 'SVODE--  MF (=I1) illegal     '
      CALL XERRWV (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
      GO TO 700
 609  MSG = 'SVODE--  ML (=I1) illegal.. .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWV (MSG, 50, 9, 1, 2, ML, NEQ, 0, ZERO, ZERO)
      GO TO 700
 610  MSG = 'SVODE--  MU (=I1) illegal.. .lt.0 or .ge.NEQ (=I2)'
      CALL XERRWV (MSG, 50, 10, 1, 2, MU, NEQ, 0, ZERO, ZERO)
      GO TO 700
 611  MSG = 'SVODE--  MAXORD (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 11, 1, 1, MAXORD, 0, 0, ZERO, ZERO)
      GO TO 700
 612  MSG = 'SVODE--  MXSTEP (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 12, 1, 1, MXSTEP, 0, 0, ZERO, ZERO)
      GO TO 700
 613  MSG = 'SVODE--  MXHNIL (=I1) .lt. 0  '
      CALL XERRWV (MSG, 30, 13, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
      GO TO 700
 614  MSG = 'SVODE--  TOUT (=R1) behind T (=R2)      '
      CALL XERRWV (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
      MSG = '      integration direction is given by H0 (=R1)  '
      CALL XERRWV (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
      GO TO 700
 615  MSG = 'SVODE--  HMAX (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
      GO TO 700
 616  MSG = 'SVODE--  HMIN (=R1) .lt. 0.0  '
      CALL XERRWV (MSG, 30, 16, 1, 0, 0, 0, 1, HMIN, ZERO)
      GO TO 700
 617  CONTINUE
      MSG='SVODE--  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWV (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
      GO TO 700
 618  CONTINUE
      MSG='SVODE--  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWV (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
      GO TO 700
 619  MSG = 'SVODE--  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
      GO TO 700
 620  MSG = 'SVODE--  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWV (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'SVODE--  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWV (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
      GO TO 700
 622  CONTINUE
      MSG='SVODE--  TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWV (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='SVODE--  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWV (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='SVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWV (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='SVODE--  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWV (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'SVODE--  At start of problem, too much accuracy   '
      CALL XERRWV (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      requested for precision of machine..  see TOLSF (=R1) '
      CALL XERRWV (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG='SVODE--  Trouble from SVINDY.  ITASK = I1, TOUT = R1.       '
      CALL XERRWV (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)

 700  CONTINUE
      ISTATE = -3
      RETURN

 800  MSG = 'SVODE--  Run aborted.. apparent infinite loop     '
      CALL XERRWV (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
      RETURN
!----------------------- End of Subroutine SVODE -----------------------
      END SUBROUTINE SVODE

      SUBROUTINE SVHIN (N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,  &
         EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)
      USE messy_clams_global, ONLY: prec
      EXTERNAL F
      REAL(PREC) T0, Y0, YDOT, RPAR, TOUT, UROUND, EWT, ATOL, Y,  &
         TEMP, H0
      INTEGER N, IPAR, ITOL, NITER, IER
      DIMENSION Y0(*), YDOT(*), EWT(*), ATOL(*), Y(*),  &
         TEMP(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
!                        EWT, ITOL, ATOL, Y, TEMP
! Call sequence output -- H0, NITER, IER
! COMMON block variables accessed -- None
!
! Subroutines called by SVHIN.. F
! Function routines called by SVHIN.. SVNORM
!-----------------------------------------------------------------------
! This routine computes the step size, H0, to be attempted on the
! first step, when the user has not supplied a value for this.
!
! First we check that TOUT - T0 differs significantly from zero.  Then
! an iteration is done to approximate the initial second derivative
! and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
! A bias factor of 1/2 is applied to the resulting h.
! The sign of H0 is inferred from the initial values of TOUT and T0.
!
! Communication with SVHIN is done with the following variables..
!
! N      = Size of ODE system, input.
! T0     = Initial value of independent variable, input.
! Y0     = Vector of initial conditions, input.
! YDOT   = Vector of initial first derivatives, input.
! F      = Name of subroutine for right-hand side f(t,y), input.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
! TOUT   = First output value of independent variable
! UROUND = Machine unit roundoff
! EWT, ITOL, ATOL = Error weights and tolerance parameters
!                   as described in the driver routine, input.
! Y, TEMP = Work arrays of length N.
! H0     = Step size to be attempted, output.
! NITER  = Number of iterations (and of f evaluations) to compute H0,
!          output.
! IER    = The error flag, returned with the value
!          IER = 0  if no trouble occurred, or
!          IER = -1 if TOUT and T0 are considered too close to proceed.
!-----------------------------------------------------------------------

! Type declarations for local variables --------------------------------

      REAL(PREC) AFI, ATOLI, DELYI, HALF, HG, HLB, HNEW, HRAT,  &
           HUB, HUN, PT1, T1, TDIST, TROUND, TWO, YDDNRM
      INTEGER I, ITER

! Type declaration for function subroutines called ---------------------

! ju_nt_20160617
!!$      REAL(PREC) SVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE HALF, HUN, PT1, TWO
      DATA HALF /0.5E0/, HUN /100.0E0/, PT1 /0.1E0/, TWO /2.0E0/

      NITER = 0
      TDIST = ABS(TOUT - T0)
      TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
      IF (TDIST .LT. TWO*TROUND) GO TO 100

! Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
      HLB = HUN*TROUND
! Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
      HUB = PT1*TDIST
      ATOLI = ATOL(1)
      DO 10 I = 1, N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        DELYI = PT1*ABS(Y0(I)) + ATOLI
        AFI = ABS(YDOT(I))
        IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
 10     CONTINUE

! Set initial guess for h as geometric mean of upper and lower bounds. -
      ITER = 0
      HG = SQRT(HLB*HUB)
! If the bounds have crossed, exit with the mean value. ----------------
      IF (HUB .LT. HLB) THEN
        H0 = HG
        GO TO 90
      ENDIF

! Looping point for iteration. -----------------------------------------
 50   CONTINUE
! Estimate the second derivative as a difference quotient in f. --------
      T1 = T0 + HG
      DO 60 I = 1, N
 60     Y(I) = Y0(I) + HG*YDOT(I)
      CALL F (N, T1, Y, TEMP, RPAR, IPAR)
      DO 70 I = 1, N
 70     TEMP(I) = (TEMP(I) - YDOT(I))/HG
      YDDNRM = SVNORM (N, TEMP, EWT)
! Get the corresponding new value of h. --------------------------------
      IF (YDDNRM*HUB*HUB .GT. TWO) THEN
        HNEW = SQRT(TWO/YDDNRM)
      ELSE
        HNEW = SQRT(HG*HUB)
      ENDIF
      ITER = ITER + 1
!-----------------------------------------------------------------------
! Test the stopping conditions.
! Stop if the new and previous h values differ by a factor of .lt. 2.
! Stop if four iterations have been done.  Also, stop with previous h
! if HNEW/HG .gt. 2 after first iteration, as this probably means that
! the second derivative value is bad because of cancellation error.
!-----------------------------------------------------------------------
      IF (ITER .GE. 4) GO TO 80
      HRAT = HNEW/HG
      IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
      IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
        HNEW = HG
        GO TO 80
      ENDIF
      HG = HNEW
      GO TO 50

! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
 80   H0 = HNEW*HALF
      IF (H0 .LT. HLB) H0 = HLB
      IF (H0 .GT. HUB) H0 = HUB
 90   H0 = SIGN(H0, TOUT - T0)
      NITER = ITER
      IER = 0
      RETURN
! Error return for TOUT - T0 too small. --------------------------------
 100  IER = -1
      RETURN
!----------------------- End of Subroutine SVHIN -----------------------
      END SUBROUTINE SVHIN

      SUBROUTINE SVINDY (T, K, YH, LDYH, DKY, IFLAG)
      USE messy_clams_global, ONLY: prec
      USE messy_clamschem_asad_blas, ONLY: sscal
      REAL(PREC) T, YH, DKY
      INTEGER K, LDYH, IFLAG
      DIMENSION YH(LDYH,*), DKY(*)
!-----------------------------------------------------------------------
! Call sequence input -- T, K, YH, LDYH
! Call sequence output -- DKY, IFLAG
! COMMON block variables accessed..
!     /SVOD01/ --  H, TN, UROUND, L, N, NQ
!     /SVOD02/ --  HU
!
! Subroutines called by SVINDY.. SSCAL, XERRWV
! Function routines called by SVINDY.. None
!-----------------------------------------------------------------------
! SVINDY computes interpolated values of the K-th derivative of the
! dependent variable vector y, and stores it in DKY.  This routine
! is called within the package with K = 0 and T = TOUT, but may
! also be called by the user for any K up to the current order.
! (See detailed instructions in the usage documentation.)
!-----------------------------------------------------------------------
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH.  This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is..
!              q
!  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
!             j=K
! where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
! communicated by COMMON.  The above sum is done in reverse order.
! IFLAG is returned negative if either K or T is out of bounds.
!
! Discussion above and comments in driver explain all variables.
!-----------------------------------------------------------------------

! Type declarations for labeled COMMON block SVOD01 --------------------

      REAL(PREC) ACNRM, CCMXJ, CONP, CRATE, DRC, EL,  &
           ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
           RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
              L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
              LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
              N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
              NSLP, NYH

! Type declarations for labeled COMMON block SVOD02 --------------------

      REAL(PREC) HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

! Type declarations for local variables --------------------------------

      REAL(PREC) C, HUN, R, S, TFUZZ, TN1, TP, ZERO
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
! op_pj_20160830+
!!$   CHARACTER*80 MSG
      CHARACTER(LEN=80) MSG
! op_pj_20160830-
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE HUN, ZERO

      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),  &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,  &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
                      NSLP, NYH
      COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

      DATA HUN /100.0E0/, ZERO /0.0E0/

      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TFUZZ = HUN*UROUND*(TN + HU)
      TP = TN - HU - TFUZZ
      TN1 = TN + TFUZZ
      IF ((T-TP)*(T-TN1) .GT. ZERO) GO TO 90

      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1, NQ
 10     IC = IC*JJ
 15   C = REAL(IC)
      DO 20 I = 1, N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1, JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1, J
 30       IC = IC*JJ
 35     C = REAL(IC)
        DO 40 I = 1, N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      CALL SSCAL (N, R, DKY, 1)
      RETURN

 80   MSG = 'SVINDY-- K (=I1) illegal      '
      CALL XERRWV (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
      IFLAG = -1
      RETURN
 90   MSG = 'SVINDY-- T (=R1) illegal      '
      CALL XERRWV (MSG, 30, 52, 1, 0, 0, 0, 1, T, ZERO)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWV (MSG, 60, 52, 1, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
!----------------------- End of Subroutine SVINDY ----------------------
      END SUBROUTINE SVINDY

      SUBROUTINE SVSTEP (Y, YH, LDYH, YH1, EWT, SAVF, VSAV, ACOR,  &
                        WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR)
      USE messy_clams_global, ONLY: prec
      USE messy_clamschem_asad_blas, ONLY: scopy, sscal, saxpy
      EXTERNAL F, JAC, PSOL, VNLS
      REAL(PREC) Y, YH, YH1, EWT, SAVF, VSAV, ACOR, WM, RPAR
      INTEGER LDYH, IWM, IPAR
      DIMENSION Y(*), YH(LDYH,*), YH1(*), EWT(*), SAVF(*), VSAV(*),  &
         ACOR(*), WM(*), IWM(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,
!                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM
! COMMON block variables accessed..
!     /SVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
!               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
!               L, LMAX, MAXORD, MITER, N, NEWQ, NQ, NQWAIT
!     /SVOD02/  HU, NCFN, NETF, NFE, NQU, NST
!
! Subroutines called by SVSTEP.. F, SAXPY, SCOPY, SSCAL,
!                               SVJUST, VNLS, SVSET
! Function routines called by SVSTEP.. SVNORM
!-----------------------------------------------------------------------
! SVSTEP performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! SVSTEP calls subroutine VNLS for the solution of the nonlinear system
! arising in the time step.  Thus it is independent of the problem
! Jacobian structure and the type of nonlinear system solution method.
! SVSTEP returns a completion flag KFLAG (in COMMON).
! A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
! consecutive failures occurred.  On a return with KFLAG negative,
! the values of TN and the YH array are as of the beginning of the last
! step, and H is the last step size attempted.
!
! Communication with SVSTEP is done with the following variables..
!
! Y      = An array of length N used for the dependent variable vector.
! YH     = An LDYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! LDYH   = A constant integer .ge. N, the first dimension of YH.
!          N is the number of ODEs in the system.
! YH1    = A one-dimensional array occupying the same space as YH.
! EWT    = An array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = An array of working storage, of length N.
!          also used for input of YH(*,MAXORD+2) when JSTART = -1
!          and MAXORD .lt. the current order NQ.
! VSAV   = A work array of length N passed to subroutine VNLS.
! ACOR   = A work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = Real and integer work arrays associated with matrix
!          operations in VNLS.
! F      = Dummy name for the user supplied subroutine for f.
! JAC    = Dummy name for the user supplied Jacobian subroutine.
! PSOL   = Dummy name for the subroutine passed to VNLS, for
!          possible use there.
! VNLS   = Dummy name for the nonlinear system solving subroutine,
!          whose real name is dependent on the method used.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
!-----------------------------------------------------------------------

! Type declarations for labeled COMMON block SVOD01 --------------------

      REAL(PREC) ACNRM, CCMXJ, CONP, CRATE, DRC, EL,  &
           ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
           RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
              L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
              LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
              N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
              NSLP, NYH

! Type declarations for labeled COMMON block SVOD02 --------------------

      REAL(PREC) HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

! Type declarations for local variables --------------------------------

      REAL(PREC) ADDON, BIAS1,BIAS2,BIAS3, CNQUOT, DDN, DSM, DUP,  &
           ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF,  &
           ETAQ, ETAQM1, ETAQP1, FLOTL, ONE, ONEPSM,  &
           R, THRESH, TOLD, ZERO
      INTEGER I, I1, I2, IBACK, J, JB, KFC, KFH, MXNCF, NCF, NFLAG

! Type declaration for function subroutines called ---------------------

! ju_nt_20160617
!!$      REAL(PREC) SVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE ADDON, BIAS1, BIAS2, BIAS3,  &
           ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF,  &
           KFC, KFH, MXNCF, ONEPSM, THRESH, ONE, ZERO
!-----------------------------------------------------------------------
      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),  &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,  &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
                      NSLP, NYH
      COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

      DATA KFC/-3/, KFH/-7/, MXNCF/10/
      DATA ADDON  /1.0E-6/,    BIAS1  /6.0E0/,     BIAS2  /6.0E0/,  &
           BIAS3  /10.0E0/,    ETACF  /0.25E0/,    ETAMIN /0.1E0/,  &
           ETAMXF /0.2E0/,     ETAMX1 /1.0E4/,     ETAMX2 /10.0E0/,  &
           ETAMX3 /10.0E0/,    ONEPSM /1.00001E0/, THRESH /1.5E0/
      DATA ONE/1.0E0/, ZERO/0.0E0/

      KFLAG = 0
      TOLD = TN
      NCF = 0
      JCUR = 0
      NFLAG = 0
      IF (JSTART .GT. 0) GO TO 20
      IF (JSTART .EQ. -1) GO TO 100
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  ETAMAX is the maximum ratio by which H can be increased
! in a single step.  It is normally 1.5, but is larger during the
! first 10 steps to compensate for the small initial H.  If a failure
! occurs (in corrector convergence or error test), ETAMAX is set to 1
! for the next increase.
!-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      NQNYH = NQ*LDYH
      TAU(1) = H
      PRL1 = ONE
      RC = ZERO
      ETAMAX = ETAMX1
      NQWAIT = 2
      HSCAL = H
      GO TO 200
!-----------------------------------------------------------------------
! Take preliminary actions on a normal continuation step (JSTART.GT.0).
! If the driver changed H, then ETA must be reset and NEWH set to 1.
! If a change of order was dictated on the previous step, then
! it is done here and appropriate adjustments in the history are made.
! On an order decrease, the history array is adjusted by SVJUST.
! On an order increase, the history array is augmented by a column.
! On a change of step size H, the history array YH is rescaled.
!-----------------------------------------------------------------------
 20   CONTINUE
      IF (KUTH .EQ. 1) THEN
        ETA = MIN(ETA,H/HSCAL)
        NEWH = 1
        ENDIF
 50   IF (NEWH .EQ. 0) GO TO 200
      IF (NEWQ .EQ. NQ) GO TO 150
      IF (NEWQ .LT. NQ) THEN
        CALL SVJUST (YH, LDYH, -1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
        ENDIF
      IF (NEWQ .GT. NQ) THEN
        CALL SVJUST (YH, LDYH, 1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
      ENDIF
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! If N was reduced, zero out part of YH to avoid undefined references.
! If MAXORD was reduced to a value less than the tentative order NEWQ,
! then NQ is set to MAXORD, and a new H ratio ETA is chosen.
! Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
! In any case, NQWAIT is reset to L = NQ + 1 to prevent further
! changes in order for that many steps.
! The new H ratio ETA is limited by the input H if KUTH = 1,
! by HMIN if KUTH = 0, and by HMXI in any case.
! Finally, the history array YH is rescaled.
!-----------------------------------------------------------------------
 100  CONTINUE
      LMAX = MAXORD + 1
      IF (N .EQ. LDYH) GO TO 120
      I1 = 1 + (NEWQ + 1)*LDYH
      I2 = (MAXORD + 1)*LDYH
      IF (I1 .GT. I2) GO TO 120
      DO 110 I = I1, I2
 110    YH1(I) = ZERO
 120  IF (NEWQ .LE. MAXORD) GO TO 140
      FLOTL = REAL(LMAX)
      IF (MAXORD .LT. NQ-1) THEN
        DDN = SVNORM (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        ENDIF
      IF (MAXORD .EQ. NQ .AND. NEWQ .EQ. NQ+1) ETA = ETAQ
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ+1) THEN
        ETA = ETAQM1
        CALL SVJUST (YH, LDYH, -1)
        ENDIF
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ) THEN
        DDN = SVNORM (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        CALL SVJUST (YH, LDYH, -1)
        ENDIF
      ETA = MIN(ETA,ONE)
      NQ = MAXORD
      L = LMAX
 140  IF (KUTH .EQ. 1) ETA = MIN(ETA,ABS(H/HSCAL))
      IF (KUTH .EQ. 0) ETA = MAX(ETA,HMIN/ABS(HSCAL))
      ETA = ETA/MAX(ONE,ABS(HSCAL)*HMXI*ETA)
      NEWH = 1
      NQWAIT = L
      IF (NEWQ .LE. MAXORD) GO TO 50
! Rescale the history array for a change in H by a factor of ETA. ------
 150  R = ONE
      DO 180 J = 2, L
        R = R*ETA
        CALL SSCAL (N, R, YH(1,J), 1 )
 180    CONTINUE
      H = HSCAL*ETA
      HSCAL = H
      RC = RC*ETA
      NQNYH = NQ*LDYH
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! SVSET is called to calculate all integration coefficients.
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
!-----------------------------------------------------------------------
 200  TN = TN + H
      I1 = NQNYH + 1
      DO 220 JB = 1, NQ
        I1 = I1 - LDYH
        DO 210 I = I1, NQNYH
 210      YH1(I) = YH1(I) + YH1(I+LDYH)
 220  CONTINUE
      CALL SVSET
      RL1 = ONE/EL(2)
      RC = RC*(RL1/PRL1)
      PRL1 = RL1

! Call the nonlinear system solver. ------------------------------------

      CALL VNLS (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,  &
                 F, JAC, PSOL, NFLAG, RPAR, IPAR)

      IF (NFLAG .EQ. 0) GO TO 450
!-----------------------------------------------------------------------
! The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
! The YH array is retracted to its values before prediction.
! The step size H is reduced and the step is retried, if possible.
! Otherwise, an error exit is taken.
!-----------------------------------------------------------------------
        NCF = NCF + 1
        NCFN = NCFN + 1
        ETAMAX = ONE
        TN = TOLD
        I1 = NQNYH + 1
        DO 430 JB = 1, NQ
          I1 = I1 - LDYH
          DO 420 I = I1, NQNYH
 420        YH1(I) = YH1(I) - YH1(I+LDYH)
 430      CONTINUE
        IF (NFLAG .LT. -1) GO TO 680
        IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 670
        IF (NCF .EQ. MXNCF) GO TO 670
        ETA = ETACF
        ETA = MAX(ETA,HMIN/ABS(H))
        NFLAG = -1
        GO TO 150
!-----------------------------------------------------------------------
! The corrector has converged (NFLAG = 0).  The local error test is
! made and control passes to statement 500 if it fails.
!-----------------------------------------------------------------------
 450  CONTINUE
      DSM = ACNRM/TQ(2)
      IF (DSM .GT. ONE) GO TO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH and TAU arrays and decrement
! NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
! for use in a possible order increase on the next step.
! If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
!-----------------------------------------------------------------------
      KFLAG = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 IBACK = 1, NQ
        I = L - IBACK
 470    TAU(I+1) = TAU(I)
      TAU(1) = H
      DO 480 J = 1, L
        CALL SAXPY (N, EL(J), ACOR, 1, YH(1,J), 1 )
 480    CONTINUE
      NQWAIT = NQWAIT - 1
      IF ((L .EQ. LMAX) .OR. (NQWAIT .NE. 1)) GO TO 490
      CALL SCOPY (N, ACOR, 1, YH(1,LMAX), 1 )
      CONP = TQ(5)
 490  IF (ETAMAX .NE. ONE) GO TO 560
      IF (NQWAIT .LT. 2) NQWAIT = 2
      NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for the
! same order.  After repeated failures, H is forced to decrease
! more rapidly.
!-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      NETF = NETF + 1
      NFLAG = -2
      TN = TOLD
      I1 = NQNYH + 1
      DO 520 JB = 1, NQ
        I1 = I1 - LDYH
        DO 510 I = I1, NQNYH
 510      YH1(I) = YH1(I) - YH1(I+LDYH)
 520  CONTINUE
      IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 660
      ETAMAX = ONE
      IF (KFLAG .LE. KFC) GO TO 530
! Compute ratio of new H to current H at the current order. ------------
      FLOTL = REAL(L)
      ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      ETA = MAX(ETA,HMIN/ABS(H),ETAMIN)
      IF ((KFLAG .LE. -2) .AND. (ETA .GT. ETAMXF)) ETA = ETAMXF
      GO TO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more consecutive failures
! have occurred.  It is assumed that the elements of the YH array
! have accumulated errors of the wrong order.  The order is reduced
! by one, if possible.  Then H is reduced by a factor of 0.1 and
! the step is retried.  After a total of 7 consecutive failures,
! an exit is taken with KFLAG = -1.
!-----------------------------------------------------------------------
 530  IF (KFLAG .EQ. KFH) GO TO 660
      IF (NQ .EQ. 1) GO TO 540
      ETA = MAX(ETAMIN,HMIN/ABS(H))
      CALL SVJUST (YH, LDYH, -1)
      L = NQ
      NQ = NQ - 1
      NQWAIT = L
      GO TO 150
 540  ETA = MAX(ETAMIN,HMIN/ABS(H))
      H = H*ETA
      HSCAL = H
      TAU(1) = H
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      DO 550 I = 1, N
 550    YH(I,2) = H*SAVF(I)
      NQWAIT = 10
      GO TO 200
!-----------------------------------------------------------------------
! If NQWAIT = 0, an increase or decrease in order by one is considered.
! Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
! be multiplied at order q, q-1, or q+1, respectively.
! The largest of these is determined, and the new order and
! step size set accordingly.
! A change of H or NQ is made only if H increases by at least a
! factor of THRESH.  If an order change is considered and rejected,
! then NQWAIT is set to 2 (reconsider it after 2 steps).
!-----------------------------------------------------------------------
! Compute ratio of new H to current H at the current order. ------------
 560  FLOTL = REAL(L)
      ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      IF (NQWAIT .NE. 0) GO TO 600
      NQWAIT = 2
      ETAQM1 = ZERO
      IF (NQ .EQ. 1) GO TO 570
! Compute ratio of new H to current H at the current order less one. ---
      DDN = SVNORM (N, YH(1,L), EWT)/TQ(1)
      ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
 570  ETAQP1 = ZERO
      IF (L .EQ. LMAX) GO TO 580
! Compute ratio of new H to current H at current order plus one. -------
      CNQUOT = (TQ(5)/CONP)*(H/TAU(2))**L
      DO 575 I = 1, N
 575    SAVF(I) = ACOR(I) - CNQUOT*YH(I,LMAX)
      DUP = SVNORM (N, SAVF, EWT)/TQ(3)
      ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
 580  IF (ETAQ .GE. ETAQP1) GO TO 590
      IF (ETAQP1 .GT. ETAQM1) GO TO 620
      GO TO 610
 590  IF (ETAQ .LT. ETAQM1) GO TO 610
 600  ETA = ETAQ
      NEWQ = NQ
      GO TO 630
 610  ETA = ETAQM1
      NEWQ = NQ - 1
      GO TO 630
 620  ETA = ETAQP1
      NEWQ = NQ + 1
      CALL SCOPY (N, ACOR, 1, YH(1,LMAX), 1)
! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
 630  IF (ETA .LT. THRESH .OR. ETAMAX .EQ. ONE) GO TO 640
      ETA = MIN(ETA,ETAMAX)
      ETA = ETA/MAX(ONE,ABS(H)*HMXI*ETA)
      NEWH = 1
      HNEW = H*ETA
      GO TO 690
 640  NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
!-----------------------------------------------------------------------
! All returns are made through this section.
! On a successful return, ETAMAX is reset and ACOR is scaled.
!-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  IF (NFLAG .EQ. -2) KFLAG = -3
      IF (NFLAG .EQ. -3) KFLAG = -4
      GO TO 720
 690  ETAMAX = ETAMX3
      IF (NST .LE. 10) ETAMAX = ETAMX2
 700  R = ONE/TQ(2)
      CALL SSCAL (N, R, ACOR, 1)
 720  JSTART = 1
      RETURN
!----------------------- End of Subroutine SVSTEP ----------------------
      END SUBROUTINE SVSTEP

      SUBROUTINE SVSET
!-----------------------------------------------------------------------
! Call sequence communication.. None
! COMMON block variables accessed..
!     /SVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
!                 METH, NQ, NQWAIT
!
! Subroutines called by SVSET.. None
! Function routines called by SVSET.. None
!-----------------------------------------------------------------------
! SVSET is called by SVSTEP and sets coefficients for use there.
!
! For each order NQ, the coefficients in EL are calculated by use of
!  the generating polynomial lambda(x), with coefficients EL(i).
!      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
! For the backward differentiation formulas,
!                                     NQ-1
!      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
!                                     i = 1
! For the Adams formulas,
!                              NQ-1
!      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
!                              i = 1
!      lambda(-1) = 0,    lambda(0) = 1,
! where c is a normalization constant.
! In both cases, xi(i) is defined by
!      H*xi(i) = t sub n  -  t sub (n-i)
!              = H + TAU(1) + TAU(2) + ... TAU(i-1).
!
!
! In addition to variables described previously, communication
! with SVSET uses the following..
!   TAU    = A vector of length 13 containing the past NQ values
!            of H.
!   EL     = A vector of length 13 in which vset stores the
!            coefficients for the corrector formula.
!   TQ     = A vector of length 5 in which vset stores constants
!            used for the convergence test, the error test, and the
!            selection of H at a new order.
!   METH   = The basic method indicator.
!   NQ     = The current order.
!   L      = NQ + 1, the length of the vector stored in EL, and
!            the number of columns of the YH array being used.
!   NQWAIT = A counter controlling the frequency of order changes.
!            An order change is about to be considered if NQWAIT = 1.
!-----------------------------------------------------------------------

! Type declarations for labeled COMMON block SVOD01 --------------------

      USE messy_clams_global, ONLY: prec
      REAL(PREC) ACNRM, CCMXJ, CONP, CRATE, DRC, EL,  &
           ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
           RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
              L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
              LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
              N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
              NSLP, NYH

! Type declarations for local variables --------------------------------

      REAL(PREC) AHATN0, ALPH0, CNQM1, CORTES, CSUM, ELP, EM,  &
           EM0, FLOTI, FLOTL, FLOTNQ, HSUM, ONE, RXI, RXIS, S, SIX,  &
           T1, T2, T3, T4, T5, T6, TWO, XI, ZERO
      INTEGER I, IBACK, J, JP1, NQM1, NQM2

      DIMENSION EM(13)
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE CORTES, ONE, SIX, TWO, ZERO

      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),  &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,  &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
                      NSLP, NYH

      DATA CORTES /0.1E0/
      DATA ONE  /1.0E0/, SIX /6.0E0/, TWO /2.0E0/, ZERO /0.0E0/

      FLOTL = REAL(L)
      NQM1 = NQ - 1
      NQM2 = NQ - 2
      GO TO (100, 200), METH

! Set coefficients for Adams methods. ----------------------------------
 100  IF (NQ .NE. 1) GO TO 110
      EL(1) = ONE
      EL(2) = ONE
      TQ(1) = ONE
      TQ(2) = TWO
      TQ(3) = SIX*TQ(2)
      TQ(5) = ONE
      GO TO 300
 110  HSUM = H
      EM(1) = ONE
      FLOTNQ = FLOTL - ONE
      DO 115 I = 2, L
 115    EM(I) = ZERO
      DO 150 J = 1, NQM1
        IF ((J .NE. NQM1) .OR. (NQWAIT .NE. 1)) GO TO 130
        S = ONE
        CSUM = ZERO
        DO 120 I = 1, NQM1
          CSUM = CSUM + S*EM(I)/REAL(I+1)
 120      S = -S
        TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
 130    RXI = H/HSUM
        DO 140 IBACK = 1, J
          I = (J + 2) - IBACK
 140      EM(I) = EM(I) + EM(I-1)*RXI
        HSUM = HSUM + TAU(J)
 150    CONTINUE
! Compute integral from -1 to 0 of polynomial and of x times it. -------
      S = ONE
      EM0 = ZERO
      CSUM = ZERO
      DO 160 I = 1, NQ
        FLOTI = REAL(I)
        EM0 = EM0 + S*EM(I)/FLOTI
        CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
 160    S = -S
! In EL, form coefficients of normalized integrated polynomial. --------
      S = ONE/EM0
      EL(1) = ONE
      DO 170 I = 1, NQ
 170    EL(I+1) = S*EM(I)/REAL(I)
      XI = HSUM/H
      TQ(2) = XI*EM0/CSUM
      TQ(5) = XI/EL(L)
      IF (NQWAIT .NE. 1) GO TO 300
! For higher order control constant, multiply polynomial by 1+x/xi(q). -
      RXI = ONE/XI
      DO 180 IBACK = 1, NQ
        I = (L + 1) - IBACK
 180    EM(I) = EM(I) + EM(I-1)*RXI
! Compute integral of polynomial. --------------------------------------
      S = ONE
      CSUM = ZERO
      DO 190 I = 1, L
        CSUM = CSUM + S*EM(I)/REAL(I+1)
 190    S = -S
      TQ(3) = FLOTL*EM0/CSUM
      GO TO 300

! Set coefficients for BDF methods. ------------------------------------
 200  DO 210 I = 3, L
 210    EL(I) = ZERO
      EL(1) = ONE
      EL(2) = ONE
      ALPH0 = -ONE
      AHATN0 = -ONE
      HSUM = H
      RXI = ONE
      RXIS = ONE
      IF (NQ .EQ. 1) GO TO 240
      DO 230 J = 1, NQM2
! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
        HSUM = HSUM + TAU(J)
        RXI = H/HSUM
        JP1 = J + 1
        ALPH0 = ALPH0 - ONE/REAL(JP1)
        DO 220 IBACK = 1, JP1
          I = (J + 3) - IBACK
 220      EL(I) = EL(I) + EL(I-1)*RXI
 230    CONTINUE
      ALPH0 = ALPH0 - ONE/REAL(NQ)
      RXIS = -EL(2) - ALPH0
      HSUM = HSUM + TAU(NQM1)
      RXI = H/HSUM
      AHATN0 = -EL(2) - RXI
      DO 235 IBACK = 1, NQ
        I = (NQ + 2) - IBACK
 235    EL(I) = EL(I) + EL(I-1)*RXIS
 240  T1 = ONE - AHATN0 + ALPH0
      T2 = ONE + REAL(NQ)*T1
      TQ(2) = ABS(ALPH0*T2/T1)
      TQ(5) = ABS(T2/(EL(L)*RXI/RXIS))
      IF (NQWAIT .NE. 1) GO TO 300
      CNQM1 = RXIS/EL(L)
      T3 = ALPH0 + ONE/REAL(NQ)
      T4 = AHATN0 + RXI
      ELP = T3/(ONE - T4 + T3)
      TQ(1) = ABS(ELP/CNQM1)
      HSUM = HSUM + TAU(NQ)
      RXI = H/HSUM
      T5 = ALPH0 - ONE/REAL(NQ+1)
      T6 = AHATN0 - RXI
      ELP = T2/(ONE - T6 + T5)
      TQ(3) = ABS(ELP*RXI*(FLOTL + ONE)*T5)
 300  TQ(4) = CORTES*TQ(2)
      RETURN
!----------------------- End of Subroutine SVSET -----------------------
      END SUBROUTINE SVSET

      SUBROUTINE SVJUST (YH, LDYH, IORD)
      USE messy_clams_global, ONLY: prec
      USE messy_clamschem_asad_blas, ONLY: saxpy
      REAL(PREC) YH
      INTEGER LDYH, IORD
      DIMENSION YH(LDYH,*)
!-----------------------------------------------------------------------
! Call sequence input -- YH, LDYH, IORD
! Call sequence output -- YH
! COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
! COMMON block variables accessed..
!     /SVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
!
! Subroutines called by SVJUST.. SAXPY
! Function routines called by SVJUST.. None
!-----------------------------------------------------------------------
! This subroutine adjusts the YH array on reduction of order,
! and also when the order is increased for the stiff option (METH = 2).
! Communication with SVJUST uses the following..
! IORD  = An integer flag used when METH = 2 to indicate an order
!         increase (IORD = +1) or an order decrease (IORD = -1).
! HSCAL = Step size H used in scaling of Nordsieck array YH.
!         (If IORD = +1, SVJUST assumes that HSCAL = TAU(1).)
! See References 1 and 2 for details.
!-----------------------------------------------------------------------

! Type declarations for labeled COMMON block SVOD01 --------------------

      REAL(PREC) ACNRM, CCMXJ, CONP, CRATE, DRC, EL,  &
           ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
           RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
              L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
              LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
              N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
              NSLP, NYH

! Type declarations for local variables --------------------------------

      REAL(PREC) ALPH0, ALPH1, HSUM, ONE, PROD, T1, XI,XIOLD, ZERO
      INTEGER I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE ONE, ZERO

      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),  &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,  &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
                      NSLP, NYH

      DATA ONE /1.0E0/, ZERO /0.0E0/

      IF ((NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN
      NQM1 = NQ - 1
      NQM2 = NQ - 2
      GO TO (100, 200), METH
!-----------------------------------------------------------------------
! Nonstiff option...
! Check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
 100  CONTINUE
      IF (IORD .EQ. 1) GO TO 180
! Order decrease. ------------------------------------------------------
      DO 110 J = 1, LMAX
 110    EL(J) = ZERO
      EL(2) = ONE
      HSUM = ZERO
      DO 130 J = 1, NQM2
! Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
        HSUM = HSUM + TAU(J)
        XI = HSUM/HSCAL
        JP1 = J + 1
        DO 120 IBACK = 1, JP1
          I = (J + 3) - IBACK
 120      EL(I) = EL(I)*XI + EL(I-1)
 130    CONTINUE
! Construct coefficients of integrated polynomial. ---------------------
      DO 140 J = 2, NQM1
 140    EL(J+1) = REAL(NQ)*EL(J)/REAL(J)
! Subtract correction terms from YH array. -----------------------------
      DO 170 J = 3, NQ
        DO 160 I = 1, N
 160      YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
 170    CONTINUE
      RETURN
! Order increase. ------------------------------------------------------
! Zero out next column in YH array. ------------------------------------
 180  CONTINUE
      LP1 = L + 1
      DO 190 I = 1, N
 190    YH(I,LP1) = ZERO
      RETURN
!-----------------------------------------------------------------------
! Stiff option...
! Check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
 200  CONTINUE
      IF (IORD .EQ. 1) GO TO 300
! Order decrease. ------------------------------------------------------
      DO 210 J = 1, LMAX
 210    EL(J) = ZERO
      EL(3) = ONE
      HSUM = ZERO
      DO 230 J = 1,NQM2
! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        HSUM = HSUM + TAU(J)
        XI = HSUM/HSCAL
        JP1 = J + 1
        DO 220 IBACK = 1, JP1
          I = (J + 4) - IBACK
 220      EL(I) = EL(I)*XI + EL(I-1)
 230    CONTINUE
! Subtract correction terms from YH array. -----------------------------
      DO 250 J = 3,NQ
        DO 240 I = 1, N
 240      YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
 250    CONTINUE
      RETURN
! Order increase. ------------------------------------------------------
 300  DO 310 J = 1, LMAX
 310    EL(J) = ZERO
      EL(3) = ONE
      ALPH0 = -ONE
      ALPH1 = ONE
      PROD = ONE
      XIOLD = ONE
      HSUM = HSCAL
      IF (NQ .EQ. 1) GO TO 340
      DO 330 J = 1, NQM1
! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        JP1 = J + 1
        HSUM = HSUM + TAU(JP1)
        XI = HSUM/HSCAL
        PROD = PROD*XI
        ALPH0 = ALPH0 - ONE/REAL(JP1)
        ALPH1 = ALPH1 + ONE/XI
        DO 320 IBACK = 1, JP1
          I = (J + 4) - IBACK
 320      EL(I) = EL(I)*XIOLD + EL(I-1)
        XIOLD = XI
 330    CONTINUE
 340  CONTINUE
      T1 = (-ALPH0 - ALPH1)/PROD
! Load column L + 1 in YH array. ---------------------------------------
      LP1 = L + 1
      DO 350 I = 1, N
 350    YH(I,LP1) = T1*YH(I,LMAX)
! Add correction terms to YH array. ------------------------------------
      NQP1 = NQ + 1
      DO 370 J = 3, NQP1
        CALL SAXPY (N, EL(J), YH(1,LP1), 1, YH(1,J), 1 )
 370  CONTINUE
      RETURN
!----------------------- End of Subroutine SVJUST ----------------------
      END SUBROUTINE SVJUST

      SUBROUTINE SVNLSD (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,  &
                       F, JAC, PDUM, NFLAG, RPAR, IPAR)
      USE messy_clams_global, ONLY: prec
      USE messy_clamschem_asad_blas, ONLY: scopy, saxpy, sscal
      EXTERNAL F, JAC, PDUM
      REAL(PREC) Y, YH, VSAV, SAVF, EWT, ACOR, WM, RPAR
      INTEGER LDYH, IWM, NFLAG, IPAR
      DIMENSION Y(*), YH(LDYH,*), VSAV(*), SAVF(*), EWT(*), ACOR(*),  &
                IWM(*), WM(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, SAVF, EWT, ACOR, IWM, WM,
!                        F, JAC, NFLAG, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM, NFLAG
! COMMON block variables accessed..
!     /SVOD01/ ACNRM, CRATE, DRC, H, RC, RL1, TQ(5), TN, ICF,
!                JCUR, METH, MITER, N, NSLP
!     /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Subroutines called by SVNLSD.. F, SAXPY, SCOPY, SSCAL, SVJAC, SVSOL
! Function routines called by SVNLSD.. SVNORM
!-----------------------------------------------------------------------
! Subroutine SVNLSD is a nonlinear system solver, which uses functional
! iteration or a chord (modified Newton) method.  For the chord method
! direct linear algebraic system solvers are used.  Subroutine SVNLSD
! then handles the corrector phase of this integration package.
!
! Communication with SVNLSD is done with the following variables. (For
! more details, please see the comments in the driver subroutine.)
!
! Y          = The dependent variable, a vector of length N, input.
! YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
!              and output.  On input, it contains predicted values.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! VSAV       = Unused work array.
! SAVF       = A work array of length N.
! EWT        = An error weight vector of length N, input.
! ACOR       = A work array of length N, used for the accumulated
!              corrections to the predicted y vector.
! WM,IWM     = Real and integer work arrays associated with matrix
!              operations in chord iteration (MITER .ne. 0).
! F          = Dummy name for user supplied routine for f.
! JAC        = Dummy name for user supplied Jacobian routine.
! PDUM       = Unused dummy subroutine name.  Included for uniformity
!              over collection of integrators.
! NFLAG      = Input/output flag, with values and meanings as follows..
!              INPUT
!                  0 first call for this time step.
!                 -1 convergence failure in previous call to SVNLSD.
!                 -2 error test failure in SVSTEP.
!              OUTPUT
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or singular matrix.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in solution (cannot occur
!                    here).
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
!
! IPUP       = Own variable flag with values and meanings as follows..
!              0,            do not update the Newton matrix.
!              MITER .ne. 0, update Newton matrix, because it is the
!                            initial step, order was changed, the error
!                            test failed, or an update is indicated by
!                            the scalar RC or step counter NST.
!
! For more details, see comments in driver subroutine.
!-----------------------------------------------------------------------
! Type declarations for labeled COMMON block SVOD01 --------------------

      REAL(PREC) ACNRM, CCMXJ, CONP, CRATE, DRC, EL,  &
           ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
           RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
              L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
              LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
              N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
              NSLP, NYH  

! Type declarations for labeled COMMON block SVOD02 --------------------

      REAL(PREC) HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

! Type declarations for local variables --------------------------------

      REAL(PREC) CCMAX, CRDOWN, CSCALE, DCON, DEL, DELP, ONE,  &
           RDIV, TWO, ZERO
      INTEGER I, IERPJ, IERSL, M, MAXCOR, MSBP

! Type declaration for function subroutines called ---------------------

! ju_nt_20160617
!!$      REAL(PREC) SVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE CCMAX, CRDOWN, MAXCOR, MSBP, RDIV, ONE, TWO, ZERO

      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),  &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,  &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
                      NSLP, NYH
      COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

      DATA CCMAX /0.3E0/, CRDOWN /0.3E0/, MAXCOR /3/, MSBP /20/,  &
           RDIV  /2.0E0/
      DATA ONE /1.0E0/, TWO /2.0E0/, ZERO /0.0E0/
!-----------------------------------------------------------------------
! On the first step, on a change of method order, or after a
! nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
! to force a Jacobian update when MITER .ne. 0.
!-----------------------------------------------------------------------
      IF (JSTART .EQ. 0) NSLP = 0
      IF (NFLAG .EQ. 0) ICF = 0
      IF (NFLAG .EQ. -2) IPUP = MITER
      IF ( (JSTART .EQ. 0) .OR. (JSTART .EQ. -1) ) IPUP = MITER
! If this is functional iteration, set CRATE .eq. 1 and drop to 220
      IF (MITER .EQ. 0) THEN
        CRATE = ONE
        GO TO 220
      ENDIF
!-----------------------------------------------------------------------
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force SVJAC to be called, if a Jacobian is involved.
! In any case, SVJAC is called at least every MSBP steps.
!-----------------------------------------------------------------------
      DRC = ABS(RC-ONE)
      IF (DRC .GT. CCMAX .OR. NST .GE. NSLP+MSBP) IPUP = MITER
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(i).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
 220  M = 0
      DELP = ZERO
      CALL SCOPY (N, YH(1,1), 1, Y, 1 )
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
!-----------------------------------------------------------------------
! If indicated, the matrix P = I - h*rl1*J is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
      CALL SVJAC (Y, YH, LDYH, EWT, ACOR, SAVF, WM, IWM, F, JAC, IERPJ,  &
                 RPAR, IPAR)
      IPUP = 0
      RC = ONE
      DRC = ZERO
      CRATE = ONE
      NSLP = NST
! If matrix is singular, take error return to force cut in step size. --
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = ZERO
! This is a looping point for the corrector iteration. -----------------
 270  IF (MITER .NE. 0) GO TO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
      DO 280 I = 1,N
 280    SAVF(I) = RL1*(H*SAVF(I) - YH(I,2))
      DO 290 I = 1,N
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = SVNORM (N, Y, EWT)
      DO 300 I = 1,N
 300    Y(I) = YH(I,1) + SAVF(I)
      CALL SCOPY (N, SAVF, 1, ACOR, 1)
      GO TO 400
!-----------------------------------------------------------------------
! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix.  The correction is scaled by the factor
! 2/(1+RC) to account for changes in h*rl1 since the last SVJAC call.
!-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = (RL1*H)*SAVF(I) - (RL1*YH(I,2) + ACOR(I))
      CALL SVSOL (WM, IWM, Y, IERSL)
      NNI = NNI + 1
      IF (IERSL .GT. 0) GO TO 410
      IF (METH .EQ. 2 .AND. RC .NE. ONE) THEN
        CSCALE = TWO/(ONE + RC)
        CALL SSCAL (N, CSCALE, Y, 1)
      ENDIF
      DEL = SVNORM (N, Y, EWT)
      CALL SAXPY (N, ONE, Y, 1, ACOR, 1)
      DO 380 I = 1,N
 380    Y(I) = YH(I,1) + ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M .gt. 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(CRDOWN*CRATE,DEL/DELP)
      DCON = DEL*MIN(ONE,CRATE)/TQ(4)
      IF (DCON .LE. ONE) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. RDIV*DELP) GO TO 410
      DELP = DEL
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      GO TO 270

 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220

 430  CONTINUE
      NFLAG = -1
      ICF = 2
      IPUP = MITER
      RETURN

! Return for successful step. ------------------------------------------
 450  NFLAG = 0
      JCUR = 0
      ICF = 0
      IF (M .EQ. 0) ACNRM = DEL
      IF (M .GT. 0) ACNRM = SVNORM (N, ACOR, EWT)
      RETURN
!----------------------- End of Subroutine SVNLSD ----------------------
      END SUBROUTINE SVNLSD

      SUBROUTINE SVJAC (Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM, F, JAC,  &
                       IERPJ, RPAR, IPAR)
      USE messy_clams_global, ONLY: prec
      USE messy_clamschem_asad_linpack, ONLY: sgefa, sgbfa
      USE messy_clamschem_asad_blas, ONLY: scopy, sscal
      EXTERNAL F, JAC
      REAL(PREC) Y, YH, EWT, FTEM, SAVF, WM, RPAR
      INTEGER LDYH, IWM, IERPJ, IPAR
      DIMENSION Y(*), YH(LDYH,*), EWT(*), FTEM(*), SAVF(*),  &
         WM(*), IWM(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, EWT, FTEM, SAVF, WM, IWM,
!                        F, JAC, RPAR, IPAR
! Call sequence output -- WM, IWM, IERPJ
! COMMON block variables accessed..
!     /SVOD01/  CCMXJ, DRC, H, RL1, TN, UROUND, ICF, JCUR, LOCJS,
!               MSBJ, NSLJ
!     /SVOD02/  NFE, NST, NJE, NLU
!
! Subroutines called by SVJAC.. F, JAC, SACOPY, SCOPY, SGBFA, SGEFA,
!                              SSCAL
! Function routines called by SVJAC.. SVNORM
!-----------------------------------------------------------------------
! SVJAC is called by SVSTEP to compute and process the matrix
! P = I - h*rl1*J , where J is an approximation to the Jacobian.
! Here J is computed by the user-supplied routine JAC if
! MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
! If MITER = 3, a diagonal approximation to J is used.
! If JSV = -1, J is computed from scratch in all cases.
! If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
! considered acceptable, then P is constructed from the saved J.
! J is stored in wm and replaced by P.  If MITER .ne. 3, P is then
! subjected to LU decomposition in preparation for later solution
! of linear systems with P as coefficient matrix. This is done
! by SGEFA if MITER = 1 or 2, and by SGBFA if MITER = 4 or 5.
!
! Communication with SVJA! is done with the following variables.  (For
! more details, please see the comments in the driver subroutine.)
! Y          = Vector containing predicted values on entry.
! YH         = The Nordsieck array, an LDYH by LMAX array, input.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! EWT        = An error weight vector of length N.
! SAVF       = Array containing f evaluated at predicted y, input.
! WM         = Real work space for matrices.  In the output, it containS
!              the inverse diagonal matrix if MITER = 3 and the LU
!              decomposition of P if MITER is 1, 2 , 4, or 5.
!              Storage of matrix elements starts at WM(3).
!              Storage of the saved Jacobian starts at WM(LOCJS).
!              WM also contains the following matrix-related data..
!              WM(1) = SQRT(UROUND), used in numerical Jacobian step.
!              WM(2) = H*RL1, saved for later use if MITER = 3.
! IWM        = Integer work space containing pivot information,
!              starting at IWM(31), if MITER is 1, 2, 4, or 5.
!              IWM also contains band parameters ML = IWM(1) and
!              MU = IWM(2) if MITER is 4 or 5.
! F          = Dummy name for the user supplied subroutine for f.
! JAC        = Dummy name for the user supplied Jacobian subroutine.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
! RL1        = 1/EL(2) (input).
! IERPJ      = Output error flag,  = 0 if no trouble, 1 if the P
!              matrix is found to be singular.
! JCUR       = Output flag to indicate whether the Jacobian matrix
!              (or approximation) is now current.
!              JCUR = 0 means J is not current.
!              JCUR = 1 means J is current.
!-----------------------------------------------------------------------

! Type declarations for labeled COMMON block SVOD01 --------------------

      REAL(PREC) ACNRM, CCMXJ, CONP, CRATE, DRC, EL,  &
           ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
           RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
              L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
              LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
              N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
              NSLP, NYH

! Type declarations for labeled COMMON block SVOD02 --------------------

      REAL(PREC) HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

! Type declarations for local variables --------------------------------

      REAL(PREC) CON, DI, FAC, HRL1, ONE, PT1, R, R0, SRUR, THOU,  &
           YI, YJ, YJJ, ZERO
      INTEGER I, I1, I2, IER, II, J, J1, JJ, JOK, LENP, MBA, MBAND,  &
              MEB1, MEBAND, ML, ML3, MU, NP1

! Type declaration for function subroutines called ---------------------

! ju_nt_20160617
!!$      REAL(PREC) SVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this subroutine.
!-----------------------------------------------------------------------
      SAVE ONE, PT1, THOU, ZERO
!-----------------------------------------------------------------------
      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),  &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,  &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
                      NSLP, NYH
      COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST

      DATA ONE /1.0E0/, THOU /1000.0E0/, ZERO /0.0E0/, PT1 /0.1E0/

      IERPJ = 0
      HRL1 = H*RL1
! See whether J should be evaluated (JOK = -1) or not (JOK = 1). -------
      JOK = JSV
      IF (JSV .EQ. 1) THEN
        IF (NST .EQ. 0 .OR. NST .GT. NSLJ+MSBJ) JOK = -1
        IF (ICF .EQ. 1 .AND. DRC .LT. CCMXJ) JOK = -1
        IF (ICF .EQ. 2) JOK = -1
      ENDIF
! End of setting JOK. --------------------------------------------------

      IF (JOK .EQ. -1 .AND. MITER .EQ. 1) THEN
! If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian. ------------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = ZERO
      CALL JAC (N, TN, Y, 0, 0, WM(3), N, RPAR, IPAR)
      IF (JSV .EQ. 1) CALL SCOPY (LENP, WM(3), 1, WM(LOCJS), 1)
      ENDIF

      IF (JOK .EQ. -1 .AND. MITER .EQ. 2) THEN
! If MITER = 2, make N calls to F to approximate the Jacobian. ---------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      FAC = SVNORM (N, SAVF, EWT)
      R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
      IF (R0 .EQ. ZERO) R0 = ONE
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = ONE/R
        CALL F (N, TN, Y, FTEM, RPAR, IPAR)
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
      LENP = N*N
      IF (JSV .EQ. 1) CALL SCOPY (LENP, WM(3), 1, WM(LOCJS), 1)
      ENDIF

      IF (JOK .EQ. 1 .AND. (MITER .EQ. 1 .OR. MITER .EQ. 2)) THEN
      JCUR = 0
      LENP = N*N
      CALL SCOPY (LENP, WM(LOCJS), 1, WM(3), 1)
      ENDIF

      IF (MITER .EQ. 1 .OR. MITER .EQ. 2) THEN
! Multiply Jacobian by scalar, add identity, and do LU decomposition. --
      CON = -HRL1
      CALL SSCAL (LENP, CON, WM(3), 1)
      J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + ONE
 250    J = J + NP1
      NLU = NLU + 1
      CALL SGEFA (WM(3), N, N, IWM(31), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
      ENDIF
! End of code block for MITER = 1 or 2. --------------------------------

      IF (MITER .EQ. 3) THEN
! If MITER = 3, construct a diagonal approximation to J and P. ---------
      NJE = NJE + 1
      JCUR = 1
      WM(2) = HRL1
      R = RL1*PT1
      DO 310 I = 1,N
 310    Y(I) = Y(I) + R*(H*SAVF(I) - YH(I,2))
      CALL F (N, TN, Y, WM(3), RPAR, IPAR)
      NFE = NFE + 1
      DO 320 I = 1,N
        R0 = H*SAVF(I) - YH(I,2)
        DI = PT1*R0 - H*(WM(I+2) - SAVF(I))
        WM(I+2) = ONE
        IF (ABS(R0) .LT. UROUND/EWT(I)) GO TO 320
        IF (ABS(DI) .EQ. ZERO) GO TO 330
        WM(I+2) = PT1*R0/DI
 320    CONTINUE
      RETURN
 330  IERPJ = 1
      RETURN
      ENDIF
! End of code block for MITER = 3. -------------------------------------

! Set constants for MITER = 4 or 5. ------------------------------------
      ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N

      IF (JOK .EQ. -1 .AND. MITER .EQ. 4) THEN
! If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian. ------------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      DO 410 I = 1,LENP
 410    WM(I+2) = ZERO
      CALL JAC (N, TN, Y, ML, MU, WM(ML3), MEBAND, RPAR, IPAR)
      IF (JSV .EQ. 1)  &
         CALL SACOPY (MBAND, N, WM(ML3), MEBAND, WM(LOCJS), MBAND)
      ENDIF

      IF (JOK .EQ. -1 .AND. MITER .EQ. 5) THEN
! If MITER = 5, make N calls to F to approximate the Jacobian. ---------
      NJE = NJE + 1
      NSLJ = NST
      JCUR = 1
      MBA = MIN(MBAND,N)
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = SVNORM (N, SAVF, EWT)
      R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
      IF (R0 .EQ. ZERO) R0 = ONE
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (N, TN, Y, FTEM, RPAR, IPAR)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = ONE/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
      IF (JSV .EQ. 1)  &
         CALL SACOPY (MBAND, N, WM(ML3), MEBAND, WM(LOCJS), MBAND)
      ENDIF

      IF (JOK .EQ. 1) THEN
      JCUR = 0
      CALL SACOPY (MBAND, N, WM(LOCJS), MBAND, WM(ML3), MEBAND)
      ENDIF

! Multiply Jacobian by scalar, add identity, and do LU decomposition.
      CON = -HRL1
      CALL SSCAL (LENP, CON, WM(3), 1 )
      II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + ONE
 580    II = II + MEBAND
      NLU = NLU + 1
      CALL SGBFA (WM(3), MEBAND, N, ML, MU, IWM(31), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
! End of code block for MITER = 4 or 5. --------------------------------

!----------------------- End of Subroutine SVJAC -----------------------
      END SUBROUTINE SVJAC

      SUBROUTINE SACOPY (NROW, NCOL, A, NROWA, B, NROWB)
      USE messy_clams_global, ONLY: prec
      USE messy_clamschem_asad_blas, ONLY: scopy
      REAL(PREC) A, B
      INTEGER NROW, NCOL, NROWA, NROWB
      DIMENSION A(NROWA,NCOL), B(NROWB,NCOL)
!-----------------------------------------------------------------------
! Call sequence input -- NROW, NCOL, A, NROWA, NROWB
! Call sequence output -- B
! COMMON block variables accessed -- None
!
! Subroutines called by SACOPY.. SCOPY
! Function routines called by SACOPY.. None
!-----------------------------------------------------------------------
! This routine copies one rectangular array, A, to another, B,
! where A and B may have different row dimensions, NROWA and NROWB.
! The data copied consists of NROW rows and NCOL columns.
!-----------------------------------------------------------------------
      INTEGER IC

      DO 20 IC = 1,NCOL
        CALL SCOPY (NROW, A(1,IC), 1, B(1,IC), 1)
 20     CONTINUE

      RETURN
!----------------------- End of Subroutine SACOPY ----------------------
      END SUBROUTINE SACOPY

      SUBROUTINE SVSOL (WM, IWM, X, IERSL)
      USE messy_clams_global, ONLY: prec
      USE messy_clamschem_asad_linpack, ONLY: sgbsl, sgesl
      REAL(PREC) WM, X
      INTEGER IWM, IERSL
      DIMENSION WM(*), IWM(*), X(*)
!-----------------------------------------------------------------------
! Call sequence input -- WM, IWM, X
! Call sequence output -- X, IERSL
! COMMON block variables accessed..
!     /SVOD01/ -- H, RL1, MITER, N
!
! Subroutines called by SVSOL.. SGESL, SGBSL
! Function routines called by SVSOL.. None
!-----------------------------------------------------------------------
! This routine manages the solution of the linear system arising from
! a chord iteration.  It is called if MITER .ne. 0.
! If MITER is 1 or 2, it calls SGESL to accomplish this.
! If MITER = 3 it updates the coefficient H*RL1 in the diagonal
! matrix, and then computes the solution.
! If MITER is 4 or 5, it calls SGBSL.
! Communication with SVSOL uses the following variables..
! WM    = Real work space containing the inverse diagonal matrix if
!         MITER = 3 and the LU decomposition of the matrix otherwise.
!         Storage of matrix elements starts at WM(3).
!         WM also contains the following matrix-related data..
!         WM(1) = SQRT(UROUND) (not used here),
!         WM(2) = HRL1, the previous value of H*RL1, used if MITER = 3.
! IWM   = Integer work space containing pivot information, starting at
!         IWM(31), if MITER is 1, 2, 4, or 5.  IWM also contains band
!         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
! X     = The right-hand side vector on input, and the solution vector
!         on output, of length N.
! IERSL = Output flag.  IERSL = 0 if no trouble occurred.
!         IERSL = 1 if a singular matrix arose with MITER = 3.
!-----------------------------------------------------------------------

! Type declarations for labeled COMMON block SVOD01 --------------------

      REAL(PREC) ACNRM, CCMXJ, CONP, CRATE, DRC, EL,  &
           ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
           RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
              L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
              LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
              N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
              NSLP, NYH

! Type declarations for local variables --------------------------------

      INTEGER I, MEBAND, ML, MU
      REAL(PREC) DI, HRL1, ONE, PHRL1, R, ZERO
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE ONE, ZERO

      COMMON /SVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),  &
                      ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,  &
                      RC, RL1, TAU(13), TQ(5), TN, UROUND,  &
                      ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,  &
                      L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,  &
                      LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,  &
                      N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,  &
                      NSLP, NYH

      DATA ONE /1.0E0/, ZERO /0.0E0/

      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
 100  CALL SGESL (WM(3), N, N, IWM(31), X, 0)
      RETURN

 300  PHRL1 = WM(2)
      HRL1 = H*RL1
      WM(2) = HRL1
      IF (HRL1 .EQ. PHRL1) GO TO 330
      R = HRL1/PHRL1
      DO 320 I = 1,N
        DI = ONE - R*(ONE - ONE/WM(I+2))
        IF (ABS(DI) .EQ. ZERO) GO TO 390
 320    WM(I+2) = ONE/DI

 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN

 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL SGBSL (WM(3), MEBAND, N, ML, MU, IWM(31), X, 0)
      RETURN
!----------------------- End of Subroutine SVSOL -----------------------
      END SUBROUTINE SVSOL

      SUBROUTINE SVSRCO (RSAV, ISAV, JOB)
      USE messy_clams_global, ONLY: prec
      REAL(PREC) RSAV
      INTEGER ISAV, JOB
      DIMENSION RSAV(*), ISAV(*)
!-----------------------------------------------------------------------
! Call sequence input -- RSAV, ISAV, JOB
! Call sequence output -- RSAV, ISAV
! COMMON block variables accessed -- All of /SVOD01/ and /SVOD02/
!
! Subroutines/functions called by SVSRCO.. None
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of the
! COMMON blocks SVOD01 and SVOD02, which are used internally by SVODE.
!
! RSAV = real array of length 49 or more.
! ISAV = integer array of length 41 or more.
! JOB  = flag indicating to save or restore the COMMON blocks..
!        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV).
!        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV).
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
      REAL(PREC) RVOD1, RVOD2
      INTEGER IVOD1, IVOD2
      INTEGER I, LENIV1, LENIV2, LENRV1, LENRV2
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE LENRV1, LENIV1, LENRV2, LENIV2

      COMMON /SVOD01/ RVOD1(48), IVOD1(33)
      COMMON /SVOD02/ RVOD2(1), IVOD2(8)
      DATA LENRV1/48/, LENIV1/33/, LENRV2/1/, LENIV2/8/

      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRV1
 10     RSAV(I) = RVOD1(I)
      DO 15 I = 1,LENRV2
 15     RSAV(LENRV1+I) = RVOD2(I)

      DO 20 I = 1,LENIV1
 20     ISAV(I) = IVOD1(I)
      DO 25 I = 1,LENIV2
 25     ISAV(LENIV1+I) = IVOD2(I)

      RETURN

 100  CONTINUE
      DO 110 I = 1,LENRV1
 110     RVOD1(I) = RSAV(I)
      DO 115 I = 1,LENRV2
 115     RVOD2(I) = RSAV(LENRV1+I)

      DO 120 I = 1,LENIV1
 120     IVOD1(I) = ISAV(I)
      DO 125 I = 1,LENIV2
 125     IVOD2(I) = ISAV(LENIV1+I)

      RETURN
!----------------------- End of Subroutine SVSRCO ----------------------
      END SUBROUTINE SVSRCO

      SUBROUTINE SEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
      USE messy_clams_global, ONLY: prec
      REAL(PREC) RTOL, ATOL, YCUR, EWT
      INTEGER N, ITOL
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
!-----------------------------------------------------------------------
! Call sequence input -- N, ITOL, RTOL, ATOL, YCUR
! Call sequence output -- EWT
! COMMON block variables accessed -- None
!
! Subroutines/functions called by SEWSET.. None
!-----------------------------------------------------------------------
! This subroutine sets the error weight vector EWT according to
!     EWT(i) = RTOL(i)*abs(YCUR(i)) + ATOL(i),  i = 1,...,N,
! with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
! depending on the value of ITOL.
!-----------------------------------------------------------------------
      INTEGER I

      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1, N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 20   CONTINUE
      DO 25 I = 1, N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      RETURN
 30   CONTINUE
      DO 35 I = 1, N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 40   CONTINUE
      DO 45 I = 1, N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      RETURN
!----------------------- End of Subroutine SEWSET ----------------------
      END SUBROUTINE SEWSET

      REAL(PREC) FUNCTION SVNORM (N, V, W)
      USE messy_clams_global, ONLY: prec
      REAL(PREC) V, W
      INTEGER N
      DIMENSION V(N), W(N)
!-----------------------------------------------------------------------
! Call sequence input -- N, V, W
! Call sequence output -- None
! COMMON block variables accessed -- None
!
! Subroutines/functions called by SVNORM.. None
!-----------------------------------------------------------------------
! This function routine computes the weighted root-mean-square norm
! of the vector of length N contained in the array V, with weights
! contained in the array W of length N..
!   SVNORM = sqrt( (1/N) * sum( V(i)*W(i) )**2 )
!-----------------------------------------------------------------------
      REAL(PREC) SUM
      INTEGER I

      SUM = 0.0E0
      DO 10 I = 1, N
 10     SUM = SUM + (V(I)*W(I))**2
      SVNORM = SQRT(SUM/REAL(N))
      RETURN
!----------------------- End of Function SVNORM ------------------------
      END FUNCTION SVNORM

      REAL(PREC) FUNCTION R1MACH (IDUM)
      USE messy_clams_global, ONLY: prec
      INTEGER IDUM
!-----------------------------------------------------------------------
! This routine computes the unit roundoff of the machine.
! This is defined as the smallest positive machine number
! u such that  1.0 + u .ne. 1.0
!
! Subroutines/functions called by R1MACH.. None
!-----------------------------------------------------------------------
      REAL(PREC) U, COMP
      U = 1.0E0
 10   U = U*0.5E0
      COMP = 1.0E0 + U
      IF (COMP .NE. 1.0E0) GO TO 10
      R1MACH = U*2.0E0
      RETURN
!----------------------- End of Function R1MACH ------------------------
      END FUNCTION R1MACH

      SUBROUTINE XERRWV (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
      USE messy_clams_global, ONLY: prec
      REAL(PREC) R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
! op_pj_20160830+
!!$   CHARACTER*1 MSG(NMES)
      CHARACTER(LEN=*) MSG
! op_pj_20160830-
!-----------------------------------------------------------------------
! Subroutines XERRWV, XSETF, XSETUN, and the two function routines
! MFLGSV and LUNSAV, as given here, constitute a simplified version of
! the SLATEC error handling package.
! Written by A. C. Hindmarsh and P. N. Brown at LLNL.
! Version of 13 April, 1989.
! This version is in single precision.
!
! All arguments are input arguments.
!
! MSG    = The message (character array).
! NMES   = The length of MSG (number of characters).
! NERR   = The error number (not used).
! LEVEL  = The error level..
!          0 or 1 means recoverable (control returns to caller).
!          2 means fatal (run is aborted--see note below).
! NI     = Number of integers (0, 1, or 2) to be printed with message.
! I1,I2  = Integers to be printed, depending on NI.
! NR     = Number of reals (0, 1, or 2) to be printed with message.
! R1,R2  = Reals to be printed, depending on NR.
!
! Note..  this routine is machine-dependent and specialized for use
! in limited context, in the following ways..
! 1. The argument MSG is assumed to be of type CHARACTER, and
!    the message is printed with a format of (1X,80A1).
! 2. The message is assumed to take only one line.
!    Multi-line messages are generated by repeated calls.
! 3. If LEVEL = 2, control passes to the statement   STOP
!    to abort the run.  This statement may be machine-dependent.
! 4. R1 and R2 are assumed to be in single precision and are printed
!    in E21.13 format.
!
! For a different default logical unit number, change the data
! statement in function routine LUNSAV.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERRWV.. None
! Function routines called by XERRWV.. MFLGSV, LUNSAV
!-----------------------------------------------------------------------

! ju_nt_20160617
!!$      INTEGER I, LUNIT, LUNSAV, MESFLG, MFLGSV
      INTEGER I, LUNIT, MESFLG

! Get message print flag and logical unit number. ----------------------
      MESFLG = MFLGSV (0,.FALSE.)
      LUNIT = LUNSAV (0,.FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
! Write the message. ---------------------------------------------------
! op_pj_20160830+
!!$   WRITE (LUNIT,10) (MSG(I),I=1,NMES)
      WRITE (LUNIT,10) (MSG(I:I),I=1,NMES)
! op_pj_20160830-
 10   FORMAT(1X,80A1)
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',E21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',E21.13,3X,'R2 =',E21.13)
! Abort the run if LEVEL = 2. ------------------------------------------
 100  IF (LEVEL .NE. 2) RETURN
      STOP
!----------------------- End of Subroutine XERRWV ----------------------
      END SUBROUTINE XERRWV

      SUBROUTINE XSETF (MFLAG)
!-----------------------------------------------------------------------
! This routine resets the print control flag MFLAG.
!
! Subroutines called by XSETF.. None
! Function routines called by XSETF.. MFLGSV
!-----------------------------------------------------------------------
! ju_nt_20160617
!!$      INTEGER MFLAG, JUNK, MFLGSV
     INTEGER MFLAG, JUNK

      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = MFLGSV (MFLAG,.TRUE.)
      RETURN
!----------------------- End of Subroutine XSETF -----------------------
      END SUBROUTINE XSETF

      SUBROUTINE XSETUN (LUN)
!-----------------------------------------------------------------------
! This routine resets the logical unit number for messages.
!
! Subroutines called by XSETUN.. None
! Function routines called by XSETUN.. LUNSAV
!-----------------------------------------------------------------------
! ju_nt_20160617
!!$   INTEGER LUN, JUNK, LUNSAV
      INTEGER LUN, JUNK

      IF (LUN .GT. 0) JUNK = LUNSAV (LUN,.TRUE.)
      RETURN
!----------------------- End of Subroutine XSETUN ----------------------
      END SUBROUTINE XSETUN

      INTEGER FUNCTION MFLGSV (IVALUE, ISET)
      LOGICAL ISET
      INTEGER IVALUE
!-----------------------------------------------------------------------
! MFLGSV saves and recalls the parameter MESFLG which controls the
! printing of the error messages.
!
! Saved local variable..
!
!   MESFLG = Print control flag..
!            1 means print all messages (the default).
!            0 means no printing.
!
! On input..
!
!   IVALUE = The value to be set for the MESFLG parameter,
!            if ISET is .TRUE. .
!
!   ISET   = Logical flag to indicate whether to read or write.
!            If ISET=.TRUE., the MESFLG parameter will be given
!            the value IVALUE.  If ISET=.FALSE., the MESFLG
!            parameter will be unchanged, and IVALUE is a dummy
!            parameter.
!
! On return..
!
!   The (old) value of the MESFLG parameter will be returned
!   in the function value, MFLGSV.
!
! This is a modification of the SLATEC library routine J4SAVE.
!
! Subroutines/functions called by MFLGSV.. None
!-----------------------------------------------------------------------
      INTEGER MESFLG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE MESFLG
      DATA MESFLG/1/

      MFLGSV = MESFLG
      IF (ISET) MESFLG = IVALUE
      RETURN
!----------------------- End of Function MFLGSV ------------------------
      END FUNCTION MFLGSV

      INTEGER FUNCTION LUNSAV (IVALUE, ISET)
      LOGICAL ISET
      INTEGER IVALUE
!-----------------------------------------------------------------------
! LUNSAV saves and recalls the parameter LUNIT which is the logical
! unit number to which error messages are printed.
!
! Saved local variable..
!
!  LUNIT   = Logical unit number for messages.
!            The default is 6 (machine-dependent).
!
! On input..
!
!   IVALUE = The value to be set for the LUNIT parameter,
!            if ISET is .TRUE. .
!
!   ISET   = Logical flag to indicate whether to read or write.
!            If ISET=.TRUE., the LUNIT parameter will be given
!            the value IVALUE.  If ISET=.FALSE., the LUNIT
!            parameter will be unchanged, and IVALUE is a dummy
!            parameter.
!
! On return..
!
!   The (old) value of the LUNIT parameter will be returned
!   in the function value, LUNSAV.
!
! This is a modification of the SLATEC library routine J4SAVE.
!
! Subroutines/functions called by LUNSAV.. None
!-----------------------------------------------------------------------
      INTEGER LUNIT
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE LUNIT
      DATA LUNIT/6/

      LUNSAV = LUNIT
      IF (ISET) LUNIT = IVALUE
      RETURN
!----------------------- End of Function LUNSAV ------------------------
      END FUNCTION LUNSAV

    End Module messy_clamschem_asad_svode

    
