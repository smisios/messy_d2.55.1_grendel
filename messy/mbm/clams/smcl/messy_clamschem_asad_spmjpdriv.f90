Module messy_clamschem_asad_spmjpdriv

contains

! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer. 
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution. 
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Driver for fully implicit ODE integrator.
!    Part of the ASAD chemical solver.
!
!   Called from ASAD_CDRIVE
!
!
!     MJPDRIV  - Driver for MJP fully implicit ODE integrator.
!
!     Michael Prather            Earth System Science
!     Oliver Wild                University of California, Irvine
!
!     ASAD: mjpdriv              Version: mjpdriv.f 1.0 04/17/97
!
!     Modification (02/2015) by J.-U.Grooss, Forschungszentrum Juelich
!     Check each cell individually for convergence. If convergence is
!     achieved, that is converged(jl)==.true. then the comical composition
!     is noch changed. This avoids the results being dependent on
!     the distribution onto compute cores, especilly reset to init values
!     for all if one of the cells in a compute core has not converged. 
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the MJP implicit integrator.
!
!     Interface
!     ---------
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
!     This routine calls the MJP integrator once for each gridpoint
!     of the one-dimensional arrays passed to it, setting 'jlst' to
!     the current level in the same way as the SVODE integrator.
!     If convergence isn't achieved, the time step is halved, and the
!     integrator called again - this is continued until either
!     convergence is achieved or the minimum time step length is
!     encountered (currently 1.E-05 seconds).
!
!     Local variables
!     ---------------
!     ncst    -  Stores number of basic chemical steps 'ncsteps'
!     ctrd    -  Stores basic chemical time step length 'ctd'
!     zf      -  Stores family concentrations at beginning of call
!     ndxraf  -  Error code from the integrator:
!                   0 = successful return
!                   1 = negatives encountered
!                   2 = convergence failure after 'nrsteps' iterations
!                   3 = convergence failure due to divergence - 'NaN's
!                   4 = convergence failure (as '2') but set debugging
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE asad_spmjpdriv(nlev,n_points)

      USE messy_clamschem_asad_mod
      USE messy_clamschem_asad_dummy, ONLY: ereport
      USE messy_clams_global, ONLY: prec
      USE messy_clamschem_asad_mod_clams,  only: theta_field_size, mype, &
                                 printstatus, prstatus_oper
      USE messy_clamschem_asad_diffun,   ONLY: ASAD_DIFFUN
      USE messy_clamschem_asad_ftoy,     ONLY: ASAD_FTOY
      USE messy_clamschem_asad_spimpmjp, ONLY: ASAD_SPIMPMJP


      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points
      INTEGER, INTENT(IN) :: nlev

! Local variables
      INTEGER, PARAMETER :: max_redo=128      ! Max times for halving TS, was 16
      LOGICAL :: ltrig
      INTEGER :: ndxraf
      INTEGER :: ncst
      INTEGER :: iredo
      INTEGER :: jl
      INTEGER :: jtr
      INTEGER :: i
      INTEGER :: nl

      INTEGER :: errcode                ! Variable passed to ereport

      CHARACTER(LEN=72) :: cmessage

      REAL(PREC) :: ctrd
      REAL(PREC) :: zf(theta_field_size,jpctr)


      COMMON /trig/ltrig
!

! initialize logical array converged (jug, 02/2015)
      converged(:) = .false.

      ncst = ncsteps
      ctrd = cdt
      ltrig=.FALSE.
!
      nl = n_points
! DEPENDS ON: asad_diffun
      CALL asad_diffun( nl )
!
      iredo = 1
      zf(1:n_points,:)=f(1:n_points,:)
!
! Start iterations here.
      i = 1
      DO WHILE (i <= iredo)
! DEPENDS ON: asad_spimpmjp
        CALL asad_spimpmjp(ndxraf, nlev, n_points)

!  Debug slow convergence systems - switch this on in 'spimpmjp'
        IF (ndxraf == 4) THEN
          IF (ltrig) THEN
            errcode=1
            cmessage='Slow-converging system, '//                       &
              'Set printstatus for Jacobian debug'
!            DO jtr=1,jpspec
!              write(6,'(a4,i6,a12,2e14.5,i12)') 'y: ',jtr,speci(jtr),     &
!                  maxval(y(:,jtr)), minval(y(:,jtr)),size(y(:,jtr))
!            END DO
            CALL ereport('ASAD_SPMJPDRIV',errcode,cmessage)
          END IF
          ltrig=.TRUE.
          do jl=1, n_points
             if (.not. converged(jl)) f(jl,:)=zf(jl,:)
          end do
! DEPENDS ON: asad_ftoy
          CALL asad_ftoy( .FALSE., nitfg, n_points )
! DEPENDS ON: asad_diffun
          CALL asad_diffun( nl )
          i = 1
        ELSE
!
!  Reset for failed convergence
          IF (ndxraf > 1) THEN
            ncsteps = ncsteps*2
            cdt = cdt/2.
            iredo = iredo*2
            IF(cdt < 1.0e-05) THEN
              errcode=2
              cmessage=' Time step now too short'
              CALL ereport('ASAD_SPMJPDRIV',errcode,cmessage)
            END IF
            do jl=1, n_points
               if (.not. converged(jl)) f(jl,:)=zf(jl,:)
            end do
! Drop out at some point - if 3 successive halvings fail
            IF (iredo >= max_redo) THEN
              IF (printstatus >= prstatus_oper) THEN  
                WRITE(6,"(' Resetting array after',i4,' iterations')")  &
                    iredo
                WRITE(6,"('NO CONVERGENCE nlev: ',i4,' pe: ',i4)")      &
                    nlev,mype
              ENDIF
              ncsteps = ncst
              cdt = ctrd
              GOTO 9999
            END IF
! DEPENDS ON: asad_ftoy
            CALL asad_ftoy( .FALSE., nitfg, n_points )
! DEPENDS ON: asad_diffun
            CALL asad_diffun( nl )
            i = 1
          ELSE
            i = i + 1
          END IF
        END IF
      END DO

      IF(iredo > 1) THEN
        IF (iredo > 2) WRITE(6,"('   No. iterations =',i2)") iredo
        ncsteps = ncst
        cdt = ctrd
      END IF

 9999 CONTINUE
      RETURN
      END SUBROUTINE asad_spmjpdriv

    End Module messy_clamschem_asad_spmjpdriv


