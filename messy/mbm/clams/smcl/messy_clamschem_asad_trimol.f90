Module messy_clamschem_asad_trimol

contains

!
!     Paul D.Brown & Glenn Carver, Centre for Atmospheric Science,
!                                  University of Cambridge.
!
!
! Purpose: Calculates trimolecular rate coefficients
!
!          Called from ASAD_CDRIVE
!
!     Method
!     ------
!     See the IUPAC reference material on their website for details on
!     calculation of termolecular rates.
!     http://www.iupac-kinetic.ch.cam.ac.uk/
!
!     Local variables
!     ---------------
!     zo         Low pressure limit to rate*density
!     zi         High pressure limit to rate
!     zr         Ratio of zo/zi
!     iho2       Reaction index for HO2+HO2+M
!     ih2o       Array index for advected tracer H2O
!     in2o5      Reaction index for N2O5+M  
!     ino2no3    Reaction index for NO2+NO3+M
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_TRIMOL(n_points)

        USE messy_clamschem_asad_mod, ONLY: rk, at, ntrkx, spt, t300, t, peps,   &
                                            tnd, f, wp, advt
        USE messy_clams_global,             ONLY: prec
        USE messy_clamschem_asad_mod_clams, ONLY: jpctr, jptk

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n_points

!       Local variables

        INTEGER, SAVE :: ih2o              ! Index for h2o in tracer array
        INTEGER       :: iho2              ! Index for ho2+ho2 in rk array
        INTEGER       :: in2o5             ! Reaction index for N2O5+M  
        INTEGER       :: ino2no3           ! Reaction index for NO2+NO3+M
! ju_nt_20140212
!   iohco:   moved here from  bimol (jug, 12/2013)    
        INTEGER       :: iohco             ! Reaction index for OH + CO.
        INTEGER       :: j                 ! Loop variable
        INTEGER       :: jl                ! Loop variable
        INTEGER       :: jtr               ! Loop variable
        INTEGER       :: jr                ! Index

        REAL(PREC) :: zo                         ! k_0
        REAL(PREC) :: zi                         ! k_infinity
        REAL(PREC) :: zfc                        ! F_c
        REAL(PREC) :: zr                         ! k_0/k_infinity
        REAL(PREC), ALLOCATABLE, SAVE :: nf(:)   ! component of broadening factor

        LOGICAL, SAVE :: first = .true.


!       1.  Calculate trimolecular rate coefficients
!           --------- ------------ ---- ------------

        iho2  = 0
        in2o5   = 0  
        ino2no3 = 0
! ju_nt_20140312: moved here from  bimol (jug, 12/2013)    
        iohco = 0

        IF (first) THEN
! Calculate the N factor to be used in broadening factor
          ALLOCATE(nf(jptk+1))
          nf(:) = 1.0
! ju_nt_20140312: broadening factor not used for JPL
!!$          WHERE (at(:,1) > 1E-3) nf(:) = 0.75 - 1.27*ALOG10(at(:,1))
          
! Check if H2O is an advected tracer
          ih2o = 0
          DO jtr = 1, jpctr
            IF ( advt(jtr)  ==  'H2O       ' ) ih2o = jtr
          END DO
          first = .false.
        END IF

        DO j = 1, jptk
          jr = ntrkx(j)

          IF ( spt(j,1) == 'HO2    '.AND. spt(j,2) == 'HO2    ' )       &
               iho2 = jr

! ju_nt_20140212
!      moved here from  bimol (jug, 12/2013)    
          if ( ( spt(j,1) == 'OH     ' .and. spt(j,2) == 'CO     ' ) .or. &
               ( spt(j,1) == 'CO     ' .and. spt(j,2) == 'OH     ' ) ) &
               iohco = jr

!         N2O5 + M  
!         reset to zero each step  

          in2o5 = 0  
          IF (spt(j,1) == 'N2O5      ') in2o5 = jr  
  
!         NO2 + NO3 + M  
!         reset to zero each step  

          ino2no3 = 0  
          IF ((spt(j,1) == 'NO2       ' .AND. spt(j,2) == 'NO3       ') &   
               .OR.(spt(j,2) == 'NO2       ' .AND.                      &  
               spt(j,1) == 'NO3       ')) ino2no3 = jr  

          DO jl = 1, n_points
            zo = at(j,2) * t300(jl)**at(j,3) *                          &
                           EXP( -at(j,4)/t(jl) ) * tnd(jl)
            zi = at(j,5) * t300(jl)**at(j,6) * EXP( -at(j,7)/t(jl) )
            IF ( zo < peps ) THEN
              rk(jl,jr) = zi
            ELSE IF ( zi < peps ) THEN
              rk(jl,jr) = zo 
            ELSE
              IF( at(j,1) <= 1.0 ) THEN
                zfc = at(j,1)
              ELSE IF ( in2o5 /= 0 .OR. ino2no3 /= 0 ) THEN ! dependent 
                                                            ! reactions  
                zfc = 2.5*EXP(-1950.0/t(jl))+0.9*EXP(-t(jl)/at(j,1))
              ELSE              ! temperature dependent Fc
                zfc = EXP( -t(jl)/at(j,1) )
! ju_nt_20140312: broadening factor not used for JPL
!!$                nf(j) = 0.75 - 1.27*ALOG10(zfc)
              END IF
              zr = zo / zi
! ju_nt_20140402: use LOG10 instead of ALOG10:
!!$              rk(jl,jr) = (zo/(1.0+zr)) *                               &
!!$                           zfc**(1.0/(1.0 + (ALOG10(zr)/nf(j))**2))
              rk(jl,jr) = (zo/(1.0+zr)) *                               &
                           zfc**(1.0/(1.0 + (LOG10(zr)/nf(j))**2))
            END IF
          END DO
        END DO              ! end of loop over jptk

!       2. Dependent reactions.
!          --------- ----------

! HO2 + HO2 [+ M]
        IF (ih2o /= 0 .AND. iho2 /= 0 ) THEN
!         h2o is an advected tracer
          DO jl = 1, n_points
            rk(jl,iho2) = rk(jl,iho2) *                                &
            ( 1.0 + 1.4E-21*f(jl,ih2o)*EXP(2200./t(jl)) )
          END DO
        ELSE IF (ih2o == 0 .AND. iho2 /= 0) THEN
!         use modelled water concentration
          DO jl = 1, n_points
            rk(jl,iho2) = rk(jl,iho2) *                                &
            ( 1.0 + 1.4E-21*wp(jl)*tnd(jl)*EXP(2200./t(jl)) )
          END DO
        END IF

! ju_nt_20140212: following lines added
! OH + CO [+ M]
        if ( iohco /= 0) then
!          assume to have reaction with product HOCO in from ratt.d
!          add reaction with product H + CO2 here  (JPL 2011)
!          needs to be changed for reaction rate update 
!          t300 == T/300
!          note that z0 is defined differently (w/o *tnd) [jug, 12/2013]
           DO jl = 1, n_points
              zo = 1.5E-13 * t300(jl)**0.6
              zi = 2.1E+09 * t300(jl)**6.1
              zr = zo * tnd(jl) / zi
! ju_nt_20140402: use LOG10 instead of ALOG10:
!!$              rk(jl,iohco) = rk(jl,iohco) + &
!!$                   zo/(1.0+zr) *  0.6**(1.0 / (1.0+((alog10(zr))**2)))
              rk(jl,iohco) = rk(jl,iohco) + &
                   zo/(1.0+zr) *  0.6**(1.0 / (1.0+((log10(zr))**2)))
           ENDDO
        endif

        RETURN
        END SUBROUTINE ASAD_TRIMOL

      End Module messy_clamschem_asad_trimol
