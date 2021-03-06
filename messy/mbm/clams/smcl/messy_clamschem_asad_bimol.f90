Module messy_clamschem_asad_bimol

contains

!**** asad_bimol - calculates bimolecular rate coefficients
!
!     Paul D.Brown & Glenn Carver, Centre for Atmospheric Science,
!                                  University of Cambridge.
!
!
! Purpose: To calculate bimolecular rate coefficients using 
!          data from clams_chem_data module
!
!
!          Called from ASAD_CDRIVE 
!
!     Method
!     ------
!     Bimolecular rate coefficients are calculated from the Ahrenius
!     expression: k(T) = k * (T/300)^a * exp(-b/T) where T is the
!     temperature (in K). The parameters k, a and b are taken from
!     the bimolecular ratefile ratb.d.
!
!     The reactions CO + OH -> H + CO2 and OH + HONO2 are pressure
!     or density dependent and therefore need to be calculated
!     separately. However, this code should never need changing even
!     if neither of these reaction are included in the chemistry.
!
!     The reactions HO2+MeCO3, MeOO+MeOO and OH+C3H8 have temperature
!     dependent branching ratios. Therefore, they are calculated
!     separately.
!
!     Local Variables
!     ---------------
!     ih2o           Tracer index for H2O if it's an advective tracer
!     iohco          Reaction index for OH + CO.
!     iohhno3        Reaction index for OH + HONO2.
!     iso3h2o        Reaction index for SO3 + H2O -> H2SO4 + H2O
!     ics2oh         Reaction index for CS2 + OH -> COS + SO2
!     z1,z3,z4       Variables used to calculate pressure
!                    dependent rate coefficients.
!     ratioa2b       Branching ratio of branch A to branch B
!
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      SUBROUTINE ASAD_BIMOL( n_points )

      USE messy_clamschem_asad_mod, ONLY: t, t300, advt, spb, ab, rk, tnd, p,    &
                                          wp, f, peps, nbrkx, jpspb
      USE messy_clams_global,                ONLY: prec
      USE messy_clamschem_asad_mod_clams,    ONLY: jpctr, jpbk
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n_points

! Local variables

      INTEGER, SAVE :: iohco = 0
      INTEGER, SAVE :: iohhno3 = 0
      INTEGER, SAVE :: iho2 = 0
      INTEGER, SAVE :: iso3h2o = 0
      INTEGER, SAVE :: ics2oh = 0
      INTEGER, SAVE :: ih2o = 0
      INTEGER, SAVE :: io2  = 0
      INTEGER, SAVE :: iho2no = 0
      INTEGER :: jtr
      INTEGER :: j
      INTEGER :: jr

      REAL(PREC), ALLOCATABLE :: z1(:)                ! Intermediate result
      REAL(PREC), ALLOCATABLE :: z3(:)                !          "
      REAL(PREC), ALLOCATABLE :: z4(:)                !          "
      REAL(PREC), ALLOCATABLE :: z5(:)                !          "
      REAL(PREC), ALLOCATABLE :: z6(:)                !          "
      REAL(PREC), ALLOCATABLE :: alpha(:)             ! Multiplication factor
      REAL(PREC), ALLOCATABLE :: ratioa2b(:)          ! Branching ratio
      REAL(PREC), ALLOCATABLE :: ratiob2total(:)      ! Branching ratio


      LOGICAL, SAVE :: first = .TRUE.


      ALLOCATE(z1(1:n_points))
      ALLOCATE(z3(1:n_points))
      ALLOCATE(z4(1:n_points))
      ALLOCATE(z5(1:n_points))
      ALLOCATE(z6(1:n_points))
      ALLOCATE(ratioa2b(1:n_points))
      ALLOCATE(ratiob2total(1:n_points))

      z1(:) = 0.0
      z3(:) = 0.0
      z4(:) = 0.0
      z5(:) = 0.0
      z6(:) = 0.0
      ratioa2b(:) = 0.0
      ratiob2total(:) = 0.0

!       1. Calculate bimolecular rate coefficients
!          --------- ----------- ---- ------------

!       Compute intermediate results

      t300(1:n_points) = t(1:n_points) / 300.0

!       Check if H2O is an advected tracer

      IF (first) then
        first = .FALSE.
        DO jtr = 1, jpctr
          IF ( advt(jtr)  ==  'H2O       ' ) ih2o = jtr
          IF ( advt(jtr)  ==  'O2        ' ) io2 = jtr
        ENDDO

!       Look for the reactions which need special treatment.

      DO j = 1, jpbk
        jr = nbrkx(j)


!        IF ( ( trim(spb(j,1)) == 'HO2' .AND. trim(spb(j,2)) == 'NO' ).OR.  &
!           (  trim(spb(j,1)) == 'NO' .AND. trim(spb(j,2)) == 'HO2' ) )     &
!          iho2no = jr

        
! ju_nt_20140212
! OH + CO treated as termoclecular reaction -- moved to trimol (jug, 12/2013)
!!$        IF ( ( spb(j,1) == 'OH     '.AND. spb(j,2) == 'CO     ' ).OR.  &
!!$           ( spb(j,1) == 'CO     '.AND. spb(j,2) == 'OH     ' ) )      &
!!$          iohco = jr
        IF ( ( spb(j,1) == 'OH     '.AND. spb(j,2) == 'HONO2  ' ).OR.  &
           (  spb(j,1) == 'HONO2  '.AND. spb(j,2) == 'OH     ' ) )     &
          iohhno3 = jr
        IF ( spb(j,1) == 'HO2    '.AND. spb(j,2) == 'HO2    ' )        &
          iho2 = jr

! code for stratospheric sulphur scheme
          IF ( spb(j,1) == 'SO3      '.AND. spb(j,2) == 'H2O   ')    &   
            iso3h2o = jr
          IF ( spb(j,1) == 'CS2      '.AND. spb(j,2) == 'OH    ')    &   
            ics2oh = jr

      END DO  ! end of loop (j) over jpbk

      END IF   ! first

!       1.2  Compute rates

      DO j = 1, jpbk
        jr = nbrkx(j)

        IF ( ABS(ab(j,2)) < peps .AND. ABS(ab(j,3)) < peps ) THEN
          rk(1:n_points,jr) = ab(j,1)
        ELSE IF ( ABS(ab(j,2)) < peps ) THEN
          rk(1:n_points,jr) = ab(j,1) * exp(-ab(j,3)/t(1:n_points))
        ELSE IF ( ABS(ab(j,3)) < peps ) THEN
          rk(1:n_points,jr) = ab(j,1) * t300(1:n_points)**ab(j,2)
        ELSE
          rk(1:n_points,jr) = ab(j,1) * t300(1:n_points)**ab(j,2) *     &
                                exp(-ab(j,3)/t(1:n_points))
        END IF
      END DO  ! end of loop (j) over jpbk

!       2. Dependent reactions.
!          --------- ----------

! OH + CO; updated with IUPAC March 2005 (Paul Young)
! GC: note the original March IUPAC summary had a mistake for this
! reaction. The values given in the datasheet are correct.
! k = k' . (1 + [N2]/4.2E19).. we use TND below instead of [N2] but
! Paul suggests the reaction would probably go with [O2] anyway.

! ju_nt_20140212
!    OH + CO treated as termoclecular reaction -- moved to trimol (jug, 12/2013)
!!$      IF ( iohco /= 0 ) THEN
!!$        rk(1:n_points,iohco)=rk(1:n_points,iohco)*                      &
!!$                        (1.0 + tnd(1:n_points)/4.2e19)
!!$      END IF

! OH + HONO2; no change with IUPAC Jan 2009 (CJ)
      IF ( iohhno3 /= 0 ) THEN
        z1(:) = 2.4e-14 * EXP(460.0/t(1:n_points))
        z3(:) = (6.5e-34 * EXP(1335.0/t(1:n_points)))*tnd(1:n_points)
        z4(:) = 2.7e-17 * EXP(2199.0/t(1:n_points))
        rk(1:n_points,iohhno3) = z1(:) + z3(:)/(1.0+z3(:)/z4(:))
      END IF

! HO2 + HO2; no change with IUPAC Nov 2003 (Paul Young)
      IF ( ih2o /= 0 .AND. iho2 /= 0) THEN
!       water is an advected tracer
        rk(1:n_points,iho2) = rk(1:n_points,iho2) *                 &
          (1.0+1.4E-21*f(1:n_points,ih2o)*exp(2200.0/t(1:n_points)))
      ELSE IF (ih2o == 0 .AND. iho2 /= 0) then
!       use model water concentration
        rk(1:n_points,iho2) = rk(1:n_points,iho2) *                 &
          (1.0+1.4E-21*wp(1:n_points)*tnd(1:n_points)*              &
           EXP(2200./t(1:n_points)))
      END IF

! HO2 + NO -> HONO2 (with extra temp and pressure dependence)
! Added by Alex 2012  
      ! IF ( iho2no /= 0 ) THEN  
      !   rk(1:n_points,iho2no)=rk(1:n_points,iho2no)*                    &  
      !    ((530.0/t(1:n_points)) + 8.53E-4*(1E-2*p(1:n_points))-1.73)/100.0  
      ! END IF  

! SO3 + H2O: 2nd H2O molecule dealt with here by multiplying rate by [H2O]
      IF ( ih2o /= 0 .AND. iso3h2o /= 0 ) THEN
!       water is an advected tracer
        rk(1:n_points,iso3h2o) = rk(1:n_points,iso3h2o) *               &
             f(1:n_points,ih2o) 
      ELSE IF (ih2o == 0 .AND. iso3h2o /= 0) THEN
!       use model water concentration
        rk(1:n_points,iso3h2o) = rk(1:n_points,iso3h2o) *               &
             wp(1:n_points)*tnd(1:n_points)
      END IF


! IUPAC, 2011
!            rk(1:n_points,idmsoh) = rk(1:n_points,idmsoh)*             &
!                                    f(1:n_points,io2)/                 &
!              (1 + 7.5E-29*EXP(5610/t(1:n_points))*f(1:n_points,io2))

      
!       3. Temperature-Dependent branching ratios
!          ----------- --------- --------- ------
!       rk above was calculated using the total rate coefficients.
!       Here, rk is reduced according to the branching ratio.


! !  MeOO + MeOO -> MeOH + HCHO       ... Branch A
! !  MeOO + MeOO -> 2HO2 + 2HCHO      ... Branch B
!       IF (imeoomeooa /= 0 .AND. imeoomeoob /=0) THEN
!         ratiob2total(1:n_points) = 1.0/(1.0+(EXP(1300.0/                &
!                                      t(1:n_points)))/33.0)
!         rk(1:n_points,imeoomeooa) = rk(1:n_points,imeoomeooa)*          &
!                                      (1.0-ratiob2total(1:n_points))
!         rk(1:n_points,imeoomeoob) = rk(1:n_points,imeoomeoob)*          &
!                                      (ratiob2total(1:n_points))
!       END IF

! !  Added from IUPAC
! !  MeOO + HO2 -> MeOOH     ... Branch A
! !  MeOO + HO2 -> HCHO      ... Branch B

!       IF ( imeooho2a /= 0 .AND. imeooho2b /= 0 ) THEN
!         ratiob2total(1:n_points) = 1.0 / (1.0 + 498.0*EXP(-1160.0/      &
!                                      t(1:n_points)))
!         rk(1:n_points,imeooho2a) = rk(1:n_points,imeooho2a)*            &
!                                      (1.0-ratiob2total(1:n_points))
!         rk(1:n_points,imeooho2b) = rk(1:n_points,imeooho2b)*            &
!                                      (ratiob2total(1:n_points))
!       END IF


      IF (ALLOCATED(z1)) DEALLOCATE(z1)
      IF (ALLOCATED(z3)) DEALLOCATE(z3)
      IF (ALLOCATED(z4)) DEALLOCATE(z4)
      IF (ALLOCATED(z5)) DEALLOCATE(z5)
      IF (ALLOCATED(z6)) DEALLOCATE(z6)
      IF (ALLOCATED(ratioa2b)) DEALLOCATE(ratioa2b)
      IF (ALLOCATED(ratiob2total)) DEALLOCATE(ratiob2total)
      IF (ALLOCATED(alpha)) DEALLOCATE(alpha)

      RETURN
      END SUBROUTINE ASAD_BIMOL

end Module messy_clamschem_asad_bimol
