MODULE messy_gwave_rayleighbg

  USE messy_main_constants_mem, ONLY: prcn=>dp, pi

  IMPLICIT NONE

  PRIVATE

  !v Rayleigh drag variables
  !REAL(prcn) :: max_rayleigh = 0.6
  !REAL(prcn) :: ht1          = 25.e3
  !REAL(prcn) :: ht2          = 80.e3

  ! Rayleigh params. 1- Rayleigh everywhere, 2-just in winter
  INTEGER      :: add_rayleigh = 2
  REAL(prcn) :: max_rayleigh = 0.2
  REAL(prcn) :: ht1          = 25.e3
  REAL(prcn) :: ht2          = 65.e3

  PUBLIC :: rayleighbackground

CONTAINS

  SUBROUTINE rayleighbackground(vy1d,h,edrag, dayno,  ht_dim, lat)

    !r Let's add a small Rayleigh drag background

    !r Inputs
    !r --------------------------------------------------------------
  !r h(ht_dim)       : 1-d array of Height (m)
    !r vy1d         : 1d array of east-west velocity (m*s-1)
    !r lat          : latitude 
    !r ht_dim       : size of vertical 1d arrays
    !r
    !r Outputs
    !r -------------------------------------------------------------
    !r
    !r edrag     : 1d array of drag (m*s-1*day-1))

    IMPLICIT NONE ! op_pj_20110202

    INTEGER,      INTENT(IN)    :: ht_dim
    REAL(prcn),   INTENT(IN)    :: vy1d(ht_dim)
    REAL(prcn),   INTENT(OUT)   :: edrag(ht_dim)
    INTEGER,      INTENT(IN)    :: dayno
    REAL(prcn),   INTENT(IN)    :: lat
    REAL(prcn),   INTENT(IN)    :: h(ht_dim)     

    ! Solstice factor, peaks at solstices
    REAL(prcn) :: solfac
    ! Which hemisphere winter solstice is
    INTEGER    :: solhem = 1

    REAL(prcn) :: rayleigh

    INTEGER :: n

    ! solhem
    solhem=2
    IF(dayno > 80 .and. dayno < 270) solhem = 1

    ! add_rayleigh==1 everywhere, 2-Only apply to winter hemisphere
    IF((solhem == 1 .AND. lat <= 0. ) .OR.   &
         (solhem == 2 .AND. lat > 0.  ) .OR.   &
         add_rayleigh == 1)  THEN

       IF(add_rayleigh == 1) THEN
          solfac=1.
       ELSE
          ! Solfac COS peak at solstices
          solfac = COS((dayno/360. + 10/360.)*4*PI)
          IF(solfac < 0) solfac = 0.
       ENDIF

       DO n = 1, ht_dim

          ! Zero
          IF(h(n) .LT. ht1) rayleigh = 0.

          ! Increase
          !IF(h(n) .GE. ht1 .AND. h(n) .LE. ht2)                          &
          !     rayleigh = max_rayleigh*(EXP(-(ht2-h(n))*(1.e-3/10.)))

          IF(h(n) .GE. ht1 .AND. h(n) .LE. ht2)                           &
               rayleigh = max_rayleigh*EXP((h(n)-ht2)/27.e3)

          ! Decrease
          IF(h(n) .GT. ht2) rayleigh = rayleigh*EXP((ht2-h(n))*1.e-3/20.)

          edrag(n) = - solfac*1.*(rayleigh*vy1d(n)/86400.)

       ENDDO

    ENDIF

  END SUBROUTINE rayleighbackground

END MODULE messy_gwave_rayleighbg
