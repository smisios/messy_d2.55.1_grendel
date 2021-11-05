MODULE messy_m7_box

! DESCRIPTION
! Interface layer for M7 box model, analog messy_m7_e5.f90
!
! AUTHOR
! Swen Metzger, Max Planck Institute for Chemistry, Mainz, Germany
! questions/suggestions: metzger@mpch-mainz.mpg.de
!
! LAST CHANGES
! 6. October 2004, S. Metzger

!*****************************************************************************

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************

  USE messy_m7,                 ONLY: modstr,modver, nmod
  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: m7_initialize
  PUBLIC :: m7_physc
  PUBLIC :: m7_output

  INTEGER, PARAMETER, PUBLIC                  :: iounit=8,  &       ! I/O unit
                                                 irunit=9,  &       ! unit to read input from file
                                                 iwunit=10, &       ! unit to write ouput to file
                                                 nci =  9,  &       ! input  array dimensions
                                                 ndo = 14,  &       ! output array dimensions
                                                 ntime = 3600*24    ! integration time (to be multiplied by ndays)
  INTEGER, PUBLIC                             :: ndays,&            ! number of input data sets (days)
                                                 mtime,&            ! integration time ntime multiplied by ndays
                                                 kproma             ! vector lenght
  INTEGER, PUBLIC                             :: COUNT1,COUNT2,             &
                                                 COUNT_RATE,COUNT_MAX,      &
                                                 icount_hour,               &
                                                 ifirst_minute, ifirst_hour,&
                                                 ifirst_day, ifirst_month,  &
                                                 ifirst_year, ilast_year,   &
                                                 ilast_minute, ilast_hour,  &
                                                 ilast_day, ilast_month,    &
                                                 imin, ihour, iday, imon, iyear

  INTEGER, DIMENSION(8), PUBLIC               :: values             ! utility array
  INTEGER, DIMENSION(12), PUBLIC              :: icalender = (/31,28,31,30,31,30,31,31,30,31,30,31/) ! number of days per month

  INTEGER, PUBLIC, PARAMETER ::                                                 &
              NS=1,     KS=2,     AS=3,     CS=4,     KI=5,     AI=6,     CI=7
  REAL(dp),PUBLIC, PARAMETER                  :: zero = 0._dp              ! zero definition
  REAL(dp),PUBLIC, PARAMETER                  :: time_step_len = 1800._dp  ! integration time step
  REAL(dp), ALLOCATABLE, PUBLIC               :: yi(:,:,:)                 ! box model input  array
  REAL(dp), ALLOCATABLE, PUBLIC               :: yo(:,:,:,:)               ! box model output array
  REAL(dp), DIMENSION(nmod)                   :: factor
  CHARACTER(LEN=8),  PUBLIC                   :: date
  CHARACTER(LEN=10), PUBLIC                   :: time
  CHARACTER(LEN=12), PUBLIC                   :: zone
  CHARACTER(LEN=150),PUBLIC                   :: ctitle
  CHARACTER(LEN=6),  PUBLIC                   :: ctime
  CHARACTER(LEN=17), PUBLIC                   :: cdate
  CHARACTER(LEN=30), PUBLIC                   :: ifile
  CHARACTER(LEN=30), PUBLIC                   :: ofile
  CHARACTER(LEN=2),  PUBLIC                   :: cdummy,      &
                                                 cmin, chour, &
                                                 cday, cmon
  CHARACTER(LEN=4),  PUBLIC                   :: cyear

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE m7_initialize

    USE messy_m7,         ONLY: m7_read_nml_ctrl, m7_initialize_core
    USE messy_main_blather, ONLY: start_message, end_message

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER               :: substr = 'm7_initialize'

    INTEGER                                   :: status ! error status
    INTEGER                                   :: i,j,iou,jmod

    CALL start_message(TRIM(modstr),'read namelist and input data', substr)

    ! read CTRL namelist
    CALL m7_read_nml_ctrl(status, 99)
    IF (status /= 0) STOP

     DO jmod = 1, nmod

     WRITE(ifile,'(A15,I1,A4)') 'input/m7_input_', jmod,'.txt'
     WRITE(ofile,'(A16,I1,A4)') 'output/m7_input_',jmod,'.dat'
     WRITE(*,*) ' read  input data from file ',ifile,' for mode ',jmod
     WRITE(*,*) ' write input data to   file ',ofile

     OPEN(unit=irunit,file=ifile,status='old',form='formatted')
     OPEN(unit=iwunit,file=ofile,status='unknown',form='formatted')
     WRITE(*,*) ' '
     READ(irunit,'(A150)')  ctitle
     READ(irunit,*)  ndays, factor(jmod)
     IF (.NOT. ALLOCATED(yi)) ALLOCATE(yi(ndays+1, nci,nmod)); yi(:,:,jmod)=zero
     WRITE(*,*)      ctitle
     WRITE(*,*)      ndays, factor(jmod)
     WRITE(iwunit,*) '#',ctitle
     WRITE(iwunit,*) '#',ndays, factor(jmod)
     DO i=1,ndays
        READ (irunit,*) cdate,ctime,(yi(i,j,jmod),j=1,nci)
        WRITE(*,*)      cdate,ctime,(yi(i,j,jmod),j=1,nci)
        WRITE(iwunit,'(3A,19E20.9)') TRIM(cdate),'.',TRIM(ctime), (yi(i,j,jmod),j=1,nci)
        IF(i == 1) THEN
           OPEN(unit=iounit,file='output/first_data_date',status='unknown',form='formatted')
           WRITE (iounit,*) TRIM(cdate),'.',TRIM(ctime)
           REWIND(iounit)
           READ  (iounit,'(A1,2(I2,A1),I4,2(A1,I2))') &
                  cdummy,ifirst_day,  cdummy, ifirst_month, cdummy, ifirst_year, cdummy, ifirst_hour, cdummy, ifirst_minute
           CLOSE(iounit)
        END IF
     END DO

     CLOSE(irunit)
     CLOSE(iwunit)

     END DO ! jmod

     OPEN(unit=iounit,file='output/last_data_date',status='unknown',form='formatted')
     WRITE (iounit,*) TRIM(cdate),'.',TRIM(ctime)
     REWIND(iounit)
     READ  (iounit,'(A1,2(I2,A1),I4,2(A1,I2))') &
            cdummy,ilast_day,  cdummy, ilast_month, cdummy, ilast_year, cdummy, ilast_hour, cdummy, ilast_minute
     CLOSE(iounit)

     WRITE(*,*) ' '
     WRITE(*,'(2(A17,I4))') 'ifirst_year   = ', ifirst_year,  ' ilast_year   = ', ilast_year
     WRITE(*,'(2(A17,I4))') 'ifirst_month  = ', ifirst_month, ' ilast_month  = ', ilast_month
     WRITE(*,'(2(A17,I4))') 'ifirst_day    = ', ifirst_day,   ' ilast_day    = ', ilast_day
     WRITE(*,'(2(A17,I4))') 'ifirst_hour   = ', ifirst_hour,  ' ilast_hour   = ', ilast_hour
     WRITE(*,'(2(A17,I4))') 'ifirst_minute = ', ifirst_minute,' ilast_minute = ', ilast_minute

     WRITE(*,*) ' '
     WRITE(*,*) ' input data successfully read in'
     WRITE(*,*) ' '

     WRITE(*,*) ' initialized m7 core'
     CALL m7_initialize_core
     WRITE(*,*) ' '

     CALL end_message(TRIM(modstr),'read namelist and input data', substr)

  END SUBROUTINE m7_initialize

  ! --------------------------------------------------------------------------
  
  SUBROUTINE m7_physc

  ! Author
  ! Swen Metzger    (metzger@mpch-mainz.mpg.de), MPI-CHEM, June 2004

  USE messy_m7, ONLY: m7_main, nmod, naermod, nsol, l_m7 => lm7,&
           iso4ns, iso4ks,iso4as,iso4cs,        & !- Sulfate
                   ibcks, ibcas, ibccs, ibcki,  & !- Black Carbon
                   iocks, iocas, ioccs, iocki,  & !- Organic Carbon
                          issas, isscs,         & !- Sea Salt
                          iduas, iducs,        iduai, iduci, &  !- Dust
           inucs,  iaits, iaccs, icoas, iaiti, iacci, icoai     !- Number

  USE messy_main_blather, ONLY: start_message, end_message

  IMPLICIT NONE

  !--- Parameter list:

  ! Local variables:

  CHARACTER(LEN=*), PARAMETER :: substr = 'm7_physc'

  ! Set the dimension of the fields to 0:
  
  INTEGER                                 :: jmod, klev, jrow, jtime

  INTEGER, PARAMETER                      :: nrow  = 1, &
                                             nlev  = 1, &
                                             kbdim = 1

  REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: zgso4, zrhum,   &
                                             zdz,   ztemp,   &
                                             zpress
                                           
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: zaerml, zaernl, &
                                             zm6rp,  zm6dry, &
                                             zrhop,  zww

  !--- Initialisations: -----------------------------------------------

  ! (A) GRID SPACE
  jrow          = nrow
  klev          = nlev
  kproma        = kbdim
  mtime         = ntime*ndays/NINT(time_step_len)

  CALL start_message(modstr,'m7 aerosol dynamics',substr)

  ALLOCATE(zgso4    (kbdim ,klev))        ; zgso4     (:,:) = 0._dp
  ALLOCATE(zdz      (kbdim ,klev))        ; zdz       (:,:) = 0._dp
  ALLOCATE(zpress   (kbdim ,klev))        ; zpress    (:,:) = 0._dp
  ALLOCATE(ztemp    (kbdim ,klev))        ; ztemp     (:,:) = 0._dp
  ALLOCATE(zrhum    (kbdim ,klev))        ; zrhum     (:,:) = 0._dp
  ALLOCATE(zaerml   (kbdim ,klev,naermod)); zaerml  (:,:,:) = 0._dp
  ALLOCATE(zaernl   (kbdim ,klev,nmod))   ; zaernl  (:,:,:) = 0._dp
  ALLOCATE(zm6rp    (kbdim ,klev,nmod))   ; zm6rp   (:,:,:) = 0._dp
  ALLOCATE(zm6dry   (kbdim ,klev,nsol))   ; zm6dry  (:,:,:) = 0._dp
  ALLOCATE(zrhop    (kbdim ,klev,nmod))   ; zrhop   (:,:,:) = 0._dp
  ALLOCATE(zww      (kbdim ,klev,nmod))   ; zww     (:,:,:) = 0._dp
  ALLOCATE(yo       (kbdim ,ndo, nmod, mtime)) ; yo (:,:,:,:)=0._dp

  ! get initial CPU time
  call SYSTEM_CLOCK(COUNT1, COUNT_RATE, COUNT_MAX)

  ! Begin integration: --------------------------------------------------------------

  WRITE(*,*) ''
  WRITE(*,*) 'integration time = ', mtime
  WRITE(*,*) 'time step length = ', time_step_len
  WRITE(*,*) ''

  imin   = ifirst_minute
  ihour  = ifirst_hour
  iday   = ifirst_day
  imon   = ifirst_month
  iyear  = ifirst_year
  ihour  = ihour - 1
  icount_hour = 1

  integration : DO jtime=1, mtime - 24

     IF(icount_hour == 3600/NINT(time_step_len)) THEN
        icount_hour  = 1
        imin = 60/(3600/NINT(time_step_len))
     ELSE
        icount_hour = icount_hour + 1
        ihour = ihour + 1
        imin  = 0
     END IF
     IF(imin > 60) THEN
        ihour = ihour + 1
        imin  = 0
     END IF
     IF(ihour > 23) THEN
        iday = iday + 1
        ihour= 0
        imin = 0
     END IF
     IF(iday > icalender(imon)) THEN
        imon = imon + 1
        iday = 1
        ihour= 0
        imin = 0
     END IF
     IF(imon > 12) THEN
        iyear = iyear + 1
        imon = 1
        iday = 1
        ihour= 0
        imin = 0
     END IF
     IF(imin < 10) THEN
        WRITE(cmin,'(A1,I1)') '0',imin
     ELSE
        WRITE(cmin,'(I2)') imin
     END IF
     IF(ihour < 10) THEN
        WRITE(chour,'(A1,I1)') '0',ihour
     ELSE
        WRITE(chour,'(I2)') ihour
     END IF
     IF(iday < 10) THEN
        WRITE(cday,'(A1,I1)') '0',iday
     ELSE
        WRITE(cday,'(I2)') iday
     END IF
     IF(imon < 10) THEN
        WRITE(cmon,'(A1,I1)') '0',imon
     ELSE
        WRITE(cmon,'(I2)') imon
     END IF
     IF(iyear < 100) THEN
        WRITE(cyear,'(A2,I2)') '20',iyear
     ELSE
        WRITE(cyear,'(I4)') iyear
     END IF

     WRITE(*,*) 'integration time step = ', cday,'.',cmon,'.',cyear,'.',chour,':',cmin

     Do jmod = 1, nmod

     !--- Temperature [K]:

     ztemp(1:kproma,1) = yi(iday, 1, jmod) + 273.15_dp ! T [dg. C] to [K]

     !--- Relative humidity [%] TO [0-1]

     zrhum(1:kproma,1) = yi(iday, 2, jmod) / 100._dp

     ! pressure at middle of box [hPa] to [Pa]

     zpress(1:kproma,1) = yi(iday, 3, jmod) * 100._dp

     !--- Layer thickness zdz=dp/(rho*g) [m]:

     zdz(1:kproma,1) = 1._dp

     !--- Sulfate mass [molecules cm-3]
     !--- Black Carbon [ug m-3 (air)]
     !--- Organic Carbon [ug m-3 (air)]
     !--- Sea Salt [ug m-3 (air)]
     !--- Dust [ug m-3 (air)]
     !--- Particle numbers [N cm-3 (air)]

     zgso4 (1:kproma,1) = yi(iday, 4, jmod) * factor(NS)

     !--------------------------------------------------------------------
     SELECT CASE(jmod)
     !--------------------------------------------------------------------
     CASE(NS) ! Nucleation mode soluble
     zaerml(1:kproma,1,iso4ns) = yi(iday, 4, jmod) * factor(jmod)
     zaernl(1:kproma,1,inucs)  = yi(iday, 9, jmod) * factor(jmod)
     !--------------------------------------------------------------------
     CASE(KS) ! Aitken mode soluble
     zaerml(1:kproma,1,iso4ks) = yi(iday, 4, jmod) * factor(jmod)
     zaerml(1:kproma,1,ibcks)  = yi(iday, 5, jmod) / factor(jmod)
     zaerml(1:kproma,1,iocks)  = yi(iday, 6, jmod) / factor(jmod)
     zaernl(1:kproma,1,iaits)  = yi(iday, 9, jmod) * factor(jmod)
     !--------------------------------------------------------------------
     CASE(AS) ! Accumulation mode
     zaerml(1:kproma,1,iso4as) = yi(iday, 4, jmod) * factor(jmod)
     zaerml(1:kproma,1,ibcas)  = yi(iday, 5, jmod) / factor(jmod)
     zaerml(1:kproma,1,iocas)  = yi(iday, 6, jmod) / factor(jmod)
     zaerml(1:kproma,1,issas)  = yi(iday, 7, jmod) / factor(jmod)
     zaerml(1:kproma,1,iduas)  = yi(iday, 8, jmod) / factor(jmod)
     zaernl(1:kproma,1,iaccs)  = yi(iday, 9, jmod) * factor(jmod)
     !--------------------------------------------------------------------
     CASE(CS) ! Coarse mode soluble
     zaerml(1:kproma,1,iso4cs) = yi(iday, 4, jmod) * factor(jmod)
     zaerml(1:kproma,1,ibccs)  = yi(iday, 5, jmod) / factor(jmod)
     zaerml(1:kproma,1,ioccs)  = yi(iday, 6, jmod) / factor(jmod)
     zaerml(1:kproma,1,isscs)  = yi(iday, 7, jmod) / factor(jmod)
     zaerml(1:kproma,1,iducs)  = yi(iday, 8, jmod) / factor(jmod)
     zaernl(1:kproma,1,icoas)  = yi(iday, 9, jmod) * factor(jmod)
     !--------------------------------------------------------------------
     CASE(KI) ! Aitken mode insoluble
     zaerml(1:kproma,1,ibcki)  = yi(iday, 5, jmod) / factor(jmod)
     zaerml(1:kproma,1,iocki)  = yi(iday, 6, jmod) / factor(jmod)
     zaernl(1:kproma,1,iaiti)  = yi(iday, 9, jmod) * factor(jmod)
     !--------------------------------------------------------------------
     CASE(AI) ! Accumulation mode insoluble
     zaerml(1:kproma,1,iduai)  = yi(iday, 8, jmod) / factor(jmod)
     zaernl(1:kproma,1,iacci)  = yi(iday, 9, jmod) * factor(jmod)
     !--------------------------------------------------------------------
     CASE(CI) ! Coarse mode insoluble
     zaerml(1:kproma,1,iduci)  = yi(iday, 8, jmod) / factor(jmod)
     zaernl(1:kproma,1,icoai)  = yi(iday, 9, jmod) * factor(jmod)
     !--------------------------------------------------------------------
     END SELECT
     !--------------------------------------------------------------------

     END DO ! jmod

     IF(l_m7) &
     CALL  m7_main(kproma, kbdim,       klev,   time_step_len,  &  ! ECHAM indices
                   zpress, zrhum,       ztemp,                  &  !   "   thermodynamics
                           zgso4,       zaerml, zaernl,         &  !  M7   tracers
                   zm6rp,  zm6dry,      zrhop,  zww             )  !   "   aerosol properties 

     Do jmod = 1, nmod

        yo(1:kproma, 1, jmod, jtime) = ztemp (1:kproma,1) - 273.15_dp ! T  [dg. C]
        yo(1:kproma, 2, jmod, jtime) = zrhum (1:kproma,1) * 100._dp   ! RH [%]
        yo(1:kproma, 3, jmod, jtime) = zpress(1:kproma,1) / 100._dp   ! P  [hPa]

        ! Scaling factor for different modes (box model only)

        !--- SO4 = Sulfate mass [molecules kg-1]
        !--- BC  = Black Carbon [ug m-3 (air)]
        !--- OC  = Organic Carbon [ug m-3 (air)]
        !--- SS  = Sea Salt [ug m-3 (air)]
        !--- DU  = Dust [ug m-3 (air)]
        !--- N   = Particle numbers [N cm-3 (air)]
        !--- zm6dry = Dry count mean radius from [cm] to [m]
        !--- zm6rp  = Ambient count median radius from [cm] to [m]
        !--- zrhop  = Mean mode density from [g/cm3] to [kg/m3]
        !--- zww    = Aerosol water [kg m-3(air)]

        yo(1:kproma, 4, jmod, jtime) = zgso4 (1:kproma,1)

        SELECT CASE(jmod)

        CASE(inucs)

        yo(1:kproma, 5, jmod, jtime) = zaerml(1:kproma,1,iso4ns)  ! SO4
        yo(1:kproma, 6, jmod, jtime) = zero                       ! BC
        yo(1:kproma, 7, jmod, jtime) = zero                       ! OC
        yo(1:kproma, 8, jmod, jtime) = zero                       ! SS
        yo(1:kproma, 9, jmod, jtime) = zero                       ! DU
        yo(1:kproma,10, jmod, jtime) = zaernl(1:kproma,1,inucs)   ! N

        CASE(iaits)

        yo(1:kproma, 5, jmod, jtime) = zaerml(1:kproma,1,iso4ks)  ! SO4
        yo(1:kproma, 6, jmod, jtime) = zaerml(1:kproma,1,ibcks)   ! BC
        yo(1:kproma, 7, jmod, jtime) = zaerml(1:kproma,1,iocks)   ! OC
        yo(1:kproma, 8, jmod, jtime) = zero                       ! SS
        yo(1:kproma, 9, jmod, jtime) = zero                       ! DU
        yo(1:kproma,10, jmod, jtime) = zaernl(1:kproma,1,iaits)   ! N

        CASE(iaccs)

        yo(1:kproma, 5, jmod, jtime) = zaerml(1:kproma,1,iso4as)  ! SO4
        yo(1:kproma, 6, jmod, jtime) = zaerml(1:kproma,1,ibcas)   ! BC
        yo(1:kproma, 7, jmod, jtime) = zaerml(1:kproma,1,iocas)   ! OC
        yo(1:kproma, 8, jmod, jtime) = zaerml(1:kproma,1,issas)   ! SS
        yo(1:kproma, 9, jmod, jtime) = zaerml(1:kproma,1,iduas)   ! DU
        yo(1:kproma,10, jmod, jtime) = zaernl(1:kproma,1,iaccs)   ! N

        CASE(icoas)

        yo(1:kproma, 5, jmod, jtime) = zaerml(1:kproma,1,iso4cs)  ! SO4
        yo(1:kproma, 6, jmod, jtime) = zaerml(1:kproma,1,ibccs)   ! BC
        yo(1:kproma, 7, jmod, jtime) = zaerml(1:kproma,1,ioccs)   ! OC
        yo(1:kproma, 8, jmod, jtime) = zaerml(1:kproma,1,isscs)   ! SS
        yo(1:kproma, 9, jmod, jtime) = zaerml(1:kproma,1,iducs)   ! DU
        yo(1:kproma,10, jmod, jtime) = zaernl(1:kproma,1,icoas)   ! N

        CASE(iaiti)

        yo(1:kproma, 5, jmod, jtime) = zero                       ! SO4
        yo(1:kproma, 6, jmod, jtime) = zaerml(1:kproma,1,ibcki)   ! BC
        yo(1:kproma, 7, jmod, jtime) = zaerml(1:kproma,1,iocki)   ! OC
        yo(1:kproma, 8, jmod, jtime) = zero                       ! SS
        yo(1:kproma, 9, jmod, jtime) = zero                       ! DU
        yo(1:kproma,10, jmod, jtime) = zaernl(1:kproma,1,iaiti)   ! N

        CASE(iacci)

        yo(1:kproma, 5, jmod, jtime) = zero                       ! SO4
        yo(1:kproma, 6, jmod, jtime) = zero                       ! BC
        yo(1:kproma, 7, jmod, jtime) = zero                       ! OC
        yo(1:kproma, 8, jmod, jtime) = zero                       ! SS
        yo(1:kproma, 9, jmod, jtime) = zaerml(1:kproma,1,iduai)   ! DU
        yo(1:kproma,10, jmod, jtime) = zaernl(1:kproma,1,iacci)   ! N

        CASE(icoai)

        yo(1:kproma, 5, jmod, jtime) = zero                       ! SO4
        yo(1:kproma, 6, jmod, jtime) = zero                       ! BC
        yo(1:kproma, 7, jmod, jtime) = zero                       ! OC
        yo(1:kproma, 8, jmod, jtime) = zero                       ! SS
        yo(1:kproma, 9, jmod, jtime) = zaerml(1:kproma,1,iduci)   ! DU
        yo(1:kproma,10, jmod, jtime) = zaernl(1:kproma,1,icoai)   ! N

        END SELECT

        IF(jmod <= nsol) THEN
        yo(1:kproma,11, jmod, jtime) = zm6dry(1:kproma,1,jmod)/100._dp
        ELSE
        yo(1:kproma,11, jmod, jtime) = zm6rp (1:kproma,1,jmod)/100._dp
        END IF
        yo(1:kproma,12, jmod, jtime) = zm6rp (1:kproma,1,jmod)/100._dp
        yo(1:kproma,13, jmod, jtime) = zrhop (1:kproma,1,jmod)*1.e3_dp
        yo(1:kproma,14, jmod, jtime) = zww   (1:kproma,1,jmod)*1.e9_dp

     END DO ! jmod

  ! Stop integration: --------------------------------------------------------------
  END DO integration

  ! get final CPU time
  CALL SYSTEM_CLOCK(COUNT2, COUNT_RATE, COUNT_MAX)

  DEALLOCATE(zgso4)
  DEALLOCATE(zdz)
  DEALLOCATE(ztemp)
  DEALLOCATE(zrhum)
  DEALLOCATE(zaerml)
  DEALLOCATE(zaernl)
  DEALLOCATE(zm6rp)
  DEALLOCATE(zm6dry)
  DEALLOCATE(zrhop)
  DEALLOCATE(zww)
  DEALLOCATE(zpress)

  WRITE(*,*) ' '

  CALL end_message(TRIM(modstr),'call m7_main', substr)

  END SUBROUTINE m7_physc

  ! --------------------------------------------------------------------------

  SUBROUTINE m7_output

     USE messy_m7,           ONLY: nmod
     USE messy_main_blather, ONLY: start_message, end_message

     IMPLICIT NONE

     CHARACTER(LEN=*), PARAMETER         :: substr = 'm7_output'

     INTEGER                             :: jc, jmod, jtime, iprint

     CALL start_message(TRIM(modstr),'write output data to file', substr)

     WRITE(*,*) ' '
     WRITE(*,*) ' Output of ',modstr,'_',modver,' at ',TRIM(date),'-', TRIM(time)
     WRITE(*,*) ' '

     DO jmod = 1, nmod

     WRITE(ofile,'(A10,I1,A4)') 'output/m7_',jmod,'.dat'
     WRITE(*,*) ' write data to file ',ofile
     OPEN(unit=iwunit,file=ofile,status='unknown',form='formatted')
     WRITE(iwunit,*) '# Output of ',modstr,'_',modver,' at ',TRIM(date),'-', TRIM(time)
     WRITE(iwunit,*) '#',mtime,' data sets'
     WRITE(iwunit,*) &
    '# Date.Time              T [C],             RH [%]          P [hPa],      H2SO4 [molec/kg(air)], SO4 [molec/kg(air)],', &
    '        BC [ug/m3],         OC [ug/m3],         SS [ug/m3],          DU [ug/m3],        N [1/cm3],         r_dry [m],', &
    '          r_wet [m],         d_wet [kg/m3],     H2O [ug/m3]'

     imin   = ifirst_minute
     ihour  = ifirst_hour
     iday   = ifirst_day
     imon   = ifirst_month
     iyear  = ifirst_year
     ihour  = ihour - 1
     icount_hour = 1
     iprint = 0

     DO jtime=1, mtime

        IF(icount_hour == 3600/NINT(time_step_len)) THEN
           icount_hour  = 1
           imin = 60/(3600/NINT(time_step_len))
        ELSE
           icount_hour = icount_hour + 1
           ihour = ihour + 1
           imin  = 0
        END IF
        IF(imin > 60) THEN
           ihour = ihour + 1
           imin  = 0
        END IF
        IF(ihour > 23) THEN
           iday = iday + 1
           ihour= 0
           imin = 0
        END IF
        IF(iday > icalender(imon)) THEN
           imon = imon + 1
           iday = 1
           ihour= 0
           imin = 0
        END IF
        IF(imon > 12) THEN
           iyear  = iyear + 1
           imon = 1
           iday = 1
           ihour= 0
           imin = 0
        END IF
        IF(imin < 10) THEN
           WRITE(cmin,'(A1,I1)') '0',imin
        ELSE
           WRITE(cmin,'(I2)') imin
        END IF
        IF(ihour < 10) THEN
           WRITE(chour,'(A1,I1)') '0',ihour
        ELSE
           WRITE(chour,'(I2)') ihour
        END IF
        IF(iday < 10) THEN
           WRITE(cday,'(A1,I1)') '0',iday
        ELSE
           WRITE(cday,'(I2)') iday
        END IF
        IF(imon < 10) THEN
           WRITE(cmon,'(A1,I1)') '0',imon
        ELSE
           WRITE(cmon,'(I2)') imon
        END IF
        IF(iyear < 100) THEN
           WRITE(cyear,'(A2,I2)') '20',iyear
        ELSE
           WRITE(cyear,'(I4)') iyear
        END IF

        IF(jtime == 1 .AND. jmod == nmod) THEN
           WRITE(*,*) ' '
           WRITE(*,*) 'first output date: ', cday,'.',cmon,'.',cyear,'.',chour,':',cmin
           OPEN(unit=iounit,file='output/first_calculation_date',status='unknown',form='formatted')
           WRITE (iounit,*) cday,'.',cmon,'.',cyear,'.',chour,':',cmin
           CLOSE(iounit)
        END IF

!       IF(ihour == ihour) THEN
        IF(ihour == 12) THEN
           iprint = iprint + 1
           WRITE(iwunit,'(9A,19E20.9)') cday,'.',cmon,'.',cyear,'.',chour,':',cmin, (yo(1:kproma,jc,jmod,jtime),jc=1,ndo)
        END IF

     END DO

     CLOSE(iwunit)

     END DO ! nmod

     WRITE(*,*) 'last  output date: ', cday,'.',cmon,'.',cyear,'.',chour,':',cmin
     WRITE(*,*) 'output steps: ', iprint
     WRITE(*,*) ' '
     OPEN(unit=iwunit,file='output/last_calculation_date',status='unknown',form='formatted')
     WRITE (iwunit,*) cday,'.',cmon,'.',cyear,'.',chour,':',cmin
     CLOSE(iwunit)

     CALL end_message(TRIM(modstr),'write output data to file', substr)

  END SUBROUTINE m7_output

  ! --------------------------------------------------------------------------

END MODULE messy_m7_box

!*****************************************************************************

PROGRAM m7

  USE messy_m7_box, ONLY: m7_initialize, m7_physc, m7_output, &
                          date, time, zone, values, zero, mtime, &
                          COUNT_MAX, COUNT_RATE, COUNT1, COUNT2, &
                          yi, yo

  IMPLICIT NONE

  REAL                 :: CPUs

  ! get date and time

  CALL DATE_AND_TIME(date, time, zone, values)
  WRITE(*,*)' '
  WRITE(*,*) ' Start calculation at: ',TRIM(date),'-', TRIM(time)
  WRITE(*,*)' '

  CALL m7_initialize ! read CTRL namelist and input data

  ! calculate aerosol dynamics
  CALL m7_physc
  ! write output data for each mode to new file
  CALL m7_output

  CALL DATE_AND_TIME(date, time, zone, values)
  WRITE(*,*)' '
  WRITE(*,*) ' Stop calculation at: ',TRIM(date),'-', TRIM(time)

  ! Calculate performance
  IF(COUNT2.ge.COUNT1) THEN
     CPUs= FLOAT(COUNT2-COUNT1)/COUNT_RATE
  ELSE  ! COUNT has reached COUNT_MAX, where it starts again
     CPUs= FLOAT(COUNT_MAX-COUNT1+COUNT2)/COUNT_RATE
  END IF
  WRITE(*,*) ' CPUs/step  = ',CPUs/FLOAT(mtime)
  WRITE(*,*) ' CPUs total = ',CPUs,' steps = ',mtime

  DEALLOCATE(yo)
  DEALLOCATE(yi)

END PROGRAM m7

!*****************************************************************************
