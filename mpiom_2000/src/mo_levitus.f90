MODULE mo_levitus

  USE mo_kind, ONLY: i8, wp, dp
  USE mo_param1, ONLY: ie, je, ke, ie_g, je_g
  USE mo_parallel, ONLY: global_max, global_min, p_pe, p_io, &
       read_slice, stop_all, p_ioff, p_joff
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1, ONLY: almzer, ddpo, dt, dti, eminpo, i3drest, lmonts, lweto, &
       relsao, crelsal, creltem, reltho, saf, sao, sicomo, sicsno, sictho, &
       taf, tho, weto, zo, iweto
  USE mo_io_config, ONLY: next_free_unit
  USE mo_planetary_constants, ONLY: rhoicwa, rhosnwa
  USE mo_units, ONLY: io_stdout

  IMPLICIT NONE

  PRIVATE

  LOGICAL, ALLOCATABLE ::  LTLEV(:,:,:), LSLEV(:,:,:)
  REAL(dp), ALLOCATABLE ::  SLEVI(:,:,:), TLEVI(:,:,:)

  INTEGER :: spongezone(6)
  REAL(wp) :: spzndamp_time

  !HH    IO_IN_SURS SURSAL  SEA SURFACE SALINITY LEVITUS
  INTEGER :: io_in_surs, io_in_init
  !> inisal  salinity levitus
  INTEGER :: io_in_inis
  !> input filehandle for surface temperature
  INTEGER :: io_in_surt
  PUBLIC :: init_levitus_2d ,init_levitus_3d, relax_ts, relax_surf, levitus_set_ts, &
       spongezone, spzndamp_time, slevi, tlevi, &
       levitus_read_3d_restore,  levitus_horizontal_stratification, &
       levitus_read_3d_stratification, levitus_read_surface_salinity, &
       levitus_read_surface_temperature, &
       levitus_per_month_setup
CONTAINS

  SUBROUTINE init_levitus_3d

    ALLOCATE(SLEVI(IE,JE,KE),TLEVI(IE,JE,KE))
    ALLOCATE(LSLEV(IE,JE,KE),LTLEV(IE,JE,KE))

    slevi(:,:,:) = 0.0_wp
    tlevi(:,:,:) = 0.0_wp
    LSLEV(:,:,:)=.FALSE.
    LTLEV(:,:,:)=.FALSE.

    IF (p_pe==p_io) THEN
#ifndef LITTLE_ENDIAN
      io_in_init = next_free_unit()
      OPEN(IO_IN_INIT,FILE='INITEM',STATUS='UNKNOWN',                 &
                     ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      io_in_inis = next_free_unit()
      OPEN(IO_IN_INIS,FILE='INISAL',STATUS='UNKNOWN',                 & 
                     ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
#else
#ifndef NOENDIANCONVERT
      io_in_init = next_free_unit()
      OPEN(IO_IN_INIT,FILE='INITEM',STATUS='UNKNOWN',                 &
                     ACCESS='SEQUENTIAL',FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
      io_in_inis = next_free_unit()
      OPEN(IO_IN_INIS,FILE='INISAL',STATUS='UNKNOWN',                 &
                     ACCESS='SEQUENTIAL',FORM='UNFORMATTED',CONVERT='BIG_ENDIAN')
#else
    ! ERROR: compiler does not support convert='big_endian'
#endif
#endif
    END IF ! p_pe == p_io

  END SUBROUTINE init_levitus_3d

  SUBROUTINE init_levitus_2d

    IF (p_pe==p_io) THEN
      io_in_surs = next_free_unit()
      OPEN(IO_IN_SURS,FILE='SURSAL',STATUS='UNKNOWN',                 &
           &          ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      IF (creltem .GT. ALMZER) THEN
        io_in_surt = next_free_unit()
        OPEN(io_in_surt, file='SURTEM', status='UNKNOWN',             &
             &            access='SEQUENTIAL', form='UNFORMATTED')
      END IF
    END IF ! p_pe == p_io

  END SUBROUTINE init_levitus_2d

  SUBROUTINE free_levitus
    DEALLOCATE(SLEVI,TLEVI)
    DEALLOCATE(LSLEV,LTLEV)
  END SUBROUTINE free_levitus

  SUBROUTINE levitus_per_month_setup(lmont)
    INTEGER, INTENT(in) :: lmont
    INTEGER :: i
#ifdef RESTORE_MON
    !:: REWIND SURSAL TO BEGIN OF YEAR
    IF(p_pe==p_io) THEN
      REWIND(IO_IN_SURS)
      IF (creltem .GT. almzer) REWIND(io_in_surt)
    END IF
    !:: READ APPROPRIATE SURSAL FOR MONTHLY RESTORING
    WRITE(IO_STDOUT,*)                                               &
         'READING MONTHLY SURFACE SALINITY IN MONTH ',LMONT
    DO i = 1, lmont
      CALL levitus_read_surface_salinity
      IF (creltem .GT. almzer)                                         &
           CALL levitus_read_surface_temperature
    END DO
#endif /*RESTORE_MON*/
    !::
    IF (I3DREST .GT. 1) THEN
      !:: REWIND INISAL/INITEM TO BEGIN OF YEAR
      IF(p_pe==p_io) THEN
        REWIND(IO_IN_INIT)
        REWIND(IO_IN_INIS)
      END IF
      !:: READ APPROPRIATE INISAL/INITEM FOR MONTHLY RESTORING
      WRITE(IO_STDOUT,*)                                               &
           'READING MONTHLY LEVITUS FIELDS IN MONTH ',LMONT
      DO i = 1, lmont
        CALL levitus_read_3d_restore
      END DO
    END IF ! I3DREST
  END SUBROUTINE levitus_per_month_setup


  SUBROUTINE relax_ts


    IF (I3DREST .GT. 0) THEN
      WHERE (ltlev(:,:,:)) tho(:,:,:)=tho(:,:,:)+dt*spzndamp_time*(tlevi(:,:,:)-tho(:,:,:))
      WHERE (lslev(:,:,:)) sao(:,:,:)=sao(:,:,:)+dt*spzndamp_time*(slevi(:,:,:)-sao(:,:,:))
    ENDIF

    CALL bounds_exch(1,'p',tho,'relax_ts 3')
    CALL bounds_exch(1,'p',sao,'relax_ts 4')

  END SUBROUTINE relax_ts



  !> surface relaxation, only for uncoupled runs
  !> computes boundary forcing on salt, temperature and zeta
  SUBROUTINE relax_surf

    INTEGER    i,j
    REAL(wp)       reiscon,oldso

    !     boundary forcing on temperature, salt and zeta
    !     relaxation time on salinity   : 1./ relsal

!$omp parallel private (i,j,reiscon,oldso)
!$omp do

    DO j=1,je
      DO i=1,ie

        eminpo(i,j) = 0._wp

        reiscon = 0._wp
        IF (sicomo(i, j) .LE. 0.01_wp) THEN
          reiscon = 1._wp
        ENDIF

#ifdef EISREST
        reiscon = 1._wp - sicomo(i,j)
#endif
        oldso=sao(i,j,1)
        sao(i,j,1)=sao(i,j,1)+dt*crelsal*reiscon*(relsao(i,j)-sao(i,j,1))
        tho(i,j,1)=tho(i,j,1)+dt*creltem*(reltho(i,j)-tho(i,j,1))

#ifdef ANOMALY_FORCING
        creltem=3.0e-6
        IF (ibek(i,j) .EQ. 1) THEN
          tho(i,j,1) = tho(i,j,1) &
               + dt * creltem * (reltho(i, j) - 3.0_wp - tho(i,j,1))
        ENDIF
#endif
        eminpo(i,j)=(ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa             &
             -sicsno(i,j)*rhosnwa)                                 &
             *(MAX(oldso, 1.e-3_wp)/MAX(sao(i,j,1), 1.e-3_wp) - 1._wp)
      ENDDO
    ENDDO
!$omp end do

!$omp do
    DO j=1,je
      DO i=1,ie
        zo(i,j)=(zo(i,j)+eminpo(i,j))*weto(i,j,1)
        eminpo(i,j)=eminpo(i,j)*dti                       !uwe eminpo in m/s
      ENDDO
    ENDDO
!$omp end do
!$omp end parallel

  END SUBROUTINE relax_surf


  SUBROUTINE levitus_set_ts
    INTEGER :: k, j, i
    DO K=1,KE
      DO J=1,JE
        DO I=1,IE
          THO(I,J,K)=TLEVI(I,J,K)
          SAO(I,J,K)=SLEVI(I,J,K)
        END DO
      END DO
    END DO
  END SUBROUTINE levitus_set_ts

  SUBROUTINE levitus_read_3d_restore
    ! 3d restoring mask
    INTEGER(kind=i8) :: ibla(4)
    INTEGER :: ist, ien, jst, jen, k, kst, ken, ierr

    WRITE(IO_STDOUT,*)                                             &
         &       '==> READING LEVITUS STRATIFICATION FOR 3D RESTORING.'
    IF ( spongezone(1) .GE. 1 .AND. spongezone(1) .LE. ie_g .AND.  &
         spongezone(4) .GE. 1 .AND. spongezone(4) .LE. ie_g .AND.  &
         spongezone(4) .GE. spongezone(1) .AND.                  &
         spongezone(2) .GE. 1 .AND. spongezone(2) .LE. je_g .AND.  &
         spongezone(5) .GE. 1 .AND. spongezone(5) .LE. je_g .AND.  &
         spongezone(5) .GE. spongezone(2) .AND.                  &
         spongezone(3) .GE. 1 .AND. spongezone(3) .LE. ke .AND.    &
         spongezone(6) .GE. 1 .AND. spongezone(6) .LE. ke .AND.    &
         spongezone(6) .GE. spongezone(3) ) THEN

      DO k=1,ke
        IF(p_pe==p_io) THEN
          READ(io_in_inis,iostat=ierr) ibla
          IF (ierr.GT.0) THEN
            WRITE(io_stdout,*)'no monthly levitus fields found in month= ',lmonts
            WRITE(io_stdout,*)                                             &
                 'check option dreidrest_mon and files inisal/tem !!! '
            CALL stop_all('no monthly levitus fields found => run aborted')
          ENDIF
        ENDIF

        CALL read_slice(io_in_inis,slevi(:,:,k))

        IF(p_pe==p_io) READ(io_in_init) ibla
        CALL read_slice(io_in_init,tlevi(:,:,k))

      ENDDO

      CALL bounds_exch(1,'p',tlevi,'levire 10')
      CALL bounds_exch(1,'p',slevi,'levire 11')


      IF (spongezone(1)-p_ioff <= ie) THEN ! start restoring zone here or left

        ist=MAX(1,spongezone(1)-p_ioff)
        ien=MIN(spongezone(4)-p_ioff,ie)
      ELSE
        ist = -1
        ien = -1
      ENDIF

      IF (spongezone(2)-p_joff <= je) THEN ! start restoring zone here or top

        jst=MAX(1,spongezone(2)-p_joff)
        jen=MIN(spongezone(5)-p_joff,je)

      ELSE

        jst= -1
        jen= -1

      ENDIF

      kst=spongezone(3)
      ken=spongezone(6)

      lslev(:,:,:)=.FALSE.
      ltlev(:,:,:)=.FALSE.

      IF  ( ist >= 1 .AND. ien >= 1 .AND. jst >= 1 .AND. jen >= 1 )   THEN

        lslev(ist:ien,jst:jen,kst:ken) = &
             ( (slevi(ist:ien,jst:jen,kst:ken) .LE. 40._wp) &
             .AND. (slevi(ist:ien,jst:jen,kst:ken) .GE. 0._wp) )

        ltlev(ist:ien,jst:jen,kst:ken) = &
             ( (tlevi(ist:ien,jst:jen,kst:ken) .LE. 40._wp) &
             .AND. (tlevi(ist:ien,jst:jen,kst:ken) .GE. -3.0_wp) )

      ENDIF

      WHERE (weto(:,:,:) < 0.5_wp) lslev(:,:,:)=.FALSE.
      WHERE (weto(:,:,:) < 0.5_wp) ltlev(:,:,:)=.FALSE.

      WRITE(io_stdout,*) &
           'levitus data read to slevi/tlevi for 3d restoring in month', &
           lmonts

    ELSE

      WRITE(io_stdout,*) 'levitus_stratification_3d_restore called', &
           ' but spongezone is invalid. the program is stopped.'
      CALL stop_all('levitus_stratification_3d_restore call fail => run aborted')
    ENDIF

  END SUBROUTINE levitus_read_3d_restore

  SUBROUTINE levitus_horizontal_stratification
    INTEGER :: i, j
    WRITE(IO_STDOUT,*)                                             &
         &       '==> HORIZONTAL INITIAL STRATIFICATION'
    DO j=1,je
      DO i=1,ie
        sao(i,j,1:ke)=saf(1:ke)
        tho(i,j,1:ke)=taf(1:ke)
      ENDDO
    ENDDO
  END SUBROUTINE levitus_horizontal_stratification

  SUBROUTINE levitus_read_3d_stratification

    REAL(wp) :: saomax(0:1), saomin(0:1), thomax(0:1), thomin(0:1)
    INTEGER(kind=i8) :: ibla(4)
    INTEGER :: i, j, k

    WRITE(io_stdout, *)                                            &
         &       '==> READING INITIAL 3D STRATIFICATION'
         
!mz_ap_20100830+
#ifndef MESSY
    IF(p_pe==p_io) REWIND(io_in_init)
    IF(p_pe==p_io) REWIND(io_in_inis)
#endif
!mz_ap_20100830-

    DO K=1,KE
      IF(p_pe==p_io) READ(io_in_inis)
      CALL read_slice(io_in_inis, sao(:,:,k))

      IF(p_pe==p_io) READ(IO_IN_INIT)
      CALL read_slice(io_in_init, tho(:,:,k))
    ENDDO

    saomax(0) = -1.0_wp
    saomin(0) = 99.0_wp
    saomax(1) = -1.0_wp
    saomin(1) = 99.0_wp
    thomax(0) = -1.0_wp
    thomin(0) = 99.0_wp
    thomax(1) = -1.0_wp
    thomin(1) = 99.0_wp

    DO K=1,KE
      DO J=1,JE
        DO I=1,IE
          saomax(iweto(i, j, k)) = MAX(saomax(iweto(i, j, k)), sao(i, j, k))
          saomin(iweto(i, j, k)) = MIN(saomin(iweto(i, j, k)), sao(i, j, k))
          thomax(iweto(i, j, k)) = MAX(thomax(iweto(i, j, k)), tho(i, j, k))
          thomin(iweto(i, j, k)) = MIN(thomin(iweto(i, j, k)), tho(i, j, k))
        ENDDO
      ENDDO
    ENDDO

    CALL global_max(saomax(0), saomax(1), thomax(0), thomax(1))
    CALL global_min(saomin(0), saomin(1), thomin(0), thomin(1))

#ifdef MESSY
         IF (p_pe==p_io) THEN
#endif
    WRITE(IO_STDOUT,*) 'MIN/MAX AT DRY CELLS (S): ', saomin(0), saomax(0)
    WRITE(IO_STDOUT,*) 'MIN/MAX AT WET CELLS (S): ', saomin(1), saomax(1)
    WRITE(IO_STDOUT,*) 'MIN/MAX AT DRY CELLS (T): ', thomin(0), thomax(0)
    WRITE(IO_STDOUT,*) 'MIN/MAX AT WET CELLS (T): ', thomin(1), thomax(1)
#ifdef MESSY
         ENDIF
#endif

  END SUBROUTINE levitus_read_3d_stratification

  SUBROUTINE levitus_read_surface_salinity

    INTEGER(kind=i8) :: ibla(4)
    INTEGER :: ierr

    WRITE(IO_STDOUT,*)                                             &
         &       '==> READING SURFACE SALINITY'
    IF(p_pe==p_io) THEN
      READ(io_in_surs,iostat=ierr) ibla
      IF (IERR.GT.0) THEN
        WRITE(IO_STDOUT,*)'NO MONTHLY SSS FIELDS FOUND IN MONTH= '  &
             &                        ,LMONTS
        WRITE(IO_STDOUT,*)                                          &
             &            'CHECK OPTION RESTORE_MON AND FILE SURSAL!!! '
        CALL stop_all('NO MONTHLY SSS FIELDS FOUND')
      ENDIF
    ENDIF
    CALL read_slice(io_in_surs, relsao)

  END SUBROUTINE levitus_read_surface_salinity

  SUBROUTINE levitus_read_surface_temperature

    INTEGER(kind=i8) :: ibla(4)
    INTEGER :: ierr

    WRITE(IO_STDOUT,*)                                             &
         &       '==> READING SURFACE TEMPERATURE'
    IF(p_pe==p_io) THEN
      READ(io_in_surt, IOSTAT=ierr) ibla
      IF (IERR.GT.0) THEN
        WRITE(IO_STDOUT,*)'NO MONTHLY SST FIELDS FOUND IN MONTH= '  &
             &                        ,LMONTS
        WRITE(IO_STDOUT,*)                                          &
             &            'CHECK OPTION RESTORE_MON AND FILE SURTEM!!! '
        CALL stop_all('NO MONTHLY SST FIELDS FOUND')
      ENDIF
    ENDIF
    CALL read_slice(io_in_surt, reltho)

  END SUBROUTINE levitus_read_surface_temperature


END MODULE mo_levitus
