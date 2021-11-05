!> module mo_restart
!> unifies subroutines used in writing/reading of restart files
!> joins former global subroutines aufr and aufw
!> @date 2009-04-29
MODULE mo_restart
  USE mo_kind, ONLY: wp
#ifndef MESSY
  USE mo_iso_c_kinds, ONLY: c_int64_t
#else
  USE mo_kind, ONLY : c_int64_t => i8, c_float => sp
#endif
  USE mo_param1, ONLY: ie, ie_g, je, je_g, ke, kep
  USE mo_parallel, ONLY: p_bcast, read_slice, write_slice, &
       p_pe, p_io
  USE mo_boundsexch, ONLY : bounds_exch
  USE mo_commo1, ONLY: uoo, uko, voe, vke, tho, sao, dvo, avo, wo, tice, &
       hibete, hibeto, hibzete, hibzeto, depto, ldays, lmonts, lyears, &
       sicomo, sicsno, sictho, sicuo, sicve, z1o, zo, istart, &
       istart_new_topo_update, lday, lmonth, tiestu, tiestw, uaccel, vaccel
!!$#ifdef MIST
!!$  USE mo_diagnosis, ONLY:  amld, psitro
!!$#endif
  USE mo_units, ONLY: io_in_z380, io_in_z370, io_stdout
#ifdef __coupled
  USE mo_fluxes1
#endif
  IMPLICIT NONE
  PRIVATE
  !> unit currently in use
  INTEGER :: iunit
  !> offset of current unit
  !> switch for respective restart files
  INTEGER :: iflag
  !
  !> set start and stop dates
  INTEGER :: ly_start, lm_start, ly_end
  !
  !> permit access via subroutines only
  PUBLIC :: setup_restart_defaults, setup_restart, aufr, aufw,  &
       ly_start, lm_start, ly_end, dump_restart_data_extra
CONTAINS
  !> setup default values to be optionally overwritten by user
  SUBROUTINE setup_restart_defaults
    !JJ    DEFAULT END OF RUN IN YEAR
    LY_END = 5000
    !JJ    DEFAULT OFFSET FOR YEAR COUNTER (NEG. VALUE == NO ACTION TAKEN)
    LY_START = -999
    LM_START = -999
  END SUBROUTINE setup_restart_defaults
  !> set logical unit for restart file; for first run write to io_in_z370
  SUBROUTINE setup_restart(iaufr)
    INTEGER, INTENT(in) :: iaufr
    IF (iaufr .EQ. 0) THEN
      IUNIT = IO_IN_Z380
      IFLAG = -1
    ELSE
      iflag = 1
    END IF
  END SUBROUTINE setup_restart


  !****************************************************************
  !
  !     AAAAAA  U    U  FFFFFF  RRRRRR
  !     A    A  U    U  F       R    R
  !     AAAAAA  U    U  FFFFF   RRRRRR
  !     A    A  U    U  F       R  RR
  !     A    A  UUUUUU  F       R   RR
  !
  !
  !*****************************************************************
  !
  !-----------------------------------------------------------------
  !
  SUBROUTINE AUFR(ldt)
    INTEGER, INTENT(inout) :: ldt

    INTEGER(KIND=C_INT64_T) IDATE,ICODE
    REAL(wp) RDT27, RDT28, RDT
    INTEGER i, j, k
    ! alternating restart information
    !HH    IO_IN_Z370 Z37000  RESTART-INFORMATION
    !HH    IO_IN_Z380 Z38000  RESTART-INFORMATION
    !
    IF (p_pe==p_io) THEN
      OPEN(IO_IN_Z370,FILE='Z37000',STATUS='UNKNOWN'                  &
           &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      OPEN(IO_IN_Z380,FILE='Z38000',STATUS='UNKNOWN'                  &
           &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      REWIND(IO_IN_Z370)
      REWIND(IO_IN_Z380)
      !
      READ(IO_IN_Z370,ERR=99,END=99) IDATE
      READ(IO_IN_Z370,ERR=99,END=99) RDT27
      GOTO 1
99    rdt27 = 0.0_wp
      WRITE(IO_STDOUT,*)'NO DATA ON RESTART FILE Z37000 !!!'
1     READ(IO_IN_Z380,ERR=88,END=88) IDATE
      READ(IO_IN_Z380,ERR=88,END=88) RDT28
      GOTO 2
88    rdt28 = 0.0_wp
      WRITE(IO_STDOUT,*)'NO DATA ON RESTART FILE Z38000 !!!'
2     WRITE(IO_STDOUT,6000)RDT27,RDT28
6000  FORMAT(' Z37000 TIMESTEPS : ',F10.0                             &
           &        ,' Z38000 TIMESTEPS : ',F10.0)
      !
      IF(RDT27.GT.RDT28)THEN
        IUNIT=IO_IN_Z370
        IFLAG=1
      ELSE
        IUNIT=IO_IN_Z380
        IFLAG=-1
      ENDIF
      !
      REWIND(IUNIT)
      !
      READ(IUNIT) IDATE
      READ(IUNIT) RDT
    ENDIF

    CALL p_bcast(IFLAG,p_io)
    CALL p_bcast(IUNIT,p_io)
    CALL p_bcast(IDATE,p_io)
    CALL p_bcast(RDT,p_io)

    LDT=NINT(RDT)
    LYEARS=INT(IDATE/10000_C_INT64_T)
    LMONTS=INT((IDATE-INT(LYEARS*10000,C_INT64_T))/100_C_INT64_T)
    LDAYS=INT(IDATE - INT(LYEARS*10000, C_INT64_T) - INT(LMONTS*100, C_INT64_T))
    IF (istart .EQ. istart_new_topo_update) lyears=-1

    !:: SET YEAR TO LY_START
    IF(LY_START .GE. 0) THEN
      LYEARS=LY_START
      LMONTS=0
      LDAY=0
      WRITE(IO_STDOUT,*)'ATTN !! SET OFFSET --> LYEARS TO: ',LYEARS
    ENDIF
    !:: SET MONTH TO LM_START
    IF(LM_START .GE. 0) THEN
      LMONTS=LM_START-1
      LDAY=0
      WRITE(IO_STDOUT,*)'ATTN !! SET OFFSET --> LMONTH TO: ',LMONTH
    ENDIF
    !::
    WRITE(IO_STDOUT,6001) IUNIT,LYEARS,LMONTS,LDAYS,LDT
6001 FORMAT(' READ FROM UNIT : ',I5,                                   &
         &       ' START YEAR,MONTH,DAY,LDT : ',4I10)

    !  READ ZONAL VELOCITY COMPONENT
    DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) IDATE
      CALL read_slice(IUNIT,UOO(:,:,K))
    ENDDO
    CALL bounds_exch(1,'u',UOO,'aufr 1')
    uko(:,:,:)=uoo(:,:,:)

    !  READ MERIDIONAL  VELOCITY COMPONENT
    DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) IDATE
      CALL read_slice(IUNIT,VOE(:,:,K))
    ENDDO
    CALL bounds_exch(1,'v',VOE,'aufr 2')
    vke(:,:,:)=voe(:,:,:)


    !:: UPDATE VELOCITY FIELDS
    CALL OCTIMF
    !  READ TEMPERATURE
    DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) IDATE
      CALL read_slice(IUNIT,THO(:,:,K))
    ENDDO
    CALL bounds_exch(1,'p',THO,'aufr 3')

    !  READ SALINITY
    DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) IDATE
      CALL read_slice(IUNIT,SAO(:,:,K))
    ENDDO
    CALL bounds_exch(1,'p',SAO,'aufr 4')

    !:: READ 2-D FIELDS
    !  READ SEA SURFACE ELEVATION
    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,ZO)
    CALL bounds_exch(1,'p',ZO,'aufr 5')

    !  READ SEA SURFACE ELEVATION CHANGE
    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,Z1O)
    CALL bounds_exch(1,'p',Z1O,'aufr 6')

    !  READ SEA ICE THICKNESS
    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,SICTHO)
    CALL bounds_exch(1,'p',SICTHO,'aufr 7')

    !  READ SEA ICE CONCENTRATION
    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,SICOMO)
    CALL bounds_exch(1,'p',SICOMO,'aufr 8')

    !  READ ZONAL SEA ICE VELOCITY
    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,SICUO)
    CALL bounds_exch(1,'u',SICUO,'aufr 9')

    !  READ MERIDIONAL SEA ICE VELOCITY
    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,SICVE)
    CALL bounds_exch(1,'v',SICVE,'aufr 10')
    !  READ SNOW
    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,SICSNO)
    CALL bounds_exch(1,'p',SICSNO,'aufr 11')
    !  READ HIBLER ETA/ZETA FIELDS
    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,HIBETE)
    CALL bounds_exch(1,'s',hibete,'aufr 12')

    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,HIBETO)
    CALL bounds_exch(1,'p',hibeto,'aufr 13')


    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,HIBZETE)
    CALL bounds_exch(1,'s',hibzete,'aufr 14')

    IF(p_pe==p_io) READ(IUNIT) IDATE
    CALL read_slice(IUNIT,HIBZETO)
    CALL bounds_exch(1,'p',hibzeto,'aufr 15')

    !  READ VERTICAL DIFFUSIVITY DVO
    DO K=1,KEP
      IF(p_pe==p_io) READ(IUNIT) IDATE
      CALL read_slice(IUNIT,DVO(:,:,K))
    ENDDO
    CALL bounds_exch(1,'p',dvo,'aufr 16')

    !  READ VERTICAL FRICTION AVO
    DO K=1,KEP
      IF(p_pe==p_io) READ(IUNIT) IDATE
      CALL read_slice(IUNIT,AVO(:,:,K))
    ENDDO
    CALL bounds_exch(1,'p',avo,'aufr 17')

    !  READ VERTICAL VELOCITY (DIAGNOSTIC)
    DO K=1,KEP
      IF(p_pe==p_io) READ(IUNIT) IDATE
      CALL read_slice(IUNIT,WO(:,:,K))
    ENDDO
    CALL bounds_exch(1,'p',wo,'aufr 18')

    DO J=1,JE
      DO I=1,IE
        dvo(i, j, 1) = 0._wp
        dvo(i, j, kep) = 0._wp
        avo(i, j, 1) = 0._wp
        avo(i, j, kep) = 0._wp
      ENDDO
    ENDDO
    !  READ ICETEMPERATURE
    IF(p_pe==p_io) READ(IUNIT) IDATE,ICODE
    CALL p_bcast(icode, p_io)
    IF (icode .EQ. INT(99, C_INT64_T)) THEN
      CALL read_slice(IUNIT,TICE)
      CALL bounds_exch(1,'p',tice,'aufr 19')
    else
      tice(:,:) = 0._wp
    endif

    !  READ MOMENTUM induced velocity
    DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) IDATE,ICODE
      CALL p_bcast(icode,p_io)
      IF (icode .EQ. 9_C_INT64_T) THEN
        CALL read_slice(IUNIT,uaccel(:,:,K))
      ELSE
        uaccel(:,:,:) = 0._wp
        vaccel(:,:,:) = 0._wp
        IF(p_pe==p_io) WRITE(0,*) &
             'Warning: READ of uaccel and vaccel has failed.'
        IF(p_pe==p_io) WRITE(0,*) ' You may safely ignore this warning.'
        goto 1468
      ENDIF
    ENDDO
    CALL bounds_exch(1,'u',uaccel,'aufr 20')

    !  READ MOMENTUM induced velocity
    DO K=1,KE
      IF(p_pe==p_io) READ(IUNIT) IDATE,ICODE
      CALL p_bcast(icode,p_io)
      IF (icode.EQ.10_C_INT64_T) THEN
        CALL read_slice(IUNIT,vaccel(:,:,K))
      ELSE
        uaccel(:,:,:) = 0._wp
        vaccel(:,:,:) = 0._wp
        IF(p_pe==p_io) WRITE(0,*) &
             'Warning: READ of uaccel and vaccel has failed.'
        IF(p_pe==p_io) WRITE(0,*) &
             ' You may safely ignore this warning.'
        goto 1468
      ENDIF
    ENDDO
    CALL bounds_exch(1,'u',vaccel,'aufr 21')


1468 CONTINUE

#if 0
    DO K=1,KE
      DO J=1,JE
        DO I=1,IE
          !     UKE(I,J,K)=0.
          !      SAO(I,J,K)=MAX(SAO(I,J,K),28.)
          !      THO(I,J,K)=MAX(THO(I,J,K),-2.)
          sao(i, j, k) = MIN(sao(i, j, k), 70._wp)
          tho(i, j, k) = MIN(tho(i, j, k), 70._wp)
          !
          !     VKO(I,J,K)=0.
          !     ZE(I,J)=0.
          !     ZO(I,J)=0.
          !
        END DO
      END DO
    END DO
#endif
    !
    IF(p_pe==p_io) THEN
      REWIND(IUNIT)
      !
      CLOSE(IO_IN_Z370)
      CLOSE(IO_IN_Z380)
    ENDIF
  END SUBROUTINE AUFR

  !****************************************************************
  !
  !     AAAAAA  U    U  FFFFFF  W    W
  !     A    A  U    U  F       W    W
  !     AAAAAA  U    U  FFFFF   W WW W
  !     A    A  U    U  F       WWWWWW
  !     A    A  UUUUUU  F       WW  WW
  !
  !
  !*****************************************************************
  ! SUBROUTINE AUFW
  !
  !     THIS SBR WRITES ROUTINELY AT CERTAIN TIME INTERVALS
  !     A RESTART FILE ON Z37000 OR Z38000
  !
  !-----------------------------------------------------------------
  SUBROUTINE AUFW(ldt)
    INTEGER, INTENT(in) :: ldt
    INTEGER(kind=c_int64_t) :: idate
    REAL(wp) rldt, rldays
    !
    IF(p_pe==p_io) THEN
      OPEN(IO_IN_Z370,FILE='Z37000',STATUS='UNKNOWN'                  &
           &                 ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      OPEN(IO_IN_Z380,FILE='Z38000',STATUS='UNKNOWN'                  &
           &                 ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      REWIND(IO_IN_Z370)
      REWIND(IO_IN_Z380)
    ENDIF
    !
    IUNIT=IUNIT+IFLAG
    WRITE(IO_STDOUT,*) 'AUFW: IUNIT IFLAG',IUNIT,IFLAG
    IF(p_pe==p_io) REWIND(IUNIT)
    WRITE(IO_STDOUT,*)' ++++++ WRITING RESTART FILE ON ',IUNIT        &
         &,' AFTER ',LDT                                                    &
         &,' TIME STEPS . CALCULATED YEAR : ',LYEARS,' MONTH : ',LMONTS
    IDATE = INT((LYEARS*10000)+(LMONTS*100)+LDAYS, C_INT64_T)
    rldt = REAL(ldt, wp)
    rldays = REAL(ldays, wp)
    CALL dump_restart_data_extra(iunit, idate, rldt, rldays)
    IFLAG=-IFLAG
    IF(p_pe==p_io) THEN
      REWIND(IUNIT)
      CLOSE(IO_IN_Z370)
      CLOSE(IO_IN_Z380)
    ENDIF
    RETURN
  END SUBROUTINE AUFW

  SUBROUTINE dump_restart_data_extra(iunit, idate, rldt, rldays)
    INTEGER, INTENT(in) :: iunit
    INTEGER(kind=c_int64_t), INTENT(in) :: idate
    REAL(wp), OPTIONAL, INTENT(in) :: rldt, rldays
    INTEGER(KIND=c_int64_t) :: icode,kcode,iedim
    INTEGER :: k
#ifdef MIST
    INTEGER :: l
    INTEGER(KIND=c_int64_t) :: idate2
#endif
    !:: WRITE ALL DATA IN EXTRA FORMAT
    IF (PRESENT(rldt) .AND. PRESENT(rldays)) THEN
      !:: WRITE TIME STEP INFORMATION
      KCODE = 1_C_INT64_T
      IEDIM = 2_C_INT64_T
      ICODE = 999_C_INT64_T
      IF(p_pe==p_io) THEN
        WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
        WRITE(IUNIT) RLDT,RLDAYS
      END IF
    END IF
    !  WRITE ZONAL VELOCITY COMPONENT
    ICODE = 3_C_INT64_T
    IEDIM = INT(IE_G, C_INT64_T) * INT(JE_G, C_INT64_T)
    DO K=1,KE
      KCODE = INT(TIESTU(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,UKO(:,:,K))
    ENDDO
    !  WRITE MERIDIONAL  VELOCITY COMPONENT
    ICODE = 4_C_INT64_T
    DO K=1,KE
      KCODE = INT(TIESTU(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,VKE(:,:,K))
    ENDDO
    !  WRITE TEMPERATURE
    ICODE = 2_C_INT64_T
    DO K=1,KE
      KCODE = INT(TIESTU(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,THO(:,:,K))
    ENDDO
    !  WRITE SALINITY
    ICODE = 5_C_INT64_T
    DO K=1,KE
      KCODE = INT(TIESTU(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,SAO(:,:,K))
    ENDDO
    !:: WRITE 2-D FIELDS
    KCODE = 0_C_INT64_T
    !  WRITE SEA SURFACE ELEVATION
    ICODE = 1_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,ZO)
    !  WRITE SEA SURFACE ELEVATION CHANGE
    ICODE = 82_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,Z1O)
    !  WRITE SEA ICE THICKNESS
    ICODE = 13_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,SICTHO)
    !  WRITE SEA ICE CONCENTRATION
    ICODE = 15_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,SICOMO)
    !  WRITE ZONAL SEA ICE VELOCITY
    ICODE = 35_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,SICUO)
    !  WRITE MERIDIONAL SEA ICE VELOCITY
    ICODE = 36_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,SICVE)
    !  WRITE SNOW
    ICODE = 141_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,SICSNO)
    !  WRITE HIBLER ETA/ZETA FIELDS
    icode = 501_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,HIBETE)
    icode = 502_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,HIBETO)
    icode = 503_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,HIBZETE)
    icode = 504_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,HIBZETO)
    !  WRITE VERTICAL DIFFUSIVITY DVO
    icode = 111_C_INT64_T
    DO K=1,KEP
      KCODE = INT(TIESTW(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,DVO(:,:,K))
    ENDDO
    !  WRITE VERTICAL FRICTION AVO
    icode = 110_C_INT64_T
    DO K=1,KEP
      KCODE = INT(TIESTW(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,AVO(:,:,K))
    ENDDO
    !  WRITE VERTICAL VELOCITY (DIAGNOSTIC)
    icode = 7_C_INT64_T
    DO K=1,KEP
      KCODE = INT(TIESTW(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,WO(:,:,K))
    ENDDO
    icode = 99_C_INT64_T
    kcode = 0_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,TICE)
    !  WRITE MOMENTUM ACCELERATION
    ICODE=9_C_INT64_T
    DO K=1,KE
      KCODE=INT(TIESTU(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,UACCEL(:,:,K))
    ENDDO
    !  WRITE MOMENTUM ACCELERATION
    ICODE=10_C_INT64_T
    DO K=1,KE
      KCODE=INT(TIESTU(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,VACCEL(:,:,K))
    ENDDO

#ifdef MIST
    !:: WRITE 2-D DIAGNOSTIC FIELDS
    kcode = 0_C_INT64_T
    !:: WRITE ZONAL ICE TRANSPORTS
!!$    icode = 142_C_INT64_T
!!$    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
!!$    CALL write_slice(IUNIT,EISTRX)
!!$    !:: WRITE MERIDIONAL ICE TRANSPORTS
!!$    icode = 143_C_INT64_T
!!$    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
!!$    CALL write_slice(IUNIT,EISTRY)
    !:: WRITE MAX. CONVECTION DEPTH
    icode = 69_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,FLOAT(KCONDEP(:,:)))
!!$    !:: WRITE HEAT FLUX
!!$    icode = 70_C_INT64_T
!!$    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
!!$    CALL write_slice(IUNIT,HFLM)
!!$    !:: WRITE P-E
!!$    icode = 79_C_INT64_T
!!$    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
!!$    CALL write_slice(IUNIT,PMEM)
    !:: WRITE ZONAL WIND STRESS
    icode = 52_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
#ifdef __coupled
    CALL write_slice(IUNIT,AOFLTXWO)
#else
    CALL write_slice(IUNIT,TXO)
#endif /*__coupled*/
    !:: WRITE MERIDIONAL WIND STRESS
    icode = 53_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
#ifdef __coupled
    CALL write_slice(IUNIT,AOFLTYWE)
#else
    CALL write_slice(IUNIT,TYE)
#endif /*__coupled*/
    !:: WRITE SHORT WAVE RAD FLUX
    icode = 176_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,QSWO)
    !:: WRITE LONG WAVE RAD FLUX
    icode = 177_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,QLWO)
    !:: WRITE SENS HEAT FLUX
    icode = 146_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,QSEO)
    !:: WRITE LATENT HEAT FLUX
    icode = 147_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,QLAO)
    !:: WRITE STREAM FUNCTION
    icode = 27_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,PSITRO)
!!$    !:: WRITE MONTHLY MIXED LAYER DEPTH
!!$    icode = 83_C_INT64_T
!!$    DO L=1,12
!!$      IDATE2= (LYEARS*10000)+(L*100)+LDAYS
!!$      IF(p_pe==p_io) WRITE(IUNIT) IDATE2,ICODE,KCODE,IEDIM
!!$      CALL write_slice(IUNIT,AMLD(:,:,L))
!!$    ENDDO
#endif
    !:: WRITE DEPTO
    icode = 84_C_INT64_T
    kcode = 0_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,DEPTO)



#ifdef MIST
    !:: WRITE DLXP
    icode = 85_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,DLXP)
    !:: WRITE DLYP
    icode = 86_C_INT64_T
    IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
    CALL write_slice(IUNIT,DLYP)
    !:: WRITE WETO
    icode = 506_C_INT64_T
    DO K=1,KE
      KCODE = INT(TIESTU(K), C_INT64_T)
      IF(p_pe==p_io) WRITE(IUNIT) IDATE,ICODE,KCODE,IEDIM
      CALL write_slice(IUNIT,WETO(:,:,K))
    ENDDO
#endif
  END SUBROUTINE dump_restart_data_extra

END MODULE mo_restart
