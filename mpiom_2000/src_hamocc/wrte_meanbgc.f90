      SUBROUTINE WRTE_MEANBGC
!****************************************************************
!     save bgc output
!**********************************************************************

      USE mo_param1, ONLY: ie,je,ke
      USE mo_commo1, ONLY: ldays,lmonts,lyears,ldtdayc,ddpo,ndtday
      USE mo_model_time, ONLY: monlen
      USE mo_bgcmean, ONLY: mean_2d_freq, mean_3d_freq,       &
      meantime_2d, meantime_3d

      IMPLICIT NONE

      INTEGER :: nanf2d,nend2d,nanf3d,nend3d,i

      nanf2d=0
      nanf3d=0
      nend2d=-1
      nend3d=-1

! for 2d output
! daily average
      IF (mean_2d_freq.EQ.1) THEN
        nanf2d=ldtdayc
        nend2d=ndtday
! monthly average
     ELSEIF (mean_2d_freq.EQ.2) THEN
        IF (ldays.EQ.1) THEN
          nanf2d=ldtdayc
        ELSEIF (ldays.GT.1) THEN
          nanf2d=ldtdayc+((ldays-1)*ndtday)
        ENDIF
        nend2d=ndtday*monlen(lmonts, lyears)
! yearly average
      ELSEIF (mean_2d_freq.EQ.3) THEN
        IF ((ldays.EQ.1).AND.(lmonts.EQ.1)) THEN
          nanf2d=ldtdayc
         ELSEIF ((ldays.GT.1).AND.(lmonts.EQ.1)) THEN
          nanf2d=ldtdayc+((ldays-1)*ndtday)
         ELSEIF (lmonts.GT.1) THEN
            nanf2d=0
            DO i=1,lmonts-1
               nanf2d=nanf2d+monlen(i, lyears)
            ENDDO
            nanf2d=nanf2d*ndtday
            nanf2d=nanf2d+(ldtdayc+((ldays-1)*ndtday))
        ENDIF
        nend2d=0
        DO i=1,12
          nend2d=nend2d+monlen(i, lyears)
        ENDDO
        nend2d=nend2d*ndtday
      ELSEIF (mean_2d_freq.EQ.4) THEN

! every timestep
        nanf2d=1
        nend2d=1

! no diagnostic output
       ELSEIF (mean_2d_freq.EQ.0) THEN
      nanf2d=0
      nend2d=-1
         CONTINUE
       ELSE
         STOP 'stop in wrte_bgcmean due to wrong parameter.'
      ENDIF



! for 3d output
  ! daily average
      IF (mean_3d_freq.EQ.1) THEN
        nanf3d=ldtdayc
        nend3d=ndtday
! monthly average
     ELSEIF (mean_3d_freq.EQ.2) THEN
        IF (ldays.EQ.1) THEN
          nanf3d=ldtdayc
        ELSEIF (ldays.GT.1) THEN
          nanf3d=ldtdayc+((ldays-1)*ndtday)
        ENDIF
        nend3d=ndtday*monlen(lmonts, lyears)
! yearly average
      ELSEIF (mean_3d_freq.EQ.3) THEN
        IF ((ldays.EQ.1).AND.(lmonts.EQ.1)) THEN
          nanf3d=ldtdayc
         ELSEIF ((ldays.GT.1).AND.(lmonts.EQ.1)) THEN
          nanf3d=ldtdayc+((ldays-1)*ndtday)
         ELSEIF (lmonts.GT.1) THEN
            nanf3d=0
            DO i=1,lmonts-1
               nanf3d=nanf3d+monlen(i, lyears)
            ENDDO
            nanf3d=nanf3d*ndtday
            nanf3d=nanf3d+(ldtdayc+((ldays-1)*ndtday))
        ENDIF
        nend3d=0
        DO i=1,12
          nend3d=nend3d+monlen(i, lyears)
        ENDDO
        nend3d=nend3d*ndtday
      ELSEIF (mean_3d_freq.EQ.4) THEN

! every timestep
        nanf3d=1
        nend3d=1

! no diagnostic output
       ELSEIF (mean_3d_freq.EQ.0) THEN
         CONTINUE
      nanf3d=0
      nend3d=-1

       ELSE
         STOP 'stop in wrte_bgcmean due to wrong parameter.'
      ENDIF


!      WRITE(0,*)'meanbgc: ',nanf2d,nend2d
!      WRITE(0,*)'MEAN_BGC',ldays,lmonts,ldtdayc,nanf2d,nend2d

      IF (nanf2d.EQ.nend2d) THEN
         meantime_2d = meantime_2d + 1
         CALL avrg_bgcmean_2d(ie,je,ke)
         CALL write_bgcmean_2d(ie,je,ke,ddpo)
         CALL write_bgcmean_bioz(ie,je,ke,ddpo)

      ENDIF

      IF (nanf3d.EQ.nend3d) THEN
         meantime_3d = meantime_3d + 1
         CALL avrg_bgcmean_3d(ie,je,ke)
         CALL write_bgcmean_3d(ie,je,ke,ddpo)
         CALL write_bgcmean_sed(ie,je,ke,ddpo)

      ENDIF

      END






