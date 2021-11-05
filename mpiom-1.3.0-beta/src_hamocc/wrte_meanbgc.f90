      SUBROUTINE WRTE_MEANBGC

!****************************************************************
!     save bgc output
!**********************************************************************

      USE mo_param1, ONLY: ie,je,ke
      USE MO_COMMO1, ONLY: ldays,lmonts,ldtdayc,dt,ddpo,monlen
      USE mo_bgcmean, ONLY: mean_2d_freq, mean_3d_freq,       &
      meantime_2d, meantime_3d

      IMPLICIT NONE

      INTEGER :: NDTDAY,nanf2d,nend2d,nanf3d,nend3d,i

      NDTDAY=NINT(86400./DT)


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
        nend2d=ndtday*monlen(lmonts)
! yearly average
      ELSEIF (mean_2d_freq.EQ.3) THEN
        IF ((ldays.EQ.1).AND.(lmonts.EQ.1)) THEN
          nanf2d=ldtdayc
        ELSE
          nanf2d=nanf2d+1
        ENDIF
        nend2d=0
        DO i=1,12
          nend2d=nend2d+monlen(i)
        ENDDO
        nend2d=nend2d*ndtday
      ELSEIF (mean_2d_freq.EQ.4) THEN

! every timestep   
        nanf2d=ldtdayc
        nend2d=ldtdayc

! no diagnostic output
       ELSEIF (mean_2d_freq.EQ.0) THEN
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
        nend3d=ndtday*monlen(lmonts)
! yearly average
      ELSEIF (mean_3d_freq.EQ.3) THEN
        IF ((ldays.EQ.1).AND.(lmonts.EQ.1)) THEN
          nanf3d=ldtdayc
        ELSE
          nanf3d=nanf3d+1
        ENDIF
        nend3d=0
        DO i=1,12
          nend3d=nend3d+monlen(i)
        ENDDO
        nend3d=nend3d*ndtday
      ELSEIF (mean_3d_freq.EQ.4) THEN

! every timestep   
        nanf3d=ldtdayc
        nend3d=ldtdayc

! no diagnostic output
       ELSEIF (mean_3d_freq.EQ.0) THEN
         CONTINUE
       ELSE
         STOP 'stop in wrte_bgcmean due to wrong parameter.'
      ENDIF

 
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






