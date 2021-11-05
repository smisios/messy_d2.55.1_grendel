!>
 !! Module for ocean river runoff. This includes the following subroutines :
 !!
 !! river_runoff_omip_ini : river runoff initialisation (default for omip forcing)
 !! river_runoff_stations_ini : river runoff initialisation
 !! glac_calv_ini : glacier calving initialisation
 !! river_runoff_stations : calculation of river runoff
 !! river_runoff_omip  : calculation of river runoff (default for omip forcing)
 !! glac_calv : : glacier calving initialisation glacier calving
 !!
 !! @author: Helmuth Haak, Max Planck Institute for Meteorology
 !!

MODULE mo_runoff

  USE mo_kind, ONLY: i8, dp, wp
  USE mo_param1, ONLY : ie,je
  USE mo_parallel, ONLY: p_pe,p_io,p_ioff,p_joff
  USE mo_commo1, ONLY: zo,sictho,sicsno,alat_g,alon_g,wetol1_g &
       ,ddpo, dlxp, dlyp, sao, tho, almzer, dt  &
       ,giriv,rivrun,ldays,lmonts,lyears
  USE mo_model_time, ONLY: monlen
  USE mo_commoau1, ONLY : tfreeze, entmel
  USE mo_planetary_constants, ONLY : rhoicwa,rhosnwa
  USE mo_commoau2, ONLY : preco, prech
  USE mo_io_config, ONLY: next_free_unit

  USE mo_units, ONLY: io_stdout

  IMPLICIT NONE

  LOGICAL :: luse_river_runoff_stations
  INTEGER :: numriv ! number of runoff stations (namelist parameter)

  LOGICAL :: luse_glac_calv
  INTEGER :: numglac ! number of glacier calving stations (namelist parameter)

  ! FIXME: dp because the values are READ from a file, one could consider
  ! conversion to wp after reading here
  REAL(dp), ALLOCATABLE :: RIVAL(:,:),RIVLON(:),RIVLAT(:)
  REAL(wp), ALLOCATABLE :: FRIV(:), DDRIV(:)
  INTEGER, ALLOCATABLE :: IRIVI(:),IRIVJ(:)

  REAL(wp), ALLOCATABLE :: GLACLAT(:),GLACLON(:),GLACVAL(:)
  INTEGER, ALLOCATABLE :: IGLAC(:),JGLAC(:)


CONTAINS

  SUBROUTINE river_runoff_omip_ini

    ALLOCATE(GIRIV(IE,JE))
    giriv(:,:) = 0._wp

  END SUBROUTINE river_runoff_omip_ini

  SUBROUTINE river_runoff_stations_ini

!  initialisation of river input locations

    USE mo_parallel, ONLY : p_bcast

    USE mo_grid, ONLY: p_suchij

    INTEGER :: n,m

    REAL(wp) :: dist

    INTEGER :: io_in_rval, io_in_rpos
    INTEGER(i8) :: ibla(4)

    ALLOCATE(RIVAL(NUMRIV,12),FRIV(NUMRIV))
    ALLOCATE(RIVLON(NUMRIV),RIVLAT(NUMRIV),DDRIV(NUMRIV))
    ALLOCATE(IRIVI(NUMRIV),IRIVJ(NUMRIV))

    FRIV(:) = 0._wp
    RIVLON(:) = 0._wp
    RIVLAT(:) = 0._wp
    DDRIV(:) = 0._wp
    IRIVI(:) = 0
    IRIVJ(:) = 0
    RIVAL(:,:) = 0._wp


!hh-------------------------------------------------------------------
!hh   river runoff data from ieee extra file: observations 12mon.clim.

    IF(p_pe==p_io) THEN
      io_in_rval = next_free_unit()
      OPEN(io_in_rval, file='runoff_obs', form='unformatted', action='read')
      io_in_rpos = next_free_unit()
      OPEN(io_in_rpos, file='runoff_pos', form='unformatted', action='read')
      DO n=1,numriv
        READ(io_in_rval)ibla
        READ(io_in_rval)(rival(n,m),m=1,12)
        READ(io_in_rpos)ibla
        READ(io_in_rpos)rivlat(n),rivlon(n)
      ENDDO
      CLOSE(io_in_rval)
      CLOSE(io_in_rpos)
    ENDIF
    CALL p_bcast(rival,p_io)
    CALL p_bcast(rivlat,p_io)
    CALL p_bcast(rivlon,p_io)


    DO n=1,numriv
      CALL p_suchij(rivlat(n), rivlon(n), 1, irivi(n), irivj(n), dist, 1._wp)

      IF(dist .GT. 1.e6_wp)THEN  ! if distance is larger than 1e6m
                            ! the river mouth found is  ignored
        irivi(n)=1
        irivj(n)=1
      ENDIF

      IF(p_pe==p_io) THEN
        WRITE(io_stdout,*)'river nr: ',rivlat(n),rivlon(n),n              &
             ,irivi(n),irivj(n)                         &
             ,dist,                                     &
             wetol1_g(irivi(n),irivj(n)),               &
             alon_g(irivi(n),irivj(n)),                 &
             alat_g(irivi(n),irivj(n))
      ENDIF



    ENDDO


  END SUBROUTINE river_runoff_stations_ini

  SUBROUTINE glac_calv_ini

!  initialisation of glacier calcing locations

    USE mo_parallel, ONLY : p_bcast
    USE mo_grid, ONLY: p_suchij


    INTEGER :: n, io_in_glac

    REAL(wp) :: dist

    ALLOCATE(GLACLAT(NUMGLAC),GLACLON(NUMGLAC),GLACVAL(NUMGLAC))
    ALLOCATE(IGLAC(NUMGLAC),JGLAC(NUMGLAC))

    IF(p_pe==p_io) THEN
      io_in_glac = next_free_unit()
      OPEN(io_in_glac, file='gletscher_5653', &
           form='unformatted', action='read')
      DO n=1,numglac
        READ(io_in_glac)glaclat(n),glaclon(n),glacval(n)
      ENDDO
      CLOSE (io_in_glac)
    ENDIF
    CALL p_bcast(glaclat,p_io)
    CALL p_bcast(glaclon,p_io)
    CALL p_bcast(glacval,p_io)

    DO n=1,numglac

      CALL p_suchij(glaclat(n), glaclon(n), 1, iglac(n), jglac(n), dist, 1._wp)

      IF(p_pe==p_io) THEN
        WRITE(io_stdout,*)'glacier nr: ',glaclat(n),glaclon(n),n          &
             ,iglac(n),jglac(n)                         &
             ,dist,                                     &
             wetol1_g(iglac(n),jglac(n)),               &
             alon_g(iglac(n),jglac(n)),                 &
             alat_g(iglac(n),jglac(n))
      ENDIF
    ENDDO


  END SUBROUTINE glac_calv_ini




  SUBROUTINE river_runoff_stations


    INTEGER i,j
    REAL(wp) zzzdz
    INTEGER :: monmon,n
    REAL(wp) :: awert, ewert

    rivrun(:,:) = 0._wp                   !  runoff diagnostic

!
!  UPDATE FORCING ONCE PER DAY

       IF(LDAYS.LE.(MONLEN(LMONTS, lyears)+1)/2)THEN
        MONMON=LMONTS-1
        IF(MONMON.LT.1)MONMON=12
        awert = 0.5_wp &
             + (REAL(ldays, wp) - 0.5_wp) / REAL(monlen(lmonts, lyears), wp)
       ELSE
        MONMON=LMONTS+1
        IF(MONMON.GT.1)MONMON=1
        awert = 0.5_wp &
             + (REAL(monlen(lmonts, lyears) - ldays, wp) + 0.5_wp) &
             & / REAL(monlen(lmonts, lyears), wp)
       ENDIF
        ewert = 1._wp - awert

       DO N=1,NUMRIV
         FRIV(N)=(AWERT*RIVAL(N,LMONTS)+EWERT*RIVAL(N,MONMON))
       ENDDO

!$omp parallel private(i,j,n,zzzdz)
!$omp do
    DO n=1,numriv
       i=irivi(n)-p_ioff
       j=irivj(n)-p_joff
       IF(i>=1 .AND. i<=ie .AND. j>=1 .AND. j<=je) THEN
          zzzdz=MAX(almzer,ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa      &
               -sicsno(i,j)*rhosnwa)
!         if(weto(i,j,1).lt.0.5)write(io_stdout,*)'alarm! river ',n,i,j
          ddriv(n)=friv(n)*dt/(dlxp(i,j)*dlyp(i,j))
!uwe      use actual layerthickness for mass/salt conservation
          sao(i,j,1)=sao(i,j,1)*zzzdz/(zzzdz+ddriv(n))
          zo(i,j)=zo(i,j)+ddriv(n)
          prech(i,j)=prech(i,j)+ddriv(n)/dt
          preco(i,j)=preco(i,j)+ddriv(n)/dt
          rivrun(i,j)=ddriv(n)/dt
!        write(io_stdout,*)' river ',n,ddriv(n),zo(i,j),sao(i,j,1)
       ENDIF
    ENDDO
!$omp end do
!$omp end parallel


  END SUBROUTINE river_runoff_stations

  SUBROUTINE river_runoff_omip

    INTEGER i,j
    REAL(wp) driv,zzzdz

    rivrun(:,:) = 0._wp                   !  runoff diagnostic


!$omp parallel private(i,j,driv,zzzdz)
!$omp do
    DO j=1,je
       DO i=1,ie
          zzzdz=MAX(almzer,ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa        &
               -sicsno(i,j)*rhosnwa)
          !       if(weto(i,j,1).lt.0.5.and.giriv(i,j).gt.0.)then
          !           write(io_stdout,*)'alarm! river ',n,i,j
          !       endif
          driv=giriv(i,j)*dt/(dlxp(i,j)*dlyp(i,j))
          !uwe      use actual layerthickness for mass/salt conservation
          sao(i,j,1)=sao(i,j,1)*zzzdz/(zzzdz+driv)
          zo(i,j)=zo(i,j)+driv
          preco(i,j)=preco(i,j)+driv/dt
          prech(i,j)=prech(i,j)+driv/dt
          rivrun(i,j)=driv/dt
       ENDDO
    ENDDO
!$omp end do
!$omp end parallel


  END SUBROUTINE river_runoff_omip


  SUBROUTINE glac_calv

    USE mo_planetary_constants, ONLY: rocp

    INTEGER :: i,j,n
    REAL(wp) :: thick, wert, tempfac

!     glacier calving
!
    tempfac = tfreeze - entmel/rocp

    !     assumption glacier calving as very cold water to conserve heat

    DO n=1,numglac
      i=iglac(n)-p_ioff
      j=jglac(n)-p_joff
      IF(i>=1 .AND. i<=ie .AND. j>=1 .AND. j<=je) THEN
        wert=glacval(n)/(dlxp(i,j)*dlyp(i,j))

        !        ice compactness minimum glacinput*0.5

        thick=ddpo(i,j,1)+zo(i,j)-sictho(i,j)*rhoicwa                  &
             &       +sicsno(i,j)*rhosnwa
        sao(i,j,1)=sao(i,j,1)*thick/(thick+wert*dt)
        tho(i,j,1)=(tho(i,j,1)*thick+wert*dt*tempfac)                  &
     &              /(thick+wert*dt)
        zo(i,j)=zo(i,j)+wert*dt
      ENDIF
    ENDDO

    !      update p-e field for diagnostics
    DO n=1,numglac
      i=iglac(n)-p_ioff
      j=jglac(n)-p_joff
      IF(i>=1 .AND. i<=ie .AND. j>=1 .AND. j<=je) THEN
         preco(i,j)=preco(i,j)+glacval(n)/(dlxp(i,j)*dlyp(i,j))
         prech(i,j)=prech(i,j)+glacval(n)/(dlxp(i,j)*dlyp(i,j))
       ENDIF
     ENDDO

   END SUBROUTINE glac_calv


END MODULE mo_runoff
