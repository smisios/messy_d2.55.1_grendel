MODULE mo_mean

  USE mo_kind
  USE mo_param1
  USE mo_commo1

  IMPLICIT NONE

  REAL(wp), ALLOCATABLE :: sum_eminpo(:,:),sum_zo(:,:)         &
                      ,sum_flum(:,:),sum_pem(:,:),sum_sictho(:,:)      &
                      ,sum_sicomo(:,:),sum_sicuo(:,:),sum_sicve(:,:)   &
                      ,sum_sictru(:,:),sum_sictrv(:,:),sum_tho(:,:,:)  &
                      ,sum_sao(:,:,:)                                  &
                      ,sum_sicsno(:,:),sum_prech(:,:),sum_psitro(:,:)  &
                      ,sum_zmld(:,:),sum_rivrun(:,:)                   &
                      ,sum_qswo(:,:),sum_qlwo(:,:),sum_qlao(:,:)       &
                      ,sum_qseo(:,:),sum_wo(:,:,:)

  REAL(wp), ALLOCATABLE :: sum_ukomfl(:,:,:),sum_vkemfl(:,:,:)

! lconvdiag
  REAL(wp), ALLOCATABLE :: tmceo(:,:),tmcdo(:,:),tmepo(:,:)                &
                       ,sum_tmceo(:,:),sum_tmcdo(:,:)

! lforcediag
  REAL(wp), ALLOCATABLE :: sum_txo(:,:),sum_tye(:,:),sum_tafo(:,:)         &
                      ,sum_fswr(:,:),sum_fclou(:,:),sum_fprec(:,:)     &
                      ,sum_ftdew(:,:),sum_wu10(:,:)

! ldiffdiag
  REAL(wp), ALLOCATABLE :: sum_avo(:,:,:),sum_dvo(:,:,:)                    &
                      ,wtmix(:,:,:),sum_wtmix(:,:,:)                    &
                      ,rinu(:,:,:),sum_rinu(:,:,:)                      &
                      ,duvdz(:,:,:),sum_duvdz(:,:,:)                    &
                      ,drhdz(:,:,:),sum_drhdz(:,:,:)

! lgmdiag
  REAL(wp), ALLOCATABLE :: sum_wgo(:,:,:),sum_bolx(:,:),sum_boly(:,:)


! lhfldiag
  REAL(wp), ALLOCATABLE ::   dqswo(:,:),dqlwo(:,:)                           &
                        ,dqseo(:,:),dqlao(:,:), dqtho(:,:)               &
                        ,dqswi(:,:),dqlwi(:,:), dqsei(:,:),dqlai(:,:)    &
                        ,dqthi(:,:),sum_dqswo(:,:),sum_dqlwo(:,:)        &
                        ,sum_dqseo(:,:),sum_dqlao(:,:),sum_dqtho(:,:)    &
                        ,sum_dqswi(:,:),sum_dqlwi(:,:),sum_dqsei(:,:)    &
                        ,sum_dqlai(:,:),sum_dqthi(:,:)                   &
                        ,dticeo(:,:),sum_dticeo(:,:)


#ifdef __coupled
  REAL(wp), ALLOCATABLE ::   sum_aoflnhw(:,:),sum_aoflshw(:,:),sum_aoflchi(:,:) &
                        ,sum_aoflfrw(:,:),sum_aoflfri(:,:),sum_aofltxw(:,:) &
                        ,sum_aofltyw(:,:),sum_aofltxi(:,:),sum_aofltyi(:,:) &
                        ,sum_aoflwsv(:,:),sum_aoflrhi(:,:)
#endif /*__coupled*/


CONTAINS

  SUBROUTINE alloc_mem_mean

    ALLOCATE( sum_eminpo(ie,je),sum_zo(ie,je)                 &
             ,sum_flum(ie,je),sum_pem(ie,je),sum_sictho(ie,je)              &
             ,sum_sicomo(ie,je),sum_sicuo(ie,je),sum_sicve(ie,je)           &
             ,sum_sictru(ie,je),sum_sictrv(ie,je),sum_tho(ie,je,ke)         &
             ,sum_sao(ie,je,ke)        &
             ,sum_sicsno(ie,je),sum_prech(ie,je),sum_psitro(ie,je)          &
             ,sum_zmld(ie,je),sum_rivrun(ie,je)                             &
             ,sum_qswo(ie,je),sum_qlwo(ie,je),sum_qlao(ie,je)               &
             ,sum_qseo(ie,je),sum_wo(ie,je,ke)                              &
             ,sum_ukomfl(ie,je,ke),sum_vkemfl(ie,je,ke) )


    sum_eminpo(:,:) = 0.0_wp
    sum_zo(:,:) = 0.0_wp
    sum_flum(:,:) = 0.0_wp
    sum_pem(:,:) = 0.0_wp
    sum_sictho(:,:) = 0.0_wp
    sum_sicomo(:,:) = 0.0_wp
    sum_sicuo(:,:) = 0.0_wp
    sum_sicve(:,:) = 0.0_wp
    sum_sictru(:,:) = 0.0_wp
    sum_sictrv(:,:) = 0.0_wp
    sum_tho(:,:,:) = 0.0_wp
    sum_sao(:,:,:) = 0.0_wp
    sum_sicsno(:,:) = 0.0_wp
    sum_prech(:,:) = 0.0_wp
    sum_psitro(:,:) = 0.0_wp
    sum_zmld(:,:) = 0.0_wp
    sum_rivrun(:,:) = 0.0_wp
    sum_qswo(:,:) = 0.0_wp
    sum_qlwo(:,:) = 0.0_wp
    sum_qlao(:,:) = 0.0_wp
    sum_qseo(:,:) = 0.0_wp
    sum_wo(:,:,:) = 0.0_wp
    sum_ukomfl(:,:,:) = 0.0_wp
    sum_vkemfl(:,:,:) = 0.0_wp


    IF (lforcediag) THEN
       ALLOCATE( sum_txo(ie,je),sum_tye(ie,je),sum_tafo(ie,je)              &
            ,sum_fswr(ie,je),sum_fclou(ie,je),sum_fprec(ie,je)              &
            ,sum_ftdew(ie,je),sum_wu10(ie,je) )

       sum_txo(:,:) = 0.0_wp
       sum_tye(:,:) = 0.0_wp
       sum_tafo(:,:) = 0.0_wp
       sum_fswr(:,:) = 0.0_wp
       sum_fclou(:,:) = 0.0_wp
       sum_fprec(:,:) = 0.0_wp
       sum_ftdew(:,:) = 0.0_wp
       sum_wu10(:,:) = 0.0_wp

    ENDIF

    IF (ldiffdiag) THEN
       ALLOCATE( sum_avo(ie,je,ke),sum_dvo(ie,je,ke)                        &
            ,wtmix(ie,je,ke),sum_wtmix(ie,je,ke)                            &
            ,rinu(ie,je,ke),sum_rinu(ie,je,ke))

       sum_avo(:,:,:) = 0.0_wp
       sum_dvo(:,:,:) = 0.0_wp
       sum_wtmix(:,:,:) = 0.0_wp
       wtmix(:,:,:) = 0.0_wp
       sum_rinu(:,:,:) = 0.0_wp
       rinu(:,:,:) = 0.0_wp

    ENDIF

    IF (lgmdiag) THEN
       ALLOCATE( sum_wgo(ie,je,ke),sum_bolx(ie,je),sum_boly(ie,je) )

       sum_wgo(:,:,:) = 0.0_wp
       sum_bolx(:,:) = 0.0_wp
       sum_boly(:,:) = 0.0_wp

    ENDIF

    IF (lconvdiag) THEN
       ALLOCATE(tmepo(ie,je),tmceo(ie,je),tmcdo(ie,je)           &
                ,sum_tmceo(ie,je),sum_tmcdo(ie,je))

       tmepo(:,:) = 0.0_wp
       tmceo(:,:) = 0.0_wp
       tmcdo(:,:) = 0.0_wp

    ENDIF

    IF (lhfldiag) THEN
       ALLOCATE( dqswo(ie,je),dqlwo(ie,je),dqseo(ie,je),dqlao(ie,je)       &
            ,dqtho(ie,je),dqswi(ie,je),dqlwi(ie,je),dqsei(ie,je)           &
            ,dqlai(ie,je),dqthi(ie,je),sum_dqswo(ie,je),sum_dqlwo(ie,je)   &
            ,sum_dqseo(ie,je),sum_dqlao(ie,je),sum_dqtho(ie,je)            &
            ,sum_dqswi(ie,je),sum_dqlwi(ie,je),sum_dqsei(ie,je)            &
            ,sum_dqlai(ie,je),sum_dqthi(ie,je),dticeo(ie,je)               &
            ,sum_dticeo(ie,je) )
    ENDIF

#ifdef __coupled
    ALLOCATE( sum_aoflnhw(ie,je),sum_aoflshw(ie,je),sum_aoflchi(ie,je)      &
             ,sum_aoflfrw(ie,je),sum_aoflfri(ie,je),sum_aofltxw(ie,je)      &
             ,sum_aofltyw(ie,je),sum_aofltxi(ie,je),sum_aofltyi(ie,je)      &
             ,sum_aoflwsv(ie,je),sum_aoflrhi(ie,je) )

    sum_aoflnhw(:,:) = 0.0_wp
    sum_aoflshw(:,:) = 0.0_wp
    sum_aoflrhi(:,:) = 0.0_wp
    sum_aoflchi(:,:) = 0.0_wp
    sum_aoflfrw(:,:) = 0.0_wp
    sum_aoflfri(:,:) = 0.0_wp
    sum_aofltxw(:,:) = 0.0_wp
    sum_aofltyw(:,:) = 0.0_wp
    sum_aofltxi(:,:) = 0.0_wp
    sum_aofltyi(:,:) = 0.0_wp
    sum_aoflwsv(:,:) = 0.0_wp

#endif /*__coupled*/


  END SUBROUTINE alloc_mem_mean


SUBROUTINE mmean2d(kdays,kmonts,kyears,mnnndt,mnnn720,iun2d,                &
                   field2,sum_field2,i4code,ivec)
  !****************************************************************
  !
  !**** *mmean2d* - average data and save.
  !
  !     *mpi-met, hh*    10.04.01
  !
  !     modified
  !     --------
  !
  !     purpose
  !     -------
  !
  !
  !     method
  !     -------
  !
  !
  !**   interface.
  !     ----------
  !
  !     *call*       *mmean2d(kdays,kmonts,kyears,mnnndt,mnnn720,
  !                           iun2d,field2,sum_field2,i4code)
  !
  !     *parameter*  *param1.h*     - grid size parameters for ocean model.
  !     *common*     *commo1.h*     - ocean/sediment tracer arrays.
  !     *common*     *units.h*      - std i/o logical units.
  !
  !**   interface to calling routine (parameter list):
  !     ----------------------------------------------
  !
  !     *integer* *kdays*    - .
  !     *integer* *kmonts*    - .
  !     *integer* *kyears*    - .
  !     *integer* *mnnndt*    - actual number of day loop.
  !     *integer* *mnnn720*   - end of the avg period
  !     *integer* *iun2d*    - .
  !     *real*    *field2*   - .
  !     *real*    *sum_field2*   - .
  !     *integer* *i4code*    - .
  !
  !
  !     externals
  !     ---------
  !     none.
  !
  !**************************************************************************

  USE mo_param1
  USE mo_mpi
  USE mo_parallel
  USE mo_commo1
  USE mo_units
  USE mo_kind

  INTEGER(kind=i4) idate,icode,klev,iedim
  INTEGER awet,i,j,ivec,kdays,kmonts,kyears,mnnndt,mnnn720,iun2d,i4code
  REAL(wp) field2(ie,je),sum_field2(ie,je)
  REAL(sp) , ALLOCATABLE :: ff_g(:,:)
  REAL(wp), ALLOCATABLE :: sum_field2_g(:,:)

  ! set extra file format parameters

  icode=i4code
  iedim=ie_g*je_g

  !h at beginning of averaging period : initialize the field

  IF (mnnndt.EQ.1) THEN
     DO i=1,ie
        DO j=1,je
           sum_field2(i,j) = 0._wp
        ENDDO
     ENDDO
  ENDIF

  !h summation and averaging of the field; check out of range for ieee

  DO i=1,ie
     DO j=1,je
        sum_field2(i,j) = sum_field2(i,j) &
             + (field2(i,j)-sum_field2(i,j))/REAL(mnnndt, wp)
     ENDDO
  ENDDO

  IF (mnnndt.EQ.mnnn720) THEN

    ! FIXME: ivec is loop-invariant and the sequence-insensitive: the
    ! loop can be written much more succinct with WHERE,
    ! also why not use -HUGE instead of -9e33?
     DO i=1,ie
        DO j=1,je
           IF (ivec.EQ.-1) awet= 1
           IF (ivec.EQ.0) awet= INT(weto(i,j,1))
           IF (ivec.EQ.1) awet= INT(amsuo(i,j,1))
           IF (ivec.EQ.2) awet= INT(amsue(i,j,1))
           IF (awet .LT. 1) THEN
              sum_field2(i,j) = -9.e33_wp
           ENDIF
!           if (abs(sum_field2(i,j)).lt.(1.e-35)) then
!              sum_field2(i,j)=0.
!           endif
        ENDDO
     ENDDO

     IF (p_pe == p_io) THEN
        ALLOCATE(sum_field2_g(ie_g,je_g), ff_g(ie_g,je_g))
     ELSE
        ALLOCATE(sum_field2_g(0,0), ff_g(0,0))
     ENDIF

     CALL gather(sum_field2,sum_field2_g,p_io)

     IF (p_pe == p_io) THEN

        !h write in extra format
        WRITE(io_stdout,*)'open unit: ',iun2d
!        open(iun2d, status='unknown',                                  &
!             &       access='sequential',                                       &
!             &       position='append',                                         &
!             &       form='unformatted')
        idate= (kyears*10000)+(kmonts*100)+kdays

        WRITE(io_stdout,*)'mmean2d: unit=',iun2d
        WRITE(io_stdout,'(1x,a11,i8)') 'yyyymmdd=',idate
        WRITE(io_stdout,*)'         code=',icode
        WRITE(io_stdout,*)'    nanf/nend=',mnnndt,mnnn720

        DO j=1,je_g
           DO i=1,ie_g
              ff_g(i,j)=REAL(sum_field2_g(i,j),sp)
           ENDDO
        ENDDO
        klev=0
        WRITE(iun2d) idate,icode,klev,iedim
        WRITE(iun2d) ff_g
!        write(iun2d)((ff_g(i,j),i=1,ie_g),j=1,je_g)
!        close(iun2d)

     ENDIF  ! p_pe==p_io

     DEALLOCATE(sum_field2_g,ff_g)

  ENDIF

END SUBROUTINE mmean2d

SUBROUTINE mmean3d(kdays,kmonts,kyears,mnnndt,mnnn720,iun3d,      &
                   field,sum_field,i4code,ivec)
  !
  !**** *mmean3d* - save mean diagnostic output.
  !
  !     chh,    *mpi-met, hh*   14.01.99
  !
  !     modified
  !     --------
  !     s.legutke,        *mpi-mad, hh*    01.10.01
  !     - separate routine extracted from ollie (main)
  !
  !     purpose
  !     -------
  !     accumulate 3d fields, average, and save.
  !
  !     method
  !     -------
  !     field is set to 0 at the first time step of each month,
  !     accumulated each step, and normalized at the end of the
  !     averaging time period. the filed is written to disk in
  !     extra format (default, or in netcdf).
  !
  !**   interface.
  !     ----------
  !
  !     *call*       *mmean3d(kdays,kmonts,kyears,mnnndt,mnnn720,iun3d,
  !                        field,sum_field,i4code)*
  !
  !     *parameter*  *mo_param1*     - grid size parameters for ocean model.
  !     *common*     *mo_commo1*     - ocean/sediment tracer arrays.
  !     *common*     *mo_units*      - std i/o logical units.
  !
  !**   interface to calling routine (parameter list):
  !     ----------------------------------------------
  !
  !     *integer* *kyears   - actual year.
  !     *integer* *kmonts*   - actual month.
  !     *integer* *kdays    - actual day.
  !
  !
  !     externals
  !     ---------
  !     none.
  !
  USE mo_param1
  USE mo_parallel
  USE mo_commo1
  USE mo_units
  USE mo_kind


  REAL(wp) field(:,:,:)
  REAL(wp) sum_field(:,:,:)

  INTEGER(kind=i4) idate,klev,icode,iedim
  INTEGER awet,i,j,k,ivec,kdays,kmonts,kyears,mnnndt,mnnn720,iun3d,i4code,kk

  REAL(kind=sp),ALLOCATABLE ::  ff_g(:,:)
  REAL(wp), ALLOCATABLE :: sum_field_g(:,:,:)

!  kk=ubound(field,3)
  kk=ke

  iedim=ie_g*je_g
  icode=i4code

  !h    initialize the field

  IF (mnnndt.EQ.1) THEN
     DO k=1,kk
        DO i=1,ie
           DO j=1,je
              sum_field(i,j,k) = 0._wp
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  !h     summation and monthly mean of the field

  DO k=1,kk
     DO j=1,je
        DO i=1,ie
           sum_field(i,j,k) = sum_field(i,j,k) &
                + (field(i,j,k)-sum_field(i,j,k))/REAL(mnnndt, wp)
        ENDDO
     ENDDO
  ENDDO

  IF (mnnndt.EQ.mnnn720) THEN


    ! FIXME: ivec is loop-invariant and the sequence-insensitive: the
    ! loop can be written much more succinct with WHERE,
    ! also why not use -HUGE instead of -9e33?
  DO k=1,kk
    DO j=1,je
       DO i=1,ie
           IF (ivec.EQ.-1) awet= 1
           IF (ivec.EQ.0) awet= INT(weto(i,j,k))
           IF (ivec.EQ.1) awet= INT(amsuo(i,j,k))
           IF (ivec.EQ.2) awet= INT(amsue(i,j,k))
           IF (awet .LT. 1) THEN
              sum_field(i,j,k) = -9.e33_wp
           ENDIF
!           if (abs(sum_field(i,j,k)).lt.(1.e-35)) then
!              sum_field(i,j,k)=0.
!           endif
        ENDDO
     ENDDO
  ENDDO

     IF (p_pe == p_io) THEN
        ALLOCATE(sum_field_g(ie_g,je_g,kk), ff_g(ie_g,je_g))
     ELSE
        ALLOCATE(sum_field_g(0,0,0), ff_g(0,0))
     ENDIF

     CALL gather(sum_field,sum_field_g,p_io)

     idate= (kyears*10000)+(kmonts*100)+kdays
     WRITE(io_stdout,*)'mmean3d: unit=',iun3d
     WRITE(io_stdout,'(1x,a11,i8)') '  yyyymmdd=',idate
     WRITE(io_stdout,*)'         code=',icode
     WRITE(io_stdout,*)'    nanf/nend=',mnnndt,mnnn720


     IF (p_pe==p_io) THEN
!        open(iun3d,    status='unknown',                            &
!                                 access='sequential',               &
!                                 position='append',                 &
!                                form='unformatted')

        DO k=1,kk
           DO j=1,je_g
              DO i=1,ie_g
                 ff_g(i,j)=REAL(sum_field_g(i,j,k),sp)
              ENDDO
           ENDDO

           IF (icode.EQ.7.OR.icode.EQ.110.OR.icode.EQ.111.OR.icode.EQ.110) THEN
              klev=NINT(tiestw(k))
           ELSE
              klev=NINT(tiestu(k))
           ENDIF
           WRITE(iun3d) idate,icode,klev,iedim
!           write(iun3d)((ff_g(i,j),i=1,ie_g),j=1,je_g)
           WRITE(iun3d) ff_g
        ENDDO
!        close(iun3d)
     ENDIF ! p_pe==p_io

     DEALLOCATE(sum_field_g,ff_g)

  ENDIF


END SUBROUTINE mmean3d

SUBROUTINE calc_avgint(km,na,ne)

  USE mo_commo1, ONLY : ldtdayc,ldays,lmonts,icontro,ndtday
  USE mo_units , ONLY : io_stdout
  USE mo_model_time, ONLY : monlen

  INTEGER, INTENT(in) :: km
  INTEGER, INTENT(inout) :: na
  INTEGER, INTENT(out) :: ne
  INTEGER :: i,nc

! set 'actual' time step and output time step for cases kmean=1/2/3:
!
!hh daily average
      IF (km.EQ.1) THEN

        na=ldtdayc
        ne=ndtday

!hh monthly average
      ELSEIF (km.EQ.2) THEN

        IF (ldays.EQ.1) THEN
          na=ldtdayc
        ELSEIF (ldays.GT.1) THEN
          na=ldtdayc+((ldays-1)*ndtday)
        ENDIF
        ne=ndtday*monlen(lmonts, lyears)

!hh yearly average
      ELSEIF (km.EQ.3) THEN

        IF ((ldays.EQ.1).AND.(lmonts.EQ.1)) THEN
          na=ldtdayc
        ELSE
          IF (lmonts.EQ.1) THEN
            na=ldtdayc+((ldays-1)*ndtday)
          ELSE
            nc=0
            DO i=1,lmonts-1
              nc=nc+ndtday*monlen(i, lyears)
            ENDDO
            na=ldtdayc+((ldays-1)*ndtday)+nc
          ENDIF
        ENDIF
        ne=0
        DO i=1,12
          ne=ne+monlen(i, lyears)
        ENDDO
        ne=ne*ndtday

      ELSEIF (km.EQ.4) THEN

!hh every timestep
        na=1
        ne=1

!hh no diagnostic output
       ELSEIF (km.EQ.0.OR.km.EQ.99) THEN
         CONTINUE
       ELSE
         STOP 'stop in calc_avgint due to wrong parameter.'
      ENDIF

      IF (icontro.NE.0)THEN
         IF (km.NE.0) THEN
            WRITE(io_stdout,*) 'na,ne',na,ne
         ENDIF
      ENDIF

END SUBROUTINE calc_avgint

END MODULE mo_mean








