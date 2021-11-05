MODULE mo_mean

  USE mo_kind
  USE mo_param1
  USE mo_commo1

  IMPLICIT NONE

  REAL, ALLOCATABLE :: flum(:,:),pem(:,:),sictru(:,:)                  &
                      ,sictrv(:,:),sum_eminpo(:,:),sum_zo(:,:)         &
                      ,sum_flum(:,:),sum_pem(:,:),sum_sictho(:,:)      &
                      ,sum_sicomo(:,:),sum_sicuo(:,:),sum_sicve(:,:)   &
                      ,sum_sictru(:,:),sum_sictrv(:,:),sum_tho(:,:,:)  &
                      ,sum_sao(:,:,:)                                  &
                      ,sum_sicsno(:,:),sum_prech(:,:),sum_psiuwe(:,:)  &
                      ,sum_zmld(:,:),sum_rivrun(:,:)                   &
                      ,sum_qswo(:,:),sum_qlwo(:,:),sum_qlao(:,:)       &
                      ,sum_qseo(:,:),sum_wo(:,:,:)                     

  REAL, ALLOCATABLE :: sum_ukomfl(:,:,:),sum_vkemfl(:,:,:)

! LCONVDIAG
  REAL, ALLOCATABLE :: TMCEO(:,:),TMCDO(:,:)                           &
                       ,SUM_TMCEO(:,:),SUM_TMCDO(:,:)

! LFORCEDIAG
  REAL, ALLOCATABLE :: sum_txo(:,:),sum_tye(:,:),sum_tafo(:,:)         &
                      ,sum_fswr(:,:),sum_fclou(:,:),sum_fprec(:,:)     &
                      ,sum_ftdew(:,:),sum_wu10(:,:) 

! LDIFFDIAG
  REAL, ALLOCATABLE :: sum_avo(:,:,:),sum_dvo(:,:,:)                    &
                      ,wtmix(:,:,:),sum_wtmix(:,:,:)                    &
                      ,rinu(:,:,:),sum_rinu(:,:,:)                      &
                      ,duvdz(:,:,:),sum_duvdz(:,:,:)                    &
                      ,drhdz(:,:,:),sum_drhdz(:,:,:)                    

! LGMDIAG      
  REAL, ALLOCATABLE :: sum_wgo(:,:,:),sum_bolx(:,:),sum_boly(:,:)       


! LHFLDIAG
  REAL, ALLOCATABLE ::   dqswo(:,:),dqlwo(:,:)                           &
                        ,dqseo(:,:),dqlao(:,:), dqtho(:,:)               &
                        ,dqswi(:,:),dqlwi(:,:), dqsei(:,:),dqlai(:,:)    &
                        ,dqthi(:,:),sum_dqswo(:,:),sum_dqlwo(:,:)        &
                        ,sum_dqseo(:,:),sum_dqlao(:,:),sum_dqtho(:,:)    &
                        ,sum_dqswi(:,:),sum_dqlwi(:,:),sum_dqsei(:,:)    &
                        ,sum_dqlai(:,:),sum_dqthi(:,:)                   &
                        ,dticeo(:,:),sum_dticeo(:,:)


#ifdef __coupled
  REAL, ALLOCATABLE ::   sum_aoflnhw(:,:),sum_aoflshw(:,:),sum_aoflchi(:,:) &
                        ,sum_aoflfrw(:,:),sum_aoflfri(:,:),sum_aofltxw(:,:) &
                        ,sum_aofltyw(:,:),sum_aofltxi(:,:),sum_aofltyi(:,:) &
                        ,sum_aoflwsv(:,:),sum_aoflrhi(:,:)
#endif /*__coupled*/


CONTAINS

  SUBROUTINE alloc_mem_mean

    ALLOCATE( flum(ie,je),pem(ie,je),sictru(ie,je)                          &
             ,sictrv(ie,je),sum_eminpo(ie,je),sum_zo(ie,je)                 &
             ,sum_flum(ie,je),sum_pem(ie,je),sum_sictho(ie,je)              &
             ,sum_sicomo(ie,je),sum_sicuo(ie,je),sum_sicve(ie,je)           &
             ,sum_sictru(ie,je),sum_sictrv(ie,je),sum_tho(ie,je,ke)         &
             ,sum_sao(ie,je,ke)        &
             ,sum_sicsno(ie,je),sum_prech(ie,je),sum_psiuwe(ie,je)          &
             ,sum_zmld(ie,je),sum_rivrun(ie,je)                             &
             ,sum_qswo(ie,je),sum_qlwo(ie,je),sum_qlao(ie,je)               &
             ,sum_qseo(ie,je),sum_wo(ie,je,ke)                              &
             ,sum_ukomfl(ie,je,ke),sum_vkemfl(ie,je,ke) )
 
    flum(:,:)=0.0
    pem(:,:)=0.0
    sictru(:,:)=0.0
    sictrv(:,:)=0.0
    sum_eminpo(:,:)=0.0
    sum_zo(:,:)=0.0
    sum_flum(:,:)=0.0
    sum_pem(:,:)=0.0
    sum_sictho(:,:)=0.0
    sum_sicomo(:,:)=0.0
    sum_sicuo(:,:)=0.0
    sum_sicve(:,:)=0.0
    sum_sictru(:,:)=0.0
    sum_sictrv(:,:)=0.0
    sum_tho(:,:,:)=0.0
    sum_sao(:,:,:)=0.0
    sum_sicsno(:,:)=0.0
    sum_prech(:,:)=0.0
    sum_psiuwe(:,:)=0.0
    sum_zmld(:,:)=0.0
    sum_rivrun(:,:)=0.0 
    sum_qswo(:,:)=0.0
    sum_qlwo(:,:)=0.0
    sum_qlao(:,:)=0.0  
    sum_qseo(:,:)=0.0
    sum_wo(:,:,:)=0.0
    sum_ukomfl(:,:,:)=0.0
    sum_vkemfl(:,:,:)=0.0


    if (LFORCEDIAG) then
       ALLOCATE( sum_txo(ie,je),sum_tye(ie,je),sum_tafo(ie,je)              &
            ,sum_fswr(ie,je),sum_fclou(ie,je),sum_fprec(ie,je)              &
            ,sum_ftdew(ie,je),sum_wu10(ie,je) )                             

       SUM_TXO(:,:)=0.0
       SUM_TYE(:,:)=0.0
       SUM_TAFO(:,:)=0.0
       SUM_FSWR(:,:)=0.0
       SUM_FCLOU(:,:)=0.0
       SUM_FPREC(:,:)=0.0
       SUM_FTDEW(:,:)=0.0
       SUM_WU10(:,:)=0.0

    endif

    if (LDIFFDIAG) then
       ALLOCATE( sum_avo(ie,je,ke),sum_dvo(ie,je,ke)                        &
            ,wtmix(ie,je,ke),sum_wtmix(ie,je,ke))                            

       SUM_AVO(:,:,:)=0.0
       SUM_DVO(:,:,:)=0.0
       SUM_WTMIX(:,:,:)=0.0
       WTMIX(:,:,:)=0.0

    endif

    if (LGMDIAG) then      
       ALLOCATE( sum_wgo(ie,je,ke),sum_bolx(ie,je),sum_boly(ie,je) )           

       SUM_WGO(:,:,:)=0.0
       SUM_BOLX(:,:)=0.0
       SUM_BOLY(:,:)=0.0

    endif

    if (LCONVDIAG) then
       ALLOCATE(tmceo(ie,je),tmcdo(ie,je),sum_tmceo(ie,je),sum_tmcdo(ie,je))

       TMCEO(:,:)=0.0
       TMCDO(:,:)=0.0

    endif

    if (LHFLDIAG) then
       ALLOCATE( dqswo(ie,je),dqlwo(ie,je),dqseo(ie,je),dqlao(ie,je)       &
            ,dqtho(ie,je),dqswi(ie,je),dqlwi(ie,je),dqsei(ie,je)           &
            ,dqlai(ie,je),dqthi(ie,je),sum_dqswo(ie,je),sum_dqlwo(ie,je)   &
            ,sum_dqseo(ie,je),sum_dqlao(ie,je),sum_dqtho(ie,je)            &
            ,sum_dqswi(ie,je),sum_dqlwi(ie,je),sum_dqsei(ie,je)            &
            ,sum_dqlai(ie,je),sum_dqthi(ie,je),dticeo(ie,je)               &
            ,sum_dticeo(ie,je) )                                    
    endif

#ifdef __coupled
    ALLOCATE( sum_aoflnhw(ie,je),sum_aoflshw(ie,je),sum_aoflchi(ie,je)      &
             ,sum_aoflfrw(ie,je),sum_aoflfri(ie,je),sum_aofltxw(ie,je)      &
             ,sum_aofltyw(ie,je),sum_aofltxi(ie,je),sum_aofltyi(ie,je)      &
             ,sum_aoflwsv(ie,je),sum_aoflrhi(ie,je) )                        

    SUM_AOFLNHW(:,:)=0.0
    SUM_AOFLSHW(:,:)=0.0
    SUM_AOFLRHI(:,:)=0.0
    SUM_AOFLCHI(:,:)=0.0
    SUM_AOFLFRW(:,:)=0.0
    SUM_AOFLFRI(:,:)=0.0
    SUM_AOFLTXW(:,:)=0.0
    SUM_AOFLTYW(:,:)=0.0
    SUM_AOFLTXI(:,:)=0.0
    SUM_AOFLTYI(:,:)=0.0
    SUM_AOFLWSV(:,:)=0.0

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
  REAL field2(ie,je),sum_field2(ie,je)
  REAL(kind=sp) ff_g(ie_g,je_g)
  REAL sum_field2_g(ie_g,je_g)

  ! set extra file format parameters

  icode=i4code
  iedim=ie_g*je_g

  !h at beginning of averaging period : initialize the field

  IF (mnnndt.EQ.1) THEN
     DO i=1,ie
        DO j=1,je
           sum_field2(i,j)=0.      
        ENDDO
     ENDDO
  ENDIF

  !h summation and averaging of the field; check out of range for ieee

  DO i=1,ie
     DO j=1,je
        sum_field2(i,j)=sum_field2(i,j)+(field2(i,j)-sum_field2(i,j))/mnnndt
     ENDDO
  ENDDO

  IF (mnnndt.EQ.mnnn720) THEN

     DO i=1,ie
        DO j=1,je
           IF (ivec.EQ.-1) awet= 1
           IF (ivec.EQ.0) awet= weto(i,j,1)
           IF (ivec.EQ.1) awet= amsuo(i,j,1)
           IF (ivec.EQ.2) awet= amsue(i,j,1)
           IF (awet.LT.0.5) THEN
              sum_field2(i,j)=-9.e33
           ENDIF
!           IF (ABS(sum_field2(i,j)).LT.(1.e-35)) THEN
!              sum_field2(i,j)=0.
!           ENDIF
        ENDDO
     ENDDO

     CALL gather_arr(sum_field2,sum_field2_g,p_io)

     IF (p_pe == p_io) THEN

        !h write in extra format
        WRITE(io_stdout,*)'open unit: ',iun2d
!        OPEN(iun2d, status='unknown',                                  &
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
!        WRITE(iun2d)((ff_g(i,j),i=1,ie_g),j=1,je_g)
!        CLOSE(iun2d)

     ENDIF  ! p_pe==p_io

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


  REAL field(ie,je,ke),sum_field(ie,je,ke)

  REAL sum_field_g(ie_g,je_g,ke)
  INTEGER(kind=i4) idate,klev,icode,iedim
  INTEGER awet,i,j,k,ivec,kdays,kmonts,kyears,mnnndt,mnnn720,iun3d,i4code

  REAL(kind=sp) ff_g(ie_g,je_g)
  iedim=ie_g*je_g
  icode=i4code

  !h    initialize the field

  IF (mnnndt.EQ.1) THEN
     DO k=1,ke
        DO i=1,ie
           DO j=1,je
              sum_field(i,j,k)=0.      
           ENDDO
        ENDDO
     ENDDO
  ENDIF

  !h     summation and monthly mean of the field

  DO k=1,ke
     DO j=1,je
        DO i=1,ie
           sum_field(i,j,k)=sum_field(i,j,k)+(field(i,j,k)-sum_field(i,j,k))/mnnndt
        ENDDO
     ENDDO
  ENDDO

  IF (mnnndt.EQ.mnnn720) THEN


  DO k=1,ke
    DO j=1,je
       DO i=1,ie
           IF (ivec.EQ.-1) awet= 1
           IF (ivec.EQ.0) awet= weto(i,j,k)
           IF (ivec.EQ.1) awet= amsuo(i,j,k)
           IF (ivec.EQ.2) awet= amsue(i,j,k)
           IF (awet.LT.0.5) THEN
              sum_field(i,j,k)=-9.e33
           ENDIF
!           IF (ABS(sum_field(i,j,k)).LT.(1.e-35)) THEN
!              sum_field(i,j,k)=0.
!           ENDIF
        ENDDO
     ENDDO
  ENDDO

     DO k=1,ke
        CALL gather_arr(sum_field(:,:,k),sum_field_g(:,:,k),p_io)
     ENDDO

     idate= (kyears*10000)+(kmonts*100)+kdays 
     WRITE(io_stdout,*)'mmean3d: unit=',iun3d
     WRITE(io_stdout,'(1x,a11,i8)') '  yyyymmdd=',idate
     WRITE(io_stdout,*)'         code=',icode
     WRITE(io_stdout,*)'    nanf/nend=',mnnndt,mnnn720


     IF (p_pe==p_io) THEN
!        OPEN(iun3d,    status='unknown',                            &
!                                 access='sequential',               &
!                                 position='append',                 &
!                                form='unformatted')

        DO k=1,ke


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
!           WRITE(iun3d)((ff_g(i,j),i=1,ie_g),j=1,je_g)
           WRITE(iun3d) ff_g
        ENDDO
!        CLOSE(iun3d)
     ENDIF ! p_pe==p_io
  ENDIF


END SUBROUTINE mmean3d

END MODULE mo_mean







