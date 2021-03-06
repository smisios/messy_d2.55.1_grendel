MODULE mo_diagnosis

  USE mo_param1
  USE mo_commo1
  USE mo_mpi
  USE mo_parallel

  IMPLICIT NONE

  INTEGER, POINTER :: ibek(:,:), ibek_g(:,:)
  REAL, POINTER :: sigh(:,:),zmld(:,:)

  REAL,ALLOCATABLE :: psiuwe_g(:,:)

  INTEGER :: grid_family

  
  INTEGER :: nmecen
  PARAMETER (nmecen=5)

  REAL,ALLOCATABLE ::                                               &
       tberi(:),tfaroe(:),tdenmar(:)                                &
       ,tdavis(:),tmedi(:),hflb(:)                                   &
       ,wflb(:),eiscb(:),eisab(:)                                    &
       ,arem(:,:),tquer(:,:),squer(:,:)                              &
       ,dvquer(:,:),avquer(:,:)                                      &
       ,tmerci(:,:)                                                  &
       ,aabw(:),rnadw(:)                                             &
       ,tvquer(:,:),svquer(:,:)                                      &
       ,tdquer(:,:)


  INTEGER iban1,jban1,iban2,jban2                                   &
       ,idra1,jdra1,idra2,jdra2                                     &
       ,iber1,jber1,iber2,jber2                                     &
       ,idberi,jdberi,iddenm,jddene,jddenw                          &
       ,idmedi,jdmedi                                               &
       ,idfarw,idfare,jdfaro,kdfaro                                 &
       ,iddavi,jddavw,jddave                                        &
       ,idfram,jdfraw,jdfrae                                        &
       ,iddav1,jddav1,iddav2,jddav2                                 &
       ,idguln,idguls,jdgulw,jdgule                                 &
       ,idkurn,idkurs,jdkurw,jdkure                                 &
       ,jdgss,isued,inord                                           &
       ,kaabw,knadw,jlink(nmecen),jrech(nmecen),imerci(nmecen)      &
       ,iddene,iddenw,jddenm                                        &
       ,jddenn,jddens                                               &
       ,jdfarn,jdfars,idfaro                                        &
       ,idgulw,idgule,jdguln,jdguls                                 &
       ,idkurw,idkure,jdkurn,jdkurs                                 &
       ,ilink(nmecen),irech(nmecen),jmerci(nmecen)


    real, allocatable :: tmerc(:,:,:),sum_tmerc(:,:,:)

CONTAINS

  SUBROUTINE alloc_mem_diag

    ! define her the grid families, i.e. take into account grid orientation
    ! for throughflows are different in gin and grob-type grids

    grid_family=1

    ! versiongin
    IF (ie_g.EQ.182.AND.je_g.EQ.84) THEN
       grid_family=2
    ENDIF


    ALLOCATE(tberi(ke),tfaroe(ke),tdenmar(ke),tdavis(ke),tmedi(ke)    &
         ,arem(ke,nbox),tquer(ke,nbox),squer(ke,nbox),dvquer(ke,nbox) &
         ,avquer(ke,nbox),tmerci(ke,nmecen),tvquer(ke,nmecen)         &
         ,svquer(ke,nmecen),tdquer(ke,nmecen),aabw(nmecen),rnadw(nmecen) &
         ,wflb(nbox),hflb(nbox),eiscb(nbox),eisab(nbox))

    ALLOCATE(ibek(ie,je),sigh(ie,je),zmld(ie,je))
    ALLOCATE(ibek_g(ie_g,je_g),psiuwe_g(ie_g,je_g))
    
    ALLOCATE(tmerc(180,2,kep),sum_tmerc(180,2,kep))

  END SUBROUTINE alloc_mem_diag

  SUBROUTINE diag_ini  

    USE mo_param1
    USE mo_mpi
    USE mo_commo1
    USE mo_commoau1
    USE mo_commoau2
    USE mo_units


    INTEGER i,j,k,n
    REAL dist

    ! initialisations
    !:: common diagnostics
    DO j=1,je
       DO i=1,ie
          ibek(i,j)=0
       ENDDO
    ENDDO


#ifdef DIAG
    ! open file for timeseries
    IF(p_pe==p_io) OPEN(io_ou_f125,file='TIMESER')
    !:: initialize diag fields
    DO k=1,ke
       tberi(k)=zero
       tfaroe(k)=zero
       tdenmar(k)=zero
       tdavis(k)=zero
       tmedi(k)=zero
    ENDDO
    DO n=1,nbox
       hflb(n)=zero
       wflb(n)=zero
       eiscb(n)=zero
       eisab(n)=zero
       DO k=1,ke
          arem(k,n)=zero
          tquer(k,n)=zero
          squer(k,n)=zero
          dvquer(k,n)=zero
          avquer(k,n)=zero
       ENDDO
    ENDDO



    DO n=1,nmecen

       IF (grid_family.EQ. 2) THEN
          jlink(n)=0
          jrech(n)=0
          imerci(n)=0
       ELSE
          ilink(n)=0
          irech(n)=0
          jmerci(n)=0
       ENDIF

       aabw(n)=zero
       rnadw(n)=zero

       DO k=1,ke
          tmerci(k,n)=zero
          tvquer(k,n)=zero
          svquer(k,n)=zero
          tdquer(k,n)=zero
       ENDDO

    ENDDO


    IF(ie_g.EQ.182.AND.je_g.EQ.84)THEN           ! versiongin

       imerci(1)=50
       imerci(2)=123
       imerci(3)=134
       imerci(4)=148
       imerci(5)=172
       jlink(1)=64
       jlink(2)=22
       jlink(3)=20
       jlink(4)=35
       jlink(5)=35
       jrech(1)=64
       jrech(2)=55
       jrech(3)=65
       jrech(4)=72
       jrech(5)=72     

    ENDIF

    IF(ie_g.EQ.130.AND.je_g.EQ.211)THEN         ! versiont43

       jmerci(1)=55
       jmerci(2)=48 
       jmerci(3)=55 
       jmerci(4)=69 
       jmerci(5)=145
       ilink(1)=7 
       ilink(2)=42
       ilink(3)=50
       ilink(4)=42
       ilink(5)=50
       irech(1)=13
       irech(2)=88
       irech(3)=76
       irech(4)=72
       irech(5)=85

    ENDIF

    IF(ie_g.EQ.60.AND.je_g.EQ.50)THEN         ! versiongr60
       jmerci(1)=19
       jmerci(2)=14 
       jmerci(3)=17 
       jmerci(4)=21 
       jmerci(5)=34 
       ilink(1)=1 
       ilink(2)=21
       ilink(3)=19
       ilink(4)=17
       ilink(5)=22
       irech(1)=3
       irech(2)=36
       irech(3)=36
       irech(4)=37
       irech(5)=36
    ENDIF


    IF(ie_g.EQ.122.AND.je_g.EQ.101)THEN         ! versiongr30

       jmerci(1)=41
       jmerci(2)=12 
       jmerci(3)=33 
       jmerci(4)=41 
       jmerci(5)=71 
       ilink(1)=120 
       ilink(2)=42
       ilink(3)=44
       ilink(4)=32
       ilink(5)=43
       irech(1)=122
       irech(2)=72
       irech(3)=74
       irech(4)=70
       irech(5)=75
    ENDIF


    IF(ie_g.EQ.256.AND.je_g.EQ.220)THEN         ! versiongr15          
       jmerci(1)=85
       jmerci(2)=34
       jmerci(3)=70
       jmerci(4)=94
       jmerci(5)=155
       ilink(1)=1
       ilink(2)=107
       ilink(3)=78
       ilink(4)=60
       ilink(5)=87
       irech(1)=15
       irech(2)=153
       irech(3)=157
       irech(4)=150
       irech(5)=160
    ENDIF

    IF(ie_g.EQ.400.AND.je_g.EQ.338)THEN         ! versiongr09
       jmerci(1)=135
       jmerci(2)=40
       jmerci(3)=109
       jmerci(4)=123
       jmerci(5)=235
       ilink(1)=398
       ilink(2)=139
       ilink(3)=145
       ilink(4)=105
       ilink(5)=142
       irech(1)=400
       irech(2)=237
       irech(3)=244
       irech(4)=231
       irech(5)=248
    ENDIF


    IF(ie_g.EQ.182.AND.je_g.EQ.84)THEN           ! versiongin

       !        bering
       idberi=50
       jdberi=63

       !        denmark
       iddenm=115
       jddene=40
       jddenw=50

       !        mediterranean outflow
       idmedi=146
       jdmedi=33

       !        faroer-bank
       jdfaro=33
       idfarw=121
       idfare=124
       kdfaro=10

       IF (ke.EQ.40)   kdfaro=15

       IF (ke.EQ.30)   kdfaro=17

       !        davis
       iddavi=92
       jddavw=67
       jddave=69

       !        ice fram
       idfram=92
       jdfraw=33
       jdfrae=45

       !        ice davis
       iddav1=15
       jddav1=64
       iddav2=36
       jddav2=70

       !        gulfstream: max of streamfunction
       idguln=130
       idguls=158
       jdgulw=62
       jdgule=84

       !        kuroshio: max of streamfunction
       idkurn=22
       idkurs=50
       jdkurw=36
       jdkure=53

       !         find gulfstream separation   
       jdgss=61
       isued=158
       inord=135

       !        aabw = max below level kaabw
       kaabw=11
       IF (ke.EQ.40)  kaabw=18
       IF (ke.EQ.30)  kaabw=16

       !        nadw = min above level knadw 
       knadw=9
       IF (ke.EQ.40) knadw=14
       IF (ke.EQ.30) knadw=13
       IF (ke.EQ.23) knadw=12

    ENDIF           ! versiongin



    IF (ie_g .EQ. 60 .AND. je_g .EQ. 50 ) THEN   ! versiongr60
       !        bering
       idberi=2
       jdberi=19
       !
       !        denmark
       iddenm=34
       jddenn=1
       jddens=6
       !
       jddenm=1
       iddenw=34
       iddene=34
       !:: mediterranean outflow
       idmedi=34    
       jdmedi=21
       !
       !        faroer-bank
       idfaro=36
       jdfarn=8
       jdfars=14
       kdfaro=12

       IF ( ke .EQ. 40 ) kdfaro=15
       IF ( ke .EQ. 30 ) kdfaro=17
       IF ( ke .EQ. 23 ) kdfaro=17

       !        davis
       iddavi=16
       jddavw=1
       jddave=9
       !
       !        ice fram
       idfram=52
       jdfraw=3
       jdfrae=11
       !
       !        ice davis
       iddav1=16
       jddav1=1
       iddav2=16
       jddav2=9
       !      gulfstream: max of streamfunction
       idgulw=35
       idgule=45
       jdguln=33
       jdguls=50

       !
       !      kuroshio: max of streamfunction
       idkurw=105
       idkure=120
       jdkurn=50
       jdkurs=60
       !
       !      find gulfstream separation   
       jdgss=45
       isued=45
       inord=35
       !
       !      aabw = max below level kaabw
       kaabw=11

       IF ( ke .EQ. 40 ) kaabw=18
       IF ( ke .EQ. 30 ) kaabw=16
       IF ( ke .EQ. 23 ) kaabw=18



       !      nadw = min above level knadw 
       knadw=9

       IF ( ke .EQ. 40 ) knadw=14
       IF ( ke .EQ. 30 ) knadw=13
       IF ( ke .EQ. 23 ) knadw=12


    ENDIF

    IF (ie_g .EQ. 122 .AND. je_g .EQ. 101 ) THEN   ! versiongr30

       !        bering
       idberi=2
       jdberi=41
       !
       !        denmark
       iddenm=69
       jddenn=2
       jddens=16
       !
       jddenm=2
       iddenw=48
       iddene=72
       !:: mediterranean outflow
       idmedi=69
       jdmedi=43
       !
       !        faroer-bank
       idfaro=73
       jdfarn=15
       jdfars=30
       kdfaro=12

       IF ( ke .EQ. 40 ) kdfaro=15
       IF ( ke .EQ. 30 ) kdfaro=17
       IF ( ke .EQ. 23 ) kdfaro=17


       !        davis
       iddavi=11
       jddavw=20
       jddave=29
       !
       !        ice fram
       idfram=101
       jdfraw=3
       jdfrae=23
       !
       !        ice davis
       iddav1=13
       jddav1=31
       iddav2=15
       jddav2=33
       !      gulfstream: max of streamfunction
       idgulw=35
       idgule=45
       jdguln=33
       jdguls=50

       !
       !      kuroshio: max of streamfunction
       idkurw=105
       idkure=120
       jdkurn=50
       jdkurs=60
       !
       !      find gulfstream separation   
       jdgss=45
       isued=45
       inord=35
       !
       !      aabw = max below level kaabw
       kaabw=11


       IF ( ke .EQ. 40 ) kaabw=18
       IF ( ke .EQ. 30 ) kaabw=16
       IF ( ke .EQ. 23 ) kaabw=18


       !      nadw = min above level knadw 
       knadw=9

       IF ( ke .EQ. 40 ) knadw=14
       IF ( ke .EQ. 30 ) knadw=13
       IF ( ke .EQ. 23 ) knadw=12

    ENDIF


    IF(ie_g.EQ.256.AND.je_g.EQ.220)THEN         ! versiongr15          

       ! bering
       idberi=8
       jdberi=85
       ! denmark
       iddenm=152
       jddenn=2
       jddens=34

       jddenm=2
       iddenw=48
       iddene=72

       ! mediterranean outflow
       idmedi=147
       jdmedi=92

       ! faroer-bank
       idfaro=157
       jdfarn=35
       jdfars=65
       kdfaro=12

       IF (ke.EQ.40) kdfaro=14
       IF (ke.EQ.30) kdfaro=17
       IF (ke.EQ.23) kdfaro=17

       ! davis
       iddavi=22
       jddavw=49
       jddave=80

       ! ice fram
       idfram=214
       jdfraw=3
       jdfrae=50

       ! ice davis
       iddav1=28
       jddav1=42
       iddav2=29
       jddav2=424

       ! gulfstream: max of streamfunction
       idgulw=74
       idgule=96
       jdguln=72
       jdguln=72
       jdguls=98

       ! kuroshio: max of streamfunction
       idkurw=105
       idkure=120
       jdkurn=230
       jdkurs=250

       ! find gulfstream separation
       jdgss=90
       isued=74
       inord=84

       ! aabw = max below level kaabw
       kaabw=11

       IF (ke.EQ.40) kaabw=18
       IF (ke.EQ.30) kaabw=16
       IF (ke.EQ.23) kaabw=18

       ! nadw = min above level knadw
       knadw=9

       IF (ke.EQ.40) knadw=14
       IF (ke.EQ.30) knadw=13
       IF (ke.EQ.23) knadw=12

    ENDIF

    IF(ie_g.EQ.400.AND.je_g.EQ.338)THEN         ! versiongr09          

       !        bering
       idberi=6
       jdberi=135
       !
       !        denmark
       iddenm=228
       jddenn=7
       jddens=52
       !
       jddenm=6
       iddenw=159
       iddene=238
       !:: mediterranean outflow
       idmedi=228
       jdmedi=142
       !
       !        faroer-bank
       idfaro=241
       jdfarn=50
       jdfars=99
       kdfaro=40

       IF (ke.EQ.40) kdfaro=15
       IF (ke.EQ.30) kdfaro=17
       IF (ke.EQ.23) kdfaro=17

       !        davis
       iddavi=37
       jddavw=66
       jddave=96
       !        ice fram
       idfram=333
       jdfraw=10
       jdfrae=76
       !
       !        ice davis
       iddav1=43
       jddav1=102
       iddav2=50
       jddav2=109
       !      gulfstream: max of streamfunction
       idgulw=35
       idgule=45
       jdguln=33
       jdguls=50
       !
       !      kuroshio: max of streamfunction
       idkurw=105
       idkure=120
       jdkurn=50
       jdkurs=60
       !
       !      find gulfstream separation
       jdgss=45
       isued=45
       inord=35
       !
       !      aabw = max below level kaabw
       kaabw=11

       IF (ke.EQ.40) kaabw=18
       IF (ke.EQ.30) kaabw=16
       IF (ke.EQ.23) kaabw=18

       !      nadw = min above level knadw
       knadw=9
       IF (ke.EQ.40) knadw=14
       IF (ke.EQ.30) knadw=13
       IF (ke.EQ.23) knadw=12

    ENDIF
    !

    IF(ie_g.EQ.130.AND.je_g.EQ.211)THEN         ! versiont43

       !         bering
       idberi=10
       jdberi=55

       !         denmark strait
       jddenm=40
       iddene=67
       iddenw=68

       iddenm=67
       jddenn=34
       jddens=40

       !         mediterranean
       idmedi=75
       jdmedi=65

       !         faroe bank
       idfaro=76
       jdfarn=45
       jdfars=50

       !         fram strait     
       idfram=101
       jdfraw=7
       jdfrae=25

       !         davis strait
       iddavi=28 
       jddavw=39
       jddave=39

       !         gulfstream: max of streamfunction
       idgulw=40 
       idgule=65 
       jdguln=50
       jdguls=80

       !         kuroshio: max of streamfunction
       idkurw=115
       idkure=125
       jdkurn=60
       jdkurs=80

       !         find gulfstream separation   
       jdgss=61
       isued=158
       inord=135

       !         aabw = max below level kaabw
       kaabw=11
       !         nadw = min above level knadw 
       knadw=9

       IF (ke.EQ.40) THEN
          kaabw=33
          knadw=23
          kdfaro=27
       ENDIF

       IF (ke.EQ.23) THEN
          kaabw=81
          knadw=12
       ENDIF

    ENDIF


#endif /*diag*/


    !:: search for streamfunction points
    CALL suchij(-20.,-60.,1,idra2,jdra2,dist,0.)
    CALL suchij(-89.,-99.,1,idra1,jdra1,dist,0.)
    CALL suchij(30.,100.,1,iban1,jban1,dist,0.)
    CALL suchij(-25.,135.,1,iban2,jban2,dist,0.)
    CALL suchij(67.,160.,1,iber1,jber1,dist,0.)
    CALL suchij(65.,-150.,1,iber2,jber2,dist,0.)
    !
    WRITE(io_stdout,*)'drake: ',idra1,jdra1,weto_g(idra1,jdra1,1)    &
         &                    ,idra2,jdra2,weto_g(idra2,jdra2,1)
    WRITE(io_stdout,*)'banda: ',iban1,jban1,weto_g(iban1,jban1,1)    &
         &                    ,iban2,jban2,weto_g(iban2,jban2,1)
    WRITE(io_stdout,*)'bering: ',iber1,jber1,weto_g(iber1,jber1,1)   &
         &                    ,iber2,jber2,weto_g(iber2,jber2,1)
    !::

  END SUBROUTINE diag_ini

  SUBROUTINE calc_mixedlayerdepth

    USE mo_param1
    USE mo_mpi
    USE mo_parallel
    USE mo_commo1
    USE mo_commoau1
    USE mo_commoau2
    USE mo_units

    REAL zzzuwe,zzzcou,zzz,sigcrit

    INTEGER i,j,k,m

    !   compute mixed-layer depth

!$OMP PARALLEL PRIVATE(i,j,k,m)

    sigcrit=0.125
    m=lmonts


    sigh(:,:)=zero
    zmld(:,:)=zero

    !$OMP DO
    DO j=1,je       
       DO i=1,ie
          sigh(i,j)=weto(i,j,1)*sigcrit
          zmld(i,j)=weto(i,j,1)*tiestu(1)
       ENDDO
    ENDDO
    !$OMP END DO

    DO k=2,ke
       !$OMP DO
       DO j=1,je
          DO i=1,ie
             IF(weto(i,j,k)*sigh(i,j).GT.1.e-6)THEN
                zzz=MIN(sigh(i,j)/(ABS(stabio(i,j,k))+almzer),dz(k))
                sigh(i,j)=MAX(0.,sigh(i,j)-zzz*stabio(i,j,k))
                zmld(i,j)=zmld(i,j)+zzz
             ELSE
                sigh(i,j)=0.
             ENDIF
          ENDDO
       ENDDO
       !$OMP END DO
    ENDDO
    !
    DO j=2,je-1
       DO i=2,ie-1
          amld(i,j,m)=MAX(amld(i,j,m),zmld(i,j))
       ENDDO
    ENDDO

    zzzuwe = SUM(zmld(2:ie-1,2:je-1)*weto(2:ie-1,2:je-1,1))
    zzzcou = SUM(weto(2:ie-1,2:je-1,1))
!$OMP SINGLE
    CALL global_sum(zzzuwe,zzzcou)
!$OMP END SINGLE
    zzzcou = zzzcou + almzer

    WRITE(io_stdout,*)'mean mld: ',zzzuwe/zzzcou
    !hh      end mixed layer depths

!$OMP END PARALLEL

  END SUBROUTINE calc_mixedlayerdepth



  SUBROUTINE calc_icecutoff
  

    USE mo_param1
    USE mo_mpi
    USE mo_parallel
    USE mo_commo1
    USE mo_commoau1
    USE mo_commoau2
    USE mo_units
#ifdef PBGC
    USE mo_carbch, only: ocetra
    USE mo_param1_bgc, only: nocetra
    INTEGER l
#endif

    REAL swmin,eismax,volice,areice,volsno,arges,zmi,zma  & 
,schwell,ddice,dlxpdlyp

    INTEGER i,j,icou,ierr
    
    swmin=0.
    icou=0
    eismax=0.
    volice=0.
    areice=0.
    volsno=0.
    arges=0.
    zmi=0.
    zma=0.

    ierr=0


    !uwe     ensure that effective thickness of uppermost waterlayer
    !        does not fall below a critical value
    !        in that case ice is converted to water
    !        not heat conserving!!!!
    !         

    DO j=1,je
       DO i=1,ie
          IF(weto(i,j,1).GT.0.5)THEN
             schwell=0.7*(dzw(1)+zo(i,j))
             IF(rhoicwa*sictho(i,j)+rhosnwa*sicsno(i,j).GT.schwell)THEN
                !-et
                WRITE(io_stdout,*)'ice cut-off at ',i,j,sictho(i,j),sao(i,j,1)
                ierr=1
                ddice=MIN((sictho(i,j)*rhoicwa+sicsno(i,j)*rhosnwa     &
                     -schwell)/rhoicwa,sictho(i,j))
                !
                sictho(i,j)=sictho(i,j)-ddice
                !uwe       salt conservation!
                sao(i,j,1)=(ddice*sice+sao(i,j,1)                      &
                     *(schwell-ddice*rhoicwa))/schwell
#ifdef PBGC
                do l=1,nocetra
                   ocetra(i,j,1,l)=(ocetra(i,j,1,l)                       &
                                *(schwell-ddice*rhoicwa))/schwell
                enddo
#endif
                WRITE(io_stdout,*)'eisneu,sneu: ',i,j,sictho(i,j),sao(i,j,1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    !
    DO j=2,je-1
       DO i=2,ie-1
          IF(weto(i,j,1).GT.0.5.AND.giph(2*i,2*j).GT.0.)THEN
             icou=icou+1
             !        wetie=wetie+dlxp(1,j)*dlyp(1,j)
             arges=arges+dlxp(i,j)*dlyp(i,j)
             swmin=MIN(swmin,fswr(i,j))
             eismax=MAX(eismax,sictho(i,j))
             zmi=MIN(zmi,zo(i,j))
             zma=MAX(zma,zo(i,j))
             volice=volice+dlxp(i,j)*dlyp(i,j)*sictho(i,j)
             volsno=volsno+dlxp(i,j)*dlyp(i,j)*sicsno(i,j)
             !-et
             dlxpdlyp=dlxp(i,j)*dlyp(i,j)
             IF(sicomo(i,j).LT.0.2)THEN
                dlxpdlyp=dlxpdlyp*MAX(0.,(sicomo(i,j)-0.1)/0.1)
             ENDIF
             areice=areice+dlxpdlyp
          ENDIF
       ENDDO
    ENDDO

    CALL global_sum(icou)
    CALL global_sum(arges,volice,volsno,areice)
    CALL global_max(eismax,zma)
    CALL global_min(swmin,zmi)

    WRITE(io_stdout,*)'eis: ',areice*1.e-12,volice*1.e-12          &
         ,volsno*1.e-12,eismax,' zeta: ',zmi,zma


  END SUBROUTINE calc_icecutoff

  SUBROUTINE calc_psi

    use mo_commo1
    use mo_parallel
    use mo_param1
    
    INTEGER I,J,K

    REAL zhilf_g(ie_g,je_g)

    zhilf_g(:,:)=0.0
    psiuwe_g(:,:)=0.0

    psiuwe(:,:)=0.

    DO j=2,je
       DO k=1,ke
          DO i=1,ie
             psiuwe(i,j)=psiuwe(i,j)+uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)
          ENDDO
       ENDDO
    ENDDO

  !  now gather and broadcast global psiuwe array, the following loop
  !  sums up globally

    !        gather global psiuwe array

    CALL gather_arr(psiuwe,psiuwe_g,p_io)
    CALL p_bcast(psiuwe_g,p_io)

    IF(p_pe==p_io) THEN
       DO j=2,je_g
          DO i=1,ie_g
             psiuwe_g(i,j)=psiuwe_g(i,j-1)+psiuwe_g(i,j)
          ENDDO
       ENDDO
       DO j=2,je_g-1
          DO i=2,ie_g-1
             zhilf_g(i,j)=0.25*(psiuwe_g(i,j)+psiuwe_g(i,j-1) &
                  +psiuwe_g(i-1,j)+psiuwe_g(i-1,j-1))
          ENDDO
       ENDDO
       DO j=1,je_g
          DO i=1,ie_g
             psiuwe_g(i,j)=zhilf_g(i,j)
          ENDDO
       ENDDO
    ENDIF


    ! for later local use:
!!    psiuwe(1:ie,1:je) = psiuwe_g(p_ioff+1:p_ioff+ie,p_joff+1:p_joff+je)

    CALL scatter_arr(psiuwe_g,psiuwe,p_io)


  END SUBROUTINE calc_psi


  SUBROUTINE diagnosis

    USE mo_param1
    USE mo_mpi
    USE mo_parallel
    USE mo_commo1
    USE mo_commoau1
    USE mo_commoau2
    USE mo_units

    INTEGER i,j,k,m,n,LDENMAR,LFAROE,LTQ1,LTQ2,LTQ3,LTQ4

    REAL kmaxlon, kmaxlat,grarad,psigolf,psiaspg            &
               ,gmaxlat,gmaxlon,smaxlat,smaxlon,psikuro,faks &
              ,rberi,RINM,RUTM,SMEDI,RIND,RUTD,RINF,RUTF,RINA &
              ,RUTA,SFRAM,SDAVIS,SDAV1,SDAV2,CPRHO

    CALL calc_mixedlayerdepth

    CALL calc_icecutoff

!    CALL calc_psi


    !******************************************************************

    IF(p_pe == p_io) THEN
       pi=4.*ATAN(1.)
       grarad=180./pi

       ! gulfstreamindex
       psigolf=0.
       !!cdir novector
       DO j=1,je_g 
          DO i=1,ie_g
             IF (ibek_g(i,j).EQ.4.OR.ibek_g(i,j).EQ.5) THEN
                psigolf=MAX(psigolf,psiuwe_g(i,j))
                IF (psigolf .EQ. psiuwe_g(i,j)) THEN
                   gmaxlat=grarad*giph_g((2*i)+1,(2*j)+1)
                   gmaxlon=grarad*gila_g((2*i)+1,(2*j)+1)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       ! atl spgyre index
       psiaspg=0.
       !!cdir novector
       DO j=1,je_g
          DO i=1,ie_g
             IF (ibek_g(i,j).EQ.4) THEN
                psiaspg=MIN(psiaspg,psiuwe_g(i,j))
                IF (psiaspg .EQ. psiuwe_g(i,j)) THEN
                   smaxlat=grarad*giph_g((2*i)+1,(2*j)+1)
                   smaxlon=grarad*gila_g((2*i)+1,(2*j)+1)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       ! kuroshioindex
       psikuro=0.
       !!cdir novector
       DO j=1,je_g-1
          DO i=1,ie_g-1
             IF (ibek_g(i,j).EQ.7) THEN
                psikuro=MAX(psikuro,psiuwe_g(i,j))
                IF (psikuro .EQ. psiuwe_g(i,j)) THEN
                   kmaxlat=grarad*giph_g((2*i)+1,(2*j)+1)
                   kmaxlon=grarad*gila_g((2*i)+1,(2*j)+1)
                ENDIF
             ENDIF
          ENDDO
       ENDDO

       WRITE(io_stdout,*)'golfstream: ', NINT(psigolf*1.e-6),                  &
            'at lat/lon ',gmaxlat,gmaxlon
       WRITE(io_stdout,*)'a-subpolar gyer: ', NINT(psiaspg*1.e-6),             &
            'at lat/lon ',smaxlat,smaxlon
       WRITE(io_stdout,*)'kuroshio: ', NINT(psikuro*1.e-6),                    &
            'at lat/lon ',kmaxlat,kmaxlon
       WRITE(io_stdout,*)'banda: '                                             &
            ,NINT((psiuwe_g(iban1,jban1)-psiuwe_g(iban2,jban2))*1.e-6)         &
            ,'  drake: '                                                       &
            ,NINT((psiuwe_g(idra1,jdra1)-psiuwe_g(idra2,jdra2))*1.e-6)         &
            ,'  bering*10: '                                                   &
            ,NINT((psiuwe_g(iber1,jber1)-psiuwe_g(iber2,jber2))*1.e-5)

    ENDIF
    !****************************************************************************


    faks=1./float(30*(lmont2+1-lmont1))

    DO j=1,je
       DO i=1,ie
          hflm(i,j)=hflm(i,j)+(qswo(i,j)+qlwo(i,j)+qlao(i,j)+qseo(i,j))*faks
          pmem(i,j)=pmem(i,j)+(prech(i,j)+eminpo(i,j))*faks
       ENDDO
    ENDDO

    DO j=1,je1 
       DO i=1,ie1
          eistrx(i,j)=eistrx(i,j)                                    &
               +sicuo(i,j)*0.5*(sictho(i,j)+sictho(i+1,j)               &
               +(sicsno(i,j)+sicsno(i+1,j))*rhosnic)*faks
          eistry(i,j)=eistry(i,j)                                    &
               +sicve(i,j)*0.5*(sictho(i,j)+sictho(i,j+1)               &
               +(sicsno(i,j)+sicsno(i,j+1))*rhosnic)*faks
       ENDDO
    ENDDO

    CALL bounds_exch('u',eistrx,'mo_diagnosis 1')
    CALL bounds_exch('v',eistry,'mo_diagnosis 2')

!====================================================================

    !
#ifdef DIAG
    !        diagnostic i
    !          transporte
    WRITE(io_stdout,*)'vor 714'
    !
    !cdir novector
    DO k=ke,1,-1
       tmedi(k)=0.0
       tberi(k)=0.
       tdavis(k)=0.
       tfaroe(k)=0.
       tdenmar(k)=0.
       IF(k.LT.ke)tmedi(k)=tmedi(k+1)
       IF(k.LT.ke)tberi(k)=tberi(k+1)
       IF(k.LT.ke)tfaroe(k)=tfaroe(k+1)
       IF(k.LT.ke)tdenmar(k)=tdenmar(k+1)
       IF(k.LT.ke)tdavis(k)=tdavis(k+1)
       !
       !hh          bering 
       IF (grid_family.EQ.2) THEN 
          i=idberi-p_ioff
          j=jdberi-p_joff
          rberi=0.
          IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
             rberi = uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)
          ENDIF

          !hh          mediterranean
          i=idmedi-p_ioff
          j=jdmedi-p_joff
          rinm=0.
          rutm=0.
          smedi=0.
          IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
             rinm=rinm+MAX(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
             rutm=rutm+MIN(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
             smedi=sao(i,j,k)*weto(i,j,k)
          ENDIF

          !hh          denmark
          i=iddenm-p_ioff
          rind=0.
          rutd=0.
          IF(i>1 .AND. i<ie) THEN
             DO j=MAX(2,jddene-p_joff),MIN(je-1,jddenw-p_joff)
                rind=rind+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                rutd=rutd+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF

          !hh          faroe k+=> 700m
          !hh          if (k.ge.kdfaro) then
          j=jdfaro-p_joff
          rinf=0.
          rutf=0.
          IF(j>1 .AND. j<je) THEN
             DO i=MAX(2,idfarw-p_ioff),MIN(ie-1,idfare-p_ioff)
                rinf=rinf+MIN(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
                rutf=rutf+MAX(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
             ENDDO
          ENDIF

          !hh          davis
          rina=0.
          ruta=0.   
          i=iddavi-p_ioff
          IF(i>1 .AND. i<ie) THEN
             DO  j=MAX(2,jddavw-p_joff),MIN(je-1,jddave-p_joff)
                rina=rina+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                ruta=ruta+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF
       ELSE
          i=idberi-p_ioff
          j=jdberi-p_joff
          rberi=0.
          IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
             rberi=vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)
          ENDIF
          !
          !hh          mediterranean
          i=idmedi-p_ioff
          j=jdmedi-p_joff
          rinm=0.
          rutm=0.
          smedi=0.
          IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
             rinm=rinm+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             rutm=rutm+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             smedi=sao(i,j,k)*weto(i,j,k)
          ENDIF
          ! 
          !hh          denmark
          j=jddenm-p_joff
          rind=0.
          rutd=0.
          IF(j>1 .AND. j<je) THEN
             DO i=MAX(2,iddene-p_ioff),MIN(ie-1,iddenw-p_ioff)
                rind=rind+MAX(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
                rutd=rutd+MIN(0.,vke(i,j,k))*dlxv(i,j)*ddue(i,j,k)
             ENDDO
          ENDIF
          i=iddenm-p_ioff
          IF(i>1 .AND. i<ie) THEN
             DO j=MAX(2,jddenn-p_joff),MIN(je-1,jddens-p_joff)
                rind=rind+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                rutd=rutd+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF
          !
          !hh          faroe k+=> 700m
          !hh          if (k.ge.kdfaro) then
          i=idfaro-p_ioff
          rinf=0.
          rutf=0.
          IF(i>1 .AND. i<ie) THEN
             DO j=MAX(2,jdfarn-p_joff),MIN(je-1,jdfars-p_joff)
                rinf=rinf+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                rutf=rutf+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF


          !hh          davis
          rina=0.
          ruta=0.   
          i=iddavi-p_ioff
          IF(i>1 .AND. i<ie) THEN
             DO  j=MAX(2,jddavw-p_joff),MIN(je-1,jddave-p_joff)
                rina=rina+MAX(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
                ruta=ruta+MIN(0.,uko(i,j,k))*dlyu(i,j)*dduo(i,j,k)
             ENDDO
          ENDIF
       ENDIF

       CALL global_sum(rberi,smedi,rinm,rutm,rind,rutd,rinf,rutf,rina,ruta)

       tberi(k)=rberi+tberi(k)
       tdavis(k)=rina+ruta+tdavis(k)
       tmedi(k)=rinm+rutm+tmedi(k)
       tfaroe(k)=rinf+rutf+tfaroe(k)
       tdenmar(k)=rind+rutd+tdenmar(k)

       !jj test
       IF (MOD(ldays,30).EQ.5) THEN
          WRITE(io_stdout,'(f5.0,a,2f6.3,a,2f6.3,a,f6.3,a,f6.3,a,2f6.3,f6.2)') &
               tiestu(k),                           &
               ' denm:  ', (rind+rutd)*1.e-6, tdenmar(k)*1.e-6, &
               ' faroe: ', (rinf+rutf)*1.e-6, tfaroe(k) *1.e-6, &
               ' ber:   ',                    tberi(k)  *1.e-6, &
               ' dav:   ',                    tdavis(k) *1.e-6, &
               ' med:   ', (rinm+rutm)*1.e-6, tmedi(k)  *1.e-6, smedi
       ENDIF

    ENDDO

    !          ice transport fram strait
    sfram=0.

    i=idfram-p_ioff
    IF(i>1 .AND. i<ie) THEN
       !cdir novector
       DO j=MAX(2,jdfraw-p_joff),MIN(je-1,jdfrae-p_joff)
          !hh             sfram=sfram+0.5*(sicuo(i,j)+sicuo(i-1,j))
          !hh     x            *(sictho(i,j)+sicsno(i,j)*rhosno/rhoice)*dlyp(i,j)
          sfram=sfram+((ABS(sicuo(i,j))+sicuo(i,j))                    &
               *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
               +((sicuo(i,j)-ABS(sicuo(i,j))))               &
               *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
               *0.5*dlyu(i,j)
          !
       ENDDO
    ENDIF
    CALL global_sum(sfram)
    sdavis=0.
    i=iddavi-p_ioff
    IF(i>1 .AND. i<ie) THEN
       !cdir novector
       DO j=MAX(2,jddavw-p_joff),MIN(je-1,jddave-p_joff)
          !hh          sdavis=sdavis+0.5*(sicuo(i,j)+sicuo(i-1,j))
          !hh  x             *(sictho(i,j)+sicsno(i,j)*rhosnic)*dlyp(i,j)
          sdavis=sdavis+((ABS(sicuo(i,j))+sicuo(i,j))               &
               *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
               +((sicuo(i,j)-ABS(sicuo(i,j))))               &
               *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
               *0.5*dlyu(i,j)
          !
       ENDDO
    ENDIF
    CALL global_sum(sdavis)

    IF (grid_family.EQ.2) THEN
       i=iddav1-p_ioff
       j=jddav1-p_joff
       sdav1=0
       IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
          sdav1=((ABS(sicuo(i,j))+sicuo(i,j))                       &
               *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
               +((sicuo(i,j)-ABS(sicuo(i,j))))               &
               *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
               *0.5*dlyu(i,j)
       ENDIF
       !
       i=iddav2-p_ioff
       j=jddav2-p_joff
       sdav2=0
       IF(i>1 .AND. i<ie .AND. j>1 .AND. j<je) THEN
          sdav2=((ABS(sicuo(i,j))+sicuo(i,j))                       &
               *(sictho(i,j)+sicsno(i,j)*rhosnic)                 &
               +((sicuo(i,j)-ABS(sicuo(i,j))))               &
               *(sictho(i+1,j)+sicsno(i+1,j)*rhosnic))            &
               *0.5*dlyu(i,j)
       ENDIF
       CALL global_sum(sdav1,sdav2)

       IF(MOD(ldays,10).EQ.5)                                       &
            WRITE(io_stdout,*)'eis fram ', NINT(sfram), ' davis ',NINT(sdavis)
    ENDIF

    DO n=1,nbox
       hflb(n)=0.
       wflb(n)=0.
       eiscb(n)=0.
       eisab(n)=0.
       DO k=1,ke
          arem(k,n)=0.
          tquer(k,n)=0.
          squer(k,n)=0.
          !
          avquer(k,n)=0.
          dvquer(k,n)=0.
       ENDDO
    ENDDO
    !hh test
    !cdir novector
    DO k=1,ke
       DO j=2,je-1
          DO i=2,ie-1
             n=ibek(i,j)
             IF(weto(i,j,k).GT.0.5.AND.n.GT.0)THEN
                arem(k,n)=arem(k,n)+dlxp(i,j)*dlyp(i,j)
                tquer(k,n)=tquer(k,n)+dlxp(i,j)*dlyp(i,j)*tho(i,j,k)
                squer(k,n)=squer(k,n)+dlxp(i,j)*dlyp(i,j)*sao(i,j,k)
                avquer(k,n)=avquer(k,n)+dlxp(i,j)*dlyp(i,j)*avo(i,j,k)
                dvquer(k,n)=dvquer(k,n)+dlxp(i,j)*dlyp(i,j)*dvo(i,j,k)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    WRITE(io_stdout,*)'check -5'
    !cdir novector
    DO j=2,je-1
       DO i=2,ie-1
          n=ibek(i,j)
          IF(weto(i,j,1).GT.0.5.AND.n.GT.0)THEN
             hflb(n)=hflb(n)+dlxp(i,j)*dlyp(i,j)*                   &
                  (qswo(i,j)+qlwo(i,j)+qlao(i,j)+qseo(i,j))
             wflb(n)=wflb(n)+dlxp(i,j)*dlyp(i,j)                    &
                  *(prech(i,j)+eminpo(i,j))
             eiscb(n)=eiscb(n)+dlxp(i,j)*dlyp(i,j)                  &
                  *(sictho(i,j)+(sicsno(i,j)*rhosnic))
             IF(sicomo(i,j).GE.0.15)                                &
                  eisab(n)=eisab(n)+dlxp(i,j)*dlyp(i,j)
          ENDIF
       ENDDO
    ENDDO
    CALL global_sum(hflb)
    CALL global_sum(wflb)
    CALL global_sum(eiscb)
    CALL global_sum(eisab)
    CALL global_sum(arem)
    arem(:,:) = arem(:,:) + 1.e-10
    CALL global_sum(tquer)
    CALL global_sum(squer)
    CALL global_sum(avquer)
    CALL global_sum(dvquer)




    DO k=1,ke
       DO n=1,nmecen
          tvquer(k,n)=0.
          svquer(k,n)=0.
       ENDDO
    ENDDO

    DO k=ke,1,-1
       DO n=1,nmecen
          IF(k.EQ.ke)THEN
             tvquer(k,n)=0.
             svquer(k,n)=0.
          ELSE
             tvquer(k,n)=tvquer(k+1,n)
             svquer(k,n)=svquer(k+1,n)
          ENDIF

          cprho=rocp
          IF (grid_family.EQ.2) THEN

             i=imerci(n)-p_ioff
             IF(i>1 .AND. i<ie) THEN
                !cdir novector
                DO j=MAX(2,jlink(n)-p_joff),MIN(je-1,jrech(n)-p_joff)
                   tvquer(k,n)=tvquer(k,n)+uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)&
                        *(tho(i,j,k)+tho(i+1,j,k))*0.5*cprho
                   svquer(k,n)=svquer(k,n)+uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)&
                        *(sao(i,j,k)+sao(i+1,j,k))*0.5
                ENDDO
             ENDIF

          ELSE 
             j=jmerci(n)-p_joff
             IF(j>1 .AND. j<je) THEN
                !cdir novector
                DO i=MAX(2,ilink(n)-p_ioff),MIN(ie-1,irech(n)-p_ioff)
                   tvquer(k,n)=tvquer(k,n)+vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)&
                        *(tho(i,j,k)+tho(i,j+1,k))*0.5*cprho
                   svquer(k,n)=svquer(k,n)+vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)&
                        *(sao(i,j,k)+sao(i,j+1,k))*0.5
                ENDDO
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    CALL global_sum(tvquer)
    CALL global_sum(svquer)

    WRITE(io_stdout,*)'atlantic overturning: '
    DO k=ke,1,-1
       DO n=1,nmecen
          IF(k.EQ.ke)THEN
             tmerci(k,n)=0.
          ELSE
             tmerci(k,n)=tmerci(k+1,n)
          ENDIF
          IF (grid_family.EQ.2) THEN
             i=imerci(n)
             IF(weto_g(i,jlink(n),k)+weto_g(i,jrech(n),k).GT.0.5)THEN
                WRITE(io_stdout,*)'alarm! jlink,jrech ',n,i,jlink(n)   &
                     ,weto_g(i,jlink(n),k),jrech(n),weto_g(i,jrech(n),k)
             ENDIF
             i=imerci(n)-p_ioff
             IF(i>1 .AND. i<ie) THEN
                !cdir novector
                DO j=MAX(2,jlink(n)-p_joff),MIN(je-1,jrech(n)-p_joff)
                   tmerci(k,n)=tmerci(k,n)+uko(i,j,k)*dlyu(i,j)*dduo(i,j,k)
                ENDDO
             ENDIF
          ELSE
             j=jmerci(n)
             IF(weto_g(ilink(n),j,k)+weto_g(irech(n),j,k).GT.0.5)THEN
                WRITE(io_stdout,*)'alarm! ilink,irech ',n,j,ilink(n)   &
                     ,irech(n),weto_g(ilink(n),j,k),weto_g(irech(n),j,k)
             ENDIF
             j=jmerci(n)-p_joff
             IF(j>1 .AND. j<je) THEN
                !cdir novector
                DO i=MAX(2,ilink(n)-p_ioff),MIN(ie-1,irech(n)-p_ioff)
                   tmerci(k,n)=tmerci(k,n)+vke(i,j,k)*dlxv(i,j)*ddue(i,j,k)
                ENDDO
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    CALL global_sum(tmerci)

    DO k=1,ke
       WRITE(io_stdout,'(f8.0,10f10.3)')tiestw(k),(tmerci(k,n)*1.e-6,n=1,nmecen)
    ENDDO

    IF (grid_family.EQ.2) THEN
       WRITE(io_stdout,'(f8.0,10f10.3)')tiestw(1)                   &      
            ,(float(imerci(n)),n=1,nmecen)
    ELSE  
       WRITE(io_stdout,'(f8.0,10f10.3)')tiestw(1)                   &
            ,(float(jmerci(n)),n=1,nmecen)
    ENDIF

    !hh        find aabw(max) and nadw(min) at 
    !hh        aabw below level kaabw
    !hh        nadw above level knadw 
    DO n=1,nmecen
       aabw(n)=0.
       rnadw(n)=0.
    ENDDO

    DO n=2,nmecen
       IF (grid_family.EQ.2) THEN
          !cdir novector
          DO k=ke,kaabw,-1
             aabw(n)=MIN(tmerci(k,n),aabw(n))
          ENDDO
       ELSE  
          !cdir novector
          DO k=ke,kaabw,-1
             aabw(n)=MAX(tmerci(k,n),aabw(n))
          ENDDO
       ENDIF

       IF (grid_family.EQ.2) THEN
          !cdir novector
          DO k=ke,knadw,-1
             rnadw(n)=MAX(tmerci(k,n),rnadw(n))
          ENDDO
       ELSE
          !cdir novector 
          DO k=ke,knadw,-1
             rnadw(n)=MIN(tmerci(k,n),rnadw(n))
          ENDDO
       ENDIF
    ENDDO
    !hh           write(io_stdout,*)'check 4'
    !cdir novector
    DO n=1,nbox-1
       eisab(nbox)=eisab(nbox)+eisab(n)
       eiscb(nbox)=eiscb(nbox)+eiscb(n)
       hflb(nbox)=hflb(nbox)+hflb(n)
       wflb(nbox)=wflb(nbox)+wflb(n)
       !cdir novector
       DO k=1,ke
          arem(k,nbox)=arem(k,nbox)+arem(k,n)
          squer(k,nbox)=squer(k,nbox)+squer(k,n)
          tquer(k,nbox)=tquer(k,nbox)+tquer(k,n)
       ENDDO
    ENDDO
    !hh           write(io_stdout,*)'check 5'
    !
    !hh        default 20 levels
    ldenmar=9
    lfaroe=11
    ltq1=1
    ltq2=7
    ltq3=11
    ltq4=16

    IF (ke.EQ.40) THEN
       ldenmar=16
       lfaroe=19
       ltq1=1
       ltq2=13
       ltq3=21
       ltq4=31
    ENDIF


    IF (ke.EQ.23) THEN
       ldenmar=11
       lfaroe=12
       ltq1=1
       ltq2=9
       ltq3=14
       ltq4=19
    ENDIF

    IF (ke.EQ.30) THEN

       ldenmar=12
       lfaroe=16
       ltq1=1
       ltq2=9
       ltq3=16
       ltq4=22
    ENDIF


    WRITE(io_ou_f125,'(10e13.6)')float(lyears)                   &
         +(float(lmonts-1)+float(ldays)/float(monlen(lmonts)))/12.    &
         ,psigolf, psikuro                                            &
         ,((psiuwe_g(iban1,jban1)-psiuwe_g(iban2,jban2))*1.e0)        &
         ,((psiuwe_g(idra1,jdra1)-psiuwe_g(idra2,jdra2))*1.e0)        &
         ,((psiuwe_g(iber1,jber1)-psiuwe_g(iber2,jber2))*1.e0)        &
         ,psiaspg                                                     &
         ,(aabw(m),rnadw(m),m=1,nmecen)                               &
         ,(tvquer(1,m),svquer(1,m),tmerci(1,m),m=1,nmecen)            &
         ,(tvquer(1,m)-tvquer(1,1),svquer(1,m)-svquer(1,1),m=2,nmecen)&
         ,tberi(1),tdenmar(ldenmar),tfaroe(lfaroe),sfram              &
         ,(eisab(n),eiscb(n),hflb(n),wflb(n)                          &
         ,tquer(ltq1,n)/arem(ltq1,n),squer(ltq1,n)/arem(ltq1,n)       &
         ,tquer(ltq2,n)/arem(ltq2,n),squer(ltq2,n)/arem(ltq2,n)       &
         ,tquer(ltq3,n)/arem(ltq3,n),squer(ltq3,n)/arem(ltq3,n)       &
         ,tquer(ltq4,n)/arem(ltq4,n),squer(ltq4,n)/arem(ltq4,n)       &
         ,n=1,nbox)
    

#endif /*diag*/


    !      **************************************************************
    !      *******************end of common diagnostic*******************
    !      **************************************************************
  END SUBROUTINE diagnosis


  subroutine calc_moc

    use mo_param1, only :ie_g,je_g,ke,kep
 
    use mo_parallel, only :gather_arr
    
    use mo_commo1, only : wo,depto,weto,      &
         tiestw,dzw,dlxp,dlyp,alat,weto_g,    &
         lyears,ldays,lmonts,ldtdayc,dt,imocdiag

    use mo_units, only :  io_ou_mocg,io_ou_moca

    implicit none

    real :: wo_g(ie_g,Je_g,kep)
    real :: depto_g(ie_g,je_g)
    real :: dlxp_g(ie_g,je_g),dlyp_g(ie_g,je_g),alat_g(ie_g,je_g)

    real :: zlat
    integer :: i,j,k,jbrei,lbrei,l,lb,i1,i2,i3,i4                   &
              ,nanf,nend,NDTDAY

    tmerc(:,:,:)=0.

    NDTDAY=NINT(86400./DT)

!    write(0,*)'in calc moc', imocdiag

!HH Daily average
      IF (IMOCDIAG.EQ.1) THEN
         NANF=LDTDAYC
         NEND=NDTDAY
!HH Monthly average
      ELSEIF (IMOCDIAG.EQ.2) THEN
         IF (LDAYS.EQ.1) THEN
            NANF=LDTDAYC
         ELSEIF (LDAYS.GT.1) THEN
            NANF=LDTDAYC+((LDAYS-1)*NDTDAY)
         ENDIF
         NEND=NDTDAY*MONLEN(LMONTS)
!HH Yearly average
      ELSEIF (IMOCDIAG.EQ.3) THEN
         IF ((LDAYS.EQ.1).AND.(LMONTS.EQ.1)) THEN
            NANF=LDTDAYC
         ELSE
            NANF=NANF+1
         ENDIF
         NEND=0
         DO I=1,12
            NEND=NEND+MONLEN(I)
         ENDDO
         NEND=NEND*NDTDAY
      ELSEIF (IMOCDIAG.EQ.4) THEN
         
!HH every timestep   
         NANF=LDTDAYC
         NEND=LDTDAYC

!HH No diagnostic output
      ELSEIF (IMOCDIAG.EQ.0) THEN
         return
      ELSE
         STOP 'Stop in calc_moc due to wrong parameter.'
      ENDIF

!      WRITE(0,*)'in calc moc :',nanf,nend

    if (IMOCDIAG.NE.0) THEN

    do k=1,kep
       call gather_arr(wo(:,:,k),wo_g(:,:,k),p_io )
    enddo

    call gather_arr(depto,depto_g,p_io)
    call gather_arr(dlxp,dlxp_g,p_io)
    call gather_arr(dlyp,dlyp_g,p_io)
    call gather_arr(alat,alat_g,p_io)

!          compute moc 

    if (p_pe==p_io) then
    jbrei=3
    do i=2,ie_g-1
       do j=1,je_g
          if(weto_g(i,j,1).eq.1.)then
             !     1   suedpol
             !     180 nordpol
             lbrei=nint(90.+alat_g(i,j))
             lbrei=max(lbrei,1)
             lbrei=min(lbrei,180)

             do k=1,ke
                zlat=min(dlxp_g(i,j),dlyp_g(i,j))/(float(2*jbrei)*111111.)
                do l=-jbrei,jbrei
                   lbrei=nint(90.+alat_g(i,j)+float(l)*zlat)
                   lbrei=max(lbrei,1)
                   lbrei=min(lbrei,180)

                   tmerc(lbrei,1,k)=tmerc(lbrei,1,k)-dlxp_g(i,j)*dlyp_g(i,j) &
                        *wo_g(i,j,k)/float(2*jbrei+1)

                   if((ibek_g(i,j).le.5))then
                      tmerc(lbrei,2,k)=tmerc(lbrei,2,k)-dlxp_g(i,j)*dlyp_g(i,j)  &
                           *wo_g(i,j,k)/float(2*jbrei+1)

                   endif
                enddo
             enddo
          endif
       enddo
    enddo

    do lb=179,1,-1
       do k=1,kep
          tmerc(lb,1,k)=tmerc(lb+1,1,k)+tmerc(lb,1,k)
          tmerc(lb,2,k)=tmerc(lb+1,2,k)+tmerc(lb,2,k)
       enddo
    enddo

    IF (nanf.EQ.1) sum_tmerc(:,:,:)=0.

    sum_tmerc(:,:,:)=sum_tmerc(:,:,:)+(tmerc(:,:,:)-sum_tmerc(:,:,:))/nanf

    IF (nanf.EQ.nend) THEN
       i1=(lyears*10000)+(lmonts*100)+ldays 
       i2=100
       i4=180*kep
       do k=1,kep
          i3=NINT(tiestw(k))   
          write(io_ou_mocg)i1,i2,i3,180
          write(io_ou_mocg) sum_tmerc(:,1,k)
       enddo
    ENDIF
    IF (nanf.EQ.nend) THEN
       i1=(lyears*10000)+(lmonts*100)+ldays 
       i2=101
       i4=180*kep
       do k=1,kep
          i3=NINT(tiestw(k))   
          write(io_ou_moca)i1,i2,i3,180
          write(io_ou_moca) sum_tmerc(:,2,k)
       enddo
    ENDIF

    endif ! p_pe==p_pio
    endif      


  end subroutine calc_moc

subroutine calc_difi(idate)

  USE mo_kind
  USE mo_commoau1
  USE mo_units


!
!  calculates 2d-fields helpful for diagnostics of ocean fluxes
!  e.g. sea ice melting from below etc.
!  snapshots at the end of months
!  interface:
!  idate parameter for header of output
!

  real hice(ie,je),tdz(ie,je),sdz(ie,je)
  integer idate 
  integer (kind=sp) i41,i42,i43,i44
  real (kind=sp) ff_4(ie_g,je_g)
  integer i,j,k
  real       ff_g(ie_g,je_g)
!
!  initialize fields

  hice(:,:)=0.
  tdz(:,:)=0.
  sdz(:,:)=0.

  do k=1,ke
   do j=1,je
   do i=1,ie
   tdz(i,j)=tdz(i,j)+ddpo(i,j,k)*tho(i,j,k)*weto(i,j,k)
   sdz(i,j)=sdz(i,j)+ddpo(i,j,k)*sao(i,j,k)*weto(i,j,k)
   enddo
  enddo
  enddo

  do j=1,je
   do i=1,ie
    hice(i,j)=(rhosnwa*sicsno(i,j)+rhoicwa*sictho(i,j))*weto(i,j,1)
    sdz(i,j)=sdz(i,j)+(sao(i,j,1)*(zo(i,j)-rhosnwa*sicsno(i,j))       &
            -(sao(i,j,1)-5.)*rhoicwa*sictho(i,j))*weto(i,j,1)
    tdz(i,j)=tdz(i,j)+tho(i,j,1)*(zo(i,j)-rhosnwa*sicsno(i,j)         &
                                -rhoicwa*sictho(i,j))*weto(i,j,1)
   enddo
  enddo


      CALL gather_arr(hice(:,:),ff_g,p_io)
      i41=idate
      i42=101
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.eq.p_io) then
       ff_4=ff_g
       write(io_ou_difi)i41,i42,i43,i44
       write(io_ou_difi)ff_4
      endif

      CALL gather_arr(tdz,ff_g,p_io)
      i41=idate
      i42=102
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.eq.p_io) then
       ff_4=ff_g
       write(io_ou_difi)i41,i42,i43,i44
       write(io_ou_difi)ff_4
      endif

      CALL gather_arr(sdz,ff_g,p_io)
      i41=idate
      i42=103
      i43=-100
      i44=ie_g*je_g


      IF(p_pe.eq.p_io) then
       ff_4=ff_g
       write(io_ou_difi)i41,i42,i43,i44
       write(io_ou_difi)ff_4
      endif



  end  subroutine calc_difi

  


END MODULE mo_diagnosis
