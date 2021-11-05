MODULE mo_units
  ! ---------------------------------------------------------------------
  !
  !*    *common* *units*      - logical units for i/o of the ocean model.
  !
  !                           - included in all subroutines.
  !
  !
  !     s. legutke          *dkrz*           07.01.00
  !
  !*    variable     type        purpose.
  !     --------     ----        --------
  !     *io_stderr*  *integer*   logical unit for ocean std eror
  !     *io_stdout*  *integer*   logical unit for ocean stdout
  !     *io_in_octl* *integer*   logical unit for ocean namelist (ocectl)
  !     ** *integer*
  !
  ! ---------------------------------------------------------------------
  !

  USE mo_io_config, ONLY: setup_io_config, add_lun_reservation, &
       default_max_lun

  IMPLICIT NONE

  INTEGER :: io_stdout, io_stderr, io_in_octl, &
       io_in_z370, io_in_z380, &
       io_in_inii

  INTEGER :: io_ou_f125, io_ou_f126, io_ou_scha, io_ou_somi, io_ou_ray,   &
       io_ou_mtho, io_ou_msao,              &
       io_ou_emip, io_ou_mmzo, io_ou_flum, io_ou_mpem,              &
       io_ou_sict, io_ou_sico, io_ou_sicu, io_ou_sicv,              &
       io_ou_sics, io_ou_tafo, io_ou_fclo, io_ou_fpre,              &
       io_ou_fswr, io_ou_ftde, io_ou_qswo, io_ou_qlwo,              &
       io_ou_qseo, io_ou_qlao, io_ou_prec, io_ou_amld, io_ou_psiu,  &
       io_ou_weto, io_ou_gila, io_ou_giph, io_ou_dept, io_ou_zmld,  &
       io_ou_sictru, io_ou_sictrv,io_ou_txo,io_ou_tye,io_ou_mwgo,   &
       io_ou_mavo,io_ou_mdvo,io_ou_mwo,io_ou_dlxp,io_ou_dlyp,       &
       io_ou_wu10,io_ou_deute,io_ou_dlxu,io_ou_dlyu,io_ou_wtmi,io_ou_rinu,&
       io_ou_dqswo,io_ou_dqlwo,io_ou_dqseo,io_ou_dqlao,io_ou_dqtho, &
       io_ou_dqswi,io_ou_dqlwi,io_ou_dqsei,io_ou_dqlai,io_ou_dqthi, &
       io_ou_dticeo,io_ou_mbolx,io_ou_mboly,io_ou_f090,io_ou_tmceo, &
       io_ou_tmcdo,io_ou_ukomfl,io_ou_vkemfl,io_ou_rivrun,          &
       io_ou_aonhw,io_ou_aoshw,io_ou_aorhi,io_ou_aochi,io_ou_aofrw, &
       io_ou_aofri,io_ou_aotxw,io_ou_aotyw,io_ou_aotxi,io_ou_aotyi, &
       io_ou_aowsv,io_ou_mocg,io_ou_moca,io_ou_difi,io_ou_tsvar,    &
       io_ou_tsdesc,io_ou_tsunit,io_ou_alat,io_ou_alatu, &
       io_ou_alatv,io_ou_alon,io_ou_alonu,io_ou_alonv,io_ou_deuto,  &
       io_ou_amsue,io_ou_amsuo,io_ou_ddpo,io_ou_dduo,io_ou_ddue,    &
       io_ou_dlxv,io_ou_dlyv,io_ou_bek,io_ou_lpv

  INTEGER ::io_ou_ut0,io_ou_us0,io_ou_uu0 &
           ,io_ou_uv0,io_ou_uw0,io_ou_vt0 &
           ,io_ou_vs0,io_ou_vu0,io_ou_vv0 &
           ,io_ou_vw0,io_ou_wt0,io_ou_ws0 &
           ,io_ou_wu0,io_ou_wv0,io_ou_ww0,io_ou_zz0

CONTAINS

  SUBROUTINE setunits

!  USE mo_commo1, only :lgmdiag,ldiffdiag                           &
!                        ,lhfldiag,lgridinfo,lconvdiag,imean

    CALL setup_io_config
    default_max_lun=400
#ifdef __coupled
    io_stdout=8
    io_stderr=0
#else
    io_stdout=6
#ifdef __stdout0
    io_stdout=0
#endif

#endif

#ifdef MESSY
    io_stdout=6  ! standard oputput stream
#endif

    io_in_octl=10

    !     restart files
    io_in_z370=27
    CALL add_lun_reservation(27)
    io_in_z380=28
    CALL add_lun_reservation(28)

    !     set the output units

    io_ou_f125=125
    CALL add_lun_reservation(125)
    io_ou_f126=126
    CALL add_lun_reservation(126)

#ifdef KONVDIAG
    io_ou_f090=73
    CALL add_lun_reservation(73)
#endif
    io_ou_amld=73
    CALL add_lun_reservation(73)

    !     amoc
    io_ou_scha=47
    CALL add_lun_reservation(47)
    io_ou_somi=83
    CALL add_lun_reservation(83)
    io_ou_ray=92
    CALL add_lun_reservation(92)

    io_ou_mocg=75
    CALL add_lun_reservation(75)
    io_ou_moca=75
    CALL add_lun_reservation(75)


!    if (IMEAN.ne.0)then
       io_ou_mtho=71
       CALL add_lun_reservation(71)
       io_ou_msao=71
       io_ou_ukomfl=71
       io_ou_vkemfl=71
       io_ou_emip=71
       io_ou_mmzo=71
       io_ou_flum=71
       io_ou_mpem=71
       io_ou_sict=71
       io_ou_sico=71
       io_ou_sicu=71
       io_ou_sicv=71
       io_ou_sics=71
       io_ou_tafo=71
       io_ou_wu10=71
       io_ou_fclo=71
       io_ou_fpre=71
       io_ou_fswr=71
       io_ou_ftde=71
       io_ou_qswo=71
       io_ou_qlwo=71
       io_ou_qseo=71
       io_ou_qlao=71
       io_ou_prec=71
       io_ou_zmld=71

       io_ou_psiu=71
       io_ou_sictru=71
       io_ou_sictrv=71
       io_ou_txo=71
       io_ou_mwo=71
       io_ou_tye=71
       io_ou_rivrun=71

!       if (LGMDIAG) then
          io_ou_mwgo=71
          io_ou_mbolx=71
          io_ou_mboly=71
!       endif

!       if (LDIFFDIAG) then
          io_ou_mavo=71
          io_ou_mdvo=71
          io_ou_wtmi=71
          io_ou_rinu=71
!       endif

!       if (LHFLDIAG) then
          io_ou_dqswo=71
          io_ou_dqlwo=71
          io_ou_dqseo=71
          io_ou_dqlao=71
          io_ou_dqtho=71
          io_ou_dqswi=71
          io_ou_dqlwi=71
          io_ou_dqsei=71
          io_ou_dqlai=71
          io_ou_dqthi=71
          io_ou_dticeo=71
!       endif

          io_ou_ut0=71
          io_ou_us0=71
          io_ou_uu0=71
          io_ou_uv0=71
          io_ou_uw0=71
          io_ou_vt0=71
          io_ou_vs0=71
          io_ou_vu0=71
          io_ou_vv0=71
          io_ou_vw0=71
          io_ou_wt0=71
          io_ou_ws0=71
          io_ou_wu0=71
          io_ou_wv0=71
          io_ou_ww0=71
          io_ou_zz0=71


!    endif


#ifdef __coupled
    io_ou_aonhw=270
    CALL add_lun_reservation(270)
    io_ou_aoshw=270
    io_ou_aorhi=270
    io_ou_aochi=270
    io_ou_aofrw=270
    io_ou_aofri=270
    io_ou_aotxw=270
    io_ou_aotyw=270
    io_ou_aotxi=270
    io_ou_aotyi=270
    io_ou_aowsv=270
#endif /*__coupled*/



!    IF (LGRIDINFO) then
       io_ou_weto=72
       CALL add_lun_reservation(72)
       io_ou_gila=73
       CALL add_lun_reservation(73)
       io_ou_giph=73
       io_ou_dept=72
       io_ou_dlxp=72
       io_ou_dlyp=72
       io_ou_deute=72
       io_ou_deuto=72
       io_ou_dlxu=72
       io_ou_dlyu=72
       io_ou_dlxv=72
       io_ou_dlyv=72
       io_ou_ddpo=72
       io_ou_dduo=72
       io_ou_ddue=72
       io_ou_amsuo=72
       io_ou_amsue=72
       io_ou_alat=72
       io_ou_alatu=72
       io_ou_alatv=72
       io_ou_alon=72
       io_ou_alonu=72
       io_ou_alonv=72
       io_ou_bek=72
!    endif

!    if (LCONVDIAG) then
       io_ou_tmceo=71
       CALL add_lun_reservation(71)
       io_ou_tmcdo=71
!    endif

      io_ou_difi=74
      CALL add_lun_reservation(74)

      io_ou_tsvar=100
      CALL add_lun_reservation(100)
      io_ou_tsdesc=101
      CALL add_lun_reservation(101)
      io_ou_tsunit=102
      CALL add_lun_reservation(102)

  END SUBROUTINE setunits

END MODULE mo_units




