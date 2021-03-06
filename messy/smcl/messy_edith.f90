! **********************************************************************
!
! Author : Alexey Vlasov, KIT-IMK, 2016
!          Stefan Versick, KIT-IMK, 2016-
!          Sabine Barthlott, KIT-IMK, 2018-
!
! **********************************************************************

! **********************************************************************
MODULE messy_edith
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP
  
  IMPLICIT NONE
  PRIVATE
  
  INTEGER :: iamu, calc_iondrag, passive_nox, ed_solvar
  LOGICAL :: use_chem
  LOGICAL :: edith_phot
  REAL(dp) :: prandtlnumber
  REAL(dp) :: ed_o2_o3p, ed_o2_o1d, ed_o2_o1s

  PUBLIC :: DP, iamu, calc_iondrag, passive_nox, use_chem, prandtlnumber, ed_solvar, edith_phot &
            , ed_o2_o1d, ed_o2_o1s, ed_o2_o3p

!ka_sb_20180419+
  REAL(dp),  DIMENSION(:), POINTER,     PUBLIC :: kp_data   => NULL() ! necessary for ed_solvar=4
!ka_sb_20180419-
  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'edith'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.3'
 

  PUBLIC :: edith_read_nml_ctrl
  PUBLIC :: edith_calc_local_amu, edith_calc_local_cp, edith_calc_local_g
  PUBLIC :: edith_calc_nocooling, edith_calc_co2cooling
  PUBLIC :: edith_vdiff_mol_temp, edith_vdiff_mol_moist, edith_vdiff_mol_trac
  PUBLIC :: edith_vdiff_mol_wind, edith_vdiff_mol_fric, edith_iondrag
!ka_sv_20180702+
  PUBLIC :: edith_photsrc, edith_rcolsph
!ka_sv_20180702-
  

CONTAINS



  ! =========================================================================
  SUBROUTINE edith_read_nml_ctrl(status, iou)

    ! ------------------------------------------------------------------
    ! This routine is used to read the CTRL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy INTERFACE
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ iamu, calc_iondrag, passive_nox, use_chem, prandtlnumber, ed_solvar, edith_phot &
                    , ed_o2_o1d, ed_o2_o1s, ed_o2_o3p

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='edith_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR
    
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    ! ### ADD HERE DIAGNOSTIC OUPUT FOR LOG-FILE
    WRITE(*,*) 'iamu:',iamu
    write(*,*) 'iondrag:',calc_iondrag
    write(*,*) 'passive NOx:',passive_nox
    write(*,*) 'Use chemistry:',use_chem
    write(*,*) 'Used Prandtl-Number:',prandtlnumber
    write(*,*) 'Solvar:',ed_solvar
    write(*,*) 'Use O2/CO2 photolysis:', edith_phot
    write(*,*) '  Branching for O2-photolysis'
    write(*,*) '    O3P: ', ed_o2_o3p
    write(*,*) '    O1D: ', ed_o2_o1d
    write(*,*) '    O1S: ', ed_o2_o1s

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR
    
  END SUBROUTINE edith_read_nml_ctrl
  ! =========================================================================

! ----------------------------------------------------------------------

SUBROUTINE edith_calc_local_amu(iamu,apm1,amu)

  !USE messy_main_tracer_mem_bi, ONLY: ntrac_gp, ti_gp
  
! Variables for I/O of the subroutine
  INTEGER, INTENT(in) :: iamu          ! switch from namelist
  REAL(dp), INTENT(in) :: apm1         ! pressure in Pascal
  REAL(dp), INTENT(out) :: amu         ! Molecular mass of Air
  
  
! Local variables


  seliamu: SELECT CASE (iamu)
  CASE (1)
    amu=28.96_dp
  CASE (2)
     IF (apm1>=0.17) amu=28.96
     IF (apm1<0.17) amu=0.02405*(log10(apm1)**3)-0.1471*(log10(apm1)**2)+ &
           0.1568*(log10(apm1))+29.17
        ! formula based on MSIS data; coefficients calculated by matlab (cubic fit)
!      ENDDO
!    ENDDO
!  CASE (3)
!#ifdef MESSY
!    amu=0._dp
!        DO jt=1,ntrac
!          IF ((ti_gp(jt)%tp%ident%fullname=='O2').OR.(ti_gp(jt)%tp%ident%fullname=='N2').OR.   &
!             (ti_gp(jt)%tp%ident%fullname=='O3P').OR.(ti_gp(jt)%tp%ident%fullname=='O1D').OR.  &
!             (ti_gp(jt)%tp%ident%fullname=='H2').OR.(ti_gp(jt)%tp%ident%fullname=='H').OR.    &
!             (ti_gp(jt)%tp%ident%fullname=='CH4').OR.(ti_gp(jt)%tp%ident%fullname=='N2O').OR.  &
!             (ti_gp(jt)%tp%ident%fullname=='CO2').OR.(ti_gp(jt)%tp%ident%fullname=='CO')) THEN
!            amu=(xtm1(jl,jk,jt,krow) + &
!                    time_step_len*xtte(jl,jk,jt,krow)) *(ti_gp(jt)%tp%meta%cask_r(R_molarmass))+amu
!          ENDIF
!        ENDDO
!      ENDDO
!      if (apm1>=0.1) amu=amu+0.008*40.          ! correction for Argon
!    ENDDO
!#else
!    write(*,*) '/----------------\'
!    write(*,*) '|     WARNING    |'
!    write(*,*) '\----------------/'
!    write(*,*) ' You have chosen iamu=3 in ECHAM-namelist'
!    write(*,*) ' This is only supported for model runs which include MESSY'
!    write(*,*) ' Please choose another iamu and restart your model run'
!    CALL finish('physc','iamu=3 not supported in this model configuration')
!#endif
!  CASE (4)
!#ifdef MESSY
!    amu=0._dp
!    
!        DO jt=1,ntrac
!          IF (ti_gp(jt)%tp%meta%cask_i(I_vdiff) == ON) THEN
!            amu=(xtm1(jl,jk,jt,krow) + &
!               time_step_len*xtte(jl,jk,jt,krow)) *(ti_gp(jt)%tp%meta%cask_r(R_molarmass))+amu
!          ENDIF
!        ENDDO
!      ENDDO
!    ENDDO
!#else
!    write(*,*) '/----------------\'
!    write(*,*) '|     WARNING    |'
!    write(*,*) '\----------------/'
!    write(*,*) ' You have chosen iamu=4 in ECHAM-namelist'
!    write(*,*) ' This is only supported for model runs which include MESSY'
!    write(*,*) ' Please choose another iamu and restart your model run'
!    CALL finish('physc','iamu=4 not supported in this model configuration')
!#endif
  CASE default
     write(*,*) 'Select a proper iamu'
  END SELECT seliamu

END SUBROUTINE edith_calc_local_amu


! ----------------------------------------------------------------------

SUBROUTINE edith_calc_local_cp(cp,amu,twoatomicvmr,oneatomicvmr)

    USE messy_main_constants_mem, only : R_gas

    REAL(dp), INTENT(out) :: cp    ! specific heat at constant pressure
    REAL(dp), INTENT(in) :: amu    ! Molecular weight of air
    REAL(dp), INTENT(in) :: twoatomicvmr, oneatomicvmr   ! volume mixing ratios of one and two atomic gases
    
    cp = (twoatomicvmr*7._dp/2._dp*R_gas/amu*1000._dp)+(oneatomicvmr*5._dp/2._dp*R_gas/amu*1000._dp) &
             / (twoatomicvmr+oneatomicvmr)
             
    
END SUBROUTINE edith_calc_local_cp

! ----------------------------------------------------------------------

SUBROUTINE edith_calc_local_g(g,lat,height)

    USE messy_main_constants_mem, ONLY: DTR
! WELMEC-Formel
    REAL(dp), INTENT(out) :: g    ! gravity acceleration
    REAL(dp), INTENT(in) :: lat, height   ! latitude (in arc) and altitude above sea level (in m)
    ! lat=sinlat
    g=(1._dp+0.0053024*sin(lat*DTR)*sin(lat*DTR)-0.0000058*sin(2._dp*lat*DTR)*sin(2._dp*lat*DTR)) &
        * 9.780318-0.000003085*height
    
END SUBROUTINE edith_calc_local_g

! ----------------------------------------------------------------------

SUBROUTINE edith_calc_nocooling(NO, O_t0, pmid, temp, cp, NOcooling)

! op_pj_20180713+: illegal use of other submodel
!!$  USE messy_mesoenergy, ONLY: vmr_to_nd
  USE messy_main_tools, ONLY: vmr_to_nd
! op_pj_20180713-

  REAL(dp), intent(in) :: NO
  REAL(dp), intent(in) :: O_t0
  REAL(dp), intent(in) :: pmid
  REAL(dp), intent(in) :: temp
  REAL(dp), intent(in) :: cp
  REAL(dp), intent(out) :: NOcooling
  REAL(dp) :: n_o

    
    call vmr_to_nd(n_o,max(O_t0,1.e-30),pmid,temp)
    n_o=n_o/1.e06;
       !tendency budget
    NOcooling = -7.965e5_dp * NO&
                      * 13.3 / ( 1._dp + 13.3_dp / 6.5e-11_dp / n_o)&
                      * exp(-2699._dp/temp)&
                      / cp

END SUBROUTINE

! ----------------------------------------------------------------------

SUBROUTINE edith_calc_co2cooling(nproma,nbdim,nlev,CO2, O3, O2, O, N2, amu, cp, &
           apm1,tm1,cossza_2d, &
           CO2nlte, CO2nir, &
           amatsave, &
           bmatsave,alsave,co2colxsave, &
           co2xsave)

! nbdim?

  USE messy_edith_rad_nlte,          ONLY: nlte_heating
  
  INTEGER, intent(in) :: nbdim
  REAL(dp), intent(in) :: CO2(:,:), O3(:,:), O2(:,:), O(:,:), N2(:,:)
  REAL(dp), intent(inout) :: tm1(:,:),apm1(:,:)
  REAL(dp), intent(in) :: amu(:,:), cp(:,:)
  REAL(dp), intent(in) :: cossza_2d(:)
  INTEGER, intent(in) :: nproma,nlev
  REAL(dp), intent(out) :: CO2nlte(:,:), CO2nir(:,:)
  REAL(dp), ALLOCATABLE :: CO2l(:,:), O3l(:,:), O2l(:,:), Ol(:,:), N2l(:,:)
  
  REAL(dp), INTENT(inout) :: amatsave(nbdim,43,9)
  REAL(dp), INTENT(inout) :: bmatsave(nbdim,43,9)
  REAL(dp), INTENT(inout) :: alsave(nbdim,17)
  REAL(dp), INTENT(inout) :: co2colxsave(nbdim,59)
  REAL(dp), INTENT(inout) :: co2xsave(nbdim,59)
  
  
  REAL(dp), ALLOCATABLE :: zhnlte(:,:), zhnir(:,:), zxm1(:,:)
  
  REAL(dp) :: zp2, zx2, zp1, zx1
  
  ALLOCATE(CO2l(nproma,nlev))
  ALLOCATE(O3l(nproma,nlev))
  ALLOCATE(O2l(nproma,nlev))
  ALLOCATE(Ol(nproma,nlev))
  ALLOCATE(N2l(nproma,nlev))
  ALLOCATE(zhnlte(nproma,nlev))
  ALLOCATE(zhnir(nproma,nlev))
  ALLOCATE(zxm1(nproma,nlev))
  
  CO2NLTE(:,:)=0._dp
  CO2NIR(:,:)=0._dp
  
  CO2l=max(CO2,1.e-30)*amu/44.
  O3l=max(O3,1.e-30)*amu/48.
  O2l=max(O2,1.e-30)*amu/32.
  Ol=max(O,1.e-30)*amu/16.
  N2l=max(N2,1.e-30)*amu/28.
  
  CALL nlte_heating(nproma,nbdim,nlev,apm1(:,:),tm1(:,:),    &
           CO2l(:,:),O3l(:,:),O2l(:,:),Ol(:,:),N2l(:,:),amu(:,:),cp(:,:), &
            cossza_2d(:), &
            amatsave(:,:,:), &
           bmatsave(:,:,:),alsave(:,:),co2colxsave(:,:), &
           co2xsave(:,:),zhnlte(:,:),zhnir(:,:))
! use Victor's scheme for p<=zp2 (Pa)
     zp2=2._dp
     zx2=log(100000._dp/zp2)

     where(apm1(1:nproma,:)<=2._dp)
          CO2NLTE(:,:)=zhnlte(:,:)
     end where
    
    ! merge between zp1 and zp2 (Pa)
     zp1=10._dp
     zx1=log(100000._dp/zp1)
     zxm1(1:nproma,:)=log(100000._dp/apm1(1:nproma,:))
     where(zp1>apm1(:,:).and.apm1(:,:)>zp2)
        CO2NLTE(:,:)=(zhnlte(:,:)*(zxm1(:,:)-zx1))/(zx2-zx1)
     end where
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! add nIR CO2 heating for p<=zpnir2
     zp2=3000._dp
     zx2=log(100000._dp/zp2)
     where(apm1(:,:)<=zp2)
        CO2NIR(:,:)=zhnir(:,:)
     end where

     ! merge between zp1 and zp2 (Pa)
     zp1=7000._dp
     zx1=log(100000._dp/zp1)
     where(zp1>apm1(:,:).and.apm1(:,:)>zp2)
        CO2NIR(:,:)=zhnir(:,:)*   &
             & (zxm1(:,:)-zx1)/(zx2-zx1)
     end where

  DEALLOCATE(CO2l)
  DEALLOCATE(O2l)
  DEALLOCATE(Ol)
  DEALLOCATE(O3l)
  DEALLOCATE(N2l)
  DEALLOCATE(zhnlte)
  DEALLOCATE(zhnir)
  DEALLOCATE(zxm1)

END SUBROUTINE



SUBROUTINE edith_vdiff_mol_temp( kproma, kbdim, klev, klevm1 &
                    , ptm1, ztmst &
                    , zalpha, zmalpha &
                    , zla, zlb, zlc &
                    , ptte, status)
! temperature part of the old vdiff_mol routine

! INPUT
INTEGER,     INTENT(IN) :: kproma, kbdim, klev, klevm1
REAL(dp),    INTENT(IN) :: ztmst
REAL(dp),    INTENT(IN) :: ptm1(kbdim,klev)
REAL(dp),    INTENT(IN) :: zalpha,zmalpha
REAL(dp),    INTENT(IN) :: zla(kbdim,klev),zlb(kbdim,klev), zlc(kbdim,klev)
! OUTPUT
INTEGER, INTENT(OUT) :: status ! error status
REAL(dp), INTENT(OUT) :: ptte(kbdim,klev)  ! temperature tendency

!INTERNAL
INTEGER :: jl, jk
INTEGER :: mess, messm
REAL(dp)    :: zut(kbdim,klev)
REAL(dp)    :: zr(kbdim,klev)
REAL(dp)    :: zbet(kbdim,klev),zgam(kbdim,klev)
REAL(dp)    :: za(kbdim,klev),zb(kbdim,klev),zc(kbdim,klev)

  status=1
!** 1: TEMPERATURE EQUATION
!   -----------------------


!** SET TRI-DIAGONAL MATRIX: zc IS UPPER  DIAGONAL
!                            zb IS MIDDLE DIAGONAL
!                            za IS LOWER  DIAGONAL

  za(1:kproma,2:klev)   =   -zalpha*zla(1:kproma,2:klev)
  zb(1:kproma,:)        = 1._dp-zalpha*zlb(1:kproma,:)
  zc(1:kproma,1:klevm1) =   -zalpha*zlc(1:kproma,1:klevm1)

!** SET RIGHT-HAND SIDE zr

  DO jk=2,klevm1
    zr(1:kproma,jk)=    zmalpha*zla(1:kproma,jk) *ptm1(1:kproma,jk-1) + &
             (1._dp+zmalpha*zlb(1:kproma,jk))*ptm1(1:kproma,jk  ) + &
                 zmalpha*zlc(1:kproma,jk) *ptm1(1:kproma,jk+1)
  ENDDO

  zr(1:kproma,1   )=(1._dp+zmalpha*zlb(1:kproma,1   ))*ptm1(1:kproma,1     ) + &
                 zmalpha*zlc(1:kproma,1   ) *ptm1(1:kproma,2     )
  zr(1:kproma,klev)=    zmalpha*zla(1:kproma,klev) *ptm1(1:kproma,klevm1) + &
             (1._dp+zmalpha*zlb(1:kproma,klev))*ptm1(1:kproma,klev  )

!** START OF MATRIX SOLVER (zut IS T(t+dt))

  mess=0
  messm=0
  DO jl=1,kproma
    IF (zb(jl,1) .EQ. 0._dp) mess=1
    zbet(jl,1)=zb(jl,1)
    zut(jl,1)=zr(jl,1)/zbet(jl,1)
!** DECOMPOSITION AND FORWARD SUBSTITUTION

    DO jk=2,klev
      zgam(jl,jk)=zc(jl,jk-1)/zbet(jl,jk-1)
      zbet(jl,jk)=zb(jl,jk)-za(jl,jk)*zgam(jl,jk)
      IF (zbet(jl,jk) .EQ. 0._dp) messm=1
      zut(jl,jk)=(zr(jl,jk)-za(jl,jk)*zut(jl,jk-1))/zbet(jl,jk)
    ENDDO

!** BACKSUBSTITUTION

    DO jk=klevm1,1,-1
      zut(jl,jk)=zut(jl,jk)-zgam(jl,jk+1)*zut(jl,jk+1)
    ENDDO
  ENDDO
!ka_av_20160105+
! changing to the messy standards, old error message is still displayed
!  if (mess.eq.1) CALL finish('vdiff_mol','(1) Solver problem')
!  if (messm.eq.1) CALL finish('vdiff_mol','(2) Solver problem')
    if (mess.eq.1) then
        write(*,*)'vdiff_mol(1) Solver problem'
        return
    end if
    if (messm.eq.1)then
        write(*,*)'vdiff_mol(2) Solver problem'
        return
    end if
!ka_av_20160105-      
  ptte=(zut-ptm1)/ztmst
  status=0

END SUBROUTINE ! edith_vdiff_mol_temp


SUBROUTINE edith_vdiff_mol_moist( kproma, kbdim, klev, klevm1 &
                    , pqm1, ztmst &
                    , zalpha, zmalpha &
                    , zla, zlb, zlc &
                    , pqte, status)
                    
! INPUT
INTEGER,     INTENT(IN) :: kproma, kbdim, klev, klevm1
REAL(dp),    INTENT(IN) :: ztmst
REAL(dp),    INTENT(IN) :: pqm1(kbdim,klev)
REAL(dp),    INTENT(IN) :: zmalpha, zalpha
REAL(dp),    INTENT(IN) :: zla(kbdim,klev),zlb(kbdim,klev), zlc(kbdim,klev)
! OUTPUT
INTEGER, INTENT(OUT) :: status ! error status
REAL(dp), INTENT(OUT) :: pqte(kbdim,klev)  ! moisture tendency

!INTERNAL
INTEGER :: jl, jk
INTEGER :: mess, messm
REAL(dp)    :: zut(kbdim,klev)
REAL(dp)    :: zr(kbdim,klev)
REAL(dp)    :: zbet(kbdim,klev),zgam(kbdim,klev)
REAL(dp)    :: za(kbdim,klev),zb(kbdim,klev),zc(kbdim,klev)

  status=1
!** 2: HUMIDITY EQUATION
!   --------------------

!** SET TRI-DIAGONAL MATRIX: zc IS UPPER  DIAGONAL
!                            zb IS MIDDLE DIAGONAL
!                            za IS LOWER  DIAGONAL

  za(1:kproma,2:klev)   =   -zalpha*zla(1:kproma,2:klev)
  zb(1:kproma,:)        = 1._dp-zalpha*zlb(1:kproma,:)
  zc(1:kproma,1:klevm1) =   -zalpha*zlc(1:kproma,1:klevm1)

!** SET RIGHT-HAND SIDE zr

  DO jk=2,klevm1
    zr(1:kproma,jk)=    zmalpha*zla(1:kproma,jk) *pqm1(1:kproma,jk-1) + &
             (1._dp+zmalpha*zlb(1:kproma,jk))*pqm1(1:kproma,jk  ) + &
                 zmalpha*zlc(1:kproma,jk) *pqm1(1:kproma,jk+1)
  ENDDO

  zr(1:kproma,1   )=(1._dp+zmalpha*zlb(1:kproma,1   ))*pqm1(1:kproma,1     ) + &
                 zmalpha*zlc(1:kproma,1   ) *pqm1(1:kproma,2     )
  zr(1:kproma,klev)=    zmalpha*zla(1:kproma,klev) *pqm1(1:kproma,klevm1) + &
             (1._dp+zmalpha*zlb(1:kproma,klev))*pqm1(1:kproma,klev  )

!** START OF MATRIX SOLVER (zut IS q(t+dt))

  mess=0
  messm=0
  DO jl=1,kproma
    IF (zb(jl,1) .EQ. 0._dp) mess=1
    zbet(jl,1)=zb(jl,1)
    zut(jl,1)=zr(jl,1)/zbet(jl,1)

!** DECOMPOSITION AND FORWARD SUBSTITUTION

    DO jk=2,klev
      zgam(jl,jk)=zc(jl,jk-1)/zbet(jl,jk-1)
      zbet(jl,jk)=zb(jl,jk)-za(jl,jk)*zgam(jl,jk)
      IF (zbet(jl,jk) .EQ. 0._dp) messm=1
      zut(jl,jk)=(zr(jl,jk)-za(jl,jk)*zut(jl,jk-1))/zbet(jl,jk)
    ENDDO

!** BACKSUBSTITUTION

    DO jk=klevm1,1,-1
      zut(jl,jk)=zut(jl,jk)-zgam(jl,jk+1)*zut(jl,jk+1)
    ENDDO
  ENDDO
!ka_av_20160105+
! changing to the messy standards, old error message is still displayed
!  if (mess.eq.1) CALL finish('vdiff_mol','(3) Solver problem')
!  if (messm.eq.1) CALL finish('vdiff_mol','(4) Solver problem')
    if (mess.eq.1)then
       write(*,*)'vdiff_mol(3) Solver problem'
       return
    end if
    if (messm.eq.1)then
       write(*,*)'vdiff_mol(4) Solver problem'
       return
    end if 
!ka_av_20160105- 
  pqte=(zut-pqm1)/ztmst
  status=0
                    
END SUBROUTINE ! edith_vdiff_mol_moist



SUBROUTINE edith_vdiff_mol_trac( kproma, kbdim, klev, klevm1 &
                    , ptm1, zpxtm1, ztmst &
                    , zalpha, zmalpha &
                    , zmah, zrhok, zgvh, zrhoh, zth, zrstar, zmm, zexpan, zhh, zrdp &
                    , zutt, status)
                    
! INPUT
INTEGER,     INTENT(IN) :: kproma, kbdim, klev, klevm1
REAL(dp),    INTENT(IN) :: ztmst
REAL(dp),    INTENT(IN) :: ptm1(kbdim,klev)
REAL(dp),    INTENT(IN) :: zmalpha,zalpha

REAL(dp), INTENT(IN)    :: zmah(kbdim,klev),zrhok(kbdim,klev)
REAL(dp), INTENT(IN)    :: zgvh(kbdim,klev),zrhoh(kbdim,klev),zth(kbdim,klev)
REAL(dp), INTENT(IN)   :: zrstar,zmm,zexpan
REAL(dp), INTENT(IN)    :: zpxtm1(kbdim,klev)
REAL(dp), INTENT(IN)    :: zhh(kbdim,klev),zrdp(kbdim,klev)
! OUTPUT
INTEGER, INTENT(OUT) :: status ! error status
REAL(dp), INTENT(OUT) :: zutt(kbdim,klev)

!INTERNAL
INTEGER :: jl, jk
INTEGER :: mess, messm
REAL(dp)    :: zrt(kbdim,klev)
REAL(dp)    :: zbet(kbdim,klev),zgam(kbdim,klev)
REAL(dp)    :: zlatr(kbdim,klev),zlbtr(kbdim,klev),zlctr(kbdim,klev)
REAL(dp)    :: zatr(kbdim,klev),zbtr(kbdim,klev),zctr(kbdim,klev)
REAL(dp)    :: zrhowd(kbdim,klev+1)
        
        
        status=1            
!** 3: TRACER EQUATION (INCLUDING DIFFUSIVE SEPARATION)
!   ---------------------------------------------------

!** DEFINE DRIFT VERTICAL VELOCITY * DENSITY

  
    DO jk=2,klev
      zrhowd(1:kproma,jk)=zrhok(1:kproma,jk)*zgvh(1:kproma,jk)/(zrstar*zth(1:kproma,jk))*        &
                      (zmah(1:kproma,jk)-zmm+zexpan*zrhoh(1:kproma,jk)*zrstar* &
                       (ptm1(1:kproma,jk-1)-ptm1(1:kproma,jk))*zhh(1:kproma,jk))
    ENDDO
  
   DO jk=2,klevm1
      zlatr(1:kproma,jk)=-zrdp(1:kproma,jk)*( zgvh(1:kproma,jk  )*zrhoh(1:kproma,jk  )*       &
         zrhok(1:kproma,jk)*zhh(1:kproma,jk  )  +0.5*zrhowd(1:kproma,jk))
      zlbtr(1:kproma,jk)=-zrdp(1:kproma,jk)*(-zgvh(1:kproma,jk  )*zrhoh(1:kproma,jk  )*       &
         zrhok(1:kproma,jk)*zhh(1:kproma,jk  )  -zgvh(1:kproma,jk+1)*zrhoh(1:kproma,jk+1)*  &
         zrhok(1:kproma,jk+1)*zhh(1:kproma,jk+1)  +0.5*zrhowd(1:kproma,jk)-     &
         0.5*zrhowd(1:kproma,jk+1))
      zlctr(1:kproma,jk)=-zrdp(1:kproma,jk)*( zgvh(1:kproma,jk+1)*zrhoh(1:kproma,jk+1)*       &
         zrhok(1:kproma,jk+1)*zhh(1:kproma,jk+1)  -0.5*zrhowd(1:kproma,jk+1))
    ENDDO



    zlbtr(1:kproma,1   )=-zrdp(1:kproma,1   )*(-zgvh(1:kproma,2   )*zrhoh(1:kproma,2   )*     &
         zrhok(1:kproma,2   )*zhh(1:kproma,2   )  -0.5*zrhowd(1:kproma,2   ))
    zlctr(1:kproma,1   )=-zrdp(1:kproma,1   )*( zgvh(1:kproma,2   )*zrhoh(1:kproma,2   )*     &
         zrhok(1:kproma,2   )*zhh(1:kproma,2   )  -0.5*zrhowd(1:kproma,2   ))
    zlatr(1:kproma,klev)=-zrdp(1:kproma,klev)*( zgvh(1:kproma,klev)*zrhoh(1:kproma,klev)*     &
         zrhok(1:kproma,klev)*zhh(1:kproma,klev)  +0.5*zrhowd(1:kproma,klev))
    zlbtr(1:kproma,klev)=-zrdp(1:kproma,klev)*(-zgvh(1:kproma,klev)*zrhoh(1:kproma,klev)*     &
         zrhok(1:kproma,klev)*zhh(1:kproma,klev)  +0.5*zrhowd(1:kproma,klev))


!** SET TRI-DIAGONAL MATRIX: zctr IS UPPER  DIAGONAL
!                            zbtr IS MIDDLE DIAGONAL
!                            zatr IS LOWER  DIAGONAL


    zatr(1:kproma,2:klev) =  -zalpha*zlatr(1:kproma,2:klev)
    zbtr(1:kproma,:) =1._dp-zalpha*zlbtr(1:kproma,:)
    zctr(1:kproma,1:klevm1) =  -zalpha*zlctr(1:kproma,1:klevm1)

!** SET RIGHT-HAND SIDE zrt

    DO jk=2,klevm1
      zrt(1:kproma,jk)=    zmalpha*zlatr(1:kproma,jk) *zpxtm1(1:kproma,jk-1) + &
                  (1._dp+zmalpha*zlbtr(1:kproma,jk))*zpxtm1(1:kproma,jk  ) + &
                      zmalpha*zlctr(1:kproma,jk) *zpxtm1(1:kproma,jk+1)
    ENDDO

    zrt(1:kproma,1   )=(1._dp+zmalpha*zlbtr(1:kproma,1   ))*zpxtm1(1:kproma,1     ) + &
                      zmalpha*zlctr(1:kproma,1   ) *zpxtm1(1:kproma,2     )
    zrt(1:kproma,klev)=    zmalpha*zlatr(1:kproma,klev) *zpxtm1(1:kproma,klevm1) + &
                  (1._dp+zmalpha*zlbtr(1:kproma,klev))*zpxtm1(1:kproma,klev  )


!** START OF MATRIX SOLVER (zutt IS tracer(t+dt))

  mess=0
  messm=0
  
    DO jl=1,kproma
      IF (zbtr(jl,1) .EQ. 0.0_dp) mess=1
      zbet(jl,1)=zbtr(jl,1)
      zutt(jl,1)=zrt(jl,1)/zbet(jl,1)

!** DECOMPOSITION AND FORWARD SUBSTITUTION

      DO jk=2,klev
        zgam(jl,jk)=zctr(jl,jk-1)/zbet(jl,jk-1)
        zbet(jl,jk)=zbtr(jl,jk)-zatr(jl,jk)*zgam(jl,jk)
        IF (zbet(jl,jk) .EQ. 0.0_dp) messm=1
        zutt(jl,jk)=(zrt(jl,jk)-zatr(jl,jk)*zutt(jl,jk-1))/zbet(jl,jk)
      ENDDO

!** BACKSUBSTITUTION

      DO jk=klevm1,1,-1
        zutt(jl,jk)=zutt(jl,jk)-zgam(jl,jk+1)*zutt(jl,jk+1)
      ENDDO
    ENDDO
  
!ka_av_20160105+
! changing to the messy standards, old error message is still displayed
!  if (mess.eq.1) CALL finish('vdiff_mol','(5) Solver problem')
!  if (messm.eq.1) CALL finish('vdiff_mol','(6) Solver problem')
    if (mess.eq.1)then
        write(*,*)'vdiff_mol(5) Solver problem'
        return
    end if
    if (messm.eq.1)then
        write(*,*)'vdiff_mol(6) Solver problem'
        return
    end if 
!ka_av_20160105- 

  ! Make sure that results ar non negative

    DO jl=1,kproma
      DO jk=1,klev
        zutt(jl,jk)=MAX(zutt(jl,jk),1.e-30_dp)
      ENDDO
    ENDDO

    status=0
 

END SUBROUTINE ! edith_vdiff_mol_trac



SUBROUTINE edith_vdiff_mol_wind( kproma, kbdim, klev, klevm1 &
                    , pum1, pvm1, ztmst &
                    , zalpha, zmalpha &
                    , zrdp, zgmurhoh, zhh &
                    , pvom, pvol, status)
                    
! INPUT
INTEGER,     INTENT(IN) :: kproma, kbdim, klev, klevm1
REAL(dp),    INTENT(IN) :: ztmst
REAL(dp),    INTENT(IN) :: pum1(kbdim,klev), pvm1(kbdim,klev)
REAL(dp),    INTENT(IN) :: zmalpha,zalpha

REAL(dp), INTENT(IN)    :: zhh(kbdim,klev),zrdp(kbdim,klev), zgmurhoh(kbdim, klev)
! OUTPUT
INTEGER, INTENT(OUT) :: status ! error status
REAL(dp), INTENT(OUT) :: pvom(kbdim,klev), pvol(kbdim, klev)

!INTERNAL
INTEGER :: jl, jk
INTEGER :: mess, messm
REAL(dp)    :: zr(kbdim,klev)
REAL(dp)    :: zbet(kbdim,klev),zgam(kbdim,klev), zut(kbdim,klev)
REAL(dp)    :: zla(kbdim,klev),zlb(kbdim,klev),zlc(kbdim,klev)
REAL(dp)    :: za(kbdim,klev),zb(kbdim,klev),zc(kbdim,klev)


  status=1
!** 4: U-WIND EQUATION
!   ------------------

  DO jk=2,klevm1
    zla(1:kproma,jk)=-zrdp(1:kproma,jk)*(zgmurhoh(1:kproma,jk  )*zhh(1:kproma,jk  ))
    zlb(1:kproma,jk)= zrdp(1:kproma,jk)*(zgmurhoh(1:kproma,jk  )*zhh(1:kproma,jk  ) + &
                           zgmurhoh(1:kproma,jk+1)*zhh(1:kproma,jk+1))
    zlc(1:kproma,jk)=-zrdp(1:kproma,jk)*(zgmurhoh(1:kproma,jk+1)*zhh(1:kproma,jk+1))
  ENDDO

  zlb(1:kproma,1   )= zrdp(1:kproma,1   )*zgmurhoh(1:kproma,2   )*zhh(1:kproma,2   )
  zlc(1:kproma,1   )=-zrdp(1:kproma,1   )*zgmurhoh(1:kproma,2   )*zhh(1:kproma,2   )

  zla(1:kproma,klev)=-zrdp(1:kproma,klev)*zgmurhoh(1:kproma,klev)*zhh(1:kproma,klev)
  zlb(1:kproma,klev)= zrdp(1:kproma,klev)*zgmurhoh(1:kproma,klev)*zhh(1:kproma,klev)

!** SET TRI-DIAGONAL MATRIX: zc IS UPPER  DIAGONAL
!                            zb IS MIDDLE DIAGONAL
!                            za IS LOWER  DIAGONAL

  za(1:kproma,2:klev)   =   -zalpha*zla(1:kproma,2:klev)
  zb(1:kproma,:)        = 1._dp-zalpha*zlb(1:kproma,:)
  zc(1:kproma,1:klevm1) =   -zalpha*zlc(1:kproma,1:klevm1)


!** SET RIGHT-HAND SIDE zr

  DO jk=2,klevm1
    zr(1:kproma,jk)=    zmalpha*zla(1:kproma,jk) *pum1(1:kproma,jk-1) + &
             (1._dp+zmalpha*zlb(1:kproma,jk))*pum1(1:kproma,jk  ) + &
                 zmalpha*zlc(1:kproma,jk) *pum1(1:kproma,jk+1)
  ENDDO

  zr(1:kproma,1   )=(1._dp+zmalpha*zlb(1:kproma,1   ))*pum1(1:kproma,1     ) + &
                 zmalpha*zlc(1:kproma,1   ) *pum1(1:kproma,2     )
  zr(1:kproma,klev)=    zmalpha*zla(1:kproma,klev) *pum1(1:kproma,klevm1) + &
             (1._dp+zmalpha*zlb(1:kproma,klev))*pum1(1:kproma,klev  )

!** START OF MATRIX SOLVER (zut IS u(t+dt))

  mess=0
  messm=0
  DO jl=1,kproma
    IF (zb(jl,1) .EQ. 0._dp) mess=1
    zbet(jl,1)=zb(jl,1)
    zut(jl,1)=zr(jl,1)/zbet(jl,1)

!** DECOMPOSITION AND FORWARD SUBSTITUTION

    DO jk=2,klev
      zgam(jl,jk)=zc(jl,jk-1)/zbet(jl,jk-1)
      zbet(jl,jk)=zb(jl,jk)-za(jl,jk)*zgam(jl,jk)
      IF (zbet(jl,jk) .EQ. 0._dp) messm=1
      zut(jl,jk)=(zr(jl,jk)-za(jl,jk)*zut(jl,jk-1))/zbet(jl,jk)
    ENDDO

!** BACKSUBSTITUTION

    DO jk=klevm1,1,-1
      zut(jl,jk)=zut(jl,jk)-zgam(jl,jk+1)*zut(jl,jk+1)
    ENDDO
  ENDDO
!ka_av_20160105+
! changing to the messy standards, old error message is still displayed
!  if (mess.eq.1) CALL finish('vdiff_mol','(7) Solver problem')
!  if (messm.eq.1) CALL finish('vdiff_mol','(8) Solver problem')
    if (mess.eq.1)then
       write(*,*)'vdiff_mol(7) Solver problem'
       return
    end if
    if (messm.eq.1)then
       write(*,*)'vdiff_mol(8) Solver problem'
       return
    end if 
!ka_av_20160105- 

  pvom=(zut-pum1)/ztmst

!** 5: V-WIND EQUATION
!   ------------------

!** SET RIGHT-HAND SIDE zr

  DO jk=2,klevm1
    zr(1:kproma,jk)=    zmalpha*zla(1:kproma,jk) *pvm1(1:kproma,jk-1) + &
             (1._dp+zmalpha*zlb(1:kproma,jk))*pvm1(1:kproma,jk  ) + &
                 zmalpha*zlc(1:kproma,jk) *pvm1(1:kproma,jk+1)
  ENDDO

  zr(1:kproma,1   )=(1._dp+zmalpha*zlb(1:kproma,1   ))*pvm1(1:kproma,1     ) + &
                 zmalpha*zlc(1:kproma,1   ) *pvm1(1:kproma,2     )
  zr(1:kproma,klev)=    zmalpha*zla(1:kproma,klev) *pvm1(1:kproma,klevm1) + &
             (1._dp+zmalpha*zlb(1:kproma,klev))*pvm1(1:kproma,klev  )

!** START OF MATRIX SOLVER (zut IS v(t+dt))
  mess=0
  messm=0
  DO jl=1,kproma
    IF (zb(jl,1) .EQ. 0._dp) mess=1
    zbet(jl,1)=zb(jl,1)
    zut(jl,1)=zr(jl,1)/zbet(jl,1)

!** DECOMPOSITION AND FORWARD SUBSTITUTION

    DO jk=2,klev
      zgam(jl,jk)=zc(jl,jk-1)/zbet(jl,jk-1)
      zbet(jl,jk)=zb(jl,jk)-za(jl,jk)*zgam(jl,jk)
      IF (zbet(jl,jk) .EQ. 0._dp) messm=1
      zut(jl,jk)=(zr(jl,jk)-za(jl,jk)*zut(jl,jk-1))/zbet(jl,jk)
    ENDDO

!** BACKSUBSTITUTION

    DO jk=klevm1,1,-1
      zut(jl,jk)=zut(jl,jk)-zgam(jl,jk+1)*zut(jl,jk+1)
    ENDDO
  ENDDO
!ka_av_20160105+
! changing to the messy standards, old error message is still displayed
!  if (mess.eq.1) CALL finish('vdiff_mol','(9) Solver problem')
!  if (messm.eq.1) CALL finish('vdiff_mol','(10) Solver problem')
    if (mess.eq.1) then
       write(*,*)'vdiff_mol(9) Solver problem'
       return
    end if
    if (messm.eq.1)then
       write(*,*)'vdiff_mol(10) Solver problem'
       return
    end if 
!ka_av_20160105- 
  pvol=(zut-pvm1)/ztmst
  status=0
  
END SUBROUTINE ! edith_vdiff_mol_wind


SUBROUTINE edith_vdiff_mol_fric( kproma, kbdim, klev, klevm1 &
                    , papm1, pum1, pvm1 &
                    , zgmurhoh, grav, cp &
                    , ptte_out, status)

! INPUT
INTEGER,     INTENT(IN) :: kproma, kbdim, klev, klevm1
REAL(dp),    INTENT(IN) :: pum1(kbdim,klev), pvm1(kbdim,klev), papm1(kbdim,klev)
REAL(dp), INTENT(IN)    :: grav(kbdim,klev),cp(kbdim,klev), zgmurhoh(kbdim, klev)

! OUTPUT
INTEGER, INTENT(OUT) :: status ! error status
REAL(dp), INTENT(OUT) :: ptte_out(kbdim,klev)

!INTERNAL
INTEGER :: jl, jk, jkk
REAL(dp)    :: zrdz,zdudz, zdvdz
REAL(dp)    :: zcoef
         
  status=1           
!** 6: FRICTIONAL HEATING
!   ---------------------

  DO jl=1,kproma
    DO jk=2,klevm1
      jkk=klev-jk+1
      zrdz=1._dp/(papm1(jl,jk-1)-papm1(jl,jk+1))
      zdudz=(pum1(jl,jk-1)-pum1(jl,jk+1))*zrdz
      zdvdz=(pvm1(jl,jk-1)-pvm1(jl,jk+1))*zrdz
      zcoef=0.5_dp*(zgmurhoh(jl,jk)+zgmurhoh(jl,jk+1))*grav(jl,jkk)/cp(jl,jk)
      ptte_out(jl,jk)=zcoef*(zdudz*zdudz+zdvdz*zdvdz)
    ENDDO
  ENDDO

  DO jl=1,kproma
    zrdz=1._dp/(papm1(jl,1)-papm1(jl,2))
    zdudz=(pum1(jl,1)-pum1(jl,2))*zrdz
    zdvdz=(pvm1(jl,1)-pvm1(jl,2))*zrdz
    zcoef=zgmurhoh(jl,2)*grav(jl,klev)/cp(jl,1)
    ptte_out(jl,1)=zcoef*(zdudz*zdudz+zdvdz*zdvdz)
  ENDDO

  DO jl=1,kproma
    zrdz=1._dp/(papm1(jl,klev)-papm1(jl,klevm1))
    zdudz=(pum1(jl,klev)-pum1(jl,klevm1))*zrdz
    zdvdz=(pvm1(jl,klev)-pvm1(jl,klevm1))*zrdz
    zcoef=zgmurhoh(jl,klev)*grav(jl,1)/cp(jl,klev)
    ptte_out(jl,klev)=zcoef*(zdudz*zdudz+zdvdz*zdvdz)
  ENDDO

   status = 0 ! NO ERROR
   
END SUBROUTINE ! edith_vdiff_mol_fric


SUBROUTINE edith_iondrag(krow, kproma, kbdim, klev, pum1, pvm1, pqm1, pgeom1, pcp, ptte,ztendu,ztendv, &
                lat, time_step_len &
                ,kp_data) !ka_sb_20180419

!**** *iondrag* - does the ion drag that is exerted on the horizontal wind
!                 fields u and v above 100 km of altitude.
!
!     Subject.
!
!       This routine computes the physical tendencies of the horizontal
!   wind fields and temperature due to ion drag above 100 km of altitude.
!   The tendencies for winds are obtained from a semi-implicit time-stepping
!   procedure that is precise to order (dtime)^2.
!
!      f(t+dtime)-f(t-dtime)    - M * ( f(t+dtime) + f(t-dtime) )
!      --------------------- =        ---------------------------
!             2*dtime                              2
!
!
!   f is is the 2D wind vector (u,v) and M is a 2x2 matrix representing the
!   drag coefficients.
!   The values for the matrix M are calculated using the method of
!   Hong and Lindzen, 1976: JAS, 33, 135-153.  for minimum solar forcing.
!   The Lorenz terms (zldrag) are from m.charron and mimic those of Hong and 
!   Lindzen (fig. 4).
!
!           _                     _
!          |                       |
!          | zdrag        zldrag   |
!      M = |                       |
!          | -zldrag   zcoef*zdrag |
!          |_                     _|
!           _   _
!          |     |
!          |  u  |
!      f = |     |
!          |  v  |
!          |_   _|
!
!
!**   Interface.
!     ----------
!
!          *iondrag* is called from *physc*.
!
!
! INPUT ARGUMENTS:
! ---------------
!
!  pum1     : zonal wind (t-dt)
!  pvm1     : meridional wind (t-dt)
!  pqm1     : humidity (t-dt)
!  pgeom1   : geopotential above surface (t-dt)
!
!
! INPUT/OUTPUT ARGUMENTS:
! ----------------------
!
!  pvol     : tendency of meridional wind
!  pvom     : tendency of zonal wind
!  ptte     : tendency of temperature
!
!
! Modifications: 
! --------------
!
!    H. Schmidt - MPI - 20020702
!       - bug fix: msis variable index counts from bottom to top
!    H. Schmidt - MPI - 20030311
!       - nproma



!USE mo_gaussgrid,      ONLY: gl_twomu
!USE mo_constants,      ONLY: g, vtmpc2
! op_pj_20180713+: illegal use of other submodel
!!$USE messy_e5vdiff,         ONLY: cvdifts
USE messy_main_constants_mem, ONLY: cvdifts
! op_pj_20180713-
!USE mo_time_control,   ONLY: time_step_len
!USE messy_main_data_bi,         ONLY: time_step_len
!USE mo_geoloc,         ONLY: ilat
! ka_av_20150825+
!!! temporali as internal parameter
!USE mo_param_switches, ONLY: solvar
! ka_av_20150825-
!USE mo_kind,           ONLY: dp
USE messy_main_constants_mem, ONLY: DTR, g, vtmpc2

IMPLICIT NONE

!----- SUBROUTINE ARGUMENTS -----

INTEGER, INTENT(IN   ) :: krow, kproma, kbdim, klev
REAL(dp),INTENT(IN   ) :: pum1  (kbdim,klev),  pvm1(kbdim,klev),  pqm1(kbdim,klev)
REAL(dp),INTENT(IN   ) :: pgeom1(kbdim,klev)
REAL(dp),INTENT(IN   ) :: pcp   (kbdim,klev)
REAL(dp),INTENT(OUT) :: ztendu  (kbdim,klev),  ztendv(kbdim,klev)
REAL(dp), INTENT(OUT) :: ptte(kbdim,klev)
REAL(dp), INTENT(IN)   :: lat(kbdim)
REAL(dp), INTENT(IN)   :: time_step_len
!ka_sb_20180419+
REAL(dp),  INTENT(IN) :: kp_data
!ka_sb_20180419-
!----- INTERNAL VARIABLES -----
! ka_av_20150825+
!INTEGER, PARAMETER              :: solvar=1
INTEGER :: solvar
! ka_av_20150825-
INTEGER                         :: jglat, jl, jk, index, jkk
REAL(dp)                        :: zlat, zcons1, zcons2, zcons4(3), ztmp(3)
REAL(dp)                        :: zalpha, zcons5, ztmst
REAL(dp), DIMENSION(kbdim)      :: zcoef, zcoef1
REAL(dp), DIMENSION(kbdim,klev) :: zdrag, zldrag, zdenum
REAL(dp), DIMENSION(kbdim,klev) :: zup1, zvp1
REAL(dp), DIMENSION(3)          :: za,zb,zc,zd,ze
REAL(dp)                        :: zal

!----- TABLE FOR ION DRAG COEFFICIENTS TAKEN FROM -----
!_____   HONG AND LINDZEN, 1976: JAS, 33, p. 152  -----

REAL(dp), PARAMETER, DIMENSION(3) :: zamin = (/  6.6E4_dp  , 1.56E5_dp ,  3.0E5_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zamax = (/ 1.15E5_dp  , 2.75E5_dp , 1.05E6_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zbmin = (/    1.4_dp  ,    1.0_dp ,   0.35_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zbmax = (/    1.4_dp  ,    1.0_dp ,   0.20_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zcmin = (/ 150.E3_dp  , 225.E3_dp , 275.E3_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zcmax = (/ 150.E3_dp  , 240.E3_dp , 300.E3_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zdmin = (/    0.2_dp  ,    0.0_dp ,    0.1_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zdmax = (/    0.2_dp  ,    0.0_dp ,    0.1_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zemin = (/   1.E3_dp  ,  42.E3_dp ,   1.E3_dp  /)
REAL(dp), PARAMETER, DIMENSION(3) :: zemax = (/   1.E3_dp  ,  52.E3_dp ,   1.E3_dp  /)
REAL(dp), PARAMETER               :: zcons3 = 5.0E-10_dp
REAL(dp), PARAMETER               :: zalmin = 37.3_dp
REAL(dp), PARAMETER               :: zalmax = 36.7_dp
!ka_sb_20180419+
REAL(dp)               :: kp_data_temp, kpweight_min, kpweight_max, kpmax, kpmin
!ka_sb_20180419-

kpmax=7._dp
kpmin=1._dp

solvar=ed_solvar

!----- SET INITIAL VALUES -----
!ka_sb_20180507+
SELECT CASE (solvar)

!IF (solvar == 1) THEN                        ! solar minimum conditions
CASE (1)                                      ! solar minimum conditions
  za=zamin
  zb=zbmin
  zc=zcmin
  zd=zdmin
  ze=zemin
  zal=zalmin
!ELSEIF (solvar == 2 .OR. solvar == 3) THEN   ! solar maximum conditions
CASE (2, 3)   ! solar maximum conditions
  za=zamax
  zb=zbmax
  zc=zcmax
  zd=zdmax
  ze=zemax
  zal=zalmax
!ka_sb_20180419+
!ELSEIF (solvar == 4) THEN! variable solar conditions
CASE (4)                                     ! variable solar conditions
! ! IF (kp_data(1)<kpmin) THEN
!!  IF (kp_data<kpmin) THEN
  IF (kp_data .LT. kpmin) THEN
    kp_data_temp=kpmin
!  ELSEIF (kp_data(1)>kpmax) THEN
!  ELSEIF (kp_data>kpmax) THEN
  ELSEIF (kp_data .GT. kpmax) THEN
    kp_data_temp=kpmax
  ELSE
!    kp_data_temp=kp_data(1)
    kp_data_temp=kp_data
  ENDIF    
kpweight_max=(kp_data_temp-kpmin)/(kpmax-kpmin)
kpweight_min=1-kpweight_max
  za=kpweight_min*zamin+kpweight_max*zamax
  zb=kpweight_min*zbmin+kpweight_max*zbmax
  zc=kpweight_min*zcmin+kpweight_max*zcmax
  zd=kpweight_min*zdmin+kpweight_max*zdmax
  ze=kpweight_min*zemin+kpweight_max*zemax
  zal=kpweight_min*zalmin+kpweight_max*zalmax
!ka_sb_20180419-
!ELSE
CASE default
  STOP 'wrong solar index in iondrag'
!ENDIF
END SELECT !(case solvar)
   !ka_sb_20180507-

ztmst=time_step_len
zalpha=cvdifts
zcons5=ztmst*zalpha
zdrag=0._dp
zldrag=0._dp
DO jl=1,kproma
!ka_sv_20180125+
  zlat=lat(jl)*DTR
!ka_sv_20180125-
  zcons1=-2._dp*TAN(zlat)
  zcons2=ATAN(zcons1)
  zcoef1(jl)=SIN(zcons2)
  zcoef(jl)=zcoef1(jl)*zcoef1(jl)
ENDDO

!ka_av_20150825+
DO jk=1,klev
!ka_av_20150825-
  DO jl=1,kproma
    IF (pgeom1(jl,jk) >= 100.E3_dp*g) THEN
      ztmp(:)=(pgeom1(jl,jk)/g-zc(:)) / (zd(:)*pgeom1(jl,jk)/g+ze(:))
      zcons4(:)=1._dp-ztmp(:)-EXP(-ztmp(:))
      DO index=1,3
        zdrag(jl,jk)=zdrag(jl,jk) + za(index)*EXP(zb(index)*zcons4(index))
      ENDDO
    ENDIF
  ENDDO
ENDDO

zdrag=zdrag*zcons3

!ka_av_20150825+
DO jk=1,klev
!ka_av_20150825-
  DO jl=1,kproma
    IF (pgeom1(jl,jk) >= 100.E3_dp*g) THEN
      zldrag(jl,jk)=zcoef1(jl)*EXP(1.9E-17_dp*(3.5_dp*(pgeom1(jl,jk)- &
             55.E3_dp*g))**3._dp/(EXP(3.5E-6_dp*(pgeom1(jl,jk)- &
             55.E3_dp*g))-1._dp)-zal)
    ENDIF
  ENDDO
ENDDO

DO jk=1,klev
  DO jl=1,kproma
    zdenum(jl,jk)=(1._dp+zcons5*zcoef(jl)*zdrag(jl,jk))*(1._dp+zcons5*zdrag(jl,jk))+   &
         zcons5*zldrag(jl,jk)*zcons5*zldrag(jl,jk)
    ztendu(jl,jk)=(-zdrag(jl,jk)*(1._dp+zcons5*zcoef(jl)*zdrag(jl,jk))*pum1(jl,jk)     &
         -zldrag(jl,jk)*(pvm1(jl,jk)+zcons5*zldrag(jl,jk)*pum1(jl,jk)))/zdenum(jl,jk)
    ztendv(jl,jk)=(-zcoef(jl)*zdrag(jl,jk)*(1._dp+zcons5*zdrag(jl,jk))*pvm1(jl,jk)         &
         +zldrag(jl,jk)*(pum1(jl,jk)-zcons5*zldrag(jl,jk)*pvm1(jl,jk)))/zdenum(jl,jk)
  ENDDO
ENDDO


zup1=pum1+ztmst*ztendu
zvp1=pvm1+ztmst*ztendv

DO jk=1,klev
  jkk=klev-jk+1
  DO jl=1,kproma

                              
    ptte(jl,jk)=( 0.5_dp*pum1(jl,jk)*pum1(jl,jk) + &
                                0.5_dp*pvm1(jl,jk)*pvm1(jl,jk) - &
                                0.5_dp*zup1(jl,jk)*zup1(jl,jk) - &
                                0.5_dp*zvp1(jl,jk)*zvp1(jl,jk))  &
                              /(ztmst*(1._dp+vtmpc2*pqm1(jl,jk))*pcp(jl,jk))

  ENDDO
ENDDO

!stop
END SUBROUTINE edith_iondrag


SUBROUTINE edith_photsrc( scolo2, sill, tgit, nk, nc, j_src, nsrc, mode )

  USE messy_main_constants_mem, ONLY : h_Planck, c_light
  
  implicit none

  INTEGER, INTENT(IN) :: nk, nc, mode, nsrc
  REAL(dp), DIMENSION(nc,nk), INTENT(IN) :: scolo2, tgit
  LOGICAL, INTENT(in), DIMENSION(nc,nk) :: sill
  REAL(dp), DIMENSION(nc,nk, nsrc), INTENT(OUT) :: j_src
  
  
!>     local
!!  variable         | description
!!  -----------------|--------------------------------------------------
!!  xsro2(nt,nbands) | O2 Schumann-Runge continuum absorption coefficient
!!        ~          | m2 at nt temperatures in nbands (averaged)
!!  xco2             | CO2 absorption coefficient at nbands
!!  wlg(3,nbands)    | wavelength, lbound, ubound (in m)
!!  sflux(nbands)    | mean solar spectral flux as W/m2/m
!!  j_src            | photolysis rate as 1/molecule/s

  integer, parameter:: nbands=6
!  real(dp), parameter:: constants_c=299792458., constants_h=6.62606876d-34

  real(dp):: xsrco2(1:nbands,1:2), wgl(1:nbands,1:3)
  real(dp):: xco2(1:nbands)
  real(dp), dimension(1:nbands):: tra, tra0, sigma, pflux, sflux, sflux0
! op_pj_20180727+: 0:1 results in error, since element (2) is used below
!!$  real(dp):: tlayer(nk), tx(0:1)
  real(dp):: tlayer(nk), tx(2) ! 1:2
! op_pj_20180727-
  integer:: i, j, k, ic
  logical:: first = .true.
!!$  =========================================================================
!!$  ini
  tx = (/ 90., 295. /)
  tra   = 1.0d0       
!!$  ini SRC cross section, wgl and some solar flux  
  xsrco2(:,1) = &
   (/ 4.1100e-23, 1.3523e-21, 1.2002e-21, 5.8720e-22, 1.5104e-22, 2.0180e-23 /)
  xsrco2(:,2) = &
   (/ 4.3338e-23, 1.4226e-21, 1.2703e-21, 6.2548e-22, 1.6849e-22, 2.2425e-23 /) 

  wgl(:,1) = &
   (/ 1.3000e-07, 1.3904e-07, 1.4808e-07, 1.5712e-07, 1.6616e-07, 1.7520e-07 /)
  wgl(:,2) = &
   (/ 1.3000e-07, 1.3452e-07, 1.4356e-07, 1.5260e-07, 1.6164e-07, 1.7068e-07 /)
  wgl(:,3) = &
   (/ 1.3452e-07, 1.4356e-07, 1.5260e-07, 1.6164e-07, 1.7068e-07, 1.7520e-07 /)

  sflux0 = &
   (/ 8.4131e+04, 4.5845e+04, 7.9565e+04, 1.8492e+05, 4.1753e+05, 7.9280e+05 /)
   
  xco2(:) = &
   (/ 5.e-23, 7.e-23, 7.e-23, 4.e-23, 1.0e-24, 5.e-25 /)

!!$  ============================================================================

  select case ( mode )
  case ( 0 )
     sflux = sflux0
  case default
     print *, 'photsrc: do not know what solar spectrum to use.'
     stop
  end select

  
!!$ calculate photon flux in interval
  pflux = sflux/h_Planck/c_light*( wgl(:,3)**2._dp - wgl(:,2)**2._dp )/2._dp

  j_src = 0.d0

  do ic=1,nc
     tlayer = tgit(ic,:)
     do k=1,nk
        if ( .not. sill(ic,k) ) cycle
        sigma = xsrco2(:,1) + ( xsrco2(:,2) - xsrco2(:,1) ) / &
             ( tx(2) - tx(1) ) * ( tlayer(k) - tx(1) )
        
        tra = exp( -sigma*scolo2(ic,k) )
        j_src(ic,k,1) = sum( tra*pflux*sigma )
        j_src(ic,k,2) = sum( tra*pflux*xco2 )
     enddo
  enddo

  return

END SUBROUTINE edith_photsrc



!c=======================================================================

SUBROUTINE edith_rcolsph( cv, nk, sza, nc, solcol, sill, press )
      
!      This routine calculates slant columns of a solar beam
!      for a spherical symmetric (dependent only on radius)
!      atmosphere. The slant column only depends on solar zenith angle
!      depth (=height). Spherical symmetry is sufficient to be
!      fulfilled only locally, that is horizontal gradients are
!      negligible. 
!      If used with pressure height as vertical coordinate
!      the slant columns calculated are only a approximation.
!      molecular mass of air is constant
!      (otherwise code has to be modified to have air density
!      as input.
! 
!      Parameters:
!      Variable      | Description
!      -----------------------------------------------------
!      cv(nk,nc)     | nc mixing ratios with nk levels 
!      z(nk)         | vertical grid
!      solcol(nk,nc) | slant column
!      sill(nk,nc)   | logical field illuminated && sza < limit     
!
!      press: needed for calculation of pressure altitude;
!             we assume that it is constant on vertical grid levels in altitudes where we need it

!c======================================================================

USE messy_main_constants_mem, ONLY : radius_earth, pi, rho_air, M_air, atm2Pa
      implicit none

!c     --- parameter
      integer, intent(in) :: nk, nc
      real(dp), intent(in) :: cv(nc,nk), sza(nc)
      real(dp), intent(out) :: solcol(nc,nk)
      real(dp) :: z(nk)
      
      logical, intent(out) :: sill(nc,nk)


!c     --- local variables
      integer :: ic, k, l, isz

!      logical, save:: first = .true.
      integer, parameter:: nsza = 101
      real(dp), parameter:: szamax = 100.
      real(dp) :: szx, gdis1, gdis0, hmin
      real(dp):: sphweight(nk,nk,nsza), szgrid(nsza)
      real(dp) :: cvprof(nk+1), zz(nk+1), rr(nk+1)
      logical:: illgrid(nk,nsza)
      
      real(dp), parameter :: H = 7000._dp       ! scale_height
      real(dp), DIMENSION(nk), intent(in) :: press
      real(dp), DIMENSION(nk) :: rou

!c     --- initialize geometrical weight
!c         fixed sza grid
!c         fixed molecular weight

         sphweight = 0.
         illgrid = .false.
         do k=1,nk
            z(k) = -H*log(press(k)/atm2Pa)
            rou(k) = exp(-z(k)/H)
            rr(k) = rou(k)
            zz(k) = z(k)
         enddo
         zz(nk) = z(nk) + H
         rr(nk+1) = rou(nk) * exp( -1.)
         do isz=1,nsza
            szgrid(isz) = szamax/(nsza-1)*( isz - 1 )
            szx = szgrid(isz)/180.*pi
            do k=1,nk
               hmin = ( radius_earth + z(k) )*sin(pi-szx) - radius_earth
               if ( szx .le. pi/2 .or. hmin .gt. 0.0 ) then
                  gdis0 = sqrt( ( radius_earth + z(k) )**2 - ( radius_earth + hmin )**2 )! 0.0
                  illgrid(k,isz) = .true.
                  do l = 1, nk
                     if ( zz(l) .gt. hmin ) then
                        gdis1 = sqrt( ( radius_earth + zz(l) )**2 - ( radius_earth + hmin )**2 )
                        if ( ( zz(l) .gt. z(k) )  &
                            .or. ( szx .gt. pi/2 ) ) then
                           sphweight(l,k,isz)  = abs( gdis1 - gdis0 )    &
                                              *( rr(l) + rr(l+1) )/2    &
                                              * rho_air / M_air
                        endif
                        if ( z(l) .le. z(k) .and. szx .gt. pi/2 ) then
                           sphweight(l,k,isz) = sphweight(l,k,isz)*2
                        endif
                        gdis0 = gdis1 ! must be inside if-clause
                     endif
                  enddo
               endif

           enddo ! nk
         enddo ! nsza
 
!c     --- column density mol/m2
      solcol = 0.d0
      do ic=1,nc
         if ( sza(ic) .ge. szamax ) then
            sill(ic,:) = .false.
            solcol(ic,:) = 0.
            cycle
         endif
!c        --- select sza index
         isz = minloc( abs( sza(ic) - szgrid ), dim=1 )
         do k=1,nk
            sill(ic,k) = illgrid(k,isz)
            if ( sill(ic,k) ) then
               do l=1,nk
                  solcol(ic,k) = solcol(ic,k)         &
                              + cv(ic,l)*sphweight(l,k,isz)
               enddo
           endif
         enddo
      enddo
         
END SUBROUTINE edith_rcolsph


! **********************************************************************
END MODULE messy_edith
! **********************************************************************

