program stratO3bud

!#############################################################################
!# calculate and save monthly mean ozone production and destruction rates 
!# from EMAC output
!# developed by Patrick Joeckel
!# implemented in Fortran by Hella Garny
!# 
!# 12.5.2011: update to reaction rates of CAABA/MECCA chemistry
!# 11.5.2016: update to work alternatively with qm1 instead of H2O
!##############################################################################
!#
!# call with stratO3bud.exe
!#
!# input: ECHAM5.nc, tracer_gp.nc, jval.nc
!# output: stratO3bud.nc (monthly mean)
!##############################################################################

implicit none

real, parameter :: Na = 6.022045E23, Rgas = 8.314409, cmfac = 1.e6

!integer*4 :: iargc
character(len=2) ::  month
character(len=4) :: year

character(len=80) :: infile, outfile
integer :: ilo, ila, ilev, it
integer :: nt, nlev, nlat, nlon
real ::  k_NO_HO2, k_NO_CH3O2, k_O3_O3P , k_O1D_H2O &
     , k_NO2_O3P, k_O3_HO2, k_OH_O3, k_HO2_O3P , k_O3_H &
     , k_OH_O3P, k_ClO_O3P , k_ClO_HO2
real :: k_BrO_ClO, k_BrO_O3P, k_ClO_ClO, press, cair !& ! op_pj_20130722
!!$     , zt_help, k0_T, kinf_T, k_ratio, nu            ! op_pj_20130722
integer :: status ! op_pj_20160511

real, dimension(:), allocatable :: hyam, hybm, dataTime, dataLat, dataLon
real, dimension(:,:,:), allocatable :: aps
real, dimension(:,:,:,:), allocatable :: O3,  NO, NO2, O1D, O3P, HO2, OH, H, CH3O2, ClO, BrO, J_O2, H2O , temp  
real, dimension(:,:), allocatable :: Prod, ProdHO2, ProdCH3O2, LossOx, LossNOx, LossHOx, LossClOx, LossBrOx 

! op_pj_20110505+
!!$
!!$call getarg(1,year)
!!$write(*,*)'year: ',year
!!$
!!$call getarg(2,month)
!!$write(*,*)'month: ',month
!!$
!!$outfile = '/work/scratch/b/b309076/StratO3Bud/TS2000_O3Budget_'//year//month//'.nc'
outfile='stratO3bud.nc'

! read data fields:
! hyam, hybm, aps(=surface pressure; i.e. nlonxnlatxnt) 
! tm1 from ECHAM5
!!$infile = '/work/bd0617/b302019/TS2000/TS2000_10/messy/TS2000_10_'//year//month//'.01_messy.nc'
infile='tr_O3_sbud.nc'
! op_pj_20110505-
CALL getDimLen(infile,nlon,nlat,nlev,nt)

allocate (temp(nlon,nlat,nlev,nt))
allocate (aps(nlon,nlat,nt))
allocate (hyam(nlev))
allocate (hybm(nlev))

CALL readData(infile,'tm1',temp,nlon,nlat,nlev,nt)
CALL readData3D(infile,'aps',aps,nlon,nlat,nt)
CALL readData1D(infile,'hyam',hyam,nlev)
CALL readData1D(infile,'hybm',hybm,nlev)


! tracers from tracer_gp 
! op_pj_20110505+
!!$infile = '/work/bd0617/b302019/TS2000/TS2000_10/tracer_gp/TS2000_10_'//year//month//'.01_tracer_gp.nc'
!CALL getDimLen(infile,nlon,nlat,nlev,nt)
! op_pj_20110505-

allocate (O3(nlon,nlat,nlev,nt))
allocate (NO(nlon,nlat,nlev,nt))
allocate (NO2(nlon,nlat,nlev,nt))
allocate (O1D(nlon,nlat,nlev,nt))
allocate (O3P(nlon,nlat,nlev,nt))
allocate (HO2(nlon,nlat,nlev,nt))
allocate (OH(nlon,nlat,nlev,nt))
allocate (H(nlon,nlat,nlev,nt))
allocate (CH3O2(nlon,nlat,nlev,nt))
allocate (ClO(nlon,nlat,nlev,nt))
allocate (BrO(nlon,nlat,nlev,nt))
allocate (H2O(nlon,nlat,nlev,nt))

CALL readData(infile,'O3',O3,nlon,nlat,nlev,nt)
CALL readData(infile,'NO',NO,nlon,nlat,nlev,nt)
CALL readData(infile,'NO2',NO2,nlon,nlat,nlev,nt)
CALL readData(infile,'O1D',O1D,nlon,nlat,nlev,nt)
CALL readData(infile,'O3P',O3P,nlon,nlat,nlev,nt)
CALL readData(infile,'HO2',HO2,nlon,nlat,nlev,nt)
CALL readData(infile,'OH',OH,nlon,nlat,nlev,nt)
CALL readData(infile,'H',H,nlon,nlat,nlev,nt)
CALL readData(infile,'CH3O2',CH3O2,nlon,nlat,nlev,nt)
CALL readData(infile,'ClO',ClO,nlon,nlat,nlev,nt)
CALL readData(infile,'BrO',BrO,nlon,nlat,nlev,nt)
! op_pj_20160511+
!!$CALL readData(infile,'H2O',H2O,nlon,nlat,nlev,nt)
!
! NOTE: If H2O is not available, it can be added (e.g. via nco based
!       preprocessing) of q(m1) as:
!       !                 q              H2O
!       ! H2O = scvmr * -----  => q = -----------
!       !               1 - q         scvmr + H2O
!       with scvmr = M_air/M_H2O = 28.970 / 18.02
CALL check_variable(infile, 'H2O', status)
IF (status == 0) then
   CALL readData(infile,'H2O',H2O,nlon,nlat,nlev,nt)
else
   write(*,*) 'variable H2O not available, trying (and converting) qm1 ...'
   CALL readData(infile,'qm1',H2O,nlon,nlat,nlev,nt)
   H2O = (28.970 / 18.02) * H2O / (1.0 - H2O)
endif
! op_pj_20160511-

! photolysis from jval_gp
! op_pj_20110505+
!!$infile = '/work/bd0617/b302019/TS2000/TS2000_10/jval_gp/TS2000_10_'//year//month//'.01_jval_gp.nc'
!!$CALL getDimLen(infile,nlon,nlat,nlev,nt)
! op_pj_20110505-

allocate (J_O2(nlon,nlat,nlev,nt))
CALL readData(infile,'J_O2',J_O2,nlon,nlat,nlev,nt)


!tracers: O3, NO, NO2, O1D, O3P, HO2, OH, H, CH3O2, ClO, BrO, J_O2, H2O 
! are in tracer_gp in  vmr =[mol/mol]

! calculate prod and loss at each gridpoint!
allocate (Prod(nlat,nlev))
allocate (ProdHO2(nlat,nlev))
allocate (ProdCH3O2(nlat,nlev))
allocate (LossOx(nlat,nlev))
allocate (LossNOx(nlat,nlev))
allocate (LossHOx(nlat,nlev))
allocate (LossClOx(nlat,nlev))
allocate (LossBrOx(nlat,nlev))

! set to 0
Prod(:,:) = 0.0
ProdHO2(:,:) = 0.0
ProdCH3O2(:,:) = 0.0
LossOx(:,:) = 0.0
LossNOx(:,:) = 0.0
LossHOx(:,:) = 0.0
LossClOx(:,:) = 0.0
LossBrOx(:,:) = 0.0

write(*,*) 'Calculate P and D values for each longitude'

do ilo=1,nlon
write(*,*) 'progress: ', ilo
  do ila=1,nlat
    do ilev=1,nlev
      do it=1,nt

        press  = hyam(ilev) + hybm(ilev)*aps(ilo,ila,it)
        cair = Na/cmfac/Rgas*press/temp(ilo,ila,ilev,it)

! convert tracers from mixing ratio to numer density!
O3(ilo,ila,ilev,it) = O3(ilo,ila,ilev,it)*cair
NO(ilo,ila,ilev,it) = NO(ilo,ila,ilev,it)*cair
NO2(ilo,ila,ilev,it) = NO2(ilo,ila,ilev,it)*cair
O1D(ilo,ila,ilev,it) = O1D(ilo,ila,ilev,it)*cair
O3P(ilo,ila,ilev,it) = O3P(ilo,ila,ilev,it)*cair
HO2(ilo,ila,ilev,it) = HO2(ilo,ila,ilev,it)*cair
OH(ilo,ila,ilev,it) = OH(ilo,ila,ilev,it)*cair
H(ilo,ila,ilev,it) = H(ilo,ila,ilev,it)*cair
CH3O2(ilo,ila,ilev,it) = CH3O2(ilo,ila,ilev,it)*cair
ClO(ilo,ila,ilev,it) = ClO(ilo,ila,ilev,it)*cair
BrO(ilo,ila,ilev,it) = BrO(ilo,ila,ilev,it)*cair
H2O(ilo,ila,ilev,it) = H2O(ilo,ila,ilev,it)*cair

! op_pj_20130722+
include 'rate_coeff.inc'
!!$! set k values:
!!$k_NO_HO2   = 3.5E-12*EXP(250./temp(ilo,ila,ilev,it))  ! <G3201>
!!$k_NO_CH3O2 = 2.8E-12*EXP(300./temp(ilo,ila,ilev,it))  ! <G4104>
!!$k_O3_O3P   = 8.E-12*EXP(-2060./temp(ilo,ila,ilev,it)) ! <G1003>
!!$k_O1D_H2O  = 1.63E-10*EXP(60./temp(ilo,ila,ilev,it))  ! <G2111>
!!$k_NO2_O3P  = 5.1e-12*exp(210./temp(ilo,ila,ilev,it))  ! <G3105>
!!$
!!$k_O3_HO2   = 1e-14*exp(-490./temp(ilo,ila,ilev,it))   ! <G2107>
!!$k_OH_O3    = 1.7E-12*EXP(-940./temp(ilo,ila,ilev,it)) ! <G2104>
!!$k_HO2_O3P  = 3.E-11*EXP(200./temp(ilo,ila,ilev,it))   ! <G2106>
!!$k_O3_H     = 1.4E-10*EXP(-470./temp(ilo,ila,ilev,it)) ! <G2101>
!!$k_OH_O3P   = 2.2E-11*EXP(120./temp(ilo,ila,ilev,it))  ! <G2103>
!!$
!!$k_ClO_O3P  = 2.5e-11*exp(110./temp(ilo,ila,ilev,it))  ! <G6101>
!!$k_ClO_HO2  = 2.2E-12*EXP(340./temp(ilo,ila,ilev,it))  ! <G6204>
!!$
!!$k_BrO_ClO = 2.9e-12*exp(220./temp(ilo,ila,ilev,it)) + 5.8e-13*exp(170./temp(ilo,ila,ilev,it)) ! <G7603b> <G7603c> (<G7603a not used !?)
!!$k_BrO_O3P = 1.9E-11*EXP(230./temp(ilo,ila,ilev,it))   ! <G7101>
!!$
!!$! k_ClO_ClO    = k_3rd_iupac(temp,cair,2.E-32,4.,1.E-11,0.,0.45) ! <G6102d>
!!$zt_help = 300./temp(ilo,ila,ilev,it)
!!$k0_T    = 2.e-32  * zt_help**4 * cair ! k_0   at current T
!!$kinf_T  = 1.e-11 * zt_help**0       ! k_inf at current T
!!$k_ratio = k0_T/kinf_T
!!$nu      = 0.75-1.27*LOG(0.45)
!!$
!!$k_ClO_ClO = k0_T/(1.+k_ratio)* &
!!$                     0.45**(1./(1.+(LOG(k_ratio)/nu)**2))
! op_pj_20130722-

! now get prod and loss rates in [molec/cm3/s]
! convert vmr to nd by *cair!


Prod(ila,ilev) = Prod(ila,ilev) + 2 * J_O2(ilo,ila,ilev,it) * 0.21 *cair


ProdHO2(ila,ilev) = ProdHO2(ila,ilev) + ( k_NO_HO2 * NO(ilo,ila,ilev,it) * HO2(ilo,ila,ilev,it) ) 


ProdCH3O2(ila,ilev) = ProdCH3O2(ila,ilev) + ( k_NO_CH3O2 * NO(ilo,ila,ilev,it) * CH3O2(ilo,ila,ilev,it) ) 


LossOx(ila,ilev) = LossOx(ila,ilev) + ( 2 * k_O3_O3P * O3(ilo,ila,ilev,it) * O3P(ilo,ila,ilev,it) + &
                   k_O1D_H2O * O1D(ilo,ila,ilev,it) * H2O(ilo,ila,ilev,it) ) 

LossNOx(ila,ilev) = LossNOx(ila,ilev) + 2 * k_NO2_O3P * NO2(ilo,ila,ilev,it) * O3P(ilo,ila,ilev,it) 

LossHOx(ila,ilev) = LossHOx(ila,ilev) +  ( k_O3_HO2 * O3(ilo,ila,ilev,it) * HO2(ilo,ila,ilev,it) + &
                    k_OH_O3 * OH(ilo,ila,ilev,it) * O3(ilo,ila,ilev,it) + &
                    k_HO2_O3P * HO2(ilo,ila,ilev,it) * O3P(ilo,ila,ilev,it) + &
                    k_O3_H * O3(ilo,ila,ilev,it) * H(ilo,ila,ilev,it) + &
                    k_OH_O3P * OH(ilo,ila,ilev,it) * O3P(ilo,ila,ilev,it) ) 

LossClOx(ila,ilev) = LossClOx(ila,ilev) +  ( 2* k_ClO_O3P * ClO(ilo,ila,ilev,it) * O3P(ilo,ila,ilev,it) + &
                     k_ClO_HO2 * ClO(ilo,ila,ilev,it) * HO2(ilo,ila,ilev,it) + &
                     2 * k_ClO_ClO * ClO(ilo,ila,ilev,it) * CLO(ilo,ila,ilev,it) ) 

LossBrOx(ila,ilev) = LossBrOx(ila,ilev) + (2 * k_BrO_ClO * BrO(ilo,ila,ilev,it) * ClO(ilo,ila,ilev,it) + &
                    2 * k_BrO_O3P * BrO(ilo,ila,ilev,it) * O3P(ilo,ila,ilev,it) ) 

      end do
    end do
  end do
end do

Prod=Prod/nlon/nt
ProdHO2=ProdHO2/nlon/nt
ProdCH3O2=ProdCH3O2/nlon/nt
LossOx=LossOx/nlon/nt
LossNOx=LossNOx/nlon/nt
LossHOx=LossHOx/nlon/nt
LossClOx=LossClOx/nlon/nt
LossBrOx=LossBrOx/nlon/nt

write(*,*) 'calculation finished'


! save monthly mean Prod and Loss rates

CALL createFile2D(outfile,nlat,nlev)

! add dimension variables
allocate(dataLat(nlat))
CALL readData1D(infile, 'lat', dataLat ,nlat)
CALL addVar1D(outfile,'lat','lat',dataLat,nlat)

CALL addVar2D(outfile,'Prod',Prod,nlat,nlev)
CALL addVar2D(outfile,'ProdHO2',ProdHO2,nlat,nlev)
CALL addVar2D(outfile,'ProdCH3O2',ProdCH3O2,nlat,nlev)
CALL addVar2D(outfile,'LossOx',LossOx,nlat,nlev)
CALL addVar2D(outfile,'LossNOx',LossNOx,nlat,nlev)
CALL addVar2D(outfile,'LossHOx',LossHOx,nlat,nlev)
CALL addVar2D(outfile,'LossClOx',LossClOx,nlat,nlev)
CALL addVar2D(outfile,'LossBrOx',LossBrOx,nlat,nlev)

deallocate(Prod)
deallocate(ProdHO2)
deallocate(ProdCH3O2)
deallocate(LossOx)
deallocate(LossNOx)
deallocate(LossHOx)
deallocate(LossClOx)
deallocate(LossBrOx)

contains

! op_pj_20130722+
! see mecca.eqn ...

  ELEMENTAL REAL FUNCTION k_3rd_iupac(temp,cair,k0_300K,n,kinf_300K,m,fc)
    ! IUPAC three body reaction formula (www.iupac-kinetic.ch.cam.ac.uk)

    INTRINSIC :: LOG10

    REAL,     INTENT(IN) :: temp      ! temperature [K]
    REAL,     INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL,     INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL,     INTENT(IN) :: n         ! exponent for low pressure limit
    REAL,     INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL,     INTENT(IN) :: m         ! exponent for high pressure limit
    REAL,     INTENT(IN) :: fc        ! broadening factor (e.g. 0.45 or 0.6...)
    REAL                 :: nu        ! N
    REAL                 :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300./temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    nu      = 0.75-1.27*LOG10(fc)
    k_3rd_iupac = k0_T/(1.+k_ratio)* &
      fc**(1./(1.+(LOG10(k_ratio)/nu)**2))

  END FUNCTION k_3rd_iupac
! op_pj_20130722-

  ! mz_rs_20180312+
  ! from kpp/util/UserRateLaws.f90 but without dp
  ELEMENTAL REAL FUNCTION k_3rd(temp,cair,k0_300K,n,kinf_300K,m,fc)

    INTRINSIC LOG10

    REAL, INTENT(IN) :: temp      ! temperature [K]
    REAL, INTENT(IN) :: cair      ! air concentration [molecules/cm3]
    REAL, INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
    REAL, INTENT(IN) :: n         ! exponent for low pressure limit
    REAL, INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
    REAL, INTENT(IN) :: m         ! exponent for high pressure limit
    REAL, INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
    REAL             :: zt_help, k0_T, kinf_T, k_ratio

    zt_help = 300./temp
    k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
    k_ratio = k0_T/kinf_T
    k_3rd   = k0_T/(1.+k_ratio)*fc**(1./(1.+LOG10(k_ratio)**2))

  END FUNCTION k_3rd
  ! mz_rs_20180312-
  
end program stratO3bud
