MODULE MICROPHYSICS

! module for original necessary microphysical variables
! Harald Rybka, 2019

use grid, only: dp, nx,ny,nzm,nz, dimx1_s,dimx2_s,dimy1_s,dimy2_s, crmvars
use micro_params

implicit none

PRIVATE :: dp

!==================================================================================

INTEGER :: micro_scheme  ! selected microphysics scheme (0=SAM1MOM, 1=M2005)
LOGICAL :: domicro_sam1mom = .FALSE.
LOGICAL :: domicro_m2005   = .FALSE.
INTEGER :: nmicro_fields ! total number of prognostic micro. variables

! microphysics prognostic variables are stored in this array:
REAL(dp), POINTER, DIMENSION(:,:,:,:) :: micro_field => NULL()

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqcl, iqci, iqr, iqs, iqg, incl, inci, inr, ins, ing

! flags and indices
INTEGER, POINTER, DIMENSION(:)        :: flag_wmass         => NULL()
INTEGER, POINTER, DIMENSION(:)        :: flag_precip        => NULL()
INTEGER, POINTER, DIMENSION(:)        :: flag_number        => NULL()
INTEGER, POINTER, DIMENSION(:)        :: flag_micro3Dout    => NULL()
INTEGER                               :: index_water_vapor  
INTEGER                               :: index_cloud_ice    

REAL(dp), POINTER, DIMENSION(:,:,:)   :: fluxbmk  => NULL() ! surface flux of tracers
REAL(dp), POINTER, DIMENSION(:,:,:)   :: fluxtmk  => NULL() ! top boundary flux of tracers
REAL(dp), POINTER, DIMENSION(:,:)     :: mkwle    => NULL() ! resolved vertical flux
REAL(dp), POINTER, DIMENSION(:,:)     :: mkwsb    => NULL() ! SGS vertical flux
REAL(dp), POINTER, DIMENSION(:,:)     :: mkadv    => NULL() ! tendency due to vertical advection
REAL(dp), POINTER, DIMENSION(:,:)     :: mklsadv  => NULL() ! tendency due to large-scale vertical advection
REAL(dp), POINTER, DIMENSION(:,:)     :: mkdiff   => NULL() ! tendency due to vertical diffusion
REAL(dp), POINTER, DIMENSION(:,:)     :: mksed    => NULL() ! sedimentation vertical flux
REAL(dp), POINTER, DIMENSION(:,:)     :: mfrac    => NULL() ! fraction of domain with microphysical quantity > 1.e-6
REAL(dp), POINTER, DIMENSION(:,:)     :: stend    => NULL() ! tendency due to sedimentation
REAL(dp), POINTER, DIMENSION(:,:)     :: mtend    => NULL() ! tendency due to micro. processes (other than sediment.)
REAL(dp), POINTER, DIMENSION(:,:)     :: trtau    => NULL() ! optical depths of various species

CHARACTER(LEN=3) , POINTER, DIMENSION(:)  :: mkname        => NULL()
CHARACTER(LEN=80), POINTER, DIMENSION(:)  :: mklongname    => NULL()
CHARACTER(LEN=10), POINTER, DIMENSION(:)  :: mkunits       => NULL()
REAL(dp)         , POINTER, DIMENSION(:)  :: mkoutputscale => NULL()

REAL(dp), POINTER, DIMENSION(:,:,:) :: qn       => NULL() ! cloud condensate (liquid + ice)
REAL(dp), POINTER, DIMENSION(:)     :: qpsrc    => NULL() ! source of precipitation microphysical processes
REAL(dp), POINTER, DIMENSION(:)     :: qpevp    => NULL() ! sink of precipitating water due to evaporation

REAL(dp), POINTER, DIMENSION(:)     :: lfac     => NULL()
REAL(dp), POINTER, DIMENSION(:,:,:) :: tauliq   => NULL()
REAL(dp), POINTER, DIMENSION(:,:,:) :: tauice   => NULL()
REAL(dp), POINTER, DIMENSION(:)     :: tmtend   => NULL()
REAL(dp), POINTER, DIMENSION(:,:,:) :: tmtend3d => NULL()

REAL(dp) :: vrain, vsnow, vgrau, crain, csnow, cgrau  ! precomputed coefs for precip terminal velocity
REAL(dp) :: sfcpcp, sfcicepcp

CONTAINS

!==================================================================================

! Supply function that computes total water in a domain:
real(dp) function total_water()

  use vars, only : nstep,nprint,adz,dz,rho
  real(dp) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water

!==================================================================================

subroutine micro_print()
  implicit none
  integer :: k

  ! print out min/max values of all microphysical variables
  do k=1,nmicro_fields
     call fminmax_print(trim(mkname(k))//':', &
          micro_field(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
  end do

end subroutine micro_print

!==================================================================================

subroutine micro_setparm()
  
  use vars, only: tlatqi
  
  use module_mp_GRAUPEL, only: dopredictNc, &
       Nc0,                &  ! initial/specified cloud droplet number conc (#/cm3)
       ccnconst, ccnexpnt, &  ! parameters for dospecifyaerosol=.false. (powerlaw CCN)
       aer_rm1, aer_rm2, &    ! two modes of aerosol for dospecifyaer...=.true.
       aer_n1, aer_n2, &      ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
       aer_sig1, aer_sig2, &  ! sig=geom standard deviation of aerosol size distn.
       dofix_pgam, pgam_fixed ! option to specify pgam (exponent of cloud water's gamma distn)

  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder

  IF (domicro_sam1mom) THEN
     ! no user-definable options in SAM1MOM microphysics.  
  ELSEIF (domicro_m2005) THEN
     Nc0 = 100.          ! default droplet number concentration
     ccnconst = 120.     ! maritime value (/cm3), adapted from Rasmussen 
     ccnexpnt = 0.4      !   et al (2002) by Hugh Morrison et al.  Values
                         !   of 1000. and 0.5 suggested for continental
     aer_rm1 = 0.052     ! two aerosol mode defaults from MPACE (from Hugh)
     aer_sig1 = 2.04
     aer_n1 = 72.2
     aer_rm2 = 1.3
     aer_sig2 = 2.5
     aer_n2 = 1.8
   
     dofix_pgam = .false.
     pgam_fixed = 5.     ! middle range value -- corresponds to radius dispersion ~ 0.4

     ! scale values of parameters for m2005micro
     aer_rm1 = 1.e-6*aer_rm1 ! convert from um to m
     aer_rm2 = 1.e-6*aer_rm2 
     aer_n1 = 1.e6*aer_n1    ! convert from #/cm3 to #/m3
     aer_n2 = 1.e6*aer_n2

     ! specify index of various quantities in micro_field array
     !  *** note that not all of these may be used if(.not.doicemicro) ***
     iqv  = 1  ! water vapor mass mixing ratio [kg H2O / kg dry air]
     iqcl = 2  ! cloud water mass mixing ratio [kg H2O / kg dry air]
     
     if(dopredictNc) then
        incl = 3  ! cloud water number mixing ratio [#/kg dry air]
        iqr = 4   ! rain mass mixing ratio [kg H2O / kg dry air]
        inr = 5   ! rain number mixing ratio [#/kg dry air]
        iqci = 6  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
        inci = 7  ! cloud ice number mixing ratio [#/kg dry air]
        iqs = 8   ! snow mass mixing ratio [kg H2O / kg dry air]
        ins = 9   ! snow number mixing ratio [#/kg dry air]
        iqg = 10  ! graupel mass mixing ratio [kg H2O / kg dry air]
        ing = 11  ! graupel number mixing ratio [#/kg dry air]
     else
        iqr = 3   ! rain mass mixing ratio [kg H2O / kg dry air]
        inr = 4   ! rain number mixing ratio [#/kg dry air]
        iqci = 5  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
        inci = 6  ! cloud ice number mixing ratio [#/kg dry air]
        iqs = 7   ! snow mass mixing ratio [kg H2O / kg dry air]
        ins = 8   ! snow number mixing ratio [#/kg dry air]
        iqg = 9   ! graupel mass mixing ratio [kg H2O / kg dry air]
        ing = 10  ! graupel number mixing ratio [#/kg dry air]
     end if

     index_water_vapor = iqv ! set SAM water vapor flag
     index_cloud_ice   = -1  ! historical variable (don't change)

     ! zero out statistics variables associated with cloud ice sedimentation
     ! in Marat's default SAM microphysics
     tlatqi = 0.

  END IF

end subroutine micro_setparm

!==================================================================================

END MODULE MICROPHYSICS

!==================================================================================

SUBROUTINE set_crm_micro_SAM1MOM(pnmicro_fields, pqcw0, pqci0,      &
                                 palphaelq, pbetaelq, pqp_threshold )

use grid,         only: dp
use micro_params, only: qcw0, qci0, alphaelq, betaelq, qp_threshold
use microphysics, only: domicro_sam1mom, domicro_m2005, nmicro_fields

IMPLICIT NONE

INTEGER :: pnmicro_fields

REAL(dp) :: pqcw0
REAL(dp) :: pqci0
REAL(dp) :: palphaelq
REAL(dp) :: pbetaelq
REAL(dp) :: pqp_threshold
 
nmicro_fields = pnmicro_fields

! autoconversion
qcw0 = pqcw0
qci0 = pqci0
alphaelq = palphaelq
betaelq  = pbetaelq
qp_threshold = pqp_threshold

domicro_sam1mom = .TRUE.
domicro_m2005   = .FALSE.

END SUBROUTINE set_crm_micro_SAM1MOM

!==================================================================================

SUBROUTINE set_crm_micro_M2005(pnmicro_fields, pdoicemicro, pdograupel, pdohail,  &
                               pdopredictNc, pdospecifyaerosol, pdosubgridw,      &
                               pdosb_warm_rain, pdoarcticicenucl, pdocloudedgeact )

use microphysics,      only: domicro_sam1mom, domicro_m2005, nmicro_fields
use module_mp_GRAUPEL, only: doicemicro, dograupel, dohail, dopredictNc,          &
                             dospecifyaerosol, dosubgridw, dosb_warm_rain,        &
                             doarcticicenucl, docloudedgeactivation

IMPLICIT NONE

INTEGER :: pnmicro_fields

LOGICAL :: pdoicemicro, pdograupel, pdohail, pdopredictNc, pdospecifyaerosol
LOGICAL :: pdosubgridw, pdosb_warm_rain, pdoarcticicenucl, pdocloudedgeact

nmicro_fields    = pnmicro_fields

doicemicro       = pdoicemicro
dograupel        = pdograupel
dohail           = pdohail
dopredictNc      = pdopredictNc
dospecifyaerosol = pdospecifyaerosol
dosubgridw       = pdosubgridw
dosb_warm_rain   = pdosb_warm_rain 
doarcticicenucl  = pdoarcticicenucl
docloudedgeactivation = pdocloudedgeact

domicro_sam1mom = .FALSE.
domicro_m2005   = .TRUE.

END SUBROUTINE set_crm_micro_M2005

!==================================================================================

SUBROUTINE CRM_ALLOCATE_MICROPHYSICS

  USE microphysics

  IMPLICIT NONE

  ALLOCATE(micro_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm, nmicro_fields)); micro_field=0.

  ALLOCATE(flag_wmass(nmicro_fields));      flag_wmass = 0 
  ALLOCATE(flag_precip(nmicro_fields));     flag_precip = 0
  ALLOCATE(flag_number(nmicro_fields));     flag_number = 0
  ALLOCATE(flag_micro3Dout(nmicro_fields)); flag_micro3Dout = 0

  ALLOCATE(fluxbmk(nx,ny,1:nmicro_fields)); fluxbmk = 0.
  ALLOCATE(fluxtmk(nx,ny,1:nmicro_fields)); fluxtmk = 0.
    
  ALLOCATE(mkwle(nz,1:nmicro_fields)); mkwle = 0.
  ALLOCATE(mkwsb(nz,1:nmicro_fields)); mkwsb = 0.
  ALLOCATE(mkadv(nz,1:nmicro_fields)); mkadv = 0.
  ALLOCATE(mklsadv(nz,1:nmicro_fields)); mklsadv = 0.
  ALLOCATE(mkdiff(nz,1:nmicro_fields)); mkdiff = 0.

  ALLOCATE(mkname(nmicro_fields))
  ALLOCATE(mklongname(nmicro_fields))
  ALLOCATE(mkunits(nmicro_fields))
  ALLOCATE(mkoutputscale(nmicro_fields)); mkoutputscale = 0.
    
  ALLOCATE(qn(nx,ny,nzm)); qn = 0.
  ALLOCATE(qpsrc(nz)); qpsrc = 0.
  ALLOCATE(qpevp(nz)); qpevp = 0.
    
  IF (domicro_sam1mom) THEN       
     ALLOCATE(accrsc(nzm)); accrsc = 0.
     ALLOCATE(accrsi(nzm)); accrsi = 0.
     ALLOCATE(accrrc(nzm)); accrrc = 0.
     ALLOCATE(coefice(nzm)); coefice= 0.
     ALLOCATE(accrgc(nzm)); accrgc = 0.
     ALLOCATE(accrgi(nzm)); accrgi = 0.
     ALLOCATE(evaps1(nzm)); evaps1 = 0.
     ALLOCATE(evaps2(nzm)); evaps2 = 0.
     ALLOCATE(evapr1(nzm)); evapr1 = 0.
     ALLOCATE(evapr2(nzm)); evapr2 = 0.
     ALLOCATE(evapg1(nzm)); evapg1 = 0.
     ALLOCATE(evapg2(nzm)); evapg2 = 0.
  ELSEIF(domicro_m2005) THEN
     ALLOCATE(lfac(1:nmicro_fields)); lfac = 0.
     ALLOCATE(tauliq(nx,ny,nzm)); tauliq = 0.
     ALLOCATE(tauice(nx,ny,nzm)); tauice = 0.

     ALLOCATE(mksed(nzm,1:nmicro_fields)); mksed = 0.
     ALLOCATE(mfrac(nzm,1:nmicro_fields)); mfrac = 0.
     ALLOCATE(stend(nzm,1:nmicro_fields)); stend = 0.
     ALLOCATE(mtend(nzm,1:nmicro_fields)); mtend = 0.
     ALLOCATE(trtau(nzm,1:nmicro_fields)); trtau = 0.
     
     ALLOCATE(tmtend(nzm)); tmtend = 0.
     ALLOCATE(tmtend3d(nx,ny,nzm)); tmtend3d = 0.
  END IF
  
END SUBROUTINE CRM_ALLOCATE_MICROPHYSICS

!==================================================================================

SUBROUTINE CRM_DEALLOCATE_MICROPHYSICS

  USE microphysics
  
  IMPLICIT NONE
    
  DEALLOCATE(micro_field)

  DEALLOCATE(flag_wmass)
  DEALLOCATE(flag_precip)
  DEALLOCATE(flag_number)
  DEALLOCATE(flag_micro3Dout)

  DEALLOCATE(fluxbmk)
  DEALLOCATE(fluxtmk)

  DEALLOCATE(mkwle)
  DEALLOCATE(mkwsb)
  DEALLOCATE(mkadv)
  DEALLOCATE(mklsadv)
  DEALLOCATE(mkdiff)
  
  DEALLOCATE(mkname)
  DEALLOCATE(mklongname)
  DEALLOCATE(mkunits)
  DEALLOCATE(mkoutputscale)
  
  DEALLOCATE(qn)
  DEALLOCATE(qpsrc)
  DEALLOCATE(qpevp)
  
  IF (domicro_sam1mom) THEN
     DEALLOCATE(accrsc)
     DEALLOCATE(accrsi)
     DEALLOCATE(accrrc)
     DEALLOCATE(coefice)
     DEALLOCATE(accrgc)
     DEALLOCATE(accrgi)
     DEALLOCATE(evaps1)
     DEALLOCATE(evaps2)
     DEALLOCATE(evapr1)
     DEALLOCATE(evapr2)
     DEALLOCATE(evapg1)
     DEALLOCATE(evapg2)
  ELSEIF (domicro_m2005) THEN
     DEALLOCATE(lfac)
     DEALLOCATE(tauliq)
     DEALLOCATE(tauice)

     DEALLOCATE(mksed)
     DEALLOCATE(mfrac)
     DEALLOCATE(stend)
     DEALLOCATE(mtend)
     DEALLOCATE(trtau)
     
     DEALLOCATE(tmtend)
     DEALLOCATE(tmtend3d)
  END IF
  
END SUBROUTINE CRM_DEALLOCATE_MICROPHYSICS
