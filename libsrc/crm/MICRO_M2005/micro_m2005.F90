module micro_m2005

! main interface to Morrison microphysics.
! original implementation by Peter Blossey, UW

#define CRM
use microphysics
use params, only: lcond, lsub, fac_cond, fac_sub, ggr, pres0

use grid,   only: dz, adz, dostatis, masterproc, &
                  docloud, doprecip,             &
                  doSAMconditionals, dosatupdnconditionals

use vars,   only: pres, rho, dtn, w, tke, t, tlatqi, condavg_mask, &
                  ncondavg, condavgname, condavglongname

use module_mp_GRAUPEL, only: GRAUPEL_INIT, M2005MICRO_GRAUPEL, &
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! use graupel
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization
      dopredictNc, &        ! prediction of cloud droplet number
      dospecifyaerosol, &   ! specify two modes of (sulfate) aerosol
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
      docloudedgeactivation,&! activate droplets at cloud edges as well as base
      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
      ccnconst, ccnexpnt, & ! parameters for dospecifyaerosol=.false. (powerlaw CCN)
      aer_rm1, aer_rm2, &   ! two modes of aerosol for dospecifyaer...=.true.
      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
      dofix_pgam, pgam_fixed ! option to specify pgam (exponent of cloud water's gamma distn)

#ifdef CRM
  use buffer, only: crmvars
#endif

implicit none

!logical :: isallocatedMICRO = .false.

!integer :: nmicro_fields ! total number of prognostic water vars

!#ifdef CRM
!  real :: micro_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,crmvars)  ! holds mphys quantities
!#else
!  real, allocatable, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities
!#endif

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
!integer :: iqv, iqcl, iqci, iqr, iqs, iqg, incl, inci, inr, ins, ing
!integer :: index_water_vapor ! separate water vapor index used by SAM

!real, allocatable, dimension(:) :: lfac
!integer, allocatable, dimension(:) :: flag_wmass, flag_precip, flag_number
!integer, allocatable, dimension(:) :: flag_micro3Dout

!integer, parameter :: index_cloud_ice = -1 ! historical variable (don't change)

!real, allocatable, dimension(:,:,:) :: fluxbmk, fluxtmk !surface/top fluxes
!real, allocatable, dimension(:,:,:) :: tauliq, tauice

!real, allocatable, dimension(:,:) :: & ! statistical arrays
!     mkwle, & ! resolved vertical flux
!     mkwsb, & ! SGS vertical flux
!     mksed, & ! sedimentation vertical flux
!     mkadv, & ! tendency due to vertical advection
!     mkdiff, &! tendency due to vertical diffusion
!     mklsadv, & ! tendency due to large-scale vertical advection
!     mfrac, & ! fraction of domain with microphysical quantity > 1.e-6
!     stend, & ! tendency due to sedimentation
!     mtend, & ! tendency due to microphysical processes (other than sedimentation)
!     trtau    ! optical depths of various species

!real, allocatable, dimension(:) :: tmtend

!real :: sfcpcp, sfcicepcp

! arrays with names/units for microphysical outputs in statistics.
!character*3, allocatable, dimension(:) :: mkname
!character*80, allocatable, dimension(:) :: mklongname
!character*10, allocatable, dimension(:) :: mkunits
!real, allocatable, dimension(:) :: mkoutputscale

! You can also have some additional, diagnostic, arrays, for example, total
! nonprecipitating cloud water, etc:

!bloss: array which holds temperature tendency due to microphysics
!real, allocatable, dimension(:,:,:), SAVE :: tmtend3d

CONTAINS

! um_hr_20190315 -> moved micro_setparm to microphysics.F90

!----------------------------------------------------------------------
!!! Read microphysical options from prm file and allocate variables
!
!subroutine micro_setparm()
!  use vars
!  implicit none

!  integer ierr, ios, ios_missing_namelist, place_holder
  
!   NAMELIST /MICRO_M2005/ &
!      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
!      dograupel, &          ! use graupel
!      dohail, &             ! graupel species has qualities of hail
!      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization in place of KK(2000)
!      dopredictNc, &        ! prediction of cloud droplet number
!      dospecifyaerosol, &   ! specify two modes of (sulfate) aerosol
!      dosubgridw, &         ! input estimate of subgrid w to microphysics
!      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
!      docloudedgeactivation,&! activate droplets at cloud edges as well as base
!      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
!      ccnconst, ccnexpnt, & ! parameters for dospecifyaerosol=.false. (powerlaw CCN)
!      aer_rm1, aer_rm2, &   ! two modes of aerosol for dospecifyaer...=.true.
!      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
!      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
!      dofix_pgam, pgam_fixed ! option to specify pgam (exponent of cloud water's gamma distn)

   !bloss: Create dummy namelist, so that we can figure out error code
   !       for a mising namelist.  This lets us differentiate between
   !       missing namelists and those with an error within the namelist.
!   NAMELIST /BNCUIODSBJCB/ place_holder

   ! define default values for namelist variables
!   doicemicro = .true.        ! use ice
!   dograupel = .true.         ! use graupel
!   dohail = .false.           ! graupel species has properties of graupel
!   dosb_warm_rain = .false.   ! use KK (2000) warm rain scheme by default
!   dopredictNc = .true.       ! prognostic cloud droplet number
!   dospecifyaerosol = .false. ! use powerlaw CCN relationship
!   dosubgridw=.false.         ! don't bother with estimating w_sgs for now
!   doarcticicenucl = .false.  ! use mid-latitude parameters
!   docloudedgeactivation  = .false. ! activate droplets at cloud base, not edges

!   Nc0 = 100. ! default droplet number concentration
   
!   ccnconst = 120.            ! maritime value (/cm3), adapted from Rasmussen 
!   ccnexpnt = 0.4             !   et al (2002) by Hugh Morrison et al.  Values
                              !   of 1000. and 0.5 suggested for continental
!   aer_rm1 = 0.052           ! two aerosol mode defaults from MPACE (from Hugh)
!   aer_sig1 = 2.04
!   aer_n1 = 72.2
!   aer_rm2 = 1.3
!   aer_sig2 = 2.5
!   aer_n2 = 1.8
   
!   dofix_pgam = .false.
!   pgam_fixed = 5. ! middle range value -- corresponds to radius dispersion ~ 0.4

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
!  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
  
  !bloss: get error code for missing namelist (by giving the name for
  !       a namelist that doesn't exist in the prm file).
!  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
!  rewind(55) !note that one must rewind before searching for new namelists

  !bloss: read in MICRO_M2005 namelist
!  read (55,MICRO_M2005,IOSTAT=ios)

!  if (ios.ne.0) then
!     !namelist error checking
!     if(ios.ne.ios_missing_namelist) then
!        write(*,*) '****** ERROR: bad specification in MICRO_M2005 namelist'
!        call task_abort()
!     elseif(masterproc) then
!        write(*,*) '****************************************************'
!        write(*,*) '****** No MICRO_M2005 namelist in prm file *********'
!        write(*,*) '****************************************************'
!     end if
!  end if
!  close(55)

!   if(dograupel.and..NOT.doicemicro) then
!      if(masterproc) write(*,*) 'doicemicro must be .true. for dograupel to be used.'
!      call task_abort()
!   end if

!   if(dohail.and..NOT.dograupel) then
!      if(masterproc) write(*,*) 'dograupel must be .true. for dohail to be used.'
!      call task_abort()
!   end if

   ! write namelist values out to file for documentation
!   if(masterproc) then
!      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.options_namelist', form='formatted', position='append')    
!      write (unit=55,nml=MICRO_M2005,IOSTAT=ios)
!      write(55,*) ' '
!      close(unit=55)
!   end if

   ! scale values of parameters for m2005micro
!   aer_rm1 = 1.e-6*aer_rm1 ! convert from um to m
!   aer_rm2 = 1.e-6*aer_rm2 
!   aer_n1 = 1.e6*aer_n1 ! convert from #/cm3 to #/m3
!   aer_n2 = 1.e6*aer_n2
   
!  nmicro_fields = 1 ! start with water vapor and cloud water mass mixing ratio
!  if(docloud) then
!     nmicro_fields = nmicro_fields + 1 ! add cloud water mixing ratio
!     if(dopredictNc) nmicro_fields = nmicro_fields + 1 ! add cloud water number concentration (if desired)
!  end if
!  if(doprecip)    nmicro_fields = nmicro_fields + 2 ! add rain mass and number (if desired)
!  if(doicemicro)  nmicro_fields = nmicro_fields + 4 ! add snow and cloud ice number and mass (if desired)
!  if(dograupel)   nmicro_fields = nmicro_fields + 2 ! add graupel mass and number (if desired)

  ! specify index of various quantities in micro_field array
  !  *** note that not all of these may be used if(.not.doicemicro) ***
!  iqv = 1   ! water vapor mass mixing ratio [kg H2O / kg dry air]
!  iqcl = 2  ! cloud water mass mixing ratio [kg H2O / kg dry air]
  
!  if(dopredictNc) then
!     incl = 3  ! cloud water number mixing ratio [#/kg dry air]
!     iqr = 4   ! rain mass mixing ratio [kg H2O / kg dry air]
!     inr = 5   ! rain number mixing ratio [#/kg dry air]
!     iqci = 6  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
!     inci = 7  ! cloud ice number mixing ratio [#/kg dry air]
!     iqs = 8   ! snow mass mixing ratio [kg H2O / kg dry air]
!     ins = 9   ! snow number mixing ratio [#/kg dry air]
!     iqg = 10  ! graupel mass mixing ratio [kg H2O / kg dry air]
!     ing = 11  ! graupel number mixing ratio [#/kg dry air]
!  else
!     iqr = 3   ! rain mass mixing ratio [kg H2O / kg dry air]
!     inr = 4   ! rain number mixing ratio [#/kg dry air]
!     iqci = 5  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
!     inci = 6  ! cloud ice number mixing ratio [#/kg dry air]
!     iqs = 7   ! snow mass mixing ratio [kg H2O / kg dry air]
!     ins = 8   ! snow number mixing ratio [#/kg dry air]
!     iqg = 9   ! graupel mass mixing ratio [kg H2O / kg dry air]
!     ing = 10  ! graupel number mixing ratio [#/kg dry air]
!  end if

  ! stop if icemicro is specified without precip -- we don't support this right now.
!  if((doicemicro).and.(.not.doprecip)) then
!     if(masterproc) write(*,*) 'Morrison 2005 Microphysics does not support both doice and .not.doprecip'
!     call task_abort()
!  end if
!  index_water_vapor = iqv ! set SAM water vapor flag

!  if(.not.isallocatedMICRO) then
!     ! allocate microphysical variables
!     allocate(&
!#ifdef CRM
!       micro_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,nmicro_fields), &
!#endif
!          fluxbmk(nx,ny,nmicro_fields), fluxtmk(nx,ny,nmicro_fields), &
!          tauliq(nx,ny,nzm), tauice(nx,ny,nzm), &
!          mkwle(nz,nmicro_fields), mkwsb(nz,nmicro_fields), &
!          mkadv(nz,nmicro_fields), mkdiff(nz,nmicro_fields), &
!          mklsadv(nz,nmicro_fields), &
!          stend(nzm,nmicro_fields), mtend(nzm,nmicro_fields), &
!          mfrac(nzm,nmicro_fields), trtau(nzm,nmicro_fields), &
!          mksed(nzm,nmicro_fields), tmtend(nzm), &
!          tmtend3d(nx,ny,nzm), flag_micro3Dout(nmicro_fields), &
!          flag_wmass(nmicro_fields), flag_precip(nmicro_fields), &
!          flag_number(nmicro_fields), lfac(nmicro_fields), &
!          mkname(nmicro_fields), mklongname(nmicro_fields), &
!          mkunits(nmicro_fields), mkoutputscale(nmicro_fields), STAT=ierr)
!     if(ierr.ne.0) then
!        write(*,*) 'Failed to allocate microphysical arrays on proc ', rank
!        call task_abort()
!     else
!        isallocatedMICRO = .true.
!     end if
!  end if

  ! zero out statistics variables associated with cloud ice sedimentation
  !   in Marat's default SAM microphysics
!  tlatqi = 0.

  ! initialize these arrays
!  micro_field = 0.
!  fluxbmk = 0.
!  fluxtmk = 0.
!  mkwle = 0.
!  mkwsb = 0.
!  mkadv = 0.
!  mkdiff = 0.
!  mklsadv = 0.

  ! initialize flag arrays to all mass, no number, no precip
!  flag_wmass = 1
!  flag_number = 0
!  flag_precip = 0
!  flag_micro3Dout = 0

!end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:
!
! this one is guaranteed to be called by SAM at the 
!   beginning of each run, initial or restart:
subroutine micro_init_2005()

  use vars
  
  implicit none
  
  real, dimension(nzm) :: qc0, qi0

  real, external :: satadj_water
  integer :: k

  ! initialize flag arrays
  if(dopredictNc) then
     ! Cloud droplet number concentration is a prognostic variable
     if(doicemicro) then
        if(dograupel) then
           flag_wmass  = (/1,1,0,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,0,1,1,0,0,1,1,1,1/)
           flag_number = (/0,0,1,0,1,0,1,0,1,0,1/)
        else
           flag_wmass  = (/1,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,0,1,1,0,0,1,1/)
           flag_number = (/0,0,1,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
           flag_wmass  = (/1,1,0,1,0/)
           flag_precip = (/0,0,0,1,1/)
           flag_number = (/0,0,1,0,1/)
        else
           flag_wmass  = (/1,1,0/)
           flag_precip = (/0,0,0/)
           flag_number = (/0,0,1/)
        end if
     end if
  else
     ! Cloud droplet number concentration is NOT a prognostic variable
     if(doicemicro) then
        if(dograupel) then
           flag_wmass  = (/1,1,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,1,1,0,0,1,1,1,1/)
           flag_number = (/0,0,0,1,0,1,0,1,0,1/)
        else
           flag_wmass  = (/1,1,1,0,1,0,1,0/)
           flag_precip = (/0,0,1,1,0,0,1,1/)
           flag_number = (/0,0,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
           flag_wmass  = (/1,1,1,0/)
           flag_precip = (/0,0,1,1/)
           flag_number = (/0,0,0,1/)
        elseif(docloud) then
           flag_wmass  = (/1,1/)
           flag_precip = (/0,0/)
           flag_number = (/0,0/)
        else
           flag_wmass  = (/1/)
           flag_precip = (/0/)
           flag_number = (/0/)
        end if
     end if
  end if

  ! output all microphysical fields to 3D output files if using more than
  !   just docloud.  Otherwise, rely on basic SAM outputs
  if(docloud.AND.(doprecip.OR.dopredictNc)) then
     flag_micro3Dout = 1
  end if

  ! initialize factor for latent heat
  lfac(:) = 1. ! use one as default for number species
  lfac(iqv) = lcond
  if(docloud) lfac(iqcl) = lcond
  if(doprecip) lfac(iqr) = lcond
  if(doicemicro) then
     lfac(iqci) = lsub
     lfac(iqs) = lsub
     if(dograupel) lfac(iqg) = lsub
  end if

  call graupel_init() ! call initialization routine within mphys module

  if(nrestart.eq.0) then

     ! condense out any water in excess of saturation
     do k = 1,nzm
        qc0(k) = 0. ! initialize cloud water and ice to zero
        qi0(k) = 0. 

        ! call saturation adjustment routine
        tabs0(k) = satadj_water(tabs0(k),q0(k),pres(k))

        ! find cloud water mixing ratio (if any)
        qc0(k) = max(0.,q0(k)-qsatw_crm(tabs0(k),pres(k)))
        qv0(k) = q0(k) - qc0(k)

        ! note that the microphysics will take care of the conversion of 
        !   the initial cloud water to ice at the appropriate temperatures.
        !   This will introduce some initial transients into the simulation, 
        !   but they should be small.
     end do

     ! initialize microphysical quantities
     do k = 1,nzm
        micro_field(:,:,k,iqv) = qv0(k)
        if(qc0(k).gt.0.) then
           micro_field(:,:,k,iqcl) = qc0(k)
           if(dopredictNc) micro_field(:,:,k,incl) = 1.e6*Nc0/rho(k) ! convert to mixing ratio
        end if
        tabs(:,:,k) = tabs0(k)
     end do

     if(docloud) call micro_diagnose_m2005()   ! leave this here

  end if

end subroutine micro_init_2005

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux_m2005()

use vars, only: fluxbq, fluxtq

fluxbmk(:,:,:) = 0. ! initialize all fluxes at surface to zero
fluxtmk(:,:,:) = 0. ! initialize all fluxes at top of domain to zero
fluxbmk(:,:,index_water_vapor) = fluxbq(:,:) ! surface qv (latent heat) flux
fluxtmk(:,:,index_water_vapor) = fluxtq(:,:) ! top of domain qv flux

end subroutine micro_flux_m2005

!----------------------------------------------------------------------
!!! compute local microphysics processes (beyond advection and SGS diffusion):
!
!  This is the place where the condensation/sublimation, accretion, coagulation, freezing,
!  melting, etc., that is  all the microphysics processes except for the spatial transport happen.

! IMPORTANT: You need to use the thermodynamic constants like specific heat, or
! specific heat of condensation, gas constant, etc, the same as in file params.f90
! Also, you should assume that the conservative thermodynamic variable during these
! proceses is the liquid/ice water static energy: t = tabs + gz - Lc (qc+qr) - Ls (qi+qs+qg) 
! It should not be changed during all of your point microphysical processes!

subroutine micro_proc_m2005()

use params, only: fac_cond, fac_sub, rgas
use grid, only: z, zi
use vars, only: docloud, t, gamaz, precsfc, precflux, qpfall, tlat, prec_xy, &
     nstep, nstatis, icycle, s_ar, dtfactor, total_water_prec


real, dimension(nzm) :: &
     tmpqcl, tmpqci, tmpqr, tmpqs, tmpqg, tmpqv, &
     tmpncl, tmpnci, tmpnr, tmpns, tmpng,  &
     tmpw, tmpwsub, tmppres, tmpdz, tmptabs, &
     tmtend1d, &
     mtendqcl, mtendqci, mtendqr, mtendqs, mtendqg, mtendqv, &
     mtendncl, mtendnci, mtendnr, mtendns, mtendng,  &
     stendqcl, stendqci, stendqr, stendqs, stendqg, stendqv, &
     stendncl, stendnci, stendnr, stendns, stendng,  &
     effg1d, effr1d, effs1d, effc1d, effi1d


real, dimension(nzm,nmicro_fields) :: stend1d, mtend1d
real :: tmpc, tmpr, tmpi, tmps, tmpg
integer :: i1, i2, j1, j2, i, j, k, m, n

double precision :: tmp_total, tmptot

call t_startf ('micro_proc')

if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
      do i=1,nx
         precsfc(i,j)=0.
      end do
   end do
   do k=1,nzm
      precflux(k) = 0.
   end do
end if

if(dostatis) then ! initialize arrays for statistics
   mfrac(:,:) = 0.
   mtend(:,:) = 0.
   trtau(:,:) = 0.
   qpfall(:)=0.
   tlat(:) = 0.
   tmtend3d(:,:,:) = 0.
end if
stend(:,:) = 0.
mksed(:,:) = 0.

!!$if(doprecip) total_water_prec = total_water_prec + total_water()
 
do j = 1,ny
   do i = 1,nx

      ! zero out mixing ratios of microphysical species
      tmpqv(:) = 0.
      tmpqcl(:) = 0.
      tmpncl(:) = 0.
      tmpqr(:) = 0.
      tmpnr(:) = 0.
      tmpqci(:) = 0.
      tmpnci(:) = 0.
      tmpqs(:) = 0.
      tmpns(:) = 0.
      tmpqg(:) = 0.
      tmpng(:) = 0.

      ! get microphysical quantities in this grid column
      tmpqv(:) = micro_field(i,j,:,iqv)
      tmpqcl(:) = micro_field(i,j,:,iqcl)
      if(dopredictNc) tmpncl(:) = micro_field(i,j,:,incl)
      if(doprecip) then
         tmpqr(:) = micro_field(i,j,:,iqr)
         tmpnr(:) = micro_field(i,j,:,inr)
      end if

      if(doicemicro) then
         tmpqci(:) = micro_field(i,j,:,iqci)
         tmpnci(:) = micro_field(i,j,:,inci)
         tmpqs(:) = micro_field(i,j,:,iqs)
         tmpns(:) = micro_field(i,j,:,ins)
         if(dograupel) then
            tmpqg(:) = micro_field(i,j,:,iqg)
            tmpng(:) = micro_field(i,j,:,ing)
         end if
      end if

      ! get absolute temperature in this column
      tmptabs(:) = t(i,j,:) &           ! liquid water-ice static energy over Cp
           - gamaz(:) &                                   ! potential energy
           + fac_cond * (tmpqcl(:) + tmpqr(:)) &          ! liquid latent energy
           + fac_sub  * (tmpqci(:) + tmpqs(:) + tmpqg(:)) ! ice latent energy

      tmpdz = adz(:)*dz
!      tmpw = 0.5*(w(i,j,1:nzm) + w(i,j,2:nz))  ! MK: changed for stretched grids 
      tmpw = ((zi(2:nz)-z(1:nzm))*w(i,j,1:nzm)+ &
             (z(1:nzm)-zi(1:nzm))*w(i,j,2:nz))/(zi(2:nz)-zi(1:nzm))
      tmpwsub = 0.

      tmppres(:) = 100.*pres(1:nzm)

      i1 = 1 ! dummy variables used by WRF convention in subroutine call
      i2 = 1
      j1 = 1
      j2 = 1

      mtendqv = 0.
      mtendqcl = 0.
      mtendqr = 0.
      mtendqci = 0.
      mtendqs = 0.
      mtendqg = 0.
      mtendncl = 0.
      mtendnr = 0.
      mtendnci = 0.
      mtendns = 0.
      mtendng = 0.

      tmtend1d = 0.

      sfcpcp = 0.
      sfcicepcp = 0.

      effc1d(:) = 10. ! default liquid and ice effective radii
      effi1d(:) = 75.

      ! explanation of variable names:
      !   mtend1d: array of 1d profiles of microphysical tendencies (w/o sed.)
      !   stend1d: array of 1d profiles of sedimentation tendencies for q*
      !   tmp**: on input, current value of **.  On output, new value of **.
      !   eff*1d: one-dim. profile of effective raduis for *
      call m2005micro_graupel(&
           mtendqcl,mtendqci,mtendqs,mtendqr, &
           mtendncl,mtendnci,mtendns,mtendnr, &
           tmpqcl,tmpqci,tmpqs,tmpqr, &
           tmpncl,tmpnci,tmpns,tmpnr, &
           tmtend1d,mtendqv, &
           tmptabs,tmpqv,tmppres,rho,tmpdz,tmpw,tmpwsub, &
           sfcpcp, sfcicepcp, &
           effc1d,effi1d,effs1d,effr1d, &
           dtn, &
           i1,i2, j1,j2, 1,nzm, i1,i2, j1,j2, 1,nzm, &
           mtendqg,mtendng,tmpqg,tmpng,effg1d,stendqg, &
           stendqr,stendqci,stendqs,stendqcl)
 
     ! update microphysical quantities in this grid column
      micro_field(i,j,:,iqv) = tmpqv(:)

      if(doprecip) then
         total_water_prec = total_water_prec + sfcpcp

         ! take care of surface precipitation
         precsfc(i,j) = precsfc(i,j) + sfcpcp/dtn/dz
         prec_xy(i,j) = prec_xy(i,j) + sfcpcp/dtn/dz

         ! update rain
         micro_field(i,j,:,iqr) = tmpqr(:)
         micro_field(i,j,:,inr) = tmpnr(:)
      else
         ! add rain to cloud
         tmpqcl(:) = tmpqcl(:) + tmpqr(:) ! add rain mass back to cloud water
         tmpncl(:) = tmpncl(:) + tmpnr(:) ! add rain number back to cloud water

         ! zero out rain 
         tmpqr(:) = 0.
         tmpnr(:) = 0.

         ! add rain tendencies to cloud
         stendqcl(:) = stendqcl(:) + stendqr(:)
         mtendqcl(:) = mtendqcl(:) + mtendqr(:)
         mtendncl(:) = mtendncl(:) + mtendnr(:)

         ! zero out rain tendencies
         stendqr(:) = 0.
         mtendqr(:) = 0.
         mtendnr(:) = 0.
      end if

      ! update cloud water
      micro_field(i,j,:,iqcl) = tmpqcl(:)
      if(dopredictNc) micro_field(i,j,:,incl) = tmpncl(:)

      ! approx optical depth so that effective radius can be
      !   reconstructed in radiation routine.
      tauliq(i,j,:) = &
           0.0018*rho(:)*dz*adz(:)*tmpqcl(:)/(1.e-20+1.e-6*effc1d(:))

      if(doicemicro) then
         micro_field(i,j,:,iqci) = tmpqci(:)
         micro_field(i,j,:,inci) = tmpnci(:)
         micro_field(i,j,:,iqs) = tmpqs(:)
         micro_field(i,j,:,ins) = tmpns(:)
         if(dograupel) then
            micro_field(i,j,:,iqg) = tmpqg(:)
            micro_field(i,j,:,ing) = tmpng(:)
         end if

         tauice(i,j,:) =&  
              0.0018*rho(:)*dz*adz(:)*tmpqci(:)/(1.e-20+1.e-6*effi1d(:)) &
              + 0.0018*rho(:)*dz*adz(:)*tmpqs(:)/(1.e-20+1.e-6*effs1d(:))
      end if

      !=====================================================
      ! update liquid-ice static energy due to precipitation
      t(i,j,:) = t(i,j,:) &
           - dtn*fac_cond*(stendqcl+stendqr) &
           - dtn*fac_sub*(stendqci+stendqs+stendqg)
      !=====================================================

      if(dostatis) then
         mtend(:,iqv) = mtend(:,iqv) + mtendqv
         mtend(:,iqcl) = mtend(:,iqcl) + mtendqcl
         if(dopredictNc) mtend(:,incl) = mtend(:,incl) + mtendncl
         if(doprecip) then
            mtend(:,iqr) = mtend(:,iqr) + mtendqr
            mtend(:,inr) = mtend(:,inr) + mtendnr
         end if

         if(doicemicro) then
            mtend(:,iqci) = mtend(:,iqci) + mtendqci
            mtend(:,inci) = mtend(:,inci) + mtendnci
            !bloss            stend(:,inci) = stend(:,inci) + stendnci

            mtend(:,iqs) = mtend(:,iqs) + mtendqs
            mtend(:,ins) = mtend(:,ins) + mtendns
            !bloss            stend(:,ins) = stend(:,ins) + stendns

            if(dograupel) then
               mtend(:,iqg) = mtend(:,iqg) + mtendqg
               mtend(:,ing) = mtend(:,ing) + mtendng
               !bloss            stend(:,ing) = stend(:,ing) + stendng
            end if
         end if

         do n = 1,nmicro_fields
            do k = 1,nzm
               if(micro_field(i,j,k,n).ge.1.e-6) mfrac(k,n) = mfrac(k,n)+1.
            end do
         end do

         ! approximate optical depth = 0.0018*lwp/effrad
         !  integrated up to level at which output
         tmpc = 0.
         tmpr = 0.
         tmpi = 0.
         tmps = 0.
         tmpg = 0.

         do k = 1,nzm
            tmpc = tmpc + 0.0018*rho(k)*dz*adz(k)*tmpqcl(k)/(1.e-20+1.e-6*effc1d(k))
            tmpr = tmpr + 0.0018*rho(k)*dz*adz(k)*tmpqr(k)/(1.e-20+1.e-6*effr1d(k))
            trtau(k,iqcl) = trtau(k,iqcl) + tmpc
            if(doprecip) trtau(k,iqr) = trtau(k,iqr) + tmpr

            if(doicemicro) then
               tmpi = tmpi + 0.0018*rho(k)*dz*adz(k)*tmpqci(k)/(1.e-20+1.e-6*effi1d(k))
               tmps = tmps + 0.0018*rho(k)*dz*adz(k)*tmpqs(k)/(1.e-20+1.e-6*effs1d(k))
               tmpg = tmpg + 0.0018*rho(k)*dz*adz(k)*tmpqg(k)/(1.e-20+1.e-6*effg1d(k))

               trtau(k,iqci) = trtau(k,iqci) + tmpi
               trtau(k,iqs) = trtau(k,iqs) + tmps
               trtau(k,iqg) = trtau(k,iqg) + tmpg
            end if
         end do

         tlat(1:nzm) = tlat(1:nzm) &
              - dtn*fac_cond*(stendqcl+stendqr) &
              - dtn*fac_sub*(stendqci+stendqs+stendqg)
         qpfall(1:nzm) = qpfall(1:nzm) + dtn*(stendqr+stendqs+stendqg)

         !bloss: temperature tendency (sensible heating) due to phase changes
         tmtend3d(i,j,1:nzm) = tmtend1d(1:nzm)

      end if ! dostatis

      stend(:,iqcl) = stend(:,iqcl) + stendqcl
      if(doprecip) then
         ! compute surface precipitation area fraction statistics
         if(sfcpcp.gt.1.e-6) s_ar=s_ar+dtfactor

         stend(:,iqr) = stend(:,iqr) + stendqr
      end if

      if(doicemicro) then
         stend(:,iqci) = stend(:,iqci) + stendqci
         stend(:,iqs) = stend(:,iqs) + stendqs
         if(dograupel) stend(:,iqg) = stend(:,iqg) + stendqg
      end if

   end do ! i = 1,nx
end do ! j = 1,ny

! back sedimentation flux out from sedimentation tendencies
tmpc = 0.
do k = 1,nzm
   m = nz-k
   tmpc = tmpc + stend(m,iqcl)*rho(m)*dz*adz(m)
   mksed(m,iqcl) = tmpc
end do
precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqcl)*dtn/dz

if(doprecip) then
   tmpr = 0.
   do k = 1,nzm
      m = nz-k
      tmpr = tmpr + stend(m,iqr)*rho(m)*dz*adz(m)
      mksed(m,iqr) = tmpr
   end do
   precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqr)*dtn/dz
end if

if(doicemicro) then
   tmpi = 0.
   tmps = 0.
   tmpg = 0.
   do k = 1,nzm
      m = nz-k
      tmpi = tmpi + stend(m,iqci)*rho(m)*dz*adz(m)
      tmps = tmps + stend(m,iqs)*rho(m)*dz*adz(m)
      tmpg = tmpg + stend(m,iqg)*rho(m)*dz*adz(m)
      mksed(m,iqci) = tmpi
      mksed(m,iqs) = tmps
      mksed(m,iqg) = tmpg
   end do
   precflux(1:nzm) = precflux(1:nzm) &
        - (mksed(:,iqci) + mksed(:,iqs) + mksed(:,iqg))*dtn/dz
end if

!!$if(doprecip) total_water_prec = total_water_prec - total_water()

if (docloud)  call micro_diagnose_m2005()   ! leave this line here

call t_stopf ('micro_proc')

end subroutine micro_proc_m2005

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and radiation:
!
!  This is the pace where the microphysics field that SAM actually cares about
!  are diagnosed.

subroutine micro_diagnose_m2005()

use vars

real omn, omp
integer i,j,k

! water vapor
qv(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqv)

! cloud liquid water
qcl(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqcl)

! rain water
if(doprecip) qpl(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqr)

! cloud ice 
if(doicemicro) then
   qci(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqci)

   if(dograupel) then
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs) &
           + micro_field(1:nx,1:ny,1:nzm,iqg)
   else
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs)
   end if
end if

end subroutine micro_diagnose_m2005

!----------------------------------------------------------------------
!!! functions to compute terminal velocity for precipitating variables:
!
! you need supply functions to compute terminal velocity for all of your 
! precipitating prognostic variables. Note that all functions should
! compute vertical velocity given two microphysics parameters var1, var2, 
! and temperature, and water vapor (single values, not arrays). Var1 and var2 
! are some microphysics variables like water content and concentration.
! Don't change the number of arguments or their meaning!

!!$real function term_vel_qr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_qr
!!$
!!$real function term_vel_Nr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_Nr
!!$
!!$real function term_vel_qs(qs,ns,tabs,rho)
!!$! .......  
!!$end function term_vel_qs

! etc.

!----------------------------------------------------------------------
!!! compute sedimentation 
!
!  The perpose of this subroutine is to prepare variables needed to call
! the precip_all() for each of the falling hydrometeor varibles
subroutine micro_precip_fall_m2005()

use vars, only : s_ar

! before calling precip_fall() for each of falling prognostic variables,
! you need to set hydro_type and omega(:,:,:) variables.
! hydro_type can have four values:
! 0 - variable is liquid water mixing ratio
! 1 - hydrometeor is ice mixing ratio
! 2 - hydrometeor is mixture-of-liquid-and-ice mixing ratio. (As in original SAM microphysics).
! 3 - variable is not mixing ratio, but, for example, rain drop concentration
! OMEGA(:,:,:) is used only for hydro_type=2, and is the fraction of liquid phase (0-1).
! for hour hypothetical case, there is no mixed hydrometeor, so omega is not actually used.

integer hydro_type
real omega(nx,ny,nzm) 

integer i,j,k

return ! do not need this routine -- sedimentation done in m2005micro.

!!$! Initialize arrays that accumulate surface precipitation flux
!!$
!!$ if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
!!$   do j=1,ny
!!$    do i=1,nx
!!$     precsfc(i,j)=0.
!!$    end do
!!$   end do
!!$   do k=1,nzm
!!$    precflux(k) = 0.
!!$   end do
!!$ end if
!!$
!!$ do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
!!$    qpfall(k)=0.
!!$    tlat(k) = 0.
!!$ end do
!!$   
!!$! Compute sedimentation of falling variables:
!!$
!!$ hydro_type=0
!!$ call precip_fall(qr, term_vel_qr, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Nr, term_vel_Nr, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qs, term_vel_qs, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ns, term_vel_Ns, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qg, term_vel_qg, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ng, term_vel_Ng, hydro_type, omega)
!!$
!!$! compute surface precipitation area fraction statistics
!!$ do j=1,ny
!!$   do i=1,nx
!!$     if((qr(i,j,1)+qs(i,j,1)+qg(i,j,1).gt.1.e-6) s_ar=s_ar+dtfactor
!!$   end do
!!$ end do

end subroutine micro_precip_fall_m2005

end module micro_m2005



