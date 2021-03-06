! Fortran code for various nucleation parameterisations available
! this is for testign before implementation in MESSy
! see the functions for explanations of parameters. We use functions here as we expect a 
! particle formation rate, i.e. one number. 
! author: sebastian ehrhart
! email: s.ehrhart@mpic.de
! date created: 11 July 2016
! modified: 15 Nov 2016

module messy_nan

  use messy_main_constants_mem, only: dp, strlen_ulong, strlen_long, strlen_medium

  implicit none
  private
  
  ! only subroutines are declared public here. Variables are set public when they are declared
  public :: read_nucleation_nml   ! read the namelist file
  public :: initialize_nucleation ! initialise and allocate variables and parameters used for the nucleation
  public :: clean_up_nucleation   ! deallocate variables
!  public :: pgasi_index_gmxe     ! assign vapour index for pgasi in gmxe
  public :: bio_ox_org_pr 
  public :: hom_pr
! ju_ak_20191107  public :: driver_NPF
  public :: nucleation_driver 
! ju_ak_20191107  public :: driver_gmxe_nucleation
!  public :: gmxe_assign_index
  ! ju_ak_20191107+
  public :: driver_vehkamaeki
  public :: GR_Nieminen2010
  public :: Anttila2010
  public :: Gordon2016
  ! ju_ak_20191107-
 

  ! GLOBAL PARAMETERS
  character(len=*), parameter, public :: modstr = 'nan'
  character(len=*), parameter, public :: modver = '1.0'

  ! fixed parameters and parameters that can be changed by readin nucleation.nml 
  character(len=strlen_ulong), dimension(1:5), parameter :: default_names=(/'H2SO4   ',&
                                                                            'NH3     ',& 
                                                                            'LAMINE  ',& 
                                                                            'HOMOH   ',& 
                                                                            'HOMO3   '/)
!                                                                            'XIODINE '/) ! XIODINE is a place holder for marine nucleation species

  integer, public, save :: nnucspec = 1 ! number of species involved in nucleation

  integer, public, save :: insa=1, & ! Index Nucleating species Sulphuric Acid (H2SO4)
                          innh3=0, & ! Index Nucleating species NH3
                          indma=0, & ! Index Nucleating species Dimethylamine, or other amines
                          inoxo1=0, & ! Index Nucleating species Biogenic oxidised organics, from Riccobono et al 2014
                          inoxo2=0, &  ! Index Nucleating species highly oxidised molecules from Kirkby et al et al 2016
                          inmar=0  ! marine nucleation species index

  character(len=strlen_medium), dimension(:), allocatable, public, save :: vapour_names ! for loops
  character(len=strlen_medium), dimension(:), allocatable, public, save :: condph_names !
  character(len=strlen_medium), public, save :: sulphuric_acid='', ammonia='', amines='', HOMOH='', HOMO3='' ! for easier input
  character(len=strlen_medium), dimension(1:5), public, save :: vapour_aero_names =(/'H2SO4 ', &
                                                                                     'NH3   ', &
                                                                                     'LAMINE', &
                                                                                     'HOMOH ',  &
                                                                                     'HOMO3 '/)
  character(len=strlen_medium), public, save :: nuclmethod='multi'

  integer, public, parameter ::     & ! ID number for nucleation mechanism of
                        idvehk =1,  & ! Vehkamaeki
                        idkul  =2,  & ! Kulmala
                        idmeri =3,  & ! Merikanto
                        iddunbn=4,  & ! Dunne binary neutral
                        iddunbi=5,  & ! Dunne binary ion induced
                        idduntn=6,  & ! Dunne ternary neutral
                        iddunti=7,  & ! Dunne ternary charged
                        idric  =8,  & ! Riccobono
                        idkirkn=9,  & ! Kirkby neutral
                        idkirki=10, & ! Kirkby ion induced
                        idalmei=11, & ! Almeida amines
                        idmar=12      ! marine nucleation ! currently unused

   character(len=strlen_medium), dimension(1:12), public, parameter :: nucchannames = &
                                                       (/'Vehkamaeki          ',&
                                                         'Kulmala             ',& 
                                                         'Merikanto           ',&
                                                         'Dunne_binary_n      ',&
                                                         'Dunne_binary_IIN    ',&
                                                         'Dunne_ternary_n     ',&
                                                         'Dunne_ternary_IIN   ',&
                                                         'Riccobono           ',&
                                                         'Kirkby_n            ',&
                                                         'Kirkby_IIN          ',&
                                                         'Almeida             ',&
                                                         'Marine              ' /)

  integer, public, parameter, dimension(3) :: alliin = (/iddunbi, iddunti, idkirki/) ! list of all ion induced nucleation channels

  ! variabels to determine via namelist or from the calling programme or just leave them as set here
  integer, public, save :: nnucchan = 1 ! number of nucleation channels
  integer, allocatable, public, save :: chid(:) ! id number for the nucleation channel, will have dimensions 1:id
  integer, allocatable, save :: iinchind(:) !index of ion induced channels in chid
  logical, allocatable, save :: lnchan(:,:)
  logical, public, save :: lselnuc(1:11) = .FALSE. ! choose the nucleation channel, index of this vector corresponds 
                                              ! to ID number of the nucleation mechanism 

  integer, allocatable, public, save :: ip4d(:) ! stores the index of nucleating species for calling aerosol model


  real(kind=dp),public, dimension(1:5), save :: kirkby_a = & 
  & (/0.04001_dp, 1.848_dp, 0.001366_dp, 1.566_dp, 0.1863_dp/) ! values from Kirkby et al 2016
  real(kind=dp), public, save :: ricco_k = 5.45d-19 ! 3.27d-21 ! prefactor for Riccobono et al 2014 parameterisation
  real(kind=dp), public, save :: dunne_pbn=3.95_dp, dunne_ubn=9.70_dp, dunne_vbn=12.6_dp, dunne_wbn=-0.00707_dp, &
  & dunne_ptn=2.89_dp, dunne_utn=182_dp, dunne_vtn=1.20_dp, dunne_wtn=-4.19_dp, dunne_pAn=8.00_dp,       &
  & dunne_an=1.6d-6, dunne_pbi=3.37_dp, dunne_ubi=-11.5_dp, dunne_vbi=25.5_dp, dunne_wbi=0.181_dp,    &
  & dunne_pti=3.14_dp, dunne_uti=-23.8_dp, dunne_vti=37.0_dp, dunne_wti=0.227_dp, dunne_pAi=3.07_dp,     &
  & dunne_ai=0.00485_dp ! fit parameters from Dunne et al SI
  real(kind=dp), public, save :: almeida_lim=2.0e8_dp, almeida_k(2)=(/2.08e-25_dp/2.5e7_dp, 1.93e-28_dp/1.798947e+32_dp/), &
  & almeida_pa(2)=(/1.0_dp, 4.36_dp/), almeida_ps(2)=(/3.7_dp, 3.7_dp/) ! Parameterisation based on data from Almeida et al
  real(kind=dp), public, save :: frhc1=1.5_dp, frhc2=0.045_dp ! RH dependence Dunne et al 2016
  real(kind=dp), public, save :: BtdOrg = 0.0_dp ! temperature dependence organic nucleation channels
  logical, public, save :: l_org_Tdep=.FALSE. ! logical value to determine if Temperature dependence for organic species is wanted
  logical, public, save :: l_Dunne_RH =.FALSE. ! RH depence in Dunne et al?

  ! the numbers below should be available in a module for all of MESSy, but GMXe uses it's own
  ! the part below is copied frome messy_gmxe.f90, these are physcial and numerical constants. 
  REAL(dp), PARAMETER :: ZERO = 0._dp ! REAL zero
  REAL(dp), PARAMETER :: bk     = 1.38e-16_dp  ! Bolzman constant [erg / K]
  REAL(dp), PARAMETER, PUBLIC :: zeps=1.e-20_dp       ! EPSILON(1._dp) (determined during initialization), ! was public #seb
  REAL(dp), PARAMETER :: cmin_epsilon   = 1.e-18_dp ! Minium value

  contains


subroutine read_nucleation_nml(status, iou)
  use messy_main_tools, only : read_nml_open, read_nml_check, read_nml_close
  implicit none
  integer, intent(out) :: status ! error status: 0 = no error, 1 = an error
  integer, intent(in)  :: iou    ! I/O unit

  NAMELIST /CTRL/ lselnuc, nnucspec, nuclmethod
  NAMELIST /NUC/ sulphuric_acid, ammonia, amines, HOMOH, HOMO3, vapour_aero_names 
  NAMELIST /PARAM/ dunne_pbn, dunne_ubn, dunne_vbn, dunne_wbn, &
  & dunne_ptn, dunne_utn, dunne_vtn, dunne_wtn, dunne_pAn, dunne_an, &
  & dunne_pbi, dunne_ubi, dunne_vbi, dunne_wbi, dunne_pti, dunne_uti, &
  & dunne_vti, dunne_wti, dunne_pAi, dunne_ai, & ! Dunne et al 2016
  & almeida_lim, almeida_k, almeida_pa, almeida_ps, &
  & kirkby_a, &  ! Kirkby et al 2016
  & ricco_k, & ! Riccobono et al 2014
  & frhc1, frhc2, l_Dunne_RH, BtdOrg, l_org_Tdep ! additional parametess
  
  ! local
  character(len=*), parameter :: substr='read_nucleation_nml'
  logical :: lex
  integer :: fstat

  status=1

  ! TODO: check if a rewind statement has to be issued to make the namelist
  ! reading independent of position

  call read_nml_open(lex, substr, iou, 'CTRL', modstr)
  if(.not.lex) return ! can't open the file

  !---------- select nucleation channels
  read(iou, nml=CTRL, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'CTRL', modstr)
  if(fstat /= 0) return ! error while reading the namelist
  

  !---------- read number of nucleating species and index
  read(iou, nml=NUC, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'NUC', modstr)
  if(fstat /= 0) return ! error while reading the namelist

  !---------- read in parameterisation
  read(iou, nml=PARAM, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'PARAM', modstr)
  if(fstat /= 0) return ! error while reading the namelist

  call read_nml_close(substr, iou, modstr)
  status = 0 ! all fine
  
end subroutine read_nucleation_nml


subroutine initialize_nucleation()
  ! define the nucleation schemes used
  implicit none
  integer :: i, k, ki
  logical, allocatable :: lcheck(:)


  insa = 0 

  nnucchan = count(lselnuc)
  write(*,*) "nnucchan= ", nnucchan
  
  allocate(lcheck(nnucspec))
  lcheck = .true. ! no index assigned yet


  ! "user friendly" approach:
  ! check if name was given by namelist, if not use default name
  ! the name can be used to inquire index numbers of tracers from the calling programme. 
  ! If this feature is not needed the user is free to choose any name, such as Alice or Bob

  allocate(vapour_names(nnucspec))
  allocate(condph_names(nnucspec))
  k = 0
  ! check for H2SO4
  if(lselnuc(idkul) .or. lselnuc(idvehk) .or. lselnuc(iddunbn) .or. lselnuc(iddunbi) &
   & .or. lselnuc(idduntn) .or. lselnuc(iddunti) .or. lselnuc(idric) .or. lselnuc(idalmei)) then
     k = k+1
     if(len_trim(sulphuric_acid) == 0) then
       vapour_names(k) = default_names(1) ! H2SO4
!       condph_names(k) = trim(default_names(1)) ! // "_ns"
     else
       vapour_names(k) = sulphuric_acid
!       condph_names(k) = trim(sulphuric_acid) !// "_ns" 
     end if
     condph_names(k) = trim(vapour_aero_names(1))
     insa = k
     write(*,*) "insa=", insa
     lcheck(k) = .false.
  end if
  ! check for NH3
  if(any((/lselnuc(idmeri), lselnuc(idduntn), lselnuc(iddunti) /))) then
    k = k+1
    if(len_trim(ammonia) == 0) then
       vapour_names(k) = default_names(2) ! NH3
    else
       vapour_names(k) = ammonia
    end if
    condph_names(k) = trim(vapour_aero_names(2))
    innh3 = k
    write(*,*) "innh3=", innh3
    lcheck(k) = .false.
  end if
  ! Amines
  if(lselnuc(idalmei)) then
    k = k + 1
    if(len_trim(amines) == 0) then 
      vapour_names(k)  = default_names(3) ! Amines
    else
      vapour_names(k) = amines
    end if
    condph_names(k) = trim(vapour_aero_names(3))
    indma = k
    write(*,*) "indma=", indma
    lcheck(k) = .false.
  end if
  ! Oxidised organics from Riccobono et al 2014
  if(any(lselnuc((/idric, idkirkn, idkirki/)))) then
    k = k+1
    if(len_trim(HOMOH) == 0) then
      vapour_names(k) = default_names(4)
    else
      vapour_names(k) = HOMOH
    end if
    condph_names(k) = trim(vapour_aero_names(4))
    inoxo1 = k
    write(*,*) "inoxo1=", inoxo1
    lcheck(k) = .false.
  end if
  ! Oxidised organics from Kirkby et al 2016
  if(any(lselnuc((/idkirkn, idkirki/)))) then
!   if(any(lselnuc((/idric, idkirkn, idkirki/)))) then
    k = k+1
    if(len_trim(HOMO3) == 0) then 
      vapour_names(k) = default_names(5)
    else
      vapour_names(k) = HOMO3
    end if
    condph_names(k) = trim(vapour_aero_names(5))
    inoxo2 = k
    write(*,*) "inoxo2=", inoxo2
    lcheck(k) = .false.
  end if 

  ! Check if any species was forgotten
  if(any(lcheck)) then
    write(*,*) "At least one nucleating species can not be assigned"
    do i=1, size(lcheck)
      if(lcheck(i)) write(*,*) "Species not found"
    end do
! abort useful if run as box model, but not if run in parallel on EMAC
!    call abort()
  end if

  deallocate(lcheck)


  ! how many ion induced channels
  k=0
  do i=1, size(alliin)
    if(lselnuc(alliin(i))) k=k+1
  end do
  write(*,*) k, " ion induced nucleation channels chosen."
  allocate(iinchind(k))

  ! create the list of chosen nucleation channels
  allocate(chid(1:nnucchan))
  k  = 1
  ki = 1
  do i=1, size(lselnuc)
    if(lselnuc(i)) then 
      chid(k) = i
      ! is it an ion induced channel? If yes put it in the list
      if(any(i == alliin)) then 
        iinchind(ki) = k
        write(*,*) "ION INDUCED CHANNEL No : ", ki, " has ID ", i, "assigned and index ", k,  &
       & "in CHID"
        ki = ki + 1
      end if
      k = k + 1
    end if
  end do 

  ! define nucleation channel matrix
  allocate(lnchan(nnucchan, nnucspec))
  lnchan = .FALSE. ! default
  do i=1,size(chid)
    select case (chid(i))
    case(idkul,idvehk,iddunbn) ! these contain H2SO4
      lnchan(i,insa)   = .TRUE.
    case(iddunbi)! H2SO4 and negative ions
      lnchan(i,insa)   = .TRUE.
    case(idmeri,idduntn) ! contains NH3 and H2SO4
      lnchan(i,insa)   = .TRUE.
      lnchan(i,innh3)  = .TRUE.
    case(iddunti) !i NH3, H2SO4 and negative ions
      lnchan(i,insa)   = .TRUE.
      lnchan(i,innh3)  = .TRUE.
    case(idalmei) ! Amines
      lnchan(i,insa)   = .TRUE.
      lnchan(i,indma)  = .TRUE.
    case(idkirkn) ! only homs
      lnchan(i,inoxo1) = .TRUE.
      lnchan(i,inoxo2) = .TRUE.
    case(idkirki) ! homs, positive and negative ions
      lnchan(i,inoxo1) = .TRUE.
      lnchan(i,inoxo2) = .TRUE.
    case(idric) ! oxidised organics only OH products Riccobono
      lnchan(i,inoxo1) = .TRUE.
      lnchan(i,insa)   = .TRUE.
    end select
  end do

  do i=1, ubound(lnchan,2)
    write(*,*) "LNCHAN= ", (lnchan(k,i), k=1,ubound(lnchan,1))
  end do

  ! define index of nucleating species
  allocate(ip4d(nnucspec))
  ip4d = 0

end subroutine initialize_nucleation


subroutine clean_up_nucleation()
  implicit none
  ! deallocate memory
  deallocate(chid)
  deallocate(lnchan)
  deallocate(iinchind)

end subroutine clean_up_nucleation



!--------------- NUCLEATION CALCULATIONS ------------------------------------------------------------------------------

! BioOxOrg production rate
real(kind=dp) elemental function bio_ox_org_pr(mt,oh,tk)
  implicit none
  real(kind=dp), intent(in) :: mt ! lumped mono terpene concentration (cm**-3)
  real(kind=dp), intent(in) :: oh ! OH concentration (cm**-3)
!  real(kind=dp), intent(in) :: cs ! condensational sink (s**-1)
  real(kind=dp), intent(in) :: tk ! temperature in Kelvin 

  bio_ox_org_pr = 1.2d-11*exp(444.0_dp/tk) * mt*oh
  
end function bio_ox_org_pr


! hom production rate
real(kind=dp) elemental function hom_pr(mt, oh, o3, tk)
  implicit none
  real(kind=dp), intent(in) :: mt ! lumped mono terpene concentration (cm**-3)
  real(kind=dp), intent(in) :: oh ! [OH] (cm**-3)
  real(kind=dp), intent(in) :: o3 ! [O3] (cm**-3)
!  real(kind=dp), intent(in) :: cs ! cond. sink (s**-1)
  real(kind=dp), intent(in) :: tk ! T in K

  !          yield  *  rate constant             
  hom_pr =   1.2d-2 * 8.05d-16*exp(-640.0_dp/tk) * mt * o3  + &
         &   0.6d-2 * 1.2d-11* exp(444.0_dp/tk)  * mt * oh 

end function hom_pr


! the simple parameterisation given in Riccobono et al 2014 10.1126/science.1243527
! this parameterisation lacks a temperature and RH dependence, it is also a 
! CLOUD parameterisations without an ion component. 
! the original publication defines BioOxOrg as products of PD + OH
! this version uses only BioOxOrg as input
real(kind=dp) elemental function riccobono2014b(sa, boxorg)! steady state formation rate of 1.7 nm particles [ #particle cm**-3 s**-1 ]
  implicit none
  real(kind=dp), intent(in) :: sa ! number concentration of H2SO4 in air [ #molecules cm**-3 ]
  real(kind=dp), intent(in) :: boxorg ! concentration of BioOxOrg in Riccobono et al  [ #molecules cm**-3 ]

  ! second order in SA first order in BioOxOrg
  riccobono2014b = ricco_k*sa*sa*boxorg

end function riccobono2014b


! parameterisation from Kirkby et al 2016 http://dx.doi.org/10.1038/nature17953 
! this parameterisation lacks temperature and humidity dependence but has a a neutral
! and ion induced term. The following function is for the neutral channel. 
! the original publication defines HOMs as production of alpha-pinene with O3 and OH
real(kind=dp) elemental function kirkby2016(hom, a, b, c) result(y)
  implicit none
  real(kind=dp), intent(in) :: hom ! Highly Oxidised Molecule, as defined in Kirkby et al [ 1e7 #molecules cm**-3 ]
  real(kind=dp), intent(in) :: a, b, c ! fitted parameter

!  if(hom > zeps) then
!    y = kirkby_a(1)*hom**(kirkby_a(2) + kirkby_a(5)/hom)
    y = a*hom**(b + c/hom)
!  else
!    y = 0.0_dp
!  end if 

end function kirkby2016


! parameterisation of Amine-sulphuric acid nucleation based on Almeida et al 2014
! The parameterisation itself is described in the SI of Dunne et al 2016
real(kind=dp) pure function almeida2014(amine,sa)
  implicit none
  real(kind=dp), intent(in) :: amine, sa ! amine and sulphuric acid concentration

  if(amine >  almeida_lim) almeida2014 = almeida_k(1) * amine**almeida_pa(1) * sa**almeida_ps(1)
  if(amine <= almeida_lim) almeida2014 = almeida_k(2) * amine**almeida_pa(2) * sa**almeida_ps(2)

end function almeida2014


!! generic nucleation function that (in the final implementation) takes the fit parameters from a namelist
!! allows changing the fit parameter easily. The paramters are variables of the function call.
!pure function pl_nucleation(sp, k, a)
!  implicit none
!  real(kind=dp), intent(in) :: sp(:)
!  real(kind=dp), intent(in) :: a(:)
!  real(kind=dp), intent(in) :: k
!  integer :: i
!  real(kind=dp) pl_nucleation
! 
!  pl_nucleation = k ! pl_nucleatio = k * sp1**a1 * sp2**a2 * ...
!  do i=1, size(sp) 
!    pl_nucleation = pl_nucleation * sp(i)**a(i)
!  end do
!
!end function pl_nucleation


! prefactors for Dunne et al 2016 parameterisation
real(kind=dp) elemental function Dunne_lnkxy(t, u, v, w)
  implicit none
  real(kind=dp), intent(in) :: t, u, v, w

  Dunne_lnkxy = u - exp(v*(t/1000.0_dp - w))

end function Dunne_lnkxy


! Dunne et al 2016 ternary concentration dependence parameterisation
real(kind=dp) elemental function Dunne_fy(sa, nh3, a, pt, pA)
  implicit none
  real(kind=dp), intent(in) :: sa, nh3, a, pt, pA
  real(kind=dp) :: sapt, nh3pA, anh3pA

  sapt = sa**pt
  nh3pA = nh3**pA
  anh3pA = a*nh3pA
  if(sapt >= anh3pA*1.0d3) then  
   Dunne_fy = nh3**(1.0_dp+pA)
  else if(sapt*1.0d3 <= anh3pA) then
   Dunne_fy = nh3*sapt/a
  else
   Dunne_fy = nh3*sapt / (a + sapt/nh3pA)
  end if
end function Dunne_fy


! binary neutral part of the Dunne et al 2016 parameterisation 
real(kind=dp) elemental function Dunne2016binary(sa,t)
  implicit none
  real(kind=dp), intent(in) :: sa, t ! [H2SO4] (cm-3), temperature (K)
  real(kind=dp) :: k

  k = exp(Dunne_lnkxy(t, dunne_ubn, dunne_vbn, dunne_wbn))
  Dunne2016binary = k * sa**dunne_pbn 

end function Dunne2016binary


! binary charged part of the Dunne et al 2016 parameterisation 
real(kind=dp) elemental function Dunne2016binaryIon(sa, t) ! this formation rate has to be multiplied by the ion concentration
  implicit none
  real(kind=dp), intent(in) :: sa, t ! [H2SO4] (cm-3), [negative ions] (cm-3), temperature (K) 
  real(kind=dp) :: k

  k = exp(Dunne_lnkxy(t, dunne_ubi, dunne_vbi, dunne_wbi))
  Dunne2016binaryIon = k * sa**dunne_pbi

end function Dunne2016binaryIon


! ternary neutral particle formation rate Dunne et al 2016
real(kind=dp) elemental function Dunne2016ternary(sa,nh3,t)
  implicit none
  real(kind=dp), intent(in) :: sa, nh3, t ! [H2SO4] (cm-3), [NH3] (cm-3), temperature (K)
  real(kind=dp) :: k
  ! out

  k = exp(Dunne_lnkxy(t, dunne_utn, dunne_vtn, dunne_wtn))
  Dunne2016ternary = k*Dunne_fy(sa, nh3, dunne_an, dunne_ptn, dunne_pAn)

end function Dunne2016ternary


! ternary neutral particle formation rate Dunne et al 2016
real(kind=dp) elemental function Dunne2016ternaryIon(sa,nh3,t) ! this formation rate has to be multiplied by the ion concentration
  implicit none
  real(kind=dp), intent(in) :: sa, nh3, t ! [H2SO4] (cm-3), [NH3] (cm-3), temperature (K)
  real(kind=dp) :: k

  k = exp(Dunne_lnkxy(t, dunne_uti, dunne_vti, dunne_wti))
  Dunne2016ternaryIon = k*Dunne_fy(sa, nh3, dunne_ai, dunne_pti, dunne_pAi)

end function Dunne2016ternaryIon


!  ! might or might not become Merikanto 2007, this parameterisation is only valid above 235 K
! ******************** ternary_parameterization.f90 *********************
! AUTHOR: J. Merikanto
! 
! J. Merikanto, I. Napari, H. Vehkamaki, T. Anttila, and M. Kulmala:
! "New parameterization of sulfuric acid-ammonia-water ternary nucleation
! rates at tropospheric conditions", Journal of Geophysical Research 
! - Atmospheres, 112, D15207, 2007.  
!
! Fortran 90 subroutine that calculates the parameterized composition 
! and nucleation rate of critical clusters in h2o-h2so4-nh3 nucleation
!
! WARNING: The fit should not be used outside its limits of validity
! (limits indicated below)
!
! IN:
! T:     temperature (K), limits 235-295 K
! rh:    relative humidity as fraction (eg. 0.5=50%) limits 0.05-0.95
! c2:    sulfuric acid concentration (molecules/cm3) 
!        limits 5x10**4 - 10**9 molecules/cm3
! c3:    ammonia mixing ratio (ppt) limits 0.1 - 1000 ppt
!
! OUT:
! J_log: logarithm of nucleation rate (1/(s cm3)) ! output changed to nucleation rate, not logarithm #seb
! ntot:  total number of molecules in the critical cluster ! removed #seb
! nacid: number of sulfuric acid molecules in the critical cluster
! namm:  number of ammonia molecules in the critical cluster
! r:     radius of the critical cluster (nm)
! ***********************************************************************

  SUBROUTINE merikanto2007(t,rh,c2,c3,J,nacid,namm)
   IMPLICIT NONE
   real(kind=dp), intent(in) :: t,rh,c2,c3
   real(kind=dp), intent(out) :: J,nacid,namm
   real(kind=dp) :: J_log ,t_onset, ntot, r
   
   t_onset=143.6002929064716 + 1.0178856665693992*rh + &
   &  10.196398812974294*Log(c2) - &
   &  0.1849879416839113*Log(c2)**2 - 17.161783213150173*Log(c3) + &
   &  (109.92469248546053*Log(c3))/Log(c2) + &
   &  0.7734119613144357*Log(c2)*Log(c3) - 0.15576469879527022*Log(c3)**2
   
   if(t_onset.gt.t) then 
   
      J_log=-12.861848898625231 + 4.905527742256349*c3 - & 
   &  358.2337705052991*rh - & 
   &  0.05463019231872484*c3*t + 4.8630382337426985*rh*t + &
   &  0.00020258394697064567*c3*t**2 - 0.02175548069741675*rh*t**2 - &
   &  2.502406532869512e-7*c3*t**3 + 0.00003212869941055865*rh*t**3 - &
   &  4.39129415725234e6/Log(c2)**2 + (56383.93843154586*t)/Log(c2)**2 -& 
   &  (239.835990963361*t**2)/Log(c2)**2 + &
   &  (0.33765136625580167*t**3)/Log(c2)**2 - &
   &  (629.7882041830943*rh)/(c3**3*Log(c2)) + &
   &  (7.772806552631709*rh*t)/(c3**3*Log(c2)) - &
   &  (0.031974053936299256*rh*t**2)/(c3**3*Log(c2)) + &
   &  (0.00004383764128775082*rh*t**3)/(c3**3*Log(c2)) + &
   &  1200.472096232311*Log(c2) - 17.37107890065621*t*Log(c2) + &
   &  0.08170681335921742*t**2*Log(c2) - &
   &  0.00012534476159729881*t**3*Log(c2) - &
   &  14.833042158178936*Log(c2)**2 + 0.2932631303555295*t*Log(c2)**2 - &
   &  0.0016497524241142845*t**2*Log(c2)**2 + &
   &  2.844074805239367e-6*t**3*Log(c2)**2 - 231375.56676032578*Log(c3) - &
   &  100.21645273730675*rh*Log(c3) + 2919.2852552424706*t*Log(c3) + &
   &  0.977886555834732*rh*t*Log(c3) - 12.286497122264588*t**2*Log(c3) - &
   &  0.0030511783284506377*rh*t**2*Log(c3) + &
   &  0.017249301826661612*t**3*Log(c3) + &
   &  2.967320346100855e-6*rh*t**3*Log(c3) + &
   &  (2.360931724951942e6*Log(c3))/Log(c2) - &
   &  (29752.130254319443*t*Log(c3))/Log(c2) + &
   &  (125.04965118142027*t**2*Log(c3))/Log(c2) - &
   &  (0.1752996881934318*t**3*Log(c3))/Log(c2) + &
   &  5599.912337254629*Log(c2)*Log(c3) - &
   &  70.70896612937771*t*Log(c2)*Log(c3) + &
   &  0.2978801613269466*t**2*Log(c2)*Log(c3) - &
   &  0.00041866525019504*t**3*Log(c2)*Log(c3) + &
   &  75061.15281456841*Log(c3)**2 - &
   &  931.8802278173565*t*Log(c3)**2 + 3.863266220840964*t**2*Log(c3)**2 - &
   &  0.005349472062284983*t**3*Log(c3)**2 - &
   &  (732006.8180571689*Log(c3)**2)/Log(c2) + &
   &  (9100.06398573816*t*Log(c3)**2)/Log(c2) - &
   &  (37.771091915932004*t**2*Log(c3)**2)/Log(c2) + &
   &  (0.05235455395566905*t**3*Log(c3)**2)/Log(c2) - &
   &  1911.0303773001353*Log(c2)*Log(c3)**2 + &
   &  23.6903969622286*t*Log(c2)*Log(c3)**2 - &
   &  0.09807872005428583*t**2*Log(c2)*Log(c3)**2 + &
   &  0.00013564560238552576*t**3*Log(c2)*Log(c3)**2 - &
   &  3180.5610833308*Log(c3)**3 + 39.08268568672095*t*Log(c3)**3 - &
   &  0.16048521066690752*t**2*Log(c3)**3 + &
   &  0.00022031380023793877*t**3*Log(c3)**3 + &
   &  (40751.075322248245*Log(c3)**3)/Log(c2) - &
   &  (501.66977622013934*t*Log(c3)**3)/Log(c2) + &
   &  (2.063469732254135*t**2*Log(c3)**3)/Log(c2) - &
   &  (0.002836873785758324*t**3*Log(c3)**3)/Log(c2) + &
   &  2.792313345723013*Log(c2)**2*Log(c3)**3 - &
   &  0.03422552111802899*t*Log(c2)**2*Log(c3)**3 + &
   &  0.00014019195277521142*t**2*Log(c2)**2*Log(c3)**3 - &
   &  1.9201227328396297e-7*t**3*Log(c2)**2*Log(c3)**3 - &
   &  980.923146020468*Log(rh) + 10.054155220444462*t*Log(rh) - &
   &  0.03306644502023841*t**2*Log(rh) + &
   &  0.000034274041225891804*t**3*Log(rh) + &
   &  (16597.75554295064*Log(rh))/Log(c2) - &
   &  (175.2365504237746*t*Log(rh))/Log(c2) + &
   &  (0.6033215603167458*t**2*Log(rh))/Log(c2) - &
   &  (0.0006731787599587544*t**3*Log(rh))/Log(c2) - &
   &  89.38961120336789*Log(c3)*Log(rh) + &
   &  1.153344219304926*t*Log(c3)*Log(rh) - &
   &  0.004954549700267233*t**2*Log(c3)*Log(rh) + &
   &  7.096309866238719e-6*t**3*Log(c3)*Log(rh) + &
   &  3.1712136610383244*Log(c3)**3*Log(rh) - &
   &  0.037822330602328806*t*Log(c3)**3*Log(rh) + &
   &  0.0001500555743561457*t**2*Log(c3)**3*Log(rh) - &
   &  1.9828365865570703e-7*t**3*Log(c3)**3*Log(rh)
   
      J=exp(J_log)
   
      ntot=57.40091052369212 - 0.2996341884645408*t + &
   &  0.0007395477768531926*t**2 - &
   &  5.090604835032423*Log(c2) + 0.011016634044531128*t*Log(c2) + &
   &  0.06750032251225707*Log(c2)**2 - 0.8102831333223962*Log(c3) + &
   &  0.015905081275952426*t*Log(c3) - &
   &  0.2044174683159531*Log(c2)*Log(c3) + &
   &  0.08918159167625832*Log(c3)**2 - &
   &  0.0004969033586666147*t*Log(c3)**2 + &
   &  0.005704394549007816*Log(c3)**3 + 3.4098703903474368*Log(J) - &
   &  0.014916956508210809*t*Log(J) + 0.08459090011666293*Log(c3)*Log(J) - &
   &  0.00014800625143907616*t*Log(c3)*Log(J) + &
   &  0.00503804694656905*Log(J)**2
    
      r=3.2888553966535506e-10 - 3.374171768439839e-12*t + &
   &  1.8347359507774313e-14*t**2 + 2.5419844298881856e-12*Log(c2) - &
   &  9.498107643050827e-14*t*Log(c2) + 7.446266520834559e-13*Log(c2)**2 + &
   &  2.4303397746137294e-11*Log(c3) + 1.589324325956633e-14*t*Log(c3) - &
   &  2.034596219775266e-12*Log(c2)*Log(c3) - &
   &  5.59303954457172e-13*Log(c3)**2 - &
   &  4.889507104645867e-16*t*Log(c3)**2 + &
   &  1.3847024107506764e-13*Log(c3)**3 + &
   &  4.141077193427042e-15*Log(J) - 2.6813110884009767e-14*t*Log(J) + &
   &  1.2879071621313094e-12*Log(c3)*Log(J) - &
   &  3.80352446061867e-15*t*Log(c3)*Log(J) - &
   &  1.8790172502456827e-14*Log(J)**2
    
      nacid=-4.7154180661803595 + 0.13436423483953885*t - & 
   &  0.00047184686478816176*t**2 - & 
   &  2.564010713640308*Log(c2) + 0.011353312899114723*t*Log(c2) + &
   &  0.0010801941974317014*Log(c2)**2 + 0.5171368624197119*Log(c3) - &
   &  0.0027882479896204665*t*Log(c3) + 0.8066971907026886*Log(c3)**2 - & 
   &  0.0031849094214409335*t*Log(c3)**2 - &
   &  0.09951184152927882*Log(c3)**3 + &
   &  0.00040072788891745513*t*Log(c3)**3 + 1.3276469271073974*Log(J) - &
   &  0.006167654171986281*t*Log(J) - 0.11061390967822708*Log(c3)*Log(J) + &
   &  0.0004367575329273496*t*Log(c3)*Log(J) + &
   &  0.000916366357266258*Log(J)**2
    
      namm=71.20073903979772 - 0.8409600103431923*t + &
   &  0.0024803006590334922*t**2 + &
   &  2.7798606841602607*Log(c2) - 0.01475023348171676*t*Log(c2) + &
   &  0.012264508212031405*Log(c2)**2 - 2.009926050440182*Log(c3) + &
   &  0.008689123511431527*t*Log(c3) - &
   &  0.009141180198955415*Log(c2)*Log(c3) + &
   &  0.1374122553905617*Log(c3)**2 - 0.0006253227821679215*t*Log(c3)**2 + &
   &  0.00009377332742098946*Log(c3)**3 + 0.5202974341687757*Log(J) - &
   &  0.002419872323052805*t*Log(J) + 0.07916392322884074*Log(c3)*Log(J) - &
   &  0.0003021586030317366*t*Log(c3)*Log(J) + &
   &   0.0046977006608603395*Log(J)**2
   
   else
   ! Nucleation rate less that 5E-6, setting j_log arbitrary small
      j_log=-300.
      J=exp(J_log)
   end if
   
   Return
  
  END SUBROUTINE merikanto2007


PURE  SUBROUTINE gnucl_kulmala(pso4g,   ptemp,  prhum,   &
                           pbnrate, oalpha, obeta, onstar )

    !  Authors:
    !  --------
    !  P. Stier, MPI-Met, Hamburg,    from the original f77 code
    !                                 in the routine gmxe_nuck        2001-2003
    !  J. Wilson, E. Vignati, JRC/EI, original source                 09/2000
    !
    !  Purpose:                                                           
    !  --------                                                           
    !  This routine calculates the instananeous nucleation rate            
    !  znucrate [molec. cm-3 s-1] from a given gas phase H2SO4 concentration      
    !  pso4g [molec. cm-3]. It also calculates the integrated change of 
    !  H2SO4 gas phase mass over one timestep due to nucleation 
    !  pa4delt(:,:,1) [molec. cm-3] as well as the number of nucleated 
    !  particles panew [1] during the timestep.
    !
    !  Interface:
    !  ----------
    !  *gnucl_kulmala* is called from *gmxe_nuck*
    !
    !  Method:
    !  -------
    !  Kulmala et al. (1998) 's formula for binary nucleation is 
    !  rewritten to take the form znucrate = exp[zalpha+ln(pso4g)*beta]. 
    !  Equation numbers are taken from Kulmala et al. (1998).
    !  After the calculation of the nucleation rate znucrate, it is 
    !  integrated in 2) analytically over one timestep, i.e.:
    !
    !  Integration of:
    ! 
    !  znucrate=d(critn*znav)/dt=exp[zalpha + zbeta*ln(znav)]
    !
    !  gives znav(t0+dt, znav(t0)) where znav(t0) is pso4g.
    !  znav is temporarily stored in zso4g_new and confined between
    !  0 and pso4g. 
    !  The number of nucleated particles is then calculated by 
    !  assuming a fixed critical mass critn of newly nucleated 
    !  particles and dividing the total mass of nucleated sulfate 
    !  by this critical mass of one nucleated particle. 
    !
    !  Externals:
    !  ----------
    !  None

    IMPLICIT NONE


!    INTEGER :: kproma, klev

    REAL(dp), INTENT(IN)  :: pso4g,        ptemp,    &
                             prhum
    REAL(dp), INTENT(OUT) :: pbnrate

    REAL(dp), INTENT(OUT), optional ::  oalpha, obeta, onstar
    REAL(dp) :: palpha, pbeta

    REAL(dp):: znwv,        zln_nac,      ztk,         zsupsat,  &
               zpeh2o,      zpeh2so4,     zra,         zxal,     &
               ztkn,        zssn,         zdelta

      !---1) Calculation of the nucleation rate: ----------------------------
      
      IF (pso4g .GT. 1e-5 .AND. prhum > ZERO) THEN  !!KP - Threshold from M7
   
       !IF (pso4g .GT. cmin_nuclmolec) THEN
       ! IF (pso4g .GT. cmin_nuclmolec .AND. prhum > ZERO) THEN
          ztk=ptemp
          zsupsat=prhum
          !
          !--- 1.1) Restrict t, and rh to limits where the parameterization ok:
          !
          ztkn = MAX(ztk, 220.0_dp)
          zssn = MIN(zsupsat, 0.90_dp)
   
          !
          !--- 1.2) Equlibrium vapour pressures (Jaeker-Mirabel (1995), JGR):
          !
          !--- H2O equlibrium vapour pressure (Tabata):
          !
          zpeh2o=0.750064*(10.**(8.42926609-1827.17843/          &
               ztkn-71208.271/ztkn/ztkn))*1333./bk/ztkn
          !
          !--- H2SO4 equlibrium vapour pressure at 360 
          !
          zpeh2so4=EXP(-10156./ztkn+16.259)*7.6e2*1333./bk/ztkn
          !
          !--- H2SO4 equlibrium vapour pressure - correction of ayers
          !    by kulmala - currently not used
          !
          !     payers=exp(-10156/360+16.259)*7.6e2
          !     zpeh2so4=exp(log(payers)+10156*(-1./ztkn+1./360.+0.38/(905-360) * &
          !              (1+log(360./ztkn)-360./ztkn)))*1333/bk/ztkn
          !
          !--- 1.3) Relative acidity (0.0 -1.0):
          ! 
          zra=pso4g/zpeh2so4
          !
          !--- 1.4) Water vapour molecule concentration [cm-3]:
          !
          znwv=zsupsat*zpeh2o
          !
          !--- 1.5) Factor delta in Eq. 22:
          ! 
          zdelta=1.0_dp+(ztkn-273.15)/273.15
          !
          !--- 1.6) Molefraction of H2SO4 in the critical cluster 
          !         minus the H2SO4(g) term in Eq. 17:
          !
          zxal = 1.2233-0.0154*zra/(zra+zssn)-0.0415*LOG(znwv)+ 0.0016*ztkn 
          !
          !--- 1.7) Exponent of the critical cluster (Eq. 18):
          !
          zln_nac = -14.5125+0.1335*ztkn-10.5462*zssn+1958.4*zssn/ztkn
          !
          !--- 1.8) Sum of all terms in Eq. 20 containing H2SO4(g):
          !
          pbeta = 25.1289 - 4890.8/ztkn + 7643.4*0.0102/ztkn - &
                         2.2479*zdelta*zssn - 1.9712*0.0102*zdelta/zssn
          ! 
          !--- 1.9) Sum all terms in Eq. 20 not containing H2SO4(g):
          !
          palpha = zln_nac*(-25.1289 + 4890.8/ztkn + 2.2479*zdelta*zssn) - &
                          1743.3/ztkn + zxal*(7643.4/ztkn - 1.9712*zdelta/zssn)
          !
          !--- 1.10) Nucleation rate [cm-3 s-1] (Kulmala et al., 1998):
          !
          pbnrate = EXP(palpha+LOG(pso4g)*pbeta)
   
       ELSE
   
          palpha =zero
          pbeta  =zero
          pbnrate=zero
   
       END IF ! pso4g > cmin_nuclmolec

       ! This allows to omit alpha and beta in cases where the rates are numerically integrated, #seb
       IF(present(oalpha)) oalpha = palpha
       IF(present(obeta))  obeta  = pbeta
       IF(present(onstar)) onstar = 100.0_dp 

  END SUBROUTINE gnucl_kulmala

!===========================================================================!

PURE  SUBROUTINE gnucl_vehkamaeki(ptemp,    prhd,  pmolecH2SO4, &  ! ECHAM5 temperature, relative humidity
                              pxtrnucr, pntot         )  ! nucleation rate, number of molecules in the
                                                         ! critical cluster
    !
    !   Authors:
    !   ---------
    !   C. TIMMRECK, MPI HAMBURG                                             2002
    !
    !   Purpose
    !   ---------
    !   Calculation of classical nucleation rate
    !               
    !   calculation of the nucleation rate after Vehkamaeki et al. (2002)
    !   The calculation of the nucrate ZKNH2SO4 is in cm^-3 s^-1
    !   and a coarse approxmation for the first class
    !
    !   Modifications:
    !   --------------
    !   R. Hommel; rewrite in f90, adopted to ECHAM5; MPI HAMBURG;      Dec. 2002
    !   P. Stier; modularisation and optimization;    MPI HAMBURG;      Jan  2003
    !
    !   H2SO4 still fixed to xxx molc/cm3, no sulfur cycle coupling yet
    !
    !   References:
    !   -----------
    !   Vehkamaeki et al. (2002), An improved parameterization for sulfuric
    !      acid/water nucleation rates for tropospheric and stratospheric
    !      conditions, J. Geophys. Res, 107, D22, 4622
    !
    !   Parameters
    !   ----------
    !   prho = prhop_neu in *sam*
    !
    !   prhd = relative humidity in sam_aeroprop & sam_nucl
    !
    !   pxtrnucr = nucleation rate in [1/m3s]
    !
    !----------------------------------------------------
    
  IMPLICIT NONE

    !----------------------------------------------------
    !
    
    REAL(dp), INTENT(IN) ::   ptemp,       &
                              prhd,        &
                              pmolecH2SO4

    REAL(dp), INTENT(OUT) ::  pxtrnucr,    &
                              pntot
    
    !----------------------------------------------------  
    ! Local Arrays
    

    REAL(dp):: zrhoa, zrh, zt, x, zjnuc, zrc, zxmole, zntot

    REAL(dp):: zrh_log, zrhoa_log
    REAL(dp):: zrh_log2, zrh_log3, zt2, zt3, zrhoa_log2, zrhoa_log3

    ! mz_ht_20081001+
    pntot   =zero
    pxtrnucr=zero   
    zjnuc = zero
    ! mz_ht_20081001-

    !----1.) Parameterization of  nucleation rate after Vehkamaeki et al. (2002)
  
          ! t: temperature in K (190.15-300.15K)                                  
          ! zrh: saturatio ratio of water (0.0001-1)                               
          ! zrhoa: sulfuric acid concentration in 1/cm3 (10^4-10^11 1/cm3)         
          ! jnuc: nucleation rate in 1/cm3s (10^-7-10^10 1/cm3s)                  
          ! ntot: total number of molecules in the critical cluster (ntot>4)      
          ! x: molefraction of H2SO4 in the critical cluster                      
          ! rc: radius of the critical cluster in nm                              
  
          ! Calculate nucleation only for valid thermodynamic conditions:
  
  !       IF( (pmolecH2SO4 > cmin_nuclmolec) .AND. & !!KP Check threshold
          IF( (pmolecH2SO4 > 1.0d4)            .AND. &  !! se_mpic_20161004 changed threshold back as it is dependent on param.
              (prhd >=1.0d-4)                            .AND. &
              (ptemp>=190.15 .OR. ptemp<=300.15_dp)    ) THEN
  
             zrhoa=MIN(pmolecH2SO4,1.e11_dp)
             zrh=MIN(prhd,1.0_dp)
             zt=ptemp
             zt2=zt*zt
             zt3=zt2*zt
             zrh_log=LOG(zrh)
             zrh_log2=zrh_log*zrh_log
             zrh_log3=zrh_log2*zrh_log
             zrhoa_log=LOG(zrhoa)
             zrhoa_log2=zrhoa_log*zrhoa_log
             zrhoa_log3=zrhoa_log2*zrhoa_log
  
  
             ! Equation (11) - molefraction of H2SO4 in the critical cluster
  
          x=0.7409967177282139_dp - 0.002663785665140117*zt   &
            + 0.002010478847383187*LOG(zrh)    &
            - 0.0001832894131464668*zt*LOG(zrh)    &
            + 0.001574072538464286*LOG(zrh)**2        &
            - 0.00001790589121766952*zt*LOG(zrh)**2    &
            + 0.0001844027436573778*LOG(zrh)**3     &
            -  1.503452308794887e-6*zt*LOG(zrh)**3    &
            - 0.003499978417957668*LOG(zrhoa)   &
            + 0.0000504021689382576*zt*LOG(zrhoa)
  
          zxmole=x !qqq
  
          ! Equation (12) - nucleation rate in 1/cm3s
  
          zjnuc =0.1430901615568665_dp + 2.219563673425199*zt -   &
                0.02739106114964264*zt**2 +     &
                0.00007228107239317088*zt**3 + 5.91822263375044/x +     &
                0.1174886643003278*LOG(zrh) + 0.4625315047693772*zt*LOG(zrh) -   &
                0.01180591129059253*zt**2*LOG(zrh) +     &
                0.0000404196487152575*zt**3*LOG(zrh) +    &
                (15.79628615047088*LOG(zrh))/x -     &
                0.215553951893509*LOG(zrh)**2 -    &
                0.0810269192332194*zt*LOG(zrh)**2 +     &
                0.001435808434184642*zt**2*LOG(zrh)**2 -    &
                4.775796947178588e-6*zt**3*LOG(zrh)**2 -     &
                (2.912974063702185*LOG(zrh)**2)/x -   &
                3.588557942822751*LOG(zrh)**3 +     &
                0.04950795302831703*zt*LOG(zrh)**3 -     &
                0.0002138195118737068*zt**2*LOG(zrh)**3 +    &
                3.108005107949533e-7*zt**3*LOG(zrh)**3 -     &
                (0.02933332747098296*LOG(zrh)**3)/x +     &
                1.145983818561277*LOG(zrhoa) -    &
                0.6007956227856778*zt*LOG(zrhoa) +    &
                0.00864244733283759*zt**2*LOG(zrhoa) -    &
                0.00002289467254710888*zt**3*LOG(zrhoa)! -    &
  
  
          zjnuc =zjnuc - &
               (8.44984513869014*LOG(zrhoa))/x +    &
                2.158548369286559*LOG(zrh)*LOG(zrhoa) +   &
                0.0808121412840917*zt*LOG(zrh)*LOG(zrhoa) -    &
                0.0004073815255395214*zt**2*LOG(zrh)*LOG(zrhoa) - &
                4.019572560156515e-7*zt**3*LOG(zrh)*LOG(zrhoa) +    &
                (0.7213255852557236*LOG(zrh)*LOG(zrhoa))/x +    &
                1.62409850488771*LOG(zrh)**2*LOG(zrhoa) -    &
                0.01601062035325362*zt*LOG(zrh)**2*LOG(zrhoa) +   &
                0.00003771238979714162*zt**2*LOG(zrh)**2*LOG(zrhoa) +    &
                3.217942606371182e-8*zt**3*LOG(zrh)**2*LOG(zrhoa) -    &
                (0.01132550810022116*LOG(zrh)**2*LOG(zrhoa))/x +    &
                9.71681713056504*LOG(zrhoa)**2 -    &
                0.1150478558347306*zt*LOG(zrhoa)**2 +    &
                0.0001570982486038294*zt**2*LOG(zrhoa)**2 +    &
                4.009144680125015e-7*zt**3*LOG(zrhoa)**2 +    &
                (0.7118597859976135*LOG(zrhoa)**2)/x -    &
                1.056105824379897*LOG(zrh)*LOG(zrhoa)**2 +    &
                0.00903377584628419*zt*LOG(zrh)*LOG(zrhoa)**2 -    &
                0.00001984167387090606*zt**2*LOG(zrh)*LOG(zrhoa)**2 +    &
                2.460478196482179e-8*zt**3*LOG(zrh)*LOG(zrhoa)**2 -    &
                (0.05790872906645181*LOG(zrh)*LOG(zrhoa)**2)/x -    &
                0.1487119673397459*LOG(zrhoa)**3 +    &
                0.002835082097822667*zt*LOG(zrhoa)**3 -    &
                9.24618825471694e-6*zt**2*LOG(zrhoa)**3 +    &
                5.004267665960894e-9*zt**3*LOG(zrhoa)**3 -    &
                (0.01270805101481648*LOG(zrhoa)**3)/x
  
          zjnuc=EXP(zjnuc)      !   add. Eq. (12) [1/(cm^3s)]
  
  
          ! Equation (13) - total number of molecules in the critical cluster
  
          zntot =-0.002954125078716302_dp - 0.0976834264241286*zt +   &
                 0.001024847927067835*zt**2 - 2.186459697726116e-6*zt**3 -    &
                 0.1017165718716887/x - 0.002050640345231486*LOG(zrh) -   &
                 0.007585041382707174*zt*LOG(zrh) +    &
                 0.0001926539658089536*zt**2*LOG(zrh) -   &
                 6.70429719683894e-7*zt**3*LOG(zrh) -    &
                 (0.2557744774673163*LOG(zrh))/x +   &
                 0.003223076552477191*LOG(zrh)**2 +   &
                 0.000852636632240633*zt*LOG(zrh)**2 -    &
                 0.00001547571354871789*zt**2*LOG(zrh)**2 +   &
                 5.666608424980593e-8*zt**3*LOG(zrh)**2 +    &
                 (0.03384437400744206*LOG(zrh)**2)/x +   &
                 0.04743226764572505*LOG(zrh)**3 -    &
                 0.0006251042204583412*zt*LOG(zrh)**3 +   &
                 2.650663328519478e-6*zt**2*LOG(zrh)**3 -    &
                 3.674710848763778e-9*zt**3*LOG(zrh)**3 -   &
                 (0.0002672510825259393*LOG(zrh)**3)/x -    &
                 0.01252108546759328*LOG(zrhoa) !+   &
  
          zntot =zntot + &
                 0.005806550506277202*zt*LOG(zrhoa) -    &
                 0.0001016735312443444*zt**2*LOG(zrhoa) +   &
                 2.881946187214505e-7*zt**3*LOG(zrhoa) +    &
                 (0.0942243379396279*LOG(zrhoa))/x -   &
                 0.0385459592773097*LOG(zrh)*LOG(zrhoa) -   &
                 0.0006723156277391984*zt*LOG(zrh)*LOG(zrhoa) + &
                 2.602884877659698e-6*zt**2*LOG(zrh)*LOG(zrhoa) +    &
                 1.194163699688297e-8*zt**3*LOG(zrh)*LOG(zrhoa) -   &
                 (0.00851515345806281*LOG(zrh)*LOG(zrhoa))/x -    &
                 0.01837488495738111*LOG(zrh)**2*LOG(zrhoa) +   &
                 0.0001720723574407498*zt*LOG(zrh)**2*LOG(zrhoa) -   &
                 3.717657974086814e-7*zt**2*LOG(zrh)**2*LOG(zrhoa) -    &
                 5.148746022615196e-10*zt**3*LOG(zrh)**2*LOG(zrhoa) +    &
                 (0.0002686602132926594*LOG(zrh)**2*LOG(zrhoa))/x -   &
                 0.06199739728812199*LOG(zrhoa)**2 +    &
                 0.000906958053583576*zt*LOG(zrhoa)**2 -   &
                 9.11727926129757e-7*zt**2*LOG(zrhoa)**2 -    &
                 5.367963396508457e-9*zt**3*LOG(zrhoa)**2 -   &
                 (0.007742343393937707*LOG(zrhoa)**2)/x +    &
                 0.0121827103101659*LOG(zrh)*LOG(zrhoa)**2 -   &
                 0.0001066499571188091*zt*LOG(zrh)*LOG(zrhoa)**2 +    &
                 2.534598655067518e-7*zt**2*LOG(zrh)*LOG(zrhoa)**2 -    &
                 3.635186504599571e-10*zt**3*LOG(zrh)*LOG(zrhoa)**2 +    &
                 (0.0006100650851863252*LOG(zrh)*LOG(zrhoa)**2)/x +   &
                 0.0003201836700403512*LOG(zrhoa)**3 -    &
                 0.0000174761713262546*zt*LOG(zrhoa)**3 +   &
                 6.065037668052182e-8*zt**2*LOG(zrhoa)**3 -    &
                 1.421771723004557e-11*zt**3*LOG(zrhoa)**3 +   &
                 (0.0001357509859501723*LOG(zrhoa)**3)/x
  
                 
             pntot=EXP(zntot)  !  add. Eq. (13)
  
            
             ! Equation (14) - radius of the critical cluster in nm
  
             zrc=EXP(-1.6524245+0.42316402*x+0.33466487*LOG(pntot))    ! [nm]
  
             ! Conversion [nm -> m]
  
  !ueberfl.  zrxc(jl)=zrc*1e-9
  
             !----1.2) Limiter
  
             IF (pntot < 4.0 ) THEN
                IF ( zt < 195.15) THEN
                   !!! zjnuc=1.e5
                   zjnuc=zero !!!KP testing - same as m7
                END IF
             END IF
  
             IF(zjnuc < 1.e-7) THEN
                zjnuc=zero
             END IF
  
             ! limitation to 1E+10 [1/cm3s]
             zjnuc=MIN(zjnuc,1.e10_dp)          !mz_kp_20080112 cmin_aernl is far too low a threshold
             !!!!zjnuc=MIN(zjnuc,cmin_aernl)
  
             pxtrnucr = zjnuc
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!KP Insered from M7  - for testing only
  !       ! mz_ak_20041102+
  !       ! convert total number of molecules in the critical cluster
  !       ! to number of sulfate molecules:
  !
  !       pntot=zntot*zxmole
  !
  !        !print*,'pntot=',pntot,zntot,zxmole
  !        ! mz_ak_20041102-
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
          ELSE ! pmolecH2SO4, ptemp , prhd out of range
  
             pntot   =zero
             pxtrnucr=zero
  
          END IF

END SUBROUTINE gnucl_vehkamaeki


pure function update_gas(vapour, nucr, pncrit, ztmst)
! this will update the gas phase concentration of condensable vapours by subtracting
! the amount that went into the particle phase due to nucleation. 
! This is a code rearrangement to allow easier maintenance and extension ot more condensable vapours
! Code originally from gmxe_nuck.
! vapour = vapour concentration to  be updated [cm-3]
! nucr   = nucleation rate [cm-3 s-1]
! pncrit = number of molecules of vapour in a critical cluster or the aerosol diameter of particle formation [#molecules]
! ztmst  = time step [s] 
implicit none
real(kind=dp), intent(in) :: vapour, nucr, ztmst, pncrit
real(kind=dp) :: update_gas ! updated gas phase concentration in [#molecules cm-3]

           IF(vapour > cmin_epsilon .AND. &
              nucr > cmin_epsilon) THEN
             update_gas = vapour - (nucr * pncrit * ztmst)  
           ELSE
!             zso4g_new(jl,jk)=MAX(pso4g(jl,jk),0.0_dp)
             update_gas = vapour
           ENDIF

end function update_gas


pure subroutine driver_vehkamaeki( ptemp, prhum, pso4g, &  ! ECHAM5 temperature, relative humidity, 
                              panew, pa4delt, &    ! new formed particles, vapours added to nucleation mode
                              nucr, ztmst) ! nucleation rate, number of molecules in the
                              
  ! Call all routines to calculate new particles due to pure binary nucleation according to the Vehkamaeki parameterisation. 
  ! It calculates also the depletion of the vapour concentration. 
  ! Author: Sebastian Ehrhart 09/2016
 
  implicit none
!  integer, intent(in) :: kproma, klev, nucm_so4m, nucm_Hp, naertot
  real(kind=dp), intent(in) :: ztmst
  real(kind=dp), intent(in), dimension(:,:) :: ptemp, prhum
  ! local
  real(kind=dp), dimension(ubound(ptemp,1),ubound(ptemp,2)) :: new_pso4g ! new gas phase [H2SO4] (#molecules cm-3)
  real(kind=dp) :: pncrit ! number of molecules in critical cluster
  real(kind=dp) :: dgas ! pso4g(jl,jk)-new_pso4g(jl,jk)
  integer :: jk, jl
  ! inout
  real(kind=dp), dimension(:,:), intent(inout) :: pso4g
  ! output
  real(kind=dp), dimension(:,:), intent(out) :: panew, pa4delt, nucr

  panew   = zero
  pa4delt = zero
  pncrit  = zero
  nucr    = zero

  do jk=1, ubound(ptemp,2)
    do jl=1, ubound(ptemp,1)
      ! 1. nucleation rate
      call gnucl_vehkamaeki(ptemp(jl,jk), prhum(jl,jk), pso4g(jl,jk), & ! ECHAM5 temperature, relative humidity
                        nucr(jl,jk), pncrit)

  ! now update for the lost H2SO4 molecules in the newly formed particles
  ! would be better to store the lower part in a subroutine. But Kulmala parameterisation 
  ! does not need the update_gas step. #seb
      if( nucr(jl,jk) > zeps  ) then  

      new_pso4g(jl,jk) = update_gas(pso4g(jl,jk), nucr(jl,jk), pncrit, ztmst)

  !--- 2.1) Security check:
      new_pso4g(jl,jk) = MAX(new_pso4g(jl,jk), zero) ! this protects against negative updated gas phase concentrations, e.g.
      ! J* ncrit * dt > H2SO4_start. It would nucleate all vapour, could happen when dt is set too large. In a multicomponent system we would need to  check all. 
      new_pso4g(jl,jk) = MIN(new_pso4g(jl,jk), pso4g(jl,jk)) ! this check protects only against negative nucleation rates 

  !--- 2.2) Calculate mass of nucleated H2SO4 (equals the net 
  !         gas phase H2SO4 loss):
  ! 
!      dgas = pso4g(jl,jk)-new_pso4g(jl,jk) ! saves one calculation and allows replacing one multiplication with addition mpic_se_20161004
      pa4delt(jl,jk) = MAX(pso4g(jl,jk)-new_pso4g(jl,jk),0.0_dp)  !mz_kp_20080112, modified mpic_se_20161004

  !
  !--- 2.3) Calculate the number of nucleated particles (nucleated mass 
  !         divided by the assumed mass of a critical cluster critn):
      if(pncrit > 0.0_dp) panew(jl,jk)=pa4delt(jl,jk)/pncrit 
           
  !--- 2.4) Calculate changes in gas phase H2SO4 due to nucleation:
!           print*, pa4delt(jl,jk,nso42m(NS)), pso4g(jl,jk), jl,jk
      pso4g(jl,jk)=pso4g(jl,jk)-pa4delt(jl,jk)

      end if
    end do
  end do

end subroutine driver_vehkamaeki


subroutine driver_kulmala(kproma, klev,        &  ! ECHAM5 dimensions
                         ptemp, prhum, pso4g, &  ! ECHAM5 temperature, relative humidity, 
                         panew, pa4delt, naertot, &  ! new formed particles, vapours added to nucleation mode
                         nucr, pncrit, ztmst, &  ! nucleation rate, number of molecules in the
                         nucm_so4m, nucm_Hp)! ! index {SO4-, H+} in pa4delt

  implicit none

  integer, intent(in) :: kproma, klev, naertot, nucm_so4m, nucm_Hp
!  real(kind=dp), intent(in) :: 
  real(kind=dp), intent(in), dimension(kproma, klev) :: ptemp, prhum
  real(kind=dp), intent(in) :: ztmst
  ! local
  real(kind=dp) :: new_pso4g ! new gas phase [H2SO4] (#molecules cm-3)
  real(kind=dp) :: zalpha, zbeta, zf1 ! used in integration of nucleated H2SO4_g
  real(kind=dp) :: dgas ! pso4g(jl,jk)-new_pso4g(jl,jk)
  integer :: jk, jl
  ! inout
  real(kind=dp), intent(inout) :: pso4g(kproma,klev)
  ! output
  real(kind=dp), intent(out) :: panew(kproma,klev), pa4delt(kproma,klev,0:naertot), & 
                                pncrit(kproma, klev), nucr(kproma, klev)

!  pncrit(:,:)=100.0_dp
!CDIR NOIEXPAND

     DO jk=1, klev
        DO jl=1, kproma

          CALL gnucl_kulmala(pso4g(jl,jk), ptemp(jl,jk), prhum(jl,jk), &
                             nucr(jl,jk), oalpha=zalpha, obeta=zbeta, onstar=pncrit(jl,jk))

     !--- 2) Analytical integration of the nucleation rate (Eq. 19) ----------------
     !       over timestep ztmst assuming no new H2SO4(g) production:
     !
     !       d(N_av/critn)/dt=exp(alpha + beta*ln(N_av)) => ... =>
     !           
     !       N_av(t0+dt)=
     !       [N_av(t0)**(1-beta) + critn*(1-beta)exp(alpha)*dt]**(1/(1-beta)
     !

           IF (nucr(jl,jk) .GT. zeps) THEN ! harmonised with vahkamaeki subroutine mpic_se_20161004
              zf1 = pso4g(jl,jk)**(1.0_dp-zbeta)-&
                          pncrit(jl,jk)*EXP(zalpha)*(1.0_dp-zbeta)*ztmst
              new_pso4g = EXP(LOG(zf1)/(1.0_dp - zbeta)) ! this will never be negative, 
              ! so we don't need the Vehkamaeki protection against negative gas phase concentrations. mpic_se_20161004
       !--- 2.1 safety check
              new_pso4g = MIN(new_pso4g, pso4g(jl,jk)) 

       !--- 2.2) Calculate mass of nucleated H2SO4 (equals the net 
       !         gas phase H2SO4 loss):
       ! 
           dgas = pso4g(jl,jk)-new_pso4g ! saves one calculation and allows replacing one multiplication with addition
           pa4delt(jl,jk,nucm_so4m) = MAX(dgas,0.0_dp)  !mz_kp_20080112, modified mpic_se_20161004
           ! adding also H+ mass to compensate for H2SO4 = 2 H+ + SO4--
           pa4delt(jl,jk,nucm_Hp)    = MAX((dgas + dgas),0.0_dp)  ! mz_ht_20090310, modified mpic_se_20161004
       !
       !--- 2.3) Calculate the number of nucleated particles (nucleated mass 
       !         divided by the assumed mass of a critical cluster critn):
           panew(jl,jk)=pa4delt(jl,jk,nucm_so4m)/pncrit(jl,jk) 
                
       !--- 2.4) Calculate changes in gas phase H2SO4 due to nucleation:
     !           print*, pa4delt(jl,jk,nso42m(NS)), pso4g(jl,jk), jl,jk
           pso4g(jl,jk)=pso4g(jl,jk)-pa4delt(jl,jk,nucm_so4m)

           ELSE
              new_pso4g = pso4g(jl,jk)
              panew(jl,jk) = zero
           END IF

        END DO
     END DO

end subroutine driver_kulmala



pure function nucleation_tendencies(lnchan, nr, nstar)
  ! calculate species tendencies due to nucleation
  implicit none
  ! input
  logical, intent(in) :: lnchan(:,:)
  real(kind=dp), intent(in) :: nr(:), nstar(:,:)
  ! local
  integer :: i, k
  ! return
  real(kind=dp), dimension(ubound(nstar,2)) :: nucleation_tendencies

  nucleation_tendencies = 0.0_dp

  do i=1, ubound(lnchan,dim=2)
    do k=1, ubound(lnchan,dim=1)
      if(lnchan(k,i)) nucleation_tendencies(i) & 
        & = nucleation_tendencies(i) - nr(k)*nstar(k,i)
    end do
  end do

end function nucleation_tendencies


pure subroutine nucleation_array(chid,temp,rh,vapour,qion,krec,kl,nr,nstar) !, iolo) 
  ! Calculate nucleation rate and number of molecules of each species for  
  implicit none
  ! input
  integer, intent(in) :: chid(:) ! id number of nucleation channel, some channels might be merged into a mechanism
  real(kind=dp), intent(in) :: vapour(:) ! vapour concentration (#molecules cm-3)
  real(kind=dp), intent(in) :: temp, rh ! temperature (K), relative humidity (%)
  real(kind=dp), intent(in) :: qion     ! ion pair production rate (i.p. cm-3 s-1)
  real(kind=dp), intent(in) :: krec     ! ion-ion recombination rate constant (cm3 s-1)
  real(kind=dp), intent(in) :: kl       ! first order loss rate of small pre-nucleation ions, e.g. to condensational sink
  ! local
  integer :: i
  real(kind=dp) :: kitot(2) ! total first order ion loss rate constant (s-1), index 1 negative ions, index 2 positive ions
  real(kind=dp) :: np, nn ! ion concentrations
  real(kind=dp) :: frac ! fraction
  real(kind=dp) :: hom
  real(kind=dp) :: RHsc, Tsc
  ! return
  real(kind=dp), intent(out) :: nr(:), nstar(:,:) ! nucleation rate in cm-3 s-1 and number of molecules in a critical cluster for each species and channel
!IONBALANCE real(kind=dp), intent(out) :: iolo(2) ! ion loss due to nucleation,  index 1 negative ions, index 2 positive ions

  kitot = kl
  nr = 0.0_dp
  nstar = 0.0_dp
!  iolo = 0.0_dp

  if(l_org_Tdep) Tsc = exp(BtdOrg*temp)
  if(l_Dunne_RH) RHsc = 1.0_dp +  frhc1*(rh-0.38_dp) + frhc2*(rh-0.38 )**3 *(temp-208.0_dp)**2

!  write(*,*) "aa"
  if(inoxo1 > 0 .and. inoxo2 > 0) hom = vapour(inoxo1) + vapour(inoxo2)
  do i=1, size(chid)
        ! call the nucleation channels if they were selected
        ! check if select case or if and integer comparison is better
    select case(chid(i))
    case(idkul)
       call gnucl_kulmala(vapour(insa), temp, rh, nr(i), onstar=nstar(i,insa))
    case(idvehk)
      call gnucl_vehkamaeki(temp, rh, vapour(insa), nr(i), nstar(i,insa) )
    case(idmeri) 
      ! important: merikanto needs NH3 as mixing ratio
!      call merikanto2007(temp, rh, vapour(insa), vapour(innh3), nr(i), nstar(i,insa), nstar(i,innh3))
    case(iddunbn) ! binary neutral
      nr(i) = Dunne2016binary(vapour(insa)/1.0e6_dp,temp)
      if(l_Dunne_RH) nr(i) =  RHsc*nr(i)
      nstar(i,insa) = 28.0_dp ! based on dry density
    case(iddunbi) ! ternary neutral
      nr(i) = Dunne2016binaryIon(vapour(insa)/1.0e6_dp, temp)!*vapour(innion)
      if(l_Dunne_RH) nr(i) =  RHsc*nr(i)
      nstar(i,insa) = 28.0_dp !  based on dry density
      kitot(1) = kitot(1) + nr(i)
    case(idduntn) ! binary iin
      nr(i) = Dunne2016ternary(vapour(insa)/1.0e6_dp, vapour(innh3)/1.0e6_dp, temp)
      if(l_Dunne_RH) nr(i) =  RHsc*nr(i)
      nstar(i,insa)  = 25.0_dp ! based on dry density
      nstar(i,innh3) = 25.0_dp ! based ond ry density
    case(iddunti) ! ternary iin
      nr(i) = Dunne2016ternaryIon(vapour(insa)/1.0e6_dp, vapour(innh3)/1.0e6_dp, temp)!*vapour(innion)
      if(l_Dunne_RH) nr(i) =  RHsc*nr(i)
      nstar(i,insa)  = 25.0_dp ! based on dry density
      nstar(i,innh3) = 25.0_dp ! based on dry density
      kitot(1) = kitot(1) + nr(i)
    case(idric) ! Riccobono et al 2014, neutral only
      nr(i) = riccobono2014b(vapour(insa), vapour(inoxo1))
      if(l_org_Tdep) nr(i) = Tsc*nr(i)
      nstar(i,insa)   = 8.0_dp ! esitmate
      nstar(i,inoxo1) = 8.0_dp ! estimate
    case(idkirkn) ! Kirkby et al 2016, neutral
      if(hom > zeps) then
        nr(i) = kirkby2016(hom/1.0e7_dp, kirkby_a(1), kirkby_a(2), kirkby_a(5))
        if(l_org_Tdep) nr(i) = Tsc*nr(i)
        ! non ideal solution the kirkby parameterisation can't be mathematically splitted into OH and O3 product contributions
        ! however the amount taken up by nucleation itself is probably small, only a few molecules per nucleated particle,
        ! since the monomers are relatively heavy (~250 g/mol and more) <
        frac = vapour(inoxo1)/hom
        nstar(i, inoxo1) = 8.0_dp*frac
        nstar(i, inoxo2) = 8.0_dp*(1-frac)
      else
        nr(i) = 0.0d0
        nstar(i, inoxo1) = 0.0_dp
        nstar(i, inoxo2) = 0.0_dp
      end if
      ! >
      ! nstar(i,inoxo2) = 4.0_dp ! guess, though justified by the exponent of organic dependence in the power law ! <<< origianl DELETE if accepted
    case(idkirki) ! Kirkby et al 2016 charged
      if(hom > zeps) then
        nr(i) = kirkby2016(hom/1.0e7_dp, kirkby_a(3), kirkby_a(4), kirkby_a(5))
        if(l_org_Tdep) nr(i) = Tsc*nr(i)
        ! non ideal solution the kirkby parameterisation can't be mathematically splitted into OH and O3 product contributions
        ! however the amount taken up by nucleation itself is probably small, only a few molecules per nucleated particle,
        ! since the monomers are relatively heavy (~250 g/mol and more) <
        frac = vapour(inoxo1)/hom
        nstar(i, inoxo1) = 8.0_dp*frac
        nstar(i, inoxo2) = 8.0_dp*(1-frac) 
      else
        nr(i) = 0.0d0
        nstar(i, inoxo1) = 0.0_dp
        nstar(i, inoxo2) = 0.0_dp
      end if
      ! >
!      nstar(i,inoxo2) = 4.0_dp ! original delete if ok
      kitot(1) = kitot(1) + nr(i)
      kitot(2) = kitot(2) + nr(i)
    case(idalmei)
      nr(i) = almeida2014(vapour(indma), vapour(insa))
      nstar(i,indma) = 20.0_dp !4.0_dp ! guess
      nstar(i,insa)  = 20.0_dp !4.0_dp ! guess
    end select
  end do

  ! loop through the ion induced chanels to include the ion contribution and limit to ion pair production rate
  ! for performance improvements this should be moved to the ion submodel. However, it can't be moved as
  nn=0.0_dp; np=0.0_dp
  nn = ( sqrt(kitot(1)**2 + 4.0_dp*krec*qion) - kitot(1) ) / (2.0_dp*krec)
  np = ( sqrt(kitot(2)**2 + 4.0_dp*krec*qion) - kitot(2) ) / (2.0_dp*krec)
  do i =1, size(iinchind)
  ! now the channels with positive and negative ion induced nucleation
  ! J = Jiin * (negative + positive) negative ions treated above
    if(chid(iinchind(i)) == idkirki) then
       nr(iinchind(i)) = nr(iinchind(i)) * (nn+np)
!IONBALANCE      iolo(1) = iolo(1) + nr(iinchind(i))
!IONBALANCE      iolo(2) = iolo(2) + nr(iinchind(i))
    else
       nr(iinchind(i)) = nr(iinchind(i)) * nn
!IONBALANCE       iolo(1) = iolo(1) + nr(iinchind(i))
    end if
  end do

end subroutine nucleation_array



pure subroutine nucleation_driver(vapour, temp, rh, tstep, qion, krec, kl, nr, particles_new, totpfr, condensed)
! drive nucleation
  implicit none
  real(kind=dp), intent(inout) :: vapour(:) ! dimension number of species
  real(kind=dp), intent(in) :: temp, rh ! temperature (K), relative humidity (%)
  real(kind=dp), intent(in) :: tstep ! sum(dt), the overall timestep of the calling modell (s)
  real(kind=dp), intent(in) :: qion     ! ion pair production rate (i.p. cm-3 s-1)
  real(kind=dp), intent(in) :: krec     ! ion-ion recombination rate constant (cm3 s-1)
  real(kind=dp), intent(in) :: kl       ! first order loss rate of small pre-nucleation ions, e.g. to condensational sink
  ! return
  real(kind=dp), intent(out) :: nr(:), &  ! on return the average new particle formation rate for each channel (# cm-3 s-1)
                                particles_new, & ! new particle concentrations (# cm-3)
                                totpfr, & ! total average particle formation rate ( #particle cm-3 s-1)
                                condensed(:) ! condensed vapour molecules per gas phase volume, dimensions like vapours
  ! local
  integer :: ntstep, nmax(1)
  integer, parameter :: maxstep=100
!IONBALANCE  real(kind=dp) :: tmpiolo(2), iolo(2)
  real(kind=dp) :: vapour_new(size(vapour))
  real(kind=dp) :: dt, trem
  real(kind=dp) :: nstar(size(nr), size(vapour)) ! number of particles in newly formed particles
  real(kind=dp) :: tmpnr(size(nr)) ! temporary storage to calculate the average
  real(kind=dp) :: td(size(vapour)) ! store tendencies
  real(kind=dp) :: tddt(size(vapour)) ! tendencies*dt
  logical :: lt(size(vapour))

  trem = tstep
  dt=tstep !/ntstep
  tmpnr = 0.0_dp
  ntstep = 0
  particles_new = 0.0_dp
  condensed = 0.0_dp
  vapour_new = vapour

  do while (trem > 0.0_dp)
    !---- create the vector of nucleation rates and assigns n*
    call nucleation_array(chid, temp, rh, vapour_new, qion, krec, kl, nr, nstar)
    !---- calculate depletion of vapours 
    td = nucleation_tendencies(lnchan, nr, nstar)! 
    lt = (vapour_new + td*dt) < 0.0d0
    if(any(lt)) then
      dt = minval(vapour_new/abs(td), mask=lt)
      dt = min(trem, dt)
    end if
!    nmax = minloc(vapour_new/abs(td))
    tddt = td*dt
    vapour_new = max(vapour_new + (tddt),0.0_dp) !nucleation_tendencies(lnchan, nr, nstar)*dt
    condensed = condensed - tddt ! condensed species 
    ! calculate number of formed particles
    particles_new = particles_new + sum(nr)*dt
!IONBALANCE    tmpiolo = tmpiolo + iolo
    tmpnr = tmpnr + nr
    trem = trem - dt
    ntstep = ntstep + 1
    if(.not. trem > 0.0_dp) exit
    dt = trem 
  end do
  ! now the average particle formation rate
  nr = tmpnr / ntstep
!IONBALANCE  totiolo = tmpiolo / nstep
  totpfr = particles_new / tstep
  vapour = vapour_new
!  condensed = 0.0_dp  !DEBUG
!  particles_new = 0.0_dp ! 1.0e4_dp !DEBUG

end subroutine nucleation_driver



!TBD subroutine gmxe_assign_index(jm)
!TBD   use messy_gmxe_mem, only : species, spec_number 
!TBD   ! store index of components in nucleation mode 
!TBD   integer, intent(in) :: jm ! nucleation mode index, should be 1 always
!TBD   ! local
!TBD   integer :: jc, jt
!TBD 
!TBD   do jc=1, spec_number
!TBD     jt = species(jc)%aermlidx(jm)
!TBD     if( insa   > 0 .and. TRIM(species(jc)%name) == vapour_names(insa) )  ip4d(insa)   = jt
!TBD     if( innh3  > 0 .and. TRIM(species(jc)%name) == vapour_names(innh3))  ip4d(innh3)  = jt
!TBD     if( indma  > 0 .and. TRIM(species(jc)%name) == vapour_names(indma))  ip4d(indma)  = jt
!TBD     if( inoxo1 > 0 .and. TRIM(species(jc)%name) == vapour_names(inoxo1)) ip4d(inoxo1) = jt
!TBD     if( inoxo2 > 0 .and. TRIM(species(jc)%name) == vapour_names(inoxo2)) ip4d(inoxo2) = jt
!TBD   end do
!TBD 
!TBD end subroutine gmxe_assign_index


! Nieminen ACP 2010 DOI: 10.5194/acp-10-9773-2010
! Growth Rate by sulphuric acid condensation
real(kind=dp) elemental function GR_Nieminen2010(sa,rh) result(y) ! growth rate in nm/h, wet radius
  implicit none
  real(kind=dp), intent(in) :: sa ! number concentration h2so4 (cm-3)
  real(kind=dp), intent(in) :: rh ! RH in %
  
  y = sa / (661.1_dp*rh*rh - 1.129e5_dp*rh + 1.549e7_dp)

end function GR_Nieminen2010


! Organic and H2SO4 growth Gordon et al
real(kind=dp) elemental function Gordon2016(sa, hom) result(y) ! growth rate in nm/h
  implicit none
  real(kind=dp), intent(in) :: sa ! number concentration h2so4 (cm-3)
  real(kind=dp), intent(in) :: hom ! number concentration oxidised organics (cm-3)
  y = 7.8e-8_dp*sa + 1.41e-7*hom

end function Gordon2016


! Anttila et al J. Aer. Scie 41 2010 621-636 http://dx.doi.org/10.1016/j.jaerosci.2010.04.008
! could be made elemental, with some modifications
real(kind=dp) elemental function Anttila2010(jnuc, coags, numnuc, gr, mnuc, d1, d2, temp) result(y)

  implicit none
  real(kind=dp), intent(in) :: jnuc ! nucleation rate (cm-3 s-1)
  real(kind=dp), intent(in) :: coags ! coagulation sink, pre-existing particles (s-1)
  real(kind=dp), intent(in) :: numnuc ! concentration of nucleation mode clusters (cm-3)
  real(kind=dp), intent(in) :: gr   ! growth rate (nm h-1)
  real(kind=dp), intent(in) :: mnuc ! molar weight of nucleating cluster (g/mol)
  real(kind=dp), intent(in) :: d1   ! initial size of nucleated clusters (nm)
  real(kind=dp), intent(in) :: d2   ! final size of nucleated cluster (nm)
  real(kind=dp), intent(in) :: temp ! temperature (K)
  
  ! 5. calculate Jrev
  y = jnuc *exp( -1.0_dp *d1  & 
                & *gam(d1,d2) &
                & *coagtot(numnuc, coags) &
                & /grtot(gr, d1, numnuc, mnuc, temp) )

  contains

  real(kind=dp) elemental function coagtot(numnuc, coags)
    ! total coagulation sink (pre-existing + new)
    implicit none
    real(kind=dp), intent(in) :: numnuc, coags
    real(kind=dp), parameter :: keff = 5e-10_dp

    coagtot = coags + keff*numnuc

  end function coagtot


  real(kind=dp) elemental function grtot(gr, d1, numnuc, mnuc, temp)
    use messy_main_constants_mem, only : Navo => N_A
    implicit none
    real(kind=dp), intent(in) :: gr, d1, numnuc, mnuc, temp
!TBD    real(kind=dp), parameter :: scpre = 1.57e-6_dp ! Growth rate prefactor self coagulation
!TBD    real(kind=dp), parameter :: lambda = 6.0_dp
    real(kind=dp), parameter :: scla = 9.42e-6_dp  ! lambda*1.57e-6_dp prefactor from Anttila

    grtot = gr + scla *d1**3.0_dp *cnuc(mnuc,temp) *numnuc *Navo

  end function grtot


  real(kind=dp) elemental function gam(d1, d2)
    ! gamma, dimensionless parameter from Kerminen Kulmala 2002
    implicit none
    real(kind=dp), intent(in) :: d1, d2
    real(kind=dp), parameter :: mp1 = -0.6_dp ! m+1 (-1.6_dp + 1.0_dp)

    gam = 1.0_dp / (mp1) *( (d2/d1)**(mp1) -1.0_dp)

  end function gam


  real(kind=dp) elemental function cnuc(m,t)
    ! mean velocity of nucleating cluster
    use messy_main_constants_mem, only : pi, Navo => N_A
    implicit none
    real(kind=dp), intent(in) :: m ! molar mass (g/mol)
    real(kind=dp), intent(in) :: t ! temperature (K) 

    cnuc = sqrt( 8.0_dp *Navo *temp / (m *pi) )

  end function cnuc

end function Anttila2010

! ju_ak_20191107+ moved to smil

! below currently unused
! real(kind=dp) pure function condsink(np,diap, nair) result(y)
!   ! CS in Kerminen-Kulmala method
!   use messy_main_constants_mem, only : pi
!   implicit none
!   real(kind=dp), intent(in) :: np(:), & ! shape 3
!                              & diap(:), & ! same as np
!                              & nair ! concentration air molecules (cm-3)
!   real(kind=dp) :: kn(size(np)) ! knudsen number
!   real(kind=dp) :: Dif ! Diffusion coefficient vapour
!   real(kind=dp), parameter :: dair2 = 0.43_dp ! collision cross section air (nm**2)
! 
!   kn = 2.0_dp *1.0_dp / (dair2 *nair) /diap
! 
!   y = 4.0_dp *pi *Dif & 
!     & *0.5_dp *sum(diap *np *( (1.0_dp +kn) /(1.0_dp +0.377_dp *kn + 1.33_dp *kn *(1.0_dp +kn)) ))
!   
! end function condsink


!!$! general nucleation scheme, particles nucleate and Grow to a predefined diameter
!!$! New Particle Formation process
!!$subroutine driver_NPF(ldebug, & 
!!$             & nuclmethod, vapour, temp, rhum, tmst, coags, d2, &
!!$             & condensed, panew, nr, nucrate, &
!!$             & qion, klion, krec)
!!$
!!$  use messy_main_constants_mem, only : pi, Navo => N_A
!!$
!!$  implicit none
!!$
!!$  character(len=strlen_medium), intent(in) :: nuclmethod
!!$
!!$  !DEBUG
!!$  logical, intent(in) :: ldebug
!!$
!!$  real(kind=dp), intent(inout) :: vapour(:,:,:) ! vapour concentration (cm-3)
!!$  real(kind=dp), intent(in) :: temp(:,:), & ! temperature (K) 
!!$                               rhum(:,:), & ! relative humidity (%)
!!$                               tmst, &      ! length of base model timestep (s)
!!$                               coags(:,:), & ! needs coupling
!!$                               d2 ! needs coupling
!!$
!!$  ! ion related variables
!!$  real(kind=dp), intent(in), optional :: qion(:,:), & ! ionisation rate (ion pairs cm-3 s-1)
!!$                                         klion(:,:), & ! first order loss rate of small ions (s-1)
!!$                                         krec(:,:)   ! small ion-ion recombination rate (cm3 s-1)
!!$
!!$  real(kind=dp), intent(out) :: condensed(:,:,:), &  ! vapour tendecy (cm-3)
!!$                                panew(:,:), & ! concentration of newly formed aerosol particles (cm-3)
!!$                                nr(:,:,:),  & ! nucleation rate for each channel
!!$                                nucrate(:,:) ! total nucleation rate [cm-3 s-1] ??? #seb
!!$
!!$  ! Local variables:
!!$  real(kind=dp), dimension(1:ubound(vapour,1), 1:ubound(vapour,2)) :: gr, mnuc, d1
!!$  integer :: jk, jl
!!$
!!$  ! local parameters
!!$  real(kind=dp), parameter :: Msa = 98.0_dp ! molar mass of H2SO4
!!$  real(kind=dp), parameter :: densa = 1.83_dp ! density of "pure" H2SO4 (g/cm-3)
!!$
!!$  nucrate = 0.0_dp
!!$  nr = 0.0_dp
!!$  panew = 0.0_dp
!!$  condensed = 0.0_dp
!!$  d1 = 0.0_dp
!!$  gr = 0.0_dp
!!$
!!$  select case (nuclmethod)
!!$  case ('vehk_gmxe')
!!$
!!$
!!$    if(ldebug) then 
!!$      write(*,*) "calculate Vehkam??ki NUCRATE"
!!$      write(*,*) "h2so4: ", minval(vapour(:,:,insa)), maxval(vapour(:,:,insa))
!!$    end if
!!$
!!$    call driver_vehkamaeki( temp, rhum*0.01_dp, vapour(:,:,insa), &  ! ECHAM5 temperature, relative humidity, 
!!$                          & panew, condensed(:,:,insa), &  ! new formed particles, vapours added to nucleation mode, aerosol modes
!!$                          & nucrate, tmst )  ! nucleation rate, number of molecules in the, length timestep
!!$    if(ldebug) then 
!!$      write(*,*) "calculate Anttila2010"
!!$      write(*,*) "nucrate: ", minval(nucrate), maxval(nucrate)
!!$      write(*,*) "panew: ", minval(panew), maxval(panew)
!!$      write(*,*) "h2so4: ", minval(vapour(:,:,insa)), maxval(vapour(:,:,insa))
!!$    end if
!!$
!!$    if(ldebug) write(*,*) "calculated Vehkam??ki NUCRATE"
!!$    do jk=1, ubound(nucrate,2)
!!$      do jl=1, ubound(nucrate,1) 
!!$        if(nucrate(jl,jk) > zeps) then
!!$          mnuc(jl,jk) = Msa*condensed(jl,jk,insa) / panew(jl,jk)
!!$          d1(jl,jk) = ( mnuc(jl,jk) / Navo / densa *6.0_dp / pi )**(1.0_dp/3.0_dp) ! check units
!!$          if(d1(jl,jk) > d2) d1(jl,jk) = d2 
!!$          ! ADD SELECT CASE for various growth models
!!$          gr(jl,jk) = 2.0_dp*GR_Nieminen2010(vapour(jl,jk,insa), rhum(jl,jk))/3600.0_dp
!!$          panew(jl,jk) = nucrate(jl,jk)*tmst ! redundant
!!$          ! update nucleation rate
!!$
!!$          nucrate(jl,jk) = Anttila2010(nucrate(jl,jk), coags(jl,jk), panew(jl,jk), gr(jl,jk),  &
!!$                         & mnuc(jl,jk), d1(jl,jk), d2, temp(jl,jk))
!!$
!!$          panew(jl,jk) = nucrate(jl,jk)*tmst
!!$      !    panew = 1.0e3_dp !DEBUG
!!$          ! update number of condensed vapour
!!$          condensed(jl,jk,insa) = d2**3.0_dp *densa *pi *Navo/(6.0_dp * Msa)*panew(jl,jk)
!!$        end if
!!$      end do
!!$    end do
!!$     
!!$    if(ldebug) then
!!$      write(*,*) "DONE Antilla2010"
!!$      write(*,*) "d1: ", d1
!!$      write(*,*) "d2: ", d2
!!$      write(*,*) "gr: ", minval(gr), maxval(gr)
!!$      write(*,*) "coags: ", minval(coags), maxval(coags)
!!$      write(*,*) "nucrate: ", minval(nucrate), maxval(nucrate)
!!$    end if
!!$
!!$  case ('multi')
!!$
!!$    if(present(qion) .and. present(klion) .and. present(krec)) then
!!$ 
!!$      do jk=1, ubound(vapour,2)
!!$        do jl=1, ubound(vapour,1)
!!$ 
!!$          call nucleation_driver(vapour(jl,jk,:), temp(jl,jk), rhum(jl,jk), tmst, &
!!$         & qion(jl,jk), krec(jl,jk), klion(jl,jk), &
!!$         & nr(jl,jk,:), panew(jl,jk), nucrate(jl,jk), condensed(jl,jk,:))
!!$ 
!!$        if(nucrate(jl,jk) > zeps) then
!!$          mnuc(jl,jk) = Msa*condensed(jl,jk,insa) / panew(jl,jk)
!!$          d1(jl,jk) = ( mnuc(jl,jk) / Navo / densa *6.0_dp / pi )**(1.0_dp/3.0_dp) ! check units
!!$          if(d1(jl,jk) > d2) d1(jl,jk) = d2 
!!$          ! TODO: options for various growth models
!!$!          gr(jl,jk) = 2.0_dp*GR_Nieminen2010(vapour(jl,jk,insa), rhum(jl,jk))/3600.0_dp
!!$          gr(jl,jk) = 2.0_dp*Gordon2016(vapour(jl,jk,insa), vapour(jl,jk,inoxo1) + vapour(jl,jk,inoxo2))/3600.0_dp
!!$          panew(jl,jk) = nucrate(jl,jk)*tmst ! redundant
!!$          ! update nucleation rate
!!$ 
!!$          nucrate(jl,jk) = Anttila2010(nucrate(jl,jk), coags(jl,jk), panew(jl,jk), gr(jl,jk), mnuc(jl,jk), d1(jl,jk), &
!!$                         & d2, temp(jl,jk)) ! define temp local nrate
!!$ 
!!$          panew(jl,jk) = nucrate(jl,jk)*tmst
!!$      !    panew = 1.0e3_dp !DEBUG
!!$          ! update number of condensed vapour
!!$          condensed(jl,jk,insa) = d2**3.0_dp *densa *pi *Navo/(6.0_dp * Msa)*panew(jl,jk) ! simplification
!!$        end if
!!$        end do
!!$      end do
!!$ 
!!$    end if
!!$  end select
!!$
!!$
!!$end subroutine driver_NPF
!!$
!!$
!!$
!!$! below copied and modified from messy_gmxe.f90 MESSy v2.52
!!$subroutine driver_gmxe_nucleation( & 
!!$             & nuclmethod, vapour, temp, rhum, tmst, &
!!$             & condensed, panew, nr, nucrate, &
!!$             & qion, klion, krec)
!!$
!!$  implicit none
!!$
!!$  character(len=strlen_medium), intent(in) :: nuclmethod
!!$
!!$  real(kind=dp), intent(inout) :: vapour(:,:,:) ! vapour concentration (cm-3)
!!$  real(kind=dp), intent(in) :: temp(:,:), &     ! 
!!$                               rhum(:,:), &
!!$                               tmst
!!$
!!$  ! ion related variables
!!$  real(kind=dp), intent(in), optional :: qion(:,:), & ! ionisation rate (ion pairs cm-3 s-1)
!!$                                         klion(:,:), & ! first order loss rate of small ions (s-1)
!!$                                         krec(:,:)   ! small ion-ion recombination rate (cm3 s-1)
!!$
!!$
!!$  real(kind=dp), intent(out) :: condensed(:,:,:), &  ! vapour tendecy (cm-3)
!!$                                panew(:,:), & ! concentration of newly formed aerosol particles (cm-3)
!!$                                nr(:,:,:),  & ! nucleation rate for each channel
!!$                                nucrate(:,:) ! total nucleation rate [cm-3 s-1] ??? #seb
!!$
!!$  ! Local variables:
!!$  integer :: i, k, j
!!$  nucrate = 0.0_dp
!!$  nr = 0.0_dp
!!$  panew = 0.0_dp
!!$  condensed = 0.0_dp
!!$
!!$  select case (nuclmethod)
!!$  case ('vehk_gmxe')      
!!$    call driver_vehkamaeki( temp, rhum*0.01_dp, vapour(:,:,insa), &  ! ECHAM5 temperature, relative humidity, 
!!$                            panew, condensed(:,:,insa), &  ! new formed particles, vapours added to nucleation mode, aerosol modes
!!$                            nucrate, tmst )  ! nucleation rate, number of molecules in the, length timestep
!!$
!!$
!!$  case ('kulm_gmxe')
!!$!TBC    call driver_kulmala(kproma, klev,        &  ! ECHAM5 dimensions
!!$!TBC                         ptemp, prhum, vapour(:,:,insa), &  ! ECHAM5 temperature, relative humidity, 
!!$!TBC                         panew, pa4delt, naertot, &  ! new formed particles, vapours added to nucleation mode
!!$!TBC                         znucrate, pncrit, ztmst, &  ! nucleation rate, number of molecules in the
!!$!TBC                         ip4d(insa), nucm_Hp)
!!$!TBC                         nucm_so4m, nucm_Hp) ! index {SO4-, H+} in pa4delt
!!$
!!$
!!$  case('multi') ! multicomponent nucleation and several nucleation channels
!!$    if(present(qion) .and. present(klion) .and. present(krec)) then
!!$      do k=1, ubound(vapour,2)
!!$        do i=1, ubound(vapour,1)
!!$
!!$          call nucleation_driver(vapour(i,k,:), temp(i,k), rhum(i,k), tmst, &
!!$         & qion(i,k), krec(i,k), klion(i,k), &
!!$         & nr(i,k,:), panew(i,k), nucrate(i,k), condensed(i,k,:))
!!$
!!$        end do
!!$      end do
!!$    end if
!!$  end select ! nnucl
!!$
!!$end subroutine driver_gmxe_nucleation 
! ju_ak_20191107-

end module messy_nan

