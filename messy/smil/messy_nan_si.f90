#include "messy_main_ppd_bi.inc" 
! defines ranks and indices for different base models, might not be needed
! Submodel interface

!!! Nucleation submodel interface layer for MESSy
!!! Author: Sebastian Ehrhart 

module messy_nan_si


  use messy_main_constants_mem, only: wp=>dp, strlen_medium
  use messy_main_channel, only : t_chaobj_cpl
  use messy_nan

  use messy_main_tools,         only: PTR_3D_ARRAY

#ifdef MESSYTENDENCY
  use messy_main_tendency_bi,   only: mtend_get_handle, mtend_register, &
                                      mtend_get_start_l, mtend_id_t, &
                                      mtend_id_q,  mtend_id_xl,      &
                                      mtend_id_xi, mtend_id_tracer,  &
                                      mtend_add_l
#endif

  implicit none
  private

  public :: nan_initialize
  public :: nan_init_memory
  public :: nan_init_coupling
  public :: nan_radiation
  public :: nan_physc
  public :: nan_free_memory


  real(kind=wp), allocatable, save :: vapour(:,:,:,:) 
  real(kind=wp), allocatable, save :: cond(:,:,:,:)

  ! for ion induced nucleation
  real(kind=wp), dimension(:,:,:), pointer, save :: total_ipr => null()
  real(kind=wp), dimension(:,:,:), pointer, save :: krec => null()
  real(kind=wp), dimension(:,:,:), pointer, save :: klion => null()

  ! for Antilla method
  real(kind=wp), dimension(:,:,:), pointer, save :: coags => null()
  real(kind=wp), dimension(:),     pointer, save :: d2 => null() 

  ! channel objects
  real(kind=wp), dimension(:,:,:), pointer, save :: panew => null() ! new particles
  real(kind=wp), dimension(:,:,:), pointer, save :: nr => null() ! nucleation rate
  real(kind=wp), dimension(:,:,:,:), pointer, save :: nrc => null() ! was pointer
  type(PTR_3D_ARRAY), dimension(:), pointer, save :: mem => null()


  integer, allocatable, save :: idx_vap(:) ! index of nucleating vapours in tracer list
  integer, allocatable, save :: idx_con(:) ! condensed phase tracer index
  integer, save :: idx_aer    ! index of nucleation mode aerosol particle number concentration
  integer, allocatable, save :: idx_LTERP(:)

  ! coupling
  logical, save :: lnucten = .FALSE. ! should this submodel calculate tendencies for nucleation mode particles?
  logical, save :: lgrow = .FALSE. ! should nucleated aerosols be grown to specific size? requires cpl_d2, cpl_coags

  character(len=strlen_medium), save :: aerosolmodel = 'GMXE' 
  character(len=strlen_medium), save :: driver_call = 'physc'

  logical, save :: b4gmxe = .TRUE.
  logical, save :: lnucmode = .TRUE. 

  type(t_chaobj_cpl), save :: cpl_d2 ! diameter to which nucleation rate is extrapolated (nm)
  type(t_chaobj_cpl), save :: cpl_coags ! coagulation sink  

  ! coupling to ion model
  type(t_chaobj_cpl), save :: cpl_ipr  ! coupling to total ion pair production rate
  type(t_chaobj_cpl), save :: cpl_krec ! coupling to recombination rate constant
  type(t_chaobj_cpl), save :: cpl_klion !coupling to ion losse to aerosol particles

  ! coupling to basemodel objects
  real(kind=wp), pointer, dimension(:,:,:) :: press_3d   => NULL()
  real(kind=wp), pointer, dimension(:,:,:) :: rhum_3d    => NULL()
  real(kind=wp), pointer, dimension(:,:,:) :: grvol      => NULL() ! m3
  real(kind=wp), pointer, dimension(:,:,:) :: grmass     => NULL() ! kg

#ifdef MESSYTENDENCY
  integer :: mt_handle
#endif


contains
!! ---------------- 1. Initialisation

subroutine nan_initialize

  use messy_main_mpi_bi,     only: p_parallel_io, p_io, p_bcast
  use messy_main_blather_bi, only: start_message_bi, end_message_bi, error_bi
  use messy_main_tools,      only: find_next_free_unit

  implicit none
  character(len=*), parameter :: substr="nan_initialize"
  integer :: status
  integer :: iou ! I/O unit
  integer :: i

  if(p_parallel_io) call start_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

  ! Read name list file
  if(p_parallel_io) then
    iou = find_next_free_unit(100,200)
    call read_nucleation_nml(status, iou)
    if(status /= 0) call error_bi(' ',substr)
  end if

  ! defines the nucleation scheme used
  ! CTRL
  do i=1, size(lselnuc)
    call p_bcast(lselnuc(i), p_io)
  end do
  call p_bcast(nnucspec, p_io) 
  call p_bcast(nuclmethod, p_io)

  ! NUC
  call p_bcast(sulphuric_acid, p_io)
  call p_bcast(ammonia,        p_io)
  call p_bcast(amines,         p_io)
  call p_bcast(HOMOH,          p_io)
  call p_bcast(HOMO3,          p_io)
  do i=1, size(vapour_aero_names)
    call p_bcast(vapour_aero_names(i), p_io)
  end do

  ! the parameterisation
  ! PARAM
  call p_bcast(dunne_pbn, p_io) 
  call p_bcast(dunne_ubn, p_io)
  call p_bcast(dunne_vbn, p_io)
  call p_bcast(dunne_wbn, p_io)
  call p_bcast(dunne_ptn, p_io)
  call p_bcast(dunne_utn, p_io)
  call p_bcast(dunne_vtn, p_io)
  call p_bcast(dunne_wtn, p_io)
  call p_bcast(dunne_pAn, p_io)
  call p_bcast(dunne_an , p_io)
  call p_bcast(dunne_pbi, p_io)
  call p_bcast(dunne_ubi, p_io)
  call p_bcast(dunne_vbi, p_io)
  call p_bcast(dunne_wbi, p_io)
  call p_bcast(dunne_pti, p_io)
  call p_bcast(dunne_uti, p_io)
  call p_bcast(dunne_vti, p_io)
  call p_bcast(dunne_wti, p_io)
  call p_bcast(dunne_pAi, p_io)
  call p_bcast(dunne_ai , p_io)
  call p_bcast(frhc1,      p_io)
  call p_bcast(frhc2,      p_io)
  call p_bcast(l_Dunne_RH, p_io)
  call p_bcast(BtdOrg,     p_io) 
  call p_bcast(l_org_Tdep, p_io)

! Organic nucleation parameterisations
  do i=1, size(kirkby_a)
    call p_bcast(kirkby_a(i), p_io) ! pure biogenic from Kirkby et al
  end do
!  do i=1, size(ricco_k)
    call p_bcast(ricco_k, p_io)  ! Riccobono et al 2014
!  end do

  if(p_parallel_io) write(*,*) "call initialize_nucleation()"

  call initialize_nucleation() ! initialise all module variables
  if(p_parallel_io) then
    write(*,*) "NUCLEATION MODEL SPECIES INDEX :"
    write(*,*) "INSA   = ", insa 
    write(*,*) "INNH3  = ", innh3
    write(*,*) "INDAM  = ", indma
    write(*,*) "INOXO1 = ", inoxo1
    write(*,*) "INOXO2 = ", inoxo2
  end if
  if(p_parallel_io) write(*,*) "AFTER call initialize_nucleation()"

  if(p_parallel_io) then
    write(*,*) "READ NUCLEATION COUPLING"
    iou = find_next_free_unit(100,200)  
    call nan_read_nml_cpl(status, iou)
    if(status /= 0) call error_bi(' ',substr)
  end if

  ! arosol ion interaction
  ! CPL
  call p_bcast(cpl_ipr%cha,  p_io)
  call p_bcast(cpl_ipr%obj,  p_io)

  call p_bcast(cpl_krec%cha, p_io)
  call p_bcast(cpl_krec%obj, p_io)

  call p_bcast(cpl_klion%cha, p_io)
  call p_bcast(cpl_klion%obj, p_io)

  call p_bcast(cpl_coags%cha, p_io)
  call p_bcast(cpl_coags%obj, p_io)

  call p_bcast(cpl_d2%cha, p_io)
  call p_bcast(cpl_d2%obj, p_io)

  call p_bcast(aerosolmodel, p_io)
  call p_bcast(lnucten, p_io)
  call p_bcast(driver_call, p_io)
  call p_bcast(b4gmxe, p_io)
  call p_bcast(lnucmode, p_io)


  if(p_parallel_io) call end_message_bi(modstr,'INITIALIZATION',substr)
 
end subroutine nan_initialize


subroutine nan_init_memory

  use messy_main_blather_bi,       only: start_message_bi, end_message_bi
  use messy_main_mpi_bi,           only: p_parallel_io
  use messy_main_grid_def_mem_bi,  only: nproma, nlev, ngpblks
  use messy_main_channel_error_bi, only: channel_halt
  use messy_main_channel_bi,       only: GP_3D_MID, DC_GP  &
                                      , DIMID_LON, DIMID_LAT, DIMID_LEV &
                                      , DC_BC &
                                      , gp_nseg, gp_start, gp_cnt &
                                      , gp_meml, gp_memu &
                                      , GP_2D_HORIZONTAL
  ! MESSy
  use messy_main_channel,    only: new_channel, new_channel_object, new_attribute
!!$  use messy_main_channel_dimensions, only: new_dimension
!!$  use messy_main_channel_repr,       only: new_representation, AUTO &
!!$                                           , set_representation_decomp &
!!$                                           , IRANK, PIOTYPE_COL

  implicit none

  ! local
  integer :: status
  integer :: i
  integer :: nmod

  character(len=*), parameter :: substr="nan_init_memory"

  nmod = nnucchan ! careful here, becomes obsolete if nmods from GMXe need to be used

  if(p_parallel_io) then
    call start_message_bi(modstr, "CHANNEL DEFINITION", substr)
    write(*,*) "NNUCSPEC =", nnucspec
  end if

  allocate(vapour(_RI_XYZN_(nproma,ngpblks,nlev,nnucspec)))
  allocate(cond(_RI_XYZN_(nproma,ngpblks,nlev,nnucspec)))
  allocate(nrc(_RI_XYZN_(nproma,ngpblks,nlev,nmod)))
  allocate(idx_vap(nnucspec))
  allocate(idx_con(nnucspec))
  allocate(mem(nnucchan))

  vapour  = 0.0_wp
  idx_vap = 0
  idx_con = 0

#ifdef MESSYTENDENCY
  call mtend_register (mt_handle,mtend_id_t)
  call mtend_register (mt_handle,mtend_id_tracer)
#endif


  ! define new channels
  call new_channel(status, modstr, reprid=GP_3D_MID)
  call channel_halt(substr, status)

  ! 1. number of freshly nucleated aerosols
  call new_channel_object(status, modstr, 'panew', p3 = panew, lrestreq=.TRUE.)
  call channel_halt(substr, status)

  ! 2. nucleation rate
  call new_channel_object(status, modstr, 'nucrate', p3 = nr, lrestreq=.TRUE.)
  call channel_halt(substr, status)

  do i=1, size(chid) 
   ! mem(i)%ptr => nrc(:,:,i,:)
    call new_channel_object(status, modstr, trim(nucchannames(chid(i))), p3=mem(i)%ptr)
    call channel_halt(substr, status)
  end do


  call end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

end subroutine nan_init_memory


!!! collect only information not provided by the calling aerosol model
subroutine nan_init_coupling

  use messy_main_blather_bi,       only: start_message_bi, end_message_bi
  use messy_main_channel,          only: get_channel_object 
  use messy_main_channel_error_bi, only: channel_halt 
  use messy_main_tracer,           only: get_tracer
  use messy_main_blather_bi,       only: error_bi, info_bi
  use messy_main_tracer_mem_bi,    only: GPTRSTR
  use messy_main_tools,            only: strcrack
  use messy_main_mpi_bi,           only: p_parallel_io
  use messy_main_data_bi,          only: basemod => modstr

  implicit none

  ! local variables here only
  character(len=*), parameter :: substr="nan_init_coupling"
  integer :: status, ierr, i, n
  character(len=strlen_medium), pointer     :: outstring(:) => null()

  call start_message_bi(modstr, 'COUPLING', substr)

  ! collect tracer IDs
  ! gas phase
  do i=1, size(vapour_names)
    if(trim(vapour_names(i)) /= '') then
      call get_tracer(ierr, GPTRSTR, trim(vapour_names(i)), idx=idx_vap(i))
      if(ierr /= 0) call error_bi("TRACER " // trim(vapour_names(i)) // " NOT FOUND!", substr)
    else
      call error_bi("VAPOUR NAME NOT DEFINED", substr)
    end if 

    if(trim(condph_names(i)) /= '') then 
      call get_tracer(ierr, GPTRSTR, trim(condph_names(i)), subname='ns', idx=idx_con(i))
      if (ierr /= 0) call error_bi("TRACER " // trim(condph_names(i)) // " NOT FOUND!", substr)
!    if(p_parallel_io) then
!       write(*,*) "SEARCHED FOR TRACER ", trim(vapour_names(i)), " TRACER ID ", idx_vap(i)
!       write(*,*) "  CORRESPONDING TRACER IN THE CONDENSED PHASE ", & !, trim(condph_names(i)), &
!                & " TRACER ID ", idx_con(i)
    else
      call error_bi("VAPOUR NAME in NUCL. MODE NOT DEFINED", substr)
    end if
  end do

  ! get nucleation mode tracer
  if(lnucten) then 
    call get_tracer(ierr, GPTRSTR, "N", subname="ns", idx=idx_aer)
!    if(p_parallel_io) write(*,*) 'IERR = ', ierr
    if(ierr /=0 ) call error_bi('TRACER N_ns not found', substr)
    if(p_parallel_io) then
      write(*,*) "NUCLEATION submodel will calculate tendencies for nucleation mode tracers."
      write(*,*) "N_NS INDEX : ", idx_aer
    end if
  end if


  ! coupling to ion model
  if(trim(cpl_ipr%cha)/='') then
    call info_bi('Coupling ion concentration: '//cpl_ipr%cha//'  '//cpl_ipr%obj)
    call get_channel_object(status, trim(cpl_ipr%cha), trim(cpl_ipr%obj), p3=total_ipr)
    call channel_halt(substr, status)
  end if

  if(trim(cpl_krec%cha)/='') then
    call info_bi('Coupling ion concentration: '//cpl_krec%cha//'  '//cpl_krec%obj)
    call get_channel_object(status,  trim(cpl_krec%cha), trim(cpl_krec%obj), p3=krec)
    call channel_halt(substr, status)
  end if

  if(trim(cpl_klion%cha)/='') then
    call info_bi('Coupling ion losses to aerosol particles: '//cpl_klion%cha//'  '//cpl_klion%obj)
    call get_channel_object(status, trim(cpl_klion%cha), trim(cpl_klion%obj), p3=klion)
    call channel_halt(substr, status)
  end if

  ! channels from GMXe
  if(trim(cpl_coags%cha)/='') then
    call info_bi('Coupling coagulation sink :'//cpl_coags%cha//'  '//cpl_coags%obj)
    call get_channel_object(status, trim(cpl_coags%cha), trim(cpl_coags%obj), p3=coags)
    call channel_halt(substr, status)
  end if

  if(trim(cpl_d2%cha)/='') then
    call info_bi('Coupling lowest aerosol diameters sink :'//cpl_d2%cha//'  '//cpl_d2%obj)
    call get_channel_object(status, trim(cpl_d2%cha), trim(cpl_d2%obj), p1=d2)
    call channel_halt(substr, status)
  end if


  ! channels from the basemodel
  call get_channel_object(status, basemod,'press', p3=press_3d)
  if(status == 1) &
    & call error_bi('channel object for press not found', substr)

  call get_channel_object(status, basemod, 'rhum', p3=rhum_3d) 
  if(status == 1) &
    & call error_bi('channel object for rhum not found', substr)

  call get_channel_object(status, 'grid_def', 'grvol', p3=grvol)
  if(status == 1) &
    & call error_bi('channel object for grvol not found', substr)

  call get_channel_object(status, 'grid_def', 'grmass', p3=grmass)
  if(status == 1) &
    & call error_bi('channel object for grmass not found', substr)

  call end_message_bi(modstr, 'COUPLING', substr)

end subroutine nan_init_coupling



subroutine nan_read_nml_cpl(status, iou)

  use messy_main_tools, only: read_nml_open, read_nml_check, read_nml_close

  implicit none
  integer, intent(in) :: iou
  integer, intent(out) :: status

  ! LOCAL
  character(len=*), parameter :: substr='nan_read_nml_cpl'
  logical                     :: lex      ! file exists ?
  integer                     :: fstat    ! file status

  namelist /CPL/ cpl_ipr,  cpl_krec, cpl_klion, cpl_coags, cpl_d2, aerosolmodel, &
               & lnucten, driver_call, b4gmxe, lnucmode  
  status = 1

  call read_nml_open(lex, substr, iou, 'CPL', modstr)
  if (.not.lex) return    ! <modstr>.nml does not exist

  read(iou, nml=CPL, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'CPL', modstr)
  if (fstat /= 0) return  ! error while reading namelist

  call read_nml_close(substr, iou, modstr)
  status = 0

end subroutine nan_read_nml_cpl

subroutine nan_radiation(flag)
  implicit none
  integer :: flag
  select case(flag)
  case(1)
     call nan_radiation_b4
  case(2)
     call nan_radiation_after
  end select
end subroutine nan_radiation

subroutine nan_radiation_b4
  ! called from basemodel where radiation processes are calculated.
  ! needed to be called with gmxe, i.e. nucleation
  ! called before gmxe
  implicit none

  if(trim(adjustl(driver_call)) == 'radiation' .and. b4gmxe) & 
    & call nan_driver_si

end subroutine nan_radiation_b4


subroutine nan_radiation_after
  ! called from basemodel where radiation processes are calculated.
  ! needed to be called with gmxe, i.e. nucleation
  ! called before gmxe
  implicit none

  if(trim(adjustl(driver_call)) == 'radiation' .and. .not. b4gmxe) call nan_driver_si

end subroutine nan_radiation_after



subroutine nan_physc
  ! called from basemodel physics submodel.
  ! needed to be called with gmxe, i.e. nucleation
  ! called before gmxe
  implicit none

  if(trim(adjustl(driver_call)) == 'physc') call nan_driver_si


end subroutine nan_physc



subroutine nan_driver_si

#ifndef MESSYTENDENCY
  use messy_main_tracer_mem_bi, only: pxtte => qxtte, pxtm1 => qxtm1
  use messy_main_data_bi,       only: tm1, tte_3d
  use messy_main_mpi_bi,        only: p_parallel_io
#endif 
  use messy_main_timer,         only: delta_time, tmst => time_step_len &
                                    , nstep => current_time_step
  use messy_main_mpi_bi,        only: p_parallel_io
  use messy_main_grid_def_mem_bi, only:   nlev       &
                                       ,nproma, kproma     &
                                       ,nrow => ngpblks    &
                                       ,jrow               
  use messy_main_constants_mem,  only: k_B, M_air

  implicit none
  ! get tracers
  character(len=*), parameter :: substr = 'nan_driver_si'

  real(kind=wp), dimension(1:kproma,1:nlev) :: temp, & ! temperature in Kelvin
                                             & conv, & ! conversion from mol/mol to molecules cm-3
                                             & rhoa    ! density of air kg/m3
  integer :: i, idtg, idta

#ifdef MESSYTENDENCY
  real(kind=wp), dimension(:,:), pointer :: xtp1_0 
#endif
   !--- Ambient properties: ----------------------------------


 ! only run this model after GMXE started
 if(nstep > 0) then

#ifndef MESSYTENDENCY    


  temp(:,:) = tm1(_RI_XYZ__(:,jrow,:))   +       &
             & tte_3d(_RI_XYZ__(:,jrow,:)) * tmst
  conv = press_3d(_RI_XYZ__(:,jrow,:))/(k_B*temp)*1.0e-6_wp ! number density in # cm-3

  rhoa = grmass(_RI_XYZ__(:,jrow,:)) / grvol(_RI_XYZ__(:,jrow,:))


  do i=1, size(idx_vap)
    idtg = idx_vap(i)
    vapour(_RI_XYZN_(:,jrow,:,i)) = max(pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idtg)) & 
                  & + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idtg)) * tmst, 0.0_wp ) &
                  & *conv
  end do


#else

  call mtend_get_start_l(mtend_id_t,  v0 = temp)

  conv = press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))/(k_B*temp)*1.0e-6_wp ! number density in # cm-3

  do i=1, size(idx_vap) 
    idtg = idx_vap(i)
!    call mtend_get_start_l(mtend_id_tracer, v0 = xtp1_0, idt=idtg) ! we don't actually need all of them
    call mtend_get_start_l(idtg, v0 = xtp1_0) ! we don't actually need all of them
    vapour(_RI_XYZN_(:,jrow,:,i)) = max(xtp1_0(:,:), 0.0_wp) *conv(:,:)
  end do 

#endif

  ! call the nucleation submodel

!  if( p_parallel_io ) write(*,*) 'CALCULATE AEROSOL NUCLEATION FOR ', aerosolmodel
  select case(aerosolmodel)
    case('GMXE')

      if(associated(total_ipr)) then
        if(lnucmode) then
          call driver_gmxe_nucleation( &
                 & nuclmethod, vapour(_RI_XYZN_(:,jrow,:,:)), temp, rhum_3d(_RI_XYZ__(:,jrow,:)), tmst, &
                 & cond(_RI_XYZN_(:,jrow,:,:)), panew(_RI_XYZ__(:,jrow,:)), nrc(_RI_XYZN_(:,jrow,:,:)), nr(_RI_XYZ__(:,jrow,:)), &
                 & qion=total_ipr(_RI_XYZ__(:,jrow,:)), klion=klion(_RI_XYZ__(:,jrow,:)), krec=krec(_RI_XYZ__(:,jrow,:)) )

        else
          call driver_NPF( p_parallel_io, & 
                 & nuclmethod, vapour(_RI_XYZN_(:,jrow,:,:)), temp, rhum_3d(_RI_XYZ__(:,jrow,:)), tmst, &
                 & coags(_RI_XYZ__(:,jrow,:)), d2(1), &
                 & cond(_RI_XYZN_(:,jrow,:,:)), panew(_RI_XYZ__(:,jrow,:)), nrc(_RI_XYZN_(:,jrow,:,:)), nr(_RI_XYZ__(:,jrow,:)), &
                 & qion=total_ipr(_RI_XYZ__(:,jrow,:)), klion=klion(_RI_XYZ__(:,jrow,:)), krec=krec(_RI_XYZ__(:,jrow,:)) )
        end if
      else
 !DEBUG       if( p_parallel_io ) write(*,*) 'TOTAL IPR NOT PRESENT'
!         call driver_gmxe_nucleation( &
!               & nuclmethod, vapour(_RI_XYZ__(:,jrow,:),:), temp, rhum_3d(_RI_XYZ__(:,jrow,:)), tmst, &
!               & cond(_RI_XYZ__(:,jrow,:),:), panew(_RI_XYZ__(:,jrow,:)), nrc(_RI_XYZ__(:,jrow,:),:), nr(_RI_XYZ__(:,jrow,:)))
      end if
    case('MADE')
    ! NOTHING YET
  end select


  do i=1, size(chid)
    mem(i)%ptr(_RI_XYZ__(:,jrow,:)) = nrc(_RI_XYZN_(:,jrow,:,i)) ! check for basemodel independence
  end do

! update tracer tendencies
!DEBUG  if( p_parallel_io ) write(*,*) 'NUCLEATION: TENDENCIES X)'
#ifndef MESSYTENDENCY
  do i=1, size(idx_vap)
    idtg = idx_vap(i)
    idta = idx_con(i) 
    pxtte(_RI_X_ZN_(1:kproma,1:nlev,idtg)) = &
                                & pxtte(_RI_X_ZN_(1:kproma,1:nlev,idtg)) - cond(_RI_XYZN_(1:kproma,jrow,1:nlev,i))/conv(:,:)/tmst ! + zxtte(:,:,i) 
    pxtte(_RI_X_ZN_(1:kproma,1:nlev,idta)) = & 
                                & pxtte(_RI_X_ZN_(1:kproma,1:nlev,idta)) + cond(_RI_XYZN_(1:kproma,jrow,1:nlev,i))/conv(:,:)/tmst
  end do 
  if(lnucten) &
    & pxtte(_RI_X_ZN_(1:kproma,1:nlev,idx_aer)) = pxtte(_RI_X_ZN_(1:kproma,1:nlev,idx_aer)) & 
    &  + panew(_RI_XYZ__(1:kproma,jrow,1:nlev))/tmst *M_air/0.001_wp/rhoa  ! needs a conversion from per cm-3 s-1 to  mol-1 s-1 
#else
  do i=1, size(idx_vap)
    idtg = idx_vap(i)
    idta = idx_con(i) 
!!$    call mtend_add_l (mt_handle, mtend_id_tracer, px = -1.0_wp* cond(:,:,jrow,i)/conv(:,:)/tmst, idt=idtg)
!!$    call mtend_add_l (mt_handle, mtend_id_tracer, px = cond(:,:,jrow,i)/conv(:,:)/tmst, idt=idta)
    call mtend_add_l (mt_handle, idtg, px = -1.0_wp* cond(_RI_XYZN_(:,jrow,:,i))/conv(:,:)/tmst)
    call mtend_add_l (mt_handle, idta, px = cond(_RI_XYZN_(:,jrow,:,i))/conv(:,:)/tmst)
  end do
  if(lnucten) &
!    & call mtend_add_l ( mt_handle, mtend_id_tracer, px=panew(:,:,jrow)/tmst* (M_air/(rhoa*0.001_wp)), idt=idx_aer ) ! needs a conversion from per cm-3 s-1 to mol-1 s-1
    & call mtend_add_l ( mt_handle, idx_aer, px=panew(_RI_XYZ__(:,jrow,:))/tmst* (M_air/(rhoa*0.001_wp))) ! needs a conversion from per cm-3 s-1 to mol-1 s-1
#endif

  end if

end subroutine nan_driver_si

!!$
!!$! below currently unused
!!$! real(kind=dp) pure function condsink(np,diap, nair) result(y)
!!$!   ! CS in Kerminen-Kulmala method
!!$!   use messy_main_constants_mem, only : pi
!!$!   implicit none
!!$!   real(kind=dp), intent(in) :: np(:), & ! shape 3
!!$!                              & diap(:), & ! same as np
!!$!                              & nair ! concentration air molecules (cm-3)
!!$!   real(kind=dp) :: kn(size(np)) ! knudsen number
!!$!   real(kind=dp) :: Dif ! Diffusion coefficient vapour
!!$!   real(kind=dp), parameter :: dair2 = 0.43_dp ! collision cross section air (nm**2)
!!$! 
!!$!   kn = 2.0_dp *1.0_dp / (dair2 *nair) /diap
!!$! 
!!$!   y = 4.0_dp *pi *Dif & 
!!$!     & *0.5_dp *sum(diap *np *( (1.0_dp +kn) /(1.0_dp +0.377_dp *kn + 1.33_dp *kn *(1.0_dp +kn)) ))
!!$!   
!!$! end function condsink
!!$
!!$
! general nucleation scheme, particles nucleate and Grow to a predefined diameter
! New Particle Formation process
subroutine driver_NPF(ldebug, & 
             & nuclmethod, vapour, temp, rhum, tmst, coags, d2, &
             & condensed, panew, nr, nucrate, &
             & qion, klion, krec)

  use messy_main_constants_mem, only : pi, Navo => N_A

  implicit none

  character(len=strlen_medium), intent(in) :: nuclmethod

  !DEBUG
  logical, intent(in) :: ldebug

  real(kind=wp), intent(inout) :: vapour(:,:,:) ! vapour concentration (cm-3)
  real(kind=wp), intent(in) :: temp(:,:), & ! temperature (K) 
                               rhum(:,:), & ! relative humidity (%)
                               tmst, &      ! length of base model timestep (s)
                               coags(:,:), & ! needs coupling
                               d2 ! needs coupling

  ! ion related variables
  real(kind=wp), intent(in), optional :: qion(:,:), & ! ionisation rate (ion pairs cm-3 s-1)
                                         klion(:,:), & ! first order loss rate of small ions (s-1)
                                         krec(:,:)   ! small ion-ion recombination rate (cm3 s-1)

  real(kind=wp), intent(out) :: condensed(:,:,:), &  ! vapour tendecy (cm-3)
                                panew(:,:), & ! concentration of newly formed aerosol particles (cm-3)
                                nr(:,:,:),  & ! nucleation rate for each channel
                                nucrate(:,:) ! total nucleation rate [cm-3 s-1] ??? #seb

  ! Local variables:
  real(kind=wp), dimension(1:ubound(vapour,1), 1:ubound(vapour,2)) :: gr, mnuc, d1
  integer :: jk, jl

  ! local parameters
  real(kind=wp), parameter :: Msa = 98.0_wp ! molar mass of H2SO4
  real(kind=wp), parameter :: densa = 1.83_wp ! density of "pure" H2SO4 (g/cm-3)

  nucrate = 0.0_wp
  nr = 0.0_wp
  panew = 0.0_wp
  condensed = 0.0_wp
  d1 = 0.0_wp
  gr = 0.0_wp

  select case (nuclmethod)
  case ('vehk_gmxe')


    if(ldebug) then 
      write(*,*) "calculate Vehkamäki NUCRATE"
      write(*,*) "h2so4: ", minval(vapour(_RI_X_ZN_(:,:,insa))), maxval(vapour(_RI_X_ZN_(:,:,insa)))
    end if

    call driver_vehkamaeki( temp, rhum*0.01_wp, vapour(_RI_X_ZN_(:,:,insa)), &  ! ECHAM5 temperature, relative humidity, 
                          & panew, condensed(_RI_X_ZN_(:,:,insa)), &  ! new formed particles, vapours added to nucleation mode, aerosol modes
                          & nucrate, tmst )  ! nucleation rate, number of molecules in the, length timestep
    if(ldebug) then 
      write(*,*) "calculate Anttila2010"
      write(*,*) "nucrate: ", minval(nucrate), maxval(nucrate)
      write(*,*) "panew: ", minval(panew), maxval(panew)
      write(*,*) "h2so4: ", minval(vapour(_RI_X_ZN_(:,:,insa))), maxval(vapour(_RI_X_ZN_(:,:,insa)))
    end if

    if(ldebug) write(*,*) "calculated Vehkamäki NUCRATE"
    do jk=1, ubound(nucrate,2)
      do jl=1, ubound(nucrate,1) 
        if(nucrate(jl,jk) > zeps) then
          mnuc(jl,jk) = Msa*condensed(_RI_X_ZN_(jl,jk,insa)) / panew(jl,jk)
          d1(jl,jk) = ( mnuc(jl,jk) / Navo / densa *6.0_wp / pi )**(1.0_wp/3.0_wp) ! check units
          if(d1(jl,jk) > d2) d1(jl,jk) = d2 
          ! ADD SELECT CASE for various growth models
          gr(jl,jk) = 2.0_wp*GR_Nieminen2010(vapour(_RI_X_ZN_(jl,jk,insa)), rhum(jl,jk))/3600.0_wp
          panew(jl,jk) = nucrate(jl,jk)*tmst ! redundant
          ! update nucleation rate

          nucrate(jl,jk) = Anttila2010(nucrate(jl,jk), coags(jl,jk), panew(jl,jk), gr(jl,jk),  &
                         & mnuc(jl,jk), d1(jl,jk), d2, temp(jl,jk))

          panew(jl,jk) = nucrate(jl,jk)*tmst
      !    panew = 1.0e3_wp !DEBUG
          ! update number of condensed vapour
          condensed(_RI_X_ZN_(jl,jk,insa)) = d2**3.0_wp *densa *pi *Navo/(6.0_wp * Msa)*panew(jl,jk)
        end if
      end do
    end do
     
    if(ldebug) then
      write(*,*) "DONE Antilla2010"
      write(*,*) "d1: ", d1
      write(*,*) "d2: ", d2
      write(*,*) "gr: ", minval(gr), maxval(gr)
      write(*,*) "coags: ", minval(coags), maxval(coags)
      write(*,*) "nucrate: ", minval(nucrate), maxval(nucrate)
    end if

  case ('multi')

    if(present(qion) .and. present(klion) .and. present(krec)) then
 
      do jk=1, ubound(vapour,_IZ_XYZ__)
        do jl=1, ubound(vapour,1)
 
          call nucleation_driver(vapour(_RI_X_ZN_(jl,jk,:)), temp(jl,jk), rhum(jl,jk), tmst, &
         & qion(jl,jk), krec(jl,jk), klion(jl,jk), &
         & nr(_RI_X_ZN_(jl,jk,:)), panew(jl,jk), nucrate(jl,jk), condensed(_RI_X_ZN_(jl,jk,:)))
 
        if(nucrate(jl,jk) > zeps) then
          mnuc(jl,jk) = Msa*condensed(jl,jk,insa) / panew(jl,jk)
          d1(jl,jk) = ( mnuc(jl,jk) / Navo / densa *6.0_wp / pi )**(1.0_wp/3.0_wp) ! check units
          if(d1(jl,jk) > d2) d1(jl,jk) = d2 
          ! TODO: options for various growth models
!          gr(jl,jk) = 2.0_wp*GR_Nieminen2010(vapour(jl,jk,insa), rhum(jl,jk))/3600.0_wp
          gr(jl,jk) = 2.0_wp*Gordon2016(vapour(_RI_X_ZN_(jl,jk,insa)), vapour(_RI_X_ZN_(jl,jk,inoxo1)) + vapour(_RI_X_ZN_(jl,jk,inoxo2)))/3600.0_wp
          panew(jl,jk) = nucrate(jl,jk)*tmst ! redundant
          ! update nucleation rate
 
          nucrate(jl,jk) = Anttila2010(nucrate(jl,jk), coags(jl,jk), panew(jl,jk), gr(jl,jk), mnuc(jl,jk), d1(jl,jk), &
                         & d2, temp(jl,jk)) ! define temp local nrate
 
          panew(jl,jk) = nucrate(jl,jk)*tmst
      !    panew = 1.0e3_wp !DEBUG
          ! update number of condensed vapour
          condensed(_RI_X_ZN_(jl,jk,insa)) = d2**3.0_wp *densa *pi *Navo/(6.0_wp * Msa)*panew(jl,jk) ! simplification
        end if
        end do
      end do
 
    end if
  end select


end subroutine driver_NPF



! below copied and modified from messy_gmxe.f90 MESSy v2.52
subroutine driver_gmxe_nucleation( & 
             & nuclmethod, vapour, temp, rhum, tmst, &
             & condensed, panew, nr, nucrate, &
             & qion, klion, krec)

  implicit none

  character(len=strlen_medium), intent(in) :: nuclmethod

  real(kind=wp), intent(inout) :: vapour(:,:,:) ! vapour concentration (cm-3)
  real(kind=wp), intent(in) :: temp(:,:), &     ! 
                               rhum(:,:), &
                               tmst

  ! ion related variables
  real(kind=wp), intent(in), optional :: qion(:,:), & ! ionisation rate (ion pairs cm-3 s-1)
                                         klion(:,:), & ! first order loss rate of small ions (s-1)
                                         krec(:,:)   ! small ion-ion recombination rate (cm3 s-1)


  real(kind=wp), intent(out) :: condensed(:,:,:), &  ! vapour tendecy (cm-3)
                                panew(:,:), & ! concentration of newly formed aerosol particles (cm-3)
                                nr(:,:,:),  & ! nucleation rate for each channel
                                nucrate(:,:) ! total nucleation rate [cm-3 s-1] ??? #seb

  ! Local variables:
  integer :: i, k, j
  nucrate = 0.0_wp
  nr = 0.0_wp
  panew = 0.0_wp
  condensed = 0.0_wp

  select case (nuclmethod)
  case ('vehk_gmxe')      
    call driver_vehkamaeki( temp, rhum*0.01_wp, vapour(_RI_X_ZN_(:,:,insa)), &  ! ECHAM5 temperature, relative humidity, 
                            panew, condensed(_RI_X_ZN_(:,:,insa)), &  ! new formed particles, vapours added to nucleation mode, aerosol modes
                            nucrate, tmst )  ! nucleation rate, number of molecules in the, length timestep


  case ('kulm_gmxe')
!TBC    call driver_kulmala(kproma, klev,        &  ! ECHAM5 dimensions
!TBC                         ptemp, prhum, vapour(:,:,insa), &  ! ECHAM5 temperature, relative humidity, 
!TBC                         panew, pa4delt, naertot, &  ! new formed particles, vapours added to nucleation mode
!TBC                         znucrate, pncrit, ztmst, &  ! nucleation rate, number of molecules in the
!TBC                         ip4d(insa), nucm_Hp)
!TBC                         nucm_so4m, nucm_Hp) ! index {SO4-, H+} in pa4delt


  case('multi') ! multicomponent nucleation and several nucleation channels
    if(present(qion) .and. present(klion) .and. present(krec)) then
      do k=1, ubound(vapour,_IZ_XYZ__)
        do i=1, ubound(vapour,1)

          call nucleation_driver(vapour(_RI_X_ZN_(i,k,:)), temp(i,k), rhum(i,k), tmst, &
         & qion(i,k), krec(i,k), klion(i,k), &
         & nr(i,k,:), panew(i,k), nucrate(i,k), condensed(_RI_X_ZN_(i,k,:)))

        end do
      end do
    end if
  end select ! nnucl

end subroutine driver_gmxe_nucleation 


subroutine nan_free_memory

  implicit none

  ! allocatables
  if(allocated(vapour))    deallocate(vapour)
  if(allocated(idx_vap))   deallocate(idx_vap)
  if(allocated(idx_con))   deallocate(idx_con)
  if(allocated(cond))      deallocate(cond)
!  ! channel objects
!  if(associated(panew))     deallocate(panew)
!  if(associated(nr))        deallocate(nr)
  if(associated(nrc))       deallocate(nrc)

end subroutine nan_free_memory

end module messy_nan_si
