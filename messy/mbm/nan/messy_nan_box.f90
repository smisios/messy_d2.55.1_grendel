! box model for messy_gmxe_nucleation
! author: sebastian ehrhart
! email: s.ehrhart@mpic.de
! date: 24. Oct 2016
! last modified:   

module box

  use messy_main_constants_mem, only: dp, strlen_ulong

  implicit none

  private

  public :: read_box_nml, filenames, number_of_lines, time_integration, recom_brasseur

  character(len=*), parameter, public :: modstr = 'box'
  character(len=*), parameter, public :: modver = '1.0'

  ! input via name list
  logical, public, save :: ltimesim, lexternal, lsteady
  character(len=10), public, save :: outsuf='dat'
  character(len=strlen_ulong), dimension(100), public, save :: exfile=''
  real(kind=dp), public, save :: temp, rh ! temperature and RH in Klevin and %
  real(kind=dp), public, save :: time=7200_dp, tstep=60_dp ! total simulation time and time step (s)
  real(kind=dp), public, save :: qion=1.84_dp, krec=1.6e-6_dp, klion=0.002 ! 

  real(kind=dp), dimension(10), public, save :: constart ! concentration at start of simulation (cm-3)
  real(kind=dp), dimension(10), public, save :: cs ! condensational sink or wall loss for each species
  real(kind=dp), dimension(10), public, save :: sol ! condensational sink or wall loss for each species
  real(kind=dp), dimension(10), public, save :: production ! production rates of species (cm-3 s-1)


  contains
  ! read namelist for the box model
subroutine read_box_nml(status, iou)

  use messy_main_tools, only : read_nml_open, read_nml_check, read_nml_close

  implicit none

  integer, intent(out) :: status ! error status: 0 = no error, 1 = an error
  integer, intent(in)  :: iou    ! I/O unit

  namelist /CTRL/ ltimesim, lexternal, lsteady !, lgmxe
  namelist /CONDI/ temp, rh, constart
  namelist /TIMER/ time, tstep
  namelist /SRCES/ production
  namelist /LOSS/ cs, sol
  namelist /IONS/ qion, krec, klion
  namelist /EXT/ exfile
  ! local
  character(len=*), parameter :: substr='box_read_nml'
  logical :: lex
  integer :: fstat

  status=1

  ! check if a rewind statement has to be issued to make the namelist
  ! reading independent of position

  call read_nml_open(lex, substr, iou, 'CTRL', modstr)
  if(.not.lex) return ! can't open the file

  read(iou, nml=CTRL, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'CTRL', modstr)
  if(fstat /= 0) return ! error while reading the namelist

!  if(ltimesim .or. .not. lexternal) then

    ! ambient conditions
    read(iou, nml=CONDI, iostat=fstat)
    call read_nml_check(fstat, substr, iou, 'CONDI', modstr)
    if(fstat /= 0) return ! error while reading the namelist

    !---------- select timestep
    read(iou, nml=TIMER, iostat=fstat)
    call read_nml_check(fstat, substr, iou, 'TIMER', modstr)
    if(fstat /= 0) return ! error while reading the namelist
  
!  end if

  !------ production terms
  read(iou, nml=SRCES, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'SRCES', modstr)
  if(fstat /= 0) return

  read(iou, nml=LOSS, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'LOSS', modstr)
  if(fstat /= 0) return

  !------ ions
  read(iou, nml=IONS, iostat=fstat)
  call read_nml_check(fstat, substr, iou, 'IONS', modstr)
  if(fstat /= 0) return

!  if(lexternal .and. lsteady) then

    read(iou, nml=EXT, iostat=fstat)
    call read_nml_check(fstat, substr, iou, 'EXT', modstr)
    if(fstat /= 0) return

!  end if

  call read_nml_close(substr, iou, modstr)

end subroutine read_box_nml


subroutine filenames(fend, fin, fout)
  implicit none
  character(len=10), intent(in) :: fend
  character(len=strlen_ulong), intent(out) :: fin(:)
  character(len=strlen_ulong), intent(out) :: fout(:)
  ! local
  integer :: i, pdot
  character(len=strlen_ulong) :: ftemp, prefix

  do i=1, size(fin)
    ftemp = exfile(i)
    pdot = index(trim(ftemp), ".", back=.true.)
    if(pdot == 0) pdot = len_trim(ftemp) ! if filename has no dot in it
    fin(i) = ftemp
    prefix = ftemp(1:pdot)
    fout(i) = trim(prefix) // trim(fend)
  end do

end subroutine filenames

function number_of_lines(iou)
  use :: iso_fortran_env
  implicit none
  integer, intent(in) :: iou
  character(len=strlen_ulong) :: s
  integer :: readstat
  integer :: number_of_lines 

  number_of_lines = 0
  do 
    read(iou,'(A)', iostat=readstat) s
    if(readstat /= 0) then
      if(readstat == iostat_end ) then
        return
      else
        write(*,*) "READ ERROR", readstat
        stop
      end if
    end if
    number_of_lines = number_of_lines + 1
  end do

end function number_of_lines


subroutine time_integration(tmax, dt, &
                            & production, cs, sol, &
                            & qion, krec, klion, &
                            & vapour, temp, rh, nr, panew, totrate, condensed)
  use messy_nan, only : nucleation_driver
  implicit none
  real(kind=dp), intent(in)    :: tmax, dt ! maximum time and delta t between steps
  real(kind=dp), intent(in)    :: temp, rh ! temeprature (K), RH (%)
  real(kind=dp), intent(in)    :: production(:), cs(:), sol(:) ! production, first and second order losses
  real(kind=dp), intent(in)    :: qion, krec, klion ! ion pair production, recombination and first order loss term (e.g. walls)
  real(kind=dp), intent(inout) :: vapour(:) ! vapour concentration
  real(kind=dp), intent(out)   :: panew, totrate, nr(:), condensed(:) ! newly produced particles, total formation rate, 
                                   ! nucleation rate, condensed vapours
  integer :: i

  do i=1, floor(tmax/dt)
     vapour = vapour + (production - cs*vapour - sol*vapour*vapour)*dt
     call nucleation_driver(vapour, temp, rh, dt, qion, krec, klion, nr, panew, totrate, condensed)
  end do 

end subroutine time_integration

elemental function recom_brasseur(t,m)
! recombination rate parameterisation given in Brasseur & Chatel 1983
  implicit none
  ! input
  real(kind=dp), intent(in) :: t ! temperature in K
  real(kind=dp), intent(in) :: m ! total concentration of air in cm-3
  ! local
  real(kind=dp) :: isct ! inverse scaled temperature, i.e. 300/t
  ! return
  real(kind=dp) :: recom_brasseur

  isct = 300.0d0/t
  recom_brasseur = 6.0d-8*sqrt(isct) + 6.0d-26*m*isct*isct*isct*isct

end function recom_brasseur
end module box


program boxmodel

  use box
  use messy_nan, only : nv => nnucspec, read_nucleation_nml, nucleation_driver, & 
 & initialize_nucleation, nc=>nnucchan
  use messy_main_constants_mem, only : dp, strlen_ulong

  implicit none

  integer :: i, i2, k, status, nt, nf, nl
  real(kind=dp), allocatable :: vapour(:), nr(:), condensed(:)
  ! below variables to test gmxe interface
  real(kind=dp), allocatable :: r3vapour(:,:,:), ptemp(:,:), r3nr(:,:,:),  & 
                               & prhum(:,:), r2panew(:,:), pa4delt(:,:,:), &
                               & r2qion(:,:), r2klion(:,:), r2krec(:,:), pncrit(:,:)
  real(kind=dp) ::  panew, totrate, ipr 
  real(kind=dp) :: current_time
  real(kind=dp) :: jexp, jsd
  character(len=strlen_ulong), allocatable :: infiles(:), outfiles(:)

  call read_nucleation_nml(status, 10)

  call initialize_nucleation()
  
  write(*,*) "allocating memory: ", " nv=", nv
! variables below have to be non allocatable arrays for g95
!  allocate(constart(nv))
!  allocate(production(nv))
!  allocate(cs(nv))
!  allocate(sol(nv))
  allocate(vapour(nv)); allocate(nr(nc)); allocate(condensed(nv))
  write(*,*) "memory allocated"

  constart=0.0_dp; production=0.0_dp; cs=0.0_dp; sol=0.0_dp
  call read_box_nml(status, 11)
  vapour = constart(1:nv)

  write(*,*) constart

  if(ltimesim) then

    nt = floor(time/tstep)

    current_time=0_dp

    ! time integration 
    do i=1, nt
!      vapour = vapour + (production(1:nv) - cs*vapour - sol*vapour*vapour)*tstep
      call time_integration(tstep, 1.0_dp, production(1:nv), cs(1:nv), sol(1:nv), &
                          &  qion, krec, klion, &
                          &  vapour, temp, rh, nr, panew, totrate, condensed)
      current_time = current_time + tstep
     write(*,*) current_time, (vapour(k), k=1, size(vapour)), &
     & (nr(k), k=1, size(nr)), panew, totrate, (condensed(k), k=1, size(condensed))
    end do


  elseif((.not. lexternal) .and. lsteady) then

    ! use constart as steady state value and provide one number, quick testing
    write(*,*) 'calculate steady state nucleation rates'
    tstep=1.0_dp ! we use 1 s so the concentration of new particles has the same value as the formation rate

    call time_integration(3600.0_dp, tstep, production(1:nv), cs(1:nv), sol(1:nv), &
                         & qion, krec, klion, &
                         &  vapour, temp, rh, nr, panew, totrate, condensed)
    write(*,*) temp, rh, (vapour(k), k=1, size(vapour)), (nr(k), k=1, size(nr)), panew, totrate, &
   & (condensed(k), k=1, size(condensed))
    write(*,*) "TOTAL PARTICLE FORMATION RATE= ", totrate 


  elseif(lexternal .and. lsteady) then

    nf = count(len_trim(exfile) > 0)
    allocate(infiles(1:nf)); allocate(outfiles(1:nf))

    call filenames(outsuf,infiles, outfiles)

    ! open file, space seperated columns, entry order:
    ! temp, rh, ipr, vapour(:), J, Jsd
    fileloop: do i=1,nf
      open(12, file=trim(infiles(i))) 
      open(22, file=trim(outfiles(i)))
      write(*,*) "opened file ", trim(infiles(i))
      nl = number_of_lines(12) ! determine number of lines

      rewind(12)

      ! read the file line by line and calculate nucleation rates in steady state
      nucloop: do i2=1, nl
        read(12, *) temp, rh, ipr, (vapour(k), k=1, size(vapour)), jexp, jsd
        production(1:nv) = vapour*cs
        qion = ipr
        krec = recom_brasseur(temp,7.337132e+21_dp/temp)
        tstep = 1.0_dp
        call time_integration(3600.0_dp, tstep, production(1:nv), cs(1:nv), sol(1:nv), &
                           &  qion, krec, klion, &
                            vapour, temp, rh, nr, panew, totrate, condensed)
          write(22,*) temp, rh, qion, (vapour(k), k=1, size(vapour)), jexp, jsd, totrate, krec, &
         & (nr(k), k=1, size(nr))
      end do nucloop

      ! clean up
      close(12)
      close(22)

    end do fileloop

  end if

!  deallocate(constart); deallocate(production); deallocate(cs); deallocate(sol)
  deallocate(vapour); deallocate(nr); deallocate(condensed)

end program boxmodel
