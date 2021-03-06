!**********************************************************************
MODULE messy_clamstraj_si

#if defined(ECHAM5) || defined(MBM_CLAMS)
!**********************************************************************
!  Submodel interface for clamstraj 
!**********************************************************************
 
  ! SMCL
  USE messy_clams_global,     ONLY: paramtype, species_type, PREC, DP, &
                                    leveldt, levelfut, &
                                    UDT,UFUT,VDT,VFUT,WDT,WFUT
  USE messy_clamstraj
  USE messy_main_timer_event, ONLY: time_event, io_time_event, event_is_active 

  IMPLICIT NONE
  PRIVATE

! op_pj_20160621+ compiler error workaround for gfortran 4.8, 4.7, ... 6.2, 7.1
!!$INTRINSIC :: NULL
#if defined(__GFORTRAN__) || defined(__G95__)
#define GF_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
! Test for GF
#if GF_VERSION > 70200
  INTRINSIC :: NULL
#endif
#else
#ifndef LF
  INTRINSIC :: NULL
#endif
#endif
! op_pj_20160621-

! op_pj_20160606+
!!$  TYPE(time_event),    PUBLIC :: trajoutevent
!!$  TYPE(io_time_event), PUBLIC :: io_trajoutevent
  TYPE(time_event),    PUBLIC, SAVE :: trajoutevent
  TYPE(io_time_event), PUBLIC, SAVE :: io_trajoutevent
! op_pj_20160606-
  TYPE(time_event),    PUBLIC, SAVE :: trajinevent
  TYPE(io_time_event), PUBLIC, SAVE :: io_trajinevent

  ! MODULE VARIABLES
  REAL(PREC), DIMENSION(:), POINTER :: LAT        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON        => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV        => NULL()
  REAL(DP),   DIMENSION(:), POINTER :: JULSEC     => NULL()

  REAL(DP),                 POINTER :: JULTIME

  TYPE(paramtype), DIMENSION(:), POINTER :: PARAM

  TYPE(species_type), DIMENSION(:), POINTER :: SPECIES
  
!!!!!
  REAL(PREC), DIMENSION(:), POINTER :: dnparts_traj, dnparts_max_traj
  REAL(PREC), DIMENSION(:), POINTER :: LAT_traj   => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LON_traj   => NULL()
  REAL(PREC), DIMENSION(:), POINTER :: LEV_traj   => NULL()
  TYPE(paramtype), DIMENSION(:), POINTER :: PARAM_traj

  PUBLIC :: clamstraj_setup
  PUBLIC :: clamstraj_initialize
  PUBLIC :: clamstraj_init_memory
  PUBLIC :: clamstraj_init_coupling
  PUBLIC :: clamstraj_global_start
  PUBLIC :: clamstraj_global_end
  PUBLIC :: clamstraj_free_memory

  PRIVATE :: get_trajfilename
!!$  PRIVATE :: get_next_dnparts
  PRIVATE :: read_trajfile

!--------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------


  SUBROUTINE clamstraj_setup

    ! Read grid information from ECMWF file

    ! BMIL
    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_mpi_bi,         ONLY: p_pe

    ! SMCL
    USE messy_main_tools,          ONLY: find_next_free_unit
    USE messy_clams_tools_utils,   ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstraj_setup'

    INTEGER :: iou, status

    IF (p_pe==0) THEN
       WRITE(*,*)
       WRITE(*,*) uppercase(substr)
    ENDIF

    ! Read namelist variables:
    iou = find_next_free_unit(100,200)

    CALL clamstraj_read_nml(status, iou) 
    IF (status /= 0) CALL error_bi('Error in clamstraj_read_nml !',substr)

   
  END SUBROUTINE clamstraj_setup
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstraj_initialize

    USE messy_main_mpi_bi,      ONLY: p_pe,  p_nprocs, p_sum, p_bcast
    USE messy_main_timer_bi,    ONLY: timer_event_init
    USE messy_main_blather_bi,  ONLY: error_bi

    USE messy_main_timer,       ONLY: delta_time
    USE messy_clams_global,     ONLY: nparams, paramnames, &
                                      dnparts, nparts, idx
    USE messy_clamstraj_global, ONLY: timestep_trajout, timestep_trajin
    USE messy_clams_tools_utils,ONLY: uppercase

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstraj_initialize'

    integer, dimension(:), pointer :: dnparts_all ! (0:ntasks-1)
    INTEGER :: i
    INTEGER :: next_dnparts

    IF (p_pe==0) THEN
       WRITE(*,*)
       WRITE(*,*) uppercase(substr)
       WRITE(*,*)
    ENDIF

    ! aus clams_si: allocate param array:
    IF (nparams > 0) THEN

       ALLOCATE (PARAM(nparams))
!!!!!
       ALLOCATE (PARAM_traj(nparams))

       DO i = 1, nparams
          PARAM(i)%name = paramnames(i)
!!!!!
          PARAM_traj(i)%name = paramnames(i)
      ENDDO

    ENDIF

!!!!!
    if (timestep_trajout > 0) then
!!$       if (mod(timestep_trajout*3600,int(delta_time)) /= 0) then
!!$          call error_bi ("TRAJ output timestep must be multiple of delta_time !!!",substr)
!!$       elseif (timestep_trajout <= 0) then
!!$          call error_bi ("TRAJ output timestep must be greater than 0 !!!",substr)
!!$       endif
       ! Define TRAJOUT event:
       io_trajoutevent%counter = timestep_trajout  
       io_trajoutevent%unit = 'hours'
       io_trajoutevent%adjustment = 'exact'
       io_trajoutevent%offset = -delta_time
       CALL timer_event_init (trajoutevent, io_trajoutevent, 'TRAJOUT_Event', 'present')
    endif

    if (timestep_trajin > 0) then
       ! Define TRAJIN event:
       io_trajinevent%counter = timestep_trajin
       io_trajinevent%unit = 'hours'
       io_trajinevent%adjustment = 'exact'
       io_trajinevent%offset = -delta_time
       CALL timer_event_init (trajinevent, io_trajinevent, 'TRAJIN_Event', 'present')
    endif

    if (timestep_trajin > 0) then
       if (timestep_trajout > 0) &
            CALL error_bi ("If trajectory input is read in (timestep_trajin>0), "//&
                            "no trajectory output can be written !!!",substr)
    endif

!!$    if (timestep_trajin > 0) then
!!$
!!$       call get_next_dnparts (modstr, next_dnparts)
!!$
!!$!!!!! dnparts and nparts modified !!!
!!$       
!!$       dnparts = next_dnparts
!!$       nparts  = p_sum(dnparts)
!!$       !write (*,*) 'in clamstraj_setup: dnparts, nparts:',p_pe, dnparts, nparts
!!$
!!$!!!!!  idx modified !!!
!!$
!!$       allocate (dnparts_all(0:p_nprocs-1))
!!$
!!$       dnparts_all(p_pe) = dnparts
!!$       do i = 0, p_nprocs-1
!!$          call p_bcast (dnparts_all(i),i)
!!$       enddo
!!$       !if (p_pe==0) write (*,*) 'dnparts_all=',dnparts_all
!!$       
!!$       idx(0,1) = 1
!!$       idx(0,2) = dnparts_all(0)
!!$       do i = 1, p_nprocs-1
!!$          idx(i,1) = idx(i-1,2)+1
!!$          idx(i,2) = idx(i-1,2)+dnparts_all(i)
!!$       enddo
!!$       !if (p_pe==0) write (*,*) idx
!!$       
!!$       deallocate (dnparts_all)
!!$       
!!$    endif
       

  END SUBROUTINE clamstraj_initialize
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstraj_init_memory

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: REPR_NTASKS, REPR_LG_CLAMS
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_blather_bi,       ONLY: error_bi

    ! SMCL
    USE messy_main_channel,      ONLY: new_channel, new_channel_object,  &
                                       new_attribute

    USE messy_clams_tools_utils, ONLY: uppercase
    USE messy_clams_global,      ONLY: nparams
    USE messy_clamstraj_global

    IMPLICIT NONE

    INTEGER :: i
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstraj_init_memory'
    INTEGER :: status

    IF (p_pe==0) write(*,*) uppercase(substr)

    CALL new_channel(status, modstr, lrestreq=.TRUE.)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, modstr//'_version', c = modver)
    CALL channel_halt(substr, status)

    ! Define channel object for dnparts
    CALL new_channel_object(status, modstr, "dnparts", p1=dnparts_traj, &
                            reprid=REPR_NTASKS, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)
    ! Define channel object for dnparts_max
    CALL new_channel_object(status, modstr, "dnparts_max", p1=dnparts_max_traj, &
                            reprid=REPR_NTASKS, lrestreq=.TRUE.)
    CALL channel_halt (substr, status)

    CALL new_channel_object(status, modstr, 'LAT', p1=LAT_traj, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LON', p1=LON_traj, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'LEV', p1=LEV_traj, reprid=REPR_LG_CLAMS)
    CALL channel_halt(substr, status)

    ! Define channel objects for all parameters 

    DO i = 1, nparams
       if (p_pe==0) &
            write (*,*) 'create channel for ',trim(PARAM(i)%name)
       CALL new_channel_object(status, modstr, trim(PARAM(i)%name), &
            p1=PARAM_traj(i)%values, reprid=REPR_LG_CLAMS)
       CALL channel_halt(substr, status)
    ENDDO

  END SUBROUTINE clamstraj_init_memory
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstraj_init_coupling

    USE messy_main_channel,          ONLY: get_channel_object
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_mpi_bi,           ONLY: p_pe

    USE messy_clams_global,          ONLY: nparams, nspec, specnames, &
                                           ldiagout
    USE messy_clamstraj_global,      ONLY: timestep_trajin

    IMPLICIT NONE


    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr =  'clamstraj_init_coupling'
    INTEGER :: status
    integer :: i
 

    ! Get arrays from CLAMS submodel:
    CALL get_channel_object(status, 'clams', 'JULTIME', p0=JULTIME)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'JULSEC', p1=JULSEC)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LAT', p1=LAT)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LON', p1=LON)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'clams', 'LEV', p1=LEV)
    CALL channel_halt(substr, status)

    DO i = 1, nparams
       if (p_pe==0 .and. ldiagout) &
            write (*,*) 'couple channel object ',trim(PARAM(i)%name)
       CALL get_channel_object (status, 'clams', &
            trim(PARAM(i)%name), p1=PARAM(i)%values)
       CALL channel_halt(substr, status)
    ENDDO

!!!!!
!!$    IF (timestep_trajin > 0) THEN
!!$       ! Get species (from CLAMS):
!!$       ALLOCATE (SPECIES(nspec))
!!$       DO i = 1, nspec
!!$          if (p_pe==0 .and. ldiagout) &
!!$               write (*,*) 'couple channel object ',trim(specnames(i))
!!$          CALL get_channel_object(status, 'clams', &
!!$               trim(specnames(i)), p1=SPECIES(i)%values)
!!$          CALL channel_halt(substr, status)
!!$       ENDDO
!!$    ENDIF
    

  END SUBROUTINE clamstraj_init_coupling
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstraj_global_start

    ! BMIL
    USE messy_main_mpi_bi,          ONLY: p_pe
    USE messy_main_blather_bi,      ONLY: error_bi

    ! SMCL: MESSy
    USE messy_main_switch,         ONLY: USE_CLAMSCIRRUS,  &
                                         USE_CLAMSMIX, USE_CLAMSBMIX
!!#D clamschem +
    USE messy_main_switch,         ONLY: USE_CLAMSCHEM
!!#D clamschem -
    USE messy_main_timer,          ONLY: lstart, delta_time,  & 
                                         YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
 
    ! SMCL: CLaMS
    USE messy_clams_global,        ONLY: mdi, eps, &
                                         dnparts_max,  &
                                         pre_year, pre_month, pre_day, pre_sec, &
                                         ldiagout
    USE messy_clamstraj_global,    ONLY: forward, timestep_trajin
!!#D clamschem +
    USE messy_clamschem_global,    ONLY: timestep_chem
!!#D clamschem -
    USE messy_clamsmix_global,     ONLY: timestep_mix
    USE messy_clamsbmix_global,    ONLY: timestep_bmix
    USE messy_clamscirrus_global,  ONLY: timestep_cirrus
    USE messy_clams_tools_utils,   ONLY: uppercase, lowercase
    USE messy_clams_tools_dateconv,ONLY: ymds2js

    USE netcdf

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstraj_global_start'

    INTEGER :: sod
    INTEGER :: i

    IF (p_pe==0 .and. ldiagout) write(*,*) uppercase(substr)


    if (p_pe==0) write(*,*) 'in clamstraj_global_start: pre_date', &
         pre_year, pre_month, pre_day, pre_sec

!!!!! ???
    IF (.not. lstart) THEN
       sod = HOUR*3600 +MINUTE*60 + SECOND
       DO i=1,dnparts_max
          IF (ABS((JULSEC(i)-mdi)/mdi)<=eps) THEN
             JULSEC(i) = mdi
          ELSE
             ! Die Trajektorien koennen zu einem spaeteren Zeitpunkt starten als zur
             ! Startzeit der Berechnung/Zeitpunkt des ersten Datensatzes:
             if (forward) then
                if (JULSEC(i) <= ymds2js(YEAR, MONTH, DAY, sod)) then
                   JULSEC(i) = ymds2js(YEAR, MONTH, DAY, sod)
                endif
             else
                if (JULSEC(i) >= ymds2js(YEAR, MONTH, DAY, sod)) then
                   JULSEC(i) = ymds2js(YEAR, MONTH, DAY, sod)
                endif
             endif
                
          ENDIF
       END DO
       JULTIME = ymds2js(YEAR, MONTH, DAY, sod)
    ENDIF ! .not. lstart

!!!!!
#ifdef MBM_CLAMS
    if (lstart) then

       if (timestep_trajin > 0) then

          IF (mod(timestep_trajin*3600,int(delta_time))/=0) &
               CALL error_bi ("Wrong timestep for trajectory input !",substr)
!!#D clamschem +
          if (USE_CLAMSCHEM .and. mod(timestep_chem,timestep_trajin*3600)/=0) &
               CALL error_bi ("Wrong timestep for trajectory input !",substr)
!!#D clamschem -
          if (USE_CLAMSMIX .and. mod(timestep_mix,timestep_trajin)/=0) &
               CALL error_bi ("Wrong timestep for trajectory input !",substr)
          
          if (USE_CLAMSBMIX .and. mod(timestep_bmix,timestep_trajin)/=0) &
               CALL error_bi ("Wrong timestep for trajectory input !",substr)

          if (USE_CLAMSCIRRUS .and. mod(timestep_cirrus,timestep_trajin)/=0) &
               CALL error_bi ("Wrong timestep for trajectory input !",substr)

       endif

    endif
#endif

    
  END SUBROUTINE clamstraj_global_start
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstraj_global_end
  ! Call of TRAJ subroutine 

    ! BMIL
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_mpi_bi,           ONLY: p_pe
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_timer_bi,         ONLY: event_state

    ! SMIL
    USE messy_clamssedi_si,       ONLY: sedievent

    ! SMCL: MESSy
    USE messy_main_channel,       ONLY: set_channel_output
    USE messy_main_timer,         ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, &
                                        current_date, next_date, delta_time, &
                                        time_days, if_equal
    USE messy_main_timer_event,   ONLY: event_next_date

    ! SMCL: CLaMS
    USE messy_clams_global,       ONLY: dnparts, dnparts_max, nparams, &
                                        mdi, eps, &
                                        lchemevent, lsedievent, &
                                        lcirrusevent, lclamsoutevent,   &
                                        ldiagout, lcoupled
    USE messy_clamstraj_global,   ONLY: forward, linterpolate, &
                                        timestep_trajin, timestep_trajout
    USE messy_clams_tools_utils,  ONLY: bubble_sort, uppercase
    USE messy_clams_tools_dateconv,ONLY: ymds2js

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstraj_global_end'

    INTEGER :: i
    INTEGER :: status
    INTEGER :: sod
    REAL(DP) :: pi
    LOGICAL :: ltrajoutevent, lsedinextstep
    LOGICAL :: ltrajinevent
    TYPE(time_days) :: sedinext_date

    CHARACTER(250) :: trajfilename

    IF (p_pe==0 .and. ldiagout) write(*,*) uppercase(substr)

    pi=4.*atan(1.)                                                        



    !----------------------------------------------------
    ! check, if trajectories should be read in from files 
    !----------------------------------------------------

    if (timestep_trajin > 0) then
       ltrajinevent = event_state (trajinevent, current_date)
    else
       ltrajinevent = .FALSE.
    endif


    !---------------------------------
    ! Set linterpolate
    !---------------------------------

    linterpolate = .FALSE.

!!!!! ACHTUNG: Wann muessen die Parameter interpoliert werden ???
    ! => kostet sehr viel Rechenzeit !!!
    
!!!!!    
    ! interpolate parameters, if CLAMS- or TRAJ-output is written 
    ! in this timestep
    if (timestep_trajout > 0) then
       ltrajoutevent = event_state(trajoutevent, current_date)
    else
       ltrajoutevent = .false.
    endif
    IF (ltrajoutevent)  linterpolate = .TRUE.
    IF (lclamsoutevent) linterpolate = .TRUE.

    ! interpolate parameters, if CIRRUS or CHEM are called in this timestep
    IF (lcirrusevent) linterpolate = .TRUE.
    IF (lchemevent)   linterpolate = .TRUE.

    ! interpolate parameters, if SEDI is called in the following timestep 
    ! (SEDI is called before TRAJ !)
    lsedinextstep = .FALSE.
    call event_next_date (sedievent, sedinext_date)  
    call if_equal (sedinext_date, next_date, lsedinextstep) 
    if (p_pe==0) write (*,*) 'lsedievent=',lsedievent
    if (p_pe==0) write (*,*) 'lsedinextstep=',lsedinextstep
    if (lsedinextstep) then
        if (p_pe==0) write (*,*) 'SEDI im naechsten Zeitschritt !!!'
        linterpolate = .TRUE.
    endif
    
    if (p_pe==0) write (*,*) 'linterpolate=',linterpolate


#ifdef MBM_CLAMS
    if (timestep_trajin==0) then
#endif

       !write (*,*) 'vor traj: p_pe, dnparts=',p_pe,dnparts
       
       if (dnparts>0) then
          CALL TRAJ (status, LAT, LON, LEV, JULSEC, PARAM,lcoupled)
          IF (status/=0) call error_bi ("Error in TRAJ !!!",substr)
       else
          CALL increase_dt(delta_time,lcoupled)
       endif
       
#ifdef MBM_CLAMS
    else ! read in trajectories

       if (p_pe==0) write (*,*) 'No trajectory calculation !!!'

       if (ltrajinevent) then

          call get_trajfilename (modstr, trajfilename)
          if (p_pe==0) write (*,*) 'trajfile:',trim(trajfilename)

          if (p_pe==0) write (*,*) 'Read trajectories'
          call read_trajfile (trajfilename)

          ! get current wind data (leveldt used for interpolation!)
          CALL increase_dt(delta_time,lcoupled)

      endif

   endif
#endif
 
    ! set julian seconds to end time of timestep
    DO i=1,dnparts_max
       IF (ABS((JULSEC(i)-mdi)/mdi)<=eps) THEN
          JULSEC(i) = mdi
       ELSE
          ! Die Trajektorien koennen zu einem spaeteren Zeitpunkt starten als zur
          ! Startzeit der Berechnung/Zeitpunkt des ersten Datensatzes:
          sod = HOUR*3600 +MINUTE*60 + SECOND
          if (forward) then
             if (JULSEC(i) <= ymds2js(YEAR, MONTH, DAY, sod)) then
                JULSEC(i) = JULSEC(i) + delta_time
             endif
          else
             if (JULSEC(i) >= ymds2js(YEAR, MONTH, DAY, sod)) then
                JULSEC(i) = JULSEC(i) - delta_time
             endif
          endif

       ENDIF
    END DO
   JULTIME = JULTIME + delta_time

!!!!!
   ! -> messy_write_output in clams_main.f90
!!$   IF (loutput_traj .and. ltrajoutevent .and. timestep_trajin==0) THEN
!!$      CALL set_channel_output(status, 'clams', .FALSE.)
!!$      CALL channel_halt(substr, status)
!!$      CALL set_channel_output(status, 'clamstraj', .TRUE.)
!!$      CALL channel_halt(substr, status)
!!$      CALL messy_write_output
!!$      if (lclamsoutevent) then
!!$         CALL set_channel_output(status, 'clams', .TRUE.)
!!$         CALL channel_halt(substr, status)
!!$      endif
!!$   END IF
   dnparts_traj = dnparts
   dnparts_max_traj = dnparts_max
   LAT_traj = LAT
   LON_traj = LON
   LEV_traj = LEV
   do i = 1, nparams
      PARAM_traj(i)%values = PARAM(i)%values
   enddo

   
  END SUBROUTINE clamstraj_global_end
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE clamstraj_free_memory
    ! Deallocates data fields (in last timestep)
    
    USE messy_clamstraj_global,  ONLY: timestep_trajin

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'clamstraj_free_memory'

    !write(*,*) substr

    DEALLOCATE (PARAM)

!!!!!    
!    IF (timestep_trajin > 0) DEALLOCATE (SPECIES)


  END SUBROUTINE clamstraj_free_memory

!*******************************************************************
!
!*******************************************************************
SUBROUTINE get_trajfilename (channelname, filename)

  USE messy_main_channel,     ONLY: EXP_NAME
  USE messy_main_timer,       ONLY: YEAR_NEXT, MONTH_NEXT, DAY_NEXT,     &
                                    HOUR_NEXT, MINUTE_NEXT, SECOND_NEXT 
  USE messy_clamstraj_global, ONLY: dir_trajin


  IMPLICIT NONE
  
  CHARACTER(*)   :: channelname
  CHARACTER(*)   :: filename

  INTEGER :: i
  INTEGER :: YEAR_DATE, MONTH_DATE, DAY_DATE     
  INTEGER :: HOUR_DATE, MINUTE_DATE, SECOND_DATE 

  filename = ''

  filename = EXP_NAME

  DO i = LEN_TRIM(EXP_NAME)+1, LEN(EXP_NAME)
     WRITE(filename(i:i),'(a1)') '_' 
  END DO

  YEAR_DATE   = YEAR_NEXT
  MONTH_DATE  = MONTH_NEXT
  DAY_DATE    = DAY_NEXT
  HOUR_DATE   = HOUR_NEXT
  MINUTE_DATE = MINUTE_NEXT
  SECOND_DATE = SECOND_NEXT

  WRITE(filename(LEN(EXP_NAME)+1:)  &
       , '(a1,i4,i2.2,i2.2,a1,i2.2,i2.2,i2.2,a1)')&
       '_', YEAR_DATE, MONTH_DATE,  DAY_DATE, '_' &
       , HOUR_DATE, MINUTE_DATE, SECOND_DATE,'_'

  filename =  TRIM(dir_trajin)//"/"//TRIM(filename)//TRIM(channelname)//".nc"

END SUBROUTINE get_trajfilename

!*******************************************************************
!
!*******************************************************************
!!$SUBROUTINE get_next_dnparts (channelname, next_dnparts)
!!$
!!$  ! BMIL
!!$  USE messy_main_blather_bi,      ONLY: error_bi
!!$  USE messy_main_mpi_bi,          ONLY: p_pe
!!$
!!$  ! SMCL
!!$  USE messy_main_channel,         ONLY: EXP_NAME
!!$  USE messy_main_timer,           ONLY: YEAR_START, MONTH_START, DAY_START,     &
!!$                                        HOUR_START, MINUTE_START, SECOND_START 
!!$  USE messy_clamstraj_global,     ONLY: dir_trajin, timestep_trajin
!!$  USE messy_clams_tools_dateconv, ONLY: incdat
!!$
!!$  USE netcdf
!!$
!!$  IMPLICIT NONE
!!$  
!!$  CHARACTER(*)   :: channelname
!!$  INTEGER        :: next_dnparts
!!$
!!$  CHARACTER(LEN=*), PARAMETER :: substr = 'get_next_dnparts'
!!$
!!$  CHARACTER(250) :: filename
!!$  INTEGER :: iyear, imonth, iday, ihour, imin, isec 
!!$  INTEGER :: i, daysec
!!$  INTEGER :: status, ncid, varid
!!$  REAL(PREC) :: dhelp(1)
!!$
!!$  filename = ''
!!$
!!$  filename = EXP_NAME
!!$
!!$  DO i = LEN_TRIM(EXP_NAME)+1, LEN(EXP_NAME)
!!$     WRITE(filename(i:i),'(a1)') '_' 
!!$  END DO
!!$
!!$  iyear   = YEAR_START
!!$  imonth  = MONTH_START
!!$  iday    = DAY_START
!!$  ihour   = HOUR_START
!!$  imin    = MINUTE_START
!!$  isec    = SECOND_START
!!$
!!$  daysec = ihour*3600+imin*60+isec
!!$  call incdat (daysec,iday,imonth,iyear,timestep_trajin*3600,0,0,0)   
!!$  ihour = daysec/3600
!!$  imin  = (daysec - ihour*3600)/60
!!$  isec  = daysec - ihour*3600 - imin*60
!!$  if (p_pe==0) write (*,*) 'ihour,imin,isec',ihour,imin,isec
!!$
!!$  WRITE(filename(LEN(EXP_NAME)+1:)  &
!!$       , '(a1,i4,i2.2,i2.2,a1,i2.2,i2.2,i2.2,a1)')&
!!$       '_', iyear, imonth,  iday, '_', ihour, imin, isec,'_'
!!$
!!$  filename =  TRIM(dir_trajin)//"/"//TRIM(filename)//TRIM(channelname)//".nc"
!!$  if (p_pe==0) write (*,*) 'read dnparts from file: ',trim(filename)
!!$
!!$  status = nf90_open (trim(filename), nf90_nowrite, ncid)
!!$  if (status /= nf90_noerr) &
!!$       call error_bi ('Can not open file '//trim(filename)//' !', substr)
!!$
!!$  status = nf90_inq_varid (ncid, "dnparts", varid)
!!$  if (status /= nf90_noerr) &
!!$       call error_bi ('Can not find variable dnparts !', substr)
!!$  status = nf90_get_var (ncid, varid, dhelp, start=(/p_pe+1,1/), count=(/1,1/))
!!$  if (status /= nf90_noerr) &
!!$       call error_bi ('Can not read variable dnparts !', substr)
!!$  next_dnparts = dhelp(1)
!!$  !write (*,*) 'in sub. get_next_dnparts: dnparts=',p_pe, next_dnparts
!!$  
!!$  status = nf90_close (ncid)
!!$
!!$END SUBROUTINE get_next_dnparts

!*******************************************************************
!
!*******************************************************************
SUBROUTINE read_trajfile (filename)

  ! BMIL
  USE messy_main_blather_bi,     ONLY: error_bi
  USE messy_main_mpi_bi,         ONLY: p_pe, p_sum

  ! SMCL
  USE messy_clams_global,       ONLY: dnparts, dnparts_max, nparts, &
                                      nparams, ldiagout

  USE netcdf

  IMPLICIT NONE

  CHARACTER(*) :: filename

  CHARACTER(LEN=*), PARAMETER :: substr='read_trajfile'
  INTEGER :: status, ncid, varid
  INTEGER :: dnparts_max_check, dnparts_check
  INTEGER :: i

  if (p_pe==0 .and. ldiagout) write (*,*) 'read from ',trim(filename)
  status = nf90_open (filename, nf90_nowrite, ncid)
  if (status /= nf90_noerr) &
       call error_bi ('Can not open file '//trim(filename)//' !', substr)

!!!!!
!!$  ! read dnparts (current number of trajectories on rank)
!!$  status = nf90_inq_varid (ncid,'dnparts',varid)
!!$  if (status /= 0) call error_bi ('Can not find variable dnparts!', substr)
!!$  status = nf90_get_var (ncid, varid, dnparts,start=(/p_pe+1/))
!!$  if (status /= 0) call error_bi ('Can not read variable dnparts!', substr)
!!$  if (ldiagout) write (*,*) 'rank, dnparts=', p_pe, dnparts
  ! check dnparts (current number of trajectories on rank)
  status = nf90_inq_varid (ncid,'dnparts',varid)
  if (status /= 0) call error_bi ('Can not find variable dnparts!', substr)
  status = nf90_get_var (ncid, varid, dnparts_check,start=(/p_pe+1/))
  if (status /= 0) call error_bi ('Can not read variable dnparts!', substr)
  if (ldiagout) write (*,*) 'rank, dnparts, dnparts_check =', p_pe, dnparts, dnparts_check
  if (dnparts_check /= dnparts) call error_bi ('dnparts has changed !!!', substr)

!!$  ! set nparts (current number of trajectories on all ranks)
!!$  nparts = p_sum (dnparts)
  if (p_pe==0 .and. ldiagout) write (*,*) 'nparts=',nparts

  ! check dnparts_max 
  status = nf90_inq_varid (ncid,'dnparts_max',varid)
  if (status /= 0) call error_bi ('Can not find variable dnparts_max!', substr)
  status = nf90_get_var (ncid, varid, dnparts_max_check)
  if (status /= 0) call error_bi ('Can not read variable dnparts_max!', substr)
  if (p_pe==0 .and. ldiagout) write (*,*) 'dnparts_max=',dnparts_max_check
  if (dnparts_max_check /= dnparts_max) call error_bi ('dnparts_max has changed !!!', substr)

  ! read latitudes
  if (p_pe==0 .and. ldiagout) write (*,*) 'read lat'
  status = nf90_inq_varid (ncid,'LAT',varid)
  if (status /= 0) call error_bi ('Can not find variable LAT!', substr)
  status = nf90_get_var (ncid,varid,LAT,start=(/p_pe*dnparts_max+1/),count=(/dnparts_max/))
  IF (status /= 0) CALL error_bi('Cannot read variable LAT ',substr)

  ! read longitudes
  if (p_pe==0 .and. ldiagout) write (*,*) 'read lon'
  status = nf90_inq_varid (ncid,'LON',varid)
  if (status /= 0) call error_bi ('Can not find variable LON!', substr)
  status = nf90_get_var (ncid,varid,LON,start=(/p_pe*dnparts_max+1/),count=(/dnparts_max/))
  IF (status /= 0) CALL error_bi('Cannot read variable LON ',substr)

  ! read level
  if (p_pe==0 .and. ldiagout) write (*,*) 'read lev'
  status = nf90_inq_varid (ncid,'LEV',varid)
  if (status /= 0) call error_bi ('Can not find variable LEV!', substr)
  status = nf90_get_var (ncid,varid,LEV,start=(/p_pe*dnparts_max+1/),count=(/dnparts_max/))
  IF (status /= 0) CALL error_bi('Cannot read variable LEV ',substr)

  ! read parameters (PARAM)
  DO i = 1, nparams

     if (p_pe==0 .and. ldiagout) write (*,*) 'read ',trim(PARAM(i)%name)
     status = nf90_inq_varid (ncid,PARAM(i)%name,varid)
     if (status /= 0) call error_bi ('Can not find variable '//trim(PARAM(i)%name), substr)
     status = nf90_get_var (ncid,varid,PARAM(i)%values,&
          start=(/p_pe*dnparts_max+1/),count=(/dnparts_max/))
     IF (status /= 0) CALL error_bi('Cannot read variable '//trim(PARAM(i)%name),substr)

  ENDDO
 
  status = nf90_close (ncid)

END SUBROUTINE read_trajfile


#endif
!**********************************************************************
END MODULE messy_clamstraj_si
!**********************************************************************
