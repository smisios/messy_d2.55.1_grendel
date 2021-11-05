MODULE mo_iolist

#ifndef NO_NEW_IO

#ifndef MESSY
  USE mo_iso_c_kinds, ONLY: c_int64_t
#else
  USE mo_kind, ONLY : c_int64_t => i8, c_float => sp
#endif

  USE mo_commo1, ONLY: ldays, lmonts, lyears
  USE mo_commo1, ONLY: amsuo, amsue, tiestu, tiestw, weto
  USE mo_param1, ONLY: ie, ie_g, je, je_g, ke
  USE mo_parallel, ONLY: p_io, p_pe
  USE mo_units, only: io_in_octl, io_stdout

  USE mo_parallel, ONLY: gather, global_sum, p_bcast, scatter, stop_all
  USE mo_util_string, ONLY: tolower

  use mo_model_time
  use mo_file_list
  use mo_average_list
  USE mo_varlist
  use mo_postprocess

#ifdef _PROFILE
  USE mo_profile,      ONLY: trace_start, trace_stop
#endif


  IMPLICIT NONE

#include "pointer_intent_macro.inc"
  INCLUDE 'cdi.inc'

  !>
  !! Global list for MPIOM data files.
  !!
  TYPE(file_list), SAVE :: ocean_iolist

  !>
  !! Global collector field for output.
  !!
  REAL(wp), ALLOCATABLE, PRIVATE :: field_g(:,:)

  !>
  !! Buffer for I/O config data.
  !!
  type(file_data), allocatable :: config_data(:)

CONTAINS

  !>
  !! Module initialization
  !!
  SUBROUTINE iolist_init

    IF(p_pe == p_io) THEN
       ALLOCATE(field_g(ie_g, je_g))
    ELSE
       ALLOCATE(field_g(0, 0))
    END IF

  END SUBROUTINE iolist_init

  !>
  !! Module finalization
  !!
  SUBROUTINE iolist_fini

    DEALLOCATE(field_g)

  END SUBROUTINE iolist_fini

  !>
  !! Read I/O configuration from name list.
  !!
  !! Read configuration from name list IOCTL as read_element
  !! using the I/O process, broadcast data to all other processes.
  !!
  !! Example:
  !! iolist(1)=99,'Z1.nc','NC',3,4,2,5,1,82,13,15,35,36,141,501,502,503,504,111,110,7,99,9,10,84
  !! iolist(2)=99,'Z2.nc','NC',3,4,2,5,1,82,13,15,35,36,141,501,502,503,504,111,110,7,99,9,10,84
  !! ...
  !! This list may be empty, but must be present.
  !!
  SUBROUTINE iolist_read_config(file_unit, ierror)

    integer, intent(in) :: file_unit
    INTEGER, intent(inout) :: ierror

    TYPE(file_data) :: iolist(100)

    NAMELIST /IOCTL/ iolist

    IF(p_pe == p_io) THEN
        READ(file_unit, IOCTL, iostat=ierror)
    ELSE
        ierror = 0
    end if

    !   broadcast each iolist element
    CALL broadcast_readelement(iolist)

    ! Store data in buffer
    allocate(config_data(COUNT(iolist%type /= 0)))
    config_data = iolist(1:size(config_data))

  END SUBROUTINE iolist_read_config

  !>
  !! Read I/O configuration from name list into internal representation.
  !!
  !! Read configuration from name list IOCTL as read_element
  !! using the I/O process, broadcast data to all other processes,
  !! and convert to internal list representation (iolist).
  !!
  SUBROUTINE iolist_create_file_list(the_file_list, the_varlist)

    ! read ioctl name list, eg.
    ! iolist(1)=99,'Z1.nc','NC',3,4,2,5,1,82,13,15,35,36,141,501,502,503,504,111,110,7,99,9,10,84
    ! iolist(2)=99,'Z2.nc','NC',3,4,2,5,1,82,13,15,35,36,141,501,502,503,504,111,110,7,99,9,10,84
    ! ....
    ! can be empty, but gives error if not present

    TYPE(file_list), INTENT(inout) :: the_file_list
    TYPE (varlist), POINTERINTENT(in)  :: the_varlist

    INTEGER :: n

    ! Check if buffer is set
    IF(.not. allocated(config_data)) then
        CALL stop_all('Internal: iolist_read_config must be called before iolist_create_file_list')
    end if

    ! Set I/O node info in list, reading from buffer.
    call file_list_set_on_io_node(the_file_list, p_pe == p_io)
    DO n = 1, size(config_data)
       ! add a readelement to the iolist; check for at least one code in list!
       IF (config_data(n)%codes(1).NE.0) &
            CALL file_list_add(the_file_list, config_data(n), the_varlist)
    ENDDO

    ! Remove buffer
    deallocate(config_data)

   !   create the iostream and generate the informatons needed by cdi (e.g. axis, varid, etc)
    IF ( p_pe == p_io ) CALL create_avg_stream(the_file_list)

  END SUBROUTINE iolist_create_file_list

  !>
  !! Finalize file list by closing open files.
  !!
  !! Currently restart files are skipped.
  !!
  subroutine iolist_close_file_list(this)
    type(file_list), intent(inout) :: this
    type(file_list) :: current

    if(this%on_io_node) then

        current = this
        do while(file_list_is_not_empty(current))

            if(all(current%head%type /= (/ file_restart, file_final, file_initial /)) .and. &
               current%head%streamid /= -1) then
#ifndef NOCDI
               call streamClose(current%head%streamid)
#endif
                current%head%streamid = -1
            end if

            current = file_list_tail(current)
        end do

    end if ! on I/O node

  end subroutine iolist_close_file_list

  !>
  !! Distribute configuration data from name list.
  !!
  SUBROUTINE broadcast_readelement(iolist)

    !   broadcast all iolist read elements

    TYPE (file_data)   :: iolist(100),packed(100)
    INTEGER            :: i,n,ncount
    INTEGER            :: bufType
    CHARACTER(len=128) :: bufName
    CHARACTER(len=3)   :: bufFormat
    INTEGER            :: bufCodes(255)

    ! remove empty entrys from iolist

    IF ( p_pe == p_io ) THEN

      packed(:)%type=0
      n=1

      DO i=1,SIZE(iolist)
        IF (iolist(i)%type /= 0 ) THEN
          packed(n)=iolist(i)
          n=n+1
        ENDIF
      ENDDO
      iolist(:)=packed(:)

    ENDIF

    IF ( p_pe == p_io ) ncount=COUNT(iolist%type.NE.0)
    CALL p_bcast(ncount,p_io)

    DO n=1,ncount

       !        PRINT*,n

       bufType=iolist(n)%type
       bufName=iolist(n)%name
       bufFormat=iolist(n)%format
       bufCodes=iolist(n)%codes

       CALL p_bcast(bufType,p_io)
       CALL p_bcast(bufName,p_io)
       CALL p_bcast(bufFormat,p_io)
       CALL p_bcast(bufCodes,p_io)

       iolist(n)%type=bufType
       iolist(n)%name=bufName
       iolist(n)%format=bufFormat
       iolist(n)%codes(:)=bufCodes(:)

    ENDDO

  END SUBROUTINE broadcast_readelement

  !>
  !! Opens stream for data files, associates variable list with streams.
  !!
  !! Does not handle write-once files
  !! (%type in { file_restart, file_final, file_initial })
  !!
  SUBROUTINE create_avg_stream(the_file_list)

    TYPE (file_list), INTENT(in) :: the_file_list

    TYPE (file_list) :: current
    INTEGER :: file_type

    ! transverse the list and print the values
    current = the_file_list ! make current as alias of list
    DO WHILE ( file_list_is_not_empty(current) )

       IF ( all(current%head%type /= (/ file_restart, file_final, file_initial /)) ) THEN

          SELECT CASE (tolower(current%head%format))
          CASE ('ext')
            file_type = filetype_ext
          CASE ('nc')
            file_type = filetype_nc
          CASE ('nc2')
            file_type = filetype_nc2
          CASE ('nc4')
            file_type = filetype_nc4
          CASE ('grb','grbsz','sz')
            file_type = filetype_grb
          CASE default
            CALL stop_all('io format in ioctl should be &
                 &ext, nc, nc2, nc4, grb, grbsz or sz')
          END SELECT

          WRITE(0,*) TRIM(current%head%name),' open for write'
#ifndef NOCDI
          current%head%streamID = &
               streamOpenWrite(TRIM(current%head%name), file_type)
          IF (current%head%streamID == cdi_elibnavail &
               .AND. file_type == filetype_nc4) THEN
            WRITE(0, *) "netCDF 4 output unavailable, retrying with netCDF 2"
            file_type = filetype_nc2
            current%head%streamID = &
                 streamOpenWrite(TRIM(current%head%name), file_type)
          END IF
          IF (tolower(current%head%format) == 'grbsz' &
               .OR.  tolower(current%head%format) == 'sz') THEN
!            WRITE(0,*)'streamDefZtype ',COMPRESS_SZIP
            CALL streamDefZtype(current%head%streamID,COMPRESS_SZIP)
          ENDIF

          IF ( current%head%streamID < 0 ) THEN
             WRITE(0,*) cdiStringError(current%head%streamID)
             WRITE(0,*) TRIM(current%head%name),' open failed '
             CALL stop_all(' ')
          ENDIF

          CALL streamDefVlist(current%head%streamID, current%head%vlistID)
#endif

!          WRITE (0, '(A,A)') &
!               "created average stream: ", &
!               trim(file_item_to_string(current%head, .true.))

       ENDIF ! type /= file_restart, file_final, file_initial

       current = file_list_tail(current)
    END DO

  END SUBROUTINE create_avg_stream

  !>
  !! Frontend routine for write averaged data files.
  !!
  SUBROUTINE iolist_poll_write_avglist(the_file_list, the_varlist, model_time)
    TYPE (file_list), INTENT(inout) :: the_file_list
    TYPE (varlist), POINTERINTENT(in) :: the_varlist
    type(time_desc), intent(in) :: model_time
    TYPE (file_list) :: current
    current = the_file_list
    DO WHILE(file_list_is_not_empty(current))
       IF( all(current%head%type /= (/ file_restart, file_final, file_initial /)) .AND. &
           is_end_of_averaging(current%head%type, model_time) ) &
            CALL iolist_element_write(current, the_varlist, model_time)
       current = file_list_tail(current)
    END DO
  END SUBROUTINE iolist_poll_write_avglist

  !>
  !! Write appropriate data into output file.
  !!
  SUBROUTINE iolist_element_write(the_file_list, the_varlist, model_time)

    TYPE (file_list), INTENT(inout) :: the_file_list
    TYPE (varlist), POINTERINTENT(in) :: the_varlist
    type(time_desc), intent(in) :: model_time

    TYPE (average_list) :: current_avglist
    TYPE (file_item), POINTER :: file
    TYPE (average_item), POINTER :: avg
    TYPE (varlist_element) :: var

    INTEGER :: level, ignored
    REAL(wp) :: scalar
    REAL(wp) :: field1(1)


#ifndef NOCDI

    file => the_file_list%head

    IF (p_pe == p_io) THEN
       WRITE (0, '("iolist_element_write: nstep = ", I0, ", name = ", A)') &
           file%nstep, TRIM(file%name)
       CALL taxisDefVdate(file%taxisID, get_date(model_time))
       CALL taxisDefVtime(file%taxisID, get_time(model_time))
       ignored = streamDefTimestep(file%streamID, file%nstep)
    END IF

    current_avglist = file%avglist
    DO WHILE(average_list_is_not_empty(current_avglist))
       avg => current_avglist%head

       var = get_var_by_code(the_varlist, avg%code) ! get data from code number

       IF (ASSOCIATED (var%scalar)) THEN

          ! Field to be collected is scalar (allready on p_io).
#ifdef _PROFILE
          CALL trace_start ('streamwritevar', 14)
#endif
          IF (p_pe == p_io) THEN
             ! Get appropriate source variable.
             if(file%snapshot) then
                scalar = var%scalar
             else
                scalar = avg%field(1,1,1)
             end if
             ! Scaling factor is only applied for non-restart files.
             if( var%factor /= 1.0_dp .and. file%type /= file_restart .and. &
                 scalar /= avg%missval ) then
                scalar = scalar * var%factor
             end if
             CALL streamWriteVar( file%streamID, avg%varID, scalar, &
                                  count((/scalar == avg%missval/)) )
          END IF
#ifdef _PROFILE
          CALL trace_stop ('streamwritevar', 14)
#endif

       ELSE IF( (associated(var%array_1d) .or. associated(var%array_2d)) .AND. &
                var%grid == grid_zonal_mean ) THEN

          ! field is a zonal mean section (allready on p_io).
#ifdef _PROFILE
          CALL trace_start ('streamwritevar', 14)
#endif
          IF (p_pe == p_io) THEN
             DO level = 1, SIZE(avg%field, 2)
                ! Scaling factor is only applied for non-restart files.
                if( var%factor /= 1.0_dp .and. file%type /= file_restart ) then
                   where(avg%field(:,level,1) /= avg%missval)
                      avg%field(:,level,1) = avg%field(:,level,1) * var%factor
                   end where
                end if
                CALL streamWriteVarSlice(file%streamID, avg%varID, level-1, &
                     avg%field(:,level,1), &
                     COUNT( avg%field == avg%missval))
             END DO
          END IF
#ifdef _PROFILE
          CALL trace_stop ('streamwritevar', 14)
#endif

       ELSE ! 2D or 3D base variable, but no zonal mean

          ! Field to be collected is scattered horizontal 2D or
          ! field to be collected is scattered "pile" of 2D fields.
          DO level = 1, var%i_shape(3)

             ! Snapshots of integration variables
             ! need to be accumulated exactly once.
             if( file%snapshot .and. &
                  (var%column_integrate .or. var%surface_integrate) ) then
                CALL avglist_element_start(current_avglist)
                CALL avglist_element_accumulate(current_avglist, var, 1)
             end if

             if(var%surface_integrate) then

                ! With surface integration we need to compute the sum.
                ! The result may have a Z axis; for now we use usual levels.

                field1(1) = avg%field(1,1,level)
                call global_sum(field1(1))
                if(p_pe == p_io) then
                   ! Scaling factor is only applied for non-restart files.
                   if( var%factor /= 1.0_dp .and. file%type /= file_restart .and. &
                        field1(1) /= avg%missval ) then
                      field1 = field1 * var%factor
                   end if
                   call streamWriteVarSlice( file%streamID, avg%varID, level-1, &
                        field1, count(field1 == avg%missval) )
                end if

             else! var%surface_integrate

                ! Actual 2d or 3d values are also written by level.
                ! The slice getter takes care of masking and scaling.

                CALL gather( &
                     iolist_element_get_slice( the_file_list, current_avglist, &
                                               var, level ), &
                     field_g, p_io)

#ifdef _PROFILE
                CALL trace_start ('streamwritevar', 14)
#endif
                IF (p_pe == p_io) THEN
                   CALL streamWriteVarSlice(file%streamID, &
                        avg%varID, level-1, field_g, &
                        COUNT(field_g == avg%missval))
                END IF
#ifdef _PROFILE
                CALL trace_stop ('streamwritevar', 14)
#endif

             end if!else var%surface_integrate

          END DO

       END IF

       current_avglist = average_list_tail(current_avglist)
    END DO ! file%avglist

#endif/*ndef NOCDI */

    ! For restart files we only ever write timestep 0,
    ! so we do not increment the number of steps.
    IF(file%type /= file_restart) &
         file%nstep = file%nstep + 1

  END SUBROUTINE iolist_element_write


  !>
  !! Write appropriate data into write-once-at-begin files.
  !!
  subroutine iolist_write_initial(this, the_varlist, model_time)
    type(file_list), intent(in) :: this
    type(varlist), POINTERINTENT(in) :: the_varlist
    type(time_desc), intent(in) :: model_time
    call iolist_write_once( this, the_varlist, model_time, &
                            (/ file_initial /) )
  end subroutine iolist_write_initial

  !>
  !! Write appropriate data into restart and write-once-at-end files.
  !!
  subroutine iolist_write_final(this, the_varlist, model_time)
    type(file_list), intent(in) :: this
    type(varlist), POINTERINTENT(in) :: the_varlist
    type(time_desc), intent(in) :: model_time
    call iolist_write_once( this, the_varlist, model_time, &
                            (/ file_restart, file_final /) )
  end subroutine iolist_write_final

  !>
  !! Write appropriate data into write-once files (restart, final, initial)
  !!
  !! These files are not pre-opened, so stream operations are done here.
  !! @see create_avg_stream
  !!
  SUBROUTINE iolist_write_once(this, the_varlist, model_time, file_types)

    TYPE (file_list), INTENT(in) :: this
    TYPE (varlist), POINTERINTENT(in) :: the_varlist
    type(time_desc), intent(in) :: model_time
    integer, intent(in) :: file_types(:)

#ifndef NOCDI
    TYPE (file_list) :: current
    INTEGER :: file_type

    current = this
    DO WHILE(file_list_is_not_empty(current))

       IF (any(current%head%type == file_types)) THEN
          ! Only handle files of one of the specified 'write-once' types,
          ! averaged fields are handled by iolist_poll_write_avglist.

          IF (p_pe == p_io) THEN

             write (0, '(A)') "write-once file: " // &
                  trim(file_item_to_string(current%head))

             SELECT CASE (tolower(current%head%format))
             CASE ('ext')
               file_type = filetype_ext
             CASE ('nc')
               file_type = filetype_nc
             CASE ('nc2')
               file_type = filetype_nc2
             CASE ('nc4')
               file_type = filetype_nc4
             CASE default
               CALL stop_all('file format in ioctl should be ext, nc, nc2, nc4')
             END SELECT

             current%head%streamID &
                  = streamOpenWrite(TRIM(current%head%name), file_type)
             IF (current%head%streamID == cdi_elibnavail &
                  .AND. file_type == filetype_nc4) THEN
               WRITE(0, *) "netCDF 4 output unavailable, retrying with netCDF 2"
               file_type = filetype_nc2
               current%head%streamID = &
                    streamOpenWrite(TRIM(current%head%name), file_type)
             END IF

             IF ( current%head%streamID < 0 ) THEN
                WRITE(0,*) cdiStringError(current%head%streamID)
                WRITE(0,*) TRIM(current%head%name),' open failed '
                CALL stop_all('in writeRestartCDI: open failed ')
             ENDIF

             CALL streamdefvlist(current%head%streamID, current%head%vlistID)

          ENDIF

          CALL iolist_element_write(current, the_varlist, model_time)

          IF (p_pe == p_io) CALL streamClose(current%head%streamID)

       END IF  ! type in file_types

       current = file_list_tail(current)
    END DO ! this
#endif
  END SUBROUTINE iolist_write_once
  !>
  !! Read data from restart files.
  !!
  !! Restart files are not pre-opened, so stream operations are done here.
  !! @see create_avg_stream
  !!
  SUBROUTINE readRestartCDI(this,var)

    TYPE (file_list), INTENT(in) :: this
    TYPE (varlist), POINTERINTENT(in) :: var

    TYPE (file_list) :: current
    TYPE (varlist_element)  :: bufferVar
    TYPE (average_list) :: currentAvg


#ifndef NOCDI
    REAL(wp), ALLOCATABLE :: field(:,:)

    INTEGER :: k,LevelID,idate,itime,status,tsID,vlistID,taxisID
    INTEGER :: VarID,code,zaxisID,nlevel,nmiss
    ! INTEGER :: nvars
    ! REAL(wp) :: rmiss

    nmiss=0

    !        IF ( p_pe == p_io ) CALL print_varlist(var)

    current = this
    DO WHILE ( file_list_is_not_empty(current) )  ! exit if null pointer

       IF (current%head%type == file_restart) THEN

          IF ( p_pe == p_io ) THEN

             write (0, '(A)') "open/read restart file: " // &
                  trim(file_item_to_string(current%head, .true.))

             current%head%streamID = streamopenread(TRIM(current%head%name))

             IF ( current%head%streamID < 0 ) THEN
                WRITE(0,*) cdiStringError(current%head%streamID)
                WRITE(0,*) TRIM(current%head%name),' open failed '
                CALL stop_all('in readRestartCDI : open failed  ')
             ENDIF

          ENDIF

          ALLOCATE(field(ie,je))

          IF (p_pe == p_io) THEN

             ! Get the variable list of the dataset
             vlistID = streamInqVlist(current%head%streamID)
             ! Get the Time axis form the variable list
             taxisID = vlistInqTaxis(vlistID)

             tsID=0
             status = streamInqTimestep(current%head%streamID, tsID)

             ! Get the verification date and time
             idate = taxisInqVdate(taxisID)
             itime = taxisInqVtime(taxisID)


             !              WRITE(0,*) 'Timestep: ', idate, itime,lyears,lmonts,ldays

             !              nvars = vlistNvars(vlistID)

          ENDIF

          !            CALL p_bcast(nvars,p_io)
          CALL p_bcast(idate,p_io)
          CALL p_bcast(itime,p_io)

          lyears=INT(idate/10000_c_int64_t)
          lmonts=INT((idate-INT(lyears*10000,c_int64_t))/100_c_int64_t)
          ldays=INT(idate - INT(lyears*10000, c_int64_t) - INT(lmonts*100, c_int64_t))



          currentAvg = current%head%avglist
          DO WHILE(average_list_is_not_empty(currentAvg))

             code = currentAvg%head%code
             bufferVar = get_var_by_code(var, code)

             IF (current%on_io_node) THEN

                varID = vlistInqVarID(vlistID,code)
                !                rmiss = vlistInqVarMissval(vlistID, varID)
                !                WRITE(0,*)'bufferVar',code,varid
                IF(varID < 0 .and. code /= 999) THEN
                   WRITE(0,*)'Warning'
                   WRITE(0,*)'CODE: ' ,currentAvg%head%code,' (', TRIM(bufferVar%name),   &
                        ') selected (ioctl) but not found in ',TRIM(current%head%name)
                END IF

             ENDIF

             ! Broadcast variable ID to all processes.
             ! So every process can check if read/scatter needs to be done.
             CALL p_bcast(varID, p_io)

             ! Skip variables that were not found in file.
             IF ( varID >= 0 ) THEN

                IF (ASSOCIATED(bufferVar%array_2d)) THEN
#ifdef _PROFILE
                   CALL trace_start ('streamreadvar', 15)
#endif
                   IF (current%on_io_node) THEN
                      CALL streamreadvar(current%head%streamID,varID,field_g,nmiss)
                   ENDIF
#ifdef _PROFILE
                   CALL trace_stop ('streamreadvar', 15)
                   CALL trace_start ('streamreadvar_scatter', 16)
#endif
                   CALL scatter(field_g,field,p_io)
                   bufferVar%array_2d=field
#ifdef _PROFILE
                   CALL trace_stop ('streamreadvar_scatter', 16)
#endif
                ENDIF ! 2D buffer associated


                IF (ASSOCIATED(bufferVar%array_3d)) THEN

                   IF (current%on_io_node) THEN
                      zaxisID = vlistInqVarZaxis(vlistID, varID)
                      nlevel = zaxisInqSize(zaxisID)
                   ENDIF

                   CALL p_bcast(nlevel,p_io)

                   DO k=1,nlevel

#ifdef _PROFILE
                      CALL trace_start ('streamreadvarslice', 17)
#endif
                      IF (current%on_io_node) THEN
                         levelID=k-1
                         CALL streamreadvarslice(current%head%streamID,varID,levelID,field_g,nmiss)
                      ENDIF
#ifdef _PROFILE
                      CALL trace_stop ('streamreadvarslice', 17)
                      CALL trace_start ('streamreadvarslice_scatter', 18)
#endif
                      CALL scatter(field_g,field,p_io)
#ifdef _PROFILE
                      CALL trace_stop ('streamreadvarslice_scatter', 18)
                      CALL trace_start ('streamreadvarslice_copy', 19)
#endif
                      bufferVar%array_3d(:,:,k)=field
#ifdef _PROFILE
                      CALL trace_stop ('streamreadvarslice_copy', 19)
#endif
                   ENDDO ! k

                ENDIF ! 3D buffer associated

             ENDIF ! varID >= 0

             currentAvg = average_list_tail(currentAvg)
          ENDDO


          IF (current%on_io_node) CALL streamclose(current%head%streamID)

       END IF ! type == file_restart

       current = file_list_tail(current)
    END DO
#endif
  END SUBROUTINE readRestartCDI

  !>
  !! Accumulate data into averaging sub-system at staggered time step.
  !!
  SUBROUTINE iolist_accumulate_staggered(this,var,model_time)

    TYPE (file_list), INTENT(in) :: this
    TYPE (varlist), POINTERINTENT(in) :: var
    type(time_desc), intent(in) :: model_time

    CALL iolist_add_avg_data(this,var,model_time,staggered=.true.)

  END SUBROUTINE iolist_accumulate_staggered

  !>
  !! Accumulate data into averaging sub-system at full time step.
  !!
  SUBROUTINE iolist_accumulate(this,var,model_time)

    TYPE (file_list), INTENT(in) :: this
    TYPE (varlist), POINTERINTENT(in) :: var
    type(time_desc), intent(in) :: model_time

    CALL iolist_add_avg_data(this,var,model_time,staggered=.false.)

  END SUBROUTINE iolist_accumulate

  !>
  !! Handle accumulation, evaluting include or exclude list for codes.
  !!
  SUBROUTINE iolist_add_avg_data(this, the_varlist, model_time, staggered)

    TYPE (file_list), INTENT(in) :: this
    TYPE (varlist), POINTERINTENT(in) :: the_varlist
    type(time_desc), intent(in) :: model_time
    logical, intent(in) :: staggered

    TYPE (file_list) :: current_iolist
    TYPE (average_list) :: current_avglist
    TYPE (varlist_element) :: var
    LOGICAL :: todo
    INTEGER :: step

    ! Initialize iolist loop
    current_iolist = this
    DO WHILE ( file_list_is_not_empty(current_iolist) )

       step = get_step_of_averaging(current_iolist%head%type, model_time)

       ! Initialize avglist loop
       current_avglist = current_iolist%head%avglist
       DO WHILE ( average_list_is_not_empty(current_avglist) )

          var = get_var_by_code(the_varlist, current_avglist%head%code)

          ! Only deal with variables at their appropriate stage
          IF (staggered .eqv. var%staggered) THEN

             ! For initialization and accumulation, skip snapshot entries.
             ! Snapshots of integration variables are computed at write time.
             IF(.not. current_iolist%head%snapshot) THEN
                IF(step == 1) THEN
                   CALL avglist_element_start(current_avglist)
                END IF
                CALL avglist_element_accumulate(current_avglist, var, step)
             END IF

          END IF

          ! Increment avglist loop
          current_avglist = average_list_tail(current_avglist)
       END DO

       ! Increment iolist loop
       current_iolist = file_list_tail(current_iolist)
    END DO

  END SUBROUTINE iolist_add_avg_data


  !>
  !! Prepare field for averaging item.
  !!
  SUBROUTINE avglist_element_start(this)
    TYPE (average_list), INTENT(inout) :: this
    this%head%field = 0._dp
  END SUBROUTINE avglist_element_start

  !>
  !! Actually accumulate values into appropriate averaging field.
  !!
  SUBROUTINE avglist_element_accumulate(this, var, step)

    TYPE (average_list), INTENT(inout) :: this
    TYPE (varlist_element), INTENT(in) :: var
    INTEGER, INTENT(in) :: step

    INTEGER :: top, bottom
    REAL(dp) :: stepi

    stepi = 1._dp/REAL(step, dp)    

    IF (ASSOCIATED (var%scalar)) THEN

       this%head%field(1,1,1) = this%head%field(1,1,1) &
            + (var%scalar - this%head%field(1,1,1))*stepi

    ELSE IF (ASSOCIATED (var%array_1d)) THEN

       this%head%field(:,1,1) = this%head%field(:,1,1) &
            + (var%array_1d(:) - this%head%field(:,1,1))*stepi

    ELSE IF (ASSOCIATED (var%array_2d)) THEN

       if(var%surface_integrate) then
          this%head%field(1,1,1) = this%head%field(1,1,1) &
               + ( surface_integral(var%array_2d(:,:)) - &
                   this%head%field(1,1,1) )*stepi
       else
          this%head%field(:,:,1) = this%head%field(:,:,1) &
               + (var%array_2d(:,:) - this%head%field(:,:,1))*stepi
       end if

    ELSE IF (ASSOCIATED (var%array_3d)) THEN

       top = merge(var%top_level_index, 1, var%top_level_index > 0)
       bottom = merge(var%level_index, size(var%array_3d, 3), var%level_index > 0)

       if(var%column_integrate .and. var%surface_integrate) then
          this%head%field(1,1,1) = this%head%field(1,1,1) &
               + ( surface_integral(column_integral_field(var%array_3d(:,:,top:bottom))) - &
                   this%head%field(1,1,1) )*stepi
       else if(var%column_integrate .and. .not. var%surface_integrate) then
          this%head%field(:,:,1) = this%head%field(:,:,1) &
               + ( column_integral_field(var%array_3d(:,:,top:bottom)) - &
                   this%head%field(:,:,1) )*stepi
       else if(.not. var%column_integrate .and. var%surface_integrate) then
          this%head%field(1,1,:) = this%head%field(1,1,:) &
               + ( level_integral_vector(var%array_3d(:,:,top:bottom)) - &
                   this%head%field(1,1,:) )*stepi
       else if(.not. var%column_integrate .and. .not. var%surface_integrate) then
          this%head%field(:,:,:) = this%head%field(:,:,:) &
               + ( var%array_3d(:,:,top:bottom) - this%head%field(:,:,:) )*stepi
       end if

    END IF

  END SUBROUTINE avglist_element_accumulate

  !>
  !! Get rid of undefined data points according to grid settings.
  !!
  subroutine mask_field(field, grid, missval, level)
    real(dp), intent(inout) :: field(:,:)
    character(*), intent(in) :: grid
    real(dp), intent(in) :: missval
    integer, intent(in) :: level

    select case(grid)
    case(grid_p, grid_p_minus, grid_s, grid_s_minus)

       if(level <= ke) then
          where(weto(:,:,level) <= 0.5_dp)
             field = missval
          end where
       else
          field = missval
       endif

    case(grid_u, grid_u_plus)

       where(amsuo(:,:,level) <= 0.5_dp)
          field = missval
       end where

    case(grid_v, grid_v_plus)

       where(amsue(:,:,level) <= 0.5_dp)
          field = missval
       end where

    end select

  end subroutine mask_field

  !>
  !! Get appropriate data for output, masked according to grid settings.
  !!
  function iolist_element_get_slice(file, file_var, var, level)
     type(file_list), intent(in) :: file
     type(average_list), intent(in) :: file_var
     type(varlist_element), intent(in) :: var
     integer, intent(in) :: level
     real(dp) :: iolist_element_get_slice(var%i_shape(1), var%i_shape(2))

     integer :: src_level

     ! Check preconditions.
     ! If not a 2d or 3d field, return undefined value.
     if(.not. (associated(var%array_2d) .or. associated(var%array_3d))) then
        iolist_element_get_slice = file_var%head%missval
        return
     end if

     src_level = &
        merge(var%top_level_index - 1, 0, var%top_level_index > 0) + level

     ! Select appropriate source field.
     if(file%head%snapshot .and. .not. var%column_integrate) then
        if(associated(var%array_2d)) then
           iolist_element_get_slice = var%array_2d
        else
           iolist_element_get_slice = var%array_3d(:,:,src_level)
        end if
     else
        iolist_element_get_slice = file_var%head%field(:,:,level)
     end if

     ! Note that restart files are unmasked by default.
     if(.not. file%head%unmasked) then
        ! Use prescribed level for masking if available.
        ! Otherwise, take level from input.
        call mask_field(iolist_element_get_slice, var%grid, &
             file_var%head%missval, &
             merge(var%mask_level_index, src_level, var%mask_level_index > 0))
     end if

     ! Scaling factor is only applied for non-restart files.
     if(file%head%type /= file_restart) then
        if(var%factor /= 1.0_dp) then
           where(iolist_element_get_slice /= file_var%head%missval)
              iolist_element_get_slice = iolist_element_get_slice * var%factor
           end where
        end if
     end if

  end function iolist_element_get_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Time and date related utilities
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !>
  !! Get current date in I/O format.
  !!
  FUNCTION get_date(model_time)
    type(time_desc), intent(in) :: model_time
    INTEGER :: get_date
    get_date = (model_time%year*100 + model_time%month)*100 + model_time%mday
  END FUNCTION get_date

  !>
  !! Get current time in I/O format.
  !!
  FUNCTION get_time(model_time)
    type(time_desc), intent(in) :: model_time
    INTEGER :: get_time
    get_time = (model_time%hour*100 + model_time%minute)*100 + model_time%second
  END FUNCTION get_time

  !>
  !! Inquire current step in averaging interval.
  !!
  function get_step_of_averaging(file_type, model_time) result(step)

    integer, intent(in) :: file_type
    type(time_desc), intent(in) :: model_time
    integer :: step

    step = -1
    select case(file_type)
    case(file_annual);   step = get_step_of_year(model_time)
    case(file_monthly);  step = get_step_of_month(model_time)
    case(file_daily);    step = get_step_of_day(model_time)
    case(file_timestep); step = 1
    case(file_12h);      step = get_step_of_partial_day(model_time, 2)
    case(file_6h);       step = get_step_of_partial_day(model_time, 4)
    case(file_3h);       step = get_step_of_partial_day(model_time, 8)
    case(file_2h);       step = get_step_of_partial_day(model_time, 12)
    case(file_1h);       step = get_step_of_partial_day(model_time, 24)
    end select

  end function

  !>
  !! Check if last step of averaging interval is reached.
  !!
  function is_end_of_averaging(file_type, model_time) result(is_end)

    integer, intent(in) :: file_type
    type(time_desc), intent(in) :: model_time
    logical :: is_end

    is_end = .false.
    select case(file_type)
    case(file_annual);   is_end = is_end_of_year(model_time)
    case(file_monthly);  is_end = is_end_of_month(model_time)
    case(file_daily);    is_end = is_end_of_day(model_time)
    case(file_timestep); is_end = .true.
    case(file_12h);      is_end = is_end_of_partial_day(model_time, 2)
    case(file_6h);       is_end = is_end_of_partial_day(model_time, 4)
    case(file_3h);       is_end = is_end_of_partial_day(model_time, 8)
    case(file_2h);       is_end = is_end_of_partial_day(model_time, 12)
    case(file_1h);       is_end = is_end_of_partial_day(model_time, 24)
    end select

  end function

#endif/*ndef NO_NEW_IO */

END MODULE mo_iolist
