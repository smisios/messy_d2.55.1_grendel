!-----------------------------------------------------------------------
! BOP
!
! !MODULE:  mo_couple.F90
!
MODULE mo_couple
  !
  ! !DESCRIPTION:
  !
  !   contains all additional code used for coupling with OASIS3
  !
  ! !USES:
  !
#if defined __coupled || __prism

#ifdef FLUSH_NEEDS_F90_UNIX_IO
  USE f90_unix_io, ONLY: flush
#endif

  USE mo_param1
  USE mo_units,           ONLY : io_in_octl, io_stdout, io_stderr
  USE mo_commo1

#ifdef __cpl_dust
  USE mo_carbch,          ONLY : dustdep
#endif
#ifdef __cpl_dms
  USE mo_carbch,          ONLY : dmsflux, dms
#endif
#ifdef __cpl_co2
  USE mo_carbch,          ONLY : co2conc, co2flux_cpl, co2trans, suppco2, atm
  USE mo_param1_bgc,      ONLY : iatmco2
#endif

  USE mo_commoau1,        ONLY : tmelt
  USE mo_planetary_constants, ONLY : rhoref_snow, rhoref_water

  USE mo_fluxes1
#ifdef FLUXCORRECT
  USE mo_commo_fluxcorr,  ONLY : fluko_hfl,fluko_frw
#endif
  USE mo_mpi
  USE mo_parallel
  USE mo_boundsexch,      ONLY : bounds_exch

  USE mo_kind

  USE mod_kinds_model,    ONLY: ip_realwp_p
  USE mod_prism_proto
  USE mod_comprism_proto
  USE mod_prism_get_proto
  USE mod_prism_put_proto
  USE mod_prism_def_partition_proto, ONLY: prism_def_partition_proto

  USE mo_rotation

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER   :: pwp = ip_realwp_p   ! working precision (psmile)

  INTEGER, ALLOCATABLE :: paral(:) ! parallel strategy
  INTEGER :: nsegments             ! no. of segments
  INTEGER :: parsize               !
  INTEGER :: nbtotproc             ! total no. of procs
  INTEGER :: nbcplproc = 1         ! no. of procs involved in data exchange
  INTEGER :: commlocal             ! local communicator
  INTEGER :: comp_id               ! model id
  INTEGER :: rank                  ! rank of processor in local communicator
  INTEGER :: info                  ! psmile info message coding
  INTEGER :: ierror                ! error codes
  INTEGER :: var_shape(2)          !
  INTEGER :: part_id               !
  INTEGER :: var_nodims(2)         !

  CHARACTER (len=6) :: modnam = 'mpiom' ! model name
#ifndef __cpl_co2
#ifdef __cpl_dust
  INTEGER, PARAMETER :: nflda2o = 12 ! no. of fields passed from atmosphere to ocean
#else
  INTEGER, PARAMETER :: nflda2o = 11 ! no. of fields passed from atmosphere to ocean
#endif
#ifdef __cpl_dms
  INTEGER, PARAMETER :: nfldo2a = 8  ! no. of fields passed from ocean to atmosphere
#else
  INTEGER, PARAMETER :: nfldo2a = 6  ! no. of fields passed from ocean to atmosphere
#endif
#else /*__cpl_co2*/
#ifdef __cpl_dust
  INTEGER, PARAMETER :: nflda2o = 14 ! no. of fields passed from atmosphere to ocean
#else
  INTEGER, PARAMETER :: nflda2o = 13 ! no. of fields passed from atmosphere to ocean
#endif
#ifdef __cpl_dms
  INTEGER, PARAMETER :: nfldo2a = 10 ! no. of fields passed from ocean to atmosphere
#else
  INTEGER, PARAMETER :: nfldo2a = 8  ! no. of fields passed from ocean to atmosphere
#endif
#endif /*__cpl_co2*/

  CHARACTER (len=8) :: clstr8a2o(nflda2o)   ! port names of exchange fields received
  CHARACTER (len=8) :: clstr8o2a(nfldo2a)   ! port names of exchange fields sent

  INTEGER    :: a2o_freq             ! exchange frequency in secs for receiving fields
  INTEGER    :: o2a_freq             ! exchange frequency in secs for sending fields
  REAL(pwp)  :: o2a_time = 0.0_pwp   ! time passed since last couple_put_o2a

  INTEGER    :: run_date_secs = 0    ! time (sec) passed since start of run

  ! should be INTEGER (kind=ip_intwp_p) rather than INTEGER*4
  INTEGER(i4)  :: portin_id (nflda2o)  ! Port IDs of exchange fields received
  INTEGER(i4)  :: portout_id(nfldo2a)  ! Port IDs of exchange fields sent

  INTEGER    :: nmseq                ! run mode : sequential/concurrent = 2/1
  INTEGER    :: num_fld_recvd = 0    ! counter for received fields
  INTEGER    :: num_fld_sent  = 0    ! counter for fields sent

  INTEGER    :: jfld

  REAL(wp), POINTER :: amsue_g_l1(:,:)
  REAL(wp), POINTER :: amsuo_g_l1(:,:)

  INTEGER       :: jpdim_ocei
  INTEGER       :: jpdim_ocej
  INTEGER       :: jpdim_oce

  ! Global fields

  REAL(wp), POINTER :: gl_sst    (:,:)
  REAL(wp), POINTER :: gl_sicomo (:,:)
  REAL(wp), POINTER :: gl_sicsno (:,:)
  REAL(wp), POINTER :: gl_sictho (:,:)
  REAL(wp), POINTER :: gl_socu   (:,:)
  REAL(wp), POINTER :: gl_socv   (:,:)
  REAL(wp), POINTER :: gl_sicu   (:,:)
  REAL(wp), POINTER :: gl_sicv   (:,:)
#ifdef __cpl_dms
  REAL(wp), POINTER :: gl_dmsflux(:,:)
  REAL(wp), POINTER :: gl_dms    (:,:)
#endif
#ifdef __cpl_co2
  REAL(wp), POINTER :: gl_co2trans(:,:)
  REAL(wp), POINTER :: gl_suppco2 (:,:)
#endif

  ! Accumulated global fields

  INTEGER, SAVE :: nacc

  REAL(wp), POINTER :: gl_sstacc (:,:)
  REAL(wp), POINTER :: gl_sitoacc(:,:)
  REAL(wp), POINTER :: gl_sicoacc(:,:)
  REAL(wp), POINTER :: gl_sntoacc(:,:)
  REAL(wp), POINTER :: gl_socuacc(:,:)
  REAL(wp), POINTER :: gl_socvacc(:,:)
#ifdef __cpl_dms
  REAL(wp), POINTER :: gl_dmsfacc(:,:)
  REAL(wp), POINTER :: gl_dmsacc (:,:)
#endif
#ifdef __cpl_co2
  REAL(wp), POINTER :: gl_co2tracc(:,:)
  REAL(wp), POINTER :: gl_co2suacc(:,:)
#endif

#ifndef __accu_by_psmile

  ! Accumulated local fields

  REAL(wp), POINTER :: sstacc (:,:)
  REAL(wp), POINTER :: sitoacc(:,:)
  REAL(wp), POINTER :: sicoacc(:,:)
  REAL(wp), POINTER :: sntoacc(:,:)
  REAL(wp), POINTER :: socuacc(:,:)
  REAL(wp), POINTER :: socvacc(:,:)
  REAL(wp), POINTER :: sntotmp(:,:)
  REAL(wp), POINTER :: sitotmp(:,:)
  REAL(wp), POINTER :: sicotmp(:,:)
#ifdef __cpl_dms
  REAL(wp), POINTER :: dmsfacc(:,:)
  REAL(wp), POINTER :: dmsacc (:,:)
#endif
#ifdef __cpl_co2
  REAL(wp), POINTER :: co2tracc(:,:)
  REAL(wp), POINTER :: co2suacc(:,:)
#endif

#endif

  ! Temporary global array without halos

  REAL(wp), POINTER :: pfield_g(:,:)

#ifdef FLUXCORRECT
  REAL(wp)          ::  fhflmax, fhflmin
#endif

  PUBLIC :: couple_prep, couple_correct_ini, couple_init,    &
       couple_end,  couple_put_o2a,     couple_get_a2o, &
       couple_calendar, &
       o2a_freq, o2a_time
  !
  ! !DESCRIPTION:
  !
  ! - contains all code related to coupling with oasis3
  ! - Note: this interface to the psmile library uses variables from
  !   the psmile modules mod_comprism_proto.
  !   In the next version (oasis4) these variables will all be received
  !   by a call of a prism_get_... routine. Then an update of the code will
  !   be needed. These variables are recognized by their prefix 'ig_'.
  !
  ! !REVISION HISTORY:
  ! 03.05.22 Stephanie Legutke, MPI-Met, M&D
  !          - created
  ! 03.08.07 Stephanie Legutke, MPI-Met, M&D
  !          - update psmile
  ! 04.10.01 Noel Keenlyside, IFM-GEOMAR
  !          - message passing version (accumulation done on global arrays)
  ! 09.12.12 Rene Redler, MPI-Met
  !          - support for different OASIS3 restart file names
  !          - Option to write OASIS3 restart files in NetCDF
  ! 10.03.09 Rene Redler, Helmuth Haak, MPI-Met
  !          - accumulation of field as local operations to avoid
  !            unneccessary and exhaustive collective communication (work in progress)
  !          - loops changed such that innermost index is treated in inner loop
  !          - unneccessay broadcasts removed
  !          - pfield_g is only allocated once in the initialisation phase
  !            rather than at each exchange
  !          - modules are only loaded once
  !          - eliminated gathering of fields that are not used anywhere
  !          - printout behind __synout followed by call flush
  !          - additional __synout output added for debugging purposes
  !          - bug fix in prep_o2a (coupling with accumulation in psmile) for
  !            tripolar grids (Note: This fix changes results.)
  ! EOP
  !-----------------------------------------------------------------------
  ! $Id: mo_couple.f90,v 1.4.2.1.10.1.4.2.2.3.4.1 2006/03/31 08:20:35 m211054 Exp $
  !-----------------------------------------------------------------------
  !

  CHARACTER(len=*), PARAMETER :: cl_version_string = &
       '$Id: mo_couple.f90,v 1.4.2.1.10.1.4.2.2.3.4.1 2006/03/31 08:20:35 m211054 Exp $'

CONTAINS

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  couple_prep
  !
  ! !INTERFACE:
  !
  SUBROUTINE couple_prep
    !
    ! !DESCRIPTION:
    !
    ! - prepare coupling.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Standard output file oceout is connected to unit ',io_stdout
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    ENDIF
#endif

  END SUBROUTINE couple_prep

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  couple_correct_ini
  !
  ! !INTERFACE:
  !
  SUBROUTINE couple_correct_ini
    !
    ! !DESCRIPTION:
    !
    !- flux correction initialization.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

#ifdef FLUXCORRECT
    CALL FLUXC_INI
#endif /*FLUXCORRECT*/

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Flux correction initialized'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    ENDIF
#endif

  END SUBROUTINE couple_correct_ini


  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  couple_init
  !
  ! !INTERFACE:
  !
  SUBROUTINE couple_init(actyear)
    !
    ! !INPUT/OUTPUT VALUE:
    !
    INTEGER, INTENT(inout) :: actyear  ! year at start of run from restart file (IN)
    !                                    year at start of run from oasis (OUT)!
    ! !DESCRIPTION:
    !
    !- Initialize communication.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------
    !    WRITE(0,*) ' couple_init 1',p_pe
#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Start couple_init'
      WRITE(io_stdout,*) ' *****************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    ENDIF
#endif

    ! Global dimensions without halos

    jpdim_ocei = ie_g-2
    jpdim_ocej = je_g

    IF (lbounds_exch_tp ) jpdim_ocej = je_g-2

    jpdim_oce  = jpdim_ocej * jpdim_ocei

    ! Allocate global memory for arrays used in coupling

    IF (p_parallel_io) THEN

      ALLOCATE(gl_sstacc (jpdim_ocei,jpdim_ocej), &
               gl_sitoacc(jpdim_ocei,jpdim_ocej), &
               gl_sicoacc(jpdim_ocei,jpdim_ocej), &
               gl_sntoacc(jpdim_ocei,jpdim_ocej))
#ifdef __cpl_dms
      ALLOCATE(gl_dmsfacc(jpdim_ocei,jpdim_ocej), &
               gl_dmsacc (jpdim_ocei,jpdim_ocej))
#endif
#ifdef __cpl_co2
      ALLOCATE(gl_co2tracc(jpdim_ocei,jpdim_ocej), &
               gl_co2suacc(jpdim_ocei,jpdim_ocej))
#endif
      ALLOCATE(amsue_g_l1(ie_g,je_g), &
               amsuo_g_l1(ie_g,je_g))

      ALLOCATE(gl_sst    (ie_g,je_g),  &
               gl_sicomo (ie_g,je_g),  &
               gl_sicsno (ie_g,je_g),  &
               gl_sictho (ie_g,je_g))

      ALLOCATE(gl_socuacc(jpdim_ocei,jpdim_ocej),  &
               gl_socvacc(jpdim_ocei,jpdim_ocej))

      ALLOCATE(gl_socu(ie_g,je_g), &
               gl_socv(ie_g,je_g), &
               gl_sicu(ie_g,je_g), &
               gl_sicv(ie_g,je_g))

#ifdef __cpl_dms
      ALLOCATE (gl_dmsflux(ie_g,je_g))
      ALLOCATE (gl_dms    (ie_g,je_g))
#endif
#ifdef __cpl_co2
      ALLOCATE (gl_co2trans(ie_g,je_g))
      ALLOCATE (gl_suppco2 (ie_g,je_g))
#endif
    ELSE
      gl_sstacc   => NULL()
      gl_sitoacc  => NULL()
      gl_sicoacc  => NULL()
      gl_sntoacc  => NULL()
      gl_socuacc  => NULL()
      gl_socvacc  => NULL()

      amsue_g_l1  => NULL()
      amsuo_g_l1  => NULL()

      gl_sst      => NULL()
      gl_sicomo   => NULL()
      gl_sicsno   => NULL()
      gl_sictho   => NULL()

      gl_socu     => NULL()
      gl_socv     => NULL()
      gl_sicu     => NULL()
      gl_sicv     => NULL()
#ifdef __cpl_dms
      gl_dmsfacc  => NULL()
      gl_dmsacc   => NULL()
      gl_dmsflux  => NULL()
      gl_dms      => NULL()
#endif
#ifdef __cpl_co2
      gl_co2tracc => NULL()
      gl_co2suacc => NULL()
      gl_co2trans => NULL()
      gl_suppco2  => NULL()
#endif
    ENDIF

#ifndef __accu_by_psmile
    ALLOCATE (sstacc (ie,je))
    ALLOCATE (sitoacc(ie,je))
    ALLOCATE (sicoacc(ie,je))
    ALLOCATE (sntoacc(ie,je))
    ALLOCATE (socuacc(ie,je))
    ALLOCATE (socvacc(ie,je))
#ifdef __cpl_dms
    ALLOCATE (dmsfacc(ie,je))
    ALLOCATE (dmsacc (ie,je))
#endif
#ifdef __cpl_co2
    ALLOCATE (co2tracc(ie,je))
    ALLOCATE (co2suacc(ie,je))
#endif
    ALLOCATE (sitotmp(ie,je))
    ALLOCATE (sicotmp(ie,je))
    ALLOCATE (sntotmp(ie,je))
#endif

    ! temporary memory for gathering and scattering

    IF (p_parallel_io) THEN
      ALLOCATE (pfield_g(ie_g,je_g))
    ELSE
      pfield_g => NULL()
    END IF

    ! write(0,*) ' couple_init 2',p_pe

    CALL gather(amsuo(:,:,1),amsuo_g_l1,p_io)
    CALL gather(amsue(:,:,1),amsue_g_l1,p_io)

#ifdef __cpl_dms
    CALL gather(dmsflux,  gl_dmsflux, p_io)
    CALL gather(dms,      gl_dms,     p_io)
#endif
#ifdef __cpl_co2
    CALL gather(co2trans, gl_co2trans, p_io)
    CALL gather(suppco2,  gl_suppco2, p_io)
#endif

#if ! (defined (use_comm_MPI1) ^ defined (use_comm_MPI2))
  "In coupled mode specification of either -Duse_comm_MP1 or -Duse_comm_MPI2 is required."
#endif
    !
    ! WRITE(0,*) ' couple_init 3',p_pe
#ifdef use_comm_MPI2
    !-- Join the communicator group.
    !
    ierror = PRISM_Ok
    CALL prism_init_comp_proto (comp_id, modnam, ierror)
    IF (ierror /= PRISM_Ok) THEN
      CALL couple_abort (modnam,'couple_init','pb init_comp',io_stderr)
    ENDIF
    !
    ! WRITE(0,*) ' couple_init 4',p_pe
#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' after prism_init_comp_proto ...'
      WRITE(io_stdout,*) ' *******************************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    ENDIF
#endif
#endif

    !
    !-- Write grids file for oasis (only one thread)
    !
    CALL grids_writing(io_stdout)
    ! WRITE(0,*) ' couple_init 5',p_pe
    !
    !-- receive local communicator
    !
#ifdef use_comm_MPI1
    commlocal = p_all_comm
#endif
#ifdef use_comm_MPI2
    ierror = PRISM_Ok
    CALL prism_get_localcomm_proto(commlocal, ierror)
    IF (ierror /= PRISM_Ok) THEN
      CALL couple_abort (modnam,'couple_init','pb get_localcomm',io_stderr)
    ENDIF
#endif
    nbtotproc = p_nprocs
    rank      = p_pe
    !
    ! Fields passed from the ocean to the atmosphere
    !
    clstr8o2a(1:nfldo2a) = cg_cnaminp(1:nfldo2a)
    !
    ! Fields passed from the atmosphere to the ocean
    !
    clstr8a2o(1:nflda2o) = cg_cnamout(nfldo2a+1:nfldo2a+nflda2o)
    !
    !-- Calendar year
    !
    actyear = ig_inidate(1)
    !
    !-- Create OASIS initial file (ocean surface conditions).
    !
    nmseq = MAXVAL (ig_clim_seq)
    ! WRITE(0,*) ' couple_init 5a',p_pe
    IF ( nmseq > 1 ) THEN
      WRITE( io_stdout,* ) 'Calling ini_wrte to create oasis3 restart files'
      CALL ini_wrte(io_in_octl,io_stdout)
    ENDIF
    !    WRITE(0,*) ' couple_init 5b',p_pe
    !WRITE(0,*) ' couple_init 5'
    !
    !-- Decomposition strategy
    !
    nsegments = 1
    parsize   = 3
    ALLOCATE(paral(parsize))
    paral ( clim_strategy ) = clim_serial
    paral ( clim_length   ) = jpdim_oce
    paral ( clim_offset   ) = 0

    IF (p_parallel_io) THEN
#ifdef __synout
      WRITE(io_stdout,*) ' calling prism_def_partition_proto ...'
      WRITE(io_stdout,*) ' *************************************'
      CALL flush(io_stdout)
#endif
      ierror   = PRISM_Ok
      CALL prism_def_partition_proto (part_id, paral, ierror)
      IF (ierror /= PRISM_Ok) THEN
        CALL couple_abort (modnam,'couple_init','pb def_partition',io_stderr)
      ENDIF

      !
      !-- Definitions for ports of incoming fields
      !
      var_nodims(1) = 1
      var_nodims(2) = 0
      var_shape(1)  = 1
      var_shape(2)  = paral (clim_length)

      DO jfld = 1,nflda2o

        ierror = PRISM_Ok
#ifdef __synout
        WRITE(io_stdout,*) ' calling prism_def_var_proto ...:', jfld
        WRITE(io_stdout,*) ' *************************************'
        CALL flush(io_stdout)
#endif
        CALL prism_def_var_proto &
             (portin_id(jfld), clstr8a2o(jfld), &
              part_id, var_nodims, PRISM_in,    &
              var_shape, PRISM_Real, ierror )
        IF (ierror /= PRISM_Ok) THEN
          WRITE (io_stdout, *) &
               'WARNING : Problem with import port ',clstr8a2o(jfld)
          WRITE (io_stdout, *) &
               'WARNING : Port ID returned is : ',portin_id(jfld)
          WRITE (io_stdout, *) '=======   Error code number = ',ierror
          CALL flush(io_stdout)
        ELSE
          WRITE (io_stdout,*)'   '
          WRITE (io_stdout,*) &
               ' couple_init: Import port ',clstr8a2o(jfld),' defined'
          WRITE (io_stdout,*) &
               ' couple_init: Port ID returned is : ',portin_id(jfld)
          WRITE (io_stdout,*)' With exchange freq.    ', &
               ig_def_freq(portin_id(jfld))
          CALL flush(io_stdout)
        ENDIF

      END DO

      !
      !-- Definitions for ports of outgoing fields
      !
      DO jfld = 1,nfldo2a

#ifdef __synout
        WRITE(io_stdout,*) ' calling prism_def_var_proto ...:', jfld
        WRITE(io_stdout,*) ' *************************************'
        CALL flush(io_stdout)
#endif
        ierror = PRISM_Ok
        CALL prism_def_var_proto ( portout_id(jfld), clstr8o2a(jfld), &
                                   part_id, var_nodims, PRISM_Out,    &
                                   var_shape, PRISM_Real, ierror )
        IF (ierror /= PRISM_Ok) THEN
          WRITE (io_stdout, *) &
               'WARNING : Problem with export port ',clstr8o2a(jfld)
          WRITE (io_stdout, *) &
               'WARNING : Port ID returned is : ',portout_id(jfld)
          WRITE (io_stdout, *) '=======   Error code number = ',ierror
          CALL flush(io_stdout)
        ELSE
          WRITE (io_stdout,*)'   '
          WRITE (io_stdout,*) &
               ' couple_init: Export port ',clstr8o2a(jfld),' defined'
          WRITE (io_stdout,*) &
               ' couple_init: Port ID returned is : ',portout_id(jfld)
          WRITE (io_stdout,*)' With exchange freq.    ', ig_def_freq(portout_id(jfld))
          CALL flush(io_stdout)
        ENDIF
      END DO
      WRITE (io_stdout, *) ' '

      ierror = PRISM_Ok
      CALL prism_enddef_proto(ierror)
      IF (ierror /= PRISM_Ok) THEN
        WRITE (io_stdout, *) 'WARNING : Problem with prism_enddef_proto'
        WRITE (io_stdout, *) '=======   Error code number = ',ierror
        CALL flush(io_stdout)
      ENDIF
      !
      !-- Check
      !   ...coupling-control parameters against those received from oasis.
      !
      CALL chck_par
      !
      !-- Display control parameters received from oasis3.
      !   ------------------------------------------------
      WRITE(io_stdout,*)' Couple_init:'
      WRITE(io_stdout,*)' -----------------------------------------------------'
      WRITE(io_stdout,*)' model name                = ',modnam
      WRITE(io_stdout,*)' model time step           = ',dt
      WRITE(io_stdout,*)' model restart-file year   = ',lyear1
      WRITE(io_stdout,*)' -----------------------------------------------------'
      WRITE(io_stdout,*)' '
      WRITE(io_stdout,*)' Parameters received from oasis3 :'
      WRITE(io_stdout,*)' mynum                         ', mynum
      WRITE(io_stdout,*)' mytid                         ', mytid
      WRITE(io_stdout,*)' model number                  ', comp_id
      WRITE(io_stdout,*)' ig_mynummod                   ', ig_mynummod
      WRITE(io_stdout,*)' Total no. of processors       ', nbtotproc
      WRITE(io_stdout,*)' No. of procs for data exchange', nbcplproc
      WRITE(io_stdout,*)' depomposition ID              ', part_id
      WRITE(io_stdout,*)' total time of the simulation  ', ig_ntime
      WRITE(io_stdout,*)' number of fields exchanged    ', ig_clim_nfield
      WRITE(io_stdout,*)' initial date                  ', ig_inidate
      WRITE(io_stdout,*)' Lag of exported fields        ', ig_def_lag
      WRITE(io_stdout,*)' coupling period of fields     ', ig_def_freq
      WRITE(io_stdout,*)' sequential index of fields    ', ig_def_seq
      WRITE(io_stdout,*)' ig_def_norstfile              ', ig_def_norstfile
      WRITE(io_stdout,*)' number of restart files       ', ig_nbr_rstfile
      WRITE(io_stdout,*)' name of restart files         ', cg_def_rstfile
      WRITE(io_stdout,*)' restart files in netcdf       ', lg_ncdfrst
      WRITE(io_stdout,*)' name of exchanged fields (in) ', cg_cnaminp
      WRITE(io_stdout,*)' name of exchanged fields (out)', cg_cnamout
      WRITE(io_stdout,*)' cg_ignout_field               ', cg_ignout_field
      WRITE(io_stdout,*)' any field going through Oasis?', lg_oasis_field
      WRITE(io_stdout,*)' seconds between 2 sent        ', o2a_freq
      WRITE(io_stdout,*)' seconds between 2 receive     ', a2o_freq
      WRITE(io_stdout,*)' calender run year             ', actyear
      IF(nmseq /= 1) THEN
        WRITE(io_stdout,*)' the models run sequentially !'
        WRITE(io_stdout,*)' no. of sequential fields: ',nmseq
      ELSE
        WRITE(io_stdout,*)' all fields have sequential no. 1 !'
        WRITE(io_stdout,*)' the models run concurrently !'
      ENDIF
#ifndef __accu_by_psmile
      WRITE (io_stdout,*) 'check namcouple: LOCTRANS = INSTANT is required.'
      WRITE (io_stdout,*) 'but mo_couple in mpiom still broken!'
      WRITE (io_stdout,*) 'activate  __accu_by_psmile'
      !     CALL couple_abort (modnam,'couple_init','unsupported code section',io_stderr)
#else
      WRITE (io_stdout,*) 'check namcouple: LOCTRANS = AVERAGE required in namcouple'
#endif
    ENDIF  ! p_parallel_io
    !
    !-- Initilize exchange fields
    !
    !   WRITE(0,*) ' couple_init 6',p_pe
    CALL ini_o2a(io_stdout)
    !   WRITE(0,*) ' couple_init 7',p_pe
    !
    CALL p_bcast(o2a_freq,p_io)
    !   CALL p_bcast(o2a_time,p_io)

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' End of couple_init'
      WRITE(io_stdout,*) ' ******************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    ENDIF
#endif

    CALL rotate2_ini

  END SUBROUTINE couple_init


  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  couple_get_a2o
  !
  ! !INTERFACE:
  !
  SUBROUTINE couple_get_a2o(ldtrun)
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(IN)   ::  ldtrun ! model run step
    !
    ! !LOCAL DECLARATIONS:
    !
    REAL(pwp)  :: exfld(jpdim_ocei,jpdim_ocej) ! buffer for receiving exchange fields
    !
    ! !DESCRIPTION:
    !
    !- Get data from coupler.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Start of couple_get_a2o'
      WRITE(io_stdout,*) ' ***********************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    ENDIF
#endif

    !
    !-- Import exchange fields; try at all dates
    !

    DO jfld = 1,nflda2o

      IF (p_parallel_io) THEN
#ifdef __synout
        WRITE(io_stdout,*)'    '
        WRITE(io_stdout,*)' ==> prism_get for ',clstr8a2o(jfld)
        WRITE(io_stdout,*) &
             ' ==> at run-step ',ldtrun,' (seconds passed=',run_date_secs,')'
        FLUSH(io_stdout)
#endif
        info = PRISM_Ok
        CALL prism_get_proto(portin_id(jfld),run_date_secs,exfld,info)

        IF ( info == PRISM_Recvd        .OR. &
             info == PRISM_FromRest     .OR. &
             info == PRISM_RecvOut      .OR. &
             info == PRISM_FromRestOut ) num_fld_recvd = num_fld_recvd + 1

        CALL digest_get_Id &
             (io_stdout,io_stderr,info,clstr8a2o(jfld),run_date_secs,num_fld_recvd)

      ENDIF
      CALL p_bcast(num_fld_recvd,p_io)
      CALL p_bcast(info, p_io)

      IF ( info == PRISM_Recvd        .OR. &
           info == PRISM_FromRest     .OR. &
           info == PRISM_RecvOut      .OR. &
           info == PRISM_FromRestOut )     &
           CALL PUT_A2O(clstr8a2o(jfld),exfld,jpdim_ocei,jpdim_ocej,io_stdout)
    ENDDO

    !
    !-- If all fields have been received ...postprocessing
    !

    IF ( num_fld_recvd == nflda2o )      &
         CALL post_a2o(num_fld_recvd, io_stdout, ldtrun)

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' End of couple_get_a2o'
      WRITE(io_stdout,*) ' *********************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    ENDIF
#endif

  END SUBROUTINE couple_get_a2o


  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  couple_put_o2a
  !
  ! !INTERFACE:
  !
  SUBROUTINE couple_put_o2a(ldtrun)
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(IN)   ::  ldtrun ! model run step
    !
    ! !LOCAL DECLARATIONS:
    !
    REAL(pwp)  :: exfld(jpdim_ocei,jpdim_ocej) ! buffer array for sending.
    !
    ! !DESCRIPTION:
    !
    !- Send data to the coupler. All o2a fields are treated here.
    !- The exchange fields is first copied to the sending buffer.
    !- Then prism_put is called with the buffer.
    !- Two options:
    !    with cpp accu_by_psmile :
    !             accumulation and normalisation is done by psmile.
    !    else:
    !             accumulation and normalisation is done by the model.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Start of couple_put_o2a'
      WRITE(io_stdout,*) ' ***********************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    END IF
#endif

#ifndef __accu_by_psmile
    !
    !-- Accumulate exchange data
    !
    CALL acc_o2a(io_stdout)
    !
    !-- At the end of every coupled time step : normalize
    !
    IF ( MOD(INT(o2a_time), o2a_freq) == 0) THEN
      CALL avg_o2a(io_stdout)
    ENDIF
#else
    !
    !-- Prepare fields for transfer to psmile
    !
    CALL prep_o2a(io_stdout)
#endif
    !
    !-- Export exchange fields; call prism_put every model step
    !
#ifndef __accu_by_psmile
#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' Before exchange ', o2a_time, o2a_freq, MOD(INT(o2a_time),o2a_freq)
      FLUSH(io_stdout)
    END IF
#endif
    IF (MOD(INT(o2a_time),o2a_freq) == 0) THEN
#endif
      DO jfld = 1, nfldo2a

        CALL get_o2a(clstr8o2a(jfld),exfld,jpdim_ocei,jpdim_ocej,io_stdout)

        IF (p_parallel_io) THEN
#ifdef __synout
          WRITE(io_stdout,*)'    '
          WRITE(io_stdout,*)' ==> prism_put for ',clstr8o2a(jfld)
          WRITE(io_stdout,*) &
               ' ==> at run-step ',ldtrun,' (seconds passed=',run_date_secs,')'
          FLUSH(io_stdout)
#endif
          info = PRISM_Ok
          CALL prism_put_proto(portout_id(jfld),run_date_secs,exfld,info)

          IF ( info == PRISM_Sent     .OR. &
               info == PRISM_ToRest   .OR. &
               info == PRISM_SentOut  .OR. &
               info == PRISM_ToRestOut ) num_fld_sent = num_fld_sent + 1

          CALL digest_put_Id &
               (io_stdout,io_stderr,info,clstr8o2a(jfld),run_date_secs,num_fld_sent)
        END IF
      ENDDO

#ifndef __accu_by_psmile

    ENDIF

    !
    !--  Reset exchange fields / accumulation interval
    !

    IF ( MOD(INT(o2a_time),o2a_freq) == 0 .AND. ldtrun*INT(dt) /= ig_ntime) THEN
      CALL ini_o2a(io_stdout)
    ENDIF
#endif
    !
    !-- Reset counter when fields are complete.
    !
    IF ( num_fld_sent == nfldo2a ) THEN

      IF (p_parallel_io) WRITE(io_stdout,*) ' All ',nfldo2a,' fields sent.'
      num_fld_sent = 0

    ENDIF

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' End of couple_put_o2a'
      WRITE(io_stdout,*) ' *********************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    END IF !p_parallel_io
#endif

    ! Update o2a_time
    ! CALL p_bcast(o2a_time,p_io)

  END SUBROUTINE couple_put_o2a

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  couple_end
  !
  ! !INTERFACE:
  !
  SUBROUTINE couple_end
    !
    ! !DESCRIPTION:
    !
    !- Stop communication with oasis.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Start of couple_end '
      WRITE(io_stdout,*) ' *******************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    END IF

    !
    !-- Deallocation of memory
    !
    ierror = PRISM_Ok
    CALL prism_terminate_proto(ierror)
    IF (ierror /= PRISM_Ok) THEN
      WRITE (io_stdout, *) 'An error occured couple_end: Error = ',ierror
    ENDIF

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' End of couple_end'
      WRITE(io_stdout,*) ' *****************'
      WRITE(io_stdout,*) ' '
      FLUSH(io_stdout)
    END IF
#endif

  END SUBROUTINE couple_end

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !ROUTINE:  couple_abort
  !
  ! !INTERFACE:
  !
  SUBROUTINE couple_abort(str1,str2,str3,kerr)

    INTEGER,          INTENT(in) :: kerr      ! std error out unit
    CHARACTER(len=*), INTENT(in) :: str1      ! name of calling model
    CHARACTER(len=*), INTENT(in) :: str2      ! name of calling routine
    CHARACTER(len=*), INTENT(in) :: str3      ! error message

    ! !DESCRIPTION:
    !
    ! - stop the MPI application with messages.
    !
    ! !REVISION HISTORY:
    ! 03.08.07  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------
    ! $Id: mo_couple.f90,v 1.4.2.1.10.1.4.2.2.3.4.1 2006/03/31 08:20:35 m211054 Exp $
    !-----------------------------------------------------------------------
    !
    CALL prism_abort_proto (comp_id, TRIM(str2), TRIM(str3))

  END SUBROUTINE couple_abort

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  grids_writing
  !
  ! !INTERFACE:

  SUBROUTINE grids_writing(kout)
    !
    ! !USES:
    !
    ! This module has to be kept here.
    ! If it is moved up some strange conflicts occur!
    !
    USE mod_prism_grids_writing

    INTEGER, INTENT(in)  :: kout         ! unit number for standard output

    INTEGER              :: gwrite       ! flag to state whether grids writing is
    ! needed or not (1 / 0)
    INTEGER              :: i,j,ii,jj,jb ! looping indicees
    INTEGER, ALLOCATABLE :: mask(:,:)    ! inverse land sea mask
    INTEGER, ALLOCATABLE :: msks(:,:)    ! land sea mask of scalar grid
    REAL(wp), ALLOCATABLE    :: lons(:,:)    ! longitudes of scalars
    REAL(wp), ALLOCATABLE    :: lats(:,:)    ! latitudes of scalars
    REAL(wp), ALLOCATABLE    :: clons(:,:,:) ! corner longitudes of scalars
    REAL(wp), ALLOCATABLE    :: clats(:,:,:) ! corner latitudes of scalars
    REAL(wp), ALLOCATABLE    :: areas(:,:)   ! grid cell area of scalar
    CHARACTER*4          :: grdacr       ! grid acronym (as used in namcouple)
    !
    ! !DESCRIPTION:
    !
    ! - Write grids and masks file for oasis.
    !
    ! !REVISION HISTORY:
    ! July 7, 2003  V. Gayler - created
    !
    ! EOP
    !-----------------------------------------------------------------------
#ifdef __synout
    WRITE(kout,*)' '
    WRITE(kout,*)' Start of grids_writing'
    WRITE(kout,*)' **********************'
    WRITE(kout,*)' '
    CALL flush(kout)
#endif

    ! WRITE(0,*)'grids_writing 3 '

    IF ( p_parallel_io ) THEN

      jj=je_g
      IF ( lbounds_exch_tp ) jj=je_g-2

      ALLOCATE ( mask (ie_g-2,jj),   msks (ie_g-2,jj) )
      ALLOCATE ( lons (ie_g-2,jj),   lats (ie_g-2,jj),   &
                 clons(ie_g-2,jj,4), clats(ie_g-2,jj,4), &
                 areas(ie_g-2,jj) )
    ENDIF

    IF ( p_parallel_io ) THEN

      ! Write grids out only on p_parallel_io
      !
      !-- start writing the arrays
      !
#ifdef __synout
      WRITE(kout,*)'prism_start_grids_writing'
      CALL flush(kout)
#endif

      ! WRITE(0,*)'grids_writing 2 '
      CALL prism_start_grids_writing(gwrite)
      ! WRITE(0,*)'grids_writing 3 '

      IF (gwrite == 1) THEN
        !
        !
        !--  extract scalar/vector grid points
        !
        !     2*ij
        !      :                              s: scalar
        !      :                              u: vector-u
        !      5   s  u  s  u  s  u           v: vector-v
        !      4   v  c  v  c  v  c           c: grid cell corners of scalars
        !      3   s  u  s  u  s  u
        !      2   v  c  v  c  v  c
        !      1   s  u  s  u  s  u           Line 0 and 1 are identical
        !          v  c  v  c  v  c
        !
        !          1  2  3  4  5  6 ... 2*ie
        !
        !
        !
        !--  create lat/lon arrays and
        !    corner arrays for SCRIP interpolation
        !
        ! jb=1
        ! IF ( lbounds_exch_tp ) jb=3
        jb = MERGE(3,1,lbounds_exch_tp)

        !
        !--     scalar
        !

        DO j = jb, je_g
          jj=j+1-jb
          DO i = 2, ie_g-1
            ii=i-1

            lons(ii,jj)    = alon_g(i,j)
            lats(ii,jj)    = alat_g(i,j)

            clons(ii,jj,1) = alonpsi_g(i  ,j-1)     ! upper right corner
            clats(ii,jj,1) = alatpsi_g(i  ,j-1)
            clons(ii,jj,2) = alonpsi_g(i-1,j-1)     ! upper left corner
            clats(ii,jj,2) = alatpsi_g(i-1,j-1)

            clons(ii,jj,3) = alonpsi_g(i-1,j)       ! lower left corner
            clats(ii,jj,3) = alatpsi_g(i-1,j)
            clons(ii,jj,4) = alonpsi_g(i  ,j)       ! lower right corner
            clats(ii,jj,4) = alatpsi_g(i  ,j)

          ENDDO
        ENDDO

        WHERE (clons(:,:,:) < 0.)
          clons(:,:,:) = clons(:,:,:) + 360.
        END WHERE
        WHERE (lons(:,:) < 0.)
          lons(:,:) = lons(:,:) + 360.
        END WHERE
        WHERE (clons(:,:,:) > 360.)
          clons(:,:,:) = clons(:,:,:) - 360.
        END WHERE
        WHERE (lons(:,:) > 360.)
          lons(:,:) = lons(:,:) - 360.
        END WHERE

        !
        !-- create mask arrays, following the OASIS convention
        !
        msks(:,:) = 0
        DO j = jb, je_g
          jj=j+1-jb

          DO i = 2, ie_g-1
            ii=i-1
            mask(ii,jj)=iwetol1_g(i,j)
          ENDDO

        ENDDO

        WHERE (mask(:,:) == 0)
          msks = 1
        END WHERE
        !
        !-- create area arrays, following the OASIS convention
        !
        DO j = jb, je_g
          jj=j+1-jb
          areas(:,jj)=dlxp_g(2:ie_g-1,j)*dlyp_g(2:ie_g-1,j)
        ENDDO
        WHERE (msks(:,:) == 1)
          areas(:,:) = 0.0
        END WHERE

        !
        !-- write grids
        !
        grdacr='oces'
#ifdef __synout
        WRITE(kout,*)'prism_write_grid: ', grdacr
        CALL flush(kout)
#endif
        jj=je_g+1-jb
        CALL prism_write_grid (grdacr, ie_g-2, jj, lons(:,:), lats(:,:))

        !-- write corners
        !
        !   writing corners is optional. If they are missing in the grids file, the
        !   corners are calculated by scrip (with good results).

        grdacr='oces'
#ifdef __synout
        WRITE(kout,*)'prism_write_corner: ', grdacr
        CALL flush(kout)
#endif
        CALL prism_write_corner (grdacr, ie_g-2, jj, 4, clons(:,:,:), clats(:,:,:))
        !
        !-- write masks
        !
        grdacr='oces'
#ifdef __synout
        WRITE(kout,*)'prism_write_mask: ', grdacr
        CALL flush(kout)
#endif
        CALL prism_write_mask (grdacr, ie_g-2, jj, msks(:,:))

        !-- write areas
        !
        grdacr='oces'
#ifdef __synout
        WRITE(kout,*)'prism_write_area: ', grdacr
        CALL flush(kout)
#endif
        CALL prism_write_area (grdacr, ie_g-2, jj, areas(:,:))
        !
        !-- all arrays are written
        !
#ifdef __synout
        WRITE(kout,*)'prism_terminate_grids_writing'
        CALL flush(kout)
#endif
        CALL prism_terminate_grids_writing

#ifdef __synout
      ELSE
        WRITE(kout,*)'  grids files are existing'
        WRITE(kout,*)'  no writing needed'
        CALL flush(kout)
#endif

      ENDIF

#ifdef __synout
      WRITE(kout,*)' '
      WRITE(kout,*)' End of grids_writing'
      WRITE(kout,*)' ********************'
      WRITE(kout,*)' '
      CALL flush(kout)
#endif
    ENDIF ! p_parallel_io

    RETURN

  END SUBROUTINE grids_writing

#ifndef __accu_by_psmile
  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  avg_o2a
  !
  ! !INTERFACE:
  !
  SUBROUTINE avg_o2a(kout)
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(in) :: kout
    !
    ! !LOCAL DECLARATIONS:
    !
    REAL(pwp)                :: fakt  ! normalizing factor
    INTEGER                  :: jx,jy ! loop indices
    INTEGER                  :: jb    ! northern halo
    !
    ! !DESCRIPTION:
    !
    !- average exchange fields at the end of a coupled time step.
    !- overlapped cells are not exchanged.
    !- transform into appropriate units if needed
    !- land points should be set to zero
    !- truncate ice thickness to max allowed value
    !  (might be needed due to truncation erros).
    !!
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    ! 10.03.10  R. Redler, H. Haak - changed from global to local
    !
    ! EOP
    !-----------------------------------------------------------------------
    !
    !     1.   Transform and average the field.
    !          -------------------------------

    ! OASIS3 averaging does not take into account time intervalls
    ! fakt = 1.0_pwp/o2a_time

    fakt = 1.0/float(nacc)

#ifdef __synout
    WRITE(io_stdout,*) ' '
    WRITE(io_stdout,*) ' Averaging of exchange fields.'
    WRITE(io_stdout,*) ' With factor ',fakt
    FLUSH(io_stdout)
#endif

    sstacc (:,:) = sstacc (:,:)*fakt ! *weto(:,:,1)+tmelt
    sitoacc(:,:) = sitoacc(:,:)*fakt ! *weto(:,:,1)
    sicoacc(:,:) = sicoacc(:,:)*fakt ! *weto(:,:,1)
    sntoacc(:,:) = sntoacc(:,:)*fakt ! *weto(:,:,1)
    socuacc(:,:) = socuacc(:,:)*fakt ! *weto(:,:,1)
    socvacc(:,:) = socvacc(:,:)*fakt ! *weto(:,:,1)

#ifdef __cpl_dms
    dmsfacc(:,:) =dmsfacc(:,:)*fakt
    dmsacc (:,:) =dmsacc (:,:)*fakt
#endif
#ifdef __cpl_co2
    co2tracc(:,:) = co2tracc(:,:)*fakt
    co2suacc(:,:) = co2suacc(:,:)*fakt
#endif


    !     2.    No computation of effectice ice thickness.
    !           Exchange sea ice in ice-meter
    !                    snow in fresh water equivalent.
    !           ----------------------------------------
    !
    !     box average ice thickness (sitoacc) divided by compactness (sicoacc)
    !     gives flow thickness which is required by ECHAM
    !     snow is converted to luiquid water equivalent
    !
    ! ... moved to acc_o2a ...
    !

!!$    DO jy=2,je-1
!!$      DO jx=2,ie-1
!!$
!!$        IF (sicoacc(jx,jy) > 0.0 ) THEN
!!$          sitoacc(jx,jy) = MIN(1000.0,(sitoacc(jx,jy)/sicoacc(jx,jy)))
!!$          sntoacc(jx,jy) = MIN(1000.0,(sntoacc(jx,jy)/sicoacc(jx,jy)))
!!$          sntoacc(jx,jy) = gl_sntoacc(jx,jy)*rhoref_snow/rhoref_water
!!$        ELSE
!!$          sitoacc(jx,jy) = 0.0
!!$          sntoacc(jx,jy) = 0.0
!!$          sicoacc(jx,jy) = 0.0
!!$        ENDIF
!!$
!!$      ENDDO
!!$    ENDDO

    jb = MERGE(3,1,lbounds_exch_tp)

    ! Fields need to be gathered on a temporary array
    ! because the gather routine requires arrays that
    ! include the halo while OASIS3 requires fields
    ! without halo.

#ifdef __synout
    IF ( p_pe == p_io ) THEN
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathering of exchange fields!', ie_g-1, jb, je_g
      FLUSH(io_stdout)
    ENDIF
#endif

    CALL gather(sstacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_sstacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Transferred exchange field 1.', MINVAL(gl_sstacc), MAXVAL(gl_sstacc)
      FLUSH(io_stdout)
#endif
    ENDIF

    CALL gather(sitoacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_sitoacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 2.', MINVAL(gl_sitoacc), MAXVAL(gl_sitoacc)
      FLUSH(io_stdout)
#endif
    ENDIF


    CALL gather(sicoacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_sicoacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 3.', MINVAL(gl_sicoacc), MAXVAL(gl_sicoacc)
      FLUSH(io_stdout)
#endif
    ENDIF


    CALL gather(sntoacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_sntoacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 4.', MINVAL(gl_sntoacc), MAXVAL(gl_sntoacc)
      FLUSH(io_stdout)
#endif
    ENDIF

    CALL gather(socuacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_socuacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 5.', MINVAL(gl_socuacc), MAXVAL(gl_socuacc)
      FLUSH(io_stdout)
#endif
    ENDIF

    CALL gather(socvacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_socvacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 6.', MINVAL(gl_socvacc), MAXVAL(gl_socvacc)
      FLUSH(io_stdout)
#endif
    ENDIF

    IF ( p_pe == p_io ) &
         CALL rotate_2_north_east(gl_socuacc,gl_socvacc,jpdim_ocei,jpdim_ocej)

#ifdef __cpl_co2
    CALL gather(co2tracc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_co2tracc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 7a.', MINVAL(gl_co2tracc), MAXVAL(gl_co2tracc)
      CALL flush(io_stdout)
#endif
    ENDIF

    CALL gather(co2suacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_co2suacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 8a.', MINVAL(gl_co2suacc), MAXVAL(gl_co2suacc)
      CALL flush(io_stdout)
#endif
    ENDIF
#endif

#ifdef __cpl_dms
    CALL gather(dmsfacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_dmsfacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 7b.', MINVAL(gl_dmsfacc), MAXVAL(gl_dmsfacc)
      CALL flush(io_stdout)
#endif
    ENDIF

    CALL gather(dmsacc, pfield_g, p_io)
    IF ( p_pe == p_io ) THEN
      DO jy=jb,je_g
        DO jx=2,ie_g-1
          gl_dmsacc(jx-1,jy-jb+1) = pfield_g(jx,jy)
        ENDDO
      ENDDO
#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Gathered exchange field 8b.', MINVAL(gl_dmsacc), MAXVAL(gl_dmsacc)
      CALL flush(io_stdout)
#endif
    ENDIF
#endif

  END SUBROUTINE avg_o2a

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  acc_o2a
  !
  ! !INTERFACE:
  !
  SUBROUTINE acc_o2a(kout)
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(IN) :: kout
    !
    ! !LOCAL DECLARATIONS:
    !
    INTEGER             :: jx,jy ! loop indices
    REAL(wp)                :: fak
    !
    ! !DESCRIPTION:
    !
    !- accumulate ocean sst, ice thickness and snow depth .
    !                   tho, sictho, sicomo, sicsno
    !
    ! !REVISION HISTORY:
    ! 03.05.22  s. legutke - created
    ! 10.03.10  R. Redler, H. Haak - changed from global to local
    !
    ! EOP
    !----------------------------------------------------------------------
    !
    !
    !     1.    No computation of effectice ice thickness.
    !           Exchange sea ice in ice-meter
    !                    snow in fresh water equivalent.
    !           ----------------------------------------
    !
    !     box average ice thickness (sitoacc) divided by compactness (sicoacc)
    !     gives flow thickness which is required by ECHAM
    !     snow is converted to luiquid water equivalent
    !
    ! halos are ignored here because they cut off later anyway
    !
    DO jy=2,je-1
      DO jx=2,ie-1
        IF (sicomo(jx,jy) > 0.0 ) THEN
          sicotmp(jx,jy) = sicomo(jx,jy)
          sitotmp(jx,jy) = MIN(1000.0,(sictho(jx,jy)/sicomo(jx,jy)))
          sntotmp(jx,jy) = MIN(1000.0,(sicsno(jx,jy)/sicomo(jx,jy)))*rhoref_snow/rhoref_water
        ELSE
          sicotmp(jx,jy) = 0.0
          sitotmp(jx,jy) = 0.0
          sntotmp(jx,jy) = 0.0
        ENDIF
      ENDDO
    ENDDO

    !
    ! 2.) accumulate local fields
    !     OASIS3 averaging does not take into account time intervalls
    !
    sstacc (:,:) = sstacc (:,:) + tho    (:,:,1)*weto (:,:,1)+tmelt  !*dt
    sitoacc(:,:) = sitoacc(:,:) + sitotmp(:,:  )*weto (:,:,1)        !*dt
    sicoacc(:,:) = sicoacc(:,:) + sicotmp(:,:  )*weto (:,:,1)        !*dt
    sntoacc(:,:) = sntoacc(:,:) + sntotmp(:,:  )*weto (:,:,1)        !*dt

#ifdef __cpl_dms
    dmsfacc(:,:) = dmsfacc(:,:) + dmsflux(:,:  )*weto (:,:,1)
    dmsacc (:,:) = dmsacc (:,:) + dms    (:,:  )*weto (:,:,1)
#endif
#ifdef __cpl_co2
    co2tracc(:,:) = co2tracc(:,:) + co2trans(:,:  )*weto (:,:,1)
    co2suacc(:,:) = co2suacc(:,:) + suppco2 (:,:  )*weto (:,:,1)
#endif

    DO jy = 2,je-1
      DO jx = 2,ie-1

        fak = MERGE(-1., 1., (lbounds_exch_tp .AND. have_g_js &
             .AND. jy == 3 .AND. jx+p_ioff >= ie_g/2))

        socuacc (jx,jy) = socuacc(jx,jy)  +  weto (jx,  jy,1) *                    &
                     ((1.-sicomo (jx,jy)) * (uko  (jx+1,jy,1)+uko  (jx,jy,1))      &
                        + sicomo (jx,jy)  * (sicuo(jx+1,jy  )+sicuo(jx,jy  )))*0.5     !*dt

        socvacc (jx,jy) = socvacc(jx,jy)  +  weto (jx,jy,1) *                      &
                     ((1.-sicomo (jx,jy)) * (vke  (jx,jy,1)+fak*vke  (jx,jy-1,1))  &
                        + sicomo (jx,jy)  * (sicve(jx,jy  )+fak*sicve(jx,jy-1  )))*0.5 !*dt
      ENDDO
    ENDDO

    o2a_time = o2a_time + dt
    nacc     = nacc + 1

#ifdef __synout
    WRITE(kout,*) ' accumulating since ', o2a_time
    WRITE(kout,*) ' sstacc:  ', MAXVAL(sstacc),  MINVAL(sstacc)
    WRITE(kout,*) ' sitoacc: ', MAXVAL(sitoacc), MINVAL(sitoacc)
    WRITE(kout,*) ' sicoacc: ', MAXVAL(sicoacc), MINVAL(sicoacc)
    WRITE(kout,*) ' sntoacc: ', MAXVAL(sntoacc), MINVAL(sntoacc)
#ifdef __cpl_dms
    WRITE(kout,*) ' dmsfacc: ' , MAXVAL(dmsfacc),  MINVAL(dmsfacc)
    WRITE(kout,*) ' dmsacc : ' , MAXVAL(dmsacc),   MINVAL(dmsacc)
#endif
#ifdef __cpl_co2
    WRITE(kout,*) ' co2tracc: ' , MAXVAL(co2tracc),  MINVAL(co2tracc)
    WRITE(kout,*) ' co2suacc: ' , MAXVAL(co2suacc),  MINVAL(co2suacc)
#endif
    FLUSH(kout)
#endif

    ! brodcast ocean time to all proccessors

    ! CALL p_bcast(o2a_time,p_io)

  END SUBROUTINE acc_o2a
#endif

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  ini_o2a
  !
  ! !INTERFACE:
  !
  SUBROUTINE ini_o2a(kout)
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(in) :: kout
    !
    ! !LOCAL DECLARATIONS:
    !
    !
    ! !DESCRIPTION:
    !
    !- initialise arrays for  accumulating ocean SST (and ice/snow
    !  thickness) at the beginning of each coupled time step.
    !!
    ! !REVISION HISTORY:
    ! 03.05.22  s. legutke - created
    !
    ! EOP
    !----------------------------------------------------------------------

    IF (p_parallel_io) THEN
      gl_sstacc (:,:) = 0.0
      gl_sitoacc(:,:) = 0.0
      gl_sicoacc(:,:) = 0.0
      gl_sntoacc(:,:) = 0.0
      gl_socuacc(:,:) = 0.0
      gl_socvacc(:,:) = 0.0
    ENDIF

#ifndef __accu_by_psmile
    sstacc (:,:) = 0.0
    sitoacc(:,:) = 0.0
    sicoacc(:,:) = 0.0
    sntoacc(:,:) = 0.0
    socuacc(:,:) = 0.0
    socvacc(:,:) = 0.0
#ifdef __cpl_dms
    dmsfacc(:,:) = 0.0
    dmsacc (:,:) = 0.0
#endif
#ifdef __cpl_co2
    co2tracc(:,:) = 0.0
    co2suacc(:,:) = 0.0
#endif

    o2a_time = 0.0_pwp
    nacc     = 0
#endif

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(kout,*) ' accumulation time and exchange fields initialized.'
      FLUSH(kout)
    END IF
#endif

  END SUBROUTINE ini_o2a

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  ini_wrte
  !
  ! !INTERFACE:
  !
  SUBROUTINE ini_wrte(kunit,kout)
    !
    INCLUDE 'netcdf.inc'
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(in)    :: kout
    INTEGER, INTENT(in)    :: kunit
    !
    ! !LOCAL DECLARATIONS:
    !
    INTEGER      :: jx,jy,jym1,jj,jb                    ! loop indices
    REAL(pwp)    :: exfld(jpdim_ocei,jpdim_ocej)
    REAL(wp)         :: fak
    !
    ! for NetCDF restart files
    !
    INTEGER      ::  stat, ncid, varid, idx, idy
    !
    ! !DESCRIPTION:
    ! - if no restart files for the interface oasis are
    ! available from a preceeding run, i.e. the run is the first of the
    ! experiment (nmseq>1), at least one model needs data to start
    ! the run and provide the data for the other models.
    ! here the ocean surface conditions file  is created using the actual
    ! sst and ice conditions of the ocean restart file or of the initial
    ! ocean conditions.
    !- change unit
    !
    ! !REVISION HISTORY:
    ! 03.05.22  s. legutke - created
    !
    ! EOP
    !----------------------------------------------------------------------

    jb = MERGE(3,1,lbounds_exch_tp)

    !
    !--   celsius --> kelvin
    !     compute effective ice thickness and initialize
    !

    ! gather fields for accumulation

#ifdef __synout
    WRITE(kout,*)'Gathering global fields'
    CALL flush(kout)
#endif
    CALL gather(tho(:,:,1),gl_sst,p_io)
    CALL gather(sictho,    gl_sictho,p_io)
    CALL gather(sicomo,    gl_sicomo,p_io)
    CALL gather(sicsno,    gl_sicsno,p_io)
    CALL gather(uko(:,:,1),gl_socu,p_io)
    CALL gather(vke(:,:,1),gl_socv,p_io)
    CALL gather(sicuo(:,:),gl_sicu,p_io)
    CALL gather(sicve(:,:),gl_sicv,p_io)

#ifdef __cpl_dms
    CALL gather(dmsflux,  gl_dmsflux, p_io)
    CALL gather(dms,      gl_dms,     p_io)
#endif
#ifdef __cpl_co2
    co2trans=0.0
    suppco2=280.0
    CALL gather(co2trans, gl_co2trans, p_io)
    CALL gather(suppco2,  gl_suppco2, p_io)
#endif


    IF (p_parallel_io) THEN

#ifdef __synout
      WRITE(kout,*)'processing global fields'
      CALL flush(kout)
#endif
      DO jy=1,jpdim_ocej
        jj=jy+jb-1
        DO jx=1,jpdim_ocei

          IF (gl_sicomo(jx+1,jj) > 0.0) THEN
            gl_sicoacc(jx,jy)=gl_sicomo(jx+1,jj)
            gl_sntoacc(jx,jy)=MIN(1000.0,gl_sicsno(jx+1,jj)/gl_sicomo(jx+1,jj))
            gl_sitoacc(jx,jy)=MIN(1000.0,gl_sictho(jx+1,jj)/gl_sicomo(jx+1,jj))
          ELSE
            gl_sitoacc(jx,jy) = 0.0
            gl_sicoacc(jx,jy) = 0.0
            gl_sntoacc(jx,jy) = 0.0
          ENDIF

          gl_sstacc (jx,jy) = gl_sst(jx+1,jj)*wetol1_g(jx+1,jj)+tmelt
          gl_sicoacc(jx,jy) = gl_sicoacc(jx,jy)*wetol1_g(jx+1,jj)
          gl_sntoacc(jx,jy) = gl_sntoacc(jx,jy)*wetol1_g(jx+1,jj)*rhoref_snow/rhoref_water
          gl_sitoacc(jx,jy) = gl_sitoacc(jx,jy)*wetol1_g(jx+1,jj)

          gl_socuacc(jx,jy) = wetol1_g(jx+1,jj)*                                 &
                        ((1.-gl_sicomo(jx+1,jj))*(gl_socu(jx+1,jj)+gl_socu(jx,jj))        &
                           + gl_sicomo(jx+1,jj) *(gl_sicu(jx+1,jj)+gl_sicu(jx,jj)))*0.5
          jym1=jj-1
          IF ( jym1 == 0) jym1 = 1
          fak=1.

          IF ( lbounds_exch_tp .AND. jym1 == 2 .AND. jx >= ie_g/2 ) fak=-1.

          gl_socvacc (jx,jy) =  wetol1_g(jx+1,jj)*                                &
                          ((1.-gl_sicomo(jx+1,jj))*(gl_socv(jx+1,jj)+fak*gl_socv(jx+1,jym1)) &
                             + gl_sicomo(jx+1,jj) *(gl_sicv(jx+1,jj)+fak*gl_sicv(jx+1,jym1)))*0.5

#ifdef __cpl_co2
          gl_co2tracc(jx,jy)=gl_co2trans(jx+1,jj)
          gl_co2suacc(jx,jy)=gl_suppco2(jx+1,jj)
#endif

#ifdef __cpl_dms
          gl_dmsfacc(jx,jy)=gl_dmsflux(jx+1,jj)
          gl_dmsacc(jx,jy)=gl_dms(jx+1,jj)
#endif




        ENDDO
      ENDDO

#ifdef __synout
      WRITE(kout,*) 'Rotate global fields'
      CALL flush(kout)
#endif
      CALL rotate_2_north_east(gl_socuacc,gl_socvacc,jpdim_ocei,jpdim_ocej)


    ENDIF

#ifndef OASIS_RESTART_FORMAT_RAW
    !
    !*    1.0     store restarts in netcdf.
    !             -------------------------------
    DO jfld = 1, nfldo2a

      CALL get_o2a(clstr8o2a(jfld),exfld,jpdim_ocei,jpdim_ocej,kout)

      IF (p_parallel_io) THEN

        !
        !*    1.1     open oasis netcdf restart file.
        !             -------------------------------
#ifdef __synout
        WRITE(kout,*)'opening file ',cg_clim_rstfile(jfld),' unit=',kunit
        CALL flush(kout)
#endif
        !rar      stat = nf_create(TRIM(cg_clim_rstfile(jfld))//'.nc',NF_CLOBBER,ncid)
        stat = nf_create(TRIM(cg_clim_rstfile(jfld)),NF_CLOBBER,ncid)
        CALL hdlerr(stat,'creating file')
        !
        !*    1.2     define netcdf header.
        !             -------------------------------
        stat = nf_def_dim(ncid, 'x', jpdim_ocei, idx)
        CALL hdlerr(stat,'define x dimension')

        stat = nf_def_dim(ncid, 'y', jpdim_ocej, idy)
        CALL hdlerr(stat,'define y dimension')

        stat=nf_def_var(ncid, TRIM(clstr8o2a(jfld)), NF_DOUBLE, 2, (/idx,idy/), varid)
        CALL hdlerr(stat,'define variable '//clstr8o2a(jfld))

        stat = nf_enddef(ncid)
        CALL hdlerr(stat,'error with enddef')
        !
        !*    1.3     write data array.
        !             -------------------------------
#ifdef __synout
        WRITE(kout,*)'ini_wrte : ',clstr8o2a(jfld),             &
             (exfld(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
        CALL flush(kout)
#endif
        stat=nf_put_var_double(ncid,varid,exfld)
        CALL hdlerr(stat,'put variable '//clstr8o2a(jfld))
        !
        !*    1.4     close file.
        !             -------------------------------
        stat = nf_close(ncid)
        CALL hdlerr(stat,'error closing file')

      END IF

    END DO
#ifdef __synout
    WRITE(kout,*)'All files written!'
    CALL flush(kout)
#endif
#else
    !
    !*    1.0     store restarts in raw binary.
    !            -----------------------------------------
    DO jfld = 1,nfldo2a
#ifdef __synout
      WRITE(kout,*)'opening file ',cg_clim_rstfile(jfld),' unit=',kunit
      CALL flush(kout)
#endif
      OPEN (kunit,status='unknown',file=cg_clim_rstfile(jfld), &
           form='unformatted',position='APPEND')
      CALL get_o2a(clstr8o2a(jfld),exfld,jpdim_ocei,jpdim_ocej,kout)
      IF (p_parallel_io) THEN
        WRITE(kunit) clstr8o2a(jfld)
        WRITE(kunit) exfld
#ifdef __synout
        WRITE(kout,*)'ini_wrte : ',clstr8o2a(jfld),             &
             (exfld(jpdim_ocei/2+1,jy),jy=3,jpdim_ocej,15)
        CALL flush(kout)
#endif
      END IF
      CLOSE (kunit)
    ENDDO

    IF (p_parallel_io) THEN
      FLUSH(kunit)
      CLOSE(kunit)
    ENDIF ! p_parallel_io
#endif

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(kout,*)
      WRITE(kout,*)' end of ini_wrte'
      WRITE(kout,*)' ***************'
      FLUSH(kout)
    ENDIF
#endif

  END SUBROUTINE ini_wrte

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  get_o2a
  !
  ! !INTERFACE:
  !
  SUBROUTINE get_o2a(clfield,pfield,jpdimi,jpdimj,kout)
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(in) :: jpdimi,jpdimj         ! field dimensions
    INTEGER, INTENT(in) :: kout                  ! unit for std output

    CHARACTER(len=8), INTENT(in) :: clfield      ! fields symbolic name
    !
    ! !RETURN VALUE:
    !
    REAL(wp), INTENT(out)   :: pfield(jpdimi,jpdimj) ! target field

    !
    ! !LOCAL DECLARATIONS:
    !
    INTEGER :: index                 ! field search index
    INTEGER :: ierr                  ! field search index
    INTEGER :: jx,jy                 ! loop indices
    !
    ! !DESCRIPTION:
    ! - moves the model field to be exchange to the array sent to OASIS.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !----------------------------------------------------------------------
    !
    ierr = 1

    index = 1

    IF ( clfield == clstr8o2a(index) ) THEN
      IF (p_parallel_io) THEN
        DO jy = 1, jpdimj
          DO jx = 1, jpdimi
            pfield(jx,jy) = gl_sstacc(jx,jy)
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      END IF
      ierr = 0
      RETURN

    ENDIF

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN
      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield(jx,jy) = gl_sitoacc(jx,jy)
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN
      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield(jx,jy) = gl_sicoacc(jx,jy)
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN
      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield(jx,jy) = gl_sntoacc(jx,jy)
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN
      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield(jx,jy) = gl_socuacc(jx,jy)
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN
      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield(jx,jy) = gl_socvacc(jx,jy)
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF

#ifdef __cpl_dms

    ! ATTN: dmsflux is only needed for diagnostics;
    !       for performace reasons it should be removed
    !       from the coupling interface

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN

      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield (jx,jy) = gl_dmsfacc(jx,jy)*32./dt ! [kmol/[m**2 timestep]] -->  [kg(S)/[m**2 s]]
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN

      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield (jx,jy) = gl_dmsacc(jx,jy)*1.e9  ![kmol/m**3 --> nanomol/l]
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF

#endif /*__cpl_dms*/

#ifdef __cpl_co2

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN

      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield (jx,jy) = gl_co2tracc(jx,jy)
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF

    index = index + 1

    IF ( clfield == clstr8o2a(index) ) THEN

      IF (p_parallel_io) THEN
        DO jy=1,jpdimj
          DO jx=1,jpdimi
            pfield (jx,jy) = gl_co2suacc(jx,jy) ![ppm CO2]
          ENDDO
        ENDDO
#ifdef __synout
        WRITE(kout,*)' get_o2a: ',clstr8o2a(index),MINVAL(pfield),MAXVAL(pfield)
#endif /*__synout*/
      ENDIF
      ierr = 0
      RETURN

    ENDIF
#endif /*__cpl_co2*/

#ifdef __synout
    FLUSH(kout)
#endif
    IF ( ierr > 0 ) &
         CALL couple_abort (modnam,'get_o2a','pb Invalid locator string: '//clfield,io_stderr)

  END SUBROUTINE get_o2a


  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE: put_a2o
  !
  ! !INTERFACE:
  !
  SUBROUTINE put_a2o(clfield,pfield,jpdimi,jpdimj,kout)
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(IN) :: jpdimi,jpdimj         ! field dimensions
    INTEGER, INTENT(IN) :: kout                  ! unit for std output
    INTEGER             :: jj, jb, il, ir,i        ! loop indices

    CHARACTER(len=8), INTENT(IN) :: clfield      ! fields symbolic name

    REAL(wp), INTENT(IN) :: pfield(jpdimi,jpdimj) ! source field
    !
    ! !LOCAL DECLARATIONS:
    !
    REAL(wp) :: pq,fak
    INTEGER :: jx,jy,j                 ! loop indeces
    !
    ! !DESCRIPTION:
    ! - only the 'inner' oceanic grid region excluding the cyclic boundary
    ! columns are exchanged by OASIS.
    ! put_a2o distributes them to the respective array defined
    ! for the complete ocean grid.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !----------------------------------------------------------------------

    !
    !JJ 01.10.99
    !JJ ATTN: FOR READING AND WIND STRESS ROTATION TEMP 2D FIELDS
    !  ARE USED
    !::      V1E===wind stress y water on u-point
    !::      Z1E===wind stress y ice on u-point
    !::      U1E===wind stress x water on v-point
    !::      B1E===wind stress x ice on v-point
    !*       2.    Forcing with fluxes only.
    !              -------------------------
    !
    !JJ 27.9.99 INTRODUCE WIND STRESS COMPONENTS ON
    !    VELOCITY POINTS, USE SURFACE FIELDS TX,TY AS TEMPORARY
    !    VARABLES BEFORE ROTATION INTO NEW GRID ORIENTATION

    jb = MERGE(3,1,lbounds_exch_tp .and. have_g_js)

    !*       2.1   Zonal wind stress over water.
    !              ----------------------------

    IF ( clfield == 'TXWOCEAS' ) THEN

      !*       2.1a  Zonal wind stress over water

      IF (p_parallel_io) THEN

        ! copy into pfield with halos
        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aofltxwo,'p','in put_a2o 1') ! x-stress on u-point

      IF (p_parallel_io) THEN

        ! copy into pfield with halos
        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(u1e,'p','in put_a2o 1') ! x-stress on v-point


#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.2   Meridional wind stress water.
    !              ----------------------------
    IF ( clfield == 'TYWOCEAS' ) THEN

      !*       2.2a  Meridional wind stress water on u-point

      IF (p_parallel_io) THEN

        ! copy into pfield with halos
        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(v1e,'p-','in put_a2o 3') ! y-stress on u-point

      IF (p_parallel_io) THEN

        ! copy into pfield with halos
        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aofltywe,'p-','in put_a2o 3') ! y-stress on v-point


#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.3   Zonal wind stress ice.
    !              ---------------------

    IF ( clfield == 'TXIOCEAS' ) THEN

      !*       2.3a  Zonal wind stress ice on u-point

      IF (p_parallel_io) THEN

        ! copy into pfield with halos
        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aofltxio,'p','in put_a2o 5') ! x-stress on u-point

      IF (p_parallel_io) THEN


         ! copy into pfield with halos
        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(b1e,'p','in put_a2o 5') ! x-stress on v-point

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.4   Meridional wind stress ice.
    !              --------------------------

    IF ( clfield == 'TYIOCEAS' ) THEN

      !*       2.4a  Meridional wind stress ice on u-point

      IF (p_parallel_io) THEN

        ! copy into pfield with halos
        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(z1e,'p-','in put_a2o 7') !y-stress on u-point

      IF (p_parallel_io) THEN


         ! copy into pfield with halos
        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aofltyie,'p-','in put_a2o 7') !y-stress on u-point


#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.5   Freshwater flux over ice.
    !              -------------------------

    IF ( clfield == 'FRIOCEAN' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb, je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aoflfrio,'p','in put_a2o 9') ! Freshwater flux over ice

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.6   Freshwater flux over water.
    !              --------------------------

    IF ( clfield == 'FRWOCEAN' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb, je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aoflfrwo,'p','in put_a2o 10') ! Freshwater flux over water

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.7   Residual heat flux over ice.
    !              ----------------------------

    IF ( clfield == 'RHIOCEAN' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb, je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aoflrhio,'p','in put_a2o 11') ! Residual heat flux over ice

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.8   Conductive heat flux through ice.
    !              ---------------------------------

    IF ( clfield == 'CHIOCEAN' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aoflchio,'p','in put_a2o 12') ! Conductive heat flux through ice

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.9   Net heat flux over water.
    !              -------------------------

    IF ( clfield == 'NHWOCEAN' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb,je_g
          jj = jy+1-jb
          DO  jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aoflnhwo,'p','in put_a2o 13') ! Net heat flux over water


#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.10   Solar heat flux.
    !               ----------------

    IF ( clfield == 'SHWOCEAN' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aoflshwo,'p','in put_a2o 14') ! Solar heat flux

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.11   Wind Stress Velocity.
    !               --------------------

    IF ( clfield == 'WSVOCEAN' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(aoflwsvo,'p','in put_a2o 15') ! Wind Stress Velocity

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

#ifdef __cpl_dust

    !*       2.12   Dust dep. flux
    !               --------------------

    IF ( clfield == 'DEPOCEAN' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(dustdep,'p','in put_a2o 16') ! Dust dep. flux

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF
#endif

#ifdef __cpl_co2

    !*       2.12   CO2 concentration
    !               --------------------

    IF ( clfield == 'CO2CONOC' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(co2conc,'p','in put_a2o 17') ! CO2 concentration

      atm(:,:,iatmco2) = co2conc(:,:)*(28.970/44.011)*1.e6

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF

    !*       2.12   CO2 flux
    !               ----------

    IF ( clfield == 'CO2FLXOC' ) THEN

      IF (p_parallel_io) THEN

        DO jy = jb,je_g
          jj = jy+1-jb
          DO jx = 2,ie_g-1
            pfield_g(jx,jy) = pfield(jx-1,jj)*wetol1_g(jx,jy)
          ENDDO
        ENDDO

      ENDIF

      CALL scatter_bounds(co2flux_cpl,'p','in put_a2o 18') ! CO2 flux

#ifdef __synout
      IF (p_parallel_io) &
           WRITE(kout,*) ' put_a2o : ',clfield, (pfield((ie_g-2)/2,jy),jy=2,je_g,14)
#endif
      RETURN

    ENDIF
#endif

    CALL couple_abort (modnam,'put_a2o','pb Invalid locator string: '//clfield,io_stderr)

6001 FORMAT(1x,4x,6(f6.3,3x))
6002 FORMAT(1x,6(f6.3,3x))

  END SUBROUTINE put_a2o


  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE: post_a2o
  !
  ! !INTERFACE:
  !
  SUBROUTINE post_a2o(k_fld_recvd, kout, kdtrun)
    !
    ! !USES:
    !
    USE mo_fluxes1, ONLY : wrte_flux_extra

    IMPLICIT NONE
    !
    ! !INPUT/OUTPUT VALUE:
    !
    INTEGER, INTENT(inout) :: k_fld_recvd    ! counter of received fields
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(in)    :: kout
    INTEGER, INTENT(in)   ::  kdtrun         ! model time step
    !
    ! !LOCAL DECLARATIONS:
    !
    INTEGER               ::  ncorrect
    INTEGER               ::  jx,jy          ! loop indices
    !
    ! !DESCRIPTION:
    ! -
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !----------------------------------------------------------------------

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(kout,*)
      WRITE(kout,*)' Start of  post_a2o'
      WRITE(kout,*)' ******************'
      WRITE(kout,*)
      FLUSH(kout)
    ENDIF
#endif
    !
    !-- Rotate received wind stress vector (JJ 01.10.99)
    !

    !average from p-point to u-point
    DO jx = 1,ie-1
      aofltxwo(jx,:)=0.5*(aofltxwo(jx,:)+aofltxwo(jx+1,:))*amsuo(jx,:,1)
      aofltxio(jx,:)=0.5*(aofltxio(jx,:)+aofltxio(jx+1,:))*amsuo(jx,:,1)
      v1e(jx,:)=0.5*(v1e(jx,:)+v1e(jx+1,:))*amsuo(jx,:,1)
      z1e(jx,:)=0.5*(z1e(jx,:)+z1e(jx+1,:))*amsuo(jx,:,1)
    ENDDO

    CALL bounds_exch(1,'u',aofltxwo)
    CALL bounds_exch(1,'u',v1e)
    CALL bounds_exch(1,'u',aofltxio)
    CALL bounds_exch(1,'u',z1e)

    CALL rotate2_u(aofltxwo,v1e,ie,je)
    CALL rotate2_u(aofltxio,z1e,ie,je)

    !average from p-point to v-point
    DO jy = 1,je-1
      aofltywe(:,jy)=0.5*(aofltywe(:,jy)+aofltywe(:,jy+1))*amsue(:,jy,1)
      aofltyie(:,jy)=0.5*(aofltyie(:,jy)+aofltyie(:,jy+1))*amsue(:,jy,1)
      u1e(:,jy)=0.5*(u1e(:,jy)+u1e(:,jy+1))*amsue(:,jy,1)
      b1e(:,jy)=0.5*(b1e(:,jy)+b1e(:,jy+1))*amsue(:,jy,1)
    ENDDO

    CALL bounds_exch(1,'vf',aofltywe)
    CALL bounds_exch(1,'vf',aofltyie)
    CALL bounds_exch(1,'vf',u1e)
    CALL bounds_exch(1,'vf',b1e)

    CALL rotate2_v(u1e,aofltywe,ie,je)
    CALL rotate2_v(b1e,aofltyie,ie,je)

    IF(icycli >= 1) THEN
      CALL bounds_exch(1,'u',aofltxwo)
      CALL bounds_exch(1,'u',aofltxio)
      CALL bounds_exch(1,'vf',aofltywe)
      CALL bounds_exch(1,'vf',aofltyie)
    ENDIF

#ifdef FLUXCORRECT
    WRITE(kout,*) 'Adding flux correction! '
    fhflmax=-999.
    fhflmin=999.
#endif /*FLUXCORRECT*/

    DO jy=1,je
      DO jx=1,ie
#ifdef FLUXCORRECT

#ifdef TEMPCORRECT
        aoflnhwo(jx,jy)=aoflnhwo(jx,jy)+fluko_hfl(jx,jy)
        aoflchio(jx,jy)=aoflchio(jx,jy)+fluko_hfl(jx,jy)
        IF (fhflmax < fluko_hfl(jx,jy)) fhflmax=fluko_hfl(jx,jy)
        IF (fhflmin > fluko_hfl(jx,jy)) fhflmin=fluko_hfl(jx,jy)
#endif /*TEMPCORRECT*/

#ifdef SALTCORRECT
        aoflfrwo(jx,jy)=aoflfrwo(jx,jy)+fluko_frw(jx,jy)
#endif /*SALTCORRECT*/

#endif /*FLUXCORRECT*/

      ENDDO
    ENDDO

#ifdef FLUXCORRECT
    WRITE(kout,*) 'max hfl fluxco ',fhflmax
    WRITE(kout,*) 'min hfl fluxco ',fhflmin
#endif /*FLUXCORRECT*/

    !-- Modify residual/conductive heat flux if negative

    ncorrect = 0
    DO jy = 1,je
      DO jx = 1,ie
        IF (aoflrhio(jx,jy) < 0.0 .AND. lweto(jx,jy,1)) THEN
          IF ( ncorrect <= 10 ) THEN
            IF ( ncorrect == 0 ) WRITE(kout,*)          &
                 &                   'Model time step of run : ',kdtrun
            WRITE(kout,*) 'Modify residual heat flux at :'&
                 &                   ,jx,jy,aoflrhio(jx,jy),aoflchio(jx,jy)
            ncorrect = ncorrect + 1
          ENDIF
          aoflchio(jx,jy) = aoflchio(jx,jy) + aoflrhio(jx,jy)
          aoflrhio(jx,jy) = 0.0
        ENDIF
      ENDDO
    ENDDO

    !
    !-- Move wind stress velocity array used in OCTHER
    !
    DO jy=1,je
      DO jx=1,ie
        fu10(jx,jy)=aoflwsvo(jx,jy)
      ENDDO
    ENDDO

    k_fld_recvd = 0
    !
#ifdef CMIP_READ_FLUX
    ! OVERREAD FLUXES BY STORED FIELDS
    CALL READ_FLUX_EXTRA
#endif
    IF ( IOASISFLUX.EQ.99 ) THEN
      CALL WRTE_FLUX_EXTRA
    ENDIF

#ifdef __synout
    WRITE(kout,*)
    WRITE(kout,*)' End of  post_a2o'
    WRITE(kout,*)' ****************'
    FLUSH(kout)
#endif
  END SUBROUTINE post_a2o


  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE: prep_o2a
  !
  ! !INTERFACE:
  !
  SUBROUTINE prep_o2a(kout)
    !
    ! !RETURN VALUE:
    !
    INTEGER, INTENT(in)    :: kout
    !
    ! !LOCAL DECLARATIONS:
    !
    INTEGER                :: jx,jy      ! loop indices
    INTEGER                :: jym1,jj,jb ! pointer to selected j-rows
    REAL(wp)                   :: fak
    !
    ! !DESCRIPTION:
    ! - prepare fields before sending to psmile.
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !----------------------------------------------------------------------
    !
    !-- Masking and unit transformation.
    !

    ! Gather arrays for accumulation
    CALL gather(tho(:,:,1),gl_sst,p_io)
    CALL gather(sictho,    gl_sictho,p_io)
    CALL gather(sicomo,    gl_sicomo,p_io)
    CALL gather(sicsno,    gl_sicsno,p_io)
    CALL gather(uko(:,:,1),gl_socu,p_io)
    CALL gather(vke(:,:,1),gl_socv,p_io)
    CALL gather(sicuo(:,:),gl_sicu,p_io)
    CALL gather(sicve(:,:),gl_sicv,p_io)

    ! jb=0
    ! IF ( lbounds_exch_tp ) jb=2

    jb = MERGE(2,0,lbounds_exch_tp)

    IF (p_parallel_io) THEN

      DO jy=1,jpdim_ocej
        jj=jy+jb
        DO jx=1,jpdim_ocei

          gl_sitoacc(jx,jy) = gl_sictho(jx+1,jj)*wetol1_g(jx+1,jj)
          gl_sicoacc(jx,jy) = gl_sicomo(jx+1,jj)*wetol1_g(jx+1,jj)
          gl_sntoacc(jx,jy) = gl_sicsno(jx+1,jj)*wetol1_g(jx+1,jj)

          gl_sstacc (jx,jy) = gl_sst(jx+1,jj)*wetol1_g(jx+1,jj)+tmelt

          gl_socuacc(jx,jy) = wetol1_g(jx+1,jj)*                                     &
                        ((1.-gl_sicomo(jx+1,jj))*(gl_socu(jx+1,jj)+gl_socu(jx,jj))   &
                           + gl_sicomo(jx+1,jj) *(gl_sicu(jx+1,jj)+gl_sicu(jx,jj)))*0.5

          jym1=jj-1
          IF ( jym1 == 0) jym1 = 1

          fak = MERGE (-1,1,(lbounds_exch_tp .AND. jy == 3 .AND. jx >= ie_g/2))

          gl_socvacc (jx,jy) =  wetol1_g(jx+1,jj)*                                           &
                          ((1.-gl_sicomo(jx+1,jj))*(gl_socv(jx+1,jj)+fak*gl_socv(jx+1,jym1)) &
                             + gl_sicomo(jx+1,jj) *(gl_sicv(jx+1,jj)+fak*gl_sicv(jx+1,jym1)))*0.5

        ENDDO
      ENDDO

      CALL rotate_2_north_east(gl_socuacc,gl_socvacc,jpdim_ocei,jpdim_ocej)

      DO jy=1,jpdim_ocej
        DO jx=1,jpdim_ocei
          IF (gl_sicoacc(jx,jy) > 0.0 ) THEN
            gl_sitoacc(jx,jy) = MIN(1000.0,(gl_sitoacc(jx,jy)/gl_sicoacc(jx,jy)))
            gl_sntoacc(jx,jy) = MIN(1000.0,(gl_sntoacc(jx,jy)/gl_sicoacc(jx,jy)))
            gl_sntoacc(jx,jy) = gl_sntoacc(jx,jy)*rhoref_snow/rhoref_water
          ELSE
            gl_sitoacc(jx,jy) = 0.0
            gl_sntoacc(jx,jy) = 0.0
            gl_sicoacc(jx,jy) = 0.0
          ENDIF
        ENDDO
      ENDDO

#ifdef __synout
      WRITE(io_stdout,*) ' '
      WRITE(io_stdout,*) ' Preparation of fields for sending done.'
      FLUSH(io_stdout)
#endif
    ENDIF ! p_parallel_io

  END SUBROUTINE prep_o2a


  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE: couple_calendar
  !
  ! !INTERFACE:
  !
  SUBROUTINE couple_calendar(pdt,kout)
    !
    ! !INPUT VALUE:
    !
    INTEGER, INTENT(in)   :: kout  ! unit std output
    REAL(wp),    INTENT(in)   :: pdt   ! model time step length (secs)
    !
    ! !DESCRIPTION:
    ! - time control of coupling aspects
    !
    ! !REVISION HISTORY:
    ! 03.05.22  S. Legutke - created
    !
    ! EOP
    !----------------------------------------------------------------------

    run_date_secs = run_date_secs + INT(pdt)

#ifdef __synout
    IF (p_parallel_io) THEN
      WRITE(kout,*) ' '
      WRITE(kout,*) ' Calendar is updated.'
      FLUSH(kout)
    ENDIF
#endif

  END SUBROUTINE couple_calendar

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  chck_par
  !
  ! !INTERFACE:
  !
  SUBROUTINE chck_par
    !
    ! !DESCRIPTION:
    !
    ! - Checks coupled model control parameter against
    !   those of the calling model mpiom.
    !
    ! !REVISION HISTORY:
    ! 03.05.15  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

    !
    !-- Check whether month of restart file is ok.
    !   ------------------------------------------
    IF(lmont1 /= ig_inidate(2)) THEN
      WRITE(io_stderr,*)'  Restart file has not the right month:'
      WRITE(io_stderr,*)'      from restart file = ',lmont1
      WRITE(io_stderr,*)'      from coupler      = ',ig_inidate(2)
      CALL couple_abort (modnam,'chck_par', &
           'Restart file has not the right month',io_stderr)
    ENDIF

    a2o_freq = 0
    DO jfld = 1,nflda2o
      a2o_freq = MAX(a2o_freq,ig_def_freq(portin_id(jfld)))
    END DO
    !
    !-- ... whether model allows for exchange algorithm.
    !
    DO jfld = 1,nflda2o
      IF(a2o_freq /= ig_def_freq(portin_id(jfld))) THEN
        CALL couple_abort (modnam,'chck_par', &
             'The algorithm is not allowed ',io_stderr)
      ENDIF
    END DO

    o2a_freq = 0
    DO jfld = 1,nfldo2a
      o2a_freq = MAX(o2a_freq,ig_def_freq(portout_id(jfld)))
    END DO

    !
    !-- ...  whether model allows for exchange algorithm.
    !
    DO jfld = 1,nfldo2a
      IF(o2a_freq /= ig_def_freq(portout_id(jfld))) THEN
        CALL couple_abort (modnam,'chck_par', &
             'Algorithm not allowed when averaging in model',io_stderr)
      ENDIF
    END DO

  END SUBROUTINE chck_par

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  digest_get_Id
  !
  ! !INTERFACE:

  SUBROUTINE digest_get_Id(kout,kerr,kinfo,clfield,kdate,kcount)
    !
    ! !INPUT VALUE:
    !
    INTEGER,            INTENT(in)  :: kout    ! unit std output
    INTEGER,            INTENT(in)  :: kerr    ! unit error output
    INTEGER,            INTENT(in)  :: kinfo   ! info ID passed from psmile
    INTEGER,            INTENT(in)  :: kdate   ! date (seconds)
    INTEGER,OPTIONAL,   INTENT(in)  :: kcount  ! exchanges field count

    CHARACTER(len=8),   INTENT(in)  :: clfield ! fields symbolic name

    !
    ! !DESCRIPTION:
    !
    ! - Print info after prism_get. Abort if error.
    !
    ! !REVISION HISTORY:
    ! 03.05.15  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

#ifdef __synout
    INTEGER              :: icount = -1

    IF (PRESENT(kcount)) THEN
      icount = kcount
    ENDIF

    WRITE(kout,*)' At date (seconds) ',kdate,clfield
    IF (icount /= -1 ) THEN
      WRITE(kout,*)' (field no. ',kcount,')'
    ENDIF
#endif

    IF ( kinfo /= PRISM_Ok          .AND. &
         kinfo /= PRISM_FromRest    .AND. &
         kinfo /= PRISM_Input       .AND. &
         kinfo /= PRISM_RecvOut     .AND. &
         kinfo /= PRISM_FromRestOut .AND. &
         kinfo /= PRISM_Recvd            )  THEN
      WRITE (kout, *)' Problem with port = ',clfield
      WRITE (kout, *)' Error code number = ',kinfo
      WRITE (kout, *)' Seconds passed    = ',kdate
      CALL couple_abort (modnam,'couple_put','pb prism_get',kerr)
#ifdef __synout
    ELSEIF (kinfo == PRISM_Recvd) THEN
      WRITE(kout,*) ' was received from another model'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_Ok) THEN
      WRITE(kout,*)' was not received; no error.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_FromRest) THEN
      WRITE(kout,*) ' was received from restart file.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_Input) THEN
      WRITE(kout,*) ' was received from input file.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_RecvOut) THEN
      WRITE(kout,*) ' was received from input file and other model.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_FromRestOut) THEN
      WRITE(kout,*) ' was received from input file ' &
           ,'and written to an output file.'
      FLUSH(kout)
#endif
    ENDIF

  END SUBROUTINE digest_get_Id

  !-----------------------------------------------------------------------
  ! BOP
  !
  ! !IROUTINE:  digest_put_Id
  !
  ! !INTERFACE:

  SUBROUTINE digest_put_Id(kout,kerr,kinfo,clfield,kdate,kcount)
    !
    ! !INPUT VALUE:
    !
    INTEGER,            INTENT(in)  :: kout    ! unit std output
    INTEGER,            INTENT(in)  :: kerr    ! unit error output
    INTEGER,            INTENT(in)  :: kinfo   ! info ID passed from psmile
    INTEGER,            INTENT(in)  :: kdate   ! date (seconds)
    INTEGER, OPTIONAL,  INTENT(in)  :: kcount  ! exchanges field count

    CHARACTER(len=8),   INTENT(in)  :: clfield ! fields symbolic name

    !
    ! !DESCRIPTION:
    !
    ! - Print info after prism_put. Abort if error.
    !
    ! !REVISION HISTORY:
    ! 03.05.15  S. Legutke - created
    !
    ! EOP
    !-----------------------------------------------------------------------

#ifdef __synout
    INTEGER              :: icount = -1

    IF (PRESENT(kcount)) THEN
      icount = kcount
    ENDIF

    WRITE(kout,*)' At date (seconds) ',kdate,clfield
    IF (icount /= -1 ) THEN
      WRITE(kout,*)' (field no. ',kcount,')'
    ENDIF
#endif

    IF ( kinfo /= PRISM_Ok        .AND. &
         kinfo /= PRISM_LocTrans  .AND. &
         kinfo /= PRISM_ToRest    .AND. &
         kinfo /= PRISM_Output    .AND. &
         kinfo /= PRISM_SentOut   .AND. &
         kinfo /= PRISM_ToRestOut .AND. &
         kinfo /= PRISM_Sent            )  THEN
      WRITE (kout, *)' Problem with port = ',clfield
      WRITE (kout, *)' Error code number = ',kinfo
      WRITE (kout, *)' Seconds passed    = ',kdate
      CALL couple_abort (modnam,'couple_put','pb prism_put',kerr)
#ifdef __synout
    ELSEIF (kinfo == PRISM_Sent) THEN
      WRITE(kout,*) ' was sent to another model'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_Ok) THEN
      WRITE(kout,*)' was not sent; no error.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_LocTrans) THEN
      WRITE(kout,*) ' was used in local transformation.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_ToRest) THEN
      WRITE(kout,*) ' was written to a restart file.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_Output) THEN
      WRITE(kout,*) ' was output to a file.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_SentOut) THEN
      WRITE(kout,*) ' was sent to another model and output to a file.'
      FLUSH(kout)
    ELSEIF (kinfo == PRISM_ToRestOut) THEN
      WRITE(kout,*) ' was sent to another model and written to a restart file.'
      FLUSH(kout)
#endif
    ENDIF

  END SUBROUTINE digest_put_Id


  SUBROUTINE hdlerr(stat,string)

    !------------------------------------------------------------------------------
    !
    !  Routine to handle netcdf errors
    !
    !------------------------------------------------------------------------------

    INCLUDE 'netcdf.inc'

    !
    ! !INPUT VALUE:
    !
    INTEGER,       INTENT(in) :: stat
    CHARACTER*(*), INTENT(in) :: string

    !------------------------------------------------------------------------------

    IF (stat /= NF_NOERR) THEN
      WRITE (6,*) '--------'
      WRITE (6,*) ' ERROR:  ', string
      WRITE (6,*) '--------'
      STOP
    END IF

  END SUBROUTINE hdlerr
  !------------------------------------------------------------------------------




SUBROUTINE scatter_bounds(field,grid,message)

  REAL(wp), INTENT(out) :: field(:,:)
  CHARACTER(*), INTENT(in) :: grid, message

  ! exchange east and west halos ; not needed in tp-case due to the
  ! following call to bounds_exchange

  IF ( p_parallel_io ) THEN
    pfield_g(1,:)= pfield_g(ie_g-1,:)
    pfield_g(ie_g,:)= pfield_g(2,:)
  END IF

  CALL scatter(pfield_g, field, p_io)         ! x-stress on u-point

  IF ( lbounds_exch_tp ) THEN
    CALL bounds_exch(1,grid,field,message)
  ENDIF

END SUBROUTINE scatter_bounds

#endif /* __coupled || __prism */

END MODULE mo_couple
