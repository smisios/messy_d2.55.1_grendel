MODULE messy_main_data_wiso_bi

  USE messy_main_tools,          ONLY: PTR_4D_ARRAY, PTR_5D_ARRAY &
                                     , PTR_3D_ARRAY, PTR_2D_ARRAY
  USE messy_main_constants_mem,  ONLY: dp
  USE messy_main_tools_wiso,     ONLY: kphase, mwiso, i_vap, i_liq, i_ice
  
  IMPLICIT NONE
  PUBLIC
  SAVE

  ! TODO: it needs to be checkes what really needs to be shared (to remain
  !       here) and what can be moved back to H2OISO ...
  
  ! ------------------------------------------------------------------
  ! This module contains the parameters, switches and pointers for
  ! the expansion of the (atmospheric) hydrological cycle by its
  ! isotopologues, i.e. HDO and H2(^18)O.
  ! These data need to be shared between the physical submodels that
  ! describe processes of the (atmospheric) hydrological cylcle.
  ! ------------------------------------------------------------------

  ! ------------------------------------------------------------------
  ! GLOBAL SWITCH: WISO is ON or OFF (set by submodel H2OISO)
  ! ------------------------------------------------------------------
  LOGICAL :: l_wiso = .FALSE.

  ! SPECIAL SWITCHES (TRANSFERRED FROM &CPL NAMELIST OF H2OISO)
  LOGICAL :: l_wiso_nocloud_dd = .FALSE.
  
  ! pointers to tracer memory
  TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER :: isotracm1 => NULL()
  TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER :: isotracte => NULL()
  
  ! tracer indices of water isotopologue tracers
  INTEGER, DIMENSION(kphase,mwiso):: idiso = 0
  
  ! new 2D diagnostic channel objects 
  INTEGER, PARAMETER :: i_ws     = 1 ! vdiff/surf
  INTEGER, PARAMETER :: i_sn     = 2 ! vdiff/surf
  INTEGER, PARAMETER :: i_wl     = 3 ! vdiff/surf
  INTEGER, PARAMETER :: i_gld    = 4 ! surf
  INTEGER, PARAMETER :: i_evapl  = 5 ! vdiff/surf
  INTEGER, PARAMETER :: i_evapot = 6 ! vdiff/surf
  INTEGER, PARAMETER :: i_rsfl   = 7 ! cloud/surf
  INTEGER, PARAMETER :: i_rsfc   = 8 ! convect/surf
  INTEGER, PARAMETER :: i_ssfl   = 9 ! cloud/surf
  INTEGER, PARAMETER :: i_ssfc   = 10! convect/surf
  INTEGER, PARAMETER :: i_snc    = 11! surf
  INTEGER, PARAMETER :: i_aprs   = 12! convect/cloud
  INTEGER, PARAMETER :: i_aprc   = 13! convect
  INTEGER, PARAMETER :: i_aprl   = 14! cloud
  INTEGER, PARAMETER :: i_qvi    = 15! cloud
  INTEGER, PARAMETER :: i_xlvi   = 16! cloud
  INTEGER, PARAMETER :: i_xivi   = 17! cloud
  INTEGER, PARAMETER :: i_alac   = 18! surf
  INTEGER, PARAMETER :: i_runoff = 19! surf
  INTEGER, PARAMETER :: i_snmel  = 20! surf
  INTEGER, PARAMETER :: i_apmegl = 21! surf
  INTEGER, PARAMETER :: i_drain  = 22! surf
  INTEGER, PARAMETER :: i_snacl  = 23! surf
  INTEGER, PARAMETER :: i_rogl   = 24! surf
  INTEGER, PARAMETER :: i_evap   = 25! vdiff
  INTEGER, PARAMETER :: i_evapw  = 26! vdiff
  INTEGER, PARAMETER :: i_evapi  = 27! vdiff
  INTEGER, PARAMETER :: i_evaplac= 28! vdiff
  INTEGER, PARAMETER :: i_evapwac= 29! vdiff
  INTEGER, PARAMETER :: i_evapiac= 30! vdiff
  INTEGER, PARAMETER :: i_snglac = 31! vdiff/surf
  INTEGER, PARAMETER :: i_daprs  = 32! global_end
  INTEGER, PARAMETER :: i_daprc  = 33! global_end
  INTEGER, PARAMETER :: i_daprl  = 34! global_end
  INTEGER, PARAMETER :: i_daprt  = 35! global_end
  INTEGER, PARAMETER :: i_2nmax  = 35 ! max. number of 2-D surf pointers
  CHARACTER(LEN=8), DIMENSION(i_2nmax), PARAMETER :: varname_2dh = &
       (/'ws      ', 'sn      ', 'wl      ', 'gld     ', 'evapl   ', 'evapot  '&
       , 'rsfl    ', 'rsfc    ', 'ssfl    ', 'ssfc    ', 'snc     ', 'aprs    '&
       , 'aprc    ', 'aprl    ', 'qvi     ', 'xlvi    ', 'xivi    ', 'alac    '&
       , 'runoff  ', 'snmel   ', 'apmegl  ', 'drain   ', 'snacl   ', 'rogl    '&
       , 'evap    ', 'evapw   ', 'evapi   ', 'evaplac ', 'evapwac ', 'evapiac '&
       , 'snglac  ', 'daprs   ', 'daprc   ', 'daprl   ', 'daprt   '/)
  TYPE(PTR_5D_ARRAY), DIMENSION(:), POINTER :: pwiso_2dh_5d => NULL()
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: pwiso_2dh => NULL()
  TYPE(PTR_2D_ARRAY), DIMENSION(i_2nmax) :: pwiso2
  
  ! new 3D diagnostic channel objects 
  INTEGER, PARAMETER :: i_xtec  = 1 ! cloud/convect
  INTEGER, PARAMETER :: i_qtec  = 2 ! cloud/convect
  INTEGER, PARAMETER :: i_numq  = 3
  INTEGER, PARAMETER :: i_numxl = 4
  INTEGER, PARAMETER :: i_numxi = 5
  INTEGER, PARAMETER :: i_dvap  = 6 ! in global_end for diagnosis 
  INTEGER, PARAMETER :: i_dliq  = 7 ! in global_end for diagnosis 
  INTEGER, PARAMETER :: i_dice  = 8 ! in global_end for diagnosis 
  INTEGER, PARAMETER :: i_dtot  = 9 ! in global_end for diagnosis 
  INTEGER, PARAMETER :: i_3nmax = 9
  CHARACTER(LEN=8), DIMENSION(i_3nmax), PARAMETER :: varname_3d = &
       (/'xtec    ','qtec    ','numq    ','numxl   ','numxi   ' &
       ,'dvap    ','dliq    ','dice    ','dtot    '/)
  TYPE(PTR_5D_ARRAY), DIMENSION(:), POINTER :: pwiso_3d_5d => NULL()
  TYPE(PTR_4D_ARRAY), DIMENSION(:), POINTER :: pwiso_3d => NULL()
  TYPE(PTR_3D_ARRAY), DIMENSION(i_3nmax)    :: pwiso3

  ! set with get_channel_object in init-coupling of H2OISO
  REAL(dp), DIMENSION(:,:), POINTER :: sn       => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: wl       => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: snc      => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: gld      => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: tsl      => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: grndcapc => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: orostd   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: pwisosw_d_HHO   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: pwisosw_d_HH18O => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: pwisosw_d_HDO   => NULL()

  ! surface reservoirs needed for steady-state restart
  REAL(dp), DIMENSION(:,:), POINTER :: ws_HHO   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: wl_HHO   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: sn_HHO   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: ws_HH18O => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: wl_HH18O => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: sn_HH18O => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: ws_HDO   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: wl_HDO   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: sn_HDO   => NULL()

  ! needed for cloud
  REAL(dp), POINTER, DIMENSION(:,:,:)   :: pqtec      => NULL()
  REAL(dp), POINTER, DIMENSION(:,:,:)   :: pvdiffp    => NULL()

  REAL(dp), POINTER, DIMENSION(:,:) :: wet_tmp => NULL()
  REAL(dp), POINTER, DIMENSION(:,:) :: evapot_2d => NULL()

  ! needed for surf but from vdiff (hence globally)
  ! because it is changed in surf
  REAL(DP), POINTER, DIMENSION(:):: pwl_surf => NULL()
  REAL(DP), POINTER, DIMENSION(:):: psn_surf => NULL()
  REAL(DP), POINTER, DIMENSION(:):: pws_surf => NULL()

  REAL(DP), POINTER, DIMENSION(:):: tsl_vdiff=> NULL()
  REAL(DP), POINTER, DIMENSION(:):: tsl_surf => NULL()
  REAL(DP), POINTER, DIMENSION(:):: grndcapc_surf => NULL()

  ! POINTERS FOR WITHIN LOCAL LOOPS
  REAL(dp), DIMENSION(:,:,:), POINTER :: pwisoqm1  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: pwisoxlm1 => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: pwisoxim1 => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: pwisoqte  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: pwisoxlte => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: pwisoxite => NULL()
  !
  REAL(dp), DIMENSION(:,:,:), POINTER :: pwiso_xtec => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: pwiso_qtec => NULL()
  !
  REAL(dp), DIMENSION(:,:),   POINTER :: pwiso_aprl  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: pwiso_qvi   => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: pwiso_xlvi  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: pwiso_xivi  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: pwiso_ssfl  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: pwiso_rsfl  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: pwiso_aprs  => NULL()
  
CONTAINS

  SUBROUTINE main_data_wiso_local_start(jrow)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN) :: jrow
    ! LOCAL
    INTEGER :: i
    
    IF (.NOT. l_wiso) RETURN

    pwisoqm1  => isotracm1(i_vap)%ptr(:,:,:,jrow)
    pwisoxlm1 => isotracm1(i_liq)%ptr(:,:,:,jrow)
    pwisoxim1 => isotracm1(i_ice)%ptr(:,:,:,jrow)
    pwisoqte  => isotracte(i_vap)%ptr(:,:,:,jrow)
    pwisoxlte => isotracte(i_liq)%ptr(:,:,:,jrow)
    pwisoxite => isotracte(i_ice)%ptr(:,:,:,jrow)

    DO i=1, i_3nmax
       pwiso3(i)%ptr => pwiso_3d(i)%ptr(:,:,:,jrow)
    END DO
    
    DO i=1, i_2nmax
       pwiso2(i)%ptr => pwiso_2dh(i)%ptr(:,:,jrow)
    END DO
                
  END SUBROUTINE main_data_wiso_local_start
  
END MODULE messy_main_data_wiso_bi
