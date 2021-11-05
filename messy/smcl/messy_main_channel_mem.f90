! **********************************************************************
MODULE messy_main_channel_mem
! **********************************************************************

  IMPLICIT NONE
  PUBLIC
  SAVE

  ! This module stores the information about the CHANNEL extension
  ! for different DOMAINS, which allows identical channel name &
  ! channel object name pairs with an additional tag, namely the
  ! DOMAIN.
  ! The defaults are set such that basemodels, which are not aware of /
  ! do not requires this extension do not need to call the
  ! corresponding setup routines.
  !
  ! Author: P. Joeckel, DLR, August 2019

  ! switch indicating, if the basemodel requires the domain extension
  ! (to be set once in the initialisation phase: main_channel_setup)
  LOGICAL :: l_dom = .FALSE.

  ! max. number of domains the basemodel can handle
  ! (to be set once in the initialisation phase: messy_setup)
  INTEGER :: n_dom_max = 1

  ! actual number of active domains
  ! (to be set once in the initialisation phase: messy_setup)
  INTEGER :: n_dom = 1

  ! abstract domain, which is not bound to a specific region or similar
  ! Note: make sure that this value is the same in
  !       messy_main_tools::split_name_domain:
  INTEGER, PARAMETER :: dom_unbound = 0

  ! id of domain within loop over domains
  ! (to be set within loop over domains: main_channel_set_domain)
  INTEGER :: dom_current = dom_unbound

  ! id of domain within loop over domains, which is however not
  ! the unbound domain ( = MAX(1, dom_current) )
  ! Note: This is required to avoid (basemodel) preprocessor
  !       directives or 'IF (l_dom)'-statements to distinguish
  !       domain-aware from non-domain-aware basemodels 
  !       in (diagnostic) submodel SMILs.
  !       
  INTEGER :: dom_curid = 1
  
! **********************************************************************
END MODULE messy_main_channel_mem
! **********************************************************************
