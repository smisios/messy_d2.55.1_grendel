! ==============================================================================
! {%CMODEL}_tag_parameters.inc
! generated: {%TIMEDATE}
!
! this module is created automatically by imtag tool - do not edit
!
! tagging inter-configuration parameters file
!
! [Gromov, MPIC, 2007-2020]
! ==============================================================================

#define TAG

! = warning: parameters set below will be active in all cfgs, unless overriden =
!#define DEBUG

! gives lots of repetitive messages
!#define DEEPDEBUG

! rare isotopologues initialized with 0 permil rel. to standard
!#define ZERO_TEST

! rare isotopologues / fraction are initialized emptied
!#define NULL_TEST

! fast first-step initialization of tagged species in E5M1/EMAC
!#define INITFAST_E5

! initialisation for clumped O2
!#define xCIO_IEXBOOST
!#define xCIO_IEXINIT

! {$CONF_LIST} [%#define @%]

