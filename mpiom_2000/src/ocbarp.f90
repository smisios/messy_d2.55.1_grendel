SUBROUTINE OCBARP
  USE MO_PARAM1
  USE MO_COMMO1
  USE MO_UNITS
  IMPLICIT NONE

  IF ( .NOT. lwith_barotropic_stokes_drift ) THEN

    UZO(:,:)=USO(:,:)
    VZE(:,:)=VSE(:,:)

  ENDIF

!--------------------------------------------------------------------
!
!     B)
!
!     WINDSTRESS INPUT ON BAROCLINIC VELOCITIES
!
!=====================================================================
!
!     C)

  VSE(:,:)=VSE(:,:)*AMSUE(:,:,1)
  VZE(:,:)=VZE(:,:)*AMSUE(:,:,1)
  USO(:,:)=USO(:,:)*AMSUO(:,:,1)
  UZO(:,:)=UZO(:,:)*AMSUO(:,:,1)


END SUBROUTINE OCBARP
