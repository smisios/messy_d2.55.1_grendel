      SUBROUTINE COMPUTE_DYN_DIFF(kpie,kpje,kpke,pddpo,psao,ptho,   &
     &                            ksave_diff,kparam)
!*******************************************************************
!
!**** *COMPUTE_DYN_DIFF* - compute the changes in the MLD.
!
!     Patrick Wetzel,        *MPIMet, HH*    24.05.04
!
!     Modified
!     --------
!     
!     Purpose
!     -------
!
!     Method
!     -------
!     -
!
!**   Interface.
!     ----------
!
!     *CALL*       *COMPUTE_DYN_DIFF(kpie,kpje,kpke)*
!
!     *MODULES*     *mo_carbch* - ocean/sediment tracer arrays.
!     *MODULES*     *mo_control_bgc*  - std I/O logical units.
!     *MODULES*     *mo_dynamic
!     *MODULES*     *mo_param1_bgc 
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st REAL :: of model grid.
!     *INTEGER* *kpje*    - 2nd REAL :: of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) REAL :: of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd REAL ::) [m].
!     *REAL*    *psao*    - salinity [psu.].
!     *REAL*    *ptho*    - potential temperature (3rd REAL).
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************


      USE mo_carbch
      USE mo_bgcmean
      USE mo_dynamic
      use mo_param1_bgc 
      USE mo_control_bgc

implicit none    

      INTEGER :: kpie,kpje,kpke,i,j,k,l
      INTEGER :: ksave_diff, kparam
      REAL :: pddpo(kpie,kpje,kpke)
      REAL psao(kpie,kpje,kpke)
      REAL ptho(kpie,kpje,kpke)

      WRITE(io_stdo_bgc,*) 'COMPUTE_DYN_DIFF mit:', ksave_diff, kparam

      if (ksave_diff.eq.kdyndiff) then
!-------------------------------------------------------------------------  
    
!$OMP parallel DO private (i,j,k)
!!$OMP+ firstprivate (bgcdyntmp)
!!$OMP+ lastprivate (bgcdyntmp)
        
        DO j=1,kpje
      DO k=1,4
        DO i=1,kpie
       if (pddpo(i,j,k).gt.0.5) then
       
        bgcdyntmp(i,j,kdphyto  ) = bgcdyntmp(i,j,kdphyto  )          &
      &                            - ocetra(i,j,k,iphy  )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdgrazer ) = bgcdyntmp(i,j,kdgrazer )          &
      &                            - ocetra(i,j,k,izoo )*pddpo(i,j,k)
        bgcdyntmp(i,j,kddoc    ) = bgcdyntmp(i,j,kddoc    )          &
      &                            - ocetra(i,j,k,idoc    )*pddpo(i,j,k)
        bgcdyntmp(i,j,kddic    ) = bgcdyntmp(i,j,kddic    )          &
      &                            - ocetra(i,j,k,isco212    )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdphosph ) = bgcdyntmp(i,j,kdphosph )          &
      &                            - ocetra(i,j,k,iphosph )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdoxygen ) = bgcdyntmp(i,j,kdoxygen )          &
      &                            - ocetra(i,j,k,ioxygen )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdiron   ) = bgcdyntmp(i,j,kdiron   )          &
      &                            - ocetra(i,j,k,iiron   )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdano3   ) = bgcdyntmp(i,j,kdano3   )          &
      &                            - ocetra(i,j,k,iano3   )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdalkali ) = bgcdyntmp(i,j,kdalkali )          &
      &                            - ocetra(i,j,k,ialkali )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdsilica ) = bgcdyntmp(i,j,kdsilica )          &
      &                            - ocetra(i,j,k,isilica )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdtemp   ) = bgcdyntmp(i,j,kdtemp   )          &
      &                            - ptho(i,j,k)           *pddpo(i,j,k)
        bgcdyntmp(i,j,kdsal    ) = bgcdyntmp(i,j,kdsal    )          &
      &                            - psao(i,j,k)           *pddpo(i,j,k)        
#ifdef PANTHROPOCO2                                                      
        bgcdyntmp(i,j,kddic_ant) = bgcdyntmp(i,j,kddic_ant)          &
      &                            - ocetra(i,j,k,isco2_ant)*pddpo(i,j,k)
        bgcdyntmp(i,j,kdalk_ant) = bgcdyntmp(i,j,kdalk_ant)          &
      &                            - ocetra(i,j,k,ialk_ant)*pddpo(i,j,k)
#endif  /*  PANTHROPOCO2 */ 
        
       endif
      ENDDO
        ENDDO
        ENDDO      
      
!$OMP parallel DO private (i,j)
        DO j=1,kpje
        DO i=1,kpie
       if (pddpo(i,j,1).gt.0.5) then      
       
        bgcdyn(i,j,kdphyto  ,kparam) = bgcdyn(i,j,kdphyto  ,kparam) &
      &                                                + bgcdyntmp(i,j,kdphyto  )
        bgcdyn(i,j,kdgrazer ,kparam) = bgcdyn(i,j,kdgrazer ,kparam) &
      &                                                + bgcdyntmp(i,j,kdgrazer )
        bgcdyn(i,j,kddoc    ,kparam) = bgcdyn(i,j,kddoc    ,kparam) &
      &                                                + bgcdyntmp(i,j,kddoc    )
        bgcdyn(i,j,kddic    ,kparam) = bgcdyn(i,j,kddic    ,kparam) &
      &                                           + bgcdyntmp(i,j,kddic    )
        bgcdyn(i,j,kdphosph ,kparam) = bgcdyn(i,j,kdphosph ,kparam) &
      &                                                + bgcdyntmp(i,j,kdphosph )
        bgcdyn(i,j,kdoxygen ,kparam) = bgcdyn(i,j,kdoxygen ,kparam) &
      &                                                + bgcdyntmp(i,j,kdoxygen )
        bgcdyn(i,j,kdiron   ,kparam) = bgcdyn(i,j,kdiron   ,kparam) &
      &                                                + bgcdyntmp(i,j,kdiron   )
        bgcdyn(i,j,kdano3   ,kparam) = bgcdyn(i,j,kdano3   ,kparam) &
      &                                                 + bgcdyntmp(i,j,kdano3   ) 
        bgcdyn(i,j,kdalkali ,kparam) = bgcdyn(i,j,kdalkali ,kparam) &
      &                                                + bgcdyntmp(i,j,kdalkali ) 
        bgcdyn(i,j,kdsilica ,kparam) = bgcdyn(i,j,kdsilica ,kparam) &
      &                                                 + bgcdyntmp(i,j,kdsilica ) 
        bgcdyn(i,j,kdtemp   ,kparam) = bgcdyn(i,j,kdtemp   ,kparam) &
      &                                                + bgcdyntmp(i,j,kdtemp   ) 
        bgcdyn(i,j,kdsal    ,kparam) = bgcdyn(i,j,kdsal    ,kparam) &
      &                                              + bgcdyntmp(i,j,kdsal    ) 
#ifdef PANTHROPOCO2                                                
        bgcdyn(i,j,kddic_ant,kparam) = bgcdyn(i,j,kddic_ant,kparam) &
      &                                              + bgcdyntmp(i,j,kddic_ant) 
        bgcdyn(i,j,kdalk_ant,kparam) = bgcdyn(i,j,kdalk_ant,kparam) &
      &                                              + bgcdyntmp(i,j,kdalk_ant) 
#endif  /*  PANTHROPOCO2 */ 
         endif
        ENDDO
        ENDDO

      endif 
      

      if (ksave_diff.eq.kdynsave) then
!-------------------------------------------------------------------------      
!$OMP parallel DO private (i,j,l)
        DO l=1,nbgcdyn
        DO j=1,kpje
        DO i=1,kpie
          bgcdyntmp(i,j,l) = 0.
        ENDDO
        ENDDO
        ENDDO      

!$OMP parallel DO private (i,j,k)
!!$OMP+ firstprivate (bgcdyntmp)
!!$OMP+ lastprivate (bgcdyntmp)

        DO j=1,kpje
      DO k=1,4
        DO i=1,kpie
       if (pddpo(i,j,k).gt.0.5) then
        
        bgcdyntmp(i,j,kdphyto  ) = bgcdyntmp(i,j,kdphyto  )          &
      &                            + ocetra(i,j,k,iphy    )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdgrazer ) = bgcdyntmp(i,j,kdgrazer )          &
      &                            + ocetra(i,j,k,izoo    )*pddpo(i,j,k)
        bgcdyntmp(i,j,kddoc    ) = bgcdyntmp(i,j,kddoc    )          &
      &                            + ocetra(i,j,k,idoc    )*pddpo(i,j,k)
        bgcdyntmp(i,j,kddic    ) = bgcdyntmp(i,j,kddic    )          &
      &                            + ocetra(i,j,k,isco212 )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdphosph ) = bgcdyntmp(i,j,kdphosph )          &
      &                            + ocetra(i,j,k,iphosph )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdoxygen ) = bgcdyntmp(i,j,kdoxygen )          &
      &                            + ocetra(i,j,k,ioxygen )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdiron   ) = bgcdyntmp(i,j,kdiron   )          &
      &                            + ocetra(i,j,k,iiron   )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdano3   ) = bgcdyntmp(i,j,kdano3   )          &
      &                            + ocetra(i,j,k,iano3   )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdalkali ) = bgcdyntmp(i,j,kdalkali )          &
      &                            + ocetra(i,j,k,ialkali )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdsilica ) = bgcdyntmp(i,j,kdsilica )          &
      &                            + ocetra(i,j,k,isilica )*pddpo(i,j,k)
        bgcdyntmp(i,j,kdtemp   ) = bgcdyntmp(i,j,kdtemp   )          &
      &                            + ptho(i,j,k)           *pddpo(i,j,k)
        bgcdyntmp(i,j,kdsal    ) = bgcdyntmp(i,j,kdsal    )          &
      &                            + psao(i,j,k)           *pddpo(i,j,k)        
#ifdef PANTHROPOCO2                                                      
        bgcdyntmp(i,j,kddic_ant) = bgcdyntmp(i,j,kddic_ant)          &
      &                            + ocetra(i,j,k,isco2_ant)*pddpo(i,j,k)
        bgcdyntmp(i,j,kdalk_ant) = bgcdyntmp(i,j,kdalk_ant)          &
      &                            + ocetra(i,j,k,ialk_ant)*pddpo(i,j,k)
#endif  /*  PANTHROPOCO2 */ 
          
       endif
      ENDDO 
        ENDDO
        ENDDO
      
      endif 
      
      RETURN
      END
