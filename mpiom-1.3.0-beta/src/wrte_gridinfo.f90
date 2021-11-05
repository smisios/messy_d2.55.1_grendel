      SUBROUTINE WRTE_GRIDINFO(KYEARS,KMONTS,KDAYS)
!
!C**** *WRTE_GRIDINFO* - save grid information.
!C
!C     CH,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    01.10.01
!     - separate routine extracted from OLLIE (MAIN)
!     - netCDF version possible (with cond.comp. PNETCDFO)
!
!     Purpose
!     -------
!     
!
!     Method
!     -------
!     
!
!**   Interface.
!     ----------
!
!     *CALL*       *WRTE_GRIDINFO(KYEARS,KMONTS,KDAYS)*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS.h*      - std I/O logical units.
!
!**   Interface to calling routine (parameter list):
!     ----------------------------------------------
!
!     *INTEGER* *KYEARS*   - actual year.
!     *INTEGER* *KMONTS*   - actual month.
!     *INTEGER* *KDAYS*    - actual day.
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************

      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      USE MO_KIND


      INTEGER(KIND=i4) I4I1,I4I2,I4I3,I4I4,I4I5
      REAL(KIND=sp) FF_G(IE_G,JE_G)


#ifdef PNETCDFO
      INCLUDE 'netcdf.inc'

      INTEGER ncid,ncstat,ncvarid,                                      &
     &        nclonid,nclatid,nclevid,nclontoid,nclattoid
#endif
      REAL dlxp_g(ie_g,je_g), dlyp_g(ie_g,je_g)
      REAL dlxu_g(ie_g,je_g), dlyu_g(ie_g,je_g)
      REAL deuto_g(ie_g,je_g), depto_g(ie_g,je_g)

      CALL gather_arr(dlxp,dlxp_g,p_io)
      CALL gather_arr(dlyp,dlyp_g,p_io)
      CALL gather_arr(dlxu,dlxu_g,p_io)
      CALL gather_arr(dlyu,dlyu_g,p_io)
      CALL gather_arr(deuto,deuto_g,p_io)
      CALL gather_arr(depto,depto_g,p_io)

      IF(p_pe/=p_io) RETURN ! Only I/O pe does the write
!
! Append to NetCDF file : NOTE: the global attribute 'date' is set when 
!                               the file is opened

#ifdef PNETCDFO
!
! Open NetCDF file
!
      ncstat = NF_OPEN('hopc.nc',NF_NOCLOBBER, ncid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     & CALL STOP_ALL('WRTE_GRIDINFO: Problem with netCDF1')
!
! Set DEFINE mode
!
      ncstat = NF_REDEF(ncid)
      IF ( ncstat .NE. NF_NOERR )                                        &
     & CALL STOP_ALL('WRTE_GRIDINFO: Problem with define mode')
!
! Get dimension IDs
!
      ncstat = NF_DEF_DIM(ncid, 'lonto', ito, nclontoid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     &   CALL STOP_ALL('AUFW: Problem with netCDF2')
      ncstat = NF_DEF_DIM(ncid, 'latto', jto, nclattoid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     &   CALL STOP_ALL('AUFW: Problem with netCDF2')
!
! Get dimension IDs
!
      ncstat = NF_INQ_DIMID(ncid,'lon',nclonid)
      ncstat = NF_INQ_DIMID(ncid,'lat',nclatid)
      ncdims(1) = nclonid
      ncdims(2) = nclatid
      ncdims(3) = nclevid
!
! Define variables
!
      ncstat = NF_DEF_VAR(ncid,'dlxp',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 1')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 & 
     &,40, 'x-distance between two scalar grid cells')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,dlxp_g(1,1) )

      ncstat = NF_DEF_VAR(ncid,'dlyp',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 2')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,40, 'y-distance between two scalar grid cells')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,dlyp_g(1,1) )

      ncstat = NF_DEF_VAR(ncid,'dlxu',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 3')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,40, 'x-distance between two vector grid cells')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,dlxu_g(1,1) )

      ncstat = NF_DEF_VAR(ncid,'dlyu',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 4')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,40, 'y-distance between two vector grid cells')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,dlyu_g(1,1) )

      ncstat = NF_DEF_VAR(ncid,'deuto',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 5')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'
     &,32, 'water depth at vector grid cells')                          &
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,deuto_g(1,1) )


      ncstat = NF_DEF_VAR(ncid,'depto',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 6')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',5, 'meter')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,32, 'water depth at scalar grid cells')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,depto_g(1,1) )
   

      ncstat = NF_DEF_VAR(ncid,'weto',NF_DOUBLE,3,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 7')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',3, '1/0')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,35, 'water/sea mask at scalar grid cells')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,weto_g(1,1,1) )
   
      ncdims(1) = nclonid
      ncdims(2) = nclatid

      ncstat = NF_DEF_VAR(ncid,'gila',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 8')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',8, 'degree N')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,25, 'latitudes of doubled grid')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,gila_g(1,1) )

      ncstat = NF_DEF_VAR(ncid,'giph',NF_DOUBLE,2,ncdims,ncvarid)
      IF (ncstat.NE.NF_NOERR)CALL STOP_ALL( 'WRTE_GRIDINFO: STOP at 9')
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'units',8, 'degree N')     
      ncstat = NF_PUT_ATT_TEXT(ncid,ncvarid,'long_name'                 &
     &,26, 'longitudes of doubled grid')
      ncstat = NF_PUT_VAR_DOUBLE(ncid,ncvarid,giph_g(1,1) )

      ncstat = NF_CLOSE(ncid)
      IF ( ncstat .NE. NF_NOERR )                                       &
     &   CALL STOP_ALL( 'AUFW: Problem with netCDF151')

#else
!
! Write to disk (EXTRA format)
!
          I4I1=((KYEARS*10000)+(KMONTS*100)+KDAYS)
          I4I3=0
          I4I4=(IE_G*JE_G)

!          OPEN(IO_OU_DLXP,STATUS='UNKNOWN',                             &
!     &                    ACCESS='SEQUENTIAL',                          &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          I4I2=85

          WRITE(IO_OU_DLXP)I4I1,I4I2,I4I3,I4I4
          WRITE(IO_OU_DLXP)real(dlxp_g,sp)
!          CLOSE(IO_OU_DLXP)

!          OPEN(IO_OU_DLYP,STATUS='UNKNOWN',                             &
!     &                    ACCESS='SEQUENTIAL',                          &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          I4I2=86
          WRITE(IO_OU_DLYP)I4I1,I4I2,I4I3,I4I4
          WRITE(IO_OU_DLYP)real(DLYP_G,sp)
!          CLOSE(IO_OU_DLYP)


!          OPEN(IO_OU_DLXU,STATUS='UNKNOWN',                             &
!     &                    ACCESS='SEQUENTIAL',                          &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          I4I2=185
          WRITE(IO_OU_DLXU)I4I1,I4I2,I4I3,I4I4
          WRITE(IO_OU_DLXU)real(DLXU_G,sp)
!          CLOSE(IO_OU_DLXU)

!          OPEN(IO_OU_DLYU,STATUS='UNKNOWN',                             &
!     &                    ACCESS='SEQUENTIAL',                          &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          I4I2=186
          WRITE(IO_OU_DLYU)I4I1,I4I2,I4I3,I4I4
          WRITE(IO_OU_DLYU)real(DLYU_G,sp)
!          CLOSE(IO_OU_DLYU)


!          OPEN(IO_OU_DEUT,STATUS='UNKNOWN',                             &
!     &                    ACCESS='SEQUENTIAL',                          &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          I4I2=184
          WRITE(IO_OU_DEUT)I4I1,I4I2,I4I3,I4I4
          WRITE(IO_OU_DEUT)real(DEUTO_G,sp)
!          CLOSE(IO_OU_DEUT)


!          OPEN(IO_OU_DEPT,STATUS='UNKNOWN',                            &
!     &                    ACCESS='SEQUENTIAL',                         &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          I4I2=84
          WRITE(IO_OU_DEPT)I4I1,I4I2,I4I3,I4I4   
          WRITE(IO_OU_DEPT)real(DEPTO_G,sp)
!          CLOSE(IO_OU_DEPT)

!          OPEN(IO_OU_GILA,STATUS='UNKNOWN',                            &
!     &                    ACCESS='SEQUENTIAL',                         &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          I4I2=54
          I4I5=2*IE_G*2*JE_G
          WRITE(IO_OU_GILA)I4I1,I4I2,I4I3,I4I5   
          WRITE(IO_OU_GILA)real(GILA_G,sp)
!          CLOSE(IO_OU_GILA)

!          OPEN(IO_OU_GIPH,STATUS='UNKNOWN',                            &
!     &                    ACCESS='SEQUENTIAL',                         &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          I4I2=55
          I4I5=2*IE_G*2*JE_G
          WRITE(IO_OU_GIPH)I4I1,I4I2,I4I3,I4I5  
          WRITE(IO_OU_GIPH)real(GIPH_G,sp)
!          CLOSE(IO_OU_GIPH)

!          OPEN(IO_OU_WETO,STATUS='UNKNOWN',                            &
!     &                    ACCESS='SEQUENTIAL',                         &
!     &                    POSITION='APPEND',                            & 
!     &                      FORM='UNFORMATTED')
          DO K=1,KE

          DO i=1,IE_G
          do j=1,je_g
            ff_G(i,j)=real(weto_g(i,j,k),sp)
          enddo
          enddo
          I4I2=172
          I4I3=K
          WRITE(IO_OU_WETO)I4I1,I4I2,I4I3,I4I4   
          WRITE(IO_OU_WETO)ff_G
          ENDDO
!          CLOSE(IO_OU_WETO)
#endif
       RETURN
       END
