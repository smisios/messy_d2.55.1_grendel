!***************************************************************************
MODULE messy_rad_long
!***************************************************************************

  !***************************************************************************
  !                                                                          *
  !                RRTM :  RAPID RADIATIVE TRANSFER MODEL                    *
  !                                                                          *
  !             ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                 *
  !                        840 MEMORIAL DRIVE                                *
  !                        CAMBRIDGE, MA 02139                               *
  !                                                                          *
  !                           ELI J. MLAWER                                  *
  !                         STEVEN J. TAUBMAN~                               *
  !                         SHEPARD A. CLOUGH                                *
  !                        ~currently at GFDL                                *
  !                       email:  mlawer@aer.com                             *
  !                                                                          *
  !        The authors wish to acknowledge the contributions of the          *
  !        following people:  Patrick D. Brown, Michael J. Iacono,           *
  !        Ronald E. Farren, Luke Chen, Robert Bergstrom.                    *
  !                                                                          *
  !***************************************************************************
  !     Reformatted for F90 by JJMorcrette, ECMWF, 980714                    *
  !***************************************************************************
  !
  ! *** This version of RRTM has been altered to interface with either
  !     the ECMWF numerical weather prediction model or the ECMWF column 
  !     radiation model (ECRT) package. 
  
  !     Revised, April, 1997;  Michael J. Iacono, AER, Inc.
  !          - initial implementation of RRTM in ECRT code
  !     Revised, June, 1999;  Michael J. Iacono and Eli J. Mlawer, AER, Inc.
  !          - to implement generalized maximum/random cloud overlap
  !     Bug-Fix, November 1999: Michael J. Iacono, AER, Inc.
  !         and  December 1999: JJMorcrette  initializations within cloud 
  !                             overlap
  !     Modified for ECHAM5, April 2000, Marco A. Giorgetta, MPI
  !     Update to ECMWF-Cy23R1, Dec2000, Marco A. Giorgetta, MPI
  !     Shift longitude loop to rrtm subroutines, Jan2001, 
  !      Marco A. Giorgetta, MPI
  
  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL(DP), parameter  :: BPADE=1._DP/0.278_DP  ! Pade constant (?)

  INTEGER, PARAMETER, PUBLIC :: JPBAND = 16  ! number of LW bands
  INTEGER, PARAMETER, PUBLIC :: JPXSEC = 4   ! number of gases in WKL
  INTEGER, PARAMETER, PUBLIC :: JPINPX = 35  ! number of tracer cross sections

  ! =========================================================================
  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA1* - RRTM COEFFICIENTS FOR INTERVAL 1
  !     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG1_loc  = 8

  REAL(DP):: FRACREFA_1(NG1_loc)  , FRACREFB_1(NG1_loc)
  REAL(DP):: ABSA_1(65,NG1_loc)
  REAL(DP):: ABSB_1(235,NG1_loc)
  REAL(DP):: SELFREF_1(10,NG1_loc), FORREF_1(NG1_loc)
!!CDIR DUPLICATE(ABSA_1,256)
!!CDIR DUPLICATE(ABSB_1,256)
!!CDIR DUPLICATE(SELFREF_1,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL
  ! FORREF  : REAL
  ! KA      : REAL     
  ! KB      : REAL     
  ! SELFREF : REAL     
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA2* - RRTM COEFFICIENTS FOR INTERVAL 2
  !     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG2_loc  = 14
  
  !     The ith set of reference fractions are from the ith reference
  !     pressure level.
  REAL(DP):: FRACREFA_2(NG2_loc,13), FRACREFB_2(NG2_loc), REFPARAM_2(13)
  REAL(DP):: ABSA_2(65,NG2_loc)
  REAL(DP):: ABSB_2(235,NG2_loc)
  REAL(DP):: SELFREF_2(10,NG2_loc), FORREF_2(NG2_loc)
!!CDIR DUPLICATE(ABSA_2,256)
!!CDIR DUPLICATE(ABSB_2,256)
!!CDIR DUPLICATE(SELFREF_2,256)
!!CDIR DUPLICATE(FRACREFA_2,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL
  ! REFPARAM: REAL
  ! KA      : REAL     
  ! KB      : REAL     
  ! SELFREF : REAL
  ! FORREF  : REAL     
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA3* - RRTM COEFFICIENTS FOR INTERVAL 3
  !     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)
  !     -----------------------------------------------------------------
  
  INTEGER, PARAMETER :: NG3_loc  = 16

  REAL(DP):: FRACREFA_3(NG3_loc,10) ,FRACREFB_3(NG3_loc,5)
  
  REAL(DP), DIMENSION(16) :: FORREF_3
  REAL(DP), DIMENSION(16) :: ABSN2OA_3
  REAL(DP), DIMENSION(16) :: ABSN2OB_3
  REAL(DP), DIMENSION(10) :: ETAREF_3
  REAL(DP), DIMENSION(59) :: H2OREF_3
  REAL(DP), DIMENSION(59) :: N2OREF_3
  REAL(DP), DIMENSION(59) :: CO2REF_3
  
  REAL(DP):: ABSA_3(650,NG3_loc)
  REAL(DP):: ABSB_3(1175,NG3_loc)
  REAL(DP):: SELFREF_3(10,NG3_loc)
  REAL(DP):: STRRAT_3

  !     ------------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14

  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! ABSN2OA : REAL
  ! ABSN2OB : REAL
  ! CO2REF  : REAL
  ! ETAREF  : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL
  ! H2OREF  : REAL
  ! KA      : REAL     
  ! KB      : REAL     
  ! N2OREF  : REAL
  ! SELFREF : REAL     
  ! STRRAT  : REAL
  !     -----------------------------------------------------------------
!!CDIR DUPLICATE(ABSA_3,256)
!!CDIR DUPLICATE(ABSB_3,256)
!!CDIR DUPLICATE(SELFREF_3,256)
!!CDIR DUPLICATE(FRACREFA_3,256)
!!CDIR DUPLICATE(FRACREFB_3,256)
!!CDIR DUPLICATE(H2OREF_3,256)
!!CDIR DUPLICATE(N2OREF_3,256)
!!CDIR DUPLICATE(CO2REF_3,256)
!!CDIR DUPLICATE(ETAREF_3,256)

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA4* - RRTM COEFFICIENTS FOR INTERVAL 5
  !     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG4_loc  = 14

  REAL(DP):: FRACREFA_4(NG4_loc,9)  ,FRACREFB_4(NG4_loc,6)
  REAL(DP):: ABSA_4(585,NG4_loc)
  REAL(DP):: ABSB_4(1410,NG4_loc)
  REAL(DP):: SELFREF_4(10,NG4_loc)
  REAL(DP):: STRRAT1_4
  REAL(DP):: STRRAT2_4
!!CDIR DUPLICATE(ABSA_4,256)
!!CDIR DUPLICATE(ABSB_4,256)
!!CDIR DUPLICATE(SELFREF_4,256)
!!CDIR DUPLICATE(FRACREFA_4,256)
!!CDIR DUPLICATE(FRACREFB_4,256)

  !     ------------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL
  ! KA      : REAL     
  ! KB      : REAL     
  ! SELFREF : REAL     
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA5* - RRTM COEFFICIENTS FOR INTERVAL 5
  !     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG5_loc  = 16

  REAL(DP):: FRACREFA_5(NG5_loc,9) ,FRACREFB_5(NG5_loc,5)

  REAL(DP), DIMENSION(NG5_loc) :: CCL4_5

  REAL(DP):: ABSA_5(585,NG5_loc)
  REAL(DP):: ABSB_5(1175,NG5_loc)
  REAL(DP):: SELFREF_5(10,NG5_loc)
  REAL(DP):: STRRAT1_5
  REAL(DP):: STRRAT2_5
!!CDIR DUPLICATE(ABSA_5,256)
!!CDIR DUPLICATE(ABSB_5,256)
!!CDIR DUPLICATE(SELFREF_5,256)
!!CDIR DUPLICATE(FRACREFA_5,256)
!!CDIR DUPLICATE(FRACREFB_5,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! CCL4    : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL
  ! KA      : REAL     
  ! KB      : REAL     
  ! SELFREF : REAL     
  ! STRRAT1 : REAL
  ! STRRAT2 : REAL    
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA6* - RRTM COEFFICIENTS FOR INTERVAL 6
  !     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG6_loc  = 8
  
  REAL(DP), DIMENSION(NG6_loc) :: FRACREFA_6
  
  REAL(DP), DIMENSION(NG6_loc) :: CFC11ADJ_6
  REAL(DP), DIMENSION(NG6_loc) :: CFC12_6
  REAL(DP), DIMENSION(NG6_loc) :: ABSCO2_6

  REAL(DP):: ABSA_6(65,NG6_loc)
  REAL(DP):: SELFREF_6(10,NG6_loc)
!!CDIR DUPLICATE(ABSA_6,256)
!!CDIR DUPLICATE(SELFREF_6,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE *
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSCO2  : REAL 
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! FRACREFA: REAL    
  ! CFC11ADJ: REAL
  ! CFC12   : REAL
  ! KA      : REAL     
  ! SELFREF : REAL     
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA7* - RRTM COEFFICIENTS FOR INTERVAL 7
  !     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG7_loc  = 12

  REAL(DP):: FRACREFA_7(NG7_loc,9)

  REAL(DP), DIMENSION(NG7_loc) :: FRACREFB_7
  REAL(DP), DIMENSION(NG7_loc) :: ABSCO2_7
  REAL(DP):: ABSA_7(585,NG7_loc)
  REAL(DP):: ABSB_7(235,NG7_loc)
  REAL(DP):: SELFREF_7(10,NG7_loc)
!!CDIR DUPLICATE(ABSA_7,256)
!!CDIR DUPLICATE(ABSB_7,256)
!!CDIR DUPLICATE(SELFREF_7,256)
!!CDIR DUPLICATE(FRACREFA_7,256)
!!CDIR DUPLICATE(FRACREFB_7,256)

  REAL(DP):: STRRAT_7

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE *
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL 
  ! ABSCO2  : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL    
  ! KA      : REAL     
  ! KB      : REAL     
  ! SELFREF : REAL  
  ! STRRAT  : REAL   
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA8* - RRTM COEFFICIENTS FOR INTERVAL 8
  !     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)
  !     -----------------------------------------------------------------
  
  INTEGER, PARAMETER :: NG8_loc  = 8

  REAL(DP), DIMENSION(NG8_loc) :: FRACREFA_8
  REAL(DP), DIMENSION(NG8_loc) :: FRACREFB_8
  REAL(DP), DIMENSION(NG8_loc) :: CFC12_8
  REAL(DP), DIMENSION(NG8_loc) :: CFC22ADJ_8
  REAL(DP), DIMENSION(NG8_loc) :: ABSCO2A_8
  REAL(DP), DIMENSION(NG8_loc) :: ABSCO2B_8
  REAL(DP), DIMENSION(NG8_loc) :: ABSN2OA_8
  REAL(DP), DIMENSION(NG8_loc) :: ABSN2OB_8
  REAL(DP), DIMENSION(59)  :: H2OREF_8
  REAL(DP), DIMENSION(59)  :: N2OREF_8
  REAL(DP), DIMENSION(59)  :: O3REF_8

  REAL(DP):: ABSA_8(35,NG8_loc)
  REAL(DP):: ABSB_8(265,NG8_loc)
  REAL(DP):: SELFREF_8(10,NG8_loc)

!!CDIR DUPLICATE(ABSA_8,256)
!!CDIR DUPLICATE(ABSB_8,256)
!!CDIR DUPLICATE(SELFREF_8,256)
!!CDIR DUPLICATE(H2OREF_8,256)
!!CDIR DUPLICATE(N2OREF_8,256)
!!CDIR DUPLICATE(O3REF_8,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE *
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! ABSCO2A : REAL     
  ! ABSCO2B : REAL     
  ! ABSN2OA : REAL     
  ! ABSN2OB : REAL 
  ! CFC12   : REAL     
  ! CFC22ADJ: REAL     
  ! FRACREFA: REAL    
  ! FRACREFB: REAL    
  ! H2OREF  : REAL    
  ! KA      : REAL     
  ! KB      : REAL     
  ! N2OREF  : REAL    
  ! O3REF   : REAL    
  ! SELFREF : REAL     
  !     -----------------------------------------------------------------
  
  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA9* - RRTM COEFFICIENTS FOR INTERVAL 9
  !     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG9_loc  = 12

  REAL(DP):: FRACREFA_9(NG9_loc,9)

  REAL(DP), DIMENSION(NG9_loc) :: FRACREFB_9
  REAL(DP), DIMENSION(13) :: N2OREF_9
  REAL(DP), DIMENSION(13) :: H2OREF_9
  REAL(DP), DIMENSION(13) :: CH4REF_9
  REAL(DP), DIMENSION(11) :: ETAREF_9
  ! 36 = 3*NG9_loc      
  REAL(DP), DIMENSION(36) :: ABSN2O_9

  REAL(DP):: ABSA_9(715,NG9_loc)
  REAL(DP):: ABSB_9(235,NG9_loc)
  REAL(DP):: SELFREF_9(10,NG9_loc)
  REAL(DP):: STRRAT_9

!!CDIR DUPLICATE(ABSA_9,256)
!!CDIR DUPLICATE(ABSB_9,256)
!!CDIR DUPLICATE(SELFREF_9,256)
!!CDIR DUPLICATE(ABSN2O_9,256)
!!CDIR DUPLICATE(FRACREFA_9,256)
!!CDIR DUPLICATE(N2OREF_9,256)
!!CDIR DUPLICATE(H2OREF_9,256)
!!CDIR DUPLICATE(CH4REF_9,256)
!!CDIR DUPLICATE(ETAREF_9,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! ABSN2O  : REAL    
  ! CH4REF  : REAL
  ! ETAREF  : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL
  ! H2OREF  : REAL
  ! N2OREF  : REAL
  ! KA      : REAL     
  ! KB      : REAL     
  ! SELFREF : REAL     
  ! STRRAT  : REAL
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA14* - RRTM COEFFICIENTS FOR INTERVAL 10
  !     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG10_loc = 6

  REAL(DP), DIMENSION(NG10_loc) :: FRACREFA_10
  REAL(DP), DIMENSION(NG10_loc) :: FRACREFB_10

  REAL(DP):: ABSA_10(65,NG10_loc)
  REAL(DP):: ABSB_10(235,NG10_loc)
!!CDIR DUPLICATE(ABSA_10,256)
!!CDIR DUPLICATE(ABSB_10,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE *
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL    
  ! KA      : REAL     
  ! KB      : REAL     
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA11* - RRTM COEFFICIENTS FOR INTERVAL 11
  !     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG11_loc = 8
  REAL(DP), DIMENSION(NG11_loc) :: FRACREFA_11
  REAL(DP), DIMENSION(NG11_loc) :: FRACREFB_11

  REAL(DP):: ABSA_11(65,NG11_loc)
  REAL(DP):: ABSB_11(235,NG11_loc)
  REAL(DP):: SELFREF_11(10,NG11_loc)
!!CDIR DUPLICATE(ABSA_11,256)
!!CDIR DUPLICATE(ABSB_11,256)
!!CDIR DUPLICATE(SELFREF_11,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE *
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! FRACREFA: REAL    
  ! FRACREFB: REAL    
  ! KA      : REAL     
  ! KB      : REAL     
  ! SELFREF : REAL     
  !     -----------------------------------------------------------------
  
  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA12* - RRTM COEFFICIENTS FOR INTERVAL 12
  !     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG12_loc = 8

  REAL(DP):: FRACREFA_12(NG12_loc,9)
  REAL(DP):: ABSA_12(585,NG12_loc)
  REAL(DP):: SELFREF_12(10,NG12_loc)
!!CDIR DUPLICATE(ABSA_12,256)
!!CDIR DUPLICATE(SELFREF_12,256)
!!CDIR DUPLICATE(FRACREFA_12,256)

  REAL(DP):: STRRAT_12

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! ABSA    : REAL
  ! ABSB    : REAL
  ! FRACREFA: REAL    
  ! KA      : REAL     
  ! SELFREF : REAL
  ! STRRAT1 : REAL     
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA13* - RRTM COEFFICIENTS FOR INTERVAL 13
  !     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG13_loc = 4

  REAL(DP):: FRACREFA_13(NG13_loc,9)

  REAL(DP):: ABSA_13(585,NG13_loc)
  REAL(DP):: SELFREF_13(10,NG13_loc)
  REAL(DP):: STRRAT_13
!!CDIR DUPLICATE(ABSA_13,256)
!!CDIR DUPLICATE(SELFREF_13,256)
!!CDIR DUPLICATE(FRACREFA_13,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! FRACREFA: REAL    
  ! KA      : REAL     
  ! SELFREF : REAL
  ! STRRAT1 : REAL     
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA14* - RRTM COEFFICIENTS FOR INTERVAL 14
  !     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG14_loc = 2

  REAL(DP), DIMENSION(NG14_loc) :: FRACREFA_14
  REAL(DP), DIMENSION(NG14_loc) :: FRACREFB_14

  REAL(DP):: ABSA_14(65,NG14_loc)
  REAL(DP):: ABSB_14(235,NG14_loc)
  REAL(DP):: SELFREF_14(10,NG14_loc)
!!CDIR DUPLICATE(ABSA_14,256)
!!CDIR DUPLICATE(ABSB_14,256)
!!CDIR DUPLICATE(SELFREF_14,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE *
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! FRACREFA: REAL    
  ! FRACREFB: REAL    
  ! KA      : REAL     
  ! KB      : REAL     
  ! SELFREF : REAL     
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA15* - RRTM COEFFICIENTS FOR INTERVAL 15
  !     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG15_loc = 2

  REAL(DP):: FRACREFA_15(NG15_loc,9)
  
  REAL(DP):: ABSA_15(585,NG15_loc)
  REAL(DP):: SELFREF_15(10,NG15_loc)
  REAL(DP):: STRRAT_15
!!CDIR DUPLICATE(ABSA_15,256)
!!CDIR DUPLICATE(SELFREF_15,256)
!!CDIR DUPLICATE(FRACREFA_15,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE *
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! FRACREFA: REAL    
  ! KA      : REAL     
  ! SELFREF : REAL 
  ! STRRAT  : REAL    
  !     -----------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTA16* - RRTM COEFFICIENTS FOR INTERVAL 16
  !     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)
  !     -----------------------------------------------------------------

  INTEGER, PARAMETER :: NG16_loc = 2

  REAL(DP):: FRACREFA_16(NG16_loc,9)

  REAL(DP):: ABSA_16(585,NG16_loc)
  REAL(DP):: SELFREF_16(10,NG16_loc)
  REAL(DP):: STRRAT_16
!!CDIR DUPLICATE(ABSA_16,256)
!!CDIR DUPLICATE(SELFREF_16,256)
!!CDIR DUPLICATE(FRACREFA_16,256)
  
  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----   : ----   : ---------------------------------------------------
  ! FRACREFA: REAL    
  ! KA      : REAL     
  ! SELFREF : REAL     
  ! STRRAT  : REAL
  !     -----------------------------------------------------------------
  
  !     ------------------------------------------------------------------
  !     Parameters relevant to AER's RRTM-LW radiation scheme
  !     980714  JJMorcrette
  !     ------------------------------------------------------------------
  
  INTEGER, PARAMETER :: JPG    = 16
  INTEGER, PARAMETER :: JPGPT  = 140

  INTEGER, PARAMETER :: NG1  = 8
  INTEGER, PARAMETER :: NG2  = 14
  INTEGER, PARAMETER :: NG3  = 16
  INTEGER, PARAMETER :: NG4  = 14
  INTEGER, PARAMETER :: NG5  = 16
  INTEGER, PARAMETER :: NG6  = 8
  INTEGER, PARAMETER :: NG7  = 12
  INTEGER, PARAMETER :: NG8  = 8
  INTEGER, PARAMETER :: NG9  = 12
  INTEGER, PARAMETER :: NG10 = 6
  INTEGER, PARAMETER :: NG11 = 8
  INTEGER, PARAMETER :: NG12 = 8
  INTEGER, PARAMETER :: NG13 = 4
  INTEGER, PARAMETER :: NG14 = 2
  INTEGER, PARAMETER :: NG15 = 2
! INTEGER, PARAMETER :: NG16 = 8

  INTEGER, PARAMETER :: NGS1  = 8
  INTEGER, PARAMETER :: NGS2  = 22
  INTEGER, PARAMETER :: NGS3  = 38
  INTEGER, PARAMETER :: NGS4  = 52
  INTEGER, PARAMETER :: NGS5  = 68
  INTEGER, PARAMETER :: NGS6  = 76
  INTEGER, PARAMETER :: NGS7  = 88
  INTEGER, PARAMETER :: NGS8  = 96
  INTEGER, PARAMETER :: NGS9  = 108
  INTEGER, PARAMETER :: NGS10 = 114
  INTEGER, PARAMETER :: NGS11 = 122
  INTEGER, PARAMETER :: NGS12 = 130
  INTEGER, PARAMETER :: NGS13 = 134
  INTEGER, PARAMETER :: NGS14 = 136
  INTEGER, PARAMETER :: NGS15 = 138

  !    -------------------------------------------------------------------

  INTEGER :: NGC(JPBAND)
  INTEGER :: NGS(JPBAND)
  INTEGER :: NGN(JPGPT)
  INTEGER :: NGB(JPGPT)

  INTEGER :: NGM(JPG*JPBAND)
  REAL(DP):: WT(JPG)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----  :  ----   : ---------------------------------------------------
  !  NGC   : INTEGER :
  !  NGS   : INTEGER :
  !  NGN   : INTEGER :
  !  NGB   : INTEGER :
  !  NGM   : INTEGER :
  !  WT    : REAL    :
  !    -------------------------------------------------------------------

! fb_mk_20150928+
#if defined(MBM_RAD)
  INTEGER, PARAMETER :: NCORR = 2048 ! op_pj_20180725
  REAL(DP):: CORR1(0:2048)
  REAL(DP):: CORR2(0:2048)
#else
  INTEGER, PARAMETER :: NCORR = 200  ! op_pj_20180725
  REAL(DP):: CORR1(0:200)
  REAL(DP):: CORR2(0:200)
#endif
! fb_mk_20150928-
!!CDIR DUPLICATE(CORR1,256)
!!CDIR DUPLICATE(CORR2,256)

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
  !  NAME     TYPE     PURPOSE
  !  ----  :  ----   : ---------------------------------------------------
  ! CORR1  :  REAL   : 
  ! CORR2  :  REAL   :
  !    -------------------------------------------------------------------

  !     -----------------------------------------------------------------
  !*    ** *MO_RRTRF* - RRTM REFERENCE ATMOSPHERE
  !     -----------------------------------------------------------------

  REAL(DP), DIMENSION(59) :: PREF
  REAL(DP), DIMENSION(59) :: PREFLOG
  REAL(DP), DIMENSION(59) :: TREF

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15
  !  NAME     TYPE     PURPOSE
  !  ----  :  ----   : ---------------------------------------------------
  ! PREF   :  REAL    
  ! PREFLOG: REAL
  ! TREF   : REAL
  !     -----------------------------------------------------------------

  INTEGER , DIMENSION(16) :: NG
  INTEGER , DIMENSION(16) :: NSPA
  INTEGER , DIMENSION(16) :: NSPB

  REAL(DP), DIMENSION(16) :: WAVENUM1
  REAL(DP), DIMENSION(16) :: WAVENUM2
  REAL(DP), DIMENSION(16) :: DELWAVE

  REAL(DP), DIMENSION(181,16) :: TOTPLNK
  REAL(DP), DIMENSION(181)    :: TOTPLK16

  !     -----------------------------------------------------------------
  !        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM LW RADIATION **
  !     J.-J. MORCRETTE       E.C.M.W.F.      98/01/15
  !  NAME     TYPE     PURPOSE
  !  ----   : ----    : -------
  !  NG     : INTEGER : Number of k-coefficients in spectral intervals
  !  NSPA   : INTEGER :
  !  NSPB   : INTEGER :
  ! WAVENUM1: REAL    : Lower wavenumber spectral limit
  ! WAVENUM2: REAL    : Higher wavenumber spectral limit
  ! DELWAVE : REAL    : Spectral interval width
  ! TOTPLNK : REAL    :
  ! TOTPLK16: REAL    :
  !     -----------------------------------------------------------------

  PUBLIC :: &
         absa_1, absb_1, fracrefa_1, fracrefb_1, forref_1, selfref_1, &
         absa_2, absb_2, fracrefa_2, fracrefb_2, forref_2, selfref_2, &
         refparam_2, &
         absa_3, absb_3, fracrefa_3, fracrefb_3, forref_3, selfref_3, &
         absn2oa_3, &
         absn2ob_3, etaref_3, h2oref_3, n2oref_3, co2ref_3, strrat_3, &
         absa_4, absb_4, fracrefa_4, fracrefb_4, selfref_4, strrat1_4, &
         strrat2_4, &
         absa_5, absb_5, ccl4_5, fracrefa_5, fracrefb_5, selfref_5, &
         strrat1_5, strrat2_5, &
         absa_6, absco2_6, cfc11adj_6, cfc12_6, fracrefa_6, selfref_6, &
         absa_7, absb_7, absco2_7, fracrefa_7, fracrefb_7, selfref_7, &
         strrat_7, &
         absa_8, absb_8, fracrefa_8, fracrefb_8, selfref_8, absco2a_8, &
         absco2b_8, &
         absn2oa_8, absn2ob_8, cfc12_8, cfc22adj_8, h2oref_8, n2oref_8, &
         o3ref_8, &
         absa_9, absb_9, fracrefa_9, fracrefb_9, selfref_9, absn2o_9, &
         ch4ref_9, &
         etaref_9, h2oref_9, n2oref_9, strrat_9, &
         absa_10, absb_10, fracrefa_10, fracrefb_10, &
         absa_11, absb_11, fracrefa_11, selfref_11, fracrefb_11, &
         absa_12, fracrefa_12, selfref_12, strrat_12, &
         absa_13, fracrefa_13, selfref_13, strrat_13, absa_14, absb_14, &
         fracrefa_14, fracrefb_14, selfref_14, &
         absa_15, fracrefa_15, selfref_15, strrat_15, &
         absa_16, fracrefa_16, selfref_16, strrat_16,        & 
         NGC, NGS, NGM, NGN, NGB, WT,                        &
         corr1, corr2,                                       &
         PREF, PREFLOG, TREF,                                &
         NG, NSPA, NSPB, WAVENUM1, WAVENUM2, DELWAVE,        &
         TOTPLNK, TOTPLK16

  ! =========================================================================

  ! SUBROUTINES
  !
  PUBLIC  :: rad_lon_RRTM_RRTM_140GP
  !PRIVATE :: rad_lon_RRTM_SETCOEF_140GP,  &
  !           rad_lon_RRTM_GASABS1A_140GP, &
  !           rad_lon_RRTM_RTRN1A_140GP,   &
  !           rad_lon_RRTM_TAUMOL1, rad_lon_RRTM_TAUMOL2, &
  !           rad_lon_RRTM_TAUMOL3, rad_lon_RRTM_TAUMOL4, &
  !           rad_lon_RRTM_TAUMOL5, rad_lon_RRTM_TAUMOL6, &
  !           rad_lon_RRTM_TAUMOL7, rad_lon_RRTM_TAUMOL8, &
  !           rad_lon_RRTM_TAUMOL9, rad_lon_RRTM_TAUMOL10, &
  !           rad_lon_RRTM_TAUMOL11, rad_lon_RRTM_TAUMOL12, &
  !           rad_lon_RRTM_TAUMOL13, rad_lon_RRTM_TAUMOL14, &
  !           rad_lon_RRTM_TAUMOL15, rad_lon_RRTM_TAUMOL16
  !
  ! CALL-TREE of the subroutines
  !
  !          rad_lon_RRTM_RRTM_140GP    
  !                  |
  !                  --> rad_lon_RRTM_SETCOEF_140GP
  !                  |
  !                  --> rad_lon_RRTM_GASABS1A_140GP
  !                  |
  !                  --> rad_lon_RRTM_RTRN1A_140GP
  !                               |
  !                               --> rad_lon_RRTM_TAUMOL1
  !                               |
  !                               --> rad_lon_RRTM_TAUMOL2
  !                               :
  !                               :
  !                               |
  !                               --> rad_lon_RRTM_TAUMOL15
  !                               |
  !                               --> rad_lon_RRTM_TAUMOL16
  !
  PUBLIC  :: rad_lw_initialize
  
CONTAINS

  ! ===========================================================================
  SUBROUTINE rad_lw_initialize
    !SUBROUTINE surrtm

    !=========================================================================
    !
    !- Description:
    !
    !   Initialisation of RRTM modules
    !
    !   M.A. Giorgetta, MPI, March 2000
    !
    !=========================================================================

    IMPLICIT NONE

    !EXTERNAL rad_lon_surrtab, rad_lon_surrtpk, rad_lon_surrtrf, &
    !         rad_lon_surrtftr, &
    !         rad_lon_surrtbg2, rad_lon_surrta

    ! init RRTM modules
    ! -----------------

    CALL rad_lon_surrtpk        ! mo_rrtwn
    CALL rad_lon_surrtrf        ! mo_rrtrf
    CALL rad_lon_surrtftr       ! mo_rrtftr
    CALL rad_lon_surrtbg2       ! mo_rrtbg2
    CALL rad_lon_surrta         ! mo_rrtaN (N=1:16) from file 

    !END SUBROUTINE surrtm
  END SUBROUTINE rad_lw_initialize
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE rad_lon_SURRTPK

    !     Adapted from Eli J. Mlawer, Atmospheric & Environmental Research.
    !     by JJMorcrette, ECMWF
    !     ------------------------------------------------------------------

    !     mag, MPI, 25 February 2000: comment added
    
    !     ------------------------------------------------------------------


    IMPLICIT NONE
    NG( :) = (/16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/)
    NSPA( :) = (/1, 1,10, 9, 9, 1, 9, 1,11, 1, 1, 9, 9, 1, 9, 9/)
    NSPB( :) = (/1, 1, 5, 6, 5, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0/)
WAVENUM1( :) = (/&
    &10._DP, 250._DP, 500._DP, 630._DP, 700._DP, 820._DP, 980._DP,1080._DP &
 &,1180._DP,1390._DP,1480._DP,1800._DP,2080._DP,2250._DP,2380._DP,2600._DP/)
WAVENUM2( :) = (/&
   &250._DP, 500._DP, 630._DP, 700._DP, 820._DP, 980._DP,1080._DP,1180._DP &
 &,1390._DP,1480._DP,1800._DP,2080._DP,2250._DP,2380._DP,2600._DP,3000._DP/)
DELWAVE( :) = (/&
   &240._DP, 250._DP, 130._DP,  70._DP, 120._DP, 160._DP, 100._DP, 100._DP &
 &, 210._DP,  90._DP, 320._DP, 280._DP, 170._DP, 130._DP, 220._DP, 400._DP/)

TOTPLNK( :, 1) = (/&
&1.13735E-06_DP,1.15150E-06_DP,1.16569E-06_DP,1.17992E-06_DP,1.19419E-06_DP,&
&1.20850E-06_DP,1.22285E-06_DP,1.23723E-06_DP,1.25164E-06_DP,1.26610E-06_DP,&
&1.28059E-06_DP,1.29511E-06_DP,1.30967E-06_DP,1.32426E-06_DP,1.33889E-06_DP,&
&1.35355E-06_DP,1.36824E-06_DP,1.38296E-06_DP,1.39772E-06_DP,1.41250E-06_DP,&
&1.42732E-06_DP,1.44217E-06_DP,1.45704E-06_DP,1.47195E-06_DP,1.48689E-06_DP,&
&1.50185E-06_DP,1.51684E-06_DP,1.53186E-06_DP,1.54691E-06_DP,1.56198E-06_DP,&
&1.57709E-06_DP,1.59222E-06_DP,1.60737E-06_DP,1.62255E-06_DP,1.63776E-06_DP,&
&1.65299E-06_DP,1.66825E-06_DP,1.68352E-06_DP,1.69883E-06_DP,1.71416E-06_DP,&
&1.72951E-06_DP,1.74488E-06_DP,1.76028E-06_DP,1.77570E-06_DP,1.79114E-06_DP,&
&1.80661E-06_DP,1.82210E-06_DP,1.83760E-06_DP,1.85313E-06_DP,1.86868E-06_DP,&
&1.88425E-06_DP,1.89985E-06_DP,1.91546E-06_DP,1.93109E-06_DP,1.94674E-06_DP,&
&1.96241E-06_DP,1.97811E-06_DP,1.99381E-06_DP,2.00954E-06_DP,2.02529E-06_DP,&
&2.04105E-06_DP,2.05684E-06_DP,2.07264E-06_DP,2.08846E-06_DP,2.10429E-06_DP,&
&2.12015E-06_DP,2.13602E-06_DP,2.15190E-06_DP,2.16781E-06_DP,2.18373E-06_DP,&
&2.19966E-06_DP,2.21562E-06_DP,2.23159E-06_DP,2.24758E-06_DP,2.26358E-06_DP,&
&2.27959E-06_DP,2.29562E-06_DP,2.31167E-06_DP,2.32773E-06_DP,2.34381E-06_DP,&
&2.35990E-06_DP,2.37601E-06_DP,2.39212E-06_DP,2.40825E-06_DP,2.42440E-06_DP,&
&2.44056E-06_DP,2.45673E-06_DP,2.47292E-06_DP,2.48912E-06_DP,2.50533E-06_DP,&
&2.52157E-06_DP,2.53781E-06_DP,2.55406E-06_DP,2.57032E-06_DP,2.58660E-06_DP,&
&2.60289E-06_DP,2.61919E-06_DP,2.63550E-06_DP,2.65183E-06_DP,2.66817E-06_DP,&
&2.68452E-06_DP,2.70088E-06_DP,2.71726E-06_DP,2.73364E-06_DP,2.75003E-06_DP,&
&2.76644E-06_DP,2.78286E-06_DP,2.79929E-06_DP,2.81572E-06_DP,2.83218E-06_DP,&
&2.84864E-06_DP,2.86510E-06_DP,2.88159E-06_DP,2.89807E-06_DP,2.91458E-06_DP,&
&2.93109E-06_DP,2.94762E-06_DP,2.96415E-06_DP,2.98068E-06_DP,2.99724E-06_DP,&
&3.01379E-06_DP,3.03036E-06_DP,3.04693E-06_DP,3.06353E-06_DP,3.08013E-06_DP,&
&3.09674E-06_DP,3.11335E-06_DP,3.12998E-06_DP,3.14661E-06_DP,3.16324E-06_DP,&
&3.17989E-06_DP,3.19656E-06_DP,3.21323E-06_DP,3.22991E-06_DP,3.24658E-06_DP,&
&3.26328E-06_DP,3.27998E-06_DP,3.29669E-06_DP,3.31341E-06_DP,3.33013E-06_DP,&
&3.34686E-06_DP,3.36360E-06_DP,3.38034E-06_DP,3.39709E-06_DP,3.41387E-06_DP,&
&3.43063E-06_DP,3.44742E-06_DP,3.46420E-06_DP,3.48099E-06_DP,3.49779E-06_DP,&
&3.51461E-06_DP,3.53141E-06_DP,3.54824E-06_DP,3.56506E-06_DP,3.58191E-06_DP,&
&3.59875E-06_DP,3.61559E-06_DP,3.63244E-06_DP,3.64931E-06_DP,3.66617E-06_DP,&
&3.68305E-06_DP,3.69992E-06_DP,3.71682E-06_DP,3.73372E-06_DP,3.75061E-06_DP,&
&3.76753E-06_DP,3.78443E-06_DP,3.80136E-06_DP,3.81829E-06_DP,3.83522E-06_DP,&
&3.85215E-06_DP,3.86910E-06_DP,3.88605E-06_DP,3.90301E-06_DP,3.91997E-06_DP,&
&3.93694E-06_DP,3.95390E-06_DP,3.97087E-06_DP,3.98788E-06_DP,4.00485E-06_DP,&
&4.02187E-06_DP/)

TOTPLNK( :, 2) = (/&
&2.13441E-06_DP,2.18076E-06_DP,2.22758E-06_DP,2.27489E-06_DP,2.32268E-06_DP,&
&2.37093E-06_DP,2.41966E-06_DP,2.46886E-06_DP,2.51852E-06_DP,2.56864E-06_DP,&
&2.61922E-06_DP,2.67026E-06_DP,2.72175E-06_DP,2.77370E-06_DP,2.82609E-06_DP,&
&2.87893E-06_DP,2.93221E-06_DP,2.98593E-06_DP,3.04008E-06_DP,3.09468E-06_DP,&
&3.14970E-06_DP,3.20515E-06_DP,3.26103E-06_DP,3.31732E-06_DP,3.37404E-06_DP,&
&3.43118E-06_DP,3.48873E-06_DP,3.54669E-06_DP,3.60506E-06_DP,3.66383E-06_DP,&
&3.72301E-06_DP,3.78259E-06_DP,3.84256E-06_DP,3.90293E-06_DP,3.96368E-06_DP,&
&4.02483E-06_DP,4.08636E-06_DP,4.14828E-06_DP,4.21057E-06_DP,4.27324E-06_DP,&
&4.33629E-06_DP,4.39971E-06_DP,4.46350E-06_DP,4.52765E-06_DP,4.59217E-06_DP,&
&4.65705E-06_DP,4.72228E-06_DP,4.78787E-06_DP,4.85382E-06_DP,4.92011E-06_DP,&
&4.98675E-06_DP,5.05374E-06_DP,5.12106E-06_DP,5.18873E-06_DP,5.25674E-06_DP,&
&5.32507E-06_DP,5.39374E-06_DP,5.46274E-06_DP,5.53207E-06_DP,5.60172E-06_DP,&
&5.67169E-06_DP,5.74198E-06_DP,5.81259E-06_DP,5.88352E-06_DP,5.95475E-06_DP,&
&6.02629E-06_DP,6.09815E-06_DP,6.17030E-06_DP,6.24276E-06_DP,6.31552E-06_DP,&
&6.38858E-06_DP,6.46192E-06_DP,6.53557E-06_DP,6.60950E-06_DP,6.68373E-06_DP,&
&6.75824E-06_DP,6.83303E-06_DP,6.90810E-06_DP,6.98346E-06_DP,7.05909E-06_DP,&
&7.13500E-06_DP,7.21117E-06_DP,7.28763E-06_DP,7.36435E-06_DP,7.44134E-06_DP,&
&7.51859E-06_DP,7.59611E-06_DP,7.67388E-06_DP,7.75192E-06_DP,7.83021E-06_DP,&
&7.90875E-06_DP,7.98755E-06_DP,8.06660E-06_DP,8.14589E-06_DP,8.22544E-06_DP,&
&8.30522E-06_DP,8.38526E-06_DP,8.46553E-06_DP,8.54604E-06_DP,8.62679E-06_DP,&
&8.70777E-06_DP,8.78899E-06_DP,8.87043E-06_DP,8.95211E-06_DP,9.03402E-06_DP,&
&9.11616E-06_DP,9.19852E-06_DP,9.28109E-06_DP,9.36390E-06_DP,9.44692E-06_DP,&
&9.53015E-06_DP,9.61361E-06_DP,9.69729E-06_DP,9.78117E-06_DP,9.86526E-06_DP,&
&9.94957E-06_DP,1.00341E-05_DP,1.01188E-05_DP,1.02037E-05_DP,1.02888E-05_DP,&
&1.03742E-05_DP,1.04597E-05_DP,1.05454E-05_DP,1.06313E-05_DP,1.07175E-05_DP,&
&1.08038E-05_DP,1.08903E-05_DP,1.09770E-05_DP,1.10639E-05_DP,1.11509E-05_DP,&
&1.12382E-05_DP,1.13257E-05_DP,1.14133E-05_DP,1.15011E-05_DP,1.15891E-05_DP,&
&1.16773E-05_DP,1.17656E-05_DP,1.18542E-05_DP,1.19429E-05_DP,1.20317E-05_DP,&
&1.21208E-05_DP,1.22100E-05_DP,1.22994E-05_DP,1.23890E-05_DP,1.24787E-05_DP,&
&1.25686E-05_DP,1.26587E-05_DP,1.27489E-05_DP,1.28393E-05_DP,1.29299E-05_DP,&
&1.30206E-05_DP,1.31115E-05_DP,1.32025E-05_DP,1.32937E-05_DP,1.33850E-05_DP,&
&1.34765E-05_DP,1.35682E-05_DP,1.36600E-05_DP,1.37520E-05_DP,1.38441E-05_DP,&
&1.39364E-05_DP,1.40288E-05_DP,1.41213E-05_DP,1.42140E-05_DP,1.43069E-05_DP,&
&1.43999E-05_DP,1.44930E-05_DP,1.45863E-05_DP,1.46797E-05_DP,1.47733E-05_DP,&
&1.48670E-05_DP,1.49608E-05_DP,1.50548E-05_DP,1.51489E-05_DP,1.52431E-05_DP,&
&1.53375E-05_DP,1.54320E-05_DP,1.55267E-05_DP,1.56214E-05_DP,1.57164E-05_DP,&
&1.58114E-05_DP/)

TOTPLNK( :, 3) = (/&
&1.34822E-06_DP,1.39134E-06_DP,1.43530E-06_DP,1.48010E-06_DP,1.52574E-06_DP,&
&1.57222E-06_DP,1.61956E-06_DP,1.66774E-06_DP,1.71678E-06_DP,1.76666E-06_DP,&
&1.81741E-06_DP,1.86901E-06_DP,1.92147E-06_DP,1.97479E-06_DP,2.02898E-06_DP,&
&2.08402E-06_DP,2.13993E-06_DP,2.19671E-06_DP,2.25435E-06_DP,2.31285E-06_DP,&
&2.37222E-06_DP,2.43246E-06_DP,2.49356E-06_DP,2.55553E-06_DP,2.61837E-06_DP,&
&2.68207E-06_DP,2.74664E-06_DP,2.81207E-06_DP,2.87837E-06_DP,2.94554E-06_DP,&
&3.01356E-06_DP,3.08245E-06_DP,3.15221E-06_DP,3.22282E-06_DP,3.29429E-06_DP,&
&3.36662E-06_DP,3.43982E-06_DP,3.51386E-06_DP,3.58876E-06_DP,3.66451E-06_DP,&
&3.74112E-06_DP,3.81857E-06_DP,3.89688E-06_DP,3.97602E-06_DP,4.05601E-06_DP,&
&4.13685E-06_DP,4.21852E-06_DP,4.30104E-06_DP,4.38438E-06_DP,4.46857E-06_DP,&
&4.55358E-06_DP,4.63943E-06_DP,4.72610E-06_DP,4.81359E-06_DP,4.90191E-06_DP,&
&4.99105E-06_DP,5.08100E-06_DP,5.17176E-06_DP,5.26335E-06_DP,5.35573E-06_DP,&
&5.44892E-06_DP,5.54292E-06_DP,5.63772E-06_DP,5.73331E-06_DP,5.82970E-06_DP,&
&5.92688E-06_DP,6.02485E-06_DP,6.12360E-06_DP,6.22314E-06_DP,6.32346E-06_DP,&
&6.42455E-06_DP,6.52641E-06_DP,6.62906E-06_DP,6.73247E-06_DP,6.83664E-06_DP,&
&6.94156E-06_DP,7.04725E-06_DP,7.15370E-06_DP,7.26089E-06_DP,7.36883E-06_DP,&
&7.47752E-06_DP,7.58695E-06_DP,7.69712E-06_DP,7.80801E-06_DP,7.91965E-06_DP,&
&8.03201E-06_DP,8.14510E-06_DP,8.25891E-06_DP,8.37343E-06_DP,8.48867E-06_DP,&
&8.60463E-06_DP,8.72128E-06_DP,8.83865E-06_DP,8.95672E-06_DP,9.07548E-06_DP,&
&9.19495E-06_DP,9.31510E-06_DP,9.43594E-06_DP,9.55745E-06_DP,9.67966E-06_DP,&
&9.80254E-06_DP,9.92609E-06_DP,1.00503E-05_DP,1.01752E-05_DP,1.03008E-05_DP,&
&1.04270E-05_DP,1.05539E-05_DP,1.06814E-05_DP,1.08096E-05_DP,1.09384E-05_DP,&
&1.10679E-05_DP,1.11980E-05_DP,1.13288E-05_DP,1.14601E-05_DP,1.15922E-05_DP,&
&1.17248E-05_DP,1.18581E-05_DP,1.19920E-05_DP,1.21265E-05_DP,1.22616E-05_DP,&
&1.23973E-05_DP,1.25337E-05_DP,1.26706E-05_DP,1.28081E-05_DP,1.29463E-05_DP,&
&1.30850E-05_DP,1.32243E-05_DP,1.33642E-05_DP,1.35047E-05_DP,1.36458E-05_DP,&
&1.37875E-05_DP,1.39297E-05_DP,1.40725E-05_DP,1.42159E-05_DP,1.43598E-05_DP,&
&1.45044E-05_DP,1.46494E-05_DP,1.47950E-05_DP,1.49412E-05_DP,1.50879E-05_DP,&
&1.52352E-05_DP,1.53830E-05_DP,1.55314E-05_DP,1.56803E-05_DP,1.58297E-05_DP,&
&1.59797E-05_DP,1.61302E-05_DP,1.62812E-05_DP,1.64327E-05_DP,1.65848E-05_DP,&
&1.67374E-05_DP,1.68904E-05_DP,1.70441E-05_DP,1.71982E-05_DP,1.73528E-05_DP,&
&1.75079E-05_DP,1.76635E-05_DP,1.78197E-05_DP,1.79763E-05_DP,1.81334E-05_DP,&
&1.82910E-05_DP,1.84491E-05_DP,1.86076E-05_DP,1.87667E-05_DP,1.89262E-05_DP,&
&1.90862E-05_DP,1.92467E-05_DP,1.94076E-05_DP,1.95690E-05_DP,1.97309E-05_DP,&
&1.98932E-05_DP,2.00560E-05_DP,2.02193E-05_DP,2.03830E-05_DP,2.05472E-05_DP,&
&2.07118E-05_DP,2.08768E-05_DP,2.10423E-05_DP,2.12083E-05_DP,2.13747E-05_DP,&
&2.15414E-05_DP/)

TOTPLNK( :, 4) = (/&
&8.90528E-07_DP,9.24222E-07_DP,9.58757E-07_DP,9.94141E-07_DP,1.03038E-06_DP,&
&1.06748E-06_DP,1.10545E-06_DP,1.14430E-06_DP,1.18403E-06_DP,1.22465E-06_DP,&
&1.26618E-06_DP,1.30860E-06_DP,1.35193E-06_DP,1.39619E-06_DP,1.44136E-06_DP,&
&1.48746E-06_DP,1.53449E-06_DP,1.58246E-06_DP,1.63138E-06_DP,1.68124E-06_DP,&
&1.73206E-06_DP,1.78383E-06_DP,1.83657E-06_DP,1.89028E-06_DP,1.94495E-06_DP,&
&2.00060E-06_DP,2.05724E-06_DP,2.11485E-06_DP,2.17344E-06_DP,2.23303E-06_DP,&
&2.29361E-06_DP,2.35519E-06_DP,2.41777E-06_DP,2.48134E-06_DP,2.54592E-06_DP,&
&2.61151E-06_DP,2.67810E-06_DP,2.74571E-06_DP,2.81433E-06_DP,2.88396E-06_DP,&
&2.95461E-06_DP,3.02628E-06_DP,3.09896E-06_DP,3.17267E-06_DP,3.24741E-06_DP,&
&3.32316E-06_DP,3.39994E-06_DP,3.47774E-06_DP,3.55657E-06_DP,3.63642E-06_DP,&
&3.71731E-06_DP,3.79922E-06_DP,3.88216E-06_DP,3.96612E-06_DP,4.05112E-06_DP,&
&4.13714E-06_DP,4.22419E-06_DP,4.31227E-06_DP,4.40137E-06_DP,4.49151E-06_DP,&
&4.58266E-06_DP,4.67485E-06_DP,4.76806E-06_DP,4.86229E-06_DP,4.95754E-06_DP,&
&5.05383E-06_DP,5.15113E-06_DP,5.24946E-06_DP,5.34879E-06_DP,5.44916E-06_DP,&
&5.55053E-06_DP,5.65292E-06_DP,5.75632E-06_DP,5.86073E-06_DP,5.96616E-06_DP,&
&6.07260E-06_DP,6.18003E-06_DP,6.28848E-06_DP,6.39794E-06_DP,6.50838E-06_DP,&
&6.61983E-06_DP,6.73229E-06_DP,6.84573E-06_DP,6.96016E-06_DP,7.07559E-06_DP,&
&7.19200E-06_DP,7.30940E-06_DP,7.42779E-06_DP,7.54715E-06_DP,7.66749E-06_DP,&
&7.78882E-06_DP,7.91110E-06_DP,8.03436E-06_DP,8.15859E-06_DP,8.28379E-06_DP,&
&8.40994E-06_DP,8.53706E-06_DP,8.66515E-06_DP,8.79418E-06_DP,8.92416E-06_DP,&
&9.05510E-06_DP,9.18697E-06_DP,9.31979E-06_DP,9.45356E-06_DP,9.58826E-06_DP,&
&9.72389E-06_DP,9.86046E-06_DP,9.99793E-06_DP,1.01364E-05_DP,1.02757E-05_DP,&
&1.04159E-05_DP,1.05571E-05_DP,1.06992E-05_DP,1.08422E-05_DP,1.09861E-05_DP,&
&1.11309E-05_DP,1.12766E-05_DP,1.14232E-05_DP,1.15707E-05_DP,1.17190E-05_DP,&
&1.18683E-05_DP,1.20184E-05_DP,1.21695E-05_DP,1.23214E-05_DP,1.24741E-05_DP,&
&1.26277E-05_DP,1.27822E-05_DP,1.29376E-05_DP,1.30939E-05_DP,1.32509E-05_DP,&
&1.34088E-05_DP,1.35676E-05_DP,1.37273E-05_DP,1.38877E-05_DP,1.40490E-05_DP,&
&1.42112E-05_DP,1.43742E-05_DP,1.45380E-05_DP,1.47026E-05_DP,1.48680E-05_DP,&
&1.50343E-05_DP,1.52014E-05_DP,1.53692E-05_DP,1.55379E-05_DP,1.57074E-05_DP,&
&1.58778E-05_DP,1.60488E-05_DP,1.62207E-05_DP,1.63934E-05_DP,1.65669E-05_DP,&
&1.67411E-05_DP,1.69162E-05_DP,1.70920E-05_DP,1.72685E-05_DP,1.74459E-05_DP,&
&1.76240E-05_DP,1.78029E-05_DP,1.79825E-05_DP,1.81629E-05_DP,1.83440E-05_DP,&
&1.85259E-05_DP,1.87086E-05_DP,1.88919E-05_DP,1.90760E-05_DP,1.92609E-05_DP,&
&1.94465E-05_DP,1.96327E-05_DP,1.98199E-05_DP,2.00076E-05_DP,2.01961E-05_DP,&
&2.03853E-05_DP,2.05752E-05_DP,2.07658E-05_DP,2.09571E-05_DP,2.11491E-05_DP,&
&2.13418E-05_DP,2.15352E-05_DP,2.17294E-05_DP,2.19241E-05_DP,2.21196E-05_DP,&
&2.23158E-05_DP/)

TOTPLNK( :, 5) = (/&
&5.70230E-07_DP,5.94788E-07_DP,6.20085E-07_DP,6.46130E-07_DP,6.72936E-07_DP,&
&7.00512E-07_DP,7.28869E-07_DP,7.58019E-07_DP,7.87971E-07_DP,8.18734E-07_DP,&
&8.50320E-07_DP,8.82738E-07_DP,9.15999E-07_DP,9.50110E-07_DP,9.85084E-07_DP,&
&1.02093E-06_DP,1.05765E-06_DP,1.09527E-06_DP,1.13378E-06_DP,1.17320E-06_DP,&
&1.21353E-06_DP,1.25479E-06_DP,1.29698E-06_DP,1.34011E-06_DP,1.38419E-06_DP,&
&1.42923E-06_DP,1.47523E-06_DP,1.52221E-06_DP,1.57016E-06_DP,1.61910E-06_DP,&
&1.66904E-06_DP,1.71997E-06_DP,1.77192E-06_DP,1.82488E-06_DP,1.87886E-06_DP,&
&1.93387E-06_DP,1.98991E-06_DP,2.04699E-06_DP,2.10512E-06_DP,2.16430E-06_DP,&
&2.22454E-06_DP,2.28584E-06_DP,2.34821E-06_DP,2.41166E-06_DP,2.47618E-06_DP,&
&2.54178E-06_DP,2.60847E-06_DP,2.67626E-06_DP,2.74514E-06_DP,2.81512E-06_DP,&
&2.88621E-06_DP,2.95841E-06_DP,3.03172E-06_DP,3.10615E-06_DP,3.18170E-06_DP,&
&3.25838E-06_DP,3.33618E-06_DP,3.41511E-06_DP,3.49518E-06_DP,3.57639E-06_DP,&
&3.65873E-06_DP,3.74221E-06_DP,3.82684E-06_DP,3.91262E-06_DP,3.99955E-06_DP,&
&4.08763E-06_DP,4.17686E-06_DP,4.26725E-06_DP,4.35880E-06_DP,4.45150E-06_DP,&
&4.54537E-06_DP,4.64039E-06_DP,4.73659E-06_DP,4.83394E-06_DP,4.93246E-06_DP,&
&5.03215E-06_DP,5.13301E-06_DP,5.23504E-06_DP,5.33823E-06_DP,5.44260E-06_DP,&
&5.54814E-06_DP,5.65484E-06_DP,5.76272E-06_DP,5.87177E-06_DP,5.98199E-06_DP,&
&6.09339E-06_DP,6.20596E-06_DP,6.31969E-06_DP,6.43460E-06_DP,6.55068E-06_DP,&
&6.66793E-06_DP,6.78636E-06_DP,6.90595E-06_DP,7.02670E-06_DP,7.14863E-06_DP,&
&7.27173E-06_DP,7.39599E-06_DP,7.52142E-06_DP,7.64802E-06_DP,7.77577E-06_DP,&
&7.90469E-06_DP,8.03477E-06_DP,8.16601E-06_DP,8.29841E-06_DP,8.43198E-06_DP,&
&8.56669E-06_DP,8.70256E-06_DP,8.83957E-06_DP,8.97775E-06_DP,9.11706E-06_DP,&
&9.25753E-06_DP,9.39915E-06_DP,9.54190E-06_DP,9.68580E-06_DP,9.83085E-06_DP,&
&9.97704E-06_DP,1.01243E-05_DP,1.02728E-05_DP,1.04224E-05_DP,1.05731E-05_DP,&
&1.07249E-05_DP,1.08779E-05_DP,1.10320E-05_DP,1.11872E-05_DP,1.13435E-05_DP,&
&1.15009E-05_DP,1.16595E-05_DP,1.18191E-05_DP,1.19799E-05_DP,1.21418E-05_DP,&
&1.23048E-05_DP,1.24688E-05_DP,1.26340E-05_DP,1.28003E-05_DP,1.29676E-05_DP,&
&1.31361E-05_DP,1.33056E-05_DP,1.34762E-05_DP,1.36479E-05_DP,1.38207E-05_DP,&
&1.39945E-05_DP,1.41694E-05_DP,1.43454E-05_DP,1.45225E-05_DP,1.47006E-05_DP,&
&1.48797E-05_DP,1.50600E-05_DP,1.52413E-05_DP,1.54236E-05_DP,1.56070E-05_DP,&
&1.57914E-05_DP,1.59768E-05_DP,1.61633E-05_DP,1.63509E-05_DP,1.65394E-05_DP,&
&1.67290E-05_DP,1.69197E-05_DP,1.71113E-05_DP,1.73040E-05_DP,1.74976E-05_DP,&
&1.76923E-05_DP,1.78880E-05_DP,1.80847E-05_DP,1.82824E-05_DP,1.84811E-05_DP,&
&1.86808E-05_DP,1.88814E-05_DP,1.90831E-05_DP,1.92857E-05_DP,1.94894E-05_DP,&
&1.96940E-05_DP,1.98996E-05_DP,2.01061E-05_DP,2.03136E-05_DP,2.05221E-05_DP,&
&2.07316E-05_DP,2.09420E-05_DP,2.11533E-05_DP,2.13657E-05_DP,2.15789E-05_DP,&
&2.17931E-05_DP/)

TOTPLNK( :, 6) = (/&
&2.73493E-07_DP,2.87408E-07_DP,3.01848E-07_DP,3.16825E-07_DP,3.32352E-07_DP,&
&3.48439E-07_DP,3.65100E-07_DP,3.82346E-07_DP,4.00189E-07_DP,4.18641E-07_DP,&
&4.37715E-07_DP,4.57422E-07_DP,4.77774E-07_DP,4.98784E-07_DP,5.20464E-07_DP,&
&5.42824E-07_DP,5.65879E-07_DP,5.89638E-07_DP,6.14115E-07_DP,6.39320E-07_DP,&
&6.65266E-07_DP,6.91965E-07_DP,7.19427E-07_DP,7.47666E-07_DP,7.76691E-07_DP,&
&8.06516E-07_DP,8.37151E-07_DP,8.68607E-07_DP,9.00896E-07_DP,9.34029E-07_DP,&
&9.68018E-07_DP,1.00287E-06_DP,1.03860E-06_DP,1.07522E-06_DP,1.11274E-06_DP,&
&1.15117E-06_DP,1.19052E-06_DP,1.23079E-06_DP,1.27201E-06_DP,1.31418E-06_DP,&
&1.35731E-06_DP,1.40141E-06_DP,1.44650E-06_DP,1.49257E-06_DP,1.53965E-06_DP,&
&1.58773E-06_DP,1.63684E-06_DP,1.68697E-06_DP,1.73815E-06_DP,1.79037E-06_DP,&
&1.84365E-06_DP,1.89799E-06_DP,1.95341E-06_DP,2.00991E-06_DP,2.06750E-06_DP,&
&2.12619E-06_DP,2.18599E-06_DP,2.24691E-06_DP,2.30895E-06_DP,2.37212E-06_DP,&
&2.43643E-06_DP,2.50189E-06_DP,2.56851E-06_DP,2.63628E-06_DP,2.70523E-06_DP,&
&2.77536E-06_DP,2.84666E-06_DP,2.91916E-06_DP,2.99286E-06_DP,3.06776E-06_DP,&
&3.14387E-06_DP,3.22120E-06_DP,3.29975E-06_DP,3.37953E-06_DP,3.46054E-06_DP,&
&3.54280E-06_DP,3.62630E-06_DP,3.71105E-06_DP,3.79707E-06_DP,3.88434E-06_DP,&
&3.97288E-06_DP,4.06270E-06_DP,4.15380E-06_DP,4.24617E-06_DP,4.33984E-06_DP,&
&4.43479E-06_DP,4.53104E-06_DP,4.62860E-06_DP,4.72746E-06_DP,4.82763E-06_DP,&
&4.92911E-06_DP,5.03191E-06_DP,5.13603E-06_DP,5.24147E-06_DP,5.34824E-06_DP,&
&5.45634E-06_DP,5.56578E-06_DP,5.67656E-06_DP,5.78867E-06_DP,5.90213E-06_DP,&
&6.01694E-06_DP,6.13309E-06_DP,6.25060E-06_DP,6.36947E-06_DP,6.48968E-06_DP,&
&6.61126E-06_DP,6.73420E-06_DP,6.85850E-06_DP,6.98417E-06_DP,7.11120E-06_DP,&
&7.23961E-06_DP,7.36938E-06_DP,7.50053E-06_DP,7.63305E-06_DP,7.76694E-06_DP,&
&7.90221E-06_DP,8.03887E-06_DP,8.17690E-06_DP,8.31632E-06_DP,8.45710E-06_DP,&
&8.59928E-06_DP,8.74282E-06_DP,8.88776E-06_DP,9.03409E-06_DP,9.18179E-06_DP,&
&9.33088E-06_DP,9.48136E-06_DP,9.63323E-06_DP,9.78648E-06_DP,9.94111E-06_DP,&
&1.00971E-05_DP,1.02545E-05_DP,1.04133E-05_DP,1.05735E-05_DP,1.07351E-05_DP,&
&1.08980E-05_DP,1.10624E-05_DP,1.12281E-05_DP,1.13952E-05_DP,1.15637E-05_DP,&
&1.17335E-05_DP,1.19048E-05_DP,1.20774E-05_DP,1.22514E-05_DP,1.24268E-05_DP,&
&1.26036E-05_DP,1.27817E-05_DP,1.29612E-05_DP,1.31421E-05_DP,1.33244E-05_DP,&
&1.35080E-05_DP,1.36930E-05_DP,1.38794E-05_DP,1.40672E-05_DP,1.42563E-05_DP,&
&1.44468E-05_DP,1.46386E-05_DP,1.48318E-05_DP,1.50264E-05_DP,1.52223E-05_DP,&
&1.54196E-05_DP,1.56182E-05_DP,1.58182E-05_DP,1.60196E-05_DP,1.62223E-05_DP,&
&1.64263E-05_DP,1.66317E-05_DP,1.68384E-05_DP,1.70465E-05_DP,1.72559E-05_DP,&
&1.74666E-05_DP,1.76787E-05_DP,1.78921E-05_DP,1.81069E-05_DP,1.83230E-05_DP,&
&1.85404E-05_DP,1.87591E-05_DP,1.89791E-05_DP,1.92005E-05_DP,1.94232E-05_DP,&
&1.96471E-05_DP/)

TOTPLNK( :, 7) = (/&
&1.25349E-07_DP,1.32735E-07_DP,1.40458E-07_DP,1.48527E-07_DP,1.56954E-07_DP,&
&1.65748E-07_DP,1.74920E-07_DP,1.84481E-07_DP,1.94443E-07_DP,2.04814E-07_DP,&
&2.15608E-07_DP,2.26835E-07_DP,2.38507E-07_DP,2.50634E-07_DP,2.63229E-07_DP,&
&2.76301E-07_DP,2.89864E-07_DP,3.03930E-07_DP,3.18508E-07_DP,3.33612E-07_DP,&
&3.49253E-07_DP,3.65443E-07_DP,3.82195E-07_DP,3.99519E-07_DP,4.17428E-07_DP,&
&4.35934E-07_DP,4.55050E-07_DP,4.74785E-07_DP,4.95155E-07_DP,5.16170E-07_DP,&
&5.37844E-07_DP,5.60186E-07_DP,5.83211E-07_DP,6.06929E-07_DP,6.31355E-07_DP,&
&6.56498E-07_DP,6.82373E-07_DP,7.08990E-07_DP,7.36362E-07_DP,7.64501E-07_DP,&
&7.93420E-07_DP,8.23130E-07_DP,8.53643E-07_DP,8.84971E-07_DP,9.17128E-07_DP,&
&9.50123E-07_DP,9.83969E-07_DP,1.01868E-06_DP,1.05426E-06_DP,1.09073E-06_DP,&
&1.12810E-06_DP,1.16638E-06_DP,1.20558E-06_DP,1.24572E-06_DP,1.28680E-06_DP,&
&1.32883E-06_DP,1.37183E-06_DP,1.41581E-06_DP,1.46078E-06_DP,1.50675E-06_DP,&
&1.55374E-06_DP,1.60174E-06_DP,1.65078E-06_DP,1.70087E-06_DP,1.75200E-06_DP,&
&1.80421E-06_DP,1.85749E-06_DP,1.91186E-06_DP,1.96732E-06_DP,2.02389E-06_DP,&
&2.08159E-06_DP,2.14040E-06_DP,2.20035E-06_DP,2.26146E-06_DP,2.32372E-06_DP,&
&2.38714E-06_DP,2.45174E-06_DP,2.51753E-06_DP,2.58451E-06_DP,2.65270E-06_DP,&
&2.72210E-06_DP,2.79272E-06_DP,2.86457E-06_DP,2.93767E-06_DP,3.01201E-06_DP,&
&3.08761E-06_DP,3.16448E-06_DP,3.24261E-06_DP,3.32204E-06_DP,3.40275E-06_DP,&
&3.48476E-06_DP,3.56808E-06_DP,3.65271E-06_DP,3.73866E-06_DP,3.82595E-06_DP,&
&3.91456E-06_DP,4.00453E-06_DP,4.09584E-06_DP,4.18851E-06_DP,4.28254E-06_DP,&
&4.37796E-06_DP,4.47475E-06_DP,4.57293E-06_DP,4.67249E-06_DP,4.77346E-06_DP,&
&4.87583E-06_DP,4.97961E-06_DP,5.08481E-06_DP,5.19143E-06_DP,5.29948E-06_DP,&
&5.40896E-06_DP,5.51989E-06_DP,5.63226E-06_DP,5.74608E-06_DP,5.86136E-06_DP,&
&5.97810E-06_DP,6.09631E-06_DP,6.21597E-06_DP,6.33713E-06_DP,6.45976E-06_DP,&
&6.58388E-06_DP,6.70950E-06_DP,6.83661E-06_DP,6.96521E-06_DP,7.09531E-06_DP,&
&7.22692E-06_DP,7.36005E-06_DP,7.49468E-06_DP,7.63084E-06_DP,7.76851E-06_DP,&
&7.90773E-06_DP,8.04846E-06_DP,8.19072E-06_DP,8.33452E-06_DP,8.47985E-06_DP,&
&8.62674E-06_DP,8.77517E-06_DP,8.92514E-06_DP,9.07666E-06_DP,9.22975E-06_DP,&
&9.38437E-06_DP,9.54057E-06_DP,9.69832E-06_DP,9.85762E-06_DP,1.00185E-05_DP,&
&1.01810E-05_DP,1.03450E-05_DP,1.05106E-05_DP,1.06777E-05_DP,1.08465E-05_DP,&
&1.10168E-05_DP,1.11887E-05_DP,1.13621E-05_DP,1.15372E-05_DP,1.17138E-05_DP,&
&1.18920E-05_DP,1.20718E-05_DP,1.22532E-05_DP,1.24362E-05_DP,1.26207E-05_DP,&
&1.28069E-05_DP,1.29946E-05_DP,1.31839E-05_DP,1.33749E-05_DP,1.35674E-05_DP,&
&1.37615E-05_DP,1.39572E-05_DP,1.41544E-05_DP,1.43533E-05_DP,1.45538E-05_DP,&
&1.47558E-05_DP,1.49595E-05_DP,1.51647E-05_DP,1.53716E-05_DP,1.55800E-05_DP,&
&1.57900E-05_DP,1.60017E-05_DP,1.62149E-05_DP,1.64296E-05_DP,1.66460E-05_DP,&
&1.68640E-05_DP/)

TOTPLNK( :, 8) = (/&
&6.74445E-08_DP,7.18176E-08_DP,7.64153E-08_DP,8.12456E-08_DP,8.63170E-08_DP,&
&9.16378E-08_DP,9.72168E-08_DP,1.03063E-07_DP,1.09184E-07_DP,1.15591E-07_DP,&
&1.22292E-07_DP,1.29296E-07_DP,1.36613E-07_DP,1.44253E-07_DP,1.52226E-07_DP,&
&1.60540E-07_DP,1.69207E-07_DP,1.78236E-07_DP,1.87637E-07_DP,1.97421E-07_DP,&
&2.07599E-07_DP,2.18181E-07_DP,2.29177E-07_DP,2.40598E-07_DP,2.52456E-07_DP,&
&2.64761E-07_DP,2.77523E-07_DP,2.90755E-07_DP,3.04468E-07_DP,3.18673E-07_DP,&
&3.33381E-07_DP,3.48603E-07_DP,3.64352E-07_DP,3.80638E-07_DP,3.97474E-07_DP,&
&4.14871E-07_DP,4.32841E-07_DP,4.51395E-07_DP,4.70547E-07_DP,4.90306E-07_DP,&
&5.10687E-07_DP,5.31699E-07_DP,5.53357E-07_DP,5.75670E-07_DP,5.98652E-07_DP,&
&6.22315E-07_DP,6.46672E-07_DP,6.71731E-07_DP,6.97511E-07_DP,7.24018E-07_DP,&
&7.51266E-07_DP,7.79269E-07_DP,8.08038E-07_DP,8.37584E-07_DP,8.67922E-07_DP,&
&8.99061E-07_DP,9.31016E-07_DP,9.63797E-07_DP,9.97417E-07_DP,1.03189E-06_DP,&
&1.06722E-06_DP,1.10343E-06_DP,1.14053E-06_DP,1.17853E-06_DP,1.21743E-06_DP,&
&1.25726E-06_DP,1.29803E-06_DP,1.33974E-06_DP,1.38241E-06_DP,1.42606E-06_DP,&
&1.47068E-06_DP,1.51630E-06_DP,1.56293E-06_DP,1.61056E-06_DP,1.65924E-06_DP,&
&1.70894E-06_DP,1.75971E-06_DP,1.81153E-06_DP,1.86443E-06_DP,1.91841E-06_DP,&
&1.97350E-06_DP,2.02968E-06_DP,2.08699E-06_DP,2.14543E-06_DP,2.20500E-06_DP,&
&2.26573E-06_DP,2.32762E-06_DP,2.39068E-06_DP,2.45492E-06_DP,2.52036E-06_DP,&
&2.58700E-06_DP,2.65485E-06_DP,2.72393E-06_DP,2.79424E-06_DP,2.86580E-06_DP,&
&2.93861E-06_DP,3.01269E-06_DP,3.08803E-06_DP,3.16467E-06_DP,3.24259E-06_DP,&
&3.32181E-06_DP,3.40235E-06_DP,3.48420E-06_DP,3.56739E-06_DP,3.65192E-06_DP,&
&3.73779E-06_DP,3.82502E-06_DP,3.91362E-06_DP,4.00359E-06_DP,4.09494E-06_DP,&
&4.18768E-06_DP,4.28182E-06_DP,4.37737E-06_DP,4.47434E-06_DP,4.57273E-06_DP,&
&4.67254E-06_DP,4.77380E-06_DP,4.87651E-06_DP,4.98067E-06_DP,5.08630E-06_DP,&
&5.19339E-06_DP,5.30196E-06_DP,5.41201E-06_DP,5.52356E-06_DP,5.63660E-06_DP,&
&5.75116E-06_DP,5.86722E-06_DP,5.98479E-06_DP,6.10390E-06_DP,6.22453E-06_DP,&
&6.34669E-06_DP,6.47042E-06_DP,6.59569E-06_DP,6.72252E-06_DP,6.85090E-06_DP,&
&6.98085E-06_DP,7.11238E-06_DP,7.24549E-06_DP,7.38019E-06_DP,7.51646E-06_DP,&
&7.65434E-06_DP,7.79382E-06_DP,7.93490E-06_DP,8.07760E-06_DP,8.22192E-06_DP,&
&8.36784E-06_DP,8.51540E-06_DP,8.66459E-06_DP,8.81542E-06_DP,8.96786E-06_DP,&
&9.12197E-06_DP,9.27772E-06_DP,9.43513E-06_DP,9.59419E-06_DP,9.75490E-06_DP,&
&9.91728E-06_DP,1.00813E-05_DP,1.02471E-05_DP,1.04144E-05_DP,1.05835E-05_DP,&
&1.07543E-05_DP,1.09267E-05_DP,1.11008E-05_DP,1.12766E-05_DP,1.14541E-05_DP,&
&1.16333E-05_DP,1.18142E-05_DP,1.19969E-05_DP,1.21812E-05_DP,1.23672E-05_DP,&
&1.25549E-05_DP,1.27443E-05_DP,1.29355E-05_DP,1.31284E-05_DP,1.33229E-05_DP,&
&1.35193E-05_DP,1.37173E-05_DP,1.39170E-05_DP,1.41185E-05_DP,1.43217E-05_DP,&
&1.45267E-05_DP/)

TOTPLNK( :, 9) = (/&
&2.61522E-08_DP,2.80613E-08_DP,3.00838E-08_DP,3.22250E-08_DP,3.44899E-08_DP,&
&3.68841E-08_DP,3.94129E-08_DP,4.20820E-08_DP,4.48973E-08_DP,4.78646E-08_DP,&
&5.09901E-08_DP,5.42799E-08_DP,5.77405E-08_DP,6.13784E-08_DP,6.52001E-08_DP,&
&6.92126E-08_DP,7.34227E-08_DP,7.78375E-08_DP,8.24643E-08_DP,8.73103E-08_DP,&
&9.23832E-08_DP,9.76905E-08_DP,1.03240E-07_DP,1.09039E-07_DP,1.15097E-07_DP,&
&1.21421E-07_DP,1.28020E-07_DP,1.34902E-07_DP,1.42075E-07_DP,1.49548E-07_DP,&
&1.57331E-07_DP,1.65432E-07_DP,1.73860E-07_DP,1.82624E-07_DP,1.91734E-07_DP,&
&2.01198E-07_DP,2.11028E-07_DP,2.21231E-07_DP,2.31818E-07_DP,2.42799E-07_DP,&
&2.54184E-07_DP,2.65983E-07_DP,2.78205E-07_DP,2.90862E-07_DP,3.03963E-07_DP,&
&3.17519E-07_DP,3.31541E-07_DP,3.46039E-07_DP,3.61024E-07_DP,3.76507E-07_DP,&
&3.92498E-07_DP,4.09008E-07_DP,4.26050E-07_DP,4.43633E-07_DP,4.61769E-07_DP,&
&4.80469E-07_DP,4.99744E-07_DP,5.19606E-07_DP,5.40067E-07_DP,5.61136E-07_DP,&
&5.82828E-07_DP,6.05152E-07_DP,6.28120E-07_DP,6.51745E-07_DP,6.76038E-07_DP,&
&7.01010E-07_DP,7.26674E-07_DP,7.53041E-07_DP,7.80124E-07_DP,8.07933E-07_DP,&
&8.36482E-07_DP,8.65781E-07_DP,8.95845E-07_DP,9.26683E-07_DP,9.58308E-07_DP,&
&9.90732E-07_DP,1.02397E-06_DP,1.05803E-06_DP,1.09292E-06_DP,1.12866E-06_DP,&
&1.16526E-06_DP,1.20274E-06_DP,1.24109E-06_DP,1.28034E-06_DP,1.32050E-06_DP,&
&1.36158E-06_DP,1.40359E-06_DP,1.44655E-06_DP,1.49046E-06_DP,1.53534E-06_DP,&
&1.58120E-06_DP,1.62805E-06_DP,1.67591E-06_DP,1.72478E-06_DP,1.77468E-06_DP,&
&1.82561E-06_DP,1.87760E-06_DP,1.93066E-06_DP,1.98479E-06_DP,2.04000E-06_DP,&
&2.09631E-06_DP,2.15373E-06_DP,2.21228E-06_DP,2.27196E-06_DP,2.33278E-06_DP,&
&2.39475E-06_DP,2.45790E-06_DP,2.52222E-06_DP,2.58773E-06_DP,2.65445E-06_DP,&
&2.72238E-06_DP,2.79152E-06_DP,2.86191E-06_DP,2.93354E-06_DP,3.00643E-06_DP,&
&3.08058E-06_DP,3.15601E-06_DP,3.23273E-06_DP,3.31075E-06_DP,3.39009E-06_DP,&
&3.47074E-06_DP,3.55272E-06_DP,3.63605E-06_DP,3.72072E-06_DP,3.80676E-06_DP,&
&3.89417E-06_DP,3.98297E-06_DP,4.07315E-06_DP,4.16474E-06_DP,4.25774E-06_DP,&
&4.35217E-06_DP,4.44802E-06_DP,4.54532E-06_DP,4.64406E-06_DP,4.74428E-06_DP,&
&4.84595E-06_DP,4.94911E-06_DP,5.05376E-06_DP,5.15990E-06_DP,5.26755E-06_DP,&
&5.37671E-06_DP,5.48741E-06_DP,5.59963E-06_DP,5.71340E-06_DP,5.82871E-06_DP,&
&5.94559E-06_DP,6.06403E-06_DP,6.18404E-06_DP,6.30565E-06_DP,6.42885E-06_DP,&
&6.55364E-06_DP,6.68004E-06_DP,6.80806E-06_DP,6.93771E-06_DP,7.06898E-06_DP,&
&7.20190E-06_DP,7.33646E-06_DP,7.47267E-06_DP,7.61056E-06_DP,7.75010E-06_DP,&
&7.89133E-06_DP,8.03423E-06_DP,8.17884E-06_DP,8.32514E-06_DP,8.47314E-06_DP,&
&8.62284E-06_DP,8.77427E-06_DP,8.92743E-06_DP,9.08231E-06_DP,9.23893E-06_DP,&
&9.39729E-06_DP,9.55741E-06_DP,9.71927E-06_DP,9.88291E-06_DP,1.00483E-05_DP,&
&1.02155E-05_DP,1.03844E-05_DP,1.05552E-05_DP,1.07277E-05_DP,1.09020E-05_DP,&
&1.10781E-05_DP/)

TOTPLNK( :,10) = (/&
&8.89300E-09_DP,9.63263E-09_DP,1.04235E-08_DP,1.12685E-08_DP,1.21703E-08_DP,&
&1.31321E-08_DP,1.41570E-08_DP,1.52482E-08_DP,1.64090E-08_DP,1.76428E-08_DP,&
&1.89533E-08_DP,2.03441E-08_DP,2.18190E-08_DP,2.33820E-08_DP,2.50370E-08_DP,&
&2.67884E-08_DP,2.86402E-08_DP,3.05969E-08_DP,3.26632E-08_DP,3.48436E-08_DP,&
&3.71429E-08_DP,3.95660E-08_DP,4.21179E-08_DP,4.48040E-08_DP,4.76294E-08_DP,&
&5.05996E-08_DP,5.37201E-08_DP,5.69966E-08_DP,6.04349E-08_DP,6.40411E-08_DP,&
&6.78211E-08_DP,7.17812E-08_DP,7.59276E-08_DP,8.02670E-08_DP,8.48059E-08_DP,&
&8.95508E-08_DP,9.45090E-08_DP,9.96873E-08_DP,1.05093E-07_DP,1.10733E-07_DP,&
&1.16614E-07_DP,1.22745E-07_DP,1.29133E-07_DP,1.35786E-07_DP,1.42711E-07_DP,&
&1.49916E-07_DP,1.57410E-07_DP,1.65202E-07_DP,1.73298E-07_DP,1.81709E-07_DP,&
&1.90441E-07_DP,1.99505E-07_DP,2.08908E-07_DP,2.18660E-07_DP,2.28770E-07_DP,&
&2.39247E-07_DP,2.50101E-07_DP,2.61340E-07_DP,2.72974E-07_DP,2.85013E-07_DP,&
&2.97467E-07_DP,3.10345E-07_DP,3.23657E-07_DP,3.37413E-07_DP,3.51623E-07_DP,&
&3.66298E-07_DP,3.81448E-07_DP,3.97082E-07_DP,4.13212E-07_DP,4.29848E-07_DP,&
&4.47000E-07_DP,4.64680E-07_DP,4.82898E-07_DP,5.01664E-07_DP,5.20991E-07_DP,&
&5.40888E-07_DP,5.61369E-07_DP,5.82440E-07_DP,6.04118E-07_DP,6.26410E-07_DP,&
&6.49329E-07_DP,6.72887E-07_DP,6.97095E-07_DP,7.21964E-07_DP,7.47506E-07_DP,&
&7.73732E-07_DP,8.00655E-07_DP,8.28287E-07_DP,8.56635E-07_DP,8.85717E-07_DP,&
&9.15542E-07_DP,9.46122E-07_DP,9.77469E-07_DP,1.00960E-06_DP,1.04251E-06_DP,&
&1.07623E-06_DP,1.11077E-06_DP,1.14613E-06_DP,1.18233E-06_DP,1.21939E-06_DP,&
&1.25730E-06_DP,1.29610E-06_DP,1.33578E-06_DP,1.37636E-06_DP,1.41785E-06_DP,&
&1.46027E-06_DP,1.50362E-06_DP,1.54792E-06_DP,1.59319E-06_DP,1.63942E-06_DP,&
&1.68665E-06_DP,1.73487E-06_DP,1.78410E-06_DP,1.83435E-06_DP,1.88564E-06_DP,&
&1.93797E-06_DP,1.99136E-06_DP,2.04582E-06_DP,2.10137E-06_DP,2.15801E-06_DP,&
&2.21576E-06_DP,2.27463E-06_DP,2.33462E-06_DP,2.39577E-06_DP,2.45806E-06_DP,&
&2.52153E-06_DP,2.58617E-06_DP,2.65201E-06_DP,2.71905E-06_DP,2.78730E-06_DP,&
&2.85678E-06_DP,2.92749E-06_DP,2.99946E-06_DP,3.07269E-06_DP,3.14720E-06_DP,&
&3.22299E-06_DP,3.30007E-06_DP,3.37847E-06_DP,3.45818E-06_DP,3.53923E-06_DP,&
&3.62161E-06_DP,3.70535E-06_DP,3.79046E-06_DP,3.87695E-06_DP,3.96481E-06_DP,&
&4.05409E-06_DP,4.14477E-06_DP,4.23687E-06_DP,4.33040E-06_DP,4.42538E-06_DP,&
&4.52180E-06_DP,4.61969E-06_DP,4.71905E-06_DP,4.81991E-06_DP,4.92226E-06_DP,&
&5.02611E-06_DP,5.13148E-06_DP,5.23839E-06_DP,5.34681E-06_DP,5.45681E-06_DP,&
&5.56835E-06_DP,5.68146E-06_DP,5.79614E-06_DP,5.91242E-06_DP,6.03030E-06_DP,&
&6.14978E-06_DP,6.27088E-06_DP,6.39360E-06_DP,6.51798E-06_DP,6.64398E-06_DP,&
&6.77165E-06_DP,6.90099E-06_DP,7.03198E-06_DP,7.16468E-06_DP,7.29906E-06_DP,&
&7.43514E-06_DP,7.57294E-06_DP,7.71244E-06_DP,7.85369E-06_DP,7.99666E-06_DP,&
&8.14138E-06_DP/)

TOTPLNK( :,11) = (/&
&2.53767E-09_DP,2.77242E-09_DP,3.02564E-09_DP,3.29851E-09_DP,3.59228E-09_DP,&
&3.90825E-09_DP,4.24777E-09_DP,4.61227E-09_DP,5.00322E-09_DP,5.42219E-09_DP,&
&5.87080E-09_DP,6.35072E-09_DP,6.86370E-09_DP,7.41159E-09_DP,7.99628E-09_DP,&
&8.61974E-09_DP,9.28404E-09_DP,9.99130E-09_DP,1.07437E-08_DP,1.15436E-08_DP,&
&1.23933E-08_DP,1.32953E-08_DP,1.42522E-08_DP,1.52665E-08_DP,1.63410E-08_DP,&
&1.74786E-08_DP,1.86820E-08_DP,1.99542E-08_DP,2.12985E-08_DP,2.27179E-08_DP,&
&2.42158E-08_DP,2.57954E-08_DP,2.74604E-08_DP,2.92141E-08_DP,3.10604E-08_DP,&
&3.30029E-08_DP,3.50457E-08_DP,3.71925E-08_DP,3.94476E-08_DP,4.18149E-08_DP,&
&4.42991E-08_DP,4.69043E-08_DP,4.96352E-08_DP,5.24961E-08_DP,5.54921E-08_DP,&
&5.86277E-08_DP,6.19081E-08_DP,6.53381E-08_DP,6.89231E-08_DP,7.26681E-08_DP,&
&7.65788E-08_DP,8.06604E-08_DP,8.49187E-08_DP,8.93591E-08_DP,9.39879E-08_DP,&
&9.88106E-08_DP,1.03834E-07_DP,1.09063E-07_DP,1.14504E-07_DP,1.20165E-07_DP,&
&1.26051E-07_DP,1.32169E-07_DP,1.38525E-07_DP,1.45128E-07_DP,1.51982E-07_DP,&
&1.59096E-07_DP,1.66477E-07_DP,1.74132E-07_DP,1.82068E-07_DP,1.90292E-07_DP,&
&1.98813E-07_DP,2.07638E-07_DP,2.16775E-07_DP,2.26231E-07_DP,2.36015E-07_DP,&
&2.46135E-07_DP,2.56599E-07_DP,2.67415E-07_DP,2.78592E-07_DP,2.90137E-07_DP,&
&3.02061E-07_DP,3.14371E-07_DP,3.27077E-07_DP,3.40186E-07_DP,3.53710E-07_DP,&
&3.67655E-07_DP,3.82031E-07_DP,3.96848E-07_DP,4.12116E-07_DP,4.27842E-07_DP,&
&4.44039E-07_DP,4.60713E-07_DP,4.77876E-07_DP,4.95537E-07_DP,5.13706E-07_DP,&
&5.32392E-07_DP,5.51608E-07_DP,5.71360E-07_DP,5.91662E-07_DP,6.12521E-07_DP,&
&6.33950E-07_DP,6.55958E-07_DP,6.78556E-07_DP,7.01753E-07_DP,7.25562E-07_DP,&
&7.49992E-07_DP,7.75055E-07_DP,8.00760E-07_DP,8.27120E-07_DP,8.54145E-07_DP,&
&8.81845E-07_DP,9.10233E-07_DP,9.39318E-07_DP,9.69113E-07_DP,9.99627E-07_DP,&
&1.03087E-06_DP,1.06286E-06_DP,1.09561E-06_DP,1.12912E-06_DP,1.16340E-06_DP,&
&1.19848E-06_DP,1.23435E-06_DP,1.27104E-06_DP,1.30855E-06_DP,1.34690E-06_DP,&
&1.38609E-06_DP,1.42614E-06_DP,1.46706E-06_DP,1.50886E-06_DP,1.55155E-06_DP,&
&1.59515E-06_DP,1.63967E-06_DP,1.68512E-06_DP,1.73150E-06_DP,1.77884E-06_DP,&
&1.82715E-06_DP,1.87643E-06_DP,1.92670E-06_DP,1.97797E-06_DP,2.03026E-06_DP,&
&2.08356E-06_DP,2.13791E-06_DP,2.19330E-06_DP,2.24975E-06_DP,2.30728E-06_DP,&
&2.36589E-06_DP,2.42560E-06_DP,2.48641E-06_DP,2.54835E-06_DP,2.61142E-06_DP,&
&2.67563E-06_DP,2.74100E-06_DP,2.80754E-06_DP,2.87526E-06_DP,2.94417E-06_DP,&
&3.01429E-06_DP,3.08562E-06_DP,3.15819E-06_DP,3.23199E-06_DP,3.30704E-06_DP,&
&3.38336E-06_DP,3.46096E-06_DP,3.53984E-06_DP,3.62002E-06_DP,3.70151E-06_DP,&
&3.78433E-06_DP,3.86848E-06_DP,3.95399E-06_DP,4.04084E-06_DP,4.12907E-06_DP,&
&4.21868E-06_DP,4.30968E-06_DP,4.40209E-06_DP,4.49592E-06_DP,4.59117E-06_DP,&
&4.68786E-06_DP,4.78600E-06_DP,4.88561E-06_DP,4.98669E-06_DP,5.08926E-06_DP,&
&5.19332E-06_DP/)

TOTPLNK( :,12) = (/&
&2.73921E-10_DP,3.04500E-10_DP,3.38056E-10_DP,3.74835E-10_DP,4.15099E-10_DP,&
&4.59126E-10_DP,5.07214E-10_DP,5.59679E-10_DP,6.16857E-10_DP,6.79103E-10_DP,&
&7.46796E-10_DP,8.20335E-10_DP,9.00144E-10_DP,9.86671E-10_DP,1.08039E-09_DP,&
&1.18180E-09_DP,1.29142E-09_DP,1.40982E-09_DP,1.53757E-09_DP,1.67529E-09_DP,&
&1.82363E-09_DP,1.98327E-09_DP,2.15492E-09_DP,2.33932E-09_DP,2.53726E-09_DP,&
&2.74957E-09_DP,2.97710E-09_DP,3.22075E-09_DP,3.48145E-09_DP,3.76020E-09_DP,&
&4.05801E-09_DP,4.37595E-09_DP,4.71513E-09_DP,5.07672E-09_DP,5.46193E-09_DP,&
&5.87201E-09_DP,6.30827E-09_DP,6.77205E-09_DP,7.26480E-09_DP,7.78794E-09_DP,&
&8.34304E-09_DP,8.93163E-09_DP,9.55537E-09_DP,1.02159E-08_DP,1.09151E-08_DP,&
&1.16547E-08_DP,1.24365E-08_DP,1.32625E-08_DP,1.41348E-08_DP,1.50554E-08_DP,&
&1.60264E-08_DP,1.70500E-08_DP,1.81285E-08_DP,1.92642E-08_DP,2.04596E-08_DP,&
&2.17171E-08_DP,2.30394E-08_DP,2.44289E-08_DP,2.58885E-08_DP,2.74209E-08_DP,&
&2.90290E-08_DP,3.07157E-08_DP,3.24841E-08_DP,3.43371E-08_DP,3.62782E-08_DP,&
&3.83103E-08_DP,4.04371E-08_DP,4.26617E-08_DP,4.49878E-08_DP,4.74190E-08_DP,&
&4.99589E-08_DP,5.26113E-08_DP,5.53801E-08_DP,5.82692E-08_DP,6.12826E-08_DP,&
&6.44245E-08_DP,6.76991E-08_DP,7.11105E-08_DP,7.46634E-08_DP,7.83621E-08_DP,&
&8.22112E-08_DP,8.62154E-08_DP,9.03795E-08_DP,9.47081E-08_DP,9.92066E-08_DP,&
&1.03879E-07_DP,1.08732E-07_DP,1.13770E-07_DP,1.18998E-07_DP,1.24422E-07_DP,&
&1.30048E-07_DP,1.35880E-07_DP,1.41924E-07_DP,1.48187E-07_DP,1.54675E-07_DP,&
&1.61392E-07_DP,1.68346E-07_DP,1.75543E-07_DP,1.82988E-07_DP,1.90688E-07_DP,&
&1.98650E-07_DP,2.06880E-07_DP,2.15385E-07_DP,2.24172E-07_DP,2.33247E-07_DP,&
&2.42617E-07_DP,2.52289E-07_DP,2.62272E-07_DP,2.72571E-07_DP,2.83193E-07_DP,&
&2.94147E-07_DP,3.05440E-07_DP,3.17080E-07_DP,3.29074E-07_DP,3.41430E-07_DP,&
&3.54155E-07_DP,3.67259E-07_DP,3.80747E-07_DP,3.94631E-07_DP,4.08916E-07_DP,&
&4.23611E-07_DP,4.38725E-07_DP,4.54267E-07_DP,4.70245E-07_DP,4.86666E-07_DP,&
&5.03541E-07_DP,5.20879E-07_DP,5.38687E-07_DP,5.56975E-07_DP,5.75751E-07_DP,&
&5.95026E-07_DP,6.14808E-07_DP,6.35107E-07_DP,6.55932E-07_DP,6.77293E-07_DP,&
&6.99197E-07_DP,7.21656E-07_DP,7.44681E-07_DP,7.68278E-07_DP,7.92460E-07_DP,&
&8.17235E-07_DP,8.42614E-07_DP,8.68606E-07_DP,8.95223E-07_DP,9.22473E-07_DP,&
&9.50366E-07_DP,9.78915E-07_DP,1.00813E-06_DP,1.03802E-06_DP,1.06859E-06_DP,&
&1.09986E-06_DP,1.13184E-06_DP,1.16453E-06_DP,1.19796E-06_DP,1.23212E-06_DP,&
&1.26703E-06_DP,1.30270E-06_DP,1.33915E-06_DP,1.37637E-06_DP,1.41440E-06_DP,&
&1.45322E-06_DP,1.49286E-06_DP,1.53333E-06_DP,1.57464E-06_DP,1.61679E-06_DP,&
&1.65981E-06_DP,1.70370E-06_DP,1.74847E-06_DP,1.79414E-06_DP,1.84071E-06_DP,&
&1.88821E-06_DP,1.93663E-06_DP,1.98599E-06_DP,2.03631E-06_DP,2.08759E-06_DP,&
&2.13985E-06_DP,2.19310E-06_DP,2.24734E-06_DP,2.30260E-06_DP,2.35888E-06_DP,&
&2.41619E-06_DP/)

TOTPLNK( :,13) = (/&
&4.53634E-11_DP,5.11435E-11_DP,5.75754E-11_DP,6.47222E-11_DP,7.26531E-11_DP,&
&8.14420E-11_DP,9.11690E-11_DP,1.01921E-10_DP,1.13790E-10_DP,1.26877E-10_DP,&
&1.41288E-10_DP,1.57140E-10_DP,1.74555E-10_DP,1.93665E-10_DP,2.14613E-10_DP,&
&2.37548E-10_DP,2.62633E-10_DP,2.90039E-10_DP,3.19948E-10_DP,3.52558E-10_DP,&
&3.88073E-10_DP,4.26716E-10_DP,4.68719E-10_DP,5.14331E-10_DP,5.63815E-10_DP,&
&6.17448E-10_DP,6.75526E-10_DP,7.38358E-10_DP,8.06277E-10_DP,8.79625E-10_DP,&
&9.58770E-10_DP,1.04410E-09_DP,1.13602E-09_DP,1.23495E-09_DP,1.34135E-09_DP,&
&1.45568E-09_DP,1.57845E-09_DP,1.71017E-09_DP,1.85139E-09_DP,2.00268E-09_DP,&
&2.16464E-09_DP,2.33789E-09_DP,2.52309E-09_DP,2.72093E-09_DP,2.93212E-09_DP,&
&3.15740E-09_DP,3.39757E-09_DP,3.65341E-09_DP,3.92579E-09_DP,4.21559E-09_DP,&
&4.52372E-09_DP,4.85115E-09_DP,5.19886E-09_DP,5.56788E-09_DP,5.95928E-09_DP,&
&6.37419E-09_DP,6.81375E-09_DP,7.27917E-09_DP,7.77168E-09_DP,8.29256E-09_DP,&
&8.84317E-09_DP,9.42487E-09_DP,1.00391E-08_DP,1.06873E-08_DP,1.13710E-08_DP,&
&1.20919E-08_DP,1.28515E-08_DP,1.36514E-08_DP,1.44935E-08_DP,1.53796E-08_DP,&
&1.63114E-08_DP,1.72909E-08_DP,1.83201E-08_DP,1.94008E-08_DP,2.05354E-08_DP,&
&2.17258E-08_DP,2.29742E-08_DP,2.42830E-08_DP,2.56545E-08_DP,2.70910E-08_DP,&
&2.85950E-08_DP,3.01689E-08_DP,3.18155E-08_DP,3.35373E-08_DP,3.53372E-08_DP,&
&3.72177E-08_DP,3.91818E-08_DP,4.12325E-08_DP,4.33727E-08_DP,4.56056E-08_DP,&
&4.79342E-08_DP,5.03617E-08_DP,5.28915E-08_DP,5.55270E-08_DP,5.82715E-08_DP,&
&6.11286E-08_DP,6.41019E-08_DP,6.71951E-08_DP,7.04119E-08_DP,7.37560E-08_DP,&
&7.72315E-08_DP,8.08424E-08_DP,8.45927E-08_DP,8.84866E-08_DP,9.25281E-08_DP,&
&9.67218E-08_DP,1.01072E-07_DP,1.05583E-07_DP,1.10260E-07_DP,1.15107E-07_DP,&
&1.20128E-07_DP,1.25330E-07_DP,1.30716E-07_DP,1.36291E-07_DP,1.42061E-07_DP,&
&1.48031E-07_DP,1.54206E-07_DP,1.60592E-07_DP,1.67192E-07_DP,1.74015E-07_DP,&
&1.81064E-07_DP,1.88345E-07_DP,1.95865E-07_DP,2.03628E-07_DP,2.11643E-07_DP,&
&2.19912E-07_DP,2.28443E-07_DP,2.37244E-07_DP,2.46318E-07_DP,2.55673E-07_DP,&
&2.65316E-07_DP,2.75252E-07_DP,2.85489E-07_DP,2.96033E-07_DP,3.06891E-07_DP,&
&3.18070E-07_DP,3.29576E-07_DP,3.41417E-07_DP,3.53600E-07_DP,3.66133E-07_DP,&
&3.79021E-07_DP,3.92274E-07_DP,4.05897E-07_DP,4.19899E-07_DP,4.34288E-07_DP,&
&4.49071E-07_DP,4.64255E-07_DP,4.79850E-07_DP,4.95863E-07_DP,5.12300E-07_DP,&
&5.29172E-07_DP,5.46486E-07_DP,5.64250E-07_DP,5.82473E-07_DP,6.01164E-07_DP,&
&6.20329E-07_DP,6.39979E-07_DP,6.60122E-07_DP,6.80767E-07_DP,7.01922E-07_DP,&
&7.23596E-07_DP,7.45800E-07_DP,7.68539E-07_DP,7.91826E-07_DP,8.15669E-07_DP,&
&8.40076E-07_DP,8.65058E-07_DP,8.90623E-07_DP,9.16783E-07_DP,9.43544E-07_DP,&
&9.70917E-07_DP,9.98912E-07_DP,1.02754E-06_DP,1.05681E-06_DP,1.08673E-06_DP,&
&1.11731E-06_DP,1.14856E-06_DP,1.18050E-06_DP,1.21312E-06_DP,1.24645E-06_DP,&
&1.28049E-06_DP/)

TOTPLNK( :,14) = (/&
&1.40113E-11_DP,1.59358E-11_DP,1.80960E-11_DP,2.05171E-11_DP,2.32266E-11_DP,&
&2.62546E-11_DP,2.96335E-11_DP,3.33990E-11_DP,3.75896E-11_DP,4.22469E-11_DP,&
&4.74164E-11_DP,5.31466E-11_DP,5.94905E-11_DP,6.65054E-11_DP,7.42522E-11_DP,&
&8.27975E-11_DP,9.22122E-11_DP,1.02573E-10_DP,1.13961E-10_DP,1.26466E-10_DP,&
&1.40181E-10_DP,1.55206E-10_DP,1.71651E-10_DP,1.89630E-10_DP,2.09265E-10_DP,&
&2.30689E-10_DP,2.54040E-10_DP,2.79467E-10_DP,3.07128E-10_DP,3.37190E-10_DP,&
&3.69833E-10_DP,4.05243E-10_DP,4.43623E-10_DP,4.85183E-10_DP,5.30149E-10_DP,&
&5.78755E-10_DP,6.31255E-10_DP,6.87910E-10_DP,7.49002E-10_DP,8.14824E-10_DP,&
&8.85687E-10_DP,9.61914E-10_DP,1.04385E-09_DP,1.13186E-09_DP,1.22631E-09_DP,&
&1.32761E-09_DP,1.43617E-09_DP,1.55243E-09_DP,1.67686E-09_DP,1.80992E-09_DP,&
&1.95212E-09_DP,2.10399E-09_DP,2.26607E-09_DP,2.43895E-09_DP,2.62321E-09_DP,&
&2.81949E-09_DP,3.02844E-09_DP,3.25073E-09_DP,3.48707E-09_DP,3.73820E-09_DP,&
&4.00490E-09_DP,4.28794E-09_DP,4.58819E-09_DP,4.90647E-09_DP,5.24371E-09_DP,&
&5.60081E-09_DP,5.97875E-09_DP,6.37854E-09_DP,6.80120E-09_DP,7.24782E-09_DP,&
&7.71950E-09_DP,8.21740E-09_DP,8.74271E-09_DP,9.29666E-09_DP,9.88054E-09_DP,&
&1.04956E-08_DP,1.11434E-08_DP,1.18251E-08_DP,1.25422E-08_DP,1.32964E-08_DP,&
&1.40890E-08_DP,1.49217E-08_DP,1.57961E-08_DP,1.67140E-08_DP,1.76771E-08_DP,&
&1.86870E-08_DP,1.97458E-08_DP,2.08553E-08_DP,2.20175E-08_DP,2.32342E-08_DP,&
&2.45077E-08_DP,2.58401E-08_DP,2.72334E-08_DP,2.86900E-08_DP,3.02122E-08_DP,&
&3.18021E-08_DP,3.34624E-08_DP,3.51954E-08_DP,3.70037E-08_DP,3.88899E-08_DP,&
&4.08568E-08_DP,4.29068E-08_DP,4.50429E-08_DP,4.72678E-08_DP,4.95847E-08_DP,&
&5.19963E-08_DP,5.45058E-08_DP,5.71161E-08_DP,5.98309E-08_DP,6.26529E-08_DP,&
&6.55857E-08_DP,6.86327E-08_DP,7.17971E-08_DP,7.50829E-08_DP,7.84933E-08_DP,&
&8.20323E-08_DP,8.57035E-08_DP,8.95105E-08_DP,9.34579E-08_DP,9.75488E-08_DP,&
&1.01788E-07_DP,1.06179E-07_DP,1.10727E-07_DP,1.15434E-07_DP,1.20307E-07_DP,&
&1.25350E-07_DP,1.30566E-07_DP,1.35961E-07_DP,1.41539E-07_DP,1.47304E-07_DP,&
&1.53263E-07_DP,1.59419E-07_DP,1.65778E-07_DP,1.72345E-07_DP,1.79124E-07_DP,&
&1.86122E-07_DP,1.93343E-07_DP,2.00792E-07_DP,2.08476E-07_DP,2.16400E-07_DP,&
&2.24568E-07_DP,2.32988E-07_DP,2.41666E-07_DP,2.50605E-07_DP,2.59813E-07_DP,&
&2.69297E-07_DP,2.79060E-07_DP,2.89111E-07_DP,2.99455E-07_DP,3.10099E-07_DP,&
&3.21049E-07_DP,3.32311E-07_DP,3.43893E-07_DP,3.55801E-07_DP,3.68041E-07_DP,&
&3.80621E-07_DP,3.93547E-07_DP,4.06826E-07_DP,4.20465E-07_DP,4.34473E-07_DP,&
&4.48856E-07_DP,4.63620E-07_DP,4.78774E-07_DP,4.94325E-07_DP,5.10280E-07_DP,&
&5.26648E-07_DP,5.43436E-07_DP,5.60652E-07_DP,5.78302E-07_DP,5.96397E-07_DP,&
&6.14943E-07_DP,6.33949E-07_DP,6.53421E-07_DP,6.73370E-07_DP,6.93803E-07_DP,&
&7.14731E-07_DP,7.36157E-07_DP,7.58095E-07_DP,7.80549E-07_DP,8.03533E-07_DP,&
&8.27050E-07_DP/)

TOTPLNK( :,15) = (/&
&3.90483E-12_DP,4.47999E-12_DP,5.13122E-12_DP,5.86739E-12_DP,6.69829E-12_DP,&
&7.63467E-12_DP,8.68833E-12_DP,9.87221E-12_DP,1.12005E-11_DP,1.26885E-11_DP,&
&1.43534E-11_DP,1.62134E-11_DP,1.82888E-11_DP,2.06012E-11_DP,2.31745E-11_DP,&
&2.60343E-11_DP,2.92087E-11_DP,3.27277E-11_DP,3.66242E-11_DP,4.09334E-11_DP,&
&4.56935E-11_DP,5.09455E-11_DP,5.67338E-11_DP,6.31057E-11_DP,7.01127E-11_DP,&
&7.78096E-11_DP,8.62554E-11_DP,9.55130E-11_DP,1.05651E-10_DP,1.16740E-10_DP,&
&1.28858E-10_DP,1.42089E-10_DP,1.56519E-10_DP,1.72243E-10_DP,1.89361E-10_DP,&
&2.07978E-10_DP,2.28209E-10_DP,2.50173E-10_DP,2.73999E-10_DP,2.99820E-10_DP,&
&3.27782E-10_DP,3.58034E-10_DP,3.90739E-10_DP,4.26067E-10_DP,4.64196E-10_DP,&
&5.05317E-10_DP,5.49631E-10_DP,5.97347E-10_DP,6.48689E-10_DP,7.03891E-10_DP,&
&7.63201E-10_DP,8.26876E-10_DP,8.95192E-10_DP,9.68430E-10_DP,1.04690E-09_DP,&
&1.13091E-09_DP,1.22079E-09_DP,1.31689E-09_DP,1.41957E-09_DP,1.52922E-09_DP,&
&1.64623E-09_DP,1.77101E-09_DP,1.90401E-09_DP,2.04567E-09_DP,2.19647E-09_DP,&
&2.35690E-09_DP,2.52749E-09_DP,2.70875E-09_DP,2.90127E-09_DP,3.10560E-09_DP,&
&3.32238E-09_DP,3.55222E-09_DP,3.79578E-09_DP,4.05375E-09_DP,4.32682E-09_DP,&
&4.61574E-09_DP,4.92128E-09_DP,5.24420E-09_DP,5.58536E-09_DP,5.94558E-09_DP,&
&6.32575E-09_DP,6.72678E-09_DP,7.14964E-09_DP,7.59526E-09_DP,8.06470E-09_DP,&
&8.55897E-09_DP,9.07916E-09_DP,9.62638E-09_DP,1.02018E-08_DP,1.08066E-08_DP,&
&1.14420E-08_DP,1.21092E-08_DP,1.28097E-08_DP,1.35446E-08_DP,1.43155E-08_DP,&
&1.51237E-08_DP,1.59708E-08_DP,1.68581E-08_DP,1.77873E-08_DP,1.87599E-08_DP,&
&1.97777E-08_DP,2.08423E-08_DP,2.19555E-08_DP,2.31190E-08_DP,2.43348E-08_DP,&
&2.56045E-08_DP,2.69302E-08_DP,2.83140E-08_DP,2.97578E-08_DP,3.12636E-08_DP,&
&3.28337E-08_DP,3.44702E-08_DP,3.61755E-08_DP,3.79516E-08_DP,3.98012E-08_DP,&
&4.17265E-08_DP,4.37300E-08_DP,4.58143E-08_DP,4.79819E-08_DP,5.02355E-08_DP,&
&5.25777E-08_DP,5.50114E-08_DP,5.75393E-08_DP,6.01644E-08_DP,6.28896E-08_DP,&
&6.57177E-08_DP,6.86521E-08_DP,7.16959E-08_DP,7.48520E-08_DP,7.81239E-08_DP,&
&8.15148E-08_DP,8.50282E-08_DP,8.86675E-08_DP,9.24362E-08_DP,9.63380E-08_DP,&
&1.00376E-07_DP,1.04555E-07_DP,1.08878E-07_DP,1.13349E-07_DP,1.17972E-07_DP,&
&1.22751E-07_DP,1.27690E-07_DP,1.32793E-07_DP,1.38064E-07_DP,1.43508E-07_DP,&
&1.49129E-07_DP,1.54931E-07_DP,1.60920E-07_DP,1.67099E-07_DP,1.73473E-07_DP,&
&1.80046E-07_DP,1.86825E-07_DP,1.93812E-07_DP,2.01014E-07_DP,2.08436E-07_DP,&
&2.16082E-07_DP,2.23957E-07_DP,2.32067E-07_DP,2.40418E-07_DP,2.49013E-07_DP,&
&2.57860E-07_DP,2.66963E-07_DP,2.76328E-07_DP,2.85961E-07_DP,2.95868E-07_DP,&
&3.06053E-07_DP,3.16524E-07_DP,3.27286E-07_DP,3.38345E-07_DP,3.49707E-07_DP,&
&3.61379E-07_DP,3.73367E-07_DP,3.85676E-07_DP,3.98315E-07_DP,4.11287E-07_DP,&
&4.24602E-07_DP,4.38265E-07_DP,4.52283E-07_DP,4.66662E-07_DP,4.81410E-07_DP,&
&4.96535E-07_DP/)

TOTPLNK( :,16) = (/&
&4.65378E-13_DP,5.41927E-13_DP,6.29913E-13_DP,7.30869E-13_DP,8.46510E-13_DP,&
&9.78750E-13_DP,1.12972E-12_DP,1.30181E-12_DP,1.49764E-12_DP,1.72016E-12_DP,&
&1.97260E-12_DP,2.25858E-12_DP,2.58206E-12_DP,2.94744E-12_DP,3.35955E-12_DP,&
&3.82372E-12_DP,4.34581E-12_DP,4.93225E-12_DP,5.59010E-12_DP,6.32711E-12_DP,&
&7.15171E-12_DP,8.07317E-12_DP,9.10159E-12_DP,1.02480E-11_DP,1.15244E-11_DP,&
&1.29438E-11_DP,1.45204E-11_DP,1.62697E-11_DP,1.82084E-11_DP,2.03545E-11_DP,&
&2.27278E-11_DP,2.53494E-11_DP,2.82424E-11_DP,3.14313E-11_DP,3.49431E-11_DP,&
&3.88064E-11_DP,4.30522E-11_DP,4.77139E-11_DP,5.28273E-11_DP,5.84308E-11_DP,&
&6.45658E-11_DP,7.12764E-11_DP,7.86103E-11_DP,8.66176E-11_DP,9.53534E-11_DP,&
&1.04875E-10_DP,1.15245E-10_DP,1.26528E-10_DP,1.38796E-10_DP,1.52123E-10_DP,&
&1.66590E-10_DP,1.82281E-10_DP,1.99287E-10_DP,2.17704E-10_DP,2.37632E-10_DP,&
&2.59182E-10_DP,2.82468E-10_DP,3.07610E-10_DP,3.34738E-10_DP,3.63988E-10_DP,&
&3.95504E-10_DP,4.29438E-10_DP,4.65951E-10_DP,5.05212E-10_DP,5.47402E-10_DP,&
&5.92707E-10_DP,6.41329E-10_DP,6.93477E-10_DP,7.49371E-10_DP,8.09242E-10_DP,&
&8.73338E-10_DP,9.41911E-10_DP,1.01524E-09_DP,1.09359E-09_DP,1.17728E-09_DP,&
&1.26660E-09_DP,1.36190E-09_DP,1.46350E-09_DP,1.57177E-09_DP,1.68709E-09_DP,&
&1.80984E-09_DP,1.94044E-09_DP,2.07932E-09_DP,2.22693E-09_DP,2.38373E-09_DP,&
&2.55021E-09_DP,2.72689E-09_DP,2.91429E-09_DP,3.11298E-09_DP,3.32353E-09_DP,&
&3.54655E-09_DP,3.78265E-09_DP,4.03251E-09_DP,4.29679E-09_DP,4.57620E-09_DP,&
&4.87148E-09_DP,5.18341E-09_DP,5.51276E-09_DP,5.86037E-09_DP,6.22708E-09_DP,&
&6.61381E-09_DP,7.02145E-09_DP,7.45097E-09_DP,7.90336E-09_DP,8.37967E-09_DP,&
&8.88092E-09_DP,9.40827E-09_DP,9.96280E-09_DP,1.05457E-08_DP,1.11583E-08_DP,&
&1.18017E-08_DP,1.24773E-08_DP,1.31865E-08_DP,1.39306E-08_DP,1.47111E-08_DP,&
&1.55295E-08_DP,1.63872E-08_DP,1.72860E-08_DP,1.82274E-08_DP,1.92132E-08_DP,&
&2.02450E-08_DP,2.13247E-08_DP,2.24541E-08_DP,2.36352E-08_DP,2.48699E-08_DP,&
&2.61602E-08_DP,2.75082E-08_DP,2.89161E-08_DP,3.03860E-08_DP,3.19203E-08_DP,&
&3.35213E-08_DP,3.51913E-08_DP,3.69330E-08_DP,3.87486E-08_DP,4.06411E-08_DP,&
&4.26129E-08_DP,4.46668E-08_DP,4.68058E-08_DP,4.90325E-08_DP,5.13502E-08_DP,&
&5.37617E-08_DP,5.62703E-08_DP,5.88791E-08_DP,6.15915E-08_DP,6.44107E-08_DP,&
&6.73404E-08_DP,7.03841E-08_DP,7.35453E-08_DP,7.68278E-08_DP,8.02355E-08_DP,&
&8.37721E-08_DP,8.74419E-08_DP,9.12486E-08_DP,9.51968E-08_DP,9.92905E-08_DP,&
&1.03534E-07_DP,1.07932E-07_DP,1.12490E-07_DP,1.17211E-07_DP,1.22100E-07_DP,&
&1.27163E-07_DP,1.32404E-07_DP,1.37829E-07_DP,1.43443E-07_DP,1.49250E-07_DP,&
&1.55257E-07_DP,1.61470E-07_DP,1.67893E-07_DP,1.74532E-07_DP,1.81394E-07_DP,&
&1.88485E-07_DP,1.95810E-07_DP,2.03375E-07_DP,2.11189E-07_DP,2.19256E-07_DP,&
&2.27583E-07_DP,2.36177E-07_DP,2.45046E-07_DP,2.54196E-07_DP,2.63634E-07_DP,&
&2.73367E-07_DP/)

TOTPLK16( :) = (/&
&4.46128E-13_DP,5.19008E-13_DP,6.02681E-13_DP,6.98580E-13_DP,8.08302E-13_DP,&
&9.33629E-13_DP,1.07654E-12_DP,1.23925E-12_DP,1.42419E-12_DP,1.63407E-12_DP,&
&1.87190E-12_DP,2.14099E-12_DP,2.44498E-12_DP,2.78793E-12_DP,3.17424E-12_DP,&
&3.60881E-12_DP,4.09698E-12_DP,4.64461E-12_DP,5.25813E-12_DP,5.94456E-12_DP,&
&6.71156E-12_DP,7.56752E-12_DP,8.52154E-12_DP,9.58357E-12_DP,1.07644E-11_DP,&
&1.20758E-11_DP,1.35304E-11_DP,1.51420E-11_DP,1.69256E-11_DP,1.88973E-11_DP,&
&2.10746E-11_DP,2.34762E-11_DP,2.61227E-11_DP,2.90356E-11_DP,3.22388E-11_DP,&
&3.57574E-11_DP,3.96187E-11_DP,4.38519E-11_DP,4.84883E-11_DP,5.35616E-11_DP,&
&5.91075E-11_DP,6.51647E-11_DP,7.17743E-11_DP,7.89797E-11_DP,8.68284E-11_DP,&
&9.53697E-11_DP,1.04658E-10_DP,1.14748E-10_DP,1.25701E-10_DP,1.37582E-10_DP,&
&1.50457E-10_DP,1.64400E-10_DP,1.79487E-10_DP,1.95799E-10_DP,2.13422E-10_DP,&
&2.32446E-10_DP,2.52970E-10_DP,2.75094E-10_DP,2.98925E-10_DP,3.24578E-10_DP,&
&3.52172E-10_DP,3.81833E-10_DP,4.13695E-10_DP,4.47897E-10_DP,4.84588E-10_DP,&
&5.23922E-10_DP,5.66063E-10_DP,6.11182E-10_DP,6.59459E-10_DP,7.11081E-10_DP,&
&7.66251E-10_DP,8.25172E-10_DP,8.88065E-10_DP,9.55155E-10_DP,1.02668E-09_DP,&
&1.10290E-09_DP,1.18406E-09_DP,1.27044E-09_DP,1.36233E-09_DP,1.46002E-09_DP,&
&1.56382E-09_DP,1.67406E-09_DP,1.79108E-09_DP,1.91522E-09_DP,2.04686E-09_DP,&
&2.18637E-09_DP,2.33416E-09_DP,2.49063E-09_DP,2.65622E-09_DP,2.83136E-09_DP,&
&3.01653E-09_DP,3.21221E-09_DP,3.41890E-09_DP,3.63712E-09_DP,3.86740E-09_DP,&
&4.11030E-09_DP,4.36641E-09_DP,4.63631E-09_DP,4.92064E-09_DP,5.22003E-09_DP,&
&5.53516E-09_DP,5.86670E-09_DP,6.21538E-09_DP,6.58191E-09_DP,6.96708E-09_DP,&
&7.37165E-09_DP,7.79645E-09_DP,8.24229E-09_DP,8.71007E-09_DP,9.20066E-09_DP,&
&9.71498E-09_DP,1.02540E-08_DP,1.08186E-08_DP,1.14100E-08_DP,1.20290E-08_DP,&
&1.26767E-08_DP,1.33544E-08_DP,1.40630E-08_DP,1.48038E-08_DP,1.55780E-08_DP,&
&1.63867E-08_DP,1.72313E-08_DP,1.81130E-08_DP,1.90332E-08_DP,1.99932E-08_DP,&
&2.09945E-08_DP,2.20385E-08_DP,2.31267E-08_DP,2.42605E-08_DP,2.54416E-08_DP,&
&2.66716E-08_DP,2.79520E-08_DP,2.92846E-08_DP,3.06711E-08_DP,3.21133E-08_DP,&
&3.36128E-08_DP,3.51717E-08_DP,3.67918E-08_DP,3.84749E-08_DP,4.02232E-08_DP,&
&4.20386E-08_DP,4.39231E-08_DP,4.58790E-08_DP,4.79083E-08_DP,5.00132E-08_DP,&
&5.21961E-08_DP,5.44592E-08_DP,5.68049E-08_DP,5.92356E-08_DP,6.17537E-08_DP,&
&6.43617E-08_DP,6.70622E-08_DP,6.98578E-08_DP,7.27511E-08_DP,7.57449E-08_DP,&
&7.88419E-08_DP,8.20449E-08_DP,8.53568E-08_DP,8.87805E-08_DP,9.23190E-08_DP,&
&9.59753E-08_DP,9.97526E-08_DP,1.03654E-07_DP,1.07682E-07_DP,1.11841E-07_DP,&
&1.16134E-07_DP,1.20564E-07_DP,1.25135E-07_DP,1.29850E-07_DP,1.34712E-07_DP,&
&1.39726E-07_DP,1.44894E-07_DP,1.50221E-07_DP,1.55711E-07_DP,1.61367E-07_DP,&
&1.67193E-07_DP,1.73193E-07_DP,1.79371E-07_DP,1.85732E-07_DP,1.92279E-07_DP,&
&1.99016E-07_DP/)
!     -----------------------------------------------------------------

   RETURN
   END SUBROUTINE rad_lon_SURRTPK
  ! ===========================================================================

  ! ===========================================================================
   SUBROUTINE rad_lon_SURRTRF

     !     Adapted from Eli J. Mlawer, Atmospheric & Environmental Research.
     !     by JJMorcrette, ECMWF
     !     ------------------------------------------------------------------

     !     mag, MPI, 25 February 2000: comment added

     !     ------------------------------------------------------------------

     !     These pressures are chosen such that the ln of the first pressure
     !     has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
     !     each subsequent ln(pressure) differs from the previous one by 0.2.

     IMPLICIT NONE

PREF( :) = (/&
    &1.05363E+03_DP,8.62642E+02_DP,7.06272E+02_DP,5.78246E+02_DP,4.73428E+02_DP,&
    &3.87610E+02_DP,3.17348E+02_DP,2.59823E+02_DP,2.12725E+02_DP,1.74164E+02_DP,&
    &1.42594E+02_DP,1.16746E+02_DP,9.55835E+01_DP,7.82571E+01_DP,6.40715E+01_DP,&
    &5.24573E+01_DP,4.29484E+01_DP,3.51632E+01_DP,2.87892E+01_DP,2.35706E+01_DP,&
    &1.92980E+01_DP,1.57998E+01_DP,1.29358E+01_DP,1.05910E+01_DP,8.67114E+00_DP,&
    &7.09933E+00_DP,5.81244E+00_DP,4.75882E+00_DP,3.89619E+00_DP,3.18993E+00_DP,&
    &2.61170E+00_DP,2.13828E+00_DP,1.75067E+00_DP,1.43333E+00_DP,1.17351E+00_DP,&
    &9.60789E-01_DP,7.86628E-01_DP,6.44036E-01_DP,5.27292E-01_DP,4.31710E-01_DP,&
    &3.53455E-01_DP,2.89384E-01_DP,2.36928E-01_DP,1.93980E-01_DP,1.58817E-01_DP,&
    &1.30029E-01_DP,1.06458E-01_DP,8.71608E-02_DP,7.13612E-02_DP,5.84256E-02_DP,&
    &4.78349E-02_DP,3.91639E-02_DP,3.20647E-02_DP,2.62523E-02_DP,2.14936E-02_DP,&
    &1.75975E-02_DP,1.44076E-02_DP,1.17959E-02_DP,9.65769E-03_DP/)

PREFLOG( :) = (/&
     &6.9600E+00_DP, 6.7600E+00_DP, 6.5600E+00_DP, 6.3600E+00_DP, 6.1600E+00_DP,&
     &5.9600E+00_DP, 5.7600E+00_DP, 5.5600E+00_DP, 5.3600E+00_DP, 5.1600E+00_DP,&
     &4.9600E+00_DP, 4.7600E+00_DP, 4.5600E+00_DP, 4.3600E+00_DP, 4.1600E+00_DP,&
     &3.9600E+00_DP, 3.7600E+00_DP, 3.5600E+00_DP, 3.3600E+00_DP, 3.1600E+00_DP,&
     &2.9600E+00_DP, 2.7600E+00_DP, 2.5600E+00_DP, 2.3600E+00_DP, 2.1600E+00_DP,&
     &1.9600E+00_DP, 1.7600E+00_DP, 1.5600E+00_DP, 1.3600E+00_DP, 1.1600E+00_DP,&
     &9.6000E-01_DP, 7.6000E-01_DP, 5.6000E-01_DP, 3.6000E-01_DP, 1.6000E-01_DP,&
    &-4.0000E-02_DP,-2.4000E-01_DP,-4.4000E-01_DP,-6.4000E-01_DP,-8.4000E-01_DP,&
    &-1.0400E+00_DP,-1.2400E+00_DP,-1.4400E+00_DP,-1.6400E+00_DP,-1.8400E+00_DP,&
    &-2.0400E+00_DP,-2.2400E+00_DP,-2.4400E+00_DP,-2.6400E+00_DP,-2.8400E+00_DP,&
    &-3.0400E+00_DP,-3.2400E+00_DP,-3.4400E+00_DP,-3.6400E+00_DP,-3.8400E+00_DP,&
    &-4.0400E+00_DP,-4.2400E+00_DP,-4.4400E+00_DP,-4.6400E+00_DP/)

!     These are the temperatures associated with the respective 
!     pressures for the MLS standard atmosphere. 
TREF( :) = (/&
     &2.9420E+02_DP, 2.8799E+02_DP, 2.7894E+02_DP, 2.6925E+02_DP, 2.5983E+02_DP,&
     &2.5017E+02_DP, 2.4077E+02_DP, 2.3179E+02_DP, 2.2306E+02_DP, 2.1578E+02_DP,&
     &2.1570E+02_DP, 2.1570E+02_DP, 2.1570E+02_DP, 2.1706E+02_DP, 2.1858E+02_DP,&
     &2.2018E+02_DP, 2.2174E+02_DP, 2.2328E+02_DP, 2.2479E+02_DP, 2.2655E+02_DP,&
     &2.2834E+02_DP, 2.3113E+02_DP, 2.3401E+02_DP, 2.3703E+02_DP, 2.4022E+02_DP,&
     &2.4371E+02_DP, 2.4726E+02_DP, 2.5085E+02_DP, 2.5457E+02_DP, 2.5832E+02_DP,&
     &2.6216E+02_DP, 2.6606E+02_DP, 2.6999E+02_DP, 2.7340E+02_DP, 2.7536E+02_DP,&
     &2.7568E+02_DP, 2.7372E+02_DP, 2.7163E+02_DP, 2.6955E+02_DP, 2.6593E+02_DP,&
     &2.6211E+02_DP, 2.5828E+02_DP, 2.5360E+02_DP, 2.4854E+02_DP, 2.4348E+02_DP,&
     &2.3809E+02_DP, 2.3206E+02_DP, 2.2603E+02_DP, 2.2000E+02_DP, 2.1435E+02_DP,&
     &2.0887E+02_DP, 2.0340E+02_DP, 1.9792E+02_DP, 1.9290E+02_DP, 1.8809E+02_DP,&
     &1.8329E+02_DP, 1.7849E+02_DP, 1.7394E+02_DP, 1.7212E+02_DP/)

!     -----------------------------------------------------------------
  RETURN
  END SUBROUTINE rad_lon_SURRTRF
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE rad_lon_SURRTFTR

    !     Adapted from Eli J. Mlawer, Atmospheric & Environmental Research.
    !     by JJMorcrette, ECMWF
    !     ------------------------------------------------------------------

    !     mag, MPI, 25 February 2000: comment added

    !     ------------------------------------------------------------------

    IMPLICIT NONE

    NGC( :) = (/8, 14, 16, 14, 16, 8, 12, 8, 12, 6, 8, 8, 4, 2, 2, 2 /)

    NGS( :) = (/&
         &8,  22,  38,  52,  68,  76,  88,  96, &
         &108, 114, 122, 130, 134, 136, 138, 140/)
    NGM( :) = (/&
           &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,             &! Band 1
           &1,2,3,4,5,6,7,8,9,10,11,12,13,13,14,14,      &! Band 2
           &1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,      &! Band 3
           &1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14,      &! Band 4
           &1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,      &! Band 5
           &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,             &! Band 6
           &1,1,2,2,3,4,5,6,7,8,9,10,11,11,12,12,        &! Band 7
           &1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,             &! Band 8
           &1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12,       &! Band 9
           &1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6,             &! Band 10
           &1,2,3,3,4,4,5,5,6,6,7,7,7,8,8,8,             &! Band 11
           &1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8,             &! Band 12
           &1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4,             &! Band 13
           &1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,             &! Band 14
           &1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,             &! Band 15
           &1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2/)           ! Band 16

    NGN( :) = (/&
           &2,2,2,2,2,2,2,2,                             &! Band 1
           &1,1,1,1,1,1,1,1,1,1,1,1,2,2,                 &! Band 2
           &1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,             &! Band 3
           &1,1,1,1,1,1,1,1,1,1,1,1,1,3,                 &! Band 4
           &1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,             &! Band 5
           &2,2,2,2,2,2,2,2,                             &! Band 6
           &2,2,1,1,1,1,1,1,1,1,2,2,                     &! Band 7
           &2,2,2,2,2,2,2,2,                             &! Band 8
           &1,1,1,1,1,1,1,1,2,2,2,2,                     &! Band 9
           &2,2,2,2,4,4,                                 &! Band 10
           &1,1,2,2,2,2,3,3,                             &! Band 11
           &1,1,1,1,2,2,4,4,                             &! Band 12
           &3,3,4,6,                                     &! Band 13
           &8,8,                                         &! Band 14
           &8,8,                                         &! Band 15
           &8,8/)                                       ! Band 16

    NGB( :) = (/&
            &1,1,1,1,1,1,1,1,                            &! Band 1
            &2,2,2,2,2,2,2,2,2,2,2,2,2,2,                &! Band 2
            &3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,            &! Band 3
            &4,4,4,4,4,4,4,4,4,4,4,4,4,4,                &! Band 4
            &5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,            &! Band 5
            &6,6,6,6,6,6,6,6,                            &! Band 6
            &7,7,7,7,7,7,7,7,7,7,7,7,                    &! Band 7
            &8,8,8,8,8,8,8,8,                            &! Band 8
            &9,9,9,9,9,9,9,9,9,9,9,9,                    &! Band 9
           &10,10,10,10,10,10,                           &! Band 10
           &11,11,11,11,11,11,11,11,                     &! Band 11
           &12,12,12,12,12,12,12,12,                     &! Band 12
           &13,13,13,13,                                 &! Band 13
           &14,14,                                       &! Band 14
           &15,15,                                       &! Band 15
           &16,16/)                                     ! Band 16

    WT( :) = (/&
         &0.1527534276_DP,0.1491729617_DP,0.1420961469_DP,0.1316886544_DP,&
         &0.1181945205_DP,0.1019300893_DP,0.0832767040_DP,0.0626720116_DP,&
         &0.0424925_DP   ,0.0046269894_DP,0.0038279891_DP,0.0030260086_DP,&
         &0.0022199750_DP,0.0014140010_DP,0.000533_DP    ,0.000075_DP    /)

    !     -----------------------------------------------------------------
    RETURN
  END SUBROUTINE rad_lon_SURRTFTR
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE rad_lon_surrtbg2

    !=========================================================================
    !
    !- Description:
    !
    !   Calculate lookup tables in mo_rrtbg2 for functions needed in routine 
    !   taumol
    !
    !   M.A. Giorgetta, MPI, June 2000
    !
    !=========================================================================

    IMPLICIT NONE
    INTRINSIC :: REAL, SQRT

    INTEGER :: i
    REAL(dp):: fp, rtfp

    corr1(0)   = 1._dp
    corr1(200) = 1._dp
    corr2(0)   = 1._dp
    corr2(200) = 1._dp

    DO i = 1,199
       fp       = 0.005_dp*REAL(i)
       rtfp     = SQRT(fp)
       corr1(i) = rtfp/fp
       corr2(i) = (1._dp-rtfp)/(1._dp-fp)
    ENDDO

!!CDIR DU_UPDATE(CORR1)
!!CDIR DU_UPDATE(CORR2)

  END SUBROUTINE rad_lon_SURRTBG2
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE rad_lon_surrta

    !=========================================================================
    !
    ! Description:
    !
    !   Reads data for mo_rrtaN (N=1:16) from formatted fortran file rrtadata
    !
    !   M.A. Giorgetta, MPI, June 2000
    !   L. Kornblueh, MPI, November 2001, fix read in parallel case
    !
    !=========================================================================

    IMPLICIT NONE

    CHARACTER(23), PARAMETER :: readform ='(/,  /,(tr1,5(es16.9)))'
    INTEGER                  :: i

    ! read data from formatted file 'rrtadata'
    ! ----------------------------------------

!  IF (p_parallel_io) THEN
    OPEN(11,file='rrtadata',status='old',form='formatted',action='read')

    DO i=1,12
       READ(11,*)
    END DO
!  END IF

    CALL rad_lon_read_rrta1
    CALL rad_lon_read_rrta2
    CALL rad_lon_read_rrta3
    CALL rad_lon_read_rrta4
    CALL rad_lon_read_rrta5
    CALL rad_lon_read_rrta6
    CALL rad_lon_read_rrta7
    CALL rad_lon_read_rrta8
    CALL rad_lon_read_rrta9
    CALL rad_lon_read_rrta10
    CALL rad_lon_read_rrta11
    CALL rad_lon_read_rrta12
    CALL rad_lon_read_rrta13
    CALL rad_lon_read_rrta14
    CALL rad_lon_read_rrta15
    CALL rad_lon_read_rrta16

    CLOSE(11)

  CONTAINS

    SUBROUTINE rad_lon_read_rrta1

      READ(11,'(//)')
      READ(11,readform) absa_1(:,:)
      READ(11,readform) absb_1(:,:)
      READ(11,readform) fracrefa_1(:)
      READ(11,readform) fracrefb_1(:)
      READ(11,readform) forref_1(:)
      READ(11,readform) selfref_1(:,:)

!!CDIR DU_UPDATE(ABSA_1)
!!CDIR DU_UPDATE(ABSB_1)
!!CDIR DU_UPDATE(SELFREF_1)

    END SUBROUTINE rad_lon_read_rrta1

    SUBROUTINE rad_lon_read_rrta2

      READ(11,'(//)')
      READ(11,readform) absa_2(:,:)
      READ(11,readform) absb_2(:,:)
      READ(11,readform) fracrefa_2(:,:)
      READ(11,readform) fracrefb_2(:)
      READ(11,readform) forref_2(:)
      READ(11,readform) selfref_2(:,:)
      READ(11,readform) refparam_2(:)
!!CDIR DU_UPDATE(ABSA_2)
!!CDIR DU_UPDATE(ABSB_2)
!!CDIR DU_UPDATE(SELFREF_2)
!!CDIR DU_UPDATE(FRACREFA_2)

    END SUBROUTINE rad_lon_read_rrta2

    SUBROUTINE rad_lon_read_rrta3

      READ(11,'(//)')
      READ(11,readform) absa_3(:,:)
      READ(11,readform) absb_3(:,:)
      READ(11,readform) fracrefa_3(:,:)
      READ(11,readform) fracrefb_3(:,:)
      READ(11,readform) forref_3(:)
      READ(11,readform) selfref_3(:,:)
      READ(11,readform) absn2oa_3(:)
      READ(11,readform) absn2ob_3(:)
      READ(11,readform) etaref_3(:)
      READ(11,readform) h2oref_3(:)
      READ(11,readform) n2oref_3(:)
      READ(11,readform) co2ref_3(:)
      READ(11,readform) strrat_3
!!CDIR DU_UPDATE(ABSA_3)
!!CDIR DU_UPDATE(ABSB_3)
!!CDIR DU_UPDATE(SELFREF_3)
!!CDIR DU_UPDATE(FRACREFA_3)
!!CDIR DU_UPDATE(FRACREFB_3)
!!CDIR DU_UPDATE(H2OREF_3)
!!CDIR DU_UPDATE(N2OREF_3)
!!CDIR DU_UPDATE(CO2REF_3)
!!CDIR DU_UPDATE(ETAREF_3)
    END SUBROUTINE rad_lon_read_rrta3

    SUBROUTINE rad_lon_read_rrta4

      READ(11,'(//)')
      READ(11,readform) absa_4(:,:)
      READ(11,readform) absb_4(:,:)
      READ(11,readform) fracrefa_4(:,:)
      READ(11,readform) fracrefb_4(:,:)
      READ(11,readform) selfref_4(:,:)
      READ(11,readform) strrat1_4
      READ(11,readform) strrat2_4
!!CDIR DU_UPDATE(ABSA_4)
!!CDIR DU_UPDATE(ABSB_4)
!!CDIR DU_UPDATE(SELFREF_4)
!!CDIR DU_UPDATE(FRACREFA_4)
!!CDIR DU_UPDATE(FRACREFB_4)

    END SUBROUTINE rad_lon_read_rrta4

    SUBROUTINE rad_lon_read_rrta5

      READ(11,'(//)')
      READ(11,readform) absa_5(:,:)
      READ(11,readform) absb_5(:,:)
      READ(11,readform) ccl4_5(:)
      READ(11,readform) fracrefa_5(:,:)
      READ(11,readform) fracrefb_5(:,:)
      READ(11,readform) selfref_5(:,:)
      READ(11,readform) strrat1_5
      READ(11,readform) strrat2_5
!!CDIR DU_UPDATE(ABSA_5)
!!CDIR DU_UPDATE(ABSB_5)
!!CDIR DU_UPDATE(SELFREF_5)
!!CDIR DU_UPDATE(FRACREFA_5)
!!CDIR DU_UPDATE(FRACREFB_5)

    END SUBROUTINE rad_lon_read_rrta5

    SUBROUTINE rad_lon_read_rrta6

      READ(11,'(//)')
      READ(11,readform) absa_6(:,:)
      READ(11,readform) absco2_6(:)
      READ(11,readform) cfc11adj_6(:)
      READ(11,readform) cfc12_6(:)
      READ(11,readform) fracrefa_6(:)
      READ(11,readform) selfref_6(:,:)

!!CDIR DU_UPDATE(ABSA_6)
!!CDIR DU_UPDATE(SELFREF_6)

    END SUBROUTINE rad_lon_read_rrta6

    SUBROUTINE rad_lon_read_rrta7

      READ(11,'(//)')
      READ(11,readform) absa_7(:,:)
      READ(11,readform) absb_7(:,:)
      READ(11,readform) absco2_7(:)
      READ(11,readform) fracrefa_7(:,:)
      READ(11,readform) fracrefb_7(:)
      READ(11,readform) selfref_7(:,:)
      READ(11,readform) strrat_7
!!CDIR DU_UPDATE(ABSA_7)
!!CDIR DU_UPDATE(ABSB_7)
!!CDIR DU_UPDATE(SELFREF_7)
!!CDIR DU_UPDATE(FRACREFA_7)
!!CDIR DU_UPDATE(FRACREFB_7)

    END SUBROUTINE rad_lon_read_rrta7

    SUBROUTINE rad_lon_read_rrta8
 
      READ(11,'(//)')
      READ(11,readform) absa_8(:,:)
      READ(11,readform) absb_8(:,:)
      READ(11,readform) fracrefa_8(:)
      READ(11,readform) fracrefb_8(:)
      READ(11,readform) selfref_8(:,:)
      READ(11,readform) absco2a_8(:)
      READ(11,readform) absco2b_8(:)
      READ(11,readform) absn2oa_8(:)
      READ(11,readform) absn2ob_8(:)
      READ(11,readform) cfc12_8(:)
      READ(11,readform) cfc22adj_8(:)
      READ(11,readform) h2oref_8(:)
      READ(11,readform) n2oref_8(:)
      READ(11,readform) o3ref_8(:)

!!CDIR DU_UPDATE(ABSA_8)
!!CDIR DU_UPDATE(ABSB_8)
!!CDIR DU_UPDATE(SELFREF_8)
!!CDIR DU_UPDATE(H2OREF_8)
!!CDIR DU_UPDATE(N2OREF_8)
!!CDIR DU_UPDATE(O3REF_8)

    END SUBROUTINE rad_lon_read_rrta8

    SUBROUTINE rad_lon_read_rrta9

      READ(11,'(//)')
      READ(11,readform) absa_9(:,:)
      READ(11,readform) absb_9(:,:)
      READ(11,readform) fracrefa_9(:,:)
      READ(11,readform) fracrefb_9(:)
      READ(11,readform) selfref_9(:,:)
      READ(11,readform) absn2o_9(:)
      READ(11,readform) ch4ref_9(:)
      READ(11,readform) etaref_9(:)
      READ(11,readform) h2oref_9(:)
      READ(11,readform) n2oref_9(:)
      READ(11,readform) strrat_9
!!CDIR DU_UPDATE(ABSA_9)
!!CDIR DU_UPDATE(ABSB_9)
!!CDIR DU_UPDATE(SELFREF_9)
!!CDIR DU_UPDATE(ABSN2O_9)
!!CDIR DU_UPDATE(FRACREFA_9)
!!CDIR DU_UPDATE(N2OREF_9)
!!CDIR DU_UPDATE(H2OREF_9)
!!CDIR DU_UPDATE(CH4REF_9)
!!CDIR DU_UPDATE(ETAREF_9)

    END SUBROUTINE rad_lon_read_rrta9

    SUBROUTINE rad_lon_read_rrta10

      READ(11,'(//)')
      READ(11,readform) absa_10(:,:)
      READ(11,readform) absb_10(:,:)
      READ(11,readform) fracrefa_10(:)
      READ(11,readform) fracrefb_10(:)

!!CDIR DU_UPDATE(ABSA_10)
!!CDIR DU_UPDATE(ABSB_10)

    END SUBROUTINE rad_lon_read_rrta10

    SUBROUTINE rad_lon_read_rrta11

      READ(11,'(//)')
      READ(11,readform) absa_11(:,:)
      READ(11,readform) absb_11(:,:)
      READ(11,readform) fracrefa_11(:)
      READ(11,readform) fracrefb_11(:)
      READ(11,readform) selfref_11(:,:)

!!CDIR DU_UPDATE(ABSA_11)
!!CDIR DU_UPDATE(ABSB_11)
!!CDIR DU_UPDATE(SELFREF_11)

    END SUBROUTINE rad_lon_read_rrta11

    SUBROUTINE rad_lon_read_rrta12

      READ(11,'(//)')
      READ(11,readform) absa_12(:,:)
      READ(11,readform) fracrefa_12(:,:)
      READ(11,readform) selfref_12(:,:)
      READ(11,readform) strrat_12

!!CDIR DU_UPDATE(ABSA_12)
!!CDIR DU_UPDATE(SELFREF_12)
!!CDIR DU_UPDATE(FRACREFA_12)

    END SUBROUTINE rad_lon_read_rrta12

    SUBROUTINE rad_lon_read_rrta13

      READ(11,'(//)')
      READ(11,readform) absa_13(:,:)
      READ(11,readform) fracrefa_13(:,:)
      READ(11,readform) selfref_13(:,:)
      READ(11,readform) strrat_13

!!CDIR DU_UPDATE(ABSA_13)
!!CDIR DU_UPDATE(SELFREF_13)
!!CDIR DU_UPDATE(FRACREFA_13)

    END SUBROUTINE rad_lon_read_rrta13

    SUBROUTINE rad_lon_read_rrta14

      READ(11,'(//)')
      READ(11,readform) absa_14(:,:)
      READ(11,readform) absb_14(:,:)
      READ(11,readform) fracrefa_14(:)
      READ(11,readform) fracrefb_14(:)
      READ(11,readform) selfref_14(:,:)

!!CDIR DU_UPDATE(ABSA_14)
!!CDIR DU_UPDATE(ABSB_14)
!!CDIR DU_UPDATE(SELFREF_14)

    END SUBROUTINE rad_lon_read_rrta14

    SUBROUTINE rad_lon_read_rrta15

      READ(11,'(//)')
      READ(11,readform) absa_15(:,:)
      READ(11,readform) fracrefa_15(:,:)
      READ(11,readform) selfref_15(:,:)
      READ(11,readform) strrat_15

!!CDIR DU_UPDATE(ABSA_15)
!!CDIR DU_UPDATE(SELFREF_15)
!!CDIR DU_UPDATE(FRACREFA_15)

    END SUBROUTINE rad_lon_read_rrta15

    SUBROUTINE rad_lon_read_rrta16

      READ(11,'(//)')
      READ(11,readform) absa_16(:,:)
      READ(11,readform) fracrefa_16(:,:)
      READ(11,readform) selfref_16(:,:)
      READ(11,readform) strrat_16

!!CDIR DU_UPDATE(ABSA_16)
!!CDIR DU_UPDATE(SELFREF_16)
!!CDIR DU_UPDATE(FRACREFA_16)

    END SUBROUTINE rad_lon_read_rrta16

  END SUBROUTINE rad_lon_surrta
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE rad_lon_RRTM_RRTM_140GP &
       ( kproma, kbdim, klev,        & ! IN
       ppave, ptave, ptl,            &
       ptbound, psemiss,             &
       pcldfrac,kcldlyr,             &
       pcoldry, pwkl, pwx,           &
       ptaucld, ptauaerl,            &
       ptotuflux, ptotufluc,         & ! OUT
       ptotdflux, ptotdfluc,         &
       psemit )

    ! *** This program is the driver for RRTM, the AER rapid model.  
    !     For each atmosphere the user wishes to analyze, this routine
    !     a) sets the atmospheric profile 
    !     b) calls SETCOEF to calculate various quantities needed for 
    !        the radiative transfer algorithm
    !     c) calls RTRN to do the radiative transfer calculation for
    !        clear or cloudy sky
    !     d) writes out the upward, downward, and net flux for each
    !        level and the heating rate for each layer

    !------------------------------Arguments--------------------------------

    IMPLICIT NONE

    ! Input arguments
    !
    INTEGER, INTENT(in)                             :: kproma    ! number of longitudes
    INTEGER, INTENT(in)                             :: kbdim     ! first dimension of 2-d arrays
    INTEGER, INTENT(in)                             :: klev      ! number of levels
    REAL(dp),INTENT(in),DIMENSION(kbdim,klev)       :: ppave     ! full level pressure [mb]
    REAL(dp),INTENT(in),DIMENSION(kbdim,klev)       :: ptave     ! full level temperature [K]
    REAL(dp),INTENT(in),DIMENSION(kbdim,klev+1)     :: ptl       ! half level temperature [K]
    REAL(dp),INTENT(in),DIMENSION(kbdim)            :: ptbound   ! surface temperature [K]
    REAL(dp),INTENT(in),DIMENSION(kbdim,jpband)     :: psemiss   ! surface emissivity in each band []
    INTEGER, INTENT(in),DIMENSION(kbdim,klev)       :: kcldlyr   ! cloud indicator, 0:clear, 1: cloudy
    REAL(dp),INTENT(in),DIMENSION(kbdim,klev)       :: pcldfrac  ! layer cloud fraction with respect to
    !                                                              cloud fraction of total column []
    REAL(dp),INTENT(in),DIMENSION(kbdim,klev)       :: pcoldry   ! number of molecules/cm2 of dry air and
    !                                                              water vapor [#/cm2]
    REAL(dp),INTENT(in),DIMENSION(kbdim,jpinpx,klev):: pwkl      ! number of molecules/cm2 of N species
    !                                                              in [#/cm2], N=JPINPX
    !                                                              1: H2O
    !                                                              2: CO2
    !                                                              3: O3
    !                                                              4: N2O
    !                                                              5: ------ empty ------
    !                                                              6: CH4
    !                                                              7... : -- empty ------
    REAL(dp),INTENT(in),DIMENSION(kbdim,jpxsec,klev):: pwx       ! number of molecules/cm2 of N species
    !                                                              in [1e20/cm2], N=JPXSEC
    !                                                              1: ------ empty ------
    !                                                              2: CFC11
    !                                                              3: CFC12
    !                                                              4... : -- empty ------
    REAL(dp),INTENT(in),DIMENSION(kbdim,klev,jpband):: ptaucld   ! optical thickness of clouds
    !                                                              in each band []
    REAL(dp),INTENT(in),DIMENSION(kbdim,klev,jpband):: ptauaerl  ! optical thickness of aerosols
    !                                                              in each band []
    !
    ! Output arguments
    !
    REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1)    :: ptotuflux ! upward flux, total sky
    REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1)    :: ptotufluc ! upward flux, clear sky
    REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1)    :: ptotdflux ! downward flux, total sky
    REAL(dp),INTENT(out),DIMENSION(kbdim,klev+1)    :: ptotdfluc ! downward flux, clear sky
    REAL(dp),INTENT(out),DIMENSION(kbdim)           :: psemit    ! surface emissivity

    !------------------------------RRTM variables---------------------------
    !
    ! These local variables have been defines in older versions in modules mo_rrtXYZ.
    ! Modules mo_rrtXYZ have been composed from fortran77 common blocks.

    !- from mo_rrtatm, COMMON INTFAC
    REAL(dp),DIMENSION(kbdim,klev)        :: FAC00
    REAL(dp),DIMENSION(kbdim,klev)        :: FAC01
    REAL(dp),DIMENSION(kbdim,klev)        :: FAC10
    REAL(dp),DIMENSION(kbdim,klev)        :: FAC11
    REAL(dp),DIMENSION(kbdim,klev)        :: FORFAC
    !
    !- from mo_rrtatm, COMMON INTIND
    INTEGER, DIMENSION(kbdim,klev)        :: JP
    INTEGER, DIMENSION(kbdim,klev)        :: JT
    INTEGER, DIMENSION(kbdim,klev)        :: JT1
    !
    !- from mo_rrtatm, COMMON PRECISE
    ! ONEMINUS is a real number just smaller than 1
    ! used in rrtm_taugbN some subroutines
    REAL(dp),PARAMETER                   :: ONEMINUS =1._dp-1.E-06_dp

    !
    !- from mo_rrtatm, COMMON PROFDATA             
    REAL(dp),DIMENSION(kbdim,klev)        :: COLH2O
    REAL(dp),DIMENSION(kbdim,klev)        :: COLCO2
    REAL(dp),DIMENSION(kbdim,klev)        :: COLO3
    REAL(dp),DIMENSION(kbdim,klev)        :: COLN2O
    REAL(dp),DIMENSION(kbdim,klev)        :: COLCH4
    REAL(dp),DIMENSION(kbdim,klev)        :: CO2MULT
    INTEGER, DIMENSION(kbdim)             :: LAYTROP
    INTEGER, DIMENSION(kbdim)             :: LAYSWTCH
    INTEGER, DIMENSION(kbdim)             :: LAYLOW
    !
    !- from mo_rrtatm, COMMON SELF             
    REAL(dp),DIMENSION(kbdim,klev)        :: SELFFAC
    REAL(dp),DIMENSION(kbdim,klev)        :: SELFFRAC
    INTEGER, DIMENSION(kbdim,klev)        :: INDSELF
    !
    !- from mo_rrtatm, COMMON SP             
    REAL(dp),DIMENSION(kbdim,jpgpt,klev)  :: PFRAC
    !
    ! former module mo_rrtctl
    ! -----------------------
    INTEGER, PARAMETER                   :: ISTART  = 1
    INTEGER, PARAMETER                   :: IEND    = JPBAND

    ! former module mo_rrtgd
    ! ----------------------
    ! - except for tau which is defined in rrtm_gasabs1a_140gp
    ! - equivalence is avoided
    REAL(dp), DIMENSION(kbdim,jpgpt*klev) :: ABSS1
    REAL(dp), DIMENSION(kbdim,jpgpt,klev) :: OD
    REAL(dp), DIMENSION(kbdim,jpgpt*klev) :: TAUSF1


    !  Calculate information needed by the radiative transfer routine
    !  that is specific to this atmosphere, especially some of the 
    !  coefficients and indices needed to compute the optical depths
    !  by interpolating data from stored reference atmospheres. 

    CALL rad_lon_RRTM_SETCOEF_140GP (KPROMA,KBDIM,KLEV,pcoldry,pwkl &
         &, FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1 &
         &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,CO2MULT &
         &, LAYTROP,LAYSWTCH,LAYLOW,ppave,ptave,SELFFAC,SELFFRAC,INDSELF)

!CDIR noiexpand
    CALL rad_lon_RRTM_GASABS1A_140GP (KPROMA,KBDIM,KLEV,ABSS1,OD,TAUSF1,pcoldry,pwx &
         &, ptauaerl,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS &
         &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,CO2MULT &
         &, LAYTROP,LAYSWTCH,LAYLOW,SELFFAC,SELFFRAC,INDSELF,PFRAC)

    !- Call the radiative transfer routine.

    !  Clear and cloudy parts of column are treated together in RTRN.
    !  Clear radiative transfer is done for clear layers and cloudy radiative
    !  transfer is done for cloudy layers as identified by icldlyr.

    ! mz_kk_20070305: compiler directive added
!CDIR noiexpand
    CALL rad_lon_RRTM_RTRN1A_140GP (KPROMA,KBDIM,KLEV,ISTART,IEND,kcldlyr,pcldfrac,ptaucld,ABSS1 &
         &, OD,TAUSF1,ptotdfluc,ptotdflux,ptotufluc,ptotuflux &
         &, ptave,ptl,ptbound,PFRAC,psemiss,psemit)

  END SUBROUTINE rad_lon_RRTM_RRTM_140GP
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE rad_lon_RRTM_SETCOEF_140GP (KPROMA,KBDIM,KLEV,COLDRY,WKL &
       &, FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1 &
       &, COLH2O,COLCO2,COLO3,COLN2O,COLCH4,CO2MULT &
       &, LAYTROP,LAYSWTCH,LAYLOW,PAVEL,TAVEL,SELFFAC,SELFFRAC,INDSELF)

    !     Reformatted for F90 by JJMorcrette, ECMWF, 980714
    !     Include longitude loop, jan2001, Marco A. Giorgetta, MPI

    !     Purpose:  For a given atmosphere, calculate the indices and
    !     fractions related to the pressure and temperature interpolations.
    !     Also calculate the values of the integrated Planck functions 
    !     for each band at the level and layer temperatures.

    IMPLICIT NONE

    !     DUMMY INTEGER SCALARS
    INTEGER :: KPROMA, KBDIM, KLEV

    REAL(DP):: COLDRY(KBDIM,KLEV)
    REAL(DP):: WKL(KBDIM,JPINPX,KLEV)

    !- from INTFAC      
    REAL(DP):: FAC00(KBDIM,KLEV)
    REAL(DP):: FAC01(KBDIM,KLEV)
    REAL(DP):: FAC10(KBDIM,KLEV)
    REAL(DP):: FAC11(KBDIM,KLEV)
    REAL(DP):: FORFAC(KBDIM,KLEV)

    !- from INTIND
    INTEGER :: JP(KBDIM,KLEV)
    INTEGER :: JT(KBDIM,KLEV)
    INTEGER :: JT1(KBDIM,KLEV)

    !- from PROFDATA             
    REAL(DP):: COLH2O(KBDIM,KLEV)
    REAL(DP):: COLCO2(KBDIM,KLEV)
    REAL(DP):: COLO3 (KBDIM,KLEV)
    REAL(DP):: COLN2O(KBDIM,KLEV)
    REAL(DP):: COLCH4(KBDIM,KLEV)
!!$REAL(DP):: COLO2 (KBDIM,KLEV)
    REAL(DP):: CO2MULT(KBDIM,KLEV)
    INTEGER :: LAYTROP(KBDIM)
    INTEGER :: LAYSWTCH(KBDIM)
    INTEGER :: LAYLOW(KBDIM)

    !- from PROFILE             
    REAL(DP):: PAVEL(KBDIM,KLEV)
    REAL(DP):: TAVEL(KBDIM,KLEV)

    !- from SELF             
    REAL(DP):: SELFFAC(KBDIM,KLEV)
    REAL(DP):: SELFFRAC(KBDIM,KLEV)
    INTEGER :: INDSELF(KBDIM,KLEV)


    ! The following source code follows that provided by Uwe Schulzweida,
    ! as given in rrtm_setcoef_140gp.f90@@/main/clean/2,
    ! except that COLO2 is not computed.
    ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    !     LOCAL INTEGER SCALARS
    INTEGER :: JP1, LAY, IPLON

    !     LOCAL REAL SCALARS
    REAL(dp) :: CO2REG, COMPFP, FACTOR, FP, FT, FT1, PLOG, SCALEFAC, STPFAC, WATER


    STPFAC = 296._dp/1013._dp

    LAYTROP(:)  = 0
    LAYSWTCH(:) = 0
    LAYLOW(:)   = 0


    DO LAY = 1, KLEV
       DO IPLON = 1, KPROMA
          !        Find the two reference pressures on either side of the
          !        layer pressure.  Store them in JP and JP1.  Store in FP the
          !        fraction of the difference (in ln(pressure)) between these
          !        two values that the layer pressure lies.
          PLOG = LOG(PAVEL(IPLON,LAY))
          JP(IPLON,LAY) = INT(36._dp - 5._DP*(PLOG+0.04_dp))
          IF (JP(IPLON,LAY)  <  1) THEN
             JP(IPLON,LAY) = 1
          ELSE IF (JP(IPLON,LAY)  >  58) THEN
             JP(IPLON,LAY) = 58
          ENDIF
          JP1 = JP(IPLON,LAY) + 1
          FP = 5._dp * (PREFLOG(JP(IPLON,LAY)) - PLOG)

          !        Determine, for each reference pressure (JP and JP1), which
          !        reference temperature (these are different for each  
          !        reference pressure) is nearest the layer temperature but does
          !        not exceed it.  Store these indices in JT and JT1, resp.
          !        Store in FT (resp. FT1) the fraction of the way between JT
          !        (JT1) and the next highest reference temperature that the 
          !        layer temperature falls.
          JT(IPLON,LAY) = INT(3._dp + (TAVEL(IPLON,LAY)-TREF(JP(IPLON,LAY)))/15._dp)
          IF (JT(IPLON,LAY)  <  1) THEN
             JT(IPLON,LAY) = 1
          ELSE IF (JT(IPLON,LAY)  >  4) THEN
             JT(IPLON,LAY) = 4
          ENDIF
          FT = ((TAVEL(IPLON,LAY)-TREF(JP(IPLON,LAY)))/15._dp) - REAL(JT(IPLON,LAY)-3, DP)
          JT1(IPLON,LAY) = INT(3._dp + (TAVEL(IPLON,LAY)-TREF(JP1))/15._dp)
          IF (JT1(IPLON,LAY)  <  1) THEN
             JT1(IPLON,LAY) = 1
          ELSE IF (JT1(IPLON,LAY)  >  4) THEN
             JT1(IPLON,LAY) = 4
          ENDIF
          FT1 = ((TAVEL(IPLON,LAY)-TREF(JP1))/15._dp) - REAL(JT1(IPLON,LAY)-3, DP)

          WATER = WKL(IPLON,1,LAY)/COLDRY(IPLON,LAY)
          SCALEFAC = PAVEL(IPLON,LAY) * STPFAC / TAVEL(IPLON,LAY)

          !        If the pressure is less than ~100mb, perform a different
          !        set of species interpolations.
          !         IF (PLOG .LE. 4.56) GO TO 5300
          !--------------------------------------         
          IF (PLOG  >  4.56_dp) THEN
             LAYTROP(IPLON) =  LAYTROP(IPLON) + 1
             !        For one band, the "switch" occurs at ~300 mb. 
             IF (PLOG  >=  5.76_dp) LAYSWTCH(IPLON) = LAYSWTCH(IPLON) + 1
             IF (PLOG  >=  6.62_dp) LAYLOW(IPLON) = LAYLOW(IPLON) + 1

             FORFAC(IPLON,LAY) = SCALEFAC / (1.0_dp+WATER)

             !        Set up factors needed to separately include the water vapor
             !        self-continuum in the calculation of absorption coefficient.
             !C           SELFFAC(IPLON,LAY) = WATER * SCALEFAC / (1.+WATER)
             SELFFAC(IPLON,LAY) = WATER * FORFAC(IPLON,LAY)
             FACTOR = (TAVEL(IPLON,LAY)-188.0_dp)/7.2_dp
             INDSELF(IPLON,LAY) = MIN(9, MAX(1, INT(FACTOR)-7))
             SELFFRAC(IPLON,LAY) = FACTOR - REAL(INDSELF(IPLON,LAY) + 7, DP)

             !        Calculate needed column amounts.
             COLH2O(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,1,LAY)
             COLCO2(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,2,LAY)
             COLO3(IPLON,LAY)  = 1.E-20_dp * WKL(IPLON,3,LAY)
             COLN2O(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,4,LAY)
             COLCH4(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,6,LAY)
!!$        COLO2(IPLON,LAY)  = 1.E-20_dp * WKL(IPLON,7,LAY)
             IF (COLCO2(IPLON,LAY)  ==  0.0_dp) COLCO2(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
             IF (COLN2O(IPLON,LAY)  ==  0.0_dp) COLN2O(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
             IF (COLCH4(IPLON,LAY)  ==  0.0_dp) COLCH4(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
             !        Using E = 1334.2 cm-1.
             CO2REG = 3.55E-24_dp * COLDRY(IPLON,LAY)
             CO2MULT(IPLON,LAY)= (COLCO2(IPLON,LAY) - CO2REG) *&
                  272.63_dp*EXP(-1919.4_dp/TAVEL(IPLON,LAY))/(8.7604E-4_dp*TAVEL(IPLON,LAY))
             !         GO TO 5400
             !------------------
          ELSE
             !        Above LAYTROP.
             ! 5300    CONTINUE

             !        Calculate needed column amounts.
             FORFAC(IPLON,LAY) = SCALEFAC / (1.0_dp+WATER)

             COLH2O(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,1,LAY)
             COLCO2(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,2,LAY)
             COLO3(IPLON,LAY)  = 1.E-20_dp * WKL(IPLON,3,LAY)
             COLN2O(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,4,LAY)
             COLCH4(IPLON,LAY) = 1.E-20_dp * WKL(IPLON,6,LAY)
!!$        COLO2(IPLON,LAY)  = 1.E-20_dp * WKL(IPLON,7,LAY)
             IF (COLCO2(IPLON,LAY)  ==  0.0_dp) COLCO2(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
             IF (COLN2O(IPLON,LAY)  ==  0.0_dp) COLN2O(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
             IF (COLCH4(IPLON,LAY)  ==  0.0_dp) COLCH4(IPLON,LAY) = 1.E-32_dp * COLDRY(IPLON,LAY)
             CO2REG = 3.55E-24_dp * COLDRY(IPLON,LAY)
             CO2MULT(IPLON,LAY)= (COLCO2(IPLON,LAY) - CO2REG) *&
                  272.63_dp*EXP(-1919.4_dp/TAVEL(IPLON,LAY))/(8.7604E-4_dp*TAVEL(IPLON,LAY))
             !----------------     
          ENDIF
          ! 5400    CONTINUE

          !        We have now isolated the layer ln pressure and temperature,
          !        between two reference pressures and two reference temperatures 
          !        (for each reference pressure).  We multiply the pressure 
          !        fraction FP with the appropriate temperature fractions to get 
          !        the factors that will be needed for the interpolation that yields
          !        the optical depths (performed in routines TAUGBn for band n).

          COMPFP = 1.0_dp - FP
          FAC10(IPLON,LAY) = COMPFP * FT
          FAC00(IPLON,LAY) = COMPFP * (1.0_dp - FT)
          FAC11(IPLON,LAY) = FP * FT1
          FAC01(IPLON,LAY) = FP * (1.0_dp - FT1)

       ENDDO
    ENDDO

    DO IPLON = 1, KPROMA
       ! MT 981104 
       !-- Set LAYLOW for profiles with surface pressure less than 750 hPa. 
       IF (LAYLOW(IPLON) == 0) LAYLOW(IPLON)=1
    ENDDO

    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! End of source code from Uwe Schulzweidas version
    ! rrtm_setcoef_140gp.f90@@/main/clean/2

    RETURN
  END SUBROUTINE rad_lon_RRTM_SETCOEF_140GP
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE rad_lon_RRTM_GASABS1A_140GP (KPROMA,KBDIM,KLEV,ABSS1,OD,TAUSF1,COLDRY,WX,&
       &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS,&
       &COLH2O,COLCO2,COLO3,COLN2O,COLCH4,CO2MULT,&
       &LAYTROP,LAYSWTCH,LAYLOW,SELFFAC,SELFFRAC,INDSELF,PFRAC)

    !     Reformatted for F90 by JJMorcrette, ECMWF, 980714
    !     Include longitude loop as in "RADOPT" of Uwe
    !         Schulzweida, Jan2001, Marco A. Giorgetta

    IMPLICIT NONE
    INTRINSIC :: INT, MAX

    ! Argument list
    ! ============
    INTEGER :: KPROMA
    INTEGER :: KBDIM
    INTEGER :: KLEV

    REAL(DP):: ABSS1(KBDIM,JPGPT*KLEV)
    REAL(DP):: TAUSF1(KBDIM,JPGPT*KLEV)

    REAL(DP):: OD(KBDIM,JPGPT,KLEV)
    REAL(DP):: COLDRY(KBDIM,KLEV)
    REAL(DP):: WX(KBDIM,JPXSEC,KLEV)

    !- from AER
    REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

    !- from INTFAC      
    REAL(DP):: FAC00(KBDIM,KLEV)
    REAL(DP):: FAC01(KBDIM,KLEV)
    REAL(DP):: FAC10(KBDIM,KLEV)
    REAL(DP):: FAC11(KBDIM,KLEV)
    REAL(DP):: FORFAC(KBDIM,KLEV)

    !- from INTIND
    INTEGER :: JP(KBDIM,KLEV)
    INTEGER :: JT(KBDIM,KLEV)
    INTEGER :: JT1(KBDIM,KLEV)

    !- from PRECISE             
    REAL(DP):: ONEMINUS

    !- from PROFDATA             
    REAL(DP):: COLH2O(KBDIM,KLEV)
    REAL(DP):: COLCO2(KBDIM,KLEV)
    REAL(DP):: COLO3 (KBDIM,KLEV)
    REAL(DP):: COLN2O(KBDIM,KLEV)
    REAL(DP):: COLCH4(KBDIM,KLEV)
    !REAL(DP):: COLO2 (KBDIM,KLEV)
    REAL(DP):: CO2MULT(KBDIM,KLEV)
    INTEGER :: LAYTROP(KBDIM)
    INTEGER :: LAYSWTCH(KBDIM)
    INTEGER :: LAYLOW(KBDIM)

    !- from SELF             
    REAL(DP):: SELFFAC(KBDIM,KLEV)
    REAL(DP):: SELFFRAC(KBDIM,KLEV)
    INTEGER :: INDSELF(KBDIM,KLEV)

    !- from SP             
    REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


    ! Local variables
    ! ===============

    ! as in ECMWF-CY23R1
    ! ------------------
    REAL(DP):: TAU   (KBDIM,JPGPT,KLEV)

    !     LOCAL INTEGER SCALARS
    INTEGER :: IPR, ITR, LAY, IR

    !     LOCAL REAL SCALARS
    REAL(DP):: ODEPTH, SECANG, TF


    ! additional source code from version "RADOPT" by Uwe Schulzweida
    ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    INTEGER :: IPLON
    INTEGER :: ICL, ICH

    !     LOCAL INTEGER ARRAYS
    INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)
    INTEGER :: IXS(KLEV), IXLOS(KBDIM,KLEV), IXHIGS(KBDIM,KLEV)

    REAL(DP), save, DIMENSION(0:5000) :: TRANS
!CDIR DUPLICATE(TRANS,256)

    logical, save             :: lfirst = .true.

    IXLOW(:,:) = 0
    IXLOS(:,:) = 0
    IXHIGH(:,:) = 0
    IXHIGS(:,:) = 0

    TAU(:,:,:)=0.0_dp

    if(lfirst)   then
       call rad_lon_SURRTAB
       lfirst = .false.
    end if

    DO LAY = 1, KLEV
       ICL = 0
       ICH = 0
       DO IPLON = 1, KPROMA
          IF ( LAY <= LAYTROP(IPLON) ) THEN
             ICL = ICL + 1
             IXLOW(ICL,LAY) = IPLON
          ELSE
             ICH = ICH + 1
             IXHIGH(ICH,LAY) = IPLON
          ENDIF
       ENDDO
       IXC(LAY) = ICL
    ENDDO

    DO LAY = 1, KLEV
       ICL = 0
       ICH = 0
       DO IPLON = 1, KPROMA
          IF ( LAY <= LAYSWTCH(IPLON) ) THEN
             ICL = ICL + 1
             IXLOS(ICL,LAY) = IPLON
          ELSE
             ICH = ICH + 1
             IXHIGS(ICH,LAY) = IPLON
          ENDIF
       ENDDO
       IXS(LAY) = ICL
    ENDDO
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! end of additional source code from "RADOPT"
    ! IXC,IXS,IXLOW,IXLOS,IXHIGH,IXHIGS are used
    ! in calls to RRTM_TAUMOLn, n=1:16, see below


    !- SECANG is equal to the secant of the diffusivity angle.
    SECANG = 1.66_DP

    ! For reasons which have to be investigated reduces inlining of the
    ! following subroutines the performance significantly.

!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL1  (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,&
         &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL2  (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,COLDRY,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,&
         &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL3  (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLCO2,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL4  (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLCO2,COLO3,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL5  (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,WX,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLCO2,COLO3,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL6  (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,WX,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
         &COLH2O,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL7  (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLO3,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL8  (KPROMA,KBDIM,KLEV,IXS,IXLOS,IXHIGS,TAU,WX,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
         &COLH2O,COLO3,COLN2O,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL9  (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLN2O,COLCH4,LAYSWTCH,LAYLOW,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL10 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
         &COLH2O,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL11 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
         &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL12 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLCO2,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL13 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL14 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
         &COLCO2,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL15 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLCO2,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)
!CDIR NOIEXPAND
    CALL rad_lon_RRTM_TAUMOL16 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
         &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
         &COLH2O,COLCH4,SELFFAC,SELFFRAC,INDSELF,PFRAC)

    !- Loop over g-channels.
    DO LAY = 1, KLEV
       DO IPR = 1, JPGPT
          DO IPLON = 1, KPROMA
             ODEPTH = SECANG * TAU(IPLON,IPR,LAY)
             ! bug-fix according to Malte Heinemann, ZMAW; (thanks to Michael Ponater)
             ODEPTH = MAX(ODEPTH, 0.0_dp)
             OD(IPLON,IPR,LAY) = ODEPTH
             TF = ODEPTH/(BPADE+ODEPTH)
             IF (ODEPTH <= 0._DP) TF = 0._DP
             ITR=INT(5.E+03_DP*TF+0.5_DP)
             IR=(LAY-1)*JPGPT+IPR
             ABSS1 (IPLON,IR) = 1._DP - TRANS(ITR)
             TAUSF1(IPLON,IR) = TF
          ENDDO
       ENDDO
    ENDDO

    RETURN

  CONTAINS

    SUBROUTINE rad_lon_SURRTAB

      !     -----------------------------------------------------------------
      !        * E.C.M.W.F. PHYSICS PACKAGE ** AER'S RRTM LW RADIATION **
      !     J.-J. MORCRETTE       E.C.M.W.F.      98/07/14
      !     -----------------------------------------------------------------
      !     mag, MPI, 25 February 2000: comment added

      IMPLICIT NONE
      INTRINSIC :: EXP, REAL

      !     LOCAL INTEGER SCALARS
      INTEGER :: ITR

      !     LOCAL REAL SCALARS
      REAL(DP) :: ZTAU, ZTFN

      TRANS(0)   =1._DP
      TRANS(5000)=0._DP
      DO ITR=1,4999
         ZTFN=REAL(ITR)/5000._DP
         ZTAU=BPADE*ZTFN/(1._DP-ZTFN)
         TRANS(ITR)=EXP(-ZTAU)
      ENDDO

!CDIR DU_UPDATE(TRANS)

      RETURN
    END SUBROUTINE rad_lon_SURRTAB

  END SUBROUTINE rad_lon_RRTM_GASABS1A_140GP
  ! ===========================================================================

  ! ===========================================================================
  SUBROUTINE rad_lon_RRTM_RTRN1A_140GP (KPROMA,KBDIM,KLEV,ISTART,IEND,ICLDLYR &
       &, CLDFRAC,TAUCLD,ABSS1 &
       &, OD,TAUSF1,TOTDFLUC,TOTDFLUX,TOTUFLUC,TOTUFLUX &
       &, TAVEL,TZ,TBOUND,PFRAC,SEMISS,SEMISLW)

    !     Reformatted for F90 by JJMorcrette, ECMWF, 980714
    !     Speed-up by D.Salmond, ECMWF, 9907
    !     Bug-fix by M.J. Iacono, AER, Inc., 9911
    !     Bug-fix by JJMorcrette, ECMWF, 991209 (RAT1, RAT2 initialization)
    !     Speed-up by D. Salmond, ECMWF, 9912
    !     Bug-fix by JJMorcrette, ECMWF, 0005 (extrapolation T<160K)
    !     Speed-up by D. Salmond, ECMWF, 000515
    !     Speed-up by Uwe Schulzweida, MPIFMET, 2001/01
    !     Change data i/o from modules to argument lists,
    !        and include new speed-up of Uwe Schulzweida
    !                     March 2001, Marco A. Giorgetta
    !     U. Schulzweida, MPI, May 2002, blocking (nproma)

    !-* This program calculates the upward fluxes, downward fluxes,
    !   and heating rates for an arbitrary atmosphere.  The input to
    !   this program is the atmospheric profile and all Planck function
    !   information.  First-order "numerical" quadrature is used for the 
    !   angle integration, i.e. only one exponential is computed per layer
    !   per g-value per band.  Cloud overlap is treated with a generalized
    !   maximum/random method in which adjacent cloud layers are treated
    !   with maximum overlap, and non-adjacent cloud groups are treated
    !   with random overlap.  For adjacent cloud layers, cloud information
    !   is carried from the previous two layers.

    IMPLICIT NONE
    INTRINSIC :: EXP, INT, MAX, MIN

    ! Argument list
    ! =============

    ! Input
    ! -----
    INTEGER :: KPROMA                    ! Number of columns
    INTEGER :: KBDIM                     ! first dimension of 2-d arrays
    INTEGER :: KLEV                      ! Number of layers
    INTEGER :: ISTART                    ! Index of 1st band for Planck emission
    INTEGER :: IEND                      ! Index of last band for Planck emission
    INTEGER :: ICLDLYR(KBDIM,KLEV)       ! Cloud indicator
    REAL(DP):: CLDFRAC(KBDIM,KLEV)       ! Cloud fraction
    REAL(DP):: TAUCLD(KBDIM,KLEV,JPBAND) ! Spectral cloud optical thickness
    REAL(DP):: ABSS1 (KBDIM,JPGPT*KLEV)  ! 
    REAL(DP):: OD    (KBDIM,JPGPT,KLEV)  ! Clear-sky optical thickness
    REAL(DP):: TAUSF1(KBDIM,JPGPT*KLEV)  ! 
    REAL(DP):: TAVEL(KBDIM,KLEV)         ! Layer temperature
    REAL(DP):: TZ(KBDIM,0:KLEV)          ! Level temperature
    REAL(DP):: TBOUND(KBDIM)             ! Surface temperature
    REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)   ! Planck function fractions
    REAL(DP):: SEMISS(KBDIM,JPBAND)      ! Surface spectral emissivity
    REAL(DP):: SEMISLW(KBDIM)            ! Surface emissivity

    ! Output
    ! ------
    REAL(DP):: TOTDFLUC(KBDIM,0:KLEV)    ! Clear-sky downward longwave flux
    REAL(DP):: TOTDFLUX(KBDIM,0:KLEV)    ! Downward longwave flux
    REAL(DP):: TOTUFLUC(KBDIM,0:KLEV)    ! Clear-sky upward longwave flux
    REAL(DP):: TOTUFLUX(KBDIM,0:KLEV)    ! Upward longwave flux


    !--------------------------------------------------------------------------
    !
    ! Maximum/Random cloud overlap variables
    ! for upward radiaitve transfer
    !  FACCLR2  fraction of clear radiance from previous layer that needs to 
    !           be switched to cloudy stream
    !  FACCLR1  fraction of the radiance that had been switched in the previous
    !           layer from cloudy to clear that needs to be switched back to
    !           cloudy in the current layer
    !  FACCLD2  fraction of cloudy radiance from previous layer that needs to 
    !           be switched to clear stream
    !  FACCLD1  fraction of the radiance that had been switched in the previous
    !           layer from clear to cloudy that needs to be switched back to
    !           clear in the current layer
    ! for downward radiaitve transfer
    !  FACCLR2D fraction of clear radiance from previous layer that needs to 
    !           be switched to cloudy stream
    !  FACCLR1D fraction of the radiance that had been switched in the previous
    !           layer from cloudy to clear that needs to be switched back to
    !           cloudy in the current layer
    !  FACCLD2D fraction of cloudy radiance from previous layer that needs to 
    !           be switched to clear stream
    !  FACCLD1D fraction of the radiance that had been switched in the previous
    !           layer from clear to cloudy that needs to be switched back to
    !           clear in the current layer
    !
    !--------------------------------------------------------------------------


    ! The following source code follows that provided by Uwe Schulzweida,
    ! as given in rrtm_rtrn1a_140gp.f90@@/main/clean/2
    ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    INTEGER :: INDLAY(KPROMA,KLEV),INDLEV(KPROMA,0:KLEV)

    REAL(dp) :: BBU1(KPROMA,JPGPT*KLEV),BBUTOT1(KPROMA,JPGPT*KLEV)
    REAL(dp) :: TLAYFRAC(KPROMA,KLEV),TLEVFRAC(KPROMA,0:KLEV)
    REAL(dp) :: BGLEV(KPROMA,JPGPT)
    REAL(dp) :: PLVL(KPROMA,JPBAND,0:KLEV),PLAY(KPROMA,JPBAND,0:KLEV),WTNUM(3)
    REAL(dp) :: SEMIS(KPROMA,JPGPT),RADUEMIT(KPROMA,JPGPT)

    REAL(dp) :: RADCLRU1(KPROMA,JPGPT) ,RADCLRD1(KPROMA,JPGPT)
    REAL(dp) :: RADLU1(KPROMA,JPGPT)   ,RADLD1(KPROMA,JPGPT)
    REAL(dp) :: TRNCLD(KPROMA,JPBAND,KLEV)
    REAL(dp) :: ABSCLDNW
    REAL(dp) :: ATOT1(KPROMA,JPGPT*KLEV)

    REAL(dp) :: SURFEMIS(KPROMA,JPBAND),PLNKEMIT(KPROMA,JPBAND)

    !******
    REAL(dp) :: clrradu(KPROMA,jpgpt),cldradu(KPROMA,jpgpt),oldcld
    REAL(dp) :: oldclr,rad(KPROMA,jpgpt)
    REAL(dp) :: faccld1(KPROMA,klev+1),faccld2(KPROMA,klev+1)
    REAL(dp) :: facclr1(KPROMA,klev+1),facclr2(KPROMA,klev+1)
    REAL(dp) :: faccmb1(KPROMA,klev+1),faccmb2(KPROMA,klev+1)
    REAL(dp) :: faccld1d(KPROMA,0:klev),faccld2d(KPROMA,0:klev)
    REAL(dp) :: facclr1d(KPROMA,0:klev),facclr2d(KPROMA,0:klev)
    REAL(dp) :: faccmb1d(KPROMA,0:klev),faccmb2d(KPROMA,0:klev)
    REAL(dp) :: clrradd(KPROMA,jpgpt),cldradd(KPROMA,jpgpt)
    INTEGER :: istcld(KPROMA,klev+1),istcldd(KPROMA,0:klev)

    !     LOCAL INTEGER SCALARS
    INTEGER :: ICCLD, ICCLDL(KPROMA), IX      
    INTEGER :: IPLON
    INTEGER :: IBAND, ICLDDN(KPROMA), IENT, INDBOUND(KPROMA), INDEX, IPR, LAY, LEV, NBI

    !     LOCAL REAL SCALARS
    REAL(dp) :: BBD(KPROMA,jpgpt),DELBGDN(KPROMA,jpgpt)
    REAL(dp) :: DELBGUP(KPROMA,jpgpt),BGLAY(KPROMA,jpgpt)
    REAL(dp) :: BBDTOT, CLDSRC, DBDTLAY, DBDTLEV,                    &
         DRAD1(KPROMA), DRADCL1(KPROMA), FACTOT1,           &
         FMAX, FMIN, GASSRC, ODSM, PLANKBND, RADCLD,                      &
         RADD, RADMOD, RAT1(KPROMA), RAT2(KPROMA), SUMPL(KPROMA),               &
         SUMPLEM(KPROMA), TBNDFRAC(KPROMA), TRNS, TTOT, URAD1(KPROMA), URADCL1(KPROMA)
    REAL(dp) :: CLDRADD_L,CLRRADD_L,RAD_L
    !******

    WTNUM(1) = 0.5_dp
    WTNUM(2) = 0.0_dp
    WTNUM(3) = 0.0_dp

    DO IPLON = 1, KPROMA
       !-start JJM_000511
       IF (TBOUND(IPLON) < 339._dp .AND. TBOUND(IPLON) >= 160._dp ) THEN
          INDBOUND(IPLON) = TBOUND(IPLON) - 159._dp
          TBNDFRAC(IPLON) = TBOUND(IPLON) - INT(TBOUND(IPLON))
       ELSE IF (TBOUND(IPLON) >= 339._dp ) THEN
          INDBOUND(IPLON) = 180
          TBNDFRAC(IPLON) = TBOUND(IPLON) - 339._dp
       ELSE IF (TBOUND(IPLON) < 160._dp ) THEN
          INDBOUND(IPLON) = 1
          TBNDFRAC(IPLON) = TBOUND(IPLON) - 160._dp
       ENDIF
       !-end JJM_000511
    END DO

    TOTUFLUC(1:KPROMA,:) = 0.0_dp
    TOTDFLUC(1:KPROMA,:) = 0.0_dp
    TOTUFLUX(1:KPROMA,:) = 0.0_dp
    TOTDFLUX(1:KPROMA,:) = 0.0_dp

    TRNCLD(1:KPROMA,:,:) = 0.0_dp

    DO LAY = 0, KLEV
       DO IPLON = 1, KPROMA
          !-start JJM_000511
          IF (TZ(IPLON,LAY) < 339._dp .AND. TZ(IPLON,LAY) >= 160._dp ) THEN
             INDLEV(IPLON,LAY) = TZ(IPLON,LAY) - 159._dp
             TLEVFRAC(IPLON,LAY) = TZ(IPLON,LAY) - INT(TZ(IPLON,LAY))
          ELSE IF (TZ(IPLON,LAY) >= 339._dp ) THEN
             INDLEV(IPLON,LAY) = 180
             TLEVFRAC(IPLON,LAY) = TZ(IPLON,LAY) - 339._dp
          ELSE IF (TZ(IPLON,LAY) < 160._dp ) THEN
             INDLEV(IPLON,LAY) = 1
             TLEVFRAC(IPLON,LAY) = TZ(IPLON,LAY) - 160._dp
          ENDIF
          !-end JJM_000511
       ENDDO
    ENDDO

    !- jjm_991209
    FACCLD1(1:KPROMA,:) = 0.0_dp
    FACCLD2(1:KPROMA,:) = 0.0_dp
    FACCLR1(1:KPROMA,:) = 0.0_dp
    FACCLR2(1:KPROMA,:) = 0.0_dp
    FACCMB1(1:KPROMA,:) = 0.0_dp
    FACCMB2(1:KPROMA,:) = 0.0_dp
    FACCLD1D(1:KPROMA,:)  = 0.0_dp
    FACCLD2D(1:KPROMA,:)  = 0.0_dp
    FACCLR1D(1:KPROMA,:)  = 0.0_dp
    FACCLR2D(1:KPROMA,:)  = 0.0_dp
    FACCMB1D(1:KPROMA,:)  = 0.0_dp
    FACCMB2D(1:KPROMA,:)  = 0.0_dp

    RAT1(1:KPROMA)    = 0.0_dp
    RAT2(1:KPROMA)    = 0.0_dp

    !- jjm_991209

    SUMPL(1:KPROMA)   = 0.0_dp
    SUMPLEM(1:KPROMA) = 0.0_dp

    ISTCLD(1:KPROMA,1)     = 1
    ISTCLDD(1:KPROMA,KLEV) = 1

    DO LEV = 1, KLEV
       DO IPLON = 1, KPROMA
          !-- DS_000515
          !-start JJM_000511
          IF (TAVEL(IPLON,LEV) < 339._dp .AND. TAVEL(IPLON,LEV) >= 160._dp ) THEN
             INDLAY(IPLON,LEV) = TAVEL(IPLON,LEV) - 159._dp
             TLAYFRAC(IPLON,LEV) = TAVEL(IPLON,LEV) - INT(TAVEL(IPLON,LEV))
          ELSE IF (TAVEL(IPLON,LEV) >= 339._dp ) THEN
             INDLAY(IPLON,LEV) = 180
             TLAYFRAC(IPLON,LEV) = TAVEL(IPLON,LEV) - 339._dp
          ELSE IF (TAVEL(IPLON,LEV) < 160._dp ) THEN
             INDLAY(IPLON,LEV) = 1
             TLAYFRAC(IPLON,LEV) = TAVEL(IPLON,LEV) - 160._dp
          ENDIF
          !-end JJM_000511
          !-- DS_000515
       END DO
    END DO

    DO LEV = 1, KLEV
       DO IPLON = 1, KPROMA
          IF (ICLDLYR(IPLON,LEV) == 1) THEN

             !mji    
             ISTCLD(IPLON,LEV+1) = 0

             IF (LEV  ==  KLEV) THEN
                FACCLD1(IPLON,LEV+1) = 0.0_dp
                FACCLD2(IPLON,LEV+1) = 0.0_dp
                FACCLR1(IPLON,LEV+1) = 0.0_dp
                FACCLR2(IPLON,LEV+1) = 0.0_dp
                FACCMB1(IPLON,LEV+1) = 0.0_dp
                FACCMB2(IPLON,LEV+1) = 0.0_dp
                !mji      ISTCLD(IPLON,LEV+1) = 0

             ELSE IF (CLDFRAC(IPLON,LEV+1)  >=  CLDFRAC(IPLON,LEV)) THEN
                FACCLD1(IPLON,LEV+1) = 0.0_dp
                FACCLD2(IPLON,LEV+1) = 0.0_dp


                IF (ISTCLD(IPLON,LEV)  ==  1) THEN
                   !mji        ISTCLD(IPLON,LEV+1) = 0
                   FACCLR1(IPLON,LEV+1) = 0.0_dp
                   !mji        
                   FACCLR2(IPLON,LEV+1) = 0.0_dp
                   IF (CLDFRAC(IPLON,LEV) < 1.0_dp) THEN
                      FACCLR2(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV+1)-CLDFRAC(IPLON,LEV)) &
                           / (1.0_dp-CLDFRAC(IPLON,LEV))
                   END IF
                ELSE
                   FMAX = MAX(CLDFRAC(IPLON,LEV),CLDFRAC(IPLON,LEV-1))
                   !mji
                   IF (CLDFRAC(IPLON,LEV+1)  >  FMAX) THEN
                      FACCLR1(IPLON,LEV+1) = RAT2(IPLON)
                      FACCLR2(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV+1)-FMAX)/(1.0_dp-FMAX)
                      !mji          
                   ELSE IF (CLDFRAC(IPLON,LEV+1) < FMAX) THEN
                      FACCLR1(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV+1)-CLDFRAC(IPLON,LEV)) &
                           / (CLDFRAC(IPLON,LEV-1)-CLDFRAC(IPLON,LEV))
                      FACCLR2(IPLON,LEV+1) = 0.0_dp
                      !mji
                   ELSE
                      FACCLR1(IPLON,LEV+1) = RAT2(IPLON)
                      FACCLR2(IPLON,LEV+1) = 0.0_dp
                   ENDIF
                ENDIF
                IF (FACCLR1(IPLON,LEV+1) > 0.0_dp .OR. FACCLR2(IPLON,LEV+1) > 0.0_dp) THEN
                   RAT1(IPLON) = 1.0_dp
                   RAT2(IPLON) = 0.0_dp
                ENDIF
             ELSE
                FACCLR1(IPLON,LEV+1) = 0.0_dp
                FACCLR2(IPLON,LEV+1) = 0.0_dp
                IF (ISTCLD(IPLON,LEV)  ==  1) THEN
                   !mji        ISTCLD(IPLON,LEV+1) = 0
                   FACCLD1(IPLON,LEV+1) = 0.0_dp
                   FACCLD2(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV)-CLDFRAC(IPLON,LEV+1)) &
                        /  CLDFRAC(IPLON,LEV)
                ELSE
                   FMIN = MIN(CLDFRAC(IPLON,LEV),CLDFRAC(IPLON,LEV-1))
                   IF (CLDFRAC(IPLON,LEV+1)  <=  FMIN) THEN
                      FACCLD1(IPLON,LEV+1) = RAT1(IPLON)
                      FACCLD2(IPLON,LEV+1) = (FMIN-CLDFRAC(IPLON,LEV+1))/FMIN
                   ELSE
                      FACCLD1(IPLON,LEV+1) = (CLDFRAC(IPLON,LEV)-CLDFRAC(IPLON,LEV+1)) &
                           / (CLDFRAC(IPLON,LEV)-FMIN)
                      FACCLD2(IPLON,LEV+1) = 0.0_dp
                   ENDIF
                ENDIF
                IF (FACCLD1(IPLON,LEV+1) > 0.0_dp .OR. FACCLD2(IPLON,LEV+1) > 0.0_dp) THEN
                   RAT1(IPLON) = 0.0_dp
                   RAT2(IPLON) = 1.0_dp
                ENDIF
             ENDIF
             !fcc
             IF (LEV == 1) THEN
                FACCMB1(IPLON,LEV+1) = 0.
                FACCMB2(IPLON,LEV+1) = FACCLD1(IPLON,LEV+1) * FACCLR2(IPLON,LEV)
             ELSE
                FACCMB1(IPLON,LEV+1) = FACCLR1(IPLON,LEV+1) * FACCLD2(IPLON,LEV) * &
                     CLDFRAC(IPLON,LEV-1)
                FACCMB2(IPLON,LEV+1) = FACCLD1(IPLON,LEV+1) * FACCLR2(IPLON,LEV) * &
                     (1.0_dp - CLDFRAC(IPLON,LEV-1)) 
             ENDIF
             !end fcc
          ELSE
             !-- DS_000515
             ISTCLD(IPLON,LEV+1) = 1
          ENDIF
       ENDDO
    ENDDO

    !- jjm_991209
    RAT1(1:KPROMA)    = 0.0_dp
    RAT2(1:KPROMA)    = 0.0_dp
    !- jjm_991209

    DO LEV = KLEV, 1, -1
       DO IPLON = 1, KPROMA
          IF (ICLDLYR(IPLON,LEV) == 1) THEN
             !mji
             ISTCLDD(IPLON,LEV-1) = 0  
             IF (LEV  ==  1) THEN
                FACCLD1D(IPLON,LEV-1) = 0.0_dp
                FACCLD2D(IPLON,LEV-1) = 0.0_dp
                FACCLR1D(IPLON,LEV-1) = 0.0_dp
                FACCLR2D(IPLON,LEV-1) = 0.0_dp
                FACCMB1D(IPLON,LEV-1) = 0.0_dp
                FACCMB2D(IPLON,LEV-1) = 0.0_dp
                !mji      ISTCLDD(IPLON,LEV-1) = 0.0_dp
             ELSEIF (CLDFRAC(IPLON,LEV-1)  >=  CLDFRAC(IPLON,LEV)) THEN
                FACCLD1D(IPLON,LEV-1) = 0.0_dp
                FACCLD2D(IPLON,LEV-1) = 0.0_dp
                IF (ISTCLDD(IPLON,LEV)  ==  1) THEN
                   !mji        ISTCLDD(IPLON,LEV-1) = 0
                   FACCLR1D(IPLON,LEV-1) = 0.0_dp
                   FACCLR2D(IPLON,LEV-1) = 0.0_dp
                   IF (CLDFRAC(IPLON,LEV) < 1.0_dp) THEN
                      FACCLR2D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV-1)-CLDFRAC(IPLON,LEV))&
                           / (1.0_dp-CLDFRAC(IPLON,LEV))
                   END IF
                ELSE
                   FMAX = MAX(CLDFRAC(IPLON,LEV),CLDFRAC(IPLON,LEV+1))
                   !mji
                   IF (CLDFRAC(IPLON,LEV-1)  >  FMAX) THEN
                      FACCLR1D(IPLON,LEV-1) = RAT2(IPLON)
                      FACCLR2D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV-1)-FMAX)/(1.0_dp-FMAX)
                      !mji
                   ELSE IF (CLDFRAC(IPLON,LEV-1) < FMAX) THEN
                      FACCLR1D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV-1)-CLDFRAC(IPLON,LEV))&
                           / (CLDFRAC(IPLON,LEV+1)-CLDFRAC(IPLON,LEV))
                      FACCLR2D(IPLON,LEV-1) = 0.0_dp
                      !mji
                   ELSE          
                      FACCLR1D(IPLON,LEV-1) = RAT2(IPLON)
                      FACCLR2D(IPLON,LEV-1) = 0.0_dp
                   ENDIF
                ENDIF
                IF (FACCLR1D(IPLON,LEV-1) > 0.0_dp .OR. FACCLR2D(IPLON,LEV-1) > 0.0_dp)THEN
                   RAT1(IPLON) = 1.0_dp
                   RAT2(IPLON) = 0.0_dp
                ENDIF
             ELSE
                FACCLR1D(IPLON,LEV-1) = 0.0_dp
                FACCLR2D(IPLON,LEV-1) = 0.0_dp
                IF (ISTCLDD(IPLON,LEV)  ==  1) THEN
                   !mji        ISTCLDD(IPLON,LEV-1) = 0
                   FACCLD1D(IPLON,LEV-1) = 0.0_dp
                   FACCLD2D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV)-CLDFRAC(IPLON,LEV-1)) &
                        /  CLDFRAC(IPLON,LEV)
                ELSE
                   FMIN = MIN(CLDFRAC(IPLON,LEV),CLDFRAC(IPLON,LEV+1))
                   IF (CLDFRAC(IPLON,LEV-1)  <=  FMIN) THEN
                      FACCLD1D(IPLON,LEV-1) = RAT1(IPLON)
                      FACCLD2D(IPLON,LEV-1) = (FMIN-CLDFRAC(IPLON,LEV-1))/FMIN
                   ELSE
                      FACCLD1D(IPLON,LEV-1) = (CLDFRAC(IPLON,LEV)-CLDFRAC(IPLON,LEV-1))&
                           / (CLDFRAC(IPLON,LEV)-FMIN)
                      FACCLD2D(IPLON,LEV-1) = 0.0_dp
                   ENDIF
                ENDIF
                IF (FACCLD1D(IPLON,LEV-1) > 0.0_dp .OR. FACCLD2D(IPLON,LEV-1) > 0.0_dp)THEN
                   RAT1(IPLON) = 0.0_dp
                   RAT2(IPLON) = 1.0_dp
                ENDIF
             ENDIF
             !mag
             IF (LEV == KLEV) THEN
                FACCMB1D(IPLON,LEV-1) = 0.
                FACCMB2D(IPLON,LEV-1) = FACCLD1D(IPLON,LEV-1) * FACCLR2D(IPLON,LEV)
             ELSE
                FACCMB1D(IPLON,LEV-1) = FACCLR1D(IPLON,LEV-1) * FACCLD2D(IPLON,LEV) &
                     * CLDFRAC(IPLON,LEV+1)
                FACCMB2D(IPLON,LEV-1) = FACCLD1D(IPLON,LEV-1) * FACCLR2D(IPLON,LEV) &
                     * (1.0_dp - CLDFRAC(IPLON,LEV+1))
             ENDIF
             !end mag
          ELSE
             ISTCLDD(IPLON,LEV-1) = 1
          ENDIF
       ENDDO
    ENDDO

    !- Loop over frequency bands.

    DO IBAND = ISTART, IEND
       DO IPLON = 1, KPROMA
          DBDTLEV  = TOTPLNK(INDBOUND(IPLON)+1,IBAND)-TOTPLNK(INDBOUND(IPLON),IBAND)
          PLANKBND = DELWAVE(IBAND) &
               * (TOTPLNK(INDBOUND(IPLON),IBAND) + TBNDFRAC(IPLON) * DBDTLEV)
          DBDTLEV  = TOTPLNK(INDLEV(IPLON,0)+1,IBAND)-TOTPLNK(INDLEV(IPLON,0),IBAND)
          PLVL(IPLON,IBAND,0) = DELWAVE(IBAND) &
               * (TOTPLNK(INDLEV(IPLON,0),IBAND) + TLEVFRAC(IPLON,0)*DBDTLEV)

          SURFEMIS(IPLON,IBAND) = SEMISS(IPLON,IBAND)
          PLNKEMIT(IPLON,IBAND) = SURFEMIS(IPLON,IBAND) * PLANKBND
          SUMPLEM(IPLON)  = SUMPLEM(IPLON) + PLNKEMIT(IPLON,IBAND)
          SUMPL(IPLON)    = SUMPL(IPLON)   + PLANKBND
          !--DS
       ENDDO
    ENDDO

    !---
    DO IBAND = ISTART, IEND
       DO LEV = 1, KLEV
          DO IPLON = 1, KPROMA
             !----              
             !- Calculate the integrated Planck functions for at the
             !  level and layer temperatures.
             !  Compute cloud transmittance for cloudy layers.
             DBDTLEV = TOTPLNK(INDLEV(IPLON,LEV)+1,IBAND) &
                  - TOTPLNK(INDLEV(IPLON,LEV),IBAND)
             DBDTLAY = TOTPLNK(INDLAY(IPLON,LEV)+1,IBAND) &
                  - TOTPLNK(INDLAY(IPLON,LEV),IBAND)
             PLAY(IPLON,IBAND,LEV) = DELWAVE(IBAND) &
                  * (TOTPLNK(INDLAY(IPLON,LEV),IBAND)+TLAYFRAC(IPLON,LEV)*DBDTLAY)
             PLVL(IPLON,IBAND,LEV) = DELWAVE(IBAND) &
                  * (TOTPLNK(INDLEV(IPLON,LEV),IBAND)+TLEVFRAC(IPLON,LEV)*DBDTLEV)
             IF (ICLDLYR(IPLON,LEV) > 0) THEN
                TRNCLD(IPLON,IBAND,LEV) = EXP(-TAUCLD(IPLON,LEV,IBAND))
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    SEMISLW(1:KPROMA) = SUMPLEM(1:KPROMA) / SUMPL(1:KPROMA)

    !- Initialize for radiative transfer.
    DO IPR = 1, JPGPT
       NBI = NGB(IPR)
       DO IPLON = 1, KPROMA
          RADCLRD1(IPLON,IPR) = 0.0_dp
          RADLD1(IPLON,IPR)   = 0.0_dp
          SEMIS(IPLON,IPR)    = SURFEMIS(IPLON,NBI)
          RADUEMIT(IPLON,IPR) = PFRAC(IPLON,IPR,1) * PLNKEMIT(IPLON,NBI)
          BGLEV(IPLON,IPR)    = PFRAC(IPLON,IPR,KLEV) * PLVL(IPLON,NBI,KLEV)
       ENDDO
    ENDDO

    !- Downward radiative transfer.
    !  *** DRAD1 holds summed radiance for total sky stream
    !  *** DRADCL1 holds summed radiance for clear sky stream

    ICLDDN(1:KPROMA) = 0

    DO LEV = KLEV, 1, -1
       DRAD1(1:KPROMA)   = 0.0_dp
       DRADCL1(1:KPROMA) = 0.0_dp

       ICCLD = 0
       DO IPLON = 1, KPROMA
          IF (ICLDLYR(IPLON,LEV) == 1) THEN
             ICCLD = ICCLD + 1
             ICCLDL(ICCLD) = IPLON
          ENDIF
       ENDDO

       IENT = JPGPT * (LEV-1)
       DO IPR = 1, JPGPT
          NBI = NGB(IPR)
          INDEX = IENT + IPR

          DO IPLON = 1, KPROMA

             BGLAY(IPLON,IPR) = PFRAC(IPLON,IPR,LEV) * PLAY(IPLON,NBI,LEV)
             !----            
             DELBGUP(IPLON,IPR)     = BGLEV(IPLON,IPR) - BGLAY(IPLON,IPR)
             BBU1(IPLON,INDEX) = BGLAY(IPLON,IPR) + TAUSF1(IPLON,INDEX) * DELBGUP(IPLON,IPR)
             !--DS            
             BGLEV(IPLON,IPR) = PFRAC(IPLON,IPR,LEV) * PLVL(IPLON,NBI,LEV-1)
             !----            
             DELBGDN(IPLON,IPR) = BGLEV(IPLON,IPR) - BGLAY(IPLON,IPR)
             BBD(IPLON,IPR) = BGLAY(IPLON,IPR) + TAUSF1(IPLON,INDEX) * DELBGDN(IPLON,IPR)
          ENDDO
       ENDDO

       DO IPR = 1, JPGPT
          NBI = NGB(IPR)
          INDEX = IENT + IPR

!CDIR NODEP
!OCL NOVREC
          DO IX = 1, ICCLD
             IPLON = ICCLDL(IX)

             !  *** Cloudy layer
             ICLDDN(IPLON) = 1

             !- total-sky downward flux          
             ODSM = OD(IPLON,IPR,LEV) + TAUCLD(IPLON,LEV,NBI)
             FACTOT1 = ODSM / (BPADE + ODSM)
             BBUTOT1(IPLON,INDEX) = BGLAY(IPLON,IPR) + FACTOT1 * DELBGUP(IPLON,IPR)
             ABSCLDNW = 1.0_dp - TRNCLD(IPLON,NBI,LEV)
             ATOT1(IPLON,INDEX) = ABSS1(IPLON,INDEX) + ABSCLDNW &
                  - ABSS1(IPLON,INDEX) * ABSCLDNW
             BBDTOT = BGLAY(IPLON,IPR) + FACTOT1 * DELBGDN(IPLON,IPR)
             GASSRC = BBD(IPLON,IPR) * ABSS1(IPLON,INDEX)
             !***
             IF (ISTCLDD(IPLON,LEV)  ==  1) THEN
                CLDRADD_L = CLDFRAC(IPLON,LEV) * RADLD1(IPLON,IPR)
                CLRRADD_L = RADLD1(IPLON,IPR) - CLDRADD_L
                RAD_L = 0.0_dp
             ELSE
                CLDRADD_L = CLDRADD(IPLON,IPR)
                CLRRADD_L = CLRRADD(IPLON,IPR)
                RAD_L = RAD(IPLON,IPR)
             ENDIF
             TTOT = 1.0_dp - ATOT1(IPLON,INDEX)
             CLDSRC = BBDTOT * ATOT1(IPLON,INDEX)

             ! Separate RT equations for clear and cloudy streams      
             CLDRADD(IPLON,IPR) = CLDRADD_L * TTOT &
                  + CLDFRAC(IPLON,LEV) * CLDSRC
             CLRRADD(IPLON,IPR) = CLRRADD_L * &
                  (1.0_dp - ABSS1(IPLON,INDEX)) + &
                  (1.0_dp - CLDFRAC(IPLON,LEV)) * GASSRC

             !  Total sky downward radiance
             RADLD1(IPLON,IPR) = CLDRADD(IPLON,IPR) + CLRRADD(IPLON,IPR)
             DRAD1(IPLON) = DRAD1(IPLON) + RADLD1(IPLON,IPR)

             !- clear-sky downward flux          
             RADCLRD1(IPLON,IPR) = RADCLRD1(IPLON,IPR) &
                  + (BBD(IPLON,IPR)-RADCLRD1(IPLON,IPR))*ABSS1(IPLON,INDEX)
             DRADCL1(IPLON) = DRADCL1(IPLON) + RADCLRD1(IPLON,IPR)

             !* Code to account for maximum/random overlap:
             !   Performs RT on the radiance most recently switched between clear and
             !   cloudy streams
             RADMOD = RAD_L * (FACCLR1D(IPLON,LEV-1) * &
                  (1.0_dp-ABSS1(IPLON,INDEX)) +                 &
                  FACCLD1D(IPLON,LEV-1) *  TTOT) -              &
                  FACCMB1D(IPLON,LEV-1) * GASSRC +              &
                  FACCMB2D(IPLON,LEV-1) * CLDSRC

             !   Computes what the clear and cloudy streams would have been had no
             !   radiance been switched       
             OLDCLD = CLDRADD(IPLON,IPR) - RADMOD
             OLDCLR = CLRRADD(IPLON,IPR) + RADMOD

             !   Computes the radiance to be switched between clear and cloudy.
             RAD(IPLON,IPR) = -RADMOD + FACCLR2D(IPLON,LEV-1)*OLDCLR - &
                  FACCLD2D(IPLON,LEV-1)*OLDCLD
             CLDRADD(IPLON,IPR) = CLDRADD(IPLON,IPR) + RAD(IPLON,IPR)
             CLRRADD(IPLON,IPR) = CLRRADD(IPLON,IPR) - RAD(IPLON,IPR)

          ENDDO
       ENDDO

       DO IPR = 1, JPGPT
          NBI = NGB(IPR)
          INDEX = IENT + IPR


          DO IPLON=1,KPROMA
             IF (ICLDLYR(IPLON,LEV) .NE. 1) THEN
                !  *** Clear layer
                !  *** DRAD1 holds summed radiance for total sky stream
                !  *** DRADCL1 holds summed radiance for clear sky stream

                !- total-sky downward flux          
                RADLD1(IPLON,IPR) = RADLD1(IPLON,IPR) &
                     + (BBD(IPLON,IPR)-RADLD1(IPLON,IPR))*ABSS1(IPLON,INDEX)
                DRAD1(IPLON) = DRAD1(IPLON) + RADLD1(IPLON,IPR)
                !- clear-sky downward flux          
                !-  Set clear sky stream to total sky stream as long as layers
                !-  remain clear.  Streams diverge when a cloud is reached.
                IF (ICLDDN(IPLON) == 1) THEN
                   RADCLRD1(IPLON,IPR) = RADCLRD1(IPLON,IPR) &
                        + (BBD(IPLON,IPR)-RADCLRD1(IPLON,IPR))*ABSS1(IPLON,INDEX)
                   DRADCL1(IPLON) = DRADCL1(IPLON) + RADCLRD1(IPLON,IPR)
                ELSE
                   RADCLRD1(IPLON,IPR) = RADLD1(IPLON,IPR)
                   DRADCL1(IPLON) = DRAD1(IPLON)
                ENDIF
             ENDIF

          ENDDO
       ENDDO

       TOTDFLUC(1:KPROMA,LEV-1) = DRADCL1(1:KPROMA) * WTNUM(1)
       TOTDFLUX(1:KPROMA,LEV-1) = DRAD1(1:KPROMA)   * WTNUM(1)

    ENDDO

    ! Spectral reflectivity and reflectance
    ! Includes the contribution of spectrally varying longwave emissivity 
    ! and reflection from the surface to the upward radiative transfer.
    ! Note: Spectral and Lambertian reflections are identical for the one
    ! angle flux integration used here.

    URAD1(1:KPROMA)   = 0.0_dp
    URADCL1(1:KPROMA) = 0.0_dp

    !IF (IREFLECT  ==  0) THEN
    !- Lambertian reflection.
    DO IPR = 1, JPGPT
       DO IPLON = 1, KPROMA
          ! Clear-sky radiance
          !      RADCLD = 2.0_dp * (RADCLRD1(IPLON,IPR) * WTNUM(1) )
          RADCLD = RADCLRD1(IPLON,IPR)
          RADCLRU1(IPLON,IPR) = RADUEMIT(IPLON,IPR) &
               + (1.0_dp - SEMIS(IPLON,IPR)) * RADCLD
          URADCL1(IPLON) = URADCL1(IPLON) + RADCLRU1(IPLON,IPR)

          ! Total sky radiance
          !      RADD = 2.0_dp * (RADLD1(IPLON,IPR) * WTNUM(1) )
          RADD = RADLD1(IPLON,IPR)
          RADLU1(IPLON,IPR) = RADUEMIT(IPLON,IPR) &
               + (1.0_dp - SEMIS(IPLON,IPR)) * RADD
          URAD1(IPLON) = URAD1(IPLON) + RADLU1(IPLON,IPR)
       ENDDO
    ENDDO
    TOTUFLUC(1:KPROMA,0) = URADCL1(1:KPROMA) * 0.5_dp
    TOTUFLUX(1:KPROMA,0) = URAD1(1:KPROMA) * 0.5_dp
    !ELSE
    !- Specular reflection.
    !  DO IPR = 1, JPGPT
    !    DO IPLON = 1, KPROMA
    !      RADCLU = RADUEMIT(IPLON,IPR)
    !      RADCLRU1(IPLON,IPR) = RADCLU + (1.0_dp - SEMIS(IPLON,IPR)) * RADCLRD1(IPLON,IPR)
    !      URADCL1(IPLON) = URADCL1(IPLON) + RADCLRU1(IPLON,IPR)
    !
    !      RADU = RADUEMIT(IPLON,IPR)
    !      RADLU1(IPLON,IPR) = RADU + (1.0_dp - SEMIS(IPLON,IPR)) * RADLD1(IPLON,IPR)
    !      URAD1(IPLON) = URAD1(IPLON) + RADLU1(IPLON,IPR)
    !    ENDDO
    !  ENDDO
    !  TOTUFLUC(1:KPROMA,0) = URADCL1(1:KPROMA) * WTNUM(1)
    !  TOTUFLUX(1:KPROMA,0) = URAD1(1:KPROMA)   * WTNUM(1)
    !ENDIF

    !- Upward radiative transfer.
    !- *** URAD1 holds the summed radiance for total sky stream
    !- *** URADCL1 holds the summed radiance for clear sky stream
    DO LEV = 1, KLEV

       URAD1(1:KPROMA)   = 0.0_dp
       URADCL1(1:KPROMA) = 0.0_dp

       ICCLD = 0
       DO IPLON = 1, KPROMA
          IF (ICLDLYR(IPLON,LEV) == 1) THEN
             ICCLD = ICCLD + 1
             ICCLDL(ICCLD) = IPLON
          ENDIF
       END DO

       IENT = JPGPT * (LEV-1)
       DO IPR = 1, JPGPT
          INDEX = IENT + IPR

!CDIR NODEP
!OCL NOVREC
          DO IX = 1, ICCLD
             IPLON = ICCLDL(IX)
             !- *** Cloudy layer
             !- total-sky upward flux          
             GASSRC = BBU1(IPLON,INDEX) * ABSS1(IPLON,INDEX)

             !- If first cloudy layer in sequence, split up radiance into clear and
             !    cloudy streams depending on cloud fraction
             IF (ISTCLD(IPLON,LEV)  ==  1) THEN
                CLDRADU(IPLON,IPR) = CLDFRAC(IPLON,LEV) * RADLU1(IPLON,IPR)
                CLRRADU(IPLON,IPR) = RADLU1(IPLON,IPR) - CLDRADU(IPLON,IPR)
                RAD(IPLON,IPR) = 0.0_dp
             ENDIF
             TTOT = 1.0_dp - ATOT1(IPLON,INDEX)
             TRNS = 1.0_dp - ABSS1(IPLON,INDEX)
             CLDSRC = BBUTOT1(IPLON,INDEX) * ATOT1(IPLON,INDEX)

             !- Separate RT equations for clear and cloudy streams      
             CLDRADU(IPLON,IPR) = CLDRADU(IPLON,IPR) * TTOT &
                  + CLDFRAC(IPLON,LEV) * CLDSRC
             CLRRADU(IPLON,IPR) = CLRRADU(IPLON,IPR) * TRNS &
                  + (1.0_dp - CLDFRAC(IPLON,LEV)) * GASSRC

             !- total sky upward flux
             RADLU1(IPLON,IPR) = CLDRADU(IPLON,IPR) + CLRRADU(IPLON,IPR)
             URAD1(IPLON) = URAD1(IPLON) + RADLU1(IPLON,IPR)

             !* Code to account for maximum/random overlap:
             !   Performs RT on the radiance most recently switched between clear and
             !   cloudy streams
             RADMOD = RAD(IPLON,IPR) * (FACCLR1(IPLON,LEV+1) * TRNS + &
                  FACCLD1(IPLON,LEV+1) *  TTOT) - &
                  FACCMB1(IPLON,LEV+1) * GASSRC + &
                  FACCMB2(IPLON,LEV+1) * CLDSRC

             !   Computes what the clear and cloudy streams would have been had no
             !   radiance been switched       
             OLDCLD = CLDRADU(IPLON,IPR) - RADMOD
             OLDCLR = CLRRADU(IPLON,IPR) + RADMOD

             !   Computes the radiance to be switched between clear and cloudy.
             RAD(IPLON,IPR) = -RADMOD + FACCLR2(IPLON,LEV+1)*OLDCLR - &
                  FACCLD2(IPLON,LEV+1)*OLDCLD
             CLDRADU(IPLON,IPR) = CLDRADU(IPLON,IPR) + RAD(IPLON,IPR)
             CLRRADU(IPLON,IPR) = CLRRADU(IPLON,IPR) - RAD(IPLON,IPR)

          ENDDO
       ENDDO

       DO IPR = 1, JPGPT
          INDEX = IENT + IPR

          DO IPLON = 1, KPROMA
             IF (ICLDLYR(IPLON,LEV) .NE. 1) THEN

                !- *** Clear layer
                !- total-sky upward flux          
                RADLU1(IPLON,IPR) = RADLU1(IPLON,IPR)+(BBU1(IPLON,INDEX) &
                     - RADLU1(IPLON,IPR))*ABSS1(IPLON,INDEX)
                URAD1(IPLON) = URAD1(IPLON) + RADLU1(IPLON,IPR)

             ENDIF

          ENDDO
       ENDDO

!CDIR UNROLL=4
       DO IPR = 1, JPGPT
          DO IPLON = 1, KPROMA
             INDEX = IENT + IPR

             !- clear-sky upward flux
             !   Upward clear and total sky streams must be separate because surface
             !   reflectance is different for each.
             RADCLRU1(IPLON,IPR) = RADCLRU1(IPLON,IPR)+(BBU1(IPLON,INDEX) &
                  - RADCLRU1(IPLON,IPR))*ABSS1(IPLON,INDEX)
             URADCL1(IPLON) = URADCL1(IPLON) + RADCLRU1(IPLON,IPR)

          ENDDO

       ENDDO

       TOTUFLUC(1:KPROMA,LEV) = URADCL1(1:KPROMA) * WTNUM(1)
       TOTUFLUX(1:KPROMA,LEV) = URAD1(1:KPROMA)   * WTNUM(1)

    END DO

    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! End of source code from Uwe Schulzweidas version
    ! rrtm_rtrn1a_140gp.f90@@/main/clean/2

    RETURN
  END SUBROUTINE rad_lon_RRTM_RTRN1A_140GP
  ! ===========================================================================

#if defined (__SX__)
#  define UNROLLNGX
#endif

!******************************************************************************
!                                                                             *
!                  Optical depths developed for the                           *
!                                                                             *
!                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
!                                                                             *
!            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
!                        840 MEMORIAL DRIVE                                   *
!                        CAMBRIDGE, MA 02139                                  *
!                                                                             *
!                           ELI J. MLAWER                                     *
!                         STEVEN J. TAUBMAN                                   *
!                         SHEPARD A. CLOUGH                                   *
!                                                                             *
!                       email:  mlawer@aer.com                                *
!                                                                             *
!        The authors wish to acknowledge the contributions of the             *
!        following people:  Patrick D. Brown, Michael J. Iacono,              *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                       *
!                                                                             *
!******************************************************************************
! Modified by:                                                                *
!      JJ Morcrette 980714 ECMWF      for use on ECMWF's Fujitsu VPP770       *
!         Reformatted for F90 by JJMorcrette, ECMWF                           * 
!         - replacing COMMONs by MODULEs                                      *
!         - changing labelled to unlabelled DO loops                          *
!         - creating set-up routines for all block data statements            *
!         - reorganizing the parameter statements                             * 
!         - passing KLEV as argument                                          *
!         - suppressing some equivalencing                                    *
!                                                                             *
!      D Salmond    9907   ECMWF      Speed-up modifications                  *
!      D Salmond    000515 ECMWF      Speed-up modifications                  *
!******************************************************************************
!     TAUMOL                                                                  *
!                                                                             *
!     This file contains the subroutines TAUGBn (where n goes from            *
!     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
!     per g-value and layer for band n.                                       *
!                                                                             *
!  Output:  optical depths (unitless)                                         *
!           fractions needed to compute Planck functions at every layer       *
!               and g-value                                                   *
!                                                                             *
!     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
!     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
!                                                                             *
!  Input                                                                      *
!                                                                             *
!     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
!     COMMON /PRECISE/  ONEMINUS                                              *
!     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
!    &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND                        *
!     COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,                              *
!    &                  COLH2O(MXLAY),COLCO2(MXLAY),                          *
!    &                  COLO3(MXLAY),COLN2O(MXLAY),COLCH4(MXLAY),             *
!    &                  COLO2(MXLAY),CO2MULT(MXLAY)                           *
!     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
!    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
!     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)                        *
!     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
!                                                                             *
!     Description:                                                            *
!     NG(IBAND) - number of g-values in band IBAND                            *
!     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
!                   atmospheres that are stored for band IBAND per            *
!                   pressure level and temperature.  Each of these            *
!                   atmospheres has different relative amounts of the         *
!                   key species for the band (i.e. different binary           *
!                   species parameters).                                      *
!     NSPB(IBAND) - same for upper atmosphere                                 *
!     ONEMINUS - since problems are caused in some cases by interpolation     *
!                parameters equal to or greater than 1, for these cases       *
!                these parameters are set to this value, slightly < 1.        *
!     PAVEL - layer pressures (mb)                                            *
!     TAVEL - layer temperatures (degrees K)                                  *
!     PZ - level pressures (mb)                                               *
!     TZ - level temperatures (degrees K)                                     *
!     LAYTROP - layer at which switch is made from one combination of         *
!               key species to another                                        *
!     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
!               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
!               respectively (molecules/cm**2)                                *
!     CO2MULT - for bands in which carbon dioxide is implemented as a         *
!               trace species, this is the factor used to multiply the        *
!               band's average CO2 absorption coefficient to get the added    *
!               contribution to the optical depth relative to 355 ppm.        *
!     FACij(LAY) - for layer LAY, these are factors that are needed to        *
!                  compute the interpolation factors that multiply the        *
!                  appropriate reference k-values.  A value of 0 (1) for      *
!                  i,j indicates that the corresponding factor multiplies     *
!                  reference k-value for the lower (higher) of the two        *
!                  appropriate temperatures, and altitudes, respectively.     *
!     JP - the index of the lower (in altitude) of the two appropriate        *
!          reference pressure levels needed for interpolation                 *
!     JT, JT1 - the indices of the lower of the two appropriate reference     *
!               temperatures needed for interpolation (for pressure           *
!               levels JP and JP+1, respectively)                             *
!     SELFFAC - scale factor needed to water vapor self-continuum, equals     *
!               (water vapor density)/(atmospheric density at 296K and        *
!               1013 mb)                                                      *
!     SELFFRAC - factor needed for temperature interpolation of reference     *
!                water vapor self-continuum data                              *
!     INDSELF - index of the lower of the two appropriate reference           *
!               temperatures needed for the self-continuum interpolation      *
!                                                                             *
!  Data input                                                                 *
!     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG) *
!        (note:  n is the band number)                                        *
!                                                                             *
!     Description:                                                            *
!     KA - k-values for low reference atmospheres (no water vapor             *
!          self-continuum) (units: cm**2/molecule)                            *
!     KB - k-values for high reference atmospheres (all sources)              *
!          (units: cm**2/molecule)                                            *
!     SELFREF - k-values for water vapor self-continuum for reference         *
!               atmospheres (used below LAYTROP)                              *
!               (units: cm**2/molecule)                                       *
!                                                                             *
!     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
!     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
!                                                                             *
!******************************************************************************

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL1 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,&
  &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     Written by Eli J. Mlawer, Atmospheric & Environmental Research.
!     Revised by Michael J. Iacono, Atmospheric & Environmental Research.

!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
 
! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)
REAL(DP):: FORFAC(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)

INTEGER :: IND0(KBDIM),IND1(KBDIM),INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(1) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)   *5+(JT1(IPLON,LAY)-1))*NSPA(1) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG1
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG1
#endif
        TAU (IPLON,IG,LAY) = COLH2O(IPLON,LAY) *              &
             (FAC00(IPLON,LAY) * ABSA_1(IND0(IPLON)  ,IG) +     &
              FAC10(IPLON,LAY) * ABSA_1(IND0(IPLON)+1,IG) +     &
              FAC01(IPLON,LAY) * ABSA_1(IND1(IPLON)  ,IG) +     &
              FAC11(IPLON,LAY) * ABSA_1(IND1(IPLON)+1,IG) +     &
            SELFFAC(IPLON,LAY) * (SELFREF_1(INDS(IPLON),IG) +   &
           SELFFRAC(IPLON,LAY) *                              &
           (SELFREF_1(INDS(IPLON)+1,IG) - SELFREF_1(INDS(IPLON),IG)))  &
             + FORFAC(IPLON,LAY) * FORREF_1(IG) )               &
             + TAUAERL(IPLON,LAY,1)
        PFRAC(IPLON,IG,LAY) = FRACREFA_1(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON)  = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(1) + 1
      IND1(IPLON)  = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(1) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG1
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG1
#endif
        TAU (IPLON,IG,LAY) = COLH2O(IPLON,LAY) *          &
             (FAC00(IPLON,LAY) * ABSB_1(IND0(IPLON)  ,IG) + &
              FAC10(IPLON,LAY) * ABSB_1(IND0(IPLON)+1,IG) + &
              FAC01(IPLON,LAY) * ABSB_1(IND1(IPLON)  ,IG) + &
              FAC11(IPLON,LAY) * ABSB_1(IND1(IPLON)+1,IG)   &
           + FORFAC(IPLON,LAY) * FORREF_1(IG) )             &
          + TAUAERL(IPLON,LAY,1)
        PFRAC(IPLON,IG,LAY) = FRACREFB_1(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL2 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,COLDRY,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,&
  &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up
!     JJMorcrette 2000-07-14 bugfix

IMPLICIT NONE
INTRINSIC :: MAX

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

REAL(DP):: COLDRY(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)
REAL(DP):: FORFAC(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


REAL(DP):: FRACINT(KBDIM)

INTEGER :: INDEX(KBDIM)

INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
INTEGER :: IFPX(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IFP, IFRAC, IG, JFRAC, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FP, H2OPARAM, WATER


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum is 
!     interpolated (in temperature) separately.

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!CDIR NODEP
!OCL NOVREC
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      WATER = 1.E20_dp * COLH2O(IPLON,LAY) / COLDRY(IPLON,LAY)
      H2OPARAM = WATER / (WATER +.002_dp)


      IF (H2OPARAM >= REFPARAM_2(2)) THEN
        INDEX(IPLON) = 2
       ELSE
!CDIR UNROLL=11
!OCL  UNROLL(11)
        DO JFRAC = 2, 12
          IF (H2OPARAM < REFPARAM_2(JFRAC)) THEN
            INDEX(IPLON) = JFRAC + 1
          END IF
        ENDDO
      ENDIF

 
     !---- JJM_000714
      IFRAC = INDEX(IPLON)
      FRACINT(IPLON) = (H2OPARAM-REFPARAM_2(IFRAC)) / &
           (REFPARAM_2(IFRAC-1)-REFPARAM_2(IFRAC))

!write(*,*)'IFRAC', IFRAC
!write(*,*)'IPLON', IPLON
!write(*,*) 'INDEX', INDEX(IPLON)
!write(*,*) 'FRACINT', FRACINT(IPLON)
!write(*,*) 'REFPARAM_2', REFPARAM_2
    ENDDO
!write(*,*)'IFRAC', IFRAC
!write(*,*)'IPLON', IPLON
!write(*,*)'REFPARAM_2', MINVAL(REFPARAM_2)



    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)

      FP = FAC11(IPLON,LAY) + FAC01(IPLON,LAY)
      IFP = 2.E2_dp*FP + 0.5_dp

      !---MI 981104        
      IFPX(IPLON) = MAX(0,IFP)

      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(2) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(2) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG2
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IFRAC = INDEX(IPLON)
      IFP   = IFPX(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=14
!OCL  UNROLL(14)
      DO IG = 1, NG2
#endif
        TAU (IPLON,NGS1+IG,LAY) = COLH2O(IPLON,LAY) *                     &
             (CORR2(IFP) * (FAC00(IPLON,LAY) * ABSA_2(IND0(IPLON)  ,IG)  +  &
                            FAC10(IPLON,LAY) * ABSA_2(IND0(IPLON)+1,IG)) +  &
              CORR1(IFP) * (FAC01(IPLON,LAY) * ABSA_2(IND1(IPLON)  ,IG)  +  &
                            FAC11(IPLON,LAY) * ABSA_2(IND1(IPLON)+1,IG)) +  &
                          SELFFAC(IPLON,LAY) * (SELFREF_2(INDS(IPLON),IG) + &
                         SELFFRAC(IPLON,LAY) *                            &
                         (SELFREF_2(INDS(IPLON)+1,IG) - SELFREF_2(INDS(IPLON),IG))) &
                       + FORFAC(IPLON,LAY) * FORREF_2(IG) )                 &
                       + TAUAERL(IPLON,LAY,2)
        PFRAC(IPLON,NGS1+IG,LAY) = FRACREFA_2(IG,IFRAC) + FRACINT(IPLON) * &
             (FRACREFA_2(IG,IFRAC-1) - FRACREFA_2(IG,IFRAC))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      FP = FAC11(IPLON,LAY) + FAC01(IPLON,LAY)
      IFP = 2.E2_dp*FP + 0.5_dp

      !---MI 981104        
      IFPX(IPLON) = MAX(0,IFP)

      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(2) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(2) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG2
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IFP   = IFPX(IPLON)
! ka_sv_20170406+
      ! avoid out of range in arrays below
      IFP = MIN(MAX(0,IFP),NCORR)
! ka_sv_20170406-
#ifdef UNROLLNGX
!CDIR UNROLL=14
!OCL  UNROLL(14)
      DO IG = 1, NG2
#endif
        TAU (IPLON,NGS1+IG,LAY) = COLH2O(IPLON,LAY) *                    &
             (CORR2(IFP) * (FAC00(IPLON,LAY) * ABSB_2(IND0(IPLON)  ,IG)  + &
                            FAC10(IPLON,LAY) * ABSB_2(IND0(IPLON)+1,IG)) + &
              CORR1(IFP) * (FAC01(IPLON,LAY) * ABSB_2(IND1(IPLON)  ,IG)  + &
                            FAC11(IPLON,LAY) * ABSB_2(IND1(IPLON)+1,IG))   &
                         + FORFAC(IPLON,LAY) * FORREF_2(IG) )              &
                        + TAUAERL(IPLON,LAY,2)
        PFRAC(IPLON,NGS1+IG,LAY) = FRACREFB_2(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL2
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL3 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,FORFAC,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)

! Modifications
!
!     D Salmond 2000-05-15 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)
REAL(DP):: FORFAC(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM), N2OMULT(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, NS, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: COLREF1, COLREF2, CURRN2O,                      &
           FP, FS, RATIO, SPECMULT, SPECPARM, WCOMB1,      &
           WCOMB2


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!CDIR NODEP
!OCL NOVREC
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT_3*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IF (JS  ==  8) THEN
        IF (FS  >=  0.9_dp) THEN
          JS = 9
          FS = 10._dp * (FS - 0.9_dp)
        ELSE
          FS = FS / 0.9_dp
        ENDIF
      ENDIF

      NS = JS + INT(FS + 0.5_dp)
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(3) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(3) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      COLREF1 = N2OREF_3(JP(IPLON,LAY))
      COLREF2 = N2OREF_3(JP(IPLON,LAY)+1)
      IF (NS  ==  10) THEN
        WCOMB1 = 1.0_dp/H2OREF_3(JP(IPLON,LAY))
        WCOMB2 = 1.0_dp/H2OREF_3(JP(IPLON,LAY)+1)
      ELSE
        WCOMB1 = (1.0_dp-ETAREF_3(NS))/(STRRAT_3 * CO2REF_3(JP(IPLON,LAY)))
        WCOMB2 = (1.0_dp-ETAREF_3(NS))/(STRRAT_3 * CO2REF_3(JP(IPLON,LAY)+1))
      ENDIF
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = SPECCOMB(IPLON) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG3
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=16
!OCL  UNROLL(16)
      DO IG = 1, NG3
#endif
        TAU (IPLON,NGS2+IG,LAY) = SPECCOMB(IPLON) *                      &
             ((1._DP - FS) *(FAC00(IPLON,LAY) * ABSA_3(IND0(IPLON)   ,IG) +   &
                          FAC10(IPLON,LAY) * ABSA_3(IND0(IPLON)+10,IG) +   &
                          FAC01(IPLON,LAY) * ABSA_3(IND1(IPLON)   ,IG) +   &
                          FAC11(IPLON,LAY) * ABSA_3(IND1(IPLON)+10,IG))+   &
                    FS * (FAC00(IPLON,LAY) * ABSA_3(IND0(IPLON)+ 1,IG) +   &
                          FAC10(IPLON,LAY) * ABSA_3(IND0(IPLON)+11,IG) +   &
                          FAC01(IPLON,LAY) * ABSA_3(IND1(IPLON)+ 1,IG) +   &
                          FAC11(IPLON,LAY) * ABSA_3(IND1(IPLON)+11,IG))) + &
                          COLH2O(IPLON,LAY) *                            &
              SELFFAC(IPLON,LAY) * (SELFREF_3(INDS(IPLON),IG) +            &
              SELFFRAC(IPLON,LAY) *                                      &
              (SELFREF_3(INDS(IPLON)+1,IG) - SELFREF_3(INDS(IPLON),IG))      &
              + FORFAC(IPLON,LAY) * FORREF_3(IG) )                         &
              + N2OMULT(IPLON) * ABSN2OA_3(IG)                             &
              + TAUAERL(IPLON,LAY,3)
        PFRAC(IPLON,NGS2+IG,LAY) = FRACREFA_3(IG,JS) + FS *         &
              (FRACREFA_3(IG,JS+1) - FRACREFA_3(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
!CDIR NODEP
!OCL NOVREC
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT_3*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 4._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      NS = JS + INT(FS + 0.5_dp)
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(3) + JS
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(3) + JS
      COLREF1 = N2OREF_3(JP(IPLON,LAY))
      COLREF2 = N2OREF_3(JP(IPLON,LAY)+1)
      IF (NS  ==  5) THEN
        WCOMB1 = 1.0_dp/H2OREF_3(JP(IPLON,LAY))
        WCOMB2 = 1.0_dp/H2OREF_3(JP(IPLON,LAY)+1)
      ELSE
        WCOMB1 = (1.0_dp-ETAREF_3(NS))/(STRRAT_3 * CO2REF_3(JP(IPLON,LAY)))
        WCOMB2 = (1.0_dp-ETAREF_3(NS))/(STRRAT_3 * CO2REF_3(JP(IPLON,LAY)+1))
      ENDIF
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = SPECCOMB(IPLON) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG3
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=16
!OCL  UNROLL(16)
      DO IG = 1, NG3
#endif
        TAU (IPLON,NGS2+IG,LAY) = SPECCOMB(IPLON) *                    &
             ((1._DP - FS) *(FAC00(IPLON,LAY) * ABSB_3(IND0(IPLON)  ,IG) +  &
                          FAC10(IPLON,LAY) * ABSB_3(IND0(IPLON)+5,IG) +  &
                          FAC01(IPLON,LAY) * ABSB_3(IND1(IPLON)  ,IG) +  &
                          FAC11(IPLON,LAY) * ABSB_3(IND1(IPLON)+5,IG))+  &
                    FS * (FAC01(IPLON,LAY) * ABSB_3(IND1(IPLON)+1,IG) +  &
                          FAC10(IPLON,LAY) * ABSB_3(IND0(IPLON)+6,IG) +  &
                          FAC00(IPLON,LAY) * ABSB_3(IND0(IPLON)+1,IG) +  &
                          FAC11(IPLON,LAY) * ABSB_3(IND1(IPLON)+6,IG)))  &
                       + COLH2O(IPLON,LAY)*FORFAC(IPLON,LAY)*FORREF_3(IG)  &
                       + N2OMULT(IPLON) * ABSN2OB_3(IG)                  &
                       + TAUAERL(IPLON,LAY,3)
        PFRAC(IPLON,NGS2+IG,LAY) = FRACREFB_3(IG,JS) + FS *       &
             (FRACREFB_3(IG,JS+1) - FRACREFB_3(IG,JS))
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL3
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL4 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLO3,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)

! Modifications
!
!     D Salmond 2000-05-15 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLO3(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT1_4*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp * (SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(4) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(4) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG4
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=14
!OCL  UNROLL(14)
      DO IG = 1, NG4
#endif
        TAU (IPLON,NGS3+IG,LAY) = SPECCOMB(IPLON) *                     &
             ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSA_4(IND0(IPLON)   ,IG)  +  &
                         FAC10(IPLON,LAY) * ABSA_4(IND0(IPLON)+ 9,IG)  +  &
                         FAC01(IPLON,LAY) * ABSA_4(IND1(IPLON)   ,IG)  +  &
                         FAC11(IPLON,LAY) * ABSA_4(IND1(IPLON)+ 9,IG)) +  &
                   FS * (FAC01(IPLON,LAY) * ABSA_4(IND1(IPLON)+ 1,IG)  +  &
                         FAC10(IPLON,LAY) * ABSA_4(IND0(IPLON)+10,IG)  +  &
                         FAC00(IPLON,LAY) * ABSA_4(IND0(IPLON)+ 1,IG)  +  &
                         FAC11(IPLON,LAY) * ABSA_4(IND1(IPLON)+10,IG))) + &
                        COLH2O(IPLON,LAY) *                             &
                       SELFFAC(IPLON,LAY) * (SELFREF_4(INDS(IPLON),IG) +  &
                       SELFFRAC(IPLON,LAY) *                            &
                       (SELFREF_4(INDS(IPLON)+1,IG) - SELFREF_4(INDS(IPLON),IG)))  &
                      + TAUAERL(IPLON,LAY,4)
        PFRAC(IPLON,NGS3+IG,LAY) = FRACREFA_4(IG,JS) + FS *        &
             (FRACREFA_4(IG,JS+1) - FRACREFA_4(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      SPECCOMB(IPLON) = COLO3(IPLON,LAY) + STRRAT2_4*COLCO2(IPLON,LAY)
      SPECPARM = COLO3(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 4._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IF (JS  >  1) THEN
        JS = JS + 1
        ELSEIF (FS  >=  0.0024_dp) THEN
        JS = 2
        FS = (FS - 0.0024_dp)/0.9976_dp
      ELSE
        JS = 1
        FS = FS/0.0024_dp
      ENDIF
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(4) + JS
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(4) + JS
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG4
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=14
!OCL  UNROLL(14)
      DO IG = 1, NG4
#endif
        TAU (IPLON,NGS3+IG,LAY) = SPECCOMB(IPLON) *                   &
             ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSB_4(IND0(IPLON)  ,IG) +  &
                         FAC10(IPLON,LAY) * ABSB_4(IND0(IPLON)+6,IG) +  &
                         FAC01(IPLON,LAY) * ABSB_4(IND1(IPLON)  ,IG) +  &
                         FAC11(IPLON,LAY) * ABSB_4(IND1(IPLON)+6,IG)) + &
                   FS * (FAC10(IPLON,LAY) * ABSB_4(IND0(IPLON)+7,IG) +  &
                         FAC01(IPLON,LAY) * ABSB_4(IND1(IPLON)+1,IG) +  &
                         FAC00(IPLON,LAY) * ABSB_4(IND0(IPLON)+1,IG) +  &
                         FAC11(IPLON,LAY) * ABSB_4(IND1(IPLON)+7,IG)))  &
                     + TAUAERL(IPLON,LAY,4)
        PFRAC(IPLON,NGS3+IG,LAY) = FRACREFB_4(IG,JS) + FS *      &
              (FRACREFB_4(IG,JS+1) - FRACREFB_4(IG,JS))
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL4
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL5 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,WX,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2, COLO3,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)

! Modifications
!
!     D Salmond 2000-05-15 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

REAL(DP):: WX(KBDIM,JPXSEC,KLEV)
!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLO3(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT1_5*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(5) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(5) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG5
#endif
!!! !OCL NOVREC
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=16
!OCL  UNROLL(16)
      DO IG = 1, NG5
#endif
        TAU (IPLON,NGS4+IG,LAY) = SPECCOMB(IPLON) *                   &
          ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSA_5(IND0(IPLON),IG) +       &
                      FAC10(IPLON,LAY) * ABSA_5(IND0(IPLON)+9,IG) +     &
                      FAC11(IPLON,LAY) * ABSA_5(IND1(IPLON)+9,IG) +     &
                      FAC01(IPLON,LAY) * ABSA_5(IND1(IPLON),IG)) +      &
                 FS* (FAC10(IPLON,LAY) * ABSA_5(IND0(IPLON)+10,IG) +    &
                      FAC01(IPLON,LAY) * ABSA_5(IND1(IPLON)+1,IG) +     &
                      FAC00(IPLON,LAY) * ABSA_5(IND0(IPLON)+1,IG) +     &
                      FAC11(IPLON,LAY) * ABSA_5(IND1(IPLON)+10,IG)) ) + &
                   COLH2O(IPLON,LAY) *                                &
                   SELFFAC(IPLON,LAY) * (SELFREF_5(INDS(IPLON),IG) +    &
                   SELFFRAC(IPLON,LAY) *                              &
                  (SELFREF_5(INDS(IPLON)+1,IG) - SELFREF_5(INDS(IPLON),IG)))     &
                + WX(IPLON,1,LAY) * CCL4_5(IG)                          &
                + TAUAERL(IPLON,LAY,5)
        PFRAC(IPLON,NGS4+IG,LAY) = FRACREFA_5(IG,JS) + FS * &
             (FRACREFA_5(IG,JS+1) - FRACREFA_5(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      SPECCOMB(IPLON) = COLO3(IPLON,LAY) + STRRAT2_5*COLCO2(IPLON,LAY)
      SPECPARM = COLO3(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 4._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(5) + JS
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(5) + JS
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG5
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=16
!OCL  UNROLL(16)
      DO IG = 1, NG5
#endif
        TAU (IPLON,NGS4+IG,LAY) = SPECCOMB(IPLON) *               &
        ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSB_5(IND0(IPLON),IG) +     &
                    FAC10(IPLON,LAY) * ABSB_5(IND0(IPLON)+5,IG) +   &
                    FAC01(IPLON,LAY) * ABSB_5(IND1(IPLON),IG) +     &
                    FAC11(IPLON,LAY) * ABSB_5(IND1(IPLON)+5,IG) ) + &
              FS * (FAC01(IPLON,LAY) * ABSB_5(IND1(IPLON)+1,IG) +   &
                    FAC10(IPLON,LAY) * ABSB_5(IND0(IPLON)+6,IG) +   &
                    FAC00(IPLON,LAY) * ABSB_5(IND0(IPLON)+1,IG) +   &
                    FAC11(IPLON,LAY) * ABSB_5(IND1(IPLON)+6,IG)))   &
             + WX(IPLON,1,LAY) * CCL4_5(IG)                         &
             + TAUAERL(IPLON,LAY,5)
        PFRAC(IPLON,NGS4+IG,LAY) = FRACREFB_5(IG,JS) + FS *  &
             (FRACREFB_5(IG,JS+1) - FRACREFB_5(IG,JS))
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL5
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL6 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,WX,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

REAL(DP):: WX(KBDIM,JPXSEC,KLEV)
!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: CO2MULT(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and
!     temperature. The water vapor self-continuum is interpolated
!     (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(6) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(6) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG6
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG6
#endif
        TAU (IPLON,NGS5+IG,LAY) = COLH2O(IPLON,LAY) *        &
             (FAC00(IPLON,LAY) * ABSA_6(IND0(IPLON)  ,IG) +    &
              FAC10(IPLON,LAY) * ABSA_6(IND0(IPLON)+1,IG) +    &
              FAC01(IPLON,LAY) * ABSA_6(IND1(IPLON)  ,IG) +    &
              FAC11(IPLON,LAY) * ABSA_6(IND1(IPLON)+1,IG) +    &
             SELFFAC(IPLON,LAY) * (SELFREF_6(INDS(IPLON),IG) + &
             SELFFRAC(IPLON,LAY)*                            &
             (SELFREF_6(INDS(IPLON)+1,IG)-SELFREF_6(INDS(IPLON),IG))))  &
             + WX(IPLON,2,LAY) * CFC11ADJ_6(IG)                &
             + WX(IPLON,3,LAY) * CFC12_6(IG)                   &
             + CO2MULT(IPLON,LAY) * ABSCO2_6(IG)               &
             + TAUAERL(IPLON,LAY,6)
        PFRAC(IPLON,NGS5+IG,LAY) = FRACREFA_6(IG)
      ENDDO
    ENDDO

    !     Nothing important goes on above LAYTROP in this band.
    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
    DO IG = 1, NG6
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG6
#endif
        TAU (IPLON,NGS5+IG,LAY) = 0.0_dp      &
             + WX(IPLON,2,LAY) * CFC11ADJ_6(IG) &
             + WX(IPLON,3,LAY) * CFC12_6(IG)    &
             + TAUAERL(IPLON,LAY,6)
        PFRAC(IPLON,NGS5+IG,LAY) = FRACREFA_6(IG)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL6
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL7 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLO3,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 7:  980-1080 cm-1 (low - H2O,O3; high - O3)

! Modifications
!
!     D Salmond 2000-05-15 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLO3(KBDIM,KLEV)
REAL(DP):: CO2MULT(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.
  
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT_7*COLO3(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*SPECPARM
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(7) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(7) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG7
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=12
!OCL  UNROLL(12)
      DO IG = 1, NG7
#endif
        TAU (IPLON,NGS6+IG,LAY) = SPECCOMB(IPLON) *                     &
             ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSA_7(IND0(IPLON),IG) +      &
                         FAC10(IPLON,LAY) * ABSA_7(IND0(IPLON)+9,IG) +    &
                         FAC01(IPLON,LAY) * ABSA_7(IND1(IPLON),IG) +      &
                         FAC11(IPLON,LAY) * ABSA_7(IND1(IPLON)+9,IG) )+   &
                   FS * (FAC01(IPLON,LAY) * ABSA_7(IND1(IPLON)+1,IG) +    &
                         FAC10(IPLON,LAY) * ABSA_7(IND0(IPLON)+10,IG) +   &
                         FAC00(IPLON,LAY) * ABSA_7(IND0(IPLON)+1,IG) +    &
                         FAC11(IPLON,LAY) * ABSA_7(IND1(IPLON)+10,IG))) + &
                        COLH2O(IPLON,LAY) *                             &
                       SELFFAC(IPLON,LAY) * (SELFREF_7(INDS(IPLON),IG) +  &
                      SELFFRAC(IPLON,LAY) *                             &
                     (SELFREF_7(INDS(IPLON)+1,IG) - SELFREF_7(INDS(IPLON),IG)))    &
                    + CO2MULT(IPLON,LAY) * ABSCO2_7(IG)                   &
                    + TAUAERL(IPLON,LAY,7)
        PFRAC(IPLON,NGS6+IG,LAY) = FRACREFA_7(IG,JS) + FS * &
             (FRACREFA_7(IG,JS+1) - FRACREFA_7(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(7) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(7) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG7
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=12
!OCL  UNROLL(12)
      DO IG = 1, NG7
#endif
        TAU (IPLON,NGS6+IG,LAY) = COLO3(IPLON,LAY) *        &
             (FAC00(IPLON,LAY) * ABSB_7(IND0(IPLON)  ,IG) +   &
              FAC10(IPLON,LAY) * ABSB_7(IND0(IPLON)+1,IG) +   &
              FAC01(IPLON,LAY) * ABSB_7(IND1(IPLON)  ,IG) +   &
              FAC11(IPLON,LAY) * ABSB_7(IND1(IPLON)+1,IG))    &
             + CO2MULT(IPLON,LAY) * ABSCO2_7(IG)              &
             + TAUAERL(IPLON,LAY,7)
        PFRAC(IPLON,NGS6+IG,LAY) = FRACREFB_7(IG)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL7
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL8 (KPROMA,KBDIM,KLEV,IXS,IXLOS,IXHIGS,TAU,WX,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,COLO3,COLN2O,CO2MULT,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 8:  1080-1180 cm-1 (low (i.e.>~300mb) - H2O; high - O3)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXS(KLEV), IXLOS(KBDIM,KLEV), IXHIGS(KBDIM,KLEV)

REAL(DP):: WX(KBDIM,JPXSEC,KLEV)
!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLO3(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)
REAL(DP):: CO2MULT(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: COLREF1, COLREF2, CURRN2O, FP, RATIO, WCOMB1, WCOMB2
REAL(DP):: N2OMULT(KBDIM)


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  

  DO LAY = 1, KLEV
    IXC0 = IXS(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOS(IXP,LAY)
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT(IPLON,LAY)-1))*NSPA(8) + 1
      IND1(IPLON) = (JP(IPLON,LAY)*5+(JT1(IPLON,LAY)-1))*NSPA(8) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
      COLREF1 = N2OREF_8(JP(IPLON,LAY))
      COLREF2 = N2OREF_8(JP(IPLON,LAY)+1)
      WCOMB1 = 1.0_dp/H2OREF_8(JP(IPLON,LAY))
      WCOMB2 = 1.0_dp/H2OREF_8(JP(IPLON,LAY)+1)
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = COLH2O(IPLON,LAY) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG8
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOS(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG8
#endif
        TAU (IPLON,NGS7+IG,LAY) = COLH2O(IPLON,LAY) *         &
             (FAC00(IPLON,LAY) * ABSA_8(IND0(IPLON)  ,IG) +     &
              FAC10(IPLON,LAY) * ABSA_8(IND0(IPLON)+1,IG) +     &
              FAC01(IPLON,LAY) * ABSA_8(IND1(IPLON)  ,IG) +     &
              FAC11(IPLON,LAY) * ABSA_8(IND1(IPLON)+1,IG) +     &
             SELFFAC(IPLON,LAY) * (SELFREF_8(INDS(IPLON),IG) +  &
             SELFFRAC(IPLON,LAY) *                            &
             (SELFREF_8(INDS(IPLON)+1,IG) - SELFREF_8(INDS(IPLON),IG)))) &
             + WX(IPLON,3,LAY) * CFC12_8(IG)                    &
             + WX(IPLON,4,LAY) * CFC22ADJ_8(IG)                 &
             + CO2MULT(IPLON,LAY) * ABSCO2A_8(IG)               &
             + N2OMULT(IPLON) * ABSN2OA_8(IG)                   &
             + TAUAERL(IPLON,LAY,8)
        PFRAC(IPLON,NGS7+IG,LAY) = FRACREFA_8(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGS(IXP,LAY)
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-7)*5+(JT (IPLON,LAY)-1))*NSPB(8) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-6)*5+(JT1(IPLON,LAY)-1))*NSPB(8) + 1
      COLREF1 = N2OREF_8(JP(IPLON,LAY))
      COLREF2 = N2OREF_8(JP(IPLON,LAY)+1)
      WCOMB1 = 1.0_dp/O3REF_8(JP(IPLON,LAY))
      WCOMB2 = 1.0_dp/O3REF_8(JP(IPLON,LAY)+1)
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = COLO3(IPLON,LAY) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG8
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGS(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG8
#endif
        TAU (IPLON,NGS7+IG,LAY) = COLO3(IPLON,LAY) *        &
             (FAC00(IPLON,LAY) * ABSB_8(IND0(IPLON)  ,IG) +   &
              FAC10(IPLON,LAY) * ABSB_8(IND0(IPLON)+1,IG) +   &
              FAC01(IPLON,LAY) * ABSB_8(IND1(IPLON)  ,IG) +   &
              FAC11(IPLON,LAY) * ABSB_8(IND1(IPLON)+1,IG))    &
             + WX(IPLON,3,LAY) * CFC12_8(IG)                  &
             + WX(IPLON,4,LAY) * CFC22ADJ_8(IG)               &
             + CO2MULT(IPLON,LAY) * ABSCO2B_8(IG)             &
             + N2OMULT(IPLON) * ABSN2OB_8(IG)                 &
             + TAUAERL(IPLON,LAY,8)
        PFRAC(IPLON,NGS7+IG,LAY) = FRACREFB_8(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL8
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL9 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLN2O,COLCH4,LAYSWTCH,LAYLOW,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 9:  1180-1390 cm-1 (low - H2O,CH4; high - CH4)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)
REAL(DP):: COLCH4(KBDIM,KLEV)
INTEGER :: LAYSWTCH(KBDIM)
INTEGER :: LAYLOW(KBDIM)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM), N2OMULT(KBDIM), FFRAC(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, IOFF(KBDIM), JS, LAY, NS, IPLON
INTEGER :: IXP, IXC0
INTEGER :: JFRAC

!     LOCAL REAL SCALARS
REAL(DP):: COLREF1, COLREF2, CURRN2O, &
           FP, FS, RATIO, SPECMULT, SPECPARM, WCOMB1, &
           WCOMB2

!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.
  

  IOFF(:) = 0

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

!CDIR NODEP
!OCL NOVREC
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT_9*COLCH4(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp * (SPECPARM)
      JS = 1 + INT(SPECMULT)
      JFRAC = JS
      FS = MOD(SPECMULT,1.0_dp)
      FFRAC(IPLON) = FS
      IF (JS  ==  8) THEN
        IF (FS .LE. 0.68_dp) THEN
          FS = FS/0.68_dp
        ELSEIF (FS  <=  0.92_dp) THEN
          JS = JS + 1
          FS = (FS-0.68_dp)/0.24_dp
        ELSE
          JS = JS + 2
          FS = (FS-0.92_dp)/0.08_dp
        ENDIF
      ELSEIF (JS == 9) THEN
        JS = 10
        FS = 1.0_dp
        JFRAC = 8
        FFRAC(IPLON) = 1.0_dp
      ENDIF
      FP = FAC01(IPLON,LAY) + FAC11(IPLON,LAY)
      NS = JS + INT(FS + 0.5_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(9) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(9) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      IF (LAY  ==  LAYLOW(IPLON)) IOFF(IPLON) = NG9
      IF (LAY  ==  LAYSWTCH(IPLON)) IOFF(IPLON) = 2*NG9
      COLREF1 = N2OREF_9(JP(IPLON,LAY))
      COLREF2 = N2OREF_9(JP(IPLON,LAY)+1)
      IF (NS  ==  11) THEN
        WCOMB1 = 1.0_dp/H2OREF_9(JP(IPLON,LAY))
        WCOMB2 = 1.0_dp/H2OREF_9(JP(IPLON,LAY)+1)
      ELSE
        WCOMB1 = (1.0_dp-ETAREF_9(NS))/(STRRAT_9 * CH4REF_9(JP(IPLON,LAY)))
        WCOMB2 = (1.0_dp-ETAREF_9(NS))/(STRRAT_9 * CH4REF_9(JP(IPLON,LAY)+1))
      ENDIF
      RATIO = (COLREF1*WCOMB1)+FP*((COLREF2*WCOMB2)-(COLREF1*WCOMB1))
      CURRN2O = SPECCOMB(IPLON) * RATIO
      N2OMULT(IPLON) = COLN2O(IPLON,LAY) - CURRN2O
      ZFS(IPLON) = FS
!    ENDDO

!#ifndef UNROLLNGX
!    DO IG = 1, NG9
!#endif
!    DO IXP = 1, IXC0
!      IPLON = IXLOW(IXP,LAY)
!      FS = ZFS(IPLON)
!#ifdef UNROLLNGX
!CDIR UNROLL=12
!OCL  UNROLL(12)
      DO IG = 1, NG9
!#endif
        TAU (IPLON,NGS8+IG,LAY) = SPECCOMB(IPLON) *                       &
           ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSA_9(IND0(IPLON)   ,IG) +       &
                       FAC10(IPLON,LAY) * ABSA_9(IND0(IPLON)+11,IG) +       &
                       FAC01(IPLON,LAY) * ABSA_9(IND1(IPLON)   ,IG) +       &
                       FAC11(IPLON,LAY) * ABSA_9(IND1(IPLON)+11,IG)) +      &
               FS* (   FAC00(IPLON,LAY) * ABSA_9(IND0(IPLON)+ 1,IG) +       &
                       FAC10(IPLON,LAY) * ABSA_9(IND0(IPLON)+12,IG) +       &
                       FAC01(IPLON,LAY) * ABSA_9(IND1(IPLON)+ 1,IG) +       &
                       FAC11(IPLON,LAY) * ABSA_9(IND1(IPLON)+12,IG))) +     &
                      COLH2O(IPLON,LAY) *                                 &
                     SELFFAC(IPLON,LAY) * (SELFREF_9(INDS(IPLON),IG) +      &
                    SELFFRAC(IPLON,LAY) *                                 &
                   (SELFREF_9(INDS(IPLON)+1,IG) - SELFREF_9(INDS(IPLON),IG))) &
                 + N2OMULT(IPLON) * ABSN2O_9(IG+IOFF(IPLON))                &
                 + TAUAERL(IPLON,LAY,9)
        PFRAC(IPLON,NGS8+IG,LAY) = FRACREFA_9(IG,JFRAC) + FFRAC(IPLON) *&
             (FRACREFA_9(IG,JFRAC+1) - FRACREFA_9(IG,JFRAC))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(9) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(9) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG9
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=12
!OCL  UNROLL(12)
      DO IG = 1, NG9
#endif
        TAU (IPLON,NGS8+IG,LAY) = COLCH4(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSB_9(IND0(IPLON)  ,IG) + &
              FAC10(IPLON,LAY) * ABSB_9(IND0(IPLON)+1,IG) + &
              FAC01(IPLON,LAY) * ABSB_9(IND1(IPLON)  ,IG) + &
              FAC11(IPLON,LAY) * ABSB_9(IND1(IPLON)+1,IG))  &
             + TAUAERL(IPLON,LAY,9)
        PFRAC(IPLON,NGS8+IG,LAY) = FRACREFB_9(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL9
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL10 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,PFRAC)

!     BAND 10:  1390-1480 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(10) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(10) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG10
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=6
!OCL  UNROLL(6)
      DO IG = 1, NG10
#endif
        TAU (IPLON,NGS9+IG,LAY) = COLH2O(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSA_10(IND0(IPLON)  ,IG) + &
              FAC10(IPLON,LAY) * ABSA_10(IND0(IPLON)+1,IG) + &
              FAC01(IPLON,LAY) * ABSA_10(IND1(IPLON)  ,IG) + &
              FAC11(IPLON,LAY) * ABSA_10(IND1(IPLON)+1,IG))  &
             + TAUAERL(IPLON,LAY,10)
        PFRAC(IPLON,NGS9+IG,LAY) = FRACREFA_10(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(10) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(10) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG10
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=6
!OCL  UNROLL(6)
      DO IG = 1, NG10
#endif
        TAU (IPLON,NGS9+IG,LAY) = COLH2O(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSB_10(IND0(IPLON)  ,IG) + &
              FAC10(IPLON,LAY) * ABSB_10(IND0(IPLON)+1,IG) + &
              FAC01(IPLON,LAY) * ABSB_10(IND1(IPLON)  ,IG) + &
              FAC11(IPLON,LAY) * ABSB_10(IND1(IPLON)+1,IG))  &
             + TAUAERL(IPLON,LAY,10)
        PFRAC(IPLON,NGS9+IG,LAY) = FRACREFB_10(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL10
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL11 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLH2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 11:  1480-1800 cm-1 (low - H2O; high - H2O)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.
  
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT( IPLON,LAY)-1))*NSPA(11) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(11) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG11
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG11
#endif
        TAU (IPLON,NGS10+IG,LAY) = COLH2O(IPLON,LAY) *        &
             (FAC00(IPLON,LAY) * ABSA_11(IND0(IPLON)  ,IG) +     &
              FAC10(IPLON,LAY) * ABSA_11(IND0(IPLON)+1,IG) +     &
              FAC01(IPLON,LAY) * ABSA_11(IND1(IPLON)  ,IG) +     &
              FAC11(IPLON,LAY) * ABSA_11(IND1(IPLON)+1,IG) +     &
             SELFFAC(IPLON,LAY) * (SELFREF_11(INDS(IPLON),IG) +  &
             SELFFRAC(IPLON,LAY) *                            &
             (SELFREF_11(INDS(IPLON)+1,IG) - SELFREF_11(INDS(IPLON),IG)))) &
             + TAUAERL(IPLON,LAY,11)
        PFRAC(IPLON,NGS10+IG,LAY) = FRACREFA_11(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(11) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(11) + 1
   ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG11
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG11
#endif
        TAU (IPLON,NGS10+IG,LAY) = COLH2O(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSB_11(IND0(IPLON)  ,IG) +  &
              FAC10(IPLON,LAY) * ABSB_11(IND0(IPLON)+1,IG) +  &
              FAC01(IPLON,LAY) * ABSB_11(IND1(IPLON)  ,IG) +  &
              FAC11(IPLON,LAY) * ABSB_11(IND1(IPLON)+1,IG))   &
             + TAUAERL(IPLON,LAY,11)
        PFRAC(IPLON,NGS10+IG,LAY) = FRACREFB_11(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL11
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL12 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 12:  1800-2080 cm-1 (low - H2O,CO2; high - nothing)

! Modifications
!
!     D Salmond   2000-05-15 speed-up
!     JJMorcrette 2000-05-17 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL (dp) :: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL (dp) :: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT_12*COLCO2(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(12) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(12) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG12
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG12
#endif
        TAU (IPLON,NGS11+IG,LAY) = SPECCOMB(IPLON) *                    &
             ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSA_12(IND0(IPLON)  ,IG) +    &
                         FAC10(IPLON,LAY) * ABSA_12(IND0(IPLON)+9,IG) +    &
                         FAC01(IPLON,LAY) * ABSA_12(IND1(IPLON)  ,IG) +    &
                         FAC11(IPLON,LAY) * ABSA_12(IND1(IPLON)+9,IG)) +   &
               FS *(     FAC01(IPLON,LAY) * ABSA_12(IND1(IPLON)+ 1,IG) +   &
                         FAC00(IPLON,LAY) * ABSA_12(IND0(IPLON)+ 1,IG) +   &
                         FAC10(IPLON,LAY) * ABSA_12(IND0(IPLON)+10,IG) +   &
                         FAC11(IPLON,LAY) * ABSA_12(IND1(IPLON)+10,IG))) + &
                        COLH2O(IPLON,LAY) *                             &
                       SELFFAC(IPLON,LAY) * (SELFREF_12(INDS(IPLON),IG) +  &
                      SELFFRAC(IPLON,LAY) *                             &
                     (SELFREF_12(INDS(IPLON)+1,IG) - SELFREF_12(INDS(IPLON),IG)))  &
                    + TAUAERL(IPLON,LAY,12)
        PFRAC(IPLON,NGS11+IG,LAY) = FRACREFA_12(IG,JS) + FS * &
             (FRACREFA_12(IG,JS+1) - FRACREFA_12(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
    DO IG = 1, NG12
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=8
!OCL  UNROLL(8)
      DO IG = 1, NG12
#endif
        TAU  (IPLON,NGS11+IG,LAY) = TAUAERL(IPLON,LAY,12)
        PFRAC(IPLON,NGS11+IG,LAY) = 0.0_dp
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL12
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL13 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 13:  2080-2250 cm-1 (low - H2O,N2O; high - nothing)

! Modifications
!
!     D Salmond 2000-05-15 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT_13*COLN2O(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(ONEMINUS,SPECPARM)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(13) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(13) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG13
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=4
!OCL  UNROLL(4)
      DO IG = 1, NG13
#endif
        TAU (IPLON,NGS12+IG,LAY) = SPECCOMB(IPLON) *                     &
             ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSA_13(IND0(IPLON),IG) +       &
                         FAC10(IPLON,LAY) * ABSA_13(IND0(IPLON)+9,IG) +     &
                         FAC01(IPLON,LAY) * ABSA_13(IND1(IPLON),IG) +       &
                         FAC11(IPLON,LAY) * ABSA_13(IND1(IPLON)+9,IG)) +    &
                 FS* (   FAC01(IPLON,LAY) * ABSA_13(IND1(IPLON)+1,IG) +     &
                         FAC10(IPLON,LAY) * ABSA_13(IND0(IPLON)+10,IG) +    &
                         FAC00(IPLON,LAY) * ABSA_13(IND0(IPLON)+1,IG) +     &
                         FAC11(IPLON,LAY) * ABSA_13(IND1(IPLON)+10,IG)) ) + &
                        COLH2O(IPLON,LAY) *                              &
                       SELFFAC(IPLON,LAY) * (SELFREF_13(INDS(IPLON),IG) +   &
                      SELFFRAC(IPLON,LAY) *                              &
                     (SELFREF_13(INDS(IPLON)+1,IG) - SELFREF_13(INDS(IPLON),IG))) &
                    + TAUAERL(IPLON,LAY,13)
        PFRAC(IPLON,NGS12+IG,LAY) = FRACREFA_13(IG,JS) + FS * &
             (FRACREFA_13(IG,JS+1) - FRACREFA_13(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
    DO IG = 1, NG13
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=4
!OCL  UNROLL(4)
      DO IG = 1, NG13
#endif
        TAU  (IPLON,NGS12+IG,LAY) = TAUAERL(IPLON,LAY,13)
        PFRAC(IPLON,NGS12+IG,LAY) = 0.0_dp
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL13
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL14 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,&
  &COLCO2,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 14:  2250-2380 cm-1 (low - CO2; high - CO2)

! Modifications
!
!     D Salmond 1999-07-14 speed-up

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PROFDATA             
REAL(DP):: COLCO2(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, LAY, IPLON
INTEGER :: IXP, IXC0


!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum 
!     is interpolated (in temperature) separately.  

  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(14) + 1
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(14) + 1
      INDS(IPLON) = INDSELF(IPLON,LAY)
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG14
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=2
!OCL  UNROLL(2)
      DO IG = 1, NG14
#endif
        TAU (IPLON,NGS13+IG,LAY) = COLCO2(IPLON,LAY) *         &
             (FAC00(IPLON,LAY) * ABSA_14(IND0(IPLON)  ,IG) +      &
              FAC10(IPLON,LAY) * ABSA_14(IND0(IPLON)+1,IG) +      &
              FAC01(IPLON,LAY) * ABSA_14(IND1(IPLON)  ,IG) +      &
              FAC11(IPLON,LAY) * ABSA_14(IND1(IPLON)+1,IG) +      &
            SELFFAC(IPLON,LAY) * (SELFREF_14(INDS(IPLON),IG) +    &
           SELFFRAC(IPLON,LAY) *                               &
           (SELFREF_14(INDS(IPLON)+1,IG) - SELFREF_14(INDS(IPLON),IG)))) &
          + TAUAERL(IPLON,LAY,14)
        PFRAC(IPLON,NGS13+IG,LAY) = FRACREFA_14(IG)
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
      IND0(IPLON) = ((JP(IPLON,LAY)-13)*5+(JT (IPLON,LAY)-1))*NSPB(14) + 1
      IND1(IPLON) = ((JP(IPLON,LAY)-12)*5+(JT1(IPLON,LAY)-1))*NSPB(14) + 1
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG14
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=2
!OCL  UNROLL(2)
      DO IG = 1, NG14
#endif
        TAU (IPLON,NGS13+IG,LAY) = COLCO2(IPLON,LAY) *     &
             (FAC00(IPLON,LAY) * ABSB_14(IND0(IPLON)  ,IG) +  &
              FAC10(IPLON,LAY) * ABSB_14(IND0(IPLON)+1,IG) +  &
              FAC01(IPLON,LAY) * ABSB_14(IND1(IPLON)  ,IG) +  &
              FAC11(IPLON,LAY) * ABSB_14(IND1(IPLON)+1,IG))   &
          + TAUAERL(IPLON,LAY,14)
        PFRAC(IPLON,NGS13+IG,LAY) = FRACREFB_14(IG)
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL14
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL15 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCO2,COLN2O,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 15:  2380-2600 cm-1 (low - N2O,CO2; high - nothing)

! Modifications
!
!     D Salmond 1999-07-14 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCO2(KBDIM,KLEV)
REAL(DP):: COLN2O(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLN2O(IPLON,LAY) + STRRAT_15*COLCO2(IPLON,LAY)
      SPECPARM = COLN2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(SPECPARM,ONEMINUS)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+(JT (IPLON,LAY)-1))*NSPA(15) + JS
      IND1(IPLON) =  (JP(IPLON,LAY)*5   +(JT1(IPLON,LAY)-1))*NSPA(15) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
    DO IG = 1, NG15
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=2
!OCL  UNROLL(2)
      DO IG = 1, NG15
#endif
        TAU (IPLON,NGS14+IG,LAY) = SPECCOMB(IPLON) *                    &
             ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSA_15(IND0(IPLON),IG) +      &
                         FAC10(IPLON,LAY) * ABSA_15(IND0(IPLON)+9,IG) +    &
                         FAC01(IPLON,LAY) * ABSA_15(IND1(IPLON),IG) +      &
                         FAC11(IPLON,LAY) * ABSA_15(IND1(IPLON)+9,IG)) +   &
                  FS *  (FAC01(IPLON,LAY) * ABSA_15(IND1(IPLON)+1,IG) +    &
                         FAC10(IPLON,LAY) * ABSA_15(IND0(IPLON)+10,IG) +   &
                         FAC00(IPLON,LAY) * ABSA_15(IND0(IPLON)+1,IG) +    &
                         FAC11(IPLON,LAY) * ABSA_15(IND1(IPLON)+10,IG))) + &
                        COLH2O(IPLON,LAY) *                             &
                       SELFFAC(IPLON,LAY) * (SELFREF_15(INDS(IPLON),IG) +  &
                      SELFFRAC(IPLON,LAY) *                             &
                     (SELFREF_15(INDS(IPLON)+1,IG) - SELFREF_15(INDS(IPLON),IG)))    &
                    + TAUAERL(IPLON,LAY,15)
        PFRAC(IPLON,NGS14+IG,LAY) = FRACREFA_15(IG,JS) + FS * &
             (FRACREFA_15(IG,JS+1) - FRACREFA_15(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
    DO IG = 1, NG15
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=2
!OCL  UNROLL(2)
      DO IG = 1, NG15
#endif
        TAU  (IPLON,NGS14+IG,LAY) = TAUAERL(IPLON,LAY,15)
        PFRAC(IPLON,NGS14+IG,LAY) = 0.0_dp
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL15
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
SUBROUTINE rad_lon_RRTM_TAUMOL16 (KPROMA,KBDIM,KLEV,IXC,IXLOW,IXHIGH,TAU,&
  &TAUAERL,FAC00,FAC01,FAC10,FAC11,JP,JT,JT1,ONEMINUS,&
  &COLH2O,COLCH4,SELFFAC,SELFFRAC,INDSELF,PFRAC)

!     BAND 16:  2600-3000 cm-1 (low - H2O,CH4; high - nothing)

! Modifications
!
!     D Salmond 1999-07-14 speed-up

IMPLICIT NONE
INTRINSIC :: INT, MIN, MOD

!     DUMMY INTEGER SCALARS
INTEGER :: KPROMA, KBDIM, KLEV
!     DUMMY INTEGER ARRAYS
INTEGER :: IXC(KLEV), IXLOW(KBDIM,KLEV), IXHIGH(KBDIM,KLEV)

!  Output
REAL(DP):: TAU(KBDIM,JPGPT,KLEV)

!- from AER
REAL(DP):: TAUAERL(KBDIM,KLEV,JPBAND)

!- from INTFAC      
REAL(DP):: FAC00(KBDIM,KLEV)
REAL(DP):: FAC01(KBDIM,KLEV)
REAL(DP):: FAC10(KBDIM,KLEV)
REAL(DP):: FAC11(KBDIM,KLEV)

!- from INTIND
INTEGER :: JP(KBDIM,KLEV)
INTEGER :: JT(KBDIM,KLEV)
INTEGER :: JT1(KBDIM,KLEV)

!- from PRECISE             
REAL(DP):: ONEMINUS

!- from PROFDATA             
REAL(DP):: COLH2O(KBDIM,KLEV)
REAL(DP):: COLCH4(KBDIM,KLEV)

!- from SELF             
REAL(DP):: SELFFAC(KBDIM,KLEV)
REAL(DP):: SELFFRAC(KBDIM,KLEV)
INTEGER :: INDSELF(KBDIM,KLEV)

!- from SP             
REAL(DP):: PFRAC(KBDIM,JPGPT,KLEV)


INTEGER :: IJS(KBDIM)
INTEGER :: IND0(KBDIM), IND1(KBDIM), INDS(KBDIM)

!     LOCAL INTEGER SCALARS
INTEGER :: IG, JS, LAY, IPLON
INTEGER :: IXP, IXC0
REAL(DP):: ZFS(KBDIM), SPECCOMB(KBDIM)

!     LOCAL REAL SCALARS
REAL(DP):: FS, SPECMULT, SPECPARM


!     Compute the optical depth by interpolating in ln(pressure), 
!     temperature, and appropriate species.  Below LAYTROP, the water
!     vapor self-continuum is interpolated (in temperature) separately. 
 
  DO LAY = 1, KLEV
    IXC0 = IXC(LAY)

    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      SPECCOMB(IPLON) = COLH2O(IPLON,LAY) + STRRAT_16*COLCH4(IPLON,LAY)
      SPECPARM = COLH2O(IPLON,LAY) / SPECCOMB(IPLON)
      SPECPARM = MIN(SPECPARM,ONEMINUS)
      SPECMULT = 8._dp*(SPECPARM)
      JS = 1 + INT(SPECMULT)
      FS = MOD(SPECMULT,1.0_dp)
      IND0(IPLON) = ((JP(IPLON,LAY)-1)*5+( JT(IPLON,LAY)-1))*NSPA(16) + JS
      IND1(IPLON) = (JP(IPLON,LAY)*5    +(JT1(IPLON,LAY)-1))*NSPA(16) + JS
      INDS(IPLON) = INDSELF(IPLON,LAY)
      ZFS(IPLON) = FS
      IJS(IPLON) = JS
    ENDDO

#ifndef UNROLLNGX
!    DO IG = 1, NG16
    DO IG = 1, 2
#endif
    DO IXP = 1, IXC0
      IPLON = IXLOW(IXP,LAY)
      FS = ZFS(IPLON)
      JS = IJS(IPLON)
#ifdef UNROLLNGX
!CDIR UNROLL=2
!OCL  UNROLL(2)
!      DO IG = 1, NG16
      DO IG = 1, 2
#endif
        TAU (IPLON,NGS15+IG,LAY) = SPECCOMB(IPLON) *                  &
           ((1._DP - FS)*(FAC00(IPLON,LAY) * ABSA_16(IND0(IPLON),IG) +      &
                       FAC10(IPLON,LAY) * ABSA_16(IND0(IPLON)+9,IG) +    &
                       FAC01(IPLON,LAY) * ABSA_16(IND1(IPLON),IG) +      &
                       FAC11(IPLON,LAY) * ABSA_16(IND1(IPLON)+9,IG)) +   &
               FS * (  FAC01(IPLON,LAY) * ABSA_16(IND1(IPLON)+1,IG) +    &
                       FAC10(IPLON,LAY) * ABSA_16(IND0(IPLON)+10,IG) +   &
                       FAC00(IPLON,LAY) * ABSA_16(IND0(IPLON)+1,IG) +    &
                       FAC11(IPLON,LAY) * ABSA_16(IND1(IPLON)+10,IG))) + &
                      COLH2O(IPLON,LAY) *                             &
                     SELFFAC(IPLON,LAY) * (SELFREF_16(INDS(IPLON),IG) +  &
                    SELFFRAC(IPLON,LAY) *                             &
                    (SELFREF_16(INDS(IPLON)+1,IG) - SELFREF_16(INDS(IPLON),IG)))   &
                  + TAUAERL(IPLON,LAY,16)
        PFRAC(IPLON,NGS15+IG,LAY) = FRACREFA_16(IG,JS) + FS * &
             (FRACREFA_16(IG,JS+1) - FRACREFA_16(IG,JS))
      ENDDO
    ENDDO

    IXC0 = KPROMA - IXC0
#ifndef UNROLLNGX
!    DO IG = 1, NG16
    DO IG = 1, 2
#endif
    DO IXP = 1, IXC0
      IPLON = IXHIGH(IXP,LAY)
#ifdef UNROLLNGX
!CDIR UNROLL=2
!OCL  UNROLL(2)
!      DO IG = 1, NG16
      DO IG = 1, 2
#endif
        TAU  (IPLON,NGS15+IG,LAY) = TAUAERL(IPLON,LAY,16)
        PFRAC(IPLON,NGS15+IG,LAY) = 0.0_dp
      ENDDO
    ENDDO

  ENDDO

  RETURN
END SUBROUTINE rad_lon_RRTM_TAUMOL16
!----------------------------------------------------------------------------

!***************************************************************************
END MODULE messy_rad_long
!***************************************************************************
