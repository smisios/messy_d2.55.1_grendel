  ?  =   k820309    ?          19.1        ??ja                                                                                                          
       messy_airtraf_gc.f90 MESSY_AIRTRAF_GC              GC_WAYPOINTS GC_AC_TRAJ WIND_EFFECT DIV_WAY D_LON D_LAT D_TIME A_LON A_LAT ADIABATIC_INDEX_AIR GAS_CONSTANT_AIR ISWITCH OUT_SWITCH VTAS_SWITCH WIND_SWITCH RARAD1 RARAD2 DCRAD1 DCRAD2 DELDC2 DELRA DELRA2 SINDIS SINDIS1 SINDIS2 THETAR PHIR GC_D_TIME FL_DIR                                                     
       DP PI DTR ONEDAY RADIUS_EARTH G                      @                              
       NN_INDEX                                                                                                                                                        
                
                 -DT?!	@        3.14159265358979323846                                                 
                   
                  9?R?Fߑ?                                                         
                
                      ?@        86400.0                                                 
                
          	           ?MXA        6371000.0                                                 
                
                 ??:?#@        9.80665#         @                                  	                    #LIST 
   #VALUE    #IDX    #IDX2    #FAC    #STATUS              
                                
                   
              &                                                     
                                      
                                                                                                                              
                                     
                                                                                                              
                
                      ??@        10668.                                                 
                
                 =
ףp=??        0.82                                                                                                   1                                                                                                   2                                                                                                   3                                                                                                   4                                                                                                   5                                                                                                   6                                                                                                   7                                                                                                   8                                                                                    	               9                                                                                    
               10                                                                                                   11                                                                                                   12                                                                                                   13                                                                                                   14                                                                                                    15                                             !                                                      16#         @                                   "                    #GC_NWAYPOINTS #   #GC_AC_ROUTES $   #GC_P2_AC_ROUTES_DESC %   #GC_PHILON &   #GC_PHILAT '   #GC_ZGL_GEOPOT_3D (   #GC_UWIND_G )   #GC_VWIND_G *   #GC_V_Z_G +   #GC_T_SCB_G ,   #GC_CPC_G -   #GC_FL_TIME .   #GC_BADA_FUEL_NOM /   #GC_VTAS_ALONG_TRAJ 0   #GC_FL_DIRECTION 1   #GC_FL_DISTANCE 2   #GC_CPC 3             
                                  #                     D                                $                   
               &                   &                                                     
                                 %                   
              &                                                     
                                 &                   
              &                                                     
                                 '                   
              &                                                     
 @                              (                   
              &                   &                   &                                                     
                                 )                   
              &                   &                   &                                                     
                                 *                   
              &                   &                   &                                                     
                                 +                   
              &                   &                   &                                                     
                                 ,                   
 	             &                   &                   &                                                     
                                 -                   
 
             &                   &                   &                                                     D                                .     
                 D                                /     
                
D @                              0                    
     p          5 ? p        r #       5 ? p        r #                               D                                 1                      D                                 2                      D @                              3     
       #         @                                   4                    #COSIS 5   #T_DIS 6   #T_DISN 7             D                                5     
                 D                                6     
                 D                                7     
       #         @                                   8                    #FLIGHT_TIME_TABLE 9   #FUEL_U ;            
                                 9                    
 '   p           5 r :   n                                       1     5 r :   n                                      1                                    D                                ;                    
 (    p           5 r :   n                                       1     5 r :   n                                      1                                      @  @                              :               ?   .      fn#fn &   ?     b   uapp(MESSY_AIRTRAF_GC )   ?  `   J  MESSY_MAIN_CONSTANTS_MEM !   =  I   J  MESSY_MAIN_TOOLS ,   ?  p       DP+MESSY_MAIN_CONSTANTS_MEM ,   ?  ?       PI+MESSY_MAIN_CONSTANTS_MEM -   |  p       DTR+MESSY_MAIN_CONSTANTS_MEM 0   ?  w       ONEDAY+MESSY_MAIN_CONSTANTS_MEM 6   c  y       RADIUS_EARTH+MESSY_MAIN_CONSTANTS_MEM +   ?  w       G+MESSY_MAIN_CONSTANTS_MEM *   S  ?       NN_INDEX+MESSY_MAIN_TOOLS /   ?  ?   a   NN_INDEX%LIST+MESSY_MAIN_TOOLS 0   d  @   a   NN_INDEX%VALUE+MESSY_MAIN_TOOLS .   ?  @   a   NN_INDEX%IDX+MESSY_MAIN_TOOLS /   ?  @   a   NN_INDEX%IDX2+MESSY_MAIN_TOOLS .   $  @   a   NN_INDEX%FAC+MESSY_MAIN_TOOLS 1   d  @   a   NN_INDEX%STATUS+MESSY_MAIN_TOOLS    ?  v       ALT      t       AC_MACH    ?  q       PROPS_LON    ?  q       PROPS_LAT    p	  q       PROPS_ALT    ?	  q       PROPS_TIME    R
  q       PROPS_AC_SPEED    ?
  q       PROPS_DIST    4  q       PROPS_FUEL_USE    ?  q       PROPS_EMIS_NOX      q       PROPS_EMIS_H2O    ?  r       PROPS_POTCOV    ?  r       PROPS_ATR20O3    k  r       PROPS_ATR20CH4    ?  r       PROPS_ATR20H2O    O  r       PROPS_ATR20CPC    ?  r       PROPS_ATR20CO2    3  r       PROPS_ATR20TOT    ?  z      GC_CALCULATE +     @   a   GC_CALCULATE%GC_NWAYPOINTS *   _  ?   a   GC_CALCULATE%GC_AC_ROUTES 2     ?   a   GC_CALCULATE%GC_P2_AC_ROUTES_DESC '   ?  ?   a   GC_CALCULATE%GC_PHILON '     ?   a   GC_CALCULATE%GC_PHILAT .   ?  ?   a   GC_CALCULATE%GC_ZGL_GEOPOT_3D (   c  ?   a   GC_CALCULATE%GC_UWIND_G (     ?   a   GC_CALCULATE%GC_VWIND_G &   ?  ?   a   GC_CALCULATE%GC_V_Z_G (   ?  ?   a   GC_CALCULATE%GC_T_SCB_G &   S  ?   a   GC_CALCULATE%GC_CPC_G (     @   a   GC_CALCULATE%GC_FL_TIME .   O  @   a   GC_CALCULATE%GC_BADA_FUEL_NOM 0   ?  ?   a   GC_CALCULATE%GC_VTAS_ALONG_TRAJ -   C  @   a   GC_CALCULATE%GC_FL_DIRECTION ,   ?  @   a   GC_CALCULATE%GC_FL_DISTANCE $   ?  @   a   GC_CALCULATE%GC_CPC      j       GC_DISTANCE "   m  @   a   GC_DISTANCE%COSIS "   ?  @   a   GC_DISTANCE%T_DIS #   ?  @   a   GC_DISTANCE%T_DISN /   -  k       TOTAL_ENERGY_MODEL_CALCULATION A   ?    a   TOTAL_ENERGY_MODEL_CALCULATION%FLIGHT_TIME_TABLE 6   ?    a   TOTAL_ENERGY_MODEL_CALCULATION%FUEL_U    ?  @      DIV_WAY 