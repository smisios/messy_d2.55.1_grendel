  �J  �   k820309    ?          19.1        ��ja                                                                                                          
       messy_attila.f90 MESSY_ATTILA       M       NCELL NGCELL AMCELL DP ATTILA_READ_NML_CTRL ATTILA_MESSAGE ATTILA_GLOBAL_INIT ATTILA_INICOM_1 ATTILA_INICOM_2 ATTILA_ALLOC ATTILA_INIPOS ATTILA_REALLOC ATTILA_COUNT_CELLS ATTILA_RESETPOS_TRAJEC ATTILA_FIRSTL ATTILA_FIRSTL_EXT ATTILA_FILL ATTILA_DRIVE ATTILA_CONVGPC ATTILA_CONVTRAJ ATTILA_SAVE_POS ATTILA_GLOBAL_EXIT MODSTR MODVER KAPPA NR1 NR2 NCHUNK CPGBAVE LLTINFO I_PBLH_METHOD ADICO LLTBLTURB LLCONV LLCAT L_LG_FMM LVDIAG I_NCELL LTRAJEC LTRAJEC_DATE LTRAJEC_SAME_DATE I_VERT PRESS_REF NGL NLON NLEV AMTOT LLTMOCA LINIPOS APOS NPOS NPOSM1 APOS_FMM PWU PWV PWW PPH SURFTEMP PPS PTPOT G_TPOT_H PTEMP SIGMA_REF_G I_REF_G PMBOX PMBOXBL KHPBL NCB NCBM1 NCBL GNCB GNCBL GNCB_MOVE UVEL VVEL WVEL APOSM1_PRESS                      @                              
       NCELL NGCELL AMCELL                                                     
      CPD CP_AIR A RADIUS_EARTH DP RD G                                                                                                                @                                                    @@                                                     @                                    
       #         @                                                       #STATUS    #IOU 	             D                                                       
  @                               	           #         @                                   
                     #         @                                                       #NGL_EX    #NLON_EX    #NLEV_EX    #DTIME_EX    #L_PM_EX    #VCT_EX    #NVCLEV_EX    #APSURF_EX    #APZERO_EX    #GL_GMU_EX    #NN_EX    #NPLVP1_EX    #NLMSGL_EX    #NPLVP2_EX    #NLMSLP_EX    #STATUS              
                                                       
                                                       
                                                       
                                      
                
                                                       
 @                                                 
 A             &                                                     
                                                       
                                      
                
                                      
                
 @                                                 
 B             &                                                     
                                                       
                                                       
                                                       
                                                       
                                                       D @                                           #         @                                                        #         @                                                       #STATUS              D                                             #         @                                                       #STATUS     #ZNCELL !   #LSPAC "   #LSPACINT #   #LSPACINTM1 $   #LSPACPRESSM1 %             D                                                         @                               !                      P                              "                   
 C              &                   &                   &                   &                   &                                                     P                              #                   
 D              &                   &                   &                   &                   &                                                     P                              $                   
 E              &                   &                   &                   &                   &                                                     P                              %                   
 F              &                   &                   &                   &                   &                                           #         @                                   &                    #HARVESTINIPOS1 '   #HARVESTINIPOS2 (   #STATUS )             
 @                              '                   
 G             &                   &                                                     
 @                              (                   
 H             &                   &                                                     D                                 )            #         @                                   *                    #FLAG +   #TMP_APOS ,   #TMP_NPOS -   #TMP_NPOSM1 .   #LSPAC /   #LSPACINT 0   #LSPACINTM1 1             
                                  +                    D                               ,                   
 O              &                   &                                                    D                               -                   
 P              &                   &                                                    D                               .                   
 Q              &                   &                                                     P                              /                   
 R              &                   &                   &                   &                   &                                                     P                              0                   
 S              &                   &                   &                   &                   &                                                     P                              1                   
 T              &                   &                   &                   &                   &                                           #         @                                   2                    #RNCB 3   #RNCBL 7   #L_INIT 8            D                                3                    
 U        p        5 r 4   p        5 r 5   p          5 r 5     5 r 4     5 r 6       5 r 5     5 r 4     5 r 6                              D                                7                    
 V      p        5 r 5   p          5 r 5     5 r 6       5 r 5     5 r 6                               
 @                               8           #         @                                   9                    #ICELL :   #LAT ;   #LON <   #PRESS =             
  @                               :                     
                                 ;     
                
                                 <     
                
  @                              =     
      #         @                                   >                    #PTSM1M ?   #PGEOM1 @   #ZTETA1 A   #ZRI B   #TPOT C             
                                 ?                   
 _             &                   &                                                     
                                 @                   
 `             &                   &                   &                                                     
                                 A                   
 a             &                   &                   &                                                     
                                 B                   
 b             &                   &                   &                                                     
                                 C                   
 c             &                   &                   &                                           #         @                                   D                    #STATUS E   #PBLH_IDX F             D                                 E                      
 @                              F                   
 h             &                   &                                           #         @                                   G                    #STATUS H             D @                               H            #         @                                   I                    #HARVESTTURB J   #HARVESTMC K   #HARVESTCAT L   #TI2 M   #LSTART N             
                                 J                   
 i             &                                                     
                                 K                   
 j             &                   &                                                     
@ @                              L                   
 k             &                                                    D P                              M                   
 l              &                   &                   &                                                     
                                  N           #         @                                   O                    #UMASSFL P   #DMASSFL Q   #TYPECONV R   #STATUS S             
                                 P                   
 r             &                   &                   &                                                     
                                 Q                   
 s             &                   &                   &                                                     
                                 R                   
 t             &                   &                                                     D                                 S            #         @                                   T                    #IUPDO U   #JTOFFSET V   #HARVESTCONV W   #ZPOS_BEG_UP X   #ZPOS_END_UP Y   #ZPOS_END_DO Z   #ZPOS_END_SUB [   #L_EXIT \             
                                  U                     
                                  V                     
                                 W                   
 u             &                                                     
D                                X                   
 v              &                                                     
D                                Y                   
 w              &                                                     
D                                Z                   
 x              &                                                     
D                                [                   
 y              &                                                     D                                 \            #         @                                   ]                    #RPOS ^   #N _             D                                ^                   
 �              &                                                     
                                  _           #         @                                   `                                                                a                                                        Cattila                                                           b                                                        C4.0a                                                            c     
                   
                  ��WöD�?                                                     d                                                      4                                             e                                                      2          @@@                               f                      @@                               g     
                 @@                                h                      @@                                i                      @@                               j                   
      p          p            p                                    @@                                k                      @@                                l                      @@                                m                                                        n                      @@                                o                      @@                                p                      @@                                q                      @@                                r                         p          p            p                                    @@                                s                      D@                                t                      D@                               u     
                 @                                6                      @@                               5                      @@                               4                       @                               v     
                  @                                w                       @                                x                     @P                              y                   
                &                   &                                                    @                              z                   
                &                   &                                                    @                              {                   
                &                   &                                                    @                              |                   
                &                   &                   &                                                     P                              }                   
                &                   &                   &                                                     P                              ~                   
                &                   &                   &                                                    @P                                                 
                &                   &                   &                                                                                   �                   
                &                   &                   &                                                                                   �                   
                &                   &                                                    @P                              �                   
                &                   &                                                     P                              �                   
                &                   &                   &                                                    @                              �                   
                &                   &                   &                                                     P                              �                   
                &                   &                   &                                                                                   �                   
                &                   &                                                                                   �                   
                &                   &                                                    @@                              �                   
                &                   &                   &                                                    @                               �                   
                &                   &                                                    @                                �                                   &                   &                                                    @                                �                                   &                   &                   &                                                    @                                �                                   &                   &                   &                                                    @                                �                                   &                   &                                                    @P                              �                   
                &                   &                   &                                                    @P                              �                   
                &                   &                                                    @P                              �                   
                &                   &                   &                                                    @                              �                   
                &                                                    @                              �                   
                &                                                    @                              �                   
                &                                                    @                              �                   
                &                                              �   &      fn#fn "   �   �  b   uapp(MESSY_ATTILA !   �  T   j  MESSY_ATTILA_MEM )   �  b   J  MESSY_MAIN_CONSTANTS_MEM ,   O  p       DP+MESSY_MAIN_CONSTANTS_MEM '   �  @       NCELL+MESSY_ATTILA_MEM (   �  @       NGCELL+MESSY_ATTILA_MEM (   ?  @       AMCELL+MESSY_ATTILA_MEM %     ]       ATTILA_READ_NML_CTRL ,   �  @   a   ATTILA_READ_NML_CTRL%STATUS )     @   a   ATTILA_READ_NML_CTRL%IOU    \  H       ATTILA_MESSAGE #   �  $      ATTILA_GLOBAL_INIT *   �  @   a   ATTILA_GLOBAL_INIT%NGL_EX +     @   a   ATTILA_GLOBAL_INIT%NLON_EX +   H  @   a   ATTILA_GLOBAL_INIT%NLEV_EX ,   �  @   a   ATTILA_GLOBAL_INIT%DTIME_EX +   �  @   a   ATTILA_GLOBAL_INIT%L_PM_EX *   	  �   a   ATTILA_GLOBAL_INIT%VCT_EX -   �	  @   a   ATTILA_GLOBAL_INIT%NVCLEV_EX -   �	  @   a   ATTILA_GLOBAL_INIT%APSURF_EX -   
  @   a   ATTILA_GLOBAL_INIT%APZERO_EX -   T
  �   a   ATTILA_GLOBAL_INIT%GL_GMU_EX )   �
  @   a   ATTILA_GLOBAL_INIT%NN_EX -      @   a   ATTILA_GLOBAL_INIT%NPLVP1_EX -   `  @   a   ATTILA_GLOBAL_INIT%NLMSGL_EX -   �  @   a   ATTILA_GLOBAL_INIT%NPLVP2_EX -   �  @   a   ATTILA_GLOBAL_INIT%NLMSLP_EX *      @   a   ATTILA_GLOBAL_INIT%STATUS     `  H       ATTILA_INICOM_1     �  T       ATTILA_INICOM_2 '   �  @   a   ATTILA_INICOM_2%STATUS    <  �       ATTILA_ALLOC $   �  @   a   ATTILA_ALLOC%STATUS $     @   a   ATTILA_ALLOC%ZNCELL #   W  �   a   ATTILA_ALLOC%LSPAC &   C  �   a   ATTILA_ALLOC%LSPACINT (   /  �   a   ATTILA_ALLOC%LSPACINTM1 *     �   a   ATTILA_ALLOC%LSPACPRESSM1      |       ATTILA_INIPOS -   �  �   a   ATTILA_INIPOS%HARVESTINIPOS1 -   '  �   a   ATTILA_INIPOS%HARVESTINIPOS2 %   �  @   a   ATTILA_INIPOS%STATUS      �       ATTILA_REALLOC $   �  @   a   ATTILA_REALLOC%FLAG (   �  �   a   ATTILA_REALLOC%TMP_APOS (   �  �   a   ATTILA_REALLOC%TMP_NPOS *   :  �   a   ATTILA_REALLOC%TMP_NPOSM1 %   �  �   a   ATTILA_REALLOC%LSPAC (   �  �   a   ATTILA_REALLOC%LSPACINT *   �  �   a   ATTILA_REALLOC%LSPACINTM1 #   �  i       ATTILA_COUNT_CELLS (       a   ATTILA_COUNT_CELLS%RNCB )     �   a   ATTILA_COUNT_CELLS%RNCBL *   �  @   a   ATTILA_COUNT_CELLS%L_INIT '   3  p       ATTILA_RESETPOS_TRAJEC -   �  @   a   ATTILA_RESETPOS_TRAJEC%ICELL +   �  @   a   ATTILA_RESETPOS_TRAJEC%LAT +   #  @   a   ATTILA_RESETPOS_TRAJEC%LON -   c  @   a   ATTILA_RESETPOS_TRAJEC%PRESS    �         ATTILA_FIRSTL %   "  �   a   ATTILA_FIRSTL%PTSM1M %   �  �   a   ATTILA_FIRSTL%PGEOM1 %   �  �   a   ATTILA_FIRSTL%ZTETA1 "   >   �   a   ATTILA_FIRSTL%ZRI #   �   �   a   ATTILA_FIRSTL%TPOT "   �!  b       ATTILA_FIRSTL_EXT )   "  @   a   ATTILA_FIRSTL_EXT%STATUS +   X"  �   a   ATTILA_FIRSTL_EXT%PBLH_IDX    �"  T       ATTILA_FILL #   P#  @   a   ATTILA_FILL%STATUS    �#  �       ATTILA_DRIVE )   $  �   a   ATTILA_DRIVE%HARVESTTURB '   �$  �   a   ATTILA_DRIVE%HARVESTMC (   M%  �   a   ATTILA_DRIVE%HARVESTCAT !   �%  �   a   ATTILA_DRIVE%TI2 $   �&  @   a   ATTILA_DRIVE%LSTART    �&  |       ATTILA_CONVGPC '   Q'  �   a   ATTILA_CONVGPC%UMASSFL '   (  �   a   ATTILA_CONVGPC%DMASSFL (   �(  �   a   ATTILA_CONVGPC%TYPECONV &   m)  @   a   ATTILA_CONVGPC%STATUS     �)  �       ATTILA_CONVTRAJ &   p*  @   a   ATTILA_CONVTRAJ%IUPDO )   �*  @   a   ATTILA_CONVTRAJ%JTOFFSET ,   �*  �   a   ATTILA_CONVTRAJ%HARVESTCONV ,   |+  �   a   ATTILA_CONVTRAJ%ZPOS_BEG_UP ,   ,  �   a   ATTILA_CONVTRAJ%ZPOS_END_UP ,   �,  �   a   ATTILA_CONVTRAJ%ZPOS_END_DO -    -  �   a   ATTILA_CONVTRAJ%ZPOS_END_SUB '   �-  @   a   ATTILA_CONVTRAJ%L_EXIT     �-  Y       ATTILA_SAVE_POS %   E.  �   a   ATTILA_SAVE_POS%RPOS "   �.  @   a   ATTILA_SAVE_POS%N #   /  H       ATTILA_GLOBAL_EXIT    Y/  �       MODSTR    �/  �       MODVER    e0  p       KAPPA    �0  q       NR1    F1  q       NR2    �1  @       NCHUNK    �1  @       CPGBAVE    72  @       LLTINFO    w2  @       I_PBLH_METHOD    �2  �       ADICO    K3  @       LLTBLTURB    �3  @       LLCONV    �3  @       LLCAT    4  @       L_LG_FMM    K4  @       LVDIAG    �4  @       I_NCELL    �4  @       LTRAJEC    5  �       LTRAJEC_DATE "   �5  @       LTRAJEC_SAME_DATE    �5  @       I_VERT    6  @       PRESS_REF    _6  @       NGL    �6  @       NLON    �6  @       NLEV    7  @       AMTOT    _7  @       LLTMOCA    �7  @       LINIPOS    �7  �       APOS    �8  �       NPOS    '9  �       NPOSM1    �9  �       APOS_FMM    �:  �       PWU    C;  �       PWV    �;  �       PWW    �<  �       PPH    w=  �       SURFTEMP    >  �       PPS    �>  �       PTPOT    {?  �       G_TPOT_H    7@  �       PTEMP    �@  �       SIGMA_REF_G    �A  �       I_REF_G    ;B  �       PMBOX    �B  �       PMBOXBL    �C  �       KHPBL    ?D  �       NCB    �D  �       NCBM1    �E  �       NCBL    [F  �       GNCB    G  �       GNCBL    �G  �       GNCB_MOVE    wH  �       UVEL    I  �       VVEL    �I  �       WVEL    J  �       APOSM1_PRESS 