  �  E   k820309    ?          19.1        ��ja                                                                                                          
       messy_contrail.f90 MESSY_CONTRAIL              DP CONTRAIL_READ_NML_CTRL CONTRAIL_POT_COV CONTRAIL_CALC CONTRAIL_CALC_DEV CONTRAIL_UV_GRAD MODSTR MODVER EI_H2O Q_FUEL ETA_AC A_SAC R_SAC                                                     
       DP                                                                                                                                                                                  	                       Ccontrail                                                                                                                   C1.1                          D@                                    
                 D@                                    
                 D@                                    
                 D@                                    
                 D@                               	     
       #         @                                   
                    #STATUS    #IOU              D                                                       
  @                                          #         H                                                      #T    #P    #Q    #ZQTE    #B_CI    #ZRHC    #DT    #B_CC    #POTCOV    #ZQSM1    #ZCONPN                                                                                                                                                                                          
  @                                   
                
                                      
                
  @                                   
                
                                      
                
                                      
                
                                      
                
                                      
                D @                                   
                 D @                                   
                 D @                                   
                 D                                     
       #         H                                                  
   #Q    #B_CI    #ZCONRAT    #POTCOV    #ZCONPN    #EMIS    #SCAL     #GBOXAREA !   #PCONCOV "   #PCONIWC #             
                                      
                
                                      
                
                                      
                
  @                                   
                
                                      
                
                                      
                
                                       
                
  @                              !     
                D @                              "     
                 D @                              #     
       #         H                                  $                    #Q %   #B_CI &   #ZCONRAT '   #POTCOV (   #GBOXAREA )   #ZCONPN *   #DT +   #EMIS ,   #SCAL -   #POTCOV_M1 .   #PCONCOV_M1 /   #PCONIWC_M1 0   #RHO 1   #ETA_DOT 2   #DUDZ 3   #DVDZ 4   #PCONCOV_SUM 5   #PCONIWC_SUM 6   #PCONCOV_NOW 7   #PCONIWC_NOW 8   #CONCOV_TE_SPREAD 9   #CONIWC_TE_SEDI :   #CONIWC_TE_POT ;             
  @                              %     
                
  @                              &     
                
  @                              '     
                
  @                              (     
                
  @                              )     
                
  @                              *     
                
                                 +     
                
  @                              ,     
                
  @                              -     
                
                                 .     
                
  @                              /     
                
  @                              0     
                
  @                              1     
                
  @                              2     
                
  @                              3     
                
  @                              4     
                D @                              5     
                 D @                              6     
                 D @                              7     
                 D @                              8     
                 D @                              9     
                 D @                              :     
                 D                                ;     
       #         @                                   <                    #U_SCB =   #V_SCB >   #SQCST_2D ?   #RHO @   #ZPF A   #DUDZ B   #DVDZ C                                     
 @                              =                   
              &                   &                                                     
                                 >                   
              &                   &                                                     
                                 ?                   
              &                                                     
                                 @                   
              &                   &                                                     
                                 A                   
              &                   &                                                     
D                                B                   
               &                   &                                                     
D                                C                   
               &                   &                                              �   *      fn#fn $   �   �   b   uapp(MESSY_CONTRAIL )   e  C   J  MESSY_MAIN_CONSTANTS_MEM ,   �  p       DP+MESSY_MAIN_CONSTANTS_MEM      �       MODSTR    �  �       MODVER    %  @       EI_H2O    e  @       Q_FUEL    �  @       ETA_AC    �  @       A_SAC    %  @       R_SAC '   e  ]       CONTRAIL_READ_NML_CTRL .   �  @   a   CONTRAIL_READ_NML_CTRL%STATUS +     @   a   CONTRAIL_READ_NML_CTRL%IOU !   B  \      CONTRAIL_POT_COV #   �  @   a   CONTRAIL_POT_COV%T #   �  @   a   CONTRAIL_POT_COV%P #     @   a   CONTRAIL_POT_COV%Q &   ^  @   a   CONTRAIL_POT_COV%ZQTE &   �  @   a   CONTRAIL_POT_COV%B_CI &   �  @   a   CONTRAIL_POT_COV%ZRHC $     @   a   CONTRAIL_POT_COV%DT &   ^  @   a   CONTRAIL_POT_COV%B_CC (   �  @   a   CONTRAIL_POT_COV%POTCOV '   �  @   a   CONTRAIL_POT_COV%ZQSM1 (   	  @   a   CONTRAIL_POT_COV%ZCONPN    ^	  �       CONTRAIL_CALC     
  @   a   CONTRAIL_CALC%Q #   X
  @   a   CONTRAIL_CALC%B_CI &   �
  @   a   CONTRAIL_CALC%ZCONRAT %   �
  @   a   CONTRAIL_CALC%POTCOV %     @   a   CONTRAIL_CALC%ZCONPN #   X  @   a   CONTRAIL_CALC%EMIS #   �  @   a   CONTRAIL_CALC%SCAL '   �  @   a   CONTRAIL_CALC%GBOXAREA &     @   a   CONTRAIL_CALC%PCONCOV &   X  @   a   CONTRAIL_CALC%PCONIWC "   �  �      CONTRAIL_CALC_DEV $     @   a   CONTRAIL_CALC_DEV%Q '   Z  @   a   CONTRAIL_CALC_DEV%B_CI *   �  @   a   CONTRAIL_CALC_DEV%ZCONRAT )   �  @   a   CONTRAIL_CALC_DEV%POTCOV +     @   a   CONTRAIL_CALC_DEV%GBOXAREA )   Z  @   a   CONTRAIL_CALC_DEV%ZCONPN %   �  @   a   CONTRAIL_CALC_DEV%DT '   �  @   a   CONTRAIL_CALC_DEV%EMIS '     @   a   CONTRAIL_CALC_DEV%SCAL ,   Z  @   a   CONTRAIL_CALC_DEV%POTCOV_M1 -   �  @   a   CONTRAIL_CALC_DEV%PCONCOV_M1 -   �  @   a   CONTRAIL_CALC_DEV%PCONIWC_M1 &     @   a   CONTRAIL_CALC_DEV%RHO *   Z  @   a   CONTRAIL_CALC_DEV%ETA_DOT '   �  @   a   CONTRAIL_CALC_DEV%DUDZ '   �  @   a   CONTRAIL_CALC_DEV%DVDZ .     @   a   CONTRAIL_CALC_DEV%PCONCOV_SUM .   Z  @   a   CONTRAIL_CALC_DEV%PCONIWC_SUM .   �  @   a   CONTRAIL_CALC_DEV%PCONCOV_NOW .   �  @   a   CONTRAIL_CALC_DEV%PCONIWC_NOW 3     @   a   CONTRAIL_CALC_DEV%CONCOV_TE_SPREAD 1   Z  @   a   CONTRAIL_CALC_DEV%CONIWC_TE_SEDI 0   �  @   a   CONTRAIL_CALC_DEV%CONIWC_TE_POT !   �  �       CONTRAIL_UV_GRAD '   �  �   a   CONTRAIL_UV_GRAD%U_SCB '   (  �   a   CONTRAIL_UV_GRAD%V_SCB *   �  �   a   CONTRAIL_UV_GRAD%SQCST_2D %   X  �   a   CONTRAIL_UV_GRAD%RHO %   �  �   a   CONTRAIL_UV_GRAD%ZPF &   �  �   a   CONTRAIL_UV_GRAD%DUDZ &   D  �   a   CONTRAIL_UV_GRAD%DVDZ 