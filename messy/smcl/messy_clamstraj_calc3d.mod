  ÛD  [   k820309    ?          19.1        µÍja                                                                                                          
       messy_clamstraj_calc3d.f90 MESSY_CLAMSTRAJ_CALC3D #         @                                                    !   #PLONS    #PLATS    #PLEVS    #CALC    #TSV    #U0    #V0    #W0    #LEVEL0    #DLEVDZ0    #UMID    #VMID    #WMID    #LEVELMID    #DLEVDZMID    #UDT    #VDT    #WDT    #LEVELDT    #DLEVDZDT    #NX 
   #NY 	   #NZ    #LONGRID    #LATGRID    #LEVELGRID    #NPARTS    #LOOP_START    #LOOP_END    #DT    #PSLAT     #LOGLEV !   #USE_DLEVDZ "                                                                   D @                                                  
     p          5  p        r        5  p        r                               D @                                                  
     p          5  p        r        5  p        r                               D @                                                  
     p          5  p        r        5  p        r                               
                                                          p          5  p        r        5  p        r                                                                                   
     p          5  p        r        5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
 	        p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
 
        p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                @                                                  
         p        5  p        r 	   p        5  p        r 
   p          5  p        r 
     5  p        r 	     5  p        r        5  p        r 
     5  p        r 	     5  p        r                                 @                               
                       @                               	                       @                                                     @                                                  
     p          & p         5  p        r 
   n                                       1     5  p        r 
   n                                      1                                     @                                                  
     p           & p          5  p        r 	   n                                       1       5  p        r 	   n                                      1p         p                                             @                                                  
     p          5  p        r        5  p        r                                D @                                                     D @                                                     D @                                                     
                                      
                
                                       
                
  @                               !                      @                               "            #         @                                  #                    #U $   #V (   #W )   #DLEVDZ *   #LEVEL +   #PLONS ,   #PLATS .   #PLEVS /   #PU 0   #PV 1   #PW 2   #PDLEVDZ 3   #NX &   #NY %   #NZ '   #LONGRID 4   #LATGRID 5   #LEVELGRID 6   #NPARTS -   #LOOP_START 7   #LOOP_END 8   #INDXPS 9   #LOGLEV :                                                                                                               
                                 $                    
 *       p        5  p        r %   p        5  p        r &   p          5  p        r &     5  p        r %     5  p        r '       5  p        r &     5  p        r %     5  p        r '                              
                                 (                    
 +       p        5  p        r %   p        5  p        r &   p          5  p        r &     5  p        r %     5  p        r '       5  p        r &     5  p        r %     5  p        r '                              
                                 )                    
 ,       p        5  p        r %   p        5  p        r &   p          5  p        r &     5  p        r %     5  p        r '       5  p        r &     5  p        r %     5  p        r '                              
                                 *                    
 -       p        5  p        r %   p        5  p        r &   p          5  p        r &     5  p        r %     5  p        r '       5  p        r &     5  p        r %     5  p        r '                              
                                 +                    
 .       p        5  p        r %   p        5  p        r &   p          5  p        r &     5  p        r %     5  p        r '       5  p        r &     5  p        r %     5  p        r '                              D @                              ,                    
 /    p          5  p        r -       5  p        r -                              D @                              .                    
 0    p          5  p        r -       5  p        r -                              
  @                              /                    
 1   p          5  p        r -       5  p        r -                              D                                0                    
 2    p          5  p        r -       5  p        r -                              D                                1                    
 3    p          5  p        r -       5  p        r -                              D                                2                    
 4    p          5  p        r -       5  p        r -                              D                                3                    
 5    p          5  p        r -       5  p        r -                               
                                  &                     
                                  %                     
  @                               '                    
                                 4                    
 6   p          & p         5  p        r &   n                                       1     5  p        r &   n                                      1                                    
                                 5                    
 7   p           & p          5  p        r %   n                                       1       5  p        r %   n                                      1p         p                                            
                                 6                    
 8   p          5  p        r '       5  p        r '                               D @                               -                      D @                               7                      D @                               8                     
                                  9                     )   p          5  p        r -       5  p        r -                               
  @                               :           #         @                                  ;                    #PLONS <   #PLATS >   #INDXPS ?   #NPARTS =   #LOOP_START @   #LOOP_END A   #PI B                                                                                  D                                <                    
 P    p          5  p        r =       5  p        r =                              D                                >                    
 Q    p          5  p        r =       5  p        r =                                                               ?                     R    p          5  p        r =       5  p        r =                                                                =                                                       @                                                       A                                                      B     
       #         @                                  C                    #PLONS D   #PLATS F   #NPARTS E   #LOOP_START G   #LOOP_END H   #PI I                                                          D                                D                    
 S    p          5  p        r E       5  p        r E                              D                                F                    
 T    p          5  p        r E       5  p        r E                                                                E                                                       G                                                       H                                                      I     
       #         @                                  J                 	   #LEVEL K   #PLEV L   #NZ M   #IX N   #IY O   #IZ P   #RZ Q   #ASCENDING R   #LOGLEV S                                   
                                 K                   
 (             &                   &                   &                                                     
  @                              L     
                
                                  M                     
                                  N                     
                                  O                     D                                 P                      D                                Q     
                 
                                  R                     
                                  S           #         @                                   T                    #ZLAT U   #ZLONG V   #IDUM W   #TIME X   #COZEN Y   #SZA Z                                                            U     
                                                 V     
                 D @                               W                    U    p          p            p                                                                    X     
                 D @                              Y     
                 D                                Z     
              :      fn#fn    Ú   ï      RUNGEK3D    É  ´   a   RUNGEK3D%PLONS    }  ´   a   RUNGEK3D%PLATS    1  ´   a   RUNGEK3D%PLEVS    å  ´   a   RUNGEK3D%CALC      ´   a   RUNGEK3D%TSV    M    a   RUNGEK3D%U0    á    a   RUNGEK3D%V0    u	    a   RUNGEK3D%W0     	    a   RUNGEK3D%LEVEL0 !       a   RUNGEK3D%DLEVDZ0    1    a   RUNGEK3D%UMID    Å    a   RUNGEK3D%VMID    Y    a   RUNGEK3D%WMID "   í    a   RUNGEK3D%LEVELMID #       a   RUNGEK3D%DLEVDZMID        a   RUNGEK3D%UDT    ©    a   RUNGEK3D%VDT    =    a   RUNGEK3D%WDT !   Ñ    a   RUNGEK3D%LEVELDT "   e    a   RUNGEK3D%DLEVDZDT    ù  @   a   RUNGEK3D%NX    9  @   a   RUNGEK3D%NY    y  @   a   RUNGEK3D%NZ !   ¹  6  a   RUNGEK3D%LONGRID !   ï  V  a   RUNGEK3D%LATGRID #   E!  ´   a   RUNGEK3D%LEVELGRID     ù!  @   a   RUNGEK3D%NPARTS $   9"  @   a   RUNGEK3D%LOOP_START "   y"  @   a   RUNGEK3D%LOOP_END    ¹"  @   a   RUNGEK3D%DT    ù"  @   a   RUNGEK3D%PSLAT     9#  @   a   RUNGEK3D%LOGLEV $   y#  @   a   RUNGEK3D%USE_DLEVDZ    ¹#         INTRUV3D    Y%    a   INTRUV3D%U    í&    a   INTRUV3D%V    (    a   INTRUV3D%W     *    a   INTRUV3D%DLEVDZ    ©+    a   INTRUV3D%LEVEL    =-  ´   a   INTRUV3D%PLONS    ñ-  ´   a   INTRUV3D%PLATS    ¥.  ´   a   INTRUV3D%PLEVS    Y/  ´   a   INTRUV3D%PU    0  ´   a   INTRUV3D%PV    Á0  ´   a   INTRUV3D%PW !   u1  ´   a   INTRUV3D%PDLEVDZ    )2  @   a   INTRUV3D%NX    i2  @   a   INTRUV3D%NY    ©2  @   a   INTRUV3D%NZ !   é2  6  a   INTRUV3D%LONGRID !   4  V  a   INTRUV3D%LATGRID #   u5  ´   a   INTRUV3D%LEVELGRID     )6  @   a   INTRUV3D%NPARTS $   i6  @   a   INTRUV3D%LOOP_START "   ©6  @   a   INTRUV3D%LOOP_END     é6  ´   a   INTRUV3D%INDXPS     7  @   a   INTRUV3D%LOGLEV    Ý7  â       PTOLL    ¿8  ´   a   PTOLL%PLONS    s9  ´   a   PTOLL%PLATS    ':  ´   a   PTOLL%INDXPS    Û:  @   a   PTOLL%NPARTS !   ;  @   a   PTOLL%LOOP_START    [;  @   a   PTOLL%LOOP_END    ;  @   a   PTOLL%PI    Û;  ¾       FIXLL    <  ´   a   FIXLL%PLONS    M=  ´   a   FIXLL%PLATS    >  @   a   FIXLL%NPARTS !   A>  @   a   FIXLL%LOOP_START    >  @   a   FIXLL%LOOP_END    Á>  @   a   FIXLL%PI    ?  ¶       GET_LEV_BOX "   ·?  ¼   a   GET_LEV_BOX%LEVEL !   s@  @   a   GET_LEV_BOX%PLEV    ³@  @   a   GET_LEV_BOX%NZ    ó@  @   a   GET_LEV_BOX%IX    3A  @   a   GET_LEV_BOX%IY    sA  @   a   GET_LEV_BOX%IZ    ³A  @   a   GET_LEV_BOX%RZ &   óA  @   a   GET_LEV_BOX%ASCENDING #   3B  @   a   GET_LEV_BOX%LOGLEV    sB         ZEN2    C  @   a   ZEN2%ZLAT    GC  @   a   ZEN2%ZLONG    C     a   ZEN2%IDUM    D  @   a   ZEN2%TIME    [D  @   a   ZEN2%COZEN    D  @   a   ZEN2%SZA 