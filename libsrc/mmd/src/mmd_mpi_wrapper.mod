  8:  ¡   k820309    ?          19.1        ªÍja                                                                                                          
       mmd_mpi_wrapper.f90 MMD_MPI_WRAPPER              gen@MMD_SEND_TO_PARENT gen@MMD_RECV_FROM_PARENT gen@MMD_SEND_TO_CHILD gen@MMD_RECV_FROM_CHILD gen@MMD_BCAST gen@MMD_INTER_BCAST                                                     
                            @                              
       M_TO_PARENT_COMM M_TO_CHILD_COMM M_MODEL_COMM M_MODEL_RANK          @          @                              
       MMD_DP                                                        u #MMD_SEND_TO_PARENT_INTEGER    #MMD_SEND_TO_PARENT_INTEGER_I2 
   #MMD_SEND_TO_PARENT_REAL_R1    #MMD_SEND_TO_PARENT_REAL_R2    #MMD_SEND_TO_PARENT_REAL_R3    #         @     @X                                                 #BUF    #N    #PARENT_RANK    #TAG    #IERR 	             
@ @                                                                &                                                     
@ @                                                    
@ @                                                    
@ @                                                    D @                               	            #         @     @X                             
                    #BUF    #N    #PARENT_RANK    #TAG    #IERR              
@ @                                                                &                   &                                                     
@ @                                                    
@ @                                                    
@ @                                                    D @                                           #         @     @X                                                 #BUF    #N    #PARENT_RANK    #TAG    #IERR              
@ @                                                 
              &                                                     
@ @                                                    
@ @                                                    
@ @                                                    D @                                           #         @     @X                                                 #BUF    #N    #PARENT_RANK    #TAG    #IERR              
@ @                                                 
              &                   &                                                     
@ @                                                    
@ @                                                    
@ @                                                    D @                                           #         @     @X                                                 #BUF    #N    #PARENT_RANK    #TAG     #IERR !             
@ @                                                 
 	             &                   &                   &                                                     
@ @                                                    
@ @                                                    
@ @                                                     D @                               !                                                                   u #MMD_RECV_FROM_PARENT_INTEGER "   #MMD_RECV_FROM_PARENT_INTEGER_I2 (   #MMD_RECV_FROM_PARENT_REAL_R1 .   #MMD_RECV_FROM_PARENT_REAL_R2 4   #MMD_RECV_FROM_PARENT_REAL_R3 :   #         @     @X                             "                    #BUF #   #N $   #PARENT_RANK %   #TAG &   #IERR '             D @                               #                                  &                                                     
@ @                               $                     
@ @                               %                     
@ @                               &                     D @                               '            #         @     @X                             (                    #BUF )   #N *   #PARENT_RANK +   #TAG ,   #IERR -             D @                               )                                  &                   &                                                     
@ @                               *                     
@ @                               +                     
@ @                               ,                     D @                               -            #         @     @X                             .                    #BUF /   #N 0   #PARENT_RANK 1   #TAG 2   #IERR 3             D @                              /                   
               &                                                     
@ @                               0                     
@ @                               1                     
@ @                               2                     D @                               3            #         @     @X                             4                    #BUF 5   #N 6   #PARENT_RANK 7   #TAG 8   #IERR 9             D @                              5                   
               &                   &                                                     
@ @                               6                     
@ @                               7                     
@ @                               8                     D @                               9            #         @     @X                             :                    #BUF ;   #N <   #PARENT_RANK =   #TAG >   #IERR ?             D @                              ;                   
 
              &                   &                   &                                                     
@ @                               <                     
@ @                               =                     
@ @                               >                     D @                               ?                                                                   u #MMD_SEND_TO_CHILD_INTEGER @   #MMD_SEND_TO_CHILD_INTEGER_I2 G   #MMD_SEND_TO_CHILD_REAL_R1 N   #MMD_SEND_TO_CHILD_REAL_R2 U   #MMD_SEND_TO_CHILD_REAL_R3 \   #         @     @X                             @                    #CHILD_ID A   #BUF B   #N C   #CHILD_RANK D   #TAG E   #IERR F             
                                  A                     
@ @                               B                                 &                                                     
@ @                               C                     
@ @                               D                     
@ @                               E                     D @                               F            #         @     @X                             G                    #CHILD_ID H   #BUF I   #N J   #CHILD_RANK K   #TAG L   #IERR M             
                                  H                     
@ @                               I                                 &                   &                                                     
@ @                               J                     
@ @                               K                     
@ @                               L                     D @                               M            #         @     @X                             N                    #CHILD_ID O   #BUF P   #N Q   #CHILD_RANK R   #TAG S   #IERR T             
                                  O                     
@ @                              P                   
              &                                                     
@ @                               Q                     
@ @                               R                     
@ @                               S                     D @                               T            #         @     @X                             U                    #CHILD_ID V   #BUF W   #N X   #CHILD_RANK Y   #TAG Z   #IERR [             
                                  V                     
@ @                              W                   
              &                   &                                                     
@ @                               X                     
@ @                               Y                     
@ @                               Z                     D @                               [            #         @     @X                             \                    #CHILD_ID ]   #BUF ^   #N _   #CHILD_RANK `   #TAG a   #IERR b             
                                  ]                     
@ @                              ^                   
              &                   &                   &                                                     
@ @                               _                     
@ @                               `                     
@ @                               a                     D @                               b                                                                   u #MMD_RECV_FROM_CHILD_INTEGER c   #MMD_RECV_FROM_CHILD_INTEGER_I2 j   #MMD_RECV_FROM_CHILD_REAL_R1 q   #MMD_RECV_FROM_CHILD_REAL_R2 x   #MMD_RECV_FROM_CHILD_REAL_R3    #         @     @X                             c                    #CHILD_ID d   #BUF e   #N f   #CHILD_RANK g   #TAG h   #IERR i             
                                  d                     
D @                               e                                  &                                                     
@ @                               f                     
@ @                               g                     
@ @                               h                     D @                               i            #         @     @X                             j                    #CHILD_ID k   #BUF l   #N m   #CHILD_RANK n   #TAG o   #IERR p             
                                  k                     
D @                               l                                  &                   &                                                     
@ @                               m                     
@ @                               n                     
@ @                               o                     D @                               p            #         @     @X                             q                    #CHILD_ID r   #BUF s   #N t   #CHILD_RANK u   #TAG v   #IERR w             
                                  r                     
D @                              s                   
               &                                                     
@ @                               t                     
@ @                               u                     
@ @                               v                     D @                               w            #         @     @X                             x                    #CHILD_ID y   #BUF z   #N {   #CHILD_RANK |   #TAG }   #IERR ~             
                                  y                     D @                              z                   
               &                   &                                                     
@ @                               {                     
@ @                               |                     
@ @                               }                     D @                               ~            #         @     @X                                                 #CHILD_ID    #BUF    #N    #CHILD_RANK    #TAG    #IERR              
                                                       D @                                                 
               &                   &                   &                                                     
@ @                                                    
@ @                                                    
@ @                                                    D @                                                                                                  u #MMD_BCAST_INTEGER    #MMD_BCAST_CHARACTER    #         @     @X                                                 #BUF    #ROOT_PE    #COMM    #IERR              
D @                                                     
@ @                                                    
 @                                                    F @                                           #         @     @X                                                 #BUF    #ROOT_PE    #COMM    #IERR              
D @                                                  1           
@ @                                                    
 @                                                    F @                                                                                                  u #MMD_INTER_BCAST_INTEGER_1    #MMD_INTER_BCAST_CHARACTER_1    #         @     @X                                                 #BUF    #SENDER    #CHILD_ID    #IERR              
D@                                                                 &                                                     
                                                       
 @                                                    F @                                           #         @     @X                                                 #BUF    #SENDER    #CHILD_ID    #IERR              
D @                                                  1           
                                                       
 @                                                    F @                                                  ,      fn#fn %   Ì      b   uapp(MMD_MPI_WRAPPER    \  @   J  MPI (     {   J  MMD_HANDLE_COMMUNICATOR      G   J  MMD_UTILITIES '   ^  ã       gen@MMD_SEND_TO_PARENT +   A  |      MMD_SEND_TO_PARENT_INTEGER /   ½     a   MMD_SEND_TO_PARENT_INTEGER%BUF -   I  @   a   MMD_SEND_TO_PARENT_INTEGER%N 7     @   a   MMD_SEND_TO_PARENT_INTEGER%PARENT_RANK /   É  @   a   MMD_SEND_TO_PARENT_INTEGER%TAG 0   	  @   a   MMD_SEND_TO_PARENT_INTEGER%IERR .   I  |      MMD_SEND_TO_PARENT_INTEGER_I2 2   Å  ¤   a   MMD_SEND_TO_PARENT_INTEGER_I2%BUF 0   i  @   a   MMD_SEND_TO_PARENT_INTEGER_I2%N :   ©  @   a   MMD_SEND_TO_PARENT_INTEGER_I2%PARENT_RANK 2   é  @   a   MMD_SEND_TO_PARENT_INTEGER_I2%TAG 3   )  @   a   MMD_SEND_TO_PARENT_INTEGER_I2%IERR +   i  |      MMD_SEND_TO_PARENT_REAL_R1 /   å     a   MMD_SEND_TO_PARENT_REAL_R1%BUF -   q  @   a   MMD_SEND_TO_PARENT_REAL_R1%N 7   ±  @   a   MMD_SEND_TO_PARENT_REAL_R1%PARENT_RANK /   ñ  @   a   MMD_SEND_TO_PARENT_REAL_R1%TAG 0   1	  @   a   MMD_SEND_TO_PARENT_REAL_R1%IERR +   q	  |      MMD_SEND_TO_PARENT_REAL_R2 /   í	  ¤   a   MMD_SEND_TO_PARENT_REAL_R2%BUF -   
  @   a   MMD_SEND_TO_PARENT_REAL_R2%N 7   Ñ
  @   a   MMD_SEND_TO_PARENT_REAL_R2%PARENT_RANK /     @   a   MMD_SEND_TO_PARENT_REAL_R2%TAG 0   Q  @   a   MMD_SEND_TO_PARENT_REAL_R2%IERR +     |      MMD_SEND_TO_PARENT_REAL_R3 /     ¼   a   MMD_SEND_TO_PARENT_REAL_R3%BUF -   É  @   a   MMD_SEND_TO_PARENT_REAL_R3%N 7   	  @   a   MMD_SEND_TO_PARENT_REAL_R3%PARENT_RANK /   I  @   a   MMD_SEND_TO_PARENT_REAL_R3%TAG 0     @   a   MMD_SEND_TO_PARENT_REAL_R3%IERR )   É  í       gen@MMD_RECV_FROM_PARENT -   ¶  |      MMD_RECV_FROM_PARENT_INTEGER 1   2     a   MMD_RECV_FROM_PARENT_INTEGER%BUF /   ¾  @   a   MMD_RECV_FROM_PARENT_INTEGER%N 9   þ  @   a   MMD_RECV_FROM_PARENT_INTEGER%PARENT_RANK 1   >  @   a   MMD_RECV_FROM_PARENT_INTEGER%TAG 2   ~  @   a   MMD_RECV_FROM_PARENT_INTEGER%IERR 0   ¾  |      MMD_RECV_FROM_PARENT_INTEGER_I2 4   :  ¤   a   MMD_RECV_FROM_PARENT_INTEGER_I2%BUF 2   Þ  @   a   MMD_RECV_FROM_PARENT_INTEGER_I2%N <     @   a   MMD_RECV_FROM_PARENT_INTEGER_I2%PARENT_RANK 4   ^  @   a   MMD_RECV_FROM_PARENT_INTEGER_I2%TAG 5     @   a   MMD_RECV_FROM_PARENT_INTEGER_I2%IERR -   Þ  |      MMD_RECV_FROM_PARENT_REAL_R1 1   Z     a   MMD_RECV_FROM_PARENT_REAL_R1%BUF /   æ  @   a   MMD_RECV_FROM_PARENT_REAL_R1%N 9   &  @   a   MMD_RECV_FROM_PARENT_REAL_R1%PARENT_RANK 1   f  @   a   MMD_RECV_FROM_PARENT_REAL_R1%TAG 2   ¦  @   a   MMD_RECV_FROM_PARENT_REAL_R1%IERR -   æ  |      MMD_RECV_FROM_PARENT_REAL_R2 1   b  ¤   a   MMD_RECV_FROM_PARENT_REAL_R2%BUF /     @   a   MMD_RECV_FROM_PARENT_REAL_R2%N 9   F  @   a   MMD_RECV_FROM_PARENT_REAL_R2%PARENT_RANK 1     @   a   MMD_RECV_FROM_PARENT_REAL_R2%TAG 2   Æ  @   a   MMD_RECV_FROM_PARENT_REAL_R2%IERR -     |      MMD_RECV_FROM_PARENT_REAL_R3 1     ¼   a   MMD_RECV_FROM_PARENT_REAL_R3%BUF /   >  @   a   MMD_RECV_FROM_PARENT_REAL_R3%N 9   ~  @   a   MMD_RECV_FROM_PARENT_REAL_R3%PARENT_RANK 1   ¾  @   a   MMD_RECV_FROM_PARENT_REAL_R3%TAG 2   þ  @   a   MMD_RECV_FROM_PARENT_REAL_R3%IERR &   >  Þ       gen@MMD_SEND_TO_CHILD *           MMD_SEND_TO_CHILD_INTEGER 3   ¥  @   a   MMD_SEND_TO_CHILD_INTEGER%CHILD_ID .   å     a   MMD_SEND_TO_CHILD_INTEGER%BUF ,   q  @   a   MMD_SEND_TO_CHILD_INTEGER%N 5   ±  @   a   MMD_SEND_TO_CHILD_INTEGER%CHILD_RANK .   ñ  @   a   MMD_SEND_TO_CHILD_INTEGER%TAG /   1  @   a   MMD_SEND_TO_CHILD_INTEGER%IERR -   q        MMD_SEND_TO_CHILD_INTEGER_I2 6   ú  @   a   MMD_SEND_TO_CHILD_INTEGER_I2%CHILD_ID 1   :  ¤   a   MMD_SEND_TO_CHILD_INTEGER_I2%BUF /   Þ  @   a   MMD_SEND_TO_CHILD_INTEGER_I2%N 8     @   a   MMD_SEND_TO_CHILD_INTEGER_I2%CHILD_RANK 1   ^  @   a   MMD_SEND_TO_CHILD_INTEGER_I2%TAG 2     @   a   MMD_SEND_TO_CHILD_INTEGER_I2%IERR *   Þ        MMD_SEND_TO_CHILD_REAL_R1 3   g  @   a   MMD_SEND_TO_CHILD_REAL_R1%CHILD_ID .   §     a   MMD_SEND_TO_CHILD_REAL_R1%BUF ,   3   @   a   MMD_SEND_TO_CHILD_REAL_R1%N 5   s   @   a   MMD_SEND_TO_CHILD_REAL_R1%CHILD_RANK .   ³   @   a   MMD_SEND_TO_CHILD_REAL_R1%TAG /   ó   @   a   MMD_SEND_TO_CHILD_REAL_R1%IERR *   3!        MMD_SEND_TO_CHILD_REAL_R2 3   ¼!  @   a   MMD_SEND_TO_CHILD_REAL_R2%CHILD_ID .   ü!  ¤   a   MMD_SEND_TO_CHILD_REAL_R2%BUF ,    "  @   a   MMD_SEND_TO_CHILD_REAL_R2%N 5   à"  @   a   MMD_SEND_TO_CHILD_REAL_R2%CHILD_RANK .    #  @   a   MMD_SEND_TO_CHILD_REAL_R2%TAG /   `#  @   a   MMD_SEND_TO_CHILD_REAL_R2%IERR *    #        MMD_SEND_TO_CHILD_REAL_R3 3   )$  @   a   MMD_SEND_TO_CHILD_REAL_R3%CHILD_ID .   i$  ¼   a   MMD_SEND_TO_CHILD_REAL_R3%BUF ,   %%  @   a   MMD_SEND_TO_CHILD_REAL_R3%N 5   e%  @   a   MMD_SEND_TO_CHILD_REAL_R3%CHILD_RANK .   ¥%  @   a   MMD_SEND_TO_CHILD_REAL_R3%TAG /   å%  @   a   MMD_SEND_TO_CHILD_REAL_R3%IERR (   %&  è       gen@MMD_RECV_FROM_CHILD ,   '        MMD_RECV_FROM_CHILD_INTEGER 5   '  @   a   MMD_RECV_FROM_CHILD_INTEGER%CHILD_ID 0   Ö'     a   MMD_RECV_FROM_CHILD_INTEGER%BUF .   b(  @   a   MMD_RECV_FROM_CHILD_INTEGER%N 7   ¢(  @   a   MMD_RECV_FROM_CHILD_INTEGER%CHILD_RANK 0   â(  @   a   MMD_RECV_FROM_CHILD_INTEGER%TAG 1   ")  @   a   MMD_RECV_FROM_CHILD_INTEGER%IERR /   b)        MMD_RECV_FROM_CHILD_INTEGER_I2 8   ë)  @   a   MMD_RECV_FROM_CHILD_INTEGER_I2%CHILD_ID 3   +*  ¤   a   MMD_RECV_FROM_CHILD_INTEGER_I2%BUF 1   Ï*  @   a   MMD_RECV_FROM_CHILD_INTEGER_I2%N :   +  @   a   MMD_RECV_FROM_CHILD_INTEGER_I2%CHILD_RANK 3   O+  @   a   MMD_RECV_FROM_CHILD_INTEGER_I2%TAG 4   +  @   a   MMD_RECV_FROM_CHILD_INTEGER_I2%IERR ,   Ï+        MMD_RECV_FROM_CHILD_REAL_R1 5   X,  @   a   MMD_RECV_FROM_CHILD_REAL_R1%CHILD_ID 0   ,     a   MMD_RECV_FROM_CHILD_REAL_R1%BUF .   $-  @   a   MMD_RECV_FROM_CHILD_REAL_R1%N 7   d-  @   a   MMD_RECV_FROM_CHILD_REAL_R1%CHILD_RANK 0   ¤-  @   a   MMD_RECV_FROM_CHILD_REAL_R1%TAG 1   ä-  @   a   MMD_RECV_FROM_CHILD_REAL_R1%IERR ,   $.        MMD_RECV_FROM_CHILD_REAL_R2 5   ­.  @   a   MMD_RECV_FROM_CHILD_REAL_R2%CHILD_ID 0   í.  ¤   a   MMD_RECV_FROM_CHILD_REAL_R2%BUF .   /  @   a   MMD_RECV_FROM_CHILD_REAL_R2%N 7   Ñ/  @   a   MMD_RECV_FROM_CHILD_REAL_R2%CHILD_RANK 0   0  @   a   MMD_RECV_FROM_CHILD_REAL_R2%TAG 1   Q0  @   a   MMD_RECV_FROM_CHILD_REAL_R2%IERR ,   0        MMD_RECV_FROM_CHILD_REAL_R3 5   1  @   a   MMD_RECV_FROM_CHILD_REAL_R3%CHILD_ID 0   Z1  ¼   a   MMD_RECV_FROM_CHILD_REAL_R3%BUF .   2  @   a   MMD_RECV_FROM_CHILD_REAL_R3%N 7   V2  @   a   MMD_RECV_FROM_CHILD_REAL_R3%CHILD_RANK 0   2  @   a   MMD_RECV_FROM_CHILD_REAL_R3%TAG 1   Ö2  @   a   MMD_RECV_FROM_CHILD_REAL_R3%IERR    3  p       gen@MMD_BCAST "   3  r      MMD_BCAST_INTEGER &   ø3  @   a   MMD_BCAST_INTEGER%BUF *   84  @   a   MMD_BCAST_INTEGER%ROOT_PE '   x4  @   a   MMD_BCAST_INTEGER%COMM '   ¸4  @   a   MMD_BCAST_INTEGER%IERR $   ø4  r      MMD_BCAST_CHARACTER (   j5  L   a   MMD_BCAST_CHARACTER%BUF ,   ¶5  @   a   MMD_BCAST_CHARACTER%ROOT_PE )   ö5  @   a   MMD_BCAST_CHARACTER%COMM )   66  @   a   MMD_BCAST_CHARACTER%IERR $   v6         gen@MMD_INTER_BCAST *   ö6  u      MMD_INTER_BCAST_INTEGER_1 .   k7     a   MMD_INTER_BCAST_INTEGER_1%BUF 1   ÷7  @   a   MMD_INTER_BCAST_INTEGER_1%SENDER 3   78  @   a   MMD_INTER_BCAST_INTEGER_1%CHILD_ID /   w8  @   a   MMD_INTER_BCAST_INTEGER_1%IERR ,   ·8  u      MMD_INTER_BCAST_CHARACTER_1 0   ,9  L   a   MMD_INTER_BCAST_CHARACTER_1%BUF 3   x9  @   a   MMD_INTER_BCAST_CHARACTER_1%SENDER 5   ¸9  @   a   MMD_INTER_BCAST_CHARACTER_1%CHILD_ID 1   ø9  @   a   MMD_INTER_BCAST_CHARACTER_1%IERR 