import numpy as np
#import string as st
#from scipy import stats
#import matplotlib.pyplot as plt
#from netCDF4 import Dataset
import re

#bins = [ 1.0E-3, 1.0E-2, 1.0E-1, 1.0E0, 1.0E1, 1.0E2] 
#edges = [ 3.2E-4, 3.2E-3, 3.2E-2, 3.2E-1, 3.2 , 3.2E1, 3.2E2 ] 
# without including part of ORACLE chemistry
#bins = [      1.0E-3,  1.0E-1, 1.0E0, 1.0E1, 1.0E2] 
#edges = [ 3.2E-20, 3.2E-2, 3.2E-1, 3.2 , 3.2E1, 3.2E2 ] 
## including part of ORACLE chemistry
bins = [      1.0E-3,   1.0E-1, 1.0E0, 1.0E1, 1.0E2, 1.0E3] 
edges = [ 3.2E-20, 3.2E-2,  3.2E-1,  3.2 , 3.2E1, 3.2E2,  3.2E3 ] 
nbins = len(bins)
#print(nbins)

def read_spc(fname):
    #OPEN 
    file = open(fname, 'r')
    #output array
    out_spc=[]
    out_spc_comp=[]
    # rewind
    file.seek(0)
    data = file.readlines()
    for line_index,line in enumerate(data):
        line=line
        if line.startswith('{'):
           continue
        info, sep, trash = line.partition(';')
        spc, sep, spc_comp = info.partition('=')
        out_spc.append(spc)
        out_spc_comp.append(spc_comp)
    file.close() 
    return out_spc,out_spc_comp

def read_eqn(fname):
    #OPEN 
    file = open(fname, 'r')
    #output array
    out_spc_eqn=[]
    # rewind
    file.seek(0)
    data = file.readlines()
    for line_index,line in enumerate(data):
        line=line
        val = line.split()
        if int(val[1]) >= 1 or int(val[2]) >= 1 :
           out_spc_eqn.append(val[0])
    file.close() 
    return out_spc_eqn

def comp_num(in_comp, search_comp):
    if search_comp in in_comp: 
       comp_sec = in_comp.split('+')    
       for comp in comp_sec:
           try:
               test = (re.search('(.+?)'+search_comp, comp)).group(1)
           except AttributeError:
               continue
           if test.strip() == '':
              num = float(1.0)
           else:
              num = float(test)
    else:
       num = float(0.0)
    return num

def calc_vol(spc_eqn, spc, spc_comp, vol_min,vol_max):
#loop over all defined species in mecca:
    #output array
    out_trac=[]
    out_mw=[]
    for index_tracer,tracer in enumerate(spc_eqn):
        #print(tracer)
        for index_spc,spc_test in enumerate(spc):
            # check if the define tracer exist in the list
            if tracer.strip() == spc_test.strip():
               # select only organics (i.e. C + H)
               if 'C ' in spc_comp[index_spc] and 'H ' in spc_comp[index_spc]:
                  if 'Br ' not in spc_comp[index_spc] and 'Cl ' not in spc_comp[index_spc] and 'I ' not in spc_comp[index_spc] and 'F ' not in spc_comp[index_spc]:
                          #                  print(spc_comp[index_spc])
                      nC = comp_num(spc_comp[index_spc],'C')
                      nH = comp_num(spc_comp[index_spc],'H')
                      nO = comp_num(spc_comp[index_spc],'O')
                      nN = comp_num(spc_comp[index_spc],'N')
                      nS = comp_num(spc_comp[index_spc],'S')
                      # Apply equation of Li et al. 2016 (Molecular corridors and parameterizations of volatility in the chemical evolution of organic aerosols, ACP, 2016)
                      # log10 C0 - 
                      n0C = 23.80
                      bC  = 0.4861
                      bO  = 0.0
                      bCO = 0.0
                      bN  = 0.0
                      bS  = 0.0
                      if 'O ' in spc_comp[index_spc] :
                         n0C = 22.66
                         bC  = 0.4481
                         bO  = 1.656
                         bCO = -0.7790
                         bN  = 0.0
                         bS  = 0.0
                      if 'N ' in spc_comp[index_spc] :
                         n0C = 24.59
                         bC  = 0.4066
                         bO  = 0.0
                         bCO = 0.0
                         bN  = 0.9619
                         bS  = 0.0
                      if 'O ' in spc_comp[index_spc] and 'N ' in spc_comp[index_spc] :
                         n0C = 24.13
                         bC  = 0.3667
                         bO  = 0.7732
                         bCO = -0.07790
                         bN  = 1.114
                         bS  = 0.0
                      if 'O ' in spc_comp[index_spc] and 'S ' in spc_comp[index_spc]:
                         n0C = 24.06
                         bC  = 0.3637
                         bO  = 1.327
                         bCO = -0.3988
                         bN  = 0.0
                         bS  = 0.7579
                      if 'O ' in spc_comp[index_spc] and 'N ' in spc_comp[index_spc] and 'S ' in spc_comp[index_spc]:
                         n0C = 28.50
                         bC  = 0.3848
                         bO  = 1.011
                         bCO = 0.2921
                         bN  = 1.053
                         bS  = 1.316
                      vol = ((n0C-nC)*bC-nO*bO-2*(nC*nO/(nC+nO))*bCO-nN*bN-nS*bS)
                      mw = 12*nC+16*nO+14*nN+32*nS
                      #print(np.log10(vol_min),vol, np.log10(vol_max))
                      if vol >=  np.log10(vol_min) and vol <= np.log10(vol_max) :
                              out_mw.append(mw)
                              out_trac.append(tracer)
               # SPECIAL CASE TO INCLUDE MISSING VOC/OXIDATION IN ORACLE               
               elif tracer.strip() == 'LaSOGv01' :
                    vol = 1.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  150.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LaSOGv02' :
                    vol = 10.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  150.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LaSOGv03' :
                    vol = 100.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  150.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LaSOGv04' :
                    vol = 1000.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  150.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LbSOGv01' :
                    vol = 1.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  180.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LbSOGv02' :
                    vol = 10.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  180.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LbSOGv03' :
                    vol = 100.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  180.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LbSOGv04' :
                    vol = 1000.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  180.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LaOSOGv01' :
                    vol = 1.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  150.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LaOSOGv02' :
                    vol = 10.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  150.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
               elif tracer.strip() == 'LaOSOGv03' :
                    vol = 100.0
                    if vol >=  (vol_min) and vol <= (vol_max) :
                       mw =  150.0 
                       out_mw.append(mw)
                       out_trac.append(tracer)
#               elif tracer.strip() == 'LbOSOGv01' :
#                    vol = 1.0
#                    if vol >=  (vol_min) and vol <= (vol_max) :
#                       mw =  180.0 
#                       out_mw.append(mw)
#                       out_trac.append(tracer)
#               elif tracer.strip() == 'LbOSOGv02' :
#                    vol = 10.0
#                    if vol >=  (vol_min) and vol <= (vol_max) :
#                       mw =  180.0 
#                       out_mw.append(mw)
#                       out_trac.append(tracer)
#               elif tracer.strip() == 'LbOSOGv03' :
#                    vol = 100.0
#                    if vol >=  (vol_min) and vol <= (vol_max) :
#                       mw =  180.0 
#                       out_mw.append(mw)
#                       out_trac.append(tracer)


    return out_trac,out_mw

 
number_total_soa = 20 + int(nbins)
# write  first part of namelist
oracle = open("oracle.nml","w")
oracle.write("! -*- f90 -*-                                                                                                                         \n")             
oracle.write("!                                                                                                                                     \n")                                         
oracle.write("&CTRL                                                                                                                                 \n")                     
oracle.write("                             ! <--- ORACLE general options                                                                            \n")                
oracle.write("NfPOA    = 5,                !number of Primary OA species from Fossil Fuel emissions                                                 \n")    
oracle.write("NbbPOA   = 4,                !number of Primary OA species from Biomass Burning emissions                                             \n")      
oracle.write("NfSOAsv  = 2,                !number of Secondary OA species from the oxidation of fPOA with csat <  1000 (NFFOA(with csat<1000) -1)  \n")      
oracle.write("NbbSOAsv = 2,                !number of Secondary OA species from the oxidation of bbPOA with csat <  1000 (NBBOA(with csat<1000) -1) \n")         
oracle.write("NfSOAiv  = 4,                !number of Secondary OA species from the oxidation of fPOA with csat >= 1000 (NFFOA-1)                   \n")     
oracle.write("NbbSOAiv = 3,                !number of Secondary OA species from the oxidation of bbPOA with csat >= 1000 (NBBOA-1)                  \n")           
oracle.write("NSOAv    = "+str(nbins)+",                !number of Secondary OA species from the oxidation of traditional biogenic VOCs                          \n")           
oracle.write("NSOAP    = "+str(number_total_soa)+",               !sum of NfPOA,NbbPOA,NfSOAsv,NbbSOAsv,NfSOAiv,NbbSOAiv,NbSOAv,NbOSOAv,NaSOAv,NaOSOAv                     \n")           
oracle.write("aermod   ='gmxe',            !aerosol module for inorganic aerosols : use 'oracle' for its own distribution                           \n")              
oracle.write("!nmode    = 1,                !number of modes used for organic aerosols (1 to 3)                                                     \n")                                                                                  
oracle.write("!tmode(1) = 3                 !type of modes (2=KS, 3=AS, 4=CS)                                                                       \n")                                                         
oracle.write("nmode    = 3,               !number of modes used for organic aerosols (1 to 3)                                                       \n")                                                                                
oracle.write("tmode(1) = 2                !type of modes (2=KS, 3=AS, 4=CS)                                                                         \n")                                                             
oracle.write("tmode(2) = 3                !type of modes (2=KS, 3=AS, 4=CS)                                                                         \n")                                                             
oracle.write("tmode(3) = 4                !type of modes (2=KS, 3=AS, 4=CS)                                                                         \n")                                                             
oracle.write("                                                                                                                                      \n")
oracle.write("mwsoap(1:5)   = 250.d0,250.d0,250.d0,250.d0,250.d0,  !mwsoap(1:NfPOA)   : molecular weights of fPOG/fPOA species (g/mol)             \n")               
oracle.write("csat(1:5)     = 1.d-2, 1.d0,  1.d2,  1.d4,  1.d6,    !csat(1:NfPOA)     : saturation concentrations of fPOG/fPOA species (ug/m3)     \n")                    
oracle.write("cstemp(1:5)   = 300.d0,300.d0,300.d0,300.d0,300.d0,  !cstemp(1:NfPOA)   : temperatures corresponding to csat                         \n")       
oracle.write("deltah(1:5)   = 112.d3,100.d3, 88.d3, 76.d3, 64.d3,  !deltah(1:NfPOA)   : enthalpy of vaporization of fPOG/fPOA species (J/mol)      \n")                       
oracle.write("flagsoap(1:5) = 1,1,1,1,1,                           !flagsoap(1:NfPOA) : set to 1 if fPOG/fPOA species forms solutions; 0 if not    \n")                    
oracle.write("                                                                                                                                     \n")                                                
oracle.write("                                               !n1=NfPOA+1 n2=NbbPOA+NfPOA                                                           \n")                                                                          
oracle.write("mwsoap(6:9)   = 250.d0,250.d0,250.d0,250.d0,  !mwsoap(n1:n2)   : molecular weights of bbPOG/bbPOA species (g/mol)                    \n")                       
oracle.write("csat(6:9)     = 1.d-2, 1.d0,  1.d2,  1.d4,    !csat(n1:n2)     : saturation concentrations of bbPOG/bbPOA species (ug/m3)            \n")                          
oracle.write("cstemp(6:9)   = 300.d0,300.d0,300.d0,300.d0,  !cstemp(n1:n2)   : temperatures corresponding to csat                                  \n")     
oracle.write("deltah(6:9)   = 93.d3, 85.d3, 77.d3, 69.d3,   !deltah(n1:n2)   : enthalpy of vaporization of bbPOG/bbPOA species (J/mol)             \n")               
oracle.write("flagsoap(6:9) = 1,1,1,1,                      !flagsoap(n1:n2) : set to 1 if bbPOG/bbPOA species forms solutions; 0 if not           \n") 
oracle.write("                                                                                                                                      \n") 
oracle.write("                                   !n3=n2+1 n4=n2+NfSOAsv                                                                             \n")
oracle.write("mwsoap(10:11)   = 250.d0,250.d0,  !mwsoap(n3:n4)   : molecular weights of fSOGsv/fSOAsv species (g/mol)                              \n")
oracle.write("csat(10:11)     = 1.d-2, 1.d0,    !csat(n3:n4)     : saturation concentrations of fSOGsv/fSOAsv species (ug/m3)                      \n")
oracle.write("cstemp(10:11)   = 300.d0,300.d0,  !cstemp(n3:n4)   : temperatures corresponding to csat                                              \n")
oracle.write("deltah(10:11)   = 112.d3,100.d3,  !deltah(n3:n4)   : enthalpy of vaporization of fSOGsv/fSOAsv species (J/mol)                       \n")
oracle.write("flagsoap(10:11) = 1,1,            !flagsoap(n3:n4) : set to 1 if fSOGsv/fSOAsv species forms solutions; 0 if not                     \n")
oracle.write("                                                                                                                                      \n")                               
oracle.write("                                   !n5=n4+1 n6=n4+NbbSOAsv                                                                            \n")
oracle.write("mwsoap(12:13)   = 250.d0,250.d0,  !mwsoap(n5:n6)   : molecular weights of bbSOGsv/bbSOAsv species (g/mol)                            \n")
oracle.write("csat(12:13)     = 1.d-2, 1.d0,    !csat(n5:n6)     : saturation concentrations of bbSOGsv/bbSOAsv species (ug/m3)                    \n")
oracle.write("cstemp(12:13)   = 300.d0,300.d0,  !cstemp(n5:n6)   : temperatures corresponding to csat                                              \n")
oracle.write("deltah(12:13)   = 93.d3, 85.d3,   !deltah(n5:n6)   : enthalpy of vaporization of bbSOGsv/bbSOAsv species (J/mol)                     \n")
oracle.write("flagsoap(12:13) = 1,1,            !flagsoap(n5:n6) : set to 1 if bbSOGsv/bbSOAsv species forms solutions; 0 if not                   \n")
oracle.write("                                                                                                                                      \n")
oracle.write("                                                 !n7=n6+1 n8=n6+NfSOAiv                                                               \n")
oracle.write("mwsoap(14:17)   = 250.d0,250.d0,250.d0,250.d0,  !mwsoap(n7:n8)   : molecular weights of fSOGiv/fSOAiv species (g/mol)                \n")
oracle.write("csat(14:17)     = 1.d-2, 1.d0,  1.d2,  1.d4,    !csat(n7:n8)     : saturation concentrations of fSOGiv/fSOAiv species (ug/m3)        \n")
oracle.write("cstemp(14:17)   = 300.d0,300.d0,300.d0,300.d0,  !cstemp(n7:n8)   : temperatures corresponding to csat                                \n")
oracle.write("deltah(14:17)   = 112.d3,100.d3, 88.d3, 76.d3,  !deltah(n7:n8)   : enthalpy of vaporization of fSOGiv/fSOAiv species (J/mol)         \n")
oracle.write("flagsoap(14:17) = 1,1,1,1,                      !flagsoap(n7:n8) : set to 1 if fSOGiv/FFiSSOA species forms solutions; 0 if not      \n")
oracle.write("                                                                                                                                      \n")
oracle.write("                                          !n9=n8+1 n10=n8+NbbSOAiv                                                                    \n")       
oracle.write("mwsoap(18:20)   = 250.d0,250.d0,250.d0,  !mwsoap(n9:n10)   : molecular weights of bbSOGiv/bbSOAiv species (g/mol)                    \n")                                       
oracle.write("csat(18:20)     = 1.d-2, 1.d0,  1.d2,    !csat(n9:n10)     : saturation concentrations of bbSOGiv/bbSOAiv species (ug/m3)            \n")       
oracle.write("cstemp(18:20)   = 300.d0,300.d0,300.d0,  !cstemp(n9:n10)   : temperatures corresponding to csat                                      \n")               
oracle.write("deltah(18:20)   = 93.d3, 85.d3, 77.d3,   !deltah(n9:n10)   : enthalpy of vaporization of bbSOGiv/bbSOAiv species (J/mol)             \n")               
oracle.write("flagsoap(18:20) = 1,1,1,                                                                                                             \n")               
oracle.write("!n11=n10+1 n12=n10+NSOAv                                                                                                              \n")       

oracle.write("mwsoap(21:"+str(number_total_soa)+") = ")
for data in bins:
        oracle.write("180.d0,")
oracle.write("!mwsoap(n11:n12)   : molecular weights of bSOGv/bSOAv species (g/mol)  \n")

oracle.write("csat(21:"+str(number_total_soa)+") = ")  
for data in bins: 
        oracle.write('%6s '%str(data)+","), 
oracle.write(" !csat(n11:n12)     : saturation concentrations of bSOGv/bSOAv species (ug/m3)     \n")                       

oracle.write("cstemp(21:"+str(number_total_soa)+") =" )
for data in bins: 
         oracle.write("300.d0,")
oracle.write(" !cstemp(n11:n12)   : temperatures corresponding to csat                           \n")                       

oracle.write("deltah(21:"+str(number_total_soa)+") =")
# this data is from Epstein et al. A semiempirical COrrelation between Enthalpy of Vaporization
# and saturation concentration for organic aerosol: Delta_H_vap = -11 *log10(C*_300) +129 (eq 13)
# must be converted from kJ to J, i.e. times 1000.0
for data in bins: 
         delta_H = (-11.0*np.log10(data)+129.0)*1000.0
         oracle.write(str(delta_H)+", ")
oracle.write("!deltah(n11:n12)   : enthalpy of vaporization of bSOGv/bSOAv species (J/mol)       \n")               

oracle.write("flagsoap(21:"+str(number_total_soa)+") = ")
for data in bins: 
         oracle.write("1,")
oracle.write("!flagsoap(n11:n12) : set to 1 if bSOGv/bSOAv species forms solutions; 0 if not     \n")       


spc_file, spc_comp_file = read_spc("../../../../mbm/caaba/mecca/mecca.spc")
#print(spc,spc_comp)
spc_eqn_file = read_eqn("../../../../mbm/caaba/mecca/check_eqns_count.log")
#print(spc_eqn_file)

#
########
# VOLATILITY LIST
########
for bins_index,bin in enumerate(bins[0:nbins]):
        print(bins[bins_index],edges[bins_index], edges[bins_index+1])
        nml_value=str(bins_index+1)
        oracle.write("SOGv0"+nml_value+" = '")                                                                     
        bin_vol, bin_mw = calc_vol(spc_eqn_file, spc_file, spc_comp_file, edges[bins_index], edges[bins_index+1])
        for trac in bin_vol: 
            print(trac)
            oracle.write(trac+";")
        oracle.write("'")
        oracle.write("\n")

########
# MOLECULAR WEIGHT LIST
########
oracle.write("! all molecular weights!                                                                                                                                                            \n")
for bins_index,bin in enumerate(bins[0:nbins]):
        print(bins[bins_index],edges[bins_index], edges[bins_index+1])
        nml_value=str(bins_index+1)
        bin_vol, bin_mw = calc_vol(spc_eqn_file, spc_file, spc_comp_file, edges[bins_index], edges[bins_index+1])
        n_mw = str(len(bin_mw))
        oracle.write("SOGv_mw("+nml_value+",1:"+n_mw+") = ")
        for trac in bin_mw: 
             #print(trac)
            oracle.write(str(trac)+",")
        oracle.write("\n")

oracle.write("/                                                                                                                                                                                  \n")  
oracle.write("                                                                                                                                                                                    \n")  

oracle.write("&CPL                                                                                                                                                                               \n")  
oracle.write("!name of emission channel                                                                                                                                                          \n")  
oracle.write("Cemis_channel = 'onemis',               ! onlem                                                                                                                                    \n")  
oracle.write("!names of organic carbon channel elements                                                                                                                                          \n")  
oracle.write("!emis_OC_insol = 'emis_oc'             ! emdep: organic carbon mass insoluble                                                                                                      \n")  
oracle.write("!emis_OC_insol = 'OC_sum_insol',        ! onlem: organic carbon mass insoluble                                                                                                     \n")  
oracle.write("!emis_OC_sol   = 'emis_oc'             ! emdep: organic carbon mass insoluble                                                                                                      \n")  
oracle.write("!emis_OC_sol   = 'OC_sum_sol',          ! onlem: organic carbon mass soluble                                                                                                       \n")  
oracle.write("!                                                                                                                                                                                  \n")  
oracle.write("! emission setup (new way)                                                                                                                                                         \n")  
oracle.write("! for each emission flux one EMIS_CASK should be filled                                                                                                                            \n")  
oracle.write("! all possible fluxes are listed in the messy_ORACLE_e5.f90                                                                                                                        \n")  
oracle.write("! in the subroutine ORACLE_emis_e5_init                                                                                                                                            \n")  
oracle.write("! certain characteristics of the fluxes are defined there as well                                                                                                                  \n")  
oracle.write("! emis_casks - array:                                                                                                                                                              \n")  
oracle.write("!                    1st entry: name (used for identification in the                                                                                                               \n")  
oracle.write("!                               list of the e5 file section emis_init)                                                                                                             \n")  
oracle.write("!                    2nd entry: total scaling factor for incoming flux                                                                                                             \n")  
oracle.write("!                    3rd entry: channel name of the emission flux                                                                                                                  \n")  
oracle.write("!                    4th entry: name of the mass emission flux object                                                                                                              \n")  
oracle.write("!                    5th entry: name of the corresponding number                                                                                                                   \n")  
oracle.write("!                               emission flux object (If it does not exist a                                                                                                       \n")  
oracle.write("!                               number is calculated from the mass)                                                                                                                \n")  
oracle.write("!                    6th entry: \";\" separated list of tracers which should                                                                                                       \n")  
oracle.write("!                               receive emissions from this flux                                                                                                                   \n")  
oracle.write("!                    7th entry: \";\" separated list of fractions of this emission                                                                                                 \n")  
oracle.write("!                               flux for all the tracers defined in entry 5                                                                                                        \n")  
oracle.write("! WARNING: In case the 1st entry from a cask is not matching any of the fluxes in the list                                                                                         \n")  
oracle.write("!          it is ignored                                                                                                                                                           \n")  
oracle.write("!                                                                                                                                                                                  \n")  
oracle.write("!-----------------                                                                                                                                                                 \n")  
oracle.write("! ORGANIC CARBON                                                                                                                                                                   \n")  
oracle.write("!-----------------                                                                                                                                                                 \n")  
oracle.write("!                                                                                                                                                                                  \n")  
oracle.write("! CCMI                                                                                                                                                                             \n")  
oracle.write("! use only one gmxe mode (as)                                                                                                                                                      \n")  
oracle.write("EMIS_CASK(1,:)= \"oc_mass_ff_ks\"  ,\"1.3\",\"offemis\" ,\"EDGAR_ANTH_OC_flux\"   ,\"\"            ,\"fPOA01_ks;fPOA02_ks;LfPOG03;LfPOG04;LfPOG05\",\"0.09;0.23;0.48;0.7;1.0\"     \n")  
oracle.write("! GFEDv3.1                                                                                                                                                                       \n")  
oracle.write("EMIS_CASK(2,:)= \"oc_mass_bb_ks\"  ,\"1.6\",\"bioburn_gp\"    ,\"OC_flux\"    ,\"\"            ,\"bbPOA01_ks;bbPOA02_ks;LbbPOG03;LbbPOG04\",\"0.2;0.2;0.3;0.3\"                    \n")  
oracle.write("!-----------------                                                                                                                                                                 \n")  
oracle.write("l_tendency = T                                                                                                                                                                     \n")  
oracle.write("/                                                                                                                                                                                  \n")  
                                                                                                                                             
