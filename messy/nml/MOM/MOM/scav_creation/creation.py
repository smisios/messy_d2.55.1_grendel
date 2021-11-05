import numpy as np
#import string as st
#from scipy import stats
#import matplotlib.pyplot as plt
#from netCDF4 import Dataset
import os
import re

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
        if line.startswith('  INTEGER, PARAMETER,PUBLIC :: ind_'):
           trash, sep, info = line.partition('::')
           spc, sep, check = info.partition('=')
           if np.int(check) > 0 :
              out_spc.append(spc[5:])
    file.close() 
    return out_spc

def read_H(fname):
    #OPEN 
    file = open(fname, 'r')
    #output array
    out_spc_H=[]
    out_val_H=[]
    out_MW=[]
    # rewind
    file.seek(0)
    data = file.readlines()
    for line_index,line in enumerate(data):
        line=line
        if line.startswith('|-'):
          continue
        if line.startswith('| <'):
          continue
        if line.startswith('\\'):
          continue
        if line in ['\n', '\r\n']:
          continue
        if line.startswith('#'):
          continue
        val = line.split("|")
        if val[14].strip() == "" :
           continue
        #print(val[1],val[14])
        out_spc_H.append(val[1])
        out_val_H.append(np.float(val[14]))
        out_MW.append(val[3])
    file.close() 
    return out_spc_H, out_val_H, out_MW


def read_proc(fname):
    #OPEN 
    file = open(fname, 'r')
    #output array
    out_spc_scav=[]
    # rewind
    file.seek(0)
    data = file.readlines()
    for line_index,line in enumerate(data):
        line=line
        if line.startswith('|-'):
          continue
        if line.startswith('#'):
          continue
        val = line.split("|")
        #print(val[1],val[13])
        if val[13].strip() == "ON" :
           out_spc_scav.append(val[1])
    file.close() 
    return out_spc_scav

def read_incl(fname):
    #OPEN 
    file = open(fname, 'r')
    #output array
    out_spc_scav=[]
    # rewind
    file.seek(0)
    data = file.readlines()
    for line_index,line in enumerate(data):
        line=line
        if line.startswith('{'):
           continue
        info, sep, trash = line.partition(';')
        spc, sep, spc_comp = info.partition('=')
        spc = spc.strip()[:-3]
        out_spc_scav.append(spc)
    file.close() 
    return out_spc_scav

def read_scav_gas(fname):
    #OPEN 
    file = open(fname, 'r')
    #output array
    out_spc_scav=[]
    # rewind
    file.seek(0)
    data = file.readlines()
    for line_index,line in enumerate(data):
        line=line
        if line.startswith('{'):
           continue
        info, sep, trash = line.partition(';')
        spc, sep, spc_comp = info.partition('=')
        spc = spc.strip()[:]
        out_spc_scav.append(spc)
    file.close() 
    return out_spc_scav



os.system("mv ../../../../mbm/scav/mechanism/aqueous.eqn ../../../../mbm/scav/mechanism/aqueous.eqn.old")
os.system("cp  aqueous_tmpl.eqn ../../../../mbm/scav/mechanism/aqueous.eqn")
os.system("mv ../../../../mbm/scav/mechanism/aqueous.spc ../../../../mbm/scav/mechanism/aqueous.spc.old")
os.system("cp  aqueous_tmpl.spc ../../../../mbm/scav/mechanism/aqueous.spc")
os.system("mv ../../../../mbm/scav/mechanism/gas.spc ../../../../mbm/scav/mechanism/gas.spc.old")
os.system("cp  gas_tmpl.spc ../../../../mbm/scav/mechanism/gas.spc")
os.system("mv ../../../../smcl/messy_cmn_gasaq.f90 ../../../../smcl/messy_cmn_gasaq.f90.old")
os.system("cp  messy_cmn_gasaq_tmpl1.f90 ./messy_cmn_gasaq.1")
os.system("cp  messy_cmn_gasaq_tmpl2.f90 ./messy_cmn_gasaq.2")
os.system("cp  messy_cmn_gasaq_tmpl3.f90 ./messy_cmn_gasaq.3")

spc_eqn = read_spc("../../../../mbm/caaba/messy_mecca_kpp.f90")
#print(spc)
spc_scav = read_proc("../../../../mbm/caaba/mecca/process_gas.tbl")
print(spc_scav)
spc_incl = read_incl("../../../../mbm/scav/mechanism/aqueous.spc")
#print(spc_incl)
spc_H,H_val,spc_MW = read_H("../../../../mbm/tracer/chemprop/messy_main_tracer_chemprop.tbl")
#print(spc_incl,H_val,spc_MW)
spc_gas_scav = read_scav_gas("../../../../mbm/scav/mechanism/gas.spc")
#print(spc_gas_scav)

scav_eqn = open("../../../../mbm/scav/mechanism/aqueous.eqn","a")
scav_spc = open("../../../../mbm/scav/mechanism/aqueous.spc","a")
scav_gspc = open("../../../../mbm/scav/mechanism/gas.spc","a")
gasaq_1 = open("./messy_cmn_gasaq.1","a")
gasaq_2 = open("./messy_cmn_gasaq.2","a")
gasaq_3 = open("./messy_cmn_gasaq.3","a")
scav_eqn.write("// !mz_ap_20190221+: extra by H species (MOM chemistry), based on work of Sergey Gromov }        \n")                                         
scav_spc.write("{ !mz_ap_20190221+: extra by H species (MOM chemistry), based on work of Sergey Gromov }        \n")                                         
scav_gspc.write("{ !mz_ap_20190221+: extra by H species (MOM chemistry), based on work of Sergey Gromov }        \n")                                         
gasaq_1.write("!mz_ap_20190221+: added MOM species!")
gasaq_2.write("!mz_ap_20190221+: added MOM species!")
gasaq_3.write("!mz_ap_20190221+: added MOM species!")

#print("List of species used in the mechanism NOT inserted in SCAV despite l_scav=ON")
index_aqueous=1
for index_eqn,eqn in enumerate(spc_eqn):
    if eqn.strip() =="O3s":
       continue
    for index_scav,scav in enumerate(spc_scav):
        if eqn.strip() == scav.strip(): 
           check = False
           for index_spc,spc in enumerate(spc_incl):
               if eqn.strip() == spc.strip():
                  check = True
           if check==False:
              #print(eqn)
              #find Henry number
              for index_H,H_spc in enumerate(spc_H):
                  if eqn.strip() == H_spc.strip():
                     Hvalue = ""     
                     for power in range(9,0,-1):
                         if H_val[index_H] >=  10**(power) : 
                            Hvalue = Hvalue+"H"+str(power)
                     #if Hvalue=="":
                     #    Hvalue = "H9H8H7H6H54H3H2H1"
                     #print(eqn,H_val[index_H],Hvalue)
                     scav_eqn.write("{#HM"+'%03d'%index_aqueous+"f##} "+eqn.strip()+" = "+eqn.strip()+"_##   : {%TrA##Mom"+Hvalue+"}       k_exf(##,KPP_"+eqn.strip()+"); {&&} \n")
                     scav_eqn.write("{#HM"+'%03d'%index_aqueous+"b##} "+eqn.strip()+"_## = "+eqn.strip()+"   : {%TrA##Mom"+Hvalue+"}       k_exb(##,KPP_"+eqn.strip()+"); {&&} \n")
                     #print("{#HM"+'%03d'%index_aqueous+"f##} "+eqn+" = "+eqn+"_##   : {%TrA##Mom"+Hvalue+"}       k_exf(##,KPP_"+eqn+"); {&&}")
                     index_aqueous=index_aqueous+1
                     scav_spc.write(eqn.strip()+"_##  = IGNORE; {@"+eqn.strip()+"\\aq}    { "+eqn.strip()+" ("+Hvalue+") }\n")
                     exist = False
                     # add to gas.spc if necessary
                     for index_gas_spc,gas_spc in enumerate(spc_gas_scav):
                         if eqn.strip() == gas_spc.strip() :
                            exist = True
                     if exist == False :
                        scav_gspc.write(eqn.strip()+"    = IGNORE; {@"+eqn.strip()+"}        { "+eqn.strip()+" }\n")
                     #molar mass
                     print(eqn,spc_MW[index_H], H_val[index_H])
                     if "." in spc_MW[index_H]:
                        gasaq_1.write("    CALL add_species('"+eqn.strip()+"',   "+spc_MW[index_H].strip()+"_dp  ) \n")
                        gasaq_2.write("    CALL add_henry('"+eqn.strip()+"',   "+str(H_val[index_H])+"_dp, 0._dp   ) \n")
                        gasaq_3.write("    CALL add_alpha('"+eqn.strip()+"',   alpha_T0, alpha_Tdep)   \n")
                     else:
                        MW = re.sub(r'([A-Z])', r'+M\1', spc_MW[index_H].strip())
                        MW = re.sub(r'([1-9][0-9]?)', r'*\1.', MW)
                        MW = re.sub(r'^\+', r' ', MW)
                        gasaq_1.write("    CALL add_species('"+eqn.strip()+"',   "+MW+"  ) \n")
                        gasaq_2.write("    CALL add_henry('"+eqn.strip()+"',   "+str(H_val[index_H])+"_dp, 0._dp   ) \n")
                        gasaq_3.write("    CALL add_alpha('"+eqn.strip()+"',   alpha_T0, alpha_Tdep)   \n")


gasaq_1.close()
gasaq_2.close()
gasaq_3.close()

os.system("cat  messy_cmn_gasaq.1 > ../../../../smcl/messy_cmn_gasaq.f90")
os.system("cat  messy_cmn_gasaq.2 >> ../../../../smcl/messy_cmn_gasaq.f90")
os.system("cat  messy_cmn_gasaq.3 >> ../../../../smcl/messy_cmn_gasaq.f90")
os.system("cat  messy_cmn_gasaq_tmpl4.f90 >> ../../../../smcl/messy_cmn_gasaq.f90")
os.system("rm ./messy_cmn_gasaq.1 ./messy_cmn_gasaq.2 ./messy_cmn_gasaq.3")
