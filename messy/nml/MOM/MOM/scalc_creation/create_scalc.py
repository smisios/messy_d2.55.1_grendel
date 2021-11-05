### Python script for creation of SCALC namelist of SOA species from oracle.nml
### oracle.nml must be created first
### reads oracle.nml and generate scalc channel objects of the form
### SOA_species * SOA_mw, sum, 2
### be aware that primary organic aerosols are hardcoded and not read in via nml
### author: Sebastian Ehrhart, MPIC
### last modified: 18.03.2019

import f90nml # not in anconda, use pip install --user f90nml to install
import json


#############################
### functions
##############################
def scalcstring_v1(OA, chanam, mw_SOA, iscal=4):
  y = "CALC("+str(iscal)+")="+chanam+",'tracer_gp:"
  for key in OA:
    y=y+key+"%"+"{:.2f}".format(mw_SOA[key])+","
  # remove last ,
  y = y[:-1] 
  y=y+"' , 'SUM','messy_global_end',"
  return y


# def store_json(x,fjson):
#   import json
#   with open(fjson, 'w') as f:
#     json.dump(x, f, indent=0, sort_keys=True)


def primSOA(istart):
  i = istart
  def inc():
    nonlocal i
    i =  i + 1
    return i
  y = []
  y.append("CALC("+str(i)+")='mass_fSOAsv_ks','tracer_gp:fSOAsv01_ks%250.00,fSOAsv02_ks%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_fSOAsv_as','tracer_gp:fSOAsv01_as%250.00,fSOAsv02_as%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_fSOAsv_cs','tracer_gp:fSOAsv01_cs%250.00,fSOAsv02_cs%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbSOAsv_ks','tracer_gp:bbSOAsv01_ks%250.00,bbSOAsv02_ks%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbSOAsv_as','tracer_gp:bbSOAsv01_as%250.00,bbSOAsv02_as%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbSOAsv_cs','tracer_gp:bbSOAsv01_cs%250.00,bbSOAsv02_cs%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_fSOAiv_ks','tracer_gp:fSOAiv01_ks%250.00,fSOAiv02_ks%250.00,fSOAiv03_ks%250.00,fSOAiv04_ks%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_fSOAiv_as','tracer_gp:fSOAiv01_as%250.00,fSOAiv02_as%250.00,fSOAiv03_as%250.00,fSOAiv04_as%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_fSOAiv_cs','tracer_gp:fSOAiv01_cs%250.00,fSOAiv02_cs%250.00,fSOAiv03_cs%250.00,fSOAiv04_cs%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbSOAiv_ks','tracer_gp:bbSOAiv01_ks%250.00,bbSOAiv02_ks%250.00,bbSOAiv03_ks%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbSOAiv_as','tracer_gp:bbSOAiv01_as%250.00,bbSOAiv02_as%250.00,bbSOAiv03_as%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbSOAiv_cs','tracer_gp:bbSOAiv01_cs%250.00,bbSOAiv02_cs%250.00,bbSOAiv03_cs%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_fPOA_ks','tracer_gp:fPOA01_ks%250.00,fPOA02_ks%250.00,fPOA03_ks%250.00,fPOA04_ks%250.00,fPOA05_ks%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_fPOA_as','tracer_gp:fPOA01_as%250.00,fPOA02_as%250.00,fPOA03_as%250.00,fPOA04_as%250.00,fPOA05_as%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_fPOA_cs','tracer_gp:fPOA01_cs%250.00,fPOA02_cs%250.00,fPOA03_cs%250.00,fPOA04_cs%250.00,fPOA05_cs%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbPOA_ks','tracer_gp:bbPOA01_ks%250.00,bbPOA02_ks%250.00,bbPOA03_ks%250.00,bbPOA04_ks%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbPOA_as','tracer_gp:bbPOA01_as%250.00,bbPOA02_as%250.00,bbPOA03_as%250.00,bbPOA04_as%250.00' , 'SUM','messy_global_end',")
  y.append("CALC("+str(inc())+")='mass_bbPOA_cs','tracer_gp:bbPOA01_cs%250.00,bbPOA02_cs%250.00,bbPOA03_cs%250.00,bbPOA04_cs%250.00' , 'SUM','messy_global_end',")
  return y


### MOM_NAN
fnscalc = "./scalc_template.nml"

meccafn = "/mnt/lustre01/pf/b/b302070/src/messy_2.53.0_MOM_NAN/messy/mbm/caaba/mecca/messy_mecca_trac_si.inc"

scalcnml = f90nml.read(fnscalc)
iscal = len(scalcnml["CPL"]["CALC"]) + 1

OAmodes = ["ks","as","cs"] 

nml = f90nml.read("../oracle.nml")
# store_json(nml,"oracle_nml.json")

## find SOA
SOA = {}
SOGv_mw = {}
SOGv = [ nml["CTRL"]["SOGv"+str(i).zfill(2)].split(";") for  i in range(1,nml["CTRL"]["NSOAv"]+1) ]
nSOG = [ sum( len(i) > 0 for i in k ) for k in SOGv ]


for key in OAmodes:
  for iSOA in range(1,nml["CTRL"]["NSOAv"]+1): 
    SOA[str(iSOA).zfill(2)+"_"+key] = [ "SOAv"+str(iSOA).zfill(2)+str(iSOG).zfill(2)+"_"+key for iSOG in range(1,nSOG[iSOA-1]+1)  ]


for key in OAmodes:
  for iSOA in range(1,nml["CTRL"]["NSOAv"]+1):
    for iSOG in range(1,nSOG[iSOA-1]+1):
      nameSOA = "SOAv"+str(iSOA).zfill(2)+str(iSOG).zfill(2)+"_"+key
      SOGv_mw[nameSOA] = nml["CTRL"]["SOGv_mw"][iSOG-1][iSOA-1]


bla = []
for key in SOA:
  bla.append(scalcstring_v1(SOA[key], "'mass_SOA_"+key+"'", SOGv_mw, iscal))
  iscal = iscal + 1



oOAnames = ["fSOAsv","bbSOAsv","fSOAiv", "bbSOAiv", "fPOA","bbPOA"]

bla.extend(primSOA(iscal))


### read scalc.nml
with open(fnscalc, "r") as f:
  orisc = f.readlines()


newsc = []
for line in orisc:
  if line == "/\n":
    for item in bla:
      newsc.append(item+"\n")
      print(len(item))
  else:
    newsc.append(line) 
newsc.append("/\n")

with open("./scalc.nml","w") as f:
  for line in newsc:
    f.write(line)

# store_json(bla,"scalc.json")


