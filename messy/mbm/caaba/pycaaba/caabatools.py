#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-03 17:21:14 sander>

# caabatools
# Rolf Sander, 2018

##############################################################################

import os
# import sys, shutil
import re # regexp
from netCDF4 import Dataset

##############################################################################

# list all KPP species in spcfile
def list_spc(spcfile):
    spclist = []
    regexp = re.compile('^[ \t]*([A-Za-z][A-Za-z0-9_#]*)[ \t]*=.*;')
    f = open(spcfile,'r')
    for line in f:
        result = regexp.search(line)
        if result:
            spclist.append(result.group(1))
    f.close()
    # for spc in spclist:
    #     print spc
    return spclist
 
##############################################################################

def split_caaba_mecca_nc(verbosity=1):

    with Dataset('caaba_mecca.nc', 'r') as ncid:
        if any('RR' in s[0:2] for s in ncid.variables):
            if (verbosity>0):
                print('\nSplitting caaba_mecca.nc into mixing ratios and reaction rates...')
            os.rename('caaba_mecca.nc', 'caaba_mecca-ori.nc')
            with Dataset('caaba_mecca-ori.nc', 'r')                          as ncid_in,   \
                 Dataset('caaba_mecca.nc',    'w', format='NETCDF3_CLASSIC') as ncid_out1, \
                 Dataset('caaba_mecca_rr.nc', 'w', format='NETCDF3_CLASSIC') as ncid_out2:
                # copy global attributes all at once via dictionary:
                ncid_out1.setncatts(ncid_in.__dict__)
                ncid_out2.setncatts(ncid_in.__dict__)
                # copy dimensions:
                for name, dimension in list(ncid_in.dimensions.items()):
                    ncid_out1.createDimension(
                        name, (len(dimension) if not dimension.isunlimited() else None))
                    ncid_out2.createDimension(
                        name, (len(dimension) if not dimension.isunlimited() else None))
                for name, variable in list(ncid_in.variables.items()):
                    # put mixing ratios and dimension variables into caaba_mecca.nc:
                    if (name[0:2]!='RR' or name in ncid_in.dimensions):
                        ncid_out1.createVariable(name, variable.datatype, variable.dimensions)
                        ncid_out1[name][:] = ncid_in[name][:]             # copy data
                        ncid_out1[name].setncatts(ncid_in[name].__dict__) # copy attributes
                    # put reaction rates and dimension variables into caaba_mecca_rr.nc:
                    if (name[0:2]=='RR' or name in ncid_in.dimensions):
                        ncid_out2.createVariable(name, variable.datatype, variable.dimensions)
                        ncid_out2[name][:] = ncid_in[name][:]             # copy data
                        ncid_out2[name].setncatts(ncid_in[name].__dict__) # copy attributes
            os.remove('caaba_mecca-ori.nc')

##############################################################################

if __name__ == '__main__':
    # only for testing
    #split_caaba_mecca_nc()
    #spclist = list_spc('mecca/gas.spc')
    pass
