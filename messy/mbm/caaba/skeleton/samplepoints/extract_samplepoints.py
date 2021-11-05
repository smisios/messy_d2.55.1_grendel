#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-07 14:03:33 sander>

# Andrea Pozzer: original code (2017)
# Rolf Sander: some changes and additions (2017-2018)

import numpy as np
from netCDF4 import Dataset

# define final data array
data = []
name = []

# This script must be run on a computer which has access to the global
# model results, e.g., mistral.dkrz.de. The results must be copied back
# into the local directory, e.g.:
#   scp extract_samplepoints.py b302023@mistral.dkrz.de:
#   (execute extract_samplepoints.py on mistral)
#   scp b302023@mistral.dkrz.de:skeleton_samplepoints_30_loworg.log .
#   scp b302023@mistral.dkrz.de:skeleton_samplepoints_30_loworg.nc  .
INPUTNCPATH = '/pf/b/b302052/work/MINOS/MOM/EXP/mom-11_i1m/data/' # on mistral at dkrz
data_tracer = Dataset(INPUTNCPATH+'mom-11_i1m_____20120701_0840_tracer_gp.nc')
data_jval   = Dataset(INPUTNCPATH+'mom-11_i1m_____20120701_0840_jval_gp.nc')
data_physc  = Dataset(INPUTNCPATH+'mom-11_i1m_____20120701_0840_ECHAM5.nc')

# choose a mask:
#MASK = None
MASK = 'loworg'
if (MASK==None):
    OUTPUT = 'skeleton_samplepoints_30'
if (MASK=='loworg'):
    OUTPUT = 'skeleton_samplepoints_30_loworg'

##############################################################################

def boxinfo(i, tracer, level, minmaxave, mylon, mylat, mytime, mytemp, mypress, myval):
    print('Sample point %3d: At level %3d, %s for tracer %s = %g mol/mol' % (i, level+1, minmaxave, tracer, myval))
    print('*** Sample point %3d: At level %3d, %s for tracer %s = %g mol/mol' % (i, level+1, minmaxave, tracer, myval), file=LOGFILE)
    print('  time step %3d = %s (YYYYMMDD)' % (mytime+1, data_physc.variables['YYYYMMDD'][mytime]), file=LOGFILE)
    levelval = data_physc.variables['hyam'][level] + data_physc.variables['hybm'][level] * data_physc.variables['aps'][mytime, mylat, mylon]
    print('  longitude %3d = %s °E' % (mylon+1, data_physc.variables['lon'][mylon]), file=LOGFILE)
    latitude = data_physc.variables['lat'][mylat]
    if (latitude>0):
        print('  latitude  %3d = %s °N' % (mylat+1, latitude), file=LOGFILE)
    else:
        print('  latitude  %3d = %s °S' % (mylat+1, abs(latitude)), file=LOGFILE)
    print('  level     %3d = %g Pa' % (level+1, levelval), file=LOGFILE)
    print('  press         = %g Pa' % (mypress), file=LOGFILE)
    print('  temperature   = %g K' % (mytemp), file=LOGFILE)

##############################################################################

if __name__ == '__main__':

    LOGFILE = open(OUTPUT+'.log','w+', 1) # 1=line-buffered

    #OPEN TARGETS FILE
    target_file = open('targets.txt')
    target_data = target_file.readlines()
    targets=[]

    ntarget=0
    for line_index,line in enumerate(target_data):
        line=line.strip()
        columns = line.split(' ')
        if not line.startswith('#'):
            #print line.rstrip()
            ntarget = ntarget+1
            targets.append(columns[0])

    target_file.close()

    print('%3d targets:' % (len(targets)), end=' ', file=LOGFILE)
    print(targets[:], file=LOGFILE)

    #OPEN LEVEL/PRESSURE FILE
    level_file = open('levels.txt')
    level_data = level_file.readlines()
    levels=[]

    for line_index,line in enumerate(level_data):
        line=line.strip()
        columns = line.split(' ')
        if not line.startswith('#'):
            #print line.rstrip()
            levels.append(int(columns[0]))

    level_file.close()

    print('%3d levels: ' % (len(levels)), end=' ', file=LOGFILE)
    print(levels[:], file=LOGFILE)

    #STORE VARIABLE NAMES
    # adding temperature name in the data
    data.append([]) # value
    name.append('nml_temp')
    # adding pressure name in the data
    data.append([]) # value
    name.append('nml_press')

    dims = data_tracer.dimensions
    for nfield,field in enumerate(data_tracer.variables):
        if field not in dims:
            if field != 'YYYYMMDD' and field !='dt' and field !='nstep' and field!='hyam' and field !='hybm' :
                name.append(field)
                #append elements
                data.append([]) # value
    dims = data_jval.dimensions
    for nfield,field in enumerate(data_jval.variables):
        if field not in dims:
            if field != 'YYYYMMDD' and field !='dt' and field !='nstep' and field!='hyam' and field !='hybm' :
                name.append(field)
                #append elements
                data.append([]) # value

    samponum = 0 # sample point number

    for nlevel,level in enumerate(levels):
        # define a mask of boxes that should be ignored:
        if (MASK=='loworg'):
            mask = np.logical_or.reduce((
                data_tracer.variables['C5H8'][:,level-1,:,:]>1e-10,
                data_tracer.variables['APINENE'][:,level-1,:,:]>1e-10,
                data_tracer.variables['TOLUENE'][:,level-1,:,:]>1e-10))
        for ntrac,target in enumerate(targets):
            index=0
            # define trac either as full array or as masked array:
            if (MASK==None):
                trac = data_tracer.variables[target][:,level-1,:,:] # order: time,lev,lat,lon
            else:
                trac = np.ma.masked_where(mask, data_tracer.variables[target][:,level-1,:,:])
            #LOCATION CALCULATION
            #print(np.unravel_index(np.argmin(trac),trac.shape))
            #print(np.argmin(trac))
            (loc_min_t,loc_min_y,loc_min_x) = np.unravel_index(np.argmin(trac),trac.shape)
            (loc_max_t,loc_max_y,loc_max_x) = np.unravel_index(np.argmax(trac),trac.shape)
            (loc_ave_t,loc_ave_y,loc_ave_x) = np.unravel_index((np.abs(np.mean(trac) - trac)).argmin(),trac.shape)
            trac_min = trac[loc_min_t,loc_min_y,loc_min_x]
            trac_ave = trac[loc_ave_t,loc_ave_y,loc_ave_x]
            trac_max = trac[loc_max_t,loc_max_y,loc_max_x]

            # extract temp and press from data
            temp_min  = data_physc.variables['tm1'][loc_min_t,level-1,loc_min_y,loc_min_x]
            temp_ave  = data_physc.variables['tm1'][loc_ave_t,level-1,loc_ave_y,loc_ave_x]
            temp_max  = data_physc.variables['tm1'][loc_max_t,level-1,loc_max_y,loc_max_x]
            press_min = data_physc.variables['press'][loc_min_t,level-1,loc_min_y,loc_min_x]
            press_ave = data_physc.variables['press'][loc_ave_t,level-1,loc_ave_y,loc_ave_x]
            press_max = data_physc.variables['press'][loc_max_t,level-1,loc_max_y,loc_max_x]
            # store temperature
            data[index].append(temp_min)
            data[index].append(temp_ave)
            data[index].append(temp_max)
            index = index+1
            # store pressure
            data[index].append(press_min)
            data[index].append(press_ave)
            data[index].append(press_max)
            index = index+1
            # store tracer_values
            dims = data_tracer.dimensions
            for nfield,field in enumerate(data_tracer.variables):
                if field not in dims:
                    if field != 'YYYYMMDD' and field !='dt' and field !='nstep' and field!='hyam' and field !='hybm' :
                        val_min=max(0.,data_tracer.variables[field][loc_min_t,level-1,loc_min_y,loc_min_x])
                        val_ave=max(0.,data_tracer.variables[field][loc_ave_t,level-1,loc_ave_y,loc_ave_x])
                        val_max=max(0.,data_tracer.variables[field][loc_max_t,level-1,loc_max_y,loc_max_x])
                        #append elements
                        data[index].append(val_min)
                        data[index].append(val_ave)
                        data[index].append(val_max)
                        index=index+1
            dims = data_jval.dimensions
            for nfield,field in enumerate(data_jval.variables):
                if field not in dims:
                    if field != 'YYYYMMDD' and field !='dt' and field !='nstep' and field!='hyam' and field !='hybm' :
                        #val=data_jval.variables[field][:,level-1,:,:]
                        val_min=data_jval.variables[field][loc_min_t,level-1,loc_min_y,loc_min_x]
                        val_ave=data_jval.variables[field][loc_ave_t,level-1,loc_ave_y,loc_ave_x]
                        val_max=data_jval.variables[field][loc_max_t,level-1,loc_max_y,loc_max_x]
                        #append elements
                        data[index].append(val_min)
                        data[index].append(val_ave)
                        data[index].append(val_max)
                        index=index+1

            samponum += 1
            boxinfo(samponum, target, level-1, 'MIN', loc_min_x, loc_min_y, loc_min_t, temp_min, press_min, trac_min)
            samponum += 1
            boxinfo(samponum, target, level-1, 'AVE', loc_ave_x, loc_ave_y, loc_ave_t, temp_ave, press_ave, trac_ave)
            samponum += 1
            boxinfo(samponum, target, level-1, 'MAX', loc_max_x, loc_max_y, loc_max_t, temp_max, press_max, trac_max)

    #CREATE NETCDF FILES
    samplepoints = Dataset(OUTPUT+'.nc', 'w', format='NETCDF3_CLASSIC')

    #create dimensions
    numtime = len(data[0][:])
    print('*** Total number of samplepoints: %s' % (numtime))
    samplepoints.createDimension('time',numtime)

    for ntrac,tracer in enumerate(name):
        #define variables
        tracer = samplepoints.createVariable(name[ntrac],'d',('time'))
        tracer[:] = data[ntrac][:]

    #close files:
    samplepoints.close()
    LOGFILE.close()

