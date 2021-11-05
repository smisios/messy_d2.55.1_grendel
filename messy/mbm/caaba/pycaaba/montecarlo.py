#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-04 14:41:55 sander>

# montecarlo.py = run Monte-Carlo simulations with CAABA
# Rolf Sander, 2018-2019

##############################################################################

import os, sys, shutil
CAABADIR = os.path.realpath(os.path.dirname(__file__)+'/..')
sys.path.append(os.path.realpath(CAABADIR+'/pycaaba'))
# import re # regexp
# import subprocess
import datetime
from glob import glob
from pyteetime import tee # from pycaaba
from rstools import HLINE, grep_i # from pycaaba
from netCDF4 import Dataset, MFDataset

LOGFILENAME = 'montecarlo.log'
DATE_NOW = datetime.datetime.now().strftime('%Y-%m-%d')
TIME_NOW = datetime.datetime.now().strftime('%H:%M:%S')

##############################################################################

if __name__ == '__main__':

    # ensure that output/montecarlo directory exists:
    outputbasedir = CAABADIR + '/output/montecarlo'
    if (not os.path.isdir(outputbasedir)):
        os.mkdir(outputbasedir)
    outputdir = 'output/montecarlo/%s-%s' % (DATE_NOW, TIME_NOW)
    os.mkdir(outputdir)
    os.mkdir(outputdir+'/runs')

    # delete output from old runs:
    list(map(os.remove, glob('caaba_*.nc')))

    # Define 'maxline', the number of Monte-Carlo runs here. Seeds are taken
    # from mcexp_seed.txt which can for example be created with the number
    # generator messy/mbm/rnd/rnd.f90-mecca. Note that maxline must not be
    # larger than the number of lines in mcexp_seed.txt.
    maxline=100
    #maxline=1000
    DifferentMaxTime = False

    LOGFILE = tee.stdout_start(LOGFILENAME, append=False) # stdout
    print('\n%s\nSTART OF MONTE-CARLO MODEL SIMULATIONS\n%s' % (HLINE, HLINE))
    SEED = open('montecarlo/mcexp_seed.txt')
    for line in range(1,maxline+1):
        # read one line from mcexp_seed.txt:
        mcexp_seed = int(SEED.readline())
        # define line number with 4 digits:
        line4 = '%4.4d' % (line)
        os.mkdir(outputdir+'/runs/'+line4)
        # test if mecca.nml is suitable for Monte-Carlo runs:
        test = grep_i('mcexp_seed', 'mecca.nml')
        if (test == ''):
            print('ERROR: mecca.nml does not contain mcexp_seed')
            sys.exit(1)
        os.rename('mecca.nml', 'tmp_mecca.nml') # move to a temporary file
        os.system("sed 's/.*mcexp_seed.*/ mcexp_seed = %d/' tmp_mecca.nml > mecca.nml" % (
            mcexp_seed))
        print('Simulation %s with seed = %10d' % (line4, mcexp_seed), end=' ')
        # run the CAABA/MECCA box model:
        os.system('./caaba.exe > caaba.log')
        # set exitstatus to the value in the file status.log:
        with open('status.log') as f: exitstatus = int(f.read())
        if (exitstatus != 0):
            print("ERROR: exit status from 'caaba.exe' is: %s" % (exitstatus))
            sys.exit(1)
        os.system('cp -p caaba.log caaba.nml mecca.nml *.nc %s/runs/%s' % (outputdir, line4))
        with Dataset('caaba_mecca.nc', 'r') as ncid:
            MaxTime = len(ncid.dimensions['time'])        
        print(' has finished (MaxTime=%d)' % (MaxTime))
        # check if MaxTime is different from previous run:
        if (line > 1):
            if (MaxTime != MaxTime0):
                DifferentMaxTime = True
        MaxTime0 = MaxTime
    SEED.close()
    tee.stdout_stop()

    ##############################################################################

    print('Creating zip file of caaba model code. Please wait...')
    os.system('gmake zip >> %s' % (LOGFILENAME))
    shutil.move('%s.zip' % (os.path.basename(CAABADIR)), outputdir)

    print('\nMerging the netcdf files:')
    olddir = os.getcwd()
    os.chdir(outputdir)
    if (DifferentMaxTime):
        fullfilenames = ['caaba_mecca_c_end.nc', 'caaba_mecca_k_end.nc']
    else:
        fullfilenames = sorted(glob('runs/0001/*.nc'))
    for fullfilename in fullfilenames:
        ncfile = os.path.basename(fullfilename)
        print('Working on %s' % (ncfile))
        print('Monte-Carlo run', end=' ')
        for mcrun in sorted(glob('runs/*')):
            basemcrun = os.path.basename(mcrun)
            print(' %s' % (basemcrun), end=' ')
            # put Monte-Carlo run number into time:
            ncid_out = Dataset(mcrun+'/'+ncfile, 'r+', format='NETCDF3_CLASSIC')
            # r+ = read and modify
            ncid_out.variables['time'][:] = basemcrun
            ncid_out.close()
        print(' done')
        # concatenate files along time:
        ncid_caaba_out = Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')
        ncid_caaba_out.createDimension('time') # unlimited dimension (no size specified)
        ncid_caaba_in = MFDataset('runs/*/'+ncfile)
        for var in ncid_caaba_in.variables:
            if (var=='lon'): continue
            if (var=='lat'): continue
            if (var=='lev'): continue
            # print var, ncid_in.variables[var][:]
            species = ncid_caaba_out.createVariable(var,'d',('time')) # f=float=sp, d=dp
            species[:] = ncid_caaba_in.variables[var][:]
        ncid_caaba_in.close()
        ncid_caaba_out.close()

    # ferret jnl files for scatter plots:
    shutil.move('%s/_scatterplot1.jnl' % (CAABADIR), 'tmp_scatterplot1.jnl')
    shutil.move('%s/_scatterplot2.jnl' % (CAABADIR), '.')
    shutil.move('%s/_histogram_k.jnl'  % (CAABADIR), '.')
    os.system('sort tmp_scatterplot1.jnl > _scatterplot1.jnl')
    os.system('ln -s %s/jnl/tools              .' % (CAABADIR))
    os.system('ln -s %s/jnl/montecarlo.jnl     .' % (CAABADIR))
    os.system('ln -s %s/jnl/scatterplot_mc.jnl .' % (CAABADIR))
    os.remove('tmp_scatterplot1.jnl')
    
    os.chdir(olddir)

    # ##############################################################################

    print('\nThe model output is now in:')
    print(outputdir)
    print('start ferret in that directory and run montecarlo.jnl to plot the results')

    # cleanup:
    list(map(os.remove, glob('tmp_*')))

    print('%s\nEND OF MONTE-CARLO MODEL SIMULATIONS\n%s' % (HLINE, HLINE))
    sys.exit(0)

##############################################################################
