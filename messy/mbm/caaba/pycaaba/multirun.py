#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-12-11 17:10:35 sander>

# multirun.py calls CAABA/MECCA for each entry in input file
# Rolf Sander, 2018

##############################################################################

import os, sys, shutil
caabadir = os.path.realpath(os.path.dirname(__file__)+'/..')
sys.path.append(os.path.abspath(caabadir+'/pycaaba'))
from netCDF4 import Dataset, MFDataset
import subprocess
import time
from glob import glob
from pyteetime import tee # from pycaaba
from rstools import HLINE, runcmd, isnumber, corename # from pycaaba/
import f90nml # pip3 install --user wheel f90nml

DONTEDIT = 'created automatically by ' + os.path.basename(__file__) + ', DO NOT EDIT!'

##############################################################################

class multirun(object):
    """ perform multple CAABA runs
    """

    @classmethod
    def get_MaxTime(cls, ncfile):
        ncid = Dataset(ncfile, 'r') # r = read only
        MaxTime = len(ncid.dimensions['time'])
        ncid.close()
        return MaxTime

    ##########################################################################

    @classmethod
    def makeruns(cls, fullncfile, outputdir, write_status=False):
        from caabatools import split_caaba_mecca_nc
        TMP_NML = 'tmp_caaba.nml'
        corencfile = corename(fullncfile)
        if (not os.path.isfile(fullncfile)):
            sys.exit('ERROR: %s does not exist' % (fullncfile))
        os.makedirs(outputdir+'/runs')
        ncid_in = Dataset(fullncfile, 'r')

        nmldict = f90nml.read('caaba.nml') # read *.nml file
        nmlfilename = os.readlink('caaba.nml') # save name of original link
        os.system('ln -fs ' + TMP_NML + ' caaba.nml')

        timedim = None
        unlimited = ''
        if ('time' in ncid_in.dimensions):
            # normally, data are listed along the 'time' axis:
            timedim = 'time'
        else:
            # otherwise, obtain name and length of the unlimited dimension:
            for dim in ncid_in.dimensions:
                if (ncid_in.dimensions[dim].isunlimited()):
                    timedim = dim
                    unlimited = 'unlimited '
        if (timedim):
            lines = len(ncid_in.dimensions[timedim])
            print('%s\nSTART OF MULTIRUN (%d runs from %s"%s" axis of %s.nc)\n%s\n' % (
                HLINE, lines, unlimited, timedim, corencfile, HLINE))
        else:
            sys.exit('ERROR: Input file %s has neither "time" axis nor an unlimited dimension'
                     % (corencfile))

        # check if input file overwrites any namelist values:
        for var in ncid_in.variables:
            if (var[0:4]=='nml_') and (var[4:] in nmldict['CAABA']):
                print('WARNING: The value for ' + var[4:] + ' in ' + nmlfilename + ' will be')
                print('         overwritten by values from ' + fullncfile)
                    
        #---------------------------------------------------------------------

        for line in range(lines):
            # define line number with 4 digits, starting with 0001:
            line4 = '%4.4d' % (line+1)
            input_ncfile = 'input_'+corencfile+'_'+line4+'.nc'
            # extract one line from fullncfile to input_ncfile:
            ncid_out = Dataset(input_ncfile, 'w', format='NETCDF3_CLASSIC')
            ncid_out.createDimension('time') # no size -> unlimited dimension
            for var in ncid_in.variables:
                #print var, ncid_in.variables[var][line]
                if (var[0:4]=='nml_'):
                    nmldict['CAABA'][var[4:]] = ncid_in.variables[var][line]
                else:
                    species = ncid_out.createVariable(var,'d',('time')) # f=float=sp, d=dp
                    species[:] = ncid_in.variables[var][line]
            ncid_out.close()

            # create nml:
            nmldict['CAABA']['init_spec']   = input_ncfile # initial conc
            nmldict['CAABA']['input_readj'] = input_ncfile # j-values
            nmldict.write(TMP_NML, force=True)
            
            print('Running simulation %s' % (line4), end=' ')
            
            # remove old *.nc files:
            list(map(os.remove, glob('caaba_*.nc')))
            # run the CAABA/MECCA box model:
            runcmd('./caaba.exe', 'caaba.log', verbosity=0)
            split_caaba_mecca_nc(verbosity=0)
            # move/copy model output to another directory:
            destdir = outputdir+'/runs/'+line4 # destination directory
            os.mkdir(destdir)
            shutil.move(input_ncfile, destdir+'/input.nc')
            for wildcard in (['caaba.log', 'caaba.nml', 'mecca.nml', '*.nc', '*.dat']):
                for filename in glob(wildcard):
                    #print 'copying %s to %s' % (filename, destdir)
                    shutil.copy2(filename, outputdir+'/runs/'+line4)
            print(': finished (MaxTime=%d)' % (cls.get_MaxTime('caaba_mecca.nc')))

        #---------------------------------------------------------------------

        if write_status:
            with open(outputdir + '/status.log', 'w+') as f:
                print('0', file=f)

        # cleanup:
        list(map(os.remove, glob('tmp_*')))
        os.system('ln -fs ' + nmlfilename + ' caaba.nml')

    ##########################################################################

    @classmethod
    def postprocessing(cls, outputdir):
        if (os.path.isfile('Makefile')):
            zipfilename = os.path.basename(os.getcwd())+'.zip'
            runcmd('gmake zip', 'gmake_zip.log')
            shutil.move(zipfilename, outputdir)

        print('\nMerging the netcdf files:')
        olddir = os.getcwd()
        os.chdir(outputdir)

        # for all runs, copy MaxTime info from caaba_mecca.nc to caaba_mecca_c_end.nc:
        for murun in sorted(glob('runs/*')):
            ncid_caaba_out = Dataset( # r+ = read and modify:
                murun+'/caaba_mecca_c_end.nc', 'r+', format='NETCDF3_CLASSIC')
            MaxTimevar = ncid_caaba_out.createVariable(
                'MaxTime','d',('time','lev','lat','lon')) # f=float=sp, d=dp
            MaxTimevar[:] = cls.get_MaxTime(murun+'/caaba_mecca.nc')
            ncid_caaba_out.close()

        fullfilenames = ['caaba_mecca_c_end.nc', 'caaba_mecca_k_end.nc', 'input.nc']
        for fullfilename in fullfilenames:
            ncfile = os.path.basename(fullfilename)
            print('Working on %s' % (ncfile))
            print('Multirun', end=' ')
            for murun in sorted(glob('runs/*')):
                basemurun = os.path.basename(murun)
                print(' %s' % (basemurun), end=' ')
                if (not fullfilename == 'input.nc'):
                    # put Multirun number into time:
                    ncid_out = Dataset(murun+'/'+ncfile, 'r+', format='NETCDF3_CLASSIC')
                    # r+ = read and modify
                    ncid_out.variables['time'][:] = basemurun
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
        
        #---------------------------------------------------------------------

        keywords = ['real', 'user', 'sys']
        timing_results = { 'real': list(), 'user': list(), 'sys': list() }
        for murun in sorted(glob('runs/*')):
            with open(murun + '/caaba.log', 'r') as logfile:
                lines = logfile.readlines()
                # scan the log file for the timings
                for keyword in keywords:
                    sindex = -1
                    for i in reversed(range(len(lines))):
                        sindex = lines[i].find(keyword)
                        if sindex >= 0:
                            line = lines[i]
                            while line[sindex].isalpha() or line[sindex].isspace():
                                sindex += 1
                            start_index = sindex
                            while isnumber(line[sindex]) or line[sindex] == '.':
                                sindex += 1
                            result_str = line[start_index : sindex]
                            timing_results[keyword].append(float(result_str))
                            break
                    assert sindex >= 0, keyword + ' time not found in log file'
        ncid_time = Dataset('timings.nc', 'w', format='NETCDF3_CLASSIC')
        ncid_time.createDimension('time')
        real_var = ncid_time.createVariable('real_time', 'f', ('time'))
        user_var = ncid_time.createVariable('user_time', 'f', ('time'))
        sys_var  = ncid_time.createVariable('sys_time',  'f', ('time'))
        real_var[:] = timing_results['real']
        user_var[:] = timing_results['user']
        sys_var[:]  = timing_results['sys']
        ncid_time.close()
        
        #---------------------------------------------------------------------

        os.chdir(olddir)

    ##########################################################################

    @classmethod
    def complete(cls, fullncfile, outputdir=False, write_status=False):
        if (not outputdir):
            outputdir = caabadir + '/output/multirun/' + corename(fullncfile)
            # if outputdir exists already, rename it, including its modification time:
            if (os.path.isdir(outputdir)):
                os.rename(outputdir, outputdir + '-' + time.strftime(
                    '%Y-%m-%d-%H:%M:%S', time.localtime(os.path.getmtime(outputdir))))
        cls.makeruns(fullncfile, outputdir, write_status)
        cls.postprocessing(outputdir)
        print('\n%s\nEND OF MULTIRUN\n%s' % (HLINE, HLINE))
        return outputdir

##############################################################################

if __name__ == '__main__':

    # ensure that multirun.py is started from main caaba/ directory:
    if (os.getcwd() != caabadir):
        sys.exit('ERROR: multirun.py must be started from caaba ' +
                 'directory, preferably via xcaaba.py')

    logfilename = 'multirun.log'
    tee.stdout_start(logfilename, append=False) # stdout

    # first command line parameter must be a netcdf input file:
    if len(sys.argv) > 1:
        fullncfile = os.path.abspath(sys.argv[1]) # convert to fullfilename
    else:
        sys.exit('ERROR: start multirun.py via xcaaba.py ' +
                 '(or provide a netcdf input file as a parameter)')

    if len(sys.argv) > 2:
        outputdir = multirun.complete(fullncfile, sys.argv[2])
    else:
        outputdir = multirun.complete(fullncfile)
    
    tee.stdout_stop()
    shutil.move(logfilename, outputdir)
