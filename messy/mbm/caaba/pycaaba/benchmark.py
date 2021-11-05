#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-27 12:52:52 sander>

# benchmark.py does multiple runs for one integrator with
# user defined parameters automatically.

# Authors:
#   Florian Strunk, 2018-2019
#   Rolf Sander,    2018-2019

# (small) manual:
"""
call ./pycaaba/benchmark.py from the main directory
use the -h command line argument for help.

Examples:
./pycaaba/benchmark.py -integ rosenbrock_posdef_h211b_qssa,_,rosenbrock_posdef,_ -icntrl 3 1,2,1,2
./pycaaba/benchmark.py -q -integ rosenbrock_posdef_h211b_qssa,_,_,rosenbrock_posdef,_,_ -icntrl 3 1,2,3,1,2,3
./pycaaba/benchmark.py -integ rosenbrock_posdef,_,kpp_radau5,_ -icntrl 3 1,2,1,2
./pycaaba/benchmark.py -integ rosenbrock_posdef_h211b_qssa,_,kpp_radau5,_,rosenbrock_posdef,_ -icntrl 3 1,2,1,2,1,2
./pycaaba/benchmark.py -integ rosenbrock_posdef,_,kpp_radau5,_ -icntrl 3 1,2,1,2 -i skeleton/samplepoints/skeleton_samplepoints_30.nc -b mom
./pycaaba/benchmark.py -icntrl 3 1,2,3
./pycaaba/benchmark.py -integ 0 -icntrl 3 1,2,3
./pycaaba/benchmark.py -p
./pycaaba/benchmark.py -p -d benchmark_renamed_directory
./pycaaba/benchmark.py -icntrl 3 range 0 5
./pycaaba/benchmark.py -rcntrl 8 range 1.0 7.0 2.0
./pycaaba/benchmark.py -icntrl 3 1,2,3 -ref
./pycaaba/benchmark.py -combine -icntrl 3 1,2,3 -rcntrl 8 1.0,2.0,3.0

--------------------------
Usage notes:
- If you want to automatically insert a reference run with default values for all settings
  using the posdef integrator:
  use -ref / --reference
- If you want a specific settings/integrator combination to be the reference, remember to put that combination
  as the first setting. This is important when plotting the runs as the timing results are plotted relative to the
  first setting.
- Use "-integ 0" argument if you don't wish to recompile the default integrator
- Ranges of setting values work for integer and floating point numbers (floating point might be less safe to use)
  syntax is '-'('icntrl' | 'rcntrl') <SETTING_INDEX> 'range' <START> <STOP> [STEP]
- When -combine is set, all possible combinations of the specificed settings are used for the runs,
  this excludes integrators
- Use the -p argument to plot the results in output/benchmark.
  You can specify the output directory to plot with -d <DIRECTORY NAME>
  (<DIRECTORY NAME> is relative to output/)

"""

# info about rosenbrock_posdef_h211b_qssa:
# ICNTRL(5)  -> 0, 1 oder 2 (posdef: never, only at end, at every substep)
# ICNTRL(6)  -> 0 oder 1...5 (QSSA 0 = off or 1...5 = on)
# RCNTRL(8)  -> 1...8 see Eqn (31) in Söderlind (2003) (4.5 is best)
# RCNTRL(9)  -> order from ICNTRL(3), e.g. 4 for Ros4
# RCNTRL(10) -> thres_tau threshold related to qssa (default 1E-3)
# RCNTRL(11) -> Iqssamax tolerance related to qssa (default 1E-2)
# RCNTRL(12) -> rel tolerance (default 1E-2 is in messy_mecca_kpp.kpp)
# RCNTRL(13) -> abs tolerance (default 10 is in messy_mecca_kpp.kpp)
# NOTE: setting RCNTRL(12) and RCNTRL(13) overwrites reltol and
#       abstol from messy_mecca_kpp.kpp

##############################################################################

# structure of output/ directory:
# benchmark[-%Y-%m-%d-%H:%M:%S]/   -> results from previous runs of benchmark.py
# benchmark/                       -> the latest results
#    benchmark.log                 -> qexec.##.[out|err]
#    integr_settings_####/         -> one directory for each integrator settings
#       benchmark.bat              -> only when integrator changes
#       xmecca.log                 -> only when integrator changes
#       gmake.log                  -> only when integrator changes
#       integr_settings.pkl        -> integrator, icntrl, rcntrl
#       caaba_mecca_c_end.nc       -> final concentrations for all runs
#       timings.nc                 -> real, user, sys for all runs
#       runs/                      -> directory for the runs
#          ####/                   -> one subdirectory for each run
#             caaba_mecca_c_end.nc -> one for each run

##############################################################################

import os, sys, shutil
caabadir = os.path.realpath(os.path.dirname(__file__)+'/..')
sys.path.append(os.path.abspath(caabadir+'/pycaaba'))
from netCDF4 import Dataset, MFDataset
import subprocess
import time
from glob import glob
from pyteetime import tee # from pycaaba
from rstools import HLINE, runcmd, isfloat, isint, drange, corename # from pycaaba
from datetime import datetime
import shutil
from ast import literal_eval
import decimal
import copy
import numpy as np
import re
import matplotlib.pyplot as plt
from viewport import viewport
import time

DONTEDIT = 'created automatically by ' + os.path.basename(__file__) + ', DO NOT EDIT!'

#@TODO define '/integr_settings' as variable

#@TODO: replace raw asserts with this. Use global debug setting
def benchmark_assert(condition, message, debug = True):
    if debug:
        assert condition, message

def run_xmecca(batchfile, integrator_name):
    # prepare_batch_file:
    print('\nmodifying batch file %s.bat' % (batchfile))
    batchdir = caabadir + '/mecca/batch/'
    with open(batchdir + batchfile + '.bat', 'r') as inputfile, \
         open(batchdir + 'benchmark.bat',   'w') as outputfile:
        outputfile.write('# '+DONTEDIT+'\n')
        for line in inputfile:
            line0 = line # backup of original line
            line = re.sub('(^\s*set\s+integr\s+=\s+).*', r'\1'+integrator_name, line)
            line = re.sub('(^\s*set\s+latex\s+=\s+).*', r'\1n', line)    # latex = n
            line = re.sub('(^\s*set\s+graphviz\s+=\s+).*', r'\1n', line) # graphviz = n
            outputfile.write(line)
            if line!=line0: 
                print('<'+line0, '>'+line) # print what was changed
    print('benchmark.bat has been created\n')
    # compile_integrator:
    olddir = os.getcwd()
    os.chdir(caabadir + '/mecca')
    runcmd('xmecca benchmark', 'xmecca.log')
    os.chdir(caabadir)
    runcmd('gmake', 'gmake.log')
    os.chdir(olddir)

def prepare_mecca_nml_file(icntrl_dict = dict(), rcntrl_dict = dict()):

    nml_default   = 'nml/mecca_default.nml'
    nml_benchmark = 'nml/mecca_benchmark.nml'

    print('\nmodifying mecca_default.nml')
    with open(nml_default,   'r') as inputfile, \
         open(nml_benchmark, 'w') as outputfile:
        outputfile.write('! '+DONTEDIT+'\n')
        for line in inputfile:
            line0 = line # backup of original line
            for key, value in icntrl_dict.items():
                line = re.sub('(^\s*&CTRL_KPP.*)', r'\1\nicntrl(%s) = %s ! ADDED'
                              % (key, value), line)
                line = re.sub('(^\s*icntrl\s*\(\s*%s\s*\)\s*=.*)' % (key),
                              r'! DEACTIVATED \1', line)
            for key, value in rcntrl_dict.items():
                line = re.sub('(^\s*&CTRL_KPP.*)', r'\1\nrcntrl(%s) = %s ! ADDED'
                              % (key, value), line)
                line = re.sub('(^\s*rcntrl\s*\(\s*%s\s*\)\s*=.*)' % (key),
                              r'! DEACTIVATED \1', line)
            outputfile.write(line)
            if line!=line0: 
                print('<'+line0, '>'+line) # print what was changed
    print('mecca_benchmark.nml has been created\n')
    if (os.path.isfile('mecca.nml')): os.remove('mecca.nml')
    os.symlink(nml_benchmark, 'mecca.nml')

def qexec_multirun(full_nc_file, multirun_outputdir, qexec=False):
    """ runs multirun synchronously or asynchronously inside an isolated copy of caabadir """

    os.mkdir(multirun_outputdir)
    # copy info about mechanism to multirun outputdir:
    shutil.copy(caabadir + '/mecca/batch/benchmark.bat', multirun_outputdir)
    shutil.copy(caabadir + '/mecca/xmecca.log',          multirun_outputdir)
    shutil.copy(caabadir + '/gmake.log',                 multirun_outputdir)

    workdir = multirun_outputdir + '/workdir'
    # copy the necessary files from the current directory (caabadir) into workdir:
    # copy some directories completely:
    shutil.copytree('nml',            workdir+'/nml',            symlinks=True)
    shutil.copytree('input/multirun', workdir+'/input/multirun', symlinks=True)
    shutil.copytree('pycaaba',        workdir+'/pycaaba',        symlinks=True)
    #@NOTE: this is needed for mom chemistry (mabye we can only copy the relevant files?)
    shutil.copytree('skeleton',       workdir+'/skeleton',       symlinks=True)
    # copy some individual files:
    shutil.copy('caaba.exe', workdir)
    shutil.copy('jval.nml',  workdir) # maybe copy more *.nml files?
    # create some links:
    os.system('ln -s ' + os.readlink('caaba.nml') + ' ' + workdir + '/caaba.nml')
    os.system('ln -s ' + os.readlink('mecca.nml') + ' ' + workdir + '/mecca.nml')
    cwd = os.getcwd()
    os.chdir(workdir)
    
    if qexec:
        command = '/bin/tcsh -c "qexec python3 ./pycaaba/multirun.py %s %s"' % (
            full_nc_file, multirun_outputdir)
        print('Executing command:\n  %s' % (command))
        os.system(command)
    else:
        from multirun import multirun
        print('Executing multirun with:\n  multirun.complete(%s, %s)' % (
            full_nc_file, multirun_outputdir))
        multirun.complete(full_nc_file, multirun_outputdir, write_status=True)

    os.chdir(cwd)

def do_runs(icntrl_dict, rcntrl_dict, batchfile, integrators = list(),
            num_runs = 1, log_verbosity = 1, 
            full_nc_file='input/multirun/benchmark_small.nc', qexec=False):

    CAABA_NML_FILE = 'nml/caaba_benchmark.nml'
    os.system('ln -fs ' + CAABA_NML_FILE + ' caaba.nml')
    
    corencfile = corename(full_nc_file) # no dir, no suffix
    if (not os.path.isfile(full_nc_file)):
        sys.exit('ERROR: %s does not exist' % (full_nc_file))

    if not os.path.exists(caabadir + '/output'):
        os.mkdir(caabadir + '/output')

    outputdir = caabadir+'/output/benchmark'

    # if outputdir exists already, rename it, including its modification time:
    if (os.path.isdir(outputdir)):
        os.rename(outputdir, outputdir + '-' + time.strftime(
            '%Y-%m-%d-%H:%M:%S', time.localtime(os.path.getmtime(outputdir))))

    os.mkdir(outputdir)
    ncid_in = Dataset(full_nc_file, 'r')
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
        #print '%s\nSTART OF MULTIRUN (%d runs from %s"%s" axis of %s.nc)\n%s\n' % (
        #    HLINE, lines, unlimited, timedim, corencfile, HLINE)
    else:
        sys.exit('ERROR: Input file %s has neither "time" axis nor an unlimited dimension'
                 % (corencfile))

    ncid_in.close()

    #---------------------------------------------------------------------

    #NOTE: this is now correct as h211b is used as the default integrator
    prev_integrator = 'rosenbrock_posdef_h211b_qssa'

    for index in range(num_runs):

        if integrators is not None:
            if len(integrators) - 1 >= index:
                if integrators[index] != '_':
                    run_xmecca(batchfile, integrators[index])

        run_name = '%4.4d' % (index+1)

        # replace all parameters in the nml file with the input parameters

        icntrl_split_dict = dict( (key, icntrl_dict[key][index])
                                  for key, value in icntrl_dict.items())
        rcntrl_split_dict = dict( (key, rcntrl_dict[key][index])
                                  for key, value in rcntrl_dict.items())

        prepare_mecca_nml_file(icntrl_split_dict, rcntrl_split_dict)

        icntrl_copy = copy.deepcopy(icntrl_split_dict)
        rcntrl_copy = copy.deepcopy(rcntrl_split_dict)
        current_integ = prev_integrator if (
            integrators is None or len(integrators) - 1 < index or
            integrators[index] == '_') else integrators[index]

        # do the runs
        multirun_outputdir = outputdir + '/integr_settings_' + run_name
        qexec_multirun(full_nc_file, multirun_outputdir, qexec)
        
        # save integrator settings for plotting
        integr_settings = dict()
        integr_settings['integrator'] = current_integ
        integr_settings['icntrl'] = icntrl_split_dict
        integr_settings['rcntrl'] = rcntrl_split_dict
        save_obj(integr_settings, multirun_outputdir+'/integr_settings')

        if integrators is not None and len(integrators) - 1 >= index:
            if integrators[index] != '_':
                prev_integrator = integrators[index]

        #@NOTE(FS): delay next run to circumvent the rate limit on gaia
        # time.sleep(2)

    return outputdir

def parse_command_line_arguments():
    import argparse

    parser = argparse.ArgumentParser(
        description='make multiple runs with different integrators')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    parser.add_argument('-n', '--num_runs', type=int, dest='N', default=0,
                        help='number of runs')
    parser.add_argument('-combine', help='use all possibe combinations?',
                        action='store_true')
    parser.add_argument('-icntrl', action='append', nargs='*', metavar=('INUM','IVAL'),
                        help='define "icntrl(INUM) = IVAL", where IVAL is a ' +
                        'comma-separated list of values or "range <X> <Y> [<Z>]" ' +
                        '(<X> = start, <Y> = end, <Z> = step)')
    parser.add_argument('-rcntrl', action='append', nargs='*', metavar=('RNUM','RVAL'),
                        help='define "rcntrl(RNUM) = RVAL", where RVAL is a ' +
                        'comma-separated list of values or "range <X> <Y> [<Z>]" ' +
                        '(<X> = start, <Y> = end, <Z> = step)')
    parser.add_argument('-integ',
                        help='INTEG is a comma seperated list of integrators to use ' +
                        'for each run. Values can also be "_" if you do not wish to ' +
                        'compile a new integrator. The number of integrators does not ' +
                        'have to equal the number of runs; this is the same as using "_" ' +
                        'for all further runs', default = 'rosenbrock_posdef_h211b_qssa')
    parser.add_argument('-p', '--plot', action='store', type=str, nargs='?',
                        const='NO2,OH,O3',
                        help='plot the latest data set. No simulation is run if neither ' +
                        '-icntrl nor -rcntrl is specified. Use -d to specify directory ' +
                        'to load from.')
    parser.add_argument('-b', '--batchfile', action='store', type=str, default='mbl',
                        help='batchfile for xmecca')
    parser.add_argument('-i', '--input', action='store', type=str,
                        default='input/multirun/benchmark_small.nc',
                        help='path to netCDF4 input file.')
    parser.add_argument('-d', '--directory', action='store', type=str, default='benchmark',
                        help='Directory to load output from. Used when the only other ' +
                        'argument is -p / --plot.')
    parser.add_argument('-q','--qexec', help='use qexec to run multirun', action='store_true')
    parser.add_argument('-r', '--reference', help='prepend posdef to list of integrators', action='store_true')
    return parser.parse_args(), parser

# expand the "range <X> <Y> (<Step>)" syntax into actually usable arguments

def expand_range_argument(args):
    assert type(args) is list
    assert len(args) == 4 or len(args) == 5, 'incorrect number of range arguments'

    start = args[2]
    end = args[3]
    step = '1'

    if len(args) == 5:
        step = args[4]

    start = float(start) if isfloat(start) else int(start)
    end = float(end) if isfloat(end) else int(end)
    step = float(step) if isfloat(step) else int(step)

    if type(start) is float and step == 1:
        step = 1.0

    assert type(start) is type(end) is type(step), 'All range arguments must be of the same type'

    expanded = ''
    rng = list(range(start, end + step, step)) if type(start) is int else drange(start, end + step, step)
    for i in rng:
        expanded += str(i) + (',' if i < end else '')

    return expanded

def scatterplots(outputdir, plotspecies, LOGSCALE=True):

    if outputdir == caabadir + '/output/benchmark' and check_run_state() != 'ENDED':
        print('Error: Run has not ended yet, or did not end without errors')
        return
        
    def make_plot(ylabel, allxdata, allydata): # inner function
        from matplotlib import colors as mcolors

        xdata_min = np.min(allxdata, axis=1)  # min time for all integrators
        xdata_avg = np.mean(allxdata, axis=1) # avg time for all integrators
        xdata_max = np.max(allxdata, axis=1)  # max time for all integrators
        ydata_min = np.min(allydata, axis=1)
        ydata_avg = np.mean(allydata, axis=1)
        ydata_max = np.max(allydata, axis=1)
        # choose which data to plot:
        xdata = copy.deepcopy(xdata_avg)
        xdata /= xdata[0]
        ydata = ydata_max

        for val, idx in enumerate(ydata):
            if val == 0:
                print('[WARNING] Setting ' + str(idx+1) + ' has no y diff to reference integrator')

        # discard reference integrator (assumption is that index(ref) == 0)
        # xdata = xdata[1:]
        # ydata = ydata[1:]
        # color the default data point differently:
        #colors = np.array([('blue' if i == 0 else 'red') for i in range(0, len(xdata))])
        colors = list()
        # colors = 'red'
        #colors.append('blue')
        ax = viewport.next()

        #@TODO: TEMPORARY
        color_map_i3 = {
            '0': 'red',
            '1': 'green',
            '2': 'blue', 
            '3': 'orange',
            '4': 'darkmagenta', 
            '5': 'cyan'
        }
        color_map_i5 = {
            '0': 'red', 
            '1': 'green', 
            '2': 'blue'
        }
        color_map_i6 = {
            '0': 'red', 
            '1': 'green', 
            '2': 'blue'
        }
        color_map_r8 = {
            '0.0': 'red', 
            '1.0': 'cyan', 
            '2.0': 'green', 
            '3.0': 'darkmagenta',
            '4.0': 'salmon', 
            '5.0': 'orange', 
            '6.0': 'fuchsia', 
            '7.0': 'sienna', 
            '8.0': 'olivedrab'
        }
        color_map_r5 = {
            '0.0': 'red',
            '1E5': 'blue'
        }
        color_map_r9 = {
            '0.0': 'red',
            '2.0': 'green',
            '3.0': 'blue', 
            '4.0': 'orange',
        }
        color_map_r10 = {
            '1E-3': 'red', 
            '1E-2': 'blue', 
            '1E-1': 'green', 
            '1': 'orange', 
            '0': 'cyan'
        }

        # write logs for plot
        print('\n%-42s%s' % (
            ylabel, 'MIN       AVG     MAX                 MIN           AVG           MAX'))
        for i in range(len(integrators) - 1):
            plt.annotate(i + 1, (xdata[i], ydata[i]), fontsize=6)
            icntrl_s = ''
            colors.append(color_map_i5[icntrl[i]['5']])
            for key, value in icntrl[i].items():
                icntrl_s += 'icntrl %s: %s ' % (key, value)
            rcntrl_s = ''
            for key, value in rcntrl[i].items():
                rcntrl_s += 'rcntrl %s: %s ' % (key, value)
            s = ('%2d) %30s, t: %6g, %6g, %6g, diff: %12g, %12g, %12g, %s, %s' % (
                i + 1, integrators[i], xdata_min[i], xdata_avg[i], xdata_max[i],
                ydata_min[i], ydata_avg[i], ydata_max[i], icntrl_s, rcntrl_s))
            print(s)

        # WIP: error bars
        # https://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html
        # https://matplotlib.org/2.1.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        # plt.errorbar(xdata, ydata, xerr=[xdata-xdata_min,xdata_max-xdata],
        #              yerr=[ydata-ydata_min,ydata_max-ydata],
        #              fmt='o', linestyle='None')

        plt.scatter(xdata, ydata, color=colors, s=200, marker='.', linewidths=0)
        if (LOGSCALE):
            plt.yscale('log')
            # do not consider ydata values of 0 for calculating the logscale plot range:
            ydata = np.ma.masked_where(ydata==0, ydata)
            yrange = np.max(ydata) / np.min(ydata)
            if (yrange>1):
                ymin = np.min(ydata) / yrange**0.1
                ymax = np.max(ydata) * yrange**0.1
            else:
                ymin = np.min(ydata) / 1.001
                ymax = np.max(ydata) * 1.001
        else:
            yrange = np.max(ydata) - np.min(ydata)
            if (yrange>0):
                ymin = np.min(ydata) - yrange * 0.1
                ymax = np.max(ydata) + yrange * 0.1
            else:
                ymin = np.min(ydata) - 1e-20
                ymax = np.max(ydata) + 1e-20
        print ('%s: LOGSCALE=%s, ymin= %g, ymax= %g, plotmin= %g, plotmax= %g' %
               (ylabel, LOGSCALE, np.min(ydata), np.max(ydata), ymin, ymax))
        plt.ylim(ymin, ymax)
        plt.xlabel('User time')
        plt.ylabel(ylabel+' difference')

    plotspecies = plotspecies.split(',')
    run_dirs = sorted(glob(outputdir + '/integr_settings*'))

    #@NOTE: change REF if the integrator setting you want to use as the reference is not
    # at index 0
    # mom2000 runs: settings 0019 and 0073 can be used as references (so index 18 and 72)
    REF = 0 # reference integrator    

    # define integrator, icntrl, rcntrl, and user_time:
    integrators = list()
    icntrl = list()
    rcntrl = list()
    user_times = list()
    for run_dir in run_dirs:
        if not os.path.isfile(run_dir + '/timings.nc'):
            print("Warning: skipping %s, did not find file timings.nc" % (run_dir))
            continue
        integr_settings = load_obj(run_dir + '/integr_settings')
        integrators.append(integr_settings['integrator'])
        icntrl.append(integr_settings['icntrl'])
        rcntrl.append(integr_settings['rcntrl'])
        ncid_time = Dataset(run_dir + '/timings.nc', 'r')
        user_times.append(ncid_time.variables['user_time'][:])
        ncid_time.close()
    user_times = np.array(user_times) # convert list to numpy array
    # calculate diffs in final concentrations (c_end) between integrators:
    viewport.init(1, 3, outputdir+'/benchmark_scatterplots.pdf', 12, 8) # open pdf
    viewport.newpage()
    #plt.rcParams.update({ 'font.size': 8 })
    for species in plotspecies:
        alldata = list()
        for run_dir in run_dirs:
            if not os.path.isfile(run_dir + '/caaba_mecca_c_end.nc'):
                print("Warning: skipping %s, did not find file caaba_mecca_c_end.nc" % (run_dir))
                continue
            ncid_c_end = Dataset(run_dir + '/caaba_mecca_c_end.nc', 'r')
            alldata.append(ncid_c_end.variables[species][:])
            ncid_c_end.close()
        alldata = np.array(alldata) # convert list to numpy array
        diffs = np.abs(alldata[:,:] - alldata[REF,:]) # diff to reference integrator
        make_plot(species, user_times, diffs)
    viewport.exit() # close pdf

# pickle serialization utils:

def save_obj(obj, path):
    import pickle

    with open(path + '.pkl', 'wb') as f:
        #@NOTE: HIGHEST_PROTOCOL for binary format/faster serialization:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(path):
    import pickle

    with open(path + '.pkl', 'rb') as f:
        return pickle.load(f)

# dictionary utils (for "-combine" cmd line argument)

def tag_and_merge_dicts(a, b, tag_a, tag_b):
    c = dict()

    for key, value in a.items():
        c[key + tag_a] = value

    for key, value in b.items():
        c[key + tag_b] = value

    return c

def product_dict(**kwargs):
    import itertools
    keys = list(kwargs.keys())
    vals = list(kwargs.values())
    for instance in itertools.product(*vals):
        yield list(zip(keys, instance))

def produce_all_input_combinations(icntrl_dict, rcntrl_dict):
    import itertools

    merged_dict = tag_and_merge_dicts(icntrl_dict, rcntrl_dict, 'i', 'r')

    new_icntrl_dict = dict()
    new_rcntrl_dict = dict()

    for param_list in product_dict(**merged_dict):
        for input_param in param_list:

            if 'i' in input_param[0]:
                key = input_param[0].replace('i', '')
                if key not in new_icntrl_dict:
                    new_icntrl_dict[key]  = list()

                new_icntrl_dict[key].append(input_param[1])

            elif 'r' in input_param[0]:
                key = input_param[0].replace('r', '')
                if key not in new_rcntrl_dict:
                    new_rcntrl_dict[key] = list()

                new_rcntrl_dict[key].append(input_param[1])

    return new_icntrl_dict, new_rcntrl_dict


def get_status_output(*args, **kwargs):
    p = subprocess.Popen(*args, **kwargs)
    stdout, stderr = p.communicate()
    return p.returncode, stdout, stderr


def check_run_state():
    if not os.path.isdir(caabadir + '/output/benchmark'):
        return

    import getpass
    try:
        status, _, _ = get_status_output('squeue')
    except:
        status = 1

    if (status == 0):
        user_name = getpass.getuser().strip()
        job_list = subprocess.check_output(['squeue', '--user', user_name])
        lines = job_list.splitlines()
        if len(lines) > 1:
            return 'RUNNING'
        else:
            #@NOTE: this does not mean that there were no errors, but if the system gets stuck
            # and thinks that there are still jobs running, this is the only way to get it unstuck
            return 'ENDED'
    else:
        print('\n[WARNING]: Unable to guarantee that jobs running in parallel have finished\n')

    dirs = glob(caabadir + '/output/benchmark/integr_settings*')
    num_end_files = 0
    num_errors = 0
    for settings_dir in dirs:
        end_files = glob(settings_dir + '/status.log')
        if len(end_files) == 1:
            with open(end_files[0]) as f:
                if int(f.readlines()[0].split()[0]) == 0:
                    num_end_files += 1
                else:
                    num_errors += 1

    if num_end_files == len(dirs):
        if num_errors == 0:
            return 'ENDED'
        else:
            return 'ENDED_WITH_ERRORS'
    else:
        return 'RUNNING'

def main():
    # ensure that benchmark.py is started from main caaba/ directory:
    if (os.getcwd() != caabadir):
        sys.exit('ERROR: benchmark.py must be started from caaba directory')

    logfilename = caabadir + '/benchmark.log'
    tee.stdout_start(logfilename, append=False) # stdout

    # command line argument parsing:
    args, parser = parse_command_line_arguments()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        return

    print('\nCommand line arguments:')
    for arg in vars(args):
        print('%-12s: %s' % (arg, getattr(args, arg)))

    if args.plot and args.icntrl is None and args.rcntrl is None:
        if (args.directory[0:1] == '/'):
            outputdir = args.directory
        else:
            outputdir = caabadir + '/output/' + args.directory
        scatterplots(outputdir, args.plot)
        tee.stdout_stop()
        shutil.move(logfilename, outputdir+'/benchmark_scatterplots.log')
        return
    elif check_run_state() == 'RUNNING':
        print("Error: Previous run has not ended yet")
        return

    icntrl_dict = dict()
    rcntrl_dict = dict()

    verbosity = args.verbose
    integrators = list()

    if args.icntrl is not None:
        for icn in args.icntrl:
            if icn[1] == 'range':
                icntrl_dict[icn[0]] = expand_range_argument(icn).split(',')
            else:
                icntrl_dict[icn[0]] = icn[1].split(',')

    if args.rcntrl is not None:
        for rcn in args.rcntrl:
            if rcn[1] == 'range':
                rcntrl_dict[rcn[0]] = expand_range_argument(rcn).split(',')
            else:
                rcntrl_dict[rcn[0]] = rcn[1].split(',')

    if args.integ != '0':
        integrators = args.integ.split(',')

    # use every combination of input parameters.
    #@NOTE(F): This may produce undesired results when using different
    #integrators for each run.

    if args.combine:
        icntrl_dict, rcntrl_dict = produce_all_input_combinations(icntrl_dict, rcntrl_dict)

    if args.reference:
        hit_icntrl_3 = False
        integrators.insert(0, 'rosenbrock_posdef')
        for icntrl_val, param_list in icntrl_dict.items():
            if icntrl_val == '3':
                param_list.insert(0, '2')
                hit_icntrl_3 = True
            else:
                param_list.insert(0, '0')
        for param_list in rcntrl_dict.values():
            param_list.insert(0, '0')

        if not hit_icntrl_3:
            icntrl_dict['3'] = ['2']

    # determine number of runs
    #@NOTE(F): setting num_runs to sys.maxsize SHOULD be safe because it
    # is tested almost immediately after.

    num_runs = sys.maxsize
    if len(icntrl_dict) > 0:
        num_runs = min(min(list(map(len, iter(icntrl_dict.values())))), num_runs)
    if len(rcntrl_dict) > 0:
        num_runs = min(min(list(map(len, iter(rcntrl_dict.values())))), num_runs)

    if num_runs == sys.maxsize:
        num_runs = 0

    if args.N < num_runs and args.N > 0:
        num_runs = args.N

    print('\nPreparing %s runs, using these parameters: ' % (str(num_runs)))
    print('icntrl:')
    print(icntrl_dict)
    print('rcntrl:')
    print(rcntrl_dict)
    print('integrators:')
    print(integrators)
    print('batchfile:')
    if (os.path.isfile(caabadir + '/mecca/batch/' + args.batchfile + '.bat')):
        print(args.batchfile + '.bat')
    else:
        sys.exit('ERROR: Cannot find ' + args.batchfile + '.bat')
    print('netCDF input file:')
    print(args.input)

    outputdir = do_runs(icntrl_dict, rcntrl_dict, args.batchfile,
                        integrators, num_runs = num_runs,
                        log_verbosity = verbosity, full_nc_file = args.input,
                        qexec=args.qexec)
    if outputdir is None:
        print("[Error]: aborting benchmark/multirun")
        return

    if args.plot:
        print()
        scatterplots(outputdir, args.plot)

    tee.stdout_stop()
    shutil.move(logfilename, outputdir)

if __name__ == '__main__':
    main()
