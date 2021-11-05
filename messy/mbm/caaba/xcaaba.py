#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-07 15:05:38 sander>

# xcaaba.py = eXecute CAABA
# Rolf Sander, 2018

##############################################################################

import os, sys, shutil
CAABADIR = os.path.realpath(os.path.dirname(__file__))
sys.path.append(os.path.realpath(CAABADIR+'/pycaaba'))
import re # regexp
import subprocess
import datetime
from glob import glob
from pyteetime import tee # from pycaaba
from rstools import HLINE, HLINE2, cat, grep_i, runcmd, tail # from pycaaba
from caabatools import split_caaba_mecca_nc
from multirun import multirun
from netCDF4 import Dataset

TMPFILE = 'tmp.txt'
LOGFILENAME_TMP = 'xcaaba.logfile' # suffix .log could be deleted by gmake clean
LOGFILENAME = os.path.splitext(LOGFILENAME_TMP)[0] + '.log'
HOME = os.environ['HOME']
HOST = os.environ['HOST']
USER = os.environ['USER']
DATE_NOW = datetime.datetime.now().strftime('%Y-%m-%d')
TIME_NOW = datetime.datetime.now().strftime('%H:%M:%S')

##############################################################################

def initial_checks():
    
    # ensure that $HOME/tmp directory exists:
    tmpdir = HOME+'/tmp'
    if (not os.path.isdir(tmpdir)):
        os.mkdir(tmpdir)

    # check that the input directory is okay:
    if (not os.path.isdir('input')): # if it is not a directory
        if (USER == 'sander'):
            os.symlink(HOME+'/e2/caaba_input', 'input')
        else:
            sys.exit('ERROR: The CAABA input directory is missing.')
    # check that the testsuite directory is okay:
    if (not os.path.isdir('testsuite')): # if it is not a directory
        if (USER == 'sander'):
            os.symlink(HOME+'/e2/caaba_testsuite', 'testsuite')
        else:
            sys.exit('ERROR: The CAABA testsuite directory is missing.')

##############################################################################

def maybe_xmecca():
    l_xmecca = False # default is not to start l_xmecca
    batfile = ''
    if len(sys.argv)>1:
        batfile = sys.argv[1]
        l_xmecca = True
    else:
        print('\nYou can create a new chemical mechanism with xmecca now.')
        print('This is necessary, if:')
        print('- you have changed any of the *.eqn, *.spc, or *.kpp files')
        print('- you have changed your *.bat batch file')
        inputstring = input('Start xmecca? [y|n|q|default=n] ')
        if (inputstring == 'q'): sys.exit(0)
        if (inputstring == 'y'):
            l_xmecca = True

    if (not l_xmecca):
        # check that kpp-generated f90 files are newer than kpp files:
        gaseqntime     = os.path.getmtime('mecca/gas.eqn')
        gasspctime     = os.path.getmtime('mecca/gas.spc')
        aqueouseqntime = os.path.getmtime('mecca/aqueous.eqn')
        aqueousspctime = os.path.getmtime('mecca/aqueous.spc')
        kpptime        = os.path.getmtime('mecca/messy_mecca_kpp.kpp')
        f90time        = os.path.getmtime('messy_mecca_kpp.f90')
        l_xmeccayes = False
        if (gaseqntime>f90time):
            print('WARNING: gas.eqn is newer than kpp-generated f90 files')
            l_xmeccayes = True
        if (gasspctime>f90time):
            print('WARNING: gas.spc is newer than kpp-generated f90 files')
            l_xmeccayes = True
        if (aqueouseqntime>f90time):
            print('WARNING: aqueous.eqn is newer than kpp-generated f90 files')
            l_xmeccayes = True
        if (aqueousspctime>f90time):
            print('WARNING: aqueous.spc is newer than kpp-generated f90 files')
            l_xmeccayes = True
        if (kpptime>f90time):
            print('WARNING: messy_mecca_kpp.kpp is newer than kpp-generated f90 files')
            l_xmeccayes = True
        if (l_xmeccayes):
            print('It is strongly suggested to create new f90 files via xmecca.')
            inputstring = input('Start xmecca? [y|n|q|default=y] ')
            if (inputstring == 'q'): sys.exit(0)
            if (inputstring != 'n'):
                l_xmecca = True

    if (l_xmecca):
        print('\nstarting xmecca.')
        olddir = os.getcwd()
        os.chdir(CAABADIR+'/mecca') # cd to MECCA directory
        cmd = './xmecca '+batfile
        exitstatus = subprocess.call(cmd, shell=True)
        os.chdir(olddir)
        if (exitstatus == 0):
            print('xmecca has finished successfully')
        else:
            print('\nStopping xcaaba here because xmecca returned exit status: %d' % (
                exitstatus))
            sys.exit(1)
        print('\n%s' % (HLINE))

    return l_xmecca

##############################################################################

def compile_code(l_xmecca):
    s_string = 'Start from scratch ("gmake clean", then compile)'
    c_string = 'Compile recently changed files with "gmake"'
    r_string = 'Run existing executable'
    print('\nChoose an option [default=c]:')
    print('s = %s' % (s_string))
    print('c = %s' % (c_string))
    if (not l_xmecca):
        print('r = %s' % (r_string))
    print('q = Quit')
    inputstring = input('')
    if (inputstring == 'q'):
        sys.exit(0)
    elif (inputstring == 'r'):
        if (l_xmecca):
            print('ERROR: You must choose option c')
            print('       after creating a new chemical mechanism with xmecca')
            sys.exit(1)
        else:
            option = 'r'
            print('You have chosen: r = %s' % (r_string))
    elif (inputstring == 's'):
        option = 's'
        print('You have chosen: s = %s' % (s_string))
    else:
        option = 'c'
        print('You have chosen: c = %s' % (c_string))
    print('\n%s' % (HLINE))
    if (option == 's'):
        print('\ngmake clean')
        subprocess.call('gmake clean', shell=True)
    if ((option == 'c') or (option == 's')):
        print('\ngmake validate')
        exitstatus = subprocess.call('gmake validate', shell=True)
        print('exit status from "gmake validate" is: %s' % (exitstatus))
        if (exitstatus != 0): sys.exit(1)
        print('\n%s' % (HLINE))
        print('\ngmake (writing output to gmake.log)')
        cmd = 'gmake -j' # the option '-j' allows simultaneous jobs
        # activate next line (instead of the one above) to show more debugging info:
        # cmd = 'gmake -j --debug=basic'
        CMDLOGFILE = open('gmake.log','w+', 1)
        exitstatus = subprocess.call(cmd, stdout=CMDLOGFILE, stderr=CMDLOGFILE, shell=True)
        CMDLOGFILE.close()
        print('exit status from "gmake" is: %s' % (exitstatus))
        if (exitstatus != 0):
            tail('gmake.log',20)
            print('\nOnly the last 20 lines of the compiler output are shown.')
            print('For further details, check gmake.log!\n')
            sys.exit(1)
        print('\n%s' % (HLINE))

##############################################################################

def select_nml():
    defaultnml = os.readlink('caaba.nml')
    allfiles = sorted(glob('nml/*.nml'))
    print('\nChoose a namelist file from the nml/ directory:')
    for i, nmlfile in enumerate(allfiles): # list all possibilities
        print('%2d) %s' % (i+1, os.path.basename(nmlfile)))
    print(' q) quit')
    inputstring = input('The default is %s (same as last time)\n' % (
        os.path.basename(defaultnml)))
    if (inputstring == 'q'): sys.exit(1)
    try:
        nmlindex = int(inputstring)-1
    except ValueError:
        nmlindex = -1
    if((nmlindex < len(allfiles)) and (nmlindex >= 0)):
        nmlfile = allfiles[nmlindex]
    else:
        nmlfile = defaultnml
    if (os.path.isfile('caaba.nml')): os.remove('caaba.nml')
    os.symlink(nmlfile, 'caaba.nml')
    print('%s\n\nThe active contents of caaba.nml is:' % (HLINE))
    os.system("sed -ne '/^&CAABA[ ]*$/,/^\//p' caaba.nml | " +
              "grep -v '^!' | grep -v '^ *$' > %s" % (TMPFILE))
    cat(TMPFILE)
    print('%s\n\nThe active contents of mecca.nml is:' % (HLINE))
    os.system("sed -ne '/^&CTRL/,/^\//p' mecca.nml | " +
              " grep -v '^!' | grep -v '^ *$' > %s" % (TMPFILE))
    cat(TMPFILE)
    print('%s\n\nBefore you continue, ensure that the selected namelist' % (HLINE))
    print(os.path.basename(nmlfile))
    print('is consistent with the selected chemistry mechanism!')
    # ----------------------
    # mecca.nml:
    if (os.path.isfile('mecca.nml')): os.remove('mecca.nml')
    os.symlink('nml/mecca_default.nml', 'mecca.nml')

##############################################################################

def check_if_monte_carlo():
    l_montecarlo = grep_i('REQ_MCFCT *= *\.TRUE\.', 'messy_mecca_kpp*.f90')
    if (l_montecarlo):
        inputstring = input(
            '\nRun Monte-Carlo simulations with CAABA? [y|n|q|default=y] ')
        if (inputstring == 'q'): sys.exit(0)
        if (inputstring != 'n'):
            print('\nStarting montecarlo.py...')
            subprocess.call('./pycaaba/montecarlo.py', shell=True)
            finalize()
            sys.exit(0)

##############################################################################

def run_caaba_exe():
    print('\nRun CAABA/MECCA?')
    print('y = yes (default)')
    print('m = multirun')
    print('n = no')
    print('q = quit')
    inputstring = input('')
    print(HLINE)
    if (inputstring == 'q'): sys.exit(0)
    if (inputstring == 'm'): # multirun
        print('\nChoose an input file from the input/multirun/ directory:')
        print('[q|number|default=loop over all files]')
        allfiles = sorted(glob('input/multirun/*.nc'))
        for i, ncfile in enumerate(allfiles): # list all possibilities
            print('%2d) %s' % (i+1, os.path.basename(ncfile)))
        inputstring = input('')
        if (inputstring == 'q'): sys.exit(1)
        try:
            ncindex = int(inputstring)-1
        except ValueError:
            ncindex = -1
        if((ncindex < len(allfiles)) and (ncindex >= 0)):
            ncfile = allfiles[ncindex]
            multirun.complete(ncfile)
        else:
            for ncfile in allfiles:
                multirun.complete(ncfile)
        finalize()
        sys.exit(0)
    if (inputstring != 'n'):
        # remove old files, if any:
        list(map(os.remove, glob('caaba_*.nc')))
        print('Starting CAABA. Please wait...\n')
        # run the CAABA/MECCA box model (script for line buffering):
        os.system("script -q -c '\\time -p ./caaba.exe' %s" % (TMPFILE))
        # convert to unix format and add to logfile:
        with open(TMPFILE, 'r') as f:
            print(f.read().replace('\r', ''), file=LOGFILE)
        # set exitstatus to the value in the file status.log:
        with open('status.log') as f: exitstatus = int(f.read())
        if (exitstatus == 0):
            print('\ncaaba.exe has finished successfully')
        else:
            print("\nERROR: exit status from 'caaba.exe' is: %s" % (exitstatus))
            sys.exit(1)

##############################################################################

def save_model_output():
    if not os.path.exists(CAABADIR + '/output'):
        os.mkdir(CAABADIR + '/output')
    defaultoutputdir = 'output/%s-%s' % (DATE_NOW, TIME_NOW)
    print('\nSave the output and model code?')
    print('Choose an option or type a directory name:')
    print('y (default) = yes, save in %s' % (defaultoutputdir))
    print('n           = no (output stays in current directory)')
    print('q           = quit')
    print('<dirname>   = save in output/<dirname>')
    inputstring = input('')
    if (inputstring == 'q'):
        sys.exit(0)
    elif (inputstring == '' or inputstring == 'y'):
        outputdir = defaultoutputdir
    elif (inputstring == 'n'):
        outputdir = False
    else:
        # remove unsuitable characters:
        myoutputdir = re.sub('[^-A-Za-z0-9:_+=.]', '', inputstring)
        if (inputstring != myoutputdir):
            print('Unsuitable characters (e.g., from pressing arrow keys) have been')
            print('removed from the name of the output directory. The new name is:')
            print(myoutputdir)
        outputdir = 'output/%s' % (myoutputdir)
        # ensure directory doesn't exist yet:
        if (os.path.isdir(outputdir)):
            print('Directory %s already exists.' % (outputdir))
            print('Using default directory instead.')
            outputdir = defaultoutputdir
        # confirmation:
        print('Name of output directory = %s' % (outputdir))
        inputstring2 = input('Confirm to save output now [y|n|q, default=y]\n')
        if (inputstring2 == 'q'): sys.exit(0)
        if (inputstring2 == 'n'): outputdir = False
    if (outputdir):
        print('Creating zip file of caaba model code. Please wait...')
        os.mkdir(outputdir)
        os.system('gmake zip > %s' % (TMPFILE))
        with open(TMPFILE, 'r') as f: print(f.read(), file=LOGFILE) # add to logfile
        shutil.move('%s.zip' % (os.path.basename(CAABADIR)), outputdir)
        os.system('cp -p *.nc *.dat pycaaba/_*.py %s' % (outputdir))
        print('Model code and output have been saved to %s:' % (outputdir))
        for thefile in sorted(os.listdir(outputdir)): print(thefile)
    else:
        print('Model code and output have not been saved.')
    return outputdir

##############################################################################

def visualize(outputdir):
    from caabaplot import makeplots_xxxg
    inputstring = input('\nPlot all species with caabaplot.py? [y|n, default=n]\n')
    if (inputstring == 'y'):
        if (outputdir):
            modelruns = [[outputdir, os.path.basename(outputdir)]]
        else:
            modelruns = [['.', 'latest run']]
        makeplots_xxxg(modelruns)

##############################################################################

def finalize(outputdir=False):
    if (os.path.isfile(TMPFILE)): os.remove(TMPFILE)
    print('\n%s\n%48s\n%s\n' % (HLINE, 'xcaaba has finished', HLINE))
    tee.stdout_stop()
    os.rename(LOGFILENAME_TMP, LOGFILENAME) # final suffix .log
    if outputdir:
        os.system('cp -p %s %s' % (LOGFILENAME, outputdir))
    print('Log output is now available in %s' % (LOGFILENAME))

##############################################################################

if __name__ == '__main__':

    LOGFILE = tee.stdout_start(LOGFILENAME_TMP, append=False) # stdout
    print('xcaaba.py was started on %s at %s by user %s on machine %s' % (
        DATE_NOW, TIME_NOW, USER, HOST), file=LOGFILE)
    initial_checks()
    l_xmecca = maybe_xmecca()
    compile_code(l_xmecca)
    select_nml()
    check_if_monte_carlo()
    # run the CAABA/MECCA model:
    run_caaba_exe()
    # check if there are any reaction rates RR* in caaba_mecca.nc:
    split_caaba_mecca_nc()
    outputdir = save_model_output()
    # plot results:
    visualize(outputdir)
    finalize(outputdir)
