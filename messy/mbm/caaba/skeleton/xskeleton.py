#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-07 14:08:19 sander>

# xskeleton: execute mechanism reduction to obtain a skeletal mechanism
# Rolf Sander, 2016-2017

##############################################################################

import os, sys, shutil
caabadir = os.path.realpath(os.path.dirname(__file__)+'/..')
sys.path.append(os.path.realpath(caabadir+'/pycaaba'))
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#plt.rcParams.update({'figure.max_open_warning': 0}) # 'More than 20 figures'
from matplotlib.backends.backend_pdf import PdfPages
import time
import datetime
import numpy as np
from glob import glob
from mecca import mecca # from pycaaba
from caabaplot import caabaplot # from pycaaba
from pyteetime import tee # from pycaaba
from rstools import HLINE, HLINE2, runcmd # from pycaaba

KPPMODE  = '// -*- kpp -*- kpp mode for emacs'
DONTEDIT = '// created automatically by %s at %s, DO NOT EDIT!' % (
    sys.argv[0], datetime.datetime.now().strftime("%Y-%m-%d-%H:%M:%S"))
# Read info about targets from file:
# python3 needs this ugly workaround for dtype. Otherwise the "bytes" type is created
targetdata = np.genfromtxt('targets.txt', dtype='U80,float,float', comments='#')
#targetdata = np.genfromtxt('targets.txt', dtype=None, comments='#')
outputdir = caabadir+'/output/skeleton'
if (not os.path.isdir(outputdir)):
    os.mkdir(outputdir)
global epslist

# run control default values (usually overwritten by ctrl/*.py):
samplepointfile   = ''   # empty, create only plots
eps0              = 5E-4 # start value for epsilon_ep
eps_increase      = 1.2  # factor for increasing epsilon_ep
plot_delta_skel   = 1    # 0=no plots, 1=all
plot_targets      = 1    # 0=no plots, 1=all
plot_samplepoints = 1    # 0=no plots, 1=all, 2=some

##############################################################################

def cleanup():
    if (os.path.isdir(outputdir)):
        os.rename(outputdir, outputdir + '-' + time.strftime(
            '%Y-%m-%d-%H:%M:%S', time.localtime(os.path.getmtime(outputdir))))
    os.mkdir(outputdir)
    
##############################################################################

def get_workdir(skelnum, mkdir=True):
    if (skelnum):
        # workdir contains skelnum with 3 digits:
        workdir = '%s/skeleton_%3.3d' % (outputdir, skelnum)
    else:
        workdir = '%s/fullmech' % (outputdir)
    if (full_calc and mkdir):
        os.mkdir(workdir)
    return workdir

##############################################################################

def show_reaction(rxn):
    reagents = ''
    products = ''
    for spc_num, spc_stoic in enumerate(StoichNum[rxn]): # loop over species
        if (spc_stoic < 0):
            reagents += ' + %g %s' % (-spc_stoic, oicdata[spc_num][1])
        if (spc_stoic > 0):
            products += ' + %g %s' % (spc_stoic, oicdata[spc_num][1])
    return '%-10s %s -> %s' % ('<'+EQN_TAGS[rxn]+'>', reagents, products)

##############################################################################

def create_skeletal_mechanism(eps):
    global del_rxn
    MECHLOGFILE = open(workdir+'/mechanism.log','w+', 1) # 1=line-buffered
    print('%s\neps = %g\n%s' % (HLINE, eps, HLINE), file=MECHLOGFILE)
    # Create empty list of reactions to delete:
    del_rxn = [False] * len(StoichNum)
    N_var_skel = 0
    # Create mechanism including all species with OIC>eps:
    for num, (oic, name) in enumerate(oicdata): # loop over species
        if (oic > eps): # keep!
            print('\nKEEP   %4d %s %15g' % (num+1, name, oic), file=MECHLOGFILE)
            N_var_skel += 1
        else: # delete!
            print('\nDELETE %4d %s %15g' % (num+1, name, oic), file=MECHLOGFILE)
            # Find rxns that contain this species:
            for rxn in range(N_rxns): # loop over reactions
                if (StoichNum[rxn][num] != 0):
                    if (StoichNum[rxn][num] < 0):
                        del_rxn[rxn] = True # mark for deletion
                        print('  '+show_reaction(rxn), file=MECHLOGFILE)
                        print('    delete because reagent %10g %15s' % (
                                StoichNum[rxn][num], name), file=MECHLOGFILE)
                    if (StoichNum[rxn][num] > 0):
                        del_rxn[rxn] = True # mark for deletion
                        print('  '+show_reaction(rxn), file=MECHLOGFILE)
                        print('    delete because product %10g %15s' % (
                                StoichNum[rxn][num], name), file=MECHLOGFILE)
    print('\n%s\nSUMMARY OF DELETED REACTIONS\n%s:' % (
        HLINE, HLINE), file=MECHLOGFILE)
    for rxn in range(N_rxns): # loop over reactions
        if (del_rxn[rxn]):
            print('DELETE: '+show_reaction(rxn), file=MECHLOGFILE)
    print('\n%s\nSUMMARY OF KEPT REACTIONS\n%s:' % (HLINE, HLINE), file=MECHLOGFILE)
    for rxn in range(N_rxns): # loop over reactions
        if (not del_rxn[rxn]):
            print('KEEP: '+show_reaction(rxn), file=MECHLOGFILE)
    MECHLOGFILE.close()
    skelinfo = 'nvar = %d/%d, nreact = %d/%d, eps = %g' % (
        N_var_skel, N_var_full, del_rxn.count(False), N_rxns, eps)
    print(skelinfo)
    epslist.append(eps)
    # Create rpl file:
    create_skeleton_rpl(del_rxn, skelinfo)
    del_rxnlist.append(del_rxn)

##############################################################################

def create_skeleton_rpl(delrxn, skelinfo=None):
    # A rpl file that deletes several reactions is used to create the
    # skeletal mechanism. Initially, an empty rpl file creates the full
    # mechanism.
    rplfilename = caabadir+'/skeleton/skeleton.rpl'
    RPLFILE = open(rplfilename,'w+')
    print(KPPMODE + '\n' + DONTEDIT, file=RPLFILE)
    print('// ctrlfile = %s' % (ctrlfile), file=RPLFILE)
    print('//\n// target               abstol     reltol', file=RPLFILE)
    for num, (target, abstol, reltol) in enumerate(targetdata):
        print('// %-15s %10G %8G ' % (target, abstol, reltol), file=RPLFILE)
    print('//\n// sample point file = %s.nc' % (samplepointfile), file=RPLFILE)
    print('// epsilon0          = %G' % (eps0), file=RPLFILE)
    print('// eps_increase      = %G' % (eps_increase), file=RPLFILE)    
    if skelinfo:
        print('// skeleton number = %s, %s' % (skelnum, skelinfo), file=RPLFILE)
    else:
        print('// full mechanism', file=RPLFILE)
    if (delrxn):
        for rxn, delete in enumerate(delrxn): # loop over reactions
            if (delete):
                print('#REPLACE %-10s' % ('<'+EQN_TAGS[rxn]+'>'), file=RPLFILE)
                print('#ENDREPLACE', file=RPLFILE)
    RPLFILE.close()
    # Save *.rpl file in workdir:
    os.system('cp -p ' + rplfilename + ' ' + workdir)
    print('replacement file skeleton.rpl is in:\n  %s' % (workdir))

##############################################################################

def caaba_multirun():
    from multirun import multirun
    if (full_calc):
        olddir = os.getcwd()
        os.chdir(caabadir+'/mecca') # cd to MECCA directory
        # Create mechanism:
        runcmd('xmecca skeleton', workdir+'/xmecca.log')
        os.chdir(caabadir) # cd to CAABA base directory
        # Use the namelist file caaba_skeleton.nml:
        os.system('ln -fs nml/caaba_skeleton.nml caaba.nml')
        # Compile caaba:
        runcmd('gmake', workdir+'/gmake.log')
        # Run CAABA for all sample points:
        multirun.complete('skeleton/samplepoints/'+samplepointfile+'.nc')
        # Move output to final destination:
        os.rename('output/multirun/' + samplepointfile, workdir+'/multirun')
        os.chdir(olddir)
    else:
        print('Skipping CAABA multirun because full_calc=0')

##############################################################################

def get_samplepointnames():
    fullsamplepointnames = sorted(glob(outputdir+'/fullmech/multirun/runs/*'))
    samplepointnames = list(map(os.path.basename, fullsamplepointnames))
    N_samplepoints = len(samplepointnames)
    #print 'samplepoints   = ', samplepointnames
    print('N_samplepoints = ', N_samplepoints)
    return samplepointnames, N_samplepoints

##############################################################################

def calc_oic():
    global oicdata, N_var_full, StoichNum, N_rxns, EQN_TAGS, EQN_NAMES
    if (full_calc):
        print('\nCalculate DICs and OICs for full mechanism.')
        # Compile skeleton:
        runcmd('gmake', workdir+'/gmake_oic.log')
        # Run oic:
        runcmd('./oic.exe', workdir+'/oic.log')
        os.system('cp -p OIC.dat StoichNum.dat EQN_*.dat ' + workdir)
    else:
        print('Skipping calculation of OIC with oic.exe because full_calc=0')

    if (os.path.basename(workdir) == 'fullmech'):
        # Load OIC data from oic.exe into a numpy structured array:
        # http://docs.scipy.org/doc/numpy/user/basics.rec.html
        oicdata = np.genfromtxt('OIC.dat', dtype=None)
        N_var_full = len(oicdata) # number of variable species
        # Load StoichNum from oic.exe:
        StoichNum = np.genfromtxt('StoichNum.dat')
        N_rxns = len(StoichNum) # number of reactions
        # Load EQN_TAGS from oic.exe:
        EQN_TAGS = np.genfromtxt('EQN_TAGS.dat', dtype='str')
        # Load EQN_NAMES from oic.exe:
        EQN_NAMES = [line.rstrip() for line in open('EQN_NAMES.dat')]

##############################################################################

def save_rates(delrxn):
    rates = np.zeros((N_samplepoints,len(delrxn)))
    keep_rxn = np.invert(np.array(delrxn))
    # sample point loop:
    for samplepointnum, samplepoint in enumerate(samplepointnames):
        ratesfile = workdir+'/multirun/runs/'+samplepoint+'/caaba_mecca_a_end.dat'
        rates0 = np.genfromtxt(ratesfile, dtype=None)
        rates[samplepointnum,keep_rxn] = rates0
    return rates

##############################################################################

def calc_error():
    olddir = os.getcwd()
    os.chdir(workdir)
    ERRORFILE = open('skel_error.dat','w+', 1) # 1=line-buffered
    # delta_samplepoints is a 2D array (list of lists) containing the
    # errors of the current skeletal mechanism compared to the full
    # mechanism for all targets and for all sample points:
    delta_samplepoints = []
    for samplepoint in samplepointnames: # sample point loop
        print('sample point: %s' % (samplepoint), file=ERRORFILE)
        filename_part2 = 'multirun/runs/'+samplepoint+'/caaba_mecca_c_end.nc'
        # ncfile_full = NetCDFFile('../fullmech/'+filename_part2)
        # ncfile_skel = NetCDFFile(filename_part2)
        ncfile_full = Dataset('../fullmech/'+filename_part2)
        ncfile_skel = Dataset(filename_part2)
        # delta_targets is a 1D array (list) containing the errors of
        # the current skeletal mechanism compared to the full mechanism
        # for all targets:
        delta_targets = [None] * len(targetdata)
        print('target           abstol      reltol    ' + \
            'mixrat_skel    mixrat_full     err/reltol', file=ERRORFILE)
        # target loop:
        for num, (target, abstol, reltol) in enumerate(targetdata):
            # read one number in a 4D array:
            mixrat_full = ncfile_full.variables[target][0][0][0][0]
            mixrat_skel = ncfile_skel.variables[target][0][0][0][0]
            delta_targets[num] = abs(max(mixrat_skel,abstol)/ \
              max(mixrat_full,abstol)-1) / reltol
            print('%-15s %10G %8G %14G %14G %14G' % (
                target, abstol, reltol, mixrat_skel, mixrat_full,
                delta_targets[num]), file=ERRORFILE)
        delta_samplepoints.append(delta_targets)
        ncfile_full.close()
        ncfile_skel.close()
    ERRORFILE.close()
    os.chdir(olddir)
    return delta_samplepoints

##############################################################################

def analyze_results():
    global all_rates
    delta_samplepoints = calc_error() # calculate error relative to full mechanism
    delta_skel = np.amax(delta_samplepoints) # max error over all targets and sample points
    # add rates in current skeletal mechanism to all_rates:
    all_rates = np.dstack((all_rates, save_rates(del_rxn)))
    if (delta_skel > 1.):
        RATESFILE = open(outputdir+'/rates.dat','w+', 1) # 1=line-buffered
        # difference of current to previous skeletal mechanism:
        diff = all_rates[:,:,-1] - all_rates[:,:,-2]
        print('\n%s\nError of skeletal mechanism has become too big.\n%s' % (HLINE, HLINE))
        print('List of sample points with delta_skel>1:')
        for samplepointnum in range(N_samplepoints): # sample point loop
            idx = np.argsort(abs(diff[samplepointnum,:]))
            # check if this is a problem sample point:
            if (max(delta_samplepoints[samplepointnum]) > 1.):
                print(HLINE+'\n\n*** sample point: %04d\n' % (samplepointnum+1), file=RATESFILE)
                for targetnum, delta_target in enumerate(delta_samplepoints[samplepointnum]):
                    if (delta_target > 1.):
                        print('Sample point: %04d, delta_skel: %10G, target: %s' % (
                            samplepointnum+1, delta_target, targetdata[targetnum][0]))
                    print('delta_skel: %10G, target: %s' % (
                        delta_target, targetdata[targetnum][0]), file=RATESFILE)
                print('\nColumn 1: Reaction rate in previous '+ \
                    '(s%03d) skeletal mechanism [cm-3 s-1]' % (skelnum-1), file=RATESFILE)
                print('Column 2: Is reaction also included in current ' \
                    '(s%03d) skeletal mechanism?' % (skelnum), file=RATESFILE)
                print('Column 3: Difference (s%03d-s%03d) of current ' \
                    'minus previous reaction rate\n' % (skelnum, skelnum-1), file=RATESFILE)
                print('        s%03d  s%03d         diff reaction' % (
                    skelnum-1, skelnum), file=RATESFILE)
                for i in range(idx.shape[0]):
                    x = idx[-i-1]
                    if (not del_rxnlist[-2][x]): # if rxn was in previous skeleton
                        print('%12G %5s %12G %-10s %s' % (
                            all_rates[samplepointnum,x,-2], not del_rxnlist[-1][x],
                            diff[samplepointnum,x], '<'+EQN_TAGS[x]+'>', EQN_NAMES[x]), file=RATESFILE)
                print(file=RATESFILE)
        print(HLINE, file=RATESFILE)
        RATESFILE.close()
        # copy previous replacement file to outputdir:
        os.system('cp -p %s/skeleton.rpl %s/' % (
            get_workdir(skelnum-1, mkdir=False), outputdir))
        # copy ctrl file and targets.txt to outputdir:
        os.system('cp -p %s/skeleton/%s %s/' % (caabadir, ctrlfile, outputdir))
        os.system('cp -p %s/skeleton/targets.txt %s/' % (caabadir, outputdir))
        
    return delta_skel, delta_samplepoints

##############################################################################

def list_species(oicdata, epslist):
    SPECFILE = open(outputdir+'/species.dat','w+', 1) # 1=line-buffered
    print(HLINE, file=SPECFILE)
    for epsnum,eps in enumerate(epslist):
        if (epsnum == 0):
            print('full', end=' ', file=SPECFILE)
        else:
            print('s%03d' % (epsnum), end=' ', file=SPECFILE)
    print('OIC          species', file=SPECFILE)
    print(HLINE, file=SPECFILE)
    for oic, species in np.sort(oicdata):
        for epsnum,eps in enumerate(epslist):
            if (oic < eps):
                print('    ', end=' ', file=SPECFILE)
            else:
                if (epsnum == 0):
                    print('full', end=' ', file=SPECFILE)
                else:
                    print('s%03d' % (epsnum), end=' ', file=SPECFILE)
        print('%12E %s' % (oic, species), file=SPECFILE)
    print(HLINE, file=SPECFILE)
    SPECFILE.close()

##############################################################################

def list_reactions(del_rxnlist):
    RXNSFILE = open(outputdir+'/reactions.dat','w+', 1) # 1=line-buffered
    print(HLINE, file=RXNSFILE)
    for rxn in range(N_rxns): # loop over reactions
        for mechnum, delrxn in enumerate(del_rxnlist):
            if (delrxn[rxn]):
                print('    ', end=' ', file=RXNSFILE)
            else:
                if (mechnum == 0):
                    print('full', end=' ', file=RXNSFILE)
                else:
                    print('s%03d' % (mechnum), end=' ', file=RXNSFILE)
        # use either EQN_TAGS+EQN_NAMES or show_reaction:
        print('%-10s %s' % ('<'+EQN_TAGS[rxn]+'>', EQN_NAMES[rxn]), file=RXNSFILE)
        #print >> RXNSFILE, show_reaction(rxn)
    print(HLINE, file=RXNSFILE)
    RXNSFILE.close()

##############################################################################

def finalize(info):
    print(HLINE)
    print('Summary of results: %d targets, %d sample points, %d skeletal mechanisms' % (
        info[2], info[1], info[0]))
    print(HLINE)
    print('Sample points                  acro samplepoints/'+samplepointfile+'.pdf')
    print('Go to output directory:        cd ../output/skeleton')
    print('Logfile:                       e xskeleton.log')
    print('List of species:               e species.dat')
    print('List of reactions:             e reactions.dat')
    print('List of rates:                 e rates.dat')
    print('Best replacement file:         e skeleton.rpl')
    print('Errors of skeletal mechanisms: e skeleton_*/skel_error.dat')
    print('Plots of errors:               acro delta_skel.pdf')
    print('Plots of target species:       acro targets.pdf')
    print('Plots of sample points:        acro samplepoint_*.pdf')
    print(HLINE)
    tee.stdout_stop()
    os.system('cp -p xskeleton.log %s/' % (outputdir)) # copy logfile to outputdir

##############################################################################

def make_target_plots(plot_targets):
    from viewport import viewport # from pycaaba
    if (not plot_targets): return
    print('\nPlotting these skeletal mechanisms:\n', all_skel, '\n')
    viewport.init(4, 4, outputdir+'/targets.pdf', 17, 8) # open pdf
    print(HLINE)
    for num, (target, abstol, reltol) in enumerate(targetdata): # target loop
        print('Plotting target %-15s' % (target))
        for samplepoint in samplepointnames: # sample point loop
            caabaplot.plot_0d(
                modelruns = [[outputdir+'/'+skel+'/multirun/runs/'+samplepoint, skel]
                               for skel in all_skel],
                species   = target,
                pagetitle = 'Target species: '+target,
                plottitle = 'Sample point: '+samplepoint,
                timeformat  = '%-Hh')
    viewport.exit() # close pdf

##############################################################################

def make_samplepoint_plots(plot_samplepoints):
    if (not plot_samplepoints): return
    # delete old plots:
    list(map(os.remove, glob(outputdir+'/samplepoint_*.pdf')))
    # define verbosity to select species for plotting:
    verbosity = 2
    plotspecies = mecca.set_species(verbosity)
    print(HLINE, '\n')
    # define plotsamplepoints:
    if (plot_samplepoints == 1):
        plotsamplepoints = samplepointnames
    else:
        plotsamplepoints = ['0003', '0027'] # selected samplepoints
    for samplepoint in plotsamplepoints: # sample point loop
        print('Plotting sample point %-15s' % (samplepoint))
        caabaplot.xxxg(
            modelruns   = [[outputdir+'/'+skel+'/multirun/runs/'+samplepoint, skel]
                           for skel in all_skel],
            plotspecies = plotspecies,
            pdffile     = outputdir+'/samplepoint_'+samplepoint,
            pagetitle   = 'sample point: '+samplepoint,
            timeformat  = '%-Hh')
        print()

##############################################################################

def make_delta_skel_plots(plot_delta_skel):
    from viewport import viewport # from pycaaba
    from cycler import cycler
    if (not plot_delta_skel): return
    linecolors = ['k', 'r', 'g', 'b', 'm', 'y', 'c']
    print(HLINE, '\nPlotting delta_skel errors')
    viewport.init(4, 4, outputdir+'/delta_skel.pdf', 17, 8) # open pdf
    for samplepointnum, samplepoint in enumerate(samplepointnames): # sample point loop
        ax = viewport.next()
        ax.set_prop_cycle(cycler('color', linecolors))
        if (viewport.current == 1):
            # on new page, start with legend on a dummy plot:
            for targetnum, (target, abstol, reltol) in enumerate(targetdata): # target loop
                lines = plt.plot([0,0], linewidth=3, label='%s (abstol=%G, reltol=%G)' % (
                    target, abstol, reltol))
            plt.axis('off')
            legend = plt.legend(loc='center',
                                mode='expand',
                                fontsize = 'small',
                                title='delta_skel errors',
                                fancybox=True,
                                shadow=True,
                                borderaxespad=0.)
            plt.setp(legend.get_title(),fontsize='large')
            ax = viewport.next()
            ax.set_prop_cycle(cycler('color', linecolors))
        for targetnum, (target, abstol, reltol) in enumerate(targetdata): # target loop
            mydata = delta_skel_all[:,samplepointnum,targetnum]
            xval = np.arange(1, len(mydata)+1, 1)
            plt.plot(xval, mydata[:], '*', linestyle='solid')
        plt.xlim(0,len(mydata)+1)
        plt.ylim(0.,1.)
        plt.title('sample point: '+samplepoint)
        plt.xlabel('skeletal mechanism number')
        plt.ylabel('delta_skel')
    viewport.exit() # close pdf

##############################################################################

def read_ctrl():
    global ctrlfile
    if len(sys.argv)>1:
        ctrlfile = 'ctrl/%s.py' % (
            os.path.splitext(os.path.basename(sys.argv[1]))[0]) # no dir, no suffix
        exec(compile(open(ctrlfile).read(), ctrlfile, 'exec'), globals()) # data in ctrl file become global variables
    else:
        ctrlfile = '(no ctrl file specified, using default values)'
    print('%s\nxskeleton.py %s\n%s\n' % (HLINE2, ctrlfile, HLINE2))
    if (samplepointfile):
        full_calc = True
        print('samplepointfile   = %s.nc' % (samplepointfile))
        print('eps0              = %10g (start value for epsilon_ep)' % (eps0))
        print('eps_increase      = %10g (factor for increasing epsilon_ep)' % (eps_increase))
    else:
        full_calc = False
    print('plot_delta_skel   = %10d (0=no plots, 1=all)' % (plot_delta_skel))
    print('plot_targets      = %10d (0=no plots, 1=all)' % (plot_targets))
    print('plot_samplepoints = %10d (0=no plots, 1=all, 2=some)' % (plot_samplepoints))
    return full_calc

##############################################################################

if __name__ == '__main__':

    LOGFILE = tee.stdout_start('xskeleton.log', append=False) # stdout
    full_calc = read_ctrl() # read and print run control variables
    eps = eps0
    if (full_calc):
        cleanup()
    print('%s\n***** Full mechanism *****\n%s\n' % (HLINE2, HLINE2))
    workdir = get_workdir(None)
    create_skeleton_rpl(None)
    caaba_multirun()
    samplepointnames, N_samplepoints = get_samplepointnames()
    calc_oic()
    del_rxn = [False] * len(StoichNum)
    all_rates = save_rates(del_rxn)
    skelnum = 1
    # list of all mechanisms (full mechanism means del_rxnlist[:]=True):
    del_rxnlist = [del_rxn]
    # list of all eps (full mech means eps=0):
    epslist = [0]
    N_var = N_var_full
    # delta_skel_all_list is a list of lists of lists containing the errors
    # of all current skeletal mechanisms compared to the full mechanism
    # for all sample points and for all targets:
    delta_skel_all_list = []
    while True: # loop over skeletal mechanisms for skelnum = 1, 2, 3, ...
        eps_last = eps
        workdir = get_workdir(skelnum)
        print('%s\n***** Skeletal mechanism %d *****\n%s\n' % (HLINE2, skelnum, HLINE2))
        # Increase epsilon_ep to include less reactions:
        eps_too_small = True
        while eps_too_small:
            eps *= eps_increase
            # Confirm that new mechanism contains less species:
            N_var_new = np.sum(oicdata['f0']>eps) # f0=oic values, f1=species
            if (N_var_new < N_var):
                print('Number of species reduced from %d to %d' % (
                    N_var, N_var_new))
                N_var = N_var_new
                eps_too_small = False
            else:
                print('Still %d species in mechanism with eps = %g.' % (
                    N_var, eps))
        if (N_var_new == 0):
            sys.exit('ERROR: N_var=0')
        # Create skeletal mechanism excluding species with OIC < eps:
        create_skeletal_mechanism(eps)
        caaba_multirun()
        # maybe calc_oic() here again?
        delta_skel, delta_samplepoints = analyze_results()
        delta_skel_all_list.append(delta_samplepoints)
        # Exit loop when error is too big:
        if (delta_skel > 1.):
            break
        skelnum += 1
    print('\n%s\nEnd of skeletal mechanism testing loop\n%s' % (HLINE2, HLINE2))
    eps = eps_last
    # convert list of lists of lists to 3D numpy array:
    delta_skel_all = np.asarray(delta_skel_all_list)
    list_species(oicdata, epslist)
    list_reactions(del_rxnlist)
    # define directories of full and all skeletal mechanisms:
    all_skel = ['fullmech']
    for skeldir in sorted(glob(outputdir+'/skeleton_*')):
        all_skel.append(os.path.basename(skeldir))
    # if there are too many skeletal mechanisms, show only the
    # last 4 and the full mechanism:
    if (len(all_skel)>5):
        all_skel = [all_skel[i] for i in [0]+list(range(-4,0))]        
    make_target_plots(plot_targets)
    make_samplepoint_plots(plot_samplepoints)
    make_delta_skel_plots(plot_delta_skel)
    finalize(delta_skel_all.shape)

    ##############################################################################
