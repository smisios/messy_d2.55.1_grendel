#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2020-09-15 17:27:29 sander>

# jvpp_plot_step1: plot results from step 1 of JVPP
# Rolf Sander, 2020

##############################################################################

import os, sys, shutil
assert sys.version_info >= (3, 6)
CAABADIR = os.path.realpath(os.path.dirname(__file__)+'/../..')
sys.path.insert(1, os.path.realpath(CAABADIR+'/pycaaba'))
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from viewport import viewport
from rstools import corename

##############################################################################

def maketwoplots(species):

    def makeplot(scale):
        ax = viewport.next()
        if (scale=='log'):
            plt.yscale('log')
            plt.ylabel('lg('+ylabel+')')
            mycolor = 'orange'
            title_suffix = '(log-scale)'
        if (scale=='lin'):
            plt.ylabel(ylabel)
            mycolor = 'red'
            title_suffix = '(linear)'
        # linestyles:
        # https://matplotlib.org/3.3.1/gallery/lines_bars_and_markers/linestyles.html
        # markers:
        # https://matplotlib.org/3.3.1/api/markers_api.html
        # plot literature values:
        plt.plot(data_lit[:,0], data_lit[:,1],
                 marker='o', linestyle=(0,(1,5)), color=mycolor, label='mylabel', linewidth=1)
        # plot values regridded to 176 bins (only the first 142 values are used):
        plt.plot(wavedata[0:141], data_176[0:141],
                 marker='+', linestyle='solid', color='black', label='mylabel', linewidth=1)
        plt.title(f'{title_prefix} {species} {title_suffix}')
        plt.xlabel('wavelength (nm)')

    if (ext=='sig'):
        title_prefix = 'sigma'
        ylabel = 'cross section (cm2)'
    if (ext=='phi'):
        title_prefix = 'phi'
        ylabel = 'quantum yield'        
    # read literature values:
    data_lit = np.genfromtxt('dat_lit/spectra/'+species+'.'+ext)
    # read values regridded to 176 bins:
    data_176 = np.genfromtxt('workdir_176/'+species+'.'+ext+'_176')
    makeplot('lin')
    makeplot('log')

##############################################################################

if __name__ == '__main__':

    viewport.init(2, 2, 'jvpp_plot_step1.pdf',17,8)

    # Load data from wave.dat into a numpy structured array:
    # http://docs.scipy.org/doc/numpy/user/basics.rec.html
    wavedata = np.genfromtxt('workdir_176/wave.dat')

    for ext in ['sig','phi']:
        for species in sorted(glob('workdir_176/*.'+ext+'_176')):
            species = corename(species)
            print(species)
            maketwoplots(species)

    viewport.exit()

##############################################################################
