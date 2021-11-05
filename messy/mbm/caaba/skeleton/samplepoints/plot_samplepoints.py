#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-07 14:04:52 sander>

# xskeleton: execute mechanism reduction to obtain a skeletal mechanism
# Rolf Sander, 2016

##############################################################################

import os, sys
sys.path.append(os.path.abspath('../../pycaaba'))
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from viewport import viewport # from pycaaba
import numpy as np

# select a file with the sample points (w/o suffix '.nc'):
#samplepointfilename = 'skeleton_samplepoints_small'
samplepointfilename = 'skeleton_samplepoints_30'
#samplepointfilename = 'skeleton_samplepoints_30_loworg'

##############################################################################

def scientificNotation(value):
    if value == 0:
        return '0'
    else:
        e = np.log10(np.abs(value))
        m = np.sign(value) * 10 ** (e - int(e))
        # formatstring = r'${:.0f} \cdot 10^{{{:d}}}$'.format(m, int(e))
        formatstring = '%gE%d' % (m, int(e))
        #print formatstring
        return formatstring

def makelinplot(species, mydata):
    xval = np.arange(1, len(mydata)+1, 1)
    ax = viewport.next()
    plt.plot(xval, mydata[:], '*', linestyle='dotted', color='r', label=species)
    plt.xlim(0,len(mydata)+1)
    plt.title(species)
    plt.xlabel('sample point number')
    plt.ylabel('mol/mol')
    # format y labels:
    # show all digits instead of using offset for y-axis:
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
    formatter = ticker.FuncFormatter(lambda x, p: scientificNotation(x))
    ax.yaxis.set_major_formatter(formatter)
    # max number of tick intervals:
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    #print plt.ylim() # data range

def makelogplot(species, mydata):
    xval = np.arange(1, len(mydata)+1, 1)
    ax = viewport.next()
    plt.plot(xval, mydata[:], '*', linestyle='dotted', color='g', label=species)
    plt.xlim(0,len(mydata)+1)
    plt.title(species)
    plt.xlabel('sample point number')
    plt.ylabel('mol/mol')
    ax.set_yscale('log')
    #ax.yaxis.set_major_locator(ticker.LogLocator(base = 1000.0))
    ax.yaxis.set_major_locator(ticker.LogLocator(numticks=6))
    
def makeplots():
    viewport.init(4, 4, samplepointfilename+'.pdf', 17, 8) # open pdf
    ncid = Dataset(samplepointfilename+'.nc')
    # species loop:
    #for species in ['HCHO', 'NO', 'OH', 'HO2', 'O3']: # small list, only for testing 
    for species in sorted(ncid.variables):
        data = ncid.variables[species]
        print('%25s MIN: %e MAX: %e' % (species, min(data), max(data)))
        makelinplot(species+' (LIN)', data)
        if (max(data)>min(data)):
            makelogplot(species+' (LOG)', (data))
        else:
            ax = viewport.next()
            plt.plot([0,1]) # dummy plot
    viewport.exit() # close pdf
    ncid.close()

##############################################################################

if __name__ == '__main__':

  makeplots()

##############################################################################
