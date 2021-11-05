#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-03 12:08:25 sander>

# cams plots for CAABA-4.0 paper
# Rolf Sander, 2018

import os, sys, shutil
# define caabadir using $PWD and abspath (not realpath) because
# testsuite/ is a symlink in the caaba/ directory:
caabadir = os.path.abspath(os.getenv('PWD')+'/../..')
sys.path.append(os.path.realpath(caabadir+'/pycaaba'))

#from caabaplot import caabaplot
from viewport import viewport
import sys
import _mecca_spc # created automatically by spc2mpl
from netCDF4 import Dataset

# select model runs:
modelruns = []
modelruns.append(['mom',        'mom'])           # black
modelruns.append(['cb05bascoe', 'cb05bascoe'])    # red
modelruns.append(['mozart',     'mozart'])        # green
modelruns.append(['mim1',       'mim1'])          # blue
modelruns.append(['mcm',        'mcm'])           # magenta
modelruns.append(['jam',        'jam'])           # cyan
timeformat = ''
#timeformat = '%-H'
#tmin=144
#tmax=217

pdffile = 'cams.pdf'

##############################################################################

def plot_0d(modelruns, species, plottitle,
            ncfilename, timeformat='', tmin=0, tmax=0):
    from cycler import cycler
    from netCDF4 import num2date
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.dates import AutoDateFormatter, AutoDateLocator, DateFormatter
    linecolors = ['k', 'r', 'g', 'b', 'm', 'c', 'y']
    ax = viewport.next()
    ax.set_prop_cycle(cycler('color', linecolors))
    baserun = modelruns[0][1]
    # loop over all model runs:
    if (tmax>0):
        fulltrange = False # plot only time[tmin:tmax]
    else:
        fulltrange = True # plot full time range
        tmin = 0
    for (modelrundir, modelrunname) in modelruns:
        ncfullfilename = modelrundir+'/'+ncfilename
        ncid = Dataset(ncfullfilename)
        # define time:
        time = ncid.variables['time']
        if (fulltrange):
            tmax = len(time) # plot the whole model run
        t = num2date(time[tmin:tmax],time.units)
        # plot data:
        #----------------------
        if (species=='LTERP'):
            plottitle = 'terpenes'
            if (modelrunname=='mom'):
                plotdata = ncid.variables['APINENE'][tmin:tmax,0,0,0] + \
                           ncid.variables['BPINENE'][tmin:tmax,0,0,0] + \
                           ncid.variables['CAMPHENE'][tmin:tmax,0,0,0] + \
                           ncid.variables['CARENE'][tmin:tmax,0,0,0] + \
                           ncid.variables['SABINENE'][tmin:tmax,0,0,0]
            elif ((modelrunname=='mcm')|(modelrunname=='jam')):
                plotdata = ncid.variables['APINENE'][tmin:tmax,0,0,0] + \
                           ncid.variables['BPINENE'][tmin:tmax,0,0,0]
            elif (modelrunname=='mim1'):
                plotdata = 0*time[tmin:tmax] # zero dummy data (no terpenes in mim1)
            else:
                mydata = ncid.variables[species]
                plotdata = mydata[tmin:tmax,0,0,0]
            lines = plt.plot(t, plotdata, label=modelrunname)
            plt.ylabel('mol/mol')
        else:
            mydata = ncid.variables[species]
            plotdata = mydata[tmin:tmax,0,0,0]
            lines = plt.plot(t, plotdata, label=modelrunname)
            plt.ylabel(mydata.units)
        if (modelrunname != baserun):
            plt.setp(lines, linestyle='dotted', linewidth=2)
        #----------------------
        ncid.close()
    plt.title(plottitle)
    plt.xlabel('date')
    ax.grid(True) # hoizontal and vertical gridlines in plot
    # x-axis:
    ax.xaxis_date() # x-axis is a date
    xtick_locator = AutoDateLocator()
    # define locations of ticks on time axis:
    ax.xaxis.set_major_locator(xtick_locator) # automatic
    #ax.xaxis.set_major_locator(plt.MaxNLocator(5)) # max number of tick intervals
    # define format of ticks on time axis:
    if (timeformat):
        xformatter = DateFormatter(timeformat)
    else:
        xformatter = ticker.FuncFormatter(viewport.timeformat) # function
    #xformatter = AutoDateFormatter(xtick_locator)
    ax.xaxis.set_major_formatter(xformatter)
    # y-axis:
    # adjust yrange:
    # print 'yrange:', plt.ylim()[0], plt.ylim()[1]
    # plt.ylim(plt.ylim()[0] * 0.9, plt.ylim()[1] * 1.1)
    # show all digits instead of using offset for y-axis:
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
    yformatter = ticker.FuncFormatter(
        lambda x, p: viewport.scientificNotation(x))
    ax.yaxis.set_major_formatter(yformatter)
    ax.yaxis.set_major_locator(plt.MaxNLocator(5)) # max number of y-tick intervals
    ax.xaxis.set_major_locator(plt.MaxNLocator(5)) # max number of x-tick intervals

##############################################################################

# plot all species:
plotspecies = []
ncid = Dataset(modelruns[0][0]+'/caaba_mecca.nc')
for var in sorted(ncid.variables):
    if (ncid.variables[var].ndim==4): # exclude lon, lat, lev, time
        plotspecies.append(var)
ncid.close()
# plot selected species:
plotspecies = ['C5H8', 'LTERP', 'PAN', 'O3', 'OH', 'NO2']

viewport.init(3, 2, pdffile, 16,6) # open pdf
viewport.newpage()
print('Plotting these species:')
spc_names = _mecca_spc.spc_names() # load dictionary
for species in plotspecies: # species loop
    print('%s' % (species), end=' ') ; sys.stdout.flush()
    plottitle = r'$\sf ' + spc_names[species] + r'$'
    plot_0d(modelruns, species, plottitle, 'caaba_mecca.nc')
viewport.exit() # close pdf
print('\nCreated the plotfile:\n  %s' % (pdffile))

##############################################################################
