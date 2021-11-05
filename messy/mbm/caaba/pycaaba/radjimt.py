#!/usr/bin/env python3
# -*- coding: utf-8 -*- Time-stamp: <2019-06-03 17:26:46 sander>

# radjimt plots for CAABA-4.0 paper
# Rolf Sander, 2017

from caabaplot import caabaplot
from viewport import viewport
import sys
import _mecca_spc # created automatically by spc2mpl
from netCDF4 import Dataset

# select model runs:
modelruns = []
modelruns.append(['../testsuite/radjimt-00', 'equator'])
modelruns.append(['../testsuite/radjimt-50', '50 degree N'])
# plotspecies = ['em', 'O1D', 'O3P', 'O2', 'O3', 'Op', 'O2p', 'H', 'H2',
#                'OH', 'HO2', 'H2O', 'H2O2', 'NOx', 'NOy', 'N', 'N2', 'N2O', 'NO', 'NO2',
#                'NO3', 'N2O5', 'HNO3', 'Np', 'N2p', 'NOp', 'CO', 'CO2',
#                'HCHO', 'CH3', 'CH3O', 'CH4', 'Cl', 'ClO']
pagetitle = '' # 'radjimt' # no plot legend if title is empty
#timeformat = ''
timeformat = '%-H'
tmin=144
tmax=217

# plot selected species:
pdffile = '../testsuite/radjimt.pdf'
plotspecies = ['em', 'O1D', 'O3', 'Op', 'N', 'NOp']
viewport.init(2, 3, pdffile, 8,8) # open pdf
viewport.newpage()
print('Plotting these species:')
spc_names = _mecca_spc.spc_names() # load dictionary
for species in plotspecies: # species loop
    print('%s' % (species), end=' ') ; sys.stdout.flush()
    plottitle = r'$\sf ' + spc_names[species] + r'$'
    caabaplot.plot_0d(modelruns, species, pagetitle, plottitle,
                      'caaba_mecca.nc', timeformat, tmin=tmin, tmax=tmax)
viewport.exit() # close pdf
print('\nCreated the plotfile:\n  %s' % (pdffile))

# plot all species:
pdffile = '../testsuite/radjimt_mixrat.pdf'
plotspecies = []
ncid = Dataset(modelruns[0][0]+'/caaba_mecca.nc')
for var in sorted(ncid.variables):
    if (ncid.variables[var].ndim==4): # exclude lon, lat, lev, time
        plotspecies.append(var)
viewport.init(4, 4, pdffile, 17, 8) # open pdf
viewport.newpage()
print('Plotting these species:')
spc_names = _mecca_spc.spc_names() # load dictionary
for species in plotspecies: # species loop
    print('%s' % (species), end=' ') ; sys.stdout.flush()
    plottitle = r'$\sf ' + spc_names[species] + r'$'
    caabaplot.plot_0d(modelruns, species, pagetitle, plottitle,
                      'caaba_mecca.nc', timeformat, tmin=tmin, tmax=tmax)
viewport.exit() # close pdf
print('\nCreated the plotfile:\n  %s' % (pdffile))

# plot all J-values:
pdffile = '../testsuite/radjimt_jvalues.pdf'
plotspecies = []
ncid = Dataset(modelruns[0][0]+'/caaba_radjimt.nc')
for var in sorted(ncid.variables):
    if (ncid.variables[var].ndim==4): # exclude lon, lat, lev, time
        plotspecies.append(var)
viewport.init(4, 4, pdffile, 17, 8) # open pdf
viewport.newpage()
print('Plotting these J-values:')
for species in plotspecies: # species loop
    print('%s' % (species), end=' ') ; sys.stdout.flush()
    plottitle = species
    caabaplot.plot_0d(modelruns, species, pagetitle, plottitle,
                'caaba_radjimt.nc', timeformat, tmin=tmin, tmax=tmax)
viewport.exit() # close pdf
print('\nCreated the plotfile:\n  %s' % (pdffile))
